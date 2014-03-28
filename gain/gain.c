/*
 * gain.c
 * Author: Jan Eysermans 
 * 2014
 */

// General libs
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdio.h>     
#include <stdlib.h> 
#include <math.h> 

// ROOT libs
#include <TROOT.h>
#include <TFile.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TMath.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoPcon.h>
#include <TGeoHalfSpace.h>
#include <TGeoMatrix.h>
#include <TGeoCompositeShape.h>
#include <TLegend.h>
#include <TBenchmark.h>

// Garfield++ libs
#include "ComponentAnsys123.hh"
#include "ViewField.hh"
#include "MediumMagboltz.hh"
#include "Sensor.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "Random.hh"
#include "Plotting.hh"
#include "DriftLineRKF.hh"
#include "TrackHeed.hh"
#include "ComponentAnalyticField.hh"
#include "ViewSignal.hh" 
#include "ViewFEMesh.hh"
#include "ViewFEMesh.hh"

using namespace Garfield;
const std::string GARFIELD = "/afs/cern.ch/user/j/jeyserma/garfieldpp/";

TString int2str(Int_t t) {
    TString str; 
    str.Form("%d", t);
    return str;
} 

int main(int argc, char * argv[]) {
    
    TStopwatch watch;

	TApplication app("app", &argc, argv);
 	plottingEngine.SetDefaultStyle();
 	gRandom = new TRandom3(0); // set random seed
 	int i, j, k, l; // loop constants
    
    const bool signal = false;
    const bool ionSignal = false;
    const bool debug = false;
    const int it = 1000;

	// Plot options
	gROOT->ProcessLine(".x rootstyle.c");
	gROOT->SetStyle("newStyle");
	
	// Working directory settings
	std::string dir = "includes/";
	dir.append(argv[1]);
	std::string geo = dir + "/geometry.txt";

	// Load GEM dimensions and list files
	FILE *fp;
	fp = fopen(geo.c_str(), "r");
    char str1[100], str2[100];
    
    fscanf(fp, "%s %s\n", str1, str2);
    const int layers = (int)atof(str2);
    fscanf(fp, "%s %s\n", str1, str2);
    const double drift = atof(str2);
    fscanf(fp, "%s %s\n", str1, str2);
    const double transfer = atof(str2);
    fscanf(fp, "%s %s\n", str1, str2);
    const double induct = atof(str2);
    fscanf(fp, "%s %s\n", str1, str2);
    const double kapton = atof(str2);
    fscanf(fp, "%s %s\n", str1, str2);
    const double metal = atof(str2);
    fscanf(fp, "%s %s\n", str1, str2);
    const double pitch = atof(str2);
    fscanf(fp, "%s %s\n", str1, str2);
    const double outdia = atof(str2);
    fscanf(fp, "%s %s\n", str1, str2);
    const double middia = atof(str2);
    fscanf(fp, "%s %s\n", str1, str2);
    const double rim = atof(str2);
   
    fclose(fp);
    
    const double layerth = kapton + 2*metal; // thickness one GEM layer
    const double totalth = (layers-1)*transfer + induct + drift + layers*layerth; // total GEM thickness
    
	// Load the field map
	ComponentAnsys123* fm = new ComponentAnsys123();
	std::string efile = dir + "/ELIST.lis";
	std::string nfile = dir + "/NLIST.lis";
	std::string mfile = dir + "/MPLIST.lis";
	std::string sfile = dir + "/PRNSOL.lis";
	std::string wfile = dir + "/WSOL.lis";
	if(!fm->Initialise(efile, nfile, mfile, sfile, "mm")) {
		std::cout << "Error while loading the ANSYS field map files." << std::endl;
	}
	fm->EnableMirrorPeriodicityX();
	fm->EnableMirrorPeriodicityY();
    fm->PrintRange();
 	fm->SetWeightingField(wfile, "readout"); //Loading the weighting field
	
	// Gas setup
	MediumMagboltz* gas = new MediumMagboltz();
  	gas->SetComposition("ar", 75., "co2", 25.);
  	gas->SetTemperature(293.15);
  	gas->SetPressure(760.0);	
    gas->SetMaxElectronEnergy(200.);
  	gas->EnableDebugging();
	gas->Initialise();
	gas->DisableDebugging();
    const double rPenning = 0.57;
    const double lambdaPenning = 0.;
    gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
	gas->LoadIonMobility(GARFIELD + "Data/IonMobility_Ar+_Ar.txt");
     
	 //Associate the gas with the corresponding field map material.
 	const int nMaterials = fm->GetNumberOfMaterials();
 	for(i=0; i<nMaterials; ++i) {
 	
   		const double eps = fm->GetPermittivity(i);
   		if(fabs(eps - 1.) < 1.e-3) fm->SetMedium(i, gas);
 	}
    fm->PrintMaterials();
    
	// Sensor setup
	Sensor* sensor = new Sensor();
	sensor->AddComponent(fm);
	sensor->SetArea(-5*pitch, -5*pitch, 0.0, 5*pitch, 5*pitch, totalth); // Volume of the sensor
	
	// Setup electron transport
	AvalancheMicroscopic* aval = new AvalancheMicroscopic();
	aval->SetSensor(sensor);
	//aval->EnableAvalancheSizeLimit(1000);
    
	if(signal) {
		sensor->AddElectrode(fm, "readout");
 		sensor->SetTimeWindow(0., 1, 100); // min, bin, max
		aval->EnableSignalCalculation();
	}

	// Setup ion transport
    AvalancheMC* iondrift = new AvalancheMC(); // Object must be initialized to avoid errors
    if(ionSignal) {
          
        iondrift->SetSensor(sensor);
        //drift->SetTimeSteps(0.05); // default: 10 ps
        iondrift->SetDistanceSteps(2e-4); // default: 10 um
        //drift->SetCollisionSteps(100); // default: 100
        if(signal) {
            iondrift->EnableSignalCalculation();
        }
    }

	// Calculate avalanche
	int ne, ni, np;
	double x0, y0, z0, x1, y1, z1, x2, y2, z2, xi1, yi1, zi1, xi2, yi2, zi2, t0, t1, t2, ti1, ti2, e0, e1, e2;
	double vx0, vy0, vz0; // initial velocity coordinates (if null vector, a random direction is generated)
	
	int status;
	int topMetalEl[layers], topMetalIon[layers]; // # electrons arrived on top metal electrode for each GEM layer
	int bottomMetalEl[layers], bottomMetalIon[layers]; // # electrons arrived on bottom metal electrode for each GEM layer
	int kaptonEl[layers], kaptonIon[layers]; // # electrons in kapton for each GEM layer
	int signalEl, signalIon; // electrons in induction region --> creation of the signal!
	int otherEl, otherIon; // lost electrons
    
    int totalEl = 0, totalIon = 0;
    int totalTopMetalEl[layers], totalTopMetalIon[layers];
    int totalBottomMetalEl[layers], totalBottomMetalIon[layers];
    int totalKaptonEl[layers], totalKaptonIon[layers];
    int totalSignalEl = 0, totalSignalIon = 0;
    int totalOtherEl = 0, totalOtherIon = 0;
    
    // All arrays to zero;
    memset(totalTopMetalEl, 0, sizeof(totalTopMetalEl));
    memset(totalBottomMetalEl, 0, sizeof(totalBottomMetalEl));
    memset(totalKaptonEl, 0, sizeof(totalKaptonEl));
	memset(totalTopMetalIon, 0, sizeof(totalTopMetalIon));
    memset(totalBottomMetalIon, 0, sizeof(totalBottomMetalIon));
    memset(totalKaptonIon, 0, sizeof(totalKaptonIon));
    
    // Setup histograms
    TH1F *hSignalEl = new TH1F("signalEl", "Electrons in induction region", 150, 0.0, 1500.);
    TH1F *hOtherEl = new TH1F("otherEl", "Other electrons", 100, 0.0, 100.);
    TH1F *hSignalIon = new TH1F("signalIon", "Ions in drift region", 150, 0.0, 1500.);
    TH1F *hOtherIon = new TH1F("otherIon", "Other ions", 100, 0.0, 100.);

    TH1F *hTopMetalEl[layers];
    TH1F *hBottomMetalEl[layers];
    TH1F *hKaptonEl[layers];
    TH1F *hTopMetalIon[layers];
    TH1F *hBottomMetalIon[layers];
    TH1F *hKaptonIon[layers];

    TString htitle, hname;
    for(l=0; l<layers; l++) {
    
        htitle = "Top metal electrons, layer " + int2str(l+1);
        hname = "tmetalEl-" + int2str(l+1);
        hTopMetalEl[l] = new TH1F(hname, htitle, 100, 0.0, 100.);
        
        htitle = "Bottom metal electrons, layer " + int2str(l+1);
        hname = "bmetalEl-" + int2str(l+1);
        hBottomMetalEl[l] = new TH1F(hname, htitle, 100, 0.0, 100.);
        
        htitle = "Kapton electrons, layer " + int2str(l+1);
        hname = "kaptonEl-" + int2str(l+1);
        hKaptonEl[l] = new TH1F(hname, htitle, 100, 0.0, 100.);

        htitle = "Top metal ions, layer " + int2str(l+1);
        hname = "tmetalIon-" + int2str(l+1);
        hTopMetalIon[l] = new TH1F(hname, htitle, 100, 0.0, 100.);
        
        htitle = "Bottom metal ions, layer " + int2str(l+1);
        hname = "bmetalIon-" + int2str(l+1);
        hBottomMetalIon[l] = new TH1F(hname, htitle, 100, 0.0, 100.);
        
        htitle = "Kapton ions, layer " + int2str(l+1);
        hname = "kaptonIon-" + int2str(l+1);
        hKaptonIon[l] = new TH1F(hname, htitle, 100, 0.0, 100.);
    }

    // Set random velocity direction
    vx0 = 0.0;
    vy0 = 0.0;
    vz0 = 0.0;
    
    // Set initial energy and time
    e0 = 30.0;
    t0 = 0.0;
    
    TRandom3 r;
    r.SetSeed(0);
    
    /*
    FILE *progress;
    progress = fopen("progress.txt", "w");
    
    fprintf(progress, "---------------------------------------------------------------\n");
    fflush(progress);
    */

    // Start iteration
	for(i=0; i<it; i++) {
        
        std::cout << i << std::endl;
        
        /*
        fprintf(progress, "%f \n", 100.*(i+1)/it);
        fflush(progress);
        */
        
        // Generate random (x,y,z) position in unit cell
        x0 = r.Uniform()*pitch/2; 
        y0 = r.Uniform()*pitch*TMath::Sqrt(3)/2;
        z0 = r.Uniform()*drift;

        // Reset variables to zero
        memset(topMetalEl, 0, sizeof(topMetalEl));
        memset(bottomMetalEl, 0, sizeof(bottomMetalEl));
        memset(kaptonEl, 0, sizeof(kaptonEl));
        memset(topMetalIon, 0, sizeof(topMetalIon));
        memset(bottomMetalIon, 0, sizeof(bottomMetalIon));
        memset(kaptonIon, 0, sizeof(kaptonIon));
		signalEl = 0;
        otherEl = 0;
        signalIon = 0;
        otherIon = 0;
        
		aval->AvalancheElectron(x0, y0, z0, t0, e0, vx0, vy0, vz0);
		aval->GetAvalancheSize(ne, ni);
		np = aval->GetNumberOfElectronEndpoints();
        totalEl = totalEl + np;
        totalIon = totalIon + np;
     
        // Loop over all electrons in avalanche
        for(k=0; k<np; k++) {

			aval->GetElectronEndpoint(k, x1, y1, z1, t1, e1, x2, y2, z2, t2, e2, status);
            
            if(ionSignal) {
                
                iondrift->DriftIon(x1, y1, z1, t1);
                iondrift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, status);   
                  
                // Trace the ion
                if(zi2 < 0.05*drift) {
                    signalIon++;
                    totalSignalIon++;
                }	
                // calculate amount of ions in each GEM element
                for(l=0; l<layers; l++) {

                    // Upper metal
                    if(zi2 > drift + l*transfer + l*layerth && zi2 < drift + l*transfer + l*layerth + metal) {
                        topMetalIon[l]++;  
                        totalTopMetalIon[l]++;  
                        break;
                    }
                    // Lower metal
                    else if(zi2 > drift + l*transfer + (l+1)*layerth - metal && zi2 < drift + l*transfer + (l+1)*layerth) {
                        bottomMetalIon[l]++;
                        totalBottomMetalIon[l]++;
                        break;
                    }
                    // Kapton
                    else if(zi2 > drift + l*transfer + l*layerth + metal && zi2 < drift + l*transfer + (l+1)*layerth +  metal) {
                        kaptonIon[l]++;
                        totalKaptonIon[l]++;
                        break;
                    }
                }          
            }

            // See whether the electron has ended in the induct region
			if(z2 > totalth - 0.95*induct) {
				signalEl++;
                totalSignalEl++;
                continue;
			}	
            
            // calculate amount of electrons in each GEM element
			for(l=0; l<layers; l++) {
			
				// Upper metal
				if(z2 > drift + l*transfer + l*layerth && z2 < drift + l*transfer + l*layerth + metal) {
					topMetalEl[l]++;  
                    totalTopMetalEl[l]++;  
                    break;
				}
				// Lower metal
				else if(z2 > drift + l*transfer + (l+1)*layerth - metal && z2 < drift + l*transfer + (l+1)*layerth) {
					bottomMetalEl[l]++;
                    totalBottomMetalEl[l]++;
                    break;
				}
				// Kapton
				else if(z2 > drift + l*transfer + l*layerth + metal && z2 < drift + l*transfer + (l+1)*layerth +  metal) {
					kaptonEl[l]++;
                    totalKaptonEl[l]++;
                    break;
				}
			}
		}
        
        // Calculate total amount of other electrons
        otherEl = np - signalEl;
        otherIon = np - signalIon;
        for(l=0; l<layers; l++) {
            otherEl = otherEl - (int)topMetalEl[l] - (int)bottomMetalEl[l] - (int)kaptonEl[l];
            otherIon = otherIon - (int)topMetalIon[l] - (int)bottomMetalIon[l] - (int)kaptonIon[l];      
        }
        totalOtherEl = totalOtherEl + otherEl;
        totalOtherIon = totalOtherIon + otherIon;

        
        // Fill histograms  
        if(signalEl != 0) hSignalEl->Fill((double)signalEl);
        if(otherEl != 0) hOtherEl->Fill((double)otherEl);
        if(signalIon != 0) hSignalIon->Fill((double)signalIon);
        if(otherIon != 0) hOtherIon->Fill((double)otherIon);
        for(l=0; l<layers; l++) {

            if(topMetalEl[l] != 0) hTopMetalEl[l]->Fill((double)topMetalEl[l]);
            if(bottomMetalEl[l] != 0) hBottomMetalEl[l]->Fill((double)bottomMetalEl[l]);
            if(kaptonEl[l] != 0) hKaptonEl[l]->Fill((double)kaptonEl[l]);
            if(topMetalIon[l] != 0) hTopMetalIon[l]->Fill((double)topMetalIon[l]);
            if(bottomMetalIon[l] != 0) hBottomMetalIon[l]->Fill((double)bottomMetalIon[l]);
            if(kaptonIon[l] != 0) hKaptonIon[l]->Fill((double)kaptonIon[l]);
        }
    
        
        if(debug) {
            
            std::cout << "---------------------" << std::endl;
            std::cout << "Signal electrons " << signalEl << std::endl;
            std::cout << "Other electrons " << otherEl << std::endl;
            std::cout << "kapton electrons " << kaptonEl[0] << std::endl;
            std::cout << "top metal electrons " << topMetalEl[0] << std::endl;
            std::cout << "bottom metal electrons " << bottomMetalEl[0] << std::endl;
            std::cout << "---------------------" << std::endl;
            std::cout << "Signal ions " << signalIon << std::endl;
            std::cout << "Other ions " << otherIon << std::endl;
            std::cout << "kapton ions " << kaptonIon[0] << std::endl;
            std::cout << "top metal ions " << topMetalIon[0] << std::endl;
            std::cout << "bottom metal ions " << bottomMetalIon[0] << std::endl;     
        } 
	}
    //fclose(progress);
    
    if(debug) {
        std::cout << "Total produced electrons " << totalEl << std::endl;
        std::cout << "Total produced ions " << totalIon << std::endl;
    }
        
    // Open file to write statistics
    ofstream out;
    TString filename = dir + "/output/results.txt";
    out.open (filename);
    out << "---------------------------------------------------------------" << std::endl;
    out << "Total produced electrons \t" << totalEl << std::endl;
    out << "Total produced ions \t\t" << totalIon << std::endl;
    out << "Total signal electrons \t\t" << totalSignalEl << " (" << 100.*totalSignalEl/totalEl << "%)" << std::endl;
    out << "Total signal ions \t\t\t" << totalSignalIon << " (" << 100.*totalSignalIon/totalIon << "%)" << std::endl;
    out << "Mean signal electrons \t\t" << hSignalEl->GetMean() << " (RMS: " << hSignalEl->GetRMS() << ")" << std::endl;
    out << "Mean signal ions \t\t\t" << hSignalIon->GetMean() << " (RMS: " << hSignalIon->GetRMS() << ")" << std::endl;

    for(l=0; l<layers; l++) {
        
        out << "----------------------------------------------------------" << std::endl;
        out << "LAYER " << l+1 << std::endl;
        out << "Top metal electrons \t\t"   << totalTopMetalEl[l] << " (" << 100.*totalTopMetalEl[l]/totalEl << "%)" << std::endl;
        out << "Top metal ions \t\t\t\t"      << totalTopMetalIon[l] << " (" << 100.*totalTopMetalIon[l]/totalIon << "%)" << std::endl;
        out << "Bottom metal electrons\t\t"    << totalBottomMetalEl[l] << " (" << 100.*totalBottomMetalEl[l]/totalEl << "%)" << std::endl;
        out << "Bottom metal ions\t\t\t"    << totalBottomMetalIon[l] << " (" << 100.*totalBottomMetalIon[l]/totalIon << "%)" << std::endl;
        out << "Kapton electrons:  \t\t\t"        << totalKaptonEl[l] << " (" << 100.*totalKaptonEl[l]/totalEl << "%)" << std::endl;
        out << "Kapton ions: \t\t\t\t" << totalKaptonIon[l] << " (" << 100.*totalKaptonIon[l]/totalIon << "%)" << std::endl;
    }

    out.close();


	// Write the histogram
	filename = dir + "/output/results.root";
	TFile *f = new TFile(filename, "update");
	f->cd();
	hSignalEl->Write();
	hOtherEl->Write();
    hSignalIon->Write();
	hOtherIon->Write();
    for(l=0; l<layers; l++) {

        hTopMetalEl[l]->Write();
        hBottomMetalEl[l]->Write();
        hKaptonEl[l]->Write();
        hTopMetalIon[l]->Write();
        hBottomMetalIon[l]->Write();
        hKaptonIon[l]->Write();
    }
	f->Close();
    
	if(signal) {  
	
		TString filename = dir + "/output/time.root";
		TCanvas* c1 = new TCanvas("c1", "Signal");
		ViewSignal* signalView = new ViewSignal();
		signalView->SetSensor(sensor);
		signalView->SetCanvas(c1);
		signalView->PlotSignal("readout");
		c1->SaveAs(filename);
	}

    std::cout << "--------------- TIMING ---------------" << std::endl;
    std::cout << watch.CpuTime() << std::endl;
	return 0;
	app.Run(); 
}
