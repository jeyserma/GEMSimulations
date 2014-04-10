/*
 * primaryIonization.c
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
#include <TTree.h>
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

using namespace Garfield;
const std::string GARFIELD = "/afs/cern.ch/user/j/jeyserma/garfieldpp/";

struct cluster {
    
    Double_t x0, y0, z0, t0, ec;
    Int_t nc;
    std::vector<Double_t> x1, y1, z1, t1, e1, vx1, vy1, vz1;
};

struct test {
    
    Double_t x0;
};


struct particle {
    
    //TString type;
    std::vector<test> clusters;
    Double_t x0, y0, z0, vx0, vy0, vz0, e0, t0;
    Int_t noClusters;
    Double_t lambda, stoppingPower;
    
};


struct avalancheI {
    
    Double_t x0, y0, z0, t0;
    Int_t ni;
    std::vector<Double_t> x1, y1, z1, x2, y2, z2, t1, e1, t2, e2;
    std::vector<Int_t> status;
};

TString int2str(Int_t t) {
    TString str; 
    str.Form("%d", t);
    return str;
} 

int main(int argc, char * argv[]) {

	TApplication app("app", &argc, argv);
 	gRandom = new TRandom3(0); // set random seed
    
    std::string workingdir = "includes/";
    workingdir.append(argv[1]);
    workingdir.append("/");

    bool debug = true;
    const int it = 5;

	// Load GEM dimensions and list files
	FILE *fp;
    std::string geo = workingdir + "geometry.txt";
    std::cout << geo << std::endl;
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
	std::string efile = workingdir + "ELIST.lis";
	std::string nfile = workingdir + "NLIST.lis";
	std::string mfile = workingdir + "MPLIST.lis";
	std::string sfile = workingdir + "PRNSOL.lis";
	if(!fm->Initialise(efile, nfile, mfile, sfile, "mm")) {
		std::cout << "Error while loading the ANSYS field map files." << std::endl;
	}
	fm->EnableMirrorPeriodicityX();
	fm->EnableMirrorPeriodicityY();
    if(debug) {
        fm->PrintRange();
    }
	
	// Gas setup
	MediumMagboltz* gas = new MediumMagboltz();
    gas->LoadGasFile("includes/config1/Ar-CO2-CF4-40-25-35.gas");
    //gas->SetMaxElectronEnergy(200.);
  	gas->EnableDebugging();
	gas->Initialise();
	gas->DisableDebugging();
	gas->LoadIonMobility(GARFIELD + "Data/IonMobility_Ar+_Ar.txt");  
     
	//Associate the gas with the corresponding field map material.
 	const int nMaterials = fm->GetNumberOfMaterials();
 	for(int i=0; i<nMaterials; ++i) {
 	
   		const double eps = fm->GetPermittivity(i);
   		if(fabs(eps - 1.) < 1.e-3) fm->SetMedium(i, gas);
 	}
    if(debug) {
        fm->PrintMaterials();
    }
    
	// Sensor setup
	Sensor* sensor = new Sensor();
	sensor->AddComponent(fm);
	sensor->SetArea(-5*pitch, -5*pitch, 0.0, 5*pitch, 5*pitch, totalth);
	
    // Setup HEED
    TrackHeed* heed = new TrackHeed();
    heed->SetSensor(sensor);
    //heed->DisableDeltaElectronTransport();
    heed->SetParticle("mu");
    heed->SetMomentum(100.e9); // 100 GeV
    if(debug) {
        heed->EnableDebugging();
    }

	// Calculate avalanche
	int nc;
	double x0, y0, z0, t0, e0, x1, y1, z1, t1, e1, x2, y2, z2, t2, e2;
	double vx0, vy0, vz0, vx2, vy2, vz2;
    double ec, extra;
	
    TRandom3 r;
    r.SetSeed(0);
    
    // Prepare trees and saving file
    TString savefile = workingdir + "output/primaryIonization.root";
    TFile *file = new TFile(savefile, "RECREATE");
    
    // Prepare tree for electrons
    particle p;
    TTree *pTree = new TTree("pTree", "Charged particle");
    pTree->Branch("test", &p.clusters);
    pTree->Branch("x0", &p.x0);
    pTree->Branch("y0", &p.y0);
    pTree->Branch("z0", &p.z0);
    pTree->Branch("vx0", &p.vx0);
    pTree->Branch("vy0", &p.vy0);
    pTree->Branch("vz0", &p.vz0);
    pTree->Branch("t0", &p.t0);
    pTree->Branch("e0", &p.e0);
    //pTree->Branch("type", &p.type);
    pTree->Branch("noClusters", &p.noClusters);
    pTree->Branch("stoppingpower", &p.stoppingPower);
    pTree->Branch("lambda", &p.lambda);
    

    
    // Start iteration
	for(int i=0; i<it; i++) {
        
        if(debug) {
            std::cout << "Progress: " << 100.*(i+1)/it << "%" << std::endl;
        }
        
        // Set random velocity direction, in positive hemisphere
        const double ctheta = 1. - r.Uniform();
        const double stheta = sqrt(1. - ctheta*ctheta);
        const double phi = TwoPi*r.Uniform();
        vx0 = cos(phi)*stheta;
        vy0 = sin(phi)*stheta;
        vz0 = ctheta;
  
        // Set initial time
        t0 = 0.0;

        // Generate random (x,y) position in unit cell
        x0 = r.Uniform()*pitch/2; 
        y0 = r.Uniform()*pitch*TMath::Sqrt(3)/2;
        z0 = 0.;
        
        p.x0 = x0;
        p.y0 = y0;
        p.z0 = z0;
        p.vx0 = vx0;
        p.vy0 = vy0;
        p.vz0 = vz0;
        p.t0 = t0;
        p.e0 = 100.e9   ;
        //p.type = "muon";
        p.noClusters = 0;
        p.lambda = 1/heed->GetClusterDensity();
        p.stoppingPower = heed->GetStoppingPower();
        
    
        heed->NewTrack(x0, y0, z0, t0, vx0, vy0, vz0);
        
        // Loop over clusters
        while(heed->GetCluster(x1, y1, z1, t1, nc, ec, extra)) {
            
            p.noClusters++;
            
            test tmpCl;
            tmpCl.x0 = x1;
            std::cout << p.clusters.size() << std::endl;
           /*
            tmpCl.y0 = y1;
            tmpCl.z0 = z1;
            tmpCl.nc = nc;
            tmpCl.ec = ec;
         
            // Skip the clusters which are not in the drift region
            if(z1 > drift) {
                continue;
            }

            for(int j=0; j<nc; j++) {

                heed->GetElectron(j, x2, y2, z2, t2, e2, vx2, vy2, vz2);
                    
                tmpCl.x1.push_back(x2);
                tmpCl.y1.push_back(y2);
                tmpCl.z1.push_back(z2);
                tmpCl.t1.push_back(t2);
                tmpCl.e1.push_back(e2);
                tmpCl.vx1.push_back(vx2);
                tmpCl.vy1.push_back(vy2);
                tmpCl.vz1.push_back(vz2);
            }
            * */
    
            p.clusters.push_back(tmpCl);
           
        }

        pTree->Fill();
        p.clusters.clear(); 
	}

    pTree->Write();
    //ITree->Write();
    file->Close();
    
    
    /*
	TCanvas* c1 = new TCanvas("c1", "Signal");
    ViewSignal* signalView = new ViewSignal();
    TString filename;
	signalView->SetSensor(sensor);
	signalView->SetCanvas(c1);

    filename = workingdir + "/output/" + argv[2] + "drift/timeEl.root";
    signalView->PlotSignal("readout");
	c1->SaveAs(filename);
    c1->Clear();
    
    filename = workingdir + "/output/" + argv[2] + "drift/timeIon.root";
	signalView->PlotSignal("ions");
	c1->SaveAs(filename);
    c1->Clear();
    
    
    
    sensor->SetTransferFunction(transferf);
    sensor->ConvoluteSignal();
    
    tau = 1.;
    sensor->SetTransferFunction(transferf);
    sensor->ConvoluteSignal();
    
    filename = workingdir + "/output/" + argv[2] + "drift/timeElTransfer1.root";
    signalView->PlotSignal("readout");
	c1->SaveAs(filename);
    c1->Clear();
    
    filename = workingdir + "/output/" + argv[2] + "drift/timeIonTransfer1.root";
	signalView->PlotSignal("ions");
	c1->SaveAs(filename);
    c1->Clear();
    
    tau = 5.;
    sensor->SetTransferFunction(transferf);
    sensor->ConvoluteSignal();
    
    filename = workingdir + "/output/" + argv[2] + "drift/timeElTransfer5.root";
    signalView->PlotSignal("readout");
	c1->SaveAs(filename);
    c1->Clear();
    
    filename = workingdir + "/output/" + argv[2] + "drift/timeIonTransfer5.root";
	signalView->PlotSignal("ions");
	c1->SaveAs(filename);
    c1->Clear();
    
    tau = 10.;
    sensor->SetTransferFunction(transferf);
    sensor->ConvoluteSignal();
    
    filename = workingdir + "/output/" + argv[2] + "drift/timeElTransfer10.root";
    signalView->PlotSignal("readout");
	c1->SaveAs(filename);
    c1->Clear();
    
    filename = workingdir + "/output/" + argv[2] + "drift/timeIonTransfer10.root";
	signalView->PlotSignal("ions");
	c1->SaveAs(filename);
    c1->Clear();
    
    tau = 25.;
    sensor->SetTransferFunction(transferf);
    sensor->ConvoluteSignal();
    
    filename = workingdir + "/output/" + argv[2] + "drift/timeElTransfer25.root";
    signalView->PlotSignal("readout");
	c1->SaveAs(filename);
    c1->Clear();
    
    filename = workingdir + "/output/" + argv[2] + "drift/timeIonTransfer25.root";
	signalView->PlotSignal("ions");
	c1->SaveAs(filename);
    c1->Clear();

     * */
	return 0;
	app.Run(); 
}


