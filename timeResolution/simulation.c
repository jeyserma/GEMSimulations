/*
 * simulation.c
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
double tau = 1.;

struct avalancheE {
    
    Double_t x0, y0, z0, t0, e0, vx0, vy0, vz0;
    Int_t ne;
    std::vector<Double_t> x1, y1, z1, x2, y2, z2, t1, e1, t2, e2;
    std::vector<Int_t> status;
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

double transferf(double t);

int main(int argc, char * argv[]) {

	TApplication app("app", &argc, argv);
 	gRandom = new TRandom3(0); // set random seed
    
    std::string workingdir = "includes/";
    workingdir.append(argv[1]);
    workingdir.append("/");

    bool debug = true;
    const int it = 1000;

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
	std::string wfile = workingdir + "WSOL.lis";
    std::string dfile = workingdir + "WSOLD.lis";
	if(!fm->Initialise(efile, nfile, mfile, sfile, "mm")) {
		std::cout << "Error while loading the ANSYS field map files." << std::endl;
	}
	fm->EnableMirrorPeriodicityX();
	fm->EnableMirrorPeriodicityY();
    if(debug) {
        fm->PrintRange();
    }
 	fm->SetWeightingField(wfile, "readout");
    fm->SetWeightingField(dfile, "ions");
	
	// Gas setup
	MediumMagboltz* gas = new MediumMagboltz();
  	gas->SetComposition("ar", 40., "co2", 25., "cf4", 35.);
  	gas->SetTemperature(293.15);
  	gas->SetPressure(760.0);	
    //gas->SetMaxElectronEnergy(200.);
  	gas->EnableDebugging();
	gas->Initialise();
	gas->DisableDebugging();
    const double rPenning = 0.57;
    const double lambdaPenning = 0.;
    gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
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
	
	// Setup electron transport
	AvalancheMicroscopic* aval = new AvalancheMicroscopic();
	aval->SetSensor(sensor);
	//aval->EnableAvalancheSizeLimit(1000);
    
	sensor->AddElectrode(fm, "readout");
    sensor->AddElectrode(fm, "ions");
    const double tMin = 0.;
    const double tMax = 50.;
    const double tStep = 0.1;
    const int nTimeBins = int((tMax - tMin)/tStep);
    sensor->SetTimeWindow(0., tStep, nTimeBins);
	aval->EnableSignalCalculation();

	// Setup ion transport
    AvalancheMC* iondrift = new AvalancheMC();
    iondrift->SetSensor(sensor);
    iondrift->EnableSignalCalculation();
    iondrift->SetDistanceSteps(2e-4);

	// Calculate avalanche
	int ne, ni, np, status;
	double x0, y0, z0, t0, e0, x1, y1, z1, t1, e1, x2, y2, z2, t2, e2;
	double vx0, vy0, vz0;
	
    TRandom3 r;
    r.SetSeed(0);
    
    // Prepare trees and saving file
    TString cmd = "mkdir " + workingdir + "output/" + argv[2] + "drift";
    system(cmd);
    TString savefile = workingdir + "output/" + argv[2] + "drift/results.root";
    TFile *file = new TFile(savefile, "RECREATE");
    
    // Prepare tree for electrons
    avalancheE aE;
    TTree *ETree = new TTree("ETree", "Avalanche electrons");
    ETree->Branch("x0", &aE.x0);
    ETree->Branch("y0", &aE.y0);
    ETree->Branch("z0", &aE.z0);
    ETree->Branch("vx0", &aE.vx0);
    ETree->Branch("vy0", &aE.vy0);
    ETree->Branch("vz0", &aE.vz0);
    ETree->Branch("t0", &aE.t0);
    ETree->Branch("e0", &aE.e0);
    ETree->Branch("ne", &aE.ne);
    ETree->Branch("x1", &aE.x1);
    ETree->Branch("y1", &aE.y1);
    ETree->Branch("z1", &aE.z1);
    ETree->Branch("t1", &aE.t1);
    ETree->Branch("e1", &aE.e1);
    ETree->Branch("x2", &aE.x2);
    ETree->Branch("y2", &aE.y2);
    ETree->Branch("z2", &aE.z2);
    ETree->Branch("t2", &aE.t2);
    ETree->Branch("e2", &aE.e2);
    ETree->Branch("status", &aE.status);
    
    // Prepare tree for ions
    avalancheI aI;
    TTree *ITree = new TTree("ITree", "Avalanche ions");
    ITree->Branch("x0", &aI.x0);
    ITree->Branch("y0", &aI.y0);
    ITree->Branch("z0", &aI.z0);
    ITree->Branch("t0", &aI.t0);
    ITree->Branch("ni", &aI.ni);
    ITree->Branch("x1", &aI.x1);
    ITree->Branch("y1", &aI.y1);
    ITree->Branch("z1", &aI.z1);
    ITree->Branch("t1", &aI.t1);
    ITree->Branch("x2", &aI.x2);
    ITree->Branch("y2", &aI.y2);
    ITree->Branch("z2", &aI.z2);
    ITree->Branch("t2", &aI.t2);
    ITree->Branch("status", &aI.status);
    
    // Start iteration
	for(int i=0; i<it; i++) {
        
        if(debug) {
            std::cout << "Progress: " << 100.*(i+1)/it << "%" << std::endl;
        }
        
        // Set random velocity direction
        vx0 = 0.0;
        vy0 = 0.0;
        vz0 = 0.0;
    
        // Set initial energy and time
        e0 = 1.0;
        t0 = 0.0;

        // Generate random (x,y,z) position in unit cell
        x0 = r.Uniform()*pitch/2; 
        y0 = r.Uniform()*pitch*TMath::Sqrt(3)/2;
        //z0 = r.Uniform()*drift;
        z0 = atof(argv[2])*drift;
        
        aval->AvalancheElectron(x0, y0, z0, t0, e0, vx0, vy0, vz0);
		aval->GetAvalancheSize(ne, ni);
		np = aval->GetNumberOfElectronEndpoints();
        
        if(debug) {
            std::cout << "np: " << np << std::endl;
        }
        
        aE.ne = ne;
        aE.x0 = x0;
        aE.y0 = y0;
        aE.z0 = z0;
        aE.vx0 = vx0;
        aE.vy0 = vy0;
        aE.vz0 = vz0;
        aE.t0 = t0;
        aE.e0 = e0;
        aI.ni = ni;
        aI.x0 = x0;
        aI.y0 = y0;
        aI.z0 = z0;
        aI.t0 = t0;
        
        // Loop over all electrons in avalanche
        for(int k=0; k<np; k++) {

			aval->GetElectronEndpoint(k, x1, y1, z1, t1, e1, x2, y2, z2, t2, e2, status);
            aE.x1.push_back(x1);
            aE.y1.push_back(y1);
            aE.z1.push_back(z1);
            aE.t1.push_back(t1);
            aE.e1.push_back(e1);
            aE.x2.push_back(x2);
            aE.y2.push_back(y2);
            aE.z2.push_back(z2);
            aE.t2.push_back(t2);
            aE.e2.push_back(e2);
            aE.status.push_back(status);

            iondrift->DriftIon(x1, y1, z1, t1);
            iondrift->GetIonEndpoint(0, x1, y1, z1, t1, x2, y2, z2, t2, status);
            aI.x1.push_back(x1);
            aI.y1.push_back(y1);
            aI.z1.push_back(z1);
            aI.t1.push_back(t1);
            aI.x2.push_back(x2);
            aI.y2.push_back(y2);
            aI.z2.push_back(z2);
            aI.t2.push_back(t2);
            aI.status.push_back(status);
        }
        
        ETree->Fill();
        ITree->Fill();
        
        // Reset vectors
        aE.x1.clear();
        aE.y1.clear();
        aE.z1.clear();
        aE.t1.clear();
        aE.e1.clear();
        aE.x2.clear();
        aE.y2.clear();
        aE.z2.clear();
        aE.t2.clear();
        aE.e2.clear();
        aE.status.clear();
        
        aI.x1.clear();
        aI.y1.clear();
        aI.z1.clear();
        aI.t1.clear();
        aI.e1.clear();
        aI.x2.clear();
        aI.y2.clear();
        aI.z2.clear();
        aI.t2.clear();
        aI.e2.clear();
        aI.status.clear();
	}

    ETree->Write();
    ITree->Write();
    file->Close();
    
    
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

	return 0;
	app.Run(); 
}

double transferf(double t) {
    
    //const double tau = 1.;
    const double G = 2.;
    return G*0.46*TMath::Power(t/tau,3)*TMath::Exp(-3*t/tau);
}
