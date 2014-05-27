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
#include <TH1D.h>
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

struct GEMconfig {
    
    Int_t GEMlayers;
    Double_t driftT, inductT, driftV, inductV, totalT, totalV;
    Double_t kapton, metal, pitch, outdia, middia, rim;
    std::vector<Double_t> VGEM, transferV, transferT;
};

void loadGEMconfig(std::string filename, GEMconfig &g);
void parseGasName(TString gasmixt[6], TString &n);
int simulation(TString gasmixt[6], TString GEMid);
int simulationParam(TString gasmixt[6], TString GEMid);

int main(int argc, char * argv[]) {
    
    TString gasmixt[6] = { "C3F8", "CO2", "", "50", "50", "" };
    simulation(gasmixt, "GEM3");
    //simulationParam(gasmixt, "GEM3m");
    

   
    return 0;
}

int simulationParam(TString gasmixt[6], TString GEMid) {
    
    TStopwatch watch;
    gRandom = new TRandom3(0); // set random seed

    // Gas setup
    TString gasname;
    parseGasName(gasmixt, gasname);
    std::string workingdir = (std::string)("includes/" + GEMid + "/");
    TString savedir = "includes/" + GEMid + "/output/" + gasname + "/";
    system("mkdir " + savedir);

	
	double Emin = 1.e9; // 5 GeV
	double Emax = 1000.e9; // 1 TeV
    std::string particleType = "mu";
    Double_t particleEnergy = 1000.e9; // default energy
    //std::string particleType = "pi-";
    //Double_t particleEnergy = 350.e6;
    Double_t thrsFrac; // threshold fraction
    bool debug = false;
    Int_t it = 10000;
    
    // Graphical styles
    gROOT->ProcessLine(".x rootstyle.c");
    gROOT->SetStyle("newStyle");
    gROOT->ForceStyle();
    
	gStyle->SetLabelOffset(0.01, "x");
	gStyle->SetLabelOffset(0.01, "y");
	gStyle->SetTitleOffset(1.3, "x");
	gStyle->SetTitleOffset(0.9, "y");

    // Load GEM dimensions
    GEMconfig g;
    loadGEMconfig(workingdir, g);
    
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
    
    // Gas setup
    MediumMagboltz* gas = new MediumMagboltz();
    gas->SetComposition((std::string)gasmixt[0], atof(gasmixt[3]), (std::string)gasmixt[1], atof(gasmixt[4]), (std::string)gasmixt[2], atof(gasmixt[5]));
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
    //gas->LoadIonMobility(GARFIELD + "Data/IonMobility_CO2+_CO2");
     
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
    sensor->SetArea(-50.*(g.pitch), -50.*(g.pitch), 0.0, 50.*(g.pitch), 50.*(g.pitch), g.totalT);
    
    // Setup HEED
    TrackHeed* heed = new TrackHeed();
    heed->SetSensor(sensor);
    //heed->DisableDeltaElectronTransport();
    heed->SetParticle(particleType);
    heed->SetMomentum(particleEnergy);
    if(debug) {
        heed->EnableDebugging();
    }
	

    // Calculate avalanche
    int nc;
    double ec, extra;
    double x0, y0, z0, t0, e0, x1, y1, z1, t1, e1;
    double vx0, vy0, vz0, vx1, vy1, vz1;
	
    TRandom3 r;
    r.SetSeed(0);

    
    TCanvas *c1 = new TCanvas("c1", "", 600, 300);
    std::vector<Double_t> cluster_z0; // z-coordinate of clusters
    std::vector<Double_t> cluster_t0; // t-coordinate of clusters
    TH1F* z_lastCluster = new TH1F("z_lastCluster", "", 100, 0., 1.05*g.driftT*10000.);
    TH1F* t_lastCluster = new TH1F("t_lastCluster", "", 200, 0., 7.); // units in ps

    // Start iteration
    std::cout << "---------------------------------------------------------------" << std::endl;
    int i=0;
    while(i<it) {
        
        if(debug) {
            std::cout << "Progress: " << 100.*(i+1)/it << "%" << std::endl;
        }
        
        particleEnergy = r.Uniform(Emin, Emax);
        heed->SetMomentum(particleEnergy);
        
        // Set random velocity direction, in positive hemisphere
        const double ctheta = 1. - r.Uniform();
        const double stheta = sqrt(1. - ctheta*ctheta);
        const double phi = TwoPi*r.Uniform();
        vx0 = cos(phi)*stheta;
        vy0 = sin(phi)*stheta;
        vz0 = ctheta;

	if(TMath::Sin(22.83*TMath::Pi()/180.) < stheta || stheta < TMath::Sin(13.96*TMath::Pi()/180.)  ) {
            
		continue;
        }

	i++;
	//std::cout << TMath::ASin(stheta)*180./TMath::Pi() <<  std::endl;
  
        // Set initial time
        t0 = 0.0;

        // Generate random (x,y) position in unit cell
        x0 = r.Uniform()*g.pitch/2; 
        y0 = r.Uniform()*g.pitch*TMath::Sqrt(3)/2;
        z0 = 0.;
        
        // Set muon perpendicular in midpoint of GEM
        x0 = 0.;
        y0 = 0.;
        //vx0 = TMath::Sin(angle);
        //vy0 = 0.;
        //vz0 = TMath::Cos(angle);
        
        heed->NewTrack(x0, y0, z0, t0, vx0, vy0, vz0); // generate particle track
       
        // Loop over clusters
        int l=0;
        while(heed->GetCluster(x0, y0, z0, t0, nc, ec, extra)) {
            
            // Skip the clusters which are not in the drift region
            if(z0 > g.driftT) {
                continue;
            }
            
           // std::cout << z0*10000. << std::endl;
            
            if(debug) {
                std::cout << "  cluster " << l << " (# electrons = " << nc << ")" << " z =" << z0*10000. << std::endl;
            } 

            cluster_z0.push_back(z0*10000.); // convert cm to µm
            cluster_t0.push_back(t0*1000.); // convert ns to ps  
            
            // Loop over electrons in cluster
            for(int j=0; j<nc; j++) {

                //heed->GetElectron(j, x1, y1, z1, t1, e1, vx1, vy1, vz1);
            }     
            l++;
        }
        
        std::sort(cluster_z0.begin(), cluster_z0.end());
        std::sort(cluster_t0.begin(), cluster_t0.end());
        //z_lastCluster->Fill(cluster_z0[cluster_z0.size()-1]);
        //t_lastCluster->Fill(cluster_t0[cluster_t0.size()-1]);
        
        
        // If vector has zero length --> no clusters created..
        if(cluster_z0.size() != 0) {
        	z_lastCluster->Fill(cluster_z0[cluster_z0.size()-1]);
        	t_lastCluster->Fill(cluster_t0[cluster_t0.size()-1]);
        }

        cluster_z0.clear();
        cluster_t0.clear();
        
    }
    
    gStyle->SetStatX(0.55);
    gStyle->SetStatY(0.85);
    gStyle->SetStatFont(42);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.40);
    gStyle->SetStatFontSize(0.06);
    gStyle->SetOptStat("emr");
    gStyle->SetStatBorderSize(1);
    
 
    c1->cd();
    z_lastCluster->SetStats(1);
    z_lastCluster->GetXaxis()->SetTitle("z-coordinate [#mum]");
    z_lastCluster->GetYaxis()->SetTitle("# last clusters / 10.5 #mum");
    z_lastCluster->GetYaxis()->SetTitleOffset(1.0);
    z_lastCluster->SetFillColor(19);
    z_lastCluster->Draw();
    c1->SaveAs((TString)(savedir + "z_lastClusterMuonEnergyCMS.pdf"));
    c1->Clear();

    
   
    c1->cd();
    t_lastCluster->SetStats(1);
    t_lastCluster->GetXaxis()->SetTitle("time [ps]");
    t_lastCluster->GetYaxis()->SetTitle("# last clusters / 0.035 ps");
    t_lastCluster->GetYaxis()->SetTitleOffset(1.0);
    t_lastCluster->SetFillColor(19);
    t_lastCluster->Draw();
    c1->SaveAs((TString)(savedir + "t_lastClusterMuonEnergyCMS.pdf"));
    c1->Clear();
    

    std::cout << "--------------- TIMING ---------------" << std::endl;
    std::cout << watch.CpuTime() << std::endl;

    return 0;
}

int simulation(TString gasmixt[6], TString GEMid) {
    
    TStopwatch watch;
    gRandom = new TRandom3(0); // set random seed

    // Gas setup
    TString gasname;
    parseGasName(gasmixt, gasname);
    std::string workingdir = (std::string)("includes/" + GEMid + "/");
    TString savedir = "includes/" + GEMid + "/output/" + gasname + "/";
    system("mkdir " + savedir);

    std::string particleType = "mu";
    Double_t particleEnergy = 100.e9;
    //std::string particleType = "pi-";
    //Double_t particleEnergy = 350.e6;
    Double_t thrsFrac; // threshold fraction
    bool debug = false;
    Int_t it = 10000;
    
    // Graphical styles
    gROOT->ProcessLine(".x rootstyle.c");
    gROOT->SetStyle("newStyle");
    gROOT->ForceStyle();
    
	gStyle->SetLabelOffset(0.01, "x");
	gStyle->SetLabelOffset(0.01, "y");
	gStyle->SetTitleOffset(1.3, "x");
	gStyle->SetTitleOffset(0.9, "y");

    // Load GEM dimensions
    GEMconfig g;
    loadGEMconfig(workingdir, g);
    
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
    
    // Gas setup
    MediumMagboltz* gas = new MediumMagboltz();
    gas->SetComposition((std::string)gasmixt[0], atof(gasmixt[3]), (std::string)gasmixt[1], atof(gasmixt[4]), (std::string)gasmixt[2], atof(gasmixt[5]));
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
    //gas->LoadIonMobility(GARFIELD + "Data/IonMobility_CO2+_CO2");
     
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
    sensor->SetArea(-5.*(g.pitch), -5.*(g.pitch), 0.0, 5.*(g.pitch), 5.*(g.pitch), g.totalT);
    
    // Setup HEED
    TrackHeed* heed = new TrackHeed();
    heed->SetSensor(sensor);
    //heed->DisableDeltaElectronTransport();
    heed->SetParticle(particleType);
    heed->SetMomentum(particleEnergy);
    if(debug) {
        heed->EnableDebugging();
    }
	

    // Calculate avalanche
    int nc;
    double ec, extra;
    double x0, y0, z0, t0, e0, x1, y1, z1, t1, e1;
    double vx0, vy0, vz0, vx1, vy1, vz1;
	
    TRandom3 r;
    r.SetSeed(0);

    
    TCanvas *c1 = new TCanvas("c1", "", 600, 300);
    std::vector<Double_t> cluster_z0; // z-coordinate of clusters
    std::vector<Double_t> cluster_t0; // t-coordinate of clusters
    TH1F* z_lastCluster = new TH1F("z_lastCluster", "", 100, 0., 1.05*g.driftT*10000.);
    TH1F* t_lastCluster = new TH1F("t_lastCluster", "", 200, 0., 7.); // units in ps

    // Start iteration
    std::cout << "---------------------------------------------------------------" << std::endl;
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
        x0 = r.Uniform()*g.pitch/2; 
        y0 = r.Uniform()*g.pitch*TMath::Sqrt(3)/2;
        z0 = 0.;
        
        // Set muon perpendicular in midpoint of GEM
        x0 = 0.;
        y0 = 0.;
        vx0 = 0.;
        vy0 = 0.;
        vz0 = 1.;
        
        heed->NewTrack(x0, y0, z0, t0, vx0, vy0, vz0); // generate particle track
       
        // Loop over clusters
        int l=0;
        while(heed->GetCluster(x0, y0, z0, t0, nc, ec, extra)) {
            
            // Skip the clusters which are not in the drift region
            if(z0 > g.driftT) {
                continue;
            }
            
            if(debug) {
                std::cout << "  cluster " << l << " (# electrons = " << nc << ")" << " z =" << z0*10000. << std::endl;
            } 

            cluster_z0.push_back(z0*10000.); // convert cm to µm
            cluster_t0.push_back(t0*1000.); // convert ns to ps  
            
            // Loop over electrons in cluster
            for(int j=0; j<nc; j++) {

                heed->GetElectron(j, x1, y1, z1, t1, e1, vx1, vy1, vz1);
            }     
            l++;
        }
        
        std::sort(cluster_z0.begin(), cluster_z0.end());
        std::sort(cluster_t0.begin(), cluster_t0.end());
        //z_lastCluster->Fill(cluster_z0[cluster_z0.size()-1]);
        //t_lastCluster->Fill(cluster_t0[cluster_t0.size()-1]);
        
        
        // If vector has zero length --> no clusters created..
        if(cluster_z0.size() != 0) {
        	z_lastCluster->Fill(cluster_z0[cluster_z0.size()-1]);
        	t_lastCluster->Fill(cluster_t0[cluster_t0.size()-1]);
        }

        cluster_z0.clear();
        cluster_t0.clear();
    }
    
    gStyle->SetStatX(0.55);
    gStyle->SetStatY(0.85);
    gStyle->SetStatFont(42);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.40);
    gStyle->SetStatFontSize(0.06);
    gStyle->SetOptStat("emr");
    gStyle->SetStatBorderSize(1);
    
 
    c1->cd();
    z_lastCluster->SetStats(1);
    z_lastCluster->GetXaxis()->SetTitle("z-coordinate [#mum]");
    z_lastCluster->GetYaxis()->SetTitle("# last clusters / 10.5 #mum");
    z_lastCluster->GetYaxis()->SetTitleOffset(1.0);
    z_lastCluster->SetFillColor(19);
    z_lastCluster->Draw();
    c1->SaveAs((TString)(savedir + "z_lastCluster.pdf"));
    c1->Clear();

    
   
    c1->cd();
    t_lastCluster->SetStats(1);
    t_lastCluster->GetXaxis()->SetTitle("time [ps]");
    t_lastCluster->GetYaxis()->SetTitle("# last clusters / 0.035 ps");
    t_lastCluster->GetYaxis()->SetTitleOffset(1.0);
    t_lastCluster->SetFillColor(19);
    t_lastCluster->Draw();
    c1->SaveAs((TString)(savedir + "t_lastCluster.pdf"));
    c1->Clear();
    
     /*
    
    c1->cd();
    z_lastCluster->SetStats(1);
    z_lastCluster->GetXaxis()->SetTitle("z-coordinate [#mum]");
    z_lastCluster->GetYaxis()->SetTitle("# last clusters / 10.5 #mum");
    z_lastCluster->SetFillColor(19);
    z_lastCluster->Draw();
    
    
    TF1 *exp = new TF1("exp", "[1]*TMath::Exp(-x*[0])", 0, 1e4*g.driftT);
    exp->SetParameter(0, heed->GetClusterDensity()*1e-4);
    exp->SetParameter(1, 1);
    //pois->SetParName(0,"#mu"); // mean
    //pois->SetParName(1,"norm. cte");
    z_lastCluster->Fit("exp", "R");
    //Double_t param = pois->GetParameter(0);
    z_lastCluster->GetFunction("exp")->SetLineColor(kRed);
    gStyle->SetOptStat("emr"); // emr = entries, mean, rms
    gStyle->SetOptFit(1);
    z_lastCluster->SetStats(1);
    c1->Update();
    c1->SaveAs((TString)workingdir + "output/" + gasname + "/cluster_total.pdf");
    c1->Clear();
    */

    std::cout << "--------------- TIMING ---------------" << std::endl;
    std::cout << watch.CpuTime() << std::endl;

    return 0;
}

double transferf(double t) {
    
    const double tau = 25.;
    const double G = 2.;
    return G*0.46*TMath::Power(t/tau,3)*TMath::Exp(-3*t/tau);
}

void loadGEMconfig(std::string workingdir, GEMconfig &g) {
    
    std::ifstream inp;
    std::string line;
    std::string file = workingdir + "geometry.txt";
    inp.open(file.c_str());
    
    TString name;
    Double_t value;
    
    g.GEMlayers = 5;
    
    while(std::getline(inp, line)) {
        
        std::stringstream ss(line);
        ss >> name >> value;
        
        if(name.Contains("layers")) g.GEMlayers = value;
        if(name.Contains("kapton")) g.kapton = value;
        if(name.Contains("metal")) g.metal = value;
        if(name.Contains("pitch")) g.pitch = value;
        if(name.Contains("rim")) g.rim = value;
        if(name.Contains("outdia")) g.outdia = value;
        if(name.Contains("middia")) g.middia = value;
        if(name.Contains("driftRegionVoltage")) g.driftV = value;
        if(name.Contains("inductRegionVoltage")) g.inductV = value;
        if(name.Contains("driftRegionThickness")) g.driftT = value;
        if(name.Contains("inductRegionThickness")) g.inductT = value;
        if(name.Contains("totalThickness")) g.totalT = value;
        if(name.Contains("totalVoltage")) g.totalV = value;
        if(name.Contains("GEMVoltage")) g.VGEM.push_back(value);
        if(name.Contains("transferVoltage")) g.transferV.push_back(value);
        if(name.Contains("transferThickness")) g.transferT.push_back(value);
    }
}

void parseGasName(TString gasmixt[6], TString &n) {
    
    if(gasmixt[1] == "") {
        n = gasmixt[0] + "-" + gasmixt[3];
    }
    else if(gasmixt[2] == "") {
        n = gasmixt[0] + "-" + gasmixt[1] + "-" + gasmixt[3] + "-" + gasmixt[4];
    }
    else {
        n = gasmixt[0] + "-" + gasmixt[1] + "-" + gasmixt[2] + "-" + gasmixt[3] + "-" + gasmixt[4] + "-" + gasmixt[5];
    }
}
