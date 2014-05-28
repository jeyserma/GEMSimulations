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

struct cluster {
    
    Double_t x0, y0, z0, t0, ec;
    Int_t nc;
    std::vector<Double_t> x1, y1, z1, t1, e1, vx1, vy1, vz1;
};

struct particle {
    
    TString type;
    Double_t x0, y0, z0, vx0, vy0, vz0, e0, t0;
    Int_t noClusters;
    Double_t lambda, stoppingPower;
    std::vector<cluster> clusters;
};

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
void loadGEMconfig(std::string filename, GEMconfig &g);

int main(int argc, char * argv[]) {
    
    TStopwatch watch;
    gRandom = new TRandom3(0); // set random seed
    gROOT->ProcessLine(".L loader.c+"); // Initialize struct objects (see dictionaries)
    
    // Gas setup
    TString gasmixt[6] = { "C5H12", "CF4", "", "60", "40", "" };
    TString output = gasmixt[0] + "-" + gasmixt[1] + "-" + gasmixt[2] + "-" + gasmixt[3] + "-" + gasmixt[4] + "-" + gasmixt[5];
    
    std::string workingdir = "includes/";
    workingdir.append("GEM5"); // Name of the working directory which contains the GEM files
    workingdir.append("/");

    std::string particleType = "mu";
    Double_t particleEnergy = 100.e9;
    bool debug = true;
    Int_t it = 100;

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
    fm->SetWeightingField(wfile, "readout");
    fm->SetWeightingField(dfile, "ions");
	
    // Gas setup
    MediumMagboltz* gas = new MediumMagboltz();
    gas->SetComposition((std::string)gasmixt[0], atof(gasmixt[3]), (std::string)gasmixt[1], atof(gasmixt[4]), (std::string)gasmixt[2], atof(gasmixt[5]));
    gas->SetTemperature(293.15);
    gas->SetPressure(760.0);	
    //gas->SetMaxElectronEnergy(200.);
    gas->EnableDebugging();
    gas->Initialise();
    gas->DisableDebugging();
    //const double rPenning = 0.57;
    //const double lambdaPenning = 0.;
    //gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
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
	
    // Setup electron transport
    AvalancheMicroscopic* aval = new AvalancheMicroscopic();
    aval->SetSensor(sensor);
    //aval->EnableAvalancheSizeLimit(1000);
    
    sensor->AddElectrode(fm, "readout");
    sensor->AddElectrode(fm, "ions");
    const double tMin = 0.;
    const double tMax = 75.;
    const double tStep = 0.2;
    const int nTimeBins = int((tMax - tMin)/tStep);
    sensor->SetTimeWindow(0., tStep, nTimeBins);
    aval->EnableSignalCalculation();
    
    ViewSignal* signalView = new ViewSignal();
    signalView->SetSensor(sensor);
    TH1D* h; // tmp storage of timing histogram

    // Setup ion transport
    AvalancheMC* iondrift = new AvalancheMC();
    iondrift->SetSensor(sensor);
    iondrift->EnableSignalCalculation();
    iondrift->SetDistanceSteps(2e-4);

    // Calculate avalanche
    int ne, ni, np, status, nc;
    double ec, extra;
    double x0, y0, z0, t0, e0, x1, y1, z1, t1, e1, x2, y2, z2, t2, e2, x3, y3, z3, t3, e3;
    double vx0, vy0, vz0, vx1, vy1, vz1;
	
    TRandom3 r;
    r.SetSeed(0);

    TString savefile = workingdir + "output/" + output + ".root";
    TFile f(savefile, "recreate");
    TDirectory *dir = f.mkdir("signals");
    dir->cd();
    
    // Prepare tree for charged particle and clusters
    particle p;
    TTree *pTree = new TTree("pTree", "Charged particle");
    pTree->Branch("x0", &p.x0);
    pTree->Branch("y0", &p.y0);
    pTree->Branch("z0", &p.z0);
    pTree->Branch("vx0", &p.vx0);
    pTree->Branch("vy0", &p.vy0);
    pTree->Branch("vz0", &p.vz0);
    pTree->Branch("t0", &p.t0);
    pTree->Branch("e0", &p.e0);
    pTree->Branch("type", "TString", &p.type);
    pTree->Branch("noClusters", &p.noClusters);
    pTree->Branch("stoppingpower", &p.stoppingPower);
    pTree->Branch("lambda", &p.lambda);
    pTree->Branch("clusters", "std::vector<cluster>", &p.clusters);
  
    // Prepare tree for electrons
    avalancheE aE;
    TTree *ETree = new TTree("eTree", "Avalanche electrons");
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
    TTree *ITree = new TTree("iTree", "Avalanche ions");
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
    std::cout << "---------------------------------------------------------------" << std::endl;
    for(int i=0; i<it; i++) {
        
        if(debug) {
            std::cout << "Progress: " << 100.*(i+1)/it << "%" << std::endl;
        }

	system("mail -s 'timeResolution10' janeysermans@gmail.com <<< 'Progress: " + int2str(i) + "' ");
        
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
        
        // Storage of TRACK coordinates and primary ionization properties
        p.x0 = x0;
        p.y0 = y0;
        p.z0 = z0;
        p.vx0 = vx0;
        p.vy0 = vy0;
        p.vz0 = vz0;
        p.t0 = t0;
        p.e0 = particleEnergy;
        p.type = particleType;
        p.noClusters = 0;
        p.lambda = 1/heed->GetClusterDensity();
        p.stoppingPower = heed->GetStoppingPower();

        // Loop over clusters
        int l=0;
        while(heed->GetCluster(x0, y0, z0, t0, nc, ec, extra)) {
            
            // Skip the clusters which are not in the drift region
            if(z0 > g.driftT) {
                continue;
            }
            
            if(debug) {
                std::cout << "  cluster " << l << " (# electrons = " << nc << ", ec = " << ec <<" )" << std::endl;
            } 

            cluster c;
            // Storage of cluster information
            p.noClusters++;
            c.x0 = x0;
            c.y0 = y0;
            c.z0 = z0;
            c.t0 = t0;
            c.nc = nc; // amount of primary electrons
            c.ec = ec; // energy transferred to gas
            
            // Loop over electrons in cluster
            for(int j=0; j<nc; j++) {

                heed->GetElectron(j, x1, y1, z1, t1, e1, vx1, vy1, vz1);
                
                c.x1.push_back(x1);
                c.y1.push_back(y1);
                c.z1.push_back(z1);
                c.t1.push_back(t1);
                c.e1.push_back(e1);
                c.vx1.push_back(vx1);
                c.vy1.push_back(vy1);
                c.vz1.push_back(vz1);

                // Calculate the drift of electrons and ions
                aval->AvalancheElectron(x1, y1, z1, t1, e1, vx1, vy1, vz1);
                aval->GetAvalancheSize(ne, ni);
                np = aval->GetNumberOfElectronEndpoints();
        
                if(debug) {
                    std::cout << "    avalanche electrons = " << np << std::endl;
                }
                
                aE.ne = ne;
                aE.x0 = x1;
                aE.y0 = y1;
                aE.z0 = z1;
                aE.vx0 = vx1;
                aE.vy0 = vy1;
                aE.vz0 = vz1;
                aE.t0 = t1;
                aE.e0 = e1;
                aI.ni = ni;
                aI.x0 = x1;
                aI.y0 = y1;
                aI.z0 = z1;
                aI.t0 = t1;
                       
                // Loop over all electrons in avalanche
                for(int k=0; k<np; k++) {

                    aval->GetElectronEndpoint(k, x2, y2, z2, t2, e2, x3, y3, z3, t3, e3, status);
                    aE.x1.push_back(x2);
                    aE.y1.push_back(y2);
                    aE.z1.push_back(z2);
                    aE.t1.push_back(t2);
                    aE.e1.push_back(e2);
                    aE.x2.push_back(x3);
                    aE.y2.push_back(y3);
                    aE.z2.push_back(z3);
                    aE.t2.push_back(t3);
                    aE.e2.push_back(e3);
                    aE.status.push_back(status);  

                    iondrift->DriftIon(x2, y2, z2, t2);
                    iondrift->GetIonEndpoint(0, x2, y2, z2, t2, x3, y3, z3, t3, status);
                    aI.x1.push_back(x2);
                    aI.y1.push_back(y2);
                    aI.z1.push_back(z2);
                    aI.t1.push_back(t2);
                    aI.x2.push_back(x3);
                    aI.y2.push_back(y3);
                    aI.z2.push_back(z3);
                    aI.t2.push_back(t3);
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
            
            p.clusters.push_back(c);
            l++;
        }

        pTree->Fill();
        p.clusters.clear();
        
        signalView->PlotSignal("readout");
        h = signalView->GetHistogram();
        h->Write("signal");

        sensor->SetTransferFunction(transferf);
        sensor->ConvoluteSignal();
        signalView->PlotSignal("readout");
        h = signalView->GetHistogram();
        h->Write("signalT");    
        sensor->ClearSignal();
    }
    
    f.cd(); // go to top directory
    ETree->Write();
    ITree->Write();
    pTree->Write();
    f.Close();
    
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
