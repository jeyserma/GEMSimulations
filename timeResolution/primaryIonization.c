/*
 * primaryIonization.c
 * Author: Jan Eysermans 
 * 2014
 */

// General libs
#include <iostream>

// ROOT libs
#include <TROOT.h>
#include <TFile.h>
#include <TApplication.h>
#include <TMath.h>
#include <TTree.h>
#include <TBranch.h>

// Garfield++ libs
#include "ComponentAnsys123.hh"
#include "MediumMagboltz.hh"
#include "Sensor.hh"
#include "Random.hh"
#include "TrackHeed.hh"


using namespace Garfield;
const std::string GARFIELD = "/afs/cern.ch/user/j/jeyserma/garfieldpp/";

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

int main(int argc, char * argv[]) {
    
    // Charged particle settings
    std::string particleType = "mu";
    Double_t particleEnergy = 100.e9; 
    
    bool debug = true;
    const int it = 100;

	TApplication app("app", &argc, argv);
 	gRandom = new TRandom3(0); // set random seed
    gROOT->ProcessLine(".L loader.c+"); // Initialize struct objects (see dictionaries)
    
    std::string workingdir = "includes/";
    workingdir.append(argv[1]);
    workingdir.append("/");

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
    heed->SetParticle(particleType);
    heed->SetMomentum(particleEnergy); // 100 GeV
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
    
    TString savefile = workingdir + "output/primaryIonization.root";
    TFile *file = new TFile(savefile, "RECREATE");
    
    particle p;
    cluster c;
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
        p.e0 = particleEnergy;
        p.type = particleType;
        p.noClusters = 0;
        p.lambda = 1/heed->GetClusterDensity();
        p.stoppingPower = heed->GetStoppingPower();
        
        heed->NewTrack(x0, y0, z0, t0, vx0, vy0, vz0);
        
        // Loop over clusters
        while(heed->GetCluster(x1, y1, z1, t1, nc, ec, extra)) {
            
            p.noClusters++;
            c.x0 = x1;
            c.y0 = y1;
            c.z0 = z1;
            c.nc = nc;
            c.ec = ec;
         
            // Skip the clusters which are not in the drift region
            if(z1 > drift) {
                continue;
            }

            for(int j=0; j<nc; j++) {

                heed->GetElectron(j, x2, y2, z2, t2, e2, vx2, vy2, vz2);
                    
                c.x1.push_back(x2);
                c.y1.push_back(y2);
                c.z1.push_back(z2);
                c.t1.push_back(t2);
                c.e1.push_back(e2);
                c.vx1.push_back(vx2);
                c.vy1.push_back(vy2);
                c.vz1.push_back(vz2);
            }     

            p.clusters.push_back(c);
        }

        pTree->Fill();
        p.clusters.clear(); 
	}
   
    pTree->Write();
    file->Close();
    
	return 0;
	app.Run(); 
}
