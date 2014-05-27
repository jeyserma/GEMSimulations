/*
 * clusterDensity.c
 * Author: Jan Eysermans 
 * 2014
 */

#include <sstream>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>

// ROOT libs
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <TPad.h>
#include <TMath.h>
#include <TLegend.h>
#include <TString.h>
#include <TGraph.h>
#include <TLatex.h>

// Garfield++ libs
#include "MediumMagboltz.hh"
#include "Sensor.hh"
#include "Random.hh"
#include "Plotting.hh"
#include "DriftLineRKF.hh"
#include "TrackHeed.hh"
#include "SolidBox.hh"
#include "ComponentConstant.hh"
#include "GeometrySimple.hh"

using namespace Garfield;
const std::string GARFIELD = "/afs/cern.ch/user/j/jeyserma/garfieldpp/";

int main(int argc, char * argv[]) {
    
    //CH4           Methane (CH_{4})
    //CD4           Deuterated methane (CD_{4})
    //C2H6          Ethane (C_{2}H_{6})
    //C3H8          Propane (C_{3}H_{8})
    //nC4H10        n-butane (C_{4}H_{10})
    //iC4H10        i-butane (C_{4}H_{10})
    //nC5H12        Pentane (C_{5}H_{12})
    //neo-C5H12     Pentane (C_{5}H_{12}
    //C2H4          Ethene (C_{2}H_{4})
    //C2H2          Acetylene (C_{2}H_{2})
    //C3H6          Propene (C_{3}H_{6})
    
    //CO2           Carbon dioxide (CO_{2})
    //CO            Carbon monoxide (CO)
   
    //CF4           Tetrafluormethane (CF_{4})
    //CHF3          Fluoroform (CHF_{3})
    //C2F6          Hexafluoroethane (C_{2}F_{6})
    //C3F8          Octafluoropropane (C_{3}F_{8})
    //SF6           Sulfur hexafluoride (SF_{6})
    //BF3           Boron trifluoride (BF_{3})
    //CF3Br         Bromotrifluoromethane (CF_{3}Br)
    
    TString gasmixt[6] = { "CF3Br", "", "", "100", "", "" };
    TString title = "Bromotrifluoromethane (CF_{3}Br)";

	// Settings
    Double_t start = TMath::Log10(0.2); // start in GeV
    Double_t end = TMath::Log10(1000); // start in GeV
    Double_t steps = 100;
    TString output = gasmixt[0] + "-" + gasmixt[1] + "-" + gasmixt[2] + "-" + gasmixt[3] + "-" + gasmixt[4] + "-" + gasmixt[5];
    
    // Graphical styles
	gROOT->ProcessLine(".x rootstyle.c");
	gROOT->SetStyle("newStyle");
    gROOT->ForceStyle();
    
	gStyle->SetLabelOffset(0.01, "x");
	gStyle->SetLabelOffset(0.01, "y");
	gStyle->SetTitleOffset(1.3, "x");
	gStyle->SetTitleOffset(0.8, "y");

	MediumMagboltz* gas = new MediumMagboltz();
    gas->SetComposition((std::string)gasmixt[0], atof(gasmixt[3]), (std::string)gasmixt[1], atof(gasmixt[4]), (std::string)gasmixt[2], atof(gasmixt[5]));
  	gas->SetTemperature(293.15);
  	gas->SetPressure(760.0);	

  	gas->EnableDebugging();
	gas->Initialise();
	gas->DisableDebugging();
    
    // Disable penning transfer
    gas->DisablePenningTransfer();
	
    // Make solid rectangular box and apply Ex = 100 V / cm
	const double width = 1.;
	SolidBox* box = new SolidBox(width / 2., 0., 0., width / 2., 10., 10.);
	GeometrySimple* geo = new GeometrySimple();
	geo->AddSolid(box, gas);
	
	ComponentConstant* cmp = new ComponentConstant();
	cmp->SetGeometry(geo);
	cmp->SetElectricField(100., 0., 0.);

	Sensor* sensor = new Sensor();
	sensor->AddComponent(cmp);

	TrackHeed* heed = new TrackHeed();
    heed->SetParticle("mu");
    heed->SetSensor(sensor);

    double x0 = 0., y0 = 0., z0 = 0., t0 = 0.;
    double dx0 = 1., dy0 = 0., dz0 = 0.; 
    
    TGraph *graph = new TGraph();
    Double_t E = start;
    Int_t j=0;
    
    Double_t inc = (end-start)/steps;
    
    while(E <= end) {

    	//heed->SetBetaGamma(i);
        heed->SetEnergy(TMath::Power(10, E)*1.0e9);
    	heed->NewTrack(x0, y0, z0, t0, dx0, dy0, dz0);
    	
    	graph->SetPoint(j, TMath::Power(10, E), heed->GetClusterDensity());
    	E = E + inc;
    	j++;
    }
    
    TCanvas *c1 = new TCanvas("c1", "Cluster density", 600, 300);
    TLegend *legend = new TLegend(0.5, 0.32, 0.9, 0.4);
	legend->AddEntry(graph, title, "l");
    legend->SetTextSize(0.05);

	c1->cd();
    c1->SetLogx();
    graph->GetXaxis()->SetTitle("E [GeV]");
	graph->GetYaxis()->SetTitle("#lambda^{-1} [cm^{-1}]");
	graph->GetXaxis()->SetLimits(TMath::Power(10, start), TMath::Power(10, end));
    graph->SetLineWidth(2);
    graph->Draw("AC");
    legend->Draw("SAME");
    c1->SaveAs("data/" + output + ".pdf");
    
	// Write to file
	TFile *f = new TFile("data/" + output + ".root", "recreate");
	f->cd();
    graph->Write();
    f->Close();
    
    return 0;
}
