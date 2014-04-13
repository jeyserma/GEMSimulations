/*
 * gasTable.c
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
#include <TBenchmark.h>

// Garfield++ libs
#include "MediumMagboltz.hh"
#include "Random.hh"
#include "TrackHeed.hh"
#include "ComponentConstant.hh"
#include "Plotting.hh"

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
     
    MediumMagboltz* gas = new MediumMagboltz();
    gas->SetComposition("ar", 50., "co2", 25., "cf4", 25);
    gas->SetTemperature(293.15);
    gas->SetPressure(760.);
  
    gas->PrintGas();
    // Calculation from 0 kV/cm to 40 kV/cm in 100 steps
    gas->SetFieldGrid(0., 40000., 50, false, 0., 0., 1, 0., 0., 1);
    
    gas->EnableDebugging();
    gas->Initialise();  
    gas->DisableDebugging();
    
    // Disable penning transfer
    //const double rPenning = 0.;
    //const double lambdaPenning = 0.;
    //gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
    gas->DisablePenningTransfer();
    
    // Generate gas table
    gas->GenerateGasTable(5, true);
    
    gas->WriteGasFile("data/Ar-CO2-CF4-50-25-25.gas");
    
    std::cout << "--------------- TIMING ---------------" << std::endl;
    std::cout << watch.CpuTime() << std::endl;

    return 0;
    app.Run(); 
}
