/*
 * plotGasTable.c
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
#include <string>

// ROOT libs
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TStyle.h>
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

Double_t extractField(TString a);
Double_t extractTownsend(TString a);
Double_t extractDriftVelocity(TString a);
Double_t extractAttachment(TString a);
Double_t extractDiffusion(TString a);

void plotGasTable() {
    
    TString gas = "Ar-CO2-CF4-45-15-40";

    // Prepare some strings
    std::ifstream inp;
    std::string line;
    std::string dir = (std::string)("data/" + gas + "/");
    std::string file = (std::string)("data/" + gas + "/" + gas + ".output");
    inp.open(file.c_str());
    
    // Graphical styles
    gROOT->ProcessLine(".x rootstyle.c");
    gROOT->SetStyle("newStyle");
    gROOT->ForceStyle();
    
    gStyle->SetLabelOffset(0.01, "x");
    gStyle->SetLabelOffset(0.01, "y");
    gStyle->SetTitleOffset(1.3, "x");
    gStyle->SetTitleOffset(0.8, "y");

    // Prepare plots
    TGraph *driftVelocity = new TGraph();
    TGraph *townsend = new TGraph();
    TGraph *diffusionTrans = new TGraph();
    TGraph *diffusionLong = new TGraph();
    TGraph *attachment = new TGraph();
   
    Double_t efield, alpha;
    TString name;   
    int i=-1;
    while(std::getline(inp, line)) {
        
        std::stringstream ss(line);
        name = (TString)line;

        if(name.Contains("ELECTRIC FIELD")) {
            
            efield = extractField(line);
            std::cout << "E-field (kV/cm) " << efield/1000. << std::endl;
            i++;
        }
            
        if(name.Contains("Z DRIFT VELOCITY")) {

            std::cout << "\tDrift velocity (um/ns) " << extractDriftVelocity(line) << std::endl;
            driftVelocity->SetPoint(i, efield/1000., extractDriftVelocity(line));
        }
        
        if(name.Contains("IONISATION RATE /CM.")) {

            townsend->SetPoint(i, efield/1000., extractTownsend(line));
            std::cout << "\tTownsend (1/cm) " << extractTownsend(line) << std::endl;
        }
        
        if(name.Contains("ATTACHMENT RATE /CM.")) {

            attachment->SetPoint(i, efield/1000., extractAttachment(line));
            std::cout << "\tAttachment (1/cm) " << extractAttachment(line) << std::endl;
        }
        
        if(name.Contains("TRANSVERSE DIFFUSION")) {

            diffusionTrans->SetPoint(i, efield/1000., extractDiffusion(line));
            std::cout << "\tTransversal diffusion (cm^2/s) " << extractDiffusion(line) << std::endl;
        }
        
        if(name.Contains("LONGITUDINAL DIFFUSION")) {

            diffusionLong->SetPoint(i, efield/1000., extractDiffusion(line));
            std::cout << "\tLongitudinal diffusion (cm^2/s) " << extractDiffusion(line) << std::endl;
        }  
    }
    
    Double_t maxEfield = efield/1000; // maximum E-field is the last entry
    TString savefile = dir + gas + ".root";;
    TFile f(savefile, "RECREATE");
    
    TCanvas *c1 = new TCanvas("c1", "", 600, 300);
    TLegend *legend = new TLegend(0.15, 0.8, 0.5, 0.90);
    legend->SetTextSize(0.05);
    c1->cd();
	
    // driftVelocity
    driftVelocity->GetXaxis()->SetTitle("E field [kV/cm]");
    driftVelocity->GetYaxis()->SetTitle("v_{d} [#mum/ns]");
    driftVelocity->GetXaxis()->SetLimits(0., maxEfield);
    driftVelocity->SetLineWidth(2);
    legend->AddEntry(driftVelocity, gas, "l");
    driftVelocity->Draw("AC");
    legend->Draw("SAME");
    savefile = dir + "driftVelocity.pdf";
    c1->SaveAs(savefile);
    driftVelocity->Write("driftvelocity");
    legend->Clear();
    c1->Clear();
    
    
    // townsendEl
    townsend->GetXaxis()->SetTitle("E field [kV/cm]");
    townsend->GetYaxis()->SetTitle("#alpha [cm^{-1}]");
    townsend->GetXaxis()->SetLimits(0., 100.);
    townsend->SetLineWidth(2);
    townsend->SetMarkerSize(5);
    townsend->SetMarkerStyle(20);
    legend->AddEntry(townsend, gas, "l");
    townsend->Draw("AC");
    legend->Draw("SAME");
    savefile = dir + "townsend.pdf";
    c1->SaveAs(savefile);
    townsend->Write("townsend");
    legend->Clear();
    c1->Clear();   
    

    // diffusionTrans
    diffusionTrans->GetXaxis()->SetTitle("E field [kV/cm]");
    diffusionTrans->GetYaxis()->SetTitle("D_{T} [#sqrt{cm}]");
    diffusionTrans->GetXaxis()->SetLimits(2., maxEfield);
    diffusionTrans->SetLineWidth(2);
    legend->AddEntry(diffusionTrans, gas, "l");
    diffusionTrans->Draw("AC");
    legend->Draw("SAME");
    savefile = dir + "diffusionTrans.pdf";
    c1->SaveAs(savefile);
    diffusionTrans->Write("diffusionTrans");
    legend->Clear();
    c1->Clear();   
    
    // diffusionLongl
    diffusionLong->GetXaxis()->SetTitle("E field [kV/cm]");
    diffusionLong->GetYaxis()->SetTitle("D_{L} [#sqrt{cm}]");
    diffusionLong->GetXaxis()->SetLimits(2., maxEfield);
    diffusionLong->SetLineWidth(2);
    legend->AddEntry(diffusionLong, gas, "l");
    diffusionLong->Draw("AC");
    legend->Draw("SAME");
    savefile = dir + "diffusionLong.pdf";
    c1->SaveAs(savefile);
    diffusionLong->Write("diffusionLong");
    legend->Clear();
    c1->Clear(); 

    
    // attachment
    attachment->GetXaxis()->SetTitle("E field [kV/cm]");
    attachment->GetYaxis()->SetTitle("k [cm^{-1}]");
    attachment->GetXaxis()->SetLimits(2., maxEfield);
    attachment->SetLineWidth(2);
    legend->AddEntry(attachment, gas, "l");
    attachment->Draw("AC");
    legend->Draw("SAME");
    savefile = dir + "attachment.pdf";
    c1->SaveAs(savefile);
    attachment->Write("attachment");
    legend->Clear();
    c1->Clear();      
    
    c1->Close();
}


Double_t extractField(TString a) {

    std::string str1 = (std::string)a;
    unsigned pos1 = str1.find("=");
    unsigned pos2 = str1.find("VOLTS/CM"); 
    std::string str2 = str1.substr(pos1+1, pos2-pos1-1);
    return atof(str2.c_str());
}

Double_t extractTownsend(TString a) {

    std::string str1 = (std::string)a;
    unsigned pos1 = str1.find("=");
    unsigned pos2 = str1.find("+/-"); 
    std::string str2 = str1.substr(pos1+1, pos2-pos1-1);
    return atof(str2.c_str());
}

Double_t extractDriftVelocity(TString a) {

    std::string str1 = (std::string)a;
    unsigned pos1 = str1.find("=");
    unsigned pos2 = str1.find("MICRONS/NANOSECOND"); 
    std::string str2 = str1.substr(pos1+1, pos2-pos1-1);
    return atof(str2.c_str());
}

Double_t extractAttachment(TString a) {

    std::string str1 = (std::string)a;
    unsigned pos1 = str1.find("=");
    unsigned pos2 = str1.find("+/-"); 
    std::string str2 = str1.substr(pos1+1, pos2-pos1-1);
    return atof(str2.c_str());
}

Double_t extractDiffusion(TString a) {

    std::string str1 = (std::string)a;
    unsigned pos1 = str1.find("=");
    unsigned pos2 = str1.find("+-"); 
    std::string str2 = str1.substr(pos1+1, pos2-pos1-1);
    return atof(str2.c_str());
}
