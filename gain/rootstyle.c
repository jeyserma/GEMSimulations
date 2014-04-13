{

    TStyle *newStyle= new TStyle("newStyle","new style");
	
	// Histogram title
	newStyle->SetTitleAlign(23) ;
	newStyle->SetTitleX(.5);
	newStyle->SetTitleBorderSize(0);
	newStyle->SetCanvasColor(0);
	
	
	// use plain black on white colors
	Int_t icol=0;
	newStyle->SetFrameBorderMode(icol);
	newStyle->SetCanvasBorderMode(icol);
	newStyle->SetPadBorderMode(icol);
	newStyle->SetPadColor(icol);
	newStyle->SetCanvasColor(icol);
	newStyle->SetStatColor(icol);
	newStyle->SetFillColor(icol);
	newStyle->SetFrameFillColor(icol);
    
    newStyle->SetPaperSize(20,26);
	newStyle->SetPadTopMargin(0.1);
	newStyle->SetPadRightMargin(-1);
	newStyle->SetPadBottomMargin(0.18);
	newStyle->SetPadLeftMargin(0.12);


	// use large fonts
	Int_t font = 42;
	Double_t tsize = 0.05;
	Double_t tsize1 = 0.05;
	newStyle->SetTextFont(font);

	newStyle->SetTitleFontSize(tsize1);
	newStyle->SetTextSize(tsize1);
	newStyle->SetLabelFont(font,"x");
	newStyle->SetTitleFont(font,"x");
	newStyle->SetLabelFont(font,"y");
	newStyle->SetTitleFont(font,"y");
	newStyle->SetLabelFont(font,"z");
	newStyle->SetTitleFont(font,"z");

	newStyle->SetLabelSize(tsize,"x");
	newStyle->SetTitleSize(tsize1,"x");
	newStyle->SetLabelSize(tsize,"y");
	newStyle->SetTitleSize(tsize1,"y");
	newStyle->SetLabelSize(tsize,"z");
	newStyle->SetTitleSize(tsize1,"z");
	newStyle->SetLabelSize(tsize);


	//use bold lines and markers
	newStyle->SetMarkerStyle(20);
	newStyle->SetMarkerSize(1.2);
	newStyle->SetHistLineWidth(1.5);
	
	newStyle->SetLabelOffset(0.02, "x");
    newStyle->SetLabelOffset(0.02, "y");
	newStyle->SetTitleOffset(1.3, "x");
    newStyle->SetTitleOffset(1.6, "y");

	//do not display any of the standard histogram decorations
	newStyle->SetOptTitle(0);
	newStyle->SetOptStat(0);
	newStyle->SetOptFit(0);


	gROOT->SetStyle("Plain");

	 
	// Legends
	newStyle->SetLegendBorderSize(0);
	
	// Statistics
	newStyle->SetPalette(1);
	newStyle->SetCanvasColor(0);
	newStyle->SetFrameFillColor(0);
	newStyle->SetFrameBorderMode(0); 
	newStyle->SetStatBorderSize(1);
	newStyle->SetStatColor(0);

	newStyle->SetStatX(0.85);
  	newStyle->SetStatY(0.85);
	newStyle->SetStatFontSize(0.04);
	newStyle->SetStatFont(42);
	newStyle->SetStatW(0.15);
  	newStyle->SetStatH(0.2);
}
