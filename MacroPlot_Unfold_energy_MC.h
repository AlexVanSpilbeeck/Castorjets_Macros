//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <TThread.h>
#include <TStopwatch.h>
#include <TLegend.h>
#include <TPad.h>
#include <TColor.h>
#include <TText.h>
#include <TLine.h>
#include <TPaletteAxis.h>

//STANDARD C++ INCLUDES
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <cstring>
#include <vector>

#include "../../../RooUnfold-1.1.1/src/RooUnfold.h"
#include "../../../RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"
#include "../../../RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "../../../RooUnfold-1.1.1/src/RooUnfold.h"
#include "../../../RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "color.h"
using namespace std;

/*****************************************
* Plot the bin-by-bin correction factors *
*****************************************/

void PlotCorrectionFactors_BinByBin(TH1D* &histo, TString filename , TString legend, TString drawoption){
   gROOT->SetStyle ("Plain");
   gStyle->SetPalette(1);
   gSystem->Load("~/RooUnfold-1.1.1/libRooUnfold.so");

   /* Open DATA and MC */

   TFile *_file0 = TFile::Open( filename, "Read");
 
   /* Extract GEN and DET (MC) and DATA distribution */

   TH1D *hDet = (TH1D*)_file0->Get("hCastorJet_energy");	hDet->Sumw2();
   TH1D *hGen = (TH1D*)_file0->Get("hGenJet_energy");	hGen->Sumw2();

   histo = (TH1D*)hDet->Clone(legend);  
   histo->SetName("Gen/Det energy ratio");

   for(int bin = 0; bin < hDet->GetNbinsX(); bin++){
     double det_value = hDet->GetBinContent(bin),
	    gen_value =  hGen->GetBinContent(bin);
     double det_error = hDet->GetBinError(bin),
            gen_error =  hGen->GetBinError(bin);
	    
     double original_ratio;
     ( det_value == 0) ? original_ratio = 0 : original_ratio = gen_value/det_value;
     
     cout << "bin gen det ratio " << bin << "\t" << gen_value << "\t" << det_value << "\t" << original_ratio << endl;

     histo->SetBinContent(bin, original_ratio);    
   }
 }
 
 
// --------------------------------------------------------------------------------------


/***************************************
* Plot the bayesian correction factors *
***************************************/

void PlotCorrectionFactors_Bayesian(TH1D* &histo, TString filename , TString legend, TString drawoption){
   gROOT->SetStyle ("Plain");
   gStyle->SetPalette(1);
   gSystem->Load("~/RooUnfold-1.1.1/libRooUnfold.so");

   /* Open DATA and MC */

   TFile *_file0 = TFile::Open( filename, "Read");
 
   if( _file0->GetListOfKeys()->Contains( "response" ) ){
 
   /* Extract GEN and DET (MC) and DATA distribution */

   RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response");

   TH1D *hDet = (TH1D*)_file0->Get("hCastorJet_energy");	hDet->Sumw2();
   TH1D *hGen = (TH1D*)_file0->Get("hGenJet_energy");	hGen->Sumw2();

   int nbins_x = hGen->GetXaxis()->GetNbins();
   double x_axisMax = hGen->GetXaxis()->GetBinUpEdge( nbins_x );
  
   /****************/
   /** Unfold MC. **/
   /****************/

   RooUnfoldBayes unfold_mc_bayes (response, hDet, 4);
     TH1D* hReco_mc_bayes= (TH1D*) unfold_mc_bayes.Hreco();

    hReco_mc_bayes->SetLineStyle(3);
     hReco_mc_bayes->SetLineColor(kGreen);
     hReco_mc_bayes->SetLineWidth(2);
     hReco_mc_bayes->SetName("hReco_bayes_1");
     
   histo = (TH1D*)hReco_mc_bayes->Clone(legend);
   
   }
}


// --------------------------------------------------------------------------------------


/***********************************************
* Plot the bin-by-bin unfolding from RooUnfold *
***********************************************/

void PlotCorrectionFactors_BinByBin_Unfolded(TH1D* &histo, TString filename , TString legend, TString drawoption){
   gROOT->SetStyle ("Plain");
   gStyle->SetPalette(1);
   gSystem->Load("~/RooUnfold-1.1.1/libRooUnfold.so");

   /* Open DATA and MC */

   TFile *_file0 = TFile::Open( filename, "Read");
   if( _file0->GetListOfKeys()->Contains( "response" ) ){
   
 
   /* Extract GEN and DET (MC) and DATA distribution */

   RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response");

   TH1D *hDet = (TH1D*)_file0->Get("hCastorJet_energy");	hDet->Sumw2();
   TH1D *hGen = (TH1D*)_file0->Get("hGenJet_energy");	hGen->Sumw2();

   int nbins_x = hGen->GetXaxis()->GetNbins();
   double x_axisMax = hGen->GetXaxis()->GetBinUpEdge( nbins_x );
  
   /****************/
   /** Unfold MC. **/
   /****************/

   RooUnfoldBinByBin unfold_bbb (response, hDet);
     TH1D* hReco_mc_bbb = (TH1D*) unfold_bbb.Hreco();

    hReco_mc_bbb->SetLineStyle(3);
     
   histo = (TH1D*)hReco_mc_bbb->Clone(legend);
   }
}



// --------------------------------------------------------------------------------------

/**********************************
* Plot the Det level distribution *
**********************************/

void Plot_DetLevel(TH1D* &histo, TString filename , TString legend, TString drawoption){
   gROOT->SetStyle ("Plain");
   gStyle->SetPalette(1);
   gSystem->Load("~/RooUnfold-1.1.1/libRooUnfold.so");

   /* Open DATA and MC */

   TFile *_file0 = TFile::Open( filename, "Read");
   if( _file0->GetListOfKeys()->Contains( "response" ) ){
   
 
   /* Extract DET (MC) distribution */

   TH1D *hDet = (TH1D*)_file0->Get("hCastorJet_energy");	hDet->Sumw2();
   histo = (TH1D*)hDet->Clone(legend);  
   }
 }
 
 
// --------------------------------------------------------------------------------------

/**********************************
* Plot the Gen level distribution *
**********************************/

void Plot_GenLevel(TH1D* &histo, TString filename , TString legend, TString drawoption){
   gROOT->SetStyle ("Plain");
   gStyle->SetPalette(1);
   gSystem->Load("~/RooUnfold-1.1.1/libRooUnfold.so");

   /* Open DATA and MC */

   TFile *_file0 = TFile::Open( filename, "Read");
   if( _file0->GetListOfKeys()->Contains( "response" ) ){
   
 
   /* Extract GEN distribution */

   TH1D *hGen = (TH1D*)_file0->Get("hGenJet_energy");	hGen->Sumw2();
   histo = (TH1D*)hGen->Clone(legend);  
   }
 } 
  

// --------------------------------------------------------------------------------------

/**********************************
* Plot the Gen level eta distribution *
**********************************/



void Plot_GenLevelEta(TH1D* &histo, TString filename , TString legend, TString drawoption){
   gROOT->SetStyle ("Plain");
   gStyle->SetPalette(1);
   gSystem->Load("~/RooUnfold-1.1.1/libRooUnfold.so");

   /* Open DATA and MC */

   TFile *_file0 = TFile::Open( filename, "Read");
 
   /* Extract GEN and DET (MC) and DATA distribution */
   
   if( _file0->GetListOfKeys()->Contains("hEtaDiff") ){ 
   
     TH1D *hGen = (TH1D*)_file0->Get("hEtaDiff");	hGen->Sumw2();
     histo = (TH1D*)hGen->Clone(legend);  
   }
 } 
 
  
// --------------------------------------------------------------------------------------

/**********************************
* Plot the JES distribution *
**********************************/  
  
void Plot_JES_vs_E(TH1D* &histo, TString filename , TString legend, TString variable){

   /* Open DATA and MC */
   TFile *_file0 = TFile::Open( filename, "Read");
   cout << "$$$ " << variable << endl;
   if( _file0->GetListOfKeys()->Contains( variable ) ){ 
   
   cout << "$$$ " << filename << "\tcontains " << variable << endl;

      /* Declare necessary values. */
      // -- fit range.
      double fitlow = -1.;
      double fithigh = 5.;
    
      // -- vectors that will become our graphs.
      std::vector<double> sigma_, mu_, chi2_, E_axis;     
  
     /* Extract 2D hist and apply fits */
     TH2D *hJER = (TH2D*)_file0->Get( variable ); 	// The original 2D.
     cout << "$$$ Rebin" << endl;
     hJER->RebinX(4);
     cout << "$$$ Rebinned" << endl;
     TH1D *hE = (TH1D*)hJER->ProjectionX("e", 1, -1,"");	// The energy values.
     cout << "$$$ hJER" << endl; 
 
     for(int bin_E = 0; bin_E < hE->GetNbinsX(); bin_E++){
      cout << "$$$ bin " << bin_E << endl;

       E_axis.push_back( hE->GetBinCenter( bin_E ) );
  
       TH1D* hJER_Slice = (TH1D*)hJER->ProjectionY("slice", bin_E, bin_E,"");
       
       TF1* f1 = new TF1("f1", "gaus",  fitlow, fithigh);
       hJER_Slice->Fit("f1", "SN0", "", fitlow, fithigh);	// Store fit result in a pointer but DO NOT DRAW!
       double mu = f1->GetParameter(1);		mu_.push_back( mu );
       double sigma = f1->GetParameter(2);   	sigma_.push_back( sigma );
     }
   
     /* Convert vectors to arrays and create graphs. */
     double * E_axis_arr= &E_axis[0],
   	    * mu_arr 	= &mu_[0],
	    * sigma_arr	= &sigma_[0],
  	    * chi2_arr	= &chi2_[0];
  	    
     cout << "graph's" << endl;
	  
     TGraph *sigma_plot 	= new TGraph( E_axis.size(), E_axis_arr, sigma_arr );
     TGraph *mu_plot 		= new TGraph( E_axis.size(), E_axis_arr, mu_arr );
   
     histo = (TH1D*)hE->Clone(legend);
   
     for(int bin_E = 0; bin_E < E_axis.size(); bin_E++){
       histo->SetBinContent( bin_E, mu_arr[bin_E]);
     }
   }
}


// --------------------------------------------------------------------------------------

/**********************************
* Plot Bayesian unfolding vs. Gen *
**********************************/  

void Plot_Bayes_vs_Gen(TH1D* &histo, TString filename , TString legend){

  TH1D* hGen;
  	Plot_GenLevel(hGen, filename, legend, "");
  TH1D* hBayes;
  	PlotCorrectionFactors_Bayesian( hBayes, filename, legend, "");
	
  histo = (TH1D*)hBayes->Clone(legend);
  histo->Divide(hGen);
} 


// --------------------------------------------------------------------------------------

/************************************
* Plot Bin-by-bin unfolding vs. Gen *
************************************/  

void Plot_BinByBin_vs_Gen(TH1D* &histo, TString filename , TString legend){

  TH1D* hGen;
  	Plot_GenLevel(hGen, filename, legend, "");
  TH1D* hUnfold;
  	PlotCorrectionFactors_BinByBin_Unfolded( hUnfold, filename, legend, "");
	
  histo = (TH1D*)hUnfold->Clone(legend);
  histo->Divide(hGen);
}



// --------------------------------------------------------------------------------------

/*********************
* Determine function *
*********************/  

void Determine_function( void(*current_function)(TH1D* &, TString, TString, TString), TString func_name ){

  if( func_name == "Plot_GenLevel"){
    current_function = &Plot_GenLevel;
  }
  else if(func_name == "Plot_DetLevel"){
    current_function = &Plot_DetLevel;    
  }  
}


