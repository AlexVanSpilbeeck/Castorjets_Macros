//STANDARD ROOT INCLUDES
#include <TROOT.h>
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

   RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response");

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
 
   /* Extract GEN and DET (MC) and DATA distribution */

   RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response");

   TH1D *hDet = (TH1D*)_file0->Get("hCastorJet_energy");	hDet->Sumw2();
   TH1D *hGen = (TH1D*)_file0->Get("hGenJet_energy");	hGen->Sumw2();

   int nbins_x = hGen->GetXaxis()->GetNbins();
   double x_axisMax = hGen->GetXaxis()->GetBinUpEdge( nbins_x );
  
   /****************/
   /** Unfold MC. **/
   /****************/

   RooUnfoldBayes unfold_mc_bayes (response, hDet, 10);
     TH1D* hReco_mc_bayes= (TH1D*) unfold_mc_bayes.Hreco();

    hReco_mc_bayes->SetLineStyle(3);
     hReco_mc_bayes->SetLineColor(kGreen);
     hReco_mc_bayes->SetLineWidth(2);
     hReco_mc_bayes->SetName("hReco_bayes_1");
     
   histo = (TH1D*)hReco_mc_bayes->Clone(legend);

  /*****************************/
  /** Get correction factors. **/
  /*****************************/
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
 
   /* Extract GEN and DET (MC) and DATA distribution */

   RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response");

   TH1D *hDet = (TH1D*)_file0->Get("hCastorJet_energy");	hDet->Sumw2();
   histo = (TH1D*)hDet->Clone(legend);  
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
 
   /* Extract GEN and DET (MC) and DATA distribution */

   RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response");

   TH1D *hGen = (TH1D*)_file0->Get("hGenJet_energy");	hGen->Sumw2();
   histo = (TH1D*)hGen->Clone(legend);  
 } 
  
