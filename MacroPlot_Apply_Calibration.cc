/*
 *The goal of this macro is to calculate the mean and spread of JER (E_gen). 
 *
 * 1. Retrieve (E_gen, JER) 2D plot from the root files.
 * 2. Fit a Gaussian to slices of constant E_gen.
 * 3. Plot the mean and spread of these fits as functions of E_gen.
 *
 * 4. Plot ratios of hadronic to electromagnetic detector jets and spread.
 */

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
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TGraph.h>

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

#include "color.h"
#include "NonZeroMinimum.h"
using namespace std;

int main(){

   gROOT->SetStyle ("Plain");
   gStyle->SetPalette(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   // Step 1 - open Root files and retrieve 2D histograms.
   
   TLegend *legend = new TLegend(0.75, 0.60, 0.99, 0.99);
     legend->SetFillColor( kWhite );

//   TString date = "20150125_had_em_jets"; 
//   TString date = "20150205_no_gen_cut";
//   TString date = "20150225_iso-cut";
   TString date = "20150305_had";

//   TString numb = "10000000";  TString suffix = "_173_1_iQO.root";
   TString numb = "12137851";
//   TString numb = "1000000"; TString suffix = "_13_1_irJ.root"; 
   //  TString suffix = "_29_1_WJW.root";	// Hardest gen jet matched to closest (phi) det jet.
   // TString numb = "100000"; TString suffix = "_129_1_D3d.root";

  int Slice_threshold = 0.;

   double E_gen_cut = 0.;
   TString distr = "gaus";

   std::vector<TString> filenames;
   std::vector<TString> legendEntries;
   std::vector<TString> plotVariables;
   std::vector<TString> jetSelection;
     std::vector<TString> plotX;
     std::vector<TString> plotY;
     std::vector<TString> plotTitle;
   std::vector<int>     colours;
   std::vector<double> plot_min;
   std::vector<TString> can_suffix;
   
   bool draw_legend = false;
 
   filenames.push_back("LoopRootFiles/20150317_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.400000_12137851_0_sectors.root");
    legendEntries.push_back("E_{I} = 0, #eta-band = 0.4");
    colours.push_back(getColor(1)); 

     /* Open files. */

     for(int file = 0; file < filenames.size(); file++){
       cout << "$$$ file $$$\t" << file << endl;

       TFile *_file0 = TFile::Open( filenames[file], "Read");

       TH2D* hResponse 	= (TH2D*)_file0->Get("hResponse");
       TH2D* hDetE 	= (TH2D*)_file0->Get("hCastorJet_energy_response");
       TH2D* hGenE	= (TH2D*)_file0->Get("hGen_fine");
       
       vector<double> mean_response;
       vector<double> mean_response_inverse;
       vector<double> mean_eDet;
       vector<double> mean_eGen;
       
       for(int bin_E = 0; bin_E < hResponse->GetNbinsX(); bin_E++){
         
	 // Get 1D projections.
         TH1D* hResponse_1D 	= (TH1D*)hResponse->ProjectionX(TString::Format("Response_1D_bin_%i", bin_E), bin_E, bin_E, "do");
	 TH1D* hDetE_1D		= (TH1D*)hDetE->ProjectionX(TString::Format("eDet_1D_bin_%i", bin_E), bin_E, bin_E, "do");
	 TH1D* hGenE_1D		= (TH1D*)hGenE->ProjectionX(TString::Format("eGen_1D_bin_%i", bin_E), bin_E, bin_E, "do");
	 
	 // Create gaussian.
	 TF1 *fit_response	= new TF1(TString::Format("Fit_Response_1D_bin_%i", bin_E), "gaus");
	 TF1 *fit_eDet		= new TF1(TString::Format("Fit_eDet_1D_bin_%i", bin_E), "gaus");
	 TF1 *fit_eGen		= new TF1(TString::Format("Fit_eGen_1D_bin_%i", bin_E), "gaus", hGenE_1D->GetBinLowEdge(bin_E), hGenE_1D->GetBinLowEdge(bin_E+1));
	 
	 // Fit Gaussian.
	 hResponse_1D	->Fit( fit_response , "S0");
	 hDetE_1D	->Fit( fit_eDet , "S0");
	 hGenE_1D	->Fit( fit_eGen , "S0");
	 
	 // Extract mean values.
	 mean_response		.push_back( fit_response->GetParameter(1) );
	 mean_response_inverse	.push_back( 1./fit_response->GetParameter(1) );	 
	 mean_eDet		.push_back( fit_eDet->GetParameter(1) );
	 mean_eGen		.push_back( fit_eGen->GetParameter(1) );	 
       }

       // Transform the vectors into arrays.       
       double * E_gen_axis	= &mean_eGen[0],
       	      * E_det_axis	= &mean_eDet[0],
	      * response_val	= &mean_response[0],
	      * response_inverse= &mean_response_inverse[0];
	      
       // Create graphs from the arrays.
       TGraph * Response_true = new TGraph( mean_eGen.size(), E_gen_axis, response_val);
       TGraph * Response_meas = new TGraph( mean_eDet.size(), E_det_axis, response_inverse);
 
       // Draw.
       TCanvas *can_true = new TCanvas("can_true", "can_true", 1.);	Response_true->Draw("AC*");
       TCanvas *can_meas = new TCanvas("can_meas", "can_meas", 1.);	Response_meas->Draw("AC*");
       
       // Save.
       can_true->SaveAs("Test_true.C");	can_true->SaveAs("Test_true.pdf");
       can_meas->SaveAs("Test_meas.C"); can_meas->SaveAs("Test_meas.pdf");
              
     } // Loop over files.
 return (0);
}
