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
#include "Function_average_histogram.h"
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
    TString calib_ = "";
/*
   filenames.push_back("LoopRootFiles/20150317_calibrated_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.400000_12137851_0_sectors.root");
    legendEntries.push_back("E_{I} = 0, #eta-band = 0.4, calibrated");
    colours.push_back(getColor(1)); 
    TString calib_  = "_calib";
*/
     /* Open files. */

     for(int file = 0; file < filenames.size(); file++){
       cout << "$$$ file $$$\t" << file << endl;

       TFile *_file0 = TFile::Open( filenames[file], "Read");
       
       TCanvas *c = new TCanvas("can", "can",  1000, 1000);
        c->SetLeftMargin(0.20);
        c->SetRightMargin(0.18);
        c->SetBottomMargin(0.20);  

       TH2D* hResponse 	= (TH2D*)_file0->Get("hResponse");			hResponse->Draw("colz");
       		c->SetLogz();	
       		c->SaveAs("hResponse" + calib_ + ".pdf");
		c->SaveAs("hResponse" + calib_ + ".C");
       
       TH2D* hDetE 	= (TH2D*)_file0->Get("hCastorJet_energy_response");	hDetE->Draw("colz");	
       		c->SetLogz();
       		c->SaveAs("hDetE" + calib_ + ".pdf");
		c->SaveAs("hDetE" + calib_ + ".C");
       TH2D* hGenE	= (TH2D*)_file0->Get("hGen_fine");			hGenE->Draw("colz");	
       		c->SetLogz();
       		c->SaveAs("hGen_fine" + calib_ + ".pdf");
		c->SaveAs("hGen_fine" + calib_ + ".C");
       
       TH1D* hGen	= (TH1D*)_file0->Get("hGenJet_energy");		
       
       vector<double> mean_response;
       vector<double> mean_response_inverse;
       vector<double> mean_eDet;
       vector<double> mean_eGen;
       vector<double> eGen_center;
       
       cout << "Bin \tResponse \t1/Response \tEdet \tEgen" << endl;
       
       int valid_fits = 0;
       
       for(int bin_E = 0; bin_E < hResponse->GetNbinsX(); bin_E++){
         
	 // Get 1D projections.
         TH1D* hResponse_1D 	= (TH1D*)hResponse->ProjectionY(TString::Format("Response_1D_bin_%i", bin_E), bin_E, bin_E, "do");
	 TH1D* hDetE_1D		= (TH1D*)hDetE->ProjectionX(TString::Format("eDet_1D_bin_%i", bin_E), bin_E, bin_E, "do");
	 TH1D* hGenE_1D		= (TH1D*)hGenE->ProjectionY(TString::Format("eGen_1D_bin_%i", bin_E), bin_E, bin_E, "do");
	 
	 // Create gaussian.
	 TF1 *fit_response	= new TF1(TString::Format("Fit_Response_1D_bin_%i", bin_E), "gaus");
	 
	 // Fit Gaussian.
	 hResponse_1D	->Fit( fit_response , "SQ0");	 
	 
	 // Extract mean values.
	 if( fit_response->GetParameter(1) != 0. ){
	   mean_response	.push_back( fit_response->GetParameter(1) );
	   mean_response_inverse.push_back( 1./fit_response->GetParameter(1) );
	   mean_eDet		.push_back( GetAverage(hDetE_1D) );
	   mean_eGen		.push_back( GetAverage(hGenE_1D) );
	   eGen_center		.push_back( hGen->GetBinCenter( bin_E ) );
	   	   
	   cout << hGen->GetBinLowEdge(bin_E) << "\t" << hGen->GetBinLowEdge(bin_E + 1) << "\t" << mean_response[valid_fits] << "\t" << mean_response_inverse[valid_fits] << "\t" << mean_eDet[valid_fits] << "\t" << mean_eGen[valid_fits] << endl << endl;

	   valid_fits++;
	   
	 }
       }
       
       for( int bin = 0; bin < mean_eDet.size(); bin++){
         cout << "lowedge.push_back(\t" << mean_eDet[bin] << "\t); muval.push_back(\t" << mean_response_inverse[bin] << ");" << endl;
       }

       // Transform the vectors into arrays.       
       double * E_gen_axis	= &mean_eGen[0],
       	      * E_det_axis	= &mean_eDet[0],
	      * response_val	= &mean_response[0],
	      * response_inverse= &mean_response_inverse[0],
	      * eGen_bin	= &eGen_center[0];
	      
       // Create graphs from the arrays.
       // -- E gen versus response
       TGraph * Response_true = new TGraph( mean_eGen.size(), E_gen_axis, response_val);
       	 Response_true	->GetXaxis()->SetTitle("E_{gen}");
	 Response_true	->GetYaxis()->SetTitle("< #frac{E_{det}}{E_{gen}} >");
	 Response_true	->GetYaxis()->SetRangeUser(0., 2.5);
	 
       // -- E det versus 1/response
       TGraph * Response_meas = new TGraph( mean_eDet.size(), E_det_axis, response_inverse);
       	 Response_meas	->GetXaxis()->SetTitle("E_{det}");
	 Response_meas	->GetYaxis()->SetTitle("< #frac{E_{gen}}{E_{det}} >");  
	 Response_meas	->GetYaxis()->SetRangeUser(0.,2.5);      
 	
       // -- E gen versus average Egen.
       TGraph * Gen_average = new TGraph( eGen_center.size(), eGen_bin, E_gen_axis);
         Gen_average	->GetXaxis()->SetTitle("E_{gen}");
	 Gen_average	->GetYaxis()->SetTitle("<E>");
	 
       // -- E gen versus average Edet.
       TGraph * Det_average = new TGraph( eGen_center.size(), eGen_bin, E_det_axis);
         Det_average	->GetXaxis()->SetTitle("E_{gen}");
	 Det_average	->GetYaxis()->SetTitle("<E>");	
 
 
       // Draw.
       TCanvas *can_true = new TCanvas("can_true" + calib_, "can_true" + calib_, 1.);	Response_true->Draw("A*");
       TCanvas *can_meas = new TCanvas("can_meas" + calib_, "can_meas" + calib_, 1.);	Response_meas->Draw("A*");
       TCanvas *can_aver = new TCanvas("can_aver" + calib_, "can_aver" + calib_, 1000, 1000);
        	can_aver->SetLeftMargin(0.20);
        	can_aver->SetRightMargin(0.18);
        	can_aver->SetBottomMargin(0.20);
		  
       		Gen_average->Draw("A*");
		
		Det_average->SetMarkerColor( kRed + 2);
		Det_average->SetMarkerStyle(22);
		Det_average->Draw("psame");
       
       // Save.
       can_true->SaveAs("Test_true" + calib_ + ".C");	can_true->SaveAs("Test_true" + calib_ + ".pdf");
       can_meas->SaveAs("Test_meas" + calib_ + ".C"); 	can_meas->SaveAs("Test_meas" + calib_ + ".pdf");
       can_aver->SaveAs("Test_aver" + calib_ + ".C");	can_aver->SaveAs("Test_aver" + calib_ + ".pdf");
              
     } // Loop over files.
 return (0);
}
