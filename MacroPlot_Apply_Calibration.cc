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
#include <TGraphErrors.h>
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

   int rebinner = 4;

   gROOT->SetStyle ("Plain");
   gStyle->SetPalette(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   // Step 1 - open Root files and retrieve 2D histograms.
   
   TLegend *legend_det = new TLegend(0.30, 0.15, 0.80, 0.38);
     legend_det->SetFillStyle( 0 );
     legend_det->SetBorderSize( 0 );
     
   TLegend *legend_gen = new TLegend(0.30, 0.70, 0.80, 0.93);
     legend_gen->SetFillStyle( 0 );
     legend_gen->SetBorderSize( 0 );   

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
  int event_threshold_ = 200.;

  TString label = "20150330_Plots_Unfolding";

   /*******************************************
   * Create a new subdirectory for the plots. *
   ********************************************/

   int 	new_dir = mkdir( ("Plots/" + label).Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  );

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
   TString calib_ = "_compare";
   TString drawoptions = "A";
   int color_index = 0;
  
/*
   filenames.push_back("LoopRootFiles/20150330_ak5ak5_JetAnalyzer_etaband_0.400000_12137851_0_sectors_had.root");
    legendEntries.push_back("E_{I} = 0, #eta-band = 0.4");
    colours.push_back(getColor(color_index++)); 
    calib_ = "";
*/
/*
   filenames.push_back("LoopRootFiles/20150330_ak5ak5_calibrated_piecewise_logform_to_slope_JetAnalyzer_etaband_0.400000_12137851_0_sectors_had.root");
    legendEntries.push_back("E_{I} = 0, #eta-band = 0.4, calibrated (ax + c)");
    colours.push_back(getColor(color_index++)); 
    calib_  = "_calib_f2";
*/
/*
 * //-- These two files work as calibration and test of calibration.  
   filenames.push_back("LoopRootFiles/AprilVis_JetAnalyzer_etaband_0.400000_12137851_0_sectors_had.root");
    legendEntries.push_back("E_{I} = 0, #eta-band = 0.4");
    colours.push_back(getColor(color_index++)); 
    calib_ = "";

   filenames.push_back("LoopRootFiles/20150402_calibrated_function_JetAnalyzer_etaband_0.400000_12137851_0_sectors_had.root");
    legendEntries.push_back("E_{I} = 0, #eta-band = 0.4");
    colours.push_back(getColor(color_index++)); 
    calib_ = "_calib_f2";
*/

//20150410_JetAnalyzer_etaband_0.400000_12137851_0_sectors_had.root

   filenames.push_back("Histograms/ak5ak5_EIcut_JetAnalyzer_etaband_0.400000_12137851_0_sectors_had.root");
    legendEntries.push_back("E_{I} = 0, #eta-band = 0.4");
    colours.push_back(getColor(color_index++)); 
    calib_ = "_calib_f2";

    if( filenames.size() > 1 ){ 
    	calib_ = "comparison";
	draw_legend = true;
    }

    TCanvas *can_true = new TCanvas("can_true" + calib_, "can_true" + calib_, 1000, 700);
      can_true	->SetLeftMargin(0.16);
      can_true	->SetTopMargin(0.07);
      can_true	->SetBottomMargin(0.20);	 
	 
    TCanvas *can_meas = new TCanvas("can_meas" + calib_, "can_meas" + calib_, 1000, 700);
      can_meas	->SetLeftMargin(0.16);
      can_meas	->SetTopMargin(0.07);
      can_meas	->SetBottomMargin(0.20); 
      
    TCanvas *can_fit = new TCanvas("can_fit" + calib_, "can_fit" + calib_, 1000, 700);
      can_fit	->SetLeftMargin(0.16);
      can_fit	->SetTopMargin(0.07);
      can_fit	->SetBottomMargin(0.20);         
	 
    TCanvas *can_aver = new TCanvas("can_aver" + calib_, "can_aver" + calib_, 1000, 700);
      can_aver	->SetLeftMargin(0.16);
      can_aver	->SetRightMargin(0.07);
      can_aver	->SetBottomMargin(0.20);
      
    TLine *line = new TLine(0.,1.,1733.,1.);
    line->SetLineWidth(2);

     /* Open files. */

     for(int file = 0; file < filenames.size(); file++){
       cout << "$$$ file $$$\t" << file << endl;

       TFile *_file0 = TFile::Open( filenames[file], "Read");
       
       TCanvas *c = new TCanvas("can", "can",  1000, 1000);
        c->SetLeftMargin(0.20);
        c->SetRightMargin(0.18);
        c->SetBottomMargin(0.20);  

       TH2D* hResponse 	= (TH2D*)_file0->Get("hResponse_gen");			hResponse->Draw("colz");
		hResponse->RebinX( rebinner );
		hResponse->RebinY( rebinner );
       		c->SetLogz();	
       		c->SaveAs("Plots/" + label +"/hResponse" + calib_ + ".pdf");
		c->SaveAs("Plots/" + label +"/hResponse" + calib_ + ".C");

/*
       for(int nx = 0; nx < hResponse->GetNbinsX(); nx++){
         for( int ny = 0; ny < hResponse->GetNbinsY(); ny++){
	
	  if( hResponse->GetBinContent( nx, ny ) > 0.) { cout << "(" << nx << ", " << ny << ")\t" << hResponse->GetBinContent( nx, ny ) << "\t+/-\t" << hResponse->GetBinError( nx, ny ) << endl; }

	 }
       }
*/       
       TH2D* hDetE 	= (TH2D*)_file0->Get("hCastorJet_energy_response");	hDetE->Draw("colz");	
		hDetE->RebinX( rebinner );
		hDetE->RebinY( rebinner );
       		c->SetLogz();
       		c->SaveAs("Plots/" + label +"/hDetE" + calib_ + ".pdf");
		c->SaveAs("Plots/" + label +"/hDetE" + calib_ + ".C");

       TH2D* hGenE	= (TH2D*)_file0->Get("hGen_fine");			hGenE->Draw("colz");	
                hGenE->RebinX( rebinner );
                hGenE->RebinY( rebinner );

       		c->SetLogz();
       		c->SaveAs("Plots/" + label +"/hGen_fine" + calib_ + ".pdf");
		c->SaveAs("Plots/" + label +"/hGen_fine" + calib_ + ".C");
       
       TH1D* hGen	= (TH1D*)_file0->Get("hGenJet_energy");		

// cout << "// -- BINS COMPARISON\t" << hGenE->GetNbinsX() << "\t" << hGenE->GetXaxis()->GetTitle() << "\t" << hGenE->GetNbinsY() << "\t" << hGenE->GetYaxis()->GetTitle() << "\t"  << hDetE->GetNbinsX() << "\t" << hDetE->GetXaxis()->GetTitle() << "\t" <<
//hDetE->GetNbinsY() << "\t" << hDetE->GetYaxis()->GetTitle() << "\t"  << hResponse->GetNbinsX() << "\t" << hResponse->GetXaxis()->GetTitle() << "\t" << hResponse->GetNbinsY() << "\t" << hResponse->GetYaxis()->GetTitle() << endl;	

       
       vector<double> mean_response;
       vector<double> mean_response_inverse;
       vector<double> mean_eDet;
       vector<double> error_eDet;
       vector<double> mean_eGen;
       vector<double> error_eGen;
       vector<double> eGen_center;
       vector<double> error_energy;
       
       cout << "Bin \tResponse \t1/Response \tEdet \tEgen" << endl;
       
       int valid_fits = 0;
       
       // Create a histogram to store slices with low statistics.
       TH1D* hSlice_storage	= (TH1D*)hResponse->ProjectionY("Storage", 1, 1, "do");
         for(int bins_store = 0; bins_store <= hSlice_storage->GetNbinsX(); bins_store++){ hSlice_storage->SetBinContent(bins_store, 0.); }

       // Do the same for the energies.
       TH1D* hEgen_storage = (TH1D*)hGenE->ProjectionY("Storage_gen", 0, 0, "do"); 
	 for(int bins_store = 0; bins_store <= hEgen_storage->GetNbinsX(); bins_store++){ hEgen_storage->SetBinContent(bins_store, 0.); }
       TH1D* hEdet_storage = (TH1D*)hDetE->ProjectionX("Storage_det", 0, 0, "do");
	 for(int bins_store = 0; bins_store <= hEdet_storage->GetNbinsX(); bins_store++){ hEdet_storage->SetBinContent(bins_store, 0.); }


       /******************
       * Loop over bins. *       
       ******************/                
      
             
       for(int bin_E = 0; bin_E < hResponse->GetNbinsX(); bin_E++){
         
	 // Get 1D projections.
         TH1D* hResponse_1D 	= (TH1D*)hResponse->ProjectionY(TString::Format("Response_1D_bin_%i", bin_E), bin_E, bin_E, "do");
	 TH1D* hDetE_1D		= (TH1D*)hDetE->ProjectionX(TString::Format("eDet_1D_bin_%i", bin_E), bin_E, bin_E, "do");
	 TH1D* hGenE_1D		= (TH1D*)hGenE->ProjectionY(TString::Format("eGen_1D_bin_%i", bin_E), bin_E, bin_E, "do");
	 
	 // Create gaussian.
	 TF1 *fit_response	= new TF1(TString::Format("Fit_Response_1D_bin_%i", bin_E), "gaus");

	 if( (hResponse_1D->Integral() + hSlice_storage->Integral() )< event_threshold_ ){
	   hSlice_storage	->Add( hResponse_1D );
	   hEdet_storage	->Add( hDetE_1D );
	   hEgen_storage	->Add( hGenE_1D );
	   continue;
	 }

         else if( ( hResponse_1D->Integral() + hSlice_storage->Integral() ) >= event_threshold_ ){
           hResponse_1D ->Add( hSlice_storage );
           hDetE_1D     ->Add( hEdet_storage );
           hGenE_1D     ->Add( hEgen_storage );
         }

	 for(int bins_store = 0; bins_store <= hSlice_storage->GetNbinsX(); bins_store++){ hSlice_storage->SetBinContent(bins_store, 0.); 	hSlice_storage->SetBinError(bins_store, 0.);}
	 for(int bins_store = 0; bins_store <= hEgen_storage->GetNbinsX(); bins_store++){ hEgen_storage->SetBinContent(bins_store, 0.); 	hEgen_storage->SetBinError(bins_store, 0.);}
	 for(int bins_store = 0; bins_store <= hEdet_storage->GetNbinsX(); bins_store++){ hEdet_storage->SetBinContent(bins_store, 0.); 	hEdet_storage->SetBinError(bins_store, 0.);}

	 TCanvas* can_slice = new TCanvas( TString::Format("Slice_%i", bin_E), TString::Format("Slice_%i", bin_E), 1);

	 // -- Fit Gaussian.
	 hResponse_1D	->Draw();
	 
	 if( file == 0 ){
  	   hResponse_1D	->Fit( fit_response , "SQ");	 
	   can_slice->SaveAs(TString::Format("Plots/" + label + "/Slice_%i.C", bin_E) );
           can_slice->SaveAs(TString::Format("Plots/" + label + "/Slice_%i.pdf", bin_E) );
	}
	else{
	  hResponse_1D ->Fit( fit_response , "SQ0");
	}
	 
	 // -- Extract mean values.
	 if( fit_response->GetParameter(1) != 0. ){
	   double R = fit_response->GetParameter(1);
	   double sR= fit_response->GetParError(1);	

	   mean_response	.push_back( R );
	   mean_response_inverse.push_back( 1./R );
	   mean_eDet		.push_back( GetAverage(hDetE_1D) );
	   mean_eGen		.push_back( GetAverage(hGenE_1D) );
	   eGen_center		.push_back( hGen->GetBinCenter( bin_E ) );

	   if( R == 0.){ 	error_eGen	.push_back( 0. ); }
	   else{ 		error_eGen	.push_back( sR ); }
           if( R == 0.){	error_eDet	.push_back( 0. ); }
	   else{ 		error_eDet	.push_back( sR/(R*R) ); }
	   error_energy		.push_back( 0. );	   
	   
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
	      * eGen_bin	= &eGen_center[0],
	      * err_R		= &error_eGen[0],
	      * err_1overR	= &error_eDet[0],
	      * err_energy	= &error_energy[0];
       // Create graphs from the arrays.
       // -- E gen versus response
       TGraphErrors * Response_true = new TGraphErrors( mean_eGen.size(), E_gen_axis, response_val, err_energy, err_R);
       	 Response_true	->GetXaxis()->SetTitle("E_{gen}");
	 Response_true	->GetXaxis()->SetTitleOffset(0.9);
	 Response_true	->GetXaxis()->SetTitleSize(0.08);
	 Response_true	->GetXaxis()->SetLabelSize(0.07);
	 Response_true	->GetXaxis()->SetNdivisions(505);
	 	 
	 Response_true	->GetYaxis()->SetTitle("< #frac{E_{det}}{E_{gen}} >");
	 Response_true	->GetYaxis()->SetTitleOffset(0.75);
	 Response_true	->GetYaxis()->SetTitleSize(0.08);
	 Response_true	->GetYaxis()->SetLabelSize(0.08);	 
	 Response_true	->GetYaxis()->SetRangeUser(0., 2.5);
	 
	 Response_true	->SetMarkerColor( colours[file] );
	 Response_true	->SetLineColor( colours[file] );
	 Response_true	->SetMarkerStyle( 20 + file );
	 
       // -- E det versus 1/response
       TH1F* hEdet_axis = new TH1F("hEdet_axis", "hEdet_axis", 100, 0, mean_eGen[ mean_eGen.size()-1] );
         hEdet_axis->GetYaxis()->SetRangeUser(0., 2.5);
       	 hEdet_axis	->GetXaxis()->SetTitle("E_{det}");
	 hEdet_axis	->GetXaxis()->SetTitleOffset(0.9);
	 hEdet_axis	->GetXaxis()->SetTitleSize(0.08);
	 hEdet_axis	->GetXaxis()->SetLabelSize(0.08);
	 hEdet_axis	->GetXaxis()->SetNdivisions(505);
	 	 
	 hEdet_axis	->GetYaxis()->SetTitle("< #frac{E_{gen}}{E_{det}} >");  
	 hEdet_axis	->GetYaxis()->SetRangeUser(0.,2.5);      
	 hEdet_axis	->GetYaxis()->SetTitleOffset(0.75);
	 hEdet_axis	->GetYaxis()->SetTitleSize(0.08);
	 hEdet_axis	->GetYaxis()->SetLabelSize(0.08); 
	       
       TGraphErrors * Response_meas = new TGraphErrors( mean_eDet.size(), E_det_axis, response_inverse,  err_energy, err_1overR);
       	 Response_meas	->GetXaxis()->SetTitle("E_{det}");
	 Response_meas	->GetXaxis()->SetTitleOffset(0.9);
	 Response_meas	->GetXaxis()->SetTitleSize(0.08);
	 Response_meas	->GetXaxis()->SetLabelSize(0.05);
	 Response_meas	->GetXaxis()->SetNdivisions(505);
	 	 
	 Response_meas	->GetYaxis()->SetTitle("< #frac{E_{gen}}{E_{det}} >");  
	 Response_meas	->GetYaxis()->SetRangeUser(0.,2.5);      
	 Response_meas	->GetYaxis()->SetTitleOffset(1.15);
	 Response_meas	->GetYaxis()->SetTitleSize(0.06);
	 Response_meas	->GetYaxis()->SetLabelSize(0.07);
	 
	 
	 Response_meas	->SetMarkerColor( colours[file] );
	 Response_meas	->SetLineColor( colours[file] );
	 Response_meas	->SetMarkerStyle( 20 + file );
	 
		  	
       // -- E gen versus average Egen.
       TGraph * Gen_average = new TGraph( eGen_center.size(), eGen_bin, E_gen_axis);
         Gen_average	->GetXaxis()->SetTitle("E_{gen}");
	 Gen_average	->GetYaxis()->SetTitle("<E>");
	 
       // -- E gen versus average Edet.
       TGraph * Det_average = new TGraph( eGen_center.size(), eGen_bin, E_det_axis);
         Det_average	->GetXaxis()->SetTitle("E_{gen}");
	 Det_average	->GetYaxis()->SetTitle("<E>");	
 
       legend_det->AddEntry( Response_meas, legendEntries[file], "lp");
       legend_gen->AddEntry( Response_meas, legendEntries[file], "lp");
 
       // Draw.

	 
       can_true->cd();	 
       Response_true->Draw("p" + drawoptions);
       line->Draw();

       can_meas->cd(); 
       if( file == 0 ){ hEdet_axis->Draw(); }    
       Response_meas->Draw("psame");
       line->Draw();
	 	
       can_aver->cd();		  
       Gen_average->Draw("p" + drawoptions);
		
	Det_average->SetMarkerColor( kRed + 2);
	Det_average->SetMarkerStyle(22);
	Det_average->Draw("psame");
       
       /*******************************************
       * Fit an analytical function to the graph. *
       *******************************************/

       if( filenames.size() == 1 || file == 0){
         can_fit->cd();
	 hEdet_axis->Draw();
	 Response_meas->Draw("pe"); 

	 double e_frac_1 = 60.;
	 double e_frac_2 = 200.;
 
         TF1* analytical = new TF1("Analytical", "( x < 60.) * ([5]*x*x + [6]*x + [7] ) + (x > 60. && x < 200.) * ( [0] + [1] * log( [2] + x) ) + (x > 200. )* ( [3] + [4] * x) ", 0., 900.);
//         TF1* analytical = new TF1("Analytical", "( x < 200.) * ( [0] + [1] * log( [2] + x) ) + (x > 200. )* ( [3] + [4] * x) ", 0., 900.);
	 analytical->SetLineColor( kRed );
         analytical->SetLineWidth( 2 );
         Response_meas->Fit( analytical, "S", "", 0., 900.);
       }
       
       drawoptions = "same";
     } // Loop over files.
     
       
     // Save.
     can_true->cd();
     legend_gen->Draw();
     can_true->SaveAs("Plots/" + label +"/CalibrationFactors_true_" + calib_ + ".C");		can_true->SaveAs("Plots/" + label +"/CalibrationFactors_true_" + calib_ + ".pdf");
     
     can_meas->cd();
     legend_det->Draw();
     can_meas->SaveAs("Plots/" + label +"/CalibrationFactors_meas_" + calib_ + ".C"); 		can_meas->SaveAs("Plots/" + label +"/CalibrationFactors_meas_" + calib_ + ".pdf");

     can_aver->cd();
     legend_gen->Draw();
     can_aver->SaveAs("Plots/" + label +"/AverageEnergyValues_per_EgenBin_" + calib_ + ".C");	can_aver->SaveAs("Plots/" + label +"/AverageEnergyValues_per_EgenBin_" + calib_ + ".pdf");  
        
    can_fit->SaveAs("Plots/" + label +"/CalibrationFactors_meas_with_fit" + calib_ + ".C");             can_fit->SaveAs("Plots/" + label +"/CalibrationFactors_meas_with_fit" + calib_ + ".pdf");
     
 return (0);
}
