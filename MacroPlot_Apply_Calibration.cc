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

/*
     filenames.push_back("LoopRootFiles/" + date + "_000_had_Output_JetAnalyzer_radii_strippedTree_GEN_ak5_DET_ak5_margin_0.500000_" + numb + "_0_sectors.root");
     legendEntries.push_back("E_{I} = 0.0 ");
     colours.push_back(getColor(8));

     filenames.push_back("LoopRootFiles/" + date + "_001_had_Output_JetAnalyzer_radii_strippedTree_GEN_ak5_DET_ak5_margin_0.500000_" + numb + "_0_sectors.root");
     legendEntries.push_back("E_{I} < 0.01 E_{det}");
     colours.push_back(getColor(4));

     filenames.push_back("LoopRootFiles/" + date + "_01_had_Output_JetAnalyzer_radii_strippedTree_GEN_ak5_DET_ak5_margin_0.500000_" + numb + "_0_sectors.root");
     legendEntries.push_back("E_{I} < 0.1 E_{det}");
     colours.push_back(getColor(7));

     filenames.push_back("LoopRootFiles/" + date + "_05_had_Output_JetAnalyzer_radii_strippedTree_GEN_ak5_DET_ak5_margin_0.500000_" + numb + "_0_sectors.root");
     legendEntries.push_back("E_{I} < 0.5 E_{det}");
     colours.push_back(getColor(5));

     filenames.push_back("LoopRootFiles/" + date + "_1_had_Output_JetAnalyzer_radii_strippedTree_GEN_ak5_DET_ak5_margin_0.500000_" + numb + "_0_sectors.root");
     legendEntries.push_back("E_{I} < 1 E_{det}");
     colours.push_back(getColor(6));
 */
 
     filenames.push_back("LoopRootFiles/20150311_function_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.400000_12137851_0_sectors.root");
     legendEntries.push_back("E_{I} = 0, #eta-band = 0.4");
     colours.push_back(getColor(1)); 

     plotVariables.push_back("hJER_per_energy");
       plotTitle.push_back("#DeltaE/E for all detector jets");
       can_suffix.push_back("all");
       plotX.push_back("E_{gen}");

     plotVariables.push_back("hJER_per_eDet");
       plotTitle.push_back("#DeltaE/E for all detector jets");
       can_suffix.push_back("all");
       plotX.push_back("E_{det}");
/*
     plotVariables.push_back("hJER_per_eGen");
       plotTitle.push_back("#DeltaE/E for all detector jets");
       can_suffix.push_back("all");
       plotX.push_back("E_{gen}");
*/

  	 for( int plot = 0; plot < plotVariables.size(); plot++ ){
     cout << "$$$ plot $$$\t" << plot << endl;
   
     double max_val_sigma, max_val_mu;
	if( plotVariables[plot].Contains("energy") ){ max_val_sigma = -10., max_val_mu = -10.; }
        if( plotVariables[plot].Contains("eDet") ){ max_val_sigma = -10., max_val_mu = 10.; }

     TH1D* originalSpread, *originalMean;
     TString variable = plotVariables[plot];
     

     TCanvas *can_sigma = new TCanvas("can_sigma_" + can_suffix[plot], "can_sigma_" + can_suffix[plot], 1.);

        can_sigma->SetLeftMargin(0.14);
        can_sigma->SetRightMargin(0.03);
        can_sigma->SetBottomMargin(0.14);

     TCanvas *can_mu  = new TCanvas("can_mu_" + can_suffix[plot], "can_mu_" + can_suffix[plot] , 1.);
	can_mu->SetRightMargin(0.03);
	can_mu->SetLeftMargin(0.14);
	can_mu->SetBottomMargin(0.14);

     TCanvas *can_JER  = new TCanvas("can_JER_" + can_suffix[plot], "can_JER_" + can_suffix[plot] , 1.);

     TString drawoptions = "E";


     /* Open files. */

     for(int file = 0; file < filenames.size(); file++){
       cout << "$$$ file $$$\t" << file << endl;

       TFile *_file0 = TFile::Open( filenames[file], "Read");

       TH2D *hJER;
       hJER = (TH2D*)_file0->Get( variable );        
       hJER->GetXaxis()->SetTitle( plotX[plot] );
       hJER->GetYaxis()->SetTitle("#DeltaE/E");
       hJER->SetTitle("#DeltaE/E");
       hJER->RebinX(4);
       
       TH1D *hCastorJet_energy;
       hCastorJet_energy = (TH1D*)_file0->Get( "hCastorJet_energy" );
  
       /* 
       // Now we loop over the bins in E_gen and plot each separate slice.
       // We fit three times and store all values.
       */

       if( variable == "hJER_per_eDet" && file == 0){
	 cout << "$$$ Plotting $$$" << endl;

	 TCanvas *can_2D = new TCanvas(date + "_" + variable, date + "_" + variable ,1.);
	 hJER->GetXaxis()->SetRangeUser(0., 1800.);
         hJER->GetYaxis()->SetRangeUser(-1,5.);
	 can_2D->SetLogz();

	 hJER->Draw("colz");
	 can_2D->SaveAs("Plots/Fits/" + date + "_" + variable + "_2D.pdf");
         can_2D->SaveAs("Plots/Fits/" + date + "_" + variable + "_2D.C");

	//
	// Loop over all slices in our 2D histogram. Each slices represents a DeltaE/E for a fixed energy.
	//
	
	std::vector<double> sigma0, sigma1, sigma2;
	std::vector<double> mu0, mu1, mu2;
	std::vector<double> chi2_0, chi2_1, chi2_2;
	std::vector<double> E_axis;
	std::vector<double> Edet, Egen, Nevents;

	// -- Loop over the 1D energy distribution to fill the Edet and Nevents (x- and y-axes).
	hCastorJet_energy->Scale(1./hCastorJet_energy->Integral() );
	for(int bin_det = 0; bin_det < hCastorJet_energy->GetNbinsX(); bin_det++){
	  Edet.push_back( hCastorJet_energy->GetBinCenter( bin_det ) );
          Egen.push_back( hCastorJet_energy->GetBinCenter( bin_det ) );
	  Nevents.push_back( hCastorJet_energy->GetBinContent( bin_det ) );
	}

	TH1D* hSlice_to_plot = (TH1D*)hJER->ProjectionY("slice_to_plot", 0, 0,"");;	// This will collect slices until it is full enough.

	TH1D* hE_gen = (TH1D*)hJER->ProjectionX("e_gen", 1, -1,"");
	for(int bin_E = 0; bin_E < hJER->GetNbinsX(); bin_E++){
	
	 

          TCanvas *can_Slice = new TCanvas(TString::Format("can_Energy_fraction_Slice_%i", bin_E) + variable , TString::Format("can_Energy_fraction_Slice_%i", bin_E) + variable , 1.);
            can_Slice->SetLeftMargin(0.07);
            can_Slice->SetRightMargin(0.01);
            can_Slice->SetBottomMargin(0.14);

          TCanvas *can_Slice_noFit = new TCanvas(TString::Format("can_Energy_fraction_Slice_%i_noFit", bin_E) + variable , TString::Format("can_Energy_fraction_Slice_%i_noFit", bin_E) + variable , 1.);
            can_Slice_noFit->SetLeftMargin(0.07);
            can_Slice_noFit->SetRightMargin(0.01);
            can_Slice_noFit->SetBottomMargin(0.14);	
	 
  	  TH1D* hJER_Slice = (TH1D*)hJER->ProjectionY("slice", bin_E, bin_E,"do");

	  hSlice_to_plot = hJER_Slice;

	  if( hSlice_to_plot->Integral() < Slice_threshold ){
	    continue;
	  }

	  hSlice_to_plot->Draw("");	  	 
	  TString slice_name = date + "_" + variable + "_2D_slice" + TString::Format("%i", bin_E); 
	  can_Slice_noFit->SaveAs( "Plots/Fits" + slice_name + ".pdf");
   
          if( hJER_Slice->Integral() > 0){
	  
 	    if( variable == "hJER_per_energy") { hSlice_to_plot->SetMarkerColor(kRed); }
	    else{ 				 hSlice_to_plot->SetMarkerColor(kBlue); }
	    hSlice_to_plot->GetXaxis()->SetRangeUser(-1., 5.);
	    hSlice_to_plot->SetMarkerColor( kRed );

	    gStyle->SetOptTitle(1);
	    double e_gen_val_min = hE_gen->GetBinLowEdge( bin_E );
            double e_gen_val_max = hE_gen->GetBinLowEdge( bin_E+1 );
	    E_axis.push_back( hE_gen->GetBinCenter( bin_E ) );

	    hSlice_to_plot->SetTitle(TString::Format("E_{gen} = %f to %f GeV", e_gen_val_min, e_gen_val_max) );
	        
	    /*
	    can_Slice_noFit->cd();
	      hSlice_to_plot->Draw( "ep" );
	      can_Slice_noFit->SaveAs("Plots/Fits/" + date + TString::Format("_can_Energy_fraction_Slice_%i_noFit_", bin_E) +  variable + ".pdf" );
              can_Slice_noFit->SaveAs("Plots/Fits/" + date + TString::Format("_can_Energy_fraction_Slice_%i_noFit_", bin_E) +  variable + ".C" );
	    */

	    can_Slice->cd();

	    //
	    // Determine range within which the fit is made.
	    //
	    
	    double fitlow, fithigh;
	    if( variable == "hJER_per_energy" ){
	      fitlow = -1;
	      fithigh = 2.5;
	    } // hJER_per_enery.

	    if( variable.Contains("hJER_per_eDet") ){
	      fitlow = E_gen_cut/hE_gen->GetBinLowEdge( bin_E ) - 1.; // Because of the lower cut on E_gen (300), (E_gen - E_det)/E_det experiences a lower cut off as well, dependent on E_det.
	      //  cout << "Lower bound is\t" << fitlow << "\tfor\t" << hE_gen->GetBinLowEdge( bin_E + 1 ) << endl;
	      fithigh = 5.;
	    } // hJER_per_eDet in variable.
	
	    TH1D* hJER_Slice_2 = (TH1D*)hSlice_to_plot->Clone("hJER_Slice_2");

	    // First Gaussian.

            TF1* f1 = new TF1("f1", distr,  fitlow, fithigh);	    
	    f1->SetLineStyle(2);
	    hSlice_to_plot->Fit("f1", "S", "", fitlow, fithigh);

	    double mean0 = f1->GetParameter(1), spread0 = f1->GetParameter(2);
	    mu0.push_back(mean0), sigma0.push_back(spread0), chi2_0.push_back( f1->GetChisquare() );   

	    // Second Gaussian
	    double fitlow_1 = max( mean0-spread0, fitlow );

            TF1* f2 = new TF1("f2", distr, fitlow_1, mean0+spread0);	    
	    hJER_Slice_2->Fit("f2", "S", "", fitlow_1, mean0+spread0);
	    f2->SetLineColor(kGreen);
	    f2->SetLineStyle(1);
	    f2->Draw("same");

            double mean1 = f2->GetParameter(1), spread1 = f2->GetParameter(2);
            mu1.push_back(mean1), sigma1.push_back(spread1), chi2_1.push_back( f2->GetChisquare() );

            double fitlow_2 = max( mean1-spread1, fitlow );

            TF1* f3 = new TF1("f3", distr, fitlow_2, mean0+spread0);
            hSlice_to_plot->Fit("f3", "S", "", fitlow_2, mean1+spread1);
            f3->SetLineColor(kBlue);
            f3->SetLineStyle(1);

            double mean2 = f3->GetParameter(1), spread2 = f3->GetParameter(2);
            mu2.push_back(mean2), sigma2.push_back(spread2), chi2_2.push_back( f3->GetChisquare() );

            f3->Draw("same"); f1->Draw("same"); f2->Draw("same");

	    TString plot_name = date + TString::Format("_can_Energy_fraction_Slice_%i_", bin_E) +  variable + "_" + distr;   
   		
	    ofstream fitPlot_tex;
   	    fitPlot_tex.open("Plots/Fits/" + date + "_" + variable + "_FitSlices_" + distr + ".tex", std::ofstream::app);
	    
	   cout << "\t\t\t+++ tex file +++" << endl;
	    fitPlot_tex	<< "\\begin{frame}\n";
            		if( variable == "hJER_per_energy" ){ fitPlot_tex << " \\frametitle{$" + variable + "$ as function of $E_{\\text{gen}}$}\n"; }
			else{				     fitPlot_tex << " \\frametitle{$" + variable + "$ as function of $E_{\\text{det}}$}\n"; }	
            fitPlot_tex	<< "  \\begin{center}\n"
			<< "   \\begin{figure}\n"
			<< "    \\includegraphics[width=.7\\linewidth]{Figures/Fits/" << plot_name << "}\n"
                        << "   \\end{figure}\n"
                        << "  \\end{center}\n"
			<< "  \\tiny\n"
			<< "  $\\mu_0 = " << mean0 << "$, $\\sigma_0 = " << spread0 << "$, $\\chi^2 = " << f1->GetChisquare() << "$, depletion below $E_{\\text{det}} = " << fitlow << "$\n\n"
                        << "  $\\mu_1 = " << mean1 << "$, $\\sigma_1 = " << spread1 << "$, $\\chi^2 = " << f2->GetChisquare() << "$, no fit below $E_{\\text{det}} = " 	<< fitlow_1 << "$\n\n"
                        << "  $\\mu_2 = " << mean2 << "$, $\\sigma_2 = " << spread2 << "$, $\\chi^2 = " << f3->GetChisquare() << "$, no fit below $E_{\\text{det}} = " 	<< fitlow_2 << "$\n\n"
			<< "\\end{frame}\n" << endl;

	    fitPlot_tex.close();
	
	    can_Slice->SaveAs( "Plots/Fits/" + plot_name + ".pdf" );
            can_Slice->SaveAs( "Plots/Fits/" + plot_name + ".C" );

	    cout << endl << endl << endl;
	  
	  } // Integral > 0.
	} // Bins.


	//
	// Create graphs from arrays.
	//
	double * E_axis_arr = &E_axis[0];
	double * sigma0_arr = &sigma0[0],
	       * sigma1_arr = &sigma1[0],
	       * sigma2_arr = &sigma2[0],
	       * mu0_arr    = &mu0[0],
	       * mu1_arr    = &mu1[0],
	       * mu2_arr    = &mu2[0],
	       * chi2_0_arr = &chi2_0[0],
	       * chi2_1_arr = &chi2_1[0],
	       * chi2_2_arr = &chi2_2[0],
	       // -- needed for Edet to Egen calibration
	       * Edet_0_arr = &Edet[0],
	       * Egen_0_arr = &Egen[0],
	       * Nevents_0_arr = &Nevents[0],
	       * Edet_1_arr = &Edet[0],
	       * Egen_1_arr = &Egen[0],
	       * Nevents_1_arr = &Nevents[0],
	       * Edet_2_arr = &Edet[0],
	       * Egen_2_arr = &Egen[0],
	       * Nevents_2_arr = &Nevents[0];

	// -- Transform Edet to Egen.
	double max_Egen_0 = 0., max_Egen_1 = 0., max_Egen_2 = 0.;
	for(int mubin = 0; mubin < E_axis.size(); mubin++){
	  Egen_0_arr[mubin] = Edet_0_arr[mubin] * (1. + mu0_arr[mubin]);
          Egen_1_arr[mubin] = Edet_1_arr[mubin] * (1. + mu1_arr[mubin]);
          Egen_2_arr[mubin] = Edet_2_arr[mubin] * (1. + mu2_arr[mubin]);

	  // -- Compact if statement.

	  if( Egen_0_arr[mubin] > max_Egen_0){  max_Egen_0 = Egen_0_arr[mubin]; }
          if( Egen_1_arr[mubin] > max_Egen_1){  max_Egen_1 = Egen_1_arr[mubin]; }
          if( Egen_2_arr[mubin] > max_Egen_2){  max_Egen_2 = Egen_2_arr[mubin]; }
	  
	  cout << "lowedge.push_back(\t" << hE_gen->GetBinLowEdge( mubin ) << "\t);\t muval.push_back(\t" << mu0_arr[mubin] << "\t);" << endl;

	}	

	TGraph *calibrated_plot_0 = new TGraph( E_axis.size(), Egen_0_arr, Nevents_0_arr);
        TGraph *calibrated_plot_1 = new TGraph( E_axis.size(), Egen_0_arr, Nevents_0_arr);
        TGraph *calibrated_plot_2 = new TGraph( E_axis.size(), Egen_2_arr, Nevents_2_arr);

	TCanvas * can_calib = new TCanvas( "can_calib", "can_calib", 1.);
	TH1D* hGenJet_energy = (TH1D*)_file0->Get("hGenJet_energy");
	TH1D* hCastorJet_energy = (TH1D*)_file0->Get("hCastorJet_energy");

	// hGenJet_energy->Scale( 1. / hGenJet_energy->Integral() );
	hGenJet_energy->Scale( 0.12 / hGenJet_energy->GetMaximum() );
	hCastorJet_energy->Scale( 0.12 / hGenJet_energy->GetMaximum() );
	  hCastorJet_energy->SetLineColor( kBlue );

	hGenJet_energy->Draw();
	hCastorJet_energy->Draw("lsame");
	calibrated_plot_0->SetMarkerColor( getColor(1) );
	calibrated_plot_0->SetMarkerStyle( 22 );
	calibrated_plot_0->Draw("*same");
	

	
	
/*
        calibrated_plot_1->SetMarkerColor( getColor(2) );
        calibrated_plot_1->SetMarkerStyle( 23 );

        calibrated_plot_1->Draw("psame");
        calibrated_plot_2->SetMarkerColor( getColor(3) );
        calibrated_plot_2->SetMarkerStyle( 24 );

	calibrated_plot_2->Draw("psame");
*/

	can_calib->SaveAs( date + "_calibrationTest.pdf");
        can_calib->SaveAs( date + "_calibrationTest.C");




	TGraph *sigma0_plot = new TGraph( E_axis.size(), E_axis_arr, sigma0_arr );
        TGraph *sigma1_plot = new TGraph( E_axis.size(), E_axis_arr, sigma1_arr );
        TGraph *sigma2_plot = new TGraph( E_axis.size(), E_axis_arr, sigma2_arr );

	TCanvas *can_sigmas = new TCanvas( date + "_Compare_JER_" +  variable, date + "_Compare_JER_" +  variable + "_" + distr);
	sigma0_plot->Draw("a*e");
	  sigma0_plot->GetXaxis()->SetTitle("E_{det}");
	  sigma0_plot->GetYaxis()->SetTitle("#sigma(#DeltaE/E)");
	  sigma0_plot->GetYaxis()->SetRangeUser(0.,10.);
	  sigma0_plot->SetTitle("JER (E_{det}): " + legendEntries[file]);

	  sigma1_plot->SetMarkerColor(kGreen+3);
	sigma1_plot->Draw("same*");
	  sigma2_plot->SetMarkerColor(kBlue);
	 sigma2_plot->Draw("same*");

	can_sigmas->SaveAs( "Plots/Fits/" + date + "_Compare_JER_" +  variable + "_" + distr + TString::Format("_%i", file) + ".pdf" );
	can_sigmas->SaveAs( "Plots/Fits/" + date + "_Compare_JER_" +  variable + "_" + distr + TString::Format("_%i", file) + ".C" );

	// Do the same for mu.

        TGraph *mu0_plot = new TGraph( E_axis.size(), E_axis_arr, mu0_arr );
        TGraph *mu1_plot = new TGraph( E_axis.size(), E_axis_arr, mu1_arr );
        TGraph *mu2_plot = new TGraph( E_axis.size(), E_axis_arr, mu2_arr );

        TCanvas *can_mus = new TCanvas( date + "_Compare_JER_" +  variable, date + "_Compare_JER_" +  variable + "_" + distr);
        mu0_plot->Draw("a*e");
          mu0_plot->GetXaxis()->SetTitle("E_{det}");
          mu0_plot->GetYaxis()->SetTitle("#mu(#DeltaE/E)");
	  mu0_plot->GetYaxis()->SetRangeUser(-1., 10.);
          mu0_plot->SetTitle("JES (E_{det}): " + legendEntries[file]);	 

          mu1_plot->SetMarkerColor(kGreen+3);
        mu1_plot->Draw("same*");
          mu2_plot->SetMarkerColor(kBlue);
         mu2_plot->Draw("same*");

        can_mus->SaveAs( "Plots/Fits/" + date + "_Compare_JES_" +  variable + "_" + distr + TString::Format("_%i", file) + ".pdf" );
        can_mus->SaveAs( "Plots/Fits/" + date + "_Compare_JES_" +  variable + "_" + distr + TString::Format("_%i", file) + ".C" );

	// Do the same for chi2.

        TGraph *chi2_0_plot = new TGraph( E_axis.size(), E_axis_arr, chi2_0_arr );
        TGraph *chi2_1_plot = new TGraph( E_axis.size(), E_axis_arr, chi2_1_arr );
        TGraph *chi2_2_plot = new TGraph( E_axis.size(), E_axis_arr, chi2_2_arr );

        TCanvas *can_chi2_s = new TCanvas( date + "_Compare_JER_" +  variable, date + "_Compare_JER_" +  variable + "_" + distr);
	can_chi2_s->SetLogy();	

        chi2_0_plot->Draw("a*e");
          chi2_0_plot->GetXaxis()->SetTitle("E_{det}");
          chi2_0_plot->GetYaxis()->SetTitle("#chi2_(#DeltaE/E)");
          chi2_0_plot->SetTitle("JES (E_{det}): " + legendEntries[file]);

          chi2_1_plot->SetMarkerColor(kGreen+3);
        chi2_1_plot->Draw("same*");
          chi2_2_plot->SetMarkerColor(kBlue);
         chi2_2_plot->Draw("same*");

        can_chi2_s->SaveAs( "Plots/Fits/" + date + "_Compare_CHI2_" +  variable + "_" + distr + TString::Format("_%i", file) + ".pdf" );
        can_chi2_s->SaveAs( "Plots/Fits/" + date + "_Compare_CHI2_" +  variable + "_" + distr + TString::Format("_%i", file) + ".C" );

      }// if plot_type, file.
     
      
	/*
	 * We now will plot the distributions of \Delta E/E versus E_gen and E_det.
	 * FitSlicesY creates a series of histograms containing the mean and spread values of each distribution, indicated with _1 and _2.
	 * These histograms will be plotted.
	 */
 
  
       
       hJER->FitSlicesY(0,0,-1,20);   
   
       can_JER->cd();
       hJER->Draw("colz");      
       can_JER->SaveAs("TEST_" + variable + ".C");  
 
       // Mu/JES of our fit.
       can_mu->cd();
       TH1D *hMean = (TH1D*)gDirectory->Get(variable + "_1");
	cout << "=====\t" << variable << "\thMean integral\t" << hMean->Integral() << "\t" << hMean->GetName() << endl;
       hMean->SetName("hMean_" + TString::Format("%i", file) );
       hMean->SetLineWidth(2);
       hMean->SetLineColor(colours[file]);
       hMean->SetLineStyle(file +1);

       hMean->GetXaxis()->SetLabelOffset(0.007);
       hMean->GetXaxis()->SetLabelSize(0.04);
       hMean->GetXaxis()->SetTitleSize(0.07);
       hMean->GetXaxis()->SetTitleOffset(0.8);
       
       hMean->GetYaxis()->SetLabelOffset(0.007);
       hMean->GetYaxis()->SetLabelSize(0.05);
       hMean->GetYaxis()->SetTitleSize(0.07);
       hMean->GetYaxis()->SetTitleOffset(0.8);
       
       hMean->GetYaxis()->SetTitle("#mu");
       hMean->GetXaxis()->SetTitle("E_{det} (GeV)");
       // hMean->GetXaxis()->SetTitle( plotX[variable] );

       hMean->GetXaxis()->SetRangeUser(0., 1000.);
       
	// -- Draw fit.
	// Suppose a function of the form 1 - log(E) 
	TF1* calib_fit = new TF1("calib_fit", "[0] + [1] * exp(x * [2])", 1., 1000.);
	
	hMean->Fit("calib_fit", "S0", "", 0., 1000.);
	
	cout 	<< "\t\t\t\t[0] = " << calib_fit->GetParameter(0) << "\n"
	  	<< "\t\t\t\t[1] = " << calib_fit->GetParameter(1) << endl;
	
	
	calib_fit->SetLineColor( kBlue );
	calib_fit->Draw("same");      
	
	// -- Fit drawn.       

	if( (file == 2 || file == 3) && (variable.Contains("eDet") )){
//	  TF1* fLine = new TF1("fLine", "[0] + [1 ]* x + [2] * x * x + [3] * x * x * x", 100., 2500.);
        TF1* fLine = new TF1("fLine", "[0] + [1] * exp(x)", 100., 1500.);	  
	  hMean->Fit("fLine", "S", "", 100.,1500.);
	  cout << endl << endl << endl << endl << "+++\t" << file << "\t" << variable << "\t" 
		<< fLine->GetParameter(0) << " + \t" 
		<< fLine->GetParameter(1) << " x + \t" 
		<< fLine->GetParameter(2) << " x^2 + \t" 
		<<fLine->GetChisquare() << endl << endl << endl;
	}
	

       if( max_val_sigma < hMean->GetMaximum() ){ max_val_sigma = hMean->GetMaximum();}
	hMean->Draw( drawoptions );

       if(plot == 0) {legend->AddEntry( hMean, legendEntries[file], "l" );}

       // Sigma/JER of our fit.
       can_sigma->cd();
       TH1D *hSpread = (TH1D*)gDirectory->Get(variable + "_2");
       if( max_val_mu < hSpread->GetMaximum() ){ max_val_mu = hSpread->GetMaximum();}
       hSpread->SetName("hSpread_" + TString::Format("%i", file) );
       hSpread->SetLineColor(colours[file]);
       hSpread->SetLineStyle(file+1);

       hSpread->SetLineWidth(2);
       hSpread->GetXaxis()->SetLabelOffset(0.007);
       hSpread->GetXaxis()->SetLabelSize(0.06);
       hSpread->GetXaxis()->SetTitleSize(0.08);
       hSpread->GetXaxis()->SetTitleOffset(0.8);
       // hSpread->GetXaxis()->SetTitle( plotX[variable] );
       hSpread->GetYaxis()->SetLabelOffset(0.007);
       hSpread->GetYaxis()->SetLabelSize(0.06);
       hSpread->GetYaxis()->SetTitleSize(0.05);
       hSpread->GetYaxis()->SetTitleOffset(1.75);
        

//       hSpread->GetXaxis()->SetRangeUser(0., 2300.);


       hSpread->Draw( drawoptions );
	can_sigma->SaveAs("Plots/" + date + "_" + variable + "_" + distr + "_sigma.C" );    
 
       drawoptions = "sameE";
       if( file == 0){
         originalMean = hMean;
         originalSpread = hSpread;
       }         
     } // Loop over files.

     can_sigma->cd();
     if( draw_legend){ legend->Draw(); }
     can_sigma->SaveAs( "Plots/" + date + "_" + variable + "_" + distr + "_sigma.pdf" );

     can_mu->cd();
     if( draw_legend){ legend->Draw(); }
     can_mu->SaveAs( "Plots/" + date + "_" + variable + "_" + distr + "_mu_fit.pdf" );
     can_mu->SaveAs( "Plots/" + date + "_" + variable + "_" + distr + "_mu.C" );
     
 } // Loop over jet types.


 return (0);
}
