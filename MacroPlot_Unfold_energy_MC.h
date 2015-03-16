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
#include "Function_make_Tex.h"
#include "Function_FirstPlot.h"

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

void PlotCorrectionFactors_Bayes(TH1D* &histo, TString filename , TString legend, TString iterations_){
   gROOT->SetStyle ("Plain");
   gStyle->SetPalette(1);
   gSystem->Load("~/RooUnfold-1.1.1/libRooUnfold.so");

   /* Open DATA and MC */

   TFile *_file0 = TFile::Open( filename, "Read");
 
   if( _file0->GetListOfKeys()->Contains( "response" ) ){
 
   /* Extract GEN and DET (MC) and DATA distribution */

   RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response");

   TH1D *hDet = (TH1D*)_file0->Get("hCastorJet_energy");	hDet->Sumw2();
  
   /****************/
   /** Unfold MC. **/
   /****************/

   int iterations;
   (iterations_ == "") ? iterations = 4 :  iterations = atoi(iterations_);

   RooUnfoldBayes unfold_mc_bayes (response, hDet, 4);
     TH1D* hReco_mc_bayes= (TH1D*) unfold_mc_bayes.Hreco();

    hReco_mc_bayes->SetLineStyle(3);
     hReco_mc_bayes->SetLineColor(kGreen);
     hReco_mc_bayes->SetLineWidth(2);
     hReco_mc_bayes->SetName("hReco_bayes_1");
     
   histo = (TH1D*)hReco_mc_bayes->Clone(legend);
   
   
   }
}

//-- Overloaded functon

void PlotCorrectionFactors_Bayes_chi2(TH1D* &histo, TString filename , double &chi2, int iterations){
   gROOT->SetStyle ("Plain");
   gStyle->SetPalette(1);
   gSystem->Load("~/RooUnfold-1.1.1/libRooUnfold.so");

   /* Open DATA and MC */

   TFile *_file0 = TFile::Open( filename, "Read");

   if( _file0->GetListOfKeys()->Contains( "response" ) ){

   /* Extract GEN and DET (MC) and DATA distribution */

   RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response");

   TH1D *hDet = (TH1D*)_file0->Get("hCastorJet_energy");        hDet->Sumw2();	//hDet->Scale( 1./hDet->Integral() );
   TH1D *hGen = (TH1D*)_file0->Get("hGenJet_energy");		hGen->Sumw2();	hGen->Scale( 1./hGen->Integral() );
   


   /****************/
   /** Unfold MC. **/
   /****************/

   RooUnfoldBayes unfold_mc_bayes (response, hDet, iterations);
     TH1D* hReco_mc_bayes= (TH1D*) unfold_mc_bayes.Hreco();

    hReco_mc_bayes->SetLineStyle(3);
     hReco_mc_bayes->SetLineColor(kGreen);
     hReco_mc_bayes->SetLineWidth(2);
     hReco_mc_bayes->SetName("hReco_bayes_1");

   histo = (TH1D*)hReco_mc_bayes->Clone(TString::Format("Bayes_%i", iterations));

   chi2 = unfold_mc_bayes.Chi2(hGen);

   }
}



// --------------------------------------------------------------------------------------

/*************************************
* Plot Bayes with several iterations *
*************************************/

void Bayes_iterations( TH1D* &histo, TString filename , TString legend, TString iterations_){

  cout << "Bayes iterations " << iterations_ << endl;

  TFile *_file0 = TFile::Open( filename, "Read");
  if( _file0->GetListOfKeys()->Contains( "response" ) ){
  
    RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response");
    TH1D *hGen = (TH1D*)_file0->Get("hGenJet_energy"); 		// hGen->Rebin(4);   
    TH1D *hDet = (TH1D*)_file0->Get("hCastorJet_energy");	// hDet->Rebin(4);    

  cout << "Integral gen " << hGen->Integral() << endl;

   int iterations;
   (iterations_ == "") ? iterations = 4 :  iterations = atoi(iterations_);
   cout << "Iterations " << iterations << endl;

   TH1D* current_bayes;
   TH1D* current_ratio;
   TH1D* original;    
   TH1D* original_com;
   vector<double> iteration_vector;
   vector<double> chi2_vector;
   double max_val=0., min_val=0.;
   double max_val_com = 0., min_val_com = 0.;

   TCanvas *can = new TCanvas("can", "can", 1.);
   TCanvas *com = new TCanvas("comparison", "comparison", 1.);
   
   TLegend *legend = new TLegend(0.65, 0.65, 0.99, 0.99);
     legend->SetFillColor( kWhite );   
   
   TString drawoptions = "";
   double chi2 = -1.;
   int colours = 1; 
   
   double gen_scale = 0.;

   for(int it = 1; it <= iterations; it++){
     PlotCorrectionFactors_Bayes_chi2(current_bayes, filename, chi2, it);
     iteration_vector.push_back(static_cast<double>(it));
     
     chi2_vector.push_back(chi2/current_bayes->GetNbinsX());
     
     current_bayes->SetLineColor( getColor(colours++) );
     

     com->cd();  
     
     current_ratio = (TH1D*)current_bayes->Clone(TString::Format("ratio_%i", it));   
     
     cout << "Integral Gen " << hGen->Integral() << "\tIntegral unfold " << current_bayes->Integral() << endl;
     
     current_ratio->Divide(hGen);
     current_ratio->SetLineStyle(colours);
     First_Plot( original_com, current_ratio, it-1, max_val_com, min_val_com); 
     current_ratio->Draw("hist" + drawoptions);     
     
     com->Update();

     current_bayes->Scale(1./current_bayes->Integral() );
     First_Plot( original, current_bayes, it-1, max_val, min_val); 

     can->cd();
     current_bayes->SetName(TString::Format("iteration%i", it));
     current_bayes->Draw("hist" + drawoptions);
     can->Update();   

     
     
     drawoptions = "same";
     legend->AddEntry( current_bayes, TString::Format("%i iterations", it), "l");
     
   }
    hGen->Scale(1./ hGen->Integral() );   
   
   can->cd();
   hGen->SetMarkerStyle(21);
   hGen->SetMarkerSize(0.6);
   hGen->Draw("samep");
   legend->AddEntry( hGen, "Generator level", "p");
   legend->Draw();
      
   filename.ReplaceAll(".","");
   filename.ReplaceAll("/","_");

   can->SaveAs("LetsTestIterations_" + filename + ".pdf"); 
   can->SaveAs("LetsTestIterations_" + filename + ".C"); 

   com->cd();
   hGen->Divide(hGen);
   hGen->Draw("samep");
   legend->Draw();

   com->SaveAs("LetsTest_Ratios_" + filename + ".pdf"); 
   com->SaveAs("LetsTest_Ratios_" + filename + ".C"); 


   //-- Plot evolution of chi2.

   double *X_axis = &iteration_vector[0];
   double *Y_axis = &chi2_vector[0];

   can->cd();

   TGraph *iteration_plot = new TGraph( iterations, X_axis, Y_axis );
   iteration_plot->SetMarkerStyle(22);
   iteration_plot->SetMarkerSize(2);
   iteration_plot->Draw("ap");
   can->SaveAs("LetsTestIterations_chi2_" + filename + ".pdf"); 
   can->SaveAs("LetsTestIterations_chi2_" + filename + ".C"); 
   
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
   if( _file0->GetListOfKeys()->Contains( "hJER_per_eDet" ) ){ 
   
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
      
     for(int bin_E = 0; bin_E < E_axis.size(); bin_E++){
       hE->SetBinContent( bin_E, mu_arr[bin_E]);
     }
     histo = (TH1D*)hE->Clone(legend);
   }
}


// --------------------------------------------------------------------------------------

/**********************************
* Plot Bayes unfolding vs. Gen *
**********************************/  

void Plot_Bayes_vs_Gen(TH1D* &histo, TString filename , TString legend, TString not_needed){

  TH1D* hGen;
  	Plot_GenLevel(hGen, filename, legend, "");
  TH1D* hBayes;
  	PlotCorrectionFactors_Bayes( hBayes, filename, legend, "");
	
  histo = (TH1D*)hBayes->Clone(legend);
  
//  histo->Scale(1./histo->Integral() );
//  hGen->Scale(1./hGen->Integral() );
  
  histo->Divide(hGen);
} 


// --------------------------------------------------------------------------------------

/************************************
* Plot Bin-by-bin unfolding vs. Gen *
************************************/  

void Plot_BinByBin_vs_Gen(TH1D* &histo, TString filename , TString legend, TString not_needed){

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




// --------------------------------------------------------------------------------------

 /**********************
 * Unfold data - Bayes *
 **********************/ 

void Unfold_data_Bayes( TH1D* &histo, TString MC_file, TString legend, TString data_file){

   gROOT->SetStyle ("Plain");
   gStyle->SetPalette(1);
   gSystem->Load("~/RooUnfold-1.1.1/libRooUnfold.so");

   /* Open DATA and MC */

   TFile *_file_MC = TFile::Open( MC_file, "Read");
   if( _file_MC->GetListOfKeys()->Contains( "response" ) ){

     /* Extract GEN and DET (MC) and DATA distribution */

     RooUnfoldResponse* response = (RooUnfoldResponse*)_file_MC->Get("response");
     TH1D *hGen = (TH1D*)_file_MC->Get("hGenJet_energy");   hGen->Sumw2();

     TFile *_file_data = TFile::Open(data_file, "Read");
     TH1D *hDet = (TH1D*)_file_data->Get("hCastorJet_energy");   hDet->Sumw2();

     int iterations = 4;

     RooUnfoldBayes unfold_mc_bayes (response, hDet, iterations);
     TH1D* hReco_mc_bayes= (TH1D*) unfold_mc_bayes.Hreco();

     hReco_mc_bayes->SetName(TString::Format("hReco_bayes_%i_iterations", iterations));

     histo = (TH1D*)hReco_mc_bayes->Clone(legend);
   }
}


// --------------------------------------------------------------------------------------

  /***************************
  * Unfold data - bin-by-bin *
  ***************************/
    


void Unfold_data_BinByBin( TH1D* &histo, TString MC_file, TString legend, TString data_file){

   gROOT->SetStyle ("Plain");
   gStyle->SetPalette(1);
   gSystem->Load("~/RooUnfold-1.1.1/libRooUnfold.so");

   /* Open DATA and MC */

   TFile *_file_MC = TFile::Open( MC_file, "Read");
   if( _file_MC->GetListOfKeys()->Contains( "response" ) ){

     /* Extract GEN and DET (MC) and DATA distribution */

     RooUnfoldResponse* response = (RooUnfoldResponse*)_file_MC->Get("response");
     TH1D *hGen = (TH1D*)_file_MC->Get("hGenJet_energy");   hGen->Sumw2();

     TFile *_file_data = TFile::Open(data_file, "Read");
     TH1D *hDet = (TH1D*)_file_data->Get("hCastorJet_energy");   hDet->Sumw2();

     int iterations = 4;

     RooUnfoldBinByBin unfold_mc_binbybin (response, hDet);
     TH1D* hReco_mc_binbybin= (TH1D*) unfold_mc_binbybin.Hreco();

     hReco_mc_binbybin->SetName(TString::Format("hReco_binbybin"));

     histo = (TH1D*)hReco_mc_binbybin->Clone(legend);
   }
}



// --------------------------------------------------------------------------------------


  /*****************************
  * Create 2D plots from file. *
  *****************************/

void Plot_2D_Energy_Response(TH1D* &histo, TString MC_file, TString legend, TString label){

  cout << "\t\t\tLabel is\t" << label << endl;

  TFile *_file0 = TFile::Open( MC_file ); 

  TString save_name;
  save_name = MC_file;
  save_name.ReplaceAll("LoopRootFiles/","");
  save_name.ReplaceAll(".root","");


  vector<TString> distributions;
   map<TString, TString> xTitle;
   map<TString, TString> yTitle;

  TString plot = "";
  
  /****************************
  * Prepare plots and titles. * 
  ****************************/	

  plot = "hJER_per_energy";
  distributions.push_back( plot );
    xTitle[plot] = "E_{gen} (GeV)";
    yTitle[plot] = "#frac{E_{det}-E_{gen}}{E_{gen}}";

  plot = "hJER_per_eDet";
  distributions.push_back( plot );
    xTitle[plot] = "E_{det} (GeV)";
    yTitle[plot] = "#frac{E_{gen}-E_{det}}{E_{det}}";

  plot = "hJER_per_eGen";
  distributions.push_back( plot );
    xTitle[plot] = "E_{gen} (GeV)";
    yTitle[plot] = "#frac{E_{gen}-E_{det}}{E_{det}}";

  plot = "hCastorJet_energy_response";
  distributions.push_back( plot );
    xTitle[plot] = "E_{det} (GeV)";
    yTitle[plot] = "E_{gen} (GeV)";

  /******************
  * Prepare Canvas. *
  ******************/  
  
  TCanvas *can = new TCanvas("can", "can", 1.);
    can->SetLeftMargin(0.20);
    can->SetRightMargin(0.18);
    can->SetBottomMargin(0.20);  
    can->SetTopMargin(0.02);

  vector<TString> plot_names;

  for(int plot_ = 0; plot_ < distributions.size(); plot_++){
 
    if( !_file0->GetListOfKeys()->Contains( distributions[plot_] ) ){ continue; }

    TH2D *hDistr = (TH2D*)_file0->Get( distributions[plot_] );
    can->cd();
    hDistr->Scale( 1./hDistr->Integral() );
    hDistr->Draw("colz");
    hDistr->GetXaxis()->SetTitle( xTitle[plot_] ); 
    hDistr->GetYaxis()->SetTitle( yTitle[plot_] );
    can->SetLogz();

    TString plot_name_ = label + "/" + distributions[plot_] + "_" + save_name;
    plot_name_.ReplaceAll("1.", "1");
    plot_name_.ReplaceAll("0.", "0");

    plot_names.push_back( plot_name_ );

    can->SaveAs( "Plots/" + plot_names[plot_] + ".pdf" );
    can->SaveAs( "Plots/" + plot_names[plot_] + ".C" );
  }

  if( plot_names.size() == 4 ){
    Four_plots( TString::Format(plot_names[0] + ".pdf"), TString::Format(plot_names[1] + ".pdf"), TString::Format(plot_names[2] + ".pdf"), TString::Format(plot_names[3] + ".pdf"), "Energy response distributions: " + legend,"2D_plots.tex");
  }


}
