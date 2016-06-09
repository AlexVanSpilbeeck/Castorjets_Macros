#include "iostream"
#include "fstream"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TColor.h"
#include "TApplication.h"
#include "stdio.h"
#include "cstring"
#include "TMath.h"
#include <sys/stat.h>

#include "../Functions/Function_average_histogram.h"
#include "../Functions/Function_CorrectSectorNumber.h"
#include "../Functions/Function_FinishCanvas.C"
#include "../Functions/Function_Prepare1Dplot.h"
#include "../Functions/Function_Prepare2Dplot.h"
#include "../Functions/Function_PrepareCanvas.h"
#include "../Functions/Function_Rebin.h"
#include "../Functions/Function_SetRangeToIncludeAll.h"
#include "../GetSubHistogram.h"
#include "../Correct_Sector.h"
//#include "Function_make_Tex.h"
#include "../Functions/Function_NonZeroMinimum.h"

#include "../../../../RooUnfold-1.1.1/src/RooUnfold.h"
#include "../../../../RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "../../../../RooUnfold-1.1.1/src/RooUnfoldBayes.h"

#include "../tdrstyle.C"
#include "../color.h"
using namespace std; 

#define pi 3.14159
#define set_title 0
#define Ethresh_ 150.


void Plot_Unfolded();



int main(){
  gStyle->SetOptTitle( 0 );
  setTDRStyle();
  gStyle->SetPalette(1);

  Plot_Unfolded();


  return 0;
}

//----------------------------------------------------------------------------------------------------//
//-- The following two functions plot the distributions extracted from the MC sample and plot them. --//
//-- They are compared with distributions from MC sample(s) with other Emin or delta phi maxes.	    --//
//----------------------------------------------------------------------------------------------------//





void Plot_Unfolded(){

  double ticksize = gStyle->GetTickLength();

  
  gStyle->SetOptTitle( 0 );
  gStyle->SetOptStat( 0 );

  int file = 0;

  vector<double> etaband;
  map<double, int> eta_markers;  
      etaband.push_back(0.0);
      etaband.push_back(0.2);
      etaband.push_back(0.5);

  vector<double> deltaPhiMax;
      deltaPhiMax.push_back(0.2);
      deltaPhiMax.push_back(0.5);

  vector<TString> matching;
      matching.push_back("_matchE");

  vector<TString> match_symbol;
      match_symbol.push_back("E");
      match_symbol.push_back("#varphi");

  vector<TString> model;
      model.push_back("");

  vector<double> Emin;
      Emin.push_back(150.);


  TString cuts_;
  std::vector<TString> vector_cuts_;
  std::map<TString, TString> legend_cuts;
  std::map<TString, TString> save_cuts;  
        
  cuts_ = "";
  vector_cuts_.push_back( cuts_ );
  legend_cuts[cuts_] = "Default selection";
  save_cuts[ cuts_ ] = "fwd11003";
  
/*
  cuts_ = "_BSC_OR";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "BSC OR";


  cuts_ = "_noVertex_HF";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "HF, No vertex";


  cuts_ = "_noVertex_BSC";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "BSC OR, No vertex";

  cuts_ = "_atLeastOneVertex_HF";

  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "HF, #geq 1 Vtx";

  cuts_ = "_atLeastOneVertex_BSC";

  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "BSC, #geq 1 Vtx";
*/
    /*
  cuts_ = "_one_none_vertex_BSCor_HFor";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "BSC OR, HF OR, #leq 1 gVtx";   
  */
    
  cuts_ = "_MostInclusive";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "BSC OR, HF OR, #leq 1 gVtx";   
  save_cuts[ cuts_ ] = "mostinclusive";         

  for(int _cuts = 0; _cuts < vector_cuts_.size(); _cuts++, file++) {
    cuts_ = vector_cuts_[_cuts];
              
    for(int _eta = 0; _eta < etaband.size(); _eta++){
      for(int _phi = 0; _phi < deltaPhiMax.size(); _phi++){
	for(int _Emin = 0; _Emin < Emin.size(); _Emin++){
   	    for(int _match = 0; _match < matching.size(); _match++){
		      
	      //== Canvas.
    	      TCanvas *can_unfoldedDistributions;
	      PrepareCanvas(can_unfoldedDistributions, "can_unfoldedDistributions");
    	      TString drawoptions = "hist";

	      //== Legend.
	      double ticksize = gStyle->GetTickLength();
	      TLegend *legend = new TLegend( 0.60, 0.4, 1.-can_unfoldedDistributions->GetRightMargin() - ticksize, 1.-can_unfoldedDistributions->GetTopMargin() - ticksize);
	      legend->SetFillStyle( 0 );
	      legend->SetBorderSize( 0 );
	      legend->SetHeader( TString::Format( "#Delta#phi_{max} = 0.%i, #eta_{acc} = 0.%i", 
	        static_cast<int>(10.*deltaPhiMax[_phi]), 
	        static_cast<int>(10.*etaband[_eta]) ) );
	      

	      int color_ = 1;

 	      TString MCname_ = TString::Format(
 	      	"/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo" 
 	      	+ cuts_ 
 	      	+ "_unfold_Emin_150.000000_deltaPhiMax_%f_etaband_%f_all_matchE.root", deltaPhiMax[_phi], etaband[_eta] );
	      
	      TFile* _MCfile= TFile::Open( MCname_, "read" );
	      if( !_MCfile ){ continue; }
	      
 	      TString dataname_ = TString::Format(
 	      	"/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5_data" + cuts_ + "_unfold_Emin_150.000000_all.root" );
	      
	      TFile* _datafile= TFile::Open( dataname_, "read" );
	      if( !_datafile ){ continue; }	      

	      //== Response object from RooUnfold.
	      can_unfoldedDistributions->cd();
	      RooUnfoldResponse* response = (RooUnfoldResponse*)_MCfile->Get("response");
	      
	      //== data distribution.
	      TH1D* hData = (TH1D*)_datafile->Get("hCastorJet_energy");
	      
      	      for(int iterations = 1; iterations <= 30; iterations += 3, color_++){
                // Determine the inrease in Bayesian iterations.

      		RooUnfoldBayes unfold_bayes(response, hData, iterations); 
      		unfold_bayes.SetVerbose(0);
      		TH1D* hUnfold = (TH1D*) unfold_bayes.Hreco( RooUnfold::kCovariance );	            		
      		hUnfold->GetYaxis()->SetTitle("d#sigma/dE [mb/GeV]");
      		hUnfold->GetXaxis()->SetTitle("E [GeV]");
      		hUnfold->SetLineColor( getColor( color_ ) );
      		hUnfold->SetLineStyle( color_ );
      		hUnfold->SetLineWidth( 2 );
      		SetDnDx( hUnfold );      
      		Prepare_1Dplot( hUnfold );
      		hUnfold->GetXaxis()->SetLabelSize( hUnfold->GetXaxis()->GetLabelSize()*0.8 );	
      		
      		legend->AddEntry( hUnfold, TString::Format("%i it.", iterations), "l");      		
      		hUnfold->Draw( drawoptions );      		
      		drawoptions = "histsame";
              }  // Iterations.
              
              TString savename = TString::Format( 
	        "can_unfoldedDistributions_phi_0%i_eta_0%i_" + save_cuts[ cuts_ ], 
	        static_cast<int>(10.*deltaPhiMax[_phi]), 
	        static_cast<int>(10.*etaband[_eta]) );
              can_unfoldedDistributions->cd();
              can_unfoldedDistributions->SetLogy();
              legend->Draw();	        
	      can_unfoldedDistributions->SaveAs( savename + ".C" );
	      can_unfoldedDistributions->SaveAs( savename + ".pdf" ); 
      	    } // Match.
	  } // Emin.
        } // Deltaphi.
      } // Eta.


  } // Event selection. 
}




