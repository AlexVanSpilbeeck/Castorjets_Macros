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

#include "../tdrstyle.C"
#include "../color.h"
using namespace std; 

#define pi 3.14159
#define set_title 0
#define Ethresh_ 150.


void PlotStartingDistributions_comparingEmin(TString distribution);
void PlotStartingDistributions_comparingSelection(TString distribution);
void PlotStartingDistributions_comparingRatiosSelection(TString distribution);
void PlotStartingDistributions_stacked(TString distribution);
void Test_misses();



int main(){
  gStyle->SetOptTitle( 0 );
  setTDRStyle();
  gStyle->SetPalette(1);

//  Test_misses();
//  PlotStartingDistributions_stacked("detector");
//  PlotStartingDistributions_stacked("generator");  
/*
  PlotStartingDistributions_comparingEmin("fake");
  PlotStartingDistributions_comparingEmin("miss");
  PlotStartingDistributions_comparingEmin("detector");
  PlotStartingDistributions_comparingEmin("generator");
  PlotStartingDistributions_comparingEmin("match_meas");
  PlotStartingDistributions_comparingEmin("match_true");
*/  
  PlotStartingDistributions_comparingSelection("purity");
  PlotStartingDistributions_comparingSelection("stability");  

  return 0;
}

//----------------------------------------------------------------------------------------------------//
//-- The following two functions plot the distributions extracted from the MC sample and plot them. --//
//-- They are compared with distributions from MC sample(s) with other Emin or delta phi maxes.	    --//
//----------------------------------------------------------------------------------------------------//

void PlotStartingDistributions_comparingSelection(TString distribution){
  cout << "Unfolder::PlotStartingDistributions_comparingEmin\t" << distribution << endl;

  double ticksize = gStyle->GetTickLength();

  
  ofstream ratios;
  ratios.open("Ratios.txt");
  
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
    TString drawoptions = "hist";
      
    TCanvas *can_startingDistributions;
    PrepareCanvas(can_startingDistributions, "can_startingDistributions");
    TLegend *legend = new TLegend( 0.55, 0.65, 1.-can_startingDistributions->GetTopMargin() - ticksize, 1.-can_startingDistributions->GetRightMargin() - ticksize);
    legend->SetFillStyle( 0 );
    legend->SetBorderSize( 0 );      
      
    for(int _eta = 0; _eta < etaband.size(); _eta++){
      for(int _phi = 0; _phi < deltaPhiMax.size(); _phi++){
	for(int _Emin = 0; _Emin < Emin.size(); _Emin++){
   	    for(int _match = 0; _match < matching.size(); _match++){

 	      TString filename_ = TString::Format(
 	      	"/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo" 
 	      	+ cuts_ 
 	      	+ "_unfold_Emin_150.000000_deltaPhiMax_%f_etaband_%f_all_matchE.root", deltaPhiMax[_phi], etaband[_eta] );
	      
	      TFile* _file= TFile::Open( filename_, "read" );
	      if( !_file ){ 
	        cout << "No\t" << filename_ << endl;
	        continue; 
	      }

	      TH1D* hDistribution;
              RooUnfoldResponse* response = (RooUnfoldResponse*)_file->Get("response");
	      TH2D* hResponse = (TH2D*)response->Hresponse();

	      TString histname = TString::Format( distribution + "_Emin_%i_deltaPhiMax_0%i_etaband_0%i" + matching[_match] + "_" , 
		static_cast<int>(Emin[_Emin]), 
		static_cast<int>(10. * deltaPhiMax[_phi]), 
		static_cast<int>(10. * etaband[_eta]));
		
	      if( distribution == "purity" ){
	      
	        TH1D* hPurity =(TH1D*)hResponse->ProjectionY();
  		hPurity->Reset();
  		hPurity->GetXaxis()->SetTitle("E_{gen} [GeV]");

  		for(int biny = 1; biny <= hResponse->GetNbinsY(); biny++){
  		  double Nx = 0.;
  		  double Ny = 0.;
 
    		  for(int binx = 1; binx <= hResponse->GetNbinsX(); binx++){
     		   Nx += hResponse->GetBinContent( binx, biny );      
    		  }

    		  Ny = hResponse->GetBinContent( biny, biny );
    		  hPurity->SetBinContent( biny, Ny/Nx );
 		}	      
 		hDistribution = hPurity;
	      }//== Purity
	      
	      else if( distribution == "stability" ){
  		TH1D* hStability =(TH1D*)hResponse->ProjectionX();
		hStability->Reset();
  		hStability->GetXaxis()->SetTitle("E_{det} [GeV]"  );

  		for(int binx = 1; binx <= hResponse->GetNbinsX(); binx++){
    		  double Nx = 0.;
    		  double Ny = 0.;
 
                  for(int biny = 1; biny <= hResponse->GetNbinsY(); biny++){
                    Ny += hResponse->GetBinContent( binx, biny );      
                  }
                  Nx = hResponse->GetBinContent( binx, binx );
                  hStability->SetBinContent( binx, Nx/Ny );
                }//== Loop over bins X	
 		hDistribution = hStability;      	     
	      }//== Stability

	      hDistribution->SetTitle( histname );
	      hDistribution->SetName( histname );
	      hDistribution->SetLineWidth( 3 );

	      hDistribution->SetLineStyle( _eta + _phi * etaband.size() + 1 ); 
	      hDistribution->SetLineColor( getColor( _eta + _phi * etaband.size() + 1 ) ); 
 
              hDistribution->SetMarkerStyle( _eta + _phi * etaband.size() + 1 );

              hDistribution->SetMarkerColor( getColor( _eta + _phi * etaband.size() + 1 ) );

	      hDistribution->SetMarkerSize( 1.8 );
	      hDistribution->GetXaxis()->SetNdivisions( 504 );

	      can_startingDistributions->cd();
	      Prepare_1Dplot( hDistribution );	      
	      hDistribution->GetYaxis()->SetRangeUser(0.,0.55);
	      hDistribution->DrawClone( drawoptions );
	      drawoptions = "phsame";
	      legend->AddEntry( hDistribution, TString::Format( "#Delta#phi_{max}=0.%i, #eta_{acc}=0.%i", static_cast<int>(10.*deltaPhiMax[_phi]), static_cast<int>(10.*etaband[_eta]) ), "lp" );
	      legend->SetFillColor( 0 );

	    } // Match.
	  } // Emin.
        } // Deltaphi.
      } // Eta.
      legend->Draw();

      can_startingDistributions->SaveAs( TString::Format( "can_" + distribution + "_" + save_cuts[ cuts_ ] + ".C" ) );
      can_startingDistributions->SaveAs( TString::Format( "can_" + distribution + "_" + save_cuts[ cuts_ ] + ".pdf" ) );     

  } // Event selection.       
  ratios.close();
}














