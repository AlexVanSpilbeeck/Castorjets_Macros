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


#include "../color.h"
using namespace std; 

#define pi 3.14159
#define set_title 0
#define Ethresh_ 150.


void PlotStartingDistributions_comparingEmin(TString distribution);
void PlotStartingDistributions_comparingSelection(TString distribution);
void PlotStartingDistributions_comparingRatiosSelection(TString distribution);
void PlotStartingDistributions_stacked(TString distribution);



int main(){

  PlotStartingDistributions_stacked("detector");
  PlotStartingDistributions_stacked("generator");  

  PlotStartingDistributions_comparingEmin("fake");
  PlotStartingDistributions_comparingEmin("miss");
  PlotStartingDistributions_comparingEmin("detector");
  PlotStartingDistributions_comparingEmin("generator");
  PlotStartingDistributions_comparingEmin("match_meas");
  PlotStartingDistributions_comparingEmin("match_true");
  
  PlotStartingDistributions_comparingSelection("fake");
  PlotStartingDistributions_comparingSelection("miss");
  PlotStartingDistributions_comparingSelection("detector");
  PlotStartingDistributions_comparingSelection("generator");
  PlotStartingDistributions_comparingSelection("match_meas");
  PlotStartingDistributions_comparingSelection("match_true");  
  
  /*
  PlotStartingDistributions_comparingRatiosSelection("fake");
  PlotStartingDistributions_comparingRatiosSelection("miss");  
*/
  return 0;
}

//----------------------------------------------------------------------------------------------------//
//-- The following two functions plot the distributions extracted from the MC sample and plot them. --//
//-- They are compared with distributions from MC sample(s) with other Emin or delta phi maxes.	    --//
//----------------------------------------------------------------------------------------------------//





void PlotStartingDistributions_comparingEmin(TString distribution){
  cout << "Unfolder::PlotStartingDistributions_comparingEmin\t" << distribution << endl;

  double ticksize = gStyle->GetTickLength();

  TCanvas *can_startingDistributions;
  PrepareCanvas(can_startingDistributions, "can_startingDistributions");
  TLegend *legend = new TLegend( 0.55, 0.65, 1.-can_startingDistributions->GetTopMargin() - ticksize, 1.-can_startingDistributions->GetRightMargin() - ticksize);
  legend->SetFillStyle( 0 );
  legend->SetBorderSize( 0 );
  
  gStyle->SetOptTitle( 0 );
  gStyle->SetOptStat( 0 );
    
  TString distributionname;
  if (distribution == "fake"){ 		distributionname = "hCastorJet_fake_all"; }
  if (distribution == "detector"){ 	distributionname = "hCastorJet_energy"; }
  if (distribution == "generator"){ 	distributionname = "hGenJet_energy"; }
  if (distribution == "miss"){ 		distributionname = "hCastorJet_miss_all"; }
  if (distribution == "match_meas"){	distributionname = "match_meas";}
  if (distribution == "match_true"){	distributionname = "match_true";}

  TString drawoptions = "phist";

  int file = 0;

  vector<double> etaband;
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

  vector<TString> model_legend;
      model_legend.push_back("p6");
      model_legend.push_back("p8");

  vector<double> Emin;
      Emin.push_back(150.);

  vector<int> markers;
      markers.push_back( 20 );
      markers.push_back( 22 );
      markers.push_back( 29 );
      markers.push_back( 21 );

    for(int _eta = 0; _eta < etaband.size(); _eta++){
      for(int _phi = 0; _phi < deltaPhiMax.size(); _phi++){
	for(int _Emin = 0; _Emin < Emin.size(); _Emin++){
   	    for(int _match = 0; _match < matching.size(); _match++, file++){

		TString filename_ = TString::Format("/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo_unfold_Emin_150.000000_deltaPhiMax_%f_etaband_%f_all_matchE.root", deltaPhiMax[_phi], etaband[_eta] );
	      cout << "\n=*=\t" << filename_ << endl;
	      TFile* _file= TFile::Open( filename_, "read" );

	      TString histname = TString::Format( distributionname + "_Emin_%i_deltaPhiMax_0%i_etaband_0%i" + matching[_match] + "_" , 
		static_cast<int>(Emin[_Emin]), 
		static_cast<int>(10. * deltaPhiMax[_phi]), 
		static_cast<int>(10. * etaband[_eta]));

    	      //== Extracting the distribution.
  	      TH1D* hDistribution = (TH1D*)_file->Get( distributionname );

	      if( !distributionname.Contains("match") ){
		hDistribution = (TH1D*)_file->Get( distributionname );
	      }
	      else{
	        RooUnfoldResponse* response = (RooUnfoldResponse*)_file->Get("response");
	        TH2D* hResponse = (TH2D*)response->Hresponse();
	        if( distributionname = "match_meas" ){
	          hDistribution = (TH1D*)hResponse->ProjectionX();
	        }
	        else if( distributionname = "match_true" ){
	          hDistribution = (TH1D*)hResponse->ProjectionY();
	        }
	      }

	      hDistribution->SetTitle( histname );
	      hDistribution->SetName( histname );
	      hDistribution->SetLineWidth( 3 );

	      hDistribution->SetLineStyle( file + 1 ); 
	      hDistribution->SetLineColor( getColor(file + 1) ); 
 
              hDistribution->SetMarkerStyle( file + 20 );

              hDistribution->SetMarkerColor( getColor( file+1 ) );

	      hDistribution->SetMarkerSize( 1.8 );
	      hDistribution->GetXaxis()->SetNdivisions( 504 );
	      hDistribution->GetXaxis()->SetTitle("E_{" + distribution + "}[GeV]");
	      hDistribution->GetYaxis()->SetTitle("d#sigma/dE [GeV]");

	      can_startingDistributions->cd();
	      Prepare_1Dplot( hDistribution );
	      hDistribution->DrawClone( drawoptions );
	      drawoptions = "phsame";
	      legend->AddEntry( hDistribution, TString::Format( " #Delta#varphi_{max}=0.%i, #eta_{acc}=0.%i", 
		static_cast<int>(10. * deltaPhiMax[_phi] ), 
		static_cast<int>(10. * etaband[_eta] ) ), "lp" );
	      legend->SetFillColor( 0 );


	    } // Match.
	  } // Emin.
        } // Deltaphi.
      } // Eta.
   
  legend->Draw();

  can_startingDistributions->SetLogy();
  can_startingDistributions->SaveAs( TString::Format( "can_startingDistribution_" + distribution + "_Emin_%i.C", static_cast<int>( Ethresh_ ) ) );
  can_startingDistributions->SaveAs( TString::Format( "can_startingDistribution_" + distribution + "_Emin_%i.pdf", static_cast<int>( Ethresh_) ) );

    for(int _eta = 0; _eta < etaband.size(); _eta++){
      for(int _phi = 0; _phi < deltaPhiMax.size(); _phi++){
	for(int _Emin = 0; _Emin < Emin.size(); _Emin++, file++){

//          TString filename_ = TString::Format( "/user/avanspil/Castor_Analysis/ak5ak5_displaced_unfold_Emin_%f_deltaPhiMax_%f_etaband_%f.root", Emin[_Emin], deltaPhiMax[_phi], etaband[_eta] ) ;
		TString filename_ = TString::Format("/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo_unfold_Emin_150.000000_deltaPhiMax_%f_etaband_%f_all_matchE.root", deltaPhiMax[_phi], etaband[_eta] );
	  TFile* _file= TFile::Open( filename_, "read" );

	  RooUnfoldResponse *resp = (RooUnfoldResponse*)_file->Get("response");
	  TH2D* hResponse = (TH2D*)resp->Hresponse();
	  TH1D* hTruth = (TH1D*)resp->Htruth();
	  TH1D* hMeasured = (TH1D*)resp->Hmeasured();

	  double rFake = 1. - hResponse->Integral()/hMeasured->Integral() ;
	  double rMiss = 1. - hResponse->Integral()/hTruth->Integral() ;

	  cout << "\t" << deltaPhiMax[ _phi ] << "\t" << etaband[ _eta ] 
		<< "\tFraction fakes\t" << rFake
		<< "\tFraction misses\t" << rMiss << "\t\t" << ( 1. - rFake)*(1. + rMiss) << endl << endl << endl;

	}
      }
    } 
}








void PlotStartingDistributions_comparingSelection(TString distribution){
  cout << "Unfolder::PlotStartingDistributions_comparingEmin\t" << distribution << endl;

  double ticksize = gStyle->GetTickLength();
  TCanvas *can_startingDistributions;
  PrepareCanvas(can_startingDistributions, "can_startingDistributions");
  TLegend *legend = new TLegend( 0.55, 0.65, 1.-can_startingDistributions->GetTopMargin() - ticksize, 1.-can_startingDistributions->GetRightMargin() - ticksize);
  legend->SetFillStyle( 0 );
  legend->SetBorderSize( 0 );
  
  gStyle->SetOptTitle( 0 );
  gStyle->SetOptStat( 0 );
    
  TString distributionname;
  if (distribution == "fake"){ 		distributionname = "hCastorJet_fake_all"; }
  if (distribution == "detector"){ 	distributionname = "hCastorJet_energy"; }
  if (distribution == "generator"){ 	distributionname = "hGenJet_energy"; }
  if (distribution == "miss"){ 		distributionname = "hCastorJet_miss_all"; }
  if (distribution == "match_meas"){	distributionname = "match_meas";}
  if (distribution == "match_true"){	distributionname = "match_true";}

  TString drawoptions = "phist";

  int file = 0;

  vector<double> etaband;
  map<double, int> eta_markers;  
      etaband.push_back(0.0);

  vector<double> deltaPhiMax;
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
        
  cuts_ = "";
  vector_cuts_.push_back( cuts_ );
  legend_cuts[cuts_] = "Default selection";
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
    
  cuts_ = "_one_none_vertex_BSCor_HFor";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "BSC OR, HF OR, #leq 1 gVtx";   
  
    
  cuts_ = "_MostInclusive";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "BSC OR, HF OR, #leq 1 gVtx";        


    for(int _eta = 0; _eta < etaband.size(); _eta++){
      for(int _phi = 0; _phi < deltaPhiMax.size(); _phi++){
	for(int _Emin = 0; _Emin < Emin.size(); _Emin++){
   	    for(int _match = 0; _match < matching.size(); _match++){
   	      for(int _cuts = 0; _cuts < vector_cuts_.size(); _cuts++, file++) {
   	        cuts_ = vector_cuts_[_cuts];

		TString filename_ = TString::Format("/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo" + cuts_ + "_unfold_Emin_150.000000_deltaPhiMax_%f_etaband_%f_all_matchE.root", deltaPhiMax[_phi], etaband[_eta] );
	      cout << "\n=*=\t" << filename_ << endl;	      
	      TFile* _file= TFile::Open( filename_, "read" );
	      if( !_file ){ 
	        cout << "No\t" << filename_ << endl;
	        continue; 
	      }

	      TString histname = TString::Format( distributionname + "_Emin_%i_deltaPhiMax_0%i_etaband_0%i" + matching[_match] + "_" , 
		static_cast<int>(Emin[_Emin]), 
		static_cast<int>(10. * deltaPhiMax[_phi]), 
		static_cast<int>(10. * etaband[_eta]));

    	      //== Extracting the distribution.
  	      TH1D* hDistribution = (TH1D*)_file->Get( distributionname );
  	      if( !distributionname ) continue; 

	      if( !distributionname.Contains("match") ){
		hDistribution = (TH1D*)_file->Get( distributionname );
	      }
	      else{
	        RooUnfoldResponse* response = (RooUnfoldResponse*)_file->Get("response");
	        TH2D* hResponse = (TH2D*)response->Hresponse();
	        if( distributionname = "match_meas" ){
	          hDistribution = (TH1D*)hResponse->ProjectionX();
	        }
	        else if( distributionname = "match_true" ){
	          hDistribution = (TH1D*)hResponse->ProjectionY();
	        }
	      }

	      hDistribution->SetTitle( histname );
	      hDistribution->SetName( histname );
	      hDistribution->SetLineWidth( 3 );

	      hDistribution->SetLineStyle( file + 1 ); 
	      hDistribution->SetLineColor( getColor(file + 1) ); 
 
              hDistribution->SetMarkerStyle( file + 20 );

              hDistribution->SetMarkerColor( getColor( file+1 ) );

	      hDistribution->SetMarkerSize( 1.8 );
	      hDistribution->GetXaxis()->SetNdivisions( 504 );
	      hDistribution->GetXaxis()->SetTitle("E_{" + distribution + "}[GeV]");
	      hDistribution->GetYaxis()->SetTitle("d#sigma/dE [GeV]");

	      can_startingDistributions->cd();
	      Prepare_1Dplot( hDistribution );
	      hDistribution->DrawClone( drawoptions );
	      drawoptions = "phsame";
//	      legend->AddEntry( hDistribution, TString::Format( " #Delta#varphi_{max}=0.%i, #eta_{acc}=0.%i, " + legend_cuts[cuts_],
	      legend->AddEntry( hDistribution, TString::Format( legend_cuts[cuts_],  
		static_cast<int>(10. * deltaPhiMax[_phi] ), 
		static_cast<int>(10. * etaband[_eta] ) ), "lp" );
	      legend->SetFillColor( 0 );

	      } // Event selection.
	    } // Match.
	  } // Emin.
        } // Deltaphi.
      } // Eta.
   
  legend->Draw();

  can_startingDistributions->SetLogy();
  can_startingDistributions->SaveAs( TString::Format( "can_startingDistribution_" + distribution + "_Emin_%i.C", static_cast<int>( Ethresh_ ) ) );
  can_startingDistributions->SaveAs( TString::Format( "can_startingDistribution_" + distribution + "_Emin_%i.pdf", static_cast<int>( Ethresh_) ) );

    for(int _eta = 0; _eta < etaband.size(); _eta++){
      for(int _phi = 0; _phi < deltaPhiMax.size(); _phi++){
	for(int _Emin = 0; _Emin < Emin.size(); _Emin++, file++){

//          TString filename_ = TString::Format( "/user/avanspil/Castor_Analysis/ak5ak5_displaced_unfold_Emin_%f_deltaPhiMax_%f_etaband_%f.root", Emin[_Emin], deltaPhiMax[_phi], etaband[_eta] ) ;
		TString filename_ = TString::Format("/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo_unfold_Emin_150.000000_deltaPhiMax_%f_etaband_%f_all_matchE.root", deltaPhiMax[_phi], etaband[_eta] );
	  TFile* _file= TFile::Open( filename_, "read" );

	  RooUnfoldResponse *resp = (RooUnfoldResponse*)_file->Get("response");
	  TH2D* hResponse = (TH2D*)resp->Hresponse();
	  TH1D* hTruth = (TH1D*)resp->Htruth();
	  TH1D* hMeasured = (TH1D*)resp->Hmeasured();

	  double rFake = 1. - hResponse->Integral()/hMeasured->Integral() ;
	  double rMiss = 1. - hResponse->Integral()/hTruth->Integral() ;

	  cout << "\t" << deltaPhiMax[ _phi ] << "\t" << etaband[ _eta ] 
		<< "\tFraction fakes\t" << rFake
		<< "\tFraction misses\t" << rMiss << "\t\t" << ( 1. - rFake)*(1. + rMiss) << endl << endl << endl;

	}
      }
    } 
}









void PlotStartingDistributions_stacked(TString distribution){
  cout << "Unfolder::PlotStartingDistributions_stacked\t" << distribution << endl;

  double ticksize = gStyle->GetTickLength();
    
  gStyle->SetOptTitle( 0 );
  gStyle->SetOptStat( 0 );
    
  TString distributionname;
  if (distribution == "detector"){ 	distributionname = "hCastorJet_energy"; }
  else if (distribution == "generator"){ 	distributionname = "hGenJet_energy"; }
  else{ return; }
  
  TString drawoptions = "phist";

  int file = 0;

  vector<double> etaband;
  map<double, int> eta_markers;  
      etaband.push_back(0.0);

  vector<double> deltaPhiMax;
      deltaPhiMax.push_back(0.5);

  vector<TString> matching;
      matching.push_back("_matchE");

  vector<TString> match_symbol;
      match_symbol.push_back("E");

  vector<TString> model;
      model.push_back("");

  vector<double> Emin;
      Emin.push_back(150.);


  TString cuts_;
  std::vector<TString> vector_cuts_;
  std::map<TString, TString> legend_cuts;
  std::map<TString, TString> save_cuts;
        
  cuts_ = "_DefaultInclusive";
  vector_cuts_.push_back( cuts_ );
  legend_cuts[cuts_] = "Default selection";
  save_cuts[ cuts_ ] = "fwd11003";
  

  cuts_ = "_MostInclusive";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "BSC OR, HF OR, #leq 1 gVtx";        
  save_cuts[ cuts_ ] = "mostinclusive";


    for(int _eta = 0; _eta < etaband.size(); _eta++){
      for(int _phi = 0; _phi < deltaPhiMax.size(); _phi++){
	for(int _Emin = 0; _Emin < Emin.size(); _Emin++){
   	    for(int _match = 0; _match < matching.size(); _match++){
   	      for(int _cuts = 0; _cuts < vector_cuts_.size(); _cuts++, file++) {
   	        cuts_ = vector_cuts_[_cuts];
   	        
   	        TCanvas *can_startingDistributions;
  		PrepareCanvas(can_startingDistributions, "can_startingDistributions");
  		TLegend *legend = new TLegend( 0.55, 0.65, 1.-can_startingDistributions->GetTopMargin() - ticksize, 1.-can_startingDistributions->GetRightMargin() - ticksize);
  		legend->SetFillStyle( 0 );
  		legend->SetBorderSize( 0 );

		TString filename_ = TString::Format("/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo" + cuts_ + "_unfold_Emin_150.000000_deltaPhiMax_%f_etaband_%f_all_matchE.root", deltaPhiMax[_phi], etaband[_eta] );
	        TFile* _file= TFile::Open( filename_, "read" );
	        if( !_file ){ 
	          cout << "No\t" << filename_ << endl;
	          continue; 
	        }

                RooUnfoldResponse* response = (RooUnfoldResponse*)_file->Get("response");
	        TH2D* hResponse = (TH2D*)response->Hresponse();
	      
	        THStack* stack = stack = new THStack( 
	      	  TString::Format( "Stack_" +  distribution ), 
	      	  TString::Format( "Stack_" +  distribution ) );

	        //== Plot detector level distributions.
	        if( distribution == "detector"){
	      
	          TH1D* hFakes = (TH1D*)response->Hfakes();
	          TH1D* hMatched = (TH1D*)hResponse->ProjectionX();
	          TH1D* hMeasured = (TH1D*)response->Hmeasured();
	        
	          //== Fill and draw stack.
	          SetDnDx( hFakes );
	          hFakes->SetFillColor( getColor( 2 ) );
	          hFakes->GetXaxis()->SetTitle( "E_{det} [GeV]" );
	          hFakes->GetYaxis()->SetTitle( "dN/dE [1/GeV]" );
	          Prepare_1Dplot( hFakes );
	          SetDnDx( hMatched );
	          hMatched->SetFillColor( getColor( 3 ) );
	        	        
	          stack->Add( hFakes );
	          stack->Add( hMatched );	        
	    	
	          can_startingDistributions->cd();       
	        
	          stack->Draw("hist");
	          stack->GetXaxis()->SetTitle( "E_{det} [GeV]" );
	          stack->GetYaxis()->SetTitle( "dN/dE [1/GeV]" );
	        
		  Prepare_1Dplot( stack );		
	          stack->GetXaxis()->SetLabelSize( stack->GetXaxis()->GetLabelSize()*0.9 );
	          stack->GetYaxis()->SetLabelSize( stack->GetYaxis()->GetLabelSize()*0.9 );
	        
		  gPad->Update();	 
		       	        
	          //== Draw measured distribution.
	          SetDnDx( hMeasured );
	          hMeasured->SetMarkerStyle( 29 );
	          hMeasured->SetMarkerSize( hMeasured->GetMarkerSize()*2 );
	          hMeasured->Draw("psame");
	      
	          legend->SetHeader("Detector level");
	          legend->AddEntry( hMeasured, "Measured", "p");
	          legend->AddEntry( hMatched, "Matched jets", "f");
	          legend->AddEntry( hFakes, "Fake jets", "f");	      	      
	        }
	      
	        //== Plot generator level distributions.
	        else{
	      	            
	          TH1D* hGeneratorLevel = (TH1D*)response->Htruth();
	          TH1D* hMatched = (TH1D*)hResponse->ProjectionY();	        
	          TH1D* hFreeMisses = (TH1D*)_file->Get("hGenJet_energy_freemiss_bins");
	          if( !hFreeMisses ){
	            hFreeMisses = (TH1D*)hMatched->Clone("freemisses");
	            hFreeMisses->Scale(0.);
	          }
	        	        
	          TH1D* hMisses = (TH1D*)hGeneratorLevel->Clone("Misses");	
	          hMisses->Add( hMatched, -1.);
	          hMisses->Add( hFreeMisses, -1.);	
	          
	          //== Fill and draw stack.	  
	          SetDnDx( hFreeMisses );
	          hFreeMisses->SetFillColor( getColor(  4 ) );      
	          SetDnDx( hMisses );
	          hMisses->SetFillColor( getColor(  2 ) );     
	          SetDnDx( hMatched );
	          hMatched->SetFillColor( getColor(  3 ) ); 
	          
	          stack->Add( hMisses );
	          stack->Add( hFreeMisses );
	          stack->Add( hMatched );
	          
	          can_startingDistributions->cd();       
	        
	          stack->Draw("hist");
	          stack->GetXaxis()->SetTitle( "E_{gen} [GeV]" );
	          stack->GetYaxis()->SetTitle( "dN/dE [1/GeV]" );
	        
		  Prepare_1Dplot( stack );		
	          stack->GetXaxis()->SetLabelSize( stack->GetXaxis()->GetLabelSize()*0.9 );
	          stack->GetYaxis()->SetLabelSize( stack->GetYaxis()->GetLabelSize()*0.9 );
	          
	          double minvalue = GetMinimumValue( stack );
	          double maxvalue = stack->GetMaximum();
	          
//	          stack->GetYaxis()->SetRangeUser( minvalue * 0.5, maxvalue * 1.1 );
		  stack->GetHistogram()->GetYaxis()->SetRangeUser( 1e-4, 1e4 );
	        
		  gPad->Update();	 
		       	        
	          //== Draw measured distribution.
	          SetDnDx( hGeneratorLevel );
	          hGeneratorLevel->SetMarkerStyle( 29 );
	          hGeneratorLevel->SetMarkerSize( hGeneratorLevel->GetMarkerSize()*2 );
	          hGeneratorLevel->Draw("psame");
	      
	          legend->SetHeader("Generator level");
	          legend->AddEntry( hGeneratorLevel, "Truth", "p");
	          legend->AddEntry( hMatched, "Matched jets", "f");
	          legend->AddEntry( hMisses, "Missed jets (passed events)", "f");	
	          legend->AddEntry( hFreeMisses, "Missed jets (unpassed events)", "f");	
	          	          	         	          	          
	                
	        }
	        	           
  		legend->Draw();

  		can_startingDistributions->SetLogy();
		can_startingDistributions->SaveAs( TString::Format( "can_startingDistribution_" + distribution + "_" + save_cuts[ cuts_ ] + "_stacked.C" ) );
  		can_startingDistributions->SaveAs( TString::Format( "can_startingDistribution_" + distribution + "_" + save_cuts[ cuts_ ] + "_stacked.pdf" ) );

	      } // Event selection.
	    } // Match.
	  } // Emin.
        } // Deltaphi.
      } // Eta.

}





























