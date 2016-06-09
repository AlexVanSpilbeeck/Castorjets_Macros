#include "iostream"
#include "fstream"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TColor.h"
#include "TMath.h"
#include "TApplication.h"

#include "stdio.h"
#include "cstring"
#include <vector>
#include <sys/stat.h>

#include "../tdrstyle.C"

#include "../color.h"

#include "../Functions/Function_PrepareCanvas.h"
#include "../Functions/Function_Prepare1Dplot.h"
#include "../Functions/Function_DivideByBinWidth.h"
//#include "../Functions/Function_SetRangeToIncludeAll.h"

using namespace std; 


















void PlotStartingDistributions(TString distribution, double deltaPhiMax, double etaacc, TString model){

  TCanvas *can_startingDistributions;
  PrepareCanvas(can_startingDistributions, "can_startingDistributions_" + distribution);


  TLegend *legend = new TLegend( 
	1. - can_startingDistributions->GetRightMargin() - 0.42,
	0.65,
	1. - can_startingDistributions->GetRightMargin(),
	1. - can_startingDistributions->GetTopMargin());

  legend->SetFillStyle( 0 );
  legend->SetBorderSize( 0 );

  TString distributionname_unmatched, distributionname_all;
  if (distribution == "fake"){	distributionname_unmatched = "hCastorJet_fake_all"; distributionname_all = "hCastorJet_energy"; }
  if (distribution == "miss"){	distributionname_unmatched = "hCastorJet_miss_all"; distributionname_all = "hGenJet_energy"; }
  if (distribution == "meas"){	distributionname_unmatched = "hCastorJet_energy"; }
  if (distribution == "true"){	distributionname_unmatched = "hGenJet_energy"; }
  

  std::map<TString, TString> yaxis;
   yaxis["fake"] 		= "dN_{fake}/dE [1/Gev]";
   yaxis["miss"] 		= "dN_{miss}/dE [1/Gev]";
   yaxis["meas"] 		= "dN_{det.}/dE [1/Gev]";
   yaxis["true"] 		= "dN_{gen.}/dE [1/Gev]";

  TString drawoptions = "phist";

  int file = 0;

  vector<TString> cuts;
    cuts.push_back("");
    cuts.push_back("_noHFcut");
    cuts.push_back("_noVertexcut");
    cuts.push_back("_noHFcut_noVertexcut");

  std::map<TString, TString> legends;
   legends[""] 			= "Vtx, HF";
   legends["_noHFcut"] 		= "Vtx";
   legends["_noVertexcut"]	= "HF";
   legends["_noHFcut_noVertexcut"] = "";

  vector<int> markers;
      markers.push_back( 20 );
      markers.push_back( 22 );
      markers.push_back( 29 );
      markers.push_back( 21 );

  vector<TH1D*> histlist;
  TH1D* hFirst;
  int setup = 0;

  for(int _cut = 0; _cut < cuts.size(); _cut++, file++){

  cout << "model\t" << model;
  cout << "cut\t" << cuts[_cut] << endl;
  cout << "DeltaPhi\t" << deltaPhiMax << endl;
  cout << "eta\t" << etaacc << endl;

	    TString filename_ = TString::Format( 
		"/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_" 
		+ model + "_NewGeo" 
		+ cuts[_cut]  + "_unfold_Emin_150.000000_deltaPhiMax_%f_etaband_%f_all_matchE.root", deltaPhiMax, etaacc ) ;
	    
	      TFile* _file= TFile::Open( filename_, "read" );

	      if( !_file ) continue;	//== Continue if file does not exist.

	      TString setup_label = TString::Format( model + "displaced_unfold_Emin_150.000000_deltaPhiMax_0%i_etaband_0%i" + cuts[_cut], 
		static_cast<int>(150), 
		static_cast<int>(10. * deltaPhiMax), 
		static_cast<int>(10. * etaacc));

	      TString histname = TString::Format( distributionname_unmatched + "_Emin_150.000000_deltaPhiMax_0%i_etaband_0%i" + cuts[_cut] +  "_" + model, 
		static_cast<int>(150), 
		static_cast<int>(10. * deltaPhiMax), 
		static_cast<int>(10. * etaacc));

    	      //== Extracting the distribution.
  	      TH1D* hUnmatched = (TH1D*)_file->Get( distributionname_unmatched );

	      //== Scaling.
              int total_events_nocuts;

	      SetDnDx( hUnmatched );
	      hUnmatched->SetTitle( histname );
	      hUnmatched->SetName( histname );
	      hUnmatched->SetLineWidth( 3 );
	      hUnmatched->SetLineColor( getColor( _cut+1 ) ); 
	      hUnmatched->SetLineStyle( file + 1 ); 
	      hUnmatched->SetLineColor( getColor(file + 1) ); 
              hUnmatched->SetMarkerStyle( file + 20 );
              hUnmatched->SetMarkerColor( getColor( file+1 ) );
	      hUnmatched->SetMarkerSize( 1.8 );
	      hUnmatched->GetXaxis()->SetNdivisions( 504 );
	      hUnmatched->GetXaxis()->SetTitle("E [GeV]");
	      hUnmatched->GetYaxis()->SetTitle( yaxis[distribution] );

	     

	      can_startingDistributions->cd();

	      Prepare_1Dplot( hUnmatched );

	      if( setup == 0 ){
		hFirst = (TH1D*)hUnmatched->Clone("First");
		hFirst->Draw("hist");
		setup++;	
	      }
	      else{
	        hUnmatched->DrawClone( "histsame" );
	      }
	      histlist.push_back( hUnmatched );

//	      SetRangeToIncludeAll( hFirst, histlist );
//	      hFirst->GetYaxis()->SetRangeUser(1e-5, 1.);
	      hFirst->GetYaxis()->SetRangeUser(
		hFirst->GetMinimum()*0.9,
		hFirst->GetMaximum()*1.1);

	      legend->AddEntry( hUnmatched, TString::Format( " #Delta#varphi_{max}=0.%i, #eta_{acc}=0.%i, " + legends[cuts[_cut]], 
		static_cast<int>(10. * deltaPhiMax ), 
		static_cast<int>(10. * etaacc ) ), "lp" );
	      legend->SetFillColor( 0 );

  } 
   
  legend->Draw();

  can_startingDistributions->SetLogy();
  can_startingDistributions->SaveAs( TString::Format(  "can_absolute_" + distribution + "_" + model + "_Emin_150.C") );
  can_startingDistributions->SaveAs( TString::Format(  "can_absolute_" + distribution + "_" + model + "_Emin_150.pdf" ) );

}

























void PlotStartingDistributions_ratio(TString distribution, double deltaPhiMax, double etaacc, TString model){

  TCanvas *can_startingDistributions;
  PrepareCanvas(can_startingDistributions, "can_startingDistributions_" + distribution);

  TLegend *legend = new TLegend( 
	1. - can_startingDistributions->GetRightMargin() - 0.42,
	0.65,
	1. - can_startingDistributions->GetRightMargin(),
	1. - can_startingDistributions->GetTopMargin());
  legend->SetFillStyle( 0 );
  legend->SetBorderSize( 0 );

  TString distributionname_unmatched, distributionname_all;
  if (distribution == "fake"){ 		distributionname_unmatched = "hCastorJet_fake_all"; distributionname_all = "hCastorJet_energy"; }
  if (distribution == "miss"){ 		distributionname_unmatched = "hCastorJet_miss_all"; distributionname_all = "hGenJet_energy"; }

  std::map<TString, TString> yaxis;
   yaxis["fake"] 		= "N_{fake}/N_{meas.}";
   yaxis["miss"] 		= "N_{miss}/N_{true}";

  TString drawoptions = "phist";

  int file = 0;

  vector<TString> cuts;
    cuts.push_back("");
    cuts.push_back("_noHFcut");
    cuts.push_back("_noVertexcut");
    cuts.push_back("_noHFcut_noVertexcut");

  std::map<TString, TString> legends;
   legends[""] 			= "Vtx, HF";
   legends["_noHFcut"] 		= "Vtx";
   legends["_noVertexcut"]	= "HF";
   legends["_noHFcut_noVertexcut"] = "";



  vector<int> markers;
      markers.push_back( 20 );
      markers.push_back( 22 );
      markers.push_back( 29 );
      markers.push_back( 21 );

  vector<TH1D*> histlist;
  TH1D* hFirst;
  int setup = 0;

  for(int _cut = 0; _cut < cuts.size(); _cut++, file++){

  cout << "model\t" << model;
  cout << "cut\t" << cuts[_cut] << endl;
  cout << "DeltaPhi\t" << deltaPhiMax << endl;
  cout << "eta\t" << etaacc << endl;

	    TString filename_ = TString::Format( 
		"/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_" 
		+ model + "_NewGeo" 
		+ cuts[_cut]  + "_unfold_Emin_150.000000_deltaPhiMax_%f_etaband_%f_all_matchE.root", deltaPhiMax, etaacc ) ;

cout << "Filename is\t" << filename_ << endl;
	    
	      TFile* _file= TFile::Open( filename_, "read" );

	      if( !_file ) continue;	//== Continue if file does not exist.

	      TString setup_label = TString::Format( model + "displaced_unfold_Emin_150.000000_deltaPhiMax_0%i_etaband_0%i" + cuts[_cut], 
		static_cast<int>(150), 
		static_cast<int>(10. * deltaPhiMax), 
		static_cast<int>(10. * etaacc));

	      TString histname = TString::Format( distributionname_unmatched + "_Emin_150.000000_deltaPhiMax_0%i_etaband_0%i" + cuts[_cut] +  "_" + model, 
		static_cast<int>(150), 
		static_cast<int>(10. * deltaPhiMax), 
		static_cast<int>(10. * etaacc));

    	      //== Extracting the distribution.
  	      TH1D* hUnmatched = (TH1D*)_file->Get( distributionname_unmatched );
  	      TH1D* hAll = (TH1D*)_file->Get( distributionname_all );

	      //== Scaling.
              int total_events_nocuts;

	      hUnmatched->SetTitle( histname );
	      hUnmatched->SetName( histname );
	      hUnmatched->SetLineWidth( 3 );

	      // Linecolor = model -- linestyle = match -- markerstyle = etaband -- markercolor = deltaphi
	      hUnmatched->SetLineColor( getColor( _cut+1 ) ); 

	      hUnmatched->SetLineStyle( file + 1 ); 
	      hUnmatched->SetLineColor( getColor(file + 1) ); 
 
              hUnmatched->SetMarkerStyle( file + 20 );

              hUnmatched->SetMarkerColor( getColor( file+1 ) );

	      hUnmatched->SetMarkerSize( 1.8 );
	      hUnmatched->GetXaxis()->SetNdivisions( 504 );
	      hUnmatched->GetXaxis()->SetTitle("E_{" + distribution + "}[GeV]");
	      hUnmatched->GetYaxis()->SetTitle( yaxis[distribution] );

	      can_startingDistributions->cd();

	      hUnmatched->Divide( hAll );
	      hUnmatched->GetYaxis()->SetRangeUser(hUnmatched->GetMinimum()*0.9, 1.);
	      Prepare_1Dplot( hUnmatched );

	      if( setup == 0 ){
		hFirst = (TH1D*)hUnmatched->Clone("First");
		hFirst->Draw("hist");
		setup++;	
	      }
	      else{
	        hUnmatched->DrawClone( "histsame" );
	      }
	      histlist.push_back( hUnmatched );

	      //SetRangeToIncludeAll( hFirst, histlist );
	      hFirst->GetYaxis()->SetRangeUser(hFirst->GetMinimum()*0.9, 1.);

	      legend->AddEntry( hUnmatched, TString::Format( " #Delta#varphi_{max}=0.%i, #eta_{acc}=0.%i, " + legends[cuts[_cut]], 
		static_cast<int>(10. * deltaPhiMax ), 
		static_cast<int>(10. * etaacc ) ), "lp" );
	      legend->SetFillColor( 0 );

  } 
   
  legend->Draw();

  can_startingDistributions->SetLogy();
  can_startingDistributions->SaveAs( TString::Format(  "can_ratio_" + distribution + "_" + model + "_Emin_150.C") );
  can_startingDistributions->SaveAs( TString::Format(  "can_ratio_" + distribution + "_" + model + "_Emin_150.pdf" ) );

}



int main(){

  gStyle->SetOptTitle(0);

  PlotStartingDistributions_ratio("fake", 0.5, 0.0, "Pythia6Z2star");
  PlotStartingDistributions_ratio("miss", 0.5, 0.0, "Pythia6Z2star");

  PlotStartingDistributions_ratio("fake", 0.5, 0.0, "Pythia84C");
  PlotStartingDistributions_ratio("miss", 0.5, 0.0, "Pythia84C");

  PlotStartingDistributions_ratio("fake", 0.5, 0.0, "EPOS");
  PlotStartingDistributions_ratio("miss", 0.5, 0.0, "EPOS");

  PlotStartingDistributions("fake", 0.5, 0.0, "Pythia6Z2star");
  PlotStartingDistributions("miss", 0.5, 0.0, "Pythia6Z2star");
  PlotStartingDistributions("meas", 0.5, 0.0, "Pythia6Z2star");
  PlotStartingDistributions("true", 0.5, 0.0, "Pythia6Z2star");

  PlotStartingDistributions("fake", 0.5, 0.0, "Pythia84C");
  PlotStartingDistributions("miss", 0.5, 0.0, "Pythia84C");
  PlotStartingDistributions("meas", 0.5, 0.0, "Pythia84C");
  PlotStartingDistributions("true", 0.5, 0.0, "Pythia84C");

  PlotStartingDistributions("fake", 0.5, 0.0, "EPOS");
  PlotStartingDistributions("miss", 0.5, 0.0, "EPOS");
  PlotStartingDistributions("meas", 0.5, 0.0, "EPOS");
  PlotStartingDistributions("true", 0.5, 0.0, "EPOS");

 return 0;

}
