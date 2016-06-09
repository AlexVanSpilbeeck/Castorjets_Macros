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
#include <TLine.h>
#include <TMath.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <THnSparse.h>
#include <TThread.h>
#include <TStopwatch.h>
#include <TLegend.h>
#include <TPad.h>
#include <TColor.h>
#include <TString.h>
#include <TText.h>
#include <TObjString.h>

#include <TLine.h>
#include <TPaletteAxis.h>

//STANDARD C++ INCLUDES
#include <sys/stat.h>
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


#include "../tdrstyle.C"
#include "../color.h"

#include "../Functions/Function_DivideByBinWidth.h"
#include "../Functions/Function_PrepareCanvas.h"
#include "../Functions/Function_Prepare1Dplot.h"
#include "../Functions/Function_Prepare2Dplot.h"

#include "../../../../RooUnfold-1.1.1/src/RooUnfold.h"
#include "../../../../RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"
#include "../../../../RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "../../../../RooUnfold-1.1.1/src/RooUnfoldResponse.h"


using namespace std;

int main(){

  gStyle->SetOptTitle( 0 );
  setTDRStyle();
  gStyle->SetPalette(1);

  //== MPI.
  TString MPI_;
  vector<TString> vector_MPI_;
  map<TString, TString> legends_MPI_;

  MPI_ = "_off";
  vector_MPI_.push_back( MPI_ );
  legends_MPI_[ MPI_ ] = "no MPI";

  MPI_ = "_on";
  vector_MPI_.push_back( MPI_ );
  legends_MPI_[ MPI_ ] = "MPI";

  //== Diffraction.
  TString diff_;
  vector<TString> vector_diff_;
  map<TString, TString> legends_diff_;

  diff_ = "_ND";
  vector_diff_.push_back( diff_ );
  legends_diff_[ diff_ ] = "ND";

  diff_ = "_SD";
  vector_diff_.push_back( diff_ );
  legends_diff_[ diff_ ] = "SD";

  diff_ = "_DD";
  vector_diff_.push_back( diff_ );
  legends_diff_[ diff_ ] = "DD";

  //== Unfolded level canvas. 
  TString canvastitle = "can_comparison_diffractive";
  TCanvas *can_;  
  PrepareCanvas_flat( can_, canvastitle );
    
  TLegend *leg = new TLegend( 0.4, 0.7, 0.7, 1. - can_->GetTopMargin() - gStyle->GetTickLength());
  leg->SetFillStyle( 0 );
  leg->SetBorderSize( 0 );    
    
  TString drawoptions = "chist";   
  double maximum_value = 0.;  
  int _color = 2;
  double y2_legend;

  for(int d_ = 0 ; d_ < vector_diff_.size(); d_++){ 
    diff_ = vector_diff_[ d_ ];

    //== Loop.
    for(int c_ = 0; c_ < vector_MPI_.size(); c_++){
      MPI_ = vector_MPI_[ c_ ];

      //== Open data, extract histogram.
      TFile* _datafile = TFile::Open( "/user/avanspil/Pythia_Fastjet_files/20160527_Tune4C_MPI" + MPI_ +  diff_  + "_softQCD_1000000events_presetTune.root", "read" );
      if( !_datafile ){ continue; }
      TH1D* _datahist = (TH1D*)_datafile->Get("hGenJet_eta" + MPI_ + diff_ );
      if( !_datahist ){ cout << "No\t" << diff_ << "\t" << MPI_ << endl; continue; }
      
      if( maximum_value < _datahist->GetMaximum() ){  maximum_value = _datahist->GetMaximum(); }

      //== Unfolded.
      _datahist->SetLineColor( getColor( d_ + 2 ) );
      _datahist->SetLineStyle( c_ +1  );
      _datahist->SetLineWidth( 2 );
      _datahist->GetXaxis()->SetTitle("#eta");
      _datahist->GetYaxis()->SetTitle("dN/d#eta");
      Prepare_1Dplot_flat( _datahist );
      leg->AddEntry( _datahist, legends_diff_[ diff_ ] + " " + legends_MPI_[ MPI_ ] , "l");
      _datahist->Draw( drawoptions );
  
      drawoptions = "csame";
      
    } 
  } //== Loop over use of same/different MC sample.

  TLine *linetracker = new TLine(-2.5, 0., -2.5, maximum_value*1.1 );
  linetracker->SetLineColor( kBlack );
  linetracker->SetLineWidth( 2 );
  linetracker->Draw();
  
  linetracker = new TLine(2.5, 0., 2.5, maximum_value*1.1 );
  linetracker->SetLineColor( kBlack );
  linetracker->SetLineWidth( 2 );
  linetracker->Draw();  

  can_->SetLogy();
  leg->Draw();
  can_->SaveAs( "can_comparison_diffractions.pdf" );

  return 0;
}
