 //STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPad.h>
#include <TColor.h>
#include <TString.h>

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

using namespace std;

int main(){

  gStyle->SetOptTitle( 0 );
  setTDRStyle();
  gStyle->SetPalette(1);
  gStyle->SetErrorX(0);

  //== Canvas.
  TString canvastitle = "can_comparison_diffractive";
  TCanvas *can_;  
  PrepareCanvas( can_, canvastitle );
  TPad* pad_abs_, *pad_ratio_;
  SplitCanvas( can_, pad_abs_, pad_ratio_ );
  can_->SetLogy();  
  double ticksize = gStyle->GetTickLength(); 
  double xwidth = ( (1. - pad_abs_->GetRightMargin() - ticksize) - 0.3)/2.;

  //== MPI.
  TString MPI_;
  vector<TString> vector_MPI_;
  map<TString, TString> draw_MPI_;
  map<TString, TString> legends_MPI_;
  map<TString, TLegend*> legendBox_MPI_;

  MPI_ = "_on";
  vector_MPI_.push_back( MPI_ );
  legends_MPI_[ MPI_ ] = "MPI";    
    TLegend *leg1 = new TLegend( 0.3, 0.5, 0.3+xwidth, 1. - pad_abs_->GetTopMargin() - ticksize);
    leg1->SetFillStyle( 0 );
    leg1->SetBorderSize( 0 ); 
  legendBox_MPI_[ MPI_ ] = leg1;
  draw_MPI_[ MPI_ ] = "hist";

  MPI_ = "_off";
  vector_MPI_.push_back( MPI_ );
  legends_MPI_[ MPI_ ] = "no MPI";
    TLegend *leg2= new TLegend(  0.3+xwidth, 0.5, 1. - pad_abs_->GetRightMargin() - ticksize, 1. - pad_abs_->GetTopMargin() - ticksize);
    leg2->SetFillStyle( 0 );
    leg2->SetBorderSize( 0 ); 
  legendBox_MPI_[ MPI_ ] = leg2;    
  draw_MPI_[ MPI_ ] = "p";

  //== Diffraction.
  TString diff_;
  vector<TString> vector_diff_;
  map<TString, TString> legends_diff_;
  map<TString, double> xsec_diff_;
  map<TString, double> max_ratio;  
  
  diff_ = "_ND";
  vector_diff_.push_back( diff_ );
  legends_diff_[ diff_ ] = "ND";
  xsec_diff_[ diff_ ] = 50.9;
  max_ratio[ diff_ ] = 1.34;

  diff_ = "_SD";
  vector_diff_.push_back( diff_ );
  legends_diff_[ diff_ ] = "SD";
  xsec_diff_[ diff_ ] = 2.*6.19; 
  max_ratio[ diff_ ] = 4.74; 

  diff_ = "_DD";
  vector_diff_.push_back( diff_ );
  legends_diff_[ diff_ ] = "DD";
  xsec_diff_[ diff_ ] = 8.10;  
  max_ratio[ diff_ ] = 7.74;

  diff_ = "_NSDD";
  vector_diff_.push_back( diff_ );
  legends_diff_[ diff_ ] = "ND/SD/DD";
  xsec_diff_[diff_] = 0.;
  max_ratio[ diff_ ] = 1.34;

  //== Plots.
  TString distr_;
  vector<TString> vector_distr_;
  map<TString, TString> legends_distr_;  
  map<TString, int> color_index;
  map<TString, THStack*> stack_distr_;
  THStack* stack;  

  
  distr_ = "hHF_energy_and_tracker_MPI_";
  vector_distr_.push_back( distr_ );
  legends_distr_[ distr_ ] = "AND + vertex";
  color_index[distr_] = 9; 
  stack = new THStack( TString::Format( "Stack" + distr_ ), "" );   
  stack_distr_[ distr_ ] = stack;  
  
  distr_ = "hHF_energy_and_MPI_";
  vector_distr_.push_back( distr_ );
  legends_distr_[ distr_ ] = "and";
  color_index[distr_] = 2;    
  stack = new THStack( TString::Format( "Stack" + distr_ ), "" );
  stack_distr_[ distr_ ] = stack;
  
  distr_ = "hHF_energy_or_MPI_";
  vector_distr_.push_back( distr_ );
  legends_distr_[ distr_ ] = "or";   
  color_index[distr_] = 1;
  stack = new THStack( TString::Format( "Stack" + distr_ ), "" );
  stack_distr_[ distr_ ] = stack;
  /*
  distr_ = "hHF_energy_minus_MPI_";
  vector_distr_.push_back( distr_ );
  legends_distr_[ distr_ ] = "HF (-)";
  color_index[distr_] = 3; 
  stack = new THStack( TString::Format( "Stack" + distr_ ), "" );
  stack_distr_[ distr_ ] = stack;
     
  distr_ = "hHF_energy_xor_minus_MPI_";
  vector_distr_.push_back( distr_ );
  legends_distr_[ distr_ ] = "HF (xor,-)"; 
  color_index[distr_] = 4;
  stack = new THStack( TString::Format( "Stack" + distr_ ), "" );
  stack_distr_[ distr_ ] = stack;
  
  distr_ = "hHF_energy_plus_MPI_";
  vector_distr_.push_back( distr_ );
  legends_distr_[ distr_ ] = "HF (+)";
  color_index[distr_] = 5;
  stack = new THStack( TString::Format( "Stack" + distr_ ), "" );
  stack_distr_[ distr_ ] = stack;
  
  distr_ = "hHF_energy_xor_plus_MPI_";
  vector_distr_.push_back( distr_ );
  legends_distr_[ distr_ ] = "HF (xor,+)";
  color_index[distr_] = 6;   
  stack = new THStack( TString::Format( "Stack" + distr_ ), "" );   
  stack_distr_[ distr_ ] = stack;
   
  distr_ = "hHF_energy_xor_MPI_";
  vector_distr_.push_back( distr_ );
  legends_distr_[ distr_ ] = "HF (xor)"; 
  color_index[distr_] = 7;     
  stack = new THStack( TString::Format( "Stack" + distr_ ), "" );
  stack_distr_[ distr_ ] = stack;
  */
  
  //== N events in file.
  vector<int> nEvents;
  map<int, TString> event_legend;  
  int nev;
  
  nev = 10000000;  
  nEvents.push_back(nev);      
  event_legend[nev] = "(10^{7})";
  
  nev = 1000000;  
  nEvents.push_back(nev);      
  event_legend[nev] = "(10^{6})";
  
  nev = 100000;  
  nEvents.push_back(nev);      
  event_legend[nev] = "(10^{5})";

  nev = 10000;  
  nEvents.push_back(nev);      
  event_legend[nev] = "(10^{4})";
 
  
  //== Loop. 

  for(int d_ = 0 ; d_ < vector_diff_.size(); d_++){ 
    diff_ = vector_diff_[ d_ ];       
    TString drawoptions = "hist";   
    double maximum_value = 0.;  
    double y2_legend;   
    TH1D* hFirst; 
    bool first = true;        
    
    //== Loop.
    for(int c_ = 0; c_ < vector_MPI_.size(); c_++){
      MPI_ = vector_MPI_[ c_ ];
      TString MPI_hist = MPI_;
      MPI_hist.ReplaceAll( "_", "");
      drawoptions.ReplaceAll("hist", draw_MPI_[MPI_]);
      
      TLegend *leg = legendBox_MPI_[ MPI_ ];      
      leg->Clear();
      leg->SetHeader( legends_diff_[ diff_ ] + ", MPI " + MPI_hist);

      //== Open data, extract histogram.
      int event_counter = 0;      

      for(int h_ = 0; h_ < vector_distr_.size(); h_++){
        distr_ = vector_distr_[h_]; 
        TH1D* _datahist;
        TFile* _datafile;        
        event_counter = 0;
        bool good_hist = false;
        
        TString dataname;
        TString histname = distr_ + MPI_hist + diff_;
        
        for( int n_ = 0; n_ < nEvents.size(); n_++){
          nev = nEvents[n_];
          dataname = TString::Format("/user/avanspil/Pythia_Fastjet_files/20160603_Tune4C_MPI" + MPI_ + diff_ + "_Tree_%ievents_presetTune.root", nev);
          _datafile = TFile::Open( dataname , "Read" );
          if( !_datafile ){ continue; }
          
          _datahist = (TH1D*)_datafile->Get( histname );
          if( _datahist ){
            if( _datahist->Integral() == _datahist->Integral() ){ break; }	//== Found the histogram.
            else{ continue; }
          }        
        }

	if( !_datafile ){ cout << "\t\tNo\t" << dataname << endl; }
        if( !_datahist ){ continue; }     
        
        SetDnDx( _datahist );        
        _datahist->Scale( 1./double(nev) );                
        
        if( maximum_value < _datahist->GetMaximum() ){  maximum_value = _datahist->GetMaximum(); }        
        if( first ){ hFirst = (TH1D*)_datahist->Clone("First");  first = false;}

	int _color = color_index[distr_];
	//== Prepare plots.
        _datahist->SetLineWidth( 2 );
        _datahist->SetMarkerColor( getColor( _color ) );
        _datahist->SetMarkerStyle( _color + 19 );        
        _datahist->SetLineColor( getColor( _color ) );        
        _datahist->SetLineStyle( _color );
        _datahist->GetXaxis()->SetTitle("E [GeV]");
        _datahist->GetYaxis()->SetTitle("1/N_{tot.} dN/dE");
        
        TString legend_filler = "l";
        if( draw_MPI_[MPI_] == "p" ){ legend_filler = "p"; }
        //leg->AddEntry( _datahist, legends_distr_[ distr_ ] + " " + event_legend[nev]  , legend_filler);
        leg->AddEntry( _datahist, legends_distr_[ distr_ ]   , legend_filler);

 	//== Absolute values.
	pad_abs_->cd();
        Prepare_1Dplot( _datahist );
        _datahist->GetYaxis()->SetRangeUser(2.e-10, 5.);        
        _datahist->DrawClone( drawoptions );
        
        //== Relative values.
        pad_ratio_->cd();
        Prepare_1Dplot_ratio( _datahist );
        _datahist->Divide( hFirst );
        _datahist->GetXaxis()->SetLabelSize( _datahist->GetXaxis()->GetLabelSize()*0.75 );
        _datahist->GetYaxis()->SetRangeUser(0.,  max_ratio[ diff_ ] );
        _datahist->GetYaxis()->SetTitle("Ratio");
        _datahist->GetYaxis()->CenterTitle();
        _datahist->DrawClone( drawoptions );
  
        drawoptions = draw_MPI_[MPI_] +"same";
        event_counter = 0;
        
        //== Prepare stack.
        if( diff_ != "NSDD" ){
          
        }
        
      }
      pad_abs_->cd();      
      leg->Draw();      
    }    
   
    pad_abs_->SetLogy();

    can_->SaveAs( "can_comparison_energy" + diff_ + ".pdf" );
     
  } //== Loop over use of same/different MC sample.


       



  return 0;
}







