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

using namespace std;

int main(){

  gStyle->SetOptTitle( 0 );
  setTDRStyle();
  gStyle->SetPalette(1);
  
  int  new_dir = mkdir( ( TString::Format("Plots_compareEventSelection") ).Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  );  


  //== Possible cuts on original file.

  TString cuts_;
  std::vector<TString> vector_cuts_;
  std::map<TString, TString> legend_cuts;
/*
  cuts_ = "";
  vector_cuts_.push_back( cuts_ );
  legend_cuts[cuts_] = "Default selection";
*/
  cuts_ = "_BSC_OR";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "BSC OR, == 1 gVtx";
  /*
  cuts_ = "_noVertex_HF";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "HF, 0 gVtx";

  cuts_ = "_noVertex_BSC";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "BSC OR, 0 gVtx";
  
  cuts_ = "_atLeastOneVertex_HF";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "HF, #geq 1 gVtx";
  */
  cuts_ = "_atLeastOneVertex_BSC";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "BSC OR, #geq 1 gVtx";
  /*
  cuts_ = "_one_none_vertex_BSC";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "BSC OR, #leq 1 gVtx";  
  */
  cuts_ = "_one_none_vertex_BSCor_HFor";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "BSC OR, #leq 1 gVtx";    


  //== Source of sample.

  TString type_;  
  std::vector<TString> vector_type_;
  std::map<TString, TString> legends_type_;
  std::map<TString, TString> match_type_;
  std::map<TString, TString> suffix_type_;
  std::map<TString, TString> save_;
  std::map<TString, TString> drawType_;

  type_ = "_data";
  vector_type_.push_back( type_ );
  legends_type_[ type_ ] = "Data";
  suffix_type_[ type_ ] = "";
  match_type_[ type_ ] = "";
  save_[ type_ ] = "Data";
  drawType_[ type_ ] = "p";

  type_ = "_Pythia6Z2star_NewGeo";
  vector_type_.push_back( type_ );
  legends_type_[ type_ ] = "Pythia6 (Z2*), New Geo.";
  match_type_[ type_ ] = "deltaPhiMax_0.500000_etaband_0.000000_";
  suffix_type_[ type_ ] = "_matchE";
  save_[ type_ ] = "Pythia6Z2star_NewGeo";
  drawType_[ type_ ] = "l";
/*
  type_ = "_Pythia84C_NewGeo";
  vector_type_.push_back( type_ );
  legends_type_[ type_ ] = "Pythia8 (4C), New Geo.";
  match_type_[ type_ ] = "deltaPhiMax_0.500000_etaband_0.000000_";
  suffix_type_[ type_ ] = "_matchE";
  save_[ type_ ] = "Pythia84C_NewGeo";
  drawType_[ type_ ] = "l";
 
  type_ = "_EPOS_NewGeo";
  vector_type_.push_back( type_ );
  legends_type_[ type_ ] = "EPOS-LHC, New Geo.";
  match_type_[ type_ ] = "deltaPhiMax_0.500000_etaband_0.000000_";
  suffix_type_[ type_ ] = "_matchE";
  save_[ type_ ] = "EPOS_NewGeo";
  drawType_[ type_ ] = "l";
*/
  //== Plots.

  TString plots;
  std::vector<TString> histograms;
  std::map<TString, TString> xtitle;
  std::map<TString, TString> ytitle;
    
  plots = "hCastorJet_energy" ;
  histograms.push_back( plots );
  xtitle[ plots ] = "E_{cal} [GeV]";
  ytitle[ plots ] = "dN/dE [1/GeV]";

  for(int t_ = 0; t_ < vector_type_.size(); t_++ ){
    type_ = vector_type_[t_];

  for( int p_ = 0; p_ < histograms.size(); p_++){
    plots = histograms[ p_ ];

    //== Canvas.
    TCanvas *can_;
    TPad* pad_abs_, *pad_ratio_;
    TString can_title = plots ;
    PrepareCanvas( can_, plots );
    SplitCanvas( can_, pad_abs_, pad_ratio_ );
    TString drawoptions = "ephist";
    
    TH1D* hFirst;
    
    TLegend* leg = new TLegend(
      0.65,
      0.4,
      1. - pad_abs_->GetRightMargin()-0.05,      
      1. - pad_abs_->GetTopMargin()-0.05
    );
    leg->SetFillStyle( 0 );
    leg->SetBorderSize( 0 );
    leg->SetHeader("Measured");

    //== Loop over cuts in sample files.
    for( int c_ = 0; c_ < vector_cuts_.size(); c_++){
      cuts_ = vector_cuts_[ c_ ];
            
      cout << "cut\t" << cuts_ << "\t" << type_ << endl;
    
      TFile* _file;
      TString filename;
      if( type_ == "_data" ){
        filename = TString::Format( "/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5_data" + cuts_ + "_unfold_Emin_150.000000_all.root" );
        _file = TFile::Open( filename , "Read" ); }
      else{
        filename = TString::Format("/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5" + type_ + cuts_ + "_unfold_Emin_150.000000_deltaPhiMax_0.500000_etaband_0.000000_all_matchE.root");
        _file = TFile::Open( filename, "Read" ); }
 
      if( !_file ){ 
      continue; }
 
      //== Extraction.
      TH1D* _h = (TH1D*)_file->Get( plots );
      if( !_h ){ 
        continue; }
      if( _h->Integral() != _h->Integral() ) continue; 
      _h->Sumw2();      
      
      SetDnDx( _h );
      
      if( c_ == 0 ){
        hFirst = (TH1D*)_h->Clone("First");
      }

      //== Draw.
      Prepare_1Dplot_ratio( _h ) ;
      _h->SetLineColor( getColor( c_ + 1 ) );
      _h->SetLineStyle(  c_ + 1  );
      _h->SetLineWidth(  2  );
      _h->SetMarkerColor( getColor( c_ + 1 ) ); 
      _h->SetMarkerStyle(  c_ + 20  );
      _h->GetXaxis()->SetRangeUser(150., 1600.);
      _h->SetMarkerSize( 1.2 * _h->GetMarkerSize() );   
      _h->GetXaxis()->SetTitle( xtitle[plots]  );
      _h->GetYaxis()->SetTitle( ytitle[plots]  );
      
      pad_abs_->cd();     
      _h->DrawClone( drawoptions );
      pad_abs_->SetLogy();           
      
      pad_ratio_->cd();
      _h->Divide( hFirst );
      _h->GetYaxis()->SetRangeUser(0.95, 1.36);
      _h->GetYaxis()->SetNdivisions( 504 );
      _h->GetYaxis()->SetTitle( "Ratio" );      
      _h->DrawClone( drawoptions );
      
      drawoptions = "ephistsame";
      
      leg->AddEntry( _h, legend_cuts[cuts_], "lp");

    } // Loop over type_.      

    pad_abs_->cd();
    leg->Draw();
    
//    pad_ratio_->SetLogy();

    //== Save
    can_->SetLogy();
    can_->SaveAs(  TString::Format( "Plots_compareEventSelection/can_comparing_eventSelection_" + plots +  "_"  + save_[type_] + ".pdf" ) );

  } // Loop over cuts.
  
  }
  
    
 return 0;
}
