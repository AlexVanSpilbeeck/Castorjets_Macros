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


int main(){

  gStyle->SetOptTitle( 0 );
   setTDRStyle();

  TString plot_;
  std::vector<TString> vector_plot_;
  std::map<TString, TString> save_plot_;
  std::map<TString, bool> dndx_;
   
  plot_ = "hCastorJet_energy";
  vector_plot_.push_back( plot_ );
  save_plot_[ plot_ ] = "Det_energy";
  dndx_[ plot_ ] = true;

  plot_ = "depth";
  vector_plot_.push_back( plot_ );
  save_plot_[ plot_ ] = "depth";
  dndx_[ plot_ ] = false;

  plot_ = "width";
  save_plot_[ plot_ ] = "width";
  vector_plot_.push_back( plot_ );
  dndx_[ plot_ ] = false;

  plot_ = "ehad";
  vector_plot_.push_back( plot_ );
  save_plot_[ plot_ ] = "ehad";
  dndx_[ plot_ ] = false;

  plot_ = "eem";
  vector_plot_.push_back( plot_ );
  save_plot_[ plot_ ] = "eem";
  dndx_[ plot_ ] = false;

  plot_ = "fem";
  vector_plot_.push_back( plot_ );
  save_plot_[ plot_ ] = "fem";
  dndx_[ plot_ ] = false;

  plot_ = "fhot";
  save_plot_[ plot_ ] = "fhot";
  vector_plot_.push_back( plot_ );
  dndx_[ plot_ ] = false;

  plot_ = "sigmaz";
  vector_plot_.push_back( plot_ );
  save_plot_[ plot_ ] = "sigmaz";
  dndx_[ plot_ ] = false;

  plot_ = "phi";
  vector_plot_.push_back( plot_ );
  save_plot_[ plot_ ] = "phi";
  dndx_[ plot_ ] = false;

  plot_ = "nTowers";
  vector_plot_.push_back( plot_ );
  save_plot_[ plot_ ] = "nTowers";
  dndx_[ plot_ ] = false;


  //== Loop over plots.
  for(int p_ = 0; p_ < vector_plot_.size(); p_++){

    plot_ = vector_plot_ [ p_ ];
    bool dndx = dndx_[ plot_ ];

    TString type_;  
    std::vector<TString> vector_type_;
    std::map<TString, TString> legends_type_;
    std::map<TString, TString> match_type_;
    std::map<TString, TString> suffix_type_;
    std::map<TString, TString> save_;

    type_ = "ak5_data_";
    vector_type_.push_back( type_ );
    legends_type_[ type_ ] = "Data";
    suffix_type_[ type_ ] = "";
    match_type_[ type_ ] = "";
    save_[ type_ ] = "Data";

    type_ = "ak5ak5_Pythia6Z2star_NewGeo_";
    vector_type_.push_back( type_ );
    legends_type_[ type_ ] = "Pythia6 (Z2*), New Geo.";
    match_type_[ type_ ] = "deltaPhiMax_0.500000_etaband_0.000000_";
    suffix_type_[ type_ ] = "_matchE";
    save_[ type_ ] = "Pythia6Z2star_NewGeo";

    type_ = "ak5ak5_Pythia6Z2star_down_";
    vector_type_.push_back( type_ );
    legends_type_[ type_ ] = "Pythia6 (Z2*), down";
    match_type_[ type_ ] = "deltaPhiMax_0.500000_etaband_0.000000_";
    suffix_type_[ type_ ] = "_matchE";
    save_[ type_ ] = "Pythia6Z2star_down";

    type_ = "ak5ak5_Pythia6Z2star_up_";
    vector_type_.push_back( type_ );
    legends_type_[ type_ ] = "Pythia6 (Z2*), up";
    match_type_[ type_ ] = "deltaPhiMax_0.500000_etaband_0.000000_";
    suffix_type_[ type_ ] = "_matchE";
    save_[ type_ ] = "Pythia6Z2star_up";

    type_ = "ak5ak5_Pythia84C_NewGeo_";
    vector_type_.push_back( type_ );
    legends_type_[ type_ ] = "Pythia8 (4C), New Geo.";
    match_type_[ type_ ] = "deltaPhiMax_0.500000_etaband_0.000000_";
    suffix_type_[ type_ ] = "_matchE";
    save_[ type_ ] = "Pythia84C_NewGeo";

    type_ = "ak5ak5_EPOS_NewGeo_";
    vector_type_.push_back( type_ );
    legends_type_[ type_ ] = "EPOS, New Geo.";
    match_type_[ type_ ] = "deltaPhiMax_0.500000_etaband_0.000000_";
    suffix_type_[ type_ ] = "_matchE";
    save_[ type_ ] = "EPOS_NewGeo";


    //== Loop over different sources (MC/Data)
    for(int t_ = 0; t_ < vector_type_.size(); t_++){

      //== Prepare names of files.
      TString file_;
      std::vector<TString> vector_;
      std::map<TString, TString> legends_;

      type_ = vector_type_[ t_ ];
      TString suffix_ = suffix_type_[ type_ ];
      TString match_ = match_type_[ type_ ];

      file_ = "/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/" + type_ + "unfold_Emin_150.000000_" + match_ + "all" + suffix_ + ".root";
      vector_.push_back( file_ );
      legends_[ file_ ] = "Original";

      file_ = "/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/" + type_ + "noHFcut_unfold_Emin_150.000000_" + match_ + "all" + suffix_ + ".root";
      vector_.push_back( file_ );
      legends_[ file_ ] = "no HF rq." ;

      file_ = "/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/" + type_ + "noVertexcut_unfold_Emin_150.000000_" + match_ + "all" + suffix_ + ".root";
      vector_.push_back( file_ );
      legends_[ file_ ] =  "no Vtx. rq." ;

      file_ = "/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/" + type_ + "noHFcut_noVertexcut_unfold_Emin_150.000000_" + match_ + "all" + suffix_ + ".root";
      vector_.push_back( file_ );
      legends_[ file_ ] =  "no Vtx./HF rq." ;

      file_ = "/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/" + type_ + "Inclusive_unfold_Emin_150.000000_" + match_ + "all" + suffix_ + ".root";
      vector_.push_back( file_ );
      legends_[ file_ ] =  "most incl." ;

      TCanvas* can_;
      TPad *pad_abs_, *pad_ratio_;
      PrepareCanvas( can_, "Comparing_data" );
      SplitCanvas( can_, pad_abs_, pad_ratio_ );

      TLegend* leg_ = new TLegend( 0.5, 0.65, 1. - pad_abs_->GetRightMargin(), 1. - pad_abs_->GetTopMargin() );
      leg_->SetFillStyle( 0 );
      leg_->SetBorderSize( 0 );
      leg_->SetHeader( legends_type_[ type_ ] );

      TString drawoptions = "hist";
      TH1D* hFirst;
      double ymax = 0.;

      //== Loop over files.
      for(int f_ = 0; f_ < vector_.size(); f_++){
 
        TFile *_file = TFile::Open( vector_[f_], "Read" );
        if( !_file ){ continue; }

        TH1D* hData = (TH1D*)_file->Get( plot_ ); 
      
	if( !hData ){ continue; }
        Prepare_1Dplot( hData );

        if( !dndx ){ 
         hData->GetYaxis()->SetTitle("N"); 
        }
        else{ 
          SetDnDx( hData );
          hData->GetYaxis()->SetTitle("dN/dE"); 
          hData->GetYaxis()->SetRangeUser( 1., 1e5);
        }

        if( f_ == 0 ){
          hFirst = (TH1D*)hData->Clone();
        }

        hData->SetLineColor( getColor( f_ + 1) );
        hData->SetLineStyle( f_ + 1 );
        hData->SetLineWidth( 2 );

        pad_abs_->cd();
        hData->DrawClone( drawoptions );

        pad_ratio_->cd();
        hData->Divide( hFirst );
        hData->GetYaxis()->SetRangeUser( 0. , 1.8 );
        hData->GetYaxis()->SetNdivisions( 504 );
//        hData->GetXaxis()->SetTitle( "E_{cal.} [GeV]" );
        hData->GetXaxis()->SetTitleOffset( 3* hData->GetXaxis()->GetTitleOffset() );
        hData->GetYaxis()->SetTitle( "ratio" );
        hData->Draw( drawoptions );

        drawoptions = "histsame";
  
        leg_->AddEntry( hData, legends_[ vector_[f_] ] , "l");
      }

      pad_abs_->cd();
      leg_->Draw();
  
      if( plot_!= "phi" ) { pad_abs_->SetLogy(); }

      if( !dndx ){can_->SaveAs( TString::Format( "can_" + save_[ type_ ] + "_" + save_plot_[ plot_ ] + "_Comparing_cuts.pdf") ); }
      else{  	can_->SaveAs( TString::Format(   "can_" + save_[ type_ ] + "_" + save_plot_[ plot_ ] + "_Comparing_cuts_dndx.pdf") ); }
    }

  }

 return 0;
}
