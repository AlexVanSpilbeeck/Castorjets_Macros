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
  
  int  new_dir = mkdir( ( TString::Format("Plots_vertices") ).Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  );  

  //== Titles of the axes.

  std::map<int, TString> axis_title;
  axis_title[0] = "good";
  axis_title[1] = "fake";
  axis_title[2] = "NDoF";
  axis_title[3] = "z";
  axis_title[4] = "rho";   
  axis_title[6] = "Nvtx";   
  axis_title[5] = "gvtx";   
  axis_title[7] = "HFact";   

  //== Possible cuts on original file.

  TString cuts_;
  std::vector<TString> vector_cuts_;

  cuts_ = "";
  vector_cuts_.push_back( cuts_ );

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

  //== Set of selection criteria.
  TString selectedEvents_;
  std::vector<TString> vector_selectedEvents_;

  std::map< TString, std::map<TString, double> > axis_min;
  std::map< TString, std::map<TString, double> > axis_max;
  std::map< TString, TString> selectionLegend_;


  //== The most inclusive one, plots all vertices.
  selectedEvents_ = "";
  vector_selectedEvents_.push_back( selectedEvents_ );
  selectionLegend_[  selectedEvents_ ] = "Inclusive";

  axis_min[ selectedEvents_ ]["good"] = -1.;
  axis_min[ selectedEvents_ ]["fake"] = -1.;
  axis_min[ selectedEvents_ ]["NDoF"] = -1.;
  axis_min[ selectedEvents_ ]["rho"] = 	-1.;
  axis_min[ selectedEvents_ ]["z"] = 	-1.;
  axis_min[ selectedEvents_ ]["Nvtx"] = -1.;
  axis_min[ selectedEvents_ ]["gvtx"] = -1.;
  axis_min[ selectedEvents_ ]["HFact"]= 0.;

  axis_max[ selectedEvents_ ]["good"] = 1.;
  axis_max[ selectedEvents_ ]["fake"] = 1.;
  axis_max[ selectedEvents_ ]["NDoF"] = 20.;
  axis_max[ selectedEvents_ ]["rho"] = 	2.;
  axis_max[ selectedEvents_ ]["z"] = 	15.;
  axis_max[ selectedEvents_ ]["Nvtx"] = 20.;
  axis_max[ selectedEvents_ ]["gvtx"] = 20.;
  axis_max[ selectedEvents_ ]["HFact"]= 1.;

  //== This plots the properties of all vertices that have been labelled "good".
  selectedEvents_ = "_onlyGoodVertices";
  vector_selectedEvents_.push_back( selectedEvents_ );
  selectionLegend_[  selectedEvents_ ] = "Good Vtx";

  axis_min[ selectedEvents_ ]["good"] = 1.;
  axis_min[ selectedEvents_ ]["fake"] = 0.;
  axis_min[ selectedEvents_ ]["NDoF"] = 4.;
  axis_min[ selectedEvents_ ]["rho"] = 0.;
  axis_min[ selectedEvents_ ]["z"] = 0.;
  axis_min[ selectedEvents_ ]["Nvtx"] = 0.;
  axis_min[ selectedEvents_ ]["gvtx"] = 0.;
  axis_min[ selectedEvents_ ]["HFact"] = 0.;

  axis_max[ selectedEvents_ ]["good"] = 1.;
  axis_max[ selectedEvents_ ]["fake"] = 0.;
  axis_max[ selectedEvents_ ]["NDoF"] = 21.;
  axis_max[ selectedEvents_ ]["rho"] = 2.;
  axis_max[ selectedEvents_ ]["z"] = 15.;
  axis_max[ selectedEvents_ ]["Nvtx"] = 20.;
  axis_max[ selectedEvents_ ]["gvtx"] = 20.;
  axis_max[ selectedEvents_ ]["HFact"] = 1.;


  //== This plots the properties of vertices in events with exactly 1 good vertex.
  selectedEvents_ = "_oneGoodVertex";
  vector_selectedEvents_.push_back( selectedEvents_ );
  selectionLegend_[  selectedEvents_ ] = "1 good Vtx/evt";

  axis_min[ selectedEvents_ ]["good"] = 1.;
  axis_min[ selectedEvents_ ]["fake"] = 0.;
  axis_min[ selectedEvents_ ]["NDoF"] = 4.;
  axis_min[ selectedEvents_ ]["rho"] = 	0.;
  axis_min[ selectedEvents_ ]["z"] = 	0.;
  axis_min[ selectedEvents_ ]["Nvtx"] = 1.;
  axis_min[ selectedEvents_ ]["gvtx"] = 1.;
  axis_min[ selectedEvents_ ]["HFact"]= 0.;

  axis_max[ selectedEvents_ ]["good"] = 1.;
  axis_max[ selectedEvents_ ]["fake"] = 0.;
  axis_max[ selectedEvents_ ]["NDoF"] = 21.;
  axis_max[ selectedEvents_ ]["rho"] = 	2.;
  axis_max[ selectedEvents_ ]["z"] = 	15.;
  axis_max[ selectedEvents_ ]["Nvtx"] = 1.;
  axis_max[ selectedEvents_ ]["gvtx"] = 1.;
  axis_max[ selectedEvents_ ]["HFact"]= 1.;

  //== This plots the properties of all vertices that have been labelled "good".
  selectedEvents_ = "_oneGoodVertices_HFactivity";
  vector_selectedEvents_.push_back( selectedEvents_ );
  selectionLegend_[  selectedEvents_ ] = "1 good Vtx/evt + HF";

  axis_min[ selectedEvents_ ]["good"] =	0.;
  axis_min[ selectedEvents_ ]["fake"] = axis_min[ "_oneGoodVertex" ]["fake"];
  axis_min[ selectedEvents_ ]["NDoF"] = axis_min[ "_oneGoodVertex" ]["NDoF"];
  axis_min[ selectedEvents_ ]["rho"] = 	axis_min[ "_oneGoodVertex" ]["rho"];
  axis_min[ selectedEvents_ ]["z"] = 	axis_min[ "_oneGoodVertex" ]["z"];
  axis_min[ selectedEvents_ ]["Nvtx"] = axis_min[ "_oneGoodVertex" ]["Nvtx"] ;
  axis_min[ selectedEvents_ ]["gvtx"] = 1.;
  axis_min[ selectedEvents_ ]["HFact"]= 1.;

  axis_max[ selectedEvents_ ]["good"] = axis_max[ "_oneGoodVertex" ]["good"];
  axis_max[ selectedEvents_ ]["fake"] = axis_max[ "_oneGoodVertex" ]["fake"];
  axis_max[ selectedEvents_ ]["NDoF"] = axis_max[ "_oneGoodVertex" ]["NDoF"];
  axis_max[ selectedEvents_ ]["rho"] = 	axis_max[ "_oneGoodVertex" ]["rho"];
  axis_max[ selectedEvents_ ]["z"] = 	axis_max[ "_oneGoodVertex" ]["z"];
  axis_max[ selectedEvents_ ]["Nvtx"] = axis_max[ "_oneGoodVertex" ]["Nvtx"] ;
  axis_max[ selectedEvents_ ]["gvtx"] = axis_max[ "_oneGoodVertex" ]["gvtx"];
  axis_max[ selectedEvents_ ]["HFact"]= 1.;

  //== This plots the properties of all vertices that have been labelled "good".
  selectedEvents_ = "_HFactivity";
  vector_selectedEvents_.push_back( selectedEvents_ );
  selectionLegend_[  selectedEvents_ ] = "HF";

  axis_min[ selectedEvents_ ]["good"] = 0.;
  axis_min[ selectedEvents_ ]["fake"] = 0.;
  axis_min[ selectedEvents_ ]["NDoF"] = 0.;
  axis_min[ selectedEvents_ ]["rho"] = 	0.;
  axis_min[ selectedEvents_ ]["z"] = 	0.;
  axis_min[ selectedEvents_ ]["Nvtx"] = 0.;
  axis_min[ selectedEvents_ ]["gvtx"] = 0.;
  axis_min[ selectedEvents_ ]["HFact"]= 1.;

  axis_max[ selectedEvents_ ]["good"] = 1.;
  axis_max[ selectedEvents_ ]["fake"] = 1.;
  axis_max[ selectedEvents_ ]["NDoF"] = 20.;
  axis_max[ selectedEvents_ ]["rho"] = 	2.;
  axis_max[ selectedEvents_ ]["z"] = 	15.;
  axis_max[ selectedEvents_ ]["Nvtx"] = 20.;
  axis_max[ selectedEvents_ ]["gvtx"] = 20.;
  axis_max[ selectedEvents_ ]["HFact"]= 1.;
  
  //== This plots the properties of all vertices that have been labelled "good".
  selectedEvents_ = "_0_goodVtx";
  vector_selectedEvents_.push_back( selectedEvents_ );
  selectionLegend_[  selectedEvents_ ] = "0 good vertex";

  axis_min[ selectedEvents_ ]["good"] = 0.;
  axis_min[ selectedEvents_ ]["fake"] = 0.;
  axis_min[ selectedEvents_ ]["NDoF"] = 0.;
  axis_min[ selectedEvents_ ]["rho"] = 	0.;
  axis_min[ selectedEvents_ ]["z"] = 	0.;
  axis_min[ selectedEvents_ ]["Nvtx"] = 0.;
  axis_min[ selectedEvents_ ]["gvtx"] = 0.;
  axis_min[ selectedEvents_ ]["HFact"]= 0.;

  axis_max[ selectedEvents_ ]["good"] = 0.;
  axis_max[ selectedEvents_ ]["fake"] = 20.;
  axis_max[ selectedEvents_ ]["NDoF"] = 20.;
  axis_max[ selectedEvents_ ]["rho"] = 	2.;
  axis_max[ selectedEvents_ ]["z"] = 	15.;
  axis_max[ selectedEvents_ ]["Nvtx"] = 20.;
  axis_max[ selectedEvents_ ]["gvtx"] = 20.;
  axis_max[ selectedEvents_ ]["HFact"]= 1.;  
  
  
  //== This plots the properties of all vertices that have been labelled "good".
  selectedEvents_ = "_exact1_badVtx";
  vector_selectedEvents_.push_back( selectedEvents_ );
  selectionLegend_[  selectedEvents_ ] = "ext. 1 bad vertex";

  axis_min[ selectedEvents_ ]["good"] = 0.;
  axis_min[ selectedEvents_ ]["fake"] = 1.;
  axis_min[ selectedEvents_ ]["NDoF"] = 0.;
  axis_min[ selectedEvents_ ]["rho"] = 	0.;
  axis_min[ selectedEvents_ ]["z"] = 	0.;
  axis_min[ selectedEvents_ ]["Nvtx"] = 0.;
  axis_min[ selectedEvents_ ]["gvtx"] = 0.;
  axis_min[ selectedEvents_ ]["HFact"]= 0.;

  axis_max[ selectedEvents_ ]["good"] = 0.;
  axis_max[ selectedEvents_ ]["fake"] = 1.;
  axis_max[ selectedEvents_ ]["NDoF"] = 20.;
  axis_max[ selectedEvents_ ]["rho"] = 	2.;
  axis_max[ selectedEvents_ ]["z"] = 	15.;
  axis_max[ selectedEvents_ ]["Nvtx"] = 20.;
  axis_max[ selectedEvents_ ]["gvtx"] = 20.;
  axis_max[ selectedEvents_ ]["HFact"]= 1.;    



  //== Plots.

  TString plots;
  std::vector<TString> histograms;
  
  plots = "hVertex_good" ;
  histograms.push_back( plots );

  plots = "hVertex_fake" ;
  histograms.push_back( plots );
  
  plots = "hVertex_ndof" ;
  histograms.push_back( plots );

  plots = "hVertex_z" ;
  histograms.push_back( plots );

  plots = "hVertex_rho" ;
  histograms.push_back( plots );

  //== Loop over cuts in sample files.
  for( int c_ = 0; c_ < vector_cuts_.size(); c_++){
    cuts_ = vector_cuts_[ c_ ];

    //== Loop over selection criteria.
    for( int s_ = 0; s_ < vector_selectedEvents_.size(); s_++){
      selectedEvents_ = vector_selectedEvents_[ s_ ];

      //== Loop over files.
      for(int t_ = 0; t_ < vector_type_.size() ; t_++) {
        type_ = vector_type_[ t_ ];

        TFile* _file = TFile::Open( "20160511_counting_vertices" + type_ + cuts_ + ".root", "Read" );
cout << "\tFILE\t" << TString::Format( "20160511_counting_vertices" + type_ + cuts_ + ".root" ) << endl;
        if( !_file ){ continue; }
        THnSparseD *sparse = (THnSparseD*)_file->Get("hVertex_info");

	for(int axis = 0; axis < axis_title.size(); axis++){
	  sparse->GetAxis( axis )->SetRangeUser( axis_min[ selectedEvents_ ][ axis_title[axis] ], axis_max[ selectedEvents_  ][ axis_title[axis] ] );
	}

        for(int axis1 = 0; axis1 < axis_title.size(); axis1++){
          TString ax1 = axis_title[axis1];

          for(int axis2 = axis1+1; axis2 < axis_title.size(); axis2++){

            TString ax2 = axis_title[axis2];

            //== Extraction.
            TH2D* _h2 = (TH2D*)sparse->Projection( axis2, axis1);

            _h2->GetXaxis()->SetTitle( ax1 );
            _h2->GetYaxis()->SetTitle( ax2 );

            //== Canvas.
            TCanvas *can_;
            TString can_title = ax1 + "_vs_" + ax2 ;
            PrepareCanvas_2D( can_, can_title );

            //== Draw.
            Prepare_2Dplot( _h2 ) ;
            _h2->Draw("colz");

            //== Save.
            can_->SetLogz();
            can_->SaveAs(  TString::Format( "can_vertex_" + ax1 + "_vs_" + ax2 + "_" + type_ + cuts_ + selectedEvents_ + ".pdf" ) );
          } //== Axis 2
        } //== Axis1
      } // Loop over type_.


      for(int axis1 = 0; axis1 < axis_title.size(); axis1++){
        TString ax1 = axis_title[axis1];      

        TCanvas *can_;
        TString can_title = ax1 + "_1D" ;
        PrepareCanvas( can_, can_title );
	TString drawoptions = "phist";
	TLegend *leg_ = new TLegend( 0.5, 0.65, 1. - can_->GetRightMargin(), 1.-can_->GetTopMargin() );
	leg_->SetBorderSize( 0 );
	leg_->SetFillStyle( 0 );
	leg_->SetHeader( selectionLegend_[  selectedEvents_ ] );

        //== Loop over files.
        for(int t_ = 0; t_ < vector_type_.size() ; t_++) {

          type_ = vector_type_[ t_ ];
          TFile* _file = TFile::Open( "20160511_counting_vertices" + type_ + cuts_ + ".root", "Read" );
          if( !_file ){ continue; }
          THnSparseD *sparse = (THnSparseD*)_file->Get("hVertex_info");
          THnSparseD *sparse_uncut = (THnSparseD*)_file->Get("hVertex_info");

	  for(int axis = 0; axis < axis_title.size(); axis++){
 	    sparse->GetAxis( axis )->SetRangeUser( axis_min[ selectedEvents_ ][ axis_title[axis] ], axis_max[ selectedEvents_  ][ axis_title[axis] ] );
cout << "\t" <<  selectedEvents_ << "\t" << axis << "\t" << axis_title[axis] << "\t" << axis_min[ selectedEvents_ ][ axis_title[axis] ] << "\tto\t" << axis_max[ selectedEvents_  ][ axis_title[axis] ] << endl;
	  }

	  TH1D* _h1 = (TH1D*)sparse->Projection( axis1 );
	  TH1D* _h1_uncut = (TH1D*)sparse_uncut->Projection( axis1 );

	  _h1->SetLineColor( getColor( t_ + 1) );
	  _h1->SetLineStyle( t_ + 1 );
	  _h1->SetLineWidth( 2 );
	  _h1->Scale( 1./_h1_uncut->Integral() );
	  Prepare_1Dplot( _h1 );
	  if( (axis_title[ axis1 ] == "fake" || axis_title[ axis1 ] == "good" || (axis_title[ axis1 ]).Contains("HF") ) ){ _h1->GetYaxis()->SetRangeUser(0.,1.); }
	  _h1->Draw( drawoptions );
	  drawoptions = "histsame";
	  leg_->AddEntry( _h1, legends_type_[ type_ ], drawType_[ type_ ]);

        } //== Loop over file.

	leg_->Draw();
	if( ! (axis_title[ axis1 ] == "fake" || axis_title[ axis1 ] == "good" || (axis_title[ axis1 ]).Contains("HF") )  ) { can_->SetLogy(); }
        can_->SaveAs(  TString::Format( "can_vertex_comparingSamples_" + ax1 + "_" + cuts_ + selectedEvents_ + ".pdf" ) );
      } //== Loop over axes.
    } // Loop over all or just selected events.
    //=========================================================

    /**********************************************************************
    * Compare evolution of distribution of 1 quantity with changing cuts. *
    **********************************************************************/ 

    //== Loop over axes.
    for(int axis1 = 0; axis1 < axis_title.size(); axis1++){
      TString ax1 = axis_title[axis1];  

      //== Loop over samples.
      for(int t_ = 0; t_ < vector_type_.size() ; t_++) {
        type_ = vector_type_[ t_ ];
        TFile* _file = TFile::Open( "20160511_counting_vertices" + type_ + cuts_ + ".root", "Read" );
        if( !_file ){ continue; }
        THnSparseD *sparse = (THnSparseD*)_file->Get("hVertex_info");
        THnSparseD *sparse_uncut = (THnSparseD*)_file->Get("hVertex_info");

        //== The canvas.
        TCanvas *can_;
        TString canvastitle = plots + cuts_;
        PrepareCanvas( can_, canvastitle );	
	TString drawoptions = "hist";
	if( type_.Contains("ata") ){ drawoptions = "phist"; }
	TLegend *leg_ = new TLegend( 0.5, 0.65, 1. - can_->GetRightMargin(), 1.-can_->GetTopMargin() );
	leg_->SetBorderSize( 0 );
	leg_->SetFillStyle( 0 );
	leg_->SetHeader( legends_type_[type_] );

        //== Loop over selection criteria.
        for( int s_ = 0; s_ < vector_selectedEvents_.size(); s_++){
          selectedEvents_ = vector_selectedEvents_[ s_ ];

	  //== Extract the histogram and cut events.
          THnSparseD *sparse = (THnSparseD*)_file->Get("hVertex_info");
	  for(int axis = 0; axis < axis_title.size(); axis++){
 	    sparse->GetAxis( axis )->SetRangeUser( axis_min[ selectedEvents_ ][ axis_title[axis] ], axis_max[ selectedEvents_  ][ axis_title[axis] ] );
	  }

	  TH1D* _h1 = (TH1D*)sparse->Projection( axis1 );
	  TH1D* _h1_uncut = (TH1D*)sparse_uncut->Projection( axis1 );

	  _h1->SetLineColor( getColor( s_ + 1) );
	  _h1->SetLineStyle( s_ + 1 );
	  _h1->SetLineWidth( 2 );
	  _h1->SetMarkerColor( getColor( s_ + 1) );
	  _h1->SetMarkerStyle( s_ + 19 );
	  _h1->SetMarkerSize( 1.5 );
	  _h1->Scale( 1./_h1_uncut->Integral() );
	  Prepare_1Dplot( _h1 );
	  if( (axis_title[ axis1 ] == "fake" || axis_title[ axis1 ] == "good" || (axis_title[ axis1 ]).Contains("HF") ) ){ _h1->GetYaxis()->SetRangeUser(0.,1.); }
	  _h1->Draw( drawoptions );
	  drawoptions = "histsame";
	  if( type_.Contains("ata") ){ drawoptions = "psame"; }
	  leg_->AddEntry( _h1 , selectionLegend_[  selectedEvents_ ], drawType_[ type_ ] );
	  

        }//== Loop over selection criteria.

        leg_->Draw();
	if( ! (axis_title[ axis1 ] == "fake" || axis_title[ axis1 ] == "good" || (axis_title[ axis1 ]).Contains("HF") )  ) { can_->SetLogy(); }
        can_->SaveAs(  TString::Format( "can_vertex_comparing_eventSelection_" + ax1 + "_" + type_ + ".pdf" ) );        

      } //== Loop over samples.
    } //== Loop over axes.
  } // Loop over cuts.
  
  
  
  
  /***********************************************/
  /*/						/*/
  /***********************************************/
 
  //== Loop over cuts in sample files.
  for( int c_ = 0; c_ < vector_cuts_.size(); c_++){
    cuts_ = vector_cuts_[ c_ ];

    //== Loop over files.
    for(int t_ = 0; t_ < vector_type_.size() ; t_++) {
      type_ = vector_type_[ t_ ];

      TFile* _file = TFile::Open( "20160517_counting_vertices" + type_ + ".root", "Read" );
      if( !_file ){ continue; }
      
      TCanvas* canvas_tower_vs_vtx;
      PrepareCanvas_2D( canvas_tower_vs_vtx, "can_nTowers_vsnVtx" );
      
      TH2I* hCastorTowers_vs_nVtx = (TH2I*)_file->Get("hCastorTowers_vs_Nvtx");
      if( !hCastorTowers_vs_nVtx ) continue;
      
//      hCastorTowers_vs_nVtx->Scale( 1./hCastorTowers_vs_nVtx->Integral() );
      hCastorTowers_vs_nVtx->GetXaxis()->SetRangeUser(0, 16 );
      hCastorTowers_vs_nVtx->Draw("surf3");
      
      canvas_tower_vs_vtx->SetLogz();
      canvas_tower_vs_vtx->SaveAs("can_nTowers_vsnVtx" + type_ + ".pdf");
      
      
    }
  }
  
 return 0;
}
