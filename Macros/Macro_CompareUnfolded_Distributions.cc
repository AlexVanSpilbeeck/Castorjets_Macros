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
#include "../Functions/Function_Empty_UOflow.h"
#include "../Functions/Function_PrepareCanvas.h"
#include "../Functions/Function_Prepare1Dplot.h"
#include "../Functions/Function_Prepare2Dplot.h"

#include "../../../../RooUnfold-1.1.1/src/RooUnfold.h"
#include "../../../../RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"
#include "../../../../RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "../../../../RooUnfold-1.1.1/src/RooUnfoldResponse.h"


using namespace std;

void PlotFakes(   vector<TString> , vector<TString> , map<TString, TString> , vector<bool> , map<bool, TString>  , map<bool, TString> );
void PlotMisses(   vector<TString> , vector<TString> , map<TString, TString> , vector<bool> , map<bool, TString>  , map<bool, TString> );


int main(){

  gStyle->SetOptTitle( 0 );
  setTDRStyle();
  gStyle->SetPalette(1);

  //== Define suffixes.
  TString cal_;
  vector<TString> vector_cal_;
  map<TString, TString> legends_;

  cal_ = "";
  vector_cal_.push_back( cal_ );
  legends_[ cal_ ] = "Had. cal. (100%)";
/*
  cal_ = "_allJetsCalibrated";
  vector_cal_.push_back( cal_ );
  legends_[ cal_ ] = "All. cal.";

  cal_ = "_90_percentHadjetsCalibrated";
  vector_cal_.push_back( cal_ );
  legends_[ cal_ ] = "Had. cal. (90%)";

  cal_ = "_80_percentHadjetsCalibrated";
  vector_cal_.push_back( cal_ );
  legends_[ cal_ ] = "Had. cal. (80%)";

  cal_ = "_70_percentHadjetsCalibrated";
  vector_cal_.push_back( cal_ );
  legends_[ cal_ ] = "Had. cal. (70%)";
*/  
  //== Event selection.
  TString sel_;
  vector<TString> vector_sel_;
  map<TString, TString> legends_sel_;
  
  sel_ = "";
  vector_sel_.push_back( sel_ );
  legends_sel_[ sel_ ] = "Default";
  
  sel_ = "_one_none_vertex_BSCor_HFor";
  vector_sel_.push_back( sel_ );
  legends_sel_[ sel_ ] = "New";  
  
  sel_ = "_MostInclusive";
  vector_sel_.push_back( sel_ );
  legends_sel_[ sel_ ] = "Most Inclusive";  
    
  //== Use same MC sample or different MC sample.
  bool same_sample_;
  vector<bool> vector_sample_;
  map<bool, TString> legends_sample_;
  map<bool, TString> save_sample_;
  
  same_sample_ = true;
  vector_sample_.push_back( same_sample_ );
  legends_sample_[ same_sample_ ] = "Uniform MC";
  save_sample_[ same_sample_ ] = "_unaltered_mc";
/*
  same_sample_ = false;
  vector_sample_.push_back( same_sample_ );
  legends_sample_[ same_sample_ ] = "Non-uniform MC";
  save_sample_[ same_sample_ ] = "_altered_mc";
*/  
  //== First hist, for ratio.
  TH1D* hFirst;
  TH1D* hFirst_cal;

  ofstream ratios;
  ratios.open("Ratios.txt");

  for(int mc_ = 0 ; mc_ < vector_sample_.size(); mc_++){  
    same_sample_ = vector_sample_[mc_];
    TString save_mc_ = save_sample_[ same_sample_ ];
    TString legend_mc_ = legends_sample_[ same_sample_ ];
    
    for(int s_ = 0; s_ < vector_sel_.size(); s_++){  
      sel_ = vector_sel_[s_];
      
      ratios << "***********\n";
      ratios << sel_ << endl;
      ratios << "***********\n";
      
      //== Unfolded level canvas. 
      TString canvastitle = "can_comparison_calibrations";
      TCanvas *can_;  
      PrepareCanvas( can_, canvastitle );
      TPad* pad_abs_, *pad_ratio_;
      SplitCanvas(can_, pad_abs_, pad_ratio_);

      //== Det level canvas.
      TCanvas *can_cal_;
      PrepareCanvas( can_cal_, "can_comparison_calibrations_detlevel" );
      TPad* cal_abs_, *cal_ratio_;
      SplitCanvas(can_cal_, cal_abs_, cal_ratio_);

      //== Response matrix canvas.
      TCanvas *can_resp_;
      PrepareCanvas_2D( can_resp_, "can_response_matrix");

      //== Fakes canvas.
      TCanvas *can_fakes_;  
      PrepareCanvas( can_fakes_, "can_comparison_calibrations_fakes" );
      TPad* pad_abs_fakes_, *pad_ratio_fakes_;
      SplitCanvas_equally( can_fakes_, pad_abs_fakes_, pad_ratio_fakes_ );      
      
      //== Legend.
      double ticksize = gStyle->GetTickLength();
    
      TLegend *leg = new TLegend( 0.7, 0.75, 1. - pad_abs_->GetRightMargin() - ticksize, 1. - pad_abs_->GetTopMargin() - ticksize);
      leg->SetFillStyle( 0 );
      leg->SetBorderSize( 0 );    
    
      TString drawoptions = "hist";  
      int file_ = 0;
 
      //== Loop.
      for(int c_ = 0; c_ < vector_cal_.size(); c_++){
        cal_ = vector_cal_[ c_ ];

        ratios << "\n\n" << cal_ << endl;

        //== Open data, extract histogram.
        TFile* _datafile = TFile::Open( "/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5_data" + sel_ + "_unfold_Emin_150.000000_all" + cal_ + ".root", "read" );
        if( !_datafile ){ continue; }
        TH1D* _datahist = (TH1D*)_datafile->Get("hCastorJet_energy");

        double lumi_ =  1.23115 * 1e5;
        lumi_ = lumi_*1/0.982;
        _datahist->Scale( 1./lumi_);
        TH1D* _datadet = (TH1D*)_datahist->Clone("DetectorLevel");    

        //== Open MC, extract RooUnfold.  
        TFile *_file0;
        if( !same_sample_) _file0 = new TFile( "/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo" + sel_ + "_unfold_Emin_150.000000_deltaPhiMax_0.500000_etaband_0.000000_all_matchE" + cal_ + ".root", "read");
        if( same_sample_) _file0 = new TFile( "/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo" + sel_ + "_unfold_Emin_150.000000_deltaPhiMax_0.500000_etaband_0.000000_all_matchE.root", "read");
        if( !_file0 ) { continue; }
        RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response");

        // Response matrix.
        TH2D* hResponse = (TH2D*)response->Hresponse();
        can_resp_->cd();
        Prepare_2Dplot( hResponse );
        hResponse->Draw("colz");
        can_resp_->SetLogz();
        can_resp_->SaveAs( TString::Format("can_response_" + cal_ + "_" + sel_ + ".pdf") );

        // Use hist as input for the unfolding.
        RooUnfoldBayes unfold_bayes(response, _datahist, 30);
        _datahist = (TH1D*) unfold_bayes.Hreco( ); 

        SetDnDx( _datahist );
        SetDnDx( _datadet );

        if( file_ == 0 ){
          hFirst = (TH1D*)_datahist->Clone( "first");
          hFirst_cal = (TH1D*)_datadet->Clone( "first");
          file_++;
        }

        //== Unfolded.
        _datahist->SetLineColor( getColor( c_ + 1 ) );
        _datahist->SetLineStyle( c_ + 1  );
        _datahist->SetLineWidth( 2 );
        _datahist->GetXaxis()->SetTitle("E [GeV]");
        _datahist->GetYaxis()->SetTitle("d#sigma/dE [mb/GeV]");
        _datahist->GetYaxis()->SetRangeUser(2*1e-5, 0.5 );
        Prepare_1Dplot( _datahist );
        leg->AddEntry( _datahist, legends_[ cal_ ] , "l");
  
        pad_abs_->cd();
        _datahist->DrawClone( drawoptions );
         
        pad_ratio_->cd();
        _datahist->Divide( hFirst );
        _datahist->GetYaxis()->SetTitle("Ratio");
        _datahist->GetYaxis()->SetNdivisions(504);
        _datahist->GetYaxis()->SetRangeUser( 0.5, 1.25);
        _datahist->GetXaxis()->SetTitleOffset( 3. * _datahist->GetXaxis()->GetTitleOffset() );
        _datahist->Draw( drawoptions );

        //== Detector level plot.
        cal_abs_->cd();
        _datadet->SetLineColor( getColor( c_ + 1 ) );
        _datadet->SetLineStyle( c_ + 1  );
        _datadet->SetLineWidth( 2 );
        _datadet->GetXaxis()->SetTitle("E_{det} [GeV]");
        _datadet->GetYaxis()->SetTitle("d#sigma/dE [mb/GeV]");
        _datadet->GetYaxis()->SetRangeUser(2*1e-5, 0.5 );
        Prepare_1Dplot( _datadet );
        _datadet->DrawClone( drawoptions );

        cal_ratio_->cd();
        _datadet->Divide( hFirst_cal );
        _datadet->GetYaxis()->SetTitle("Ratio");
        _datadet->GetYaxis()->SetNdivisions(504);
        _datadet->GetYaxis()->SetRangeUser( 0.5, 1.25);
        _datadet->GetXaxis()->SetTitleOffset( 3. * _datadet->GetXaxis()->GetTitleOffset() );
        _datadet->Draw( drawoptions );
 
        drawoptions = "histsame";
      }    
 
      pad_abs_->cd();
      pad_abs_->SetLogy();
      leg->SetHeader("Unfolded " + legend_mc_);
      leg->Draw();
      can_->SaveAs( "can_comparison_calibrations" + save_mc_ + "_" + sel_ + ".pdf" );

      cal_abs_->cd();
      cal_abs_->SetLogy();
      leg->SetHeader( "Detector level");
      leg->Draw();
      can_cal_->SaveAs( "can_comparison_calibrations_detlevel" + save_mc_ + "_" + sel_ + ".pdf" );


    }
  } //== Loop over use of same/different MC sample.

  PlotFakes( vector_cal_, vector_sel_ , legends_sel_, vector_sample_, legends_sample_ , save_sample_);
  PlotMisses( vector_cal_, vector_sel_ , legends_sel_, vector_sample_, legends_sample_ , save_sample_);

  return 0;
}



void PlotMisses(
	vector<TString> vector_cal_, 
	vector<TString> vector_sel_ , map<TString, TString> legends_sel_, 
	vector<bool> vector_sample_,  map<bool, TString> legends_sample_ , map<bool, TString> save_sample_){

  int  new_dir = mkdir( ( TString::Format("Plots_distributions/") ).Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  );

  ofstream ratios;
  ratios.open("Ratios.txt", ios::app);
  
  //== Loop.
  for(int mc_ = 0 ; mc_ < vector_sample_.size(); mc_++){  
    bool same_sample_ = vector_sample_[mc_];
    TString save_mc_ = save_sample_[ same_sample_ ];
    TString legend_mc_ = legends_sample_[ same_sample_ ];

    TString drawoptions = "hist";
        
    //== Misses canvas.
    TCanvas *can_misses_;  
    PrepareCanvas( can_misses_, "can_comparison_calibrations_misses" );
    TPad* pad_abs_, *pad_ratio_;
    SplitCanvas_equally( can_misses_, pad_abs_, pad_ratio_ );
      
    //== Legend.
    double ticksize = gStyle->GetTickLength();
    
    TLegend *leg = new TLegend( 0.65, 0.35, 1. - pad_abs_->GetRightMargin() - ticksize, 1. - pad_abs_->GetTopMargin() - ticksize);
    leg->SetFillStyle( 0 );
    leg->SetBorderSize( 0 );   

    int distribution_nr = 1;

    for(int s_ = 0; s_ < vector_sel_.size(); s_++){  
      TString sel_ = vector_sel_[s_];             
      
      for(int c_ = 0; c_ < vector_cal_.size(); c_++){
        TString cal_ = vector_cal_[ c_ ];


        //== Open MC, extract RooUnfold.  
        TFile *_file0;
        if( !same_sample_) _file0 = new TFile( "/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo" + sel_ + "_unfold_Emin_150.000000_deltaPhiMax_0.500000_etaband_0.000000_all_matchE" + cal_ + ".root", "read");
        if( same_sample_) _file0 = new TFile( "/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo" + sel_ + "_unfold_Emin_150.000000_deltaPhiMax_0.500000_etaband_0.000000_all_matchE.root", "read");
        if( !_file0 ) { continue; }

        RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response");
        TH2D* _hResponse = (TH2D*)response->Hresponse();
      
        //== Misses.
        can_misses_->cd();

        TH1D* _hGen 	= (TH1D*)_file0->Get("hGenJet_energy"); Empty_UOflow( _hGen );
        TH1D* _hMatched = (TH1D*)_hResponse->ProjectionY();	Empty_UOflow( _hMatched );
        TH1D* _hTruth 	= (TH1D*)response->Htruth();  		Empty_UOflow( _hTruth );
        TH1D* _hMiss 	= (TH1D*)_hTruth->Clone( "miss" );	Empty_UOflow( _hMiss );
        _hMiss->Add( _hMatched, -1. );  
  
  	TH1D* _hFreeGenJet = (TH1D*)_file0->Get("hGenJet_energy_noDet"); 
  	if( _hFreeGenJet ){ 
  	  if( _hFreeGenJet->Integral() > 0.){
  	    Empty_UOflow( _hFreeGenJet  ); SetDnDx( _hFreeGenJet ); 
  	  }
  	}
  
  cout << "Filename\t" << _file0->GetName() << endl;
  
  	ratios << "\n\t" << _file0->GetName() << endl;
        ratios << "Misses\t" << 	_hMiss->Integral()/_hGen->Integral() << endl;  
        ratios << "Misses\t" << 	_hMiss->Integral()/_hTruth->Integral() << endl;    
        ratios << "Misses\t" << 1. - 	_hMatched->Integral()/_hGen->Integral() << endl;
        ratios << "Misses\t" << 1. - 	_hMatched->Integral()/_hTruth->Integral() << endl;
        ratios << "Misses\t" << 1. - 	_hResponse->Integral()/_hTruth->Integral() << endl;        
         
        SetDnDx( _hMiss );
        SetDnDx( _hTruth );
        SetDnDx( _hMatched );  
        SetD2nDxDy( _hResponse );

        _hMiss->GetXaxis()->SetTitle("E [GeV]");
        _hMiss->GetXaxis()->SetRangeUser(150., 1600.);
        _hMiss->GetYaxis()->SetTitle("dN/dE [1/GeV]");
        _hMiss->GetYaxis()->SetRangeUser( 5e-1, 5e4 );
        _hMiss->SetLineColor( getColor( distribution_nr  ) );
        _hMiss->SetLineStyle( distribution_nr++  );
        _hMiss->SetLineWidth( 2 );
        Prepare_1Dplot( _hMiss );
        Prepare_1Dplot( _hTruth );                        
        
        pad_abs_->cd();
        _hMiss->DrawClone( drawoptions );
        
        _hTruth->SetLineColor( getColor( distribution_nr ) );
        _hTruth->SetLineStyle( distribution_nr++  );        
        _hTruth->Draw("histsame");        

        pad_ratio_->cd();
        _hMiss->Divide( _hTruth );
        _hMiss->GetXaxis()->SetTitleOffset( _hMiss->GetXaxis()->GetTitleOffset() * 2. );
        _hMiss->GetYaxis()->SetTitle("Misses/Gen. jets");
        _hMiss->GetYaxis()->SetRangeUser(1e-3, 1.2);
        _hMiss->Draw( drawoptions ); 
        
	if( _hFreeGenJet ){
  	  if( _hFreeGenJet->Integral() > 0.){	
	    pad_abs_->cd();
	    _hFreeGenJet->SetMarkerColor( getColor(distribution_nr++ )) ;
	    _hFreeGenJet->SetMarkerStyle( 20 + distribution_nr++ );
	    _hFreeGenJet->DrawClone("psame");
            leg->AddEntry( _hFreeGenJet, "Free miss: " + legends_sel_[ sel_ ] , "p"); 	  
          
            pad_ratio_->cd();
            _hFreeGenJet->Divide( _hTruth );
            _hFreeGenJet->DrawClone( "psame" );
          }
	}        
        
        leg->AddEntry( _hMiss, "Misses: " + legends_sel_[ sel_ ] , "l"); 
        leg->AddEntry( _hTruth, "True: " + legends_sel_[ sel_ ] , "l");         
        
        drawoptions = "histsame";

      }    
    }
    can_misses_->cd();
    pad_abs_->SetLogy();
    pad_abs_->cd();
//    leg->SetHeader("Gen. level");
    leg->Draw();
    pad_ratio_->SetLogy();
    can_misses_->SaveAs( "Plots_distributions/can_comparison_calibrations_misses" + save_mc_ + ".pdf" );        
          
  }
}








void PlotFakes(
	vector<TString> vector_cal_, 
	vector<TString> vector_sel_ , map<TString, TString> legends_sel_, 
	vector<bool> vector_sample_,  map<bool, TString> legends_sample_ , map<bool, TString> save_sample_){

  int  new_dir = mkdir( ( TString::Format("Plots_distributions/") ).Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  );

  ofstream ratios;
  ratios.open("Ratios.txt", ios::app);
  
  //== Loop.
  

  for(int mc_ = 0 ; mc_ < vector_sample_.size(); mc_++){  
    bool same_sample_ = vector_sample_[mc_];
    TString save_mc_ = save_sample_[ same_sample_ ];
    TString legend_mc_ = legends_sample_[ same_sample_ ];

    TString drawoptions = "hist";
        
    //== Misses canvas.
    TCanvas *can_;  
    PrepareCanvas( can_, "can_comparison_calibrations_fakes" );
    TPad* pad_abs_, *pad_ratio_;
    SplitCanvas_equally( can_, pad_abs_, pad_ratio_ ); 
      
    //== Legend.
    double ticksize = gStyle->GetTickLength();
    
    TLegend *leg = new TLegend( pad_abs_->GetLeftMargin() + ticksize, ticksize, 0.5, 0.55);
    leg->SetFillStyle( 0 );
    leg->SetBorderSize( 0 );  
      
    int distribution_nr = 1;
      
    for(int s_ = 0; s_ < vector_sel_.size(); s_++){  
      TString sel_ = vector_sel_[s_];              
      
      for(int c_ = 0; c_ < vector_cal_.size(); c_++){
        TString cal_ = vector_cal_[ c_ ];


        //== Open MC, extract RooUnfold.  
        TFile *_file0;
        if( !same_sample_) _file0 = new TFile( "/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo" + sel_ + "_unfold_Emin_150.000000_deltaPhiMax_0.500000_etaband_0.000000_all_matchE" + cal_ + ".root", "read");
        if( same_sample_) _file0 = new TFile( "/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo" + sel_ + "_unfold_Emin_150.000000_deltaPhiMax_0.500000_etaband_0.000000_all_matchE.root", "read");
        if( !_file0 ) { continue; }

        RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response");
        TH2D* _hResponse = (TH2D*)response->Hresponse();
      
        //== Fakes.
        can_->cd();
        TH1D* _hFake = (TH1D*)_file0->Get("hCastorJet_fake_all");	Empty_UOflow( _hFake );
        TH1D* _hFakeres = (TH1D*) response->Hfakes();			Empty_UOflow( _hFakeres );
        TH1D* _hMeas = (TH1D*)response->Hmeasured();    		Empty_UOflow( _hMeas );
        TH1D* _hDet = (TH1D*)_file0->Get("hCastorJet_energy");   	Empty_UOflow( _hDet );
               
  	ratios << "\n\t" << _file0->GetName() << endl;        
  	        
        ratios << "Number of fakes\t" << _hFake->Integral() << endl;
        ratios << "Number of dets\t" << _hDet->Integral() << endl;
       
        ratios << "Fakes\t" << _hFake->Integral()/_hDet->Integral() << endl;             
  /*
        SetDnDx( _hDet );
        SetDnDx( _hMeas );
        SetDnDx( _hFake );
  */
  
        _hFake->GetXaxis()->SetTitle("E [GeV]");	
        _hFake->GetXaxis()->SetRangeUser(150, 1600.);
        _hFake->GetYaxis()->SetTitle("dN/dE [1/GeV]");
        _hFake->GetYaxis()->SetRangeUser( 2, 2e5);
        _hFake->SetLineColor( getColor( distribution_nr ) );
        _hFake->SetLineStyle( distribution_nr++  );
        _hFake->SetLineWidth( 2 );
        Prepare_1Dplot( _hFake );
        leg->AddEntry( _hFake, "Fakes: " + legends_sel_[ sel_ ] , "l");        
        pad_abs_->cd();
        _hFake->DrawClone( drawoptions );

        _hMeas->SetLineStyle( distribution_nr  );        
        _hMeas->SetLineColor( getColor( distribution_nr++ ) );
        _hMeas->Draw("histsame");     
        leg->AddEntry( _hMeas, "Meas.: " + legends_sel_[ sel_ ] , "l");         
      
        pad_ratio_->cd();
        _hFake->Divide( _hDet );
        _hFake->GetXaxis()->SetTitleOffset( _hFake->GetXaxis()->GetTitleOffset() * 2. );
        _hFake->GetYaxis()->SetTitle("Fakes/Det. jets"); 
        _hFake->GetYaxis()->SetRangeUser(1e-3, 1.2);
        _hFake->Draw( drawoptions );
        
        drawoptions = "histsame";

      }
    }
    can_->cd();
    pad_abs_->SetLogy();
    pad_abs_->cd();
//    leg->SetHeader("Det. level");
    leg->Draw();
    pad_ratio_->SetLogy();
    can_->SaveAs( "Plots_distributions/can_comparison_calibrations_fakes" + save_mc_ + ".pdf" );                 
  }
}
