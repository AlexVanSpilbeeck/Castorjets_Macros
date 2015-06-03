//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TGraph.h>
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
#include <TLine.h>
#include <TPaletteAxis.h>

//STANDARD C++ INCLUDES
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
#include <vector>

#include "../../../RooUnfold-1.1.1/src/RooUnfold.h"
#include "../../../RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"
#include "../../../RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "../../../RooUnfold-1.1.1/src/RooUnfold.h"
#include "../../../RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "color.h"
#include "Function_Rebin.h"
//#include "Function_make_Tex.h"
//#include "Function_FirstPlot.h"

//#define normalise_1 true

// Define the plotter class.
class Unfolder{
  public:
   Unfolder(vector<TString> MC_files, TString datafile, TString folder, int normalise); // list of MC files - datafile - map to store histograms.
   ~Unfolder();

   // Get histograms from files.
   void Get_DetEnergy(int file_, TH1D* &hist_);
   void Get_DetEnergy_lead(int file_, TH1D* &hist_);
   void Get_Distribution(int file_, TH1D* &hist_, TString variable);
   void Hist_DetLevel();
   void Hist_DetLevel_lead();
   void Get_GenEnergy(int file_, TH1D* &hist_);
   void Get_GenEnergy_lead(int file_, TH1D* &hist_); 
   void Hist_GenLevel();   
   void Hist_GenLevel(TPad* &pad_, bool isFirst);

   void Hist_getTruth();

   // Unfolding.
   void Plot_Unfolded();
   void Plot_Unfolded_Ratio();
//   void Get_DetUnfolded(int file_, TH1D* &hist_, int iterations = 4);
//   void Get_GenSmeared(int file_, TH1D* &hist_, TH1D* hGen);

   void Get_DetUnfolded(int file_, TH1D* &hist_, int iterations = 4, TString variable = "all");
   void Get_GenSmeared(int file_, TH1D* &hist_, TH1D* hGen, TString variable = "all");

   void Hist_DetLevel_DataWithSmear();
   void Hist_DetLevel_MCWithSmear();

   // Helper functions to prepare canvas, titles, labels, ...
   void PrepareCanvas( TCanvas* &can_, TString label);
   void PrepareLegend( map<TString, TString> entries );
   void PrepareTitles( map<TString, TString> xtitle ,  map<TString, TString> ytitle ,  map<TString, TString> htitle);
   void LabelPlots( TString label );
 
   void DoublePaddedComparison_unfolding(TString variable);
   void Plot_Unfolded(TPad* &pad_, TString variable);
   void Plot_Unfolded_Ratio(TPad* &pad_, TString variable);

   void SplitCanvas(TCanvas* &can_, TPad* &pad_abs_, TPad* &pad_ratio_);

   void DoublePaddedComparison(TString variable);
   void Plot_Absolute(TPad* &pad_, TString variable);
   void Plot_Ratio(TPad* &pad_, TString variable);
   void Plot_Absolute(TPad* &pad_, vector<TH1D*>);
   void Plot_Ratio(TPad* &pad_, vector<TH1D*>);

   void Systematics_CompareDetLevel();
   void Systematics_CompareGenLevel();
   void CompareGenLevel();

   void ClosureTest_data(TString variable = "lead");
   void ClosureTest_MC(TString variable = "lead");
   void ClosureTest_MC_detLevel(TString variable = "lead");

  private:
   vector<TString> MC_files_;
   TString datafile_;
   TString folder_; 
   TString label_;
   bool normalise_1;
   map<TString, TString> legend_info_;
   map<TString, TString> xtitle_;
   map<TString, TString> ytitle_;
   map<TString, TString> htitle_;

   void MakeDoublePaddedComparison(TCanvas * &can, vector<TH1D*>, TLegend *leg);
   double Chi2_test( TH1D* hist_data, TH1D* hist_MC);
   double Det_to_gen_scale(int file);

};

Unfolder::Unfolder( vector<TString> MC_files, TString datafile, TString folder, int normalise )
{
  gStyle->SetStatStyle(0);
  gStyle->SetTitleStyle(0);

  MC_files_ = MC_files;
  datafile_ = datafile;
  folder_ = "Plots/" + folder;
  
  int  new_dir = mkdir( ("Plots/" + folder).Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  );

  normalise_1 = false;
  if (normalise == 1 ){ normalise_1 = true; }
}
Unfolder::~Unfolder(){
}


/*****************************************
* Plot the Det level energy distribution *
*****************************************/

// -- Plot
void Unfolder::Get_DetEnergy(int file_, TH1D* &hist_){
   cout << "\n\n\n// -- Detector level energy - file selection" << endl;

   TH1D *hDet;
   TFile *_file0;
   TString drawoptions = "";
   
   if( file_ < MC_files_.size() ){	_file0 = new TFile( MC_files_[file_], "Read"); 	drawoptions = "hist";}
   else if( file_ == -1 ){		_file0 = new TFile( datafile_, "Read");		drawoptions = "data";}
         
   hDet = (TH1D*)_file0->Get("hCastorJet_energy");	hDet->Sumw2();
   hDet->GetXaxis()->SetTitle("E_{det} (GeV)");
   hDet->GetYaxis()->SetTitle("#frac{dN}{dE}");
  
   if( normalise_1 ) hDet->Scale( 1./hDet->Integral() );
 
   hist_ = hDet;
}



void Unfolder::Get_DetEnergy_lead(int file_, TH1D* &hist_){
   TH1D *hDet;
   TFile *_file0;
   TString drawoptions = "";

   if( file_ < MC_files_.size() ){      _file0 = new TFile( MC_files_[file_], "Read");  drawoptions = "hist";}
   else if( file_ == -1 ){              _file0 = new TFile( datafile_, "Read");         drawoptions = "data";}

   hDet = (TH1D*)_file0->Get("hCastorJet_energy_lead");      hDet->Sumw2();
   hDet->GetXaxis()->SetTitle("E_{det} (GeV)");
   hDet->GetYaxis()->SetTitle("#frac{dN}{dE}");

   if( normalise_1 ) hDet->Scale( 1./hDet->Integral() );

   hist_ = hDet;
}



void Unfolder::Get_Distribution(int file_, TH1D* &hist_, TString variable){

   TH1D *hVar;
   TFile *_file0;
  
   if( file_ < MC_files_.size() ){      _file0 = TFile::Open( MC_files_[file_], "Read"); }
   else if( file_ == -1 ){              _file0 = TFile::Open( datafile_, "Read");        }   

   hVar = (TH1D*)_file0->Get( variable );      hVar->Sumw2();

   hist_ = (TH1D*)hVar->Clone(variable);
}

////////////////////////////////////
// Plot energy at detector level. //
////////////////////////////////////

void Unfolder::Hist_DetLevel_lead(){

  TH1D* hist_;
  TString drawoptions = "hist", legendoptions = "l";
  TCanvas *can; 
  PrepareCanvas( can, "Detector_Level_Energy_lead" );
  vector<TH1D*> histos;
  TLegend*leg = new TLegend( 0.65, 0.75, 0.95, .95);
  bool done_data = false;

  double max_ = 0.;

  for( int file = 0; file <= MC_files_.size(); file++ ){
    int file_ = file;
    if( file == MC_files_.size() ){
      file_ = -1;
      drawoptions = "datasame";
      legendoptions = "p";
      done_data = true;
    }

    Get_DetEnergy_lead( file_ , hist_ );
    hist_->SetLineColor( getColor( file_ +1) );
    hist_->SetLineWidth( 3 );
    hist_->SetMarkerColor( getColor( file_ +1) );
    hist_->SetMarkerSize( 1 );

    if( file_ != -1) {histos.push_back( hist_ ); }
    if( file_ == -1) {histos.insert( histos.begin(), hist_) ; }

    TString legend_info;
    if( file_ == -1 ){ legend_info = legend_info_[datafile_]; }
    else{ legend_info = legend_info_[MC_files_[file_]];  }

    leg->AddEntry( hist_, legend_info, legendoptions );
    drawoptions = "histsame";
    if( done_data ){ break; }
  }
   MakeDoublePaddedComparison(can, histos, leg );
 
  can->SaveAs(folder_ + "/DetectorLevelEnergy_lead" + label_ + ".C");
  can->SaveAs(folder_ + "/DetectorLevelEnergy_lead" + label_ + ".pdf");
}

////////////////////////
// Leading jets only. //
////////////////////////

void Unfolder::Hist_DetLevel(){

  TH1D* hist_;
  TString drawoptions = "hist", legendoptions = "l";
  TCanvas *can; 
  PrepareCanvas( can, "Detector_Level_Energy" );
  vector<TH1D*> histos;
  TLegend*leg = new TLegend( 0.65, 0.75, 0.95, .95);
  bool done_data = false;

  double max_ = 0.;

  for( int file = 0; file <= MC_files_.size(); file++ ){
    int file_ = file;
    if( file == MC_files_.size() ){
      file_ = -1;
      drawoptions = "datasame";
      legendoptions = "p";
      done_data = true;
    }

    Get_DetEnergy( file_ , hist_ );
    hist_->SetLineColor( getColor( file_ +1) );
    hist_->SetLineWidth( 3 );
    hist_->SetMarkerColor( getColor( file_ +1) );
    hist_->SetMarkerSize( 1 );

    if( file_ != -1) {histos.push_back( hist_ ); }
    if( file_ == -1) {histos.insert( histos.begin(), hist_) ; }

    TString legend_info;
    if( file_ == -1 ){ legend_info = legend_info_[datafile_]; }
    else{ legend_info = legend_info_[MC_files_[file_]];  }

    leg->AddEntry( hist_, legend_info, legendoptions );
    drawoptions = "histsame";
    if( done_data ){ break; }
  }
   MakeDoublePaddedComparison(can, histos, leg );
 
  can->SaveAs(folder_ + "/DetectorLevelEnergy" + label_ + ".C");
  can->SaveAs(folder_ + "/DetectorLevelEnergy" + label_ + ".pdf");
}







void Unfolder::Hist_getTruth(){
  int file_ = 0;
  TFile* _file0 = TFile::Open( MC_files_[0], "Read");
  RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response");

  TH1D* hTruth;
  hTruth = (TH1D*) response->Htruth();
  cout << "Truth\t" << hTruth->Integral() << endl;

}



/****************************************************************************
* Compare detector level data with unfolded -> smeared data to get closure. *
****************************************************************************/

void Unfolder::Hist_DetLevel_DataWithSmear(){
  cout << "\n// -- Detector level energy on a single plot" << endl;

  int file_ = -1;

  TH1D* hist_, *hist_unfold_, *hist_smeared_;
  TString drawoptions = "hist";
  TCanvas *can;
  PrepareCanvas( can, "Detector_Level_Energy" );

  Get_DetEnergy( file_ , hist_ );
  hist_->SetLineColor( getColor( file_ +1) );
  hist_->SetLineWidth( 3 );
  hist_->SetMarkerColor( getColor( file_ +1) );
  hist_->SetMarkerSize( 1 );
  hist_->Draw(drawoptions);

  drawoptions = "histsame";

  Get_DetUnfolded( file_, hist_unfold_ );
  
  Get_GenSmeared( file_, hist_smeared_ , hist_unfold_);
  hist_smeared_->SetLineColor( getColor( file_+2) );
  hist_smeared_->SetLineWidth(3);
  hist_smeared_->SetMarkerColor( getColor( file_ +2) );
  hist_smeared_->SetMarkerSize( 1 );

  hist_smeared_->Draw(drawoptions);

  can->SaveAs(folder_ + "/DetectorLevelData" + label_ + ".C");
  can->SaveAs(folder_ + "/DetectorLevelData" + label_ + ".pdf");
}




/*************************************************************************
* Compare detector level MC with -> smeared gen level MC to get closure. *
*************************************************************************/

void Unfolder::Hist_DetLevel_MCWithSmear(){

  int file_ = 0;

  TH1D* hist_, *hist_unfold_, *hist_smeared_;
  TString drawoptions = "hist";
  TCanvas *can; 
  PrepareCanvas( can, "Detector_Level_Energy" );

  Get_DetEnergy( file_ , hist_ );
  hist_->SetLineColor( getColor( file_ +1) );
  hist_->SetLineWidth( 3 );
  hist_->SetMarkerColor( getColor( file_ +1) );
  hist_->SetMarkerSize( 1 );
  hist_->Draw(drawoptions);

  drawoptions = "histsame";

  Get_GenEnergy( file_, hist_unfold_ );
  
  Get_GenSmeared( file_, hist_smeared_ , hist_unfold_);


  hist_smeared_->SetLineColor( getColor( file_+2) );
  hist_smeared_->SetLineWidth(3);
  hist_smeared_->SetMarkerColor( getColor( file_ +2) );
  hist_smeared_->SetMarkerSize( 1 );

  hist_smeared_->Draw(drawoptions);

  can->SaveAs(folder_ + "/DetectorLevelMC" + label_ + ".C");
  can->SaveAs(folder_ + "/DetectorLevelMC" + label_ + ".pdf");
}





/*****************************************
* Plot the Gen level energy distribution *
*****************************************/

// -- Plot
void Unfolder::Get_GenEnergy(int file_, TH1D* &hist_){
   TH1D *hGen;
   TFile *_file0;
   TString drawoptions = "";
   
   if( file_ < MC_files_.size() ){	_file0 = new TFile( MC_files_[file_], "Read"); 	drawoptions = "hist";}
   
   hGen = (TH1D*)_file0->Get("hGenJet_energy");	hGen->Sumw2();
   hGen->GetXaxis()->SetTitle("E_{det} (GeV)");
   hGen->GetYaxis()->SetTitle("#frac{dN}{dE}");

   if( normalise_1 ) hGen->Scale( 1./hGen->Integral() );
   
   hist_ = hGen;
}

// -- Plot
void Unfolder::Get_GenEnergy_lead(int file_, TH1D* &hist_){
   TH1D *hGen;
   TFile *_file0;
   TString drawoptions = "";
   
   if( file_ < MC_files_.size() ){	_file0 = new TFile( MC_files_[file_], "Read"); 	drawoptions = "hist";}
   
   hGen = (TH1D*)_file0->Get("hGenJet_energy_lead");	hGen->Sumw2();
   hGen->GetXaxis()->SetTitle("E_{det} (GeV)");
   hGen->GetYaxis()->SetTitle("#frac{dN}{dE}");

   if( normalise_1 ) hGen->Scale( 1./hGen->Integral() );
   
   hist_ = hGen;
}
void Unfolder::Hist_GenLevel(){

  TH1D* hist_;
  TString drawoptions = "hist";
  TCanvas *can; 
  PrepareCanvas( can, "Generator_Level_Energy" );
  TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.95);

  for( int file = 0; file < MC_files_.size(); file++ ){
    int file_ = file;
    Get_GenEnergy( file_ , hist_ );
    hist_->SetLineColor( getColor( file_ +1) );
    hist_->SetLineWidth( 3 );
    hist_->SetMarkerColor( getColor( file_ +1) );
    hist_->SetMarkerSize( 1 );
    hist_->Scale( 1./hist_->Integral() );
    
    hist_->Draw(drawoptions);
    drawoptions = "histsame";  

    leg->AddEntry( hist_, TString::Format("MC %i", file), "l");
  }

  leg->Draw();  

  can->SaveAs(folder_ + "/GeneratorLevelEnergy" + label_ + ".C");
  can->SaveAs(folder_ + "/GeneratorLevelEnergy" + label_ + ".pdf");
  
}

void Unfolder::Hist_GenLevel(TPad* &pad_, bool isFirst){
  pad_->cd();

  cout << "\t// -- GEN LEVEL" << endl;

  TH1D* hist_;
  TString drawoptions = "hist";	if( !isFirst ) drawoptions ="histsame";

  for( int file = 0; file < MC_files_.size(); file++ ){
    int file_ = file;
    Get_GenEnergy( file_ , hist_ );
    hist_->SetLineColor( getColor( file_ +1) );
    hist_->SetLineWidth( 3 );
    hist_->SetMarkerColor( getColor( file_ +1) );
    hist_->SetMarkerSize( 1 );
    hist_->Scale( 1./hist_->Integral() );

    hist_->Draw(drawoptions);
    drawoptions = "histsame";
  }
}

/************************************************************************
* Plot generator level, unfolded detector and data next to one another. *
************************************************************************/

void Unfolder::Plot_Unfolded(){

  int color_ = 1;

  TH1D* hGen;
    Get_GenEnergy(0, hGen);
    hGen->SetLineColor( getColor(color_++) );
    hGen->SetLineStyle( 1 );
  
  TH1D* hDet;
    Get_DetEnergy( 0 , hDet );
    Get_DetUnfolded(0, hDet, 4);
    hDet->SetLineColor( getColor(color_++) );
    hDet->SetLineStyle( 2 );

  
  TH1D* hData;
    Get_DetEnergy( -1 , hData );
    Get_DetUnfolded(-1, hData, 4);
    hData->SetLineColor( getColor( 0 ) );
    hData->SetMarkerColor( getColor(0) );
  
  TCanvas *can;
  PrepareCanvas( can , "Comparison_Energy_Unfolded");

  TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.95);  

  double min_val, max_val;
  
  hGen->	Draw("hist");
    hGen->	Scale( 1./hGen->Integral() );
    First_Plot( hGen, hGen, 0, min_val, max_val);
    leg->	AddEntry( hGen, "MC (Gen)", "l");

  hDet->	Draw("histsame");
    hDet->	Scale( 1./hDet->Integral() );
    First_Plot( hGen, hDet, 1, min_val, max_val);
    leg->	AddEntry( hDet, "MC (Det)", "l");
    
  hData->	Draw("edatasame");
    hData->	Scale( 1./hData->Integral() );
    First_Plot( hGen, hData, 2, min_val, max_val);
    leg->	AddEntry( hData, "Data", "p");
    
  leg->Draw();
  can->SetLogy();
  can->Update();
  
  can->SaveAs(folder_ + "/Compare_UnfoldedEnergy" + label_ + ".C");
  can->SaveAs(folder_ + "/Compare_UnfoldedEnergy" + label_ + ".pdf");    
}



void Unfolder::Plot_Unfolded_Ratio(){

  int color_ = 1;

  TH1D* hGen;
    Get_GenEnergy(0, hGen);
    hGen->SetLineColor( getColor(color_++) );
    hGen->SetLineStyle( 1 );
  
  TH1D* hDet;
    Get_DetEnergy(0, hDet);
    Get_DetUnfolded(0, hDet, 4);
    hDet->SetLineColor( getColor(color_++) );
    hDet->SetLineStyle( 2 );
  
  TH1D* hData;
    Get_DetEnergy(-1, hData);
    Get_DetUnfolded(-1, hData, 4);
    hData->SetLineColor( getColor( 0 ) );
    hData->SetMarkerColor( getColor(0) );
  
  TCanvas *can;
  PrepareCanvas( can , "Comparison_Energy_Unfolded");
  double min_val, max_val;
  
  hGen->	Scale( 1./hGen->Integral() );    
    hGen->	Draw("hist");
    First_Plot( hGen, hGen, 0, min_val, max_val);
    
  hDet->	Scale( 1./hDet->Integral() );
    hDet->	Divide( hGen );
    hDet->	Draw("histsame");
    First_Plot( hGen, hDet, 1, min_val, max_val);

  hData->	Scale( 1./hData->Integral() );
    hData->	Divide( hGen );
    hData->	Draw("edatasame");
    First_Plot( hGen, hData, 2, min_val, max_val);
    
    hGen->	Divide( hGen );
    First_Plot( hGen, hGen, 0, min_val, max_val);    
    
    cout << "Min and max\t" << min_val << "\t" << max_val << endl;
    
    can->Update();
  
  can->SaveAs(folder_ + "/Compare_UnfoldedEnergy_Ratio" + label_ + ".C");
  can->SaveAs(folder_ + "/Compare_UnfoldedEnergy_Ratio" + label_ + ".pdf");    
}






void Unfolder::Plot_Absolute( TPad* & pad_, TString variable){
  cout << "\tPlot_Absolute\t" << variable << endl;

  double min_val, max_val;
  TH1D* hDistr, *hFirst;
  bool firstDistr = true, passed_data = false;
  int nDistr = 0;
  TString drawoptions = "phist";
  TString legendoptions = "l";

  TLegend *leg = new TLegend(0.6, 0.7, 0.95, 0.95  );

  cout << "Prepare for loop" << endl;
  for(int file_ = 0; file_ <= MC_files_.size(); file_++){
    cout << "Loop iteration\t" << file_ << endl;

    TString filename;
    if( (file_-1) < MC_files_.size() && file_ > 0){ filename = MC_files_[ (file_-1) ]; }
    if( (file_-1) == -1 ){ filename = datafile_; legendoptions = "p";}
   
    TFile* _file = TFile::Open( filename, "read" );   
    if(! _file->GetListOfKeys()->Contains( variable ) ){ continue; }

    Get_Distribution( (file_-1) , hDistr, variable );
    hDistr->Rebin( 4 );
    hDistr->GetXaxis()->SetTitle( xtitle_[variable] );
    hDistr->GetYaxis()->SetTitle( ytitle_[variable] );
    hDistr->GetYaxis()->SetTitleOffset( 0.85 );
    hDistr->SetTitle( htitle_[variable] );
    hDistr->GetXaxis()->SetNdivisions( 504 );
 
    if( firstDistr ){ 
      hFirst = (TH1D*)hDistr->Clone(variable + "_first");
      hDistr = hFirst;
      firstDistr = false;
    }

    hDistr->Scale( 1./hDistr->Integral() );
    First_Plot( hFirst, hDistr, nDistr, max_val, min_val);
    
    hDistr->GetXaxis()->SetRangeUser(0., 2000.);
    leg->AddEntry( hDistr, legend_info_[filename], legendoptions );

    hDistr->SetLineColor( getColor( nDistr+1 ) );
    hDistr->SetLineStyle( nDistr+1 );
    hDistr->SetLineWidth( 3 );
    hDistr->SetMarkerColor(  getColor( nDistr+1 ) );
  
    pad_->cd();
    pad_->SetLogy();
    hDistr->Draw( drawoptions );
   
    drawoptions = "histsame";
    nDistr++;
//    leg->Draw();
    
    if( nDistr > 0 ) pad_->Update();
    if( file_ == MC_files_.size() ){
      leg->Draw();
      cout << "\t\t\t\t\t\t\tLegend drawn" << endl;
      break;
    }
  }
}




void Unfolder::Plot_Ratio( TPad* & pad_, TString variable){
  double min_val, max_val;
  TH1D* hDistr, *hFirst, *hFirst_abs;
  bool firstDistr = true, passed_data = false;
  int nDistr = 0;
  TString drawoptions = "phist";

  for(int file_ = 0; file_ <= MC_files_.size(); file_++){

    TString filename;
    if( (file_-1) < MC_files_.size() ) filename = MC_files_[ (file_-1) ];
    if( (file_-1) == -1){ filename = datafile_; }

    TFile* _file = TFile::Open( filename, "read" );
    if(! _file->GetListOfKeys()->Contains( variable ) ){ continue; }

    Get_Distribution( (file_-1) , hDistr, variable );
    hDistr->Rebin( 4 );
    hDistr->GetXaxis()->SetTitle( xtitle_[variable] );
    hDistr->GetXaxis()->SetNdivisions( 504 );
    hDistr->Scale( 1./hDistr->Integral() );

    if( firstDistr ){
      hFirst = (TH1D*)hDistr->Clone(variable + "_first");
      hFirst_abs = (TH1D*)hFirst->Clone(variable + "_first_absolute");
      hFirst->Divide( hFirst );
      firstDistr = false;
    }
   

    hDistr->Divide( hFirst_abs );
    First_Plot( hFirst, hDistr, nDistr, max_val, min_val);

    hDistr->GetXaxis()->SetRangeUser(0., 2000.);
    hDistr->GetYaxis()->SetRangeUser(0., 3.);

    hDistr->SetLineColor( getColor( file_+1 ) );
    hDistr->SetLineWidth( 3 );
    hDistr->SetLineStyle( file_+1 );
    hDistr->SetMarkerColor(  getColor( file_+1 ) );


    pad_->cd();
//    pad_->SetLogy();
    hDistr->Draw( drawoptions );
    hFirst->GetYaxis()->SetRangeUser(0., 3.);
 
    drawoptions = "histsame";
    nDistr++;

    if( nDistr > 0 ) pad_->Update();
    if( passed_data ) break;
  }


}


void Unfolder::Get_DetUnfolded(int file_, TH1D* &hist_, int iterations, TString variable){

  TFile *_file0 = new TFile( MC_files_[0], "read");
  RooUnfoldResponse* response;
  if( variable == "all"){ response = (RooUnfoldResponse*)_file0->Get("response"); }
  if( variable == "lead"){ response = (RooUnfoldResponse*)_file0->Get("response_lead"); }
  
  // Use hist as input for the unfolding.
  RooUnfoldBayes unfold_bayes(response, hist_, iterations);
  // the original hist has become obsolete, transform it into the unfolded histogram.
  hist_ = (TH1D*) unfold_bayes.Hreco(); 
}




void Unfolder::Get_GenSmeared(int file_, TH1D* &hist_, TH1D* hGen, TString variable){

  TFile *_file0 = new TFile( MC_files_[0], "read");
  RooUnfoldResponse* response;
  if( variable == "all"){ response = (RooUnfoldResponse*)_file0->Get("response"); }
  if( variable == "lead"){ response = (RooUnfoldResponse*)_file0->Get("response_lead"); }


  hist_ = (TH1D*) response->ApplyToTruth( hGen, "Smearing");
}




/**********************************
* Set the right canvas properties *
**********************************/

void Unfolder::PrepareCanvas( TCanvas* &can_, TString label){

  cout << "\n\n\n// --  Prepare canvas" << endl;

  can_ = new TCanvas( label, label, 900, 900 );

  can_->SetLeftMargin(0.20);
  can_->SetTopMargin(0.07);
  can_->SetBottomMargin(0.14);
  can_->SetRightMargin(0.05);


}

void Unfolder::PrepareLegend( map<TString, TString> legend_info ){

  legend_info_ = legend_info;

}

void Unfolder::PrepareTitles( map<TString, TString> xtitle, map<TString, TString> ytitle, map<TString, TString> htitle ){

  xtitle_ = xtitle;
  ytitle_ = ytitle;
  htitle_ = htitle;

}


/***********************************************
* -- Compare unfolded plots on a split canvas. *
***********************************************/



void Unfolder::DoublePaddedComparison_unfolding(TString variable){

  TCanvas* can_;
  TPad* pad_abs_, *pad_ratio_;
  PrepareCanvas(can_, "Comparison_UnfoldedEnergy_Data_and_MC");

  SplitCanvas(can_, pad_abs_, pad_ratio_);

  Hist_GenLevel(pad_abs_, true);
  Plot_Unfolded(pad_abs_, variable);
  Plot_Unfolded_Ratio(pad_ratio_, variable);

  can_->SaveAs( folder_ + "/UnfoldedEnergy_Data_and_MC_" + variable + "_" + label_ + ".pdf");
  can_->SaveAs( folder_ + "/UnfoldedEnergy_Data_and_MC_" + variable + "_" + label_ + ".C");
}

//
// -- Absolute distributions.
//


void Unfolder::Plot_Unfolded(TPad* & pad_, TString variable){
  pad_->cd();

  TLegend *leg = new TLegend(0.6, 0.50, 0.95, 0.95);

  int color_ = 1;
  double min_val, max_val;

  TH1D* hGen;
    if( variable == "all") { Get_GenEnergy(0, hGen); }
    if( variable == "lead") { Get_GenEnergy_lead(0, hGen); }
    hGen->SetLineColor( 1  );
    hGen->SetLineStyle( 1 );
    hGen->SetLineWidth( 2 );
    leg->AddEntry( hGen, "MC (Gen)", "l");

  hGen->        Draw("hist");
    First_Plot( hGen, hGen, 0, min_val, max_val);


  // Determine the normalisation of the detector level distribution: MC(Gen)/MC(Det)
  // For this we need the unscaled MC distributions.
  double det_to_gen = Det_to_gen_scale(0);

  TH1D* hDet;
    if( variable == "all") { Get_DetEnergy(0, hDet);}
    if( variable == "lead"){ Get_DetEnergy_lead(0, hDet);}
    hDet->Scale( det_to_gen);
    Get_DetUnfolded(0, hDet, 4, variable);
    hDet->SetLineColor( getColor(color_++) );
    hDet->SetLineStyle( 2 );
    hDet->SetLineWidth( 2 );
    leg->AddEntry( hDet, "MC (Det)", "l");
    hDet->        Draw("histsame");
    First_Plot( hGen, hDet, 1, min_val, max_val);

  TH1D* hData;

  for(int iterations_ = 1; iterations_ < 16; iterations_++){
    if( variable == "all") { Get_DetEnergy(-1, hData); }
    if( variable == "lead"){ Get_DetEnergy_lead(-1, hData); }
    hData->Scale( det_to_gen);
    Get_DetUnfolded(-1, hData, iterations_, variable);
    hData->SetLineColor( getColor( color_) );
    hData->SetMarkerColor( getColor(color_++) );
    hData->SetMarkerStyle( 24 );
    leg->AddEntry( hData, TString::Format("Data %i it.", iterations_), "p");
    hData->       DrawClone("edatasame");

    First_Plot( hGen, hData, 2, min_val, max_val);
  }
  
  leg->Draw();

  pad_->SetLogy();
  pad_->Update();
}

//
// -- Relative distributions.
//

void Unfolder::Plot_Unfolded_Ratio(TPad* & pad_, TString variable){
  pad_->cd();

  int color_ = 1;



  TH1D* hGen;
    if( variable == "all") { Get_GenEnergy(0, hGen); }
    if( variable == "lead") { Get_GenEnergy_lead(0, hGen); }

    hGen->SetLineColor( 1  );
    hGen->SetLineStyle( 1 );
    hGen->SetLineWidth( 2 );

  TH1D* hFirst = (TH1D*)hGen->Clone("Reference_histogram");
  hGen->Divide(hFirst);
  hGen->        Draw("hist");

  // Determine the normalisation of the detector level distribution: MC(Gen)/MC(Det)
  // For this we need the unscaled MC distributions.
  double det_to_gen = Det_to_gen_scale(0);

  TH1D* hDet;
    if( variable == "all") { Get_DetEnergy(0, hDet);}
    if( variable == "lead"){ Get_DetEnergy_lead(0, hDet);}
    hDet->Scale( det_to_gen);
    Get_DetUnfolded(0, hDet, 4, variable);
    hDet->SetLineColor( getColor(color_++) );
    hDet->SetLineStyle( 2 );
    hDet->SetLineWidth( 2 );
    hDet->Divide( hFirst);
    hDet->        Draw("histsame");

  TH1D* hData;

  for(int iterations_ = 1; iterations_ < 16; iterations_++){
    if( variable == "all") { Get_DetEnergy(-1, hData); }
    if( variable == "lead"){ Get_DetEnergy_lead(-1, hData); }
    hData->Scale( det_to_gen);
    Get_DetUnfolded(-1, hData, iterations_, variable);
    hData->SetLineColor( getColor( color_) );
    hData->SetMarkerColor( getColor(color_++) );
    hData->SetMarkerStyle( 24 );
    hData->Divide( hFirst);
    hData->       DrawClone("edatasame");

  }

  hGen->GetYaxis()->SetRangeUser(0., 3.);
  pad_->Update();
}




////////////////////////////////////////////////////
// Absolute and relative values for any variable. //
////////////////////////////////////////////////////

void Unfolder::DoublePaddedComparison(TString variable){
  cout << "--- Double padded comparison" << endl;

  TCanvas* can_;
  TPad* pad_abs_, *pad_ratio_;
  PrepareCanvas(can_, "Comparison_UnfoldedEnergy_Data_and_MC");

  SplitCanvas(can_, pad_abs_, pad_ratio_);

  Plot_Absolute(pad_abs_, variable);
  Plot_Ratio(pad_ratio_, variable);

  can_->SaveAs( folder_ + "/Comparison_" + variable + "_data_and_MC" + label_ + ".pdf");
  can_->SaveAs( folder_ + "/Comparison_" + variable + "_data_and_MC" + label_ + ".C");
}




void Unfolder::SplitCanvas(TCanvas* &can_, TPad* &pad_abs_, TPad* &pad_ratio_){

  can_->cd();

  double left_ 	= 0.20;
  double top_ 	= 0.05;
  double right_ = 0.05;
  double bottom_= 0.15;

  pad_abs_ = new TPad("Absolute_Values", "Absolute_Values", 0., 0.4,1.,1.);
  pad_abs_->SetLeftMargin(left_);
  pad_abs_->SetRightMargin(right_);
  pad_abs_->SetBottomMargin(0.);
  pad_abs_->SetTopMargin(top_);
  pad_abs_->Draw();

  can_->cd();
  pad_ratio_ = new TPad("Ratios", "Ratios", 0., 0., 1., 0.4);
  pad_ratio_->SetLeftMargin(left_);
  pad_ratio_->SetRightMargin(right_);
  pad_ratio_->SetTopMargin(0);
  pad_ratio_->SetBottomMargin(bottom_);
  pad_ratio_->Draw();
  
}



void Unfolder::Systematics_CompareDetLevel(){

  TH1D* hist_;
  TString drawoptions = "hist", legendoptions = "l";
  TCanvas *can; 
  PrepareCanvas( can, "Systematics_Detector_Level_Energy" );
  vector<TH1D*> histos;
  TLegend *leg = new TLegend(0.65, 0.7, 0.95, 0.95);
  bool done_data = false;

  for( int file = 0; file <= MC_files_.size(); file++ ){
    int file_ = file;
    if( file == MC_files_.size() ){
      file_ = -1;
      drawoptions = "datasame";
      legendoptions = "p";
      done_data = true;
    }

    Get_DetEnergy( file_ , hist_ );
    Rebin_to( hist_, 50 ); 
    hist_->GetXaxis()->SetTitle( xtitle_["hCastorJet_energy"] );
    hist_->GetYaxis()->SetTitle( ytitle_["hCastorJet_energy"] );
    hist_->GetYaxis()->SetTitleOffset( 0.85 );
    hist_->GetYaxis()->SetTitleSize( 0.09 );
    hist_->GetYaxis()->SetLabelSize(0.07);

    hist_->GetXaxis()->SetNdivisions( 504 );
    hist_->SetLineColor( getColor( file_ +1) );
    hist_->SetLineWidth( 3 );
    hist_->SetMarkerColor( getColor( file_ +1) );
    hist_->SetMarkerSize( 1 );

    if( file_ != -1) {histos.push_back( hist_ ); }
    if( file_ == -1) {histos.insert( histos.begin(), hist_) ; }

    TString legend_info;
    if( file_ == -1 ){ legend_info = legend_info_[datafile_]; }
    else{ legend_info = legend_info_[MC_files_[file_]];  }

    leg->AddEntry( hist_, legend_info, legendoptions );
    drawoptions = "histsame";
    if( done_data ){ break; }
  }

  MakeDoublePaddedComparison(can, histos, leg );

  can->SaveAs(folder_ + "/Systematics_DetectorLevelEnergy" + label_ + ".C");
  can->SaveAs(folder_ + "/Systematics_DetectorLevelEnergy" + label_ + ".pdf");
}


void Unfolder::Systematics_CompareGenLevel(){

  TH1D* hist_;
  TString drawoptions = "hist", legendoptions = "l";
  TCanvas *can; 
  PrepareCanvas( can, "Systematics_Generator_Level_Energy" );
  vector<TH1D*> histos;
  TLegend *leg = new TLegend(0.65, 0.7, 0.95, 0.95);
  bool done_data = false;

  double Det_to_gen = Det_to_gen_scale(0);

  for( int file = 0; file <= MC_files_.size(); file++ ){
    int file_ = file;
    if( file == MC_files_.size() ){
      file_ = -1;
      drawoptions = "datasame";
      legendoptions = "p";
      done_data = true;
    }

    Get_DetEnergy( file_ , hist_ );
    Rebin_to( hist_, 50 ); 
    hist_->Scale( Det_to_gen );

    Get_DetUnfolded( file, hist_, 4);

    hist_->GetXaxis()->SetTitle( xtitle_["hGenJet_energy"] );
    hist_->GetYaxis()->SetTitle( ytitle_["hGenJet_energy"] );
    hist_->GetYaxis()->SetTitleOffset( 0.85 );
    hist_->GetYaxis()->SetTitleSize( 0.09 );
    hist_->GetYaxis()->SetLabelSize(0.07);

    hist_->GetXaxis()->SetNdivisions( 504 );
    hist_->SetLineColor( getColor( file_ +1) );
    hist_->SetLineWidth( 3 );
    hist_->SetMarkerColor( getColor( file_ +1) );
    hist_->SetMarkerSize( 1 );

    if( file_ != -1) {histos.push_back( hist_ ); }
    if( file_ == -1) {histos.insert( histos.begin(), hist_) ; }

    TString legend_info;
    if( file_ == -1 ){ legend_info = legend_info_[datafile_]; }
    else{ legend_info = legend_info_[MC_files_[file_]];  }

    leg->AddEntry( hist_, legend_info, legendoptions );
    drawoptions = "histsame";
    if( done_data ){ break; }
  }

  MakeDoublePaddedComparison(can, histos, leg );

  can->SaveAs(folder_ + "/Systematics_GeneratorLevelEnergy" + label_ + ".C");
  can->SaveAs(folder_ + "/Systematics_GeneratorLevelEnergy" + label_ + ".pdf");
}






void Unfolder::CompareGenLevel(){

  TH1D* hist_;
  TString drawoptions = "hist", legendoptions = "l";
  TCanvas *can; 
  PrepareCanvas( can, "Systematics_Generator_Level_Energy" );
  vector<TH1D*> histos;
  TLegend *leg = new TLegend(0.65, 0.7, 0.95, 0.95);
  bool done_data = false;

  double Det_to_gen = Det_to_gen_scale(0);

  for( int file = 0; file <= 1; file++ ){
    int file_ = 1;
    if( file == MC_files_.size() ){
      file_ = -1;
      drawoptions = "datasame";
      legendoptions = "p";
      done_data = true;
    }

    Get_DetEnergy( file_ , hist_ );
    Rebin_to( hist_, 50 ); 
    hist_->Scale( Det_to_gen );

    Get_DetUnfolded( file, hist_, 4);

    hist_->GetXaxis()->SetTitle( xtitle_["hGenJet_energy"] );
    hist_->GetYaxis()->SetTitle( ytitle_["hGenJet_energy"] );
    hist_->GetYaxis()->SetTitleOffset( 0.85 );
    hist_->GetYaxis()->SetTitleSize( 0.09 );
    hist_->GetYaxis()->SetLabelSize(0.07);

    hist_->GetXaxis()->SetNdivisions( 504 );
    hist_->SetLineColor( getColor( file_ +1) );
    hist_->SetLineWidth( 3 );
    hist_->SetMarkerColor( getColor( file_ +1) );
    hist_->SetMarkerSize( 1 );

    if( file_ != -1) {histos.push_back( hist_ ); }
    if( file_ == -1) {histos.insert( histos.begin(), hist_) ; }

    TString legend_info;
    if( file_ == -1 ){ legend_info = legend_info_[datafile_]; }
    else{ legend_info = legend_info_[MC_files_[file_]];  }

    leg->AddEntry( hist_, legend_info, legendoptions );
    drawoptions = "histsame";
    if( done_data ){ break; }
  }

  MakeDoublePaddedComparison(can, histos, leg );

  can->SaveAs(folder_ + "/GeneratorLevelEnergy" + label_ + ".C");
  can->SaveAs(folder_ + "/GeneratorLevelEnergy" + label_ + ".pdf");
}






void Unfolder::ClosureTest_data(TString variable){
  // Prepare variables and objects.
  TH1D* hist_original;
  TH1D* hist_result;  
  TString htitle;
  TLegend *leg = new TLegend(0.65, 0.45, 0.95, 0.95);
  int iterations_ = 15;
  double xaxisgraph[iterations_], yaxisgraph[iterations_];

  vector<TH1D*> histos;
  
  // Get the detector level distribution.
  if( variable == "lead"){	Get_DetEnergy_lead(-1, hist_original );	htitle = "E_{det} spectrum (leading)";}
  else{				Get_DetEnergy(-1, hist_original );	htitle = "E_{det} spectrum (all)";	}

  if( normalise_1 ){ hist_original->Scale( 1./hist_original->Integral() ); }

  // Set detector level properties and legend.
  hist_original->SetLineColor( getColor( 1 ) );
  	hist_original->SetLineStyle( 1 );
  	hist_original->SetLineWidth( 3 );
  	hist_original->SetMarkerColor(  getColor( 1 ) );

  leg->AddEntry( hist_original, "Actual data", "p");
  histos.push_back( hist_original );

  // Loop over number of Bayesian iterations, and unfold->smear detector level distribution.
  // Calculate chi2 as we go.
  for(int iterations = 1; iterations <= iterations_; iterations++){

    cout << "\t\t\tClosure test data - iteration " << iterations << " of variable " << variable << endl;

    hist_result = (TH1D*)hist_original->Clone(TString::Format("Closure_%i_iterations", iterations) );
    hist_result->SetName( TString::Format("Closure_%i_iterations", iterations) );
    hist_result->SetTitle( TString::Format("Closure_%i_iterations", iterations) );
    Get_DetUnfolded( -1, hist_result, iterations, variable );
    Get_GenSmeared( -1, hist_result , hist_result, variable); 
    //hist_int = (TH1D*)hist_result->Clone(TString::Format("it_unfolded_%i", iterations) );
    //list_of_iterations.push_back( hist_int);
    hist_result->SetLineColor( getColor( iterations + 1 ) );
    hist_result->SetLineStyle( iterations + 1 );
    hist_result->SetLineWidth( 3 );
    hist_result->SetMarkerColor(  getColor( iterations+1 ) );
    leg->AddEntry( hist_result, TString::Format("Treated data, %i it.", iterations) , "p");
    histos.push_back( hist_result );


    // Chi2 test.
    double chi2 = Chi2_test( hist_original, hist_result );
    xaxisgraph[iterations] = iterations;
    yaxisgraph[iterations] = chi2;
  }

  // Make sure to set histogram axes.
  hist_original->GetXaxis()->SetTitle("E_{det} [GeV]");
  hist_original->GetYaxis()->SetTitle("#frac{1}{N}.#frac{dN}{dE}");

  // Plot the distributions and their ratios.
  TCanvas *can;
  PrepareCanvas( can, "ClosureTest_data" + variable);
  MakeDoublePaddedComparison(can, histos, leg );
  can->SaveAs("ClosureTest_data_" + variable + label_ + ".C");
  can->SaveAs("ClosureTest_data_" + variable + label_ + ".pdf");

  for( int i = 0; i < histos.size(); i++){

    cout << "histo\t" << i << "\t" << (histos[i])->GetNbinsX() << endl;

  }

  // Finish chi2 study.
  TCanvas *can_chi2;
  PrepareCanvas( can_chi2, "CHI2_Test_" + variable);
  TGraph* chi2_evolution = new TGraph(iterations_, xaxisgraph, yaxisgraph);
  chi2_evolution->GetXaxis()->SetTitle("N_{it.}");
  chi2_evolution->GetYaxis()->SetTitle("#chi^{2}");
  chi2_evolution->Draw("A*");
  can_chi2->SaveAs("Chi2_Test_data_" + variable + label_ + ".C");
  can_chi2->SaveAs("Chi2_Test_data_" + variable + label_ + ".pdf");
}




void Unfolder::ClosureTest_MC_detLevel(TString variable){

  // Prepare variables and objects.
  TH1D* hist_original;
  TH1D* hist_result;  
  TString htitle;
  TLegend *leg = new TLegend(0.65, 0.45, 0.95, 0.95);
  int iterations_ = 15;
  double xaxisgraph[iterations_], yaxisgraph[iterations_];

  vector<TH1D*> histos;
  
  // Get the detector level distribution.
  if( variable == "lead"){	Get_DetEnergy_lead(0, hist_original );	htitle = "E_{det} spectrum (leading)";}
  else{				Get_DetEnergy(0, hist_original );	htitle = "E_{det} spectrum (all)";	}

  if( normalise_1 ){ hist_original->Scale( 1./hist_original->Integral() ); }

  // Set detector level properties and legend.
  hist_original->SetLineColor( getColor( 1 ) );
  	hist_original->SetLineStyle( 1 );
  	hist_original->SetLineWidth( 3 );
  	hist_original->SetMarkerColor(  getColor( 1 ) );

  leg->AddEntry( hist_original, "Actual data", "p");
  histos.push_back( hist_original );

  // Loop over number of Bayesian iterations, and unfold->smear detector level distribution.
  // Calculate chi2 as we go.
  for(int iterations = 1; iterations <= iterations_; iterations++){

    hist_result = (TH1D*)hist_original->Clone(TString::Format("Closure_%i_iterations", iterations) );
    hist_result->SetName( TString::Format("Closure_%i_iterations", iterations) );
    hist_result->SetTitle( TString::Format("Closure_%i_iterations", iterations) );
    Get_DetUnfolded( 0, hist_result, iterations, variable );
    Get_GenSmeared( 0, hist_result , hist_result, variable); 
    //hist_int = (TH1D*)hist_result->Clone(TString::Format("it_unfolded_%i", iterations) );
    //list_of_iterations.push_back( hist_int);
    hist_result->SetLineColor( getColor( iterations + 1 ) );
    hist_result->SetLineStyle( iterations + 1 );
    hist_result->SetLineWidth( 3 );
    hist_result->SetMarkerColor(  getColor( iterations+1 ) );

    histos.push_back( hist_result );
    leg->AddEntry( hist_result, TString::Format("Treated data, %i it.", iterations) , "p");

    // Chi2 test.
    double chi2 = Chi2_test( hist_original, hist_result );
cout << "CHI2 is " << chi2 << endl;
    xaxisgraph[iterations] = iterations;
    yaxisgraph[iterations] = chi2;
  }

  // Make sure to set histogram axes.
  hist_original->GetXaxis()->SetTitle("E_{det} [GeV]");
  hist_original->GetYaxis()->SetTitle("#frac{1}{N}.#frac{dN}{dE}");

  // Plot the distributions and their ratios.
  TCanvas *can;
  PrepareCanvas( can, "ClosureTest_MC_detLevel" + variable);
  MakeDoublePaddedComparison(can, histos, leg );
  can->SaveAs("ClosureTest_MC_detLevel_" + variable + label_ + ".C");
  can->SaveAs("ClosureTest_MC_detLevel_" + variable + label_ + ".pdf");

  // Finish chi2 study.
  TCanvas *can_chi2;
  PrepareCanvas( can_chi2, "CHI2_Test_" + variable);
  TGraph* chi2_evolution = new TGraph(iterations_, xaxisgraph, yaxisgraph);
  chi2_evolution->GetXaxis()->SetTitle("N_{it.}");
  chi2_evolution->GetYaxis()->SetTitle("#chi^{2}");
  chi2_evolution->Draw("A*");
  can_chi2->SaveAs("Chi2_Test_MC_det_" + variable + label_ + ".C");
  can_chi2->SaveAs("Chi2_Test_MC_det_" + variable + label_ + ".pdf");
}




void Unfolder::ClosureTest_MC(TString variable){

  // Prepare variables and objects.
  TH1D* hist_original;
  TH1D* hist_result;  
  TString htitle;
  TLegend *leg = new TLegend(0.65, 0.45, 0.95, 0.95);
  int iterations_ = 15;
  double xaxisgraph[iterations_], yaxisgraph[iterations_];

  vector<TH1D*> histos;
  
  // Get the detector level distribution.
  if( variable == "lead"){	Get_GenEnergy_lead(0, hist_original );	htitle = "E_{det} spectrum (leading)";}
  else{				Get_GenEnergy(0, hist_original );	htitle = "E_{det} spectrum (all)";	}

  if( normalise_1 ){ hist_original->Scale( 1./hist_original->Integral() ); }

  // Set detector level properties and legend.
  hist_original->SetLineColor( getColor( 1 ) );
  	hist_original->SetLineStyle( 1 );
  	hist_original->SetLineWidth( 3 );
  	hist_original->SetMarkerColor(  getColor( 1 ) );

  leg->AddEntry( hist_original, "Generator level", "p");
  histos.push_back( hist_original );


  // Loop over number of Bayesian iterations, and unfold->smear detector level distribution.
  // Calculate chi2 as we go.
  for(int iterations = 1; iterations <= iterations_; iterations++){

    hist_result = (TH1D*)hist_original->Clone(TString::Format("Closure_%i_iterations", iterations) );
    hist_result->SetName( TString::Format("Closure_%i_iterations", iterations) );
    hist_result->SetTitle( TString::Format("Closure_%i_iterations", iterations) );
    Get_GenSmeared( 0, hist_result , hist_result, variable ); 
    Get_DetUnfolded( 0, hist_result, iterations, variable );

    hist_result->SetLineColor( getColor( iterations + 1 ) );
    hist_result->SetLineStyle( iterations + 1 );
    hist_result->SetLineWidth( 3 );
    hist_result->SetMarkerColor(  getColor( iterations+1 ) );

    histos.push_back( hist_result );
    leg->AddEntry( hist_result, TString::Format("Treated data, %i it.", iterations) , "p");

    // Chi2 test.
    double chi2 = Chi2_test( hist_original, hist_result );
    cout << "CHI2 is " << chi2 << endl;
    xaxisgraph[iterations] = iterations;
    yaxisgraph[iterations] = chi2;
  }

  // Make sure to set histogram axes.
  hist_original->GetXaxis()->SetTitle("E_{det} [GeV]");
  hist_original->GetYaxis()->SetTitle("#frac{1}{N}.#frac{dN}{dE}");

  // Plot the distributions and their ratios.
  TCanvas *can;
  PrepareCanvas( can, "ClosureTest_MC" + variable);
  MakeDoublePaddedComparison(can, histos, leg );
  can->SaveAs("ClosureTest_MC_" + variable + label_ + ".C");
  can->SaveAs("ClosureTest_MC_" + variable + label_ + ".pdf");

  // Finish chi2 study.
  TCanvas *can_chi2;
  PrepareCanvas( can_chi2, "CHI2_Test_" + variable);
  TGraph* chi2_evolution = new TGraph(iterations_, xaxisgraph, yaxisgraph);
  chi2_evolution->GetXaxis()->SetTitle("N_{it.}");
  chi2_evolution->GetYaxis()->SetTitle("#chi^{2}");
  chi2_evolution->Draw("A*");
  can_chi2->SaveAs("Chi2_Test_MC_" + variable + label_ + ".C");
  can_chi2->SaveAs("Chi2_Test_MC_" + variable + label_ + ".pdf");
}



/*************************************************************************************************
* -- The functions below split a canvas into two parts, for absolute and relative distributions. *
*************************************************************************************************/


void Unfolder::MakeDoublePaddedComparison(TCanvas * &can_, vector<TH1D*> histos, TLegend *legend){

  TPad* pad_abs_, *pad_ratio_;
  PrepareCanvas(can_, "Comparison_UnfoldedEnergy_Data_and_MC");

  SplitCanvas(can_, pad_abs_, pad_ratio_);

//  Hist_GenLevel(pad_abs_, true);
  Plot_Absolute(pad_abs_, histos);
  Plot_Ratio(pad_ratio_, histos);

  pad_abs_->cd();
  legend->Draw();
  pad_abs_->Update();
}


void Unfolder::Plot_Absolute( TPad* & pad_, vector<TH1D*> histos){
  cout << "----- Plot Absolute Value ---" << endl;

  double min_val, max_val;
  TH1D* hDistr, *hFirst;
  bool firstDistr = true, passed_data = false;
  int nDistr = 0;
  TString drawoptions = "phist";
  TString legendoptions = "l";

  for(int hist_ = 0; hist_ < histos.size(); hist_++){

    TH1D* hDistr = histos[hist_];

    Rebin_to( hDistr, 50);
    hDistr->GetYaxis()->SetTitleOffset( 0.85 );
    hDistr->GetYaxis()->SetTitleSize(0.09);
    hDistr->GetYaxis()->SetLabelSize(0.07);
    hDistr->GetXaxis()->SetNdivisions( 504 );
    hDistr->SetName( TString::Format("Absolute_%i", hist_) );

    if( firstDistr ){ 
      hFirst = (TH1D*)hDistr->Clone("First_absolute");
      hDistr = hFirst;
      firstDistr = false;
    }

    First_Plot( hFirst, hDistr, nDistr, max_val, min_val);
 
    pad_->cd();
    pad_->SetLogy();
    hDistr->Draw( drawoptions );
   
    drawoptions = "histsame";
    nDistr++;
    
    if( nDistr > 0 ) pad_->Update();
  }
}




void Unfolder::Plot_Ratio( TPad* & pad_, vector<TH1D*> histos){
  cout << "----- Plot Relative Value ---" << endl;

  double min_val, max_val;
  TH1D* hDistr, *hFirst, *hFirst_abs;
  bool firstDistr = true, passed_data = false;
  int nDistr = 0;
  TString drawoptions = "phist";

  for(int hist_ = 0; hist_ < histos.size(); hist_++){

    hDistr = (TH1D*) histos[hist_]->Clone( TString::Format("Ratio_%i", hist_) );

    Rebin_to( hDistr, 50);
    hDistr->GetXaxis()->SetNdivisions( 504 );

    if( firstDistr ){
      hFirst = (TH1D*)hDistr->Clone("First_ratio");
      hFirst_abs = (TH1D*)hFirst->Clone("_first_absolute");
      hFirst->Divide( hFirst );
      hFirst->GetYaxis()->SetTitle("Ratio");
      firstDistr = false;
    }

    hDistr->Divide( hFirst_abs );
    First_Plot( hFirst, hDistr, nDistr, max_val, min_val);

    pad_->cd();
    hDistr->Draw( drawoptions );
    hFirst->GetYaxis()->SetRangeUser(0., 3.);
 
    drawoptions = "histsame";
    nDistr++;

    if( nDistr > 0 ) pad_->Update();
    if( passed_data ) break;
  }
}

/**************************************************
* -- Prepare a label to be attached to all plots. *
**************************************************/


void Unfolder::LabelPlots( TString label){

  label_ = "_" + label;
}

/***********************
* -- CHI2 calculation. *
***********************/

double Unfolder::Chi2_test( TH1D* hist_data, TH1D* hist_MC){

  int ndf = hist_data->GetNbinsX();
  double chi2 = 0.;

  for( int bin = 0; bin <= ndf; bin++){

   double h1 = hist_data->GetBinContent( bin );
   double h2 = hist_MC->GetBinContent( bin );

   if( h1 != 0.){ 
     chi2 += (h1-h2)*(h1-h2)/h1;
    }
  }
  return chi2/ndf;
}

/*
*
*/

double Unfolder::Det_to_gen_scale(int file){
  TFile* file_ = TFile::Open(MC_files_[file], "read");
  TH1D* hDet = (TH1D*)file_->Get("hCastorJet_energy");
  TH1D* hGen = (TH1D*)file_->Get("hGenJet_energy");

  double scale_ratio = static_cast<double>( hDet->Integral() )/static_cast<double>( hGen->Integral() );
  return scale_ratio;
}	
