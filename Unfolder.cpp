//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>
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
#include "Function_average_histogram.h"
//#include "Function_make_Tex.h"


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
   void PrepareLegend( map<TString, TString> entries, map<TString, TString> printLabel );
   void PrepareTitles( map<TString, TString> xtitle ,  map<TString, TString> ytitle ,  map<TString, TString> htitle);
   void LabelPlots( TString label );
 
   void DoublePaddedComparison_unfolding(TString variable, int iterations = 0);
   void Plot_Unfolded(TPad* &pad_, TString variable, int iterations = 0);
   void Plot_Unfolded_Ratio(TPad* &pad_, TString variable, int iterations = 0);

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

   void Chi2_comparison_data(TString variable = "lead");
   void Chi2_comparison_MC(TString variable = "lead");
   void Chi2_comparison_MC_det(TString variable = "lead");
  
   double Chi2_testII( TH1D* hist_data, TH1D* hist_MC);
   double Chi2_testIII( TH1D* hist_data, TH1D* hist_MC);

   void CalibrationFunction(int phi_first, int phi_last, TGraphErrors* &gre);
   void CalibrationFunction_workingsectors(int first_sector, int last_sector, TGraphErrors* &gre);
   void CalibrationFunction_sectors();
   void CalculateSystematics(TString setup = "separate", int sector = 0);
   void CalculateSystematics_comparison();

   void Plot_Calibrated_functions();

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
   map<TString, TString> printLabel_;
   map<TString, vector<double> > calibration_parameters_;

   void MakeDoublePaddedComparison(TCanvas * &can, vector<TH1D*>, TLegend *leg);
   double Chi2_test( TH1D* hist_data, TH1D* hist_MC);
   double Det_to_gen_scale(int file);

   void Histogram_settings_ratio(TH1* hist);
   void Histogram_settings_absolute(TH1* hist);

   void Analyze_response(TH2D* hResponse_selection, TH2D* hGenE_selection, TH2D* hDetE_selection, TH2D* hGenE, TString label, TH1D* hGen, int file, TGraphErrors* &gre_meas, TGraphErrors* &gre_true );

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

  can_->SetLeftMargin(0.25);
  can_->SetTopMargin(0.07);
  can_->SetBottomMargin(0.14);
  can_->SetRightMargin(0.05);


}

void Unfolder::PrepareLegend( map<TString, TString> legend_info, map<TString, TString> printLabel ){

  legend_info_ = legend_info;
  printLabel_ = printLabel;

}

void Unfolder::PrepareTitles( map<TString, TString> xtitle, map<TString, TString> ytitle, map<TString, TString> htitle ){

  xtitle_ = xtitle;
  ytitle_ = ytitle;
  htitle_ = htitle;

}


/***********************************************
* -- Compare unfolded plots on a split canvas. *
***********************************************/



void Unfolder::DoublePaddedComparison_unfolding(TString variable, int iterations){

  TCanvas* can_;
  TPad* pad_abs_, *pad_ratio_;
  PrepareCanvas(can_, "Comparison_UnfoldedEnergy_Data_and_MC");

  SplitCanvas(can_, pad_abs_, pad_ratio_);

  Hist_GenLevel(pad_abs_, true);
  Plot_Unfolded(pad_abs_, variable, iterations);
  Plot_Unfolded_Ratio(pad_ratio_, variable, iterations);

  can_->SaveAs( folder_ + "/UnfoldedEnergy_Data_and_MC_" + variable + "_" + label_ + ".pdf");
  can_->SaveAs( folder_ + "/UnfoldedEnergy_Data_and_MC_" + variable + "_" + label_ + ".C");
}

//
// -- Absolute distributions.
//


void Unfolder::Plot_Unfolded(TPad* & pad_, TString variable, int iterations){
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

  int iterations_min = 1;
  int iterations_max = 15;

  if( iterations > 0 ){
    iterations_min = iterations;
    iterations_max = iterations;

  }
  for(int iterations_ = iterations_min; iterations_ < iterations_max+1; iterations_++){
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

void Unfolder::Plot_Unfolded_Ratio(TPad* & pad_, TString variable, int iterations){
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

  int iterations_min = 1;
  int iterations_max = 15;

  if( iterations > 0 ){
    iterations_min = iterations;
    iterations_max = iterations;

  }
  for(int iterations_ = iterations_min; iterations_ < iterations_max+1; iterations_++){

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
  int iterations_ = 16;
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
  double chi2 = hist_original->Chi2Test( hist_original, "CHI2/NDF");
  double chi2_prev = chi2;

  for(int iterations = 1; iterations <= iterations_; iterations++){

    hist_result = (TH1D*)hist_original->Clone(TString::Format("Closure_%i_iterations", iterations) );
    hist_result->SetName( TString::Format("Closure_%i_iterations", iterations) );
    hist_result->SetTitle( TString::Format("Closure_%i_iterations", iterations) );
    Get_DetUnfolded( -1, hist_result, iterations, variable );
    Get_GenSmeared( -1, hist_result , hist_result, variable); 

    hist_result->SetLineColor( getColor( iterations + 1) );
    hist_result->SetLineStyle( iterations + 1 );
    hist_result->SetLineWidth( 3 );
    hist_result->SetMarkerColor(  getColor( iterations + 1) );

    leg->AddEntry( hist_result, TString::Format("Treated data, %i it.", iterations) , "p");
    histos.push_back( hist_result );

    // Chi2 test.
    // chi2 = hist_original->Chi2Test( hist_result, "CHI2/NDF");
    chi2 = Chi2_test( hist_original, hist_result);
    xaxisgraph[iterations] = iterations;
    yaxisgraph[iterations] = chi2_prev - chi2;
    chi2_prev = chi2;
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

  // Finish chi2 study.
  TCanvas *can_chi2;
  PrepareCanvas( can_chi2, "CHI2_Test_" + variable);
  TGraph* chi2_evolution = new TGraph(iterations_, xaxisgraph, yaxisgraph);
  chi2_evolution->GetXaxis()->SetTitle("N_{it.}");
  chi2_evolution->GetYaxis()->SetTitle("#chi^{2}/NDF");
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
//    double chi2 = Chi2_test( hist_original, hist_result );
    double chi2 = hist_original->Chi2Test( hist_result, "CHI2/NDF");
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
  chi2_evolution->GetYaxis()->SetTitle("#chi^{2}/NDF");
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
//    double chi2 = Chi2_test( hist_original, hist_result );
    double chi2 = hist_original->Chi2Test( hist_result, "CHI2/NDF");
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
  chi2_evolution->GetYaxis()->SetTitle("#chi^{2}/NDF");
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
   double H1 = hist_data->Integral( );

   if( h1 != 0.){ 
     chi2 += (h1-h2)*(h1-h2)/h2;
    }
  }
  return chi2/ndf;
}

/**************************
* -- CHI2 calculation II. *
**************************/

double Unfolder::Chi2_testII( TH1D* hist_data, TH1D* hist_MC){

  int ndf = hist_data->GetNbinsX();
  double chi2 = 0.;

  for( int bin = 0; bin <= ndf; bin++){

   double h1 = hist_data->GetBinContent( bin );
   double h2 = hist_MC->GetBinContent( bin );
   double H1 = hist_data->Integral( );

   if( h1 != 0.){ 
     chi2 += (h1-h2)*(h1-h2)/(h1 + h2);
    }
  }
  return chi2/ndf;
}


/***************************
* -- CHI2 calculation III. *
***************************/

double Unfolder::Chi2_testIII( TH1D* hist_data, TH1D* hist_MC){

  int ndf = hist_data->GetNbinsX();
  double chi2 = 0.;

  for( int bin = 0; bin <= ndf; bin++){

   double h1 = hist_data->GetBinContent( bin );
   double h2 = hist_MC->GetBinContent( bin );
   double H1 = hist_data->Integral( );

   if( h1 != 0.){ 
     chi2 += (h1-h2)*(h1-h2)/H1;
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







/****************************************
* Compare different chi2 evolutions.
*/







void Unfolder::Chi2_comparison_data(TString variable){
  // Prepare variables and objects.
  TH1D* hist_original;
  TH1D* hist_result;  
  TString htitle;

  int iterations_ = 20;
  double xaxisgraph[iterations_], yaxisgraph_chi0[iterations_], yaxisgraph_chi1[iterations_], yaxisgraph_chi2[iterations_];
  
  // Get the detector level distribution.
  if( variable == "lead"){	Get_DetEnergy_lead(-1, hist_original );	htitle = "E_{det} spectrum (leading)";}
  else{				Get_DetEnergy(-1, hist_original );	htitle = "E_{det} spectrum (all)";	}

  if( normalise_1 ){ hist_original->Scale( 1./hist_original->Integral() ); }

  // Loop over number of Bayesian iterations, and unfold->smear detector level distribution.
  // Calculate chi2 as we go.
  double chi2 = hist_original->Chi2Test( hist_original, "CHI2/NDF");
  double chi2_prev = chi2;

  for(int iterations = 1; iterations <= iterations_; iterations++){

    hist_result = (TH1D*)hist_original->Clone(TString::Format("Closure_%i_iterations", iterations) );
    hist_result->SetName( TString::Format("Closure_%i_iterations", iterations) );
    hist_result->SetTitle( TString::Format("Closure_%i_iterations", iterations) );
    Get_DetUnfolded( -1, hist_result, iterations, variable );
    Get_GenSmeared( -1, hist_result , hist_result, variable); 


    // Chi2 test.
//    chi2 = hist_original->Chi2Test( hist_result, "CHI2/NDF");
    xaxisgraph[iterations] 	= iterations;
    yaxisgraph_chi0[iterations] = Chi2_test( hist_original, hist_result);
    yaxisgraph_chi1[iterations] = Chi2_testII( hist_original, hist_result);
    yaxisgraph_chi2[iterations] = Chi2_testIII( hist_original, hist_result);
  }


  // Finish chi2 study.
  TCanvas *can_chi2;
  PrepareCanvas( can_chi2, "CHI2_Test_" + variable);

  TGraph* chi2_evolution_0 = new TGraph(iterations_, xaxisgraph, yaxisgraph_chi0);
  chi2_evolution_0->GetXaxis()->SetTitle("N_{it.}");
  chi2_evolution_0->GetYaxis()->SetTitle("#chi^{2}/NDF");
  chi2_evolution_0->Draw("A*");
  /*
  TGraph* chi2_evolution_1 = new TGraph(iterations_, xaxisgraph, yaxisgraph_chi1);
  chi2_evolution_1->SetMarkerColor( kRed );
  chi2_evolution_1->SetMarkerStyle( 21 );
  chi2_evolution_1->Draw("psame");

  TGraph* chi2_evolution_2 = new TGraph(iterations_, xaxisgraph, yaxisgraph_chi2);
  chi2_evolution_2->SetMarkerColor( kGreen );
  chi2_evolution_2->SetMarkerStyle( 22 );
  chi2_evolution_2->Draw("psame");
  */

  can_chi2->SaveAs("Chi2_comparison_data_" + variable + label_ + ".C");
  can_chi2->SaveAs("Chi2_comparison_data_" + variable + label_ + ".pdf");
}










void Unfolder::Chi2_comparison_MC_det(TString variable){
  // Prepare variables and objects.
  TH1D* hist_original;
  TH1D* hist_result;  
  TString htitle;

  int iterations_ = 20;
  double xaxisgraph[iterations_], yaxisgraph_chi0[iterations_], yaxisgraph_chi1[iterations_], yaxisgraph_chi2[iterations_];
  
  // Get the detector level distribution.
  if( variable == "lead"){	Get_DetEnergy_lead(0, hist_original );	htitle = "E_{det} spectrum (leading)";}
  else{				Get_DetEnergy(0, hist_original );	htitle = "E_{det} spectrum (all)";	}

  if( normalise_1 ){ hist_original->Scale( 1./hist_original->Integral() ); }

  // Loop over number of Bayesian iterations, and unfold->smear detector level distribution.
  // Calculate chi2 as we go.
  double chi2 = hist_original->Chi2Test( hist_original, "CHI2/NDF");
  double chi2_prev = chi2;

  for(int iterations = 1; iterations <= iterations_; iterations++){

    hist_result = (TH1D*)hist_original->Clone(TString::Format("Closure_%i_iterations", iterations) );
    hist_result->SetName( TString::Format("Closure_%i_iterations", iterations) );
    hist_result->SetTitle( TString::Format("Closure_%i_iterations", iterations) );
    Get_DetUnfolded( 0, hist_result, iterations, variable );
    Get_GenSmeared( 0, hist_result , hist_result, variable); 


    // Chi2 test.
//    chi2 = hist_original->Chi2Test( hist_result, "CHI2/NDF");
    xaxisgraph[iterations] 	= iterations;
    yaxisgraph_chi0[iterations] = Chi2_test( hist_original, hist_result);
    yaxisgraph_chi1[iterations] = Chi2_testII( hist_original, hist_result);
    yaxisgraph_chi2[iterations] = Chi2_testIII( hist_original, hist_result);
  }


  // Finish chi2 study.
  TCanvas *can_chi2;
  PrepareCanvas( can_chi2, "CHI2_Test_" + variable);

  TGraph* chi2_evolution_0 = new TGraph(iterations_, xaxisgraph, yaxisgraph_chi0);
  chi2_evolution_0->GetXaxis()->SetTitle("N_{it.}");
  chi2_evolution_0->GetYaxis()->SetTitle("#chi^{2}/NDF");
  chi2_evolution_0->Draw("A*");

  /*
  TGraph* chi2_evolution_1 = new TGraph(iterations_, xaxisgraph, yaxisgraph_chi1);
  chi2_evolution_1->SetMarkerColor( kRed );
  chi2_evolution_1->SetMarkerStyle( 21 );
  chi2_evolution_1->Draw("psame");

  TGraph* chi2_evolution_2 = new TGraph(iterations_, xaxisgraph, yaxisgraph_chi2);
  chi2_evolution_2->SetMarkerColor( kGreen );
  chi2_evolution_2->SetMarkerStyle( 22 );
  chi2_evolution_2->Draw("psame");
  */

  can_chi2->SaveAs("Chi2_comparison_MC_det_" + variable + label_ + ".C");
  can_chi2->SaveAs("Chi2_comparison_MC_det_" + variable + label_ + ".pdf");
}
















void Unfolder::Chi2_comparison_MC(TString variable){
  // Prepare variables and objects.
  TH1D* hist_original;
  TH1D* hist_result;  
  TString htitle;

  int iterations_ = 20;
  double xaxisgraph[iterations_], yaxisgraph_chi0[iterations_], yaxisgraph_chi1[iterations_], yaxisgraph_chi2[iterations_];
  
  // Get the detector level distribution.
  if( variable == "lead"){	Get_GenEnergy_lead(0, hist_original );	htitle = "E_{det} spectrum (leading)";}
  else{				Get_GenEnergy(0, hist_original );	htitle = "E_{det} spectrum (all)";	}

  if( normalise_1 ){ hist_original->Scale( 1./hist_original->Integral() ); }

  // Loop over number of Bayesian iterations, and unfold->smear detector level distribution.
  // Calculate chi2 as we go.
  double chi2 = hist_original->Chi2Test( hist_original, "CHI2/NDF");
  double chi2_prev = chi2;

  for(int iterations = 1; iterations <= iterations_; iterations++){

    hist_result = (TH1D*)hist_original->Clone(TString::Format("Closure_%i_iterations", iterations) );
    hist_result->SetName( TString::Format("Closure_%i_iterations", iterations) );
    hist_result->SetTitle( TString::Format("Closure_%i_iterations", iterations) );
    Get_GenSmeared( 0, hist_result , hist_result, variable); 
    Get_DetUnfolded( 0, hist_result, iterations, variable );



    // Chi2 test.
//    chi2 = hist_original->Chi2Test( hist_result, "CHI2/NDF");
    xaxisgraph[iterations] 	= iterations;
    yaxisgraph_chi0[iterations] = Chi2_test( hist_original, hist_result);
    yaxisgraph_chi1[iterations] = Chi2_testII( hist_original, hist_result);
    yaxisgraph_chi2[iterations] = Chi2_testIII( hist_original, hist_result);
  }


  // Finish chi2 study.
  TCanvas *can_chi2;
  PrepareCanvas( can_chi2, "CHI2_Test_" + variable);
  TGraph* chi2_evolution_0 = new TGraph(iterations_, xaxisgraph, yaxisgraph_chi0);
  TGraph* chi2_evolution_1 = new TGraph(iterations_, xaxisgraph, yaxisgraph_chi1);
  TGraph* chi2_evolution_2 = new TGraph(iterations_, xaxisgraph, yaxisgraph_chi2);
  chi2_evolution_0->GetXaxis()->SetTitle("N_{it.}");
  chi2_evolution_0->GetYaxis()->SetTitle("#chi^{2}/NDF");
  chi2_evolution_0->Draw("A*");
  chi2_evolution_1->SetMarkerColor( kRed );
  chi2_evolution_1->SetMarkerStyle( 21 );
  chi2_evolution_1->Draw("psame");
  chi2_evolution_2->SetMarkerColor( kGreen );
  chi2_evolution_2->SetMarkerStyle( 22 );
  chi2_evolution_2->Draw("psame");
  can_chi2->SaveAs("Chi2_comparison_MC_" + variable + label_ + ".C");
  can_chi2->SaveAs("Chi2_comparison_MC_" + variable + label_ + ".pdf");
}








/************************************************************
* Code to extract a calibration function from the MC files. *
************************************************************/

void Unfolder::CalibrationFunction(int first_sector, int last_sector, TGraphErrors* &gre){

  ofstream calibrating_values;
  calibrating_values.open("Calibrating_values.txt", ios::out | ios::app | ios::binary);

  int rebinner = 1;

  // Step 1 - open Root files and retrieve 2D histograms.
   
  TLegend *legend_det = new TLegend(0.25, 0.50, 0.50, 0.95);
    legend_det->SetFillStyle( 0 );
    legend_det->SetBorderSize( 0 );
   
  TLegend *legend_gen = new TLegend(0.25, 0.50, 0.50, 0.95);
    legend_gen->SetFillStyle( 0 );
    legend_gen->SetBorderSize( 0 );   
 
  TString date = "20150305_had";

  TString numb = "12137851";

  int Slice_threshold = 0.;
  int event_threshold_ = 50.;


  /*******************************************
  * Create a new subdirectory for the plots. *
  ********************************************/
  double E_gen_cut = 0.;
  TString distr = "gaus";

  std::vector<TString> legendEntries;
  std::vector<TString> plotVariables;
  std::vector<TString> jetSelection;
    std::vector<TString> plotX;
    std::vector<TString> plotY;
    std::vector<TString> plotTitle;
  std::vector<int>     colours;
  std::vector<double> plot_min;
  std::vector<TString> can_suffix;
   
  bool draw_legend = false;
  TString calib_ = "_compare";
  TString drawoptions = "A";
  int color_index = 0;

  TCanvas *can_true = new TCanvas("can_true" + calib_, "can_true" + calib_, 1000, 1000);
    can_true	->SetLeftMargin(0.25);
    can_true	->SetTopMargin(0.05);
    can_true	->SetBottomMargin(0.14);	 
	 
  TCanvas *can_meas = new TCanvas("can_meas" + calib_, "can_meas" + calib_, 1000, 1000);
    can_meas	->SetLeftMargin(0.25);
    can_meas	->SetTopMargin(0.05);
    can_meas	->SetBottomMargin(0.14); 
      
  TCanvas *can_fit = new TCanvas("can_fit" + calib_, "can_fit" + calib_, 1000, 1000);
    can_fit	->SetLeftMargin(0.25);
    can_fit	->SetTopMargin(0.05);
    can_fit	->SetBottomMargin(0.14);         
	 
  TCanvas *can_aver = new TCanvas("can_aver" + calib_, "can_aver" + calib_, 1000, 1000);
    can_aver	->SetLeftMargin(0.25);
    can_aver	->SetRightMargin(0.05);
    can_aver	->SetBottomMargin(0.14);
      
  TLine *line = new TLine(0.,1.,1733.,1.);
  line->SetLineWidth(2);

//  TF1* analytical = new TF1("Analytical", " ( x < 200.) * ( [0] + [1] * log( [2] + x) ) + (x > 200. )* ( [3] + [4] * x) ", 0., 900.);
   TF1* analytical = new TF1("Analytical", " ( x > 100.) * ( [0] + [1] * log( [2] + x) )", 200., 900.);
    analytical->SetParLimits(2, -199., 1000.);

  /* Open files. */
  for(int file = 0; file < MC_files_.size(); file++){

    TFile *_file0 = TFile::Open( MC_files_[file], "Read");
    TString label = TString::Format( printLabel_[ MC_files_[file] ] + "_sector_%i_to_%i", first_sector, last_sector);


    label = "Test_this_" + label;
    int	new_dir = mkdir( ("Plots/" + label).Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  );
       
    TCanvas *c = new TCanvas("can", "can",  1000, 1000);
      c->SetLeftMargin(0.20);
      c->SetRightMargin(0.18);
      c->SetBottomMargin(0.20);  

    // -- Extract the response, Edet and Egen distributions from the files.

    THnSparseD* hResponse 	= (THnSparseD*)_file0->Get("hResponse_gen_phi");			//hResponse->Draw("colz");
    THnSparse* hGenE_fine_phi   = (THnSparse*)_file0->Get("hGen_fine_phi");
    THnSparse* hDetE_fine_phi = (THnSparse*)_file0->Get("hDet_fine_phi");
      
    TH2D* hDetE 	= (TH2D*)_file0->Get("hCastorJet_energy_response");//	hDetE->Draw("colz");	
      hDetE->RebinX( rebinner );
      hDetE->RebinY( rebinner );

    TH2D* hGenE	= (TH2D*)_file0->Get("hGen_fine");//			hGenE->Draw("colz");	
      hGenE->RebinX( rebinner );*
      hGenE->RebinY( rebinner );
       
    TH1D* hGen	= (TH1D*)_file0->Get("hGenJet_energy");		

    // Templates for 1D projections of THnSparse.
    TH1D* hSlice_storage;
    TH1D* hEgen_storage;
    TH1D* hEdet_storage;

    double mean_alpha_nom = 0., mean_alpha_denom = 0.;
    double mean_beta_nom = 0., mean_beta_denom = 0.;
    double mean_gamma_nom = 0., mean_gamma_denom = 0.;

    for( int bin_phi = first_sector; bin_phi <=last_sector; bin_phi++){//  hResponse->GetNbinsZ(); bin_phi++){
  //     if( bin_phi > 11 && bin_phi < 16 ){ continue; }       



       vector<double> mean_response;
       vector<double> mean_response_inverse;
       vector<double> mean_eDet;
       vector<double> error_eDet;
       vector<double> mean_eGen;
       vector<double> error_eGen;
       vector<double> eGen_center;
       vector<double> error_energy;
      
       int valid_fits = 0;

       THnSparse* hGenE_sparse = 	(THnSparse*)hGenE_fine_phi->Clone("");//TString::Format("GenE_bin_%i", bin_phi) );
       THnSparse* hDetE_sparse = 	(THnSparse*)hDetE_fine_phi->Clone("");//TString::Format("DetE_bin_%i", bin_phi) );
       THnSparse* hResponse_sparse =	(THnSparse*)hResponse	->Clone("");//TString::Format("Response_bin_%i", bin_phi) );

       hGenE_sparse ->GetAxis(2)->SetRange(bin_phi, bin_phi);
       hDetE_sparse ->GetAxis(2)->SetRange(bin_phi, bin_phi);
       hResponse_sparse      ->GetAxis(2)->SetRange(bin_phi, bin_phi);

       TH2D* hGenE_selection  = (TH2D*)hGenE_sparse->Projection(1,0);
         hGenE_selection->RebinX( rebinner );
	 hGenE_selection->RebinY( rebinner );
       TH2D* hDetE_selection  = (TH2D*)hDetE_sparse->Projection(1,0);
         hDetE_selection->RebinX( rebinner );
	 hDetE_selection->RebinY( rebinner );
       TH2D* hResponse_selection = (TH2D*)hResponse_sparse->Projection(1,0);
	 hResponse_selection->RebinX( rebinner );
	 hResponse_selection->RebinY( rebinner );
	 
       TGraphErrors * Response_true;
       TGraphErrors * Response_meas;

       Analyze_response( hResponse_selection, hGenE_selection, hDetE_selection, hGenE, label, hGen, file, Response_meas, Response_true);



       /*		  	
       // -- E gen versus average Egen.
       TGraph * Gen_average = new TGraph( eGen_center.size(), eGen_bin, E_gen_axis);
         Gen_average	->GetXaxis()->SetTitle("E_{gen}");
	 Gen_average	->GetYaxis()->SetTitle("<E>");
	 
       // -- E gen versus average Edet.
       TGraph * Det_average = new TGraph( eGen_center.size(), eGen_bin, E_det_axis);
         Det_average	->GetXaxis()->SetTitle("E_{gen}");
	 Det_average	->GetYaxis()->SetTitle("<E>");	
       */ 

       legend_det->AddEntry( Response_meas, TString::Format("Sector %i", bin_phi), "lp");
//       legend_gen->AddEntry( Response_meas, legendEntries[file], "lp");
 
       // Draw.
       can_true->cd();	 
       Response_true->Draw("p" + drawoptions);
       line->Draw();

	can_true->SaveAs(TString::Format("Plots/" + label + "/CalibrationFactors_true_calib_phi_%i_" + label + ".C", bin_phi) );
        can_true->SaveAs(TString::Format("Plots/" + label + "/CalibrationFactors_true_calib_phi_%i_.pdf", bin_phi) );

	// Draw.

       can_fit->cd();
       Response_meas->Draw("p" + drawoptions);
       drawoptions="same";	

       can_meas->cd(); 
       Response_meas->Draw("ape");// + drawoptions);
       line->Draw();

       cout << "Sector\t" << bin_phi << endl;
       analytical->SetLineColor( getColor(bin_phi) );
       Response_meas->Fit( analytical, "0", "", 200., 900.); 
       if( first_sector >= 12 && last_sector <= 15){ can_fit->cd(); analytical->DrawCopy("lsame"); }

       can_fit->SaveAs(TString::Format("Plots/" + label + "/CalibrationFactors_fit_calib_phi_%i_" + label + ".C", bin_phi) );
       can_fit->SaveAs(TString::Format("Plots/" + label + "/CalibrationFactors_fit_calib_phi_%i_" + label + ".pdf", bin_phi) );
         
       if( bin_phi < 12 || bin_phi > 15){
         double alpha = analytical->GetParameter( 0);
         double salpha = analytical->GetParError( 0);
         mean_alpha_nom += alpha/(salpha*salpha);
         mean_alpha_denom += 1./(salpha*salpha);       

         double beta = analytical->GetParameter( 1);
         double sbeta = analytical->GetParError( 1);
         mean_beta_nom += beta/(sbeta*sbeta);
         mean_beta_denom += 1./(sbeta*sbeta);   

         double gamma = analytical->GetParameter( 2);
         double sgamma = analytical->GetParError( 2);
         mean_gamma_nom += gamma/(sgamma*sgamma);
         mean_gamma_denom += 1./(sgamma*sgamma);   
       }

	can_meas->SaveAs(TString::Format("Plots/" + label + "/CalibrationFactors_meas_calib_phi_%i_" + label + ".C", bin_phi) );
        can_meas->SaveAs(TString::Format("Plots/" + label + "/CalibrationFactors_meas_calib_phi_%i_" + label + ".pdf", bin_phi) );


	// Delete the histograms to avoid memory leak.
       hGenE_selection->~TH2();
       hDetE_selection->~TH2();
       hResponse_selection->~TH2();

       hGenE_sparse->~THnSparse();
       hDetE_sparse->~THnSparse();
       hResponse_sparse->~THnSparse();
     } // Loop over phi.

     calibrating_values << legend_info_[ MC_files_[file] ] << endl;
     calibrating_values << " --- The Parameters --- " << endl;
     calibrating_values << "Mean alpha\t" << mean_alpha_nom/mean_alpha_denom << "\t+/-\t" << sqrt( 1./ mean_alpha_denom) << endl;
     calibrating_values << "Mean beta\t" << mean_beta_nom/mean_beta_denom << "\t+/-\t" << sqrt( 1./ mean_beta_denom) << endl;
     calibrating_values << "Mean gamma\t" << mean_gamma_nom/mean_gamma_denom << "\t+/-\t" << sqrt( 1./ mean_gamma_denom) << endl << endl;

     std::vector<double> parameters_calibration;
       parameters_calibration.push_back( mean_alpha_nom/mean_alpha_denom );
       parameters_calibration.push_back( mean_beta_nom/mean_beta_denom );
       parameters_calibration.push_back( mean_gamma_nom/mean_gamma_denom );

     calibration_parameters_[ MC_files_[file] ] = parameters_calibration;

     // Save.
     can_true->cd();
     legend_gen->Draw();
     	can_true->SaveAs("Plots/" + label +"/CalibrationFactors_true_" + calib_ + "_" + label + ".C");		
	can_true->SaveAs("Plots/" + label +"/CalibrationFactors_true_" + calib_ + "_" + label + ".pdf");
     
     can_meas->cd();
     legend_det->Draw();
     	can_meas->SaveAs("Plots/" + label +"/CalibrationFactors_meas_" + calib_ + "_" + label + ".C"); 		
	can_meas->SaveAs("Plots/" + label +"/CalibrationFactors_meas_" + calib_ + "_" + label + ".pdf");

     can_fit->cd();
     legend_det->Draw();
//        analytical->SetParameters(mean_alpha_nom/mean_alpha_denom, mean_beta_nom/mean_beta_denom, mean_gamma_nom/mean_gamma_denom);
	TF1* analytical_mean = new TF1("Analytical_mean", "( [0] + [1] * log( [2] + x) )", 200., 900.);
        analytical_mean->SetParameters(mean_alpha_nom/mean_alpha_denom, mean_beta_nom/mean_beta_denom, mean_gamma_nom/mean_gamma_denom); 
	analytical_mean->SetLineColor( TColor::GetColor("#FC00E7") );

	can_fit->SaveAs("Plots/" + label +"/CalibrationFactors_meas_with_fit" + calib_ + "_" + label + ".C");             
	can_fit->SaveAs("Plots/" + label +"/CalibrationFactors_meas_with_fit" + calib_ + "_" + label + ".pdf");

    drawoptions = "A";

    can_fit->Clear();
    can_meas->Clear();
    can_true->Clear();
    legend_det->Clear();

  } // Loop over files.
}




/************************************************************
* Code to extract a calibration function from the MC files. *
************************************************************/

void Unfolder::CalibrationFunction_workingsectors(int first_sector, int last_sector, TGraphErrors* &gre){

  ofstream calibrating_values;
  calibrating_values.open("Calibrating_values.txt", ios::out | ios::app | ios::binary);

  int rebinner = 1;

  // Step 1 - open Root files and retrieve 2D histograms.
   
  TLegend *legend_det = new TLegend(0.25, 0.60, 0.50, 0.95);
    legend_det->SetFillStyle( 0 );
    legend_det->SetBorderSize( 0 );
   
  TLegend *legend_gen = new TLegend(0.25, 0.60, 0.50, 0.95);
    legend_gen->SetFillStyle( 0 );
    legend_gen->SetBorderSize( 0 );  

  int Slice_threshold = 0.;
  int event_threshold_ = 50.;


  /*******************************************
  * Create a new subdirectory for the plots. *
  ********************************************/
  double E_gen_cut = 0.;
  TString distr = "gaus";

  std::vector<TString> legendEntries;
  std::vector<TString> plotVariables;
  std::vector<TString> jetSelection;
    std::vector<TString> plotX;
    std::vector<TString> plotY;
    std::vector<TString> plotTitle;
  std::vector<int>     colours;
  std::vector<double> plot_min;
  std::vector<TString> can_suffix;
   
  bool draw_legend = false;
  TString calib_ = "_compare";
  TString drawoptions = "A";
  int color_index = 0;

  TCanvas *can_true = new TCanvas("can_true" + calib_, "can_true" + calib_, 1000, 1000);
    can_true	->SetLeftMargin(0.25);
    can_true	->SetTopMargin(0.05);
    can_true	->SetBottomMargin(0.14);	 
	 
  TCanvas *can_meas = new TCanvas("can_meas" + calib_, "can_meas" + calib_, 1000, 1000);
    can_meas	->SetLeftMargin(0.25);
    can_meas	->SetTopMargin(0.05);
    can_meas	->SetBottomMargin(0.14); 
      
  TCanvas *can_fit = new TCanvas("can_fit" + calib_, "can_fit" + calib_, 1000, 1000);
    can_fit	->SetLeftMargin(0.25);
    can_fit	->SetTopMargin(0.05);
    can_fit	->SetBottomMargin(0.14);         
	 
  TCanvas *can_aver = new TCanvas("can_aver" + calib_, "can_aver" + calib_, 1000, 1000);
    can_aver	->SetLeftMargin(0.25);
    can_aver	->SetRightMargin(0.05);
    can_aver	->SetBottomMargin(0.14);
      
  TLine *line = new TLine(0.,1.,1733.,1.);
  line->SetLineWidth(2);

//  TF1* analytical = new TF1("Analytical", " ( x < 200.) * ( [0] + [1] * log( [2] + x) ) + (x > 200. )* ( [3] + [4] * x) ", 0., 900.);
  TF1* analytical = new TF1("Analytical", " ( [0] + [1] * log( [2] + x) )", 200., 900.);
    analytical->SetParLimits(2, -199., 1000.);

  /* Open files. */

  for(int file = 0; file < MC_files_.size(); file++){

    TFile *_file0 = TFile::Open( MC_files_[file], "Read");
    TString label;
    if( first_sector == 1 && last_sector == 16){ label = TString::Format( printLabel_[ MC_files_[file] ] + "_allsectors"); }
    else if( first_sector != last_sector ){ label = TString::Format( printLabel_[ MC_files_[file] ] + "_sectors_%i_to_%i", first_sector, last_sector); }
    else{ label = TString::Format( printLabel_[ MC_files_[file] ] + "_sector_%i", first_sector); }

    label = "O7O7_" + label;
    int	new_dir = mkdir( ("Plots/" + label).Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  );
       
    TCanvas *c = new TCanvas("can", "can",  1000, 1000);
      c->SetLeftMargin(0.20);
      c->SetRightMargin(0.18);
      c->SetBottomMargin(0.20);  

    // -- Extract the response, Edet and Egen distributions from the files.

    THnSparseD* hResponse 	= (THnSparseD*)_file0->Get("hResponse_gen_phi");			//hResponse->Draw("colz");
    THnSparse* hGenE_fine_phi   = (THnSparse*)_file0->Get("hGen_fine_phi");
    THnSparse* hDetE_fine_phi = (THnSparse*)_file0->Get("hDet_fine_phi");
      
    TH2D* hDetE 	= (TH2D*)_file0->Get("hCastorJet_energy_response");//	hDetE->Draw("colz");	
      hDetE->RebinX( rebinner );
      hDetE->RebinY( rebinner );

    TH2D* hGenE	= (TH2D*)_file0->Get("hGen_fine");//			hGenE->Draw("colz");	
      hGenE->RebinX( rebinner );*
      hGenE->RebinY( rebinner );
       
    TH1D* hGen	= (TH1D*)_file0->Get("hGenJet_energy");		

    // Templates for 1D projections of THnSparse.
    TH1D* hSlice_storage;
    TH1D* hEgen_storage;
    TH1D* hEdet_storage;

    // Prepare the histograms;
    // 1. Project the 3D histogram into its desired axes.
    // 2. Rebin the histogram to its desired number of bins per axis.
    // 3. Reset the contents etc.
    TH2D* hGenE_selection  = (TH2D*)hGenE_fine_phi->Projection(1,0);
       hGenE_selection->RebinX( rebinner );
       hGenE_selection->RebinY( rebinner );
       hGenE_selection->Reset();
    TH2D* hDetE_selection  = (TH2D*)hDetE_fine_phi->Projection(1,0);
       hDetE_selection->RebinX( rebinner );
       hDetE_selection->RebinY( rebinner );
       hDetE_selection->Reset();
    TH2D* hResponse_selection = (TH2D*)hResponse->Projection(1,0);
       hResponse_selection->RebinX( rebinner );
       hResponse_selection->RebinY( rebinner );
       hResponse_selection->Reset();

    for( int bin_phi = first_sector; bin_phi <= last_sector; bin_phi++){//  hResponse->GetNbinsZ(); bin_phi++){
       if( first_sector != last_sector && ( bin_phi > 11 && bin_phi < 16) ){ continue; }

       THnSparse* hGenE_sparse = 	(THnSparse*)hGenE_fine_phi->Clone("");//TString::Format("GenE_bin_%i", bin_phi) );
       THnSparse* hDetE_sparse = 	(THnSparse*)hDetE_fine_phi->Clone("");//TString::Format("DetE_bin_%i", bin_phi) );
       THnSparse* hResponse_sparse =	(THnSparse*)hResponse	->Clone("");//TString::Format("Response_bin_%i", bin_phi) );

       hGenE_sparse 	->GetAxis(2)->SetRange(bin_phi, bin_phi);
       hDetE_sparse 	->GetAxis(2)->SetRange(bin_phi, bin_phi);
       hResponse_sparse ->GetAxis(2)->SetRange(bin_phi, bin_phi);

       TH2D* hGenE_selection_phi  = (TH2D*)hGenE_sparse->Projection(1,0);
         hGenE_selection_phi	->RebinX( rebinner );
	 hGenE_selection_phi	->RebinY( rebinner );
       TH2D* hDetE_selection_phi  = (TH2D*)hDetE_sparse->Projection(1,0);
         hDetE_selection_phi	->RebinX( rebinner );
	 hDetE_selection_phi	->RebinY( rebinner );
       TH2D* hResponse_selection_phi = (TH2D*)hResponse_sparse->Projection(1,0);
	 hResponse_selection_phi->RebinX( rebinner );
	 hResponse_selection_phi->RebinY( rebinner );

       hGenE_selection		->Add( hGenE_selection_phi );
       hDetE_selection		->Add( hDetE_selection_phi );
       hResponse_selection	->Add( hResponse_selection_phi );

       // Delete the histograms to avoid memory leak.
       hGenE_selection_phi->~TH2();
       hDetE_selection_phi->~TH2();
       hResponse_selection_phi->~TH2();

       hGenE_sparse->~THnSparse();
       hDetE_sparse->~THnSparse();
       hResponse_sparse->~THnSparse();

/*	 
cout << "// -- BINS COMPARISON\t" << hGenE_selection->GetNbinsX() << "\t" << hGenE_selection->GetXaxis()->GetTitle() << "\t" << hGenE_selection->GetNbinsY() << "\t" << hGenE_selection->GetYaxis()->GetTitle() << "\t"
				  << hDetE_selection->GetNbinsX() << "\t" << hDetE_selection->GetXaxis()->GetTitle() << "\t" << hDetE_selection->GetNbinsY() << "\t" << hDetE_selection->GetYaxis()->GetTitle() << "\t"
				  << hResponse_selection->GetNbinsX() << "\t" << hResponse_selection->GetXaxis()->GetTitle() << "\t" << hResponse_selection->GetNbinsY() << "\t" << hResponse_selection->GetYaxis()->GetTitle() << endl;	
*/

     } // Loop over phi.


     TGraphErrors * Response_true;
     TGraphErrors * Response_meas;

     Analyze_response( hResponse_selection, hGenE_selection, hDetE_selection, hGenE, label, hGen, file, Response_meas, Response_true);
		  	


       // -- E gen versus average Egen.
       //TGraph * Gen_average = new TGraph( eGen_center.size(), eGen_bin, E_gen_axis);
       //  Gen_average	->GetXaxis()->SetTitle("E_{gen}");
       //  Gen_average	->GetYaxis()->SetTitle("<E>");
	 
       // -- E gen versus average Edet.
       //TGraph * Det_average = new TGraph( eGen_center.size(), eGen_bin, E_det_axis);
       //  Det_average	->GetXaxis()->SetTitle("E_{gen}");
       //  Det_average	->GetYaxis()->SetTitle("<E>");	
 
     legend_det->AddEntry( Response_meas, TString::Format("Sectors 1-11, 16"), "lp");
       //  legend_gen->AddEntry( Response_meas, legendEntries[file], "lp");
 
     // Draw.
     can_true->cd();	 
     Response_true->Draw("p" + drawoptions);
     line->Draw();

     can_true->SaveAs(TString::Format("Plots/" + label + "/CalibrationFactors_true_calib_" + label + ".C") );
     can_true->SaveAs(TString::Format("Plots/" + label + "/CalibrationFactors_true_calib_" + label + ".pdf") );

     // Draw.

     can_fit->cd();
     Response_meas->Draw("p" + drawoptions);
     
     can_fit->SaveAs(TString::Format("Plots/" + label + "/CalibrationFactors_fit_calib_" + label + ".C") );
     can_fit->SaveAs(TString::Format("Plots/" + label + "/CalibrationFactors_fit_calib_" + label + ".pdf") );
	
     drawoptions="same";	

     can_meas->cd(); 
     Response_meas->Draw("ape");// + drawoptions);
     line->Draw();

     analytical->SetLineColor( getColor(file+1) );
     Response_meas->Fit( analytical, "", "", 200., 900.);

     can_meas->SaveAs(TString::Format("Plots/" + label + "/CalibrationFactors_meas_calib_" + label + ".C") );
     can_meas->SaveAs(TString::Format("Plots/" + label + "/CalibrationFactors_meas_calib_" + label + ".pdf") );

     std::vector<double> parameters_calibration;
       parameters_calibration.push_back( analytical->GetParameter(0) );
       parameters_calibration.push_back( analytical->GetParameter(1) );
       parameters_calibration.push_back( analytical->GetParameter(2) );

     calibration_parameters_[ MC_files_[file] ] = parameters_calibration;

     // Save.
     can_true->cd();
     legend_gen->Draw();
     	can_true->SaveAs("Plots/" + label +"/CalibrationFactors_true_" + calib_ + "_" + label + ".C");		
	can_true->SaveAs("Plots/" + label +"/CalibrationFactors_true_" + calib_ + "_" + label + ".pdf");
     
     can_meas->cd();
     legend_det->Draw();
     	can_meas->SaveAs("Plots/" + label +"/CalibrationFactors_meas_" + calib_ + "_" + label + ".C"); 		
	can_meas->SaveAs("Plots/" + label +"/CalibrationFactors_meas_" + calib_ + "_" + label + ".pdf");

     can_fit->cd();
     legend_det->Draw();

	can_fit->SaveAs("Plots/" + label +"/CalibrationFactors_meas_with_fit" + calib_ + "_" + label + ".C");             
	can_fit->SaveAs("Plots/" + label +"/CalibrationFactors_meas_with_fit" + calib_ + "_" + label + ".pdf");

    drawoptions = "A";

    can_fit->Clear();
    can_meas->Clear();
    can_true->Clear();
    legend_det->Clear();

     cout << "Finished file\t" << file << endl;

  } // Loop over files.
  cout << "Finished all files" << endl;
}





void Unfolder::Plot_Calibrated_functions(){

  TGraphErrors *graph_sector;
  TString drawoptions = "ape";

  ofstream parameters;
  parameters.open("Parameters.txt");

  // Create a canvas to plot.
  TCanvas *can;
  PrepareCanvas( can, "Systematics_comparison");

  TLegend *leg = new TLegend(0.25, 0.7, 0.5, 0.9);
    leg->SetFillStyle( 0 );
    leg->SetBorderSize( 0 );
  
for(int file_ = 0; file_ < MC_files_.size(); file_++){

  for(int sector = 11; sector <= 15; sector++){
    if( sector == 11 ){ CalibrationFunction_workingsectors(1,16, graph_sector); } 
    else{ CalibrationFunction_workingsectors(sector, sector, graph_sector); } 

    can->cd();

    if( sector != 11){ 
      graph_sector->SetMarkerColor( getColor( sector ) ); 
      graph_sector->SetLineColor( getColor( sector ));
    }

    graph_sector->Draw(drawoptions);
    drawoptions = "psame";

    TF1* analytical = new TF1("Analytical", " ( [0] + [1] * log( [2] + x) )", 200., 900.);
    analytical->SetParLimits(2, -199., 1000.);   

    analytical->SetLineColor( graph_sector->GetLineColor() );
    graph_sector->Fit( analytical, "0", "", 200., 900.);  

    bool line = false;  
/*
    if( analytical->GetParameter(2) <= -200. ){ 
	analytical = new TF1("Analytical", " ( [0] + [1] *  x )", 200., 900.); 
	graph_sector->Fit( analytical, "0", "", 200., 900.);
	line = true;
    }
*/
    analytical->Draw("lsame"); 

    TString legend_entry;

    if( sector != 11){ legend_entry = TString::Format("Sector %i", sector); }
    else{ legend_entry = TString::Format("Sectors 1-11,16");}

    if( file_ == 0){ leg->AddEntry(graph_sector, legend_entry, "lp"); }

    parameters << legend_entry 		<< "	& $"
	 << analytical->GetParameter(0) << "$	& $\\pm "
	 << analytical->GetParError(0) 	<< "$	& $"
	 << analytical->GetParameter(1) << "$	& $\\pm "
	 << analytical->GetParError(1) 	<< "$	& $";
    parameters << analytical->GetParameter(2) << "$	& $\\pm " << analytical->GetParError(2) << "$	\\\\ \n"; 
//    else{       	parameters << "$&\\\\\n"; }
  }

}
  parameters << "\\hline ";
  parameters.close();

  leg->Draw();
  can->SaveAs("DoesThisWork.C");
  can->SaveAs("DoesThisWork.pdf");
  
}






/************************************************************
* Code to extract a calibration function from the MC files. *
************************************************************/

void Unfolder::CalibrationFunction_sectors(){

  ofstream calibrating_values;
  calibrating_values.open("Calibrating_values.txt", ios::out | ios::app | ios::binary);

  int rebinner = 1;

  // Step 1 - open Root files and retrieve 2D histograms.
   
  TLegend *legend_det = new TLegend(0.25, 0.60, 0.75, 0.95);
    legend_det->SetFillStyle( 0 );
    legend_det->SetBorderSize( 0 );
   
  TLegend *legend_gen = new TLegend(0.25, 0.60, 0.50, 0.95);
    legend_gen->SetFillStyle( 0 );
    legend_gen->SetBorderSize( 0 );  

  int Slice_threshold = 0.;
  int event_threshold_ = 50.;

  /*******************************************
  * Create a new subdirectory for the plots. *
  ********************************************/
  double E_gen_cut = 0.;
  TString distr = "gaus";

  std::vector<TString> legendEntries;
  std::vector<TString> plotVariables;
  std::vector<TString> jetSelection;
    std::vector<TString> plotX;
    std::vector<TString> plotY;
    std::vector<TString> plotTitle;
  std::vector<int>     colours;
  std::vector<double> plot_min;
  std::vector<TString> can_suffix;
   
  bool draw_legend = false;
  TString calib_ = "_compare";
  TString drawoptions = "A";
  int color_index = 0;

  TCanvas *can_true = new TCanvas("can_true" + calib_, "can_true" + calib_, 1000, 1000);
    can_true	->SetLeftMargin(0.25);
    can_true	->SetTopMargin(0.05);
    can_true	->SetBottomMargin(0.14);	 
	 
  TCanvas *can_meas = new TCanvas("can_meas" + calib_, "can_meas" + calib_, 1000, 1000);
    can_meas	->SetLeftMargin(0.25);
    can_meas	->SetTopMargin(0.05);
    can_meas	->SetBottomMargin(0.14); 
      
  TCanvas *can_fit = new TCanvas("can_fit" + calib_, "can_fit" + calib_, 1000, 1000);
    can_fit	->SetLeftMargin(0.25);
    can_fit	->SetTopMargin(0.05);
    can_fit	->SetBottomMargin(0.14);         
	 
  TCanvas *can_aver = new TCanvas("can_aver" + calib_, "can_aver" + calib_, 1000, 1000);
    can_aver	->SetLeftMargin(0.25);
    can_aver	->SetRightMargin(0.05);
    can_aver	->SetBottomMargin(0.14);
      
  TLine *line = new TLine(0.,1.,1733.,1.);
  line->SetLineWidth(2);

//  TF1* analytical = new TF1("Analytical", " ( x < 200.) * ( [0] + [1] * log( [2] + x) ) + (x > 200. )* ( [3] + [4] * x) ", 0., 900.);
  TF1* analytical = new TF1("Analytical", " ( [0] + [1] * log( [2] + x) )", 200., 900.);
    analytical->SetParLimits(2, -199., 1000.);

  /* Open files. */


  for(int sector = 11; sector <= 15; sector++){
    TString label;
    if( sector == 11 )	{ label = TString::Format( "Systematics_comparison_allsectors"); }
    else 		{ label = TString::Format( "Systematics_comparison_sector_%i", sector); }

    int	new_dir = mkdir( ("Plots/" + label).Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  );

    for(int file = 0; file < MC_files_.size(); file++){

      TFile *_file0 = TFile::Open( MC_files_[file], "Read");
       
      TCanvas *c = new TCanvas("can", "can",  1000, 1000);
        c->SetLeftMargin(0.20);
        c->SetRightMargin(0.18);
        c->SetBottomMargin(0.20);  

       // -- Extract the response, Edet and Egen distributions from the files.
  
      THnSparseD* hResponse 	= (THnSparseD*)_file0->Get("hResponse_gen_phi");			//hResponse->Draw("colz");
      THnSparse* hGenE_fine_phi   = (THnSparse*)_file0->Get("hGen_fine_phi");
      THnSparse* hDetE_fine_phi = (THnSparse*)_file0->Get("hDet_fine_phi");
      
      TH2D* hDetE 	= (TH2D*)_file0->Get("hCastorJet_energy_response");//	hDetE->Draw("colz");	
        hDetE->RebinX( rebinner );
        hDetE->RebinY( rebinner );

      TH2D* hGenE	= (TH2D*)_file0->Get("hGen_fine");//			hGenE->Draw("colz");	
        hGenE->RebinX( rebinner );*
        hGenE->RebinY( rebinner );
       
      TH1D* hGen	= (TH1D*)_file0->Get("hGenJet_energy");		

      // Templates for 1D projections of THnSparse.
      TH1D* hSlice_storage;
      TH1D* hEgen_storage;
      TH1D* hEdet_storage;

      // Prepare the histograms;
      // 1. Project the 3D histogram into its desired axes.
      // 2. Rebin the histogram to its desired number of bins per axis.
      // 3. Reset the contents etc.
      TH2D* hGenE_selection  = (TH2D*)hGenE_fine_phi->Projection(1,0);
         hGenE_selection->RebinX( rebinner );
         hGenE_selection->RebinY( rebinner );
         hGenE_selection->Reset();
      TH2D* hDetE_selection  = (TH2D*)hDetE_fine_phi->Projection(1,0);
         hDetE_selection->RebinX( rebinner );
         hDetE_selection->RebinY( rebinner );
         hDetE_selection->Reset();
      TH2D* hResponse_selection = (TH2D*)hResponse->Projection(1,0);
         hResponse_selection->RebinX( rebinner );
         hResponse_selection->RebinY( rebinner );
         hResponse_selection->Reset();



      for( int bin_phi = 1; bin_phi <= 16; bin_phi++){
       if( sector == 11 && (bin_phi > 11 && bin_phi < 16)){ continue; }
       else if( sector != 11 && sector != bin_phi ){ continue; }

       THnSparse* hGenE_sparse = 	(THnSparse*)hGenE_fine_phi->Clone("");	//TString::Format("GenE_bin_%i", bin_phi) );
       THnSparse* hDetE_sparse = 	(THnSparse*)hDetE_fine_phi->Clone("");	//TString::Format("DetE_bin_%i", bin_phi) );
       THnSparse* hResponse_sparse =	(THnSparse*)hResponse	->Clone("");	//TString::Format("Response_bin_%i", bin_phi) );

       hGenE_sparse 	->GetAxis(2)->SetRange(bin_phi, bin_phi);
       hDetE_sparse 	->GetAxis(2)->SetRange(bin_phi, bin_phi);
       hResponse_sparse ->GetAxis(2)->SetRange(bin_phi, bin_phi);

       TH2D* hGenE_selection_phi  = (TH2D*)hGenE_sparse->Projection(1,0);
         hGenE_selection_phi	->RebinX( rebinner );
	 hGenE_selection_phi	->RebinY( rebinner );
       TH2D* hDetE_selection_phi  = (TH2D*)hDetE_sparse->Projection(1,0);
         hDetE_selection_phi	->RebinX( rebinner );
	 hDetE_selection_phi	->RebinY( rebinner );
       TH2D* hResponse_selection_phi = (TH2D*)hResponse_sparse->Projection(1,0);
	 hResponse_selection_phi->RebinX( rebinner );
	 hResponse_selection_phi->RebinY( rebinner );

       hGenE_selection		->Add( hGenE_selection_phi );
       hDetE_selection		->Add( hDetE_selection_phi );
       hResponse_selection	->Add( hResponse_selection_phi );

       // Delete the histograms to avoid memory leak.
       hGenE_selection_phi->~TH2();
       hDetE_selection_phi->~TH2();
       hResponse_selection_phi->~TH2();

       hGenE_sparse->~THnSparse();
       hDetE_sparse->~THnSparse();
       hResponse_sparse->~THnSparse();
      } // Loop over phi.

      TGraphErrors * Response_true;
      TGraphErrors * Response_meas;

      Analyze_response( hResponse_selection, hGenE_selection, hDetE_selection, hGenE, label, hGen, file, Response_meas, Response_true);

       Response_meas->SetName( TString::Format("%s_" + printLabel_[ MC_files_[file] ], Response_meas->GetName()) );
       Response_true->SetName( TString::Format("%c_" + printLabel_[ MC_files_[file] ], Response_true->GetName()) );
		  	
      if( sector == 11){legend_det->AddEntry( Response_meas, TString::Format(legend_info_[ MC_files_[file] ] + " (Sectors 1-11, 16)"), "lp"); }
      else{		legend_det->AddEntry( Response_meas, TString::Format(legend_info_[ MC_files_[file] ] + " (Sector %i)", sector), "lp"); }
 
       // Draw.
       can_true->cd();	 
       Response_true->Draw("p" + drawoptions);
       line->Draw();

       can_true->SaveAs(TString::Format("Plots/" + label + "/CalibrationFactors_true_calib_" + label + ".C") );
       can_true->SaveAs(TString::Format("Plots/" + label + "/CalibrationFactors_true_calib_" + label + ".pdf") );

       // Draw.
 
       can_fit->cd();
       Response_meas->Draw("p" + drawoptions);
     
       can_fit->SaveAs(TString::Format("Plots/" + label + "/CalibrationFactors_fit_calib_" + label + ".C") );
       can_fit->SaveAs(TString::Format("Plots/" + label + "/CalibrationFactors_fit_calib_" + label + ".pdf") );
	
       drawoptions="same";	

       can_meas->cd(); 
       Response_meas->Draw("ape");// + drawoptions);
       line->Draw();

       analytical->SetLineColor( getColor(file+1) );
       Response_meas->Fit( analytical, "", "", 200., 900.);

       can_meas->SaveAs(TString::Format("Plots/" + label + "/CalibrationFactors_meas_calib_" + label + ".C") );
       can_meas->SaveAs(TString::Format("Plots/" + label + "/CalibrationFactors_meas_calib_" + label + ".pdf") );

       std::vector<double> parameters_calibration;
         parameters_calibration.push_back( analytical->GetParameter(0) );
         parameters_calibration.push_back( analytical->GetParameter(1) );
         parameters_calibration.push_back( analytical->GetParameter(2) );

       calibration_parameters_[ MC_files_[file] ] = parameters_calibration;

     }
     cout << "Finished all files" << endl;

       // Save.
       can_true->cd();
       legend_gen->Draw();
       	can_true->SaveAs("Plots/" + label +"/CalibrationFactors_true_" + calib_ + "_" + label + ".C");		
	can_true->SaveAs("Plots/" + label +"/CalibrationFactors_true_" + calib_ + "_" + label + ".pdf");
     
       can_meas->cd();
       legend_det->Draw();
     	can_meas->SaveAs("Plots/" + label +"/CalibrationFactors_meas_" + calib_ + "_" + label + ".C"); 		
	can_meas->SaveAs("Plots/" + label +"/CalibrationFactors_meas_" + calib_ + "_" + label + ".pdf");

       can_fit->cd();
       legend_det->Draw();

	can_fit->SaveAs("Plots/" + label +"/CalibrationFactors_meas_with_fit" + calib_ + "_" + label + ".C");             
	can_fit->SaveAs("Plots/" + label +"/CalibrationFactors_meas_with_fit" + calib_ + "_" + label + ".pdf");

      drawoptions = "A";

      can_fit->Clear();
      can_meas->Clear();
      can_true->Clear();
      legend_det->Clear();


  } // Loop over sectors sects.
}













void Unfolder::CalculateSystematics(TString setup, int sector){
  // Start by executing the calibration determination to obtain the necessary parameters.
  TGraphErrors* gre;
  if( setup == "all" ){ CalibrationFunction_workingsectors(1,16, gre); }
  else{ CalibrationFunction_workingsectors(sector, sector, gre); }

  // Create a canvas to plot.
  TCanvas *can;
  PrepareCanvas( can, "Systematics_comparison");

  TPad *pad1, *pad2;
  SplitCanvas(can, pad1, pad2);

  TLegend *legend_det = new TLegend(0.25, 0.0, 0.50, 0.50);
    legend_det->SetFillStyle( 0 );
    legend_det->SetBorderSize( 0 );

  // 1. Loop over files.
  // 2. Get parameters and use in function.
  // 3. Plot function as Graph.
  // 4. Profit.

  TString drawoptions = "";
  TH1F* original_calibration;

  // Needed for the range of the y-axis.
  double min_val, max_val;

  for(int file = 0; file < MC_files_.size(); file++){

    double alpha = (calibration_parameters_[ MC_files_[file] ])[0];
    double beta = (calibration_parameters_[ MC_files_[file] ])[1];
    double gamma = (calibration_parameters_[ MC_files_[file] ])[2];

    TF1* calibration_function = new TF1("calibration_function_", "[0] + [1] * log( [2] + x)", 200., 900.);
    calibration_function->SetParameters( alpha, beta, gamma );
    TH1* calibration_histogram = calibration_function->GetHistogram();

    if( file == 0 ){
      original_calibration = (TH1F*)calibration_histogram->Clone("Original_calibration_function_" + setup);
      original_calibration->GetYaxis()->SetRangeUser(0., 2.);
    }

    legend_det->AddEntry( calibration_histogram, legend_info_[ MC_files_[file]], "l");

    pad1->cd();

    calibration_histogram->SetLineColor( getColor( file+1 ) );
    calibration_histogram->GetXaxis()->SetTitle("E_{det}");
    calibration_histogram->GetYaxis()->SetRangeUser(0., 4.);
    calibration_histogram->GetYaxis()->SetTitle("<#frac{E_{gen}}{E_{det}}>");
    calibration_histogram->DrawCopy("hist" + drawoptions);

    pad1->Update();

    pad2->cd();
    calibration_histogram->Divide( original_calibration );
    calibration_histogram->GetYaxis()->SetTitle( "Syst./MC");
    calibration_histogram->GetYaxis()->SetRangeUser(0.75,1.2);

    calibration_histogram->Draw("hist" + drawoptions);

    drawoptions = "same"; 
  }
  pad1->cd();
  legend_det->Draw();

  can->SaveAs(TString::Format( "Systematics_comparison_" + setup + "_%i.C", sector));
  can->SaveAs(TString::Format( "Systematics_comparison_" + setup + "_%i.pdf", sector) );
}







void Unfolder::Histogram_settings_absolute(TH1* hist){

  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.07);
  hist->GetYaxis()->SetTitleOffset( 1.25);  

}

void Unfolder::Histogram_settings_ratio(TH1* hist){

  hist->GetYaxis()->SetLabelSize(0.065);
  hist->GetYaxis()->SetTitleSize(0.08);
  hist->GetYaxis()->SetTitleOffset( 1.00);  

  hist->GetXaxis()->SetLabelSize(0.06);
  hist->GetXaxis()->SetTitleSize(0.08);
  hist->GetXaxis()->SetTitleOffset( 0.8); 

}


void Unfolder::Analyze_response(TH2D* hResponse_selection, TH2D* hGenE_selection, TH2D* hDetE_selection, TH2D* hGenE, TString label, TH1D* hGen, int file, TGraphErrors* &gre_meas, TGraphErrors* &gre_true){

     int event_threshold_ = 50;
     int rebinner = 1;
     int valid_fits = 0;

     vector<double> mean_response;
     vector<double> mean_response_inverse;
     vector<double> mean_eDet;
     vector<double> error_eDet;
     vector<double> mean_eGen;
     vector<double> error_eGen;
     vector<double> eGen_center;
     vector<double> error_energy;

     TH1D* hSlice_storage;
     TH1D* hEgen_storage;
     TH1D* hEdet_storage;

     // Create a histogram to store slices with low statistics.
     hSlice_storage     = (TH1D*)hResponse_selection->ProjectionY("Storage", 1, 1, "do");
       for(int bins_store = 0; bins_store <= hSlice_storage->GetNbinsX(); bins_store++){ hSlice_storage->SetBinContent(bins_store, 0.); }

     // Do the same for the energies.
     hEgen_storage = (TH1D*)hGenE_selection->ProjectionY("Storage_gen", 0, 0, "do"); hEgen_storage->Rebin( rebinner*100 );
       for(int bins_store = 0; bins_store <= hEgen_storage->GetNbinsX(); bins_store++){ hEgen_storage->SetBinContent(bins_store, 0.); }

     hEdet_storage = (TH1D*)hDetE_selection->ProjectionY("Storage_det", 0, 0, "do"); hEdet_storage->Rebin( rebinner*100 );
       for(int bins_store = 0; bins_store <= hEdet_storage->GetNbinsX(); bins_store++){ hEdet_storage->SetBinContent(bins_store, 0.); }

     /******************
     * Loop over bins. *       
     ******************/                 
                    
     for(int bin_E = 1; bin_E < hGenE->GetNbinsX() ; bin_E++){

       // Get 1D projections.
       TH1D* hResponse_1D 	= (TH1D*)hResponse_selection->ProjectionY(TString::Format("Response_1D_bin_%i", bin_E), 	bin_E, bin_E, "do");

       TH1D* hGenE_1D  = (TH1D*)hGenE_selection->ProjectionY(TString::Format("eGen_1D_bin_%i", bin_E), bin_E, bin_E, "o");
	 hGenE_1D->Rebin( rebinner*100 );

       TH1D* hDetE_1D  = (TH1D*)hDetE_selection->ProjectionY(TString::Format("eDet_1D_bin_%i", bin_E), bin_E, bin_E, "o");
	 hDetE_1D->Rebin( rebinner*100 );
 
       // Create gaussian.
       TF1 *fit_response	= new TF1(TString::Format("Fit_Response_1D_bin_%i", bin_E), "gaus");

       // If our current histogram has too few events, we save it for later.
       if( (hResponse_1D->Integral() + hSlice_storage->Integral() )< event_threshold_ ){
 	 hSlice_storage	->Add( hResponse_1D );	
	 hEgen_storage	->Add( hGenE_1D );	
	 hEdet_storage	->Add( hDetE_1D );
	 continue;
       }

       // Our current histogram, combined with the stored histogram, has enough statistics. Add them and go on.
       else if( ( hResponse_1D->Integral() + hSlice_storage->Integral() ) >= event_threshold_ ){
         hResponse_1D ->Add( hSlice_storage );
         hGenE_1D     ->Add( hEgen_storage );
	 hDetE_1D	->Add( hEdet_storage );
       }

       // Empty the bins of our storage histogram for a next run.
       for(int bins_store = 0; bins_store <= hSlice_storage->GetNbinsX(); bins_store++){ hSlice_storage->SetBinContent(bins_store, 0.); 	hSlice_storage->SetBinError(bins_store, 0.);}
       for(int bins_store = 0; bins_store <= hEgen_storage->GetNbinsX(); bins_store++){ hEgen_storage->SetBinContent(bins_store, 0.); 	hEgen_storage->SetBinError(bins_store, 0.);}
       for(int bins_store = 0; bins_store <= hEdet_storage->GetNbinsX(); bins_store++){ hEdet_storage->SetBinContent(bins_store, 0.); 	hEdet_storage->SetBinError(bins_store, 0.);}

	 TCanvas* can_slice = new TCanvas( TString::Format("Slice_%i", bin_E), TString::Format("Slice_%i", bin_E), 1);

	 // -- Fit Gaussian.
	 hResponse_1D	->Draw();
	 
	 if( file == 0 ){
  	   hResponse_1D	->Fit( fit_response , "Q0");	 
	   can_slice->SaveAs(TString::Format("Plots/" + label + "/Slice_%i.C", bin_E) );
           can_slice->SaveAs(TString::Format("Plots/" + label + "/Slice_%i.pdf", bin_E) );
	}
	else{
	  hResponse_1D ->Fit( fit_response , "Q0");
	}
	 
	 // -- Extract mean values.
	 if( fit_response->GetParameter(1) != 0. ){
	   double R = fit_response->GetParameter(1);
	   double sR= fit_response->GetParError(1);	

	   mean_response	.push_back( R );
	   mean_response_inverse.push_back( 1./R );
	   mean_eDet		.push_back( GetAverage(hDetE_1D) );
	   mean_eGen		.push_back( GetAverage(hGenE_1D) );
	   eGen_center		.push_back( hGen->GetBinCenter( bin_E ) );

	   if( R == 0.){ 	error_eGen	.push_back( 0. ); }
	   else{ 		error_eGen	.push_back( sR ); }
           if( R == 0.){	error_eDet	.push_back( 0. ); }
	   else{ 		error_eDet	.push_back( sR/(R*R) ); }
	   error_energy		.push_back( 0. );	   
	   
	   valid_fits++;
	   
	 }
       } // Loop over energybins.
       
       // Print out our values.
       //for( int bin = 0; bin < mean_eDet.size(); bin++){
         //cout << "lowedge.push_back(\t" << mean_eGen[bin] << "\t); muval.push_back(\t" << mean_response[bin] << ");" << endl;
       //}

       // Transform the vectors into arrays.       
     double * E_gen_axis	= &mean_eGen[0],
       	      * E_det_axis	= &mean_eDet[0],
	      * response_val	= &mean_response[0],
	      * response_inverse= &mean_response_inverse[0],
	      * eGen_bin	= &eGen_center[0],
	      * err_R		= &error_eGen[0],
	      * err_1overR	= &error_eDet[0],
	      * err_energy	= &error_energy[0];
       // Create graphs from the arrays.

       // -- E gen versus response
     TGraphErrors * Response_true = new TGraphErrors( mean_eGen.size(), E_gen_axis, response_val, err_energy, err_R);
     Response_true	->GetXaxis()->SetTitle("E_{gen}");
     Response_true	->GetXaxis()->SetTitleOffset(0.9);
     Response_true	->GetXaxis()->SetTitleSize(0.06);
     Response_true	->GetXaxis()->SetLabelSize(0.05);
     Response_true	->GetXaxis()->SetNdivisions(505);
	 	 
     Response_true	->GetYaxis()->SetTitle("< #frac{E_{det}}{E_{gen}} >");
     Response_true	->GetYaxis()->SetTitleOffset(1.8);
     Response_true	->GetYaxis()->SetTitleSize(0.06);
     Response_true	->GetYaxis()->SetLabelSize(0.07);	 
     Response_true	->GetYaxis()->SetRangeUser(0., 4.);
     
     Response_true	->SetMarkerColor( getColor(1) );
     Response_true	->SetLineColor( getColor(1) );
     Response_true	->SetMarkerStyle( 23 );

     Response_true	->SetName( TString::Format("Truth_response") );

     gre_true = (TGraphErrors*)Response_true->Clone( "TGraphErrors_Calibration_response_true_" + label );
	 
       // -- E det versus 1/response
     TH1F* hEdet_axis = new TH1F("hEdet_axis", "hEdet_axis", 100, 0, mean_eGen[ mean_eGen.size()-1] );
     hEdet_axis->GetYaxis()->SetRangeUser(0., 4.0);
     hEdet_axis	->GetXaxis()->SetTitle("E_{det}");
     hEdet_axis	->GetXaxis()->SetTitleOffset(0.9);
     hEdet_axis	->GetXaxis()->SetTitleSize(0.06);
     hEdet_axis	->GetXaxis()->SetLabelSize(0.05);
     hEdet_axis	->GetXaxis()->SetNdivisions(505);
	 	 
     hEdet_axis	->GetYaxis()->SetTitle("< #frac{E_{gen}}{E_{det}} >");  
     hEdet_axis	->GetYaxis()->SetRangeUser(0.,4.0);      
     hEdet_axis	->GetYaxis()->SetTitleOffset(1.8);
     hEdet_axis	->GetYaxis()->SetTitleSize(0.06);
     hEdet_axis	->GetYaxis()->SetLabelSize(0.07); 
	       
     TGraphErrors * Response_meas = new TGraphErrors( mean_eDet.size(), E_det_axis, response_inverse,  err_energy, err_1overR);
     Response_meas	->GetXaxis()->SetTitle("E_{det}");
     Response_meas	->GetXaxis()->SetTitleOffset(0.9);
     Response_meas	->GetXaxis()->SetTitleSize(0.06);
     Response_meas	->GetXaxis()->SetLabelSize(0.05);
     Response_meas	->GetXaxis()->SetNdivisions(505);
	 	 
     Response_meas	->GetYaxis()->SetTitle("< #frac{E_{gen}}{E_{det}} >");  
     Response_meas	->GetYaxis()->SetRangeUser(0.,4.0);      
     Response_meas	->GetYaxis()->SetTitleOffset(1.8);
     Response_meas	->GetYaxis()->SetTitleSize(0.06);
     Response_meas	->GetYaxis()->SetLabelSize(0.07);
	 
	 
     Response_meas	->SetMarkerColor( getColor(file+1) );
     Response_meas	->SetLineColor( getColor(file+1) );
     Response_meas	->SetMarkerStyle( 20 + file );

     gre_meas = (TGraphErrors*)Response_meas->Clone( "TGraphErrors_Calibration_response_meas_" + label );
}
