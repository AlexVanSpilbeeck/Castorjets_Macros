//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TColor.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <TTree.h>
#include <TText.h>
#include <TThread.h>
#include <TVectorD.h>


//STANDARD C++ INCLUDES
#include <sstream>
#include <iostream>
#include <iomanip>
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
#include <ctime>

#include "../../../RooUnfold-1.1.1/src/RooUnfold.h"
#include "../../../RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"
#include "../../../RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "../../../RooUnfold-1.1.1/src/RooUnfold.h"
#include "../../../RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "color.h"
#include "Function_Rebin.h"
#include "Function_average_histogram.h"
#include "GetSubHistogram.h"
//#include "Function_make_Tex.h"
#include "NonZeroMinimum.h"
	
//#define normalise_1 true

// Define the plotter class.
class Unfolder{
  public:
   Unfolder(vector<TString> MC_files, TString datafile, std::map< TString, std::map<TString, TString> >, double Eplotmin, double Ethresh, TString folder, int normalise); // list of MC files - datafile - map to store histograms.
   ~Unfolder();

   // Get histograms from files.

   void Chi2_comparison_data(TString variable = "lead");
   void Chi2_comparison_MC(TString variable = "lead");
   void Chi2_comparison_MC_det(TString variable = "lead");
   double Chi2_testII( TH1D* hist_data, TH1D* hist_MC);
   double Chi2_testIII( TH1D* hist_data, TH1D* hist_MC);
   void ClosureTest_data(TString variable = "lead", TString file = "Displaced", int method = 1);
   void ClosureTest_MC(TString variable = "lead");
   void ClosureTest_MC_detLevel(TString variable = "lead");
   void CompareGenLevel();
   void DoublePaddedComparison(TString variable);
   void DoublePaddedComparison_unfolding(TString variable, int iterations = 0);
   void Get_DetEnergy(int file_, TH1D* &hist_);
   void Get_DetEnergy_lead(int file_, TH1D* &hist_);
   void Get_DetUnfolded(int file_, int MC_, TH1D* &hist_, TMatrixD& covariance_m, int iterations = 4, TString variable = "all");
   void Get_DetUnfolded(int file_, TH1D* &hist_, int iterations = 4, TString variable = "all");
   void Get_Distribution(int file_, TH1D* &hist_, TString variable);
   void Get_GenEnergy(int file_, TH1D* &hist_);
   void Get_GenEnergy_lead(int file_, TH1D* &hist_); 
   void Get_GenSmeared(int file_, TH1D* &hist_, TH1D* hGen, TString variable = "all");
   void Get_GenSmeared(int file_, int MC_, TH1D* &hist_, TH1D* hGen, TString variable = "all");
   void Get_Misses(int file_, TH1D* &hist_);
   void Get_ResponseMatrix(int file_, TH2D* &hist_);
   void Get_ResponseMatrixTHn(int file_, THnSparseD* & hist_);
   void Hist_DetLevel();
   void Hist_DetLevel_DataWithSmear();
   void Hist_DetLevel_lead();
   void Hist_DetLevel_MCWithSmear();
   void Hist_GenLevel();   
   void Hist_GenLevel(TString variable);   
   void Hist_GenLevel(TPad* &pad_, bool isFirst);
   void Hist_getFake(int file_, TH1D* &hist_);
   void Hist_getMiss(int file_, TH1D* &hist_);
   void Hist_getTruth();
   void LabelPlots( TString label );
   void Plot_Absolute(TPad* &pad_, TString variable);
   void Plot_Absolute(TPad* &pad_, vector<TH1D*>);
   void Plot_Ratio(TPad* &pad_, TString variable);
   void Plot_Ratio(TPad* &pad_, vector<TH1D*>);
   void Plot_Unfolded();
   void Plot_Unfolded(TPad* &pad_, TString variable, int iterations = 0);
   void Plot_Unfolded_Ratio();
   void Plot_Unfolded_Ratio(TPad* &pad_, TString variable, int iterations = 0);
   void Prepare_2Dplot( TH2D* &hist_);
   void PrepareCanvas( TCanvas* &can_, TString label);
   void PrepareCanvas_2D( TCanvas* &can_, TString label);
   void PrepareLegend( map<TString, TString> entries, map<TString, TString> printLabel );
   void PrepareTitles( map<TString, TString> xtitle ,  map<TString, TString> ytitle ,  map<TString, TString> htitle);
   void SetAddLabel( TString label );
   void SetCastorJetEnergy_norm( double renorm );
   void SetSubhistogram_cut(double Ecut);
   void SetFitDraw( bool fitnotdrawn );
   void SplitCanvas(TCanvas* &can_, TPad* &pad_abs_, TPad* &pad_ratio_);
   void Systematics_CompareDetLevel();
   void Systematics_CompareGenLevel();

   void Unfolding_data(TString variable = "lead", int iterations_ = 10);


   void CalibrationFunction(int phi_first, int phi_last, TGraphErrors* &gre);
   void CalibrationFunction_workingsectors(int first_sector, int last_sector, TGraphErrors* &gre);
   void CalibrationFunction_sectors(int first_sector, int last_sector, int whichfile, TGraphErrors* &gre);
   void CalculateSystematics(TString setup = "separate", int first_sector = 1, int last_sector = 16);
   void CalculateSystematics_comparison();

   void Plot_Calibrated_functions();
   void PlotResponseMatrix(int file_ );

   void PlotFromResponseMatrix(int file_, TString axis_1, TString axis_2);
   void PlotFromResponseMatrix(int file_, TString axis_1);

   void Dissect_ResponseObject();

   void CalibrationFactors_oneCanvas(bool draw_functions);

   double Calculate_chi2(TH1D* hist_ref, TH1D* hist_res);

   void Get_UnfoldSmearError( TH1D* vUnfold, TH1D* &vSmeared, int file_, TString variable, int iterations );
   void Unfold_and_smear(int file_, TH1D* &hist_, int MC_, int iterations, TString variable, int method , double &chi2);

   void Calculate_smearedBackError(int file_, int MC_, int iterations, TH1D* &hActual_hist, TH1D* hUnfold);
   double Calculate_smearedBackError_covariance( TH1D* hData, TH1D* hUnfold, RooUnfoldResponse* response, int iterations );

   void CovarianceMatrix( TH1D* hUnfold, TH1D* hSmeared, TMatrixD& cov_);

  private:
   vector<TString> MC_files_;
   TString datafile_;
   TString folder_; 
   TString label_;
   TString setup_calibration_;
   TString addLabel_;

   bool normalise_1;
   double Eplotmin_;
   double Ethresh_; 
   double fitting_threshold_;

   TH2D* hError_hist;

   map<TString, TString> legend_info_;
   map<TString, TString> xtitle_;
   map<TString, TString> ytitle_;
   map<TString, TString> htitle_;
   map<TString, TString> printLabel_;
   map<TString, vector<double> > calibration_parameters_;

   map<TString, TGraphErrors*> calibration_graphs_;

   map<TString, map<TString, TString> > set_of_tags_;

   void MakeDoublePaddedComparison(TCanvas * &can, vector<TH1D*>, TLegend *leg);
   double Chi2_test( TH1D* hist_data, TH1D* hist_MC);
   double Det_to_gen_scale(int file);
   bool BadNumerals( TH1D* hist );


   double renorm_;

   double Ecut_;

   bool fitnotdrawn_;

   void Histogram_settings_ratio(TH1* hist);
   void Histogram_settings_absolute(TH1* hist);

   void Analyze_response(TH2D* hResponse_selection, TH2D* hGenE_selection, TH2D* hDetE_selection, TH2D* hGenE, TString label, TH1D* hGen, int file, TGraphErrors* &gre_meas, TGraphErrors* &gre_true );
   int Extract_2D_energy_distributions(int file_, TH2D* &hGenE_selection, TH2D* &hDetE_selection, TH2D* &hResponse_selection, TH2D* &hDetE, TH2D* &hGenE, TH1D* &hGen, TString setup);

   TGraphErrors* gre_;
   int SetAxisTHnSparse(TString axis);

   void FillAnew_1D(TH1D* hOld, TH1D* &hNew, TRandom3* &rand);
   void FillAnew_2D(TH2D* hOld, TH2D* &hNew, TRandom3* &rand);
   void FillAnew_1D(TH1D* hOld, TH1D* hMin, TH1D* &hNew, TRandom3* &rand);


};

Unfolder::Unfolder( vector<TString> MC_files, TString datafile, std::map< TString, std::map<TString, TString> > set_of_tags, double Eplotmin, double Ethresh, TString folder, int normalise )
{
  gStyle->SetStatStyle(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetOptStat(0);

  MC_files_ = MC_files;
  datafile_ = datafile;
  folder_ = "Plots/" + folder + "/";
  set_of_tags_ = set_of_tags;
  Eplotmin_ = Eplotmin;
  Ethresh_ = Ethresh;
  addLabel_ = "";  

  int  new_dir = mkdir( ("Plots/" + folder).Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  );

  normalise_1 = false;
  if (normalise == 1 ){ normalise_1 = true; }
  renorm_ = 1;

}
Unfolder::~Unfolder(){
}

void Unfolder::SetFitDraw(bool fitnotdrawn){
  fitnotdrawn_ = !fitnotdrawn;
}


/*****************************************
* Plot the Det level energy distribution *
*****************************************/
// -- Plot
void Unfolder::Get_DetEnergy(int file_, TH1D* &hist_){
   cout << "\n\n\n// -- Detector level energy - file selection\t" << file_ << endl;

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





void Unfolder::Get_ResponseMatrix(int file_, TH2D* &hist_){

   TFile *_file0;
   TString drawoptions = "";

   if( file_ < MC_files_.size() ){      _file0 = new TFile( MC_files_[file_], "Read");  drawoptions = "hist";}
   else if( file_ == -1 ){              _file0 = new TFile( datafile_, "Read");         drawoptions = "data";}

   TH2D* hRes = (TH2D*)_file0->Get("hCastorJet_energy_response_fine");      hRes->Sumw2();
   hRes->GetXaxis()->SetTitle("E_{det} (GeV)");
   hRes->GetYaxis()->SetTitle("E_{gen} (GeV)");

   if( normalise_1 ) hRes->Scale( 1./hRes->Integral() );

   hist_ = hRes;
}




void Unfolder::Get_ResponseMatrixTHn(int file_, THnSparseD* & hist_){

   cout << "---\tGet_ResponseMatrixTHn" << endl;
   TFile *_file0;

   if( file_ < MC_files_.size() ){      _file0 = new TFile( MC_files_[file_], "Read");  }
   else if( file_ == -1 ){              _file0 = new TFile( datafile_, "Read");         }

   cout << "---\tOpened file " << file_ << "\t" << MC_files_[file_] << endl;
   THnSparseD* hRes = (THnSparseD*)_file0->Get("hResponse_leading");      hRes->Sumw2();
   cout << "---\tGot THnSparse" << endl;

   hist_ = hRes;

   for( int nD = 0; nD < hist_->GetNdimensions(); nD++){

     cout << "\tAxis\t" << nD << "\t" << hist_->GetAxis( nD )->GetTitle() << endl;

   }
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



void Unfolder::Hist_getMiss(int file_, TH1D* &hMiss){
  TFile* _file0 = TFile::Open( MC_files_[file_], "Read");
  RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response");

//  hMiss= (TH1D*) response->Hmisses();
//  cout << "Misses\t" << hMiss->Integral() << endl;
}



void Unfolder::Hist_getFake(int file_, TH1D* &hFake){
  TFile* _file0 = TFile::Open( MC_files_[file_], "Read");
  RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response");

  hFake = (TH1D*) response->Hfakes();
  cout << "Fakes\t" << hFake->Integral() << endl;

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
   cout << "\n\n\n// -- Generator level energy - file selection\t" << file_ << endl;
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



void Unfolder::Get_Misses(int file_, TH1D* &hist_){
   TH1D *hGen;
   TFile *_file0;
   TString drawoptions = "";
   
   if( file_ < MC_files_.size() ){	_file0 = new TFile( MC_files_[file_], "Read"); 	drawoptions = "hist";}
   
   hGen = (TH1D*)_file0->Get("hCastorJet_miss_all");	hGen->Sumw2();
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

  TH1D* hist_, *hist_misses_;
  TString drawoptions = "hist";
  TCanvas *can; 
  PrepareCanvas( can, "Generator_Level_Energy" );
  can->SetLogy();
  TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.95);
  TString scalefactor_data = (set_of_tags_[ "scalefactors" ])[ datafile_ ];
  double scalefactors_data_ = scalefactor_data.Atof();

  for( int file = 0; file < MC_files_.size(); file++ ){
    int file_ = file;
    TString scalefactor_MC = (set_of_tags_[ "scalefactors" ])[ MC_files_[file_] ];
    double scalefactors_MC_ = scalefactor_MC.Atof() ;

    Get_GenEnergy( file_ , hist_ );
    Get_Misses( file_, hist_misses_ );
    hist_->Add( hist_misses_ );
    hist_->SetLineColor( getColor( file_ +1) );
    hist_->SetLineWidth( 3 );
    hist_->SetMarkerColor( getColor( file_ +1) );
    hist_->SetMarkerSize( 1 );
    hist_->GetXaxis()->SetNdivisions(504);

    hist_->Scale( scalefactors_data_/scalefactors_MC_ );
    
    hist_->Draw(drawoptions);
    drawoptions = "histsame";  

    hist_misses_->Scale( scalefactors_data_/scalefactors_MC_ );
    hist_misses_->Draw( drawoptions );

    leg->AddEntry( hist_, legend_info_[ MC_files_[file_] ], "l");
  }

  TFile* _file0 = new TFile( MC_files_[0], "Read");
  TH1D* hMeas;	Get_DetEnergy(-1, hMeas );
  RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response_all");
  RooUnfoldBayes unfold_bayes(response, hMeas, 60);
  TH1D* hUnfold = (TH1D*) unfold_bayes.Hreco( RooUnfold::kCovariance ); 
  hUnfold->SetMarkerStyle( 22 );
  hUnfold->SetMarkerColor( kRed );
  hUnfold->Draw("edatasame");

  leg->AddEntry( hUnfold, TString::Format("Data") , "p");

  leg->Draw();  

  can->SaveAs(folder_ + "/GeneratorLevelEnergy" + label_ + "_lead.C");
  can->SaveAs(folder_ + "/GeneratorLevelEnergy" + label_ + "_lead.pdf");
  
}







void Unfolder::Hist_GenLevel(TString variable){

  //-- Preamble for detector level distributions.
  TFile* _file0 = new TFile( MC_files_[0], "Read");
  RooUnfoldResponse* response;

  TH1D* hist_, *hist_misses_, *hMeas;
  TString drawoptions = "hist";
  TCanvas *can; 
  PrepareCanvas( can, "Generator_Level_Energy" );
  can->SetLogy();
  TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.95);

  if( variable == "all" ){ 	response = (RooUnfoldResponse*)_file0->Get("response"); 	
    Get_DetEnergy(-1, hMeas ); }
  else if( variable == "lead" ){response = (RooUnfoldResponse*)_file0->Get("response_lead"); 	
    Get_DetEnergy_lead(-1, hMeas ); }
  else{ return; }

  TH1D* hMeasured_response = (TH1D*) response->Hmeasured();
  double scalefactors_MC_ = hMeasured_response->Integral();
  double scalefactors_data_ = hMeas->Integral();

  //-- Get generator level distribution.
  for( int file = 0; file < MC_files_.size(); file++ ){
    //-- Preamble.
    int file_ = file;  

    //-- Generator level jet energy distribution.
    hist_ = (TH1D*) response->Htruth();



    hist_->SetLineColor( getColor( file_ +1) );
    hist_->SetLineWidth( 3 );
    hist_->SetMarkerColor( getColor( file_ +1) );
    hist_->SetMarkerSize( 1 );
    hist_->GetXaxis()->SetNdivisions(504);
cout << "Hist_ before\t" << hist_->Integral() << endl;	
    hist_->Scale( scalefactors_data_/scalefactors_MC_ );
cout << "Unfolder - " << scalefactors_data_/scalefactors_MC_ << endl;
cout << "Hist_ after\t" << hist_->Integral() << endl;
    hist_->DrawCopy(drawoptions);
    drawoptions = "histsame";  

    leg->AddEntry( hist_, "MC (Gen), with misses", "l");
  }

  //-- Unfold data.
//  hMeas->Scale( scalefactors_MC_/scalefactors_data_ );

  RooUnfoldBayes unfold_bayes(response, hMeas, 40);
  unfold_bayes.SetVerbose(-1);
  TH1D* hUnfold = (TH1D*) unfold_bayes.Hreco( RooUnfold::kCovariance ); 
  hUnfold->SetMarkerStyle( 22 );
  hUnfold->SetMarkerColor( kRed );
  hUnfold->SetMarkerSize( 2 );
  hUnfold->DrawCopy("edatasame");
  leg->AddEntry( hUnfold, TString::Format("Data unfolded") , "p");

  //-- Unfold the MC.
  TH1D* hMeas_MC;
  if( variable == "all" ){ 	Get_DetEnergy(0, hMeas_MC); 	}
  else if( variable == "lead" ){Get_DetEnergy_lead(0, hMeas_MC);}

  RooUnfoldBayes unfold_bayes_MC(response, hMeas_MC, 40);
	
  hMeas_MC->Scale( scalefactors_data_ / scalefactors_MC_ );

  unfold_bayes_MC.SetVerbose(-1);
  TH1D* hUnfold_MC = (TH1D*) unfold_bayes_MC.Hreco( RooUnfold::kCovariance ); 
  hUnfold_MC->SetMarkerStyle( 20 );
  hUnfold_MC->SetMarkerColor( kGreen );
  hUnfold_MC->SetMarkerSize( 2 );
  hUnfold_MC->DrawCopy("edatasame");
  leg->AddEntry( hUnfold_MC, TString::Format("MC unfolded") , "p");  

  leg->Draw();  
  can->SaveAs( TString::Format(folder_ + "/Unfolding_genLevel_" + variable + "-spectrum_data-and-MC_%iGeV.C", (int)Ethresh_ ) );
  can->SaveAs( TString::Format(folder_ + "/Unfolding_genLevel_" + variable + "-spectrum_data-and-MC_%iGeV.pdf", (int)Ethresh_ ) );
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
    
    // cout << "Min and max\t" << min_val << "\t" << max_val << endl;
    
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
    cout << "Loop iteration\t" << file_-1 << endl;

    TString filename;
    if( (file_-1) < MC_files_.size() && file_ > 0){ filename = MC_files_[ (file_-1) ]; }
    if( (file_-1) == -1 ){ filename = datafile_; legendoptions = "p";}
 
    /*  
    TFile* _file = TFile::Open( filename, "read" );   
    if(! _file->GetListOfKeys()->Contains( variable ) ){ continue; }
    Get_Distribution( (file_-1) , hDistr, variable );
    hDistr->Rebin( 4 );
    */

    if( variable == "all"){ Get_DetEnergy( file_-1 , hDistr ); }
    else{ Get_DetEnergy_lead( file_-1 , hDistr ); }

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


    if( variable == "all"){ Get_DetEnergy( file_-1 , hDistr ); }
    else{ Get_DetEnergy_lead( file_-1 , hDistr ); }


//    Get_Distribution( (file_-1) , hDistr, variable );
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


void Unfolder::Get_DetUnfolded(int file_, int MC_, TH1D* &hist_, TMatrixD& covariance_m, int iterations, TString variable){
  cout << "Unfolder::Get_DetUnfolded(" << file_ << ", " << MC_ << ", "  << iterations << ", " << variable << ") with Matrix" << endl;

  TFile *_file0 = new TFile( MC_files_[MC_], "read");
  RooUnfoldResponse* response;
  if( variable == "all"){ response = (RooUnfoldResponse*)_file0->Get("response"); }
  if( variable == "lead"){ response = (RooUnfoldResponse*)_file0->Get("response_lead"); }

  // Use hist as input for the unfolding.
  RooUnfoldBayes unfold_bayes(response, hist_, iterations);
  // the original hist has become obsolete, transform it into the unfolded histogram.
  hist_ = (TH1D*) unfold_bayes.Hreco( RooUnfold::kCovariance ); 

  RooUnfold::ErrorTreatment et = RooUnfold::kCovariance;

//  covariance_m = (TMatrixD) unfold_bayes.Ereco( et );

}




void Unfolder::Get_DetUnfolded(int file_, TH1D* &hist_, int iterations, TString variable){
  cout << "Unfolder::Get_DetUnfolded(" << file_ << ", "  << iterations << ", " << variable << ") no Matrix" << endl;

  TFile *_file0 = new TFile( MC_files_[0], "read");
  RooUnfoldResponse* response;
  if( variable == "all"){ response = (RooUnfoldResponse*)_file0->Get("response"); }
  if( variable == "lead"){ response = (RooUnfoldResponse*)_file0->Get("response_lead"); }
  
  // Use hist as input for the unfolding.
  RooUnfoldBayes unfold_bayes(response, hist_, iterations);
  // the original hist has become obsolete, transform it into the unfolded histogram.
  hist_ = (TH1D*) unfold_bayes.Hreco(); 
}




void Unfolder::Get_GenSmeared(int file_, int MC_, TH1D* &hist_, TH1D* hGen, TString variable){
  cout << "Unfolder::Get_GenSmeared(" << file_ << ", " << MC_ << ", " << variable << ")" << endl;
  TFile *_file0 = new TFile( MC_files_[ MC_ ], "read");
  RooUnfoldResponse* response;
  if( variable == "all"){ response = (RooUnfoldResponse*)_file0->Get("response"); }
  if( variable == "lead"){ response = (RooUnfoldResponse*)_file0->Get("response_lead"); }

  hist_ = (TH1D*) response->ApplyToTruth( hGen, "Smearing");
}



void Unfolder::Get_UnfoldSmearError( TH1D* vUnfold, TH1D* &vSmeared, int MC_, TString variable, int iterations ){
  cout << "Get_UnfoldSmearError" << endl;

  TFile *_file0 = new TFile( MC_files_[ MC_ ], "read");
  RooUnfoldResponse* response;
  if( variable == "all"){ response = (RooUnfoldResponse*)_file0->Get("response"); }
  if( variable == "lead"){ response = (RooUnfoldResponse*)_file0->Get("response_lead"); }

  // -- The response matrix.
  TMatrixD& responsematrix = (TMatrixD&) response->Mresponse();

  // -- Construct the RooUnfolfBayes object to obtain the covariance matrix.
//  TMatrixD& covariancematrix = (TMatrixD&) response->Ereco();
  
  for(int bin_smear = 0; bin_smear < vSmeared->GetNbinsX(); bin_smear++){
//    cout << "USE---\t" << bin_smear << "/" << vSmeared->GetNbinsX() << endl;

    double err = 0;
    for(int col = 0; col < responsematrix.GetNcols(); col++){
//      cout << "USE---\t\t" << col << "/" << responsematrix.GetNcols() << endl;
      double res = responsematrix(bin_smear, col);
      double unf = vUnfold->GetBinContent(col);
      double s_u = vUnfold->GetBinError(col);
      double s_r = responsematrix(bin_smear, col);

      err += unf*unf*s_r*s_r + res*res*s_u*s_u;

    }// Loop over truth bins.
    vSmeared->SetBinError( bin_smear, sqrt(err) );   
  }// Loop over smear vector.
}



/////////////////////////////////////
// Unfold-and-smearing on its own. //
/////////////////////////////////////

void Unfolder::Unfold_and_smear(int file_, TH1D* &hist_, int MC_, int iterations, TString variable, int method, double &chi2 ){
  cout << "Unfolder::Unfold-and-smear\t" << MC_files_[ MC_ ] << endl;
  variable = "all";

  ofstream smearing_matrix;
  smearing_matrix.open("Smearing_matrix.txt");

  TH1D* theHist = (TH1D*)hist_->Clone("theHist");

  TFile *_file0 = new TFile( MC_files_[ MC_ ], "read");
  TFile *_file_det = new TFile( datafile_, "read");
  RooUnfoldResponse* response;
  
  if( method == 1 ){
    if( variable == "all"){ response = (RooUnfoldResponse*)_file0->Get("response"); }
    else if( variable == "lead"){ response = (RooUnfoldResponse*)_file0->Get("response_lead"); }
 }
 else if( method == 2 ){
    if( variable == "all"){ response = (RooUnfoldResponse*)_file0->Get("response_match"); }
    else if( variable == "lead"){ response = (RooUnfoldResponse*)_file0->Get("response_lead"); }
 }

  //-- Extract histograms.
  TH1D* hMiss = (TH1D*)_file0->Get("hCastorJet_miss_all");
  TH1D* hFake = (TH1D*)_file0->Get("hCastorJet_fake_all");	
  TH2D* hResponse 	= (TH2D*)response->Hresponse(); 
  TH1D* hMeasured	= (TH1D*)response->Hmeasured();
  TH1D* hTruth		= (TH1D*)response->Htruth();
  TH1D* hSmear;

  // -- S = R * (U - M) + F

  //-- Method 1.
  if( method == 1 ){  
    // -- Construct the RooUnfolfBayes object to obtain the covariance matrix and the unfolded distribution.
    RooUnfoldBayes unfold_bayes(response, theHist, iterations); 
    
    unfold_bayes.SetVerbose(0);
    TH1D* hUnfold = (TH1D*) unfold_bayes.Hreco( RooUnfold::kCovariance ); 

    TCanvas *can___ = new TCanvas("can___", "can___", 1.);
    hUnfold->Draw("hist");
    TH1D* hGen = (TH1D*)response->Htruth();
    hGen->SetLineColor( kRed );
    hGen->Draw("histsame");
    can___->SaveAs("tessssssst.C");

 
//    hMiss->Scale( 1./renorm_ );
    hFake->Scale( 1./renorm_ );   
 
    hSmear = (TH1D*) response->ApplyToTruth( hUnfold );
    hSmear->Add( hFake );    

//  Calculate_smearedBackError(-1, 0, iterations, hSmear, hUnfold);
    chi2 = Calculate_smearedBackError_covariance( theHist, hUnfold, response, iterations );
  }


  // Method 2.
  //--- Alternative way of calculating.
  //--- Response matrix/object only contain matched jets, treat misses and fakes separately.

  if( method == 2 ){

    hMiss->Scale( 1./renorm_ );
    hFake->Scale( 1./renorm_ );

    // Substract fakes (scaled MC) from data to obtain matched jet pair energies.
    theHist->Add( hFake, -1. );	
    // Unfold the spectrum.
    RooUnfoldBayes unfold_bayes(response, theHist, iterations); 
    unfold_bayes.SetVerbose(0);
    TH1D* hUnfold = (TH1D*) unfold_bayes.Hreco( RooUnfold::kCovariance ); 

    hResponse->Scale( 1./renorm_ );

    hSmear = (TH1D*) response->ApplyToTruth( hUnfold );

    hSmear->Add( hFake );  
  }


  hist_ = hSmear;



  cout << "Done smearing" << endl;
  cout << "\t\thTruth\t" << hTruth->GetXaxis()->GetTitle() << endl;
}









void Unfolder::Get_GenSmeared(int file_, TH1D* &hist_, TH1D* hGen, TString variable){
  cout << "Unfolder::Get_GenSmeared(" << file_ << ", " << variable << ")" << endl;

  TFile *_file0 = new TFile( MC_files_[ 0 ], "read");
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

void Unfolder::PrepareCanvas_2D( TCanvas* &can_, TString label){

  cout << "\n\n\n// --  Prepare canvas 2D" << endl;

  can_ = new TCanvas( label, label, 900, 900 );

  can_->SetLeftMargin(0.17);
  can_->SetTopMargin(0.09);
  can_->SetBottomMargin(0.14);
  can_->SetRightMargin(0.15);
}


void Unfolder::Prepare_2Dplot(TH2D* &hist_){

  hist_->GetYaxis()->SetTitleOffset(1.30);
  hist_->GetYaxis()->SetTitleSize(0.06);
  hist_->GetYaxis()->SetLabelSize(0.04);

  hist_->GetXaxis()->SetLabelSize(0.04);
  hist_->GetXaxis()->SetTitleSize(0.06);
  hist_->GetXaxis()->SetNdivisions(304);

  double min = GetMinimumValue( hist_ );
  double max = hist_->GetMaximum();

  hist_->GetZaxis()->SetRangeUser( min*0.9, max*1.1);

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





void Unfolder::Unfolding_data(TString variable, int iterations_){
  // Prepare variables and objects.
  TH1D* hist_original;
  TH1D* hist_reference;
  TH1D* hist_result;  
  TString htitle;
  TLegend *leg = new TLegend(0.65, 0.45, 0.95, 0.95);
  iterations_ = 16;
  double xaxisgraph[iterations_], yaxisgraph[iterations_];

  vector<TH1D*> histos;
  
  // Get the detector level distribution.
  if( variable == "lead"){	Get_DetEnergy_lead(-1, hist_original );	htitle = "E_{det} spectrum (leading)";}
  else{				Get_DetEnergy(-1, hist_original );	htitle = "E_{det} spectrum (all)";	}

  if( variable == "lead"){	Get_GenEnergy_lead(0, hist_reference );htitle = "E_{det} spectrum (leading)";}
  else{				Get_GenEnergy(0, hist_reference );	htitle = "E_{det} spectrum (all)";	}

  hist_reference->Scale( hist_original->Integral() / hist_reference->Integral() );

  if( normalise_1 ){ hist_original->Scale( 1./hist_original->Integral() ); }
  if( normalise_1 ){ hist_reference->Scale( 1./hist_reference->Integral() ); }

  // Set detector level properties and legend.
  hist_original->SetLineColor( getColor( 1 ) );
  	hist_original->SetLineStyle( 1 );
  	hist_original->SetLineWidth( 3 );
  	hist_original->SetMarkerColor(  getColor( 1 ) );

  hist_reference->SetMarkerColor( kBlack );

  leg->AddEntry( hist_reference, "Pythia6 (Z2*)", "p");
  histos.push_back( hist_reference );

  for(int iterations = 1; iterations <= iterations_; iterations++){

    hist_result = (TH1D*)hist_original->Clone(TString::Format("Unfolding_%i_iterations", iterations) );
    hist_result->SetName( TString::Format("Unfolding_%i_iterations", iterations) );
    hist_result->SetTitle( TString::Format("Unfolding_%i_iterations", iterations) );
    Get_DetUnfolded( -1, hist_result, iterations, variable );

    hist_result->SetLineColor( getColor( iterations + 1) );
    hist_result->SetLineStyle( iterations + 1 );
    hist_result->SetLineWidth( 3 );
    hist_result->SetMarkerColor(  getColor( iterations + 1) );

    leg->AddEntry( hist_result, TString::Format("Treated data, %i it.", iterations) , "p");
    histos.push_back( hist_result );

  }

  // Make sure to set histogram axes.
  hist_original->GetXaxis()->SetTitle("E_{det} [GeV]");
  hist_original->GetYaxis()->SetTitle("#frac{1}{N}.#frac{dN}{dE}");

  // Plot the distributions and their ratios.
  TCanvas *can;
  PrepareCanvas( can, "Unfolding_data" + variable);
  MakeDoublePaddedComparison(can, histos, leg );
  can->SaveAs(folder_ + "/Unfolding_data_" + variable + label_ + ".C");
  can->SaveAs(folder_ + "/Unfolding_data_" + variable + label_ + ".pdf");
}








void Unfolder::ClosureTest_data(TString variable, TString file, int method){
  ofstream iterations_and_errors;
  iterations_and_errors.open("Iterations_and_errors.txt");

  // Prepare variables and objects.
  TH1D* hist_original, *hist_MCdet;
  TH1D* hist_result; 
  TH1D* hist_reference; 
  TString htitle;
  TLegend *leg = new TLegend(0.65, 0.45, 0.95, 0.95);
    leg->SetFillColor(0);
  int iterations_ = 5;
  int actual_iterations = 0;


  double xaxisgraph[iterations_], yaxisgraph[iterations_];

  TString norm = "notNorm";

  vector<TH1D*> histos;
  TMatrixD cov_m;

  int current_color = 1;
  
  // Get the detector level distribution.
  if( variable == "lead"){	Get_DetEnergy_lead(-1, hist_original );	htitle = "E_{det} spectrum (leading)";}
  else{				Get_DetEnergy(-1, hist_original );	htitle = "E_{det} spectrum (all)";	}

  if( variable == "lead"){	Get_DetEnergy_lead(0, hist_MCdet );	htitle = "E_{det} spectrum (leading)";}
  else{				Get_DetEnergy(0, hist_MCdet );		htitle = "E_{det} spectrum (all)";	}  
 
  SetCastorJetEnergy_norm( static_cast<double>(547922)/ static_cast<double>(4661641) ); 

  // Hist_reference is the same as the original, but with modified lower edge.
  GetSubHistogram( hist_original, hist_reference, Eplotmin_, 2000.);


  hError_hist = new TH2D("hError_hist", "hError_hist;E;iteration", hist_reference->GetNbinsX(), hist_reference->GetBinLowEdge( 1 ), hist_reference->GetXaxis()->GetBinUpEdge(  hist_reference->GetNbinsX() ), iterations_, 0., iterations_ );

   TCanvas* can_ = new TCanvas("can_", "can_", 1.);
   hist_original->Draw("hist");
   hist_reference->SetLineColor( kRed );
   hist_reference->SetLineStyle( 2 );
   hist_reference->Draw("histsame");
   can_->SaveAs("compare_subhistograms.C");

  // Set detector level properties and legend.
  hist_reference->SetLineColor( getColor( 1 ) );
  	hist_reference->SetLineStyle( 1 );
  	hist_reference->SetLineWidth( 3 );
  	hist_reference->SetMarkerColor(  getColor( 1 ) );
	hist_reference->SetMarkerStyle( 20 );


  leg->AddEntry( hist_reference, "Actual data", "p");
  histos.push_back( hist_reference );

  // Loop over number of Bayesian iterations, and unfold->smear detector level distribution.
  // Calculate chi2 as we go.
  double chi2 = hist_original->Chi2Test( hist_reference, "CHI2/NDF");
  double chi2_prev = chi2;

  //---------------------------------------------------------------//
  //-- Loop over the files, if there are several to be unfolded. --//
  //---------------------------------------------------------------//

  for(int file_ = 0; file_ < MC_files_.size(); file_++){
    if( file != printLabel_[ MC_files_[file_] ] ){ continue; }

    //--------------------------------------------------//
    //-- Loop over the number of Bayesian iterations. --//
    //--------------------------------------------------//

    int increase_iterations = 3;
    for(int iterations = 1; iterations <= iterations_; iterations+=increase_iterations){
      // Determine the inrease in Bayesian iterations.
      if(iterations >= 10 ){ increase_iterations = 5; }
      if(iterations >= 20 ){ increase_iterations = 15; }
      if(iterations >= 100 ){ increase_iterations = 20; }

      hist_result = (TH1D*)hist_original->Clone(TString::Format("Closure_%i_iterations", iterations) );
      hist_result->SetName( TString::Format("Closure_%i_iterations", iterations) );
      hist_result->SetTitle( TString::Format("Closure_%i_iterations", iterations) );

      double chi2;
      Unfold_and_smear( -1, hist_result, file_, iterations, variable, method, chi2 );

      can_->cd();
      hist_reference->Draw("hist");
      hist_result->SetLineStyle( 5 );
      hist_result->DrawCopy("histsame");

      GetSubHistogram( hist_result, hist_result, Eplotmin_, 2000.);

      hist_result->SetLineColor( getColor( current_color ) );
      hist_result->SetLineStyle( iterations + 1 );
      hist_result->SetLineWidth( 3 );
      hist_result->SetMarkerColor(  getColor( current_color ) );

      hist_result->DrawCopy("histsame");
      can_->SaveAs("compare_2.C");

      for(int binx = 0; binx <= hist_result->GetNbinsX(); binx++){

        double bin_err = hist_result->GetBinContent( binx );
        hError_hist->SetBinContent( binx, iterations, bin_err );

      }

      leg->AddEntry( hist_result, TString::Format("Treated data, %i it.", iterations) , "l");
      histos.push_back( hist_result );

      //-- Chi2 test.
//      chi2 = Chi2_test( hist_reference, hist_result);
//      chi2 = hist_reference->Chi2Test( hist_result, "CHI2/NDF");
      xaxisgraph[actual_iterations] = iterations;
      yaxisgraph[actual_iterations] = chi2;    
      chi2_prev = chi2;

      current_color++;
      actual_iterations++;

    }

    // Make sure to set histogram axes.
    hist_original->GetXaxis()->SetTitle("E_{det} [GeV]");
    hist_original->GetYaxis()->SetTitle("#frac{1}{N}.#frac{dN}{dE}");

    // Plot the distributions and their ratios.
    TCanvas *can;
    PrepareCanvas( can, "ClosureTest_data" + variable);
    MakeDoublePaddedComparison(can, histos, leg );
    can->SaveAs( TString::Format(folder_ + "ClosureTest_data_method_%i_" + variable + label_ + "_%i_iterations_" + norm + ".C", method, iterations_) );
    can->SaveAs( TString::Format(folder_ + "ClosureTest_data_method_%i_" + variable + label_ + "_%i_iterations" + norm + ".pdf", method, iterations_) );

    // Finish chi2 study.
    TCanvas *can_chi2;
    PrepareCanvas( can_chi2, "CHI2_Test_" + variable);
    TGraph* chi2_evolution = new TGraph(actual_iterations, xaxisgraph, yaxisgraph);
    chi2_evolution->GetXaxis()->SetTitle("N_{it.}");
    chi2_evolution->GetYaxis()->SetTitle("#chi^{2}/NDF");
    chi2_evolution->Draw("A*");
    can_chi2->SaveAs( TString::Format(folder_ + "Chi2_Test_data_method_%i_" + variable + label_ + "_%i_iterations_" + norm + ".C", method, iterations_) );
    can_chi2->SaveAs( TString::Format(folder_ + "Chi2_Test_data_method_%i_" + variable + label_ + "_%i_iterations" + norm + ".pdf", method, iterations_) );

    iterations_and_errors << "€€€ We plot from\t" << Eplotmin_ << "\tto\t" << hist_original->GetBinLowEdge( hist_original->GetNbinsX()+1 ) << "\toriginal hist has\t" << hist_original->GetNbinsX() << "\t bins" << "\nFrom file\t" << datafile_ << endl;
  }

  iterations_and_errors.close();

/*
  for(int bin_ = 1; bin_ <= hist_original->GetNbinsX(); bin_++){ cout << "bin_\t" << bin_ << "\t" << hist_original->GetBinContent( bin_ ) << "\t" << hist_original->GetBinError( bin_ ) << endl;  }
  for(int bin_ = 1; bin_ <= hist_reference->GetNbinsX(); bin_++){ cout << "bin_\t" << bin_ << "\t:" << hist_reference->GetBinContent( bin_ ) << "\t" << hist_reference->GetBinError( bin_ ) << endl; }
*/

  TCanvas* canerror = new TCanvas("can", "can", 1.);
  hError_hist->Draw("colz");
  canerror->SaveAs("error_plot.C");
}













void Unfolder::ClosureTest_MC_detLevel(TString variable){
cout << "++++++++++++++++++++++++++++++++++++++This is the function++++++++++++++++++++++++++++++++++++++" << endl;
  ofstream iterations_and_errors;
  iterations_and_errors.open("Iterations_and_errors.txt");
  TString file = "Displaced";
  int method = 1;

  // Prepare variables and objects.
  TH1D* hist_original, *hist_MCdet;
  TH1D* hist_result; 
  TH1D* hist_reference; 
  TString htitle;
  TLegend *leg = new TLegend(0.65, 0.45, 0.95, 0.95);
    leg->SetFillColor(0);
  int iterations_ = 80;
  double xaxisgraph[iterations_], yaxisgraph[iterations_];

  TString norm = "notNorm";

  vector<TH1D*> histos;
  TMatrixD cov_m;

  int current_color = 1;
  
  // Get the detector level distribution.
  if( variable == "lead"){	Get_DetEnergy_lead(0, hist_original );	htitle = "E_{det} spectrum (leading)";}
  else{				Get_DetEnergy(0, hist_original );	htitle = "E_{det} spectrum (all)";	}
 
  SetCastorJetEnergy_norm( 1. ); 
 
  // Hist_reference is the same as the original, but with modified lower edge.
  GetSubHistogram( hist_original, hist_reference, Eplotmin_, 2000.);

  TCanvas* can_ = new TCanvas("can_", "can_", 1.);
  hist_original->Draw("hist");
  hist_reference->SetLineColor( kRed );
  hist_reference->SetLineStyle( 2 );
  hist_reference->Draw("histsame");
  can_->SaveAs("compare_subhistograms.C");

  //  Rebin_to( hist_reference, 40);

  // Set detector level properties and legend.
  hist_reference->SetLineColor( getColor( 1 ) );
  	hist_reference->SetLineStyle( 1 );
  	hist_reference->SetLineWidth( 3 );
  	hist_reference->SetMarkerColor(  getColor( 1 ) );
	hist_reference->SetMarkerStyle( 20 );

  leg->AddEntry( hist_reference, "Actual data", "p");
  histos.push_back( hist_reference );

  // Loop over number of Bayesian iterations, and unfold->smear detector level distribution.
  // Calculate chi2 as we go.
  double chi2 = hist_original->Chi2Test( hist_reference, "CHI2/NDF");
  double chi2_prev = chi2;

  for(int file_ = 0; file_ < MC_files_.size(); file_++){
    if( file != printLabel_[ MC_files_[file_] ] ){ continue; }

    int increase_iterations = 3;
    for(int iterations = 1; iterations <= iterations_; iterations+=increase_iterations){
      if(iterations >= 10 ){ increase_iterations = 5; }

      hist_result = (TH1D*)hist_original->Clone(TString::Format("Closure_%i_iterations", iterations) );
      hist_result->SetName( TString::Format("Closure_%i_iterations", iterations) );
      hist_result->SetTitle( TString::Format("Closure_%i_iterations", iterations) );

      double chi2;
      Unfold_and_smear( -1, hist_result, file_, iterations, variable, method, chi2 );

      can_->cd();
      hist_reference->Draw("hist");
      hist_result->SetLineStyle( 5 );
      hist_result->DrawCopy("histsame");

      GetSubHistogram( hist_result, hist_result, Eplotmin_, 2000.);

      iterations_and_errors << "Unfolder::ClosureTest_MC_detLevel\titeration\t" << iterations << "\t" << hist_result->GetBinLowEdge( 1 ) << endl;

      hist_result->SetLineColor( getColor( current_color ) );
      hist_result->SetLineStyle( iterations + 1 );
      hist_result->SetLineWidth( 3 );
      hist_result->SetMarkerColor(  getColor( current_color ) );

      hist_result->DrawCopy("histsame");
      can_->SaveAs("compare_2.C");

      leg->AddEntry( hist_result, TString::Format("Treated data, %i it.", iterations) , "l");
      histos.push_back( hist_result );

      //-- Chi2 test.
//      chi2 = Chi2_test( hist_reference, hist_result);
      chi2 = hist_reference->Chi2Test( hist_result, "CHI2/NDF");
      xaxisgraph[iterations] = iterations;
      yaxisgraph[iterations] = chi2;    
      chi2_prev = chi2;

      current_color++;

    }

    // Make sure to set histogram axes.
    hist_original->GetXaxis()->SetTitle("E_{det} [GeV]");
    hist_original->GetYaxis()->SetTitle("#frac{1}{N}.#frac{dN}{dE}");

    // Plot the distributions and their ratios.
    TCanvas *can;
    PrepareCanvas( can, "ClosureTest_data" + variable);
    MakeDoublePaddedComparison(can, histos, leg );
    can->SaveAs( TString::Format(folder_ + "ClosureTest_MC_detLevel_method_%i_" + variable + label_ + "_%i_iterations_" + norm + ".C", method, iterations_) );
    can->SaveAs( TString::Format(folder_ + "ClosureTest_MC_detLevel_method_%i_" + variable + label_ + "_%i_iterations" + norm + ".pdf", method, iterations_) );

    // Finish chi2 study.
    TCanvas *can_chi2;
    PrepareCanvas( can_chi2, "CHI2_Test_" + variable);
    TGraph* chi2_evolution = new TGraph(iterations_, xaxisgraph, yaxisgraph);
    chi2_evolution->GetXaxis()->SetTitle("N_{it.}");
    chi2_evolution->GetYaxis()->SetTitle("#chi^{2}/NDF");
    chi2_evolution->Draw("A*");
    can_chi2->SaveAs( TString::Format(folder_ + "Chi2_Test_MC_detLevel_method_%i_" + variable + label_ + "_%i_iterations_" + norm + ".C", method, iterations_) );
    can_chi2->SaveAs( TString::Format(folder_ + "Chi2_Test_MC_detLevel_method_%i_" + variable + label_ + "_%i_iterations" + norm + ".pdf", method, iterations_) );

    iterations_and_errors << "€€€ We plot from\t" << Eplotmin_ << "\tto\t" << hist_original->GetBinLowEdge( hist_original->GetNbinsX()+1 ) << "\toriginal hist has\t" << hist_original->GetNbinsX() << "\t bins" << "\nFrom file\t" << datafile_ << endl;
  }

  iterations_and_errors.close();

  for(int bin_ = 1; bin_ <= hist_original->GetNbinsX(); bin_++){

    cout << "bin_\t" << bin_ << "\t" << hist_original->GetBinContent( bin_ ) << "\t" << hist_original->GetBinError( bin_ ) << endl;

  }
  for(int bin_ = 1; bin_ <= hist_reference->GetNbinsX(); bin_++){
    cout << "bin_\t" << bin_ << "\t" << hist_reference->GetBinContent( bin_ ) << "\t" << hist_reference->GetBinError( bin_ ) << endl;

  }
}





























/*
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
//    Get_DetUnfolded( 0, hist_result, iterations, variable );
//    Get_GenSmeared( 0, hist_result , hist_result, variable); 

    Unfold_and_smear( -1, hist_result, 0, iterations, variable, 1 );

    //hist_int = (TH1D*)hist_result->Clone(TString::Format("it_unfolded_%i", iterations) );
    //list_of_iterations.push_back( hist_int);
    hist_result->SetLineColor( getColor( iterations + 1 ) );
    hist_result->SetLineStyle( iterations + 1 );
    hist_result->SetLineWidth( 3 );
    hist_result->SetMarkerColor(  getColor( iterations+1 ) );

    histos.push_back( hist_result );
    leg->AddEntry( hist_result, TString::Format("Treated data, %i it.", iterations) , "l");

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
  can->SaveAs(folder_ + "ClosureTest_MC_detLevel_" + variable + label_ + ".C");
  can->SaveAs(folder_ + "ClosureTest_MC_detLevel_" + variable + label_ + ".pdf");

  // Finish chi2 study.
  TCanvas *can_chi2;
  PrepareCanvas( can_chi2, "CHI2_Test_" + variable);
  TGraph* chi2_evolution = new TGraph(iterations_, xaxisgraph, yaxisgraph);
  chi2_evolution->GetXaxis()->SetTitle("N_{it.}");
  chi2_evolution->GetYaxis()->SetTitle("#chi^{2}/NDF");
  chi2_evolution->Draw("A*");
  can_chi2->SaveAs(folder_ + "Chi2_Test_MC_det_" + variable + label_ + ".C");
  can_chi2->SaveAs(folder_ + "Chi2_Test_MC_det_" + variable + label_ + ".pdf");
}
*/



void Unfolder::ClosureTest_MC(TString variable){

  // Prepare variables and objects.
  TH1D* hist_original;
  TH1D* hist_result;  
  TString htitle;
  TLegend *leg = new TLegend(0.65, 0.45, 0.95, 0.95);
  int iterations_ = 5;
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
  can->SaveAs(folder_ + "ClosureTest_MC_" + variable + label_ + ".C");
  can->SaveAs(folder_ + "ClosureTest_MC_" + variable + label_ + ".pdf");

  // Finish chi2 study.
  TCanvas *can_chi2;
  PrepareCanvas( can_chi2, "CHI2_Test_" + variable);
  TGraph* chi2_evolution = new TGraph(iterations_, xaxisgraph, yaxisgraph);
  chi2_evolution->GetXaxis()->SetTitle("N_{it.}");
  chi2_evolution->GetYaxis()->SetTitle("#chi^{2}/NDF");
  chi2_evolution->Draw("A*");
  can_chi2->SaveAs(folder_ + "Chi2_Test_MC_" + variable + label_ + ".C");
  can_chi2->SaveAs(folder_ + "Chi2_Test_MC_" + variable + label_ + ".pdf");
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
  TString drawoptions = "pehist";
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
   
    drawoptions = "ehistsame";
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
  cout << "My Chi2" << endl;


  ofstream chi2_printout;
  chi2_printout.open("chi2_printout.txt", ios::out | ios::app | ios::binary);
  
  chi2_printout << "bin \th1 BinCenter \th2 BinCenter \tH1 \th1 \te1 \tH2 \th2 \te2 \t\tsigma \tdelta \tchi2" << endl;

  int ndf = hist_data->GetNbinsX();
  double chi2 = 0.;

  double H1 = hist_data->Integral( );
  double H2 = hist_MC->Integral();

  for( int bin = 1; bin <= ndf; bin++){
   double h1 = hist_data->GetBinContent( bin );
   double h2 = hist_MC->GetBinContent( bin );

   if( h1 == 0. && h2 == 0. ){ --ndf; continue; }

   double e1 = hist_data->GetBinError( bin );
   double e2 = hist_MC->GetBinError( bin );

   double sigma = H1*H1 * e2*e2 + H2*H2 * e1*e1;
   double delta = H1*h2 - H2*h1;

   chi2 += delta*delta/sigma;

   chi2_printout << bin << "\t" << hist_data->GetBinCenter( bin ) << "\t" << hist_MC->GetXaxis()->GetBinCenter( bin ) << "\t" << H1 << "\t" << h1 << "\t" << e1 << "\t" << H2 << "\t" << h2 << "\t" << e2 << "\t\t" << sigma << "\t" << delta << "\t" << chi2 << endl;
  }

  chi2_printout << "\n\tCHI2/NDF\t" << chi2/ndf << endl << endl << endl;

  chi2_printout.close();

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

  can_chi2->SaveAs(folder_ + "Chi2_comparison_data_" + variable + label_ + ".C");
  can_chi2->SaveAs(folder_ + "Chi2_comparison_data_" + variable + label_ + ".pdf");
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

  can_chi2->SaveAs(folder_ + "Chi2_comparison_MC_det_" + variable + label_ + ".C");
  can_chi2->SaveAs(folder_ + "Chi2_comparison_MC_det_" + variable + label_ + ".pdf");
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
  can_chi2->SaveAs(folder_ + "Chi2_comparison_MC_" + variable + label_ + ".C");
  can_chi2->SaveAs(folder_ + "Chi2_comparison_MC_" + variable + label_ + ".pdf");
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

  int Slice_threshold = 0.;
  int event_threshold_ = 50.;


  // -- Prepare variables, canvasses, ...
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

  // -- The fit function.

   TF1* analytical = new TF1("Analytical", "  ( [0] + [1] * log( [2] + x) )", fitting_threshold_, 900.);
    analytical->SetParLimits(2, fitting_threshold_+1., 1000.);

  /* Open files. */
  for( int bin_phi = first_sector; bin_phi <=last_sector; bin_phi++){//  hResponse->GetNbinsZ(); bin_phi++){

    for(int file_ = 0; file_ < MC_files_.size(); file_++){

      TFile *_file0 = TFile::Open( MC_files_[file_], "Read");
      TString label = TString::Format( printLabel_[ MC_files_[file_] ] + "_sector_%i_to_%i", first_sector, last_sector);

      label = "Test_this_" + label;
      TString label_short = "Test_this";
      int	new_dir = mkdir( (folder_).Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  );
      new_dir = mkdir( (folder_).Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  );

      TH2D* hGenE_selection;
      TH2D* hDetE_selection;
      TH2D* hResponse_selection;
      TH2D* hDetE;
      TH2D* hGenE;
      TH1D* hGen;
      TString setup = "good_sectors";

      if( setup_calibration_ == "all_sectors") 	setup = "all_sectors";
      else if( setup_calibration_ == "separate_sectors") setup = TString::Format("%i", bin_phi); 

       // -- Extract the response, Edet and Egen distributions from the files.
      Extract_2D_energy_distributions(file_, hGenE_selection, hDetE_selection, hResponse_selection, hDetE, hGenE, hGen, setup);
            	 
      TGraphErrors * Response_true;
      TGraphErrors * Response_meas;

      cout << "CALIBRATIONFUNCTION - Analyze_response\tfile\t" << file_ << endl;

      Analyze_response( hResponse_selection, hGenE_selection, hDetE_selection, hGenE, label, hGen, file_, Response_meas, Response_true);

      cout << "CALIBRATIONFUNCTION - Analyze_response done" << endl;

      TString legend_info;
      if( setup_calibration_ == "all_sectors") { legend_info = "All sectors"; }
      else if( setup_calibration_ == "good_sectors"){ legend_info = "Good sectors"; }
      else if( setup_calibration_ == "separate_sectors"){ legend_info = TString::Format("Sector %i", bin_phi); }
       legend_det->AddEntry( Response_meas, legend_info, "lp");
 
       // Draw.
       can_true->cd();	 
       Response_true->Draw("p" + drawoptions);
       line->Draw();

	can_true->SaveAs(TString::Format(folder_ + "/CalibrationFactors_true_calib_phi_%i_" + label + ".C", bin_phi) );
        can_true->SaveAs(TString::Format(folder_ + "/CalibrationFactors_true_calib_phi_%i_.pdf", bin_phi) );

	// Draw.

       can_fit->cd();
       Response_meas->Draw("p" + drawoptions);
       drawoptions="same";	

       can_meas->cd(); 
       Response_meas->Draw("ape");// + drawoptions);
       line->Draw();

       cout << "Sector\t" << bin_phi << endl;
       analytical->SetLineColor( getColor(bin_phi) );

       analytical->SetLineColor( getColor( file_ + 1) );
       Response_meas->Fit( analytical, "", "", fitting_threshold_, 900.); 
       if( first_sector >= 12 && last_sector <= 15){ can_fit->cd(); analytical->DrawCopy("lsame"); }

       can_fit->SaveAs(TString::Format(folder_ + "/CalibrationFactors_fit_calib_phi_%i_" + label + ".C", bin_phi) );
       can_fit->SaveAs(TString::Format(folder_ + "/CalibrationFactors_fit_calib_phi_%i_" + label + ".pdf", bin_phi) );

      double mean_alpha_nom = 0., mean_alpha_denom = 0.;
      double mean_beta_nom = 0., mean_beta_denom = 0.;
      double mean_gamma_nom = 0., mean_gamma_denom = 0.;         
//       if( bin_phi < 12 || bin_phi > 15){
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
//       }

       can_meas->SaveAs(TString::Format(folder_ + "/CalibrationFactors_meas_calib_phi_%i_" + label + ".C", bin_phi) );
       can_meas->SaveAs(TString::Format( folder_ + "/CalibrationFactors_meas_calib_phi_%i_" + label + ".pdf", bin_phi) );

	// Delete the histograms to avoid memory leak.
       hGenE_selection->~TH2();
       hDetE_selection->~TH2();
       hResponse_selection->~TH2();

       std::vector<double> parameters_calibration;
         parameters_calibration.push_back( mean_alpha_nom/mean_alpha_denom );
         parameters_calibration.push_back( mean_beta_nom/mean_beta_denom );
         parameters_calibration.push_back( mean_gamma_nom/mean_gamma_denom );

       calibration_parameters_[ MC_files_[file_] ] = parameters_calibration;
       calibration_graphs_[ MC_files_[file_] ] = Response_meas;

       // Save.
       can_true->cd();
       legend_gen->Draw();
       	can_true->SaveAs( folder_ +"/CalibrationFactors_true_" + calib_ + "_" + label + ".C");		
	can_true->SaveAs( folder_ +"/CalibrationFactors_true_" + calib_ + "_" + label + ".pdf");
     
       can_meas->cd();
       legend_det->Draw();
     	 can_meas->SaveAs( folder_ +"/CalibrationFactors_meas_" + calib_ + "_" + label + ".C"); 		
	 can_meas->SaveAs( folder_ +"/CalibrationFactors_meas_" + calib_ + "_" + label + ".pdf");

       can_fit->cd();
       legend_det->Draw();

       TF1* analytical_mean = new TF1("Analytical_mean", "( [0] + [1] * log( [2] + x) )", fitting_threshold_, 900.);
         analytical_mean->SetParameters(mean_alpha_nom/mean_alpha_denom, mean_beta_nom/mean_beta_denom, mean_gamma_nom/mean_gamma_denom); 
         analytical_mean->SetLineColor( TColor::GetColor("#FC00E7") );
 
       can_fit->SaveAs( folder_ + "/CalibrationFactors_meas_with_fit" + calib_ + "_" + label + ".C");             
       can_fit->SaveAs( folder_ +"/CalibrationFactors_meas_with_fit" + calib_ + "_" + label + ".pdf");

       drawoptions = "A";

       can_fit->Clear();
       can_meas->Clear();
       can_true->Clear();
       legend_det->Clear();

      TString calibration_tag = (set_of_tags_[ "calibration_tag" ])[ MC_files_[file_] ];

      int Ethresh_int_ = (int)round(Ethresh_);
      TString energy_threshold = TString::Format("%i", Ethresh_int_);
      
      if( (setup == "good_sectors" || setup == "all_sectors") && bin_phi == 1){
        if( calibration_tag == "MC" || calibration_tag == "data"){
          ofstream calibrating_values;
          calibrating_values.open( "Calibrating_values.h", ios::out | ios::app | ios::binary);

          calibrating_values << "double Energy_" << setup << "_" << calibration_tag << "_" << energy_threshold <<"( double edet ){" << endl;
          calibrating_values << "  double alpha = " << analytical->GetParameter(0) << ";" << endl;
          calibrating_values << "  double beta = "  << analytical->GetParameter(1) << ";" << endl;
          calibrating_values << "  double gamma = " << analytical->GetParameter(2) << ";" << endl;
          calibrating_values << "return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; \n}\n" << endl; 

	  calibrating_values.close();
        }
      }
    } // Loop over phi.
  } // Loop over files.
}



// Calculates calibration function for every single sector.
void Unfolder::CalibrationFunction_sectors(int first_sector, int last_sector, int which_file, TGraphErrors* &gre){


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

  double lowerbound_fit = Ethresh_;
  TF1* analytical = new TF1("Analytical", "  [0] + [1] * log( [2] + x)  ", lowerbound_fit, 900.);
    analytical->SetParLimits(2, 1.-lowerbound_fit, 1000.);


  // Create a map of canvasses: each canvas contains the plots from one file.
  map<TString, TCanvas*> canvasses;
  for(int file_ = 0; file_ < MC_files_.size(); file_++){
    TCanvas* can;
    PrepareCanvas( can, "CalibrationFactors_" + printLabel_[ MC_files_[file_] ] );

    canvasses[ printLabel_[ MC_files_[file_] ] ] = can;
  }
  TString drawoptions_sector ="A";
 
  /* Open files. */


  for(int sector = first_sector; sector <= last_sector; sector++){
    cout << "// -- Sector\t" << sector << endl;
    TString label;


    TCanvas *can_true = new TCanvas(TString::Format("can_true" + calib_ + "_sector_%i", sector), TString::Format("can_true" + calib_ + "_sector_%i", sector), 1000, 1000);
    can_true	->SetLeftMargin(0.25);
    can_true	->SetTopMargin(0.05);
    can_true	->SetBottomMargin(0.14);	 
	 
    TCanvas *can_meas = new TCanvas(TString::Format("can_meas" + calib_ + "_sector_%i", sector), TString::Format("can_meas" + calib_ + "_sector_%i", sector), 1000, 1000);
    can_meas	->SetLeftMargin(0.25);
    can_meas	->SetTopMargin(0.05);
    can_meas	->SetBottomMargin(0.14); 
      
    TCanvas *can_fit = new TCanvas(TString::Format("can_fit" + calib_ + "_sector_%i", sector), TString::Format("can_fit" + calib_ + "_sector_%i", sector), 1000, 1000);
    can_fit	->SetLeftMargin(0.25);
    can_fit	->SetTopMargin(0.05);
    can_fit	->SetBottomMargin(0.14);         
	 
    TCanvas *can_aver = new TCanvas(TString::Format("can_aver" + calib_ + "_sector_%i", sector), TString::Format("can_aver" + calib_ + "_sector_%i", sector), 1000, 1000);
    can_aver	->SetLeftMargin(0.25);
    can_aver	->SetRightMargin(0.05);
    can_aver	->SetBottomMargin(0.14);
      
    TLine *line = new TLine(0.,1.,1733.,1.);
    line->SetLineWidth(2);


    label = TString::Format( "Systematics_comparison_sector_%i", sector);

    int	new_dir = mkdir( TString::Format( folder_), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  );

    // -- Determine which files.
    int first_file, last_file;

    if( which_file = -1){
      first_file = 0;
      last_file = MC_files_.size();
    }
    else{
      first_file = which_file;
      last_file = which_file;
    }


    for(int file_ = 0; file_ < MC_files_.size(); file_++){

      TFile *_file0 = TFile::Open( MC_files_[file_], "Read");
       
      TCanvas *c = new TCanvas("can", "can",  1000, 1000);
        c->SetLeftMargin(0.20);
        c->SetRightMargin(0.18);
        c->SetBottomMargin(0.20);  

      TH2D* hGenE_selection;
      TH2D* hDetE_selection;
      TH2D* hResponse_selection;
      TH2D* hDetE;
      TH2D* hGenE;
      TH1D* hGen;
      TString setup = TString::Format("%i", sector);

       // -- Extract the response, Edet and Egen distributions from the files.
      Extract_2D_energy_distributions(file_, hGenE_selection, hDetE_selection, hResponse_selection, hDetE, hGenE, hGen, setup);

      TGraphErrors * Response_true;
      TGraphErrors * Response_meas;

      Analyze_response( hResponse_selection, hGenE_selection, hDetE_selection, hGenE, label, hGen, file_, Response_meas, Response_true);

      Response_meas->SetName( TString::Format("%s_" + printLabel_[ MC_files_[file_] ], Response_meas->GetName()) );
      Response_true->SetName( TString::Format("%c_" + printLabel_[ MC_files_[file_] ], Response_true->GetName()) );
      Response_meas->SetLineWidth(2);	
	  	
      // cout << "Names\t" << Response_meas->GetName() << "\t\t" << Response_true->GetName() << endl;

      legend_det->AddEntry( Response_meas, TString::Format(legend_info_[ MC_files_[file_] ] + " (Sector %i)", sector), "lp"); 
 
      // Draw.
      can_true->cd();	 
      Response_true->Draw("p" + drawoptions);

//      Response_true->Fit( analytical_2, "", "", 100., 1000.);
      line->Draw();

      can_true->SaveAs(TString::Format( folder_ + "/CalibrationFactors_true_calib_" + label + ".C") );
      can_true->SaveAs(TString::Format( folder_ + "/CalibrationFactors_true_calib_" + label + ".pdf") );

      // Draw.
 
      can_fit->cd();
      Response_meas->Draw("p" + drawoptions);
     
      can_fit->SaveAs(TString::Format( folder_ + "/CalibrationFactors_fit_calib_" + label + ".C") );
      can_fit->SaveAs(TString::Format( folder_ + "/CalibrationFactors_fit_calib_" + label + ".pdf") );
	
      drawoptions="same";	

      can_meas->cd(); 
      Response_meas->Draw("ape");// + drawoptions);
      line->Draw();

      analytical->SetLineColor( getColor(file_+1) );

      TString fitoptions = "";
      if( fitnotdrawn_ ) fitoptions = "0";
      Response_meas->Fit( analytical, fitoptions, "", lowerbound_fit, 900.);

      can_meas->SaveAs(TString::Format( folder_ + "/CalibrationFactors_meas_calib_" + label + ".C") );
      can_meas->SaveAs(TString::Format( folder_ + "/CalibrationFactors_meas_calib_" + label + ".pdf") );

      std::vector<double> parameters_calibration;
        parameters_calibration.push_back( analytical->GetParameter(0) );
        parameters_calibration.push_back( analytical->GetParameter(1) );
        parameters_calibration.push_back( analytical->GetParameter(2) );

      calibration_parameters_[ MC_files_[file_] ] = parameters_calibration;

       // -- Draw the calibration factors on the file-separate canvas.
      (canvasses[ printLabel_[MC_files_[file_]]])->cd();

       TGraphErrors* gre_new = (TGraphErrors*)Response_meas->Clone( TString::Format("response_measured_%i", sector) );
       gre_new->SetMarkerColor( getColor( sector ) );
       gre_new->SetLineColor( getColor( sector ) );
       gre_new->Draw("ap" + drawoptions_sector);



       (canvasses[ printLabel_[MC_files_[file_]]])->SaveAs("AndWhatDoesThisDo_" + printLabel_[MC_files_[file_]] + ".C");


      cout << "XXX---XXX\t" << Response_meas->GetName() << endl;

      gre = (TGraphErrors*)Response_meas->Clone( TString::Format("response_measured_%i", sector) );


      TString calibration_tag = (set_of_tags_[ "calibration_tag" ])[ MC_files_[file_] ];

      cout << "Calibration_tag is\t" << calibration_tag << endl;

      int Ethresh_int_ = (int)round(Ethresh_);
      TString energy_threshold = TString::Format("%i", Ethresh_int_);
      
      if( calibration_tag == "MC" || calibration_tag == "data"){
        ofstream calibrating_values;
        calibrating_values.open( "Calibrating_values.h", ios::out | ios::app | ios::binary);

        calibrating_values << "double Energy_sector_" << sector - 1 << "_" << calibration_tag << "_" << energy_threshold <<"( double edet ){" << endl;
        calibrating_values << "  double alpha = " << analytical->GetParameter(0) << ";" << endl;
        calibrating_values << "  double beta = "  << analytical->GetParameter(1) << ";" << endl;
        calibrating_values << "  double gamma = " << analytical->GetParameter(2) << ";" << endl;
        calibrating_values << "return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; \n}\n" << endl; 

	calibrating_values.close();
      }

       drawoptions_sector = "same";
    } // Loop over files.

    cout << "Finished all files" << endl;

    // Save.
    can_true->cd();
    legend_det->Draw();
    can_true->SaveAs( TString::Format( folder_ + "/CalibrationFactors_true_" + calib_ + "_" + label + "_%iGeV.C", static_cast<int>(Ethresh_) ) );		
    can_true->SaveAs( TString::Format( folder_ + "/CalibrationFactors_true_" + calib_ + "_" + label + "_%iGeV.pdf", static_cast<int>(Ethresh_) ) );
     
    can_meas->cd();
    legend_det->Draw();
    can_meas->SaveAs( TString::Format( folder_ + "/CalibrationFactors_meas_" + calib_ + "_" + label + "_%iGeV.C", static_cast<int>(Ethresh_) )); 		
    can_meas->SaveAs( TString::Format( folder_ + "/CalibrationFactors_meas_" + calib_ + "_" + label + "_%iGeV.pdf", static_cast<int>(Ethresh_) ));

    can_fit->cd();
    legend_det->Draw();

    can_fit->SaveAs( TString::Format( folder_ + "/CalibrationFactors_meas_with_fit" + calib_ + "_" + label + "_%iGeV.C", static_cast<int>(Ethresh_) ));             
    can_fit->SaveAs( TString::Format( folder_ + "/CalibrationFactors_meas_with_fit" + calib_ + "_" + label + "_%iGeV.pdf", static_cast<int>(Ethresh_) ) );

    drawoptions = "A";

    can_fit->Clear();
    can_meas->Clear();
    can_true->Clear();
    legend_det->Clear();
  
    cout << "YYY---YYY\t" << gre->GetName() << endl;

  } // Loop over sectors sects.
    cout << "ZZZ---ZZZ\t" << gre->GetName() << endl;
//  calibrating_values.close();
}

/************************************************************
* Code to extract a calibration function from the MC files. *
************************************************************/

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
//    if( sector == 11 ){ CalibrationFunction_sectors(1,16, graph_sector); } 
//    else{ CalibrationFunction_sectors(sector, sector, graph_sector); } 

    can->cd();

    if( sector != 11){ 
      graph_sector->SetMarkerColor( getColor( sector ) ); 
      graph_sector->SetLineColor( getColor( sector ));
    }

    graph_sector->Draw(drawoptions);
    drawoptions = "psame";

    TF1* analytical = new TF1("Analytical", " ( [0] + [1] * log( [2] + x) )", 100., 900.);
    analytical->SetParLimits(2, -99., 1000.);   

    analytical->SetLineColor( graph_sector->GetLineColor() );
    graph_sector->Fit( analytical, "0", "", fitting_threshold_, 900.);  

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














// Calculates systematics for all good sectors in one function.
void Unfolder::CalculateSystematics(TString setup, int first_sector, int last_sector){
  // Start by executing the calibration determination to obtain the necessary parameters.

  TCanvas *can;		// -- Needed for the functions.
  TCanvas *can_graph;	// -- Needed for the graphs.

  fitting_threshold_ = 0.;

  TString plots_label;
  TString legend_info;

  setup_calibration_ = setup;


  if( setup != "separate_sectors" ){
    first_sector = 1;
    last_sector = 1;
  }
  for( int sector = first_sector; sector <= last_sector; sector++){
    TGraphErrors* gre_meas;

    if( setup == "all_sectors"){
      plots_label = "all_sectors";
      legend_info = "All sectors";
    }
    else if( setup == "separate_sectors"){
      plots_label = TString::Format("sector_%i", sector);
      legend_info = TString::Format("Sector %i", sector);
    }
    else{
      plots_label = "good_sectors";
      legend_info = "Good sectors";
    }

    // Calculate the systematics by getting the behaviour of all well behaved sectors in a single function.
    CalibrationFunction(sector,sector, gre_meas); 

    // Create a canvas to plot.
    PrepareCanvas( can, TString::Format("Systematics_comparison_%i", sector) );
    PrepareCanvas( can_graph, TString::Format("Systematics_graphs_%i", sector) );

    TPad *pad1, *pad2;
    SplitCanvas(can, pad1, pad2);

    TLegend *legend_det = new TLegend(0.25, 0.45, 0.75, 0.95);
      legend_det->SetFillStyle( 0 );
      legend_det->SetBorderSize( 0 );

   TLegend* legend_graph = new TLegend(0.25, 0.60, 0.75, 0.90);
      legend_graph->SetFillStyle( 0 );
      legend_graph->SetBorderSize( 0 );

    // 1. Loop over files.
    // 2. Get parameters and use in function.
    // 3. Plot function as Graph.
    // 4. Profit.

    TString drawoptions = "";
    TString drawoptions_graph = "ape";
    TH1F* original_calibration;

    // Needed for the range of the y-axis.
    double min_val, max_val;

    cout << "Before loop over files" << endl;
    for(int file = 0; file < MC_files_.size(); file++){

      can->cd();

      // -- Extract the calibration function's parameters.
      double alpha = (calibration_parameters_[ MC_files_[file] ])[0];
      double beta = (calibration_parameters_[ MC_files_[file] ])[1];
      double gamma = (calibration_parameters_[ MC_files_[file] ])[2];
 
      TString calibration_tag = ""; 
      cout << "QQQ---\t" << calibration_tag << "\t" << setup << endl;

      if( (calibration_tag == "MC" || calibration_tag == "data") && (setup == "good_sectors" || setup == "all_sectors" )  ){
        ofstream calibrating_values;
        calibrating_values.open( "Calibrating_" + setup + "_values.h", ios::out | ios::app | ios::binary);

        calibrating_values << "double Energy_" << setup << "_" << calibration_tag << "( double edet ){" << endl;
        calibrating_values << "  double alpha = " << alpha << ";" << endl;
        calibrating_values << "  double beta = "  << beta << ";" << endl;
        calibrating_values << "  double gamma = " << gamma << ";" << endl;
        calibrating_values << "return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; \n}\n" << endl; 

	calibrating_values.close();
      }


      // -- Calculate and draw the calibration function.
      TF1* calibration_function = new TF1(TString::Format("calibration_function_sector_%i_" + printLabel_[ MC_files_[file] ], sector), "[0] + [1] * log( [2] + x)", fitting_threshold_, 900.);
      calibration_function->SetParameters( alpha, beta, gamma );
      TH1* calibration_histogram = calibration_function->GetHistogram();

      // -- Prepare the first MC sample to serve as a reference.
      if( file == 0 ){
        original_calibration = (TH1F*)calibration_histogram->Clone("Original_calibration_function");
        original_calibration->GetYaxis()->SetRangeUser(0., 2.);
      }

      legend_det->AddEntry( calibration_histogram, legend_info_[ MC_files_[file]] + " (" + legend_info + ")", "l");
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

      pad1->cd();
      legend_det->Draw();

      drawoptions = "same";

      can->SaveAs(TString::Format( folder_ + "/Systematics_function_" + plots_label + "_%iGeV.C", static_cast<int> (Ethresh_ ) ) );
      can->SaveAs(TString::Format( folder_ + "/Systematics_function_" + plots_label + "_%iGeV.pdf", static_cast<int> (Ethresh_ ) ) );

      can_graph->cd();
cout << "ZZZ\tto extract\t" << MC_files_[file] << endl;
      TGraphErrors* graph_errs = calibration_graphs_[ MC_files_[file] ];
cout << "ZZZ\textracted\t" << MC_files_[file] << endl;
      graph_errs->Draw(drawoptions_graph);
      drawoptions_graph = "pesame";

      cout << "Ethresh_\t" << Ethresh_ << "\t" << static_cast<int> (Ethresh_ ) << endl;

      legend_graph->AddEntry( graph_errs, legend_info_[ MC_files_[file]] + " (" + legend_info + ")", "p");

      legend_graph->Draw();
      
      can_graph->SaveAs(TString::Format( folder_ + "/Systematics_graph_" + plots_label + "_%iGeV.C", static_cast<int> (Ethresh_ ) ) );
      can_graph->SaveAs(TString::Format( folder_ + "/Systematics_graph_" + plots_label + "_%iGeV.pdf", static_cast<int> (Ethresh_ ) ) );
    } // Loop over files.
    
    calibration_parameters_.clear();
  } // Loop over sectors.
}













void Unfolder::CalibrationFactors_oneCanvas(bool draw_functions){
  // We do not wish to draw fits in intermediate functions.
  SetFitDraw( false );

  // Start by executing the calibration determination to obtain the necessary parameters.

  TCanvas *can;

  for(int file = 0; file < MC_files_.size(); file++){

      // Create a canvas to plot.

      PrepareCanvas( can, TString::Format("CalibrationFactors_oneCanvas_" + printLabel_[ MC_files_[file] ] ) );

      TLegend *legend_det = new TLegend(0.65, 0.58, 0.95, 0.93);
        legend_det->SetFillStyle( 0 );
        legend_det->SetBorderSize( 0 );

      TString drawoptions_gre = "APE";
      TString drawoptions = "";

    for( int sector = 1; sector <= 16; sector++){
      TGraphErrors* gre_meas;
      CalibrationFunction_sectors(sector,sector, file, gre_meas); 

      cout << "!!!---!!!\t" << gre_meas->GetName() << endl;

      // 1. Loop over files.
      // 2. Get parameters and use in function.
      // 3. Plot function as Graph.
      // 4. Profit.

      // Needed for the range of the y-axis.
      double min_val, max_val;
      can->cd();
      double alpha = (calibration_parameters_[ MC_files_[file] ])[0];
      double beta = (calibration_parameters_[ MC_files_[file] ])[1];
      double gamma = (calibration_parameters_[ MC_files_[file] ])[2];
  
      TF1* calibration_function = new TF1(TString::Format("calibration_function_sector_%i_" + printLabel_[ MC_files_[file] ], sector), "[0] + [1] * log( [2] + x)", fitting_threshold_, 900.);
      calibration_function->SetParameters( alpha, beta, gamma );
      TH1* calibration_histogram_orig = calibration_function->GetHistogram();
      TH1D* calibration_histogram = (TH1D*)calibration_histogram_orig->Clone(TString::Format("Calibration_values_sector_%i", sector) );
      
      legend_det->AddEntry( calibration_histogram, TString::Format( legend_info_[ MC_files_[file]] + " (sector %i)", sector), "l");

      can->cd();

      calibration_histogram->SetLineColor( getColor( sector ) );
      calibration_histogram->GetXaxis()->SetTitle("E_{det}");
      calibration_histogram->GetYaxis()->SetRangeUser(0., 4.);
      calibration_histogram->GetYaxis()->SetTitle("<#frac{E_{gen}}{E_{det}}>");

      if( draw_functions ){
        calibration_histogram->DrawCopy("hist" + drawoptions);
      }

      gre_meas->SetMarkerColor( getColor(sector) );
      gre_meas->SetLineColor( getColor(sector) );
      gre_meas->GetXaxis()->SetTitle("E_{det}");
      gre_meas->GetYaxis()->SetRangeUser(0., 5.5);
      gre_meas->GetYaxis()->SetTitle("<#frac{E_{gen}}{E_{det}}>");

      gre_meas->Draw( drawoptions_gre );
      drawoptions_gre = "pesame";
      drawoptions = "same";

    } // Loop over sectors.

      legend_det->Draw();

    TString with_fit = "_fit";
    if( !draw_functions ){ with_fit = "_no_fit"; }
    TString savenameC = "Plots/Systematics_Test/CalibrationFactors_" + printLabel_[ MC_files_[file] ] + with_fit + ".C";
    TString savenamepdf = "Plots/Systematics_Test/CalibrationFactors_" + printLabel_[ MC_files_[file] ] + with_fit + ".pdf";

    can->SaveAs( savenameC);
    can->SaveAs( savenamepdf );
    calibration_parameters_.clear();
  } // Loop over files.
  SetFitDraw( true );
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


/*******************************************************************
* Take 2D histograms and extract from them the calibration values. *
*******************************************************************/


void Unfolder::Analyze_response(TH2D* hResponse_selection, TH2D* hGenE_selection, TH2D* hDetE_selection, TH2D* hGenE, TString label, TH1D* hGen, int file, TGraphErrors* &gre_meas, TGraphErrors* &gre_true){

     cout << "ANALYZE_RESPONSE\t" << file << endl;
     cout << hResponse_selection->Integral() << "\t" << hGenE_selection->Integral() << "\t" << hDetE_selection->Integral() << "\t" << hGenE->Integral() << "\t" << endl;
     int event_threshold_ = 25;
//     int event_threshold_ = 50;
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
         // cout << "\t@@@Ebin " << bin_E << ":\t" << hGenE_1D->Integral() << "\t" << hDetE_1D->Integral() << "\t" << hResponse_1D->Integral() << endl;
	 continue;
       }

       // Our current histogram, combined with the stored histogram, has enough statistics. Add them and go on.
       else if( ( hResponse_1D->Integral() + hSlice_storage->Integral() ) >= event_threshold_ ){
         hResponse_1D ->Add( hSlice_storage );
         hGenE_1D     ->Add( hEgen_storage );
	 hDetE_1D	->Add( hEdet_storage );
         // cout << "\t@@@Before fitting:\t" << hGenE_1D->Integral() << "\t" << hDetE_1D->Integral() << "\t" << hResponse_1D->Integral() << endl;
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
	 if( fit_response->GetParameter(1) != 0. && fit_response->GetParameter(1) == fit_response->GetParameter(1) ){
	   double R = fit_response->GetParameter(1);
	   double sR= fit_response->GetParError(1);

	   cout << "\t€€€ Averages\t" << GetAverage(hDetE_1D) << "\t" << GetAverage(hGenE_1D) << endl;

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
     
     Response_true	->SetMarkerColor( getColor(file + 1) );
     Response_true	->SetLineColor( getColor(file + 1) );
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


/****************************
* Plot the response matrix. *
****************************/


void Unfolder::PlotResponseMatrix(int file_ ){

  TCanvas *can;
  PrepareCanvas_2D(can, "Response_matrix_" + printLabel_[ MC_files_[file_] ] );

  TH2D* hRes;
  Get_ResponseMatrix( file_, hRes );
  Prepare_2Dplot( hRes );

  hRes->Draw("colz");

  can->SetLogz();

  can->SaveAs( TString::Format("Plots/ResponseMatrix/Response_matrix_" + label_ + "_" + printLabel_[ MC_files_[file_] ] + ".C") );
  can->SaveAs( TString::Format("Plots/ResponseMatrix/Response_matrix_" + label_ + "_" + printLabel_[ MC_files_[file_] ] + ".pdf") );

  cout << "€€€\t" << hRes->Integral() << "\tevents" << endl;
}



void Unfolder::PlotFromResponseMatrix(int file_, TString axis_1, TString axis_2){

  int axis_1_, axis_2_;
  cout << "Prepare axes" << endl;

  axis_1_ = SetAxisTHnSparse( axis_1 );
  axis_2_ = SetAxisTHnSparse( axis_2 );
  
  cout << "Axes\t" << axis_1_ << "\t" << axis_2_ << endl;

  THnSparseD* hSparse;
  
  cout << "Prepared";
  Get_ResponseMatrixTHn(file_, hSparse);
  cout << "\tdone" << endl;

  TH2D* hResponseMatrix_selection  = (TH2D*)hSparse->Projection(axis_2_, axis_1_); 

  TCanvas *can;
  PrepareCanvas_2D(can, "Plot_" + axis_1 + "_vs_" + axis_2 + "_" + printLabel_[ MC_files_[file_] ]);
  Prepare_2Dplot( hResponseMatrix_selection );

  hResponseMatrix_selection->Draw("colz");

  can->SetLogz();
  can->Update();
  


  can->SaveAs("Plots/ResponseMatrix_extended/" +  axis_1 + "_vs_" + axis_2 + "_" + label_ + "_" + printLabel_[ MC_files_[file_] ] + ".C" );
  can->SaveAs("Plots/ResponseMatrix_extended/" +  axis_1 + "_vs_" + axis_2 + "_" + label_ + "_" + printLabel_[ MC_files_[file_] ] + ".pdf" );

}



void Unfolder::PlotFromResponseMatrix(int file_, TString axis_1){

  int axis_1_;
  cout << "Prepare axes" << endl;

  axis_1_ = SetAxisTHnSparse( axis_1 );
  cout << "Axis\t" << axis_1_ << endl;

  THnSparseD* hSparse;
  
  cout << "Prepared";
  Get_ResponseMatrixTHn(file_, hSparse);
  cout << "\tdone" << endl;

  TH1D* hResponseMatrix_selection  = (TH1D*)hSparse->Projection(axis_1_); 

  TCanvas *can;
  PrepareCanvas(can, "Plot_" + axis_1 + "_" + printLabel_[ MC_files_[file_] ]);

  hResponseMatrix_selection->Draw("hist");

//  can->SetLogy();
  can->Update();
  
  can->SaveAs("Plots/ResponseMatrix_extended/" +  axis_1 + "_" + label_ + "_" + printLabel_[ MC_files_[file_] ] + ".C" );
  can->SaveAs("Plots/ResponseMatrix_extended/" +  axis_1 + "_" + label_ + "_" + printLabel_[ MC_files_[file_] ] + ".pdf" );

}
int Unfolder::SetAxisTHnSparse(TString axis){

  if(axis == "Edet")  
	{ return 0; }
  else if( axis == "Egen")
	{ return 1; }
  else if( axis == "nMatch")
	{ return 2; }
  else if( axis == "nFake")
	{ return 3; }
  else if( axis == "nMiss")
	{ return 4; }
  else if( axis == "nDet")
	{ return 5; }
  else if( axis == "nGen")
	{ return 6;}

  else{ return -1; }
}


void Unfolder::Dissect_ResponseObject(){

  for(int file_ = 0; file_ < MC_files_.size(); file_++){
    TString file_label_ = printLabel_[ MC_files_[file_] ];

    TH1D* hMiss, *hFake;

    TCanvas *can_fake;
    PrepareCanvas( can_fake, TString::Format("Fakes_" + file_label_) );
    Hist_getFake(file_, hFake);    
    hFake->Draw("hist");
    can_fake->SaveAs("Plots/ResponseMatrix_extended/Fakes_" + file_label_ + ".C");
    can_fake->SaveAs("Plots/ResponseMatrix_extended/Fakes_" + file_label_ + ".pdf");
 /* 
    TCanvas *can_miss;
    PrepareCanvas( can_miss, TString::Format("Misses_" + file_label_) );
    Hist_getMiss(file_, hMiss);    
    hMiss->Draw("hist");
    can_fake->SaveAs("Plots/ResponseMatrix/Misses_" + file_label_ + ".C");
    can_fake->SaveAs("Plots/ResponseMatrix/Misses_" + file_label_ + ".pdf");

    TCanvas *can_resp;
    PrepareCanvas( can_resp, TString::Format("Response_" + file_label_) );
*/
  }
}



double Unfolder::Calculate_chi2(TH1D* hist_reference, TH1D* hist_result  ){

  double sum_ref = 0., sum_res = 0.;
  
  for(int i = 0; i <= hist_reference->GetNbinsX(); i++){
    sum_ref += hist_reference->GetBinContent( i );
  }

  double chi2 = 0.;
  for(int i = 0; i <= hist_reference->GetNbinsX(); i++){
    double bin_ref = hist_reference->GetBinContent( i );
    double bin_res = hist_reference->GetBinContent( i );

    double err_ref = sqrt(bin_ref);
    double err_res = hist_reference->GetBinError( i );

    double nom = pow(  ( sum_res * bin_ref - sum_ref * bin_res), 2. );
    double denom = pow(sum_ref * err_res, 2.) + pow(sum_res*err_ref, 2.);

    chi2 += nom/denom;
  }
  
  return chi2;
}





int Unfolder::Extract_2D_energy_distributions(int file_, TH2D* &hGenE_selection, TH2D* &hDetE_selection, TH2D* &hResponse_selection, TH2D* &hDetE, TH2D* &hGenE, TH1D* &hGen, TString setup){
      cout << "\%\%\%" << "Extract 2D\t" << setup <<  endl;

      TFile* _file0 = TFile::Open(MC_files_[file_], "read");
      int rebinner = 1;
      int first_sector, last_sector;
      vector<int> bad_sectors;

      if( setup == "all_sectors" ){ 
	first_sector = 1;
	last_sector = 16;
      }
      else if( setup == "good_sectors" ){
	first_sector = 1;
	last_sector = 16;
	bad_sectors.push_back(13);
	bad_sectors.push_back(14);
      }
      else if( setup.IsDigit() ){
	if( setup.Atoi() >= 1 && setup.Atoi() <= 16 ){
	  first_sector = setup.Atoi();
	  last_sector = setup.Atoi();
	}
	else return 0;

      }

       // -- Extract the response, Edet and Egen distributions from the files.
  
      THnSparseD* hResponse 	= (THnSparseD*)_file0->Get("hResponse_gen_phi");			//hResponse->Draw("colz");
      THnSparse* hGenE_fine_phi   = (THnSparse*)_file0->Get("hGen_fine_phi");
      THnSparse* hDetE_fine_phi = (THnSparse*)_file0->Get("hDet_fine_phi");
      
      cout << "EXTRACT\t" << hResponse->GetEntries() << "\t" << hGenE_fine_phi->GetEntries() << "\t" << hDetE_fine_phi->GetEntries() << endl;

      hDetE 	= (TH2D*)_file0->Get("hCastorJet_energy_response");//	hDetE->Draw("colz");	
        hDetE->RebinX( rebinner );
        hDetE->RebinY( rebinner );

      hGenE	= (TH2D*)_file0->Get("hGen_fine");//			hGenE->Draw("colz");	
        hGenE->RebinX( rebinner );*
        hGenE->RebinY( rebinner );
       
      hGen	= (TH1D*)_file0->Get("hGenJet_energy");		

      // Templates for 1D projections of THnSparse.
      TH1D* hSlice_storage;
      TH1D* hEgen_storage;
      TH1D* hEdet_storage;

      // Prepare the histograms;
      // 1. Project the 3D histogram into its desired axes.
      // 2. Rebin the histogram to its desired number of bins per axis.
      // 3. Reset the contents etc.
      hGenE_selection  = (TH2D*)hGenE_fine_phi->Projection(1,0);
         hGenE_selection->RebinX( rebinner );
         hGenE_selection->RebinY( rebinner );
         hGenE_selection->Reset();
      hDetE_selection  = (TH2D*)hDetE_fine_phi->Projection(1,0);
         hDetE_selection->RebinX( rebinner );
         hDetE_selection->RebinY( rebinner );
         hDetE_selection->Reset();
      hResponse_selection = (TH2D*)hResponse->Projection(1,0);
         hResponse_selection->RebinX( rebinner );
         hResponse_selection->RebinY( rebinner );
         hResponse_selection->Reset();

      for(int bin_phi = first_sector; bin_phi <= last_sector; bin_phi++){
	cout << "Sector\t" << bin_phi;
	if( setup == "good_sectors" ){
	  if( binary_search(bad_sectors.begin(), bad_sectors.end(), bin_phi ) ){
	    cout << "\tbad" << endl;
	    continue;
	  }
	  cout << "Sector\t" << bin_phi << "\tgood";
	}
	cout << endl;

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
  cout << "Done extractin'" << endl;

  return 0;
}



void Unfolder::Calculate_smearedBackError(int file_, int MC_, int iterations, TH1D* &hActual_hist, TH1D* hUnfold){



  cout << "Unfolder::Calculate_smearedBackError" << endl;
  int nPoisson = 10000;

  //-- Procedure.
  /* 
  1. Extract DET distribution from file_; RooUnfoldResponse, misses and fakes from MC_.
  2. Copy 
  */
  
  // -- Open files.
  TFile* _file_det, *_file_unfold;
  if( file_ < MC_files_.size() ){	_file_det = new TFile( MC_files_[file_], "Read");}
  else if( file_ == -1 ){		_file_det = new TFile( datafile_, "Read");	}

  _file_unfold = new TFile( MC_files_[MC_], "Read"); 

  // (1) Extract the histograms.
  // cout << "Unfolder::Calculate_smearedBackError (1)" << endl;

  TH1D* hDet = (TH1D*)_file_det->Get("hCastorJet_energy");	//hDet->Scale( renorm_ );

  TH1D* hMiss = (TH1D*)_file_unfold->Get("hCastorJet_miss_all");	hMiss->Scale( 1./renorm_ );
  TH1D* hFake = (TH1D*)_file_unfold->Get("hCastorJet_fake_all");	hFake->Scale( 1./renorm_ );
  RooUnfoldResponse* response 	= (RooUnfoldResponse*)_file_unfold->Get("response");
    TH2D* hResponse 	= (TH2D*)response->Hresponse(); 		hResponse->Scale( 1./renorm_ );
    TH1D* hTruth 	= (TH1D*)response->Htruth();			hTruth->Scale( 1./renorm_ );
    TH1D* hMeasured	= (TH1D*)response->Hmeasured();			hMeasured->Scale( 1./renorm_ );

  // (2) Create empty copies.
  // cout << "Unfolder::Calculate_smearedBackError (2)" << endl;

  TH1D* hSmear = (TH1D*)hDet->Clone("hSmear");
  TH1D* hDet_new = (TH1D*)hDet->Clone("hDet_new");			
  TH1D* hMiss_new = (TH1D*)hMiss->Clone("hMiss_new");			
  TH1D* hFake_new = (TH1D*)hFake->Clone("hFake_new");
  TH1D* hTruth_new = (TH1D*)hTruth->Clone("hTruth_new");		
  TH1D* hMeasured_new = (TH1D*)hMeasured->Clone("hMeasured_new");		
  TH2D* hResponse_new = (TH2D*)hResponse->Clone("hResponse_new");	

  // (2b) Create THnSparse to store unfold and smeared distributions.
  int bins_sparseUnf[2] = { hTruth->GetNbinsX(), nPoisson };
  double mins_sparseUnf[2] = {0., 0.};
  double maxs_sparseUnf[2] = {hTruth->GetBinLowEdge( 1 ), nPoisson};
  THnSparseD sparseUnf("Sparse_unf", "Sparse_unf", 2, bins_sparseUnf, mins_sparseUnf, maxs_sparseUnf);

  int bins_sparseSm[2] = { hMeasured->GetNbinsX(), nPoisson };
  double mins_sparseSm[2] = {0., 0.};
  double maxs_sparseSm[2] = {hMeasured->GetBinLowEdge( 1 ), nPoisson};
  THnSparseD sparseSm("Sparse_sm", "Sparse_sm", 2, bins_sparseSm, mins_sparseSm, maxs_sparseSm);


//  TH1D* hUnfold = (TH1D*)hMiss->Clone("hUnfold");			

  // cout << "Unfolder::Calculate_smearedBackError - Clone spread" << endl;
  int binsx = hDet_new->GetNbinsX();
  double minx = hDet_new->GetXaxis()->GetBinLowEdge( 1 ), maxx = hDet_new->GetXaxis()->GetBinUpEdge( binsx );
  int binsy = 1000;
  double miny = 0., maxy = 250000.;

  TH2D* hSmear_spread = new TH2D("hSmear_spread", "Smear spread", binsx, minx, maxx, binsy, miny, maxy);
    hSmear_spread->Reset();
//    hSmear_spread->GetYaxis()->Set( 1000, 0., 100000. );
    hSmear_spread->GetYaxis()->SetTitle("Distribution of jets in this E-bin");
  
  TRandom3* rand = new TRandom3();

  // -- Repeat the following algorithm N times.
  for(int n_spread = 0; n_spread < nPoisson; n_spread++){
    if( n_spread%100 == 0){  cout << "Iterations\t" << iterations << "\tIteration\t" << n_spread << endl; }

    hDet_new->Reset();
    hMiss_new->Reset();
    hFake_new->Reset();
    hTruth_new->Reset();
    hMeasured_new->Reset();
    hResponse_new->Reset(); 
    hUnfold->Reset();

    // (3) Fill with Poissonian distribution.
    // cout << "Unfolder::Calculate_smearedBackError (3)" << endl;

    FillAnew_1D( hDet, hDet_new, rand);
    FillAnew_1D( hFake, hFake_new, rand);
    FillAnew_2D( hResponse, hResponse_new, rand);

    // (4) Unfold-and-smear.
    //cout << "Unfolder::Calculate_smearedBackError (4)" << endl; 
    RooUnfoldResponse *response_new = new RooUnfoldResponse( hMeasured, hTruth, hResponse_new );

    // -- S = R * U + F		-- Smearing back applied to all generator level jets.

    hSmear = (TH1D*) response_new->ApplyToTruth( hUnfold );

    hSmear->Add( hFake );

    for(int bin_smear = 0; bin_smear <= hSmear_spread->GetNbinsX(); bin_smear++){
      double smear_center = hSmear_spread->GetXaxis()->GetBinCenter( bin_smear );   
      double smear_value = hSmear->GetBinContent( bin_smear );
      hSmear_spread->Fill( smear_center, smear_value );
    }

  } // Loop over algorithm.

  // (5) Extract the spread from the distribution.
  // -- Loop over bins.
  for(int bin_smear = 0; bin_smear <= hSmear_spread->GetNbinsX(); bin_smear++){
    
    TH1D* hCurrent_smear_bin = (TH1D*)hSmear_spread->ProjectionX(TString::Format("hSmear_spread_1D_bin_%i", bin_smear), bin_smear, bin_smear, "do");
    TF1 *fit_response	= new TF1(TString::Format("Fit_Response_1D_bin_%i", bin_smear), "gaus");
    if( hCurrent_smear_bin->Integral() > 0 ){ 
      hCurrent_smear_bin->Fit(fit_response);
      double spread = fit_response->GetParameter(2);
      hActual_hist->SetBinError( bin_smear, spread );
    }
    else{ hActual_hist->SetBinError( bin_smear, 0. ); }    
  }


  TCanvas* can = new TCanvas(TString::Format("canvas_spread_%i", iterations), TString::Format("canvas_spread_%i", iterations), 1.);
  hSmear_spread->Draw("colz");
  can->SaveAs("Plotted_smear.pdf");
  can->SaveAs("Plotted_smear.C");
}


void Unfolder::FillAnew_1D(TH1D* hOld, TH1D* &hNew, TRandom3* &rand){
   

  for(int bin = 1; bin <= hOld->GetNbinsX(); bin++){
    // Measured distribution.
    int nevents = hOld->GetBinContent( bin );
    int nevents_new;

    if( nevents < 1000 ){ 
      nevents_new = rand->Poisson(  nevents );
    }
    else{ 
      nevents_new = rand->Gaus(  nevents,  sqrt(nevents) );
      //if( nevents_new < 0) nevents_new = -1. * nevents_new;
      while( nevents_new < 0. ) nevents_new = rand->Gaus(  nevents, sqrt(nevents));      
    }

    hNew->SetBinContent( bin, static_cast<double>(nevents_new) );


  } 
}



void Unfolder::FillAnew_1D(TH1D* hOld, TH1D* hMin, TH1D* &hNew, TRandom3* &rand){
   

  for(int bin = 1; bin <= hOld->GetNbinsX(); bin++){
    // Measured distribution.
    int nevents = hOld->GetBinContent( bin );
    int nevents_new;
    int nevents_min = hMin->GetBinContent( bin );

    if( nevents < 1000 ){ 
      nevents_new = rand->Poisson(  nevents );
      while( nevents_new < nevents_min ) nevents_new = rand->Poisson(  nevents );
    }
    else{ 
      nevents_new = rand->Gaus(  nevents,  sqrt(nevents) );
      //if( nevents_new < 0) nevents_new = -1. * nevents_new;
      while( nevents_new < nevents_min ) nevents_new = rand->Gaus(  nevents, sqrt(nevents));      
    }

    hNew->SetBinContent( bin, static_cast<double>(nevents_new) );


  } 
}


void Unfolder::FillAnew_2D(TH2D* hOld, TH2D* &hNew, TRandom3* &rand){
 
  ofstream testing_covariance, testing_cov2;
  testing_covariance.open("Testing_filling.txt", ios::app);  

  for(int binx = 1; binx <= hOld->GetNbinsX(); binx++){
    for(int biny = 1; biny <= hOld->GetNbinsY(); biny++){ 
      // Measured distribution.
      int nevents = hOld->GetBinContent( binx, biny );
      int nevents_new;

      if( nevents < 1000 ){ 
        nevents_new = rand->Poisson( nevents );
      }
      else{ 
        nevents_new = rand->Gaus( nevents, sqrt(nevents) );
        while( nevents_new < 0.  ) nevents_new = rand->Gaus(   nevents, sqrt(nevents) );
      }
      hNew->SetBinContent( binx, biny, static_cast<double>(nevents_new) );
    }
  }
  testing_covariance.close(); 
}





void Unfolder::SetCastorJetEnergy_norm( double renorm ){
  renorm_ = renorm;
}



void Unfolder::CovarianceMatrix( TH1D* hUnfold, TH1D* hSmeared, TMatrixD& cov_){

  cout << "Unfolder:CovarianceMatrix" << endl;
  int _nUnf = hUnfold->GetNbinsX();
  int _nSm = hSmeared->GetNbinsX();

  int _Nunf = hUnfold->Integral();
  int _Nsm = hSmeared->Integral();

  cout << "Unfolder:CovarianceMatrix\tResize matrix" << endl;
  cov_.ResizeTo(_nSm+1, _nUnf+1 );

  // Average in unfolded distribution.
  double _avUnf = 0.;
  for(int bin = 1; bin <= _nUnf; bin++){
    double bincenter = hUnfold->GetBinCenter( bin );
    double binvalue = hUnfold->GetBinContent( bin );
    _avUnf += bincenter*binvalue;
  }
  _avUnf = _avUnf/_Nunf;

  
  // Average in smeared distribution.
  double _avSm = 0.;
  for(int bin = 1; bin <= _nUnf; bin++){
    double bincenter = hSmeared->GetBinCenter( bin );
    double binvalue = hSmeared->GetBinContent( bin );
    _avSm += bincenter*binvalue;
  }
  _avSm = _avSm/_Nsm;
  

  for(int bin_unf = 1; bin_unf <= _nUnf; bin_unf++){
    for(int bin_sm = 1; bin_sm <= _nSm; bin_sm++){

      double val_unf = hUnfold->GetBinContent( bin_unf );
      double val_sm = hSmeared->GetBinContent( bin_sm );

      cov_(bin_sm, bin_unf) = (val_unf - _avUnf) * (val_sm - _avSm ) * val_unf/_Nunf * val_sm/_Nsm;
    }
  }  
}



double Unfolder::Calculate_smearedBackError_covariance(TH1D* hData, TH1D* hUnfold, RooUnfoldResponse* response, int iterations){

  cout << "Unfolder::Calculate_smearedBackError_covariance" << endl;
  int nPoisson = 100;
  int file_ = -1;
  int MC_ = 0;

  ofstream testing_covariance;
  testing_covariance.open("Testing_covariance.txt", ios::app);

  cout << "***\t" << iterations << "\t" << hUnfold->Integral() << endl;

  // -- Open files.
  TFile* _file_det, *_file_unfold;
  if( file_ < MC_files_.size() ){	_file_det = new TFile( MC_files_[file_], "Read");}
  else if( file_ == -1 ){		_file_det = new TFile( datafile_, "Read");	}

  _file_unfold = new TFile( MC_files_[MC_], "Read"); 

  // (1) Extract the histograms.
  // cout << "Unfolder::Calculate_smearedBackError (1)" << endl;

  TH1D* hMiss 		= (TH1D*)_file_unfold->Get("hCastorJet_miss_all");	hMiss->Scale( 1./renorm_ );	//GetSubHistogram( hMiss, hMiss, Ecut_, 2000);
  TH1D* hFake 		= (TH1D*)_file_unfold->Get("hCastorJet_fake_all");	hFake->Scale( 1./renorm_ );	//GetSubHistogram( hFake, hFake, Ecut_, 2000);
  TH2D* hResponse_proj 	= (TH2D*)response->Hresponse(); 			hResponse_proj->Scale( 1./renorm_ );	//GetSubHistogram( hResponse, hResponse, Ecut_, 2000);
  TH1D* hTruth 		= (TH1D*)response->Htruth();				hTruth->Scale( 1./renorm_ );	//GetSubHistogram( hTruth, hTruth, Ecut_, 2000);
  TH1D* hMeasured	= (TH1D*)response->Hmeasured();				hMeasured->Scale( 1./renorm_ );	//GetSubHistogram( hMeasured, hMeasured, Ecut_, 2000);


  TCanvas *can_resp; 
  PrepareCanvas_2D(can_resp, "Response");
  can_resp->SetLogz();

  TH2D* hResponse = (TH2D*)hResponse_proj->Clone("Response_2D");
  hResponse->SetName("Response_2D");
  hResponse->Draw("colz");
  can_resp->SaveAs(folder_ + "Response.pdf");
  can_resp->SaveAs(folder_ + "Response.C");


  cout << "$$$$Mis\t" << hMiss->Integral() << endl;

  TCanvas *canData = new TCanvas("canData", "canData", 1.);
  hData->Draw("ehist");

  TCanvas *canMiss = new TCanvas("canMiss", "canMiss", 1.);
  hMiss->Draw("ehist");

  TCanvas *canFake = new TCanvas("canFake", "canFake", 1.);
  hFake->Draw("ehist");

  TCanvas *canMeasured = new TCanvas("canMeasured", "canMeasured", 1.);
  hMeasured->Draw("ehist");

  TCanvas *canTruth = new TCanvas("canTruth", "canTruth", 1.);
  hTruth->Draw("ehist");
  
  // (2) Create empty copies.
  // cout << "Unfolder::Calculate_smearedBackError (2)" << endl;

  TH1D* hSmear = (TH1D*)hData->Clone("hSmear");
  TH1D* hData_new = (TH1D*)hData->Clone("hData_new");			
  TH1D* hMiss_new = (TH1D*)hMiss->Clone("hMiss_new");			
  TH1D* hFake_new = (TH1D*)hFake->Clone("hFake_new");
  TH1D* hTruth_new = (TH1D*)hTruth->Clone("hTruth_new");		
  TH1D* hMeasured_new = (TH1D*)hMeasured->Clone("hMeasured_new");		
  TH2D* hResponse_new = (TH2D*)hResponse->Clone("hResponse_new");	

  // (2b) Create THnSparse to store data and smeared distributions.
  int bins_sparseData[2] = { hData->GetNbinsX(), nPoisson };
  double mins_sparseData[2] = {hData->GetXaxis()->GetBinLowEdge(1), 0.};
  double maxs_sparseData[2] = {hData->GetXaxis()->GetBinUpEdge( hData->GetNbinsX() ), nPoisson};
  THnSparseD sparseData("Sparse_data", "Sparse_data;E_{data};n_{poisson};", 2, bins_sparseData, mins_sparseData, maxs_sparseData);

  int bins_sparseSm[2] = { hMeasured->GetNbinsX(), nPoisson };
  double mins_sparseSm[2] = { hMeasured->GetXaxis()->GetBinLowEdge(1), 0.};
  double maxs_sparseSm[2] = {hMeasured->GetXaxis()->GetBinUpEdge( hMeasured->GetNbinsX() ), nPoisson};
  THnSparseD sparseSmear("Sparse_smear", "Sparse_smear;E_{smear};n_{poisson};", 2, bins_sparseSm, mins_sparseSm, maxs_sparseSm);
  THnSparseD sparseSmear_sq("Sparse_smear_squared", "(E-#mu)^2;E_{smear};n_{poisson};", 2, bins_sparseSm, mins_sparseSm, maxs_sparseSm);

  cout << "Unfolder::Calculate_smearedBackError_covariance - Clone spread" << endl;
 
  TRandom3* rand = new TRandom3();
    rand->SetSeed( iterations  );

    //-- Get a time dependent seed.
    time_t timer;
    struct tm y2k = {0};
    double seconds;

    y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
    y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;    
    seconds = difftime(timer,mktime(&y2k));

    rand->SetSeed( seconds );

  //----------------------------------------------//
  // -- Repeat the following algorithm N times. --//
  //----------------------------------------------//

  for(int n_spread = 0; n_spread < nPoisson; n_spread++){
    if( n_spread%100 == 0){  cout << "Iterations\t" << iterations << "\tIteration\t" << n_spread << endl; }

    hFake_new->Reset();
    hMiss_new->Reset();
    hTruth_new->Reset();
    hMeasured_new->Reset();
    hResponse_new->Reset(); 
    //hUnfold->Reset();

    // (3) Fill with Poissonian distribution.
    //-- MC distributions.

    FillAnew_1D( hFake, hFake_new, rand);
    FillAnew_1D( hMiss, hMiss_new, rand);
    FillAnew_2D( hResponse, hResponse_new, rand);
//    FillAnew_1D( hMeasured, hResponse_new->ProjectionX(), hMeasured_new, rand);
//    FillAnew_1D( hTruth, hResponse_new->ProjectionY(), hTruth_new, rand);

    //-- The code below draws the varied distributions and compares them to the actual fakes/measured/truth distributions.
    if( n_spread < 10 ){

      canMiss->cd();
      hMiss_new->SetLineColor( getColor(n_spread) );
      hMiss_new->SetLineStyle( n_spread );
      hMiss_new->DrawClone("histsame");

      canFake->cd();
      hFake_new->SetLineColor( getColor(n_spread) );
      hFake_new->SetLineStyle( n_spread );
      hFake_new->DrawClone("histsame");

      canMeasured->cd();
      hMeasured_new->SetLineColor( getColor(n_spread) );
      hMeasured_new->SetLineStyle( n_spread );
      hMeasured_new->DrawClone("histsame");

      canTruth->cd();
      hTruth_new->SetLineColor( getColor(n_spread) );
      hTruth_new->SetLineStyle( n_spread );
      hTruth_new->DrawClone("histsame");
    }

    //-----------------------//
    // (4) Unfold-and-smear. //
    //-----------------------//

    //-- Determine the number of fakes as: fakes = measured - matched(det)
    TH1D* hRes_X = (TH1D*)hResponse_new->ProjectionX();
    hMeasured_new = (TH1D*)hFake_new->Clone("hMeasured_new");
    hMeasured_new->Add( hRes_X, 1.); 

    TH1D* hRes_Y = (TH1D*)hResponse_new->ProjectionY();
    hTruth_new = (TH1D*)hMiss_new->Clone("hTruth_new");
    hTruth_new->Add( hRes_Y, 1.); 

    testing_covariance << "Variation\tMiss\tFake\tMeasured\tTruth\tResponse\t" <<
	n_spread << "\t" << 
	hMiss_new->Integral() << "\t" <<
	hFake_new->Integral() << "\t" <<
	hMeasured_new->Integral() << "\t" <<
	hTruth_new->Integral() << "\t" <<
	hResponse_new->Integral() << endl;

    //-- Create a new RooUnfold response object.
    RooUnfoldResponse *response_new = new RooUnfoldResponse( hMeasured_new, hTruth_new, hResponse_new );

    //-- Make sure that the unfolded distribution has only real, possible numbers.
    if( !BadNumerals(hUnfold) ){ cout << "This is not good\t" << n_spread << endl; n_spread--; continue; }

    //-----------------//
    //-- S = R * U + F //
    //-----------------//
    hSmear = (TH1D*) response_new->ApplyToTruth( hUnfold );
    hSmear->Add( hFake_new );

    //-- Draw the smeared distribution.
    if( n_spread < 10){
      canData->cd();
      hSmear->SetLineColor( getColor( n_spread ) );
      hSmear->SetLineStyle( n_spread );
      hSmear->Draw("histsame");
    }

    //-- Store histogram in THnSparse, bin-by-bin.
    for(int bin_smear = 1; bin_smear <= hSmear->GetNbinsX(); bin_smear++){
      //-- THnSparse for the original distribution.
      double smear_value = hSmear->GetBinContent( bin_smear );
      double data_value = hData->GetBinContent( bin_smear );
      int sparse_entry_smear[2] = { bin_smear, n_spread };

      //--Look at difference between data and smeared.
      //smear_value -= data_value;

      sparseSmear.SetBinContent( sparse_entry_smear, smear_value );
      if( n_spread < 10){
	//testing_covariance << setprecision(8) <<  smear_value << endl;
      }
    }
  } // Loop over algorithm.

  //-----------------------------------------------------//
  //-- (CHECK) Save canvasses with varied distributions. //
  //-----------------------------------------------------//

  canData->SaveAs( TString::Format( folder_ + "canData_%i.C", iterations) );  
  canFake->SaveAs( TString::Format( folder_ + "canFake_%i.C", iterations) );
  canTruth->SaveAs( TString::Format( folder_ + "canTruth_%i.C", iterations) );
  canMeasured->SaveAs( TString::Format( folder_ + "canMeasured_%i.C", iterations) );

  //-------------------------------------------------//
  // (5) Extract the average from the distributions. //
  //-------------------------------------------------//

  cout << "Unfolder::Calculate_smearedBackError_covariance - (5): Averages	" << endl;
  //-- Project 2nd axis on first axis and divide by nPoisson to get average value.
  TH1D* hAverage_smear = (TH1D*)sparseSmear.Projection( 0 );
  hAverage_smear->Scale( 1./nPoisson );

  TH1D* hAverage_distribution = (TH1D*)hAverage_smear->Clone("hAverage_distribution");

  TCanvas* can_avg;
  PrepareCanvas(can_avg, "Averages");
  hAverage_distribution->Draw("hist");
  can_avg->SaveAs(folder_ + "Average_smeared.C");

  hAverage_smear->Reset();

  //---------------------------------------------------------------------------------------------------------------//
  //-- (6) Loop over all iterations (i.e. all smeared out distributions) and substract the average distribution. --//
  //---------------------------------------------------------------------------------------------------------------//

  // <(E-mu)^2>
  TH1D* hAverage_smear_sq = (TH1D*)sparseSmear_sq.Projection( 0 );
  hAverage_smear_sq->Reset(); 

  TCanvas *can_SmearSpread = new TCanvas("can_smearSpread", "can_smearSpread", 1.) ;
  TString drawoptions = "hist";  
  int color_ = 1;

  for(int bin_smear = 1; bin_smear <= hAverage_smear->GetNbinsX(); bin_smear++){
    //-- Substract average - smeared.
    const double average_smear = hAverage_distribution->GetBinContent( bin_smear );
    double average_e_minus_mu2 = 0.;
    double average_e_minus_mu = 0.;

    TH1D* hSmear_spread = new TH1D(TString::Format("hSmear_spread_%i", bin_smear), TString::Format("hSmear_spread_%i", bin_smear), 100., -3000.,3000.);

    for(int i = 0; i < nPoisson; i++){
      int bin_sparse[2] = { bin_smear, i };
 
      //-- Calculate (E-mu) and (E-mu)² distribution.
      double content_smear = sparseSmear.GetBinContent( bin_sparse );
      content_smear -= average_smear;
      hSmear_spread->Fill( content_smear );

      sparseSmear.SetBinContent( bin_sparse, content_smear );
      sparseSmear_sq.SetBinContent( bin_sparse, content_smear*content_smear );

      average_e_minus_mu2 += content_smear*content_smear;
      average_e_minus_mu += fabs(content_smear);	// This serves only as a check, not as a true quantity to measure.
    }

    //-- Calculate <(E-mu)> and <(E-mu)²>
    average_e_minus_mu2 = average_e_minus_mu2/nPoisson;
    average_e_minus_mu = average_e_minus_mu/nPoisson;

    hAverage_smear_sq->SetBinContent( bin_smear, average_e_minus_mu2 );
    hAverage_smear->SetBinContent( bin_smear, average_e_minus_mu );

    //-- (CHECK) Draw the spread per bin.
    can_SmearSpread->cd();
    if( bin_smear%3 == 0 && bin_smear > 0){ color_++; }
    hAverage_smear_sq->SetLineColor( getColor( color_ ) );
    hAverage_smear_sq->SetLineStyle( color_ );
    hAverage_smear_sq->DrawClone(drawoptions);

    drawoptions = "histsame";
    can_SmearSpread->SaveAs(folder_ + TString::Format("Smear_spread_%i.C", bin_smear) );
  } 

  //-------------------------------------------------//
  //-- Calculate denominator for correlation matrix. //
  //-------------------------------------------------//

  TH2D* hCorrelation_denom = new TH2D("hCorrelation_denom", "hCorrelation_denom", hAverage_smear_sq->GetNbinsX(), 0, 26, hAverage_smear_sq->GetNbinsX(), 0, 26 );
  for(int row = 1; row <= hAverage_smear_sq->GetNbinsX(); row++){
    double nrow = hAverage_smear_sq->GetBinContent( row );

    for(int col = 1; col <= hAverage_smear_sq->GetNbinsX(); col++){
      double ncol = hAverage_smear_sq->GetBinContent( col );

      hCorrelation_denom->SetBinContent( row, col, sqrt( ncol*nrow ) );
      if( sqrt(ncol*nrow) < 1. ){ cout << "Small corr.\t" << col << "\t" << row << "\t" << sqrt( ncol*nrow )<< "\t" << ncol << "\t" << nrow << endl; }

//      //testing_covariance <<"correlation\t(row, col)\t" << row << "\t" << col << "\t" << sqrt( ncol*nrow ) << endl; 
    }
  }

  //-- Poisson-per-poisson iteration 2D histogram.
  TH2D* currentCov = new TH2D("currentCov", "currentCov", hAverage_smear->GetNbinsX(), 0, 26, hAverage_smear->GetNbinsX(), 0, 26);
  TCanvas *can_currentCov;
  PrepareCanvas(can_currentCov, "can_2D");

  hCorrelation_denom->Draw("colz");
  can_currentCov->SaveAs(folder_ +  TString::Format("mCorrelationDenom_%i.C", iterations) );
  can_currentCov->SaveAs(folder_ +  TString::Format("mCorrelationDenom_%i.pdf", iterations) );

  //---------------------------------------------------------------------------------------------------------//
  //-- (7) Loop over all iterations and multiply the different (E-mu_E)_DAT x (E-mu_E)_SM with each other. --//
  //---------------------------------------------------------------------------------------------------------//

  //testing_covariance <<"Unfolder::Calculate_smearedBackError_covariance - (7): Multiply bin-per-bin" << endl;
  //-- We need a THnSparse with nPoisson TH2 matrices in it.
  int bins_3D[3] = {nPoisson, hAverage_smear->GetNbinsX(), hAverage_smear->GetNbinsX() };

  double min_smear = hAverage_smear->GetXaxis()->GetBinLowEdge(1);
  double mins_3D[3] = {0, min_smear, min_smear };

  double max_smear = hAverage_smear->GetXaxis()->GetBinUpEdge( hAverage_smear->GetNbinsX() );
  double maxs_3D[3] = { nPoisson, max_smear, max_smear };

  THnSparseD covariance_3D( "hCovariance_3D", "hCovariance_3D;n_{poisson};E_{smear,i}-#mu_;E_{smear};", 3, bins_3D, mins_3D, maxs_3D );
  THnSparseD correlation_3D( "hCorrelation_3D", "hCorrelation_3D;n_{poisson};E_{smear,i}-#mu_;E_{smear};", 3, bins_3D, mins_3D, maxs_3D );  

  TLegend *leg_smear = new TLegend(0.7, 0.7, 0.95, 0.95);

  hAverage_distribution->SetLineColor( kGreen );
  hAverage_distribution->SetLineWidth( 2 );
  hAverage_distribution->Draw("hist");

  leg_smear->AddEntry(hAverage_distribution, "#mu_{E}", "l");

  can_currentCov->SaveAs(folder_ + "Average_minusSmeared.C");
  can_currentCov->SaveAs(folder_ + "Average_minusSmeared.pdf");

  TCanvas *can_e_mu; 	PrepareCanvas(can_e_mu, "E_mu");  
  TString drawoptions_emu = "hist";
  TH1D* hDistr;

  PrepareCanvas_2D(can_currentCov, "can_2D");

  //-- Loop over variations of response matrix.
  for(int i = 0; i < nPoisson; i++){

   //(CHECK) What does the distribution look like?
   /*
   if( i < 10 ){
     hDistr = (TH1D*)hAverage_distribution->Clone(TString::Format("hE_mu_%i", iterations) );
     hDistr->Reset();
   }
   */

    if( i%500 == 0 ){ cout << "Unfolder::Calculate_smearedBackError_covariance - (7)\t" << i <<  endl; }

    //-- Loop over rows of smear.
    for(int bin_smear_row = 1; bin_smear_row <= hAverage_smear->GetNbinsX(); bin_smear_row++){

      int bin_sparse_smear_row[2] = { bin_smear_row, i };
      double e_muE_smear_row = sparseSmear.GetBinContent( bin_sparse_smear_row );

      //(CHECK) Fill the distribution.
      /*
      if( i < 10 ){
	hDistr->SetBinContent(bin_smear_row, e_muE_smear_row);
 	hDistr->SetLineColor( getColor(i) );
      }
      */

      //cout << "\t" << bin_smear_row << "\t" << i << endl;

      //-- Loop over cols of smear.
      for(int bin_smear_col = 1; bin_smear_col <= hAverage_smear->GetNbinsX(); bin_smear_col++){

	int bin_sparse_smear_col[2] = { bin_smear_col, i };
	double e_muE_smear_col = sparseSmear.GetBinContent( bin_sparse_smear_col );
      
        //-- Calculate value of "covariance".
        double covariance_value = e_muE_smear_row * e_muE_smear_col;
        testing_covariance <<"nPoisson\t" << i << "\t(row, col)\t(" << bin_smear_row << " ,\t" << bin_smear_col << ")\t" << covariance_value << endl;
        //-- Bin.
        int bin_3D[3] = { i, bin_smear_row, bin_smear_col };
        covariance_3D.SetBinContent( bin_3D, covariance_value );

//        currentCov->SetBinContent( bin_smear_row, bin_smear_col, covariance_value);
      }
    }

    //(CHECK) Save the 1D and 2D distributions.
    if( i < 0 ){ 

      cout << "\t" << i << "th average" << endl;

      can_currentCov->cd();
      covariance_3D.GetAxis(0)->SetRange(i,i);
      TH2D* hCov = (TH2D*)covariance_3D.Projection(2, 1);
      hCov->Scale( hCov->Integral()/(i+1.) );

      double maximumPlot = hCov->GetMaximum();
      double minimumPlot = GetMinimumValue( hCov );

      hCov->GetZaxis()->SetRangeUser(minimumPlot * 0.9, maximumPlot * 1.1);
      can_currentCov->SetLogz();
      hCov->Draw("colz");
      can_currentCov->Print( folder_ + "covariance_evolution.gif+10" );
      can_currentCov->SaveAs( folder_ + "covariance_evolution.C" );
    }
  }

  //--------------------------------------------------------------------------------//
  //-- (8) Average out over all histograms by projecting and scaling over nPoisson. //
  //--------------------------------------------------------------------------------//

  //-- Cut off the unneeded parts of vectors and matrices.
  TH1D* hGaugeLengthVectors = (TH1D*)hData->Clone("GaugeLength");
  GetSubHistogram( hData, hGaugeLengthVectors, Ecut_, 2000.);
  int nbins_newlength = hGaugeLengthVectors->GetNbinsX();

  //cout << "***\t" << nbins_newlength << endl;

  cout << "Unfolder::Calculate_smearedBackError_covariance - (8)" << endl;
  TH2D* hCovariance_matrix = (TH2D*) covariance_3D.Projection( 2, 1);
  hCovariance_matrix->Reset();

  //-- Prepare correlation matrix.
  TH2D* hCorrelation_matrix = (TH2D*)hCovariance_matrix->Clone("hCorrelation_matrix");

 
  for(int row = 1; row <= hCovariance_matrix->GetNbinsX(); row++){
    for(int col = 1; col <= hCovariance_matrix->GetNbinsY(); col++){
      double x[nPoisson], y[nPoisson];

      //testing_covariance  << "\t(row, col)\t(" << row << " ,\t" << col << ")\t"; 
      double binval = 0.;
      int poisson = 0;
      for(; poisson < nPoisson; poisson++){
	int bin_sparse_smear_col[3] = { poisson, row, col };        
        double sparse_val = covariance_3D.GetBinContent( bin_sparse_smear_col );
        binval += sparse_val;

	//if( binval == sparse_val && poisson > 0 ){ cout << "\t\t\t---\t" << binval << "\t" << poisson << "\t" << row << "\t" << col << endl; }

	//testing_covariance <<binval/( static_cast<double>(poisson) + 1. ) << "\t";
	x[poisson] = poisson+1;
	y[poisson] = binval/( static_cast<double>(poisson) + 1.);
      }

      //-- Add data to the diagonal of the matrix.
      if( row == col ){ binval += hData->GetBinContent( row ); }

      binval = binval/static_cast<double>(nPoisson);
      //testing_covariance <<binval << endl;
      //-- Set the element in the covariance matrix.

      hCovariance_matrix->SetBinContent( row, col, binval );

      double correlation_denom = hCorrelation_denom->GetBinContent( row, col );
      hCorrelation_matrix->SetBinContent( row, col, binval/correlation_denom);
    }
  }


  TCanvas *can_corr;
  PrepareCanvas_2D( can_corr, TString::Format("canvas_correlation_%i", iterations) );
//  hCorrelation_matrix->Draw("colz");
  TMatrixD *mCorrelation = RooUnfoldResponse::H2M( hCorrelation_matrix, hCorrelation_matrix->GetNbinsX(),hCorrelation_matrix->GetNbinsY()  );
  mCorrelation->ResizeTo( hCovariance_matrix->GetNbinsX() - nbins_newlength, hCovariance_matrix->GetNbinsX()-1, hCovariance_matrix->GetNbinsY() - nbins_newlength, hCovariance_matrix->GetNbinsY()-1 );

  mCorrelation->Draw("colz");
  can_corr->SetLogz();
  can_corr->SaveAs(folder_ + TString::Format("mCor_%i.C", iterations) );
  can_corr->SaveAs(folder_ + TString::Format("mCor_%i.pdf", iterations) );

  TH1D* hDiff = (TH1D*)hData->Clone("hDiff");

  hSmear = (TH1D*) response->ApplyToTruth( hUnfold );
  hSmear->Add( hFake );

  hDiff->Add( hSmear, -1.);	

  for(int bin_data = 1; bin_data <= hDiff->GetNbinsX(); bin_data++){
    //testing_covariance << iterations << "\t" << bin_data << "\t"  << hData->GetBinContent( bin_data ) << "\t" << hSmear->GetBinContent( bin_data ) << "\t"  << hDiff->GetBinContent( bin_data ) << endl;
  }

  TMatrixD *mCovariance = RooUnfoldResponse::H2M( hCovariance_matrix, hCovariance_matrix->GetNbinsX(), hCovariance_matrix->GetNbinsY() );

  //-- Save a copy of the true covariance matrix.
  TCanvas *can_cov = new TCanvas(TString::Format("can_cov_%i", iterations), TString::Format("can_cov_%i", iterations), 1. );
  PrepareCanvas_2D( can_cov, TString::Format("canvas_covariance_%i", iterations) );

  mCovariance->Draw("colz");
  can_cov->SetLogz();
  can_cov->SaveAs(TString::Format("Backsmearing/CovarianceMatrices/mCov_%i.C", iterations) );
  can_cov->SaveAs(TString::Format("Backsmearing/CovarianceMatrices/mCov_%i.pdf", iterations) );

  hCovariance_matrix->Draw("colz");
  can_cov->SetLogz();
  can_cov->SaveAs(folder_ + TString::Format("hCov_%i.C", iterations) );
  can_cov->SaveAs(folder_ + TString::Format("hCov_%i.pdf", iterations) );

  mCovariance->Draw("colz");
  can_cov->SetLogz(),
  can_cov->SaveAs(folder_ +  TString::Format("mCov_%i_Data" + addLabel_ + ".C", iterations) );  
  can_cov->SaveAs(folder_ +  TString::Format("mCov_%i_Data" + addLabel_ + ".pdf", iterations) );  

  //-- Save a copy of the covariance matrix.
  mCovariance->Draw("colz");
  can_cov->SetLogz();
  can_cov->SaveAs( folder_ + TString::Format("mCov_%i_Data" + addLabel_ + ".C", iterations) );  
  can_cov->SaveAs( folder_ + TString::Format("mCov_%i_Data" + addLabel_ + ".pdf", iterations) );  

  TVectorD *vDiff = RooUnfoldResponse::H2V( hDiff, hDiff->GetNbinsX());

  cout << vDiff->GetNoElements() << "\t" << hDiff->GetNbinsX() - nbins_newlength << endl;
    vDiff->ResizeTo( hDiff->GetNbinsX() - nbins_newlength, hDiff->GetNbinsX()-1 );

  //-- Save a copy of the difference between measured and smeared data.
  vDiff->Draw("hist");
  can_cov->SaveAs(folder_ + TString::Format("/vDiff_%i.C", iterations ));
  can_cov->SaveAs(folder_ + TString::Format("/vDiff_%i.pdf", iterations ));

  cout << "\n\n\n" << folder_ + TString::Format("/vDiff_%i.C", iterations ) << endl << endl << endl;

  //-- Invert the matrix.
  TMatrixD mInvertCovariance = (*mCovariance);
  mInvertCovariance.SetTol(1.e-23);
  mInvertCovariance.Invert();

  mInvertCovariance.ResizeTo(hCovariance_matrix->GetNbinsX() - nbins_newlength, hCovariance_matrix->GetNbinsX()-1, hCovariance_matrix->GetNbinsY() - nbins_newlength, hCovariance_matrix->GetNbinsY()-1 );

  TCanvas *can_inv;
  PrepareCanvas_2D( can_inv, "Canvas_inverseMatrices" );
  mInvertCovariance.Draw("colz");
  can_inv->SetLogz();
  can_inv->SaveAs(folder_ + TString::Format("Inverted_C_%i" + addLabel_ + ".C", iterations));
  can_inv->SaveAs(folder_ + TString::Format("Inverted_C_%i" + addLabel_ + ".pdf", iterations));

  TVectorD vTemp = (*vDiff);
  vTemp *= mInvertCovariance;
  double chi2 = (*vDiff) * vTemp;

  can_cov->SetLogz();

  //-- Test if inverse is good inverse.
  TMatrixD mUnity = mInvertCovariance;
  mUnity *= (*mCovariance); //mCovariance_unmodded;
  TCanvas *can_unity = new TCanvas("Unity", "Unity", 1.);
  mUnity.Draw("colz");
  can_unity->SaveAs(folder_ + TString::Format("Unity_%i" + addLabel_ + ".C", iterations));
  can_unity->SaveAs(folder_ + TString::Format("Unity_%i" + addLabel_ + ".pdf", iterations));

  cout << "CHI2 for " << iterations << " iterations is " << chi2 << "�\t" << MC_files_[0] << endl;

  for(int idata = 1; idata <= hData->GetNbinsX(); idata++){
    cout << idata << "\tdata\t" << hData->GetBinContent( idata ) << endl;
  }

  return chi2/ ( vDiff->GetNoElements() ) ;
}


void Unfolder::SetSubhistogram_cut(double Ecut){

  Ecut_ = Ecut;
}


bool Unfolder::BadNumerals( TH1D* hist ){
  bool goodNumerals= true;
  for(int bin = 0; bin <= hist->GetNbinsX(); bin++){
    if( hist->GetBinContent(bin) != hist->GetBinContent(bin) ){ 
      cout << "Shiiit\t" <<  hist->GetBinContent(bin) << endl;
      goodNumerals = false; 
      break; 
    }
  }
  return goodNumerals;
}


void Unfolder::SetAddLabel( TString label ){
  addLabel_ = "_" + label + "_";
}

