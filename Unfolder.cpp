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
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TLatex.h>
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
#include <TString.h>
#include <TObjString.h>
#include <TVectorD.h>

//#include <TPythia8.h>


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
#include "../../../RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "color.h"

#include "CMS_lumi.C"

#include "Functions/Function_average_histogram.h"
#include "Functions/Function_changeBinEdges.h"
#include "Functions/Function_CorrectSectorNumber.h"
#include "Functions/Function_FinishCanvas.C"
#include "Functions/Function_FirstPlot.h"
#include "Functions/Function_PlotSystematics.h"
#include "Functions/Function_Prepare1Dplot.h"
#include "Functions/Function_Prepare2Dplot.h"
#include "Functions/Function_PrepareCanvas.h"
#include "Functions/Function_Rebin.h"
#include "Functions/Function_SetRangeToIncludeAll.h"
#include "GetSubHistogram.h"
#include "Correct_Sector.h"
//#include "Function_make_Tex.h"
#include "Functions/Function_NonZeroMinimum.h"
	
//#define normalise_1 true

// Define the plotter class.
class Unfolder{
  public:
   Unfolder(vector<TString> MC_files, TString datafile, std::map< TString, std::map<TString, TString> >, double Eplotmin, double Ethresh, double deltaPhiMax, double etawidth, TString folder, int normalise); // list of MC files - datafile - map to store histograms.
   ~Unfolder();

   // Get histograms from files.

   void Chi2_comparison_data(TString variable = "lead");
   void Chi2_comparison_MC(TString variable = "lead");
   void Chi2_comparison_MC_det(TString variable = "lead");
   double Chi2_testII( TH1D* hist_data, TH1D* hist_MC);
   double Chi2_testIII( TH1D* hist_data, TH1D* hist_MC);

   void CompareGenLevel();
   void Cut_efficiency(); 	// Counts the number of events in Data: how efficient is each cut?
   void DoublePaddedComparison(TString variable);
   void DoublePaddedComparison_statistics(TString variable, int iterations = 0);
   void DoublePaddedComparison_unfolding(TString variable, int iterations = 0);

   int Get_DetEnergy(int file_, TH1D* &hist_);
   int Get_DetEnergy_unCalib(int file_, TH1D* &hist_);

   int Get_DetPt(int file_, TH1D* &hist_);

   void Get_DetEnergy_lead(int file_, TH1D* &hist_);
   void Get_DetEnergy_JESup(int file_, TH1D* &hist_);
   void Get_DetEnergy_JESdown(int file_, TH1D* &hist_);
   void Get_DetPt_JESup(int file_, TH1D* &hist_);
   void Get_DetPt_JESdown(int file_, TH1D* &hist_);
   void Get_DetUnfolded(int file_, int MC_, TH1D* &hist_, TMatrixD& covariance_m, int iterations = 4, TString variable = "all");
   void Get_DetUnfolded(int file_, TH1D* &hist_, int iterations = 4, TString variable = "all");
   void Get_DetUnfolded_pt(int file_, TH1D* &hist_, int iterations, TString variable);
   void Get_Distribution(int file_, TH1D* &hist_, TString variable);
   void Get_GenEnergy(int file_, TH1D* &hist_);
   void Get_GenPt(int file_, TH1D* &hist_);
   void Get_GenEnergy_lead(int file_, TH1D* &hist_); 
   void Get_GenEnergy_response(int file_, TH1D* &hist_);
   void Get_GenPt_response(int file_, TH1D* &hist_);
   void Get_PureGenEnergy(int file_, TH1D* &hist_);
   void Get_PureGenPt(int file_, TH1D* &hist_);
   void Get_GenSmeared(int file_, TH1D* &hist_, TH1D* hGen, TString variable = "all");
   void Get_GenSmeared(int file_, int MC_, TH1D* &hist_, TH1D* hGen, TString variable = "all");
   void Get_Misses(int file_, TH1D* &hist_);
   void Get_ResponseMatrix(int file_, TH2D* &hist_);
   void Get_ResponseMatrix(TString file_, TH2D* &hist_);
   void Get_ResponseMatrixTHn(int file_, THnSparseD* & hist_);
   void Get_ResponseMatrix_RooUnfold(int file_, TH2D* &hist_);
   void Get_EgenEtaMatrix(int file_, TH2D* &hist_);

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
   void Plot_covariance(int file_, int iterations);

   void Plot_DistributionsResponseObject(int files);
   void Plot_Ratio(TPad* &pad_, TString variable);
   void Plot_Ratio(TPad* &pad_, vector<TH1D*>);
   void Plot_Unfolded();
   void Plot_Unfolded_xf();
   void Plot_Unfolded(TPad* &pad_, TString variable, int iterations = 0);
   void Plot_Unfolded_Ratio();
   void Plot_Unfolded_Ratio(TPad* &pad_, TString variable, int iterations = 0);
   void Plot_Unfolded_Ratio_statistics(TPad* &pad_, TString variable, int iterations = 0);
   void Plot_Unfolded_Ratio_allSystematics(TCanvas* can_, TString variable, int iterations);
   void Plot_Unfolded_Ratio_allSystematics_normData(TCanvas* can_, TString variable, int iterations);
   void Plot_Unfolded_Ratio_allSystematics_pT(TCanvas* can_, TString variable, int iterations, TString plot_as, TString which_models);
   void Plot_Unfolded_Ratio_allSystematics_separate(TCanvas* can_, TString variable, int iterations, TString plot_as, TString which_models);

   //== Systematic uncertainties: plotted on a double padded canvas (abs. and ratio), and the two functions.
   void DoublePaddedComparison_positionDependence(TString variable, int iterations = 0);
   void Plot_Unfolded_positionDependence(TPad* & pad_, TString variable, int iterations);
   void Plot_Unfolded_Ratio_positionDependence(TPad* & pad_, TString variable, int iterations);

   void DoublePaddedComparison_modelDependence(TString variable, int iterations = 0);
   void Plot_Unfolded_modelDependence(TPad* & pad_, TString variable, int iterations);
   void Plot_Unfolded_Ratio_modelDependence(TPad* & pad_, TString variable, int iterations);

   void DoublePaddedComparison_JESDependence(TString variable, int iterations = 0);
   void Plot_Unfolded_JESDependence(TPad* & pad_, TString variable, int iterations);
   void Plot_Unfolded_Ratio_JESDependence(TPad* & pad_, TString variable, int iterations);

   void DoublePaddedComparison_JESmeasured(TString variable, int iterations = 0);
   void Plot_measured_JESDependence(TPad* & pad_, TString variable, int iterations);
   void Plot_measured_Ratio_JESDependence(TPad* & pad_, TString variable, int iterations);
   
   //== Calculate systematic uncertainties.
   void Calculate_jetID_systematic(	int iterations, double *x, 	double *y  );
   void Calculate_JES_systematic( 	int iterations, double *x_up, 	double *y_up, double *x_down, double *y_down );
   void Calculate_model_systematic( 	int iterations, double *x, 	double *y );
   void Calculate_position_systematic( 	int iterations, double *x_up, 	double *y_up, double *x_down, double *y_down );
   void Calculate_lumi_systematic( 	int iterations, double *x, 	double *y);
   void PlotSystematics(int iterations);
   
   //== Plot simple functions.
   void Plot_EtaDiff();
   void Plot_JER();
   void Plot_Isolation();
   void Plot_CastorJetID();   
   void Plot_Unfolded_JES(int iterations_);
   void Plot_Measured_JES();
   void Plot_NjetsMatrix(int file_ );
   void Plot_EgenEtaMatrix(int file_ );

   void Plot_Egen_pt();
  
   //== Plot and compare distributions stored in files for several setups.
   void PlotStartingDistributions();
   void PlotStartingDistributions(TString distribution);
   void PlotStartingDistributions_comparingEmin(TString distribution);
   void PlotStartingDistributions_MCfiles(TString distribution);
   void PlotStartingDistributions_ratio(TString distribution);

   //== Plot stability and purity of unfolded distributions.
   void Plot_Stability( int file_);
   void Plot_Purity( int file_);

   int Scale_to_xsec(int file_, TH1D* &hist_);

   void ScaleToData( int MC, double &scale );
   void Scale_to_Ntotal(int file_, TH1D* &hist_);
   void SetAddLabel( TString label );
   void SetCastorJetEnergy_norm( double renorm );
   void SetScaleFactorData( double scalefactors_Data );
   void SetScaleToData( int yes_or_no);
   void SetSubhistogram_cut(double Ecut);
   void SetSubhistogram_max(double Ecut);
   void SetFitDraw( bool fitnotdrawn );

   void Systematics_CompareDetLevel();
   void Systematics_CompareGenLevel();

   void Get_Njets(int file_, TH2D* &hist_);

   //== Conversion between energy and pT.
   void Convert_E_to_pt(TH1D* &hist_);
   void Calculate_averageEta_perEgen(int file_, TH1D* &hist_);
   void Calculate_averageEgen_perpT(int file_, TGraph* &hist_);
   void Convert_E_to_xF(TH1D* &hist_);
   void Calculate_averageEta_perpt(int file_, TH1D* &hist_);

   void Unfolding_data(TString variable = "lead", int iterations_ = 10);


   //== Calibration functions.
   void CalibrationFactors_oneCanvas(bool draw_functions);
   void CalibrationFunction(int phi_first, int phi_last, TGraphErrors* &gre);
   void CalibrationFunction_workingsectors(int first_sector, int last_sector, TGraphErrors* &gre);
   void CalibrationFunction_sectors(int first_sector, int last_sector, int whichfile, TGraphErrors* &gre);

   void Calibration_SetParLimits(TF1*& analytical);

   void CalculateSystematics(TString setup = "separate", int first_sector = 1, int last_sector = 16);
   void CalculateSystematics_comparison();

   void Plot_Calibrated_functions();
   void Fit_unfolded_distribution( TF1* &analytical);

   //== Response matrix proper.
   void PlotResponseMatrix(int file_, TString setup );
   void PlotResponseMatrix(TString file_, TString setup );
   
   //== Plot average eta vs. egen
   void Plot_averageEta_perEgen(int file_);

   //== Derived from response matrix.
   void PlotFromResponseMatrix(int file_, TString axis_1, TString axis_2);
   void PlotFromResponseMatrix(int file_, TString axis_1);
   void Dissect_ResponseObject();

   //== Chi2 and covariance matrix. Closure test
   void ClosureTest_data(TString variable = "lead", TString file = "Displaced", int method = 1);
   void ClosureTest_MC(TString variable = "lead");
   void ClosureTest_MC_detLevel(TString variable = "lead");
   void Chi2diff_test_data(TString variable, TString file, int method);

   double Calculate_chi2(TH1D* hist_ref, TH1D* hist_res);
   double Calculate_smearedBackError_covariance( TH1D* hData, TH1D* hUnfold, RooUnfoldResponse* response, int iterations );
   void CovarianceMatrix( TH1D* hUnfold, TH1D* hSmeared, TMatrixD& cov_);

   void Smear_gen(int file);

   void Get_UnfoldSmearError( TH1D* vUnfold, TH1D* &vSmeared, int file_, TString variable, int iterations );
   void Unfold_and_smear(int file_, TH1D* &hist_, int MC_, int iterations, TString variable, int method , double &chi2);
   void Unfold_and_smear(int file_, TH1D* &hist_, TFile* MC_, int iterations, TString variable, int method , double &chi2);   

   //== Initialization functions.
   void Initialize_xsec( std::map<TString, double> );
   void Trigger_efficiency( double trigger_eff );
   void PrepareLegend( map<TString, TString> legend_info, map<TString, TString> printLabel, map<TString, TString> legend_info_gen_ );
   void PrepareTitles( map<TString, TString> xtitle ,  map<TString, TString> ytitle ,  map<TString, TString> htitle);


   //Concerning resolution.
   void Determine_resolution( TString setup );
   void Determine_resolution_eDet( TString setup );

   //== Set string to contain cutsetup.
   void SetCutSetup( TString cutsetup );

   void Set_phimax_etaacc( double phimax, double etaacc);


// Stupid checking functions.
void Plot_GenLevels_();

void Calculate_averagepT_perEgen(int file, TH1D* &hist_);

void DetermineScaling();







  private:
  

   
   double lumi_;
   double eBeam_;
   double trigger_efficiency_correction_;

   double norm_cut_;	//== This is used to normalize the backsmearing fakes.


   vector<TString> MC_files_;
   TString datafile_;
   TString folder_; 
   TString cutsetup_;
   TString label_;
   TString setup_calibration_;
   TString addLabel_;

   bool normalise_1;
   double Eplotmin_;
   double Ethresh_; 
   double fitting_threshold_;
   double deltaPhiMax_;
   double etawidth_;
   bool set_xf_;

   double scalefactors_Data_;

   TH2D* hError_hist;

   map<TString, TString> legend_info_;	// Legend entries.
   map<TString, TString> legend_info_gen_;	// Legend entries for plots without any det. level influence.
   map<TString, TString> xtitle_;
   map<TString, TString> ytitle_;
   map<TString, TString> htitle_;
   map<TString, TString> printLabel_;
   map<TString, vector<double> > calibration_parameters_;

   map<TString, TGraphErrors*> calibration_graphs_;

   map<TString, map<TString, TString> > set_of_tags_;

   map<TString, double> xsec_;

   void MakeDoublePaddedComparison(TCanvas * &can, vector<TH1D*>, TLegend *leg);
   void MakeDoublePaddedComparison(TCanvas * &can_, TPad* &pad_abs_, TPad* &pad_ratio_, vector<TH1D*> histos, TLegend *legend);
   
   

   double Chi2_test( TH1D* hist_data, TH1D* hist_MC);
   double Det_to_gen_scale(int file);
   bool BadNumerals( TH1D* hist );

   double renorm_;

   double Ecut_;
   double Emax_;

   bool fitnotdrawn_;
   bool scaletodata_;

   void Histogram_settings_ratio(TH1* hist);
   void Histogram_settings_absolute(TH1* hist);

   void Analyze_response(TH2D* hResponse_selection, TH2D* hGenE_selection, TH2D* hDetE_selection, TH2D* hGenE, TString label, TH1D* hGen, int file, TGraphErrors* &gre_meas, TGraphErrors* &gre_true );
   int Extract_2D_energy_distributions(int file_, TH2D* &hGenE_selection, TH2D* &hDetE_selection, TH2D* &hResponse_selection, TH2D* &hInverseResponse_selection,TH2D* &hDetE, TH2D* &hGenE, TH1D* &hGen, TString setup);

   TGraphErrors* gre_;
   int SetAxisTHnSparse(TString axis);

   void FillAnew_1D(TH1D* hOld, TH1D* &hNew, TRandom3* &rand);
   void FillAnew_2D(TH2D* hOld, TH2D* &hNew, TRandom3* &rand);
   void FillAnew_1D(TH1D* hOld, TH1D* hMin, TH1D* &hNew, TRandom3* &rand);
   
   void Convert_to_xf( TH1D* &hist );
   void Set_to_xf( bool set );
   
   double phimax_global_;
   double etaacc_global_;

};

Unfolder::Unfolder( vector<TString> MC_files, TString datafile, std::map< TString, std::map<TString, TString> > set_of_tags, double Eplotmin, double Ethresh, double deltaPhiMax, double etawidth, TString folder, int normalise )
{
  gStyle->SetStatStyle(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetPalette(1);

  MC_files_ = MC_files;
  datafile_ = datafile;
  folder_ = "Plots/" + folder + "/";
  set_of_tags_ = set_of_tags;
  Eplotmin_ = Eplotmin;
  Ethresh_ = Ethresh;
  deltaPhiMax_ = deltaPhiMax;
  etawidth_ = etawidth;
  addLabel_ = "";  

  lumi_ =  1.23115 * 1e5;
  lumi_ = lumi_*1/0.982;
  eBeam_ = 3500.;


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



void Unfolder::Initialize_xsec( std::map<TString, double> xsec ){  
  xsec_ = xsec;
}

void Unfolder::Trigger_efficiency( double trigger_eff ){
  trigger_efficiency_correction_ = trigger_eff;
}

void Unfolder::SetCutSetup( TString cutsetup ){
  cutsetup_ = cutsetup;
}

void Unfolder::Set_to_xf(bool set){

  set_xf_ = set;

}

void Unfolder::Set_phimax_etaacc( double phimax, double etaacc){

  phimax_global_ = phimax;
  etaacc_global_ = etaacc;

}





/*****************************************
* Plot the Det level energy distribution *
*****************************************/
// -- Plot
int Unfolder::Get_DetEnergy(int file_, TH1D* &hist_){

   cout << "=== GET DET ===" << endl;

   TH1D *hDet;
   TFile *_file0;
   TString drawoptions = "";
   
   if( file_ == -1 ){			_file0 = new TFile( datafile_, "Read");		drawoptions = "data";}
   else if( file_ < MC_files_.size() ){	_file0 = new TFile( MC_files_[file_], "Read"); 	drawoptions = "hist";}
 
   //== Check whether we need energy or xf.
   if( !set_xf_ ){ hDet = (TH1D*)_file0->Get("hCastorJet_energy"); }
   else{ hDet = (TH1D*)_file0->Get("hCastorJet_xf"); }   

   cout << "=== FILE ===\t" << _file0->GetName() << endl;
   cout << "=== xf ===\t"; if( set_xf_ ){ cout << "YES" << endl; } else{ cout << "NO" << endl; }
   cout << "=== DET ===\t" <<  hDet << endl;

   if( !hDet ) return 0;

   hDet->GetXaxis()->SetTitle("E_{det} [GeV]");
   hDet->GetYaxis()->SetTitle("#frac{dN}{dE}");
  
   if( normalise_1 ) hDet->Scale( 1./hDet->Integral() );
 
   hist_ = hDet;
   cout << "Got det" << endl; 
}


int Unfolder::Get_DetPt(int file_, TH1D* &hist_){

   TH1D *hDet;
   TFile *_file0;
   TString drawoptions = "";
   
   if( file_ == -1 ){			_file0 = new TFile( datafile_, "Read");		drawoptions = "data";}
   else if( file_ < MC_files_.size() ){	_file0 = new TFile( MC_files_[file_], "Read"); 	drawoptions = "hist";}

   hDet = (TH1D*)_file0->Get("hCastorJet_pT");

   if( !hDet ) return 0;

//   hDet->Sumw2();
   hDet->GetXaxis()->SetTitle("E_{det} [GeV]");
   hDet->GetYaxis()->SetTitle("#frac{dN}{dE}");
  
   if( normalise_1 ) hDet->Scale( 1./hDet->Integral() );
 
   hist_ = hDet;
   cout << "\t===// -- Detector level energy - end\t" << file_ << endl;   
}


int Unfolder::Get_DetEnergy_unCalib(int file_, TH1D* &hist_){
   cout << "\n\t===// -- Detector level uncalibrated energy - file selection\t" << file_ << endl;

   TH1D *hDet;
   TFile *_file0;
   TString drawoptions = "";
   
   cout << "\t===//\tUnfolder::Get_DetEnergy_unCalib\topening files" << endl;
   if( file_ == -1 ){			_file0 = new TFile( datafile_, "Read");		drawoptions = "data";}
   else if( file_ < MC_files_.size() ){	_file0 = new TFile( MC_files_[file_], "Read"); 	drawoptions = "hist";}

   cout << "\tUnfolder::Get_DetEnergy\topened file" << endl;
         
   hDet = (TH1D*)_file0->Get("hCastorJet_energy_unCalib");

   if( !hDet ) return 0;

   hDet->Sumw2();
   hDet->GetXaxis()->SetTitle("E_{det} [GeV]");
   hDet->GetYaxis()->SetTitle("#frac{dN}{dE}");
 
 
   hist_ = hDet;
   cout << "\t===// -- Detector level uncalibrated energy - end\t" << file_ << endl;   
}



void Unfolder::Get_DetEnergy_lead(int file_, TH1D* &hist_){
   TH1D *hDet;
   TFile *_file0;
   TString drawoptions = "";

   if( file_ == -1 ){              _file0 = new TFile( datafile_, "Read");         drawoptions = "data";}
   else if( file_ < MC_files_.size() ){      _file0 = new TFile( MC_files_[file_], "Read");  drawoptions = "hist";}


   hDet = (TH1D*)_file0->Get("hCastorJet_energy_lead");      hDet->Sumw2();
   hDet->GetXaxis()->SetTitle("E_{det} [GeV]");
   hDet->GetYaxis()->SetTitle("#frac{dN}{dE}");

   if( normalise_1 ) hDet->Scale( 1./hDet->Integral() );

   hist_ = hDet;
}

void Unfolder::SetScaleFactorData( double scalefactors_Data ){
 scalefactors_Data_ = scalefactors_Data; 
}



void Unfolder::Get_DetEnergy_JESup(int file_, TH1D* &hist_){
   cout << "\n\n\n// -- Detector level energy - JESup\t" << file_ << endl;

   TH1D *hDet;
   TFile *_file0;
   TString drawoptions = "";
   
   if( file_ == -1 ){			_file0 = new TFile( datafile_, "Read");		drawoptions = "data";}
   else if( file_ < MC_files_.size() ){	_file0 = new TFile( MC_files_[file_], "Read"); 	drawoptions = "hist";}

         
   if( !set_xf_ ){ hDet = (TH1D*)_file0->Get("hCastorJet_energy_JES_up");	hDet->Sumw2(); }
   if( set_xf_ ){ hDet = (TH1D*)_file0->Get("hCastorJet_xf_JES_up");	hDet->Sumw2(); }
      
   hDet->GetXaxis()->SetTitle("E_{det} [GeV]");
   hDet->GetYaxis()->SetTitle("#frac{dN}{dE}");
  
   if( normalise_1 ) hDet->Scale( 1./hDet->Integral() );
 
   hist_ = hDet;
}




void Unfolder::Get_DetEnergy_JESdown(int file_, TH1D* &hist_){
   cout << "\n\n\n// -- Detector level energy - JESup\t" << file_ << endl;

   TH1D *hDet;
   TFile *_file0;
   TString drawoptions = "";
   
   if( file_ == -1 ){			_file0 = new TFile( datafile_, "Read");		drawoptions = "data";}
   else if( file_ < MC_files_.size() ){	_file0 = new TFile( MC_files_[file_], "Read"); 	drawoptions = "hist";}

         
   if( !set_xf_ ){ hDet = (TH1D*)_file0->Get("hCastorJet_energy_JES_down");	hDet->Sumw2(); }
   if( set_xf_ ){ hDet = (TH1D*)_file0->Get("hCastorJet_xf_JES_down");	hDet->Sumw2(); }
   
   hDet->GetXaxis()->SetTitle("E_{det} [GeV]");
   hDet->GetYaxis()->SetTitle("#frac{dN}{dE}");
  
   if( normalise_1 ) hDet->Scale( 1./hDet->Integral() );
 
   hist_ = hDet;
}





void Unfolder::Get_DetPt_JESup(int file_, TH1D* &hist_){

   TH1D *hDet;
   TFile *_file0;
   TString drawoptions = "";
   
   if( file_ == -1 ){			_file0 = new TFile( datafile_, "Read");		drawoptions = "data";}
   else if( file_ < MC_files_.size() ){	_file0 = new TFile( MC_files_[file_], "Read"); 	drawoptions = "hist";}
        
   hDet = (TH1D*)_file0->Get("hCastorJet_pT_JES_up");	hDet->Sumw2();
   hDet->GetXaxis()->SetTitle("p_{T,det} [GeV]");
   hDet->GetYaxis()->SetTitle("#frac{dN}{dp_{T}}");
  
   if( normalise_1 ) hDet->Scale( 1./hDet->Integral() );
 
   hist_ = hDet;

}

void Unfolder::Get_DetPt_JESdown(int file_, TH1D* &hist_){

   TH1D *hDet;
   TFile *_file0;
   TString drawoptions = "";
   
   if( file_ == -1 ){			_file0 = new TFile( datafile_, "Read");		drawoptions = "data";}
   else if( file_ < MC_files_.size() ){	_file0 = new TFile( MC_files_[file_], "Read"); 	drawoptions = "hist";}

         
   hDet = (TH1D*)_file0->Get("hCastorJet_pT_JES_down");	hDet->Sumw2();
   hDet->GetXaxis()->SetTitle("p_{T,det} [GeV]");
   hDet->GetYaxis()->SetTitle("#frac{dN}{dp_{T}}");
  
   if( normalise_1 ) hDet->Scale( 1./hDet->Integral() );
 
   hist_ = hDet;
}






void Unfolder::Get_ResponseMatrix(int file_, TH2D* &hist_){

   TFile *_file0;
   TString drawoptions = "";

   if( file_ < MC_files_.size() ){      _file0 = new TFile( MC_files_[file_], "Read");  drawoptions = "hist";}
   else if( file_ == -1 ){              _file0 = new TFile( datafile_, "Read");         drawoptions = "data";}

   TH2D* hRes = (TH2D*)_file0->Get("hCastorJet_energy_response_fine");      hRes->Sumw2();
   hRes->GetXaxis()->SetTitle("E_{det} [GeV]");
   hRes->GetYaxis()->SetTitle("E_{gen} [GeV]");

   if( normalise_1 ) hRes->Scale( 1./hRes->Integral() );

   hist_ = hRes;
}






void Unfolder::Get_ResponseMatrix(TString file_, TH2D* &hist_){

   TFile *_file0;
   TString drawoptions = "";

   _file0 = new TFile( file_, "Read");  

   TH2D* hRes = (TH2D*)_file0->Get("hCastorJet_energy_response_fine");      hRes->Sumw2();
   hRes->GetXaxis()->SetTitle("E_{det} [GeV]");
   hRes->GetYaxis()->SetTitle("E_{gen} [GeV]");

   if( normalise_1 ) hRes->Scale( 1./hRes->Integral() );

   hist_ = hRes;
}





void Unfolder::Get_ResponseMatrix_RooUnfold(int file_, TH2D* &hist_){

   TFile *_file0;
   TString drawoptions = "";

   _file0 = new TFile( MC_files_[file_], "Read");  

   TH2D* hRes = (TH2D*)_file0->Get("hCastorJet_energy_response_fine");      hRes->Sumw2();
   hRes->GetXaxis()->SetTitle("E_{det} [GeV]");
   hRes->GetYaxis()->SetTitle("E_{gen} [GeV]");

   RooUnfoldResponse *response = (RooUnfoldResponse*)_file0->Get("response");
   hist_ = (TH2D*)response->Hresponse();
}








void Unfolder::Get_Njets(int file_, TH2D* &hist_){

   TFile *_file0;
   TString drawoptions = "";

   if( file_ < MC_files_.size() ){      _file0 = new TFile( MC_files_[file_], "Read");  drawoptions = "hist";}
   else if( file_ == -1 ){              _file0 = new TFile( datafile_, "Read");         drawoptions = "data";}

   TH2D* hRes = (TH2D*)_file0->Get("hNumber_of_match_jets");      hRes->Sumw2();
   hRes->GetXaxis()->SetTitle("N_{det}");
   hRes->GetYaxis()->SetTitle("N_{gen}");

   if( normalise_1 ) hRes->Scale( 1./hRes->Integral() );

   hist_ = hRes;
}


void Unfolder::Get_EgenEtaMatrix(int file_, TH2D* &hist_){

   TFile *_file0;
   TString drawoptions = "";

   if( file_ < MC_files_.size() ){      _file0 = new TFile( MC_files_[file_], "Read");  drawoptions = "hist";}
   else if( file_ == -1 ){              _file0 = new TFile( datafile_, "Read");         drawoptions = "data";}

   TH2D* hRes = (TH2D*)_file0->Get("hGenJet_energy_vs_eta");      hRes->Sumw2();
   hRes->GetXaxis()->SetTitle("E_{gen} [GeV]");
   hRes->GetYaxis()->SetTitle("#eta");

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

    //-- Data.
    if( file == MC_files_.size() ){
      file_ = -1;
      drawoptions = "edatasame";
      legendoptions = "p";
      done_data = true;
    }
 
    //-- MC.
    else{
      TString scalefactor_MC = (set_of_tags_[ "scalefactors" ])[ MC_files_[file_] ];
      double scalefactors_MC_ = scalefactor_MC.Atof() ;
      SetCastorJetEnergy_norm( scalefactors_MC_ / scalefactors_Data_ ); 
    }

    //-- Extract detector level distribution and set the properties of the histogram.
    Get_DetEnergy_lead( file_ , hist_ );
    hist_->SetLineColor( getColor( file_ +1) );
    hist_->SetLineWidth( 3 );
    hist_->SetMarkerColor( getColor( file_ +1) );
    hist_->SetMarkerSize( 1 );

    if( file_ != -1) {
      hist_->Scale(1./renorm_);
      histos.push_back( hist_ ); 
    }
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
  cout << "\n\nUnfolder::Hist_DetLevel()\n\n";

  //-- Prepare legend etc.
  TH1D* hist_;
  TString drawoptions = "hist", legendoptions = "l";
  TCanvas *can; 
  PrepareCanvas( can, "Detector_Level_Energy" );
  vector<TH1D*> histos;

  TLegend*leg = new TLegend( 
	1. - can->GetRightMargin()-0.3, 1. - can->GetTopMargin() - 0.4,
	1.-can->GetRightMargin(), 	1.-can->GetTopMargin()*4./3.);
  leg->SetFillColor(0);

  bool done_data = false;

  //-- Loop over files.
  for( int file = 0; file <= MC_files_.size(); file++ ){
    int file_ = file;
   
    //-- Data.
    if( file == MC_files_.size() ){
      file_ = -1;
      drawoptions = "epsame";
      legendoptions = "p";
      done_data = true;
    }

    else{
      if( (set_of_tags_["mc_type"])[MC_files_[file_]] != "model" && 
	  (set_of_tags_["mc_type"])[MC_files_[file_]] != "actual" &&  
	  (set_of_tags_["mc_type"])[MC_files_[file_]] != "model_"   ) { continue; }
    }

//    Get_DetEnergy_unCalib( file_ , hist_ );     hist_->GetXaxis()->SetRangeUser( 150., 2100./1.5 );
    Get_DetEnergy( file_ , hist_ );


    if( !hist_) continue;
    if( hist_->Integral() != hist_->Integral() ) continue;

    //== Scale MC to number of events.
    if( file_ != -1 ){  Scale_to_Ntotal(file, hist_); }

    if( hist_->Integral() == 0. ) continue;	

    hist_->SetLineColor( getColor( file_ +1) );
    hist_->SetLineWidth( 3 );
    hist_->SetMarkerColor( getColor( file_ +1) );
    hist_->SetMarkerSize( 1.5 );
    hist_->GetYaxis()->SetTitle( "N" );

//    SetDnDx( hist_ );
    Prepare_1Dplot( hist_ );

    if( file_ != -1) {histos.push_back( hist_ ); }
    if( file_ == -1) {
      hist_->SetLineColor( 0 );
      histos.insert( histos.begin(), hist_) ; }

    TString legend_info;
    if( file_ == -1 ){ legend_info = "Data"; }
    else{ legend_info = legend_info_gen_[MC_files_[file_]];  }

    leg->AddEntry( hist_, legend_info, legendoptions );
    drawoptions = "histsame";
    if( done_data ){ break; }
  }

  TPad* pad_abs_, *pad_ratio_;
  MakeDoublePaddedComparison(can, pad_abs_, pad_ratio_, histos, leg );

  can->cd();

  Finish_canvas( pad_abs_ );
 
  can->SaveAs( TString::Format( folder_ + "DetectorLevelEnergy" + label_ + ".C") );
  can->SaveAs( TString::Format( folder_ + "DetectorLevelEnergy" + label_ + ".pdf") );
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
int Unfolder::Scale_to_xsec(int file_, TH1D* &hist_){
   cout << "\t===// -- Scale to x-section - file selection (file, hist)\t" << file_ << endl;

   TH1D *hGen;
   TFile *_file0;
   TString drawoptions = "";
   
   if( file_ < MC_files_.size() ){	_file0 = new TFile( MC_files_[file_], "Read"); 	drawoptions = "hist";}

   int total_events_nocuts;
   TTree *tree_numbers = (TTree*)_file0->Get("useful_numbers");
   TString branch_str = (set_of_tags_[ "scalefactors" ])[ MC_files_[file_] ];
   tree_numbers->SetBranchAddress(branch_str, &total_events_nocuts);

   tree_numbers->GetEntry( 0 );

   double xsec = 1.;
   xsec = xsec_[MC_files_[file_]];
   cout << "\t==xsec\t" << xsec << endl;

   if( xsec == 0. ) return 0;

   hist_->Scale( xsec/ total_events_nocuts);

   cout << MC_files_[file_] << "\tNormalize to\t" << total_events_nocuts << "\t" << xsec << endl;

   cout << "\nReached the end of Unfolder::Scale_to_xsec" << endl;
   cout << "==Returning histogram\t" << hist_->Integral()  << endl << endl << endl;
}





// -- Plot
void Unfolder::Scale_to_Ntotal(int file_, TH1D* &hist_){
   

   cout << "\t===// -- Scale to Nevents - file selection (file, hist)\t" << file_ << endl;

   TH1D *hGen;
   TFile *_file0;
   TString drawoptions = "";
   TString branch_str = "castor_150GeV_events";   
   TString branch_data = "castor_150GeV_events";   

   if( (set_of_tags_["mc_type"])[MC_files_[file_]] == "shift_MPI_or_Tune" ){
     branch_str = "castor_150GeV_events_MC";
   }

   //Extract number of events from MC file.
   if( file_ < MC_files_.size() ){	_file0 = new TFile( MC_files_[file_], "Read"); }
   int total_events_MC;
   TTree *tree_numbers = (TTree*)_file0->Get("useful_numbers");
   tree_numbers->SetBranchAddress(branch_str, &total_events_MC);
   tree_numbers->GetEntry( 0 );

   //Extract number of events from data file.
   TFile* _fileData = new TFile( datafile_, "Read");
   int total_events_data;
   TTree *tree_numbers_data = (TTree*)_fileData->Get("useful_numbers");
   tree_numbers_data->SetBranchAddress(branch_data, &total_events_data);
   tree_numbers_data->GetEntry( 0 );

   hist_->Scale( double(total_events_data)/double(total_events_MC));

   cout << MC_files_[file_] << "\tNormalize to\t" << total_events_data/total_events_MC <<  endl;

   cout << "\nReached the end of Unfolder::Scale_to_Ntotal" << endl;
   cout << "==Returning histogram\t" << hist_->Integral()  << endl << endl << endl;
}






// -- Plot
void Unfolder::Get_GenEnergy(int file_, TH1D* &hist_){
   cout << "\t\t// -- Generator level energy - file selection (file, hist)\t" << file_ << endl;

   TH1D *hGen;
   TFile *_file0;
   TString drawoptions = "";
   
   if( file_ < MC_files_.size() ){	_file0 = new TFile( MC_files_[file_], "Read"); 	drawoptions = "hist";}

   if( !set_xf_ ){
     hGen = (TH1D*)_file0->Get("hGenJet_energy");
     if( !hGen ) hGen = (TH1D*)_file0->Get("hGenJet_Pythia84C_MPI_off");	
     if( !hGen ) hGen = (TH1D*)_file0->Get("hGenJet_Pythia84C_MPI_on");	
   }
   else{ 
     hGen = (TH1D*)_file0->Get("hGenJet_xf");
     cout << "Get gen xf" << endl;
   }

   if( hGen ){
   cout << "\nExtracted GenJetEnergy" << endl;

   hGen->Sumw2();
   hGen->GetXaxis()->SetTitle("E_{det} [GeV]");
   hGen->GetYaxis()->SetTitle("#frac{dN}{dE}");

   if( scaletodata_  ){
      int total_events_nocuts;
      TTree *tree_numbers = (TTree*)_file0->Get("useful_numbers");
      TString branch_str = (set_of_tags_[ "scalefactors" ])[ MC_files_[file_] ];
      tree_numbers->SetBranchAddress(branch_str, &total_events_nocuts);     

      tree_numbers->GetEntry( 0 );

      double xsec = 1.;
      xsec = xsec_[MC_files_[file_]];
      //hGen->Scale( xsec/ total_events_nocuts);

      cout << "\t***\tSuccesful event\t" << (set_of_tags_[ "scalefactors" ])[ MC_files_[file_] ] << "\t" << total_events_nocuts << endl;
      cout << "\t***\t" << MC_files_[file_] << "\t" << hGen->Integral() << "\t" << endl;

      cout << MC_files_[file_] << "\tNormalize to\t" << total_events_nocuts << endl;
   }


   if( normalise_1 ) hGen->Scale( 1./hGen->Integral() );
   
   hist_ = hGen;

   cout << "\nReached the end of Unfolder::Get_GenEnergy" << endl;
   cout << "==Returning histogram\t" << hGen->Integral()  << endl << endl << endl;
   
   }
}






void Unfolder::Get_GenPt(int file_, TH1D* &hist_){
   TH1D *hGen;
   TFile *_file0;
   TString drawoptions = "";
   
   if( file_ < MC_files_.size() ){	_file0 = new TFile( MC_files_[file_], "Read"); 	drawoptions = "hist";}

   hGen = (TH1D*)_file0->Get("hGenJet_pT");
   if( !hGen ) hGen = (TH1D*)_file0->Get("hGenJet_Pythia84C_MPI_off");	
   if( !hGen ) hGen = (TH1D*)_file0->Get("hGenJet_Pythia84C_MPI_on");	

   hGen->Sumw2();
   hGen->GetXaxis()->SetTitle("p_{T,det} [GeV]");
   hGen->GetYaxis()->SetTitle("#frac{dN}{dp_{T}}");

   if( scaletodata_ && !MC_files_[file_].Contains("HEJ") ){
      int total_events_nocuts;
      TTree *tree_numbers = (TTree*)_file0->Get("useful_numbers");
      TString branch_str = (set_of_tags_[ "scalefactors" ])[ MC_files_[file_] ];
      tree_numbers->SetBranchAddress(branch_str, &total_events_nocuts);     

      tree_numbers->GetEntry( 0 );

      double xsec = 1.;
      xsec = xsec_[MC_files_[file_]];
      hGen->Scale( xsec/ total_events_nocuts);
    }

   if( normalise_1 ) hGen->Scale( 1./hGen->Integral() );
   
   hist_ = hGen;
}





void Unfolder::Get_GenEnergy_response(int file_, TH1D* &hist_){
   cout << "\n\n\n// -- Generator level energy response - file selection (file, hist)\t" << file_ << endl;
   TH1D *hGen;
   TFile *_file0;
   TString drawoptions = "";
   
   if( file_ < MC_files_.size() ){	_file0 = new TFile( MC_files_[file_], "Read"); 	drawoptions = "hist";}
   

   RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response");

   hGen = (TH1D*)response->Htruth();
   hGen->GetXaxis()->SetTitle("E_{det} [GeV]");
   hGen->GetYaxis()->SetTitle("#frac{dN}{dE}");

   if( normalise_1 ) hGen->Scale( 1./hGen->Integral() );
   
   hist_ = hGen;
}






void Unfolder::Get_GenPt_response(int file_, TH1D* &hist_){
   cout << "\n\n\n// -- Generator level energy response - file selection (file, hist)\t" << file_ << endl;
   TH1D *hGen;
   TFile *_file0;
   TString drawoptions = "";
   
   if( file_ < MC_files_.size() ){	_file0 = new TFile( MC_files_[file_], "Read"); 	drawoptions = "hist";}
   

   RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response_pT");

   hGen = (TH1D*)response->Htruth();
   hGen->GetXaxis()->SetTitle("E_{det} [GeV]");
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
   hGen->GetXaxis()->SetTitle("E_{det} [GeV]");
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
   hGen->GetXaxis()->SetTitle("E_{det} [GeV]");
   hGen->GetYaxis()->SetTitle("#frac{dN}{dE}");



   if( normalise_1 ) hGen->Scale( 1./hGen->Integral() );
   
   hist_ = hGen;
}





void Unfolder::Get_PureGenEnergy(int file_, TH1D* &hist_){
   cout << "\t\t// -- Generator level energy - file selection\t" << file_ << endl;
   TH1D *hGen;
   TFile *_file0;
   TString drawoptions = "";
   
   if( file_ < MC_files_.size() ){	_file0 = new TFile( MC_files_[file_], "Read"); 	drawoptions = "hist";}
  
   if( !set_xf_ ){ hGen = (TH1D*)_file0->Get("hGenJet_energy");	 }
   if( set_xf_ ){ hGen = (TH1D*)_file0->Get("hGenJet_xf");	 }
   
   //-- Continue if the histogram exists.
   if( hGen ){
     cout << "Unfolder--Get_PureGenEnergy\thGen exists" << endl;
          
     //-- Continue if the histogram must be scaled to its number of events.
     //-- 
     if( scaletodata_ ){

      int total_events_nocuts;
      TTree *tree_numbers = (TTree*)_file0->Get("useful_numbers");
      TString branch_str = (set_of_tags_[ "scalefactors" ])[ MC_files_[file_] ];
      tree_numbers->SetBranchAddress(branch_str, &total_events_nocuts);
      tree_numbers->GetEntry( 0 );

      double xsec = 1.;
      xsec = xsec_[MC_files_[file_]];
      hGen->Scale( xsec/ total_events_nocuts);
      cout << MC_files_[file_] << "\tNormalize to\t" << total_events_nocuts << endl;
    }
    hGen->Sumw2();
    hGen->GetXaxis()->SetTitle("E_{det} [GeV]");
    hGen->GetYaxis()->SetTitle("#frac{dN}{dE}");
    
  }

  hist_ = hGen;

   cout << "\nReached the end of Unfolder--Get_PureGenEnergy" << endl;
   cout << "==Returning histogram\t" << hGen->Integral() << endl << endl << endl;
}



void Unfolder::Get_PureGenPt(int file_, TH1D* &hist_){
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

  can->SaveAs(folder_ + "GeneratorLevelEnergy" + label_ + "_lead.C");
  can->SaveAs(folder_ + "GeneratorLevelEnergy" + label_ + "_lead.pdf");
  
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
  cout << "\n\n\tUnfolder::Plot_Unfolded" << endl;

  //-- Loop over all MC files.
  for(int file_ = 0; file_ < 1; file_++){ //MC_files_.size(); file_++){

    int color_ = 1;
    TH1D* hGen, *hDet, *hFirst, *hData_ratio;
    TString legend_file = legend_info_gen_[MC_files_[file_]];
    TString print_label = printLabel_[MC_files_[file_]];
  
    TCanvas *can;
    TPad *pad_abs_, *pad_ratio_;
    PrepareCanvas( can , "Comparison_Energy_Unfolded");
    TLegend* leg = new TLegend(0.6, 0.48, 1. - can->GetRightMargin(), 1. - can->GetTopMargin());  
    leg->SetFillColor(0);

    double min_val, max_val;
    TString drawoptions = "hist";
    bool first_iteration = true;
    int it_incr = 3;
    int it_max = 30;

    for( int iterations = 1; iterations <= it_max; iterations += it_incr, color_++){
          

      Get_DetEnergy(-1, hDet);
      Get_DetUnfolded(file_, hDet, iterations);
      hDet->Scale( 1./(1.2*1e5) );      
      hDet->SetName( TString::Format( "Iterations_%i", iterations) );

      SetDnDx( hDet );
      hDet->SetLineColor( getColor(color_) );
      hDet->SetLineStyle( color_ );
      hDet->SetLineWidth( 2 );
      hDet->GetXaxis()->SetNdivisions(504);
      hDet->SetXTitle("E [GeV]");
      hDet->SetYTitle("#frac{d#sigma}{dE} [mb/GeV]");

      can->cd();
      Prepare_1Dplot( hDet );
      hDet->DrawClone( drawoptions );
      drawoptions = "histsame";


      if( first_iteration){
	hFirst = (TH1D*)hDet->Clone("hFirst");
	first_iteration = false;
      }

      First_Plot( hFirst, hDet, iterations, min_val, max_val);

//      leg->AddEntry( hDet, TString::Format("%i it., #Delta#varphi = 0.%i, #eta = 0.%i", iterations, static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) ), "l");
      leg->AddEntry( hDet, TString::Format("%i it., " + legend_info_gen_[MC_files_[0]], iterations), "l");


    }
    leg->Draw();
    can->SetLogy();
    can->Update();
  
    TString variable = "energy";
    if( set_xf_ ) variable = "xf"; 
  
    can->SaveAs(folder_ + "/Compare_Unfolded_" + variable + "_" + print_label + ".C");
    can->SaveAs(folder_ + "/Compare_Unfolded_" + variable + "_" + print_label + ".pdf");  

    can->SaveAs( TString::Format(folder_ + "Compare_UnfoldedData_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i.C", 
	it_max, static_cast<int>(Ethresh_), static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) ) );
    can->SaveAs( TString::Format(folder_ + "Compare_UnfoldedData_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i.pdf", 
	it_max, static_cast<int>(Ethresh_), static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) ) );
  }  
}

void Unfolder::Plot_Unfolded_xf(){

  set_xf_ = true;
  Plot_Unfolded();
  set_xf_ = false;


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
      Prepare_1Dplot( hFirst );
    }

    Prepare_1Dplot( hDistr );

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
    hDistr->GetYaxis()->SetNdivisions(205);
    hDistr->Scale( 1./hDistr->Integral() );

    if( firstDistr ){
      hFirst = (TH1D*)hDistr->Clone(variable + "_first");
      hFirst_abs = (TH1D*)hFirst->Clone(variable + "_first_absolute");
      hFirst->Divide( hFirst );
      hFirst->GetYaxis()->SetTitle("Ratio");
      firstDistr = false;
      hFirst->GetYaxis()->SetLabelFont(42);
      hFirst->GetYaxis()->SetLabelSize( (hFirst->GetYaxis()->GetLabelSize())*1.5 );
      hFirst->GetYaxis()->SetLabelSize( (hFirst->GetYaxis()->GetLabelSize())*1.5 );
      hFirst->GetYaxis()->SetLabelOffset( hFirst->GetYaxis()->GetLabelOffset() );


      hFirst->GetXaxis()->SetLabelFont(42);
      hFirst->GetXaxis()->SetLabelSize( ( hFirst->GetXaxis()->GetLabelSize())*1.5 );
      hFirst->GetXaxis()->SetTitleSize( ( hFirst->GetXaxis()->GetTitleSize())*1.5 );
      hFirst->GetXaxis()->SetLabelOffset( hFirst->GetXaxis()->GetLabelOffset()*1.5 );
    }
   
    hDistr->Divide( hFirst_abs );
    First_Plot( hFirst, hDistr, nDistr, max_val, min_val);

    hDistr->GetXaxis()->SetRangeUser(0., 2000.);
    hDistr->GetYaxis()->SetRangeUser(0., 3.);

    hDistr->SetLineColor( getColor( file_+1 ) );
    hDistr->SetLineWidth( 3 );
    hDistr->SetLineStyle( file_+1 );
    hDistr->SetMarkerColor(  getColor( file_+1 ) );

    Prepare_1Dplot_ratio( hDistr );



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
  TFile *_file0 = new TFile( MC_files_[MC_], "read");
  RooUnfoldResponse* response;
  
  if( !set_xf_ ){
    if( variable == "all"){ response = (RooUnfoldResponse*)_file0->Get("response"); }
    if( variable == "lead"){ response = (RooUnfoldResponse*)_file0->Get("response_lead"); }
  }
  else{
    response = (RooUnfoldResponse*)_file0->Get("response_xf"); 
    if( !response ){ cout << "No such response" << endl; }
  }

  cout << "=== RESPONSE ===\t" << response->GetName() << endl;

  // Use hist as input for the unfolding.
  RooUnfoldBayes unfold_bayes(response, hist_, iterations);
  // the original hist has become obsolete, transform it into the unfolded histogram.
  hist_ = (TH1D*) unfold_bayes.Hreco( RooUnfold::kCovariance ); 

  RooUnfold::ErrorTreatment et = RooUnfold::kCovariance;

//  covariance_m = (TMatrixD) unfold_bayes.Ereco( et );

}




void Unfolder::Get_DetUnfolded(int file_, TH1D* &hist_, int iterations, TString variable){
  cout << "\t\tUnfolder--Get_DetUnfolded(" << file_ << ", "  << iterations << ", " << variable << ") no Matrix" << endl;

  if( file_ >= 0 ){
  TFile *_file0 = new TFile( MC_files_[file_], "read");
  cout << "\tFile is\t" << MC_files_[file_] << endl;
  cout << "\tHist is\t" << hist_->GetTitle() << "\t" << hist_->Integral() << endl;

  RooUnfoldResponse* response;
  if( !set_xf_ ){
    if( variable == "all"){ response = (RooUnfoldResponse*)_file0->Get("response"); }
    if( variable == "lead"){ response = (RooUnfoldResponse*)_file0->Get("response_lead"); }
  }
  else{
    response = (RooUnfoldResponse*)_file0->Get("response_xf"); 
  }

  cout << "=== RESPONSE ===\t" << response->GetName() << endl;

  // Use hist as input for the unfolding.
  RooUnfoldBayes unfold_bayes(response, hist_, iterations);
  // the original hist has become obsolete, transform it into the unfolded histogram.
  hist_ = (TH1D*) unfold_bayes.Hreco( RooUnfold::kCovariance ); 

  //== The following lines extract the unfolding matrix and plot it as a histogram.
  TCanvas * can_unfold = new TCanvas("unfolded", "unfodled", 1.);
  PrepareCanvas_2D( can_unfold, TString::Format("Unfolding_matrix_%i_iterations", iterations) );
  TMatrixD& unfold_m = (TMatrixD&)unfold_bayes.UnfoldingMatrix();
  TH2D* unfold_h = new TH2D( unfold_m );   
  double min_unfold = GetMinimumValue( unfold_h );
  double max_unfold = unfold_h->GetMaximum();
  unfold_h->GetZaxis()->SetRangeUser( min_unfold * 0.9, max_unfold * 1.1);
  Prepare_2Dplot( unfold_h );

  unfold_h->GetXaxis()->SetTitle("E_{det} [GeV]");
  unfold_h->GetYaxis()->SetTitle("E_{gen} [GeV]");

//  unfold_h->Draw("colz");
  
  can_unfold->SetLogz();
  //== Done drawing.


  //== This draw the covariance matrix.
  TCanvas * can_cov = new TCanvas("covariance", "covariance", 1.);
  PrepareCanvas_2D( can_cov, TString::Format("Covariance_unfolding_matrix_%i_iterations", iterations) );
  TMatrixD& covariance_m = (TMatrixD&)unfold_bayes.GetMeasuredCov();
  TH2D* covariance_h = new TH2D( covariance_m );   
  double min_covariance = GetMinimumValue( covariance_h );
  double max_covariance = covariance_h->GetMaximum();
  covariance_h->GetZaxis()->SetRangeUser( min_covariance * 0.9, max_covariance * 1.1);
  Prepare_2Dplot( covariance_h );

  covariance_h->GetXaxis()->SetTitle("E_{det} [GeV]");
  covariance_h->GetYaxis()->SetTitle("E_{gen} [GeV]");
  
  can_cov->SetLogz();
  

  //== Done with the covariance matrix.

  cout << "\t\tUnfolder::Get_DetUnfolded after\t" << hist_->Integral() << endl;

  }
}






void Unfolder::Get_DetUnfolded_pt(int file_, TH1D* &hist_, int iterations, TString variable){

  if( file_ >= 0 ){
  TFile *_file0 = new TFile( MC_files_[file_], "read");

  RooUnfoldResponse* response;
  response = (RooUnfoldResponse*)_file0->Get("response_pT"); 
  
  // Use hist as input for the unfolding.
  RooUnfoldBayes unfold_bayes(response, hist_, iterations);
  // the original hist has become obsolete, transform it into the unfolded histogram.
  hist_ = (TH1D*) unfold_bayes.Hreco( RooUnfold::kCovariance ); 

  //== The following lines extract the unfolding matrix and plot it as a histogram.
  TCanvas * can_unfold = new TCanvas("unfolded", "unfodled", 1.);
  PrepareCanvas_2D( can_unfold, TString::Format("Unfolding_matrix_%i_iterations", iterations) );
  TMatrixD& unfold_m = (TMatrixD&)unfold_bayes.UnfoldingMatrix();
  TH2D* unfold_h = new TH2D( unfold_m );   
  double min_unfold = GetMinimumValue( unfold_h );
  double max_unfold = unfold_h->GetMaximum();
  unfold_h->GetZaxis()->SetRangeUser( min_unfold * 0.9, max_unfold * 1.1);
  Prepare_2Dplot( unfold_h );

  unfold_h->GetXaxis()->SetTitle("p_{T,det} [GeV]");
  unfold_h->GetYaxis()->SetTitle("p_{T,gen} [GeV]");
  can_unfold->SetLogz();

  //== Done drawing.


  //== This draw the covariance matrix.
  TCanvas * can_cov = new TCanvas("covariance", "covariance", 1.);
  PrepareCanvas_2D( can_cov, TString::Format("Covariance_unfolding_matrix_%i_iterations", iterations) );
  TMatrixD& covariance_m = (TMatrixD&)unfold_bayes.GetMeasuredCov();
  TH2D* covariance_h = new TH2D( covariance_m );   
  double min_covariance = GetMinimumValue( covariance_h );
  double max_covariance = covariance_h->GetMaximum();
  covariance_h->GetZaxis()->SetRangeUser( min_covariance * 0.9, max_covariance * 1.1);
  Prepare_2Dplot( covariance_h );

  covariance_h->GetXaxis()->SetTitle("E_{det} [GeV]");
  covariance_h->GetYaxis()->SetTitle("E_{gen} [GeV]");
  
  can_cov->SetLogz();


  }
}





void Unfolder::Plot_covariance(int file_, int iterations){

  if( file_ >= 0 ){
  TFile *_file0 = new TFile( MC_files_[file_], "read");

  TH1D* hist_;
  Get_DetEnergy( 0, hist_ );
  hist_->Scale(1./lumi_);
  RooUnfoldResponse* response;
  response = (RooUnfoldResponse*)_file0->Get("response"); 

  // Use hist as input for the unfolding.
  RooUnfoldBayes unfold_bayes(response, hist_, iterations);

 //== This draw the covariance matrix.
  TCanvas * can_cov = new TCanvas("covariance", "covariance", 1.);
  PrepareCanvas_2D( can_cov, TString::Format("Covariance_unfolding_matrix_%i_iterations", iterations) );
  TMatrixD covariance_m = (TMatrixD)unfold_bayes.Ereco( RooUnfold::kCovariance);
  TH2D* covariance_h = new TH2D( covariance_m );   
  double min_covariance = GetMinimumValue( covariance_h );
  double max_covariance = covariance_h->GetMaximum();
  covariance_h->GetZaxis()->SetRangeUser( min_covariance * 0.9, max_covariance * 1.1);
  Prepare_2Dplot( covariance_h );


  //== We need to change the bins of the covariance matrix. 
  //== We use the response matrix as template.
  TH2D* hTemplate = (TH2D*)response->Hresponse();
  ChangeBinEdges( covariance_h, covariance_h, hTemplate);

  Prepare_2Dplot( covariance_h );
  covariance_h->GetXaxis()->SetTitle("E [GeV]");
  covariance_h->GetYaxis()->SetTitle("E [GeV]");

  DrawWithNegativeLog( covariance_h, can_cov );
  
  can_cov->SetLogz();

  //Finish_canvas(can_cov, "rightish");

   TLatex *   tex = new TLatex(0.12,0.92,"Pythia6 (Z2*)");
   tex->SetTextAlign(13);
   tex->SetTextFont(43);
   tex->SetTextSize(35);
   tex->SetLineWidth(2);
//   tex->Draw();

   tex = new TLatex(0.80,0.935,"7 TeV");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.05454545);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.5622857,0.8500625,"CMS");
tex->SetNDC();
   tex->SetTextFont(61);
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.5622857,0.80,"Simulation");
tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetTextFont(52);
   tex->SetTextSize(0.0456);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.12,0.935,"anti-k_{t} (R=0.5) (-6.6<#eta<-5.2)");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->Draw();
  
  can_cov->SaveAs( folder_ + TString::Format("Covariance_UnfoldMatrix_%i.pdf", iterations) );
  can_cov->SaveAs( folder_ + TString::Format("Covariance_UnfoldMatrix_%i.C", iterations) );  


  //== Done with the covariance matrix.
 } 
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
  cout << "\n\t===Unfolder::Unfold-and-smear===\t" << MC_files_[ MC_ ] << "\tmethod\t" << method << endl;
  variable = "all";

  ofstream smearing_matrix;
  smearing_matrix.open("Smearing_matrix.txt");

  TH1D* theHist = (TH1D*)hist_->Clone("theHist");

  TFile *_file0 = new TFile( MC_files_[ MC_ ], "read");
  TFile *_file_det = new TFile( datafile_, "read");
  RooUnfoldResponse* response;
  
  response = (RooUnfoldResponse*)_file0->Get("response"); 

  //-- Extract histograms.
  TH1D* hMiss = (TH1D*)_file0->Get("hCastorJet_miss_all");
  TH1D* hFake 		= (TH1D*)response->Hfakes();
  TH2D* hResponse 	= (TH2D*)response->Hresponse(); 
  TH1D* hMeasured	= (TH1D*)response->Hmeasured();
  TH1D* hTruth		= (TH1D*)response->Htruth();
  TH1D* hSmear;

  // -- S = R * (U - M) + F

  // -- Construct the RooUnfolfBayes object to obtain the covariance matrix and the unfolded distribution.
  RooUnfoldBayes unfold_bayes(response, theHist, iterations); 
   
  unfold_bayes.SetVerbose(0);
  TH1D* hUnfold = (TH1D*) unfold_bayes.Hreco( RooUnfold::kCovariance );

  int total_events_nocuts, total_events_nocuts_data;

  //== Scale with the number measured jets in data and MC.
  double scale_fakes = theHist->Integral()/hMeasured->Integral();
  hFake->Scale( scale_fakes );

  hSmear = (TH1D*) response->ApplyToTruth( hUnfold );
  hSmear->Add( hFake );  
  // Method 1: calculation of chi² through toy model.

  if( method == 1 ){ 
    chi2 = Calculate_smearedBackError_covariance( theHist, hUnfold, response, iterations );
  }
  else{ chi2 = 0; }

  hist_ = hSmear;
  cout << "Done smearing" << endl;
}




//== Pass a custom MC file.
void Unfolder::Unfold_and_smear(int file_, TH1D* &hist_, TFile* MC_, int iterations, TString variable, int method, double &chi2 ){
  cout << "\n\t===Unfolder::Unfold-and-smear===\t" << MC_->GetName() << "\tmethod\t" << method << endl;
  variable = "all";

  ofstream smearing_matrix;
  smearing_matrix.open("Smearing_matrix.txt");

  TH1D* theHist = (TH1D*)hist_->Clone("theHist");

  TFile *_file0 = MC_; //new TFile( MC_files_[ MC_ ], "read");
  TFile *_file_det = new TFile( datafile_, "read");
  RooUnfoldResponse* response;
  
  response = (RooUnfoldResponse*)_file0->Get("response"); 

  //-- Extract histograms.
  TH1D* hMiss = (TH1D*)_file0->Get("hCastorJet_miss_all");
  TH1D* hFake = (TH1D*)_file0->Get("hCastorJet_fake_all");	
  TH2D* hResponse 	= (TH2D*)response->Hresponse(); 
  TH1D* hMeasured	= (TH1D*)response->Hmeasured();
  TH1D* hTruth		= (TH1D*)response->Htruth();
  TH1D* hSmear;

  // -- S = R * (U - M) + F

  // -- Construct the RooUnfolfBayes object to obtain the covariance matrix and the unfolded distribution.
  RooUnfoldBayes unfold_bayes(response, theHist, iterations); 
   
  unfold_bayes.SetVerbose(0);
  TH1D* hUnfold = (TH1D*) unfold_bayes.Hreco( RooUnfold::kCovariance );

  int total_events_nocuts, total_events_nocuts_data;

  //== Scale with the number measured jets in data and MC.
  double scale_fakes = theHist->Integral()/hMeasured->Integral();
  hFake->Scale( scale_fakes );

  hSmear = (TH1D*) response->ApplyToTruth( hUnfold );
  hSmear->Add( hFake );  
  // Method 1: calculation of chi² through toy model.

  if( method == 1 ){ 
    chi2 = Calculate_smearedBackError_covariance( theHist, hUnfold, response, iterations );
  }
  else{ chi2 = 0; }

  hist_ = hSmear;
  cout << "Done smearing" << endl;
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

void Unfolder::PrepareLegend( map<TString, TString> legend_info, map<TString, TString> printLabel, map<TString, TString> legend_info_gen ){

  legend_info_gen_ = legend_info_gen;
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

//  Hist_GenLevel(pad_abs_, true);
  Plot_Unfolded(pad_abs_, variable, iterations);
  Plot_Unfolded_Ratio(pad_ratio_, variable, iterations);

  can_->SaveAs( TString::Format( folder_ + "/UnfoldedEnergy_Data_and_MC_" + variable + "_" + label_ + "_%iterations.pdf", iterations) );
  can_->SaveAs( TString::Format( folder_ + "/UnfoldedEnergy_Data_and_MC_" + variable + "_" + label_ + "_%iterations.C", iterations) );
}



void Unfolder::DoublePaddedComparison_statistics(TString variable, int iterations){
  cout << "\nUnfolder::DoublePaddedComparison_statistics" << endl;

  TCanvas* can_;
  TPad* pad_abs_, *pad_ratio_;
  PrepareCanvas(can_, "Comparison_UnfoldedEnergy_Data_and_MC");

  SplitCanvas(can_, pad_abs_, pad_ratio_);

  Hist_GenLevel(pad_abs_, true);
  Plot_Unfolded(pad_abs_, variable, iterations);
  Plot_Unfolded_Ratio_statistics(pad_ratio_, variable, iterations);

  can_->SaveAs( TString::Format( folder_ + "/UnfoldedEnergy_statistics_Data_and_MC_" + variable + "_" + label_ + "_%iterations.pdf", iterations) );
  can_->SaveAs( TString::Format( folder_ + "/UnfoldedEnergy_statistics_Data_and_MC_" + variable + "_" + label_ + "_%iterations.C", iterations) );
}




void Unfolder::DoublePaddedComparison_modelDependence(TString variable, int iterations){
  cout << "\nUnfolder::DoublePaddedComparison_modelDependence" << endl;

  TCanvas* can_;
  TPad* pad_abs_, *pad_ratio_;
  PrepareCanvas(can_, "Comparison_UnfoldedEnergy_Data_and_MC");

  scaletodata_ = false;

  SplitCanvas(can_, pad_abs_, pad_ratio_);

//  Hist_GenLevel(pad_abs_, true);
  Plot_Unfolded_modelDependence(pad_abs_, variable, iterations);
  Plot_Unfolded_Ratio_modelDependence(pad_ratio_, variable, iterations);

  can_->SaveAs( TString::Format( folder_ + "/Systematics_UnfoldedEnergy_modelDependence_%iterations_deltaphi_0%i_etawidth_0%i" + cutsetup_ + ".pdf", iterations, static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_)) );
  can_->SaveAs( TString::Format( folder_ + "/Systematics_UnfoldedEnergy_modelDependence_%iterations_deltaphi_0%i_etawidth_0%i" + cutsetup_ + ".C", iterations, static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_)) );
}




void Unfolder::DoublePaddedComparison_positionDependence(TString variable, int iterations){
  cout << "\nUnfolder::DoublePaddedComparison_statistics" << endl;

  TCanvas* can_;
  TPad* pad_abs_, *pad_ratio_;
  PrepareCanvas(can_, "Comparison_UnfoldedEnergy_Data_and_MC");

  SplitCanvas(can_, pad_abs_, pad_ratio_);

//  Hist_GenLevel(pad_abs_, true);
  Plot_Unfolded_positionDependence(pad_abs_, variable, iterations);
  Plot_Unfolded_Ratio_positionDependence(pad_ratio_, variable, iterations);

  can_->SaveAs( TString::Format( folder_ + "/Systematics_UnfoldedEnergy_positionDependence_%iterations_deltaphi_0%i_etawidth_0%i" + cutsetup_ + ".pdf", iterations, static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_)) );
  can_->SaveAs( TString::Format( folder_ + "/Systematics_UnfoldedEnergy_positionDependence_%iterations_deltaphi_0%i_etawidth_0%i" + cutsetup_ + ".C", iterations, static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_)) );
}






void Unfolder::DoublePaddedComparison_JESDependence(TString variable, int iterations){
  cout << "\nUnfolder::DoublePaddedComparison_JESDependence" << endl;

  TCanvas* can_;
  TPad* pad_abs_, *pad_ratio_;
  PrepareCanvas(can_, "Comparison_UnfoldedEnergy_Data_and_MC");

  SplitCanvas(can_, pad_abs_, pad_ratio_);

//  Hist_GenLevel(pad_abs_, true);
  Plot_Unfolded_JESDependence(pad_abs_, variable, iterations);
  Plot_Unfolded_Ratio_JESDependence(pad_ratio_, variable, iterations);

  can_->SaveAs( TString::Format( folder_ + "/Compare_JES_%iterations_deltaphi_0%i_etawidth_0%i" + cutsetup_ + ".pdf", iterations, static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_)) );
  can_->SaveAs( TString::Format( folder_ + "/Compare_JES_%iterations_deltaphi_0%i_etawidth_0%i" + cutsetup_ + ".C", iterations, static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_)) );
}




void Unfolder::DoublePaddedComparison_JESmeasured(TString variable, int iterations){
  cout << "\nUnfolder::DoublePaddedComparison_JESDependence" << endl;

  TCanvas* can_;
  TPad* pad_abs_, *pad_ratio_;
  PrepareCanvas(can_, "Comparison_UnfoldedEnergy_Data_and_MC");

  SplitCanvas(can_, pad_abs_, pad_ratio_);

  Plot_measured_JESDependence(pad_abs_, variable, iterations);
  Plot_measured_Ratio_JESDependence(pad_ratio_, variable, iterations);

  can_->SaveAs( TString::Format( folder_ + "/Compare_JESmeasured_%iterations_deltaphi_0%i_etawidth_0%i" + cutsetup_ + ".pdf", iterations, static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_)) );
  can_->SaveAs( TString::Format( folder_ + "/Compare_JESmeasured_%iterations_deltaphi_0%i_etawidth_0%i" + cutsetup_ + ".C", iterations, static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_)) );
}






//
// -- Absolute distributions.
//


void Unfolder::Plot_Unfolded(TPad* & pad_, TString variable, int iterations){
  pad_->cd();

  cout << "Unfolder::Plot_Unfolded\t" << iterations << "\titerations" << endl;

  // Set up legend.
  TLegend *leg = new TLegend(0.6, 0.50, 0.95, 0.95);

  int color_ = 1;
  double min_val, max_val;

  TString drawoptions = "hist";
  int plotted_hist = 0;

  TCanvas* can_unfolded;
  PrepareCanvas(can_unfolded, "can_unfolded");

  for(int MC_ = 0; MC_ < 1; MC_++){ //MC_files_.size(); MC_++){
    cout << "Unfolder::Plot_Unfolded\tMC file \t" << MC_ << endl;

    //-- Unfold det. (Data)
    TH1D* hData, *hFirst;

    if( variable == "all") { Get_DetEnergy(-1, hData); }
    if( variable == "lead"){ Get_DetEnergy_lead(-1, hData); }


    for(int it = 1; it <= iterations; it++){ 

      Get_DetUnfolded( MC_, hData, it, variable);
      if( it == 0){ First_Plot( hFirst, hData, plotted_hist++, min_val, max_val); }
    
      hData->SetLineColor( getColor( MC_+1) );
      hData->SetMarkerColor( getColor(MC_+1) );
      hData->SetMarkerStyle( 24 + MC_);
      leg->AddEntry( hData, TString::Format("Data %i it.", iterations), "p");

      pad_->cd();
      hData->     DrawClone(drawoptions);

      can_unfolded->cd();
      hData->	DrawClone(drawoptions);

      drawoptions = "histsame";

      cout << "\tUnfolder::Plot_Unfolded\t" << MC_ << "\t" << hData->Integral() << endl;
  
      First_Plot( hFirst, hData, plotted_hist++, min_val, max_val);
    }
  }
 
  cout << "Unfolder::Plot_Unfolded - The End\t" << endl;

  leg->Draw();

  pad_->SetLogy();
  pad_->Update();

  can_unfolded->SetLogy();
  can_unfolded->Update();
  can_unfolded->SaveAs(folder_ + "Unfolded_data_30_iterations.C");
  can_unfolded->SaveAs(folder_ + "Unfolded_data_30_iterations.pdf");
}

//
// -- Relative distributions.
//

void Unfolder::Plot_Unfolded_Ratio(TPad* & pad_, TString variable, int iterations){
  pad_->cd();

  cout << "Unfolder::Plot_Unfolded_Ratio\t" << iterations << "\titerations" << endl;

  TH1D* hFirst;
  TString drawoptions = "hist";

  for(int MC_ = 0; MC_ < 1 ; MC_++){
    //-- Generator level.

    TH1D* hGen, *hDet;
    if( variable == "all") { Get_GenEnergy(MC_, hGen); Get_DetEnergy( MC_, hDet); }
    if( variable == "lead") { Get_GenEnergy_lead(MC_, hGen); }

    hGen->SetLineColor( getColor( MC_+1 ) );
    hGen->SetLineStyle( 1 );
    hGen->SetLineWidth( 2 );
    hGen->GetYaxis()->SetTitle("Ratio");

    TString scalefactor_MC = (set_of_tags_[ "scalefactors" ])[ MC_files_[MC_] ];
    double scalefactors_MC_ = scalefactor_MC.Atof() ;
  SetCastorJetEnergy_norm( scalefactors_MC_ / scalefactors_Data_ ); 

//    hGen->Scale( 1./renorm_ );
    hGen->Scale( 3425279.  / hDet->Integral() );
    cout << "\tUnfolder::Plot_Unfolded_Ratio\t" << MC_ << "\t" << hGen->Integral() << "\t" << scalefactors_MC_ << endl;

    if( MC_ == 0){ 
      hFirst = (TH1D*)hGen->Clone("Reference_histogram"); 
      hFirst->Draw("hist");
    }

    hGen->Divide(hFirst);
    hGen-> Draw( "histsame" );
    drawoptions = "histsame";

    // Determine the normalisation of the detector level distribution: MC(Gen)/MC(Det)
    // For this we need the unscaled MC distributions.

    /*
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
    */

    //-- Data unfolded.
    TH1D* hData;

    if( variable == "all") { Get_DetEnergy( -1 , hData); }
    if( variable == "lead"){ Get_DetEnergy_lead( -1 , hData); }


    Get_DetUnfolded( MC_ , hData, iterations, variable);
    hData->SetLineColor( getColor( MC_+1) );
    hData->SetMarkerColor( getColor( MC_+1 ) );
    hData->SetMarkerStyle( 24 + MC_ );
    cout << "\tUnfolder::Plot_Unfolded_Ratio\t" << MC_ << "\t" << hData->Integral() << endl;
    hData->Divide( hFirst);
    hData->       DrawClone("edatasame");

  }

  hFirst->GetYaxis()->SetRangeUser(0., 3.);
  pad_->Update();
}







void Unfolder::Plot_Unfolded_Ratio_statistics(TPad* & pad_, TString variable, int iterations){
  pad_->cd();

  cout << "Unfolder::Plot_Unfolded_Ratio_statistics\t" << iterations << "\titerations" << endl;

  TH1D* hFirst;
  TString drawoptions = "ephist";

  for(int MC_ = 0; MC_ < MC_files_.size() ; MC_++){
    //-- Generator level.


    //-- Data unfolded.
    TH1D* hData;

    if( variable == "all") { Get_DetEnergy( -1 , hData); }
    if( variable == "lead"){ Get_DetEnergy_lead( -1 , hData); }

    Get_DetUnfolded( MC_ , hData, iterations, variable);

    if( MC_ == 0){ 
      hFirst = (TH1D*)hData->Clone("Reference_histogram"); 
      hFirst->GetYaxis()->SetTitle("Ratio");
    }

    hData->GetYaxis()->SetTitle("Ratio");
    hData->SetLineColor( getColor( MC_+1) );
    hData->SetMarkerColor( getColor( MC_+1 ) );
    hData->SetMarkerStyle( 24 + MC_ );
    cout << "\tUnfolder::Plot_Unfolded_Ratio\t" << MC_ << "\t" << hData->Integral() << endl;

    hData->Divide( hFirst);
    cout << "\t\t\t" << hFirst->Integral() << endl;
    hData->       DrawClone( drawoptions );
    drawoptions = "edatasame";

  }

  hFirst->GetYaxis()->SetRangeUser(0., 3.);
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

  int iterations;

  TH1D* hist_;
  TString drawoptions = "hist", legendoptions = "l";
  TCanvas *can; 
  PrepareCanvas( can, "Systematics_Generator_Level_Energy" );
  vector<TH1D*> histos;
  TLegend *leg = new TLegend(0.65, 0.7, 0.95, 0.95);
  bool done_data = false;

  double Det_to_gen = Det_to_gen_scale(0);

  for( int file = 0; file < MC_files_.size(); file++ ){

    if( ( (set_of_tags_["mc_type"])[MC_files_[file]] == "shift_MPI_or_Tune" ) ) { 
	cout << "Bad type\t" << MC_files_[file]	<< endl;
	continue; }

/*
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
*/


    int file_ = file;
    Get_DetEnergy( -1 , hist_ );
    Get_DetUnfolded( file, hist_, iterations);
//    Get_GenEnergy(file_, hist_);    

    hist_->Scale( 1./hist_->Integral() );

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
  cout << "Done loop" << endl;

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
  can->SaveAs(folder_ + "Unfolding_data_" + variable + label_ + ".C");
  can->SaveAs(folder_ + "Unfolding_data_" + variable + label_ + ".pdf");
}








void Unfolder::ClosureTest_data(TString variable, TString file, int method){
  //== Method is now an obsolete parameter.

  cout << "\n\n\n\t************************************" << endl;
  cout << "\t************************************" << endl;
  cout << "\t* ClosureTest_data******************" << endl;
  cout << "\t************************************" << endl;
  cout << "\t************************************" << endl;

  ofstream iterations_and_errors;
  iterations_and_errors.open( TString::Format("Iterations_and_errors_%i_%i.txt", static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) ) );

  // Prepare variables and objects.
  TH1D* hist_original;
  TH1D* hist_result; 
  TH1D* hist_reference, *hist_previous; 
  TString htitle;
  int iterations_ = 30, iterations_start = 1;

  double ticksize = gStyle->GetTickLength();
  
  int file_to_unfold = -1;

  TString norm = "notNorm";

  vector<TH1D*> histos_unfolded;
   
  // Get the detector level distribution.
  Get_DetEnergy(file_to_unfold, hist_original );
  htitle = "E_{det} spectrum (all)";	

  // Hist_reference is the same as the original, but with modified lower edge.
  GetSubHistogram( hist_original, hist_reference, Eplotmin_, 2100.);

  hError_hist = new TH2D("hError_hist",  "hError_hist;E;iteration", 
		hist_reference->GetNbinsX(), 
		hist_reference->GetBinLowEdge( 1 ), 
		hist_reference->GetXaxis()->GetBinUpEdge(  hist_reference->GetNbinsX() ), 
		iterations_, 
		0., iterations_ );

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
	hist_reference->GetXaxis()->SetTitle("E [GeV]");
	hist_reference->GetYaxis()->SetTitle("#frac{dN}{dE}");

  hist_reference->Scale( 1./ (1.2 * 1E5) );  
  SetDnDx( hist_reference );

  // Loop over number of Bayesian iterations, and unfold->smear detector level distribution.
  // Calculate chi2 as we go.
  double chi2 = hist_original->Chi2Test( hist_reference, "CHI2/NDF");

  vector<double> etaband;
  map<double, int> eta_markers;  
      etaband.push_back(0.0);
      etaband.push_back(0.2);      
      etaband.push_back(0.5);      

  vector<double> deltaPhiMax;
      deltaPhiMax.push_back(0.2);  	
      deltaPhiMax.push_back(0.5);

  //==============================//
  //== Prepare canvas for chi2. ==//
  //==============================//

  TCanvas *can_chi2;
  PrepareCanvas( can_chi2, "CHI2_Test_" + variable);
  TString drawoptiongraph = "AP";
  
  TLegend *leg_chi2 = new TLegend( 0.55, 0.55, 1. - can_chi2->GetRightMargin()-ticksize , 1. - can_chi2->GetTopMargin()-ticksize);
  leg_chi2->SetFillStyle(0);
  leg_chi2->SetBorderSize(0);


  //---------------------------------------------------------------//
  //-- Loop over the files, if there are several to be unfolded. --//
  //---------------------------------------------------------------//

  int file_ = 1;
  for(int _eta = 0; _eta < etaband.size(); _eta++){
     for(int _phi = 0; _phi < deltaPhiMax.size(); _phi++, file_++){
     
       vector<TH1D*> histos;     
       histos.push_back( hist_reference );
       histos_unfolded.push_back(hist_reference);     
          
       int current_color = 1;
  
       // Plot the smeared distributions and their ratios.
       TCanvas *can;
       PrepareCanvas( can, "ClosureTest_data" + variable);
    
       TPad* dummyPad1, *dummyPad2; SplitCanvas( can, dummyPad1, dummyPad2);
    
       TLegend *leg = new TLegend(0.55, 0.45, 1. - dummyPad1->GetRightMargin()-ticksize , 1. - dummyPad1->GetTopMargin()-ticksize);
         leg->SetFillStyle(0);    
         leg->SetBorderSize(0); 
         
       leg->AddEntry( hist_reference, "Actual data", "p");
       double xaxisgraph[iterations_], yaxisgraph[iterations_], chi2diff[iterations_];
         int actual_iterations = 0;
       //--------------------------------------------------//
       //-- Loop over the number of Bayesian iterations. --//
       //--------------------------------------------------//

       cout << "=== Bayesian iterations ===" << endl;
       int increase_iterations = 1;
       increase_iterations = 3;

       for(int iterations = iterations_start; iterations <= iterations_; iterations+=increase_iterations){
         //TFile *_file0 = new TFile( MC_files_[file_], "read");
         TString MCfilename = TString::Format("/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo" + cutsetup_ + "_unfold_Emin_150.000000_deltaPhiMax_%f_etaband_%f_all_matchE.root", deltaPhiMax[_phi],  etaband[_eta] );
         TFile* _file0 = new TFile( MCfilename  , "READ");
         cout << "Filename\t" << MCfilename << "\t" << endl;
         

         //== Determine the inrease in Bayesian iterations.
         if(iterations >= 10 ){ increase_iterations = 3; }
         if(iterations >= 20 ){ increase_iterations = 3; }
         if(iterations >= 100 ){ increase_iterations = 1; }    

         //== Prepare the histogram to send away.
         TString hist_name =  TString::Format("Closure_%i_iterations_phidiff_0%i_etaband_0%i", 
		iterations, 
	   	static_cast<int>(10.* deltaPhiMax[_phi] ), 
		static_cast<int>(10. * etaband[_eta]) );

         hist_result = (TH1D*)hist_original->Clone( hist_name );
         hist_result->SetName( hist_name );
         hist_result->SetTitle( hist_name );	
   
         //== Send the histogram away for unfolding.
         double chi2 =0.;
         Unfold_and_smear( file_to_unfold, hist_result, _file0, iterations, variable, method, chi2 );
         GetSubHistogram( hist_result, hist_result, Eplotmin_, 2100.);
   
         //== Set colors.
         hist_result->SetLineColor( getColor( current_color ) );
         hist_result->SetLineStyle( current_color  );
         hist_result->SetLineWidth( 2 );
         hist_result->SetMarkerColor(  getColor( current_color ) );
         leg->AddEntry( hist_result, TString::Format("%i it., #Delta#varphi = 0.%i, #eta = 0.%i", iterations, static_cast<int>(10. * deltaPhiMax[_phi]), static_cast<int>(10. * etaband[_eta]) ), "l");

         //== Scale and divide by binsize/
         SetDnDx( hist_result );
         hist_result->Scale( 1./ (1.2 * 1E5) );  

         //== Store.
         histos.push_back( hist_result );
   
         //-- Chi2 test.
         xaxisgraph[actual_iterations] = iterations;
         yaxisgraph[actual_iterations] = chi2; 

         iterations_and_errors << "iterations\t" << iterations << "\t" << chi2 << endl;

         current_color++;
         actual_iterations++;
       }
       // Make sure to set histogram axes.
       hist_original->GetXaxis()->SetTitle("E_{det} [GeV]");
       hist_original->GetYaxis()->SetTitle("#frac{1}{N}.#frac{dN}{dE}");
      
       MakeDoublePaddedComparison(can, histos, leg );
       
       can->SaveAs( TString::Format(folder_ + "ClosureTest_data_" + variable + "_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i_" + norm + ".C", 
   	   iterations_, static_cast<int>(Ethresh_), static_cast<int>(10. * deltaPhiMax[_phi]), static_cast<int>(10. * etaband[_eta]) ) );
       can->SaveAs( TString::Format(folder_ + "ClosureTest_data_" + variable + "_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i_" + norm + ".pdf", 
	   iterations_, static_cast<int>(Ethresh_), static_cast<int>(10. * deltaPhiMax[_phi]), static_cast<int>(10. * etaband[_eta]) ) );

       // Finish chi2 study.
       can_chi2->cd();
       
       TGraph* chi2_evolution = new TGraph(actual_iterations, xaxisgraph, yaxisgraph);
       chi2_evolution->GetXaxis()->SetTitle("N_{it.}");
       chi2_evolution->GetYaxis()->SetTitle("#chi^{2}/NDF");
       chi2_evolution->GetYaxis()->SetRangeUser(0.5, 5000.);
       chi2_evolution->SetName( TString::Format("CHI2_deltaPhiMax_0%i_etaband_0%i", static_cast<int>(10. * deltaPhiMax[_phi]), static_cast<int>(10. * etaband[_eta]) ) );
       int color_ = file_;
       chi2_evolution->SetMarkerColor( getColor( color_ ) );
       chi2_evolution->SetMarkerStyle( color_ + 20 );
       chi2_evolution->SetMarkerSize( 1.25 );
       Prepare_1Dplot( chi2_evolution );
       chi2_evolution->DrawClone( drawoptiongraph );
       drawoptiongraph = "psame";
       leg_chi2->AddEntry( chi2_evolution, TString::Format("#Delta#varphi = 0.%i, #eta = 0.%i", static_cast<int>(10. * deltaPhiMax[_phi]), static_cast<int>(10. * etaband[_eta])) , "p");
     } //== Loop over phi.
  } //== Loop over eta

  iterations_and_errors.close();

  TCanvas* canerror = new TCanvas("can", "can", 1.);
  hError_hist->Draw("colz");
  canerror->SaveAs("error_plot.C");
  
  //== Save chi2.  
     
  can_chi2->cd();
  can_chi2->SetLogy();
  leg_chi2->Draw();
  can_chi2->SaveAs( TString::Format(folder_ + "Chi2_Test_data_" + variable + "_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i_" + norm + ".C", 
    iterations_, static_cast<int>(Ethresh_), static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) ) );
  can_chi2->SaveAs( TString::Format(folder_ + "Chi2_Test_data_" + variable + "_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i_" + norm + ".pdf",
    iterations_, static_cast<int>(Ethresh_), static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) ) );    
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
  GetSubHistogram( hist_original, hist_reference, Eplotmin_, 2100.);

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

      GetSubHistogram( hist_result, hist_result, Eplotmin_, 2100.);

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
  cout << "===\tUnfolder::MakeDoublePaddedComparison\t===" << endl;

  TPad* pad_abs_, *pad_ratio_;
  PrepareCanvas(can_, "Comparison_UnfoldedEnergy_Data_and_MC");

  SplitCanvas(can_, pad_abs_, pad_ratio_);

  Plot_Absolute(pad_abs_, histos);
  Plot_Ratio(pad_ratio_, histos);

  pad_abs_->cd();
  legend->Draw();
  pad_abs_->Update();
}



void Unfolder::MakeDoublePaddedComparison(TCanvas * &can_, TPad* &pad_abs_, TPad* &pad_ratio_, vector<TH1D*> histos, TLegend *legend){
  cout << "===\tUnfolder::MakeDoublePaddedComparison\t===" << endl;


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
  cout << "----- Plot Absolute Value ---\t\tPlot_Absolute( TPad* & pad_, vector<TH1D*> histos)" << endl;

  double min_val, max_val;
  TH1D* hDistr, *hFirst;
  bool firstDistr = true, passed_data = false;
  int nDistr = 0;
  TString drawoptions = "pehist";
  TString legendoptions = "l";

  TH1D* hData = histos[0];  
  for(int bin = 1; bin <=hData->GetNbinsX(); bin++){
    cout << "%%%\tBin\t" << bin << "\t" << hData->GetBinError( bin ) << "\t" << hData->GetBinContent( bin ) << "\t" << hData->GetBinError( bin )/hData->GetBinContent( bin ) << endl;

  }

  for(int hist_ = 0; hist_ < histos.size(); hist_++){

    TH1D* hDistr = (TH1D*)histos[hist_]->Clone( TString::Format("Absolute_%i", hist_) );

    Rebin_to( hDistr, 50);
    hDistr->GetYaxis()->SetLabelFont(42);
    hDistr->GetYaxis()->SetTitleOffset( 0.85 );
    hDistr->GetYaxis()->SetTitleSize(0.09);
    hDistr->GetYaxis()->SetLabelSize(0.07);
    hDistr->GetXaxis()->SetNdivisions( 504 );
    hDistr->SetName( TString::Format("Absolute_%i", hist_) );

    if( firstDistr ){ 
      hFirst = (TH1D*)hDistr->Clone("First_absolute");
      hDistr = hFirst;
      firstDistr = false;
      Prepare_1Dplot( hFirst );
    }

    First_Plot( hFirst, hDistr, nDistr, max_val, min_val);
 
    pad_->cd();
    pad_->SetLogy();

    Prepare_1Dplot( hDistr );

    hDistr->Draw( drawoptions );
   
    drawoptions = "histsame";
    nDistr++;
    
    if( nDistr > 0 ) pad_->Update();
  }
}




void Unfolder::Plot_Ratio( TPad* & pad_, vector<TH1D*> histos){
  cout << "----- Plot Relative Value ---\t\tPlot_Ratio( TPad* & pad_, vector<TH1D*> histos)" << endl;

  double min_val, max_val;
  TH1D* hDistr, *hFirst, *hFirst_abs;
  bool firstDistr = true, passed_data = false;
  int nDistr = 0;
  TString drawoptions = "pehist";

  for(int hist_ = 0; hist_ < histos.size(); hist_++){

    hDistr = (TH1D*) histos[hist_]->Clone( TString::Format("Ratio_%i", hist_) );

    Rebin_to( hDistr, 50);
    hDistr->GetXaxis()->SetNdivisions( 510 );
    hDistr->GetYaxis()->SetTitle("Ratio");

    hDistr->GetYaxis()->SetLabelSize( (hDistr->GetYaxis()->GetLabelSize())*1.5 );
    hDistr->GetYaxis()->SetLabelOffset( (hDistr->GetYaxis()->GetLabelOffset())*1.5 );
    hDistr->GetYaxis()->SetTitleSize( hDistr->GetYaxis()->GetTitleSize()*1.5 );
    hDistr->GetYaxis()->SetNdivisions(205);

    hDistr->GetXaxis()->SetLabelSize(
      hDistr->GetXaxis()->GetLabelSize() *
	( double(pad_->GetWw()) / double(pad_->GetWh()) ) * 1.5
    ); 
    cout << "***\n"
	 << "***\t" << pad_->GetWw() << "\t" << pad_->GetWh() << "\t" << pad_->GetWw() / pad_->GetWh() << endl;
	


    if( firstDistr ){
      hFirst = (TH1D*)hDistr->Clone("First_ratio");
      hFirst_abs = (TH1D*)hFirst->Clone("_first_absolute");
      hFirst->Divide( hFirst );
      hFirst->GetYaxis()->SetTitle("Ratio");
      firstDistr = false;
      hFirst->GetYaxis()->SetLabelFont(42);
      hFirst->GetYaxis()->SetLabelSize( (hFirst->GetYaxis()->GetLabelSize())*1.5 );
      hFirst->GetYaxis()->SetLabelSize( (hFirst->GetYaxis()->GetLabelSize())*1.5 );
      hFirst->GetYaxis()->SetLabelOffset( hFirst->GetYaxis()->GetLabelOffset() );

      hFirst->GetXaxis()->SetLabelFont(42);
      hFirst->GetXaxis()->SetLabelSize( ( hFirst->GetXaxis()->GetLabelSize())*1.5 );
      hFirst->GetXaxis()->SetTitleSize( ( hFirst->GetXaxis()->GetTitleSize())*1.5 );
      hFirst->GetXaxis()->SetLabelOffset( hFirst->GetXaxis()->GetLabelOffset()*1.5 );
    }

    hDistr->Divide( hFirst_abs );
    First_Plot( hFirst, hDistr, nDistr, max_val, min_val);

    Prepare_1Dplot_ratio( hDistr );

    pad_->cd();
    hDistr->Draw( drawoptions );
    hFirst->GetYaxis()->SetRangeUser(0.75, 1.25);

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


void Unfolder::Calibration_SetParLimits(TF1*& analytical){

  analytical = new TF1("Analytical", "  ( [0] + [1] * log( [2] + x) )", fitting_threshold_, 900.);
  analytical->SetParLimits(2, -80., 1000.);

}



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
  int event_threshold_ = 20.;


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

  TF1* analytical;
  Calibration_SetParLimits( analytical);

  /* Open files. */
  for( int bin_phi = first_sector; bin_phi <=last_sector; bin_phi++){//  hResponse->GetNbinsZ(); bin_phi++){

    for(int file_ = 0; file_ < MC_files_.size(); file_++){

      TFile *_file0 = TFile::Open( MC_files_[file_], "Read");
      TString label = TString::Format( printLabel_[ MC_files_[file_] ] + "_sector_%i_to_%i", first_sector, last_sector);

      label = folder_ + label;
      TString label_short = "Test_this";
      int	new_dir = mkdir( (folder_).Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  );
      new_dir = mkdir( (folder_).Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  );

      TH2D* hGenE_selection;
      TH2D* hDetE_selection;
      TH2D* hResponse_selection;
      TH2D* hInverseResponse_selection;
      TH2D* hDetE;
      TH2D* hGenE;
      TH1D* hGen;
      TString setup = "good_sectors";

      if( setup_calibration_ == "all_sectors") 	setup = "all_sectors";
      if( setup_calibration_ == "bad_sectors") 	setup = "bad_sectors";
      else if( setup_calibration_ == "separate_sectors") setup = TString::Format("%i", bin_phi); 

       // -- Extract the response, Edet and Egen distributions from the files.
      Extract_2D_energy_distributions(file_, hGenE_selection, hDetE_selection, hResponse_selection, hInverseResponse_selection, hDetE, hGenE, hGen, setup);
            	 
      TGraphErrors * Response_true;
      TGraphErrors * Response_meas;

      cout << "CALIBRATIONFUNCTION - Analyze_response\tfile\t" << file_ << endl;

      Analyze_response( hResponse_selection, hGenE_selection, hDetE_selection, hGenE, label, hGen, file_, Response_meas, Response_true);

      cout << "CALIBRATIONFUNCTION - Analyze_response done" << endl;

      TString legend_info;
      if( setup_calibration_ == "all_sectors") { legend_info = "All sectors"; }
      else if( setup_calibration_ == "good_sectors"){ legend_info = "Good sectors"; }
      else if( setup_calibration_ == "bad_sectors"){ legend_info = "Bad sectors"; }
      else if( setup_calibration_ == "separate_sectors"){ legend_info = TString::Format("Sector %i", CorrectSectorNumber(bin_phi)); }
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
      
      if( ((setup == "good_sectors" || setup == "all_sectors") && bin_phi == 1) || (setup == "bad_sectors" )){
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
  cout << "Unfolder::CalibrationFunction_sectors(int " << first_sector << ", int " << last_sector << ", int " << which_file << ", TGraphErrors* &gre)" << endl;

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
  TF1* analytical;
  
 /*= new TF1("Analytical", "  [0] + [1] * log( [2] + x)  ", lowerbound_fit, 900.);
//    analytical->SetParLimits(2, 1.-lowerbound_fit, 1000.);
    analytical->SetParLimits(2, -80., 1000.);
//    analytical->SetParLimits(1, -1000., 0.);
//    analytical->SetParLimits(0, 0., 1000.); 
*/
  Calibration_SetParLimits( analytical);


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

    if( which_file == -1){
      first_file = 0;
      last_file = MC_files_.size();
    }
    else{
      first_file = which_file;
      last_file = which_file;
    }


//    for(int file_ = 0; file_ < MC_files_.size(); file_++){
    for(int file_ = which_file; file_ <= which_file; file_++){
      cout << "Which_file\t" << file_ << endl;

      TFile *_file0 = TFile::Open( MC_files_[file_], "Read");
       
      TCanvas *c = new TCanvas("can", "can",  1000, 1000);
        c->SetLeftMargin(0.20);
        c->SetRightMargin(0.18);
        c->SetBottomMargin(0.20);  

      TH2D* hGenE_selection;
      TH2D* hDetE_selection;
      TH2D* hResponse_selection;
      TH2D* hInverseResponse_selection;
      TH2D* hDetE;
      TH2D* hGenE;
      TH1D* hGen;
      TString setup = TString::Format("%i", sector);

       // -- Extract the response, Edet and Egen distributions from the files.
      Extract_2D_energy_distributions(file_, hGenE_selection, hDetE_selection, hResponse_selection, hInverseResponse_selection, hDetE, hGenE, hGen, setup);

      TGraphErrors * Response_true;
      TGraphErrors * Response_meas;

      Analyze_response( hResponse_selection, hGenE_selection, hDetE_selection, hGenE, label, hGen, file_, Response_meas, Response_true);

      Response_meas->SetName( TString::Format("%s_" + printLabel_[ MC_files_[file_] ], Response_meas->GetName()) );
      Response_true->SetName( TString::Format("%c_" + printLabel_[ MC_files_[file_] ], Response_true->GetName()) );
      Response_meas->SetLineWidth(2);	
	  	
      // cout << "Names\t" << Response_meas->GetName() << "\t\t" << Response_true->GetName() << endl;

      legend_det->AddEntry( Response_meas, TString::Format(legend_info_[ MC_files_[file_] ] + " (Sector %i)", CorrectSectorNumber(sector)), "lp"); 
 
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

    TF1* analytical ;

/*= new TF1("Analytical", " ( [0] + [1] * log( [2] + x) )", 100., 900.);
    analytical->SetParLimits(2, -80., 1000.); 
//    analytical->SetParLimits(1, -1000., 0.);
//    analytical->SetParLimits(0, 0., 1000.);  */
  Calibration_SetParLimits( analytical);

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
      plots_label = TString::Format("sector_%i", CorrectSectorNumber(sector));
      legend_info = TString::Format("Sector %i", CorrectSectorNumber(sector));
    }
    else if( setup == "bad_sectors"){
      plots_label = "bad_sectors";
      legend_info = "Bad sectors";
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

   TLegend* legend_graph = new TLegend( can_graph->GetLeftMargin(), 0.60, 0.45, 0.75);
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

    for(int file = 0; file < MC_files_.size(); file++){

      can->cd();

      // -- Extract the calibration function's parameters.
      //== The calibration parameters have been stored in a global variable. 
      double alpha = (calibration_parameters_[ MC_files_[file] ])[0];
      double beta = (calibration_parameters_[ MC_files_[file] ])[1];
      double gamma = (calibration_parameters_[ MC_files_[file] ])[2];
 
      TString calibration_tag = ""; 
      if( MC_files_[file].Contains("ShowerLibrary")){ calibration_tag = "Data"; }
      if( MC_files_[file].Contains("FullSimulation")){ calibration_tag = "MC"; }

      //== Print out fit values.
        ofstream calibrating_values;
        calibrating_values.open( "Calibrating_" + setup + "_values.txt", ios::out | ios::app | ios::binary);

        calibrating_values << "double Energy_" << calibration_tag << "_" << sector << "( double edet ){" << endl;
        calibrating_values << "  double alpha = " << alpha << ";" << endl;
        calibrating_values << "  double beta = "  << beta << ";" << endl;
        calibrating_values << "  double gamma = " << gamma << ";" << endl;
        calibrating_values << "return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; \n}\n" << endl; 

	calibrating_values.close();   
      //== Printed out fit values.


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
      calibration_histogram->GetYaxis()->SetTitle("1/<E_{det}/E_{gen}>");
      pad1->Update();

      pad2->cd();
      calibration_histogram->Divide( original_calibration );
      calibration_histogram->GetYaxis()->SetTitle( "Syst./MC");
      calibration_histogram->GetYaxis()->SetRangeUser(0.75,1.2);

      pad1->cd();
      legend_det->Draw();

      //== Save canvas.
      can->SaveAs(TString::Format( folder_ + "/Systematics_function_" + plots_label + "_%iGeV.C", static_cast<int> (Ethresh_ ) ) );
      can->SaveAs(TString::Format( folder_ + "/Systematics_function_" + plots_label + "_%iGeV.pdf", static_cast<int> (Ethresh_ ) ) );

      can_graph->cd();
      TGraphErrors* graph_errs = calibration_graphs_[ MC_files_[file] ];
      graph_errs->SetMarkerSize( 2. );
      Prepare_1Dplot( graph_errs ) ;

      graph_errs->Draw(drawoptions_graph);
      graph_errs->GetHistogram()->GetYaxis()->SetRangeUser(0., 3.5);
      graph_errs->GetHistogram()->GetYaxis()->SetTitle("1/<E_{det}/E_{gen}>");
      drawoptions_graph = "pesame";

      legend_graph->AddEntry( graph_errs, legend_info_[ MC_files_[file]] + " (" + legend_info + ")", "p");

    } // Loop over files.

    legend_graph->SetX1( can_graph->GetLeftMargin() );
    legend_graph->SetX2( 1. - can_graph->GetRightMargin() );
    legend_graph->SetY1( 1. -can_graph->GetTopMargin() - 0.05 - 0.18);
    legend_graph->SetY2( 1. -can_graph->GetTopMargin() - 0.05 );

    legend_graph->Draw();

    TString setup_;
    if( MC_files_[0].Contains("isolated") ) setup_ = "isolated";
    else if( MC_files_[0].Contains("calibrated") ) setup_ = "calibrated";

    //== Add CMS touch to the plot.
//      Finish_canvas_narrow( can_graph );
     
    can_graph->SaveAs(TString::Format( folder_ + "/CalibrationFactors_graph_" + plots_label + "_%iGeV.C", static_cast<int> (Ethresh_ ) ) );
    can_graph->SaveAs(TString::Format( folder_ + "/CalibrationFactors_graph_" + plots_label + "_%iGeV.pdf", static_cast<int> (Ethresh_ ) ) );
    
    calibration_parameters_.clear();
  } // Loop over sectors.
}













void Unfolder::CalibrationFactors_oneCanvas(bool draw_functions){
  cout << "Unfolder::CalibrationFactors_oneCanvas" << endl;

  // We do not wish to draw fits in intermediate functions.
  SetFitDraw( false );

  // Start by executing the calibration determination to obtain the necessary parameters.

  TCanvas *can;

  for(int file = 0; file < MC_files_.size(); file++){

      // Create a canvas to plot.

      PrepareCanvas( can, TString::Format("CalibrationFactors_oneCanvas_" + printLabel_[ MC_files_[file] ] ) );
      cout << "Prepared canvas" << endl;

      // Two legends for legibility.
      double common_legend_border = ((1. - can->GetRightMargin()) + 0.65)/2.;
      double top_legend_border = 1. - can->GetTopMargin();

      TLegend *legend_det = new TLegend(
	0.65, 
	0.58, 
	common_legend_border, 
	top_legend_border);
        legend_det->SetFillStyle( 0 );
        legend_det->SetBorderSize( 0 );	

      TLegend *legend_det_2 = new TLegend(
	common_legend_border, 
	0.58, 
	1. - can->GetRightMargin(), 
	top_legend_border);
        legend_det_2->SetFillStyle( 0 );
        legend_det_2->SetBorderSize( 0 );	


      TString drawoptions_gre = "ape";
      TString drawoptions = "hist";

    for( int sector = 1; sector <= 16; sector++){
      cout << "Unfolder::CalibrationFactors_oneCanvas\tsector\t" << sector << endl;
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
  
      TF1* calibration_function = new TF1(TString::Format("calibration_function_sector_%i_" + printLabel_[ MC_files_[file] ], sector), "[0] + [1] * log( [2] + x)", Ethresh_, 900.);
      calibration_function->SetParameters( alpha, beta, gamma );
      calibration_function->SetLineWidth( 1 );
      TH1* calibration_histogram_orig = calibration_function->GetHistogram();
      TH1D* calibration_histogram = (TH1D*)calibration_histogram_orig->Clone(TString::Format("Calibration_values_sector_%i", sector) );
      
//      legend_det->AddEntry( calibration_histogram, TString::Format( legend_info_[ MC_files_[file]] + " (sector %i)", Correct_Sector(sector) ), "l");
      if( sector%2 == 1. ){ legend_det->AddEntry( calibration_histogram, TString::Format( "Sec. %i", Correct_Sector(sector) ), "l"); }
      if( sector%2 == 0. ){ legend_det_2->AddEntry( calibration_histogram, TString::Format( "Sec. %i", Correct_Sector(sector) ), "l"); }

      can->cd();

      calibration_histogram->SetLineColor( getColor( sector ) );
      calibration_histogram->GetXaxis()->SetTitle("E_{det} [GeV]");
      calibration_histogram->GetYaxis()->SetRangeUser(0., 5.);
      calibration_histogram->GetYaxis()->SetTitle("1/<E_{det}/E_{gen}>");

      if( draw_functions ){
        Prepare_1Dplot( calibration_histogram ) ;
        calibration_histogram->GetYaxis()->SetRangeUser(0., 5.);
        calibration_histogram->DrawCopy(drawoptions);
        drawoptions = "histsame";
        drawoptions_gre = "pesame";
      }

      gre_meas->SetMarkerColor( getColor(sector) );
      gre_meas->SetLineColor( getColor(sector) );
      gre_meas->GetXaxis()->SetTitle("E_{det} [GeV]");
      gre_meas->GetYaxis()->SetTitle("1/<E_{det}/E_{gen}>");

      Prepare_1Dplot( gre_meas ) ;
      gre_meas->GetYaxis()->SetRangeUser(0., 5.);
      gre_meas->Draw( drawoptions_gre );
      drawoptions_gre = "pesame";

    } // Loop over sectors.

    legend_det->Draw();
    legend_det_2->Draw();

    TString with_fit = "_fit";
    if( !draw_functions ){ with_fit = "_no_fit"; }

    TString setup;
    if( MC_files_[file].Contains("isolated") ) setup = "isolated";
    else if( MC_files_[file].Contains("calibrated") ) setup = "calibrated";

    TLatex Tl; Tl.SetTextFont(43); Tl.SetTextSize(35);
    Tl.DrawText( 180.874,4., legend_info_[ MC_files_[file]] );


    TString savenameC = TString::Format( folder_ + "CalibrationFactors_oneCanvas_" + printLabel_[ MC_files_[file] ] + with_fit + "_" + setup + ".C") ;
    TString savenamepdf = TString::Format( folder_ + "CalibrationFactors_oneCanvas_" + printLabel_[ MC_files_[file] ] + with_fit + "_" + setup + ".pdf") ;

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
	   can_slice->SaveAs(TString::Format( folder_ + "/Slice_%i.C", bin_E) );
           can_slice->SaveAs(TString::Format( folder_ + "/Slice_%i.pdf", bin_E) );
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
	 	 
     hEdet_axis	->GetYaxis()->SetTitle("1/< #frac{E_{det}}{E_{gen}} >");  
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
	 	 
     Response_meas	->GetYaxis()->SetTitle("1/< #frac{E_{det}}{E_{gen}} >");  
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


void Unfolder::PlotResponseMatrix(int file_, TString setup ){

  cout << "\tUnfolder::PlotResponseMatrix(int " << file_ << ", TString " << setup << ")"  << endl;

  TCanvas *can;
  PrepareCanvas_2D(can, "Response_matrix_" + printLabel_[ MC_files_[file_] ] );

  cout << "\tUnfolder::PlotResponseMatrix\tPrepped canvas" << endl;

  TH2D* hRes;
  Get_ResponseMatrix( file_, hRes );

  cout << "\tUnfolder::PlotResponseMatrix\tGot matrix" << endl;

  Prepare_2Dplot( hRes );

  cout << "\tUnfolder::PlotResponseMatrix\tPrepped matrix" << endl;

  hRes->GetXaxis()->SetTitle("E_{det} [GeV]");
  hRes->GetYaxis()->SetTitle("E_{gen} [GeV]");

  hRes->Draw("colz");
  
  TLatex Tl; 
  Tl.SetTextFont(43); 
  Tl.SetTextSize(35);
  Tl.SetTextAlign( 13 );
  Tl.DrawLatex( 
	can->GetLeftMargin(),
	1. - can->GetTopMargin(), legend_info_gen_[ MC_files_[file_]] );

  //== Add MC details to plot.
  if( MC_files_[file_].Contains("ated") ){  Finish_canvas_2D( can, "narrow", "rightish" ); }
  else if( MC_files_[file_].Contains("raw_wide") ) { Finish_canvas_2D( can, "raw_wide", "rightish" ); }
  else if( MC_files_[file_].Contains("raw") ) { Finish_canvas_2D( can, "raw", "rightish" ); }

  TString MCtype = "Pythia6 (Z2*)";
  TString showertype = "Shower Library";
  if( set_of_tags_["legends"][MC_files_[file_]].Contains("ull") ){ showertype = "Full Simulation"; }

  TLatex *   tex;
  tex = new TLatex(0.5622857,0.25, MCtype);
  tex->SetNDC();
  tex->SetTextAlign(11);
  tex->SetTextFont(43);
  tex->SetTextSize(35);
  tex->SetLineWidth(2);
  tex->Draw();

  tex = new TLatex(0.5622857,0.20,showertype);
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->SetTextFont(43);
  tex->SetTextSize(35);
  tex->SetLineWidth(2);
  tex->Draw();

  //== added details to plot. Save.

  can->SetLogz();

  can->SaveAs( TString::Format(folder_ + "Response_matrix_"  + "_" + printLabel_[ MC_files_[file_] ] + "_" + setup + ".C") );
  can->SaveAs( TString::Format(folder_ + "Response_matrix_"  + "_" + printLabel_[ MC_files_[file_] ] + "_" + setup +  ".pdf") );

  //== Repeat with the response matrix as stored in RooUnfold if possible.
  if( setup == "unfold"){

    cout << "\tUnfolder::PlotResponseMatrix\tUNFOLD" << endl;

    TFile* _file0 = TFile::Open( MC_files_[file_], "Read");
    RooUnfoldResponse* response = (RooUnfoldResponse*)_file0->Get("response");
    hRes = (TH2D*)response->Hresponse();
    Prepare_2Dplot( hRes );

    hRes->GetXaxis()->SetTitle("E_{det} [GeV]");
    hRes->GetYaxis()->SetTitle("E_{gen} [GeV]");
//  SetD2nDxDy( hRes );

    hRes->Draw("colz");
  
    TLatex Tl; 
    Tl.SetTextFont(43); 
    Tl.SetTextSize(35);
    Tl.SetTextAlign( 13 );
    Tl.DrawLatex( 
	can->GetLeftMargin(),
	1. - can->GetTopMargin(), legend_info_gen_[ MC_files_[file_]] );

//    Finish_canvas( can, "rightish");

   TLatex *   tex = new TLatex(0.12,0.92,"Pythia6 (Z2*)");
   tex->SetTextAlign(13);
   tex->SetTextFont(43);
   tex->SetTextSize(35);
   tex->SetLineWidth(2);
//   tex->Draw();

   tex = new TLatex(0.80,0.935,"7 TeV");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.05454545);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.5622857,0.8500625,"CMS");
tex->SetNDC();
   tex->SetTextFont(61);
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.5622857,0.80,"Simulation");
tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetTextFont(52);
   tex->SetTextSize(0.0456);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.12,0.935,"anti-k_{t} (R=0.5) (-6.6<#eta<-5.2)");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->Draw();
 
    can->SetLogz();

    can->SaveAs( TString::Format(folder_ + "RooUnfold_Response_matrix_"  + "_" + printLabel_[ MC_files_[file_] ] + "_" + setup + ".C") );
    can->SaveAs( TString::Format(folder_ + "RooUnfold_Response_matrix_"  + "_" + printLabel_[ MC_files_[file_] ] + "_" + setup +  ".pdf") );
  }
}







void Unfolder::PlotResponseMatrix(TString file_, TString setup ){

  cout << "\n===\tPlotResponseMatrix(TString " << file_ << ", TString " << setup << ")" << endl;

  TCanvas *can;
  PrepareCanvas_2D(can, "Response_matrix_" + file_ );

  TH2D* hRes;
  cout << "\n===\tPlotResponseMatrix" << endl;
  Get_ResponseMatrix( file_, hRes );
  cout << "\n===\tPlotResponseMatrix" << endl;
  Prepare_2Dplot( hRes );
  cout << "\n===\tPlotResponseMatrix" << endl;
//  SetD2nDxDy( hRes );
  cout << "\n===\tPlotResponseMatrix" << endl;

  hRes->Draw("colz");



  can->SetLogz();

  can->SaveAs( TString::Format(folder_ + "Response_matrix_" + setup +  ".C") );
  can->SaveAs( TString::Format(folder_ + "Response_matrix_" + setup +  ".pdf") );

  cout << "€€€\t" << hRes->Integral() << "\tevents" << endl;
}






void Unfolder::Plot_NjetsMatrix(int file_ ){

  TCanvas *can;
  PrepareCanvas_2D(can, "NumberOfJets_" + printLabel_[ MC_files_[file_] ] );

  TH2D* hRes;
  Get_Njets( file_, hRes );
  Prepare_2Dplot( hRes );
  SetD2nDxDy( hRes );

  hRes->Draw("colz");

  can->SetLogz();

  can->SaveAs( TString::Format(folder_ + "NumberOfJets_" + label_ + "_" + printLabel_[ MC_files_[file_] ] + ".C") );
  can->SaveAs( TString::Format(folder_ + "NumberOfJets_" + label_ + "_" + printLabel_[ MC_files_[file_] ] + ".pdf") );

  cout << "€€€\t" << hRes->Integral() << "\tevents" << endl;
}


void Unfolder::Plot_EgenEtaMatrix(int file_ ){

  TCanvas *can;
  PrepareCanvas_2D(can, "Egen_vs_eta_matrix_" + printLabel_[ MC_files_[file_] ] );

  TH2D* hRes;
  Get_EgenEtaMatrix( file_, hRes );
  Prepare_2Dplot( hRes );
  SetD2nDxDy( hRes );

  hRes->Draw("colz");

  can->SetLogz();

  can->SaveAs( TString::Format(folder_ + "Egen_vs_eta_matrix_" + label_ + "_" + printLabel_[ MC_files_[file_] ] + ".C") );
  can->SaveAs( TString::Format(folder_ + "Egen_vs_eta_matrix_" + label_ + "_" + printLabel_[ MC_files_[file_] ] + ".pdf") );

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





int Unfolder::Extract_2D_energy_distributions(int file_, TH2D* &hGenE_selection, TH2D* &hDetE_selection, TH2D* &hResponse_selection, TH2D* &hInverseResponse_selection, TH2D* &hDetE, TH2D* &hGenE, TH1D* &hGen, TString setup){

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
      else if( setup == "bad_sectors" ){
	first_sector = 13;
	last_sector = 14;
      }
      else if( setup.IsDigit() ){
	if( setup.Atoi() >= 1 && setup.Atoi() <= 16 ){
	  first_sector = setup.Atoi();
	  last_sector = setup.Atoi();
	}
	else return 0;

      }

       // -- Extract the response, Edet and Egen distributions from the files.
  
      THnSparseD* hResponse 		= (THnSparseD*)_file0->Get("hResponse_gen_phi");			//hResponse->Draw("colz");
      THnSparse* hGenE_fine_phi   	= (THnSparse*)_file0->Get("hGen_fine_phi");
      THnSparse* hDetE_fine_phi 	= (THnSparse*)_file0->Get("hDet_fine_phi");
cout << "(0) \t Yes" << endl;
      THnSparseD* hInverseResponse 	= (THnSparseD*)_file0->Get("hInverseResponse_gen_phi");

cout << "(1) \t Yes" << endl;

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

      if( hInverseResponse ){
        hInverseResponse_selection = (TH2D*)hInverseResponse->Projection(1,0);
          hInverseResponse_selection->RebinX( rebinner );
          hInverseResponse_selection->RebinY( rebinner );
          hInverseResponse_selection->Reset();
      }

cout << "(2) \t Yes" << endl;

      for(int bin_phi = first_sector; bin_phi <= last_sector; bin_phi++){
	if( setup == "good_sectors" ){
	  if( binary_search(bad_sectors.begin(), bad_sectors.end(), bin_phi ) ){
	    continue;
	  }
	}

        THnSparse* hGenE_sparse = 		(THnSparse*)hGenE_fine_phi	->Clone("");	//TString::Format("GenE_bin_%i", bin_phi) );
        THnSparse* hDetE_sparse =	 	(THnSparse*)hDetE_fine_phi	->Clone("");	//TString::Format("DetE_bin_%i", bin_phi) );
        THnSparse* hResponse_sparse =		(THnSparse*)hResponse		->Clone("");	//TString::Format("Response_bin_%i", bin_phi) );
        THnSparse* hInverseResponse_sparse;
	if(hInverseResponse ){ hInverseResponse_sparse = (THnSparse*)hInverseResponse	->Clone("");	} 

        hGenE_sparse 	->GetAxis(2)->SetRange(bin_phi, bin_phi);
        hDetE_sparse 	->GetAxis(2)->SetRange(bin_phi, bin_phi);
        hResponse_sparse ->GetAxis(2)->SetRange(bin_phi, bin_phi);
        if(hInverseResponse ){ hInverseResponse_sparse ->GetAxis(2)->SetRange(bin_phi, bin_phi); }

	//== Extract the distributions from their Sparse distributions and rebin if necessary.
        TH2D* hGenE_selection_phi  = (TH2D*)hGenE_sparse->Projection(1,0);
          hGenE_selection_phi	->RebinX( rebinner );
	  hGenE_selection_phi	->RebinY( rebinner );
        TH2D* hDetE_selection_phi  = (TH2D*)hDetE_sparse->Projection(1,0);
          hDetE_selection_phi	->RebinX( rebinner );
	  hDetE_selection_phi	->RebinY( rebinner );
        TH2D* hResponse_selection_phi = (TH2D*)hResponse_sparse->Projection(1,0);
	  hResponse_selection_phi->RebinX( rebinner );
	  hResponse_selection_phi->RebinY( rebinner );
        TH2D* hInverseResponse_selection_phi;
	  if( hInverseResponse ){ hInverseResponse_selection_phi = (TH2D*)hInverseResponse_sparse->Projection(1,0);
	  hInverseResponse_selection_phi->RebinX( rebinner );
	  hInverseResponse_selection_phi->RebinY( rebinner ); }

	//== Add the current sector to the total.
        hGenE_selection			->Add( hGenE_selection_phi );
        hDetE_selection			->Add( hDetE_selection_phi );
        hResponse_selection		->Add( hResponse_selection_phi );
        if( hInverseResponse ){ hInverseResponse_selection	->Add( hInverseResponse_selection_phi ); }

        // Delete the histograms to avoid memory leak.
        hGenE_selection_phi->~TH2();
        hDetE_selection_phi->~TH2();
        hResponse_selection_phi->~TH2();
        if( hInverseResponse ){ hInverseResponse_selection_phi->~TH2();}

        hGenE_sparse->~THnSparse();
        hDetE_sparse->~THnSparse();
        hResponse_sparse->~THnSparse();
        if( hInverseResponse ){ hInverseResponse_sparse->~THnSparse(); }
      } // Loop over phi.
  cout << "Done extractin'" << endl;

  return 0;
}


void Unfolder::FillAnew_1D(TH1D* hOld, TH1D* &hNew, TRandom3* &rand){
   

  for(int bin = 0; bin <= hOld->GetNbinsX()+1; bin++){
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
   

  for(int bin = 0; bin <= hOld->GetNbinsX()+1; bin++){
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

  for(int binx = 0; binx <= hOld->GetNbinsX()+1; binx++){
    for(int biny = 0; biny <= hOld->GetNbinsY()+1; biny++){ 
      // Measured distribution.
      int nevents = hOld->GetBinContent( binx, biny );
      int nevents_new;

      if( nevents < 10000 ){ 
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
  cout << "\t===\tUnfolder::SetCastorJetEnergy_norm\t===\t" << renorm << endl;
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

  cout << "\n\t\t===Unfolder--Calculate_smearedBackError_covariance\titerations\t" << iterations << endl;
  int nPoisson = 5000;
  int file_ = -1;
  int MC_ = 0;

  cout << "===Renorm===\t"	<< renorm_ << endl;
  cout << "===Response===\t" 	<< response->GetTitle() << endl;
  cout << "===hData===\t" 	<< hData->Integral() << endl;
  cout << "===hUnfold===\t" 	<< hUnfold->Integral() << endl;

  TCanvas* can_evolution_of_data;
  PrepareCanvas( can_evolution_of_data, "can_evolution_of_data");

  hData->Draw("hist");
  hUnfold->SetLineColor( kGreen - 3);
  hUnfold->Draw("histsame");

  ofstream testing_covariance;
  testing_covariance.open("Testing_covariance.txt", ios::app);

  cout << "***\t" << iterations << "\t" << hUnfold->Integral() << "\t" << Ethresh_ << endl;

  // -- Open files.
  TFile* _file_det, *_file_unfold;
  if( file_ < MC_files_.size() ){	_file_det = new TFile( MC_files_[file_], "Read");}
  else if( file_ == -1 ){		_file_det = new TFile( datafile_, "Read");	}

  _file_unfold = new TFile( MC_files_[MC_], "Read"); 


  // (1) Extract the histograms.
  // cout << "Unfolder--Calculate_smearedBackError (1)" << endl;

  TH1D* hMiss 		= (TH1D*)_file_unfold->Get("hCastorJet_miss_all");
  //TH1D* hFake 		= (TH1D*)_file_unfold->Get("hCastorJet_fake_all");	
  //TH1D* hFake_bis 	= (TH1D*)response->Hfakes();	
  TH1D* hFake		= (TH1D*)response->Hfakes();
  TH2D* hResponse_proj 	= (TH2D*)response->Hresponse(); 
  TH1D* hTruth 		= (TH1D*)response->Htruth();
  TH1D* hMeasured	= (TH1D*)response->Hmeasured();	

  renorm_ = hData->Integral()/hMeasured->Integral(); 

  cout << "===hFake===\t" << hFake->Integral() << endl;

  cout << "\n\t*Truth\t"	<< hTruth->GetNbinsX() 		<< "\t" << hTruth->GetBinLowEdge( 1 ) 
				<< "\t" << hTruth->GetXaxis()->GetBinUpEdge( hTruth->GetNbinsX() ) 		<< "\t" << hTruth->Integral() << endl;
  cout << "\t*Measured\t" 	<< hMeasured->GetNbinsX() 	<< "\t" << hMeasured->GetBinLowEdge( 1 ) 
				<< "\t" << hMeasured->GetXaxis()->GetBinUpEdge( hMeasured->GetNbinsX() ) 	<< "\t" << hMeasured->Integral() << endl;
  cout << "\t*Fake\t" 		<< hFake->GetNbinsX() 		<< "\t" << hFake->GetBinLowEdge( 1 ) 
				<< "\t" << hFake->GetXaxis()->GetBinUpEdge( hFake->GetNbinsX() ) 		<< "\t" << hFake->Integral() << endl;
  cout << "\t*Miss\t" 		<< hMiss->GetNbinsX() 		<< "\t" << hMiss->GetBinLowEdge( 1 ) 
				<< "\t" << hMiss->GetXaxis()->GetBinUpEdge( hMiss->GetNbinsX() ) 		<< "\t" << hMiss->Integral() << endl;

  //-- Cut off the unneeded parts of vectors and matrices.
  TH1D* hLengthVectorsEmin = (TH1D*)hData->Clone("LengthEmin");
 
  TCanvas *can_resp; 
  PrepareCanvas_2D(can_resp, "Response");
  can_resp->SetLogz();

  TH2D* hResponse = (TH2D*)hResponse_proj->Clone("Response_2D");
  hResponse->SetName("Response_2D");
  hResponse->Draw("colz");


  TCanvas *canData = new TCanvas("canData", "canData", 1.);
    TPad* pad_abs_, *pad_ratio_;
    PrepareCanvas(canData, "Pad_absolute");
    SplitCanvas(canData, pad_abs_, pad_ratio_);  
  pad_abs_->cd();
  hData->GetXaxis()->SetNdivisions(504);
  hData->GetXaxis()->SetTitle("E [GeV]");
  hData->GetYaxis()->SetTitle("dN/dE");

  TH1D* hData_subhisto;
  GetSubHistogram( hData, hData_subhisto, Eplotmin_, 2100.);
  SetDnDx( hData_subhisto );
  hData_subhisto->Scale( 1./ (1.2* 1E5) );
  hData_subhisto->GetXaxis()->SetNdivisions(504);
  hData_subhisto->GetXaxis()->SetTitle("E [GeV]");
  hData_subhisto->GetYaxis()->SetTitle("d#sigma/dE [mb/GeV]");
  hData_subhisto->Draw("ehist");
  pad_abs_->SetLogy();

    pad_ratio_->cd();
  TH1D* hData_copy = (TH1D*)hData_subhisto->Clone("copy");
  hData_copy->GetYaxis()->SetTitle("ratio");
  hData_copy->Divide( hData_subhisto );
  hData_copy->GetYaxis()->SetRangeUser(0.75, 1.25);
  hData_copy->GetYaxis()->SetTitleSize( hData_copy->GetYaxis()->GetTitleSize()*1.5 );
  hData_copy->GetYaxis()->SetLabelSize( hData_copy->GetYaxis()->GetLabelSize()*1.5 );
  hData_copy->GetYaxis()->SetTitleOffset( hData_copy->GetYaxis()->GetTitleOffset()* 0.7 );
  hData_copy->GetXaxis()->SetTitleSize( hData_copy->GetXaxis()->GetTitleSize()*1.5 );
  hData_copy->GetXaxis()->SetLabelSize( hData_copy->GetXaxis()->GetLabelSize()*1.5 );
  hData_copy->Draw("ehist");

  TLegend *leg_data = new TLegend(0.6,.6, 1- canData->GetRightMargin(), 1 - pad_abs_->GetTopMargin() );
  leg_data->SetFillColor( kWhite );
  leg_data->AddEntry( hData_copy, "Measured data", "l" );

  //==

  TCanvas *canMiss = new TCanvas("canMiss", "canMiss", 1.);
  hMiss->GetXaxis()->SetNdivisions(504);
  hMiss->GetXaxis()->SetTitle("E_{miss} [GeV]");
  hMiss->GetYaxis()->SetTitle("dN/dE");
  hMiss->Draw("ehist");

  TCanvas *canFake = new TCanvas("canFake", "canFake", 1.);
  hFake->GetXaxis()->SetNdivisions(504);
  hFake->GetXaxis()->SetTitle("E_{fake} [GeV]");
  hFake->GetYaxis()->SetTitle("dN/dE");
  hFake->Draw("ehist");

  TCanvas *canRes_meas = new TCanvas("canRes_meas", "canRes_meas", 1.);
    TH1D* hRes_X_old = (TH1D*)hResponse->ProjectionX();
  hRes_X_old->GetXaxis()->SetNdivisions(504);
  hRes_X_old->GetXaxis()->SetTitle("E_{resp. meas} [GeV]");
  hRes_X_old->GetYaxis()->SetTitle("dN/dE");
  hRes_X_old->Draw("ehist");  

  TCanvas *canRes_true = new TCanvas("canRes_true", "canRes_true", 1.);
    TH1D* hRes_Y_old = (TH1D*)hResponse->ProjectionY();
  hRes_Y_old->GetXaxis()->SetNdivisions(504);
  hRes_Y_old->GetXaxis()->SetTitle("E_{resp. true} [GeV]");
  hRes_Y_old->GetYaxis()->SetTitle("dN/dE");
  hRes_Y_old->Draw("ehist");


  TCanvas *canMeasured = new TCanvas("canMeasured", "canMeasured", 1.);
  hMeasured->GetXaxis()->SetNdivisions(504);
  hMeasured->GetXaxis()->SetTitle("E_{meas} [GeV]");
  hMeasured->GetYaxis()->SetTitle("dN/dE");
  hMeasured->Draw("ehist");

  TCanvas *canTruth = new TCanvas("canTruth", "canTruth", 1.);
  hTruth->GetXaxis()->SetNdivisions(504);
  hTruth->GetXaxis()->SetTitle("E_{true} [GeV]");
  hTruth->GetYaxis()->SetTitle("dN/dE");
  hTruth->Draw("ehist");
  
  // (2) Create empty copies.
  // cout << "Unfolder--Calculate_smearedBackError (2)" << endl;

  TH1D* hSmear; //		= (TH1D*)hData->Clone("hSmear");
  TH1D* hData_new 	= (TH1D*)hData->Clone("hData_new");			
  TH1D* hMiss_new 	= (TH1D*)hMiss->Clone("hMiss_new");			
  TH1D* hFake_new 	= (TH1D*)hFake->Clone("hFake_new");
  TH1D* hTruth_new 	= (TH1D*)hTruth->Clone("hTruth_new");		
  TH1D* hMeasured_new 	= (TH1D*)hMeasured->Clone("hMeasured_new");		
  TH2D* hResponse_new 	= (TH2D*)hResponse->Clone("hResponse_new");	

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

  cout << "Unfolder--Calculate_smearedBackError_covariance - Clone spread" << endl;
 
  pad_abs_->cd();
  hSmear = (TH1D*)response->ApplyToTruth( hUnfold );
  hSmear->Add( hFake, renorm_ );
  hSmear->SetLineStyle( 7 );
  hSmear->SetLineColor( kBlue );
  TH1D* hSmear_first_copy = (TH1D*) hSmear->Clone("Smearing");  
  
  TRandom3* rand = new TRandom3();
    rand->SetSeed( iterations  );

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

    // (3) Fill with Poissonian distribution.
    //-- MC distributions.
    FillAnew_2D( hResponse, hResponse_new, rand);

    //--------------------------------------------------------------------------------------//
    // OPTION - Vary the fake distribution, calculate the measured and truth distributions. //
    //--------------------------------------------------------------------------------------//

    FillAnew_1D( hFake, hFake_new, rand);
    FillAnew_1D( hMiss, hMiss_new, rand);

    //-- Determine the number of measured jets as: measured  = fakes + matched(det)
    TH1D* hRes_X = (TH1D*)hResponse_new->ProjectionX();
    hMeasured_new = (TH1D*)hFake_new->Clone("hMeasured_new");
    hMeasured_new->Add( hRes_X, 1.); 
    
    TH1D* hRes_Y = (TH1D*)hResponse_new->ProjectionY();
    hTruth_new = (TH1D*)hMiss_new->Clone("hTruth_new");
    hTruth_new->Add( hRes_Y, 1.); 

    //-----------------------//
    // (4) Unfold-and-smear. //
    //-----------------------//

    //-- The code below draws the varied distributions and compares them to the actual fakes/measured/truth distributions.
    if( n_spread < 1000 ){

      canMiss->cd();
      hMiss_new->SetLineColor( kYellow - 3 );
      hMiss_new->SetLineStyle( 2 );
      hMiss_new->DrawClone("histsame");

      canFake->cd();
      hFake_new->SetLineColor( kYellow - 3 );
      hFake_new->SetLineStyle( 2 );
      hFake_new->DrawClone("histsame");

      canMeasured->cd();
      hMeasured_new->SetLineColor( kYellow - 3 );
      hMeasured_new->SetLineStyle( 2 );
      hMeasured_new->DrawClone("histsame");

      canTruth->cd();
      hTruth_new->SetLineColor( kYellow - 3 );
      hTruth_new->SetLineStyle( 2 );
      hTruth_new->DrawClone("histsame");

      canRes_meas->cd();
      hRes_X->SetLineColor( kYellow - 3 );
      hRes_X->SetLineStyle( 2 );
      hRes_X->DrawClone("histsame");

      canRes_true->cd();
      hRes_Y->SetLineColor( kYellow - 3 );
      hRes_Y->SetLineStyle( 2 );
      hRes_Y->DrawClone("histsame");
    }

    //-- Create a new RooUnfold response object.
    RooUnfoldResponse *response_new = new RooUnfoldResponse( hMeasured_new, hTruth_new, hResponse_new );

    //-- Make sure that the unfolded distribution has only real, possible numbers.
    if( !BadNumerals(hUnfold) ){ cout << "This is not good\t" << n_spread << endl; n_spread--; continue; }

    //-----------------//
    //-- S = R * U + F //
    //-----------------//
    renorm_ = hData->Integral()/hMeasured_new->Integral(); 
    hSmear = (TH1D*) response_new->ApplyToTruth( hUnfold ); 
    hSmear->Add( hFake_new, renorm_ );

    //-- Draw the smeared distribution.
    if( n_spread < 100){
      canData->cd();

      pad_abs_->cd();

      TH1D* hSmear_copy;
      GetSubHistogram( hSmear, hSmear_copy, Eplotmin_, 2100.);  
      SetDnDx( hSmear_copy );
      hSmear_copy->Scale( 1./ (1.2* 1E5) );

      hSmear_copy->SetName( TString::Format("hSmear_%i", n_spread) );
      hSmear_copy->SetTitle( TString::Format("hSmear_%i", n_spread) );
      hSmear_copy->SetLineColor( kGreen - 3 );
      hSmear_copy->SetLineStyle( 2 );
      hSmear_copy->DrawCopy("histsame");

      hSmear_copy->Divide(hData_subhisto);
      pad_ratio_->cd();
      hSmear_copy->Draw("histsame");

      if( n_spread == 0 ){ leg_data->AddEntry( hSmear_copy, "Backsmeared (var.)", "l" ); }
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
      /*
      if( n_spread <= 10){
	testing_covariance << setprecision(8) <<  smear_value << endl;
      }
      */
    }
  } // Loop over algorithm.

  

  //-----------------------------------------------------//
  //-- (CHECK) Save canvasses with varied distributions. //
  //-----------------------------------------------------//

  canData->SetLogy();
  pad_abs_->cd();


  GetSubHistogram( hSmear_first_copy, hSmear_first_copy, Eplotmin_, 2100.);
  SetDnDx( hSmear_first_copy );
  hSmear_first_copy->SetName("Smearing");
  hSmear_first_copy->SetLineColor( kRed );
  hSmear_first_copy->SetMarkerColor( kRed );
  hSmear_first_copy->SetMarkerStyle( 26 );
  hSmear_first_copy->Scale( 1./ (1.2* 1E5) );  
  hSmear_first_copy->DrawCopy("histsame");

  leg_data->AddEntry( hSmear_first_copy, "Backsmeared data", "l" );
  leg_data->Draw();

  pad_ratio_->cd();
  hSmear_first_copy->Divide( hData_subhisto );  
  hSmear_first_copy->Draw("histsame");

  canData->SaveAs( TString::Format( folder_ + "canData_%i_iterations_deltaPhiMax_0%i_etaband_0%i.C", 	iterations, static_cast<int>( 10. * deltaPhiMax_), static_cast<int>(10. * etawidth_ )  ) );  
  canData->SaveAs( TString::Format( folder_ + "canData_%i_iterations_deltaPhiMax_0%i_etaband_0%i.pdf", 	iterations, static_cast<int>( 10. * deltaPhiMax_), static_cast<int>(10. * etawidth_ ) ) );  

  canFake->SetLogy();
//  hFake->DrawCopy("phistsame");
  canFake->SaveAs( TString::Format( folder_ + "canFake_%i_iterations_deltaPhiMax_0%i_etaband_0%i.C", 	iterations, static_cast<int>( 10. * deltaPhiMax_), static_cast<int>(10. * etawidth_ ) ) );
  canFake->SaveAs( TString::Format( folder_ + "canFake_%i_iterations_deltaPhiMax_0%i_etaband_0%i.pdf", 	iterations, static_cast<int>( 10. * deltaPhiMax_), static_cast<int>(10. * etawidth_ ) ) );

  canTruth->SetLogy();
  hTruth->DrawCopy("phistsame");
  canTruth->SaveAs( TString::Format( folder_ + "canTruth_%i_iterations_deltaPhiMax_0%i_etaband_0%i.C", 	iterations, static_cast<int>( 10. * deltaPhiMax_), static_cast<int>(10. * etawidth_ ) ) );
  canTruth->SaveAs( TString::Format( folder_ + "canTruth_%i_iterations_deltaPhiMax_0%i_etaband_0%i.pdf", 	iterations, static_cast<int>( 10. * deltaPhiMax_), static_cast<int>(10. * etawidth_ ) ) );

  canMiss->SetLogy();
//  hFake->DrawCopy("phistsame");
  canMiss->SaveAs( TString::Format( folder_ + "canMiss_%i_iterations_deltaPhiMax_0%i_etaband_0%i.C", 	iterations, static_cast<int>( 10. * deltaPhiMax_), static_cast<int>(10. * etawidth_ ) ) );
  canMiss->SaveAs( TString::Format( folder_ + "canMiss_%i_iterations_deltaPhiMax_0%i_etaband_0%i.pdf", 	iterations, static_cast<int>( 10. * deltaPhiMax_), static_cast<int>(10. * etawidth_ ) ) );

  canRes_meas->SetLogy();
//  hFake->DrawCopy("phistsame");
  canRes_meas->SaveAs( TString::Format( folder_ + "canRes_meas_%i_iterations_deltaPhiMax_0%i_etaband_0%i.C", 	iterations, static_cast<int>( 10. * deltaPhiMax_), static_cast<int>(10. * etawidth_ ) ) );
  canRes_meas->SaveAs( TString::Format( folder_ + "canRes_meas_%i_iterations_deltaPhiMax_0%i_etaband_0%i.pdf", 	iterations, static_cast<int>( 10. * deltaPhiMax_), static_cast<int>(10. * etawidth_ ) ) );

  canRes_true->SetLogy();
//  hFake->DrawCopy("phistsame");
  canRes_true->SaveAs( TString::Format( folder_ + "canRes_true_%i_iterations_deltaPhiMax_0%i_etaband_0%i.C", 	iterations, static_cast<int>( 10. * deltaPhiMax_), static_cast<int>(10. * etawidth_ ) ) );
  canRes_true->SaveAs( TString::Format( folder_ + "canRes_true_%i_iterations_deltaPhiMax_0%i_etaband_0%i.pdf", 	iterations, static_cast<int>( 10. * deltaPhiMax_), static_cast<int>(10. * etawidth_ ) ) );


  canMeasured->SetLogy();
  hMeasured->SetLineColor( kBlack );
  hMeasured->SetLineStyle( 1 );
  hMeasured->SetLineWidth( 3 );
  hMeasured->SetMarkerStyle( 22 );
  hMeasured->DrawCopy("phistsame");
  canMeasured->SaveAs( TString::Format( folder_ + "canMeasured_%i_iterations_deltaPhiMax_0%i_etaband_0%i.C", 	iterations, static_cast<int>( 10. * deltaPhiMax_), static_cast<int>(10. * etawidth_ ) ) );
  canMeasured->SaveAs( TString::Format( folder_ + "canMeasured_%i_iterations_deltaPhiMax_0%i_etaband_0%i.pdf", 	iterations, static_cast<int>( 10. * deltaPhiMax_), static_cast<int>(10. * etawidth_ ) ) );

  can_evolution_of_data->SetLogy();
  can_evolution_of_data->SaveAs( folder_ + "/Evolution_smeared.pdf");

  //-------------------------------------------------//
  // (5) Extract the average from the distributions. //
  //-------------------------------------------------//

  cout << "Unfolder--Calculate_smearedBackError_covariance - (5): Averages	" << endl;
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
//      testing_covariance << "i\t" << i << "\t" << bin_smear << "\t" << content_smear << "\t" << content_smear * content_smear << endl;
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

  //-- A rather unelegant but temporary-variable-less way to create TH2D.
  TH2D* hCorrelation_denom = new TH2D("hCorrelation_denom", "hCorrelation_denom", 
		hAverage_smear_sq->GetNbinsX(), 
		hAverage_smear_sq->GetXaxis()->GetBinLowEdge( 1 ), 
		hAverage_smear_sq->GetXaxis()->GetBinUpEdge( hAverage_smear_sq->GetNbinsX() ),
		hAverage_smear_sq->GetNbinsX(), 
		hAverage_smear_sq->GetXaxis()->GetBinLowEdge( 1 ), 
		hAverage_smear_sq->GetXaxis()->GetBinUpEdge( hAverage_smear_sq->GetNbinsX() ) );

  for(int row = 1; row <= hAverage_smear_sq->GetNbinsX(); row++){
    double nrow = hAverage_smear_sq->GetBinContent( row );

    for(int col = 1; col <= hAverage_smear_sq->GetNbinsX(); col++){
      double ncol = hAverage_smear_sq->GetBinContent( col );

      hCorrelation_denom->SetBinContent( row, col, sqrt( ncol*nrow ) );
    }
  }

  //-- Poisson-per-poisson iteration 2D histogram.
  TH2D* currentCov = new TH2D("currentCov", "currentCov", hAverage_smear->GetNbinsX(), 0, 26, hAverage_smear->GetNbinsX(), 0, 26);

  //---------------------------------------------------------------------------------------------------------//
  //-- (7) Loop over all iterations and multiply the different (E-mu_E)_DAT x (E-mu_E)_SM with each other. --//
  //---------------------------------------------------------------------------------------------------------//

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
/*
  can_currentCov->SaveAs(folder_ + "Average_minusSmeared.C");
  can_currentCov->SaveAs(folder_ + "Average_minusSmeared.pdf");
*/
  TCanvas *can_e_mu; 	PrepareCanvas(can_e_mu, "E_mu");  
  TString drawoptions_emu = "hist";
  TH1D* hDistr;

//  PrepareCanvas_2D(can_currentCov, "can_2D");

  //-- Loop over variations of response matrix.
  for(int i = 0; i < nPoisson; i++){

   //(CHECK) What does the distribution look like?
   
   if( i < 10 ){
     hDistr = (TH1D*)hAverage_distribution->Clone(TString::Format("hE_mu_%i", iterations) );
     hDistr->Reset();
   }
   

    if( i%500 == 0 ){ cout << "Unfolder--Calculate_smearedBackError_covariance - (7)\t" << i <<  endl; }

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
//        testing_covariance <<"nPoisson\t" << i << "\t(row, col)\t(" << bin_smear_row << " ,\t" << bin_smear_col << ")\t" << covariance_value << endl;
        //-- Bin.
        int bin_3D[3] = { i, bin_smear_row, bin_smear_col };
        covariance_3D.SetBinContent( bin_3D, covariance_value );

//        currentCov->SetBinContent( bin_smear_row, bin_smear_col, covariance_value);
      }
    }

    
  }

  //--------------------------------------------------------------------------------//
  //-- (8) Average out over all histograms by projecting and scaling over nPoisson. //
  //--------------------------------------------------------------------------------//

  //-- Cut off the unneeded parts of vectors and matrices.
  TH1D* hGaugeLengthVectors = (TH1D*)hData->Clone("GaugeLength");


  Ecut_ = 280.;
  GetSubHistogram( hData, hGaugeLengthVectors, Ecut_, 2100.);
  int nbins_newlength = hGaugeLengthVectors->GetNbinsX();

  cout << "Unfolder--Calculate_smearedBackError_covariance - (8)" << endl;
  TH2D* hCovariance_matrix = (TH2D*) covariance_3D.Projection( 2, 1);

  //-- Prepare correlation matrix.
  TH2D* hCorrelation_matrix = (TH2D*)hCovariance_matrix->Clone("hCorrelation_matrix");

  for(int row = 1; row <= hCovariance_matrix->GetNbinsX(); row++){
    for(int col = 1; col <= hCovariance_matrix->GetNbinsY(); col++){
      double x[nPoisson], y[nPoisson];


      double binval = 0.;
      int poisson = 0;
      for(; poisson < nPoisson; poisson++){
	int bin_sparse_smear_col[3] = { poisson, row, col };        
        double sparse_val = covariance_3D.GetBinContent( bin_sparse_smear_col );
        binval += sparse_val;

	//if( binval == sparse_val && poisson > 0 ){ cout << "\t\t\t---\t" << binval << "\t" << poisson << "\t" << row << "\t" << col << endl; }

//	testing_covariance <<binval/( static_cast<double>(poisson) + 1. ) << "\t";
	x[poisson] = poisson+1;
	y[poisson] = binval/( static_cast<double>(poisson) + 1.);
      }

//      testing_covariance  << "\t(row, col)\t(" << row << " ,\t" << col << ")\t" << binval << endl; 

      //-- Add data to the diagonal of the matrix.
      if( row == col ){ binval += hData->GetBinContent( row ); }

      binval = binval/static_cast<double>(nPoisson);

      hCovariance_matrix->SetBinContent( row, col, binval );

      double correlation_denom = hCorrelation_denom->GetBinContent( row, col );
      hCorrelation_matrix->SetBinContent( row, col, binval/correlation_denom);
    }
  }

  TCanvas *can_corr;
  PrepareCanvas_2D( can_corr, TString::Format("canvas_correlation_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i", 
	iterations, static_cast<int>( Ethresh_), static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) ) );
  TMatrixD *mCorrelation = RooUnfoldResponse::H2M( hCorrelation_matrix, hCorrelation_matrix->GetNbinsX(),hCorrelation_matrix->GetNbinsY()  );

  mCorrelation->Draw("colz");
  can_corr->SetLogz();
  can_corr->SaveAs(folder_ + TString::Format("mCor_uncut_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i.C", 
	iterations, static_cast<int>( Ethresh_), static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) ) );
  can_corr->SaveAs(folder_ + TString::Format("mCor_uncut_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i.pdf", 
	iterations, static_cast<int>( Ethresh_), static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) ) );

  mCorrelation->ResizeTo( 
	hCovariance_matrix->GetNbinsX() - nbins_newlength, hCovariance_matrix->GetNbinsX()-1, 
	hCovariance_matrix->GetNbinsY() - nbins_newlength, hCovariance_matrix->GetNbinsY()-1 );

  mCorrelation->Draw("colz");
  can_corr->SetLogz();
  can_corr->SaveAs(folder_ + TString::Format("mCor_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i.C", 
	iterations, static_cast<int>( Ethresh_), static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) ) );
  can_corr->SaveAs(folder_ + TString::Format("mCor_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i.pdf", 
	iterations, static_cast<int>( Ethresh_), static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) ) );

  TH1D* hDiff = (TH1D*)hData->Clone("hDiff");

  renorm_ = hData->Integral()/hMeasured->Integral(); 
  hSmear = (TH1D*) response->ApplyToTruth( hUnfold );
  hSmear->Add( hFake, renorm_ );

  hDiff->Add( hSmear, -1.);	

  TMatrixD *mCovariance = RooUnfoldResponse::H2M( hCovariance_matrix, hCovariance_matrix->GetNbinsX(), hCovariance_matrix->GetNbinsY() );

  //-- Save a copy of the true covariance matrix.
  TCanvas *can_cov = new TCanvas(TString::Format("can_cov_%i", iterations), TString::Format("can_cov_%i", iterations), 1. );
  PrepareCanvas_2D( can_cov, TString::Format("canvas_covariance_%i", iterations) );

  ChangeBinEdges( hCovariance_matrix, hCovariance_matrix, hResponse_proj);
  Prepare_2Dplot( hCovariance_matrix, "E [GeV]", "E [GeV]", 504, 504);

  TH2D* hCovariance_matrix_drawn = (TH2D*)hCovariance_matrix->Clone("hCovariance");
  SetD2nDxDy( hCovariance_matrix_drawn );

  can_cov->SetLogz();

  GetSubHistogram( hCovariance_matrix_drawn, hCovariance_matrix_drawn, Ethresh_, 2100.);
  Prepare_2Dplot( hCovariance_matrix_drawn, "E [GeV]", "E [GeV]", 504, 504);
  Prepare_2Dplot( hCovariance_matrix, "E [GeV]", "E [GeV]", 504, 504);

  hCovariance_matrix->GetXaxis()->LabelsOption("h");
  hCovariance_matrix->GetYaxis()->LabelsOption("v");

  hCovariance_matrix_drawn->Draw("colz");
  hCovariance_matrix_drawn->GetZaxis()->SetRangeUser(
    GetMinimumValue( hCovariance_matrix_drawn ) * 0.9,
    hCovariance_matrix_drawn->GetMaximum()*1.1 );

  DrawWithNegativeLog( hCovariance_matrix_drawn, can_cov );
  Finish_canvas_2D(can_cov);

  can_cov->SaveAs(folder_ + TString::Format("hCov_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i.pdf", 
	iterations, static_cast<int>( Ethresh_), static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) ) );
  can_cov->SaveAs(folder_ + TString::Format("hCov_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i.C", 
	iterations, static_cast<int>( Ethresh_), static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) ) );

  //------------------------------------//
  //-- Measured - smeared difference. --//
  //------------------------------------//

//  GetSubHistogram( hDiff, hDiff, 280., 2100. );
  hDiff->Draw("hist");  
  TVectorD *vDiff = RooUnfoldResponse::H2V( hDiff, hDiff->GetNbinsX());
  vDiff->ResizeTo( hDiff->GetNbinsX() - nbins_newlength , hDiff->GetNbinsX() - 1);
/*
  can_cov->SaveAs(folder_ + TString::Format("/vDiff_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i.C", 
	iterations, static_cast<int>( Ethresh_ ) , static_cast<int>( 10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) ) );
  can_cov->SaveAs(folder_ + TString::Format("/vDiff_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i.pdf", 
	iterations, static_cast<int>( Ethresh_ ) , static_cast<int>( 10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) ) );
*/
  TH1D* hDiff_rel = (TH1D*)hDiff->Clone( "relative_difference");
  for(int bin = 0; bin <= hDiff_rel->GetNbinsX(); bin++){

    hDiff_rel->SetBinContent( bin, 
	hDiff_rel->GetBinContent(bin)/hData->GetBinContent(bin)
	);
  }

  //------------------------//
  //-- Invert the matrix. --//
  //------------------------//
  TMatrixD mInvertCovariance = (*mCovariance);
  mInvertCovariance.SetTol(1.e-23);
  mInvertCovariance.Invert();

  TCanvas *can_inv;
  PrepareCanvas_2D( can_inv, "Canvas_inverseMatrices" );

  mCovariance->Draw("colz");
  can_inv->SaveAs("notinverted_stuff.C");
  mInvertCovariance.Draw("colz");
  can_inv->SaveAs("inverted_stuff.C");


  TH2D* hInvertCovariance = (TH2D*)hCovariance_matrix_drawn->Clone("Inverted_covariance");
  ChangeBinEdges( hInvertCovariance, mInvertCovariance, hCovariance_matrix_drawn);

  hInvertCovariance->GetXaxis()->SetTitle("E [GeV]");
  hInvertCovariance->GetYaxis()->SetTitle("E [GeV]");

  SetD2nDxDy_inverted( hInvertCovariance );
  DrawWithNegativeLog( hInvertCovariance, can_inv );
  Finish_canvas_2D(can_inv);  
  	
  can_inv->SetLogz();
  can_inv->SaveAs(folder_ + TString::Format("Inverted_uncut_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i.C", 
	iterations, static_cast<int>( Ethresh_ ) , static_cast<int>( 10. * deltaPhiMax_ ), static_cast<int>(10. * etawidth_) ) );
  can_inv->SaveAs(folder_ + TString::Format("Inverted_uncut_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i.pdf", 
	iterations, static_cast<int>( Ethresh_ ) , static_cast<int>( 10. * deltaPhiMax_ ), static_cast<int>(10. * etawidth_) ) );

  //-------------------//
  //-- Unity matrix. --//
  //-------------------//  
  TMatrixD mUnity = mInvertCovariance;
  mUnity *= (*mCovariance); //mCovariance_unmodded;
  TCanvas *can_unity = new TCanvas("Unity", "Unity", 1.);
  mUnity.Draw("colz");
  can_unity->SaveAs(folder_ + TString::Format("Unity_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i.C", 
	iterations, static_cast<int>( Ethresh_ ) , static_cast<int>( 10. * deltaPhiMax_ ), static_cast<int>(10. * etawidth_) ) );
  can_unity->SaveAs(folder_ + TString::Format("Unity_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i.pdf", 
	iterations, static_cast<int>( Ethresh_ ) , static_cast<int>( 10. * deltaPhiMax_ ), static_cast<int>(10. * etawidth_) ) );

  mInvertCovariance.ResizeTo(
	hCovariance_matrix->GetNbinsX() - nbins_newlength, hCovariance_matrix->GetNbinsX() -1, 
	hCovariance_matrix->GetNbinsY() - nbins_newlength, hCovariance_matrix->GetNbinsY() -1);

  mInvertCovariance.Draw("colz");
  can_inv->SetLogz();
  can_inv->SaveAs(folder_ + TString::Format("Inverted_C_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i.C", 
	iterations, static_cast<int>( Ethresh_ ) , static_cast<int>( 10. * deltaPhiMax_ ), static_cast<int>(10. * etawidth_) ) );
  can_inv->SaveAs(folder_ + TString::Format("Inverted_C_%i_iterations_%i_GeV_deltaPhiMax_0%i_etaband_0%i.pdf", 
	iterations, static_cast<int>( Ethresh_ ) , static_cast<int>( 10. * deltaPhiMax_ ), static_cast<int>(10. * etawidth_) ) );

  TVectorD vTemp = (*vDiff);
  vTemp *= mInvertCovariance;
  double chi2 = (*vDiff) * vTemp;

  can_cov->SetLogz();

  cout << "\nCHI2/NDF for " << iterations << " iterations is " << chi2 << "/" << vDiff->GetNoElements() << " = " <<  chi2/vDiff->GetNoElements() << "\t" << MC_files_[0] << endl;
  cout << "Newbinslength\t" << nbins_newlength << "\t" << Ecut_ << endl;
  cout << "renorm_\t" << renorm_ << endl;

  for(int idata = 1; idata <= hData->GetNbinsX(); idata++){
    cout << idata << "\tdata\t" << hData->GetBinContent( idata ) << endl;
  }




  return chi2/ ( vDiff->GetNoElements() ) ;
}




void Unfolder::SetSubhistogram_cut(double Ecut){

  Ecut_ = Ecut;
}

void Unfolder::SetSubhistogram_max(double Ecut){

  Emax_ = Ecut;
}


bool Unfolder::BadNumerals( TH1D* hist ){
  bool goodNumerals= true;
  for(int bin = 0; bin <= hist->GetNbinsX(); bin++){
    if( hist->GetBinContent(bin) != hist->GetBinContent(bin) ){ 
      cout << "Shiiit\t" <<  bin << "\t" << hist->GetBinContent(bin) << endl;
      goodNumerals = false; 
      break; 
    }
  }
  return goodNumerals;
}


void Unfolder::SetAddLabel( TString label ){
  addLabel_ = "_" + label + "_";
}











void Unfolder::Smear_gen(int file_){

  TFile* _file_det;
  if( file_ < MC_files_.size() ){
    _file_det = new TFile( MC_files_[file_], "Read");
  
    RooUnfoldResponse *response = (RooUnfoldResponse*)_file_det->Get("response");

    TH1D* hDet = (TH1D*)_file_det->Get("hCastorJet_energy");	//hDet->Scale( renorm_ );

    TH1D* hMiss = (TH1D*)_file_det->Get("hCastorJet_miss_all");	
    TH1D* hFake = (TH1D*)_file_det->Get("hCastorJet_fake_all");	
    TH2D* hResponse 	= (TH2D*)response->Hresponse(); 		
    TH1D* hTruth 	= (TH1D*)response->Htruth();			
    TH1D* hMeasured	= (TH1D*)response->Hmeasured();			

    TCanvas *can;
    PrepareCanvas(can, "Smeared_MC");

    TH1D *hSmeared = (TH1D*)response->ApplyToTruth( hTruth );

//    hSmeared->Add( hFake );
    hSmeared->SetLineColor( kRed );
    hSmeared->SetLineWidth( 3 );
    hSmeared->Draw("hist");

    hDet->SetLineStyle( 3 );
    hDet->SetLineWidth( 3 );
    hDet->Draw("histsame");

    can->SetLogy();

    can->SaveAs(TString::Format("Plots/SmearingBackGen_%i.C", static_cast<int>(Ethresh_)) );
    can->SaveAs(TString::Format("Plots/SmearingBackGen_%i.pdf", static_cast<int>(Ethresh_)) );
  }
}





void Unfolder::Chi2diff_test_data(TString variable, TString file, int method){

  // Prepare variables and objects.
  TH1D* hist_original, *hist_MCdet;
  TH1D* hist_result; 
  TH1D* hist_reference, *hist_previous; 
  TString htitle;
  int iterations_ = 50, iterations_start = 1;
  int actual_iterations = 0;

  double xaxisgraph[ (iterations_ - iterations_start) - 1], chi2diff[(iterations_ - iterations_start) - 1];

  TString norm = "notNorm";

  vector<TH1D*> histos;
  TMatrixD cov_m;

  int current_color = 1;
  
  // Get the detector level distribution.
  if( variable == "lead"){	Get_DetEnergy_lead(-1, hist_original );	htitle = "E_{det} spectrum (leading)";}
  else{				Get_DetEnergy(-1, hist_original );	htitle = "E_{det} spectrum (all)";	}

  TString scalefactor_MC = (set_of_tags_[ "scalefactors" ])[ MC_files_[0] ];
  double scalefactors_MC_ = scalefactor_MC.Atof() ;

  SetCastorJetEnergy_norm( scalefactors_MC_/ scalefactors_Data_ ); 

  TCanvas *can;
  PrepareCanvas( can, "Compare_unfolded_histograms");
  
  double tickx = gStyle->GetTickLength();
  TLegend *legend = new TLegend(0.55, 0.65, 1.-can->GetRightMargin()-tickx, 1.-can->GetTopMargin()-tickx);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0); 

  histos.push_back( hist_reference );

  // Loop over number of Bayesian iterations, and unfold->smear detector level distribution.
  // Calculate chi2 as we go.

  //---------------------------------------------------------------//
  //-- Loop over the files, if there are several to be unfolded. --//
  //---------------------------------------------------------------//


    TCanvas *can_chi2diff;
    PrepareCanvas( can_chi2diff, "CHI2_Diff_" + variable);
    TString drawoptions = "apc";

    vector<double> etaband;
      etaband.push_back(0.0);
      etaband.push_back(0.2);
      etaband.push_back(0.5);

    vector<double> deltaPhiMax;
      deltaPhiMax.push_back(0.2);
      deltaPhiMax.push_back(0.5);

    vector<int> marker_color;

      marker_color.push_back( 4 );

    vector<TString> model;
      model.push_back("");

    vector<TString> matching;
      matching.push_back("_matchE");

    vector<double> Emin;
      Emin.push_back(150.);

    vector<TString> match_symbol;
      match_symbol.push_back("E");

    vector<TString> model_legend;
      model_legend.push_back("p6");
      model_legend.push_back("p8");


    vector<int> markers;
//      markers.push_back( 20 );
//      markers.push_back( 22 );
      markers.push_back( 29 );
      markers.push_back( 21 );

    int file_ = 0;

    for(int _model = 0; _model < model.size(); _model++){
      for(int _eta = 0; _eta < etaband.size(); _eta++){
        for(int _phi = 0; _phi < deltaPhiMax.size(); _phi++){
	  for(int _Emin = 0; _Emin < Emin.size(); _Emin++){
   	    for(int _match = 0; _match < matching.size(); _match++, file_++){

//              TString filename_ = TString::Format( "/user/avanspil/Castor_Analysis/ak5ak5_" + model[_model] + "displaced_unfold_Emin_%f_deltaPhiMax_%f_etaband_%f" + matching[_match] + ".root", Emin[_Emin], deltaPhiMax[_phi], etaband[_eta] ) ;
              TString filename_ = TString::Format( "/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo"+ cutsetup_ + "_unfold_Emin_%f_deltaPhiMax_%f_etaband_%f_all" + matching[_match] + ".root", Emin[_Emin], deltaPhiMax[_phi], etaband[_eta] ) ;

      	      //-- Unfolding for delta chi2.
      	      TFile *_file0 = new TFile( filename_, "read");

      	      RooUnfoldResponse* response;  
      	      if( variable == "all"){ 		response = (RooUnfoldResponse*)_file0->Get("response"); }
      	      else if( variable == "lead"){ 	response = (RooUnfoldResponse*)_file0->Get("response_lead"); }

      	      //--------------------------------------------------//
      	      //-- Loop over the number of Bayesian iterations. --//
      	      //--------------------------------------------------//

      	      for(int iterations = iterations_start; iterations <= iterations_; iterations++){
                // Determine the inrease in Bayesian iterations.

      		RooUnfoldBayes unfold_bayes(response, hist_original, iterations); 
      		unfold_bayes.SetVerbose(0);
      		TH1D* hUnfold = (TH1D*) unfold_bayes.Hreco( RooUnfold::kCovariance );

      		cout << "\tUnfolded\t" << hist_original->Integral() << "\t" << endl;

      		if( iterations == iterations_start ){       
		  hist_previous = (TH1D*)hUnfold->Clone(); 
      		}

      		double chi2diff_ = 0.;

		GetSubHistogram( hist_previous, hist_previous, Eplotmin_, 2100.);	
     		GetSubHistogram( hUnfold, hUnfold, Eplotmin_, 2100.);	

      		chi2diff_ = hUnfold->Chi2Test( hist_previous, "CHI2/NDF"); 


		if( iterations > iterations_start ){
	
		  cout << "Enter\t" << iterations - iterations_start << "\t" << iterations << "\t" << chi2diff_ << endl;

 		  chi2diff[iterations - iterations_start -1] = chi2diff_;
      		  xaxisgraph[iterations - iterations_start -1 ] = iterations;
		}

      		can->cd();

      		cout << "Iterations\t" << iterations << endl;

      		if( iterations == iterations_start ){ hUnfold->Draw("hist"); }
      		else{ hUnfold->SetLineColor( iterations ); hUnfold->Draw("histsame"); }

      		current_color++;
      		actual_iterations++;

      		hist_previous = (TH1D*)hUnfold->Clone();

    	      } // Iterations

			

	      can_chi2diff->cd();

              TGraph* chi2_diff = new TGraph( (iterations_ - iterations_start) - 1, xaxisgraph, chi2diff );
	      chi2_diff->GetHistogram()->GetYaxis()->SetRangeUser(1e-3, 1e3);
	      chi2_diff->GetHistogram()->GetYaxis()->SetTitleOffset(1.25);

	      TString histname = TString::Format( "Chi2NDF_Emin_%i_deltaPhiMax_0%i_etaband_0%i" + matching[_match] + "_" + model[_model], 
		static_cast<int>(Emin[_Emin]), 
		static_cast<int>(10. * deltaPhiMax[_phi]), 
		static_cast<int>(10. * etaband[_eta]));

	      chi2_diff->SetTitle( histname );
	      chi2_diff->SetName( histname );
	      chi2_diff->SetLineWidth( 3 );

	      // Linecolor = model -- linestyle = Emin -- markerstyle = etaband -- markercolor = deltaphi
	      chi2_diff->SetLineColor( getColor( file_+1 ) ); 
	      chi2_diff->SetLineStyle( file_ + 1 );  
              chi2_diff->SetMarkerStyle( file_ + 20 );
              chi2_diff->SetMarkerColor( getColor( file_+1 ) );
	      chi2_diff ->SetMarkerColor( getColor( file_+1 ) );
	      chi2_diff->SetMarkerColor( getColor( file_+1 ) );
	      chi2_diff->SetMarkerColor( getColor( file_+1 ) );
	      chi2_diff->SetMarkerSize( 1.8 );

	      chi2_diff->SetMarkerSize( 1.25 );
	      chi2_diff->GetXaxis()->SetNdivisions( 504 );
	      chi2_diff->GetXaxis()->SetTitle("E [GeV]");

              chi2_diff->GetXaxis()->SetTitle("N_{it.}");
              chi2_diff->GetYaxis()->SetTitle("#Delta#chi^{2}/NDF");
              Prepare_1Dplot( chi2_diff );
              chi2_diff->Draw(drawoptions);
	      drawoptions = "pcsame";

	      legend->AddEntry( chi2_diff, TString::Format( " #Delta#varphi_{max}=0.%i, #eta_{acc}=0.%i", 
		static_cast<int>(10. * deltaPhiMax[_phi] ), 
		static_cast<int>(10. * etaband[_eta] ) ), "lp" );
	      legend->SetFillColor( 0 );
  	    }
	  }
	}
      }
    }

    legend->Draw();
    legend->SetFillColor( 0 );  

    can_chi2diff->SetLogy();

    can_chi2diff->SaveAs( TString::Format(folder_ + "Chi2_diff_data_" + variable + "_%i_iterations_deltaPhiMax_0%i_etaband_0%i_" + norm + ".C", 
 method, iterations_, static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) ) );
    can_chi2diff->SaveAs( TString::Format(folder_ + "Chi2_diff_data_" + variable + "_%i_iterations_deltaPhiMax_0%i_etaband_0%i_" + norm + ".pdf", 
	method, iterations_, static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) ) );
/*
    can->SetLogy();
    can->SaveAs( TString::Format( "Plots/Compare_unfolded_data.C" ) );
    can->SaveAs( TString::Format( "Plots/Compare_unfolded_data.pdf" ) );
*/
    // Finish delta chi2.

}



























//
//-- Model dependence.
//


void Unfolder::Plot_Unfolded_modelDependence(TPad* & pad_, TString variable, int iterations){
  pad_->cd();

  TLegend* legend = new TLegend( 0.60, 0.5, 1. - pad_->GetRightMargin(), 1. - pad_->GetTopMargin()  );
  legend->SetFillStyle(  0);
  legend->SetBorderSize( 0 );


  TH1D* hAverage, *hFirst;
  TString drawoptions = "hist";
  double models = 0.;

  TH1D* hGen, *hDet;

  for(int MC_ = 0; MC_ < MC_files_.size() ; MC_++){
    //-- Generator level.
    if( (set_of_tags_["mc_type"])[MC_files_[MC_]] != "model" && (set_of_tags_["mc_type"])[MC_files_[MC_]] != "actual"){ continue; }
    models++;

    cout << "Working on\t" << MC_files_[MC_] << endl;
    
    //-- Data unfolded.
    TH1D* hData;

    if( variable == "all") { Get_DetEnergy( -1 , hData); }
    if( variable == "lead"){ Get_DetEnergy_lead( -1 , hData); }
    hData->GetYaxis()->SetTitle("#frac{dN}{dE} [1/GeV]");

    Get_DetUnfolded( MC_ , hData, iterations, variable);
    SetDnDx( hData );
    if( MC_ == 0 ){ 
      hAverage = (TH1D*)hData->Clone("Average");
      hDet = (TH1D*)hData->Clone("hModeldep");
    }
    else{
      hAverage->Add( hData );
    }
  }

  hAverage->Scale( 1./models );
  cout << "\t===Average\t" << hAverage->Integral() << endl;

  hFirst = (TH1D*)hAverage->Clone("First");

  hAverage->SetLineWidth(2);
  hAverage->SetMarkerColor( 1 );
  hAverage->SetLineColor( 1 );
  hAverage->GetXaxis()->SetTitle("E [GeV]");
  hAverage->GetYaxis()->SetTitle("#frac{d#sigma}{dE}[mb/GeV]");

  pad_->cd();
  Prepare_1Dplot( hAverage );	// Set title and label size.
  hAverage->Draw("ephist");
  legend->AddEntry( hAverage, "Model average", "ep"  );
  
    drawoptions = "psame";

  for(int bin = 0; bin <= hDet->GetNbinsX(); bin++){ cout << bin << "\t" << hDet->GetBinContent( bin ) << endl; }

  TH1D* hModel_dep = (TH1D*)hDet->Clone("hModelDep");
//  hModel_dep->Draw("histsame");
  cout << "Model_dep\t" << hModel_dep->Integral();
  hModel_dep->Add( hAverage, -1. );

  //
  //-- 
  // 
  for(int bin = 0; bin <= hModel_dep->GetNbinsX(); bin++){
    if( hModel_dep->GetBinContent( bin ) < 0. ){
       hModel_dep->SetBinContent( bin, -1. * hModel_dep->GetBinContent( bin ) );
    }
  }  


  //-- Model dependence.
  for(int MC_ = 0; MC_ < MC_files_.size() ; MC_++){
    //-- Generator level.
    if( (set_of_tags_["mc_type"])[MC_files_[MC_]] != "model" && (set_of_tags_["mc_type"])[MC_files_[MC_]] != "actual"){ continue; }

    //-- Generator level.

    TH1D* hGen, *hDet, *hData;

    /*
    //-- Data unfolded.
    TH1D* hData = (TH1D*)hModel_dep->Clone("model_dep");
    if( MC_ == 0 ){ hData->Scale( -1. ); }
    else{ hData->Scale( 1. ); }
    hData->Add( hAverage );
    */

    Get_DetEnergy( -1 , hData); 

    Get_DetUnfolded( MC_ , hData, iterations, variable);
    SetDnDx( hData );

    cout << "MC\t" << hData->Integral() << "\t" << drawoptions << endl;

    pad_->cd();
    hData->SetLineColor( getColor( MC_+2) );
    hData->SetMarkerColor( getColor( MC_+2 ) );
    hData->SetMarkerStyle( 24 + MC_ );
    hData->SetLineWidth( 2 );
    hData->Draw( drawoptions );

    legend->AddEntry( hData, legend_info_gen_[ MC_files_[MC_] ], "ep"  );

    for(int bin = 0; bin <= hData->GetNbinsX(); bin++){ cout << bin << "\t" << hData->GetBinContent( bin ) << endl; }
  }


//  hAverage->Draw( drawoptions );
//  hAverage->GetYaxis()->SetRangeUser(0., 2.);
  legend->Draw();
  pad_->SetLogy();
  pad_->Update();
}



















void Unfolder::Plot_Unfolded_Ratio_modelDependence(TPad* & pad_, TString variable, int iterations){
  pad_->cd();

  TH1D* hAverage, *hFirst;
  TString drawoptions = "hist";
  double models = 0.;

  TH1D* hGen, *hDet;

  for(int MC_ = 0; MC_ < MC_files_.size() ; MC_++){
    //-- Generator level.
    if( (set_of_tags_["mc_type"])[MC_files_[MC_]] != "model" && (set_of_tags_["mc_type"])[MC_files_[MC_]] != "actual"){ continue; }
    models++;

    cout << "Working on\t" << MC_files_[MC_] << endl;
    
    //-- Data unfolded.
    TH1D* hData;

    if( variable == "all") { Get_DetEnergy( -1 , hData); }
    if( variable == "lead"){ Get_DetEnergy_lead( -1 , hData); }

    Get_DetUnfolded( MC_ , hData, iterations, variable);
    SetDnDx( hData );
    if( MC_ == 0 ){ 
      hAverage = (TH1D*)hData->Clone("Average");
      hDet = (TH1D*)hData->Clone("hModeldep");
    }
    else{
      hAverage->Add( hData );
    }
  }

  hAverage->Scale( 1./models );
  cout << "\t===Average\t" << hAverage->Integral() << endl;

  hFirst = (TH1D*)hAverage->Clone("First");

  hAverage->SetLineWidth(2);
  hAverage->SetMarkerColor( 1 );
  hAverage->SetLineColor( 1 );
  hAverage->GetXaxis()->SetTitle("E [GeV]");

  pad_->cd();
  Prepare_1Dplot( hAverage );	// Set title and label size.
  hAverage->GetXaxis()->SetTitleOffset(
	hAverage->GetXaxis()->GetTitleOffset() * 2 );
  hAverage->Divide( hFirst );

  Prepare_1Dplot( hAverage, pad_ );
  hAverage->GetYaxis()->SetRangeUser(0., 3.);
  hAverage->GetYaxis()->SetNdivisions(504);
  hAverage->Draw("ephist");
  
    drawoptions = "psame";

  for(int bin = 0; bin <= hDet->GetNbinsX(); bin++){ cout << bin << "\t" << hDet->GetBinContent( bin ) << endl; }

  TH1D* hModel_dep = (TH1D*)hDet->Clone("hModelDep");
//  hModel_dep->Draw("histsame");
  cout << "Model_dep\t" << hModel_dep->Integral();
  hModel_dep->Add( hAverage, -1. );

  //
  //-- 
  // 
  for(int bin = 0; bin <= hModel_dep->GetNbinsX(); bin++){
    if( hModel_dep->GetBinContent( bin ) < 0. ){
       hModel_dep->SetBinContent( bin, -1. * hModel_dep->GetBinContent( bin ) );
    }
  }  


  //-- Model dependence.
  for(int MC_ = 0; MC_ < MC_files_.size() ; MC_++){
    //-- Generator level.
    if( (set_of_tags_["mc_type"])[MC_files_[MC_]] != "model" && (set_of_tags_["mc_type"])[MC_files_[MC_]] != "actual"){ continue; }

    //-- Generator level.

    TH1D* hGen, *hDet, *hData;

    /*
    //-- Data unfolded.
    TH1D* hData = (TH1D*)hModel_dep->Clone("model_dep");
    if( MC_ == 0 ){ hData->Scale( -1. ); }
    else{ hData->Scale( 1. ); }
    hData->Add( hAverage );
    */

    Get_DetEnergy( -1 , hData); 

    Get_DetUnfolded( MC_ , hData, iterations, variable);
    SetDnDx( hData );

    cout << "MC\t" << hData->Integral() << "\t" << drawoptions << endl;

    pad_->cd();
    hData->SetLineColor( getColor( MC_+2) );
    hData->SetMarkerColor( getColor( MC_+2 ) );
    hData->SetMarkerStyle( 24 + MC_ );
    hData->SetLineWidth( 2 );
    hData->Divide( hFirst );
    hData->Draw( drawoptions );
    for(int bin = 0; bin <= hData->GetNbinsX(); bin++){ cout << bin << "\t" << hData->GetBinContent( bin ) << endl; }
  }


//  hAverage->Draw( drawoptions );
  hAverage->GetYaxis()->SetRangeUser(0., 3.);
  pad_->Update();
}



//
//-- Systematics uncertainty on position.
//








void Unfolder::Plot_Unfolded_positionDependence(TPad* & pad_, TString variable, int iterations){
  pad_->cd();

  TLegend* legend = new TLegend( 0.60, 0.5, 1. - pad_->GetRightMargin(), 1. - pad_->GetTopMargin()  );

  legend->SetFillStyle(  0);
  legend->SetBorderSize( 0 );


  TH1D *hFirst, *hData;
  TString drawoptions = "ephist";

  //-- Detector level energy extraction.
  if( variable == "all") { Get_DetEnergy( -1 , hData); }
  if( variable == "lead"){ Get_DetEnergy_lead( -1 , hData); }   

  //-- Unfolding
  Get_DetUnfolded( 0 , hData, iterations, variable);
  cout << "Unfolding position\t" << MC_files_[0] << "\t" << hData->Integral() << endl;

  SetDnDx( hData );

  hData->SetMarkerColor( kBlack );
  hData->SetLineColor( kBlack );
  hData->GetXaxis()->SetTitle("E [GeV]");
  hData->GetYaxis()->SetTitle("#frac{dN}{dE} [1/GeV]");
  pad_->cd();

  Prepare_1Dplot( hData );
  hData->Draw( drawoptions );
  drawoptions = "psame";
  legend->AddEntry( hData, legend_info_gen_[ MC_files_[0] ], "ep"  );  

  //-- Model dependence.
  for(int MC_ = 0; MC_ < MC_files_.size() ; MC_++){
    //-- Generator level.
    if( (set_of_tags_["mc_type"])[MC_files_[MC_]] != "position" ){ continue; }
    //-- Generator level.


    TH1D* hGen, *hDet;

    //-- Data unfolded.
    TH1D* hData;

    Get_DetEnergy( -1 , hData); 
    Get_DetUnfolded( MC_ , hData, iterations, variable);
    SetDnDx( hData );

    hData->SetLineColor( getColor( MC_+2) );
    hData->SetMarkerColor( getColor( MC_+2 ) );
    hData->SetMarkerStyle( 24 + MC_ );
    hData->SetLineWidth( 2 );
   
    pad_->cd();
    hData->       Draw( drawoptions );
    drawoptions = "psame";
    legend->AddEntry( hData, legend_info_gen_[ MC_files_[MC_] ], "ep"  );
  }

  legend->Draw();
  pad_->SetLogy();
  pad_->Update();
}




void Unfolder::PlotStartingDistributions_ratio(TString distribution){
  cout << "Unfolder::PlotStartingDistributions_comparingEmin\t" << distribution << endl;

  TCanvas *can_startingDistributions;
  PrepareCanvas(can_startingDistributions, "can_startingDistributions_" + distribution);
  TLegend *legend = new TLegend( 
	can_startingDistributions->GetLeftMargin(),
	can_startingDistributions->GetBottomMargin(),
	0.5,
	0.45);
  legend->SetFillStyle( 0 );
  legend->SetBorderSize( 0 );

  TString distributionname_unmatched, distributionname_all;
  if (distribution == "fake"){ 		distributionname_unmatched = "hCastorJet_fake_all"; distributionname_all = "hCastorJet_energy"; }
  if (distribution == "miss"){ 		distributionname_unmatched = "hCastorJet_miss_all"; distributionname_all = "hGenJet_energy"; }

  TString drawoptions = "phist";

  int file = 0;

  vector<double> etaband;
      etaband.push_back(0.0);
      etaband.push_back(0.2);
      etaband.push_back(0.5);


  vector<double> deltaPhiMax;
      deltaPhiMax.push_back(0.2);
//      deltaPhiMax.push_back(0.4);
      deltaPhiMax.push_back(0.5);

  vector<TString> matching;
      matching.push_back("_matchE");
//      matching.push_back("_matchPhi");

  vector<TString> match_symbol;
      match_symbol.push_back("E");
      match_symbol.push_back("#varphi");

  vector<TString> model;
      model.push_back("Pythia6Z2star");
//      model.push_back("Pythia84C_");

  vector<TString> model_legend;
      model_legend.push_back("p6");
      model_legend.push_back("p8");

  vector<int> model_events;
      model_events.push_back( 547922 );
      model_events.push_back( 670215 );

  vector<double> Emin;
//      Emin.push_back(0.);
      Emin.push_back(150.);

  vector<int> markers;
      markers.push_back( 20 );
      markers.push_back( 22 );
      markers.push_back( 29 );
      markers.push_back( 21 );

  vector<TH1D*> histlist;
  TH1D* hFirst;
  int setup = 0;

  for(int _model = 0; _model < model.size(); _model++){
    for(int _eta = 0; _eta < etaband.size(); _eta++){
      for(int _phi = 0; _phi < deltaPhiMax.size(); _phi++){
	for(int _Emin = 0; _Emin < Emin.size(); _Emin++){
   	    for(int _match = 0; _match < matching.size(); _match++, file++){

              TString filename_ = TString::Format( "/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_" + model[_model] + "_NewGeo_unfold_Emin_%f_deltaPhiMax_%f_etaband_%f_all" + matching[_match] + ".root", Emin[_Emin], deltaPhiMax[_phi], etaband[_eta] ) ;
	    
	      cout << "\n=*=\t" << filename_ << endl;
	      TFile* _file= TFile::Open( filename_, "read" );

	      TString setup_label = TString::Format( model[_model] + "displaced_unfold_Emin_%i_deltaPhiMax_0%i_etaband_0%i" + matching[_match] , 
		static_cast<int>(Emin[_Emin]), 
		static_cast<int>(10. * deltaPhiMax[_phi]), 
		static_cast<int>(10. * etaband[_eta]));
	      PlotResponseMatrix( filename_, setup_label );

	      TString histname = TString::Format( distributionname_unmatched + "_Emin_%i_deltaPhiMax_0%i_etaband_0%i" + matching[_match] + "_" + model[_model], 
		static_cast<int>(Emin[_Emin]), 
		static_cast<int>(10. * deltaPhiMax[_phi]), 
		static_cast<int>(10. * etaband[_eta]));

    	      //== Extracting the distribution.
  	      TH1D* hUnmatched = (TH1D*)_file->Get( distributionname_unmatched );
  	      TH1D* hAll = (TH1D*)_file->Get( distributionname_all );

	      //== Scaling.
              int total_events_nocuts;
	      TTree *tree_numbers = (TTree*)_file->Get("useful_numbers");
	      TString branch_str = (set_of_tags_[ "scalefactors" ])[ MC_files_[0] ];

	      hUnmatched->SetTitle( histname );
	      hUnmatched->SetName( histname );
	      hUnmatched->SetLineWidth( 3 );

	      // Linecolor = model -- linestyle = match -- markerstyle = etaband -- markercolor = deltaphi
	      hUnmatched->SetLineColor( getColor( _model+1 ) ); 

	      hUnmatched->SetLineStyle( file + 1 ); 
	      hUnmatched->SetLineColor( getColor(file + 1) ); 
 
              hUnmatched->SetMarkerStyle( file + 20 );

              hUnmatched->SetMarkerColor( getColor( file+1 ) );

	      hUnmatched->SetMarkerSize( 1.8 );
	      hUnmatched->GetXaxis()->SetNdivisions( 504 );
	      hUnmatched->GetXaxis()->SetTitle("E_{" + distribution + "}[GeV]");
	      hUnmatched->GetYaxis()->SetTitle("Ratio unmatched");

	      can_startingDistributions->cd();

	      hUnmatched->Divide( hAll );
	      hUnmatched->GetYaxis()->SetRangeUser(0.9*1e-3, 1.);
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
	      hFirst->GetYaxis()->SetRangeUser(1e-5, 1.);

	      //legend->AddEntry( hUnmatched, TString::Format( distribution + "#Delta #varphi_{max} = 0.%i, E_{min} = %i GeV, #eta = 0.%i" + matching[_match],
	      legend->AddEntry( hUnmatched, TString::Format( " #Delta#varphi_{max}=0.%i, #eta_{acc}=0.%i", 
		static_cast<int>(10. * deltaPhiMax[_phi] ), 
		static_cast<int>(10. * etaband[_eta] ) ), "lp" );
	      legend->SetFillStyle( 0 );


	    } // Match.
	  } // Emin.
        } // Deltaphi.
      } // Eta.
    } // Model.
   
  legend->Draw();

  can_startingDistributions->SetLogy();
  can_startingDistributions->SaveAs( TString::Format(folder_ + "can_ratio_" + distribution + "_Emin_%i.C", static_cast<int>( Ethresh_ ) ) );
  can_startingDistributions->SaveAs( TString::Format(folder_ + "can_ratio_" + distribution + "_Emin_%i.pdf", static_cast<int>( Ethresh_) ) );

}





void Unfolder::Plot_Unfolded_JESDependence(TPad* & pad_, TString variable, int iterations){
  pad_->cd();

  TLegend* legend = new TLegend( 0.60, 0.7, 0.95, 0.9  );

  TH1D *hFirst, *hData, *hJESup, *hJESdown;
  TString drawoptions = "ephist";

  //-- Detector level energy extraction.
  if( variable == "all") { Get_DetEnergy( -1 , hData); }
  Get_DetEnergy_JESup( -1 , hJESup);
  Get_DetEnergy_JESdown( -1 , hJESdown);


  //-- Unfolding
  Get_DetUnfolded( 0 , hData, iterations, variable);

  SetDnDx( hData );

  pad_->cd();
  hFirst = (TH1D*)hData->Clone("First");
  hData->SetMarkerColor( kBlack );
  hData->GetXaxis()->SetTitle("E [GeV]");
  hData->GetYaxis()->SetTitle("#frac{dN}{dE} [1/GeV]");
  Prepare_1Dplot( hData );
  hData->GetXaxis()->SetRangeUser(300., 1800.);
  hData->Draw( drawoptions );
  drawoptions = "psame";
  legend->AddEntry( hData, legend_info_[ MC_files_[0] ], "ep"  );  

  Get_DetUnfolded( 0 , hJESup, iterations, variable); 
    SetDnDx( hJESup );
    hJESup->SetLineColor( getColor( 2) );
    hJESup->SetMarkerColor( getColor( 2 ) );
    hJESup->SetMarkerStyle( 25  );
    hJESup->SetLineWidth( 2 );


  Get_DetUnfolded( 0 , hJESdown, iterations, variable); 
    SetDnDx( hJESdown );
    hJESdown->SetLineColor( getColor( 3) );
    hJESdown->SetMarkerColor( getColor( 3 ) );
    hJESdown->SetMarkerStyle( 26 );
    hJESdown->SetLineWidth( 2 );
 
  pad_->cd();
  hJESup->Draw( drawoptions );
  hJESdown->Draw( drawoptions );

  pad_->SetLogy();
  pad_->Update();
}









void Unfolder::Plot_Unfolded_Ratio_JESDependence(TPad* & pad_, TString variable, int iterations){
  pad_->cd();

  TLegend* legend = new TLegend( 0.60, 0.7, 0.95, 0.9  );

  TH1D *hFirst, *hData, *hJESup, *hJESdown;
  TString drawoptions = "ephist";

  //-- Detector level energy extraction.
  if( variable == "all") { Get_DetEnergy( -1 , hData); }
  Get_DetEnergy_JESup( -1 , hJESup);
  Get_DetEnergy_JESdown( -1 , hJESdown);


  //-- Unfolding
  Get_DetUnfolded( 0 , hData, iterations, variable);

  pad_->cd();
  hFirst = (TH1D*)hData->Clone("First");
  hData->SetMarkerColor( kBlack );
  hData->GetXaxis()->SetTitle("E [GeV]");
  hData->GetYaxis()->SetTitle("Ratio");

  hData->Divide( hFirst );
  Prepare_1Dplot( hData, pad_ );
  hData->GetXaxis()->SetRangeUser(300., 1800.);
  hData->GetYaxis()->SetRangeUser(0., 4.);
  hData->GetXaxis()->SetTitleOffset( hData->GetXaxis()->GetTitleOffset() * 2 );

  Prepare_1Dplot( hData, pad_ );
  hData->GetYaxis()->SetRangeUser(0., 3.7);
  hData->GetYaxis()->SetNdivisions(504);
  hData->Draw( drawoptions );
  drawoptions = "psame";
  legend->AddEntry( hData, legend_info_[ MC_files_[0] ], "ep"  );  

  Get_DetUnfolded( 0 , hJESup, iterations, variable); 
    hJESup->Divide( hFirst );
    hJESup->SetLineColor( getColor( 2) );
    hJESup->SetMarkerColor( getColor( 2 ) );
    hJESup->SetMarkerStyle( 25  );
    hJESup->SetLineWidth( 2 );


  Get_DetUnfolded( 0 , hJESdown, iterations, variable); 
    hJESdown->Divide( hFirst );
    hJESdown->SetLineColor( getColor( 3) );
    hJESdown->SetMarkerColor( getColor( 3 ) );
    hJESdown->SetMarkerStyle( 26 );
    hJESdown->SetLineWidth( 2 );
 
  pad_->cd();
  hJESup->Draw( drawoptions );
  hJESdown->Draw( drawoptions );

  pad_->Update();
}





void Unfolder::Plot_measured_JESDependence(TPad* & pad_, TString variable, int iterations){
  pad_->cd();

  TLegend* legend = new TLegend( 0.60, 0.7, 0.95, 0.9  );

  TH1D *hFirst, *hData, *hJESup, *hJESdown;
  TString drawoptions = "ephist";

  //-- Detector level energy extraction.
  if( variable == "all") { Get_DetEnergy( -1 , hData); }
  Get_DetEnergy_JESup( -1 , hJESup);
  Get_DetEnergy_JESdown( -1 , hJESdown);


  //-- Unfolding

  SetDnDx( hData );

  hFirst = (TH1D*)hData->Clone("First");

  pad_->cd();
  hData->SetMarkerColor( kBlack );
  hData->SetLineColor( kBlack );
  hData->GetXaxis()->SetTitle("E [GeV]");
  hData->GetYaxis()->SetTitle("#frac{dN}{dE} [1/GeV]");

  Prepare_1Dplot( hData, pad_ );
  hData->GetXaxis()->SetRangeUser(300., 1800.);
  hData->Draw( drawoptions );
  drawoptions = "psame";
  legend->AddEntry( hData, legend_info_[ MC_files_[0] ], "ep"  );  

    SetDnDx( hJESup );
    hJESup->SetLineColor( getColor( 2) );
    hJESup->SetMarkerColor( getColor( 2 ) );
    hJESup->SetMarkerStyle( 25  );
    hJESup->SetLineWidth( 2 );


    SetDnDx( hJESdown );
    hJESdown->SetLineColor( getColor( 3) );
    hJESdown->SetMarkerColor( getColor( 3 ) );
    hJESdown->SetMarkerStyle( 26 );
    hJESdown->SetLineWidth( 2 );
 
  pad_->cd();
  hJESup->Draw( drawoptions );
  hJESdown->Draw( drawoptions );


  pad_->SetLogy();
  pad_->Update();
}






void Unfolder::Plot_measured_Ratio_JESDependence(TPad* & pad_, TString variable, int iterations){
  pad_->cd();

  TLegend* legend = new TLegend( 0.60, 0.7, 0.95, 0.9  );

  TH1D *hFirst, *hData, *hJESup, *hJESdown;
  TString drawoptions = "ephist";

  //-- Detector level energy extraction.
  if( variable == "all") { Get_DetEnergy( -1 , hData); }
  Get_DetEnergy_JESup( -1 , hJESup);
  Get_DetEnergy_JESdown( -1 , hJESdown);

  SetDnDx( hData );

  hData->SetMarkerColor( kBlack );
  hData->GetXaxis()->SetTitle("E [GeV]");
  hData->GetYaxis()->SetTitle("Ratio");
  hData->GetYaxis()->SetRangeUser(0., 3.);
  hData->GetXaxis()->SetRangeUser(300., 1800.);

  pad_->cd();
//  hFirst->GetXaxis()->SetRangeUser(150., 1800.);
//  hFirst->GetYaxis()->SetRangeUser(0., 3.);
  hFirst = (TH1D*)hData->Clone("First");
  hFirst->Divide( hData );
  Prepare_1Dplot( hFirst );
  hFirst->GetXaxis()->SetTitleOffset( hFirst->GetXaxis()->GetTitleOffset() * 2 );
  Prepare_1Dplot( hFirst, pad_); 
  hFirst->GetYaxis()->SetRangeUser(0., 3.7);
  hFirst->GetYaxis()->SetNdivisions(504);
  hFirst->Draw( drawoptions );

  drawoptions = "psame";
  legend->AddEntry( hData, legend_info_[ MC_files_[0] ], "ep"  );  

  SetDnDx( hJESup );
  hJESup->SetLineColor( getColor( 2) );
  hJESup->SetMarkerColor( getColor( 2 ) );
  hJESup->SetMarkerStyle( 25  );
  hJESup->SetLineWidth( 2 );


  SetDnDx( hJESdown );
  hJESdown->SetLineColor( getColor( 3) );
  hJESdown->SetMarkerColor( getColor( 3 ) );
  hJESdown->SetMarkerStyle( 26 );
  hJESdown->SetLineWidth( 2 );
 
  pad_->cd();
  hJESup->Divide( hData );
  hJESup->Draw( drawoptions );

  hJESdown->Divide( hData );
  hJESdown->Draw( drawoptions );


  pad_->Update();
}





void Unfolder::Plot_Unfolded_Ratio_positionDependence(TPad* & pad_, TString variable, int iterations){
  pad_->cd();

  TH1D *hFirst, *hData;
  TString drawoptions = "ephist";

  //-- Detector level energy extraction.
  Get_DetEnergy( -1 , hData);   

  //-- Unfolding
  Get_DetUnfolded( 0 , hData, iterations, variable);

  hFirst = (TH1D*)hData->Clone("First");

  hData->Divide( hFirst );
  hData->GetYaxis()->SetRangeUser(0., 3.);
  hData->SetMarkerColor( kBlack );
  hData->SetLineColor( kBlack );
  hData->GetXaxis()->SetTitle("E [GeV]");
  hData->GetYaxis()->SetTitle("Ratio");
  pad_->cd();

  Prepare_1Dplot( hData, pad_ );
  hData->GetYaxis()->SetRangeUser(0., 3.4);
  hData->GetYaxis()->SetNdivisions(504);

  hData->Draw( drawoptions );
  drawoptions = "psame";

  //-- Model dependence.
  for(int MC_ = 0; MC_ < MC_files_.size() ; MC_++){
    //-- Generator level.
    if( (set_of_tags_["mc_type"])[MC_files_[MC_]] != "position"){ continue; }
    //-- Generator level.

    TH1D* hGen, *hDet;

    //-- Data unfolded.
    TH1D* hData;

    Get_DetEnergy( -1 , hData); 
    Get_DetUnfolded( MC_ , hData, iterations, variable);
    hData->SetLineColor( getColor( MC_+2) );
    hData->SetMarkerColor( getColor( MC_+2 ) );
    hData->SetMarkerStyle( 24 + MC_ );
   
    pad_->cd();
    hData->Divide( hFirst);
    hData->       Draw( drawoptions );
    drawoptions = "psame";
  }
 
  hFirst->GetYaxis()->SetRangeUser(0., 3.);
  pad_->Update();
}









void Unfolder::Plot_DistributionsResponseObject(int files){

  TCanvas *can_;
  TPad *pad_abs_, *pad_ratio_;
  PrepareCanvas(can_, "DistributionsResponseObject");
  SplitCanvas(can_, pad_abs_, pad_ratio_);


  TString drawoptions_abs = "hist";
  TString drawoptions_ratio = "hist";

  TLegend* legend = new TLegend( 0.5, 0.65, 0.93, .95);

  int compared_files;
  TString compared_files_label;
  if( files == 1 ){ 
	compared_files = 1;
	compared_files_label = "leadingSample"; }
  else{ compared_files = MC_files_.size();
	compared_files_label = "allSamples"; }


  for(int file_ = 0; file_ < compared_files; file_++){

    TString filename = MC_files_[ (file_) ];
    TFile* _file = TFile::Open( filename, "read");
    RooUnfoldResponse* response = (RooUnfoldResponse*)_file->Get("response");

    TH1D* hMeasured = (TH1D*)response->Hmeasured();
	SetDnDx( hMeasured );
	hMeasured->GetYaxis()->SetTitle("#frac{dN}{dE}");
    TH1D* hTruth = (TH1D*)response->Htruth();
	SetDnDx( hTruth );
	hTruth->GetYaxis()->SetTitle("#frac{dN}{dE}");

    legend->AddEntry( hMeasured, legend_info_[MC_files_[file_]], "l");

    //-- Absolute values.
    pad_abs_->cd();

    hTruth->SetLineColor( getColor( file_ + 1) );
    hTruth->SetLineStyle( 2 );
    hTruth->SetLineWidth( 2 );
    hMeasured->SetLineColor( getColor( file_ + 1) );
    hMeasured->SetLineWidth( 2 );

    hTruth->Draw( drawoptions_abs );
    drawoptions_abs = "histsame";
    hMeasured->DrawClone( drawoptions_abs );

    pad_abs_->SetLogy();

    //-- Measured to truth ratio.
    pad_ratio_->cd(); 
    hMeasured->Divide( hTruth );
    hMeasured->GetYaxis()->SetRangeUser(0., 2.);
    hMeasured->Draw( drawoptions_ratio );
    hMeasured->GetXaxis()->SetTitle("E [GeV]");
    hMeasured->GetYaxis()->SetTitle("Ratio");
    drawoptions_ratio = "histsame";
 
  }

  pad_abs_->cd();
  legend->Draw();



  can_->SaveAs( TString::Format(folder_ + "/Compare_truth_" + compared_files_label + "_measured_deltaphi_0%i_etaband_0%i.C", static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) )  );
  can_->SaveAs( TString::Format(folder_ + "/Compare_truth_" + compared_files_label + "_measured_deltaphi_0%i_etaband_0%i.pdf", static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_) ) );


}



//-- Compare the etaDiff distributions between the different MC samples.
void Unfolder::Plot_EtaDiff(){

 TString canvasname = TString::Format("Compare_eta_diff_Emin_%i_deltaPhiMax_0%i_etaband_0%i", static_cast<int>( Ethresh_ ), static_cast<int>( 10. * deltaPhiMax_ ), static_cast<int>( 10. * etawidth_ ) );

 TCanvas* can_;
 PrepareCanvas( can_, canvasname );

 TString drawoptions = "hist";

 double min_val, max_val;
 TH1D* hFirst;

 for(int file_ = 0; file_ < MC_files_.size(); file_++){

   // Open file.
   TString filename = MC_files_[ (file_) ];
   TFile* _file = TFile::Open( filename, "read");

   // Prepare scaling.
   TString scalefactor_MC = (set_of_tags_[ "scalefactors" ])[ MC_files_[file_] ];
   double scalefactors_MC_ = scalefactor_MC.Atof() ;

   // Extract histogram, scale and prepare.
   TH1D* hGen = (TH1D*)_file->Get("hGenJet_eta");
   hGen->Scale( 1./scalefactors_MC_ );   

   hGen->SetLineColor( getColor( file_+1 ) );
   hGen->SetLineWidth( 2 );
   hGen->SetLineStyle( file_+1);

   if( file_ == 0 ){
     hFirst = (TH1D*)hGen->Clone("hFirst");
     hFirst->GetXaxis()->SetTitle("#eta_{gen}");
     hFirst->Draw( drawoptions );
   }
   else{
     hGen->Draw( drawoptions );
   }
   First_Plot( hFirst, hGen, file_, max_val, min_val);

   hGen->Draw( drawoptions );

   drawoptions = "histsame";
 }

 hFirst->GetYaxis()->SetRangeUser(0.9 * min_val, 1.1 * max_val);
 can_->SetLogy();

 can_->SaveAs( TString::Format( folder_ + canvasname + ".C") );
 can_->SaveAs( TString::Format( folder_ + canvasname + ".pdf") );



}


















void Unfolder::Plot_Unfolded_JES(int iterations_){
  cout << "\n\n\tUnfolder::Plot_Unfolded_JES" << endl;

  int iterations = iterations_;

  //-- Loop over all MC files.
  for(int file_ = 0; file_ < 1; file_++){


    int color_ = 2;
    TH1D* hGen, *hDet, *hData, *hData_ratio, *hJESup, *hJESup_ratio, *hJESdown, *hJESdown_ratio;
    TString legend_file = legend_info_[MC_files_[file_]];
    TString print_label = printLabel_[MC_files_[file_]];
  
    TCanvas *can;
    TPad *pad_abs_, *pad_ratio_;
    PrepareCanvas( can , "Comparison_Energy_UnfoldedJES");
    SplitCanvas( can, pad_abs_, pad_ratio_ );
    TLegend* leg = new TLegend(0.6, 0.55, 0.95, 0.95);  
    leg->SetFillColor(0);

    double min_val, max_val;

    //-- Absolute value of data.
    pad_abs_->cd();
    Get_DetEnergy(-1, hData);
    Get_DetUnfolded(0, hData, iterations);
    hDet = (TH1D*)hData->Clone( TString::Format( "Unfolded") );

    SetDnDx( hData );
    hData->GetXaxis()->SetTitle("E [GeV]");
    hData->GetYaxis()->SetTitle("#frac{dN}{dE}");
    Prepare_1Dplot( hData, pad_abs_);
    hData->	Draw("ephist");
    hData->GetXaxis()->SetNdivisions(504);
    First_Plot( hData, hData, 0, min_val, max_val);
    leg->	AddEntry( hData, TString::Format("Unfolded Data, %i it.", iterations), "p");

    //-- Ratio of data.
    pad_ratio_->cd();
    hData_ratio = (TH1D*)hData->Clone("hData_ratio");
    Prepare_1Dplot( hData_ratio, pad_ratio_);
    hData_ratio->Divide( hData );
    hData_ratio->GetYaxis()->SetRangeUser(0.,2.);
    hData_ratio->GetYaxis()->SetTitle("Ratio");
    hData_ratio->Draw("ephist");

    //-- JES up
    Get_DetEnergy_JESup(-1, hJESup);
    Get_DetUnfolded(0, hJESup, iterations);
    hJESup->SetName( "JESup" );

    SetDnDx( hJESup );
    hJESup->SetLineColor( getColor(color_) );
    hJESup->SetLineStyle( color_++ );
    hJESup->SetLineWidth( 2 );

    pad_abs_->cd();
    hJESup->DrawClone("histsame");
    First_Plot( hData, hJESup, 1, min_val, max_val);
    leg->AddEntry( hJESup, TString::Format("Unfolded JES+, %i it.", iterations), "lp" );

    pad_ratio_->cd();
    hJESup->Divide( hData );
    hJESup->Draw("histsame");

    //-- JES up
    Get_DetEnergy_JESdown(-1, hJESdown);
    Get_DetUnfolded(0, hJESdown, iterations);
    hJESdown->SetName( "hJESup" );

    SetDnDx( hJESdown );
    hJESdown->SetLineColor( getColor(color_) );
    hJESdown->SetLineStyle( color_++ );
    hJESdown->SetLineWidth( 2 );

    pad_abs_->cd();
    hJESdown->DrawClone("histsame");
    First_Plot( hData, hJESdown, 1, min_val, max_val);
    leg->AddEntry( hJESdown, TString::Format("Unfolded JES-, %i it.", iterations), "lp" );

    pad_ratio_->cd();
    hJESdown->Divide( hData );
    hJESdown->Draw("histsame");
    

    pad_abs_->cd();
    leg->Draw();
    pad_abs_->SetLogy();
    can->Update();
  
    can->SaveAs(folder_ + "/Compare_JES" + print_label + cutsetup_ + ".C");
    can->SaveAs(folder_ + "/Compare_JES" + print_label + cutsetup_ + ".pdf");  
  }  
}





void Unfolder::Plot_Measured_JES(){
  cout << "\n\n\tUnfolder::Plot_Unfolded_JES" << endl;

  //-- Loop over all MC files.
  for(int file_ = 0; file_ < 1; file_++){


    int color_ = 2;
    TH1D* hGen, *hDet, *hData, *hData_ratio, *hJESup, *hJESup_ratio, *hJESdown, *hJESdown_ratio;
    TString legend_file = legend_info_[MC_files_[file_]];
    TString print_label = printLabel_[MC_files_[file_]];
  
    TCanvas *can;
    TPad *pad_abs_, *pad_ratio_;
    PrepareCanvas( can , "Comparison_Energy_UnfoldedJES");
    SplitCanvas( can, pad_abs_, pad_ratio_ );
    TLegend* leg = new TLegend(0.6, 0.55, 0.95, 0.95);  
    leg->SetFillColor(0);

    double min_val, max_val;

    //-- Absolute value of data.
    pad_abs_->cd();
    Get_DetEnergy(-1, hData);
    hDet = (TH1D*)hData->Clone( TString::Format( "Unfolded") );

    SetDnDx( hData );
    hData->GetXaxis()->SetTitle("E [GeV]");
    hData->GetYaxis()->SetTitle("#frac{dN}{dE}");
    hData->	Draw("ephist");
    hData->GetXaxis()->SetNdivisions(504);
    First_Plot( hData, hData, 0, min_val, max_val);
    leg->	AddEntry( hData, "Measured Data", "p");

    //-- Ratio of data.
    pad_ratio_->cd();
    hData_ratio = (TH1D*)hData->Clone("hData_ratio");
    hData_ratio->Divide( hData );
    hData_ratio->GetYaxis()->SetRangeUser(0.,2.);
    hData_ratio->GetYaxis()->SetTitle("Ratio");
    hData_ratio->Draw("ephist");

    //-- JES up
    Get_DetEnergy_JESup(-1, hJESup);
//    Get_DetUnfolded(file_, hDet, iterations);
    hJESup->SetName( "JESup" );

    SetDnDx( hJESup );
    hJESup->SetLineColor( getColor(color_) );
    hJESup->SetLineStyle( color_++ );
    hJESup->SetLineWidth( 2 );

    pad_abs_->cd();
    hJESup->DrawClone("histsame");
    First_Plot( hData, hJESup, 1, min_val, max_val);
    leg->AddEntry( hJESup, "JES +", "lp" );

    pad_ratio_->cd();
    hJESup->Divide( hData );
    hJESup->Draw("histsame");

    //-- JES down
    Get_DetEnergy_JESdown(-1, hJESdown);
    hJESdown->SetName( "hJESup" );

    SetDnDx( hJESdown );
    hJESdown->SetLineColor( getColor(color_) );
    hJESdown->SetLineStyle( color_++ );
    hJESdown->SetLineWidth( 2 );

    pad_abs_->cd();
    hJESdown->DrawClone("histsame");
    First_Plot( hData, hJESdown, 1, min_val, max_val);
    leg->AddEntry( hJESdown, "JES -", "lp" );

    pad_ratio_->cd();
    hJESdown->Divide( hData );
    hJESdown->Draw("histsame");
    

    pad_abs_->cd();
    leg->Draw();
    pad_abs_->SetLogy();
    can->Update();
  
    can->SaveAs(folder_ + "/Compare_JES_measured" + print_label + cutsetup_ + ".C");
    can->SaveAs(folder_ + "/Compare_JES_measured" + print_label + cutsetup_ + ".pdf");  
  }  
}




void Unfolder::Fit_unfolded_distribution( TF1* &analytical){

 TH1D* hData;

 Get_DetEnergy(-1, hData);
 Get_DetUnfolded( 0, hData);

 SetDnDx( hData );

 TCanvas* can_fit;
 PrepareCanvas(can_fit, "Fit_unfolded_data");

 hData->Draw("phist");
 hData->GetXaxis()->SetNdivisions(504);

 double dataint = hData->Integral();

// analytical = new TF1("Analytical", "  [0] * exp( [1] * x) * (x < 1000.) + [2] * exp( [3] * x) * (x > 1000.)", 150., 2000.);

 //--- 150-1000.
 // 1.22009e+00
 // -7.38385e-03
 //-- 
 // 1.59985e-01
 // -4.68686e-03
 //--
 analytical = new TF1("Analytical", "  [0] * exp( [1] * x) * (x < 894.07) + [2] * exp( [3] * x) * (x > 894.07)", 150., 2000.);

 analytical->SetLineColor( getColor( 2) );

 analytical->SetParLimits(0, 1.2*dataint, 1.3*dataint);
 analytical->SetParLimits(2, 0.15*dataint, 0.17*dataint);

 analytical->SetParLimits(1, -8.0e-03, -7.0e-03);
 analytical->SetParLimits(3, -5.0e-03, -4.0e-03);

 hData->Fit( analytical, "", "", 150., 2000. );

 can_fit->SetLogy();
 
 can_fit->SaveAs(folder_ + "Fit_to_unfolded_data.C");
 can_fit->SaveAs(folder_ + "Fit_to_unfolded_data.pdf");
}






void Unfolder::Plot_JER(){
 cout << "Unfolder::Plot_JER" << endl;

 for(int file = 0; file < MC_files_.size(); file++){

   TFile* _file0 = TFile::Open( MC_files_[file] , "READ");
   TCanvas* can = new TCanvas( "canvas_JER", "canvas_JER", 800., 800. );

   TH1D* hJER = (TH1D*)_file0->Get("hJER");
     hJER->SetLineWidth( 2 );
     hJER->GetXaxis()->SetTitle("JES");
     hJER->GetXaxis()->SetRangeUser(-1., 1.5);

   TH1D* hJER_had = (TH1D*)_file0->Get("hJER_had");
     hJER_had->SetLineWidth( 2 );
     hJER_had->SetLineStyle( 2 );
     hJER_had->SetLineColor( getColor( 2 ) );

   TH1D* hJER_em = (TH1D*)_file0->Get("hJER_em");
     hJER_em->SetLineWidth( 2 );
     hJER_em->SetLineStyle( 3);
     hJER_em->SetLineColor( getColor( 3 ) );

   TH1D* hJER_other = (TH1D*)_file0->Get("hJER_other");
     hJER_other->SetLineWidth( 2 );
     hJER_other->SetLineStyle( 4 );
     hJER_other->SetLineColor( getColor( 4 ) );

   hJER->Draw("hist");
   hJER_had->Draw("histsame");
   hJER_em->Draw("histsame");
   hJER_other->Draw("histsame");

    TLegend *leg = new TLegend( 
	1. - can->GetRightMargin() - 0.25,
	1. - can->GetTopMargin() - 0.2,
	1. - can->GetRightMargin(),
	1. - can->GetTopMargin() );

    leg->SetFillStyle( 0 );
    leg->SetBorderSize( 0 );

    leg->AddEntry( hJER, "All jets", "l");
    leg->AddEntry( hJER_had, "Had. jets", "l");
    leg->AddEntry( hJER_em, "EM. jet", "l");
    leg->AddEntry( hJER_other, "other", "l");

    leg->Draw();

    TLatex Tl; Tl.SetTextFont(43); Tl.SetTextSize(35);
    Tl.DrawText( 
	0.075,
	hJER->GetMaximum() * 0.5, 
	legend_info_[ MC_files_[file]] );

    TString setup;
    if( MC_files_[file].Contains("isolated") ) setup = "isolated";
    else if( MC_files_[file].Contains("calibrated") ) setup = "calibrated";



    TString savenameC = TString::Format( folder_ + "JES_" + printLabel_[ MC_files_[file] ] + "_" + setup + ".C") ;
    TString savenamepdf = TString::Format( folder_ + "JES_" + printLabel_[ MC_files_[file] ] + "_" + setup + ".pdf") ;

    can->SaveAs( savenameC );
    can->SaveAs( savenamepdf);

  }


}



void Unfolder::Plot_Isolation(){


 for(int file = 0; file < MC_files_.size(); file++){

   TFile* _file0 = TFile::Open( MC_files_[file] , "READ");

  //-- Plot 2D EI versus nTowers.

  TCanvas *can_2D;
  PrepareCanvas_2D(can_2D, "EI_vs_nTowers");

  TH2D* EI_vs_nTowers = (TH2D*)_file0->Get("hIsolationEnergy");
  EI_vs_nTowers->GetXaxis()->SetRangeUser(0., 4.);
  EI_vs_nTowers->GetXaxis()->SetNdivisions(405);
  EI_vs_nTowers->GetXaxis()->SetTitle("N_{tower}");

  EI_vs_nTowers->GetYaxis()->SetRangeUser(0., 2000.);
  EI_vs_nTowers->GetYaxis()->SetTitle("E_{I} [GeV]");
  EI_vs_nTowers->GetYaxis()->SetLabelSize( 0.75 * EI_vs_nTowers->GetYaxis()->GetLabelSize() );
  EI_vs_nTowers->GetYaxis()->SetTitleOffset( 1.2 * EI_vs_nTowers->GetYaxis()->GetTitleOffset() );

  EI_vs_nTowers->GetZaxis()->SetRangeUser( 0.9 * GetMinimumValue(EI_vs_nTowers), 1.  );

  EI_vs_nTowers->Scale(1./EI_vs_nTowers->Integral() );
  EI_vs_nTowers->Draw("colz");

  can_2D->SetLogz();
  can_2D->SaveAs(folder_ + "/EI_vs_nTowers.C");
  can_2D->SaveAs(folder_ + "/EI_vs_nTowers.pdf");


  //--- Compare distributions per number of tower.
  PrepareCanvas(can_2D, "EI_vs_nTowers_projected");

  TH1D* h3T = (TH1D*)EI_vs_nTowers->ProjectionY("nTower_3", 4, 4);
  h3T->SetLineColor( getColor( 3 ) );
  h3T->SetLineStyle( 3 );
  h3T->SetLineWidth( 2 );
  h3T->GetXaxis()->SetTitleOffset( h3T->GetXaxis()->GetTitleOffset() * 1.1 );
  h3T->GetYaxis()->SetTitleOffset( h3T->GetYaxis()->GetTitleOffset() * 1.6 );
  h3T->GetXaxis()->SetTitleSize( h3T->GetXaxis()->GetTitleSize() * 0.8 );
  h3T->GetYaxis()->SetTitle("#frac{1}{N_{tot.}}.#frac{dN}{dE}");
  h3T->GetXaxis()->SetNdivisions(504);
  h3T->Draw("hist");

  TH1D* h2T = (TH1D*)EI_vs_nTowers->ProjectionY("nTower_2", 3, 3);
  h2T->SetLineColor( getColor( 2 ) );
  h2T->SetLineStyle( 2 );
  h2T->SetLineWidth( 2 );
  h2T->Draw("histsame");

  TH1D* h1T = (TH1D*)EI_vs_nTowers->ProjectionY("nTower_1", 2, 2);
  h1T->Draw("histsame");

  TLegend *leg_proj = new TLegend( 0.7, 0.8, 0.9, 0.9);
  leg_proj->SetFillColor( 0 );
  leg_proj->AddEntry( h1T, "1 sector", "l");
  leg_proj->AddEntry( h2T, "2 sector", "l");
  leg_proj->AddEntry( h3T, "3 sector", "l");
  leg_proj->Draw();
  
  can_2D->SetLogy();
  can_2D->SaveAs(folder_ + "EI_per_nTower.C");
  can_2D->SaveAs(folder_ + "EI_per_nTower.pdf"); 


  //-- Plot 1D.

  TCanvas *can_1D_relative = new TCanvas( "EI_as_function_of_energy", "EI_as_function_of_energy",800, 500 );
  PrepareCanvas(can_1D_relative, "EI_as_function_of_energy");
  can_1D_relative->SetLogx();
  can_1D_relative->SetLogy();

  TLegend *leg = new TLegend( 0.7, 0.8, 0.9, 0.9);
  leg->SetFillColor( 0 );

  //-- Versus detector level energy.
  TH1D* EI_vs_Edet = (TH1D*)_file0->Get("hIsolationEnergy_Edet");
  EI_vs_Edet->Scale( 1./EI_vs_Edet->Integral() );
  EI_vs_Edet->GetXaxis()->SetNdivisions(405);
  EI_vs_Edet->GetXaxis()->SetTitle("E_{I}/E");
  EI_vs_Edet->GetXaxis()->SetTitleSize( EI_vs_Edet->GetXaxis()->GetTitleSize() * 0.8 );
  EI_vs_Edet->GetXaxis()->SetTitleOffset( EI_vs_Edet->GetXaxis()->GetTitleOffset() * 1.1 );

  EI_vs_Edet->GetYaxis()->SetRangeUser( GetMinimumValue( EI_vs_Edet )*0.9, 1.);
  EI_vs_Edet->GetYaxis()->SetTitle("#frac{1}{N}.#frac{dN}{dE}");
  EI_vs_nTowers->GetYaxis()->SetTitleOffset( 0.7 * EI_vs_nTowers->GetYaxis()->GetTitleOffset() );
  EI_vs_Edet->Draw("hist");  

  
  //-- Versus generator level energy.
  TH1D* EI_vs_Egen = (TH1D*)_file0->Get("hIsolationEnergy_Egen");
  EI_vs_Egen->Scale( 1./EI_vs_Egen->Integral() );
  EI_vs_Egen->GetXaxis()->SetNdivisions(405);
  EI_vs_Egen->GetXaxis()->SetTitle("E_{I/E}");
  EI_vs_Egen->GetYaxis()->SetRangeUser( GetMinimumValue( EI_vs_Egen )*0.9, 1.);
  EI_vs_Egen->GetYaxis()->SetTitle("N");
  EI_vs_Egen->SetLineColor( kRed );
  EI_vs_Egen->SetLineStyle( 2 );
  EI_vs_Egen->Draw("histsame");  

  leg->AddEntry( EI_vs_Edet, "E_{I}/E_{det}", "l" );
  leg->AddEntry( EI_vs_Egen, "E_{I}/E_{gen}", "l" );
  leg->Draw();

  can_1D_relative->SaveAs(folder_ + "EI_relative.C");
  can_1D_relative->SaveAs(folder_ + "EI_relative.pdf");

  //-- Absolute value.

  TCanvas *can_1D_absolute = new TCanvas( "EI_absolute", "EI_absolute", 800, 500 );
  PrepareCanvas(can_1D_absolute, "EI_absolute");
  can_1D_absolute->SetLogx();
  can_1D_absolute->SetLogy();

  TH1D* EI = (TH1D*)_file0->Get("hIsolationEnergy_1D_log");

  EI->Scale( 1./EI->Integral() );
  EI->GetXaxis()->SetNdivisions(405);
  EI->GetXaxis()->SetTitle("E_{I} [GeV]");
  EI->GetYaxis()->SetRangeUser( GetMinimumValue( EI )*0.9, 1.);
  EI->GetYaxis()->SetTitle("#frac{1}{N}.#frac{dN}{dE}");
  EI->Draw("hist");  

  TLatex Tl; Tl.SetTextFont(43); Tl.SetTextSize(35);
  Tl.DrawText( 
	1e-2,
	0.5 , 
	legend_info_[ MC_files_[file]] );

  TString setup;
  if( MC_files_[file].Contains("isolated") ) setup = "isolated";
  else if( MC_files_[file].Contains("calibrated") ) setup = "calibrated";

  TString savenameC = TString::Format( folder_ + "EI_" + printLabel_[ MC_files_[file] ] + "_" + setup + ".C") ;
  TString savenamepdf = TString::Format( folder_ + "EI_" + printLabel_[ MC_files_[file] ] + "_" + setup + ".pdf") ;

  can_1D_absolute->SetLogx();

  can_1D_absolute->SaveAs( savenameC );
  can_1D_absolute->SaveAs( savenamepdf);

  
  //-- Relative value: EGEN
  cout << "//-- Relative value: EGEN" << endl;

  can_1D_relative = new TCanvas( "EI_relative", "EI_relative", 800, 500 );
  PrepareCanvas(can_1D_relative, "EI_relative");
  can_1D_relative->SetLogx();
  can_1D_relative->SetLogy();

  TH1D* EI_egen = (TH1D*)_file0->Get("hIsolationEnergy_Egen_log");

  EI_egen->Scale( 1./EI_egen->Integral() );
  EI_egen->GetXaxis()->SetNdivisions(405);
  EI_egen->GetXaxis()->SetTitle("E_{I}/E_{gen}");
  EI_egen->GetYaxis()->SetRangeUser( GetMinimumValue( EI )*0.9, 1.);
  EI_egen->GetYaxis()->SetTitle("#frac{1}{N}.#frac{dN}{dE}");
  EI_egen->Draw("hist");  

  Tl.SetTextFont(43); Tl.SetTextSize(35);
  Tl.DrawText( 
	1e-2,
	0.5 , 
	legend_info_[ MC_files_[file]] );

  savenameC = TString::Format( folder_ + "EI_gen_" + printLabel_[ MC_files_[file] ] + "_" + setup + ".C") ;
  savenamepdf = TString::Format( folder_ + "EI_gen_" + printLabel_[ MC_files_[file] ] + "_" + setup + ".pdf") ;

  can_1D_relative->SetLogx();

  can_1D_relative->SaveAs( savenameC );
  can_1D_relative->SaveAs( savenamepdf);

  //-- Relative value: EDET
  cout << "//-- Relative value: EDET" << endl;

  can_1D_relative = new TCanvas( "EI_relative", "EI_relative", 800, 500 );
  PrepareCanvas(can_1D_relative, "EI_relative");
  can_1D_relative->SetLogx();
  can_1D_relative->SetLogy();

  TH1D* EI_edet = (TH1D*)_file0->Get("hIsolationEnergy_Edet_log");

  EI_edet->Scale( 1./EI_edet->Integral() );
  EI_edet->GetXaxis()->SetNdivisions(405);
  EI_edet->GetXaxis()->SetTitle("E_{I}/E_{det}");
  EI_edet->GetYaxis()->SetRangeUser( GetMinimumValue( EI )*0.9, 1.);
  EI_edet->GetYaxis()->SetTitle("#frac{1}{N}.#frac{dN}{dE}");
  EI_edet->Draw("hist");  

  Tl.SetTextFont(43); Tl.SetTextSize(35);
  Tl.DrawText( 
	1e-2,
	0.5 , 
	legend_info_[ MC_files_[file]] );

  savenameC = TString::Format( folder_ + "EI_det_" + printLabel_[ MC_files_[file] ] + "_" + setup + ".C") ;
  savenamepdf = TString::Format( folder_ + "EI_det_" + printLabel_[ MC_files_[file] ] + "_" + setup + ".pdf") ;

  can_1D_relative->SetLogx();

  can_1D_relative->SaveAs( savenameC );
  can_1D_relative->SaveAs( savenamepdf);
  }
}




void Unfolder::ScaleToData( int MC, double &scale ){
  cout << "\n\nUnfolder::ScaleToData(" << MC << ")\t" << MC_files_[MC] << endl;

  TString norm_ = "150";

  TFile *_file0 = TFile::Open( MC_files_[MC], "Read"); 
  
  TTree *tMC = (TTree*)_file0->Get("useful_numbers");
  int castor_150GeV_events_MC, castor_events_MC, total_events_MC;

  if( (set_of_tags_["mc_type"])[MC_files_[MC]] != "shift_MPI_or_Tune" ){
    tMC->SetBranchAddress("castor_" + norm_ + "GeV_events_nocuts",&castor_150GeV_events_MC);
    if( !castor_150GeV_events_MC ){ tMC->SetBranchAddress("castor_" + norm_ + "GeV_events",&castor_150GeV_events_MC); }
    tMC->SetBranchAddress("castor_events_nocuts",&castor_events_MC);
    if( !castor_events_MC ){ tMC->SetBranchAddress("castor_events",&castor_events_MC); }
    tMC->SetBranchAddress("total_events", &total_events_MC);
  }

  else{
    tMC->SetBranchAddress("castor_" + norm_ + "GeV_events",&castor_150GeV_events_MC);
    tMC->SetBranchAddress("castor_events",&castor_events_MC);
    tMC->SetBranchAddress("succesful_events", &total_events_MC);
  }

  tMC->GetEntry(0);

  TFile *_fileData = TFile::Open( datafile_, "Read"); 
  
  TTree *tData = (TTree*)_fileData->Get("useful_numbers");
  int castor_150GeV_events_data,  castor_events_data, total_events_data;
  tData->SetBranchAddress("castor_" + norm_ + "GeV_events",&castor_150GeV_events_data);
  tData->SetBranchAddress("castor_events",&castor_events_data);
  tData->SetBranchAddress("total_events", &total_events_data);
  tData->GetEntry(0);


//  scale = double(castor_events_data)  /  double(castor_events_MC) ;
//  if( MC_files_[MC].Contains("Pythia8") ){ scale = double(castor_events_data)  /  double(castor_events_MC) ; }

  scale = double(total_events_data)  /  double(total_events_MC) ;
  scale = double(castor_events_data)  /  double(castor_events_MC) ;
  scale = double(castor_150GeV_events_data)  /  double(castor_150GeV_events_MC) ;

  cout << "\nscale is " << scale << "\tfile\t" << MC_files_[MC] << "\t\n" 
	<< "\t" << castor_150GeV_events_data << "\t" << double(castor_150GeV_events_data) 
	<< "\t" << castor_events_data << "\t" << double(castor_events_data) << endl 
	<< "\t" << castor_150GeV_events_MC << "\t" << double(castor_150GeV_events_MC) 
	<< "\t" << castor_events_MC << "\t" << double(castor_events_MC) 
	<< endl << endl << endl;
}




















void Unfolder::Plot_GenLevels_(){
  cout << "\nUnfolder::Plot_GenLevels_()" << endl;

  TCanvas *can = new TCanvas();

  gStyle->SetLineWidth(2);


cout << "Opening" << endl;
  TFile *_file0 = TFile::Open("/user/avanspil/Castor_Analysis/ak5ak5_NewGeo_unfold_Emin_150.000000_deltaPhiMax_0.500000_etaband_0.000000_all_matchE.root", "READ");

  TFile *_file1 = TFile::Open("/user/avanspil/Castor_Analysis/ak5ak5_displaced_pythia6Z2star_unfold_Emin_150.000000_deltaPhiMax_0.500000_etaband_0.000000_all_matchE.root", "READ");
  TFile *_file2 = TFile::Open("/user/avanspil/Castor_Analysis/ak5ak5_displaced_pythia6Z2star_unfold_Emin_150.000000_deltaPhiMax_0.500000_etaband_0.000000_all_matchE.root", "READ");

  TFile *_file3 = TFile::Open("/user/avanspil/Castor_Analysis/ak5ak5_Pythia6Z2star_up_unfold_Emin_150.000000_deltaPhiMax_0.500000_etaband_0.000000_all_matchE.root", "READ");
  TFile *_file4 = TFile::Open("/user/avanspil/Castor_Analysis/ak5ak5_Pythia6Z2star_down_unfold_Emin_150.000000_deltaPhiMax_0.500000_etaband_0.000000_all_matchE.root", "READ");


  TFile *_file5 = TFile::Open("/user/avanspil/Castor_Analysis/ak5ak5_displaced_EPOS_unfold_Emin_150.000000_deltaPhiMax_0.500000_etaband_0.000000_all_matchE.root", "Read");

  TFile *_file6 = TFile::Open("/user/avanspil/Castor_Analysis/ak5ak5_Pythia84C_nonDisplaced_unfold_Emin_150.000000_deltaPhiMax_0.500000_etaband_0.000000_all_matchE.root", "Read");

  TFile *_file7 = TFile::Open("/user/avanspil/Castor_Analysis/ak5ak5_Pythia84C_displaced_unfold_Emin_150.000000_deltaPhiMax_0.500000_etaband_0.000000_all_matchE.root", "Read");

  TH1D* h0, *h1, *h2, *h3, *h4, *h5, *h6, *h7;
  TLegend* leg = new TLegend(0, 0.6, 0.4, 1.);

  h0 = (TH1D*)_file0->Get("hGenJet_energy_noCuts");	leg->AddEntry( h0, "NewGeo", "l");
  h1 = (TH1D*)_file1->Get("hGenJet_energy_noCuts");	leg->AddEntry( h1, "displace p6", "l");
  h2 = (TH1D*)_file2->Get("hGenJet_energy_noCuts");	leg->AddEntry( h2, "displaced", "l");
  h3 = (TH1D*)_file3->Get("hGenJet_energy_noCuts");	leg->AddEntry( h3, "up", "l");
  h4 = (TH1D*)_file4->Get("hGenJet_energy_noCuts");	leg->AddEntry( h4, "down", "l");

  h5 = (TH1D*)_file5->Get("hGenJet_energy_noCuts");	leg->AddEntry( h5, "EPOS", "l");
  h6 = (TH1D*)_file6->Get("hGenJet_energy_noCuts");	leg->AddEntry( h6, "nom. p8", "l");
  h7 = (TH1D*)_file7->Get("hGenJet_energy_noCuts");	leg->AddEntry( h7, "dis. p8", "l");

  int events0, events1, events2, events3, events4, events5, events6, events7;
  TTree *tree0, *tree1, *tree2, *tree3, *tree4, *tree5, *tree6, *tree7;

  tree0 = (TTree*)_file0->Get("useful_numbers"); if( !tree0 ) cout << "Nothin" << endl;
  tree1 = (TTree*)_file1->Get("useful_numbers");
  tree2 = (TTree*)_file2->Get("useful_numbers");
  tree3 = (TTree*)_file3->Get("useful_numbers");
  tree4 = (TTree*)_file4->Get("useful_numbers");
  tree5 = (TTree*)_file5->Get("useful_numbers");
  tree6 = (TTree*)_file6->Get("useful_numbers");
  tree7 = (TTree*)_file7->Get("useful_numbers");


      tree0->SetBranchAddress("total_events_nocuts", &events0);
      tree1->SetBranchAddress("total_events_nocuts", &events1);
      tree2->SetBranchAddress("total_events_nocuts", &events2);
      tree3->SetBranchAddress("total_events_nocuts", &events3);
      tree4->SetBranchAddress("total_events_nocuts", &events4);
      tree5->SetBranchAddress("total_events_nocuts", &events5);
      tree6->SetBranchAddress("total_events_nocuts", &events6);
      tree7->SetBranchAddress("total_events_nocuts", &events7);

tree0->GetEntry(0);
tree1->GetEntry(0);
tree2->GetEntry(0);
tree3->GetEntry(0);
tree4->GetEntry(0);
tree5->GetEntry(0);
tree6->GetEntry(0);
tree7->GetEntry(0);


cout	<< "\n0\t" << events0
	<< "\n1\t" << events1
        << "\n2\t" << events2
        << "\n3\t" << events3
        << "\n4\t" << events4
        << "\n5\t" << events5
        << "\n6\t" << events6
        << "\n7\t" << events7 << endl;


  h0->Scale( 1./events0 );
  h1->Scale( 1./events1 );
  h2->Scale( 1./events2 );
  h3->Scale( 1./events3 );
  h4->Scale( 1./events4 );
  h5->Scale( 1./events5 );
  h6->Scale( 1./events6 );
  h7->Scale( 1./events7 );

  h0->Divide( h2 );
  h0->GetYaxis()->SetRangeUser(1e-4, 2.);
  h1->SetLineWidth( 3 );
  h0->Draw("hist");

  h1->Divide( h2 );
  h1->SetLineColor( getColor(2) );
  h1->SetLineStyle( 2 );
  h1->SetLineWidth( 3 );
  h1->Draw("histsame");
/*
  h2->Divide( h2 );
  h2->SetLineColor( getColor(3) );
  h2->SetLineStyle( 3 );
  h2->Draw("histsame");
*/
  h3->Divide( h2 );
  h3->SetLineColor( getColor(4) );
  h3->SetLineStyle( 4 );
  h3->SetLineWidth( 3 );
  h3->Draw("histsame");

  h4->Divide( h2 );
  h4->SetLineColor( getColor(5) );
  h4->SetLineStyle( 2 );
  h4->SetLineWidth( 3 );
  h4->Draw("histsame");

  h6->Divide( h2 );
  h6->SetLineColor( getColor(6) );
  h6->SetLineStyle( 2 );
  h6->SetLineWidth( 3 );
  h6->Draw("histsame");

  h7->Divide( h2 );
  h7->SetLineColor( getColor(7) );
  h7->SetLineStyle( 2 );
  h7->SetLineWidth( 3 );
  h7->Draw("histsame");

  h5->Divide( h2 );
  h5->SetLineColor( getColor(8) );
  h5->SetLineStyle( 3 );
  h5->SetLineWidth( 3 );
  h5->Draw("histsame");

  h2->Divide( h2 );
  h2->SetLineColor( getColor(3) );
  h2->SetLineStyle( 3 );
  h2->SetLineWidth( 3 );
  h2->Draw("histsame");


  leg->Draw();
//  can->SetLogy();
  can->SaveAs("Three_genPlots.C");
  can->SaveAs("Three_genPlots.pdf");

}








void Unfolder::Plot_CastorJetID(){
cout << "Plot jet id" << endl;

  vector<TString> jetIDs;
  jetIDs.push_back("ehad");
  jetIDs.push_back("eem");
  jetIDs.push_back("nTowers");
  jetIDs.push_back("sigmaz");
  jetIDs.push_back("width");
  jetIDs.push_back("depth");
  jetIDs.push_back("fhot");
  jetIDs.push_back("fem");
  jetIDs.push_back("phi");
  jetIDs.push_back("hCastorJet_energy");

  vector<TString> files;
  
  files.push_back( datafile_ );
  
  for(int file_ = 0 ; file_ < MC_files_.size(); file_++){
    if(	(set_of_tags_["mc_type"])[MC_files_[file_]] != "model" && 
	(set_of_tags_["mc_type"])[MC_files_[file_]] != "actual"){ continue; }
    files.push_back( MC_files_[file_] );
  }

  for(vector<TString>::iterator var = jetIDs.begin(); var != jetIDs.end(); ++var){
    TString variable = *var;
    TCanvas *can;
    PrepareCanvas( can, TString::Format( "can_" + variable ) );
    TString drawoptions = "";
    int color = 1;
    double ticksize = gStyle->GetTickLength();

    double xmin, ymin;  
    
    if( (*var) == "depth" || (*var) == "hCastorJet_energy"  || (*var) == "fem" ){ 
      xmin = 0.30; 
      ymin = can->GetBottomMargin() + ticksize;
    }

    else if( (*var) != "ehad" && (*var) != "eem"  ){ 
       xmin = 0.35; 
       ymin = can->GetBottomMargin();
    }
    
    else{
      xmin = 1. - can->GetRightMargin() - 0.3;
      ymin = 1. - can->GetTopMargin() -0.3;
    }
    TLegend *leg = new TLegend( xmin, ymin, xmin + 0.3, ymin + 0.3);
    leg->SetFillStyle( 0 );

      for( int file_ = 0; file_ < files.size(); file_++){

	TFile* _file = TFile::Open( files[file_], "Read");
	TH1D* hist = (TH1D*)_file->Get( variable );

	hist->Scale( 1./hist->Integral() );
	hist->SetLineColor( getColor( color ) );
	
	double lowestedge = hist->GetXaxis()->GetBinLowEdge( 1 );
	double highestedge = hist->GetXaxis()->GetBinUpEdge( hist->GetNbinsX() );
	if( (*var) == "depth" ){ hist->GetXaxis()->SetNdivisions( 204  ); }
	hist->SetLineStyle ( color++ );
	hist->SetLineWidth ( 2 );
	Prepare_1Dplot( hist );
	
	if( (*var) == "phi" ){ hist->GetYaxis()->SetRangeUser(0., 0.14); }
	if( ( files[file_] ).Contains("ata") ) { hist->Draw( "p" + drawoptions ); }
	else{ hist->Draw( "hist" + drawoptions ); }
	drawoptions = "same";

	if( ( files[file_] ).Contains("Nom") ){ leg->AddEntry( hist, "Pythia6 (Z2*) [SL]" , "l" ); }
	else if( files[file_] == datafile_ ){ leg->AddEntry( hist, "Data" , "p" ); }
	else{ 	leg->AddEntry( hist, legend_info_gen_[ files[file_] ] , "l" ); }	
	
	if( (*var) == "fem" ){ 
	  hist->GetXaxis()->SetTitle( "E_{EM}/E_{tot}" );
	  hist->GetYaxis()->SetTitle( "1/N_{tot} . N/#Delta f_{EM}" );
	}
	else if( (*var) == "depth" ){ 
	  hist->GetXaxis()->SetTitle( "<z> [mm]" );
	  hist->GetYaxis()->SetTitle( "1/N_{tot} . N/#Delta <z>" );
	}

     }

     leg->Draw();

     if( (*var) != "phi" ){ can->SetLogy(); }

     can->SaveAs( TString::Format( folder_ + "CastorJetID_" + variable + ".pdf" ) );
     can->SaveAs( TString::Format( folder_ + "CastorJetID_" + variable + ".C" ) );
  }
}




















/*
To lazy to build a separate makefile/executable for just a single counting function.
*/
void Unfolder::Cut_efficiency(){

   TFile* _file0 = TFile::Open("/localgrid/avanspil/Castor_Analysis/CMSSW_4_2_10_patch2/src/UACastor/CastorTree/Analysis/20160318_ak5_Data_STRIPPED_TREE.root", "Read");

   int total_events, survive_cleaning, survive_BSC, survive_HF, survive_vertex;
   int count_total = 0, count_cleaning = 0, count_BSC = 0, count_HF = 0, count_vertex = 0;

   TTree *tree_numbers = (TTree*)_file0->Get("useful_numbers");

   //== All events.
   TString branch_total = "total_events";
   tree_numbers->SetBranchAddress(branch_total, &total_events);  

   //== Passed cleaning cuts.
   TString branch_cleaning = "survive_cleaning";
   tree_numbers->SetBranchAddress(branch_cleaning, &survive_cleaning);   

   //== Passed BSC.
   TString branch_BSC = "survive_BSC";
   tree_numbers->SetBranchAddress(branch_BSC, &survive_BSC);
  
   //== Passed HF.
   TString branch_HF = "survive_HF";
   tree_numbers->SetBranchAddress(branch_HF, &survive_HF);
  
   //== Passed vertex.
   TString branch_vertex = "survive_vertex";
   tree_numbers->SetBranchAddress(branch_vertex, &survive_vertex);


   for( int entry = 0; entry < tree_numbers->GetEntries(); entry++){
     tree_numbers->GetEntry( entry );

     count_total 	+= total_events;
     count_cleaning	+= survive_cleaning;
     count_BSC 		+= survive_BSC;
     count_HF 		+= survive_HF;
     count_vertex 	+= survive_vertex;
   }

   ofstream cut_efficiency;
   cut_efficiency.open("Cut_efficiency.txt");

   cut_efficiency << std::setprecision(4) << "Total\t" << count_total 		<< "\t" << double(count_total)/double(count_total) * 100. 	<<  " \%" << endl;
   cut_efficiency << std::setprecision(4) << "Cleanign\t" << count_cleaning 	<< "\t" << double(count_cleaning)/double(count_total) * 100. 	<< " \%" <<"\t" << double(count_cleaning)/double(count_total) * 100. 	<<  " \%" << endl;
   cut_efficiency << std::setprecision(4) << "BSC\t" << count_BSC 		<< "\t" << double(count_BSC)/double(count_total) * 100. 	<< " \%" <<"\t" << double(count_BSC)/double(count_cleaning) * 100.  	<<  " \%" << endl;
   cut_efficiency << std::setprecision(4) << "HF\t" << count_HF 		<< "\t" << double(count_HF)/double(count_total) * 100. 		<< " \%" <<"\t" << double(count_HF)/double(count_BSC) * 100.  		<<  " \%" << endl;
   cut_efficiency << std::setprecision(4) << "Vertex\t" << count_vertex 	<< "\t" << double(count_vertex)/double(count_total) * 100. 	<< " \%" <<"\t" << double(count_vertex)/double(count_HF) * 100.  	<<  " \%" << endl;

  cut_efficiency.close();
}





void Unfolder::Plot_Stability(int file_){

  TH2D* response_;
  Get_ResponseMatrix_RooUnfold(file_, response_);
  TH1D* hStability =(TH1D*)response_->ProjectionY();
  hStability->Reset();
  hStability->GetXaxis()->SetTitle("E_{det} [GeV]"  );

  for(int binx = 1; binx <= response_->GetNbinsX(); binx++){
    double Nx = 0.;
    double Ny = 0.;
 
    for(int biny = 1; biny <= response_->GetNbinsY(); biny++){
      Ny += response_->GetBinContent( binx, biny );      
    }

    Nx = response_->GetBinContent( binx, binx );

    hStability->SetBinContent( binx, Nx/Ny );
  }

  TCanvas* can_;
  PrepareCanvas( can_, "Stability" );
  hStability->Draw("hist");

  Finish_canvas( can_ );

  can_->SaveAs(folder_ + "Stability.C");
  can_->SaveAs(folder_ + "Stability.pdf"); 
}


void Unfolder::Plot_Purity(int file_){

  TH2D* response_;
  Get_ResponseMatrix_RooUnfold(file_, response_);
  TH1D* hPurity =(TH1D*)response_->ProjectionY();
  hPurity->Reset();
  hPurity->GetXaxis()->SetTitle("E_{gen} [GeV]");

  for(int biny = 1; biny <= response_->GetNbinsY(); biny++){
    double Nx = 0.;
    double Ny = 0.;
 
    for(int binx = 1; binx <= response_->GetNbinsX(); binx++){
      Nx += response_->GetBinContent( binx, biny );      
    }

    Ny = response_->GetBinContent( biny, biny );

    hPurity->SetBinContent( biny, Ny/Nx );
  }

  TCanvas* can_;
  PrepareCanvas( can_, "Purity" );
  hPurity->Draw("hist");

  Finish_canvas( can_ );

  can_->SaveAs(folder_ + "Purity.C");
  can_->SaveAs(folder_ + "Purity.pdf"); 
}





void Unfolder::Plot_averageEta_perEgen(int file_){

  TH3D* h3D_;
  TFile* _file = TFile::Open( MC_files_[file_], "READ");
  gStyle->SetPalette(1);
  TCanvas* can;   PrepareCanvas_2D(can, "can_3D");
  TCanvas* can2;  PrepareCanvas(can2, "can_2D");

  //== Matched gen. jets only.

  // 3D plot.
  h3D_ = (TH3D*)_file->Get("hEgen_Edet_eta");
  h3D_->Draw("box");
  can->SaveAs(folder_ + "3D_plot_box.C");
  can->SaveAs(folder_ + "3D_plot_box.pdf");

  h3D_->Draw("iso");
  can->SaveAs(folder_ + "3D_plot_iso.C");
  can->SaveAs(folder_ + "3D_plot_iso.pdf");

  // 2D plot: egen - eta.
  TH2D* hXZ = (TH2D*)h3D_->Project3D("xz");
  Prepare_2Dplot( hXZ );
  hXZ->Draw("colz");
  can->SetLogz();
  can->SaveAs(folder_ + "Egen_eta_matched.C");
  can->SaveAs(folder_ + "Egen_eta_matched.pdf");

  // 1D plot: <eta> vs. egen
  TH1D* hZ = (TH1D*)hXZ->ProjectionY();

  for(int bin_egen = 1; bin_egen <= hXZ->GetNbinsY(); bin_egen++){
    double eta_total = 0.;
    double events_total = 0;
    for(int bin_eta = 1; bin_eta <= hXZ->GetNbinsX(); bin_eta++){
      double nevents = hXZ->GetBinContent( bin_eta, bin_egen );
      double eta_value = hXZ->GetXaxis()->GetBinCenter( bin_eta  );

      eta_total += nevents * eta_value;
      events_total += nevents;
    }
    double eta_average = eta_total/events_total;

    hZ->SetBinContent( bin_egen, eta_average );
  }

  can2->cd();
  Prepare_1Dplot( hZ );
  hZ->GetYaxis()->SetRangeUser(-6.6, -5.2);
  hZ->DrawClone("hist");
  can2->SaveAs(folder_ + "Average_eta_per_Egen_matched.C");
  can2->SaveAs(folder_ + "Average_eta_per_Egen_matched.pdf");

  //=================
  //== All gen. jets.
  //=================
  can->cd();
  TH2D* hEgen_eta = (TH2D*)_file->Get("hEgen_eta");
  hEgen_eta->Draw("colz");
  can->SaveAs(folder_ + "Egen_eta_all.C");
  can->SaveAs(folder_ + "Egen_eta_all.pdf");

  for(int bin_egen = 1; bin_egen <= hEgen_eta->GetNbinsX(); bin_egen++){
    double eta_total = 0.;
    double events_total = 0;
    for(int bin_eta = 1; bin_eta <= hEgen_eta->GetNbinsY(); bin_eta++){
      double nevents = hEgen_eta->GetBinContent( bin_egen, bin_eta );
      double eta_value = hEgen_eta->GetYaxis()->GetBinCenter( bin_eta  );

      eta_total += nevents * eta_value;
      events_total += nevents;
    }
    double eta_average = eta_total/events_total;

    hZ->SetBinContent( bin_egen, eta_average );
  }

  can2->cd();
  Prepare_1Dplot( hZ );
  hZ->SetLineColor( kRed );
  hZ->SetLineStyle( 2 );
  hZ->Draw("histsame");
  can2->SaveAs(folder_ + "Average_eta_per_Egen_all.C");
  can2->SaveAs(folder_ + "Average_eta_per_Egen_all.pdf");
}













void Unfolder::Calculate_averageEta_perEgen(int file_, TH1D* &hist_){


  cout << "\n//===//Calculate eta per E" << endl;

  TFile* _file = TFile::Open( MC_files_[file_], "READ");

  //=================
  //== All gen. jets.
  //=================

  TH2D* hEgen_eta = (TH2D*)_file->Get("hEgen_eta");
  TH1D* hZ = (TH1D*)hEgen_eta->ProjectionY();
 

  for(int bin_egen = 1; bin_egen <= hEgen_eta->GetNbinsX(); bin_egen++){
    cout << endl;
    double eta_total = 0.;
    double events_total = 0;
    for(int bin_eta = 1; bin_eta <= hEgen_eta->GetNbinsY(); bin_eta++){
      cout << "(" << hEgen_eta->GetXaxis()->GetBinCenter(bin_egen) << "," << hEgen_eta->GetYaxis()->GetBinCenter(bin_eta) << ")\t" << endl;
      double nevents = hEgen_eta->GetBinContent( bin_egen, bin_eta );
      double eta_value = hEgen_eta->GetYaxis()->GetBinCenter( bin_eta  );

      eta_total += nevents * eta_value;
      events_total += nevents;
    }
    double eta_average = eta_total/events_total;

    hZ->SetBinContent( bin_egen, eta_average );
  }

  hist_ = hZ;

  cout << "\n//===//Calculated eta per E" << endl;
}







void Unfolder::Convert_E_to_pt(TH1D* &hist_){
  cout << "\n//===//Convert to pT" << endl;


  //== The histogram contains <eta> as function of E.
  TH1D* hEta_per_E, *hNewAxis;
  Calculate_averageEta_perEgen( 0, hEta_per_E );

  vector<double> binEdges;

  for(int bin = 1; bin <= hist_->GetNbinsX(); bin++){
  
    double edge_in_E = hist_->GetBinLowEdge( bin );
    double edge_in_pt = edge_in_E / cosh( hEta_per_E->GetBinContent(bin) );

    binEdges.push_back( edge_in_pt );

    cout << "E\t" << edge_in_E << "\tto\t" << edge_in_pt << "\t\teta\t" << hEta_per_E->GetBinContent( bin ) << "\t" <<  cosh(hEta_per_E->GetBinContent( bin )) << endl;
  }
  //== Don't forget upper edge of last bin.
  {
    double edge_in_E = hist_->GetXaxis()->GetBinUpEdge( hist_->GetNbinsX() );
    double edge_in_pt = edge_in_E / cosh(hEta_per_E->GetBinContent( hist_->GetNbinsX()));    
    binEdges.push_back( edge_in_pt );

    cout << "E\t" << edge_in_E << "\tto\t" << edge_in_pt << "\t\teta\t" << hEta_per_E->GetBinContent( hist_->GetNbinsX()) << "\t" <<  cosh(hEta_per_E->GetBinContent( hist_->GetNbinsX())) << endl;

  }

  //== Convert vector to array.
  double* pTbins = &binEdges[0];

  //== Numbers of bins.
  int nBins = binEdges.size()-1;

  //== New histogram.
  hNewAxis = new TH1D(	hist_->GetName(),	hist_->GetTitle(),	nBins, 	pTbins	);

  //== Fill new histogram.
  for(int bin = 1; bin <= hist_->GetNbinsX(); bin++){
    double value = hist_->GetBinContent( bin );
    double error = hist_->GetBinError( bin );
    hNewAxis->SetBinContent( bin, value );
    hNewAxis->SetBinError( bin, error );
  }

  hNewAxis->GetXaxis()->SetTitle("p_{T} [GeV]");
  hist_ = hNewAxis;


  cout << "\n//===//Converted to pT" << endl;
}













void Unfolder::Convert_E_to_xF(TH1D* &hist_){
  cout << "\n//===//Convert to xF" << endl;


  //== The histogram contains <eta> as function of E.
  TH1D* hEta_per_E, *hNewAxis;

  vector<double> binEdges;

  for(int bin = 1; bin <= hist_->GetNbinsX(); bin++){
  
    double edge_in_E = hist_->GetBinLowEdge( bin );

    binEdges.push_back( edge_in_E/eBeam_ );
  }
  //== Don't forget upper edge of last bin.
  {
    double edge_in_E = hist_->GetXaxis()->GetBinUpEdge( hist_->GetNbinsX() );
    binEdges.push_back( edge_in_E/eBeam_ );
  }

  //== Convert vector to array.
  double* xFbins = &binEdges[0];

  //== Numbers of bins.
  int nBins = binEdges.size()-1;

  //== New histogram.
  hNewAxis = new TH1D(	hist_->GetName(),	hist_->GetTitle(),	nBins, 	xFbins	);

  //== Fill new histogram.
  for(int bin = 1; bin <= hist_->GetNbinsX(); bin++){
    double value = hist_->GetBinContent( bin );
    double error = hist_->GetBinError( bin );
    hNewAxis->SetBinContent( bin, value );
    hNewAxis->SetBinError( bin, error );
  }

  hNewAxis->GetXaxis()->SetTitle("x_{F}");
  hist_ = hNewAxis;


  cout << "\n//===//Converted to xF" << endl;
}










/*=================================================================//
	
	1.	Plot the unfolded data
	2.	Calculate and draw the systematic uncertainties
		I.	Model + position
		II.	Model + position + JES
	3.	Draw Gen. level distribution from MC samples
	4.	Draw Gen. level distribution from standalone sample

//=================================================================*/



void Unfolder::Plot_Unfolded_Ratio_allSystematics_pT(TCanvas* can_, TString variable, int iterations, TString plot_as, TString which_models){

  //== Write down systematics.
  ofstream systematics_txt;
  systematics_txt.open("systematics.txt");

  //== Write down x-sections.
  ofstream xsec_txt;
  xsec_txt.open( TString::Format( "xsec_" + which_models  + ".txt" ) );
  
  //== Set highest energy shown.
  Emax_ = 1500.;

  //== Prepare plots to use for correlation between systematics.
  TFile *_fSystematics = TFile::Open("_fSystematics_7TeV.root", "Update");
  
  TH1D *hReference_benoit,
	*hReference_benoit_dndx,
	*hLumi_up_benoit_dndx,
	*hLumi_down_benoit_dndx,
	*hModel_up_benoit,
	*hModel_up_benoit_dndx,
	*hModel_down_benoit,
	*hModel_down_benoit_dndx,
	*hPos_up_benoit,
	*hPos_up_benoit_dndx,
	*hPos_down_benoit,
	*hPos_down_benoit_dndx,
	*hJES_up_benoit,
	*hJES_up_benoit_dndx,
	*hJES_down_benoit,
	*hJES_down_benoit_dndx;

  //== Prepare canvas and pads.
  TPad* pad_abs_, *pad_ratio_;
  PrepareCanvas(can_, "Systematics");
  SplitCanvas(can_, pad_abs_, pad_ratio_);

  //== Prepare histograms.
  TH1D* hAverage, *hFirst;
  TString drawoptions = "phist";
  TH1D* hGen, *hGen_ratio, *hDet, *hReference, *hMC, *hMC_ratio;

  //== Reference histogram.
  double E_to_xf_;
  
  cout << "Before we start: " << set_xf_ << endl;
  
  if( plot_as == "E"){ 		E_to_xf_ = 1.; 		Set_to_xf(false); }
  else if( plot_as == "xf") { 	E_to_xf_ = 3500.; 	Set_to_xf(true);}
  else{ return; }  
  
  cout << "After we start: " << set_xf_ << endl;  
  
  Get_DetEnergy( -1 , hReference); 
  
  cout << "\%\%\%\\ Got det energy" << endl;
  if( !hReference ){ cout << "Not" << endl; }
  cout << "\%\%\%\\ " << hReference->GetNbinsX() << "\t" << hReference->Integral() << endl;

  //-- Normalize hReference (DATA) to xsec by dividing by its luminosity.
  hReference->Scale( 1./ lumi_ );  
  hReference->GetYaxis()->SetTitle( "d#sigma/dE [mb/GeV]" );
  
  //== Unfold the data.
  Get_DetUnfolded( 0 , hReference, iterations, variable);	// Unfold hReference with the Response from MC sample 0.

  //== dN/dx
  SetDnDx( hReference );
  hReference->SetName("UnfoldedData");

  //== Benoit.
  {
    _fSystematics->cd();
        
    hReference_benoit_dndx = (TH1D*)hReference->Clone("Data_dsigmadE");
    hLumi_up_benoit_dndx = (TH1D*)hReference_benoit_dndx->Clone("hLumi_up");
    hLumi_up_benoit_dndx->Scale( 1.036);
    //hLumi_up_benoit_dndx->Divide(hReference_benoit_dndx);
    hLumi_up_benoit_dndx->Write();
    
    hLumi_down_benoit_dndx = (TH1D*)hReference_benoit_dndx->Clone("hLumi_down");
    hLumi_down_benoit_dndx->Scale( 0.964);
    //hLumi_down_benoit_dndx->Divide(hReference_benoit_dndx);
    hLumi_down_benoit_dndx->Write();
    
    hReference_benoit_dndx->Write();
    
    _fSystematics->Close();
    
  }
  
  

  //== Determine the first bin above 300 GeV.
  int bin_300 = 0;
  double ElowBin = 300.;
  if( plot_as == "xf" ){ ElowBin = 280./3500.;   }
  while( hReference->GetXaxis()->GetBinLowEdge( bin_300 ) < ElowBin ){ bin_300++; }
  bin_300--;
  double Eplot_lowest = hReference->GetXaxis()->GetBinLowEdge( bin_300 );
  double Emax = hReference->GetXaxis()->GetBinLowEdge( hReference->GetNbinsX() );

  //== The last bin of the histogram is not the last bin of the plot: cut it away.
  SetSubhistogram_max( Emax );
  
  //== Get the integral of the area plotted.
  TH1D* hReference_copy = (TH1D*)hReference->Clone("copy");

  //== Legend (1): for models.
  double legx, legy;
  if( which_models == "Pythia"){ legx = 1. - pad_abs_->GetRightMargin()-0.47; legy = 1. - pad_abs_->GetTopMargin() - 0.45; }
  else{ legx = 1. - pad_abs_->GetRightMargin()-0.35; legy = 1. - pad_abs_->GetTopMargin() - 0.5; }

  TLegend* legend = new TLegend( 
	legx, 	1. - pad_abs_->GetTopMargin() - 0.45, 
	1. - pad_abs_->GetRightMargin(), 	1. - pad_abs_->GetTopMargin() );

  //== Legend (2):  for systematics.
  TLegend* legend_syst	 = new TLegend( 
	pad_abs_->GetLeftMargin(), 	pad_abs_->GetBottomMargin(), 
	pad_abs_->GetLeftMargin()+0.35, 	pad_abs_->GetBottomMargin() +0.5);
  legend_syst->AddEntry( hReference, TString::Format("Data"), "lp");

  //== Luminosity uncertainty: fixed 3.6%
  double lumi_dep = 0.036;

  /******************************************************************************************
  * Begin by averaging the models and taking the difference between average and each model. *
  ******************************************************************************************/

  //== Loop over MC files.

  TH1D* hReference_ratio = (TH1D*)hReference->Clone("hReference_ratio");
  hReference_ratio->Divide( hReference );

  /*************************************************************************
  The following is needed to calculate the proper position of the datapoints.
  See Eq. (7) in  Physics Research A 355 (1995) 541-547.

  We have a spectrum with two distinct slopes: the first slope is called low, the second is called high.
  **************************************************************************/
 
  cout << "===WHERE TO STICK YOUR DATA POINTS===" << endl;

  double alow = 185509.8;
  double blow = 0.007383852;   

  double ahigh = 24325.05;
  double bhigh = 0.004686864;

  double a = alow, b = blow;

  const Int_t n = hReference_ratio->GetNbinsX();
  // x-axis and the values.
  Double_t x[2*n], y[2*n];

  // Systematic uncertainty lumi.
  Double_t ex_lumi[2*n], ey_lumi[2*n];
  Calculate_lumi_systematic( iterations,  ex_lumi, ey_lumi);
  
  // Systematic uncertainty model.
  Double_t ex_model[2*n], ey_model[2*n];  
  Calculate_model_systematic(  iterations, ex_model, ey_model);
  
  // Systematic uncertainty position.
  Double_t ex_position_up[2*n], ey_position_up[2*n]; 
  Double_t ex_position_down[2*n], ey_position_down[2*n];  
  Calculate_position_systematic(iterations,  ex_position_up, ey_position_up, ex_position_down, ey_position_down );   
  
  // Systematic uncertainty jetID.
  Double_t ex_id[2*n], ey_id[2*n]; 
  Calculate_jetID_systematic( iterations, ex_id, ey_id);
  
  // Systematic uncertainty JES.
  Double_t ex_jes_up[2*n], ey_jes_up[2*n]; 
  Double_t ex_jes_down[2*n], ey_jes_down[2*n];   
  Calculate_JES_systematic(  iterations, ex_jes_up, ey_jes_up, ex_jes_down, ey_jes_down );  
  
  // SUM: Systematic uncertainty lumi or model.
  Double_t exl_LM[2*n], eyl_LM[2*n], exh_LM[2*n], eyh_LM[2*n];

  // SUM: Systematic uncertainty lumi, model, pos.
  Double_t exl_LMI[2*n], eyl_LMI[2*n], exh_LMI[2*n], eyh_LMI[2*n];
  
  // SUM: Systematic uncertainty lumi, model, pos, ID
  Double_t exl_LMIP[2*n], eyl_LMIP[2*n], exh_LMIP[2*n], eyh_LMIP[2*n];
  
  // SUM: Systematic uncertainty lumi, model, pos, ID
  Double_t exl_LMIPJ[2*n], eyl_LMIPJ[2*n], exh_LMIPJ[2*n], eyh_LMIPJ[2*n];   
 
  double x_data[n], y_data[n], y_ratio[n];
  double exl_data[n],exh_data[n];
  double eyl_data[n],eyh_data[n];
  double eyl_ratio[n],eyh_ratio[n];

  // Remember: hist bins start at 1, graph bins at 0.
  for(int bin = 0; bin < hReference->GetNbinsX(); bin++){
    double bin_lowedge = hReference->GetXaxis()->GetBinLowEdge( bin+1 );
    double bin_highedge = hReference->GetXaxis()->GetBinUpEdge( bin+1 );

    if( bin_lowedge > 900. ){ a = ahigh; b = bhigh; }

    double binwidth = hReference->GetXaxis()->GetBinWidth( bin+1);

    double xlw = bin_lowedge 
		+ 1./b* (log(b * binwidth) )
		- 1./b* (log( 1. - exp( -1. * b * binwidth) ));

    int bin_shift = 1;
    if (bin == n-1){ bin_shift = 0; }	//== This ensures the last ratio is drawn neatly.

    x_data[bin] = xlw; //E_to_xf_;
    y_data[bin] = hReference->GetBinContent( bin+bin_shift ); // * E_to_xf_;
    y_ratio[bin] = 1.;
    exl_data[bin] = (xlw - bin_lowedge); ///E_to_xf_;
    exh_data[bin] = (bin_highedge - xlw); ///E_to_xf_;
    eyl_data[bin] = hReference->GetBinError( bin+bin_shift );
    eyh_data[bin] = hReference->GetBinError( bin+bin_shift );
    eyl_ratio[bin] = hReference->GetBinError( bin+bin_shift )/hReference->GetBinContent( bin+bin_shift  );
    eyh_ratio[bin] = hReference->GetBinError( bin+bin_shift )/hReference->GetBinContent( bin+bin_shift  );
  }

   TGraphAsymmErrors *gr_data = new TGraphAsymmErrors(n,x_data,y_data,exl_data,exh_data,eyl_data,eyh_data);
     gr_data->SetLineWidth( 1 );
     gr_data->SetMarkerSize( 1.25 );
     gr_data->SetMarkerStyle( 20 );
     gr_data->SetTitle("Graph_data");
     gr_data->SetName("Graph_data");


   TGraphAsymmErrors *gr_ratio = new TGraphAsymmErrors(n,x_data,y_ratio,exl_data,exh_data,eyl_ratio,eyh_ratio);
     gr_ratio->SetLineWidth( 1 );
     gr_ratio->SetMarkerSize( 1.25 );
     gr_ratio->SetMarkerStyle( 20 );
     gr_ratio->SetTitle("Graph_ratio");
     gr_ratio->SetName("Graph_ratio");


  for(int bin = 0; bin < n; bin++){
    int bin_shift = 1;
    //if (bin == n-1){ bin_shift = 0; }	//== This ensures the last ratio is drawn neatly.
    //if (bin == 2){ bin_shift = 2; }	//== This ensures the last ratio is drawn neatly.

    x[2*bin] = hReference->GetXaxis()->GetBinLowEdge( bin+1 ); ///E_to_xf_;
    y[2*bin] = hReference_ratio->GetBinContent( bin+1 );
    x[2*bin+1] = hReference->GetXaxis()->GetBinUpEdge( bin+1 ); ///E_to_xf_;
    y[2*bin+1] = hReference_ratio->GetBinContent( bin+1);

    cout << "Bin\t" << bin << "\tfrom\t" << x[2*bin] << "\tto\t" << x[2*bin+1] << endl;

    //== LUMI + MODEL.
    eyl_LM[2*bin] = sqrt( ey_lumi[2*bin]*ey_lumi[2*bin] + ey_model[2*bin]*ey_model[2*bin]  );
    eyh_LM[2*bin] = sqrt( ey_lumi[2*bin]*ey_lumi[2*bin] + ey_model[2*bin]*ey_model[2*bin]  );
    eyl_LM[2*bin+1] = sqrt( ey_lumi[2*bin+1]*ey_lumi[2*bin+1] + ey_model[2*bin+1]*ey_model[2*bin+1]  );
    eyh_LM[2*bin+1] = sqrt( ey_lumi[2*bin+1]*ey_lumi[2*bin+1] + ey_model[2*bin+1]*ey_model[2*bin+1]  );
    
    //== LUMI + MODEL + JetID.
    eyl_LMI[2*bin] = sqrt( eyl_LM[2*bin]*eyl_LM[2*bin] + ey_id[2*bin]*ey_id[2*bin]  );
    eyh_LMI[2*bin] = sqrt( eyh_LM[2*bin]*eyh_LM[2*bin] + ey_id[2*bin]*ey_id[2*bin]  );   
    eyl_LMI[2*bin+1] = sqrt( eyl_LM[2*bin+1]*eyl_LM[2*bin+1] + ey_id[2*bin+1]*ey_id[2*bin+1]  );
    eyh_LMI[2*bin+1] = sqrt( eyh_LM[2*bin+1]*eyh_LM[2*bin+1] + ey_id[2*bin+1]*ey_id[2*bin+1]  );    
    
    //== LUMI + MODEL + JetID + POSITION.
    eyl_LMIP[2*bin] = sqrt( eyl_LMI[2*bin]*eyl_LMI[2*bin] + ey_position_down[2*bin]*ey_position_down[2*bin]  );
    eyh_LMIP[2*bin] = sqrt( eyh_LMI[2*bin]*eyh_LMI[2*bin] + ey_position_up[2*bin]*ey_position_up[2*bin]  );    
    eyl_LMIP[2*bin+1] = sqrt( eyl_LMI[2*bin+1]*eyl_LMI[2*bin+1] + ey_position_down[2*bin+1]*ey_position_down[2*bin+1]  );
    eyh_LMIP[2*bin+1] = sqrt( eyh_LMI[2*bin+1]*eyh_LMI[2*bin+1] + ey_position_up[2*bin+1]*ey_position_up[2*bin+1]  );   
    
    //== LUMI + POSITION + MODEL + JetID + JES.
    eyl_LMIPJ[2*bin] = sqrt( eyl_LMIP[2*bin]*eyl_LMIP[2*bin] + ey_jes_down[2*bin]*ey_jes_down[2*bin]  );
    eyh_LMIPJ[2*bin] = sqrt( eyh_LMIP[2*bin]*eyh_LMIP[2*bin] + ey_jes_up[2*bin]*ey_jes_up[2*bin]  ); 
    exl_LMIPJ[2*bin] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;
    exh_LMIPJ[2*bin] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;
        
    eyl_LMIPJ[2*bin+1] = sqrt( eyl_LMIP[2*bin+1]*eyl_LMIP[2*bin+1] + ey_jes_down[2*bin+1]*ey_jes_down[2*bin+1]  );
    eyh_LMIPJ[2*bin+1] = sqrt( eyh_LMIP[2*bin+1]*eyh_LMIP[2*bin+1] + ey_jes_up[2*bin+1]*ey_jes_up[2*bin+1]  );
    exl_LMIPJ[2*bin+1] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;
    exh_LMIPJ[2*bin+1] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;      
  }
  
  pad_ratio_->cd();
  TH1D* hReference__ = (TH1D*)hReference->Clone("hNeededToProperlySetAxes");
  hReference__->Divide( hReference__ );
  if( plot_as == "E" ){ hReference__->GetXaxis()->SetRangeUser(Eplot_lowest, 1500.); hReference__->GetXaxis()->SetTitle("E [GeV]"); }
  if( plot_as == "xf" ){ hReference__->GetXaxis()->SetRangeUser(0.1, 0.38); hReference__->GetXaxis()->SetTitle("x_{F}"); }
  Prepare_1Dplot_ratio( hReference__ );
  hReference__->GetYaxis()->SetRangeUser(0.,2.75);
  hReference__->GetYaxis()->SetNdivisions( 504 );  
  hReference__->GetYaxis()->SetTitle("Ratio");
  hReference__->Draw("hist");

  TGraphAsymmErrors *gr_all_ratio, *gr_id_ratio, *gr_mod_ratio, *gr_pos_ratio, *gr_lumi_ratio;


  //== All uncertainty.
  TGraphAsymmErrors *gr = new TGraphAsymmErrors(2*n,x,y,exl_data,exh_data,eyl_LMIPJ,eyh_LMIPJ);
  gr->SetTitle("TGraphAsymmErrors All");
  int ci_all = TColor::GetColor("#FFFF00");
  gr->SetMarkerSize(0);
  gr->SetFillColor( ci_all );
  gr->SetMarkerStyle(21);

  gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
  if( plot_as == "xf" ){    gr->GetHistogram()->GetXaxis()->SetTitle("x_{F}"); }
  if( plot_as == "pT" ){    gr->GetHistogram()->GetXaxis()->SetTitle("p_{T} [GeV]"); }
  gr->GetHistogram()->GetYaxis()->SetTitle("Ratio");
  gr->GetHistogram()->GetYaxis()->SetNdivisions(205);
  gr->GetHistogram()->GetYaxis()->CenterTitle();
  gr->GetHistogram()->GetYaxis()->SetRangeUser(0, 3.35);  
  Prepare_1Dplot_ratio( gr );
  gr->SetTitle("all_ratio");
  gr->SetName("all_ratio");
  gr_all_ratio = (TGraphAsymmErrors*)gr->Clone("all_ratio");
  gr_all_ratio->Draw("E3same");
  
  pad_ratio_->RedrawAxis();
  pad_ratio_->Update();  
  
  //== Lumi + position + models + jetID
  gr = new TGraphAsymmErrors(2*n,x,y,exl_data,exh_data,eyl_LMIP,eyh_LMIP);
  int ci_id = TColor::GetColor("#000099");  
  gr->SetFillColor( ci_id );
  gr->GetHistogram()->GetYaxis()->SetRangeUser( 0.1*GetMinimumValue( hReference ), hReference->GetMaximum() * 1.1 );
  gr->GetHistogram()->GetYaxis()->SetTitle("#frac{dN}{dE}");
  gr_pos_ratio = (TGraphAsymmErrors*)gr->Clone("pos_ratio");
  gr_pos_ratio->Draw("E3same");

  //== Model + position + luminosity uncertainty. 
  gr = new TGraphAsymmErrors(2*n,x,y,ex_lumi,ex_lumi,eyl_LMI,eyh_LMI);
  gr->SetTitle("TGraphAsymmErrors No JES");
  int ci = TColor::GetColor("#C0C0C0");
  gr->SetMarkerSize(0);
  gr->SetFillColor( ci );
  gr->SetMarkerStyle(21);
  gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
  gr->GetHistogram()->GetYaxis()->SetRangeUser(0., 2.);
  gr_id_ratio = (TGraphAsymmErrors*)gr->Clone("JetID_ratio");
  gr_id_ratio->Draw("PE3same");

  //== Position + luminosity uncertainty.
  gr = new TGraphAsymmErrors(2*n,x,y,ex_lumi,ex_lumi,eyl_LM,eyh_LM);
  gr->SetTitle("TGraphAsymmErrors No JES, No model");
  int ci_noModel = TColor::GetColor("#33FF99");
  gr->SetMarkerSize(0);
  gr->SetFillColor( ci_noModel );
  gr->SetMarkerStyle(21);
  gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
  gr->GetHistogram()->GetYaxis()->SetRangeUser(0., 2.);
  gr_mod_ratio = (TGraphAsymmErrors*)gr->Clone("Model_ratio");
  gr_mod_ratio->Draw("PE3same");

  //== Luminosity uncertainty.
  gr = new TGraphAsymmErrors(2*n,x,y,ex_lumi,ex_lumi,ey_lumi,ey_lumi);
  gr->SetTitle("TGraphAsymmErrors Lumi");
  int ci_lumi =   TColor::GetColor("#FFCCCC");
  gr->SetMarkerSize(0);
  gr->SetFillColor( ci_lumi );
  gr->SetMarkerStyle(21);
  gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
  gr->GetHistogram()->GetYaxis()->SetRangeUser(0., 2.);
  gr_lumi_ratio = (TGraphAsymmErrors*)gr->Clone("Lumi_ratio");
  gr_lumi_ratio->Draw("PE3same");
  hReference_ratio->SetMarkerColor( kBlack );
  gr_ratio->Draw("psame");

  //---------------------//
  //-- Absolute value. --//
  //---------------------//

  cout << "===ABSOLUTE===" << endl;
  pad_abs_->cd();

  for(int bin = 0; bin < n; bin++){

    int bin_shift = 1;
    //if (bin == n-1){ bin_shift = 0; }	//== This ensures the last ratio is drawn neatly.
    //if (bin == 2){ bin_shift = 2; }	//== This ensures the last ratio is drawn neatly.

    x[2*bin] = hReference->GetXaxis()->GetBinLowEdge( bin+1 ); ///E_to_xf_;
    y[2*bin] = hReference->GetBinContent( bin+bin_shift ); // * E_to_xf_;
    x[2*bin+1] = hReference->GetXaxis()->GetBinUpEdge( bin+1 ); ///E_to_xf_;
    y[2*bin+1] = hReference->GetBinContent( bin+bin_shift ); //*E_to_xf_;

    //-- Error band lumi.
    ey_lumi[2*bin] = 0.036 * y[2*bin] ;
    ex_lumi[2*bin] = hReference_ratio->GetBinContent( bin+1 ); ///E_to_xf_;

    ey_lumi[2*bin+1] = 0.036 * y[2*bin] ;
    ex_lumi[2*bin+1] = hReference_ratio->GetBinWidth( bin+1 )/2; ///E_to_xf_;

    //-- Error band lumi + position.
    eyl_LM[2*bin] = hReference->GetBinContent( bin+bin_shift ) *  eyl_LM[2*bin]; // * E_to_xf_;
    eyh_LM[2*bin] = hReference->GetBinContent( bin+bin_shift ) *  eyh_LM[2*bin]; // * E_to_xf_;
    exl_LM[2*bin] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;
    exh_LM[2*bin] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;

    eyl_LM[2*bin+1] = hReference->GetBinContent( bin+bin_shift ) *  eyl_LM[2*bin+1]; // * E_to_xf_;
    eyh_LM[2*bin+1] = hReference->GetBinContent( bin+bin_shift ) *  eyh_LM[2*bin+1]; // * E_to_xf_;
    exl_LM[2*bin+1] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;
    exh_LM[2*bin+1] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;

    //-- Error band lumi + position + model.
    eyl_LMI[2*bin] = hReference->GetBinContent( bin+bin_shift ) *  eyl_LMI[2*bin] ; //* E_to_xf_;
    eyh_LMI[2*bin] = hReference->GetBinContent( bin+bin_shift ) *  eyh_LMI[2*bin] ; //* E_to_xf_;
    exl_LMI[2*bin] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;
    exh_LMI[2*bin] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;

    eyl_LMI[2*bin+1] = hReference->GetBinContent( bin+bin_shift ) *  eyl_LMI[2*bin+1]; // * E_to_xf_;
    eyh_LMI[2*bin+1] = hReference->GetBinContent( bin+bin_shift ) *  eyh_LMI[2*bin+1]; // * E_to_xf_;
    exl_LMI[2*bin+1] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;
    exh_LMI[2*bin+1] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;
    
    //-- Error band lumi + position + model + jetID. 

    eyl_LMIP[2*bin] = hReference->GetBinContent( bin+bin_shift ) *  eyl_LMIP[2*bin]; // * E_to_xf_;
    eyh_LMIP[2*bin] = hReference->GetBinContent( bin+bin_shift ) *  eyh_LMIP[2*bin]; // * E_to_xf_;
    exl_LMIP[2*bin] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;
    exh_LMIP[2*bin] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;     

    eyl_LMIP[2*bin+1] = hReference->GetBinContent( bin+bin_shift ) *  eyl_LMIP[2*bin+1]; // * E_to_xf_;   
    eyh_LMIP[2*bin+1] = hReference->GetBinContent( bin+bin_shift ) *  eyh_LMIP[2*bin+1]; // * E_to_xf_;
    exl_LMIP[2*bin+1] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;
    exh_LMIP[2*bin+1] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;  

    //-- Error band with JES.
    eyl_LMIPJ[2*bin] = hReference->GetBinContent( bin+bin_shift ) * eyl_LMIPJ[2*bin]; // * E_to_xf_;
    eyh_LMIPJ[2*bin] = hReference->GetBinContent( bin+bin_shift ) * eyh_LMIPJ[2*bin]; // * E_to_xf_;
    exl_LMIPJ[2*bin] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;
    exh_LMIPJ[2*bin] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;

    eyl_LMIPJ[2*bin+1] = hReference->GetBinContent( bin+bin_shift ) * eyl_LMIPJ[2*bin+1]; // * E_to_xf_;
    eyh_LMIPJ[2*bin+1] = hReference->GetBinContent( bin+bin_shift ) * eyh_LMIPJ[2*bin+1]; // * E_to_xf_;
    exl_LMIPJ[2*bin+1] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;
    exh_LMIPJ[2*bin+1] = hReference->GetBinWidth( bin+bin_shift )/2; ///E_to_xf_;
  }
   
  pad_abs_->SetLogy();
  Prepare_1Dplot( hReference );
  if( plot_as == "E"){ hReference->GetYaxis()->SetTitle("d#sigma/dE [mb/GeV]"); }
  if( plot_as == "xf"){ hReference->GetYaxis()->SetTitle("d#sigma/dx_{f} [mb/x_{f}]"); }  
  hReference->Draw("hist");
  if( plot_as == "E" ){ hReference->GetXaxis()->SetRangeUser(Eplot_lowest, 1500.); }
  if( plot_as == "xf" ){ 
    hReference->GetXaxis()->SetRangeUser(0.1, 0.38);
    hReference->GetXaxis()->SetTitle("x_{f}"); 
  }
  
  TFile* _file = TFile::Open("Alex_7TeV.root", "Update");
  TH1D* hReference_written = (TH1D*)hReference->Clone( TString::Format("UnfoldedData_" + plot_as) );
  hReference_written->Write();
  _file->Close();

  TGraphAsymmErrors *gr_all, *gr_id, *gr_mod, *gr_pos, *gr_lumi;

  //-- Draw graph with smaller errors on top of graph with bigger errors.
  gr = new TGraphAsymmErrors(2*n,x,y,exl_data,exh_data,eyl_LMIPJ,eyh_LMIPJ);
  gr->SetFillColor(ci_all);
  gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
 
  double ylow =   0.5*GetMinimumValue( hReference ) ; // *E_to_xf_ ;
  double yhigh = hReference->GetMaximum() * 1.1; // *E_to_xf_ ;

  ylow = 2*1e-5*E_to_xf_ ;
  yhigh = 0.5*E_to_xf_ ;
  gr->GetHistogram()->GetYaxis()->SetRangeUser(ylow, yhigh );

  if( plot_as == "E"){ gr->GetHistogram()->GetYaxis()->SetTitle("d#sigma/dE [mb/GeV]"); }
  if( plot_as == "xf"){ gr->GetHistogram()->GetYaxis()->SetTitle("d#sigma/dx_{F} [mb]"); }


  //== Only Y-axis needed.
  Prepare_1Dplot( gr );
  gr->Draw("E3same");
  gr_all = (TGraphAsymmErrors*)gr->Clone("All_absolute");

  //legend_syst->AddEntry( gr, "Syst. errors (All)", "f");
  
  //== Lumi + position + models + jetID
  gr = new TGraphAsymmErrors(2*n,x,y,exl_LMIP,exh_LMIP,eyl_LMIP,eyh_LMIP);
  gr->SetFillColor( ci_id );
  if( plot_as == "E"){ 
    gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
  }
  gr->GetHistogram()->GetYaxis()->SetRangeUser( 0.1*GetMinimumValue( hReference )*E_to_xf_, hReference->GetMaximum() * 1.1*E_to_xf_ );
  gr->GetHistogram()->GetYaxis()->SetTitle("#frac{dN}{dE}");
  gr->Draw("E3same");
  gr_id = (TGraphAsymmErrors*)gr->Clone("Position_absolute"); 

  //== Lumi + position + models.
  gr = new TGraphAsymmErrors(2*n,x,y,exl_LMI,exh_LMI,eyl_LMI,eyh_LMI);
  gr->SetFillColor(ci);
  if( plot_as == "E"){ 
    gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
  }
  gr->GetHistogram()->GetYaxis()->SetRangeUser( 0.1*GetMinimumValue( hReference ), hReference->GetMaximum() * 1.1 );
  gr->GetHistogram()->GetYaxis()->SetTitle("#frac{dN}{dE}");
  gr->Draw("E3same");
  gr_mod = (TGraphAsymmErrors*)gr->Clone("JetID_absolute");

  //== Lumi + position.
  gr = new TGraphAsymmErrors(2*n,x,y,exl_LM,exh_LM,eyl_LM,eyh_LM);
  gr->SetFillColor(ci_noModel);
  if( plot_as == "E" ){
    gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
  }
  gr->GetHistogram()->GetYaxis()->SetRangeUser( 0.1*GetMinimumValue( hReference ), hReference->GetMaximum() * 1.1 );
  gr->GetHistogram()->GetYaxis()->SetTitle("#frac{dN}{dE}");
  gr->Draw("E3same");
  gr_pos = (TGraphAsymmErrors*)gr->Clone("Model_absolute");

  //== Lumi.
  gr = new TGraphAsymmErrors(2*n,x,y,exl_data,exh_data,ey_lumi,ey_lumi);
  gr->SetFillColor(ci_lumi);
  if( plot_as == "E" ){
    gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
  }
  gr->GetHistogram()->GetYaxis()->SetRangeUser( 0.1*GetMinimumValue( hReference ), hReference->GetMaximum() * 1.1 );
  gr->GetHistogram()->GetYaxis()->SetTitle("#frac{dN}{dE}");
  gr->Draw("E3same");
  gr_lumi = (TGraphAsymmErrors*)gr->Clone("Lumi_absolute");

  hReference->SetMarkerColor( kBlack );

  gr_data->Draw("psame");

  legend_syst->AddEntry( gr_lumi, "Syst. errors (lumi.)", "f");
  legend_syst->AddEntry( gr_pos, "Syst. errors (+ mod.)", "f");
  legend_syst->AddEntry( gr_mod, "Syst. errors (+ id)", "f");  
  legend_syst->AddEntry( gr_id, "Syst. errors (+ pos.)", "f");   
  legend_syst->AddEntry( gr_all, "Syst. errors (+ JES)", "f");
  
  //-- Scale and draw generator level distribution.
  int color_ = 0;
  int line_ = 1;
  drawoptions = "][histsame";
  ofstream binwidths;
  binwidths.open("binwidths.txt");

  TString split_MCs = "";
  int count_MCs = 0;
  
  cout << "===MCs===" << endl;
  
  for(int MC_ = 0; MC_ < MC_files_.size(); MC_++){
  
    if( plot_as == "pT" ){ continue; }
    //-- Skip colors if needed.
    color_++ ;
    bool goodcolor = false;
    while( !goodcolor ){
      if( color_ == 1 || color_ == 8 || color_ == 9 || color_ == 10 || color_ == 11) color_++;
      else{ goodcolor = true; }
    }

    TH1D* hGen;
    //-- Only continue if the MC is one of the following (excludes position samples).
    if( ((set_of_tags_["mc_type"])[MC_files_[MC_]] != "model" && 
	(set_of_tags_["mc_type"])[MC_files_[MC_]] != "actual" &&  
	(set_of_tags_["mc_type"])[MC_files_[MC_]] != "shift_MPI_or_Tune" &&
	(set_of_tags_["mc_type"])[MC_files_[MC_]] != "model_") ){ continue; }
	
    cout << "Continue with file\t" << MC_files_[MC_] << "\t" << which_models << endl;

    //-- Split up the canvas.
    //== "Pythia" will plot only models from Pythia.
    //== "noPythia" will plot models without Pythia.
    //== "All" plots all models.
    //== Any other entry skips the addition of models.
    
    if( which_models == "Pythia" ){
      if( !MC_files_[MC_].Contains( "ythia") ){ continue; }
    }
    else if( which_models == "noPythia"){
      cout << "No Pythia" << endl;
      if( MC_files_[MC_].Contains( "ythia") ){ continue; }
    }
    else if( which_models != "All"){ break; }

    cout << "Continue with model" << endl;

    //-- Check for xsec = 0 or NaN.
    if( xsec_[ MC_files_[MC_] ] !=  xsec_[ MC_files_[MC_] ] ){ continue; }
    if( xsec_[ MC_files_[MC_] ] ==  0. ){ continue; }
   
    //-- If scaletodata_ is true, the gen. energy spectrum will be scaled with numbers taken from the file before being passed back.
    scaletodata_ = true;

    // There is only one gen. distribution, since there are no cuts on detector level.
    if( (set_of_tags_["mc_type"])[MC_files_[MC_]] == "shift_MPI_or_Tune" ){	Get_GenEnergy( MC_, hGen ); }
    else{ Get_PureGenEnergy( MC_, hGen );  }
    scaletodata_ = false;
    if( !hGen ){ continue; }
    if( hGen->Integral() != hGen->Integral() ){ continue; }

    hGen->SetTitle( TString::Format("hGenJet_energy_%i", MC_) );
    hGen->SetName( TString::Format("hGenJet_energy_%i", MC_) );

    SetDnDx( hGen );

    if( (set_of_tags_["mc_type"])[MC_files_[MC_]] == "shift_MPI_or_Tune" ){
      drawoptions = "][histsame";
    }

    //-- Set line properties.
    hGen->SetLineColor( getColor(color_) );
    hGen->SetLineStyle( (line_ != 1)*((line_++)%7+2) );
    hGen->SetLineWidth( 3 );
    
    pad_abs_->cd();
    //hGen->GetXaxis()->SetRangeUser(Eplot_lowest/E_to_xf_, Emax_/E_to_xf_); 
    hGen->DrawClone( drawoptions);    

    //-- Set ratio of distribution.
    pad_ratio_->cd();    
    hGen_ratio = (TH1D*)hGen->Clone( TString::Format("hGen_ratio_%i", MC_) );
    hGen_ratio->Divide( hReference );
    hGen_ratio->Draw( drawoptions);

    TString legend_file = legend_info_gen_[MC_files_[MC_]];
    
    if( (set_of_tags_["mc_type"])[MC_files_[MC_]] == "shift_MPI_or_Tune" ){
      legend->AddEntry( hGen, TString::Format( legend_file  ), "l");
    }
    else{
      legend->AddEntry( hGen, legend_file , "l");
    }
    count_MCs++;
  }

   pad_abs_->cd();

   //== Adjust height of legend.
   legend->SetY1( 1. - pad_abs_->GetTopMargin() - .09 * count_MCs   );
   legend->SetX2( 1. - pad_abs_->GetRightMargin() - gStyle->GetTickLength()   );
   legend->SetFillStyle(0);
   legend->SetLineWidth(0);
   legend->Draw();

   pad_abs_->cd();
   pad_abs_->RedrawAxis();
   pad_abs_->Update();
   legend_syst->SetHeader( TString::Format("%i iterations", iterations) );
   legend_syst->SetFillStyle(0);
   legend_syst->SetLineWidth(0);
   legend_syst->SetTextSize( legend->GetTextSize() );   
   legend_syst->Draw();

   pad_ratio_->RedrawAxis();
   pad_ratio_->Update();

   can_->cd();

   //== Text.
   Finish_canvas( pad_abs_ , "leftish", etawidth_);

   can_->SaveAs( TString::Format( folder_ + "Totaldependence_" + plot_as + "_" + which_models + "_%iit_deltaphi_0%i_etaband_0%i" + cutsetup_ + ".C", iterations, static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_)) );
   can_->SaveAs( TString::Format( folder_ + "Totaldependence_" + plot_as + "_" + which_models + "_%iit_deltaphi_0%i_etaband_0%i" + cutsetup_ + ".pdf", iterations, static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_)) ); 

   systematics_txt.close();
   xsec_txt.close();
   
   //PlotSystematics();
   Set_to_xf(false); 
}
















void Unfolder::Plot_Unfolded_Ratio_allSystematics_separate(TCanvas* can_, TString variable, int iterations, TString plot_as, TString which_models){

  //== Write down systematics.
  ofstream systematics_txt;
  systematics_txt.open("systematics.txt");

  //== Write down x-sections.
  ofstream xsec_txt;
  xsec_txt.open( TString::Format( "xsec_" + which_models  + ".txt" ) );

  //== Prepare plots to use for correlation between systematics.
  TH1D *hReference_benoit,
	*hReference_benoit_dndx,
	*hLumi_up_benoit_dndx,
	*hLumi_down_benoit_dndx,
	*hModel_up_benoit,
	*hModel_up_benoit_dndx,
	*hModel_down_benoit,
	*hModel_down_benoit_dndx,
	*hPos_up_benoit,
	*hPos_up_benoit_dndx,
	*hPos_down_benoit,
	*hPos_down_benoit_dndx,
	*hJES_up_benoit,
	*hJES_up_benoit_dndx,
	*hJES_down_benoit,
	*hJES_down_benoit_dndx;

  //== Prepare canvas and pads.
  TPad* pad_abs_, *pad_ratio_;
  PrepareCanvas(can_, "Systematics");
  SplitCanvas(can_, pad_abs_, pad_ratio_);

  //== Prepare histograms.
  TH1D* hAverage, *hFirst;
  TString drawoptions = "phist";
  TH1D* hGen, *hGen_ratio, *hDet, *hReference, *hMC, *hMC_ratio;

  //== Reference histogram.
  if( plot_as == "E") { 
    Get_DetEnergy( -1 , hReference); 
    Get_GenEnergy_response( 0, hGen);
    Get_DetEnergy( 0, hMC );
  }
  else if( plot_as == "pT") { 
    Get_DetPt( -1 , hReference); 
    Get_GenPt_response( 0, hGen);
    Get_DetPt( 0, hMC );
  }
  else{
   return;
  }

  //-- Normalize hReference (DATA) to xsec by dividing by its luminosity.
  hReference->Scale( 1./ lumi_ );  
  hReference->GetYaxis()->SetTitle( "d#sigma/dE [mb/GeV]" );
  
  //== Unfold the data.
  if( plot_as == "E" ){ 	Get_DetUnfolded( 0 , hReference, iterations, variable); }	// Unfold hReference with the Response from MC sample 0.
  else if( plot_as == "pT"){ 	Get_DetUnfolded_pt( 0 , hReference, iterations, variable); }

  //== dN/dx
  SetDnDx( hReference );

  //== Benoit.
  {
    hReference_benoit_dndx = (TH1D*)hReference->Clone("Data_dsigmadE");
    hLumi_up_benoit_dndx = (TH1D*)hReference_benoit_dndx->Clone("hLumi_up");
    hLumi_up_benoit_dndx->Scale( 1.036);
    hLumi_down_benoit_dndx = (TH1D*)hReference_benoit_dndx->Clone("hLumi_down");
    hLumi_down_benoit_dndx->Scale( 0.964);
  }

  //== Determine the first bin above 300 GeV.
  int bin_300 = 0;
  while( hReference->GetXaxis()->GetBinLowEdge( bin_300 ) < 300. ){ bin_300++; }
  bin_300--;
  
  
  double Eplot_lowest = hReference->GetXaxis()->GetBinLowEdge( bin_300 );
  double Emax = hReference->GetXaxis()->GetBinLowEdge( hReference->GetNbinsX() );

  if( plot_as == "pT" ){ bin_300 = 1; }

  //== The last bin of the histogram is not the last bin of the plot: cut it away.
  SetSubhistogram_max( Emax );
  
  //== Get the integral of the area plotted.
  TH1D* hReference_copy = (TH1D*)hReference->Clone("copy");

  //== Legend (1): for models.
  double legx, legy;
  if( which_models == "Pythia"){ legx = 1. - pad_abs_->GetRightMargin()-0.47; legy = 1. - pad_abs_->GetTopMargin() - 0.45; }
  else{ legx = 1. - pad_abs_->GetRightMargin()-0.35; legy = 1. - pad_abs_->GetTopMargin() - 0.5; }

  TLegend* legend = new TLegend( 
	legx, 	1. - pad_abs_->GetTopMargin() - 0.45, 
	1. - pad_abs_->GetRightMargin(), 	1. - pad_abs_->GetTopMargin() );

  //== Legend (2):  for systematics.
  TLegend* legend_syst	 = new TLegend( 
	pad_abs_->GetLeftMargin(), 	pad_abs_->GetBottomMargin(), 
	pad_abs_->GetLeftMargin()+0.35, 	pad_abs_->GetBottomMargin() +0.5);
  legend_syst->AddEntry( hReference, TString::Format("Data"), "lp");

  //== Luminosity uncertainty: fixed 3.6%
  double lumi_dep = 0.036;

  /******************************************************************************************
  * Begin by averaging the models and taking the difference between average and each model. *
  ******************************************************************************************/

  //== Loop over MC files.

  TH1D* hReference_ratio = (TH1D*)hReference->Clone("hReference_ratio");
  hReference_ratio->Divide( hReference );


  /*************************************************************************
  The following is needed to calculate the proper position of the datapoints.
  See Eq. (7) in  Physics Research A 355 (1995) 541-547.

  We have a spectrum with two distinct slopes: the first slope is called low, the second is called high.
  **************************************************************************/
 
  cout << "===WHERE TO STICK YOUR DATA POINTS===" << endl;

  double alow = 185509.8;
  double blow = 0.007383852;   

  double ahigh = 24325.05;
  double bhigh = 0.004686864;

  double a = alow, b = blow;

  const Int_t n = hReference_ratio->GetNbinsX();
  // x-axis and the values.
  Double_t x[2*n], y[2*n];

  // Systematic uncertainty lumi.
  Double_t ex_lumi[2*n], ey_lumi[2*n];
  Calculate_lumi_systematic( iterations, ex_lumi, ey_lumi);
  
  // Systematic uncertainty model.
  Double_t ex_model[2*n], ey_model[2*n];  
  Calculate_model_systematic( iterations, ex_model, ey_model);
  
  // Systematic uncertainty position.
  Double_t ex_position_up[2*n], ey_position_up[2*n]; 
  Double_t ex_position_down[2*n], ey_position_down[2*n];  
  Calculate_position_systematic( iterations, ex_position_up, ey_position_up, ex_position_down, ey_position_down );   
  
  // Systematic uncertainty jetID.
  Double_t ex_id[2*n], ey_id[2*n]; 
  Calculate_jetID_systematic( iterations, ex_id, ey_id);
  
  // Systematic uncertainty JES.
  Double_t ex_jes_up[2*n], ey_jes_up[2*n]; 
  Double_t ex_jes_down[2*n], ey_jes_down[2*n];   
  Calculate_JES_systematic( iterations, ex_jes_up, ey_jes_up, ex_jes_down, ey_jes_down );  
  
  // SUM: Systematic uncertainty lumi or model.
  Double_t exl_M[2*n], eyl_M[2*n], exh_M[2*n], eyh_M[2*n];

  // SUM: Systematic uncertainty lumi, model, pos.
  Double_t exl_I[2*n], eyl_I[2*n], exh_I[2*n], eyh_I[2*n];
  
  // SUM: Systematic uncertainty lumi, model, pos, ID
  Double_t exl_P[2*n], eyl_P[2*n], exh_P[2*n], eyh_P[2*n];
  
  // SUM: Systematic uncertainty lumi, model, pos, ID
  Double_t exl_J[2*n], eyl_J[2*n], exh_J[2*n], eyh_J[2*n];   
 
  double x_data[n], y_data[n], y_ratio[n];
  double exl_data[n],exh_data[n];
  double eyl_data[n],eyh_data[n];
  double eyl_ratio[n],eyh_ratio[n];

  // Remember: hist bins start at 1, graph bins at 0.
  for(int bin = 0; bin < hReference->GetNbinsX(); bin++){
    double bin_lowedge = hReference->GetXaxis()->GetBinLowEdge( bin+1 );
    double bin_highedge = hReference->GetXaxis()->GetBinUpEdge( bin+1 );

    if( bin_lowedge > 900. ){ a = ahigh; b = bhigh; }

    double binwidth = hReference->GetXaxis()->GetBinWidth( bin+1);

    double xlw = bin_lowedge 
		+ 1./b* (log(b * binwidth) )
		- 1./b* (log( 1. - exp( -1. * b * binwidth) ));

    int bin_shift = 1;
    if (bin == n-1){ bin_shift = 0; }	//== This ensures the last ratio is drawn neatly.

    x_data[bin] = xlw;
    y_data[bin] = hReference->GetBinContent( bin+bin_shift );
    y_ratio[bin] = 1.;
    exl_data[bin] = xlw - bin_lowedge;
    exh_data[bin] = bin_highedge - xlw;
    eyl_data[bin] = hReference->GetBinError( bin+bin_shift );
    eyh_data[bin] = hReference->GetBinError( bin+bin_shift );
    eyl_ratio[bin] = hReference->GetBinError( bin+bin_shift )/hReference->GetBinContent( bin+bin_shift  );
    eyh_ratio[bin] = hReference->GetBinError( bin+bin_shift )/hReference->GetBinContent( bin+bin_shift  );
  }

   TGraphAsymmErrors *gr_data = new TGraphAsymmErrors(n,x_data,y_data,exl_data,exh_data,eyl_data,eyh_data);
     gr_data->SetLineWidth( 1 );
     gr_data->SetMarkerSize( 1.25 );
     gr_data->SetMarkerStyle( 20 );
     gr_data->SetTitle("Graph_data");
     gr_data->SetName("Graph_data");


   TGraphAsymmErrors *gr_ratio = new TGraphAsymmErrors(n,x_data,y_ratio,exl_data,exh_data,eyl_ratio,eyh_ratio);
     gr_ratio->SetLineWidth( 1 );
     gr_ratio->SetMarkerSize( 1.25 );
     gr_ratio->SetMarkerStyle( 20 );
     gr_ratio->SetTitle("Graph_ratio");
     gr_ratio->SetName("Graph_ratio");


  for(int bin = 0; bin < n; bin++){
    int bin_shift = 1;
    if (bin == n-1){ bin_shift = 0; }	//== This ensures the last ratio is drawn neatly.
    if (bin == 2){ bin_shift = 2; }	//== This ensures the last ratio is drawn neatly.

    x[2*bin] = hReference_ratio->GetXaxis()->GetBinLowEdge( bin+1 );
    y[2*bin] = hReference_ratio->GetBinContent( bin+1 );
    x[2*bin+1] = hReference_ratio->GetXaxis()->GetBinUpEdge( bin+1 );
    y[2*bin+1] = hReference_ratio->GetBinContent( bin+1);

    //== MODEL.
    eyl_M[2*bin] = sqrt(  ey_model[2*bin]*ey_model[2*bin]  );
    eyh_M[2*bin] = sqrt(  ey_model[2*bin]*ey_model[2*bin]  );
    eyl_M[2*bin+1] = sqrt(  ey_model[2*bin+1]*ey_model[2*bin+1]  );
    eyh_M[2*bin+1] = sqrt(  ey_model[2*bin+1]*ey_model[2*bin+1]  );
    
    //== JetID.
    eyl_I[2*bin] = sqrt( ey_id[2*bin]*ey_id[2*bin]  );
    eyh_I[2*bin] = sqrt( ey_id[2*bin]*ey_id[2*bin]  );   
    eyl_I[2*bin+1] = sqrt( ey_id[2*bin+1]*ey_id[2*bin+1]  );
    eyh_I[2*bin+1] = sqrt( ey_id[2*bin+1]*ey_id[2*bin+1]  );    
    
    //==POSITION.
    eyl_P[2*bin] = sqrt(ey_position_down[2*bin]*ey_position_down[2*bin]  );
    eyh_P[2*bin] = sqrt( ey_position_up[2*bin]*ey_position_up[2*bin]  );    
    eyl_P[2*bin+1] = sqrt( ey_position_down[2*bin+1]*ey_position_down[2*bin+1]  );
    eyh_P[2*bin+1] = sqrt( ey_position_up[2*bin+1]*ey_position_up[2*bin+1]  );   
    
    //== JES.
    eyl_J[2*bin] = sqrt( ey_jes_down[2*bin]*ey_jes_down[2*bin]  );
    eyh_J[2*bin] = sqrt(  ey_jes_up[2*bin]*ey_jes_up[2*bin]  ); 
    eyl_J[2*bin+1] = sqrt( ey_jes_down[2*bin+1]*ey_jes_down[2*bin+1]  );
    eyh_J[2*bin+1] = sqrt(  ey_jes_up[2*bin+1]*ey_jes_up[2*bin+1]  );  
    
    cout << "Bin\t" << bin << "\tDOWN\t" << ey_jes_down[2*bin] << "\tUP\t" << ey_jes_up[2*bin] << endl;
  }
  
  pad_ratio_->cd();

  int style_ = 2;

  //== All uncertainty.
  TGraphAsymmErrors *gr = new TGraphAsymmErrors(2*n,x,y,exl_data,exh_data,eyl_J,eyh_J);
  gr->SetTitle("TGraphAsymmErrors All");
  int ci_all = TColor::GetColor("#FFFF00");
  gr->SetMarkerSize(0);
  gr->SetFillColor( kWhite );
  gr->SetLineColor( ci_all );
  gr->SetLineWidth( 3 );  
  gr->SetLineStyle( style_++);
  gr->SetMarkerStyle(21);

  gr->GetHistogram()->GetXaxis()->SetRangeUser(300., Emax_);
  gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
  if( plot_as == "xf" ){    gr->GetHistogram()->GetXaxis()->SetTitle("x_{F}"); }
  if( plot_as == "pT" ){    gr->GetHistogram()->GetXaxis()->SetTitle("p_{T} [GeV]"); }
  gr->GetHistogram()->GetYaxis()->SetTitle("Ratio");
  gr->GetHistogram()->GetYaxis()->SetNdivisions(205);
  gr->GetHistogram()->GetYaxis()->SetRangeUser(0, 3.35);

  //Prepare_1Dplot( gr, pad_ratio_ );
  Prepare_1Dplot_ratio( gr );
  gr->GetHistogram()->GetXaxis()->SetRangeUser(300., Emax_);
  gr->Draw("AE3");
  gr->Draw("E1same");

  
  //== Lumi + position + models + jetID
  gr = new TGraphAsymmErrors(2*n,x,y,exl_data,exh_data,eyl_P,eyh_P);
  int ci_id = TColor::GetColor("#000099");  
  gr->SetFillColor( ci_id );
  if( plot_as == "E"){ 
//    gr->GetHistogram()->GetXaxis()->SetRangeUser(Eplot_lowest, Emax_);
    gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
  }
  gr->GetHistogram()->GetYaxis()->SetRangeUser( 0.1*GetMinimumValue( hReference ), hReference->GetMaximum() * 1.1 );
  gr->GetHistogram()->GetYaxis()->SetTitle("#frac{dN}{dE}");
  gr->Draw("E1same");

  //== Model + position + luminosity uncertainty. 
  gr = new TGraphAsymmErrors(2*n,x,y,ex_lumi,ex_lumi,eyl_I,eyh_I);
  gr->SetTitle("TGraphAsymmErrors No JES");
  int ci = TColor::GetColor("#C0C0C0");
  gr->SetMarkerSize(0);
  gr->SetLineColor( ci );
  
  gr->SetMarkerStyle(21);

  gr->GetHistogram()->GetXaxis()->SetRangeUser(Eplot_lowest, Emax_);
  gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
  gr->GetHistogram()->GetYaxis()->SetRangeUser(0., 2.);
  gr->Draw("PE1same");

  //== Position + luminosity uncertainty.
  gr = new TGraphAsymmErrors(2*n,x,y,ex_lumi,ex_lumi,eyl_M,eyh_M);
  gr->SetTitle("TGraphAsymmErrors No JES, No model");
  int ci_noModel = TColor::GetColor("#33FF99");
  gr->SetMarkerSize(0);
  gr->SetLineColor( ci_noModel );
  gr->SetMarkerStyle(21);

  gr->GetHistogram()->GetXaxis()->SetRangeUser(Eplot_lowest, Emax_);
  gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
  gr->GetHistogram()->GetYaxis()->SetRangeUser(0., 2.);
  gr->Draw("PE1same");

  //== Luminosity uncertainty.
  gr = new TGraphAsymmErrors(2*n,x,y,ex_lumi,ex_lumi,ey_lumi,ey_lumi);
  gr->SetTitle("TGraphAsymmErrors Lumi");
  int ci_lumi =   TColor::GetColor("#FFCCCC");

  gr->SetMarkerSize(0);
  gr->SetLineColor( ci_lumi );
  gr->SetMarkerStyle(21);

  gr->GetHistogram()->GetXaxis()->SetRangeUser(Eplot_lowest, Emax_);
  gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
  gr->GetHistogram()->GetYaxis()->SetRangeUser(0., 2.);
  gr->Draw("PE1same");

  hReference_ratio->SetMarkerColor( kBlack );
  gr_ratio->GetHistogram()->GetXaxis()->SetRangeUser(Eplot_lowest, Emax_);
  gr_ratio->Draw("psame");

  //-- Absolute value.

  cout << "===ABSOLUTE===" << endl;
  pad_abs_->cd();

  for(int bin = 0; bin < n; bin++){

    int bin_shift = 1;
    if (bin == n-1){ bin_shift = 0; }	//== This ensures the last ratio is drawn neatly.
    if (bin == 2){ bin_shift = 2; }	//== This ensures the last ratio is drawn neatly.

    x[2*bin] = hReference->GetXaxis()->GetBinLowEdge( bin+1 );
    y[2*bin] = hReference->GetBinContent( bin+bin_shift );
    x[2*bin+1] = hReference->GetXaxis()->GetBinUpEdge( bin+1 );
    y[2*bin+1] = hReference->GetBinContent( bin+bin_shift );

    //-- Error band lumi.
    ey_lumi[2*bin] = 0.036 * y[2*bin] ;
    ex_lumi[2*bin] = hReference_ratio->GetBinContent( bin+1 );

    ey_lumi[2*bin+1] = 0.036 * y[2*bin] ;
    ex_lumi[2*bin+1] = hReference_ratio->GetBinWidth( bin+1 )/2;

    //-- Error band lumi + position.
    eyl_M[2*bin] = hReference->GetBinContent( bin+bin_shift ) *  eyl_M[2*bin];
    eyh_M[2*bin] = hReference->GetBinContent( bin+bin_shift ) *  eyh_M[2*bin];
    exl_M[2*bin] = hReference->GetBinWidth( bin+bin_shift )/2;
    exh_M[2*bin] = hReference->GetBinWidth( bin+bin_shift )/2;

    eyl_M[2*bin+1] = hReference->GetBinContent( bin+bin_shift ) *  eyl_M[2*bin+1];
    eyh_M[2*bin+1] = hReference->GetBinContent( bin+bin_shift ) *  eyh_M[2*bin+1];
    exl_M[2*bin+1] = hReference->GetBinWidth( bin+bin_shift )/2;
    exh_M[2*bin+1] = hReference->GetBinWidth( bin+bin_shift )/2;

    //-- Error band lumi + position + model.
    eyl_I[2*bin] = hReference->GetBinContent( bin+bin_shift ) *  eyl_I[2*bin];
    eyh_I[2*bin] = hReference->GetBinContent( bin+bin_shift ) *  eyh_I[2*bin];
    exl_I[2*bin] = hReference->GetBinWidth( bin+bin_shift )/2;
    exh_I[2*bin] = hReference->GetBinWidth( bin+bin_shift )/2;

    eyl_I[2*bin+1] = hReference->GetBinContent( bin+bin_shift ) *  eyl_I[2*bin+1];
    eyh_I[2*bin+1] = hReference->GetBinContent( bin+bin_shift ) *  eyh_I[2*bin+1];
    exl_I[2*bin+1] = hReference->GetBinWidth( bin+bin_shift )/2;
    exh_I[2*bin+1] = hReference->GetBinWidth( bin+bin_shift )/2;
    
    //-- Error band lumi + position + model + jetID. 

    eyl_P[2*bin] = hReference->GetBinContent( bin+bin_shift ) *  eyl_P[2*bin];
    eyh_P[2*bin] = hReference->GetBinContent( bin+bin_shift ) *  eyh_P[2*bin];
    exl_P[2*bin] = hReference->GetBinWidth( bin+bin_shift )/2;
    exh_P[2*bin] = hReference->GetBinWidth( bin+bin_shift )/2;     

    eyl_P[2*bin+1] = hReference->GetBinContent( bin+bin_shift ) *  eyl_P[2*bin+1];   
    eyh_P[2*bin+1] = hReference->GetBinContent( bin+bin_shift ) *  eyh_P[2*bin+1];
    exl_P[2*bin+1] = hReference->GetBinWidth( bin+bin_shift )/2;
    exh_P[2*bin+1] = hReference->GetBinWidth( bin+bin_shift )/2;  

    //-- Error band with JES.
    eyl_J[2*bin] = hReference->GetBinContent( bin+bin_shift ) * eyl_J[2*bin];
    eyh_J[2*bin] = hReference->GetBinContent( bin+bin_shift ) * eyh_J[2*bin];
    exl_J[2*bin] = hReference->GetBinWidth( bin+bin_shift )/2;
    exh_J[2*bin] = hReference->GetBinWidth( bin+bin_shift )/2;

    eyl_J[2*bin+1] = hReference->GetBinContent( bin+bin_shift ) * eyl_J[2*bin+1];
    eyh_J[2*bin+1] = hReference->GetBinContent( bin+bin_shift ) * eyh_J[2*bin+1];
    exl_J[2*bin+1] = hReference->GetBinWidth( bin+bin_shift )/2;
    exh_J[2*bin+1] = hReference->GetBinWidth( bin+bin_shift )/2;
  }
   
  pad_abs_->SetLogy();

  TGraphAsymmErrors *gr_all, *gr_id, *gr_mod, *gr_pos, *gr_lumi;

  //-- Draw graph with smaller errors on top of graph with bigger errors.
  gr = new TGraphAsymmErrors(2*n,x,y,exl_J,exh_J,eyl_J,eyh_J);
  gr->SetFillColor(ci_all);
  gr->GetHistogram()->GetXaxis()->SetRangeUser(300., Emax_);
  gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
 
  double ylow =   0.5*GetMinimumValue( hReference );
  double yhigh = hReference->GetMaximum() * 1.1;

  ylow = 2*1e-5;
  yhigh = 0.5;
  gr->GetHistogram()->GetYaxis()->SetRangeUser(ylow, yhigh );

  if( plot_as == "pT"){ gr->GetHistogram()->GetYaxis()->SetTitle("d#sigma/dp_{T} [mb/GeV]"); }
  if( plot_as == "xf"){ gr->GetHistogram()->GetYaxis()->SetTitle("d#sigma/dx_{F} [mb]"); }
  if( plot_as == "E"){ gr->GetHistogram()->GetYaxis()->SetTitle("d#sigma/dE [mb/GeV]"); }

  //== Only Y-axis needed.
  Prepare_1Dplot( gr );
  gr->GetHistogram()->GetXaxis()->SetRangeUser(300., Emax_);
  gr->Draw("AE1");
  gr_all = (TGraphAsymmErrors*)gr->Clone();

  //legend_syst->AddEntry( gr, "Syst. errors (All)", "f");
  
  //== Lumi + position + models + jetID
  gr = new TGraphAsymmErrors(2*n,x,y,exl_P,exh_P,eyl_P,eyh_P);
  gr->SetFillColor( ci_id );
  if( plot_as == "E"){ 
    gr->GetHistogram()->GetXaxis()->SetRangeUser(Eplot_lowest, Emax_);
    gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
  }
  gr->GetHistogram()->GetYaxis()->SetRangeUser( 0.1*GetMinimumValue( hReference ), hReference->GetMaximum() * 1.1 );
  gr->GetHistogram()->GetYaxis()->SetTitle("#frac{dN}{dE}");
  gr->Draw("E1same");
  gr_id = (TGraphAsymmErrors*)gr->Clone();
  //legend_syst->AddEntry( gr, "Syst. errors (lumi., pos., mod., jetID)", "f");  

  //== Lumi + position + models/
  gr = new TGraphAsymmErrors(2*n,x,y,exl_I,exh_I,eyl_I,eyh_I);
  gr->SetFillColor(ci);
  if( plot_as == "E"){ 
    gr->GetHistogram()->GetXaxis()->SetRangeUser(Eplot_lowest, Emax_);
    gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
  }
  gr->GetHistogram()->GetYaxis()->SetRangeUser( 0.1*GetMinimumValue( hReference ), hReference->GetMaximum() * 1.1 );
  gr->GetHistogram()->GetYaxis()->SetTitle("#frac{dN}{dE}");
  gr->Draw("E1same");
  gr_mod = (TGraphAsymmErrors*)gr->Clone();
  //legend_syst->AddEntry( gr, "Syst. errors (lumi., pos., mod.)", "f");

  //== Lumi + position.
  gr = new TGraphAsymmErrors(2*n,x,y,exl_M,exh_M,eyl_M,eyh_M);
  gr->SetFillColor(ci_noModel);
  if( plot_as == "E" ){
    gr->GetHistogram()->GetXaxis()->SetRangeUser(Eplot_lowest, Emax_);
    gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
  }
  gr->GetHistogram()->GetYaxis()->SetRangeUser( 0.1*GetMinimumValue( hReference ), hReference->GetMaximum() * 1.1 );
  gr->GetHistogram()->GetYaxis()->SetTitle("#frac{dN}{dE}");
  gr->Draw("E1same");
  gr_pos = (TGraphAsymmErrors*)gr->Clone();

  //legend_syst->AddEntry( gr, "Syst. errors (lumi. + pos.)", "f");

  //-- Draw graph with smaller errors on top of graph with bigger errors.
  gr = new TGraphAsymmErrors(2*n,x,y,exl_data,exh_data,ey_lumi,ey_lumi);
  gr->SetFillColor(ci_lumi);
  if( plot_as == "E" ){
    gr->GetHistogram()->GetXaxis()->SetRangeUser(Eplot_lowest, Emax_);
    gr->GetHistogram()->GetXaxis()->SetTitle("E [GeV]");
  }
  gr->GetHistogram()->GetYaxis()->SetRangeUser( 0.1*GetMinimumValue( hReference ), hReference->GetMaximum() * 1.1 );
  gr->GetHistogram()->GetYaxis()->SetTitle("#frac{dN}{dE}");
  gr->Draw("E1same");
  gr_lumi = (TGraphAsymmErrors*)gr->Clone();

  //legend_syst->AddEntry( gr, "Syst. errors (lumi.)", "f");

  hReference->SetMarkerColor( kBlack );

  gr_data->GetHistogram()->GetXaxis()->SetRangeUser(Eplot_lowest, Emax_);
  gr_data->Draw("psame");


  legend_syst->AddEntry( gr_lumi, "Syst. errors (lumi.)", "f");
  legend_syst->AddEntry( gr_pos, "Syst. errors (+ mod.)", "f");
  legend_syst->AddEntry( gr_mod, "Syst. errors (+ id)", "f");  
  legend_syst->AddEntry( gr_id, "Syst. errors (+ pos.)", "f");   
  legend_syst->AddEntry( gr_all, "Syst. errors (+ JES)", "f");
  
  //-- Scale and draw generator level distribution.
  int color_ = 0;
  int line_ = 1;
  drawoptions = "][histsame";
  ofstream binwidths;
  binwidths.open("binwidths.txt");

  TString split_MCs = "";
  int count_MCs = 0;

  cout << "===MCs===" << endl;
  for(int MC_ = 0; MC_ < MC_files_.size(); MC_++){

    if( plot_as == "pT" ){ continue; }
    //-- Skip colors if needed.
    color_++ ;
    bool goodcolor = false;
    while( !goodcolor ){
      if( color_ == 1 || color_ == 8 || color_ == 9 || color_ == 10 || color_ == 11) color_++;
      else{ goodcolor = true; }
    }

    TH1D* hGen;
    //-- Only continue if the MC is one of the following (excludes position samples).
    if( ((set_of_tags_["mc_type"])[MC_files_[MC_]] != "model" && 
	(set_of_tags_["mc_type"])[MC_files_[MC_]] != "actual" &&  
	(set_of_tags_["mc_type"])[MC_files_[MC_]] != "shift_MPI_or_Tune" &&
	(set_of_tags_["mc_type"])[MC_files_[MC_]] != "model_") ){ continue; }

    //-- Split up the canvas.
    //== "Pythia" will plot only models from Pythia.
    //== "noPythia" will plot models without Pythia.
    //== "All" plots all models.
    //== Any other entry skips the addition of models.
    if( which_models == "Pythia" ){
      if( !MC_files_[MC_].Contains( "ythia") ){ continue; }
    }
    else if( which_models == "noPythia"){
      if( MC_files_[MC_].Contains( "ythia") ){ continue; }
    }
    else if( which_models != "All"){ break; }

    //-- Check for xsec = 0 or NaN.
    if( xsec_[ MC_files_[MC_] ] !=  xsec_[ MC_files_[MC_] ] ){ continue; }
    if( xsec_[ MC_files_[MC_] ] ==  0. ){ continue; }
   
    //-- If scaletodata_ is true, the gen. energy spectrum will be scaled with numbers taken from the file before being passed back.
    scaletodata_ = true;

     // There is only one gen. distribution, since there are no cuts on detector level.
    if( (set_of_tags_["mc_type"])[MC_files_[MC_]] == "shift_MPI_or_Tune" ){ 
      if( plot_as == "E" ){		Get_GenEnergy( MC_, hGen ); }
      else if( plot_as == "pT" ){ 	Get_GenPt( MC_, hGen ); }
    }

    // There are several gen. distributions, take the one uninfluenced by the cuts on detector level.
    else{ 
      if( plot_as == "E" ){ 		Get_PureGenEnergy( MC_, hGen ); }
      else if( plot_as == "pT" ){ 	Get_PureGenPt( MC_, hGen ); }
    }
    scaletodata_ = false;
    if( !hGen ){ continue; }

    hGen->SetTitle( TString::Format("hGenJet_energy_%i", MC_) );
    hGen->SetName( TString::Format("hGenJet_energy_%i", MC_) );

    //-- Skip to next file if histogram does not exist.
    if( hGen->Integral() != hGen->Integral() ){ continue; }
 
    //== Get the integral of the area plotted.
    TH1D* hGen_copy;  
    if( plot_as == "E" ){ GetSubHistogram( hGen, hGen_copy, Eplot_lowest, Emax_); }

    xsec_txt << "Integral\t" << legend_info_gen_[MC_files_[MC_]] << "\t" << hGen_copy->Integral() << "\t" <<  endl;   

    SetDnDx( hGen );

    if( (set_of_tags_["mc_type"])[MC_files_[MC_]] == "shift_MPI_or_Tune" ){
      drawoptions = "][histsame";
    }

    //-- Multiply the distribution to the total number of measured data.
    //-- Set underflow bin to the same value to avoid dive to zero at low bin.
    hGen->GetXaxis()->SetRangeUser( Eplot_lowest, Emax);

    //-- Set line properties.
    hGen->SetLineColor( getColor(color_) );
    hGen->SetLineStyle( (line_ != 1)*((line_++)%7+2) );
    hGen->SetLineWidth( 3 );

    //-- Set ratio of distribution.
    hGen_ratio = (TH1D*)hGen->Clone( TString::Format("hGen_ratio_%i", MC_) );
    hGen_ratio->Divide( hReference );

    pad_abs_->cd();
    if( plot_as == "E" ){ hGen->GetXaxis()->SetRangeUser(Eplot_lowest, Emax_); }
    hGen->Draw( drawoptions);
    
    pad_ratio_->cd();
    if( plot_as == "E" ){ hGen_ratio->GetXaxis()->SetRangeUser(Eplot_lowest, Emax_); }
    hGen_ratio->Draw( drawoptions);

    TString legend_file = legend_info_gen_[MC_files_[MC_]];

    if( (set_of_tags_["mc_type"])[MC_files_[MC_]] == "shift_MPI_or_Tune" ){
      legend->AddEntry( hGen, TString::Format( legend_file  ), "l");
    }
    else{
      legend->AddEntry( hGen, legend_file , "l");
    }
    count_MCs++;
  }

   pad_abs_->cd();

   //== Adjust height of legend.
   legend->SetY1( 1. - pad_abs_->GetTopMargin() - .09 * count_MCs   );
   legend->SetX2( 1. - pad_abs_->GetRightMargin() - .35 + (.07 * count_MCs)   );
   legend->SetFillStyle(0);
   legend->SetLineWidth(0);
   legend->Draw();
  

   pad_abs_->cd();
   pad_abs_->RedrawAxis("g");
   legend_syst->SetFillStyle(0);
   legend_syst->SetLineWidth(0);
   legend_syst->SetTextSize( legend->GetTextSize() );   
   legend_syst->Draw();

   pad_ratio_->RedrawAxis("g");

    can_->cd();

    //== Text.
    Finish_canvas( pad_abs_ , "leftish");

   can_->SaveAs( TString::Format( folder_ + "Totaldependence_separate_" + plot_as + "_" + which_models + "_%iit_deltaphi_0%i_etaband_0%i" + cutsetup_ + ".C", iterations, static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_)) );
   can_->SaveAs( TString::Format( folder_ + "Totaldependence_separate_" + plot_as + "_" + which_models + "_%iit_deltaphi_0%i_etaband_0%i" + cutsetup_ + ".pdf", iterations, static_cast<int>(10. * deltaPhiMax_), static_cast<int>(10. * etawidth_)) ); 

   systematics_txt.close();
   xsec_txt.close();
   
   PlotSystematics(iterations);
}

























void Unfolder::Determine_resolution( TString setup ){

   //== Prepare canvas.
   TCanvas *can_1D, *can_2D;
   PrepareCanvas(can_1D, "fit");
   PrepareCanvas_2D(can_2D, "can_2D");

   // Take the necessary file.
   TFile *_file0 = TFile::Open("/user/avanspil/public/CastorJets/ak5ak5_Pythia6Z2star_ShowerLibrary_NomGeo_" + setup + "_Emin_150.000000_deltaPhiMax_0.200000_etaband_0.000000_all_matchE.root");     
   
   // Extract the 2D histogram.
   TH2D *hJER_per_energy = (TH2D*)_file0->Get("hJER_per_energy");
   	int Ebins = hJER_per_energy->GetNbinsY();
   	double Emin = 0., Emax = 3000.;
     
   // This fits a Gaussian to all slices with at least 20 non-zero bins.   
   hJER_per_energy->FitSlicesY(0,0,-1);
   TH1D* hJER_per_energy_2= (TH1D*)gDirectory->Get("hJER_per_energy_2");   
   TH1D *BinWidths = (TH1D*)hJER_per_energy_2->Clone("BinWidths");
   
   for(int i = 0; i < hJER_per_energy_2->GetNbinsX(); i++){
     double sigma = hJER_per_energy_2->GetBinContent(i);
     double egen = hJER_per_energy_2->GetBinCenter(i);
     
     BinWidths->SetBinContent(i, sigma*egen);     
   }
   
   can_1D->cd();
   BinWidths->Draw("hist");  
   Prepare_1Dplot( BinWidths );
   BinWidths->GetXaxis()->SetTitle("E_{gen} [GeV]");
   BinWidths->GetYaxis()->SetTitle("#sigma (E_{gen})");
   
   TF1 *f_linear = new TF1("f_linear", "[0] * x + [1]", 0.,1600.);

   gStyle->SetStatStyle(1);
   TFitResultPtr r = BinWidths->Fit("f_linear", "RS");   

   //== round parameters to two decimals.
   double r0 = r->Parameter(0);
   int r0_100 = static_cast<int>(r0 * 100.);
   int r0_dec = r0_100%100;
   int r0_int = (r0_100 - r0_dec)/100;
   TString par0 = TString::Format( "%i.%i", r0_int, r0_dec);

   double r1 = r->Parameter(1);
   int r1_100 = static_cast<int>(r1 * 100.);
   int r1_dec = r1_100%100;
   int r1_int = (r1_100 - r1_dec)/100;
   TString par1 = TString::Format( "%i.%i", r1_int, r1_dec);


   TLatex *   tex = new TLatex( can_1D->GetLeftMargin(),1. - can_1D->GetTopMargin(),"7 TeV (pp)");
//   tex->SetNDC()
   tex->SetTextAlign(11);
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->DrawLatex( 
	1200., 
	300.,
	TString::Format( "#sigma = " + par0 + "*E_{gen} + " + par1 ) );  

   can_1D->SaveAs( folder_ + "Resolution_egen_"  + setup +  ".pdf");
   can_1D->SaveAs( folder_ + "Resolution_egen_"  + setup +  ".C");

   can_2D->cd();
   Prepare_2Dplot( hJER_per_energy );
   hJER_per_energy->GetYaxis()->SetRangeUser( -1., 2.5);
   hJER_per_energy->GetYaxis()->SetTitle( "JES");
   hJER_per_energy->Draw("colz");
   can_2D->SetLogz();
   can_2D->SaveAs( folder_ + "JES_vs_Egen_"  + setup +  ".pdf");
   can_2D->SaveAs( folder_ + "JES_vs_Egen_"  + setup +  ".pdf");
}






void Unfolder::Determine_resolution_eDet( TString setup ){

   //== Prepare canvas.
   TCanvas *can_1D, *can_2D;
   PrepareCanvas(can_1D, "fit");
   PrepareCanvas_2D(can_2D, "can_2D");

   // Take the necessary file.
   TFile *_file0 = TFile::Open("/user/avanspil/public/CastorJets/ak5ak5_Pythia6Z2star_ShowerLibrary_NomGeo_" + setup + "_Emin_150.000000_deltaPhiMax_0.200000_etaband_0.000000_all_matchE.root");     
   
   // Extract the 2D histogram.
   TH2D *hJER_per_energy = (TH2D*)_file0->Get("hJER_per_eDet");
   	int Ebins = hJER_per_energy->GetNbinsY();
   	double Emin = 0., Emax = 3000.;
     
   /*
   // This fits a Gaussian to all slices with at least 20 non-zero bins.   
   hJER_per_energy->FitSlicesY(0,0,-1);
   TH1D* hJER_per_energy_2= (TH1D*)gDirectory->Get("hJER_per_eDet_2");   
   */
   TH1D *BinWidths = (TH1D*)hJER_per_energy->ProjectionX( "Binwidths",0, 0 );
   BinWidths->Reset();
   TH1D *BinWidths2 = (TH1D*)hJER_per_energy->ProjectionX( "Binwidths2",0, 0 );
   BinWidths2->Reset();

   TCanvas *can = new TCanvas("can", "can", 1.);
   TString drawoptions = "hist";   

   for(int i = 0; i < hJER_per_energy->GetNbinsX(); i++){

     TH1D* hJER = (TH1D*)hJER_per_energy->ProjectionY( TString::Format("bin_%i", i), i, i );
     hJER->SetLineColor( getColor(i) );
     hJER->Draw( drawoptions );
     drawoptions = "histsame";

     double sigma = hJER->GetRMS(1);
     double egen = hJER_per_energy->GetXaxis()->GetBinCenter(i);


     BinWidths->SetBinContent(i, sigma);     
     BinWidths2->SetBinContent(i, sigma*egen);     
   }

   can->SaveAs( folder_ + "Projections.pdf");
   
   can_1D->cd();
   BinWidths->SetLineColor( kRed );
   BinWidths->SetLineWidth(3);
   BinWidths->Draw("p");  

   BinWidths->GetXaxis()->SetTitle("E_{det} [GeV]");
   BinWidths->GetYaxis()->SetTitle("#sigma (E_{det})");

   can_1D->SaveAs( folder_ + "Resolution_edet_"  + setup +  ".pdf");
   can_1D->SaveAs( folder_ + "Resolution_edet_"  + setup +  ".C");

   BinWidths2->SetLineColor( kRed );
   BinWidths2->SetLineWidth(3);
   BinWidths2->Draw("p");  
   BinWidths2->GetXaxis()->SetTitle("E_{det} [GeV]");
   BinWidths2->GetYaxis()->SetTitle("#sigma (E_{det})");

   can_1D->SaveAs( folder_ + "Resolution_abs_"  + setup +  ".pdf");
   can_1D->SaveAs( folder_ + "Resolution_abs_"  + setup +  ".C");



   can_2D->cd();
   Prepare_2Dplot( hJER_per_energy );
   hJER_per_energy->GetYaxis()->SetRangeUser( -1., 2.5);
   hJER_per_energy->GetYaxis()->SetTitle( "JES");
   hJER_per_energy->Draw("colz");
   can_2D->SetLogz();
   can_2D->SaveAs( folder_ + "JES_vs_Edet_"  + setup +  ".pdf");
   can_2D->SaveAs( folder_ + "JES_vs_Edet_"  + setup +  ".pdf");
}





void Unfolder::Plot_Egen_pt(){

  TFile* _file = TFile::Open( TString::Format("/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo_unfold_Emin_150.000000_deltaPhiMax_%f_etaband_%f_all_matchE.root", phimax_global_, etaacc_global_) , "Read");

  TCanvas *can;
  PrepareCanvas_2D( can, "pt_egen");

  // Extract 2D pt-E distribution from file.

  TH2D* hEgen_pt = (TH2D*)_file->Get("hEgen_pt");

  hEgen_pt->Draw("colz");

  // Plot pt = E/ cosh(<eta>)

  TH1D* hist_avgEta;
  Calculate_averageEta_perEgen( 0 , hist_avgEta );

  TH1D* hPtFromE = (TH1D*)hist_avgEta->Clone("hPtFromE");

  for(int bin = 1; bin <= hPtFromE->GetNbinsX(); bin++){

    hPtFromE->SetBinContent( bin,
		 hEgen_pt->GetXaxis()->GetBinCenter(bin) /
		 cosh( hPtFromE->GetBinContent( bin ) )
);
    cout << bin << "\t" << hEgen_pt->GetXaxis()->GetBinCenter(bin) <<  hist_avgEta->GetBinContent( bin ) << "\t" << hPtFromE->GetBinContent( bin ) << endl;

  }

  hPtFromE->SetMarkerSize( hPtFromE->GetMarkerSize() * 2 );
  hPtFromE->SetMarkerStyle( 20 );
  
  hPtFromE->Draw("histsame");


  //== <pt> per Egen.

  TH1D* hist_avgpT;
  Calculate_averagepT_perEgen( 0, hist_avgpT);

  hist_avgpT->SetMarkerSize( hist_avgpT->GetMarkerSize() * 2);
  hist_avgpT->SetMarkerStyle(24);
  hist_avgpT->SetMarkerColor( kRed );

  hist_avgpT->Draw("histsame");


  // Plot E = pT/ cosh(<eta>)

  TH2D* hEta_pt = (TH2D*)_file->Get("hEta_pt");

  TH1D* hist_avgEta_pt;
  Calculate_averageEta_perpt( 0 , hist_avgEta_pt );

	cout << "\tPrepare graph" << endl;     

  TGraph* hEfromPt;
  double pT_calc[hist_avgEta_pt->GetNbinsX()];
  double E_calc[hist_avgEta_pt->GetNbinsX()];

	if( !hEta_pt ){ cout << "Schiiiiiiiiiiit!" << endl; }

	cout << "\tPrepared graph" << "\t\t" << hEta_pt->GetNbinsY() << endl;     

  for(int bin = 1; bin <= hist_avgEta_pt->GetNbinsX(); bin++){

	cout << "\t" << bin << "\t" << hEta_pt->GetYaxis()->GetBinCenter(bin) << "\t" << cosh( hist_avgEta_pt->GetBinContent( bin ) ) << endl;

    pT_calc[bin-1] = hEta_pt->GetYaxis()->GetBinCenter(bin);
    E_calc[bin-1] = pT_calc[bin-1] * cosh( hist_avgEta_pt->GetBinContent( bin ) );
  }

  hEfromPt = new TGraph( hist_avgEta_pt->GetNbinsX(), E_calc , pT_calc );

  hEfromPt->SetMarkerStyle( 24 );
  hEfromPt->SetMarkerColor( kOrange );
  hEfromPt->SetMarkerSize(2);

  hEfromPt->Draw("psame");


  //== <E> per pT

  TGraph* graph_avgE;
  Calculate_averageEgen_perpT( 0, graph_avgE);
 
  graph_avgE->SetMarkerStyle(21);
  graph_avgE->SetMarkerSize(2);
  graph_avgE->SetMarkerColor( kGreen + 1 );

  graph_avgE->Draw("psame");


  can->SetLogz();
  can->SaveAs(folder_ + "pT_vs_Egen.pdf");

  can->SaveAs(folder_ + "pT_vs_Egen.C");

}





void Unfolder::Calculate_averagepT_perEgen(int file_, TH1D* &hist_){

  cout << "\n//===//Calculate pt per E" << endl;

  TFile* _file = TFile::Open( TString::Format("/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo_unfold_Emin_150.000000_deltaPhiMax_%f_etaband_%F_all_matchE.root", phimax_global_, etaacc_global_), "Read");




  //=================
  //== All gen. jets.
  //=================

  TH2D* hEgen_eta = (TH2D*)_file->Get("hEgen_pt");
  TH1D* hZ = (TH1D*)hEgen_eta->ProjectionY();
 

  for(int bin_egen = 1; bin_egen <= hEgen_eta->GetNbinsX(); bin_egen++){
    double eta_total = 0.;
    double events_total = 0;
    for(int bin_eta = 1; bin_eta <= hEgen_eta->GetNbinsY(); bin_eta++){
//      cout << "(" << hEgen_eta->GetXaxis()->GetBinCenter(bin_egen) << "," << hEgen_eta->GetYaxis()->GetBinCenter(bin_eta) << ")\t" << endl;
      double nevents = hEgen_eta->GetBinContent( bin_egen, bin_eta );
      double eta_value = hEgen_eta->GetYaxis()->GetBinCenter( bin_eta  );

      eta_total += nevents * eta_value;
      events_total += nevents;
    }
    double eta_average = eta_total/events_total;
cout << hEgen_eta->GetXaxis()->GetBinCenter(bin_egen) << "\t" << eta_average << endl;
    hZ->SetBinContent( bin_egen, eta_average );
  }

  hist_ = hZ;

  cout << "\n//===//Calculated pt per E" << endl;
}




void Unfolder::Calculate_averageEgen_perpT(int file_, TGraph* &hist_){

  cout << "\n//===//Calculate <E> per pT" << endl;

  TFile* _file = TFile::Open( TString::Format("/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo_unfold_Emin_150.000000_deltaPhiMax_%f_etaband_%f_all_matchE.root", phimax_global_, etaacc_global_), "Read");


  //=================
  //== All gen. jets.
  //=================

  TH2D* hEgen_eta = (TH2D*)_file->Get("hEgen_pt");
 
  double pT_[ hEgen_eta->GetNbinsY() ];
  double E_[ hEgen_eta->GetNbinsY() ];


  for(int bin_egen = 1; bin_egen <= hEgen_eta->GetNbinsY(); bin_egen++){
    double eta_total = 0.;
    double events_total = 0;
    for(int bin_eta = 1; bin_eta <= hEgen_eta->GetNbinsX(); bin_eta++){
//      cout << "(" << hEgen_eta->GetXaxis()->GetBinCenter(bin_egen) << "," << hEgen_eta->GetYaxis()->GetBinCenter(bin_eta) << ")\t" << endl;
      double nevents = hEgen_eta->GetBinContent( bin_eta, bin_egen );
      double eta_value = hEgen_eta->GetXaxis()->GetBinCenter( bin_eta  );

      eta_total += nevents * eta_value;
      events_total += nevents;
    }
    double eta_average = eta_total/events_total;
cout << bin_egen << "\t" << hEgen_eta->GetYaxis()->GetBinCenter(bin_egen) << "\t" << eta_average << endl;
    pT_[bin_egen-1] = hEgen_eta->GetYaxis()->GetBinCenter(bin_egen);
    E_[bin_egen-1] =  eta_average ;
  }

  hist_ = new TGraph( hEgen_eta->GetNbinsY(), E_, pT_ ); 

  cout << "\n//===//Calculated <E> per pT" << endl;
}






void Unfolder::Calculate_averageEta_perpt(int file_, TH1D* &hist_){


  cout << "\n//===//Calculate eta per pT" << endl;

  TFile* _file = TFile::Open( MC_files_[file_], "READ");

  //=================
  //== All gen. jets.
  //=================

  TH2D* hEta_pt = (TH2D*)_file->Get("hEta_pt");
  TH1D* hZ = (TH1D*)hEta_pt->ProjectionX();
 
  cout << "hEta_pt\t" << hEta_pt->GetNbinsX() << "\t" << hEta_pt->GetNbinsY() << endl;

  for(int bin_egen = 1; bin_egen <= hEta_pt->GetNbinsY(); bin_egen++){
    double eta_total = 0.;
    double events_total = 0;
    for(int bin_eta = 1; bin_eta <= hEta_pt->GetNbinsX(); bin_eta++){
      cout << "(" << hEta_pt->GetXaxis()->GetBinCenter(bin_egen) << "," << hEta_pt->GetYaxis()->GetBinCenter(bin_eta) << ")\t" << endl;
      double nevents = hEta_pt->GetBinContent( bin_egen, bin_eta );
      double eta_value = hEta_pt->GetXaxis()->GetBinCenter( bin_eta  );

      eta_total += nevents * eta_value;
      events_total += nevents;
    }
    double eta_average = eta_total/events_total;

    hZ->SetBinContent( bin_egen, eta_average );
  }

  cout << "hZ\t" << hZ->GetNbinsX() << "\t" << hZ->Integral() << endl;

  hist_ = hZ;

  cout << "\n//===//Calculated eta per pT" << endl;
}






void Unfolder::DetermineScaling(){

  TFile* _fileMC = TFile::Open( TString::Format("/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo_unfold_Emin_150.000000_deltaPhiMax_%f_etaband_%f_all_matchE.root", phimax_global_, etaacc_global_ ), "READ" );
  TFile* _fileData = TFile::Open( "/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5_data_unfold_Emin_150.000000_all.root", "Read" );

  TH1D* hData = (TH1D*)_fileData->Get("hCastorJet_energy");

  TCanvas *can_;
  PrepareCanvas( can_, "canvas" );

  RooUnfoldResponse *response_all = (RooUnfoldResponse*)_fileMC->Get("response");
    RooUnfoldBayes unfold_bayes(response_all, hData, 1);
    // the original hist has become obsolete, transform it into the unfolded histogram.
    TH1D * hData_all = (TH1D*) unfold_bayes.Hreco( RooUnfold::kCovariance ); 

    TH1D* hMiss_orig = (TH1D*)_fileMC->Get("hCastorJet_miss_all");

  RooUnfoldResponse *response_noMiss = (RooUnfoldResponse*)_fileMC->Get("response_noMiss");
    RooUnfoldBayes unfold_bayes_noMiss(response_noMiss, hData, 1);
    // the original hist has become obsolete, transform it into the unfolded histogram.
    TH1D * hData_noMiss = (TH1D*) unfold_bayes_noMiss.Hreco( RooUnfold::kCovariance ); 

  hData_all->GetYaxis()->SetRangeUser( 1., hData_all->GetMaximum() * 1.1 );
  hData_all->Draw("hist");

  

  hData_noMiss->SetLineStyle( 2 );
  hData_noMiss->SetLineColor( kRed );
  hData_noMiss->Draw("histsame");

  cout << "\tµµµ\t" << hData_all->Integral() << "\tµµµ\t" << hData_noMiss->Integral() << endl;
  
  TH1D* hMisses = (TH1D*)hData_all->Clone("Scaled_misses");
  hMisses->Add( hData_noMiss, -1);

  hMisses->SetLineColor( kBlue );
  hMisses->SetLineStyle( 3 ) ;
  hMisses->Draw( "histsame" );

  hMiss_orig->SetLineColor( kGreen );
  hMiss_orig->SetLineStyle( 4 ) ;
  hMiss_orig->Draw("histsame");

  cout << "\tRatio Misses\t" << hMisses->Integral()/hMiss_orig->Integral() << endl;

  //== Calculate ratio of measured to data distributions.
  TH1D* hMeas_response = (TH1D*)response_all->Hmeasured();

    cout << "\tRatio Measured\t" << hData->Integral()/hMeas_response->Integral() << endl;

  

  TH1D* hRatio = (TH1D*)hMisses->Clone("Ratio");
  hRatio->Divide( hMiss_orig );
  hRatio->Draw("psame");
  

  can_->SetLogy();

  can_->SaveAs( folder_ + "Comparison.pdf");
  can_->SaveAs( folder_ + "Comparison.C");

  

}





void Unfolder::Calculate_jetID_systematic( int iterations, double *x, double *y ){
  cout << "=== Jet ID ===" << endl;
  
  
  ofstream uncertainties;
  uncertainties.open("Uncertainties.txt", ios::app);
  uncertainties << "\n\nID" <<  endl; 

  //== Unfold data as usual for reference.
    
  TH1D* hData; 
  Get_DetEnergy( -1 , hData);
  hData->Scale( 1./ lumi_ );
  Get_DetUnfolded( 0 , hData );      
  TH1D* hReference = (TH1D*)hData->Clone("hReference");
  SetDnDx( hData );
  
  //== Unfold data with all jets calibrated.
  TH1D* hData_overcal;
  TFile* _file_overcal = TFile::Open("/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5_data" + cutsetup_ + "_unfold_Emin_150.000000_all_allJetsCalibrated.root"); 
  hData_overcal = (TH1D*)_file_overcal->Get("hCastorJet_energy"); 
  hData_overcal->Scale( 1./ lumi_ ); 

  
  //== MC sample.
  TFile* _MC_overcal = TFile::Open( TString::Format("/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo" + cutsetup_ + "_unfold_Emin_150.000000_deltaPhiMax_%f_etaband_%f_all_matchE_allJetsCalibrated.root", phimax_global_, etaacc_global_) );  
  				   
  RooUnfoldResponse *response = (RooUnfoldResponse*)_MC_overcal->Get("response");
  
  //== Unfolding.
  RooUnfoldBayes unfold_bayes(response, hData_overcal, iterations);
  hData_overcal = (TH1D*) unfold_bayes.Hreco( RooUnfold::kCovariance ); 
  SetDnDx( hData_overcal );    
  hData_overcal->SetTitle("JetID");
  
  TH1D* overcal_copy = (TH1D*)hData_overcal->Clone("JetID");
  vector<TH1D*> vSyst;
  vSyst.push_back( overcal_copy );    
  
  //== Determine difference between the two.
  
  TCanvas *can_;
  PrepareCanvas( can_, "can_");
  
  hData->Draw("hist");
  hData_overcal->SetLineColor( kRed );
  hData_overcal->DrawClone("histsame");  
   
  hData_overcal->Add( hData, -1. ); 
  hData_overcal->Divide( hData );
  
  hData_overcal->SetLineColor( kBlue );
  hData_overcal->DrawClone("histsame");  
     
  can_->SetLogy();  
  can_->SaveAs( folder_ + "Overcalibration.pdf");    

  int n = hData->GetNbinsX();
  for(int bin = 0; bin < n; bin++){

    int bin_shift = 1;
    if (bin == n-1){ bin_shift = 0; }	//== This ensures the last ratio is drawn neatly.
    if (bin == 2){ bin_shift = 2; }	//== This ensures the last ratio is drawn neatly.

    double lowedge = hReference->GetXaxis()->GetBinLowEdge( bin+1 );
    double upedge = hReference->GetXaxis()->GetBinUpEdge( bin+1 );
    
    x[2*bin] = hReference->GetXaxis()->GetBinLowEdge( bin+1 );
    y[2*bin] = hData_overcal->GetBinContent( bin+bin_shift );
    x[2*bin+1] = hReference->GetXaxis()->GetBinUpEdge( bin+1 );
    y[2*bin+1] = hData_overcal->GetBinContent( bin+bin_shift );
    
    if( lowedge <= 300. && upedge >= 300. ){ uncertainties << "\t" << hData_overcal->GetBinContent( bin+bin_shift ) << endl; }
    if( lowedge <= 1100. && upedge >= 1000. ){ uncertainties << "\t" << hData_overcal->GetBinContent( bin+bin_shift ) << endl; }
    if( lowedge <= 1600. && upedge >= 1600. ){ uncertainties << "\t" << hData_overcal->GetBinContent( bin+bin_shift ) << endl; }     

  }//== End loop over bins
  cout << "=== Jet ID DONE===" << endl;
  
  Plot_Systematics( hData, vSyst, "JetID");    
}



void Unfolder::Calculate_JES_systematic( int iterations,  double *x_up, double *y_up, double *x_down, double *y_down){
  cout << "=== JES ===" << endl;
  
  ofstream uncertainties;
  uncertainties.open("Uncertainties.txt", ios::app);
  uncertainties << "\n\nJES" <<  endl; 

  //== Unfold data as usual for reference.
    
  TH1D* hData; 
  Get_DetEnergy( -1 , hData);
  hData->Scale( 1./ lumi_ );
  Get_DetUnfolded( 0 , hData, iterations );      
  TH1D* hReference = (TH1D*)hData->Clone("hReference");
  SetDnDx( hReference ); 	// Divide by binwidth.  
  
  TCanvas* can_;
  PrepareCanvas( can_, "JES_systematic" );
  TPad* pad_abs_, *pad_ratio_;
  SplitCanvas( can_, pad_abs_, pad_ratio_ );
  pad_abs_->cd();
  hReference->Draw("hist"); 
  
  TH1D* up_clone, *down_clone;
  
  vector<TH1D*> vSyst;
  
  //== JES: +15%
  TH1D *hJESup, *hJESdown;
  Get_DetEnergy_JESup( -1 , hJESup);	// Extract data distribution.
  hJESup->Scale( 1./ lumi_ );		// Scale to xsection.
  Get_DetUnfolded( 0 , hJESup, iterations, "all");
  SetDnDx( hJESup ); 			// Divide by binwidth.
    up_clone = (TH1D*)hJESup->Clone("hJESup");
  hJESup->Add( hReference, -1);		// Difference with unfolded data.
  hJESup->Divide( hReference );		// dN/N

  //== JES: -15%
  Get_DetEnergy_JESdown( -1 , hJESdown);// Extract data distribution.
  hJESdown->Scale( 1./ lumi_ );	// Scale to xsection.
  Get_DetUnfolded( 0 , hJESdown, iterations, "all");
  SetDnDx( hJESdown ); 			// Divide by binwidth.
    down_clone = (TH1D*)hJESdown->Clone("hJESdown");
  hJESdown->Add( hReference, -1);		// Difference with unfolded data.
  hJESdown->Scale( -1 );			// Difference with unfolded data.  
  hJESdown->Divide( hReference );		// dN/N
  
  TFile *_fSystematics = TFile::Open("_fSystematics_7TeV.root", "Update");
  
  up_clone->SetName("hJESup");
  up_clone->Write();
  down_clone->SetName("hJESdown");
  down_clone->Write();

  _fSystematics->Close();  

  int n = hData->GetNbinsX();
  for(int bin = 0; bin < n; bin++){

    int bin_shift = 1;
    if (bin == n-1){ bin_shift = 0; }	//== This ensures the last ratio is drawn neatly.
    if (bin == 2){ bin_shift = 2; }	//== This ensures the last ratio is drawn neatly.
    
    double lowedge = hReference->GetXaxis()->GetBinLowEdge( bin+1 );
    double upedge = hReference->GetXaxis()->GetBinUpEdge( bin+1 );    

    x_up[2*bin] = hReference->GetXaxis()->GetBinLowEdge( bin+1 );
    y_up[2*bin] = hJESup->GetBinContent( bin+bin_shift );
    x_up[2*bin+1] = hReference->GetXaxis()->GetBinUpEdge( bin+1 );
    y_up[2*bin+1] = hJESup->GetBinContent( bin+bin_shift );

    x_down[2*bin] = hReference->GetXaxis()->GetBinLowEdge( bin+1 );
    y_down[2*bin] = hJESdown->GetBinContent( bin+bin_shift );
    x_down[2*bin+1] = hReference->GetXaxis()->GetBinUpEdge( bin+1 );
    y_down[2*bin+1] = hJESdown->GetBinContent( bin+bin_shift );

    if( lowedge <= 300. && upedge >= 300. ){ uncertainties << "\t+" << hJESup->GetBinContent( bin+bin_shift ) << "\t-" <<  hJESdown->GetBinContent( bin+bin_shift ) << endl;}
    if( lowedge <= 1100. && upedge >= 1000. ){ uncertainties << "\t+" << hJESup->GetBinContent( bin+bin_shift ) << "\t-" <<  hJESdown->GetBinContent( bin+bin_shift ) << endl; }
    if( lowedge <= 1600. && upedge >= 1600. ){ uncertainties << "\t+" << hJESup->GetBinContent( bin+bin_shift ) << "\t-" << hJESdown->GetBinContent( bin+bin_shift ) << endl; }     
    
    
  }//== End loop over bins

  //== Drawing.  
  up_clone->SetLineColor( kRed );
  up_clone->SetLineStyle( 2 );
  up_clone->SetLineWidth( 2 );
  pad_abs_->cd();
  up_clone->DrawClone("histsame");
  pad_ratio_->cd();
  up_clone->Divide( hReference );
  up_clone->Draw("hist");
  
  down_clone->SetLineColor( kBlue );
  down_clone->SetLineStyle( 3 );
  down_clone->SetLineWidth( 2 );
  pad_abs_->cd();
  down_clone->DrawClone("histsame");
  pad_ratio_->cd();
  down_clone->Divide( hReference );
  down_clone->Draw("histsame");  
  
  TH1D*  hFirst = (TH1D*)hReference->Clone();
  hFirst->Divide(hReference) ;
  hFirst->DrawClone("histsame");
  
  pad_abs_->SetLogy();  
  
  can_->SaveAs( TString::Format( folder_ + "/Can_systematics_JES.pdf" ) );
  cout << TString::Format( folder_ + "/Can_systematics_JES.pdf" ) << endl;
  cout << "=== JES DONE===" << endl;
}


void Unfolder::Calculate_model_systematic( int iterations, double *x, double *y){
  cout << "=== MODEL ===" << endl;
  
  ofstream uncertainties;
  uncertainties.open("Uncertainties.txt", ios::app);
  uncertainties << "\n\nModel" <<  endl; 
    
  TH1D* hData; 
  Get_DetEnergy( -1 , hData);
  hData->Scale( 1./ lumi_ );
  Get_DetUnfolded( 0 , hData, iterations );      
  TH1D* hReference = (TH1D*)hData->Clone("hReference");
  SetDnDx( hReference ); 	// Divide by binwidth.  
  
  //== Loop over models.
  // 1. calculate unfolded spectrum per model.
  // 2. Loop over bins and calculate average.
  float models = 0.;
  std::vector<TH1D*> vModels;
  TH1D* hAverage;
      
  TFile* _file = TFile::Open("Unfolded_7TeV.root","Recreate");  
  
  //-- Calculate average of model dependence.
  for(int MC_ = 0; MC_ < MC_files_.size() ; MC_++){
    //-- Generator level.
    if( (set_of_tags_["mc_type"])[MC_files_[MC_]] != "model" && (set_of_tags_["mc_type"])[MC_files_[MC_]] != "actual"){      cout << "=== AFILE ===\t" << MC_files_[MC_] << endl; continue; }
    models++;
    
    cout << "=== MFILE ===\t" << MC_files_[MC_] << endl;
    
    //-- Data unfolded.
    TH1D* hModel;

    //-- Extract and properly scale data.
    Get_DetEnergy( -1 , hModel);
    hModel->Scale( 1./ lumi_ );	
    Get_DetUnfolded( MC_ , hModel, iterations, "all");    
    
    //== Data unfolded.
    TH1D* hist = (TH1D*)hModel->Clone( printLabel_[ MC_files_[MC_] ] );
    SetDnDx( hModel );
    _file->cd();
    hist->Write();
        		
    //== MC.
    TH1D* hist_MC;
    Get_PureGenEnergy( MC_, hist_MC );        		
    SetDnDx( hist_MC );
    hist_MC = (TH1D*)hist_MC->Clone( TString::Format(printLabel_[ MC_files_[MC_] ] + "_gen" ) );
    _file->cd();
    hist_MC->Write();

    //== Add the models to the summing histogram.
    if( models == 1 ){ 
      hAverage = (TH1D*)hModel->Clone("Average");
    }
    else{
      hAverage->Add( hModel );
    }
    
    if( (MC_files_[MC_]).Contains("ythia6") ) {hModel->SetTitle("Pythia6"); }
    if( (MC_files_[MC_]).Contains("ythia8") ) {hModel->SetTitle("Pythia8"); }
    if( (MC_files_[MC_]).Contains("EPOS") ) {hModel->SetTitle("EPOS"); }
    vModels.push_back( hModel );    
  }
  _file->Write();

  hAverage->Scale( 1./models );
  
  //== over all bins and check which model returns the biggest difference with the average (above and below). 
  TH1D* hModel_dep_high = (TH1D*)hAverage->Clone("hModel_dependence_up");
  TH1D* hModel_dep_low = (TH1D*)hAverage->Clone("hModel_dependence_down");
  
  TH1D* hSyst_up = (TH1D*)hAverage->Clone("hModel_dependence_up");
  TH1D* hSyst_down = (TH1D*)hAverage->Clone("hModel_dependence_down");  

  //== We loop over each bin.
  for(int bin = 0; bin <= hAverage->GetNbinsX(); bin++){
    double current_bin = hAverage->GetBinContent( bin );
    double bin_min = current_bin, bin_max = current_bin;

    //== We loop over each model in this bin and check which models has the highest and lowest value.
    for( int model = 0; model < vModels.size(); model++){
      TH1D* hMod = vModels[model];
      double model_bin = hMod->GetBinContent( bin );
      if( model_bin > bin_max ){ bin_max = model_bin; }
      if( model_bin < bin_min ){ bin_min = model_bin; }
    } 

    //== Store the difference between the highest/lowest distributions and the average as the model uncertainty in ratio!
    //== Calculate the difference between the furthest model and the average, and store the difference relative to the average.

    if( current_bin != 0.){
      hModel_dep_high->SetBinContent( bin, 	(bin_max - current_bin)/current_bin );
      hModel_dep_low->SetBinContent( bin, 	(current_bin - bin_min)/current_bin );
      
      hSyst_up->SetBinContent( bin, 	( 1. + (bin_max - current_bin)/current_bin ) * hReference->GetBinContent(bin) );
      hSyst_down->SetBinContent( bin, 	( 1. - (current_bin - bin_min)/current_bin ) * hReference->GetBinContent(bin) );
    }
    else{	//== If bin content is zero, uncertainty is zero.
      hModel_dep_high->SetBinContent( bin, bin_max );
      hModel_dep_low->SetBinContent( bin, bin_min );
    }
  }

  TFile *_fSystematics = TFile::Open("_fSystematics_7TeV.root", "Update");
  
  hSyst_up->Write();
  hSyst_down->Write();

  _fSystematics->Close();

  int n = hData->GetNbinsX();
  for(int bin = 0; bin < n; bin++){

    int bin_shift = 1;
    if (bin == n-1){ bin_shift = 0; }	//== This ensures the first ratio is drawn neatly.
    if (bin == 2){ bin_shift = 2; }	//== This ensures the last ratio is drawn neatly.

    double lowedge = hReference->GetXaxis()->GetBinLowEdge( bin+1 );
    double upedge = hReference->GetXaxis()->GetBinUpEdge( bin+1 );

    x[2*bin] = hReference->GetXaxis()->GetBinLowEdge( bin+1 );
    y[2*bin] = hModel_dep_high->GetBinContent( bin+bin_shift );
    x[2*bin+1] = hReference->GetXaxis()->GetBinUpEdge( bin+1 );
    y[2*bin+1] = hModel_dep_high->GetBinContent( bin+bin_shift );
    
    if( lowedge <= 300. && upedge >= 300. ){ uncertainties << "\t" << hModel_dep_high->GetBinContent( bin+bin_shift ) << endl; }
    if( lowedge <= 1100. && upedge >= 1000. ){ uncertainties << "\t" << hModel_dep_high->GetBinContent( bin+bin_shift ) << endl; }
    if( lowedge <= 1600. && upedge >= 1600. ){ uncertainties << "\t" << hModel_dep_high->GetBinContent( bin+bin_shift ) << endl; }     

  }//== End loop over bins
  cout << "=== MODEL DONE===" << endl;
  
  Plot_Systematics( hReference, vModels, "Model");  
}





void Unfolder::Calculate_position_systematic( int iterations, double *x_up, double *y_up, double *x_down, double *y_down ){
  cout << "=== POSITION ===" << endl;
  //== Unfold data as usual for reference.
  
  ofstream uncertainties;
  uncertainties.open("Uncertainties.txt", ios::app);
  uncertainties << "\n\nPosition" <<  endl; 
  
  vector<TH1D*> hSyst;
    
  TH1D* hData; 
  Get_DetEnergy( -1 , hData);
  hData->Scale( 1./ lumi_ );
  Get_DetUnfolded( 0 , hData, iterations );      
  TH1D* hReference = (TH1D*)hData->Clone("hReference");
  SetDnDx( hReference ); 	// Divide by binwidth.  
  
  TH1D* hPosUp, *hPosDown;  
  TH1D* hSyst_up, *hSyst_down;    
  
  //-- Position dependence.
  for(int MC_ = 0; MC_ < MC_files_.size() ; MC_++){
    //-- Generator level.
    if( (set_of_tags_["mc_type"])[MC_files_[MC_]] != "position" && (set_of_tags_["mc_type"])[MC_files_[MC_]] != "actual"){ continue; }
    //-- Generator level.

    //-- Data unfolded.
    TH1D* hPos;

    Get_DetEnergy( -1 , hPos); 
    hPos->Scale( 1./lumi_ );
    Get_DetUnfolded( MC_ , hPos, iterations, "all");
    SetDnDx( hPos );
   
    TH1D* hPos_copy = (TH1D*)hPos->Clone("copy");

    if( (MC_files_[MC_]).Contains("up") ) { 
      hPos_copy->SetTitle("PositionUp");
      hSyst.push_back( hPos_copy );
      hSyst_up = (TH1D*)hPos_copy->Clone("hPos_up");
          
      hPos->Add( hReference, -1. );
      hPos->Divide( hReference );
      hPosUp = (TH1D*)hPos->Clone("hPosUp");

    }
    else if( (MC_files_[MC_]).Contains("down") ){
      hPos_copy->SetTitle("PositionDown");
      hSyst.push_back( hPos_copy );
      hSyst_down = (TH1D*)hPos_copy->Clone("hPos_down");
    
      hPos->Add( hReference, -1. );
      hPos->Scale( -1. );
      hPos->Divide( hReference );
      hPosDown = (TH1D*)hPos->Clone("hPosDown");

    }
  } //== Loop over MC.
  
  TFile *_fSystematics = TFile::Open("_fSystematics_7TeV.root", "Update");
  
  hSyst_up->Write();
  hSyst_down->Write();

  _fSystematics->Close();

  int n = hData->GetNbinsX();
  for(int bin = 0; bin < n; bin++){

    int bin_shift = 1;
    if (bin == n-1){ bin_shift = 0; }	//== This ensures the last ratio is drawn neatly.
    if (bin == 2){ bin_shift = 2; }	//== This ensures the last ratio is drawn neatly.
    
    double shiftup = hPosUp->GetBinContent( bin+bin_shift );
    double shiftdown = hPosDown->GetBinContent( bin+bin_shift );  
    double posshift = max( shiftup, shiftdown);  
    double lowedge = hReference->GetXaxis()->GetBinLowEdge( bin+1 );
    double upedge = hReference->GetXaxis()->GetBinUpEdge( bin+1 );

    x_up[2*bin] = lowedge;
    y_up[2*bin] = posshift;
    x_up[2*bin+1] = upedge;
    y_up[2*bin+1] = posshift;

    x_down[2*bin] = lowedge;
    y_down[2*bin] = posshift;
    x_down[2*bin+1] = upedge;
    y_down[2*bin+1] = posshift;
    
    if( lowedge <= 300. && upedge >= 300. ){ uncertainties << "\t" << posshift << endl; }
    if( lowedge <= 1100. && upedge >= 1000. ){ uncertainties << "\t" << posshift << endl; }
    if( lowedge <= 1600. && upedge >= 1600. ){ uncertainties << "\t" << posshift << endl; }        
    
  }//== End loop over bins
  cout << "=== POSITION DONE===" << endl;
  
  uncertainties.close();
  
  Plot_Systematics( hReference, hSyst, "Position");
}





void Unfolder::Calculate_lumi_systematic( int iterations, double *x, double *y){
  cout << "=== LUMI ===" << endl;
  //== Unfold data as usual for reference.
    
  TH1D* hData; 
  Get_DetEnergy( -1 , hData);		cout << "=== GOT DET ===" << endl;
 

  hData->Scale( 1./ lumi_ );		cout << "=== SCALED ===" << endl;
  Get_DetUnfolded( 0 , hData, iterations );      	cout << "=== UNFOLDED ===" << endl;
  TH1D* hReference = (TH1D*)hData->Clone("hReference");	cout << "=== REFERENCE ===" << endl;
  SetDnDx( hReference ); 	cout << "=== DNDX ===" << endl; // Divide by binwidth.  
  
  
  TH1D* lumi_up = (TH1D*)hReference->Clone("LumiUp");
  lumi_up->SetTitle("LumiUp");
  lumi_up->Scale( 1.036 );
  TH1D* lumi_down = (TH1D*)hReference->Clone("LumiDown");
  lumi_down->SetTitle("LumiDown");
  lumi_down->Scale( 0.964 );
  vector<TH1D*> hSyst;  
  hSyst.push_back( lumi_up );
  hSyst.push_back( lumi_down );
  
  int n = hData->GetNbinsX();
  for(int bin = 0; bin < n; bin++){

    int bin_shift = 1;
    if (bin == n-1){ bin_shift = 0; }	//== This ensures the last ratio is drawn neatly.
    if (bin == 2){ bin_shift = 2; }	//== This ensures the last ratio is drawn neatly.

    x[2*bin] = hReference->GetXaxis()->GetBinLowEdge( bin+1 );
    y[2*bin] = 0.036;
    x[2*bin+1] = hReference->GetXaxis()->GetBinUpEdge( bin+1 );
    y[2*bin+1] = 0.036;
  }//== End loop over bins 
  
  cout << "=== LUMI DONE ===" << endl;
  
  Plot_Systematics( hReference, hSyst, "Luminosity"); 
  
  cout << "=== SYST DONE ===" << endl;
}





void Unfolder::PlotSystematics( int iterations ){
  
  //== Data.  
  TH1D* hData; 
  Get_DetEnergy( -1 , hData);
  hData->Scale( 1./ lumi_ );
  Get_DetUnfolded( 0 , hData, iterations );      
  TH1D* hReference = (TH1D*)hData->Clone("hReference");
  SetDnDx( hReference ); 	// Divide by binwidth.  
  TH1D* hFirst = (TH1D*)hReference->Clone("First");
  
  //== Luminosity up.
  TH1D* hLumi_up = (TH1D*)hReference->Clone("LumiUp");
  hLumi_up->Scale(1.036);
  
  //== Luminosity down.
  TH1D* hLumi_down = (TH1D*)hReference->Clone("LumiDown");
  hLumi_down->Scale(0.964);
  
  //== Position dependence.
  TH1D* hPosUp, *hPosDown;    
  for(int MC_ = 0; MC_ < MC_files_.size() ; MC_++){
    if( (set_of_tags_["mc_type"])[MC_files_[MC_]] != "position" && (set_of_tags_["mc_type"])[MC_files_[MC_]] != "actual"){ continue; }

    //-- Data unfolded.
    TH1D* hPos;

    Get_DetEnergy( -1 , hPos); 
    hPos->Scale( 1./lumi_ );
    Get_DetUnfolded( MC_ , hPos, iterations, "all");
    SetDnDx( hPos );

    if( (MC_files_[MC_]).Contains("up") ) { 
      hPosUp = (TH1D*)hPos->Clone("hPosUp");
    }
    else if( (MC_files_[MC_]).Contains("down") ){

      hPosDown = (TH1D*)hPos->Clone("hPosDown");
    }
  } //== Loop over MC.


  //== Model dependence.
  // 1. calculate unfolded spectrum per model.
  // 2. Loop over bins and calculate average.
  float models = 0.;
  std::vector<TH1D*> vModels;
  TH1D* hAverage;
  //-- Calculate average of model dependence.
  for(int MC_ = 0; MC_ < MC_files_.size() ; MC_++){
    //-- Generator level.
    if( (set_of_tags_["mc_type"])[MC_files_[MC_]] != "model" && (set_of_tags_["mc_type"])[MC_files_[MC_]] != "actual"){ continue; }
    models++;
    
    //-- Data unfolded.
    TH1D* hModel;

    //-- Extract and properly scale data.
    Get_DetEnergy( -1 , hModel);
    hModel->Scale( 1./ lumi_ );
    Get_DetUnfolded( MC_ , hModel, iterations, "all");
    SetDnDx( hModel );	

    //== Add the models to the summing histogram.
    if( MC_ == 0 ){ 
      hAverage = (TH1D*)hModel->Clone("Average");
    }
    else{
      hAverage->Add( hModel );
    }
    vModels.push_back( hModel );
  }

  hAverage->Scale( 1./models );

  //==========================================//
  //== Unfold data with all jets calibrated.==//
  //==========================================//
  TH1D* hData_overcal;
  TFile* _file_overcal = TFile::Open("/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5_data" + cutsetup_ + "_unfold_Emin_150.000000_all_allJetsCalibrated.root"); 
  hData_overcal = (TH1D*)_file_overcal->Get("hCastorJet_energy"); 
  hData_overcal->Scale( 1./ lumi_ );  
  
  //== MC sample.
  TFile* _MC_overcal = TFile::Open( TString::Format("/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_NewGeo" + cutsetup_ + "_unfold_Emin_150.000000_deltaPhiMax_%f_etaband_%f_all_matchE_allJetsCalibrated.root", phimax_global_, etaacc_global_ ), "READ"  );
  				     				   
  RooUnfoldResponse *response = (RooUnfoldResponse*)_MC_overcal->Get("response");
  RooUnfoldBayes unfold_bayes(response, hData_overcal, iterations);
  hData_overcal = (TH1D*) unfold_bayes.Hreco( RooUnfold::kCovariance ); 
  SetDnDx( hData_overcal );    
  
  
  //=========//
  //== JES.==//
  //=========//
  TH1D* up_clone, *down_clone;
  
  //== JES: +15%
  TH1D *hJESup, *hJESdown;
  Get_DetEnergy_JESup( -1 , hJESup);	// Extract data distribution.
  hJESup->Scale( 1./ lumi_ );		// Scale to xsection.
  Get_DetUnfolded( 0 , hJESup, iterations, "all");
  SetDnDx( hJESup ); 		

  //== JES: -15%
  Get_DetEnergy_JESdown( -1 , hJESdown);// Extract data distribution.
  hJESdown->Scale( 1./ lumi_ );	// Scale to xsection.
  Get_DetUnfolded( 0 , hJESdown, iterations, "all");
  SetDnDx( hJESdown ); 			// Divide by binwidth.
    down_clone = (TH1D*)hJESdown->Clone();
   
  //== Canvas.
  TCanvas *can_;
  PrepareCanvas( can_, "Systematics");
  TPad *pad_abs_, *pad_ratio_;
  SplitCanvas( can_, pad_abs_, pad_ratio_ );

  //== TLegend.
  TLegend* leg = new TLegend(
  	0.75, 0.5, 
  	1. - pad_abs_->GetRightMargin() - hData->GetYaxis()->GetTickLength(),
  	1. - pad_abs_->GetTopMargin() - hData->GetXaxis()->GetTickLength() );
  leg->SetBorderSize( 0 );
  leg->SetFillStyle( 0 );  
  
  //== Reference.
  pad_abs_->cd();
  hReference->GetXaxis()->SetRangeUser(300., 1600.);
  hReference->SetLineWidth( 2 );
  hReference->GetYaxis()->SetTitle("dN/dE [1/GeV]");
  Prepare_1Dplot( hReference );
  hReference->DrawClone("hist");
  pad_ratio_->cd();
  hReference->Divide( hFirst );
  Prepare_1Dplot_ratio( hReference );
  hReference->GetXaxis()->SetTitle("E [GeV]");
  hReference->GetYaxis()->SetTitle("Ratio");
  hReference->GetYaxis()->SetRangeUser(0.,2.2);  
  hReference->GetYaxis()->CenterTitle();
  hReference->GetYaxis()->SetNdivisions(504);    
  hReference->Draw("hist");
  leg->AddEntry( hReference, "Data", "l");
  
  int ci = 2;
  int ls = 2; 
  int marker = 33; 
  //== Lumi up.
  pad_abs_->cd();
  hLumi_up->SetLineColor( getColor( ci ) );
  hLumi_up->SetLineStyle( ls );
  hLumi_up->SetLineWidth( 2 );
  hLumi_up->DrawClone("histsame");
  pad_ratio_->cd();
  hLumi_up->Divide( hFirst );
  hLumi_up->Draw("histsame");  
//  leg->AddEntry( hLumi_up, "Lumi up", "l");
  
  //== Lumi down.
  pad_abs_->cd();
  hLumi_down->SetLineColor( getColor( ci++ ) );
  hLumi_down->SetLineStyle( ls++ );
  hLumi_down->SetLineWidth( 2 );
  hLumi_down->DrawClone("histsame");
  pad_ratio_->cd();
  hLumi_down->Divide( hFirst );
  hLumi_down->Draw("histsame"); 
  leg->AddEntry( hLumi_down, "Lumi (+/-)", "l"); 
  
  //== Pos up.
  pad_abs_->cd();
  hPosUp->SetLineColor( getColor( ci ) );
  hPosUp->SetLineStyle( ls );
  hPosUp->SetLineWidth( 2 );
  hPosUp->SetMarkerSize(1.2);
  hPosUp->SetMarkerColor( getColor( ci ) );  
  hPosUp->SetMarkerStyle( marker++ );
  hPosUp->DrawClone("plsame");
  pad_ratio_->cd();
  hPosUp->Divide( hFirst );
  hPosUp->Draw("plsame");
  leg->AddEntry( hPosUp, "Pos. up", "pl"); 
  
  //== Pos down.
  pad_abs_->cd();
  hPosDown->SetMarkerSize(1.2);
  hPosDown->SetMarkerColor( getColor( ci ) );  
  hPosDown->SetMarkerStyle( marker++ );
  hPosDown->SetLineColor( getColor( ci++ ) );
  hPosDown->SetLineStyle( ls++ );
  hPosDown->SetLineWidth( 2 );
  hPosDown->DrawClone("plsame");
  pad_ratio_->cd();
  hPosDown->Divide( hFirst );
  hPosDown->Draw("plsame"); 
  leg->AddEntry( hPosDown, "Pos. down", "pl");   
   
  //== Average model.
  pad_abs_->cd();
  hAverage->SetLineColor( getColor( ci++ ) );
  hAverage->SetLineStyle( ci++ );
  hAverage->SetLineWidth( 2 );
  hAverage->DrawClone("histsame");
  pad_ratio_->cd();
  hAverage->Divide( hFirst );
  hAverage->Draw("histsame"); 
  leg->AddEntry( hAverage, "Model", "l");  
  
  //== JetID.
  pad_abs_->cd();
  hData_overcal->SetLineColor( getColor( ci++ ) );
  hData_overcal->SetLineStyle( ci++ );
  hData_overcal->SetLineWidth( 2 );
  hData_overcal->DrawClone("histsame");
  pad_ratio_->cd();
  hData_overcal->Divide( hFirst );
  hData_overcal->Draw("histsame"); 
  leg->AddEntry( hData_overcal, "JetID", "l");    
     
  //== JES up.
  pad_abs_->cd();
  hJESup->SetLineColor( getColor( ci ) );
  hJESup->SetLineStyle( ls );
  hJESup->SetLineWidth( 2 );
  hJESup->DrawClone("histsame");
  pad_ratio_->cd();
  hJESup->Divide( hFirst );
  hJESup->Draw("histsame");  
  leg->AddEntry( hJESup, "JES+", "l");
  
  //== JES down.
  pad_abs_->cd();
  hJESdown->SetLineColor( getColor( ci ) );
  hJESdown->SetLineStyle( ls++ );
  hJESdown->SetLineWidth( 2 );
  hJESdown->DrawClone("histsame");
  pad_ratio_->cd();
  hJESdown->Divide( hFirst );
  hJESdown->Draw("histsame");  
  leg->AddEntry( hJESdown, "JES-", "l");  
  
  pad_abs_->cd();
  
  if( cutsetup_ == "_one_none_vertex_BSCor_HFor" ){ leg->SetHeader("New selection"); }
  if( cutsetup_ == "" ){ leg->SetHeader("Old selection"); }  
  
  leg->Draw();
  pad_abs_->SetLogy();
  can_->SaveAs( TString::Format( folder_ + "/Can_systematics" + label_ + ".pdf") );
}




void Unfolder::Convert_to_xf( TH1D* &hist ){

  vector<double> binedges;
  vector<double> bincontent;
  
  for(int bin = 1; bin <= hist->GetNbinsX() + 1; bin++ ){  
    binedges.push_back( hist->GetBinLowEdge( bin )/3500. );  
    bincontent.push_back( hist->GetBinContent( bin ) );
  }
  
  double* binedges_array = &binedges[0];
  TH1D* hnew = new TH1D("energy", "energy",  binedges.size()-1, binedges_array   );
  
  for(int bin = 1; bin <= hist->GetNbinsX() + 1; bin++ ){  
    hnew->SetBinContent( bin, hist->GetBinContent(bin) );
  }  
  
  hist = hnew;
}

