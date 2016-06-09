#include "iostream"
#include "fstream"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TColor.h"
#include "TApplication.h"
#include "stdio.h"
#include "cstring"
#include "TMath.h"
#include <sys/stat.h>

#include "../tdrstyle.C"

#include "../Functions/Function_PrepareCanvas.h"

using namespace std; 

#define pi 3.14159
#define set_title 0

double Energy_Data_good( double edet ){
  double alpha = 2.38935;
  double beta = -0.0811714;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_MC_good( double edet ){
  double alpha = 2.06553;
  double beta = -0.0761115;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}


double Energy_Data_bad( double edet ){
  double alpha = 0.335274;
  double beta = 0.413631;
  double gamma = -80;
  return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_MC_bad( double edet ){
  double alpha = -0.452917;
  double beta = 0.443985;
  double gamma = -26.3024;
  return  edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}
void Prepare_1Dplot(TGraph* &gr_){

  gr_->GetHistogram()->GetXaxis()->SetTitleFont(43);
  gr_->GetHistogram()->GetXaxis()->SetTitleSize(32);
  gr_->GetHistogram()->GetXaxis()->SetTitleOffset(1);
  gr_->GetHistogram()->GetXaxis()->SetLabelFont(43);
  gr_->GetHistogram()->GetXaxis()->SetLabelSize(26);

  gr_->GetHistogram()->GetYaxis()->SetTitleFont(43);
  gr_->GetHistogram()->GetYaxis()->SetTitleSize(32);
  gr_->GetHistogram()->GetYaxis()->SetTitleOffset(1.3);
  gr_->GetHistogram()->GetYaxis()->SetLabelFont(43);
  gr_->GetHistogram()->GetYaxis()->SetLabelSize(26);

  double min = GetMinimumValue( gr_->GetHistogram() );
  double max = gr_->GetHistogram()->GetMaximum();

  gr_->GetHistogram()->GetYaxis()->SetRangeUser( min*0.9, max*1.1);
}

int main(int argc, char *argv[]) {

   setTDRStyle();

  int n = 100;
  double Edet[n], Ecal[n], Cal[n];

  void (*current_calibration)(double) = NULL;



  //=================
  //== Shower Library
  //=================

  //== Good sectors.
  for( int i = 0; i < n; i++){
    Edet[i] = 100. + i * 20.;
    Ecal[i] = Energy_Data_good( Edet[i] );
    Cal[i] = Ecal[i]/Edet[i];
  }
  TGraph* gr_SL_good = new TGraph( n, Edet, Ecal);
  TGraph* gr_SL_good_factor = new TGraph( n, Edet, Cal);

  //== Sector BAD.
  for( int i = 0; i < n; i++){
    Edet[i] = 100. + i * 20.;
    Ecal[i] = Energy_Data_bad( Edet[i] );
    Cal[i] = Ecal[i]/Edet[i];
  }
  TGraph* gr_SL_bad = new TGraph( n, Edet, Ecal);
  TGraph* gr_SL_bad_factor = new TGraph( n, Edet, Cal);

  //==================
  //== Full Simulation
  //==================

  //== Good sectors.
  for( int i = 0; i < n; i++){
    Edet[i] = 100. + i * 20.;
    Ecal[i] = Energy_MC_good( Edet[i] );
    Cal[i] = Ecal[i]/Edet[i];
  }
  TGraph* gr_FS_good = new TGraph( n, Edet, Ecal);
  TGraph* gr_FS_good_factor = new TGraph( n, Edet, Cal);

  //== Sector BAD.
  for( int i = 0; i < n; i++){
    Edet[i] = 100. + i * 20.;
    Ecal[i] = Energy_MC_bad( Edet[i] );
    Cal[i] = Ecal[i]/Edet[i];
  }
  TGraph* gr_FS_bad = new TGraph( n, Edet, Ecal);
  TGraph* gr_FS_bad_factor = new TGraph( n, Edet, Cal);

  //==============
  //== 1-to-1 line
  //==============

  for( int i = 0; i < n; i++){
    Edet[i] = 50. + n * 20.;
    Ecal[i] = Edet[i] ;
  }
  TGraph* one_to_one = new TGraph( n, Edet, Ecal);


  //========
  //== DRAW.
  //========

  TCanvas* can_ = new TCanvas("canvas", "canvas", 800, 756);
  PrepareCanvas( can_, "Calibrated_energy" );

  gr_SL_good->SetLineStyle( 1 );
   gr_SL_good->SetLineColor( kBlue );
   gr_SL_good->GetHistogram()->GetXaxis()->SetTitle("E_{det} [GeV]");
   gr_SL_good->GetHistogram()->GetYaxis()->SetTitle("E_{cal} [GeV]");
   Prepare_1Dplot( gr_SL_good );
   gr_SL_good->GetXaxis()->SetRangeUser( 0., 1500.);
   gr_SL_good->GetYaxis()->SetRangeUser( 0., 3500.);
   gr_SL_good->Draw("al");

  one_to_one->SetLineStyle( 3 );
    one_to_one->Draw("lsame");
    gr_SL_good->Draw("lsame");


  gr_SL_bad->SetLineStyle( 1 );
   gr_SL_bad->SetLineColor( kRed );
   gr_SL_bad->Draw("lsame");

  gr_FS_good->SetLineStyle( 2 );
   gr_FS_good->SetLineColor( kBlue );
   gr_FS_good->Draw("lsame");

  gr_FS_bad->SetLineStyle( 2 );
   gr_FS_bad->SetLineColor( kRed );
   gr_FS_bad->Draw("lsame");

  TLegend* leg = new TLegend(can_->GetLeftMargin(), 0.67, 0.60, 1. - can_->GetTopMargin() - 0.05 );
    leg->AddEntry( gr_SL_good, "Shower Library [good sectors]", "l");
    leg->AddEntry( gr_SL_bad, "Shower Library [bad sectors]", "l");


    leg->AddEntry( gr_FS_good, "Full Simulation [good sectors]", "l");
    leg->AddEntry( gr_FS_bad, "Full Simulation [bad sectors]", "l");

  leg->SetFillStyle( 0 );
  leg->SetBorderSize(0);

  leg->Draw();

  can_->SaveAs("Calibrated_energy.pdf");
  can_->SaveAs("Calibrated_energy.C");


  //== 
  //== Draw the calibration factors.
  //==

  gr_SL_good_factor->SetLineStyle( 1 );
   gr_SL_good_factor->SetLineColor( kBlue );
   gr_SL_good_factor->GetHistogram()->GetXaxis()->SetTitle("E_{det} [GeV]");
   gr_SL_good_factor->GetHistogram()->GetYaxis()->SetTitle("E_{cal}/E_{det} [GeV]");
   Prepare_1Dplot( gr_SL_good_factor );
   gr_SL_good_factor->GetXaxis()->SetRangeUser( 0., 1500.);
   gr_SL_good_factor->GetYaxis()->SetRangeUser( 0., 5.);
   gr_SL_good_factor->Draw("al");

  gr_SL_bad_factor->SetLineStyle( 1 );
   gr_SL_bad_factor->SetLineColor( kRed );
   gr_SL_bad_factor->Draw("lsame");

  gr_FS_good_factor->SetLineStyle( 2 );
   gr_FS_good_factor->SetLineColor( kBlue );
   gr_FS_good_factor->Draw("lsame");

  gr_FS_bad_factor->SetLineStyle( 2 );
   gr_FS_bad_factor->SetLineColor( kRed );
   gr_FS_bad_factor->Draw("lsame");

  leg = new TLegend(can_->GetLeftMargin() ,0.67, 0.60, 1. - can_->GetTopMargin() - 0.05 );
    leg->AddEntry( gr_SL_good_factor, "Shower Library [good sectors]", "l");
    leg->AddEntry( gr_SL_bad_factor, "Shower Library [bad sectors]", "l");


    leg->AddEntry( gr_FS_good_factor, "Full Simulation [good sectors]", "l");
    leg->AddEntry( gr_FS_bad_factor, "Full Simulation [bad sectors]", "l");

  leg->SetFillStyle( 0 );
  leg->SetBorderSize(0);

  leg->Draw();

  can_->SaveAs("Calibrated_factors.pdf");
  can_->SaveAs("Calibrated_factors.C");
}
