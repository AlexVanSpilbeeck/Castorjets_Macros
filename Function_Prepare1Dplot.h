#include "NonZeroMinimum.h"

#include "TH1.h"
#include "TGraph.h"
#include "TF1.h"

void Prepare_1Dplot(TH1D* &hist_){

  hist_->GetXaxis()->SetTitleFont(43);
  hist_->GetXaxis()->SetTitleSize(20);
  hist_->GetXaxis()->SetTitleOffset(2);
  hist_->GetXaxis()->SetLabelFont(43);
  hist_->GetXaxis()->SetLabelSize(20);

  hist_->GetYaxis()->SetTitleFont(43);
  hist_->GetYaxis()->SetTitleSize(35);
  hist_->GetYaxis()->SetTitleOffset(1.5);
  hist_->GetYaxis()->SetLabelFont(43);
  hist_->GetYaxis()->SetLabelSize(35);

  double min = GetMinimumValue( hist_ );
  double max = hist_->GetMaximum();

  hist_->GetZaxis()->SetRangeUser( min*0.9, max*1.1);
}

void Prepare_1Dplot(TGraph* &gr_){

  gr_->GetHistogram()->GetXaxis()->SetTitleFont(43);
  gr_->GetHistogram()->GetXaxis()->SetTitleSize(20);
  gr_->GetHistogram()->GetXaxis()->SetTitleOffset(2);
  gr_->GetHistogram()->GetXaxis()->SetLabelFont(43);
  gr_->GetHistogram()->GetXaxis()->SetLabelSize(20);

  gr_->GetHistogram()->GetYaxis()->SetTitleFont(43);
  gr_->GetHistogram()->GetYaxis()->SetTitleSize(35);
  gr_->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
  gr_->GetHistogram()->GetYaxis()->SetLabelFont(43);
  gr_->GetHistogram()->GetYaxis()->SetLabelSize(35);

  double min = GetMinimumValue( gr_->GetHistogram() );
  double max = gr_->GetHistogram()->GetMaximum();

  gr_->GetHistogram()->GetZaxis()->SetRangeUser( min*0.9, max*1.1);
}


void Prepare_1Dplot(TGraphErrors* &gr_){

  gr_->GetHistogram()->GetXaxis()->SetTitleFont(43);
  gr_->GetHistogram()->GetXaxis()->SetTitleSize(45);
  gr_->GetHistogram()->GetXaxis()->SetTitleOffset(2);
  gr_->GetHistogram()->GetXaxis()->SetLabelFont(43);
  gr_->GetHistogram()->GetXaxis()->SetLabelSize(45);

  gr_->GetHistogram()->GetYaxis()->SetTitleFont(43);
  gr_->GetHistogram()->GetYaxis()->SetTitleSize(35);
  gr_->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
  gr_->GetHistogram()->GetYaxis()->SetLabelFont(43);
  gr_->GetHistogram()->GetYaxis()->SetLabelSize(35);

  double min = GetMinimumValue( gr_->GetHistogram() );
  double max = gr_->GetHistogram()->GetMaximum();

  gr_->GetHistogram()->GetZaxis()->SetRangeUser( min*0.9, max*1.1);
}

void Prepare_1Dplot(TGraphAsymmErrors* &gr_){

  gr_->GetHistogram()->GetXaxis()->SetTitleFont(43);
  gr_->GetHistogram()->GetXaxis()->SetTitleSize(20);
  gr_->GetHistogram()->GetXaxis()->SetTitleOffset(2);
  gr_->GetHistogram()->GetXaxis()->SetLabelFont(43);
  gr_->GetHistogram()->GetXaxis()->SetLabelSize(20);

  gr_->GetHistogram()->GetYaxis()->SetTitleFont(43);
  gr_->GetHistogram()->GetYaxis()->SetTitleSize(35);
  gr_->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
  gr_->GetHistogram()->GetYaxis()->SetLabelFont(43);
  gr_->GetHistogram()->GetYaxis()->SetLabelSize(35);

  double min = GetMinimumValue( gr_->GetHistogram() );
  double max = gr_->GetHistogram()->GetMaximum();

  gr_->GetHistogram()->GetZaxis()->SetRangeUser( min*0.9, max*1.1);
}


void Prepare_1Dplot(TF1* &hist_){

  hist_->GetXaxis()->SetTitleFont(43);
  hist_->GetXaxis()->SetTitleSize(20);
  hist_->GetXaxis()->SetTitleOffset(2);
  hist_->GetXaxis()->SetLabelFont(43);
  hist_->GetXaxis()->SetLabelSize(20);

  hist_->GetYaxis()->SetTitleFont(43);
  hist_->GetYaxis()->SetTitleSize(35);
  hist_->GetYaxis()->SetTitleOffset(1.5);
  hist_->GetYaxis()->SetLabelFont(43);
  hist_->GetYaxis()->SetLabelSize(35);


//  hist_->GetYaxis()->SetRangeUser( min*0.9, max*1.1);
}
