#include "Function_NonZeroMinimum.h"


#include "TF1.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TPad.h"

//== TH1D.

void Prepare_1Dplot(TH1D* &hist_, TPad* &pad_){

  double textsize = 32, charheight;
  double titleoffsetX = 1., titleoffsetY = 1.4;

  double pad_width  = pad_->XtoPixel(pad_->GetX2());
  double pad_height = pad_->YtoPixel(pad_->GetY1());
  charheight = textsize*pad_width;

  hist_->GetXaxis()->SetTitleFont(43);
  hist_->GetXaxis()->SetTitleSize(32);
  hist_->GetXaxis()->SetTitleOffset(1.);
  hist_->GetXaxis()->SetLabelFont(43);
  hist_->GetXaxis()->SetLabelSize(32);

  hist_->GetYaxis()->SetTitleFont(43);
  hist_->GetYaxis()->SetTitleSize(32);
  hist_->GetYaxis()->SetTitleOffset(1.35);
  hist_->GetYaxis()->SetLabelFont(43);
  hist_->GetYaxis()->SetLabelSize(32);

  double min = GetMinimumValue( hist_ );
  double max = hist_->GetMaximum();

  hist_->GetYaxis()->SetRangeUser( min*0.9, max*1.1);
}

void Prepare_1Dplot(TH1D* &hist_){

 

  hist_->GetXaxis()->SetTitleFont(43);
  hist_->GetXaxis()->SetTitleSize(32);
  hist_->GetXaxis()->SetTitleOffset(1.);
  hist_->GetXaxis()->SetLabelFont(43);
  hist_->GetXaxis()->SetLabelSize(32);

  hist_->GetYaxis()->SetTitleFont(43);
  hist_->GetYaxis()->SetTitleSize(32);
  hist_->GetYaxis()->SetTitleOffset(1.35);
  hist_->GetYaxis()->SetLabelFont(43);
  hist_->GetYaxis()->SetLabelSize(32);

  double min = GetMinimumValue( hist_ );
  double max = hist_->GetMaximum();

  hist_->GetYaxis()->SetRangeUser( min*0.9, max*1.1);
}

//== TGraph.

void Prepare_1Dplot(TGraph* &gr_){

  gr_->GetHistogram()->GetXaxis()->SetTitleFont(43);
  gr_->GetHistogram()->GetXaxis()->SetTitleSize(32);
  gr_->GetHistogram()->GetXaxis()->SetTitleOffset(1);
  gr_->GetHistogram()->GetXaxis()->SetLabelFont(43);
  gr_->GetHistogram()->GetXaxis()->SetLabelSize(29);

  gr_->GetHistogram()->GetYaxis()->SetTitleFont(43);
  gr_->GetHistogram()->GetYaxis()->SetTitleSize(32);
  gr_->GetHistogram()->GetYaxis()->SetTitleOffset(1.3);
  gr_->GetHistogram()->GetYaxis()->SetLabelFont(43);
  gr_->GetHistogram()->GetYaxis()->SetLabelSize(29);

  double min = GetMinimumValue( gr_->GetHistogram() );
  double max = gr_->GetHistogram()->GetMaximum();

  gr_->GetHistogram()->GetYaxis()->SetRangeUser( min*0.9, max*1.1);
}

//== TGraphErrors.

void Prepare_1Dplot(TGraphErrors* &gr_){

  gr_->GetHistogram()->GetXaxis()->SetTitleFont(43);
  gr_->GetHistogram()->GetXaxis()->SetTitleSize(32);
  gr_->GetHistogram()->GetXaxis()->SetTitleOffset(1);
  gr_->GetHistogram()->GetXaxis()->SetLabelFont(43);
  gr_->GetHistogram()->GetXaxis()->SetLabelSize(29);


  gr_->GetHistogram()->GetYaxis()->SetTitleFont(43);
  gr_->GetHistogram()->GetYaxis()->SetTitleSize(32);
  gr_->GetHistogram()->GetYaxis()->SetTitleOffset(1.3);
  gr_->GetHistogram()->GetYaxis()->SetLabelFont(43);
  gr_->GetHistogram()->GetYaxis()->SetLabelSize(29);

  double min = GetMinimumValue( gr_->GetHistogram() );
  double max = gr_->GetHistogram()->GetMaximum();

  gr_->GetHistogram()->GetYaxis()->SetRangeUser( min*0.9, max*1.1);
}

//== TGraphAsymmErrors.

void Prepare_1Dplot(TGraphAsymmErrors* &gr_, TPad* &pad_){

  double textsize = 32, charheight;
  double titleoffsetX = 1., titleoffsetY = 1.4;

  double pad_width  = pad_->XtoPixel(pad_->GetX2());
  double pad_height = pad_->YtoPixel(pad_->GetY1());
  charheight = textsize*pad_width;

  gr_->GetHistogram()->GetXaxis()->SetTitleFont(43);
  gr_->GetHistogram()->GetXaxis()->SetTitleSize(32);
  gr_->GetHistogram()->GetXaxis()->SetTitleOffset(1. * (756./pad_height) );
  gr_->GetHistogram()->GetXaxis()->SetLabelFont(43);
  gr_->GetHistogram()->GetXaxis()->SetLabelSize(32);

  gr_->GetHistogram()->GetYaxis()->SetTitleFont(43);
  gr_->GetHistogram()->GetYaxis()->SetTitleSize(32);
  gr_->GetHistogram()->GetYaxis()->SetTitleOffset(1.35);
  gr_->GetHistogram()->GetYaxis()->SetLabelFont(43);
  gr_->GetHistogram()->GetYaxis()->SetLabelSize(32);

  double min = GetMinimumValue( gr_->GetHistogram() );
  double max = gr_->GetHistogram()->GetMaximum();

//  gr_->GetHistogram()->GetYaxis()->SetRangeUser( min*0.9, max*1.1);
}

void Prepare_1Dplot(TGraphAsymmErrors* &gr_){

  gr_->GetHistogram()->GetXaxis()->SetTitleFont(43);
  gr_->GetHistogram()->GetXaxis()->SetTitleSize(32);
  gr_->GetHistogram()->GetXaxis()->SetTitleOffset(1. );
  gr_->GetHistogram()->GetXaxis()->SetLabelFont(43);
  gr_->GetHistogram()->GetXaxis()->SetLabelSize(29);


  gr_->GetHistogram()->GetYaxis()->SetTitleFont(43);
  gr_->GetHistogram()->GetYaxis()->SetTitleSize(32);
  gr_->GetHistogram()->GetYaxis()->SetTitleOffset(1.3);
  gr_->GetHistogram()->GetYaxis()->SetLabelFont(43);
  gr_->GetHistogram()->GetYaxis()->SetLabelSize(29);

  double min = GetMinimumValue( gr_->GetHistogram() );
  double max = gr_->GetHistogram()->GetMaximum();

  gr_->GetHistogram()->GetZaxis()->SetRangeUser( min*0.9, max*1.1);
}

//== TF1.

void Prepare_1Dplot(TF1* &hist_){

  hist_->GetXaxis()->SetTitleFont(43);
  hist_->GetXaxis()->SetTitleSize(20);
  hist_->GetXaxis()->SetTitleOffset(1);
  hist_->GetXaxis()->SetLabelFont(43);
  hist_->GetXaxis()->SetLabelSize(20);

  hist_->GetYaxis()->SetTitleFont(43);
  hist_->GetYaxis()->SetTitleSize(32);
  hist_->GetYaxis()->SetTitleOffset(1.5);
  hist_->GetYaxis()->SetLabelFont(43);
  hist_->GetYaxis()->SetLabelSize(29);


//  hist_->GetYaxis()->SetRangeUser( min*0.9, max*1.1);
}
