#include "NonZeroMinimum.h"

#include "TH1.h"
#include "TGraph.h"


void Prepare_2Dplot(TH2D* &hist_){

  hist_->GetXaxis()->SetNdivisions(304);
  hist_->GetXaxis()->SetTitleFont(43);
  hist_->GetXaxis()->SetTitleSize(35);
  hist_->GetXaxis()->SetTitleOffset(1.5);
  hist_->GetXaxis()->SetLabelFont(43);
  hist_->GetXaxis()->SetLabelSize(30);

  hist_->GetYaxis()->SetNdivisions(304);
  hist_->GetYaxis()->SetTitleFont(43);
  hist_->GetYaxis()->SetTitleSize(35);
  hist_->GetYaxis()->SetTitleOffset(1.5);
  hist_->GetYaxis()->SetLabelFont(43);
  hist_->GetYaxis()->SetLabelSize(30);

  double min = GetMinimumValue( hist_ );
  double max = hist_->GetMaximum();

  hist_->GetZaxis()->SetRangeUser( min*0.9, max*1.1);
}




void Prepare_2Dplot(TH2D* &hist_, TString xaxis, TString yaxis, int divisionsx, int divisionsy){

  hist_->GetYaxis()->SetTitle( yaxis );

  hist_->GetXaxis()->SetTitle( xaxis );

  hist_->GetXaxis()->SetNdivisions(divisionsx);
  hist_->GetXaxis()->SetTitleFont(43);
  hist_->GetXaxis()->SetTitleSize(35);
  hist_->GetXaxis()->SetTitleOffset(1.);
  hist_->GetXaxis()->SetLabelFont(43);
  hist_->GetXaxis()->SetLabelSize(30);

  hist_->GetYaxis()->SetNdivisions(divisionsy);
  hist_->GetYaxis()->SetTitleFont(43);
  hist_->GetYaxis()->SetTitleSize(35);
  hist_->GetYaxis()->SetTitleOffset(1.5);
  hist_->GetYaxis()->SetLabelFont(43);
  hist_->GetYaxis()->SetLabelSize(30);


  double min = GetMinimumValue( hist_ );
  double max = hist_->GetMaximum();

  hist_->GetZaxis()->SetRangeUser( min*0.9, max*1.1);
}


