#include "Function_NonZeroMinimum.h"

#include "TCanvas.h"
#include "TString.h"

//#include "tdrstyle.C"

void PrepareCanvas( TCanvas* &can_, TString label){

  int W = 800;
  int H = 600;
  int H_ref = 600; 
  int W_ref = 800; 
  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;


/*
  can_ = new TCanvas( label, label, 1200, 1000 );
  can_->SetLeftMargin(0.12);
  can_->SetRightMargin(0.04);
  can_->SetTopMargin(0.12);
  can_->SetBottomMargin(0.12);  
*/
  can_ = new TCanvas( label, label, 50,50,W,H);
  can_->SetFillColor(0);
  can_->SetBorderMode(0);
  can_->SetFrameFillStyle(0);
  can_->SetFrameBorderMode(0);
  can_->SetLeftMargin( L/W );
  can_->SetRightMargin( R/W );
  can_->SetTopMargin( T/H );
  can_->SetBottomMargin( B/H );
  can_->SetTickx(0);
  can_->SetTicky(0);
}

void PrepareCanvas_2D( TCanvas* &can_, TString label){

  double leftmargin = 0.12, topmargin = 0.12, bottommargin = 0.12;
/*
  int can_ww = 1050, can_wh = 900;
  can_ = new TCanvas( label, label, can_ww, can_wh );

  double rightmargin = (1. - leftmargin) - (1.-topmargin-bottommargin)*double(can_wh)/double(can_ww);

   can_->SetLeftMargin(		leftmargin);
   can_->SetRightMargin(	rightmargin);
   can_->SetTopMargin(		topmargin);
   can_->SetBottomMargin(	bottommargin);  
*/

  int W = 1050;
  int H = 900;
  int H_ref = 900; 
  int W_ref = 1050; 
  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  //float R = 0.04*W_ref;
  float R = (1. - L/W) - (1.-T/H-B/H)*double(H)/double(W);
  R = R*W_ref;

  cout << "\n\tPrepare canvas\t(L= " << L/W << ", R = " << R/W << ", T = " << T/H << ", B = " << B/H << ")\t" << endl;

  can_ = new TCanvas( label, label, 50,50,W,H);
  can_->SetLeftMargin(	L/W);
  can_->SetRightMargin(	R/W);
  can_->SetTopMargin(	T/H);
  can_->SetBottomMargin(B/H);  

  cout << "\n\tPrepared canvas\t" << endl;

}

void SplitCanvas(TCanvas* &can_, TPad* &pad_abs_, TPad* &pad_ratio_, double left_ = 0.20,   double top_ = 0.05, double right_ = 0.05, double bottom_ = 0.25){

  can_->cd();

  can_->SetLeftMargin(0);
  can_->SetRightMargin(0);
  can_->SetTopMargin(0);
  can_->SetBottomMargin(0);

    

  pad_abs_ = new TPad("Absolute_Values", "Absolute_Values", 0., 0.4,1.,1. - can_->GetTopMargin());
  pad_abs_->SetLeftMargin(0.16);
  pad_abs_->SetRightMargin(0.02);
  pad_abs_->SetBottomMargin(0.);
  pad_abs_->SetTopMargin( 0.05 * 1/0.6 );
  pad_abs_->Draw();

  can_->cd();
  pad_ratio_ = new TPad("Ratios", "Ratios", 0., 0., 1., 0.4);
  pad_ratio_->SetLeftMargin(0.16);
  pad_ratio_->SetRightMargin(0.02);
  pad_ratio_->SetTopMargin(0);
  pad_ratio_->SetBottomMargin( 0.13 * 1./.4 );
  pad_ratio_->Draw();
  
}


