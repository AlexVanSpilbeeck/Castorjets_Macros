#include "TH1D.h"
#include "TH2D.h"
/******************************************************************
* Divide the bin content of each bin in a histogram by bin width. *
*******************************************************************/

void SetDnDx(TH1D* &hist_){

  for(int bin = 1; bin <= hist_->GetNbinsX(); bin++){
    double bincontent = hist_->GetBinContent( bin );
    double binerror = hist_->GetBinError( bin );
    double binwidth = hist_->GetBinWidth( bin );
    hist_->SetBinContent( bin, bincontent/binwidth );
    hist_->SetBinError( bin, binerror/binwidth );
  }
}



void InvertDnDx(TH1D* &hist_){

  for(int bin = 1; bin <= hist_->GetNbinsX(); bin++){
    double bincontent = hist_->GetBinContent( bin );
    double binerror = hist_->GetBinError( bin );
    double binwidth = hist_->GetBinWidth( bin );
    hist_->SetBinContent( bin, bincontent*binwidth );
    hist_->SetBinError( bin, binerror*binwidth );
  }
}


void SetD2nDxDy(TH2D* &hist_){

  for(int binx = 1; binx <= hist_->GetNbinsX(); binx++){
    double binwidthx = hist_->GetXaxis()->GetBinWidth( binx );

    for(int biny = 1; biny <= hist_->GetNbinsY(); biny++){
      double binwidthy = hist_->GetYaxis()->GetBinWidth( biny );
 
      double bincontent = hist_->GetBinContent( binx, biny );    
      double binerror = hist_->GetBinError( binx, biny );    

      hist_->SetBinContent( binx, biny, bincontent/(binwidthx * binwidthy) );
      hist_->SetBinError( binx, biny, binerror/(binwidthx * binwidthy) );
      if( binerror != binerror ){ hist_->SetBinError( binx, biny, 0. ); }
    }
  }
}
