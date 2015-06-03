//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>
#include <TF1.h>
#include <vector>

#include "NonZeroMinimum.h"


void First_Plot(TH1D* &original, TH1D* &histo, int fileN, double &max_val, double& min_val){

  //if( histo->Integral() > 0. ) {  histo->Scale( 1./histo->Integral() ); }

  if( fileN == 0 ){ 
    original = histo; 
    max_val = original->GetMaximum();
    min_val = GetMinimumValue(histo);
  }
  
  if( histo->GetMaximum() > max_val) { 
    max_val = histo->GetMaximum(); 
  }
  
  if( GetMinimumValue(histo) <= min_val) { 
    min_val = GetMinimumValue(histo); 
  }
  
  else if( histo->GetMinimum() < 0. && histo->GetMinimum() <= min_val ){
    min_val = histo->GetMinimum();
  }
  
  if( min_val <= 0 && max_val <= 0)		{ original->GetYaxis()->SetRangeUser( min_val * 1.1, max_val * 0.9); }
  else if( min_val <= 0 && max_val > 0)	{ original->GetYaxis()->SetRangeUser( min_val * 1.1, max_val * 1.1);  }
  else if( min_val >= 0 && max_val >= 0)	{ original->GetYaxis()->SetRangeUser( min_val * 0.9, max_val * 1.1); }   

  std::cout << "Min and max\t" << min_val << "\t" << max_val << endl;
   
//  cout << "Y-range " << min_val << "\t" << max_val << "\tafter\t" << histo->GetName() << "\tfor\t" <<
//  original->GetName() << endl;
}

void Level_2D_plots( std::vector<TH2D*> &hist_vector, std::vector<TPad*> &pad_vector){

  double minval = 1.;
  double maxval = 0.;

  std::cout << "Hist\tminval\tmin_\tmaxval\tmax_" << std::endl;

  for(int h = 0; h < hist_vector.size(); h++){
    TH2D* &histo = hist_vector[h];
    
    double min_ = GetMinimumValue(histo);
    double max_ = histo->GetMaximum();

    std::cout << h << "\t" << minval << "\t" << min_ << "\t" << maxval << "\t" << max_ << std::endl;

    if( min_ < minval ) minval = min_;
    if( max_ > maxval ) maxval = max_;  
  }

  for(int h = 0; h < hist_vector.size(); h++){
    TH2D* &histo = hist_vector[h];
    histo->GetZaxis()->SetRangeUser( minval * 0.9, maxval * 1.1);
    TPad* &pad = pad_vector[h];
    pad->Update();
  }

  



}

