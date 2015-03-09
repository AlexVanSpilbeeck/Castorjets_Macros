//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>




void First_Plot(TH1D* &original, TH1D* &histo, int fileN, double &max_val, double& min_val){

  histo->Scale( 1./histo->Integral() );    

  if( fileN == 0 ){ 
    original = histo; 
    max_val = original->GetMaximum();
    min_val = original->GetMinimum();
  }
  
  if( histo->GetMaximum() > max_val) { max_val = histo->GetMaximum(); }
  
   if( histo->GetMinimum() < min_val) { min_val = histo->GetMinimum(); }
  
   if( min_val <= 0 && max_val <= 0)		{ original->GetYaxis()->SetRangeUser( min_val * 1.1, max_val * 0.9); }
   else if( min_val < 0 && max_val >= 0)	{ original->GetYaxis()->SetRangeUser( min_val * 1.1, max_val * 1.1); }
   else if( min_val >= 0 && max_val >= 0)	{ original->GetYaxis()->SetRangeUser( min_val * 0.9, max_val * 1.1); }
   
   
   
}
