#include <TH1.h>
#include <vector>

#include "Function_NonZeroMinimum.h"

//== This function is used to set the y-range of a histogram so that it covers the ranges of all histograms in histlist.

void SetRangeToIncludeAll( TH1D* &hist_, vector<TH1D*> histlist ){

  if( histlist.size() == 0 ) return;

  vector<TH1D*>::iterator it = histlist.begin();
  double min = GetMinimumValue( *it );
  double max = (*it)->GetMaximum();

   for(; it != histlist.end(); ++it){
    double curr_min = GetMinimumValue( *it );
    if( curr_min < min ){ min = curr_min; }

    double curr_max = ( *it )->GetMaximum();
    if( curr_max > max ){ max = curr_max; }

   }

  hist_->GetYaxis()->SetRangeUser( 0.9 * min, 1.1 * max);
}
