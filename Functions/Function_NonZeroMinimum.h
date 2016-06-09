// Code by Alex Van Spilbeeck.

#ifndef NONZEROMINIMUM_H
#define NONZEROMINIMUM_H

#include <TH1.h>
#include <THStack.h>

#include <iostream>

double GetMinimumValue(TH1 * h){

  double minimum = h->GetMaximum();

  if( minimum > 0.){
    for(int bin = 1; bin <= h->GetNbinsX(); bin++){
      double value = h->GetBinContent(bin);
      if( value < minimum && value > 0.) minimum = value;
    }
  }
  return minimum;
} 

double GetMinimumValue(TH2D * h){

  double minimum = h->GetMaximum();

  if( minimum > 0.){
    for(int binX = 1; binX <= h->GetNbinsX(); binX++){
      for(int binY = 1; binY <= h->GetNbinsY(); binY++){
        double value = h->GetBinContent(binX, binY);
        if( value < minimum && value > 0.) minimum = value;
      }
    }
  }

  return minimum;
}

//== THStack

double GetMinimumValue(THStack * h){

  double minimum = h->GetHistogram()->GetMaximum();

  if( minimum > 0.){
    for(int bin = 1; bin <= h->GetHistogram()->GetNbinsX(); bin++){
      double value = h->GetHistogram()->GetBinContent(bin);
      if( value < minimum && value > 0.) minimum = value;
    }
  }
  return minimum;
} 



#endif
