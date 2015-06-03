// Code by Alex Van Spilbeeck.

#ifndef NONZEROMINIMUM_H
#define NONZEROMINIMUM_H

#include <TH1.h>

#include <iostream>

double GetMinimumValue(TH1D * h){

  double minimum = h->GetMaximum();

  if( minimum > 0.){
    for(int bin = 0; bin < h->GetNbinsX(); bin++){
      double value = h->GetBinContent(bin);
      if( value < minimum && value > 0.) minimum = value;
    }
  }
  return minimum;
} 

double GetMinimumValue(TH2D * h){

  double minimum = h->GetMaximum();

  if( minimum > 0.){
    for(int binX = 0; binX < h->GetNbinsX(); binX++){
      for(int binY = 0; binY < h->GetNbinsY(); binY++){
        double value = h->GetBinContent(binX, binY);
        if( value < minimum && value > 0.) minimum = value;
      }
    }
  }

  return minimum;
}

#endif
