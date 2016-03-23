 //STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <TThread.h>
#include <TStopwatch.h>
#include <TLegend.h>
#include <TPad.h>
#include <TColor.h>
#include <TText.h>
#include <TLine.h>
#include <TPaletteAxis.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TGraph.h>

//STANDARD C++ INCLUDES
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <cstring>


double GetAverage(TH1D *histo){

  double sum = 0., count = 0.;
  
  for(int bin = 1; bin <= histo->GetNbinsX(); bin++){
  
    double bin_val = histo->GetBinCenter( bin );
    double value = histo->GetBinContent( bin );
    
    sum += bin_val * value;
    count += value;  
  }
  
  if( count != 0) sum = sum/count;

}
