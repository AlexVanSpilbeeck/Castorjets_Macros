 //STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>

//STANDARD C++ INCLUDES
#include <sys/stat.h>
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

void Empty_UOflow( TH1* hist ){

  hist->SetBinContent( 0, 0.);
  int nbins = hist->GetNbinsX();
  hist->SetBinContent( nbins+1, 0.);
  

}

