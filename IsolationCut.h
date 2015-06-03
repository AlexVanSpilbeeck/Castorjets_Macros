#ifndef IsolationCut_h
#define IsolationCut_h

#include <TMath.h>
#include <iostream>

bool IsolationCut( double phi_jet, int nTowers, double phi_gen_part){

  double PI = 3.14;
  double phi_width = 16./(2.*PI);
  bool remove_energy = false;

  if( nTowers == 1 ){
    if( fabs( phi_jet - phi_gen_part ) < phi_width ){ 
      remove_energy = true;
    }
  } // nTowers == 1.


 return remove_energy;
}

#endif
