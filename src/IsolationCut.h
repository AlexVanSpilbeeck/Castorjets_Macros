#ifndef IsolationCut_h
#define IsolationCut_h

#include <TMath.h>
#include <iostream>

bool IsolationCut( double phi_jet, int nTowers, double phi_gen_part){

  double PI = 3.1415926536;
  
  double phi_width = 2.*PI/16.;
  
  bool remove_energy = false;

  // -- Only one sector: straightforward determination of phi range.

  if( nTowers == 1 ){
    if( fabs( phi_jet - phi_gen_part ) < phi_width ){ 
      remove_energy = true;
    }
  } // nTowers == 1.

  // -- Two sectors: determine the closest sector center.
  
  if( nTowers == 2){
    int det_sec = -10;
    double dphi_min = 10.;
  
    for(int sec = -8; sec < 8; sec++){
        
      double phi_sec = (sec + 0.5) * phi_width ;
            
      if( (phi_jet - phi_sec) < phi_width/2. ){
	det_sec = sec;
        break;
      }  
    } // nTowers == 2.
      
    double phi_sec = (det_sec + 0.5) * phi_width;
    
    // Initialize the allowed range for a jet below the phi-center of the sector.
    double sec_min = ( det_sec - 1. ) * phi_width;
    double sec_max = ( det_sec ) * phi_width;

    if( phi_sec < phi_jet ){ // Adjust if the jet lies above the phi-center.
       sec_min = (det_sec ) * phi_width;
       sec_max = (det_sec + 1. ) * phi_width;    
    }        

    // Our jet lies in the range? Bingo!
    if( phi_gen_part > sec_min && phi_gen_part < sec_max ){ remove_energy = true; }
    
    // Our jet lies on the jump from negative to positive?
    // -- Jet center lies in sector 7.
    else if( phi_jet > 7.5 * phi_width ){ // Lies in upper half of highest sector. Jet lies in sector 7 and -8.
      if( phi_gen_part < -7. * phi_width ){ remove_energy = true; }
      if( phi_gen_part > 7. * phi_width ) { remove_energy = true; }              
    }   
    else if( phi_jet < 7.5 * phi_width ){ // Lies in lower half of upper sector. Jet lies in sector 6 and 7.
      if( phi_gen_part > 6. * phi_width ) { remove_energy = true; }
    }

    // -- Jet center lies in sector -8.
    else if( phi_jet < -7.5 * phi_width ){// Jet lies in sector -8 and 7.
      if( phi_gen_part < -7. * phi_width ){ remove_energy = true; }
      if( phi_gen_part > 7. * phi_width ) { remove_energy = true; }
    }   
    else if( phi_jet < -7. * phi_width ){
      if( phi_gen_part < -6. * phi_width ){ remove_energy = true; }
      if( phi_gen_part > 7. * phi_width ) { remove_energy = true; }
    }
  } // nTowers == 2.


  // -- 3 Towers.
  if( nTowers == 3){
    int det_sec; 
    for(int sec = -8; sec < 8; sec++){
      double phi_sec = (sec + 0.5) * phi_width ;

      if( (phi_jet - phi_sec) < phi_width/2. ){
        det_sec = sec;
	break;
      } 
    }
    
    // Remove anything in the sector before and the sector after.
    if( phi_gen_part > (det_sec - 1) * phi_width && phi_gen_part < phi_width * (det_sec + 1) ){ remove_energy = true; }
    else if( det_sec == -8 && ( phi_gen_part < -6 * phi_width || phi_gen_part > 7 * phi_width ) ){ remove_energy = true; }
    else if( det_sec == 7 &&  ( phi_gen_part < -7 * phi_width || phi_gen_part > 6 * phi_width ) ){ remove_energy = true; }
  } // nTowers == 3.

 return remove_energy;
}

#endif
