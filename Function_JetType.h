//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>
#include <TF1.h>
#include <TString.h>
#include "../src/MyCastorJet.h"

TString GetJetType(MyCastorJet castorjet){

  // Analyze jet type.

  double depth_jet = castorjet.depth;	// depth
  double fhot_jet = castorjet.fhot;		// fhot
  double fem_jet = castorjet.fem;
  double sigmaz_jet = castorjet.sigmaz;
  double width_jet = castorjet.width;
  double det_energy = castorjet.energy;				  
  double det_phi = castorjet.phi;

  TString detjettype = "other";
  // Count pions.
  if( ! (depth_jet > -14450. && det_energy < 175.) ){ // Most likely not a pion.
    if( ! (depth_jet > -14460. && det_energy > 175.) ){ // Most likely not a pion.
      if( ! (fem_jet > 0.95) ){ // Most likely no pion.
 //       hPion_energy->Fill( det_energy );
        detjettype = "had";
      } // Not a pion.
    } // Not a pion.
  } // Not a pion.
	
  // Count electrons.
  if( ! (fhot_jet < 0.45) ){
    if( ! (fem_jet < 0.9) ){
      if( ! ( sigmaz_jet > 30. && det_energy < 75.) ){
        if( ! ( sigmaz_jet > 40. && det_energy > 75. ) ){
          if( ! ( width_jet > 11.5) ){
	    if( ! ( depth_jet < -14450. && det_energy < 125.) ){
 	      if( ! ( depth_jet < -14460. && det_energy > 125.) ){						    
//		hElectron_energy->Fill( det_energy );
		detjettype = "em";
	      }
            }				
          }
	}      
      }
    }
  } // Not an electron.
  return detjettype;
}





