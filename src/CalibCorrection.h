#ifndef CalibCorrection_h
#define CalibCorrection_h

#include <TMath.h>
#include <iostream>
#include <vector>

double FillCorrectionVectors( std::vector<double>& lowedge, std::vector<double>& muval){

//-- Calibration with method described by Kostas.
/*
The procedure is very simple: in each bin of the true energy you make a distribution of the response Ereco/Etrue, of Etrue, and of Ereco. 
Then, you measure the peak position (after a Gaussian fit) of the response R = <Ereco/Etrue>, the average <Etrue>, and the average <Ereco>.
Finally, you make a graph with:

	a. R vs <Etrue>, to show the response as a function of Etrue (remember: the response is a property of the TRUE jet)
	b. 1/R vs <Ereco>, to get the correction as a function of Ereco (remember: the correction is a property of the MEASURED jet)

After you apply the correction of step b above to the measured jets, you can repeat step a to demonstrate the closure.

The above procedure ensures that you won't get biased results.
*/

/**************
* Wider bins. *
**************/

lowedge.push_back(	40.3823	); muval.push_back(	1.30966);
lowedge.push_back(	56.0522	); muval.push_back(	1.75545);
lowedge.push_back(	93.095	); muval.push_back(	1.75001);
lowedge.push_back(	122.596	); muval.push_back(	1.76632);
lowedge.push_back(	157.459	); muval.push_back(	1.78458);
lowedge.push_back(	188.561	); muval.push_back(	1.81156);
lowedge.push_back(	221.581	); muval.push_back(	1.81864);
lowedge.push_back(	254.099	); muval.push_back(	1.81921);
lowedge.push_back(	286.822	); muval.push_back(	1.81963);
lowedge.push_back(	318.67	); muval.push_back(	1.82801);
lowedge.push_back(	353.544	); muval.push_back(	1.81966);
lowedge.push_back(	386.179	); muval.push_back(	1.82516);
lowedge.push_back(	423.175	); muval.push_back(	1.8073);
lowedge.push_back(	454.399	); muval.push_back(	1.80166);
lowedge.push_back(	485.288	); muval.push_back(	1.82158);
lowedge.push_back(	514.378	); muval.push_back(	1.8326);
lowedge.push_back(	553.58	); muval.push_back(	1.78887);
lowedge.push_back(	583.282	); muval.push_back(	1.79078);
lowedge.push_back(	609.254	); muval.push_back(	1.82241);
lowedge.push_back(	655.368	); muval.push_back(	1.79466);
lowedge.push_back(	677.054	); muval.push_back(	1.78651);
lowedge.push_back(	713.051	); muval.push_back(	1.75962);
lowedge.push_back(	741.064	); muval.push_back(	1.86337);
lowedge.push_back(	808.318	); muval.push_back(	1.70804);
lowedge.push_back(	817.619	); muval.push_back(	1.83185);
lowedge.push_back(	835.263	); muval.push_back(	1.84964);
lowedge.push_back(	838.846	); muval.push_back(	1.87311);



  return 0;

}




double CalibratedDet( std::vector<double> lowedge, std::vector<double> muval, double Edet){

 // -- Loop over energies until we found the bin within which our energy lies.
 
 double Egen = Edet * 1.3;
 
 for( int bin = 0; bin < lowedge.size() - 1; bin++){
   if( Edet > lowedge[bin] && Edet < lowedge[bin + 1]){
     Egen = Edet * muval[bin];
     break;
   }
 }
 if( Edet > Edet > lowedge[ lowedge.size() ] ){
    Egen = Edet * muval[ lowedge.size() ];
 }

 return Egen;
}

// Overloading calibration functions depending on form of calibration function.


// -- Piece wise: log form & constant.
double CalibratedDet( double Edet, double par0, double par1, double par2, double par3){

 if( Edet > 200.)  return Edet * (par0 + par1 * log( Edet + par2 ));
 else return Edet * par3;
}


// -- Piece wise: log form & slope.
double CalibratedDet( double Edet, double par0, double par1, double par2, double par3, double par4){

 if( Edet > 200.)  return Edet * (par0 + par1 * log( Edet + par2 ));
 else return Edet * (par3 + par4 * Edet);
}

// -- Piece wise: log form & slope.
 double CalibratedDet( double Edet ){

  // -- First iteration.
  double I_0 = -2.51;
  double I_1 = .60;
  double I_2 = 1160.;
  double I_3 = 1.7895;
  double I_4 = 0.000023;

  double Egen = Edet;

  if( Edet < 200.)  Egen = Edet * (I_0 + I_1 * log( Edet + I_2 ));
  else Egen = Edet * (I_3 + I_4 * Edet);
/*
  // -- Second iteration.
  double II_0 = -2.51;
  double II_1 = .60;
  double II_2 = 1160.;
  double II_3 = 1.7895;
  double II_4 = 0.000023;
  double II_2 = 1160.;
  double II_3 = 1.7895;
  double II_4 = 0.000023;
*/
}


#endif
