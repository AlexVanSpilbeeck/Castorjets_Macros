#ifndef CalibCorrection_h
#define CalibCorrection_h

#include <TMath.h>
#include <iostream>
#include <vector>
#include <map>

//#include "Calibrating_Values.h"

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

lowedge.push_back(	39.7911	); muval.push_back(	1.84707);
lowedge.push_back(	53.0241	); muval.push_back(	1.80058);
lowedge.push_back(	89.703	); muval.push_back(	1.79169);
lowedge.push_back(	118.547	); muval.push_back(	1.80295);
lowedge.push_back(	152.467	); muval.push_back(	1.81477);
lowedge.push_back(	182.949	); muval.push_back(	1.83476);
lowedge.push_back(	215.637	); muval.push_back(	1.83846);
lowedge.push_back(	245.995	); muval.push_back(	1.85344);
lowedge.push_back(	277.452	); muval.push_back(	1.85631);
lowedge.push_back(	310.139	); muval.push_back(	1.85408);
lowedge.push_back(	341.752	); muval.push_back(	1.85137);
lowedge.push_back(	373.07	); muval.push_back(	1.85792);
lowedge.push_back(	402.119	); muval.push_back(	1.86838);
lowedge.push_back(	436.161	); muval.push_back(	1.8614);
lowedge.push_back(	464.443	); muval.push_back(	1.88176);
lowedge.push_back(	496.009	); muval.push_back(	1.88265);
lowedge.push_back(	523.989	); muval.push_back(	1.89504);
lowedge.push_back(	565.654	); muval.push_back(	1.84234);
lowedge.push_back(	586.786	); muval.push_back(	1.89034);
lowedge.push_back(	632.186	); muval.push_back(	1.86768);
lowedge.push_back(	660.973	); muval.push_back(	1.89525);
lowedge.push_back(	761.852	); muval.push_back(	1.83531);
lowedge.push_back(	864.03	); muval.push_back(	1.91204);




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
  double a_2 = -1.4738 	* 1E-4; 
  double a_1 =  9.83	* 1E-3;
  double a_0 =  1.684	* 1E0;
  double b_0 = -1.224	* 1E0;
  double b_1 =  4.4	* 1E-1;
  double b_2 =  9.1	* 1E2;
  double c_0 =  1.82560 * 1E0;
  double c_1 =  9.0	* 1E-5;


  double Egen = Edet;

  if( Edet < 60. )  	Egen = Edet * (a_2 * Edet*Edet+ a_1 * Edet + a_0 );
  else if( Edet < 200.)	Egen = Edet * (b_0 + b_1 * log( Edet + b_2 ));
  else 			Egen = Edet * (c_0 + c_1 * Edet);
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


/********************************
* Sector dependent calibration. *
********************************/

double FillCorrectionVectors( std::map<int, std::vector<double> >& lowedge, std::map<int, std::vector<double> >& muval){

  std::cout << "// -- FillCorrectionVectors " << std::endl;

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

std::vector<double> lowedge_1, lowedge_2, lowedge_3, lowedge_4, lowedge_5, lowedge_6, lowedge_7, lowedge_8, lowedge_9, lowedge_10, lowedge_11, lowedge_12, lowedge_13, lowedge_14, lowedge_15, lowedge_16;
std::vector<double> muval_1, muval_2, muval_3, muval_4, muval_5, muval_6, muval_7, muval_8, muval_9, muval_10, muval_11, muval_12, muval_13, muval_14, muval_15, muval_16;

// -- SECTOR	1
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- STORAGE
lowedge_1.push_back(	84.33	); muval_1.push_back(	0.569654);
lowedge_1.push_back(	89.9131	); muval_1.push_back(	1.16137);
lowedge_1.push_back(	118.463	); muval_1.push_back(	1.44261);
lowedge_1.push_back(	153.52	); muval_1.push_back(	1.4937);
lowedge_1.push_back(	189.685	); muval_1.push_back(	1.50834);
lowedge_1.push_back(	226.633	); muval_1.push_back(	1.52118);
lowedge_1.push_back(	258.778	); muval_1.push_back(	1.55122);
lowedge_1.push_back(	294.164	); muval_1.push_back(	1.56387);
lowedge_1.push_back(	332.831	); muval_1.push_back(	1.55977);
lowedge_1.push_back(	362.869	); muval_1.push_back(	1.59274);
lowedge_1.push_back(	399.974	); muval_1.push_back(	1.58769);
lowedge_1.push_back(	439.012	); muval_1.push_back(	1.5865);
lowedge_1.push_back(	468.635	); muval_1.push_back(	1.59697);
lowedge_1.push_back(	502.127	); muval_1.push_back(	1.61475);
lowedge_1.push_back(	558.757	); muval_1.push_back(	1.56658);
lowedge_1.push_back(	579.792	); muval_1.push_back(	1.59785);
lowedge_1.push_back(	609.424	); muval_1.push_back(	1.62258);
lowedge_1.push_back(	665.904	); muval_1.push_back(	1.62107);
lowedge_1.push_back(	766.128	); muval_1.push_back(	1.60862);
// -- SECTOR	2
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- STORAGE
lowedge_2.push_back(	78.2182	); muval_2.push_back(	0.616023);
lowedge_2.push_back(	85.5789	); muval_2.push_back(	1.46756);
lowedge_2.push_back(	115.014	); muval_2.push_back(	1.49292);
lowedge_2.push_back(	149.865	); muval_2.push_back(	1.51857);
lowedge_2.push_back(	186.454	); muval_2.push_back(	1.53393);
lowedge_2.push_back(	220.248	); muval_2.push_back(	1.55771);
lowedge_2.push_back(	255.545	); muval_2.push_back(	1.57033);
lowedge_2.push_back(	290.585	); muval_2.push_back(	1.56507);
lowedge_2.push_back(	325.423	); muval_2.push_back(	1.59652);
lowedge_2.push_back(	361.479	); muval_2.push_back(	1.59437);
lowedge_2.push_back(	398.639	); muval_2.push_back(	1.60412);
lowedge_2.push_back(	432.346	); muval_2.push_back(	1.61852);
lowedge_2.push_back(	466.786	); muval_2.push_back(	1.62458);
lowedge_2.push_back(	501.064	); muval_2.push_back(	1.61855);
lowedge_2.push_back(	536.85	); muval_2.push_back(	1.62099);
lowedge_2.push_back(	573.143	); muval_2.push_back(	1.62368);
lowedge_2.push_back(	587.357	); muval_2.push_back(	1.708);
lowedge_2.push_back(	652.402	); muval_2.push_back(	1.66138);
lowedge_2.push_back(	746.162	); muval_2.push_back(	1.6531);
// -- SECTOR	3
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- STORAGE
lowedge_3.push_back(	75.6443	); muval_3.push_back(	0.661765);
lowedge_3.push_back(	84.5462	); muval_3.push_back(	1.25604);
lowedge_3.push_back(	113.053	); muval_3.push_back(	1.49843);
lowedge_3.push_back(	146.875	); muval_3.push_back(	1.53504);
lowedge_3.push_back(	181.832	); muval_3.push_back(	1.56172);
lowedge_3.push_back(	217.414	); muval_3.push_back(	1.57401);
lowedge_3.push_back(	253.013	); muval_3.push_back(	1.57874);
lowedge_3.push_back(	285.914	); muval_3.push_back(	1.60221);
lowedge_3.push_back(	319.931	); muval_3.push_back(	1.62806);
lowedge_3.push_back(	354.199	); muval_3.push_back(	1.62014);
lowedge_3.push_back(	389.47	); muval_3.push_back(	1.63232);
lowedge_3.push_back(	422.228	); muval_3.push_back(	1.65653);
lowedge_3.push_back(	458.961	); muval_3.push_back(	1.63932);
lowedge_3.push_back(	486.142	); muval_3.push_back(	1.67184);
lowedge_3.push_back(	530.576	); muval_3.push_back(	1.6315);
lowedge_3.push_back(	554.137	); muval_3.push_back(	1.66932);
lowedge_3.push_back(	590.769	); muval_3.push_back(	1.69195);
lowedge_3.push_back(	638.198	); muval_3.push_back(	1.71151);
lowedge_3.push_back(	745.778	); muval_3.push_back(	1.65654);
// -- SECTOR	4
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- STORAGE
lowedge_4.push_back(	74.5898	); muval_4.push_back(	0.635569);
lowedge_4.push_back(	83.6214	); muval_4.push_back(	1.21798);
lowedge_4.push_back(	111.65	); muval_4.push_back(	1.50502);
lowedge_4.push_back(	145.134	); muval_4.push_back(	1.55292);
lowedge_4.push_back(	181.663	); muval_4.push_back(	1.56774);
lowedge_4.push_back(	214.638	); muval_4.push_back(	1.59616);
lowedge_4.push_back(	250.039	); muval_4.push_back(	1.59906);
lowedge_4.push_back(	281.023	); muval_4.push_back(	1.6335);
lowedge_4.push_back(	316.868	); muval_4.push_back(	1.63043);
lowedge_4.push_back(	352.789	); muval_4.push_back(	1.64186);
lowedge_4.push_back(	388.4	); muval_4.push_back(	1.63112);
lowedge_4.push_back(	418.497	); muval_4.push_back(	1.65966);
lowedge_4.push_back(	449.383	); muval_4.push_back(	1.67542);
lowedge_4.push_back(	490.078	); muval_4.push_back(	1.66571);
lowedge_4.push_back(	520.403	); muval_4.push_back(	1.67621);
lowedge_4.push_back(	549.425	); muval_4.push_back(	1.69467);
lowedge_4.push_back(	602.338	); muval_4.push_back(	1.63586);
lowedge_4.push_back(	636.926	); muval_4.push_back(	1.68649);
lowedge_4.push_back(	728	); muval_4.push_back(	1.68637);
// -- SECTOR	5
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- STORAGE
lowedge_5.push_back(	78.4314	); muval_5.push_back(	0.614764);
lowedge_5.push_back(	85.2307	); muval_5.push_back(	1.41223);
lowedge_5.push_back(	112.65	); muval_5.push_back(	1.52563);
lowedge_5.push_back(	146.171	); muval_5.push_back(	1.55865);
lowedge_5.push_back(	180.582	); muval_5.push_back(	1.5757);
lowedge_5.push_back(	216.811	); muval_5.push_back(	1.58503);
lowedge_5.push_back(	250.356	); muval_5.push_back(	1.5971);
lowedge_5.push_back(	282.474	); muval_5.push_back(	1.62664);
lowedge_5.push_back(	314.875	); muval_5.push_back(	1.64436);
lowedge_5.push_back(	349.222	); muval_5.push_back(	1.64844);
lowedge_5.push_back(	383.529	); muval_5.push_back(	1.64681);
lowedge_5.push_back(	417.198	); muval_5.push_back(	1.66743);
lowedge_5.push_back(	447.429	); muval_5.push_back(	1.68719);
lowedge_5.push_back(	487.16	); muval_5.push_back(	1.69209);
lowedge_5.push_back(	525.215	); muval_5.push_back(	1.65675);
lowedge_5.push_back(	555.168	); muval_5.push_back(	1.68839);
lowedge_5.push_back(	608.174	); muval_5.push_back(	1.62975);
lowedge_5.push_back(	629.003	); muval_5.push_back(	1.71846);
lowedge_5.push_back(	715.534	); muval_5.push_back(	1.70197);
// -- SECTOR	6
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- STORAGE
lowedge_6.push_back(	84.0552	); muval_6.push_back(	0.57015);
lowedge_6.push_back(	88.2685	); muval_6.push_back(	1.15409);
lowedge_6.push_back(	116.253	); muval_6.push_back(	1.49652);
lowedge_6.push_back(	149.279	); muval_6.push_back(	1.53476);
lowedge_6.push_back(	184.218	); muval_6.push_back(	1.55769);
lowedge_6.push_back(	218.206	); muval_6.push_back(	1.57772);
lowedge_6.push_back(	251.843	); muval_6.push_back(	1.6009);
lowedge_6.push_back(	290.466	); muval_6.push_back(	1.59805);
lowedge_6.push_back(	322.017	); muval_6.push_back(	1.61352);
lowedge_6.push_back(	354.519	); muval_6.push_back(	1.63613);
lowedge_6.push_back(	390.816	); muval_6.push_back(	1.62192);
lowedge_6.push_back(	423.01	); muval_6.push_back(	1.63802);
lowedge_6.push_back(	454.795	); muval_6.push_back(	1.66897);
lowedge_6.push_back(	487.123	); muval_6.push_back(	1.68248);
lowedge_6.push_back(	520.923	); muval_6.push_back(	1.65761);
lowedge_6.push_back(	556.17	); muval_6.push_back(	1.67882);
lowedge_6.push_back(	602.986	); muval_6.push_back(	1.65506);
lowedge_6.push_back(	646.091	); muval_6.push_back(	1.64751);
lowedge_6.push_back(	712.381	); muval_6.push_back(	1.73051);
// -- SECTOR	7
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- STORAGE
lowedge_7.push_back(	88.5585	); muval_7.push_back(	0.531672);
lowedge_7.push_back(	92.6129	); muval_7.push_back(	1.12846);
lowedge_7.push_back(	118.137	); muval_7.push_back(	1.47311);
lowedge_7.push_back(	152.373	); muval_7.push_back(	1.51674);
lowedge_7.push_back(	187.169	); muval_7.push_back(	1.54209);
lowedge_7.push_back(	222.874	); muval_7.push_back(	1.54917);
lowedge_7.push_back(	255.437	); muval_7.push_back(	1.58694);
lowedge_7.push_back(	289.693	); muval_7.push_back(	1.59844);
lowedge_7.push_back(	320.442	); muval_7.push_back(	1.62107);
lowedge_7.push_back(	357.28	); muval_7.push_back(	1.62681);
lowedge_7.push_back(	392.891	); muval_7.push_back(	1.64489);
lowedge_7.push_back(	416.681	); muval_7.push_back(	1.66138);
lowedge_7.push_back(	459.294	); muval_7.push_back(	1.64462);
lowedge_7.push_back(	492.103	); muval_7.push_back(	1.64595);
lowedge_7.push_back(	521.711	); muval_7.push_back(	1.67809);
lowedge_7.push_back(	567.871	); muval_7.push_back(	1.63994);
lowedge_7.push_back(	579.443	); muval_7.push_back(	1.71006);
lowedge_7.push_back(	647.463	); muval_7.push_back(	1.65436);
lowedge_7.push_back(	731.065	); muval_7.push_back(	1.68211);
// -- SECTOR	8
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- STORAGE
lowedge_8.push_back(	93.5252	); muval_8.push_back(	0.577358);
lowedge_8.push_back(	94.1948	); muval_8.push_back(	1.04585);
lowedge_8.push_back(	120.702	); muval_8.push_back(	1.40604);
lowedge_8.push_back(	155.544	); muval_8.push_back(	1.48604);
lowedge_8.push_back(	189.121	); muval_8.push_back(	1.52474);
lowedge_8.push_back(	224.504	); muval_8.push_back(	1.54089);
lowedge_8.push_back(	256.833	); muval_8.push_back(	1.57683);
lowedge_8.push_back(	291.564	); muval_8.push_back(	1.58888);
lowedge_8.push_back(	324.82	); muval_8.push_back(	1.60647);
lowedge_8.push_back(	357.661	); muval_8.push_back(	1.61935);
lowedge_8.push_back(	389.662	); muval_8.push_back(	1.64442);
lowedge_8.push_back(	421.318	); muval_8.push_back(	1.663);
lowedge_8.push_back(	456.038	); muval_8.push_back(	1.64994);
lowedge_8.push_back(	488.665	); muval_8.push_back(	1.67001);
lowedge_8.push_back(	525.094	); muval_8.push_back(	1.66159);
lowedge_8.push_back(	568.619	); muval_8.push_back(	1.62384);
lowedge_8.push_back(	567.948	); muval_8.push_back(	1.73258);
lowedge_8.push_back(	622.804	); muval_8.push_back(	1.74336);
lowedge_8.push_back(	716.071	); muval_8.push_back(	1.71608);
// -- SECTOR	9
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- STORAGE
lowedge_9.push_back(	97.7887	); muval_9.push_back(	0.485989);
lowedge_9.push_back(	98.8052	); muval_9.push_back(	1.08652);
lowedge_9.push_back(	125.435	); muval_9.push_back(	1.40164);
lowedge_9.push_back(	159.959	); muval_9.push_back(	1.45209);
lowedge_9.push_back(	194.458	); muval_9.push_back(	1.49799);
lowedge_9.push_back(	229.051	); muval_9.push_back(	1.51322);
lowedge_9.push_back(	267.909	); muval_9.push_back(	1.52981);
lowedge_9.push_back(	301.42	); muval_9.push_back(	1.54772);
lowedge_9.push_back(	335.937	); muval_9.push_back(	1.55624);
lowedge_9.push_back(	372.052	); muval_9.push_back(	1.57327);
lowedge_9.push_back(	403.132	); muval_9.push_back(	1.58151);
lowedge_9.push_back(	440.645	); muval_9.push_back(	1.57824);
lowedge_9.push_back(	479.033	); muval_9.push_back(	1.58124);
lowedge_9.push_back(	514.243	); muval_9.push_back(	1.57744);
lowedge_9.push_back(	549.256	); muval_9.push_back(	1.5985);
lowedge_9.push_back(	586.491	); muval_9.push_back(	1.58998);
lowedge_9.push_back(	616.078	); muval_9.push_back(	1.59361);
lowedge_9.push_back(	664.249	); muval_9.push_back(	1.60859);
lowedge_9.push_back(	744.143	); muval_9.push_back(	1.62948);
// -- SECTOR	10
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- STORAGE
lowedge_10.push_back(	98.8174	); muval_10.push_back(	1.02214);
lowedge_10.push_back(	100.722	); muval_10.push_back(	0.518177);
lowedge_10.push_back(	125.124	); muval_10.push_back(	1.37903);
lowedge_10.push_back(	161.28	); muval_10.push_back(	1.44993);
lowedge_10.push_back(	197.631	); muval_10.push_back(	1.47938);
lowedge_10.push_back(	234.126	); muval_10.push_back(	1.48552);
lowedge_10.push_back(	270.471	); muval_10.push_back(	1.49283);
lowedge_10.push_back(	306	); muval_10.push_back(	1.52618);
lowedge_10.push_back(	336.656	); muval_10.push_back(	1.55305);
lowedge_10.push_back(	378.199	); muval_10.push_back(	1.53496);
lowedge_10.push_back(	409.957	); muval_10.push_back(	1.55541);
lowedge_10.push_back(	440.47	); muval_10.push_back(	1.59139);
lowedge_10.push_back(	478.051	); muval_10.push_back(	1.57718);
lowedge_10.push_back(	510.271	); muval_10.push_back(	1.60021);
lowedge_10.push_back(	553.374	); muval_10.push_back(	1.57504);
lowedge_10.push_back(	574.225	); muval_10.push_back(	1.61692);
lowedge_10.push_back(	596.542	); muval_10.push_back(	1.6676);
lowedge_10.push_back(	670.16	); muval_10.push_back(	1.59037);
lowedge_10.push_back(	746.13	); muval_10.push_back(	1.62201);
lowedge_10.push_back(	910.985	); muval_10.push_back(	1.7101);
// -- SECTOR	11
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- STORAGE
lowedge_11.push_back(	99.7662	); muval_11.push_back(	1.01165);
lowedge_11.push_back(	100.627	); muval_11.push_back(	0.518893);
lowedge_11.push_back(	126.088	); muval_11.push_back(	1.37597);
lowedge_11.push_back(	159.662	); muval_11.push_back(	1.4635);
lowedge_11.push_back(	196.866	); muval_11.push_back(	1.4803);
lowedge_11.push_back(	231.465	); muval_11.push_back(	1.51034);
lowedge_11.push_back(	269.06	); muval_11.push_back(	1.52142);
lowedge_11.push_back(	299.573	); muval_11.push_back(	1.55194);
lowedge_11.push_back(	333.283	); muval_11.push_back(	1.57018);
lowedge_11.push_back(	365.369	); muval_11.push_back(	1.58292);
lowedge_11.push_back(	398.612	); muval_11.push_back(	1.59738);
lowedge_11.push_back(	437.16	); muval_11.push_back(	1.58634);
lowedge_11.push_back(	472.159	); muval_11.push_back(	1.60174);
lowedge_11.push_back(	495.972	); muval_11.push_back(	1.6389);
lowedge_11.push_back(	539.228	); muval_11.push_back(	1.61799);
lowedge_11.push_back(	574.932	); muval_11.push_back(	1.61803);
lowedge_11.push_back(	605.701	); muval_11.push_back(	1.63513);
lowedge_11.push_back(	667.247	); muval_11.push_back(	1.60173);
lowedge_11.push_back(	752.817	); muval_11.push_back(	1.62029);
// -- SECTOR	12
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- STORAGE
lowedge_12.push_back(	96.5428	); muval_12.push_back(	1.03607);
lowedge_12.push_back(	97.4012	); muval_12.push_back(	0.520889);
lowedge_12.push_back(	121.324	); muval_12.push_back(	1.45347);
lowedge_12.push_back(	154.238	); muval_12.push_back(	1.53071);
lowedge_12.push_back(	188.333	); muval_12.push_back(	1.5446);
lowedge_12.push_back(	222.941	); muval_12.push_back(	1.57222);
lowedge_12.push_back(	252.68	); muval_12.push_back(	1.6116);
lowedge_12.push_back(	285.657	); muval_12.push_back(	1.62549);
lowedge_12.push_back(	314.38	); muval_12.push_back(	1.65887);
lowedge_12.push_back(	351.175	); muval_12.push_back(	1.66387);
lowedge_12.push_back(	374.974	); muval_12.push_back(	1.69755);
lowedge_12.push_back(	411.357	); muval_12.push_back(	1.70086);
lowedge_12.push_back(	445.155	); muval_12.push_back(	1.70417);
lowedge_12.push_back(	478.201	); muval_12.push_back(	1.70287);
lowedge_12.push_back(	504.543	); muval_12.push_back(	1.74083);
lowedge_12.push_back(	522.559	); muval_12.push_back(	1.77986);
lowedge_12.push_back(	576.26	); muval_12.push_back(	1.72765);
lowedge_12.push_back(	601.143	); muval_12.push_back(	1.76417);
lowedge_12.push_back(	621.682	); muval_12.push_back(	1.81153);
lowedge_12.push_back(	732.43	); muval_12.push_back(	1.78758);
// -- SECTOR	13
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- STORAGE
lowedge_13.push_back(	73.8435	); muval_13.push_back(	1.77465);
lowedge_13.push_back(	94.4823	); muval_13.push_back(	1.84307);
lowedge_13.push_back(	119.302	); muval_13.push_back(	1.91361);
lowedge_13.push_back(	145.973	); muval_13.push_back(	1.9554);
lowedge_13.push_back(	171.778	); muval_13.push_back(	2.00205);
lowedge_13.push_back(	192.353	); muval_13.push_back(	2.06824);
lowedge_13.push_back(	217.221	); muval_13.push_back(	2.09877);
lowedge_13.push_back(	241.718	); muval_13.push_back(	2.13425);
lowedge_13.push_back(	260.29	); muval_13.push_back(	2.20978);
lowedge_13.push_back(	279.981	); muval_13.push_back(	2.24197);
lowedge_13.push_back(	307.441	); muval_13.push_back(	2.23091);
lowedge_13.push_back(	324.681	); muval_13.push_back(	2.32246);
lowedge_13.push_back(	350.448	); muval_13.push_back(	2.30895);
lowedge_13.push_back(	377.193	); muval_13.push_back(	2.31633);
lowedge_13.push_back(	389.026	); muval_13.push_back(	2.37018);
lowedge_13.push_back(	422.591	); muval_13.push_back(	2.38074);
lowedge_13.push_back(	464.371	); muval_13.push_back(	2.42362);
lowedge_13.push_back(	539.774	); muval_13.push_back(	2.47611);
// -- SECTOR	14
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- STORAGE
lowedge_14.push_back(	69.9252	); muval_14.push_back(	1.76792);
lowedge_14.push_back(	92.9296	); muval_14.push_back(	1.84329);
lowedge_14.push_back(	120.048	); muval_14.push_back(	1.87646);
lowedge_14.push_back(	147.728	); muval_14.push_back(	1.91658);
lowedge_14.push_back(	172.086	); muval_14.push_back(	1.97258);
lowedge_14.push_back(	196.969	); muval_14.push_back(	2.00773);
lowedge_14.push_back(	219.26	); muval_14.push_back(	2.07366);
lowedge_14.push_back(	243.164	); muval_14.push_back(	2.13255);
lowedge_14.push_back(	270.45	); muval_14.push_back(	2.12617);
lowedge_14.push_back(	291.518	); muval_14.push_back(	2.15743);
lowedge_14.push_back(	316.298	); muval_14.push_back(	2.18548);
lowedge_14.push_back(	333.43	); muval_14.push_back(	2.25985);
lowedge_14.push_back(	354.857	); muval_14.push_back(	2.26327);
lowedge_14.push_back(	391.476	); muval_14.push_back(	2.19561);
lowedge_14.push_back(	414.727	); muval_14.push_back(	2.22735);
lowedge_14.push_back(	434.685	); muval_14.push_back(	2.27387);
lowedge_14.push_back(	451.254	); muval_14.push_back(	2.38028);
lowedge_14.push_back(	504.724	); muval_14.push_back(	2.40624);
// -- SECTOR	15
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- STORAGE
lowedge_15.push_back(	89.6137	); muval_15.push_back(	0.567557);
lowedge_15.push_back(	90.9882	); muval_15.push_back(	1.14256);
lowedge_15.push_back(	117.875	); muval_15.push_back(	1.47212);
lowedge_15.push_back(	150.258	); muval_15.push_back(	1.53756);
lowedge_15.push_back(	183.868	); muval_15.push_back(	1.55863);
lowedge_15.push_back(	217.153	); muval_15.push_back(	1.58988);
lowedge_15.push_back(	249.765	); muval_15.push_back(	1.61734);
lowedge_15.push_back(	281.399	); muval_15.push_back(	1.64724);
lowedge_15.push_back(	314.503	); muval_15.push_back(	1.6437);
lowedge_15.push_back(	350.036	); muval_15.push_back(	1.65535);
lowedge_15.push_back(	381.342	); muval_15.push_back(	1.66869);
lowedge_15.push_back(	410.103	); muval_15.push_back(	1.70137);
lowedge_15.push_back(	443.009	); muval_15.push_back(	1.69044);
lowedge_15.push_back(	476.622	); muval_15.push_back(	1.7061);
lowedge_15.push_back(	510.302	); muval_15.push_back(	1.73274);
lowedge_15.push_back(	545.721	); muval_15.push_back(	1.73162);
lowedge_15.push_back(	569.571	); muval_15.push_back(	1.76056);
lowedge_15.push_back(	624.459	); muval_15.push_back(	1.73116);
lowedge_15.push_back(	690.72	); muval_15.push_back(	1.75616);
lowedge_15.push_back(	792.475	); muval_15.push_back(	1.7488);
// -- SECTOR	16
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- EXTRACT 2D FROM THN
// -- STORAGE
lowedge_16.push_back(	85.4992	); muval_16.push_back(	0.551595);
lowedge_16.push_back(	91.0137	); muval_16.push_back(	1.10764);
lowedge_16.push_back(	119.628	); muval_16.push_back(	1.44599);
lowedge_16.push_back(	154.558	); muval_16.push_back(	1.48675);
lowedge_16.push_back(	191.035	); muval_16.push_back(	1.50968);
lowedge_16.push_back(	227.682	); muval_16.push_back(	1.50439);
lowedge_16.push_back(	262.621	); muval_16.push_back(	1.53654);
lowedge_16.push_back(	295.892	); muval_16.push_back(	1.54831);
lowedge_16.push_back(	331	); muval_16.push_back(	1.57298);
lowedge_16.push_back(	365.51	); muval_16.push_back(	1.58206);
lowedge_16.push_back(	401.449	); muval_16.push_back(	1.59442);
lowedge_16.push_back(	435.82	); muval_16.push_back(	1.60316);
lowedge_16.push_back(	475.457	); muval_16.push_back(	1.59125);
lowedge_16.push_back(	516.136	); muval_16.push_back(	1.58169);
lowedge_16.push_back(	546.032	); muval_16.push_back(	1.61612);
lowedge_16.push_back(	580.577	); muval_16.push_back(	1.61723);
lowedge_16.push_back(	611.3	); muval_16.push_back(	1.62481);
lowedge_16.push_back(	673.784	); muval_16.push_back(	1.59481);
lowedge_16.push_back(	771.506	); muval_16.push_back(	1.61321);

  lowedge[0] = lowedge_1;
  lowedge[1] = lowedge_2;
  lowedge[2] = lowedge_3;
  lowedge[3] = lowedge_4;
  lowedge[4] = lowedge_5;
  lowedge[5] = lowedge_6;
  lowedge[6] = lowedge_7;
  lowedge[7] = lowedge_8;
  lowedge[8] = lowedge_9;
  lowedge[9] = lowedge_10;
  lowedge[10] = lowedge_11;
  lowedge[11] = lowedge_12;
  lowedge[12] = lowedge_13;
  lowedge[13] = lowedge_14;
  lowedge[14] = lowedge_15;
  lowedge[15] = lowedge_16;

  muval[0] = muval_1;
  muval[1] = muval_2;
  muval[2] = muval_3;
  muval[3] = muval_4;
  muval[4] = muval_5;
  muval[5] = muval_6;
  muval[6] = muval_7;
  muval[7] = muval_8;
  muval[8] = muval_9;
  muval[9] = muval_10;
  muval[10] = muval_11;
  muval[11] = muval_12;
  muval[12] = muval_13;
  muval[13] = muval_14;
  muval[14] = muval_15;
  muval[15] = muval_16;


    return 0;

}




double CalibratedDet( std::map<int, std::vector<double> >& lowedge_map, std::map<int, std::vector<double> >& muval_map, double Edet, int sector){

 // -- Loop over energies until we found the bin within which our energy lies.
 
 std::vector<double> lowedge = lowedge_map[sector];
 std::vector<double> muval	= muval_map[sector];
 
 double Egen = Edet * muval[0];

 for( int bin = 0; bin < lowedge.size() - 1; bin++){
   if( Edet > lowedge[bin] && Edet < lowedge[bin + 1] ){
     Egen = Edet * muval[bin];
     break;
   }
 }
 if( Edet > lowedge[ lowedge.size() -1  ]) {
    Egen = Edet * muval[ lowedge.size()-1 ];
 }

 return Egen;
}

/*
double CalibratedDet_function( std::map<int, std::vector<double> >& lowedge_map, std::map<int, std::vector<double> >& muval_map, double Edet, int sector){

 // -- Loop over energies until we found the bin within which our energy lies.
 
 // double (*sector_calibration)(double edet) = NULL;
 if( sector == 0 ){		return Energy_sector_0( Edet ); 	}
 else if( sector == 1 ){     	return Energy_sector_1( Edet );         }
 else if( sector == 2 ){        return Energy_sector_2( Edet );         }
 else if( sector == 3 ){        return Energy_sector_3( Edet );         }
 else if( sector == 4 ){        return Energy_sector_4( Edet );         }
 else if( sector == 5 ){        return Energy_sector_5( Edet );         }
 else if( sector == 6 ){        return Energy_sector_6( Edet );         }
 else if( sector == 7 ){        return Energy_sector_7( Edet );         }
 else if( sector == 8 ){        return Energy_sector_8( Edet );         }
 else if( sector == 9 ){        return Energy_sector_9( Edet );         }
 else if( sector == 10 ){       return Energy_sector_10( Edet );        }
 else if( sector == 11 ){       return Energy_sector_11( Edet );        }
 else if( sector == 12 ){       return Energy_sector_12( Edet );        }
 else if( sector == 13 ){       return Energy_sector_13( Edet );        }
 else if( sector == 14 ){       return Energy_sector_14( Edet );        }
 else if( sector == 15 ){       return Energy_sector_15( Edet );        }

 return Egen;
}
*/

/*************************************************
* Here are the calibration functions per sector. * 
*************************************************/

double Energy_sector_0_data( double edet ){
  double alpha = 3.33381;
  double beta = -0.21675;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_0_MC( double edet ){
  double alpha = 1.3117;
  double beta = 0.0513503;
  double gamma = -135.795;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_1_data( double edet ){
  double alpha = 2.4401;
  double beta = -0.0960558;
  double gamma = 462.468;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_1_MC( double edet ){
  double alpha = 1.49311;
  double beta = 0.0223115;
  double gamma = -149;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_2_data( double edet ){
  double alpha = 1.86306;
  double beta = -0.00689278;
  double gamma = -149;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_2_MC( double edet ){
  double alpha = 6.44598;
  double beta = -0.69526;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_3_data( double edet ){
  double alpha = 2.57384;
  double beta = -0.099085;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_3_MC( double edet ){
  double alpha = 1.27835;
  double beta = 0.0860718;
  double gamma = -146.597;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_4_data( double edet ){
  double alpha = 2.30363;
  double beta = -0.0609511;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_4_MC( double edet ){
  double alpha = 1.65794;
  double beta = -0.0265372;
  double gamma = -149;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_5_data( double edet ){
  double alpha = 1.8892;
  double beta = -0.00856925;
  double gamma = -149;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_5_MC( double edet ){
  double alpha = 1.3039;
  double beta = 0.074386;
  double gamma = -149;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_6_data( double edet ){
  double alpha = 2.50709;
  double beta = -0.0919485;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_6_MC( double edet ){
  double alpha = 2.09381;
  double beta = -0.075838;
  double gamma = 999.94;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_7_data( double edet ){
  double alpha = 2.8995;
  double beta = -0.145498;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_7_MC( double edet ){
  double alpha = 7.87802;
  double beta = -0.885372;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_8_data( double edet ){
  double alpha = 2.146;
  double beta = -0.0554186;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_8_MC( double edet ){
  double alpha = 1.51814;
  double beta = -0.0093327;
  double gamma = -149;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_9_data( double edet ){
  double alpha = 0.744816;
  double beta = 0.134036;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_9_MC( double edet ){
  double alpha = 1.31492;
  double beta = 0.0305184;
  double gamma = -149;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_10_data( double edet ){
  double alpha = 1.34621;
  double beta = 0.0730863;
  double gamma = -48.2154;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_10_MC( double edet ){
  double alpha = 1.28456;
  double beta = 0.0385405;
  double gamma = -149;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_11_data( double edet ){
  double alpha = 1.68625;
  double beta = 0.0225201;
  double gamma = -149;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_11_MC( double edet ){
  double alpha = 1.47187;
  double beta = 0.0102384;
  double gamma = -149;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_12_data( double edet ){
  double alpha = 1.34951;
  double beta = 0.233206;
  double gamma = -149;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_12_MC( double edet ){
  double alpha = 0.90146;
  double beta = 0.253139;
  double gamma = -149;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_13_data( double edet ){
  double alpha = 1.43492;
  double beta = 0.214868;
  double gamma = -149;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_13_MC( double edet ){
  double alpha = -5.80292;
  double beta = 1.31249;
  double gamma = 144.004;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_14_data( double edet ){
  double alpha = 1.62625;
  double beta = 0.0284726;
  double gamma = -149;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_14_MC( double edet ){
  double alpha = 1.57924;
  double beta = -0.0256388;
  double gamma = -149;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_15_data( double edet ){
  double alpha = 0.992968;
  double beta = 0.117341;
  double gamma = 217.136;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}

double Energy_sector_15_MC( double edet ){
  double alpha = -1.29666;
  double beta = 0.383205;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) ) ) ; 
}




double CalibratedDet( double Edet, int sector, TString fileLabel_){

 // -- Loop over energies until we found the bin within which our energy lies.

 double Ecal = 0.;
 
 // Calibrate Data.
 if( fileLabel_ == "data" || fileLabel_ == "Pythia6Z2star" ){ 
   if( sector == 0 ){             Ecal = Energy_sector_0_data( Edet );         }
   else if( sector == 1 ){        Ecal = Energy_sector_1_data( Edet );         }
   else if( sector == 2 ){        Ecal = Energy_sector_2_data( Edet );         }
   else if( sector == 3 ){        Ecal = Energy_sector_3_data( Edet );         }
   else if( sector == 4 ){        Ecal = Energy_sector_4_data( Edet );         }
   else if( sector == 5 ){        Ecal = Energy_sector_5_data( Edet );         }
   else if( sector == 6 ){        Ecal = Energy_sector_6_data( Edet );         }
   else if( sector == 7 ){        Ecal = Energy_sector_7_data( Edet );         }
   else if( sector == 8 ){        Ecal = Energy_sector_8_data( Edet );         }
   else if( sector == 9 ){        Ecal = Energy_sector_9_data( Edet );         }
   else if( sector == 10 ){       Ecal = Energy_sector_10_data( Edet );        }
   else if( sector == 11 ){       Ecal = Energy_sector_11_data( Edet );        }
   else if( sector == 12 ){       Ecal = Energy_sector_12_data( Edet );        }
   else if( sector == 13 ){       Ecal = Energy_sector_13_data( Edet );        }
   else if( sector == 14 ){       Ecal = Energy_sector_14_data( Edet );        }
   else if( sector == 15 ){       Ecal = Energy_sector_15_data( Edet );        }
 }

 // Calibrate MC.
 else{
   if( sector == 0 ){             Ecal = Energy_sector_0_MC( Edet );         }
   else if( sector == 1 ){        Ecal = Energy_sector_1_MC( Edet );         }
   else if( sector == 2 ){        Ecal = Energy_sector_2_MC( Edet );         }
   else if( sector == 3 ){        Ecal = Energy_sector_3_MC( Edet );         }
   else if( sector == 4 ){        Ecal = Energy_sector_4_MC( Edet );         }
   else if( sector == 5 ){        Ecal = Energy_sector_5_MC( Edet );         }
   else if( sector == 6 ){        Ecal = Energy_sector_6_MC( Edet );         }
   else if( sector == 7 ){        Ecal = Energy_sector_7_MC( Edet );         }
   else if( sector == 8 ){        Ecal = Energy_sector_8_MC( Edet );         }
   else if( sector == 9 ){        Ecal = Energy_sector_9_MC( Edet );         }
   else if( sector == 10 ){       Ecal = Energy_sector_10_MC( Edet );        }
   else if( sector == 11 ){       Ecal = Energy_sector_11_MC( Edet );        }
   else if( sector == 12 ){       Ecal = Energy_sector_12_MC( Edet );        }
   else if( sector == 13 ){       Ecal = Energy_sector_13_MC( Edet );        }
   else if( sector == 14 ){       Ecal = Energy_sector_14_MC( Edet );        }
   else if( sector == 15 ){       Ecal = Energy_sector_15_MC( Edet );        }
 }

// std::cout << "***\tEdet to Ecal\t" << Edet << "\t" << Ecal << "\tin sector\t" << sector << endl;

 return Ecal;

}


#endif
