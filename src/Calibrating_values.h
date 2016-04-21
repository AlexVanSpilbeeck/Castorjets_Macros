double Energy_sector_0_data( double edet ){
  double alpha = 1.92075;
  double beta = -0.0253481;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_0_MC( double edet ){
  double alpha = -3.2699;
  double beta = 0.680766;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_1_data( double edet ){
  double alpha = 2.01785;
  double beta = -0.0402767;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_1_MC( double edet ){
  double alpha = 1.07037;
  double beta = 0.0823896;
  double gamma = 130.603;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_2_data( double edet ){
  double alpha = 1.98282;
  double beta = -0.0294956;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_2_MC( double edet ){
  double alpha = 1.48333;
  double beta = 0.0139331;
  double gamma = 1.21439;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_3_data( double edet ){
  double alpha = 2.05586;
  double beta = -0.0371076;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_3_MC( double edet ){
  double alpha = -0.462967;
  double beta = 0.292432;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_4_data( double edet ){
  double alpha = 2.00012;
  double beta = -0.028161;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_4_MC( double edet ){
  double alpha = 3.26417;
  double beta = -0.240636;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_5_data( double edet ){
  double alpha = 2.12578;
  double beta = -0.0552012;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_5_MC( double edet ){
  double alpha = -0.142138;
  double beta = 0.242701;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_6_data( double edet ){
  double alpha = 1.98523;
  double beta = -0.028304;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_6_MC( double edet ){
  double alpha = 1.0714;
  double beta = 0.0801237;
  double gamma = 322.285;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_7_data( double edet ){
  double alpha = 2.15677;
  double beta = -0.0589826;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_7_MC( double edet ){
  double alpha = 1.44197;
  double beta = 0.0307552;
  double gamma = 49.023;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_8_data( double edet ){
  double alpha = 1.8424;
  double beta = -0.0176457;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_8_MC( double edet ){
  double alpha = 1.26168;
  double beta = 0.0452016;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_9_data( double edet ){
  double alpha = 1.80048;
  double beta = -0.0185742;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_9_MC( double edet ){
  double alpha = 1.33757;
  double beta = 0.0198719;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_10_data( double edet ){
  double alpha = 1.47354;
  double beta = 0.0505798;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_10_MC( double edet ){
  double alpha = 0.95256;
  double beta = 0.085281;
  double gamma = 134.648;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_11_data( double edet ){
  double alpha = 1.26475;
  double beta = 0.102908;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_11_MC( double edet ){
  double alpha = 1.13365;
  double beta = 0.0716389;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_12_data( double edet ){
  double alpha = -11.6994;
  double beta = 1.9918;
  double gamma = 948.012;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_12_MC( double edet ){
  double alpha = -15.5166;
  double beta = 2.4626;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_13_data( double edet ){
  double alpha = -7.62093;
  double beta = 1.47218;
  double gamma = 661.499;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_13_MC( double edet ){
  double alpha = -14.7953;
  double beta = 2.35331;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_14_data( double edet ){
  double alpha = 0.648146;
  double beta = 0.155665;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_14_MC( double edet ){
  double alpha = 0.294397;
  double beta = 0.167487;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_15_data( double edet ){
  double alpha = 1.83213;
  double beta = -0.023728;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_15_MC( double edet ){
  double alpha = -2.95464;
  double beta = 0.622508;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_0_data( double edet ){
  double alpha = 1.92075;
  double beta = -0.0253481;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_0_MC( double edet ){
  double alpha = -3.2699;
  double beta = 0.680766;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_1_data( double edet ){
  double alpha = 2.01785;
  double beta = -0.0402767;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_1_MC( double edet ){
  double alpha = 1.07037;
  double beta = 0.0823896;
  double gamma = 130.603;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_2_data( double edet ){
  double alpha = 1.98282;
  double beta = -0.0294956;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_2_MC( double edet ){
  double alpha = 1.48333;
  double beta = 0.0139331;
  double gamma = 1.21439;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_3_data( double edet ){
  double alpha = 2.05586;
  double beta = -0.0371076;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_3_MC( double edet ){
  double alpha = -0.462967;
  double beta = 0.292432;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_4_data( double edet ){
  double alpha = 2.00012;
  double beta = -0.028161;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_4_MC( double edet ){
  double alpha = 3.26417;
  double beta = -0.240636;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_5_data( double edet ){
  double alpha = 2.12578;
  double beta = -0.0552012;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_5_MC( double edet ){
  double alpha = -0.142138;
  double beta = 0.242701;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_6_data( double edet ){
  double alpha = 1.98523;
  double beta = -0.028304;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_6_MC( double edet ){
  double alpha = 1.0714;
  double beta = 0.0801237;
  double gamma = 322.285;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_7_data( double edet ){
  double alpha = 2.15677;
  double beta = -0.0589826;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_7_MC( double edet ){
  double alpha = 1.44197;
  double beta = 0.0307552;
  double gamma = 49.023;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_8_data( double edet ){
  double alpha = 1.8424;
  double beta = -0.0176457;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_8_MC( double edet ){
  double alpha = 1.26168;
  double beta = 0.0452016;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_9_data( double edet ){
  double alpha = 1.80048;
  double beta = -0.0185742;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_9_MC( double edet ){
  double alpha = 1.33757;
  double beta = 0.0198719;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_10_data( double edet ){
  double alpha = 1.47354;
  double beta = 0.0505798;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_10_MC( double edet ){
  double alpha = 0.95256;
  double beta = 0.085281;
  double gamma = 134.648;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_11_data( double edet ){
  double alpha = 1.26475;
  double beta = 0.102908;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_11_MC( double edet ){
  double alpha = 1.13365;
  double beta = 0.0716389;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_12_data( double edet ){
  double alpha = -11.6994;
  double beta = 1.9918;
  double gamma = 948.012;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_12_MC( double edet ){
  double alpha = -15.5166;
  double beta = 2.4626;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_13_data( double edet ){
  double alpha = -7.62093;
  double beta = 1.47218;
  double gamma = 661.499;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_13_MC( double edet ){
  double alpha = -14.7953;
  double beta = 2.35331;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_14_data( double edet ){
  double alpha = 0.648146;
  double beta = 0.155665;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_14_MC( double edet ){
  double alpha = 0.294397;
  double beta = 0.167487;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_15_data( double edet ){
  double alpha = 1.83213;
  double beta = -0.023728;
  double gamma = 1;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

double Energy_sector_15_MC( double edet ){
  double alpha = -2.95464;
  double beta = 0.622508;
  double gamma = 1000;
return edet * ( (alpha + beta * log( edet + gamma ) )*(edet > 100) + (edet < 100)*1 ); 
}

