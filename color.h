// Code by  Tom Cornelis.

#ifndef COLOR_H
#define COLOR_H

#include <map>
#include <TColor.h>
/*
int getColor(TString color){
  std::map<TString, int> colors;
  colors["kWhite"]         = 0;
  colors["kBlack"]         = 1;
  colors["kGray"]         = 920;
  colors["kRed"]         = 632;
  colors["kGreen"]         = 416;
  colors["kBlue"]         = 600;
  colors["kYellow"]        = 400;
  colors["kMagenta"]        = 616;
  colors["kCyan"]         = 432;
  colors["kOrange"]        = 800;
  colors["kSpring"]        = 820;
  colors["kTeal"]        = 840;
  colors["kAzure"]        = 860;
  colors["kViolet"]        = 880;
  colors["kPink"]        = 900;
  colors["kUABlue"]        = TColor::GetColor("#003D64");
  colors["kUABlue50"]        = TColor::GetColor("#60A6BF");
  colors["kUARed"]        = TColor::GetColor("#A10040");
  colors["kUARed50"]        = TColor::GetColor("#CF678F");
  colors["kUABrown"]        = TColor::GetColor("#D79A46");
  colors["kUABrown50"]        = TColor::GetColor("#EBCDA3");

  TString knownColor = colors.find(color);
  if(knownColor != colors.end()) return knownColor->second;
  if(color.IsDigit()) return color.Atoi();
  return 1;
}
*/
int getColor(int color){
  std::map<int, int> colors;
  colors[0]         = 1;
//  colors[1]         = TColor::GetColor("#ff9933");
  colors[1]  	    = 1;
  colors[2]	    = 2; // kGreen + 2
  colors[3]         = 418;
  colors[4]         = 801;
  colors[5]         = 600;
  colors[6]        = 617;
  colors[7]        = 932; // kMagenta + 1
  colors[8]         = 433; // kCyan + 1
  colors[9]        = 401; // kYellow + 1
  colors[10]        = 820;
  colors[11]        = 840;
  colors[12]        = 860;
  colors[13]        = 880;
  colors[14]        = 900;
  colors[15]        = TColor::GetColor("#003D64");
  colors[16]        = TColor::GetColor("#60A6BF");
  colors[17]        = TColor::GetColor("#A10040");
  colors[18]        = TColor::GetColor("#CF678F");
  colors[19]        = TColor::GetColor("#D79A46");
  colors[20]        = TColor::GetColor("#EBCDA3");
  colors[21]	    = TColor::GetColor("#882D61");
  colors[22]	    = TColor::GetColor("#AA8439");


  std::map<int,int>::iterator knownColor = colors.find(color);
//  if(knownColor != colors.end()) return knownColor->second;
return knownColor->second;  
//  if(color.IsDigit()) return color.Atoi();
//  return 1;
}


#endif
