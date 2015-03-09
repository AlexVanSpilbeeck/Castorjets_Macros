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

#include "../../../RooUnfold-1.1.1/src/RooUnfold.h"

#include "color.h"
#include "NonZeroMinimum.h"
#include "MacroPlot_Unfold_energy_MC.h"
#include "Function_FirstPlot.h"
using namespace std;

int main(){
   gROOT->SetStyle ("Plain");
   gStyle->SetPalette(1);

//   TString date = "20141215_03"; 		// Gen-Det jets matched hardest to hardest.
//   TString date = "20141212_05"; 		// Gen jet matched to closest det jet.
//   TString date = "Match_phi_02"; 
//   TString date = "20150107_had_em_jets"; 	// Hardest gen jet matched to closest (phi) det jet. 
//   TString date = "20150199_had_em_jets"; 
//   TString date = "20150112_had_em_jets"; 
//   TString date = "20150119_had_em_jets"; 

    TString date = "20150302";
    TString numb = "10000000"; TString suffix = "_173_1_iQO.root";
//    TString numb = "1000000"; TString suffix = "_13_1_irJ.root";
//    TString numb = "100000"; TString suffix = "_129_1_D3d.root";
    
   /// Draw all plots on one canvas. ///

   TCanvas *can_comparison = new TCanvas("Canvas_comparison", "Canvas_comparison", 1.);   
     can_comparison->SetLeftMargin(0.18);
     can_comparison->SetRightMargin(0.01);
     can_comparison->SetBottomMargin(0.20);

   TLegend *legend = new TLegend(0.65, 0.75, 0.99, 0.99);
     legend->SetFillColor( kWhite );

   std::vector<TString> filenames;
   std::vector<TString> legendEntries;
   std::vector<TString> plotVariables;
   std::vector<TString> jetSelection;
     std::vector<TString> plotX;
     std::vector<TString> plotY;
     std::vector<TString> plotTitle;
     std::vector<bool> project_2D;
   std::vector<int>     colours;
   
   std::map<TString, TString> axis_of_interest;

   cout << "Files" << endl;

 
   filenames.push_back("LoopRootFiles/20150306_test_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.400000_10000000_0_sectors.root");
     legendEntries.push_back("#eta band: 0.4");
     colours.push_back(getColor(1));
     
   filenames.push_back("LoopRootFiles/20150306_test_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.500000_10000000_0_sectors.root");
     legendEntries.push_back("#eta band: 0.5");
     colours.push_back(getColor(2));
     
   filenames.push_back("LoopRootFiles/20150306_test_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.600000_10000000_0_sectors.root");
     legendEntries.push_back("#eta band: 0.6");
     colours.push_back(getColor(3));
     
   filenames.push_back("LoopRootFiles/20150306_test_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.700000_10000000_0_sectors.root");
     legendEntries.push_back("#eta band: 0.7");
     colours.push_back(getColor(4));
     
   filenames.push_back("LoopRootFiles/20150306_test_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.800000_10000000_0_sectors.root");
     legendEntries.push_back("#eta band: 0.8");
     colours.push_back(getColor(5));
     
   filenames.push_back("LoopRootFiles/20150306_test_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.900000_10000000_0_sectors.root");
     legendEntries.push_back("#eta band: 0.9");
     colours.push_back(getColor(6));                    

   filenames.push_back("LoopRootFiles/20150306_test_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.000000_10000000_0_sectors.root");
     legendEntries.push_back("#eta band: 1.0");
     colours.push_back(getColor(7));   
     
   filenames.push_back("LoopRootFiles/20150306_test_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.100000_10000000_0_sectors.root");
     legendEntries.push_back("#eta band: 1.1");
     colours.push_back(getColor(8));  

   filenames.push_back("LoopRootFiles/20150306_test_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.200000_10000000_0_sectors.root");
     legendEntries.push_back("#eta band: 1.2");
     colours.push_back(getColor(9));   
     
   filenames.push_back("LoopRootFiles/20150306_test_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.300000_10000000_0_sectors.root");
     legendEntries.push_back("#eta band: 1.3");
     colours.push_back(getColor(10));   
     
   filenames.push_back("LoopRootFiles/20150306_test_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.400000_10000000_0_sectors.root");
     legendEntries.push_back("#eta band: 1.4");
     colours.push_back(getColor(11));                     


   TCanvas *can = new TCanvas("Test", "Test", 1.);
   TH1D* histo;
   TH1D* original;  
   double max_val, min_val;   
   TString drawoptions="";
   
   // -- Gen over Det ratio.
 
   for( int file = 0; file < filenames.size(); file++){ 
     PlotCorrectionFactors_BinByBin( histo, filenames[file] , legendEntries[file], "");  
     First_Plot( original, histo, file, max_val, min_val);     
         
     histo->SetLineColor( colours[file] );	 
     legend->AddEntry( histo, legendEntries[file], "l");
     can->cd();
     histo->Draw( "hist" + drawoptions );
     drawoptions = "same";
   }
   legend->Draw();     
   can->SaveAs("Gen-over-det_factors.pdf");
   
   // -- Bayesian unfolding.
   
   drawoptions = "";
   max_val = 0, min_val = 0;
   for( int file = 0; file < filenames.size(); file++){ 
     PlotCorrectionFactors_Bayesian( histo, filenames[file] , legendEntries[file], "");  
     First_Plot( original, histo, file, max_val, min_val);     
         
     histo->SetLineColor( colours[file] );	 
     can->cd();
     histo->Draw( "hist" + drawoptions );
     drawoptions = "same";
   }   
   legend->Draw();
   can->SaveAs("Bayesian_factors.pdf");
   
   // -- Gen level distribution.
   
   drawoptions = "";
   max_val = 0, min_val = 0;
   for( int file = 0; file < filenames.size(); file++){ 
     Plot_GenLevel( histo, filenames[file] , legendEntries[file], "");  
     First_Plot( original, histo, file, max_val, min_val);     
         
     histo->SetLineColor( colours[file] );	 
     can->cd();
     histo->Draw( "hist" + drawoptions );
     drawoptions = "same";
   }   
   legend->Draw();
   can->SaveAs("Gen_level_energy.pdf");   
   
   // -- Det level distribution.
   
   drawoptions = "";
   max_val = 0, min_val = 0;
   for( int file = 0; file < filenames.size(); file++){ 
     Plot_DetLevel( histo, filenames[file] , legendEntries[file], "");  
     First_Plot( original, histo, file, max_val, min_val);     
         
     histo->SetLineColor( colours[file] );	 
     can->cd();
     histo->Draw( "hist" + drawoptions );
     drawoptions = "same";
   }   
   legend->Draw();
   can->SaveAs("Det_level_energy.pdf");    
   
     
   
   // -- Bin-by-bin unfolded (RooUnfold).
   
   drawoptions = "";
   max_val = 0, min_val = 0;
   for( int file = 0; file < filenames.size(); file++){ 
     PlotCorrectionFactors_BinByBin_Unfolded( histo, filenames[file] , legendEntries[file], "");  
     First_Plot( original, histo, file, max_val, min_val);     
         
     histo->SetLineColor( colours[file] );	 
     can->cd();
     histo->Draw( "hist" + drawoptions );
     drawoptions = "same";
   }   
   legend->Draw();
   can->SaveAs("RooUnfold_binbybin.pdf");   
   
 
   
     
 
  return(0); 
}
