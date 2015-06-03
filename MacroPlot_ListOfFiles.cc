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

#include "../../../RooUnfold-1.1.1/src/RooUnfold.h"

#include "color.h"
#include "NonZeroMinimum.h"
#include "MacroPlot_Unfold_energy_MC.h"

#include "Unfolder.cpp"


//#include "Function_Rebin.h"

//#ifndef FUNCTION_FIRSTPLOT_H
//#include "Function_FirstPlot.h"
//#endif
//#include "Function_make_Tex.h"

using namespace std;

int main(){
   gROOT->SetStyle ("Plain");
   gStyle->SetPalette(1);
   gStyle->SetOptTitle(0);

//   TString label = "20150323_Plots_Calibrating_ak5ak5";
   TString label = "20150330_Plots_unfolding";

   // -- Prepare the necessary vectors.

   std::vector<TString> filenames;
      std::vector<TString> fileLabel;
      std::vector<TString> legendEntries;
      std::vector<int>     colours;
      std::vector<int>	linestyle;
   
   std::vector<TString> data_files;
   
   std::vector<TString> plotVariables;
     std::vector<TString> plotX;
     std::vector<TString> plotY;
     std::vector<TString> plotTitle;    
     std::vector<bool> project_2D;                
   
   std::vector<TString> Plot_list;
   
   int color_index = 0;
   
   int style_of_line = 1; 

   TString setup;  
   TString Ethresh;
 
   std::map<TString, TString> axis_of_interest;
   std::map<TString, TString> legends;
   std::map<TString, TString> xtitle, ytitle, htitle;

   bool vary_eta 	= false;
   bool vary_eI 	= false;
   bool compare_radii 	= false;
   bool after_calibration = false;
   bool after_all 	= false;
   bool compare_had 	= false;
   bool compare_em 	= false;
   bool determine_fit 	= false;
   bool systematics_ 	= true;

   cout << "Files" << endl;

     cout << "What is the setup?\t";
     cin >> setup;
     cout << "What is energy threshold?\t";
     cin >> Ethresh;

   /********
   * DATA. * 
   ********/ 

   TString datafile = "Histograms/ak5_data_" + setup + "_Emin_" + Ethresh + ".000000.root";

   if( systematics_ ){

     datafile = "Histograms/ak5_data_" + setup + "_Emin_" + Ethresh + ".000000.root";
     legends["Histograms/ak5_data_" + setup + "_Emin_" + Ethresh + ".000000.root"] = "Data";
     
     filenames.push_back("Histograms/ak5ak5_" + setup + "_Emin_" + Ethresh + ".000000.root");
       legendEntries.push_back("(ak5, ak5)");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("Actual");    
       legends["Histograms/ak5ak5_" + setup + "_Emin_" + Ethresh + ".000000.root"] = "Pythia6 Z2*";
 
     filenames.push_back("Histograms/ak5ak5_FullSimulation_" + setup + "_Emin_" + Ethresh + ".000000.root");
       legendEntries.push_back("(ak5, ak5) - Full Simulation");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("Full simulation");    
       legends["Histograms/ak5ak5_FullSimulation_" + setup + "_Emin_" + Ethresh + ".000000.root"] = "Full simulation";
   
     filenames.push_back("Histograms/ak5ak5_displaced_" + setup + "_Emin_" + Ethresh + ".000000.root");
       legendEntries.push_back("(ak5, ak5) - Nominal");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("(0,0)");    
       legends["Histograms/ak5ak5_displaced_" + setup + "_Emin_" + Ethresh + ".000000.root"] = "(0,0)";

     filenames.push_back("Histograms/ak5ak5_Pythia84C_" + setup + "_Emin_" + Ethresh + ".000000.root");
       legendEntries.push_back("(ak5, ak5) - Pythia8 (4C)");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("Pythia8 (4C)");
       legends["Histograms/ak5ak5_Pythia84C_" + setup + "_Emin_" + Ethresh + ".000000.root"] = "Pythia8 (4C)";

   }

   if( after_all ){ 
/*     
     filenames.push_back("LoopRootFiles/20150402_unfolding_JetAnalyzer_etaband_0.400000_12137851_0_sectors_had.root");
       legendEntries.push_back("(ak5, ak5) - unfolding");
       colours.push_back(getColor(color_index++)); 
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("ak5ak5_response");      
*/       

     filenames.push_back("ak5ak5_raw_calibrated_sector_E200_JetAnalyzer_etaband_1.400000_12137851_0_sectors_had.root");
       legendEntries.push_back("(ak5, ak5) - unfolding isolated");
       colours.push_back(getColor(color_index++)); 
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("ak5ak5_isolated");         
       
   } // After calibration


   if( compare_radii ){
     //-- All files in this section have E_I = 0, eta_band = 0.4

     if( compare_had ){
     filenames.push_back("");
       legendEntries.push_back("(ak5, ak5)");
       colours.push_back(getColor(color_index++)); 
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("ak5ak5");  
   
     filenames.push_back("");
       legendEntries.push_back("(ak5, ak7)");
       colours.push_back(getColor(color_index++));            
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("ak5ak7");         

     filenames.push_back("");
       legendEntries.push_back("(ak3, ak3)");
       colours.push_back(getColor(color_index++)); 
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("ak3ak3");

     filenames.push_back("");
       legendEntries.push_back("(ak3, ak5)");
       colours.push_back(getColor(color_index++)); 
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("ak3ak5");                         

     filenames.push_back("");
       legendEntries.push_back("(ak3, ak7)");
       colours.push_back(getColor(color_index++)); 
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("ak3ak7");   
       
     filenames.push_back("");
       legendEntries.push_back("(ak7, ak7)");
       colours.push_back(getColor(color_index++));            
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("ak7ak7");  
       
     } // Compare had.

     if( compare_em ){
     filenames.push_back("LoopRootFiles/20150317_ak5ak5_JetAnalyzer_etaband_0.400000_12137851_0_sectors_em.root");
       legendEntries.push_back("(ak5, ak5)");
       colours.push_back(getColor(color_index++)); 
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("ak5ak5");  
    
     filenames.push_back("LoopRootFiles/20150317_ak5ak7_JetAnalyzer_etaband_0.400000_12137851_0_sectors_em.root");
       legendEntries.push_back("(ak5, ak7)");
       colours.push_back(getColor(color_index++));            
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("ak5ak7");         

     filenames.push_back("LoopRootFiles/20150317_ak3ak3_JetAnalyzer_etaband_0.400000_12139019_0_sectors_em.root");
       legendEntries.push_back("(ak3, ak3)");
       colours.push_back(getColor(color_index++)); 
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("ak3ak3");

     filenames.push_back("LoopRootFiles/20150317_ak3ak5_JetAnalyzer_etaband_0.400000_12139019_0_sectors_em.root");
       legendEntries.push_back("(ak3, ak5)");
       colours.push_back(getColor(color_index++)); 
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("ak3ak5");                         

     filenames.push_back("LoopRootFiles/20150317_ak3ak7_JetAnalyzer_etaband_0.400000_12139019_0_sectors_em.root");
       legendEntries.push_back("(ak3, ak7)");
       colours.push_back(getColor(color_index++)); 
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("ak3ak7");   
       
      filenames.push_back("LoopRootFiles/20150317_ak7ak7_JetAnalyzer_etaband_0.400000_12135484_0_sectors_em.root");
       legendEntries.push_back("(ak7, ak7)");
       colours.push_back(getColor(color_index++));            
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("ak7ak7");  
       
     } // Compare had.                   
   }

   if( vary_eI && vary_eta){

     filenames.push_back("LoopRootFiles/20150311_Calib_Iso0_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.400000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.4");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line ); 
       
     filenames.push_back("LoopRootFiles/20150311_Calib_Iso0_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.500000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.5");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );
     
     filenames.push_back("LoopRootFiles/20150311_Calib_Iso0_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.600000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.6");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );
     
     filenames.push_back("LoopRootFiles/20150311_Calib_Iso0_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.700000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.7");
       linestyle.push_back( style_of_line );
       colours.push_back(getColor(color_index++));
     
     filenames.push_back("LoopRootFiles/20150311_Calib_Iso0_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.800000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.8");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );
     
     filenames.push_back("LoopRootFiles/20150311_Calib_Iso0_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.900000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.9");
       colours.push_back(getColor(color_index++));                    
       linestyle.push_back( style_of_line );
  
     filenames.push_back("LoopRootFiles/20150311_Calib_Iso0_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.000000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 1.0");
       colours.push_back(getColor(color_index++));   
       linestyle.push_back( style_of_line );
     
     filenames.push_back("LoopRootFiles/20150311_Calib_Iso0_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.100000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 1.1");
       colours.push_back(getColor(color_index++));  
       linestyle.push_back( style_of_line );

     filenames.push_back("LoopRootFiles/20150311_Calib_Iso0_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.200000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 1.2");
       colours.push_back(getColor(color_index++));   
       linestyle.push_back( style_of_line );
     
     filenames.push_back("LoopRootFiles/20150311_Calib_Iso0_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.300000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 1.3");
       colours.push_back(getColor(color_index++));   
       linestyle.push_back( style_of_line );
     
     filenames.push_back("LoopRootFiles/20150311_Calib_Iso0_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.400000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 1.4");
       linestyle.push_back( style_of_line );
       colours.push_back(getColor(color_index++));                     
   } // Vary eta.

   /**************************
   * MC with varied E_I cut. *
   **************************/ 

   else if(  !vary_eta && vary_eI){
     cout << "Vary EI" << endl;

     filenames.push_back("LoopRootFiles/20150225_iso-cut_000_had_Output_JetAnalyzer_radii_strippedTree_GEN_ak5_DET_ak5_margin_0.500000_9082063_0_sectors.root");
       legendEntries.push_back("E_{I} = 0 GeV");
      colours.push_back(getColor(color_index++));   
      
     filenames.push_back("LoopRootFiles/20150225_iso-cut_001_had_Output_JetAnalyzer_radii_strippedTree_GEN_ak5_DET_ak5_margin_0.500000_9082063_0_sectors.root");
       legendEntries.push_back("E_{I} < 0.01 E_{I}");
      colours.push_back(getColor(color_index++));   
      
     filenames.push_back("LoopRootFiles/20150225_iso-cut_01_had_Output_JetAnalyzer_radii_strippedTree_GEN_ak5_DET_ak5_margin_0.500000_9082063_0_sectors.root");
       legendEntries.push_back("E_{I} < 0.1 E_{I}");
      colours.push_back(getColor(color_index++));   
      
     filenames.push_back("LoopRootFiles/20150225_iso-cut_05_had_Output_JetAnalyzer_radii_strippedTree_GEN_ak5_DET_ak5_margin_0.500000_9082063_0_sectors.root");
       legendEntries.push_back("E_{I} < 0.5 E_{I}");
      colours.push_back(getColor(color_index++));   
      
     filenames.push_back("LoopRootFiles/20150225_iso-cut_1_had_Output_JetAnalyzer_radii_strippedTree_GEN_ak5_DET_ak5_margin_0.500000_9082063_0_sectors.root");
       legendEntries.push_back("E_{I} < 1 E_{I}");
      colours.push_back(getColor(color_index++));          
   }

   /********************************************
   * MC without E_I cut and variable eta band. *
   ********************************************/

   else if( vary_eta && !vary_eI){
     int style_of_line = 1;

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.400000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.4, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );
       fileLabel.push_back("ak5ak5_eta04");

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.500000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.5, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );
       fileLabel.push_back("ak5ak5_eta05");

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.600000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.6, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );
       fileLabel.push_back("ak5ak5_eta06");

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.700000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.7, no E_{I} cut");
       linestyle.push_back( style_of_line );
       colours.push_back(getColor(color_index++));
       fileLabel.push_back("ak5ak5_eta07");

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.800000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.8, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );
       fileLabel.push_back("ak5ak5_eta08");

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.900000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.9, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );
       fileLabel.push_back("ak5ak5_eta09");

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.000000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 1.0, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );
       fileLabel.push_back("ak5ak5_eta10");

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.100000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 1.1, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );
       fileLabel.push_back("ak5ak5_eta11");

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.200000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 1.2, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );
       fileLabel.push_back("ak5ak5_eta12");

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.300000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 1.3, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );
       fileLabel.push_back("ak5ak5_eta13");

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.400000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 1.4, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );
       fileLabel.push_back("ak5ak5_eta14");
   }   

   /**************************************
   * Prepare the plots and their labels. * 
   **************************************/	 

   std::map<TString, TString> xTitle;
   std::map<TString, TString> yTitle; 
   std::map<TString, TString> save; 
   std::map<TString, TString> labels;
   std::map<TString, TString> draw;
   
   std::map<TString, bool> logx;
   std::map<TString, bool> logy;
   std::map<TString, bool> normalise;

   /******************
   *  Control plots. *
   ******************/

/*
   Plot_list.push_back("Plot_PhiDiff");
     xTitle["Plot_PhiDiff"] = "#Delta#varphi";
     yTitle["Plot_PhiDiff"] = "#frac{dN}{d#Delta#varphi}";
     labels["Plot_PhiDiff"] = label;
     save["Plot_PhiDiff"] = "PhiDiff";    
     logx["Plot_PhiDiff"] = false;
     logy["Plot_PhiDiff"] = false;
     normalise["Plot_PhiDiff"]= true;
     
   Plot_list.push_back("Plot_EtaDiff");
     xTitle["Plot_EtaDiff"] = "#Delta#eta";
     yTitle["Plot_EtaDiff"] = "#frac{dN}{d#Delta#eta}";
     labels["Plot_EtaDiff"] = label;
     save["Plot_EtaDiff"] = "EtaDiff"; 
     logx["Plot_EtaDiff"] = false;
     logy["Plot_EtaDiff"] = false;
     normalise["Plot_Etadiff"]= true;
     
   Plot_list.push_back("Plot_Distance");
     xTitle["Plot_Distance"] = "#DeltaR";
     yTitle["Plot_Distance"] = "#frac{dN}{d#DeltaE}";
     labels["Plot_Distance"] = label;
     save["Plot_Distance"] = "Distance";
     logx["Plot_Distance"] = false;
     logy["Plot_Distance"] = false;
     normalise["Plot_Distance"]= true;
     
   Plot_list.push_back("Plot_EnergyIsolation");
     xTitle["Plot_EnergyIsolation"] = "E_{I} (GeV)";
     yTitle["Plot_EnergyIsolation"] = "#frac{dN}{dE_{I}}";
     labels["Plot_EnergyIsolation"] = label;
     save["Plot_EnergyIsolation"] = "EnergyIsolation";      
     logx["Plot_EnergyIsolation"] = true;
     logy["Plot_EnergyIsolation"] = true;
     normalise["Plot_EnergyIsolation"]= true;

   Plot_list.push_back("Plot_GenLevel");
     xTitle["Plot_GenLevel"] = "E_{gen} (GeV)";
     yTitle["Plot_GenLevel"] = "#frac{dN}{dE_{gen}}";
     labels["Plot_GenLevel"] = label;
     save["Plot_GenLevel"] = "Gen_energy"; 
     normalise["Plot_GenLevel"]= true;
     
   Plot_list.push_back("Plot_DetLevel");
     xTitle["Plot_DetLevel"] = "E_{det} (GeV)";
     yTitle["Plot_DetLevel"] = "#frac{dN}{dE_{det}}";
     labels["Plot_GenLevel"] = label;
     save["Plot_DetLevel"] = "Det_energy";
     normalise["Plot_DetLevel"]= true;

   Plot_list.push_back("Plot_GenLevelEta");
     xTitle["Plot_GenLevelEta"] = "#Delta#eta";
     yTitle["Plot_GenLevelEta"] = "#frac{dN}{#Delta#Eta}";
     labels["Plot_GenLevelEta"] = label;
     save["Plot_GenLevelEta"] = "Gen_eta";  
     normalise["Plot_GenLevelEta"]= true;
*/    
    
   /***************************************
   * Calibration - response and all that. *
   ***************************************/
/*   
   Plot_list.push_back("Plot_2D_Energy_Response");
     xTitle["Plot_2D_Energy_response"] = "";
     yTitle["Plot_2D_Energy_response"] = "";
     labels["Plot_2D_Energy_response"] = label;
     save["Plot_2D_Energy_response"] = "";
     draw["Plot_2D_Energy_response"] = "";   
    
   Plot_list.push_back("Plot_JetResponse");
     xTitle["Plot_JetResponse"] = "(E_{det}-E_{gen})/E_{gen}";
     yTitle["Plot_JetResponse"] = "#frac{1}{N}.#frac{E_{det}-E_{gen}}{E_{gen}}";
     labels["Plot_JetResponse"] = label;
     save["Plot_JetResponse"] = "Jet_Response";   
     logy["Plot_JetResponse"] = false;
     normalise["Plot_JetResponse"] = true;
*/
   /************
   * Unfolding *
   ************/

/*
   Plot_list.push_back("Plot_Bayes_vs_Gen");
     xTitle["Plot_Bayes_vs_Gen"] = "E";
     yTitle["Plot_Bayes_vs_Gen"] = "#frac{Bayes}{Gen}";
     labels["Plot_Bayes_vs_Gen"] = label;
     save["Plot_Bayes_vs_Gen"] = "Bayes_vs_Gen";
     draw["Plot_Bayes_vs_Gen"] = "e";

   Plot_list.push_back("PlotCorrectionFactors_Bayes");
     xTitle["PlotCorrectionFactors_Bayes"] = "E (GeV)";
     yTitle["PlotCorrectionFactors_Bayes"] = "#frac{dN}{dE}";
     labels["PlotCorrectionFactors_Bayes"] = label;
     save["PlotCorrectionFactors_Bayes"] = "Bayes";

   Plot_list.push_back("Bayes_iterations");
     xTitle["Bayes_iterations"] = "";
     yTitle["Bayes_iterations"] = "#chi^{2}";
     labels["Bayes_iterations"] = "10";
     save["Bayes_iterations"] = "";
 
   Plot_list.push_back("Bayes_iterations_data");
     xTitle["Bayes_iterations_data"] = "";
     yTitle["Bayes_iterations_data"] = "#chi^{2}";
     labels["Bayes_iterations_data"] = datafile;
     save["Bayes_iterations_data"] = "";     
              
   Plot_list.push_back("Plot_BinByBin_vs_Gen");
     xTitle["Plot_BinByBin_vs_Gen"] = "E";
     yTitle["Plot_BinByBin_vs_Gen"] = "#frac{Bin-by-bin}{Gen}";
     labels["Plot_BinByBin_vs_Gen"] = label;
     save["Plot_BinByBin_vs_Gen"] = "BinByBin_vs_Gen";

   Plot_list.push_back("PlotCorrectionFactors_BinByBin");
     xTitle["PlotCorrectionFactors_BinByBin"] = "E";
     yTitle["PlotCorrectionFactors_BinByBin"] = "#frac{dN}{dE}";
     labels["PlotCorrectionFactors_BinByBin"] = label;
     save["PlotCorrectionFactors_BinByBin"] = "BinByBin";

   Plot_list.push_back("Plot_JES_vs_E");
     xTitle["Plot_JES_vs_E"] = "E_{det} (GeV)";
     yTitle["Plot_JES_vs_E"] = "JES";
     labels["Plot_JES_vs_E"] = label;
     save["Plot_JES_vs_E"] = "JES";
*/  

/*
   Plot_list.push_back("Unfold_data_Bayes");
     xTitle["Unfold_data_Bayes"] = "E (GeV)";
     yTitle["Unfold_data_Bayes"] = "#frac{dN}{dE}";
     labels["Unfold_data_Bayes"] = datafile;
     save["Unfold_data_Bayes"] = "Unfolded_data_Bayes";

   Plot_list.push_back("Unfold_data_BinByBin");
     xTitle["Unfold_data_BinByBin"] = "E (GeV)";
     yTitle["Unfold_data_BinByBin"] = "#frac{dN}{dE}";
     labels["Unfold_data_BinByBin"] = datafile;
     save["Unfold_data_BinByBin"] = "Unfolded_data_BinByBin";     
*/   
   Plot_list.push_back("Compare_unfolded_data");
     xTitle["Compare_unfolded_data"] = "E (GeV)";
     yTitle["Compare_unfolded_data"] = "#frac{dN}{dE}";
     labels["Compare_unfolded_data"] = datafile;
     save["Compare_unfolded_data"] = "Comparison";      

/*
   Plot_list.push_back("Compare_unfolded_data_all");
     xTitle["Compare_unfolded_data_all"] = "E (GeV)";
     yTitle["Compare_unfolded_data_all"] = "#frac{dN}{dE}";
     labels["Compare_unfolded_data_all"] = datafile;
     save["Compare_unfolded_data_all"] = "Comparison";
     
   Plot_list.push_back("Compare_det_level");
     xTitle["Compare_det_level"] = "E (GeV)";
     yTitle["Compare_det_level"] = "#frac{dN}{dE}";
     labels["Compare_det_level"] = datafile;
     save["Compare_det_level"] = "Comparison";  
*/

  /**********************
  * Unfold data itself. *
  **********************/
   


  cout << "Create dir" << endl;

   /*******************************************
   * Create a new subdirectory for the plots. *
   *******************************************/

   int 	new_dir = mkdir( ("Plots/" + label).Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  );

   /************************************************
   * Prepare canvas, legend, and reused variables. *
   ************************************************/

   TCanvas *can = new TCanvas("Test", "Test", 1500, 1000);   
     can->SetLeftMargin(0.25);
     can->SetTopMargin(0.07);
     can->SetBottomMargin(0.14);
       
   TLegend *legend = new TLegend(0.65, 0.65, 0.99, 0.99);
     legend->SetFillColor( kWhite );

   TH1D* histo;
   TH1D* original;  
   double max_val, min_val;   
   TString drawoptions="";  
   
   // Fun
   void (*current_function)(TH1D*&, TString, TString, TString) = NULL;

      
   /***********************
   * Loop over the plots. *
   ***********************/  

   for( int plot = 0; plot < Plot_list.size(); plot++){  
   
     cout << "// -- FUNCTION\t" <<  Plot_list[plot] << endl;
   
     //-- 1. Assign function to function pointer.

     if	( Plot_list[plot] == "Plot_GenLevel"){ 	
     	current_function = &Plot_GenLevel; }
	
     else if( Plot_list[plot] == "Plot_DetLevel"){ 	
     	current_function = &Plot_DetLevel; }

     else if( Plot_list[plot] == "Plot_GenLevelEta"){	
     	current_function = &Plot_GenLevelEta; }
	
     else if( Plot_list[plot] == "Plot_EtaDiff"){	
     	current_function = &Plot_EtaDiff; }
	
     else if( Plot_list[plot] == "Plot_PhiDiff"){	
     	current_function = &Plot_PhiDiff; }
		
     else if( Plot_list[plot] == "Plot_Distance"){	
     	current_function = &Plot_Distance; }
	
     else if( Plot_list[plot] == "Plot_EnergyIsolation"){
         current_function = &Plot_EnergyIsolation;}
	 
     else if( Plot_list[plot] == "Plot_JetResponse"){
         current_function = &Plot_JetResponse;}

     else if( Plot_list[plot] == "PlotCorrectionFactors_BinByBin"){ 	
     	current_function = &PlotCorrectionFactors_BinByBin; }

     else if( Plot_list[plot] == "PlotCorrectionFactors_Bayes"){	
     	current_function = &PlotCorrectionFactors_Bayes; }  
	
     else if( Plot_list[plot] == "Plot_JES_vs_E"){	
     	current_function = &Plot_JES_vs_E; }  

     else if( Plot_list[plot] == "Plot_BinByBin_vs_Gen"){
        current_function = &Plot_BinByBin_vs_Gen; }

     else if( Plot_list[plot] == "Plot_Bayes_vs_Gen"){
        current_function = &Plot_Bayes_vs_Gen; }
	
     else if( Plot_list[plot] == "Unfold_data_Bayes"){
        current_function = &Unfold_data_Bayes; }     
	
     else if( Plot_list[plot] == "Unfold_data_BinByBin"){
        current_function = &Unfold_data_BinByBin; }	
	
     
	     

     else{
        current_function = NULL; }

 
     //-- Reset variables.
      
     drawoptions =  draw[Plot_list[plot]] ;
     max_val = 0, min_val = 1;
     if( current_function != NULL ){
       for( int file = 0; file < filenames.size(); file++){ 
       
         cout << "// -- FUNCTION - FILE\t" << Plot_list[plot] << "\t" << filenames[file] << endl;
	 cout << "// -- \t\t\t" << filenames[file] << "\t" << fileLabel[file] << "\t" << labels[Plot_list[plot]]
	 << endl;

	 //-- Obtain histogram from file.
         (*current_function)( histo, filenames[file] , fileLabel[file], labels[Plot_list[plot]]);  

	

	 if( normalise[ Plot_list[plot] ] && histo->Integral() != 0.){ histo->Scale( 1./histo->Integral() ); }

	 cout << "// -- FUNCTION - FILE - done" << endl;

	 //-- Adjust min. and max values of the canvas.
         First_Plot( original, histo, file, max_val, min_val);  

	 //-- Color histogram, draw.

         histo->SetLineColor( colours[file] );
	 
         histo->GetXaxis()->SetTitle( xTitle[Plot_list[plot]] );
	 histo->GetXaxis()->SetTitleOffset(0.9);
	 histo->GetXaxis()->SetTitleSize(0.07);
	 histo->GetXaxis()->SetLabelSize(0.07);	 
	 
         histo->GetYaxis()->SetTitle( yTitle[Plot_list[plot]] );
	 histo->GetYaxis()->SetTitleOffset(2.);
	 histo->GetYaxis()->SetTitleSize(0.06);
	 histo->GetYaxis()->SetLabelSize(0.07);
	 
         histo->SetLineStyle( linestyle[file] );	 
	 histo->SetLineWidth(4);
	 
 	 histo->SetMarkerStyle( linestyle[file] );
	 histo->SetMarkerColor( colours[file] );

	 
         can->cd();
	 if( logx[Plot_list[plot] ] ){ can->SetLogx(); }
	 else{	can->SetLogx(0); }
	 if( logy[Plot_list[plot] ] ){ can->SetLogy(); }
	 else{	can->SetLogy(0); }
	 
         histo->Draw( "hist" + drawoptions );
         drawoptions += "same";      
         if( plot == 0 ){ legend->AddEntry( histo, legendEntries[file], "l"); }
	 	 
       }       

       if( original->Integral() != 0. ){
         //-- Finish canvas, save.
	 cout << "Draw this" << endl;
         legend->Draw();

         can->SaveAs( "Plots/" + label + "/" + save[Plot_list[plot]] + ".pdf" );
         can->SaveAs( "Plots/" + label + "/" + save[Plot_list[plot]] + ".C" );
	  cout << "Drawn in Plots/" << label << endl;

         current_function == NULL;
        }
	if( original->Integral() == 0.) { cout << "Shit" << endl; }
       } // Loop over files.
       
      
      /*****************************************************************************
      * The functions that do not produce a series of 1D plots on the same canvas. *
      *****************************************************************************/ 
       
     if( Plot_list[plot] == "Plot_2D_Energy_Response"){
        current_function = &Plot_2D_Energy_Response; 
	cout << "2D responses" << endl;}

     else if( Plot_list[plot] == "Bayes_iterations"){
        current_function = &Bayes_iterations; } 
	
     else if( Plot_list[plot] == "Bayes_iterations_data"){
        current_function = &Bayes_iterations_data; } 	
		
     else if( Plot_list[plot] == "Compare_det_level"){
     	current_function = &Compare_det_level; }  	  

     else if( Plot_list[plot] == "Compare_det_level_all"){
        current_function = &Compare_det_level_all; }	

     //-- Reset variables.
      
     drawoptions = "";
     max_val = 0, min_val = 1.;
     if( current_function != NULL ){
       for( int file = 0; file < filenames.size(); file++){ 

	 //-- Obtain histogram from file.
	 cout << "\t\t\tLabel before\t" << labels[Plot_list[plot]] << endl; 
         (*current_function)( histo, filenames[file] , label, labels[Plot_list[plot]]);  
	 	 
       }       

     } // Loop over files.
     current_function == NULL;      
     
      
   } // Loop over plots.
 

   xtitle["hCastorJet_energy"] = "E_{Det} [GeV]";
   ytitle["hCastorJet_energy"] = "#frac{1}{N}.#frac{dN}{dE_{Det}}";
   htitle["hCastorJet_energy"] = "Detector level jet energy (all)";




   Unfolder my_first_unfolder(filenames, datafile, "20150529", 1);
   my_first_unfolder.LabelPlots( setup + "_Emin_" + Ethresh );
   my_first_unfolder.PrepareLegend( legends );
   my_first_unfolder.PrepareTitles( xtitle, ytitle, htitle );

   my_first_unfolder.Hist_DetLevel();
   my_first_unfolder.Hist_DetLevel_lead();
//   my_first_unfolder.Hist_GenLevel();
//   my_first_unfolder.Hist_getTruth();
//   my_first_unfolder.Plot_Unfolded();
//   my_first_unfolder.Plot_Unfolded_Ratio();
   if( setup == "unfold" )   my_first_unfolder.DoublePaddedComparison("all");
   if( setup == "unfold" )   my_first_unfolder.DoublePaddedComparison("lead");
//   my_first_unfolder.DoublePaddedComparison("hCastorJet_energy");
//   cout << "--- DoublePaddedComparison - hCastorJet_energy" << endl;
//   my_first_unfolder.DoublePaddedComparison("hCastorJet_energy_lead");

//   my_first_unfolder.Hist_DetLevel_DataWithSmear();
//   my_first_unfolder.Hist_DetLevel_MCWithSmear();
   my_first_unfolder.Systematics_CompareDetLevel();
   my_first_unfolder.Systematics_CompareGenLevel();
   my_first_unfolder.CompareGenLevel();
   cout << "--- Systematics" << endl;
 	
   if( setup == "unfold" )   my_first_unfolder.ClosureTest_data("lead");
   if( setup == "unfold" )   my_first_unfolder.ClosureTest_data("all");
   if( setup == "unfold" )   my_first_unfolder.ClosureTest_MC("lead");
   if( setup == "unfold" )   my_first_unfolder.ClosureTest_MC("all");
   if( setup == "unfold" )   my_first_unfolder.ClosureTest_MC_detLevel("lead");
   if( setup == "unfold" )   my_first_unfolder.ClosureTest_MC_detLevel("all");
   if( setup == "unfold" )   cout << "--- ClosureTest" << endl;
       
  return(0); 

}
