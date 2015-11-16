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
   double Eplotmin;
   double deltaPhiMax;
   double etawidth;
   TString match ="";
   TString model_ = "";
   TString modelfilename = "";
 
   std::map<TString, TString> axis_of_interest;
   std::map<TString, TString> legends;
   std::map<TString, TString> printLabel;
   std::map<TString, TString> xtitle, ytitle, htitle;

   std::map<TString, std::map<TString, TString> > set_of_tags;
   std::map<TString, TString> calibration_tag;
   std::map<TString, TString> scalefactors;
   std::map<TString, TString> mc_type;	// Is our MC a model or position dependent source?
  

   bool vary_eta 	= false;
   bool vary_eI 	= false;
   bool compare_radii 	= false;
   bool after_calibration = false;
   bool after_all 	= false;
   bool compare_had 	= false;
   bool compare_em 	= false;
   bool determine_fit 	= false;
   bool systematics_ 	= true;
   double scalefactors_Data_;

   cout << "\n\n\n" << endl;
   cout << "Welcome to the show!" << endl << endl << endl;
 
     cout << "What is the setup?\t";
     cin >> setup;





   /********
   * DATA. * 
   ********/ 

//   TString datafile = "/user/avanspil/Castor_Analysis/ak5_data_" + setup + "_Emin_" + Ethresh + ".000000.root";

   TString MCfile;
   TString datafile;

   if( systematics_ ){

     //datafile = "Histograms/ak5_data_" + setup + "_Emin_" + Ethresh + ".000000.root";
     legends[datafile] = "Data";
     printLabel[datafile] = "data";
     scalefactors_Data_ = 4661641;
     scalefactors[datafile] = TString::Format("%i", static_cast<int>(scalefactors_Data_) );


   /******
   * MC. * 
   ******/ 
   if( setup == "raw" ){
     MCfile = TString::Format("/user/avanspil/Castor_Analysis/ak5ak5_displaced_" + setup + "_Emin_" + Ethresh + ".000000_deltaPhiMax_%f_etaband_%f_" + match + ".root", deltaPhiMax, etawidth );
     filenames.push_back( MCfile );
     legendEntries.push_back("(ak5, ak5) - Displaced");
     colours.push_back(getColor(color_index++));
     linestyle.push_back( style_of_line++);
     fileLabel.push_back("Displaced");    
     printLabel[MCfile] = "Displaced_Pythia6Z2star";
     legends[MCfile] = "Displaced (Pythia6)";
     calibration_tag[MCfile] = "MC";	
     scalefactors[MCfile] = "547922";
   }

    // MCfile for unfolding.
   MCfile = TString::Format("/user/avanspil/Castor_Analysis/ak5ak5_" + modelfilename + "displaced_" + setup + "_Emin_" + Ethresh + ".000000_deltaPhiMax_%f_etaband_%f.root", deltaPhiMax, etawidth );

     //-- Smearing: only 1 file.	 
   if( setup == "smear" ){
     cout << "What is energy threshold?\t";
     cin >> Ethresh;
     cout << "Plot from which Edet value?\t";
     cin >> Eplotmin;
     cout << "Use what delta phi max?\t";
     cin >> deltaPhiMax;
     cout << "Eta band outside of CASTOR?\t";
     cin >> etawidth;
     cout << "What model?\t";
     cin >> model_;

     if( model_ == "p6" ){
       modelfilename = "";
     }
     else if( model_ == "p8" ){
       modelfilename = "Pythia84C_";
     }

     cout << "What matching procedure?\t";
     cin >> match;

     cout << 	"Ethresh\t" << Ethresh <<
		"\nEplotmin\t" << Eplotmin << 
		"\ndeltaPhiMax\t" << deltaPhiMax << 
		"\netawidth\t" << etawidth <<
		"\nmodel_\t" << model_ << 
		"\nmatch\t" << match << endl;

     MCfile = TString::Format("/user/avanspil/Castor_Analysis/ak5ak5_" + modelfilename + "displaced_unfold_Emin_" + Ethresh + ".000000_deltaPhiMax_%f_etaband_%f_" + match + ".root", deltaPhiMax, etawidth );

     if( model_ == "p6" && setup == "smear"){
       cout << "MC file 1" << endl;

       filenames.push_back( MCfile );
       legendEntries.push_back("(ak5, ak5) - Displaced");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("Displaced");    
       printLabel[MCfile] = "Displaced_Pythia6Z2star";
       legends[MCfile] = "Displaced (Pythia6)";
       calibration_tag[MCfile] = "MC";		
	scalefactors[MCfile] = "547922";
     }


     if( model_ == "p8" && setup == "smear"){
       cout << "MC file 2" << endl;

       filenames.push_back( MCfile );
       legendEntries.push_back("(ak5, ak5) - Displaced, Pythia8");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("Displaced");    
       printLabel[MCfile] = "Displaced_Pythia84C";
       legends[MCfile] = "Displaced (Pythia8)";
       calibration_tag[MCfile] = "MC";	
	scalefactors[MCfile] = "670215"; 
     }

     if( model_ == "epos" && setup == "smear"){
       cout << "MC file 3" << endl;

       filenames.push_back( MCfile );
       legendEntries.push_back("(ak5, ak5) - Displaced, EPOS");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("Displaced");    
       printLabel[MCfile] = "Displaced_EPOS";
       legends[MCfile] = "Displaced (Epos)";
       calibration_tag[MCfile] = "MC";	
	scalefactors[MCfile] = "632346"; 
     }

       datafile = "/user/avanspil/Castor_Analysis/ak5_data_unfold_Emin_" + Ethresh + ".000000.root";

   }

   

   //-- Unfolding: multiple files.
   if( setup == "unfold" ){

     cout << "What is energy threshold?\t";
     cin >> Ethresh;
     cout << "Plot from which Edet value?\t";
     cin >> Eplotmin;
     cout << "Use what delta phi max?\t";
     cin >> deltaPhiMax;
     cout << "Eta band outside of CASTOR?\t";
     cin >> etawidth;
     cout << "What matching procedure?\t";
     cin >> match;

     MCfile = TString::Format("/user/avanspil/Castor_Analysis/ak5ak5_displaced_" + setup + "_Emin_" + Ethresh + ".000000_deltaPhiMax_%f_etaband_%f_" + match + ".root", deltaPhiMax, etawidth );
     filenames.push_back( MCfile );
     legendEntries.push_back("(ak5, ak5) - Displaced");
     colours.push_back(getColor(color_index++));
     linestyle.push_back( style_of_line++);
     fileLabel.push_back("Displaced");    
     printLabel[MCfile] = "Displaced_Pythia6Z2star";
     legends[MCfile] = "Displaced (Pythia6)";
     calibration_tag[MCfile] = "MC";	
     scalefactors[MCfile] = "547922";
     mc_type[MCfile] = "actual";

     MCfile = TString::Format("/user/avanspil/Castor_Analysis/ak5ak5_Pythia84C_displaced_" + setup + "_Emin_" + Ethresh + ".000000_deltaPhiMax_%f_etaband_%f_" + match + ".root", deltaPhiMax, etawidth );
     filenames.push_back( MCfile );
     legendEntries.push_back("(ak5, ak5) - Displaced, Pythia8");
     colours.push_back(getColor(color_index++));
     linestyle.push_back( style_of_line++);
     fileLabel.push_back("Displaced");    
     printLabel[MCfile] = "Displaced_Pythia84C";
     legends[MCfile] = "Displaced (Pythia8)";
     calibration_tag[MCfile] = "MC";	
     scalefactors[MCfile] = "670215";
     mc_type[MCfile] = "model";   

     MCfile = TString::Format("/user/avanspil/Castor_Analysis/ak5ak5_EPOS_displaced_" + setup + "_Emin_" + Ethresh + ".000000_deltaPhiMax_%f_etaband_%f_" + match + ".root", deltaPhiMax, etawidth );
     filenames.push_back( MCfile );
     legendEntries.push_back("(ak5, ak5) - Displaced, EPOS");
     colours.push_back(getColor(color_index++));
     linestyle.push_back( style_of_line++);
     fileLabel.push_back("Displaced");    
     printLabel[MCfile] = "Displaced_EPOS";
     legends[MCfile] = "Displaced (EPOS)";
     calibration_tag[MCfile] = "MC";	
     scalefactors[MCfile] = "632346";
     mc_type[MCfile] = "model";   

     MCfile = TString::Format("/user/avanspil/Castor_Analysis/ak5ak5_displaced_down_" + setup + "_Emin_" + Ethresh + ".000000_deltaPhiMax_%f_etaband_%f_" + match + ".root", deltaPhiMax, etawidth );
     filenames.push_back( MCfile );
     legendEntries.push_back("(ak5, ak5) - Displaced down");
     colours.push_back(getColor(color_index++));
     linestyle.push_back( style_of_line++);
     fileLabel.push_back("Displaced_down");    
     printLabel[MCfile] = "Displaced_down";
     legends[MCfile] = "Displaced (down)";
     calibration_tag[MCfile] = "MC";	
     scalefactors[MCfile] = "608839";
     mc_type[MCfile] = "position";   

     MCfile = TString::Format("/user/avanspil/Castor_Analysis/ak5ak5_displaced_up_" + setup + "_Emin_" + Ethresh + ".000000_deltaPhiMax_%f_etaband_%f_" + match + ".root", deltaPhiMax, etawidth );
     filenames.push_back( MCfile );
     legendEntries.push_back("(ak5, ak5) - Displaced up");
     colours.push_back(getColor(color_index++));
     linestyle.push_back( style_of_line++);
     fileLabel.push_back("Displaced_up");    
     printLabel[MCfile] = "Displaced_up";
     legends[MCfile] = "Displaced (up)";
     calibration_tag[MCfile] = "MC";	
     scalefactors[MCfile] = "592822";
     mc_type[MCfile] = "position";   

     datafile = "/user/avanspil/Castor_Analysis/ak5_data_unfold_Emin_" + Ethresh + ".000000.root";
   }

   set_of_tags["calibration_tag"] = calibration_tag;
   set_of_tags["scalefactors"] = scalefactors;
   set_of_tags["legends"] = legends;
   set_of_tags["mc_type"] = mc_type;

    // MCfile for calibration of data.
//     MCfile = TString::Format("/user/avanspil/Castor_Analysis/ak5ak5_" + setup + "_Emin_" + Ethresh + ".000000_deltaPhiMax_%f_etaband_%f.root", deltaPhiMax, etawidth );

    // MC file for calibration of MC.
//     MCfile = TString::Format("/user/avanspil/Castor_Analysis/ak5ak5_FullSimulation_" + setup + "_Emin_" + Ethresh + ".000000_deltaPhiMax_%f_etaband_%f.root", deltaPhiMax, etawidth );



cout << "MCfile\t" << MCfile << endl;

     
/*
    MCfile = "/user/avanspil/Castor_Analysis/ak5ak5_" + setup + "_Emin_" + Ethresh + ".000000.root";
    filenames.push_back( MCfile );
       legendEntries.push_back("(ak5, ak5)");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("Actual");    
       printLabel[MCfile] = "Pythia6Z2star";
       legends[MCfile] = "Pythia6 Z2*";
      calibration_tag[ MCfile ] = "data";		// This tag will be appended to the file containing the calibration factors for data.       
       scalefactors[MCfile] = "12140281";

/*

     MCfile = "/user/avanspil/Castor_Analysis/ak5ak5_FullSimulation_" + setup + "_Emin_" + Ethresh + ".000000.root";
     filenames.push_back( MCfile );
       legendEntries.push_back("(ak5, ak5) - Full Simulation");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("Full simulation");    
       printLabel[MCfile] = "FullSimulation";
       legends[MCfile] = "Full simulation";
       calibration_tag[MCfile] = "MC";     


  

		// This tag will be appended to the file containing the calibation factors for MC.

     MCfile = "/user/avanspil/Castor_Analysis/ak5ak5_displaced_down_" + setup + "_Emin_" + Ethresh + ".000000.root"; 
     filenames.push_back( MCfile );
       legendEntries.push_back("(ak5, ak5) - Disp. (down)");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("Displaced_down");    
       printLabel[MCfile] = "Displaced_down";
       legends[MCfile] = "Disp. (down)";
       calibration_tag[MCfile] = "down";

     MCfile = "/user/avanspil/Castor_Analysis/ak5ak5_displaced_up_" + setup + "_Emin_" + Ethresh + ".000000.root"; 
     filenames.push_back( MCfile );
       legendEntries.push_back("(ak5, ak5) - Disp. (down)");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("Displaced_up");    
       printLabel[MCfile] = "Displaced_up";
       legends[MCfile] = "Disp. (up)";
       calibration_tag[MCfile] = "up";

     MCfile = "/user/avanspil/Castor_Analysis/ak5ak5_Pythia84C_" + setup + "_Emin_" + Ethresh + ".000000.root";
     filenames.push_back( MCfile );	
       legendEntries.push_back("(ak5, ak5) - Pythia8 (4C)");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("Pythia8 (4C)");
       printLabel[MCfile] = "Pythia84C";
       legends[MCfile] = "Pythia8 (4C)";

     MCfile = "/user/avanspil/Castor_Analysis/ak5ak5_Pythia84C_displaced_" + setup + "_Emin_" + Ethresh + ".000000.root";
     filenames.push_back( MCfile );
       legendEntries.push_back("(ak5, ak5) - Pythia8 (4C) - Disp.");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line++);
       fileLabel.push_back("Pythia8 (4C) - Disp.");
       printLabel[MCfile] = "Pythia84C_Displaced";
       legends[MCfile] = "Pythia8 (4C) - Disp.";
*/
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

   Plot_list.push_back("Compare_unfolded_data");
     xTitle["Compare_unfolded_data"] = "E (GeV)";
     yTitle["Compare_unfolded_data"] = "#frac{dN}{dE}";
     labels["Compare_unfolded_data"] = datafile;
     save["Compare_unfolded_data"] = "Comparison";      


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

   time_t theTime = time(NULL);
   struct tm *aTime = localtime(&theTime);

   int day = aTime->tm_mday;
   int month = aTime->tm_mon + 1; // Month is 0 - 11, add 1 to get a jan-dec 1-12 concept
   int year = aTime->tm_year + 1900;

   TString labeltje = TString::Format("%i%i%i_" + setup + "_" + Ethresh + "_%i_deltaPhiMax_0%i_etaband_0%i_" + modelfilename, year, month, day, static_cast<int>(Eplotmin), static_cast<int>(10. * deltaPhiMax), static_cast<int>(10. * etawidth) );

   if( setup == "unfold" ){
     labeltje += "_" + match;
   }
   labeltje += "_" + model_;

   Unfolder my_first_unfolder(	filenames, 
				datafile, 
				set_of_tags, 
				Eplotmin, atof(Ethresh ), 
				deltaPhiMax, 
				etawidth, 
				labeltje , 
				0);
   my_first_unfolder.LabelPlots( 
	TString::Format( setup + "_" + Ethresh + "_%i_deltaPhiMax_0%i_etaband_0%i_" + modelfilename, 
	static_cast<int>(Eplotmin), 
	static_cast<int>(10. * deltaPhiMax), 
	static_cast<int>(10. * etawidth) ) );

   my_first_unfolder.PrepareLegend( legends, printLabel );
   my_first_unfolder.PrepareTitles( xtitle, ytitle, htitle );

   for(int file_ = 0; file_ < filenames.size(); file_++){
     my_first_unfolder.PlotResponseMatrix( file_ );
   }

   my_first_unfolder.SetScaleFactorData( scalefactors_Data_ );


   my_first_unfolder.Hist_DetLevel();
//   my_first_unfolder.Hist_DetLevel_lead();
/*
   my_first_unfolder.Hist_GenLevel();
   my_first_unfolder.Hist_getTruth();
   my_first_unfolder.Plot_Unfolded();
   my_first_unfolder.Plot_Unfolded_Ratio();
   my_first_unfolder.DoublePaddedComparison("hCastorJet_energy_lead");
*/

  // ISOLATED
   if( setup == "isolated"){ 

//     my_first_unfolder.CalculateSystematics("20150729_test", 1, 16);
//     my_first_unfolder.CalibrationFactors_oneCanvas(true);
//     my_first_unfolder.CalibrationFactors_oneCanvas(false);
/*
     my_first_unfolder.CalculateSystematics("all"); 

     my_first_unfolder.CalculateSystematics("one", 12);
     my_first_unfolder.CalculateSystematics("one", 13);
     my_first_unfolder.CalculateSystematics("one", 14);
     my_first_unfolder.CalculateSystematics("one", 15);
*/
//     my_first_unfolder.Plot_Calibrated_functions();
     // --  Creates a canvas with the calibration functions in absolute and ratio format.

/*
     my_first_unfolder.CalculateSystematics("all_sectors",1,1);
*/
     my_first_unfolder.CalculateSystematics("good_sectors",1,1);
     my_first_unfolder.CalculateSystematics("separate_sectors", 13, 14);



/*     TGraphErrors* gre;

     my_first_unfolder.SetFitDraw( true );
     my_first_unfolder.CalibrationFunction(1, 1, gre); 

     my_first_unfolder.CalibrationFunction_sectors(13, 14, -1, gre);   
*/
   } // ISOLATED.

   
   if( setup == "smear" ){

//     my_first_unfolder.SetCastorJetEnergy_norm( static_cast<double>(547922)/ static_cast<double>(4661641) );

//     my_first_unfolder.Hist_GenLevel("lead");
//     my_first_unfolder.Hist_GenLevel("all");
     my_first_unfolder.Hist_DetLevel();

/*
     for(int file_ = 0; file_ < filenames.size(); file_++){ 
	// -- Response matrix
	my_first_unfolder.PlotFromResponseMatrix(file_, "Edet", "Egen"); 
	// -- Edet vs. fakes
	my_first_unfolder.PlotFromResponseMatrix(file_, "Edet", "nFake"); 
	// -- Edet vs. misses
	my_first_unfolder.PlotFromResponseMatrix(file_, "Edet", "nMiss"); 
	// -- Egen vs. fakes
	my_first_unfolder.PlotFromResponseMatrix(file_, "Egen", "nFake"); 
	// -- Egen vs. misses.
	my_first_unfolder.PlotFromResponseMatrix(file_, "Egen", "nMiss"); 
	// -- Misses vs. fakes.
	my_first_unfolder.PlotFromResponseMatrix(file_, "nMiss", "nFake"); 
	// -- nDet vs. fakes.
	my_first_unfolder.PlotFromResponseMatrix(file_, "nDet", "nFake"); 
	// -- nGen vs. misses.
	my_first_unfolder.PlotFromResponseMatrix(file_, "nGen", "nMiss"); 

	my_first_unfolder.PlotFromResponseMatrix(file_, "Egen", "nMatch"); 

	my_first_unfolder.PlotFromResponseMatrix(file_, "Edet", "nMatch"); 

	my_first_unfolder.PlotFromResponseMatrix(file_, "nMiss", "nMatch"); 

	my_first_unfolder.PlotFromResponseMatrix(file_, "nFake", "nMatch"); 

	my_first_unfolder.PlotFromResponseMatrix(file_, "nFake");

	my_first_unfolder.PlotFromResponseMatrix(file_, "nMiss");

	my_first_unfolder.PlotFromResponseMatrix(file_, "Edet");

	my_first_unfolder.PlotFromResponseMatrix(file_, "Egen");
	
    }
*/

/*
     my_first_unfolder.DoublePaddedComparison("all");
     my_first_unfolder.DoublePaddedComparison("lead");
     my_first_unfolder.DoublePaddedComparison_unfolding("all", 6);
*/
//     my_first_unfolder.DoublePaddedComparison_unfolding("lead", 6);


//     my_first_unfolder.Unfolding_data("all", 20);


     // -- Closure test on data and MC.	

//     my_first_unfolder.ClosureTest_data("lead", "Displaced");

     // Closure Test - method 1.
     my_first_unfolder.SetSubhistogram_cut( Eplotmin );

//     my_first_unfolder.PlotStartingDistributions();

//     my_first_unfolder.Smear_gen(0);

     //my_first_unfolder.SetAddLabel( "findingSolutions" )
     my_first_unfolder.ClosureTest_data("all", "Displaced", 1);


//     my_first_unfolder.ClosureTest_MC_detLevel("all");

     // Closure Test - method 2.
//     my_first_unfolder.ClosureTest_data("all", "Displaced", 2);
/*
     my_first_unfolder.PlotStartingDistributions_comparingEmin("fake");
     my_first_unfolder.PlotStartingDistributions_comparingEmin("miss");
     my_first_unfolder.PlotStartingDistributions_comparingEmin("detector");
*/
//     my_first_unfolder.PlotStartingDistributions_comparingEmin("generator");

/*
     my_first_unfolder.Dissect_ResponseObject();


//     my_first_unfolder.Calculate_smearedBackError(-1, 0);



/*
     my_first_unfolder.ClosureTest_MC("all");
     my_first_unfolder.ClosureTest_MC_detLevel("all");

     // -- CHI2 tests.

     my_first_unfolder.Chi2_comparison_data("all");
     my_first_unfolder.Chi2_comparison_MC("all");
     my_first_unfolder.Chi2_comparison_MC_det("all");
*/
   }


   if( setup == "unfold"){

     //-- Comparison of unfolded distributions.


     	//-- Extract the response object per file and plot detector and generator level distributions.
//    my_first_unfolder.Plot_DistributionsResponseObject(1);	
//    my_first_unfolder.Plot_DistributionsResponseObject(-1);	

	//-- Unfold with a number of Bayesian iterations and check the chi2 between two consecutive iterations.
	//-- This is Delta Chi2.
//     my_first_unfolder.Chi2diff_test_data("all", "Displaced", 2);


	//-- Plot the different model samples, calculate the average and plot the absolutes and ratios versus the average.
//     my_first_unfolder.DoublePaddedComparison_modelDependence("all", 30);
//     my_first_unfolder.DoublePaddedComparison_modelDependence("all", 50);

	//-- Plot the actual sample and the samples containing the uncertainty on the position.
	//-- Plot the absolute and plot ratio versus the actual sample.
//     my_first_unfolder.DoublePaddedComparison_positionDependence("all", 30);
//     my_first_unfolder.DoublePaddedComparison_positionDependence("all", 50);
/*
	//-- Plot all systematics in absolute and ratio.
     TCanvas* can_30it;
     my_first_unfolder.Plot_Unfolded_Ratio_allSystematics( can_30it, "all", 30);

     TCanvas* can_50it;
     my_first_unfolder.Plot_Unfolded_Ratio_allSystematics( can_50it, "all", 50);
*/

//     my_first_unfolder.DoublePaddedComparison_unfolding("all", 30);	//
//     my_first_unfolder.DoublePaddedComparison_unfolding("all", 50);	/

///     my_first_unfolder.Plot_Unfolded();


/*
     my_first_unfolder.DoublePaddedComparison_statistics("all", 12);
     my_first_unfolder.DoublePaddedComparison_statistics("all", 50);
/*
     my_first_unfolder.Hist_DetLevel();
*/
/*
     my_first_unfolder.PlotStartingDistributions_MCfiles("fake");
     my_first_unfolder.PlotStartingDistributions_MCfiles("miss");
     my_first_unfolder.PlotStartingDistributions_MCfiles("detector");

     my_first_unfolder.PlotStartingDistributions_comparingEmin("fake");
     my_first_unfolder.PlotStartingDistributions_comparingEmin("miss");
     my_first_unfolder.PlotStartingDistributions_comparingEmin("detector");
*/
/*
     TCanvas *can_modeldependence_12it = new TCanvas("Model_dependence_12its", "Model_dependence_12its", 1.);
     TPad* pad = new TPad("pad", "pad", 0.,0.,1.,1.);
     pad->Draw();
     my_first_unfolder.Plot_Unfolded_Ratio_modelDependence( pad , "all", 12);
     can_modeldependence_12it->SaveAs("Modeldependence_12it.C");
     can_modeldependence_12it->SaveAs("Modeldependence_12it.pdf");

     my_first_unfolder.Plot_Unfolded_Ratio_positionDependence( pad , "all", 12);
     can_modeldependence_12it->SaveAs("Positiondependence_12it.C");
     can_modeldependence_12it->SaveAs("Positiondependence_12it.pdf");

     TCanvas* can_12it;

//     can_12it->SaveAs("Totaldependence_12it.C");
//     can_12it->SaveAs("Totaldependence_12it.pdf");

     cout << "Off to 50" << endl;

     TCanvas* can_50it;   
     my_first_unfolder.Plot_Unfolded_Ratio_allSystematics( can_50it , "all", 50);
*/
   }

  return(0); 

}

