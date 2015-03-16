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

//#ifndef FUNCTION_FIRSTPLOT_H
//#include "Function_FirstPlot.h"
//#endif
//#include "Function_make_Tex.h"

using namespace std;

#define determine_fit false

int main(){
   gROOT->SetStyle ("Plain");
   gStyle->SetPalette(1);

   TString label = "20150312_test";

   std::vector<TString> filenames;
   std::vector<TString> data_files;

   std::vector<TString> legendEntries;
   std::vector<TString> plotVariables;
   std::vector<TString> jetSelection;
     std::vector<TString> plotX;
     std::vector<TString> plotY;
     std::vector<TString> plotTitle;
     std::vector<bool> project_2D;
   std::vector<int>     colours;
   std::vector<int>	linestyle;
   std::vector<TString> Plot_list;
   
   int color_index = 2;
   
   std::map<TString, TString> axis_of_interest;

   bool vary_eta = true;
   bool vary_eI = false;

   cout << "Files" << endl;

   /********
   * DATA. * 
   ********/ 
/*
   data_files.push_back("LoopRootFiles/20150312_stripped_tree_DATA.root");
     legendEntries.push_back("Data");
*/   
   if( determine_fit ){
   
     filenames.push_back("LoopRootFiles/20150305_had_Output_JetAnalyzer_radii_strippedTree_GEN_ak5_DET_ak5_margin_0.500000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.4");
       colours.push_back(getColor(color_index++));      
   }

   else if( vary_eI && vary_eta){
     int style_of_line = 1; 

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
     int style_of_line = 2;

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.400000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.4, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.500000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.5, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.600000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.6, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.700000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.7, no E_{I} cut");
       linestyle.push_back( style_of_line );
       colours.push_back(getColor(color_index++));

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.800000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.8, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_0.900000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 0.9, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.000000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 1.0, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.100000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 1.1, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.200000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 1.2, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.300000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 1.3, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );

     filenames.push_back("LoopRootFiles/20150311_Calib_noIso_had_Output_JetAnalyzer_radii_strippedTree_etaband_1.400000_12137851_0_sectors.root");
       legendEntries.push_back("#eta band: 1.4, no E_{I} cut");
       colours.push_back(getColor(color_index++));
       linestyle.push_back( style_of_line );
   }   

   /**************************************
   * Prepare the plots and their labels. * 
   **************************************/	 

   std::map<TString, TString> xTitle;
   std::map<TString, TString> yTitle; 
   std::map<TString, TString> save; 
   std::map<TString, TString> labels;
/*
   Plot_list.push_back("Plot_GenLevel");
     xTitle["Plot_GenLevel"] = "E_{gen} (GeV)";
     yTitle["Plot_GenLevel"] = "#frac{dN}{dE_{gen}}";
     labels["Plot_GenLevel"] = label;
     save["Plot_GenLevel"] = "Gen_energy";    

   Plot_list.push_back("Plot_DetLevel");
     xTitle["Plot_DetLevel"] = "E_{det} (GeV)";
     yTitle["Plot_DetLevel"] = "#frac{dN}{dE_{det}}";
     labels["Plot_GenLevel"] = label;
     save["Plot_DetLevel"] = "Det_energy";

   Plot_list.push_back("Plot_GenLevelEta");
     xTitle["Plot_GenLevelEta"] = "#Delta#eta";
     yTitle["Plot_GenLevelEta"] = "#frac{dN}{#Delta#Eta}";
     labels["Plot_GenLevelEta"] = label;
     save["Plot_GenLevelEta"] = "Gen_eta";

   Plot_list.push_back("PlotCorrectionFactors_BinByBin");
     xTitle["PlotCorrectionFactors_BinByBin"] = "E";
     yTitle["PlotCorrectionFactors_BinByBin"] = "#frac{dN}{dE}";
     labels["PlotCorrectionFactors_BinByBin"] = label;
     save["PlotCorrectionFactors_BinByBin"] = "BinByBin";
*/
   Plot_list.push_back("PlotCorrectionFactors_Bayes");
     xTitle["PlotCorrectionFactors_Bayes"] = "E (GeV)";
     yTitle["PlotCorrectionFactors_Bayes"] = "#frac{dN}{dE}";
     labels["PlotCorrectionFactors_Bayes"] = label;
     save["PlotCorrectionFactors_Bayes"] = "Bayes";

   Plot_list.push_back("Plot_Bayes_vs_Gen");
     xTitle["Plot_Bayes_vs_Gen"] = "E";
     yTitle["Plot_Bayes_vs_Gen"] = "#frac{Bayes}{Gen}";
     labels["Plot_Bayes_vs_Gen"] = label;
     save["Plot_Bayes_vs_Gen"] = "Bayes_vs_Gen";
     
/*
   Plot_list.push_back("Plot_BinByBin_vs_Gen");
     xTitle["Plot_BinByBin_vs_Gen"] = "E";
     yTitle["Plot_BinByBin_vs_Gen"] = "#frac{Bin-by-bin}{Gen}";
     labels["Plot_BinByBin_vs_Gen"] = label;
     save["Plot_BinByBin_vs_Gen"] = "BinByBin_vs_Gen";

//   Plot_list.push_back("Plot_JES_vs_E");
     xTitle["Plot_JES_vs_E"] = "E_{det} (GeV)";
     yTitle["Plot_JES_vs_E"] = "JES";
     labels["Plot_JES_vs_E"] = label;
     save["Plot_JES_vs_E"] = "JES_vs_E";

   Plot_list.push_back("Plot_2D_Energy_Response");
     xTitle["Plot_2D_Energy_response"] = "";
     yTitle["Plot_2D_Energy_response"] = "";
     labels["Plot_2D_Energy_response"] = label;
     save["Plot_2D_Energy_response"] = "";
*/
   Plot_list.push_back("Bayes_iterations");
     xTitle["Bayes_iterations"] = "";
     yTitle["Bayes_iterations"] = "#chi^{2}";
     labels["Bayes_iterations"] = "15";
     save["Bayes_iterations"] = "";



   /*******************************************
   * Create a new subdirectory for the plots. *
   *******************************************/

   int 	new_dir = mkdir( ("Plots/" + label).Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  );

   /************************************************
   * Prepare canvas, legend, and reused variables. *
   ************************************************/

   TCanvas *can = new TCanvas("Test", "Test", 1.);   
       
   TLegend *legend = new TLegend(0.65, 0.65, 0.99, 0.99);
     legend->SetFillColor( kWhite );

   TH1D* histo;
   TH1D* original;  
   double max_val, min_val;   
   TString drawoptions="";  
   
   void (*current_function)(TH1D*&, TString, TString, TString) = NULL;   
      
   /***********************
   * Loop over the plots. *
   ***********************/  

   for( int plot = 0; plot < Plot_list.size(); plot++){   
   
     //-- 1. Assign function to function pointer.

     if	( Plot_list[plot] == "Plot_GenLevel"){ 	
     	current_function = &Plot_GenLevel; }
	
     else if( Plot_list[plot] == "Plot_DetLevel"){ 	
     	current_function = &Plot_DetLevel; }

     else if( Plot_list[plot] == "Plot_GenLevelEta"){	
     	current_function = &Plot_GenLevelEta; }

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

     else{
        current_function = NULL; }
 
     //-- Reset variables.
      
     drawoptions = "";
     max_val = 0, min_val = 0;
     if( current_function != NULL ){
       for( int file = 0; file < filenames.size(); file++){ 

	 //-- Obtain histogram from file.
	 cout << "\t\t\tLabel before\t" << labels[Plot_list[plot]] << endl; 
         (*current_function)( histo, filenames[file] , legendEntries[file], labels[Plot_list[plot]]);  

	 //-- Adjust min. and max values of the canvas.

         First_Plot( original, histo, file, max_val, min_val);  

	 //-- Color histogram, draw.

         histo->SetLineColor( colours[file] );	 
         histo->GetXaxis()->SetTitle( xTitle[Plot_list[plot]] );
         histo->GetYaxis()->SetTitle( yTitle[Plot_list[plot]] );
         histo->SetLineStyle( linestyle[file] );

         can->cd();
         histo->Draw( "hist" + drawoptions );
         drawoptions = "same";      
         if( plot == 0 ){ legend->AddEntry( histo, legendEntries[file], "l"); }
	 
       }       

       if( original->Integral() != 0. ){
         //-- Finish canvas, save.
         legend->Draw();

         can->SaveAs( "Plots/" + label + "/" + save[Plot_list[plot]] + ".pdf" );
         can->SaveAs( "Plots/" + label + "/" + save[Plot_list[plot]] + ".C" );

         current_function == NULL;
        }
       } // Loop over files.
       
      
      /*****************************************************************************
      * The functions that do not produce a series of 1D plots on the same canvas. *
      *****************************************************************************/ 
       
     if( Plot_list[plot] == "Plot_2D_Energy_Response"){
        current_function = &Plot_2D_Energy_Response; }

     else if( Plot_list[plot] == "Bayes_iterations"){
        current_function = &Bayes_iterations; }    
	
     //-- Reset variables.
      
     drawoptions = "";
     max_val = 0, min_val = 0;
     if( current_function != NULL ){
       for( int file = 0; file < filenames.size(); file++){ 

	 //-- Obtain histogram from file.
	 cout << "\t\t\tLabel before\t" << labels[Plot_list[plot]] << endl; 
         (*current_function)( histo, filenames[file] , legendEntries[file], labels[Plot_list[plot]]);  
	 	 
       }       

     } // Loop over files.
     current_function == NULL;       
   } // Loop over plots.
       
  return(0); 
}
