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
using namespace std;

int main(){
   gROOT->SetStyle ("Plain");
   gStyle->SetPalette(1);
   
   std::vector<TString> filenames;
   std::vector<TString> legendEntries;
   std::vector<TString> plotVariables;
     std::vector<TString> plotX;
     std::vector<TString> plotY;
     std::vector<TString> plotTitle;
   std::vector<TString> colours;
   std::vector<double> nEvents;
     filenames.push_back("LoopRootFiles/20141114_Output_JetAnalyzer_GEN_ak5_DET_ak5_unsortedE_100000_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_129_1_D3d.root");
     legendEntries.push_back("MC, No corrections.");
     colours.push_back(getColor(1));
     nEvents.push_back(100000.);

     filenames.push_back("LoopRootFiles/20141114_Output_SystematicsMin_GEN_ak5_DET_ak5_unsortedE_1000000_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_13_1_irJ.root");
     legendEntries.push_back("MC, -22%");
     colours.push_back(getColor(2));
     nEvents.push_back(1000000.);

     filenames.push_back("LoopRootFiles/20141114_Output_SystematicsMax_GEN_ak5_DET_ak5_unsortedE_1000000_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_13_1_irJ.root");
     legendEntries.push_back("MC, +22%");
     colours.push_back(getColor(3));
     nEvents.push_back(1000000.);
     
     // Data. //
     
     filenames.push_back("LoopRootFiles/20141114_Output_JetAnalyzer_Data_ak5_Nevents_419996_CastorTree_data_MinBias_Commissioning10-May19ReReco-v1_7TeV_53XRECOwithCorrector_8_2_cEo.root");
     legendEntries.push_back("Data, no corrections.");
     colours.push_back(getColor(4));
     nEvents.push_back(419996.);   
     
     filenames.push_back("LoopRootFiles/20141114_Output_SystematicsMin_Data_ak5_Nevents_419996_CastorTree_data_MinBias_Commissioning10-May19ReReco-v1_7TeV_53XRECOwithCorrector_8_2_cEo.root");
     legendEntries.push_back("Data, -22%");
     colours.push_back(getColor(5));
     nEvents.push_back(419996.);      
     
     filenames.push_back("LoopRootFiles/20141114_Output_SystematicsMax_Data_ak5_Nevents_419996_CastorTree_data_MinBias_Commissioning10-May19ReReco-v1_7TeV_53XRECOwithCorrector_8_2_cEo.root");
     legendEntries.push_back("Data, +22%");
     colours.push_back(getColor(6));
     nEvents.push_back(419996.);             

     plotVariables.push_back("hCastorJet_energy");
       plotX.push_back("E (GeV)");
       plotY.push_back("#frac{dN}{dE}");
       plotTitle.push_back("Jets per energy (detector level)");

     plotVariables.push_back("hCastorJet_energy_gen");
       plotX.push_back("E (GeV)");
       plotY.push_back("#frac{dN}{dE}");
       plotTitle.push_back("Jets per energy (generator level)");


   double max_val = 0.;
   TH1D* original;
   TString drawoptions = "hist";
   int first_file = 1;

   for(int plot = 0; plot < plotVariables.size(); plot++){

    TString variable = plotVariables[ plot ],
	   htitle =   plotTitle[ plot ],
	   xtitle =   plotX[ plot ],
	   ytitle =   plotY[ plot ];

   /// Draw all plots on one canvas. ///

   TCanvas *can_comparison = new TCanvas("Canvas_" + variable, "Canvas_" + variable, 1.);   
     can_comparison->SetLeftMargin(0.14);
     can_comparison->SetBottomMargin(0.14);

   TLegend *legend = new TLegend(0.6, 0.7, 0.99, 0.99);
     legend->SetFillColor( kWhite );   
   
     for(int file = 0; file < filenames.size(); file++){
       /* Open DATA and MC */

       TFile *_file0 = TFile::Open( filenames[file], "Read");

       /* Extract GEN and DET (MC) and DATA distribution */
//       if( _file0->GetListOfKeys()->Contains( variable ) ){
       if( !(filenames[file].Contains( "Data" ) && variable.Contains("_gen") )){
         TH1D *hJER = (TH1D*)_file0->Get( variable );	hJER->Sumw2();	//hJER	 ->Scale( 1./hJER  ->Integral() );

         if( first_file == 1 ){ original = hJER; first_file = 0;}

         hJER->GetXaxis()->SetTitle( xtitle );
           hJER->GetYaxis()->SetTitle( ytitle );
           hJER->SetTitle( htitle );
           hJER->SetLineWidth(3);
           hJER->SetLineColor(file + 1);
           hJER->SetLineStyle(file + 1);

         hJER->GetXaxis()->SetLabelOffset(0.007);
	   hJER->GetXaxis()->SetLabelSize(0.03);
	   hJER->GetXaxis()->SetTitleSize(0.03);

         hJER->GetYaxis()->SetLabelOffset(0.007);
           hJER->GetYaxis()->SetLabelSize(0.03);
           hJER->GetYaxis()->SetTitleSize(0.03);
           hJER->GetYaxis()->SetTitleOffset(2);

         hJER->Scale(1./ nEvents[file] );
	 
         for(int bin = 0; bin <= hJER->GetNbinsX(); bin++){
           hJER->SetBinContent( bin,
	   			hJER->GetBinContent( bin )/hJER->GetBinWidth( bin ) );
         }

	 

  	 hJER->Draw( drawoptions );     
	 legend->AddEntry( hJER, legendEntries[file], "l" ); 
	 if( max_val < hJER->GetMaximum() ){ max_val = hJER->GetMaximum(); }
	  drawoptions = "histsame";
       }

       
	
      }
      original->GetYaxis()->SetRangeUser(0.9, max_val*1.1);
      original->GetXaxis()->SetRangeUser(200., 2000.);
      can_comparison->Update();
      can_comparison->SetLogy();
      
      drawoptions = "hist";
      legend->Draw();
	  
      can_comparison->SaveAs("Plots/20141126_Systematics" + variable + ".pdf"); 
      can_comparison->SaveAs("Plots/20141126_Systematics" + variable + ".C");
	   
    } // Loop over variables. 
    
   return(0);
}
