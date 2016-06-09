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
   gStyle->SetOptTitle(0);   

   /// Draw all plots on one canvas. ///

   TCanvas *can_comparison = new TCanvas("Canvas_comparison", "Canvas_comparison", 1.);   
     can_comparison->SetLeftMargin(0.14);
     can_comparison->SetBottomMargin(0.14);
   TString drawoptions = "hist";
   TLegend *legend = new TLegend(0.15, 0.7, 0.4, 0.99);
     legend->SetFillColor( kWhite );

   std::vector<TString> filenames;
   std::vector<TString> legendEntries;
   std::vector<int> 	colours;

     filenames.push_back("LoopRootFiles/20141024_Output_JetAnalyzer_deltaR_GEN_ak3_DET_ak3_unsortedE_10054400_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_95_1_EfU.root");
     legendEntries.push_back("Gen (ak3) - Det (ak3)");
     colours.push_back( getColor(1));

     filenames.push_back("LoopRootFiles/20141024_Output_JetAnalyzer_deltaR_GEN_ak3_DET_ak5_unsortedE_1000000_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_13_1_irJ.root");
     legendEntries.push_back("Gen (ak3) - Det (ak5)");
     colours.push_back(getColor(2));

     filenames.push_back("LoopRootFiles/20141024_Output_JetAnalyzer_deltaR_GEN_ak3_DET_ak7_unsortedE_1000000_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_13_1_irJ.root");
     legendEntries.push_back("Gen (ak3) - Det (ak7)");
     colours.push_back(getColor(3));

     filenames.push_back("LoopRootFiles/20141024_Output_JetAnalyzer_deltaR_GEN_ak5_DET_ak3_unsortedE_10000000_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_173_1_iQO.root");
     legendEntries.push_back("Gen (ak5) - Det (ak3)");
     colours.push_back(getColor(4));

     filenames.push_back("LoopRootFiles/20141101_Output_JetAnalyzer_deltaR_GEN_ak5_DET_ak5_unsortedE_10054400_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_102_1_g4j.root");
     legendEntries.push_back("Gen (ak5) - Det (ak5)");
     colours.push_back(getColor(5));

     filenames.push_back("LoopRootFiles/20141024_Output_JetAnalyzer_deltaR_GEN_ak5_DET_ak7_unsortedE_1000000_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_13_1_irJ.root");
     legendEntries.push_back("Gen (ak5) - Det (ak7)");
     colours.push_back(getColor(6));

     filenames.push_back("LoopRootFiles/20141101_Output_JetAnalyzer_deltaR_GEN_ak7_DET_ak3_unsortedE_10054400_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_54_1_f1Z.root");
     legendEntries.push_back("Gen (ak7) - Det (ak3)");
     colours.push_back(getColor(7));

     filenames.push_back("LoopRootFiles/20141101_Output_JetAnalyzer_deltaR_GEN_ak7_DET_ak5_unsortedE_10054400_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_125_1_FEG.root");
     legendEntries.push_back("Gen (ak7) - Det (ak5)");
     colours.push_back(getColor(8));

     filenames.push_back("LoopRootFiles/20141024_Output_JetAnalyzer_deltaR_GEN_ak7_DET_ak7_unsortedE_1000000_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_13_1_irJ.root");
     legendEntries.push_back("Gen (ak7) - Det (ak7)");
     colours.push_back(getColor(9));

   double max_val = 0.;
   TH1D* original;

   for(int file = 0; file < filenames.size(); file++){
     /* Open DATA and MC */

     TFile *_file0 = TFile::Open( filenames[file], "Read");

     /* Extract GEN and DET (MC) and DATA distribution */

     TH1D *hGen = (TH1D*)_file0->Get("hCastorJet_energy_gen");	hGen->Sumw2();	//hGen	 ->Scale( 1./hGen  ->Integral() );
     TH1D *hDet = (TH1D*)_file0->Get("hCastorJet_energy");		hDet->Sumw2();		//hDet->Scale( 1./hDet->Integral() );

     TH1D *hRatio = (TH1D*)hGen->Clone("hRatio");
     if( file == 0 ){ original = hRatio; }
       hRatio->GetXaxis()->SetTitle("E_{jet} (GeV)");
       hRatio->GetYaxis()->SetTitle("Gen/Det ratio");
       hRatio->SetTitle("E_{jet} correction factors");
       hRatio->Divide(hDet);
       hRatio->SetLineWidth(3);
       hRatio->SetLineColor(colours[file]);
//       hRatio->SetLineStyle(file + 1);
       hRatio->Draw( drawoptions );     
       legend->AddEntry( hRatio, legendEntries[file], "l" );       
 
     if( max_val < hRatio->GetMaximum() ){ max_val = hRatio->GetMaximum(); }

     drawoptions = "histsame";
   }

   original->GetYaxis()->SetRangeUser(0.1, max_val*1.1);
   original->GetXaxis()->SetRangeUser(200., 1700.);
   can_comparison->SetLogy();
   legend->Draw();
  
   can_comparison->SaveAs("Plots/20141126_CorrectionFactors.pdf"); 
   can_comparison->SaveAs("Plots/20141126_CorrectionFactors.C");
   
   return(0);
}
