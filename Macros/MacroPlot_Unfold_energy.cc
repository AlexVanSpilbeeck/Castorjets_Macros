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
#include "../../../RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"
#include "../../../RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "../../../RooUnfold-1.1.1/src/RooUnfold.h"

#include "color.h"
using namespace std;

int main(){
   gROOT->SetStyle ("Plain");
   gStyle->SetPalette(1);
   gSystem->Load("~/RooUnfold-1.1.1/libRooUnfold.so");

   TString Date = "20141126";
   TString det_radius = "5", gen_radius = "5";

   /* Open DATA and MC */

   TFile *_file0 = TFile::Open( "LoopRootFiles/20141101_Output_JetAnalyzer_deltaR_GEN_ak5_DET_ak5_unsortedE_10054400_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_124_1_Ilx.root", "Read");
   TFile *_file1 = TFile::Open("LoopRootFiles/20141101_Output_JetAnalyzer_Data_ak7_Nevents_419996_CastorTree_data_MinBias_Commissioning10-May19ReReco-v1_7TeV_53XRECOwithCorrector_8_2_cEo.root", "Read");
 
   /* Extract GEN and DET (MC) and DATA distribution */

   RooUnfoldResponse response = _file0.Get("response");
   RooUnfoldResponse response_onlyMatches = (RooUnfoldResponse)_file0->Get("response_onlyMatches");
   TH1D *hCas = (TH1D*)_file1->Get("hCastorJet_energy");	hCas->Sumw2();
   TH1D *hDet = (TH1D*)_file0->Get("hCastorJet_energy");	hDet->Sumw2();
   TH1D *hGen = (TH1D*)_file0->Get("hCastorJet_energy_gen");	hGen->Sumw2();

   int it = 10;

   for(int bin = 0; bin < hDet->GetNbinsX(); bin++){
     double det_value = hDet->GetBinContent(bin),
	    gen_value =  hGen->GetBinContent(bin),
	    data_value = hCas->GetBinContent(bin);
     double det_error = hDet->GetBinError(bin),
            gen_error =  hGen->GetBinError(bin),
            data_error = hCas->GetBinError(bin);

     double original_ratio = gen_value/det_value;

     cout << "Bin\t" << bin << "\t" << original_ratio << endl;

     double binwidth = hDet->GetBinWidth(bin);
/*
     hDet->SetBinContent(bin, det_value/binwidth);	hDet->SetBinError(bin, det_error/binwidth);
     hGen->SetBinContent(bin, gen_value/binwidth);	hGen->SetBinError(bin, gen_error/binwidth);
     hCas->SetBinContent(bin, data_value/binwidth);	hCas->SetBinError(bin, data_error/binwidth);
*/
   }

   hDet->Scale(1./ hDet->Integral() );	hDet->Sumw2();
   hCas->Scale(1./ hCas->Integral() );  hCas->Sumw2();
   hGen->Scale(1./ hGen->Integral() );  hGen->Sumw2();

   int nbins_x = hGen->GetXaxis()->GetNbins();
   double x_axisMax = hGen->GetXaxis()->GetBinUpEdge( nbins_x );
  
   /******************/
   /** Unfold data. **/   
   /******************/

   TCanvas *can = new TCanvas("Generator_level", "Data unfolded", 1.);
   
   TLegend *leg_gen = new TLegend(0.5, 0.6, 0.9, 0.9);
     
   /* Generator level */
     hGen->GetXaxis()->SetTitle("Energy");
     hGen->GetYaxis()->SetTitle("#frac{1}{N}.#frac{dN}{dE}");
     hGen->GetXaxis()->SetRangeUser(200., 3000.);
     hGen->Draw("hist");
     leg_gen->AddEntry(hGen, "Pythia6 Z2*", "p");

   /* Bin-by-bin unfolding. */
   RooUnfoldBinByBin unfold_bin (&response, hCas);
     TH1D* hReco_bin= (TH1D*) unfold_bin.Hreco();
     hReco_bin->SetLineColor(kBlue);
     hReco_bin->SetMarkerColor(kBlue);
     hReco_bin->SetName("Unfolded_binbybin");
     hReco_bin->Draw("datasame");
     leg_gen->AddEntry(hReco_bin, "Data (Bin-by-bin)", "p");

   RooUnfoldBayes unfold_bayes (&response, hCas, it);
   TH1D* hReco_bayes_1= (TH1D*) unfold_bayes.Hreco();
     hReco_bayes_1->SetLineColor(3);
     hReco_bayes_1->SetMarkerColor(3);
     hReco_bayes_1->SetName("Unfolded_bayes");
     hReco_bayes_1->Draw("datasame");


   leg_gen->AddEntry(hReco_bayes_1, TString::Format("Data (Bayes, %i iterations)", it), "p");        
 
   leg_gen->Draw();

   can->SetLogy();

   can->SaveAs("Plots/" + Date + "_Data_unfolded_Det_" + det_radius + "_Gen_" + gen_radius + ".C");
   can->SaveAs("Plots/" + Date + "_Data_unfolded_Det_" + det_radius + "_Gen_" + gen_radius + ".pdf");
  
   /****************/
   /** Unfold MC. **/
   /****************/
  
   TCanvas *can1 = new TCanvas("DET_unfolded", "DET unfolded", 1.);
   
   TLegend *leg_check = new TLegend(0.5, 0.6, 0.9, 0.9);
   
   /* GEN distribution. */
   hGen->SetLineWidth(2);
     hGen->Draw("hist");
     hGen->GetXaxis()->SetRangeUser(200., 3000.);
     leg_check->AddEntry(hGen, "Pythia6 Z2*", "l");
   
   RooUnfoldBinByBin unfold_mc (&response, hDet);
     TH1D* hReco_mc= (TH1D*) unfold_mc.Hreco();
     hReco_mc->SetLineStyle(2);
     hReco_mc->SetLineColor(kBlue);
     hReco_mc->SetLineWidth(2);
     hReco_mc->SetName("Response_binbybin");
     hReco_mc->Draw("histsame");
     leg_check->AddEntry(hReco_mc, "DET --> GEN MC (bin-by-bin)", "l");
   
   RooUnfoldBayes unfold_mc_bayes (&response, hDet, 10);
     TH1D* hReco_mc_bayes= (TH1D*) unfold_mc_bayes.Hreco();
     hReco_mc_bayes->SetLineStyle(3);
     hReco_mc_bayes->SetLineColor(kGreen);
     hReco_mc_bayes->SetLineWidth(2);
     hReco_mc_bayes->SetName("hReco_bayes_1");
     hReco_mc_bayes->Draw("histsame");
     leg_check->AddEntry(hReco_mc_bayes, TString::Format("DET--> GEN MC (Bayes, %i it.)", it), "l");  
   
   leg_check->Draw(); 
   can1->SetLogy();   
   can1->SaveAs("Plots/" + Date + "_MC_unfolded_Det_" + det_radius + "_Gen_" + gen_radius + ".C");
   can1->SaveAs("Plots/" + Date + "_MC_unfolded_Det_" + det_radius + "_Gen_" + gen_radius + ".pdf");



   /******************************/      
   /** Compare data and DET MC. **/
   /******************************/
 
   TCanvas *can2 = new TCanvas("Detector_level", "Detector_level", 1.);   
   TLegend *leg_det = new TLegend(0.5, 0.6, 0.9, 0.9);
   
   /* DET. */
     hDet->GetXaxis()->SetRangeUser(200., 3000.);
//     hDet->GetYaxis()->SetTitle("#frac{1}{N}.#frac{dN}{dE}");
     hDet->GetXaxis()->SetTitle("Energy (GeV)");
     hDet->SetLineColor(kRed);
     hDet->SetMarkerColor(kRed);
     hDet->Draw("hist");
   /* Data. */
     hCas->GetXaxis()->SetTitle("Energy");
//     hCas->GetYaxis()->SetTitle("#frac{1}{N}.#frac{dN}{dE}");
     hCas->Draw("datasame");
     hCas->GetXaxis()->SetRangeUser(200.,3000.);
//     leg_det->AddEntry(hCas, "Data", "l");
 
   leg_det->AddEntry(hCas, "Data", "p");
   leg_det->AddEntry( hDet, "Pythia6 Z2*", "l");
  
   leg_det->Draw();
   can2->SetLogy();

   can2->SaveAs("Plots/" + Date + "_detectorLevel_Det_" + det_radius + "_Gen_" + gen_radius + ".C");
   can2->SaveAs("Plots/" + Date + "_detectorLevel_Det_" + det_radius + "_Gen_" + gen_radius + ".pdf");  


  /*****************************/
  /** Get correction factors. **/
  /*****************************/

  for(int bin = 0; bin < hDet->GetNbinsX(); bin++){

    double detbin = hDet->GetBinContent( bin ),
	   casbin = hCas->GetBinContent( bin ),
           bbbdet = hReco_mc->GetBinContent( bin ),
	   bbbcas = hReco_bin->GetBinContent( bin ),
	   bayesdet = hReco_mc_bayes->GetBinContent( bin ),
	   bayescas_1 = hReco_bayes_1->GetBinContent( bin );
    double detratio_bbb = bbbdet/detbin,
	   casratio_bbb = bbbcas/casbin,
	   detratio_bayes = bayesdet/detbin,
	   casratio_bayes_1 = bayescas_1/casbin;
    cout << "Bin\t" << bin << "\t" << detratio_bbb << "\t" << casratio_bbb << "\t\t\t" << detratio_bayes << "\t" << casratio_bayes_1 << endl; 


  }

  TH1D* htruth = (TH1D*) response.Htruth();
  TCanvas *cat = new TCanvas("canvas_truth", "canvas_truth", 1);
  htruth->Scale( 1./htruth->Integral() );
  htruth->Draw("hist");
  TH1D* hmeas = (TH1D*) response.Hmeasured();
  hmeas->Scale( 1./hmeas->Integral() );
  hmeas->SetLineColor( kGreen + 2);
  hmeas->Draw("histsame");

  hGen->SetLineColor( kRed );
  hGen->SetLineStyle( 2 );
  hGen->Draw("histsame");
  cat->SaveAs("Truth.C");

  return(0);
}
