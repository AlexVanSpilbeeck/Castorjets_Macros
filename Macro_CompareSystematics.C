{
   gStyle->SetHistLineWidth(4);
   TString Radius = "";		// Radius of matching criterium. 
   TString Setup = ""; 	// All Castor jets.
   TString Date = "20141124";
//   TString Setup = "allCastorJets";

   /* Open DATA and MC */

   TFile *_file0 = TFile::Open("20141114_Output_JetAnalyzer_GEN_ak5_DET_ak5_unsortedE_100000_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_129_1_D3d.root", "Read");
   TFile *_file1 = TFile::Open("20141114_Output_SystematicsMin_GEN_ak5_DET_ak5_unsortedE_100000_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_129_1_D3d.root", "Read");
   TFile *_file2 =
  
TFile::Open("20141114_Output_SystematicsMax_GEN_ak5_DET_ak5_unsortedE_100000_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_129_1_D3d.root", "Read");
 
   /* Extract GEN and DET (MC) and DATA distribution */

   TH1D *hSysNul = (TH1D*)_file0->Get("hCastorJet_energy");	hSysNul->Sumw2();
   TH1D *hSysMin = (TH1D*)_file1->Get("hCastorJet_energy");	hSysMin->Sumw2();
   TH1D *hSysMax = (TH1D*)_file2->Get("hCastorJet_energy");	hSysMax->Sumw2();
   
   double maxval = 0.;

   /** Canvas for Generator level plots. **/   

   TCanvas *can = new TCanvas("Systematics_comparison", "Systematics_comparison", 800, 1000);
     can->SetLeftMargin(.15);
    
   hSysNul->SetLineWidth(2);
     hSysNul->Draw("hist");
     maxval = hSysNul->GetMaximum();
     
     
   hSysMin->SetLineWidth(2);
     hSysMin->SetLineColor(kRed);
     hSysMin->SetLineStyle(2);
     hSysMin->Draw("histsame");
     if( hSysMin->GetMaximum() > maxval ){ maxval = hSysMin->GetMaximum(); }
     
   hSysMax->SetLineWidth(2);
     hSysMax->SetLineColor(kGreen + 2);
     hSysMax->SetLineStyle(3);
     hSysMax->Draw("histsame");
     if( hSysMax->GetMaximum() > maxval ){ maxval = hSysMax->GetMaximum(); }     
     
   hSysNul->GetYaxis()->SetRangeUser(0.9, maxval*1.1);
   can->Update();

   can->SaveAs("Plots/" + Date + "_Systematics.C");
   can->SaveAs("Plots/" + Date + "_Systematics.pdf");
  

}
