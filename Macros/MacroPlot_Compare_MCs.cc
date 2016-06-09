{
  TCanvas *can = new TCanvas("canvas", "canvas", 800, 600);
    can->SetLeftMargin(0.14);
    can->SetBottomMargin(0.14);

  TPad *pad1 = new TPad("Pad1", "Pad1", 0., 0.3, 0.57, 1.);
    pad1->SetLeftMargin(0.25);
    pad1->SetRightMargin(0.);
    pad1->SetTopMargin(0.05);
    pad1->SetBottomMargin(0.);
    pad1->Draw();

    can->cd();

  TPad *pad1_ratio = new TPad("Pad1_ratio", "Pad1_ratio", 0., 0., 0.57, .3);
    pad1_ratio->SetLeftMargin(0.25);
    pad1_ratio->SetRightMargin(0.);
    pad1_ratio->SetTopMargin(0.);
    pad1_ratio->SetBottomMargin(0.2);
    pad1_ratio->Range(-50,-0.1852941,5451.852,2.5);
    pad1_ratio->Draw();
  
    can->cd();

  TPad *pad2 = new TPad("Pad2", "Pad2", 0.57, 0.3, 1., 1.);
    pad2->SetRightMargin(0.05);
    pad2->SetLeftMargin(0.);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.);
    pad2->Draw();

    can->cd();

  TPad *pad2_ratio = new TPad("Pad2_ratio", "Pad2_ratio", 0.57, 0., 1., .3);
    pad2_ratio->SetLeftMargin(0.);
    pad2_ratio->SetRightMargin(0.05);
    pad2_ratio->SetTopMargin(0.);
    pad2_ratio->SetBottomMargin(0.2);
    pad2_ratio->Range(-50,-0.1852941,5451.852,2.5);
    pad2_ratio->Draw();

    can->cd();

//  TPad *padLeg = new TPad("Pad3", "Pad3", 0.9,0.6, 1., 1.);
//    padLeg->Draw();

    can->cd();

   TLegend *legend = new TLegend(0.7, 0.7, 1., .9);
     legend->SetFillColor( kWhite );

   std::vector<TString> filenames;
   std::vector<TString> legendEntries;
   std::vector<TString> colours;

   legendEntries.push_back("Data");
   legendEntries.push_back("Pythia6 Z2*");
/*
   filenames.push_back("LoopRootFiles/20141101_Output_JetAnalyzer_deltaR_GEN_ak5_DET_ak5_unsortedE_10054400_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_102_1_g4j.root");
     legendEntries.push_back("Gen (ak5) - Det (ak5)");
     colours.push_back(kBlack);
*/
   filenames.push_back("Rivet_files/epos_500k.root");
     legendEntries.push_back("epos");
     colours.push_back(kGreen+2);

   filenames.push_back("Rivet_files/herwig_qcd_ee3c_1M.root");
     legendEntries.push_back("Herwig QCD");
     colours.push_back(kBlue);

   filenames.push_back("Rivet_files/py6_cuet_5M.root");
     legendEntries.push_back("py6_cuet");
     colours.push_back(kCyan);

   filenames.push_back("Rivet_files/py6_z2star_1M.root");
     legendEntries.push_back("py6_Z2*");
     colours.push_back(kYellow+4);

   filenames.push_back("Rivet_files/py8_4c_1M.root");
     legendEntries.push_back("py8_4c");
     colours.push_back(kOrange+2);

   filenames.push_back("Rivet_files/py8_cuet_1M.root");
     legendEntries.push_back("py8_cuet");
     colours.push_back(kGray);
   
   filenames.push_back("Rivet_files/py8_monash_1M.root");
     legendEntries.push_back("py8_monash");
     colours.push_back(kPink+9);

   filenames.push_back("Rivet_files/qgsjet_1M.root");
     legendEntries.push_back("qgsjet");
     colours.push_back(kMagenta);

   /* Unfold data. */

   TFile *fMC = new TFile("LoopRootFiles/20141101_Output_JetAnalyzer_deltaR_GEN_ak5_DET_ak5_unsortedE_10054400_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_102_1_g4j.root", "Read");
  
   TFile *fData = new TFile("LoopRootFiles/20141101_Output_JetAnalyzer_Data_ak7_Nevents_419996_CastorTree_data_MinBias_Commissioning10-May19ReReco-v1_7TeV_53XRECOwithCorrector_8_2_cEo.root", "Read");

   RooUnfoldResponse response = (RooUnfoldResponse)fMC.Get("response");
   TH1D *hGen = (TH1D*)fMC->Get("hCastorJet_energy_gen");	hGen->Sumw2();	hGen->Scale(1./hGen->Integral() );
   TH1D *hCas = (TH1D*)fData->Get("hCastorJet_energy");		hCas->Sumw2();	hCas->Scale(1./hCas->Integral() );

   for(int bins = 0; bins < hGen->GetNbinsX(); bins++){
     double binold_gen = hGen->GetBinContent( bins );
     double binwdt_gen = hGen->GetBinWidth( bins );
     double binnew_gen = binold_gen/binwdt_gen;
     hGen->SetBinContent( bins, binnew_gen );

     double binold_det = hCas->GetBinContent( bins );
     double binwdt_det = hCas->GetBinWidth( bins );
     double binnew_det = binold_det/binwdt_det;
     hCas->SetBinContent( bins, binnew_det );
   }

    hGen->Sumw2();  hGen->Scale(1./hGen->Integral() );
    hCas->Sumw2();  hCas->Scale(1./hCas->Integral() );


   RooUnfoldBayes unfold_bayes (&response, hCas, 4);
     TH1D* hData= (TH1D*) unfold_bayes.Hreco();
     hGen->GetYaxis()->SetRangeUser(1e-8, 1.);

   legend->AddEntry( hData, legendEntries[0] + " - Bayes (4)", "p");
   legend->AddEntry( hGen, legendEntries[1], "l");

   TH1D *hData_clone = (TH1D*)hData->Clone("Divide");
   hData_clone->Divide(hData);

   Pad1_ratio->cd();
     hData_clone->GetXaxis()->SetTitleSize(0.1);
     hData_clone->GetXaxis()->SetTitleOffset(0.9);
     hData_clone->GetXaxis()->SetTitle("E (GeV)");
     hData_clone->GetXaxis()->SetLabelSize(0.1);
     hData_clone->GetYaxis()->SetTitle("Ratio");
     hData_clone->GetYaxis()->CenterTitle(true);
     hData_clone->GetYaxis()->SetTitleSize(0.08);
     hData_clone->GetYaxis()->SetLabelSize(0.1);
     hData_clone->Draw("hist");
//     hData_clone->Draw("adata");

   Pad2_ratio->cd();
     hData_clone->Draw("hist");
//     hData_clone->Draw("data");

   // Make sure to draw them on both pads.

   Pad1->cd();
     hGen->SetLineColor(kBlack);
     hGen->SetLineWidth(2);
     hGen->GetXaxis()->SetTitle("E (GeV)");
     hGen->GetXaxis()->SetTitleOffset(1.);
//     hGen->GetXaxis()->SetTitleSize(1.);
     hGen->GetYaxis()->SetTitle("#frac{1}{N}.#frac{dN}{dE}");  
     hGen->GetYaxis()->SetTitleSize(0.04);
     hGen->GetYaxis()->SetTitleOffset(1.5);
     hGen->DrawCopy("hist");

     hData->SetMarkerColor(kRed);
     hData->SetLineColor(kRed);
     hData->DrawCopy("datasame");

   Pad2->cd();
   hGen->DrawCopy("hist");
   hData->DrawCopy("datasame");

   /* Data unfolded. */
   can->cd();

   for(int file = 2; file < filenames.size(); file++){
     cout << "First or second pad\t" << static_cast<double>(filenames.size())/2. << "\t" << file;

     if( static_cast<double>(filenames.size())/2. >= file ){ Pad1->cd(); cout << "\tP1";}
     else{ Pad2->cd(); cout << "\tP2";}

     cout << "\t" << filenames[file] << "\t" << legendEntries[file] << endl;

     TFile *_file0 = TFile::Open( filenames[file], "Read");
     TH1D* hRivet = (TH1D*)_file0->Get("hCastorJet_energy_Rivet");
     hRivet->Scale(1./hRivet->Integral() );
     if(file != 5){ hRivet->SetLineColor( file ); } else{ hRivet->SetLineColor( kOrange ); }
     hRivet->SetLineStyle( file );
     hRivet->SetLineWidth( 3 );
     hRivet->DrawCopy("histsame");
     legend->AddEntry(hRivet, legendEntries[file], "l");

     if( static_cast<double>(filenames.size())/2. >= file ){ Pad1_ratio->cd(); }
     else{ Pad2_ratio->cd(); }

     hRivet->Divide(hData);
     hRivet->DrawCopy("histsame");
  }   

   Pad1->cd();
   legend->Draw();

   Pad1->SetLogy();
   Pad1_ratio->SetLogy();
   pad1_ratio->Range(-50.,0.01,5451,4.);
   pad1_ratio->Update();

   Pad2->SetLogy();
   Pad2_ratio->SetLogy();
//   pad2_ratio->Range(-50,0.01, 5451.852,4.);


   can->SaveAs("Plots/20141126_CompareMCs.pdf");
   can->SaveAs("Plots/20141126_CompareMCs.C");


   TCanvas *canz2 = new TCanvas("canz2", "canz2", 1.);

     TFile *_file0 = TFile::Open( "Rivet_files/py6_z2star_1M.root", "Read");
     TH1D* hRivet = (TH1D*)_file0->Get("hCastorJet_energy_Rivet");
     hRivet->Scale( 1./hRivet->Integral() );
     hGen->Divide( hRivet );
     hGen->Draw("hist");
}
