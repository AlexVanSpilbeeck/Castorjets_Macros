{
   std::vector<TString> filenames;
   std::vector<TString> legendEntries;
   std::vector<TString> colours;
     filenames.push_back("LoopRootFiles/20141024_Output_JetAnalyzer_deltaR_GEN_ak3_DET_ak3_unsortedE_10054400_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_95_1_EfU.root");
     legendEntries.push_back("Gen (ak3) - Det (ak3)");
     colours.push_back(kBlack);

     filenames.push_back("LoopRootFiles/20141024_Output_JetAnalyzer_deltaR_GEN_ak3_DET_ak5_unsortedE_1000000_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_13_1_irJ.root");
     legendEntries.push_back("Gen (ak3) - Det (ak5)");
     colours.push_back(kRed);

     filenames.push_back("LoopRootFiles/20141024_Output_JetAnalyzer_deltaR_GEN_ak3_DET_ak7_unsortedE_1000000_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_13_1_irJ.root");
     legendEntries.push_back("Gen (ak3) - Det (ak7)");
     colours.push_back(kBlue);

     filenames.push_back("LoopRootFiles/20141024_Output_JetAnalyzer_deltaR_GEN_ak5_DET_ak3_unsortedE_10000000_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_173_1_iQO.root");
     legendEntries.push_back("Gen (ak5) - Det (ak3)");
     colours.push_back(kYellow+2);

     filenames.push_back("LoopRootFiles/20141101_Output_JetAnalyzer_deltaR_GEN_ak5_DET_ak5_unsortedE_10054400_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_102_1_g4j.root");
     legendEntries.push_back("Gen (ak5) - Det (ak5)");
     colours.push_back(kGreen+2);

     filenames.push_back("LoopRootFiles/20141024_Output_JetAnalyzer_deltaR_GEN_ak5_DET_ak7_unsortedE_1000000_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_13_1_irJ.root");
     legendEntries.push_back("Gen (ak5) - Det (ak7)");
     colours.push_back(kOrange+7);

     filenames.push_back("LoopRootFiles/20141101_Output_JetAnalyzer_deltaR_GEN_ak7_DET_ak3_unsortedE_10054400_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_54_1_f1Z.root");
     legendEntries.push_back("Gen (ak7) - Det (ak3)");
     colours.push_back(kCyan+2);

     filenames.push_back("LoopRootFiles/20141101_Output_JetAnalyzer_deltaR_GEN_ak7_DET_ak5_unsortedE_10054400_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_125_1_FEG.root");
     legendEntries.push_back("Gen (ak7) - Det (ak5)");
     colours.push_back(kViolet - 2);

     filenames.push_back("LoopRootFiles/20141024_Output_JetAnalyzer_deltaR_GEN_ak7_DET_ak7_unsortedE_1000000_Pythia6_Z2star_Default_varyR_CastorTree_MC_7TeV_42X_53XRECOwithCorrector_13_1_irJ.root");
     legendEntries.push_back("Gen (ak7) - Det (ak7)");
     colours.push_back(kGray+2);
  
   TCanvas *can = new TCanvas("can", "can", 1000, 1000);
   
   TPad *superpad = new TPad("pad", "pad", 0.05, 0.05, 0.9, 0.9);
     superpad->Draw();

   TPad *palettePad = new TPad("palette", "palette", 0.90, 0.05, 0.99, 0.95);
     palettePad->Draw();
 
   TH2D* hMatrix; 
///   can->Divide(3,3);

   for(int file = 0; file < filenames.size(); file++){

     int column = file % 3;
     cout << "File " << file << " Column " << column;
     
     int row = (file - column)/3;
     cout << "Row " << row << endl;

     double xscale = 1.;
     double yscale = 1.;
    
     if( row == 2){ ymin = 0.,	ymax = 0.4; yscale = 3./4.;}
     if( row == 1){ ymin = 0.4,	ymax = 0.7; }
     if( row == 0){ ymin = 0.7,	ymax = 1.0; }
     
     if( column == 0 ){ xmin = 0., xmax = 0.4; xscale = 3./4.;}
     if( column == 1 ){ xmin = 0.4, xmax = 0.7; }
     if( column == 2 ){ xmin = 0.7, xmax = 0.98; }

 
     superpad->cd();
//     can->cd( file+1 );
//     TPad *currentPad = can->GetPad( file+1 );
     TString padname = TString::Format("pad_%i", file);    

     TPad *currentPad = new TPad( padname, padname, 
			xmin, 
			ymin, 
			xmax, 
			ymax );
      
       currentPad->Draw();
       currentPad->cd();
       currentPad->SetLogz();
       currentPad->SetLeftMargin( 0 ); if( column == 0 ){ currentPad->SetLeftMargin( 0.25 ); }
       currentPad->SetRightMargin(0);
       currentPad->SetTopMargin(0);
       currentPad->SetBottomMargin( 0 ); if( row == 2 ){ currentPad->SetBottomMargin( 0.25 ); }
       currentPad->Update();

     TFile *_file0 = TFile::Open( filenames[file], "Read");

     TString histname = TString::Format("hCastorJet_responseMatrix_%i", file) ;
     hMatrix = (TH2D*)_file0->Get("hCastorJet_responseMatrix");
       hMatrix->SetName(histname);
       hMatrix->SetTitle( legendEntries[file] );
       hMatrix->Scale( 1./hMatrix->Integral() );
       hMatrix->SetTitle( legendEntries[file] );

       hMatrix->GetXaxis()->SetTitle("E_{det} (GeV)");
       hMatrix->GetXaxis()->SetTitleOffset(1.0 /xscale);
       hMatrix->GetXaxis()->SetTitleSize(0.08 * xscale);
       hMatrix->GetXaxis()->SetLabelSize(0.06 * xscale);
       hMatrix->GetXaxis()->SetRangeUser(-100., 3000.);
 
       hMatrix->GetYaxis()->SetTitle("E_{gen} (GeV)");
       hMatrix->GetYaxis()->SetTitleOffset(1.2 /yscale);
       hMatrix->GetYaxis()->SetTitleSize(0.08 * yscale);
       hMatrix->GetYaxis()->SetLabelSize(0.06 * yscale);
       hMatrix->GetYaxis()->SetRangeUser(-100., 3000.);

       hMatrix->GetZaxis()->SetRangeUser(1e-8,1.);
       hMatrix->Draw("col");
       
       palette_hist = (TH2D*)hMatrix->Clone("PaletteHistogram");   
//     hMatrix->GetXaxis()->SetLogz();
   }

   superpad->cd();
      
   TLine *line = new TLine(.981,0.10,.981,1.);
   line->Draw();

   can->cd();

   // Add labels about the jet radii to the plot.
   
   TText *text = new TText(0.04,0.25,"Gen (ak7)");
   text->SetTextSize(0.02);
   text->SetTextAngle(90);
   text->Draw();
   text = new TText(0.04,0.5,"Gen (ak5)");
   text->SetTextSize(0.02);
   text->SetTextAngle(90);
   text->Draw();
   text = new TText(0.04,0.75,"Gen (ak3)");
   text->SetTextSize(0.02);
   text->SetTextAngle(90);
   text->Draw();

   text = new TText(0.25,0.97,"Det (ak3)");
   text->SetTextSize(0.02);
   text->Draw();
   text = new TText(0.5,0.97,"Det (ak5)");
   text->SetTextSize(0.02);
   text->Draw();
   text = new TText(0.75,0.97,"Det (ak7)");
   text->SetTextSize(0.02);
   text->Draw();

   // Add a legenda for the colours.

   palettePad->cd();
   palettePad->Range(0.,0.,1.,1.);   
   palettePad->SetLogz();
   TPaletteAxis *palette = new TPaletteAxis(0.10,0.095,0.50, 0.94, hMatrix);
//   TPaletteAxis *palette = (TPaletteAxis*)hMatrix->GetListOfFunctions()->FindObject("palette"); 
   palette->SetLabelSize(0.02);
   palette->SetLabelOffset(0.007);
   palette->SetLabelSize(0.25);
   palette->SetTitleOffset(1);
   palette->SetTitleSize(0.30);

   palette->Draw();

   can->SaveAs("Plots/20141030_Matrix_2.C");
   can->SaveAs("Plots/20141030_Matrix_2.pdf");


}
