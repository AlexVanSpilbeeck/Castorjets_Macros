{
  gROOT->ProcessLine(".L Function_PrepareCanvas.h");
  gROOT->ProcessLine(".L Function_Prepare1Dplot.h");
  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();

  TF1* gamma = new TF1("gamma", "1./sqrt( 1. - x)", 0., 0.999999999999) ;
  TCanvas *can;
  PrepareCanvas(can, "Gamma_variable" );

  can->cd();
  gamma->GetXaxis()->SetTitle("v/c");
  gamma->GetYaxis()->SetTitle("#gamma");
  gamma->SetLineWidth(3);
  Prepare_1Dplot( gamma );

  gamma->Draw("");

  can->SaveAs("Relativity_gamma.pdf");


  TF1* gamma_speed = new TF1("gamma_speed", "x/pow(sqrt( 1. - x),3.)", 0., 0.999999999999) ;
  gamma_speed ->GetXaxis()->SetTitle("v/c");
  gamma_speed ->GetYaxis()->SetTitle("#frac{d#gamma}{dv}");
  gamma_speed ->SetLineWidth(3);
  Prepare_1Dplot( gamma_speed );  
  
  gamma_speed->Draw("");
  can->SetLogy();
  can->SaveAs("Relativity_gamma_change.pdf");
  
}
