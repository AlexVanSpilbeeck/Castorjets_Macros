

 //STANDARD ROOT INCLUDES
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>
#include <TColor.h>
#include <TString.h>

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


#include "../color.h"

#include "Function_PrepareCanvas.h"
#include "Function_Prepare1Dplot.h"



void Plot_Systematics( TH1D* hData, vector<TH1D*> vSyst, TString name ){

  //== Prepare the canvas.
  TCanvas* can_;
  PrepareCanvas( can_, "Systematics_" + name );
  TPad *pad_abs_, *pad_ratio_;
  SplitCanvas( can_, pad_abs_, pad_ratio_ );
  
  //== Legend
  double ticksize = gStyle->GetTickLength();
  TLegend *leg = new TLegend( 0.5, 0.7, 1. - pad_abs_->GetRightMargin() - ticksize, 1. - pad_abs_->GetTopMargin() - ticksize  );
  leg->SetFillStyle( 0 ) ;
  leg->SetBorderSize( 0 );
  
  //== Prepare the reference.
  TH1D* hReference = (TH1D*)hData->Clone("Reference");
  TH1D* hReference_ratio = (TH1D*)hData->Clone("Reference_ratio");  
  hReference_ratio->Divide( hReference );
  
  for( int bin = 0; bin <= hData->GetNbinsX()+1; bin++){
    hReference->SetBinError( bin, 0. );
    hReference_ratio->SetBinError( bin, 0. );
  }
  
  //== Draw the data.
  pad_abs_->cd();
  Prepare_1Dplot_ratio( hData );  
  hData->GetXaxis()->SetRangeUser(150.,1700.);
  hData->GetYaxis()->SetTitle("d#sigma/dE [mb/GeV]");
  hData->DrawClone("phist");
  
  pad_ratio_->cd();
  hData->Divide( hReference );
  
  hData->GetXaxis()->SetTitle("E [GeV]");  
  hData->GetXaxis()->SetLabelSize( hData->GetXaxis()->GetLabelSize()*0.85 );

  
  hData->GetYaxis()->SetTitle("Ratio");
  hData->GetYaxis()->CenterTitle();
  hData->GetYaxis()->SetNdivisions( 504 );
  hData->GetYaxis()->SetRangeUser( 0. , 2.3 );
  
  hData->DrawClone("phist");    
  
  //== Loop over the histograms.
  for( int _syst = 0; _syst < vSyst.size(); _syst++){
  std::cout << "systematic\t" << _syst << std::endl;
  
    TH1D* hSyst = vSyst[_syst];
    
    for(int bin = 0; bin <= hData->GetNbinsX()+1; bin++){
      double syst_content = hSyst->GetBinContent( bin );
      double hist_content = hReference->GetBinContent( bin );
      double hist_error = hReference->GetBinError( bin );
      double new_error = max( hist_error, fabs(hist_content-syst_content) );
      hReference->SetBinError( bin, new_error );
      hReference_ratio->SetBinError( bin, new_error/hist_content );
    }
    
    if( TString::Format(hSyst->GetTitle()).Contains("Pythia8") ){ leg->AddEntry( hSyst, "Pythia8 (4C)", "l" ); }
    if( TString::Format(hSyst->GetTitle()).Contains("Pythia6") ){ leg->AddEntry( hSyst, "Pythia6 (Z2*)", "l" ); }
    if( TString::Format(hSyst->GetTitle()).Contains("LumiUp") ){ leg->AddEntry( hSyst, "Luminosity (+3.6%)", "l" ); }            
    if( TString::Format(hSyst->GetTitle()).Contains("LumiDown") ){ leg->AddEntry( hSyst, "Luminosity (-3.6%)", "l" ); }            
    if( TString::Format(hSyst->GetTitle()).Contains("EPOS") ){ leg->AddEntry( hSyst, "EPOS", "l" ); }        
    if( TString::Format(hSyst->GetTitle()).Contains("JetID") ){ leg->AddEntry( hSyst, "Overcalibrated data", "l" ); }
    if( TString::Format(hSyst->GetTitle()).Contains("PositionDown") ){ leg->AddEntry( hSyst, "Position Down", "l" ); }
    if( TString::Format(hSyst->GetTitle()).Contains("PositionUp") ){ leg->AddEntry( hSyst, "Position Up", "l" ); }        

        
  }//== Loop over histograms.

  int fillcolor_ = 0;
  
  std::cout << "Name is\t" << name << endl;
  
  if( name == "Model" ){	fillcolor_ = TColor::GetColor("#C0C0C0"); }
  if( name == "Position" ){	fillcolor_ = TColor::GetColor("#33FF99"); }
  if( name == "JetID" ){	fillcolor_ = TColor::GetColor("#000099"); }
  if( name == "Luminosity" ){	fillcolor_ = TColor::GetColor("#FFCCCC"); }
  if( name == "JES" ){		fillcolor_ = TColor::GetColor("#FFFF00"); }

  pad_abs_->cd();
  hReference->SetLineWidth(2);
  hReference->DrawClone("histsame");
  hReference->SetFillColor( fillcolor_  );
  hReference->DrawClone("E2same"); 
  
  pad_ratio_->cd();
  hReference_ratio->DrawClone("histsame");
  hReference_ratio->SetFillColor(  fillcolor_ );  
  hReference_ratio->DrawClone("E2same");

  for( int _syst = 0; _syst < vSyst.size(); _syst++){
  
    TH1D* hSyst = vSyst[_syst];    
    hSyst->SetLineColor( getColor( _syst + 2 ) );
    hSyst->SetLineStyle( _syst + 2  );    
    hSyst->SetLineWidth( 2 );
    
    pad_abs_->cd();    
    hSyst->DrawClone("histsame");
    
    pad_ratio_->cd();
    hSyst->Divide( hReference );
    hSyst->DrawClone("histsame");    
  }
  
  if( leg->GetNRows() > 0 ){
    leg->SetHeader( name );
    pad_abs_->cd();
    leg->Draw();
  }


  pad_abs_->SetLogy();  
  can_->SaveAs( TString::Format("Canvas_uncertainty_" + name + ".pdf") );
}
