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
#include <THnSparse.h>
#include <TThread.h>
#include <TStopwatch.h>
#include <TLegend.h>
#include <TPad.h>
#include <TColor.h>
#include <TString.h>
#include <TText.h>
#include <TObjString.h>

#include <TLine.h>
#include <TPaletteAxis.h>

//STANDARD C++ INCLUDES
#include <sstream>
#include <iostream>
#include <iomanip>
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
#include <vector>
#include <ctime>



#include "../tdrstyle.C"
#include "../color.h"

#include "../Functions/Function_DivideByBinWidth.h"
#include "../Functions/Function_PrepareCanvas.h"
#include "../Functions/Function_Prepare1Dplot.h"
#include "../Functions/Function_Prepare2Dplot.h"

using namespace std;

int main(){


  gStyle->SetOptTitle( 0 );
  setTDRStyle();
  gStyle->SetPalette(1);
  
  //== Type of event selections.
  TString cuts_;
  std::vector<TString> vector_cuts_;
  std::map<TString, TString> legend_cuts;

  cuts_ = "";
  vector_cuts_.push_back( cuts_ );
  legend_cuts[cuts_] = "Default selection";

  cuts_ = "_BSC_OR";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "BSC OR, == 1 gVtx";
  
  cuts_ = "_noVertex_HF";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "HF, 0 gVtx";

  cuts_ = "_noVertex_BSC";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "BSC OR, 0 gVtx";
  
  cuts_ = "_atLeastOneVertex_HF";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "HF, #geq 1 gVtx";
  
  cuts_ = "_atLeastOneVertex_BSC";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "BSC OR, #geq 1 gVtx";  
  
  cuts_ = "_one_none_vertex_BSC";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "BSC OR, #leq 1 gVtx";  
    
  cuts_ = "_one_none_vertex_BSCor_HFor";
  vector_cuts_.push_back( cuts_ );  
  legend_cuts[cuts_] = "BSC OR, HF OR, #leq 1 gVtx";       
  
  

  for( int c_ = 0; c_ < vector_cuts_.size(); c_++){
    cuts_ = vector_cuts_[ c_ ];

    TFile* _file = TFile::Open("/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5_data" + cuts_ + "_unfold_Emin_150.000000_all.root", "Read");  
    TH1D* hIstlumibx = (TH1D*)_file->Get("hIstlumibx");
    if( !hIstlumibx ){ continue; }
  
    TCanvas *can_;
    PrepareCanvas(can_, "Comparing_lumi_BX");
    TPad* pad_abs_, *pad_ratio_;
    SplitCanvas( can_, pad_abs_, pad_ratio_ );
    TString drawoptions = "ephist";
    
    double ticksize = gStyle->GetTickLength();
    
    TLegend *leg = new TLegend( 0.65, 0.35, 1. - pad_abs_->GetRightMargin()-ticksize, 1. - pad_abs_->GetTopMargin()-ticksize);
    leg->SetHeader( legend_cuts[cuts_] );
    leg->SetFillStyle( 0 );
    leg->SetBorderSize( 0 );
    
    TH1D* hFirst;
  
    for(int current_lumi_bin = 1; current_lumi_bin <= hIstlumibx->GetNbinsX(); current_lumi_bin++){
  
      TString histname_ = TString::Format("hCastorJet_energy_istlumi_bin_%i", current_lumi_bin);
      TH1D *hEflow = (TH1D*)_file->Get( histname_ );
      hEflow->Sumw2();

      double norm = hIstlumibx->GetBinContent( current_lumi_bin );
      hEflow->Scale( 1./norm );
      SetDnDx( hEflow ); 
      hEflow->GetXaxis()->SetRangeUser(150., 1700.);           
      
      if( current_lumi_bin == 1 ){      
        hFirst = (TH1D*)hEflow->Clone("hFirst");
      }    
    
      hEflow->SetLineColor( getColor( current_lumi_bin ) );
      hEflow->SetLineStyle( current_lumi_bin%5 + 2 );
      hEflow->SetMarkerColor( getColor( current_lumi_bin ) );
      hEflow->SetMarkerStyle( 19 + current_lumi_bin );      
      hEflow->GetYaxis()->SetTitle("dN/dE [1/GeV]");
      hEflow->GetXaxis()->SetTitle("E_{cal} [GeV]");
      hEflow->GetYaxis()->SetRangeUser( 3.e-6, 	hEflow->GetMaximum()*1.5 );      	
      Prepare_1Dplot_ratio( hEflow );
      
      pad_abs_->cd();
      hEflow->DrawClone( drawoptions );
      
      pad_ratio_->cd();
      hEflow->Divide( hFirst );
      hEflow->GetYaxis()->SetTitle("Ratio");
      hEflow->GetYaxis()->SetRangeUser( 0.79, 1.24 );
      hEflow->GetYaxis()->SetNdivisions( 504 );
      hEflow->GetYaxis()->CenterTitle();
      hEflow->DrawClone( drawoptions );
      
      drawoptions = "ephistsame";
      
      int lowerlumi = 20;
      int lumiwidth = 2;
      leg->AddEntry( hEflow, TString::Format("L = 0.00%i - 0.00%i", 
      	lowerlumi + lumiwidth * (current_lumi_bin - 1 ),
      	lowerlumi + lumiwidth * current_lumi_bin  ), "p");
    }

    pad_abs_->cd();
    leg->Draw();
    pad_abs_->SetLogy();    
    can_->SaveAs( "eflow" + cuts_ + ".pdf");

  }

  return 0;
}
