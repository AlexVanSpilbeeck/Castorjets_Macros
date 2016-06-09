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

   TCanvas *CHI2_Test_all = new TCanvas("CHI2_Test_all", "CHI2_Test_all",0,0,800,756);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   CHI2_Test_all->SetHighLightColor(2);
   CHI2_Test_all->Range(0,0,1,1);
   CHI2_Test_all->SetFillColor(0);
   CHI2_Test_all->SetBorderMode(0);
   CHI2_Test_all->SetBorderSize(2);
   CHI2_Test_all->SetLogy();
   CHI2_Test_all->SetTickx(1);
   CHI2_Test_all->SetTicky(1);
   CHI2_Test_all->SetLeftMargin(0.12);
   CHI2_Test_all->SetRightMargin(0.04);
   CHI2_Test_all->SetTopMargin(0.06315);
   CHI2_Test_all->SetFrameFillStyle(0);
   CHI2_Test_all->SetFrameBorderMode(0);
   
   PrepareCanvas( CHI2_Test_all, "CHi2_test");

   TLegend* leg = new TLegend(
	0.65, 
	0.75, 
	1. - CHI2_Test_all->GetRightMargin() - 0.02, 
	1. - CHI2_Test_all->GetTopMargin() - 0.02);
   leg->SetFillStyle( 0 );
   leg->SetBorderSize( 0 );


   //=======================

   TGraph *graph = new TGraph(6);
   graph->SetName("Graph");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);
   graph->SetMarkerStyle(3);
   graph->SetPoint(0,1,399.214);
   graph->SetPoint(1,4,13.94366);
   graph->SetPoint(2,7,7.185133);
   graph->SetPoint(3,10,7.527302);
   graph->SetPoint(4,20,7.369105);
   graph->SetPoint(5,30,6.888159);
   
   TH1F *Graph1 = new TH1F("Graph1","Graph1",100,0,32.9);
   Graph1->SetMinimum(2.);
   Graph1->SetMaximum(500.);
   Graph1->SetDirectory(0);
   Graph1->SetStats(0);
   Graph1->SetLineStyle(0);
   Graph1->SetMarkerStyle(20);
   Graph1->GetXaxis()->SetTitle("N_{it.}");
   Graph1->GetXaxis()->SetLabelFont(42);
   Graph1->GetXaxis()->SetLabelOffset(0.007);
   Graph1->GetXaxis()->SetLabelSize(0.05);
   Graph1->GetXaxis()->SetTitleSize(0.06);
   Graph1->GetXaxis()->SetTitleOffset(0.9);
   Graph1->GetXaxis()->SetTitleFont(42);
   Graph1->GetYaxis()->SetTitle("#chi^{2}/NDF");
   Graph1->GetYaxis()->SetLabelFont(42);
   Graph1->GetYaxis()->SetLabelOffset(0.007);
   Graph1->GetYaxis()->SetLabelSize(0.05);
   Graph1->GetYaxis()->SetTitleSize(0.06);
   Graph1->GetYaxis()->SetTitleOffset(1.25);
   Graph1->GetYaxis()->SetTitleFont(42);
   Graph1->GetYaxis()->SetRangeUser(2.,500.);
   Graph1->GetZaxis()->SetLabelFont(42);
   Graph1->GetZaxis()->SetLabelOffset(0.007);
   Graph1->GetZaxis()->SetLabelSize(0.05);
   Graph1->GetZaxis()->SetTitleSize(0.06);
   Graph1->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph1);
   Prepare_1Dplot( graph );   
   graph->GetHistogram()->GetYaxis()->SetRangeUser(2., 500.);

   graph->SetMarkerStyle(21);
   graph->SetMarkerColor( kBlack );
   graph->Draw("ap");

   leg->AddEntry( graph, "Had. cal.", "p");

   //=======================

 
   TGraph *graph_allJets = new TGraph(6);
   graph_allJets->SetName("Graph");
   graph_allJets->SetTitle("Graph");
   graph_allJets->SetFillColor(1);
   graph_allJets->SetMarkerStyle(3);
   graph_allJets->SetPoint(0,1,454.6779);
   graph_allJets->SetPoint(1,4,21.47268);
   graph_allJets->SetPoint(2,7,8.187766);
   graph_allJets->SetPoint(3,10,6.556972);
   graph_allJets->SetPoint(4,20,5.828577);
   graph_allJets->SetPoint(5,30,5.861616);
   
   TH1F *Graph_allJets = new TH1F("Graph_allJets","Graph_allJets",100,0,32.9);
   Graph_allJets->SetMinimum(1.);
   Graph_allJets->SetMaximum(500.);
   Graph_allJets->SetDirectory(0);
   Graph_allJets->SetStats(0);
   Graph_allJets->SetLineStyle(0);
   Graph_allJets->SetMarkerStyle(20);
   Graph_allJets->GetXaxis()->SetTitle("N_{it.}");
   Graph_allJets->GetXaxis()->SetLabelFont(42);
   Graph_allJets->GetXaxis()->SetLabelOffset(0.007);
   Graph_allJets->GetXaxis()->SetLabelSize(0.05);
   Graph_allJets->GetXaxis()->SetTitleSize(0.06);
   Graph_allJets->GetXaxis()->SetTitleOffset(0.9);
   Graph_allJets->GetXaxis()->SetTitleFont(42);
   Graph_allJets->GetYaxis()->SetTitle("#chi^{2}/NDF");
   Graph_allJets->GetYaxis()->SetLabelFont(42);
   Graph_allJets->GetYaxis()->SetLabelOffset(0.007);
   Graph_allJets->GetYaxis()->SetLabelSize(0.05);
   Graph_allJets->GetYaxis()->SetTitleSize(0.06);
   Graph_allJets->GetYaxis()->SetTitleOffset(1.25);
   Graph_allJets->GetYaxis()->SetTitleFont(42);
   Graph_allJets->GetZaxis()->SetLabelFont(42);
   Graph_allJets->GetZaxis()->SetLabelOffset(0.007);
   Graph_allJets->GetZaxis()->SetLabelSize(0.05);
   Graph_allJets->GetZaxis()->SetTitleSize(0.06);
   Graph_allJets->GetZaxis()->SetTitleFont(42);
   graph_allJets->SetHistogram(Graph_allJets);
   
   graph_allJets->SetMarkerStyle(20 );
   graph_allJets->SetMarkerColor( kRed );
   graph_allJets->Draw("psame");

   leg->AddEntry( graph_allJets, "All cal.", "p");

   //===========


   CHI2_Test_all->Modified();
   CHI2_Test_all->SetLogy();
   CHI2_Test_all->cd();
   leg->Draw();
   CHI2_Test_all->SetSelected(CHI2_Test_all);
   CHI2_Test_all->SaveAs("can_compare_Chi2_calibration.pdf");
   CHI2_Test_all->SaveAs("can_compare_Chi2_calibration.C");
   

  return 0;
}
