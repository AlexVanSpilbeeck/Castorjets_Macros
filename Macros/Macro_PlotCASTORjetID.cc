#include "iostream"
#include "fstream"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TColor.h"
#include "TApplication.h"
#include "stdio.h"
#include "cstring"
#include "TMath.h"

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

#include "../Functions/Function_average_histogram.h"
#include "../Functions/Function_CorrectSectorNumber.h"
#include "../Functions/Function_FinishCanvas.C"
#include "../Functions/Function_Prepare1Dplot.h"
#include "../Functions/Function_Prepare2Dplot.h"
#include "../Functions/Function_PrepareCanvas.h"
#include "../Functions/Function_Rebin.h"
#include "../Functions/Function_SetRangeToIncludeAll.h"
#include "../GetSubHistogram.h"
#include "../Correct_Sector.h"
//#include "Function_make_Tex.h"
#include "../Functions/Function_NonZeroMinimum.h"

#include "../color.h"
using namespace std; 


int main(){
  gStyle->SetOptTitle(0);

  vector<TString> jetIDs;
  std::map<TString, TString> yaxis;
  jetIDs.push_back("ehad");
  yaxis["ehad"] = "1/N_{tot} dN/dE_{had} [1/GeV]";
  jetIDs.push_back("eem");
  yaxis["eem"] = "1/N_{tot} dN/dE_{EM} [1/GeV]";
  jetIDs.push_back("nTowers");
  yaxis["nTowers"] = "1/N_{tot} dN/dN_{tower}";
  jetIDs.push_back("sigmaz");
  yaxis["sigmaz"] = "1/N_{tot} dN/d#sigma_{z}";
  jetIDs.push_back("width");
  yaxis["width"] = "1/N_{tot} dN/dwidth";
  jetIDs.push_back("depth");
  yaxis["depth"] = "1/N_{tot} dN/d<z> [1/mm]";
  jetIDs.push_back("fhot");
  yaxis["fhot"] = "1/N_{tot} dN/df_{hot}";
  jetIDs.push_back("fem");
  yaxis["fem"] = "1/N_{tot} dN/df_{EM}";
  jetIDs.push_back("phi");
  yaxis["phi"] = "1/N_{tot} dN/d#phi";


  vector<TString> files;
  std::map<TString, TString> legend_entries;

  TString file_ = "/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_ShowerLibrary_NomGeo_calibrated_Emin_150.000000_deltaPhiMax_0.200000_etaband_0.000000_had_matchE.root";
  files.push_back( file_ );
  legend_entries[ file_ ] = "Pythia6 (Z2*) [SL]";

  file_ =  "/user/avanspil/Castor_Analysis/Stripped_trees_histo_files/ak5ak5_Pythia6Z2star_FullSimulation_NomGeo_calibrated_Emin_150.000000_deltaPhiMax_0.200000_etaband_0.000000_had_matchE.root";
  files.push_back( file_);
  legend_entries[ file_ ] = "Pythia6 (Z2*) [FS]";

  //== Make a new directory for the plots/
   int 	new_dir = mkdir( (TString::Format("CASTORjetIDs")).Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  );

  for(vector<TString>::iterator var = jetIDs.begin(); var != jetIDs.end(); ++var){

    TString variable = *var;
    cout << "Variable\t" << variable << endl;
    TCanvas *can;
    PrepareCanvas( can, TString::Format( "can_" + variable ) );

    TString drawoptions = "hist";
    int color = 1;

    double xmin, ymin;  

    if( (*var) != "ehad" && (*var) != "eem" &&  (*var) != "depth" ){ 
       xmin = 0.35; 
       ymin = can->GetBottomMargin();
    }
    else if( (*var) == "depth" ){ 
      xmin = can->GetLeftMargin(); 
      ymin = 1. - can->GetTopMargin() - 0.3;
    }
    else{
      xmin = 1. - can->GetRightMargin() - 0.3;
      ymin = 1. - can->GetTopMargin() - 0.3;
    }
    TLegend *leg = new TLegend( xmin, ymin, xmin + 0.3, ymin + 0.3);

    leg->SetFillStyle( 0 );
    leg->SetBorderSize( 0 );

      for(vector<TString>::iterator file_ = files.begin(); file_ < files.end(); ++file_){

	TFile* _file = TFile::Open( *file_, "Read");
	TH1D* hist = (TH1D*)_file->Get( variable );

	hist->Scale( 1./hist->Integral() );
	hist->SetLineColor( getColor( color ) );
	hist->SetLineStyle ( color++ );
	hist->SetLineWidth ( 2 );
	Prepare_1Dplot( hist );
        hist->GetYaxis()->SetRangeUser( 0.9*GetMinimumValue( hist), 1.1 * hist->GetMaximum() );
        hist->GetYaxis()->SetTitle( yaxis[*var] );
	if( *var == "depth" ) hist->GetXaxis()->SetNdivisions(504);
	hist->Draw( drawoptions );
	drawoptions = "histsame";

	leg->AddEntry( hist, legend_entries[ *file_ ] , "l" );
       } //== Loop over files

       leg->Draw();  
  
       if( *var != "phi" ) { can->SetLogy(); }

       can->SaveAs( TString::Format( "CASTORjetIDs/CastorJetID_" + variable + ".pdf" ) );
       can->SaveAs( TString::Format( "CASTORjetIDs/CastorJetID_" + variable + ".C" ) );
     } //== Loop over variables.

  return 0;
}
