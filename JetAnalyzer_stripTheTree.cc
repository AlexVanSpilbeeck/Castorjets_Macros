////////////////////////////////////////
//////// New CMSSW_4_2_X version ///////
////////////////////////////////////////

/* Updated by Alex Van Spilbeeck.

 This code takes the original trees and removes the unnecessary events. 
*/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "JetAnalyzer_stripTheTree.h"
#include "HistoRetriever.h"

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
#include <TBranchElement.h>

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

// own classes includes
#include "../src/MyCastorRecHit.h"
#include "../src/MyCastorDigi.h"
#include "../src/MyCastorTower.h"
#include "../src/MyCastorJet.h"
#include "../src/MyEvtId.h"
#include "../src/MyGenKin.h"
#include "../src/MyDiJet.h"
#include "../src/MyVertex.h"
#include "../src/MyHLTrig.h"
#include "../src/MyL1Trig.h"
#include "../src/MyJet.h"
#include "../src/MyBeamSpot.h"
#include "../src/MyGenPart.h"
#include "../src/MyCaloTower.h"
#include "../src/MyGenJet.h"
#include "../src/MyTrackJet.h"

//Fastjet
//#include "fastjet/PseudoJet.hh"
//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/PseudoJet.hh"
//using namespace fastjet;

#define jetPtThreshold 35.
#define jetEThreshold_gen 0. 
#define jetEThreshold_det 0. 
#define EbinWidth 5.
#define EbinWidth_rel 1.4
#define phi_diff_max 0.1
#define fileSteps 49

#define jet_distance 10
#define jet_distance_string "JER_allPairs__JetSorted_ak5"
//#define gen_radius_ "ak7"
//#define det_radius_ "ak7"
#define GenJetContained -0.7
//#define totalEvents_ 10000
#define PI 3.14159265359


#include "../../../RooUnfold-1.1.1/src/RooUnfold.h"
//#include "../../../RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#endif



TFile *JetAnalyzer_stripTheTree::currentStaticTFile_ = new TFile();

JetAnalyzer_stripTheTree::JetAnalyzer_stripTheTree(TString inputdir, TObjArray* filelist, bool isData, const char* outputname, TString gen_radius, TString det_radius, int totalEvents, TString date, int startFile) {
    
	std::cout << "constructing JetAnalyzer_stripTheTree class..." << std::endl;
	
    inputdir_ = inputdir;
    filelist_ = filelist;
    isData_ = isData;
    outputname_ = outputname;
    gen_radius_ = gen_radius;
    det_radius_ = det_radius;
    totalEvents_ = totalEvents;
    date_ = date;	
    startFile_ = startFile;
	std::cout << "initializing basic variables..." << std::endl;
    
    // initialize basic variables

    
    LoopOutputFile_ = "";

	currentfile_ = "";
	currentTFile_ = new TFile();
	
	std::cout << "all initialisations done, class constructed!" << std::endl;
    
}

JetAnalyzer_stripTheTree::~JetAnalyzer_stripTheTree() { }

void JetAnalyzer_stripTheTree::Loop() {

#ifdef __CINT__
  gSystem->Load("../../../RooUnfold-1.1.1/libRooUnfold.so");
#endif
	
	
	std::cout << " JetAnalyzer_stripTheTree Loop function is started " << std::endl;
	
	TString tstring = outputname_;
	std::cout << " TString outputname_ = " << tstring << std::endl;
	TString string_tag = jet_distance_string;
	TString string_det_radius = det_radius_;
        TString string_gen_radius;	if (!isData_){ string_gen_radius = gen_radius_; }
 
	
    // reweight the MC in this case
    bool reweightMC = false;
    if (tstring.Contains("Reweighted")) {
        std::cout << "We will reweight the MC now !!!" << std::endl;
        reweightMC = true;
    }
    
	
	using namespace std;
	int it = 0;
	int totalevents = 0;
	
	/////////////////////////////////////
	// Define all histograms
	/////////////////////////////////////


	int Ebins_fix = 100.;
	double Emin_fix = -100.;
	double Emax_fix = 3000.;
	
        // AVS - number of bins.
	/* We get binwidth from peak JER. */
//	double Emin = jetEThreshold_det + 0., Emax = 3000.;
//	int Ebins = static_cast<int> ( (Emax - Emin)/EbinWidth );
	double Emin = 0., Emax = 3000.;
	int Ebins = 200.;

	// AVS - count number of events with MC jets in CASTOR.
	int total_events_MC = 0, castor_events_MC = 0, castor_150GeV_events_MC = 0, castor_300GeV_events_MC = 0;


	/* If we want variable Ebins. */
/*
	double current_lowE = 200.;
	int count_bincenters = 3;

	vector<double> binEdges;
	binEdges.push_back( -50.);
	binEdges.push_back( 0. );
	binEdges.push_back( 50. );
	binEdges.push_back( jetEThreshold_det );
	binEdges.push_back( current_lowE );
	while( current_lowE/EbinWidth_rel < Emax ){

	  binEdges.push_back( current_lowE * EbinWidth_rel );

 	  cout << "\tLower edge\t" << count_bincenters << "\t" << current_lowE << "\tBinwidth\t" << current_lowE *0.2 << endl;

	  current_lowE += current_lowE*0.2;
	  count_bincenters++;
	}

//	Ebins = binEdges.size();
        float *Ebins_var = new float[Ebins];
//float Ebins_var[Ebins];
	for(int i = 0; i < binEdges.size(); i++){
	  Ebins_var[i] = binEdges[i];

	  cout << "\t\tBin " << i << "\tout of\t" << sizeof(Ebins_var) << " or " << binEdges.size() << " or " << Ebins << "\t" << Ebins_var[i] << "\t" << binEdges[i] << endl;
	}

	// We need an extended bin range for our matrix + hits & misses
	// Let us take bins 0 - 50, 50 - 100
	float *Ebins_ext = new float[Ebins+2];
	Ebins_ext[0] = 0.;
	Ebins_ext[1] = 50.;
	for(int i = 0; i < binEdges.size(); i++){
	  Ebins_ext[i+2] = binEdges[i];
	}
*/
	//-------------------------------------------//
	//-- We want bin width 20% of high bin edge. //
	//-------------------------------------------// 
	double current_lowE = 150.;
	vector<double> binEdges;
	binEdges.push_back( 150. );
	double binwidth = 75.;
	Ebins = 1;

	while( current_lowE < 2100. ){
          // Binwidth = 25% bin_low_edge;
	  binwidth = .25 * current_lowE;	
	  current_lowE += binwidth;
	  binEdges.push_back( current_lowE );
	}

	double* Ebins_var = &binEdges[0];

	for(int i = 0; i < binEdges.size(); i++){
	  Ebins_var[i] = binEdges[i];
	  cout << "bin\t" << i << "\t" << binEdges[i] << "\t" << Ebins << "\tOr\t" << binEdges.size() << endl;
	}

	for(int i = 0; i < binEdges.size(); i++){
	  cout << "Ebins_var\t" << i << "\t" << Ebins_var[i] << endl;

	}
	cout << "Done!" << endl;

	Ebins = binEdges.size() - 1;

	cout << "Number of bins is\t" << Ebins << endl;

	//-------------------------------------------//
	//-- We wanted bin width 20% of high bin edge. //
	//-------------------------------------------// 




        /* End variable ebins. */

cout << "Variable bins done" << endl;
	// detector level histograms
	
	char name [100];
	char title [100];
		
	TH1D *hCASTORTowerMulti = new TH1D("hCASTORTowerMulti","CASTOR Tower Multiplicity (N above threshold) distribution",17,0,17);
	hCASTORTowerMulti->Sumw2();
	
	// default energy flow histos - using 5 modules
	TH1D *hCASTOReflow = new TH1D("hCASTOReflow","Total CASTOR energy flow in first 5 modules",252,-30,3750);
	hCASTOReflow->Sumw2();
	
	TH2D *h2CASTOReflow_grid = new TH2D("h2CASTOReflow_grid","CASTOR energy weighted module vs sector distribution",16,1,17,14,1,15);
	
	TH1D *hCASTOReflow_channel[224];
	for (int i=0;i<224;i++) {
		sprintf(name,"hCASTOReflow_channel_%d",i+1);
		sprintf(title,"CASTOR Energy distribution for channel %d",i+1);
		hCASTOReflow_channel[i] = new TH1D(name,title,100,-10,1800);
		hCASTOReflow_channel[i]->Sumw2();
	}
	
	// Castor jet histograms
//	TH1D *hCastorJet_energy = new TH1D("hCastorJet_energy","CastorJet energy distribution",Ebins, Emin, Emax);
        TH1D *hCastorJet_energy = new TH1D("hCastorJet_energy","CastorJet energy distribution",	Ebins, Emin, Emax);
        TH1D *hGenJet_energy = new TH1D("hGenJet_energy","GenJet energy distribution",		Ebins, Emin, Emax);
        TH1D *hGenJet_energy_noCuts = new TH1D("hGenJet_energy_noCuts","GenJet energy distribution",		Ebins, Ebins_var);
	TH1D *hCastorJet_pt = new TH1D("hCastorJet_pt","CastorJet pt distribution",30,0,30);
	TH1D *hCastorJet_em = new TH1D("hCastorJet_em","CastorJet EM energy distribution",150,0,1500);
	TH1D *hCastorJet_had = new TH1D("hCastorJet_had","CastorJet HAD energy distribution",150,0,1500);
	TH1D *hCastorJet_fem = new TH1D("hCastorJet_fem","CastorJet EM/(EM+HAD) distribution",60,-0.1,1.1);
	TH1D *hCastorJet_fhot = new TH1D("hCastorJet_fhot","CastorJet Fhot distribution",60,-0.1,1.1);
	TH1D *hCastorJet_width = new TH1D("hCastorJet_width","CastorJet width distribution",100,0,1);
	TH1D *hCastorJet_depth = new TH1D("hCastorJet_depth","CastorJet depth distribution",100,-16000,-14000);
	TH1D *hCastorJet_sigmaz = new TH1D("hCastorJet_sigmaz","CastorJet sigmaz distribution",100,0,500);
	TH1D *hCastorJet_ntower = new TH1D("hCastorJet_ntower","CastorJet ntower distribution",16,1,17);
	TH1D *hCastorJet_eta = new TH1D("hCastorJet_eta","CastorJet eta distribution",14,-6.6,-5.2);
	TH1D *hCastorJet_phi = new TH1D("hCastorJet_phi","CastorJet phi distribution",16,-M_PI,+M_PI);
	TH1D *hCastorJet_multi = new TH1D("hCastorJet_multi","CastorJet multiplicity distribution",17,0,17);

	// Central-forward configuration response matrix.
        TH1D *hCastorJet_cf_energy = new TH1D("hCastorJet_cf_energy", "CastorJet energy in with leading central jet",					Ebins+10,Emin-10*EbinWidth,Emax);
        TH1D *hCastorJet_cf_energy_gen = new TH1D("hCastorJet_cf_energy_gen", "CastorJet energy in with leading central jet",				Ebins+10,Emin-10*EbinWidth,Emax);
	TH2D *hCastorJet_cf_energy_response = new TH2D("hCastorJet_cf_energy_response", "CastorJet energy response matrix for leading central jet",	Ebins+10,Emin-10*EbinWidth,Emax,Ebins+10,Emin-10*EbinWidth,Emax);
        TH1D *hCastorJet_cf_energy_fakes = new TH1D("hCastorJet_cf_energy_fakes", "CastorJet energy in with leading central jet - fakes",Ebins+10,Emin-10*EbinWidth,Emax);
        TH1D *hCastorJet_cf_energy_misses = new TH1D("hCastorJet_cf_energy_misses", "CastorJet energy in with leading central jet - misses",Ebins+10,Emin-10*EbinWidth,Emax);

	// All Casto jets response matrix.
        TH2D *hCastorJet_energy_response = new TH2D("hCastorJet_energy_response","CastorJet energy distribution - Response", 	Ebins, Emin, Emax, Ebins, Emin, Emax);	// ,Ebins+10,Emin-10*EbinWidth,Emax,Ebins+10,Emin-10*EbinWidth,Emax);
        TH2D *hCastorJet_energy_fakes 	 = new TH2D("hCastorJet_energy_fakes",	 "CastorJet energy distribution - Fakes", 	Ebins, Emin, Emax, Ebins, Emin, Emax);	// ,Ebins+10,Emin-10*EbinWidth,Emax);
        TH2D *hCastorJet_energy_misses	 = new TH2D("hCastorJet_energy_misses",	 "CastorJet energy distribution - Misses", 	Ebins, Emin, Emax, Ebins, Emin, Emax);	// ,Ebins+10,Emin-10*EbinWidth,Emax);
	
	TH2D *hCastorJet_energy_complete_response = new TH2D("hCastorJet_energy_complete_response","CastorJet energy distribution - Complete Response", Ebins_fix, Emin_fix, Emax_fix, Ebins_fix, Emin_fix, Emax_fix);

	// Response with variable binning.
/*
        TH2D *hCastorJet_energy_response = new TH2D("hCastorJet_energy_response","CastorJet energy distribution  - Response",Ebins-1, Ebins_var,Ebins-1, Ebins_var);
        TH1D *hCastorJet_energy_fakes = new TH1D("hCastorJet_energy_fakes","CastorJet energy distribution - Fakes", Ebins-1, Ebins_var);
        TH1D *hCastorJet_energy_misses = new TH1D("hCastorJet_energy_misses","CastorJet energy distribution - Misses",Ebins-1, Ebins_var);
*/
//	TH2D *hCastorJet_energy_ratio = new TH2D("hCastorJet_energy_ratio", "Ratio of generator to detector", 1000, 0., 10.,Ebins-1, Ebins_var);
        TH2D *hCastorJet_energy_ratio = new TH2D("hCastorJet_energy_ratio", "Ratio of generator to detector", 1000, 0., 10.,Ebins, Emin, Emax);

	TH2D *hCastorJet_cf_Matrix = new TH2D("hCastorJet_cf_responseMatrix", "Central-forward jet Response Matrix",Ebins+10,Emin-10*EbinWidth,Emax,Ebins+10,Emin-10*EbinWidth,Emax);
//        TH2D *hCastorJet_Matrix = new TH2D("hCastorJet_responseMatrix", "Castor jet Response Matrix",Ebins+10,Emin-10*EbinWidth,Emax,Ebins+10,Emin-10*EbinWidth,Emax);
        TH2D *hCastorJet_Matrix = new TH2D("hCastorJet_responseMatrix", "Castor jet Response Matrix",Ebins-1, Ebins_var,Ebins-1, Ebins_var);

	  TH2D *hCastorJet_Matrix_had_pi = new TH2D("hCastorJet_responseMatrix_had_pi", "Castor jet Response Matrix (had-had)"	,Ebins, Emin, Emax, Ebins, Emin, Emax);
          TH2D *hCastorJet_Matrix_had_e = new TH2D("hCastorJet_responseMatrix_had_e", 	"Castor jet Response Matrix (had-em)"	,Ebins, Emin, Emax, Ebins, Emin, Emax);
          TH2D *hCastorJet_Matrix_em_pi = new TH2D("hCastorJet_responseMatrix_em_pi", 	"Castor jet Response Matrix (em-had)"	,Ebins, Emin, Emax, Ebins, Emin, Emax);
          TH2D *hCastorJet_Matrix_em_e 	= new TH2D("hCastorJet_responseMatrix_em_e", 	"Castor jet Response Matrix (em-em)"	,Ebins, Emin, Emax, Ebins, Emin, Emax);
//          TH2D *hCastorJet_Matrix_had_pi = new TH2D("hCastorJet_responseMatrix_had_pi", "Castor jet Response Matrix (had-had)",Ebins-1, Ebins_var,Ebins-1, Ebins_var);
//          TH2D *hCastorJet_Matrix_had_pi = new TH2D("hCastorJet_responseMatrix_had_pi", "Castor jet Response Matrix (had-had)",Ebins-1, Ebins_var,Ebins-1, Ebins_var);


	// Number of trackjets.
	TH2D *hTrackjets_2D_number = new TH2D("hTrackjets_2D_number", "Trackjets vs. Charged Gen Jet", 17,0.,17., 17,0.,17.);
        TH2D *hTrackjets_2D_pt = new TH2D("hTrackjets_2D_pt", "Trackjets vs. Charged Gen Jet", 30,0.,30., 30,0.,30.);

        TH1D *hCastorJet_energy_gen = new TH1D("hCastorJet_energy_gen","CastorJet energy distribution",Ebins-1, Ebins_var);
        TH1D *hCastorJet_pt_gen = new TH1D("hCastorJet_pt_gen","CastorJet pt distribution",500,0,10.);

	// Efficiency control.
	TH1D *hJER = new TH1D("hJER", "Jet Energy Resolution", 200, -5, 5);
        TH2D *hJER_per_energy = new TH2D("hJER_per_energy", "#DeltaE/E for fixed energies;E_{gen};JER", 			Ebins, Emin, Emax, 200, -5, 5);
	     TH2D *hJER_per_energy_had_pi = new TH2D("hJER_per_energy_had_pi", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",	Ebins, Emin, Emax, 200, -5, 5);
             TH2D *hJER_per_energy_had_e  = new TH2D("hJER_per_energy_had_e",  "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",     	Ebins, Emin, Emax, 200, -5, 5);
             TH2D *hJER_per_energy_em_pi  = new TH2D("hJER_per_energy_em_pi",  "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",     	Ebins, Emin, Emax, 200, -5, 5);
             TH2D *hJER_per_energy_em_e	  = new TH2D("hJER_per_energy_em_e",   "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",     	Ebins, Emin, Emax, 200, -5, 5);

	     TH2D *hJER_per_energy_had_det= new TH2D("hJER_per_energy_had_det", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);
             TH2D *hJER_per_energy_em_det= new TH2D("hJER_per_energy_em_det", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);
             TH2D *hJER_per_energy_none_det= new TH2D("hJER_per_energy_none_det", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);

        TH2D *hJER_per_eDet = new TH2D("hJER_per_eDet", "#DeltaE/E for fixed energies;E_{det};#DeltaE/E",     Ebins, Emin, Emax, 200, -5, 5);
             TH2D *hJER_per_eDet_had_pi = new TH2D("hJER_per_eDet_had_pi", "#DeltaE/E for fixed energies;E_{det};#DeltaE/E",     Ebins, Emin, Emax, 200, -5, 5);
             TH2D *hJER_per_eDet_had_e  = new TH2D("hJER_per_eDet_had_e",  "#DeltaE/E for fixed energies;E_{det};#DeltaE/E",     Ebins, Emin, Emax, 200, -5, 5);
             TH2D *hJER_per_eDet_em_pi  = new TH2D("hJER_per_eDet_em_pi",  "#DeltaE/E for fixed energies;E_{det};#DeltaE/E",     Ebins, Emin, Emax, 200, -5, 5);
             TH2D *hJER_per_eDet_em_e   = new TH2D("hJER_per_eDet_em_e",   "#DeltaE/E for fixed energies;E_{det};#DeltaE/E",     Ebins, Emin, Emax, 200, -5, 5);

             TH2D *hJER_per_eDet_had_det= new TH2D("hJER_per_eDet_had_det", "#DeltaE/E for fixed energies;E_{det};#DeltaE/E",   Ebins, Emin, Emax, 200, -5, 5);
             TH2D *hJER_per_eDet_em_det= new TH2D("hJER_per_eDet_em_det", "#DeltaE/E for fixed energies;E_{det};#DeltaE/E",  Ebins, Emin, Emax, 200, -5, 5);
             TH2D *hJER_per_eDet_none_det= new TH2D("hJER_per_eDet_none_det", "#DeltaE/E for fixed energies;E_{det};#DeltaE/E",   Ebins, Emin, Emax, 200, -5, 5);



        TH2D *hJER_per_distance = new TH2D("hJER_per_distance", "#DeltaE/E for distance;#DeltaR;#DeltaE/E", 	200, 0., 6.5, 
															200, -5, 5);
        TH2D *hJER_per_eta = new TH2D("hJER_per_eta", "#DeltaE/E for distance;#DeltaR;#DeltaE/E",               14,-6.6,-5.2,
                                                                                                                        200, -5, 5);
        TH2D *hEnergy_per_eta = new TH2D("hEnergy_per_eta", "E_{gen} vs. #eta of leading jet", Ebins, Emin, Emax,  14,-6.6,-5.2);

        // Matches.
        TH1D *hMatched = new TH1D("hMatched", "Castor and Gen jets match", 10, -0.5, 9.5);
        TH1D *hUnmatched = new TH1D("hUnmatched", "Castor and Gen jets don't match", Ebins, Emin, Emax);
        TH1D *hJRE = new TH1D("hJRE", "Jet Reconstruction Efficiency", Ebins, Emin, Emax);

        // JER
        TH1D *hJER_all = new TH1D("hJER_all", TString::Format("Jet Energy Resolution, " + string_gen_radius + " " + string_det_radius ), 200, -5, 5);
        TH2D *hJER_per_energy_all = new TH2D("hJER_per_energy_all", "#DeltaE/E for fixed energies;E_{gen};JER",	200.,100.,3000., 	200, -5, 5);
        TH2D *hJER_per_distance_all = new TH2D("hJER_per_distance_all", "#DeltaE/E for distance;#DeltaR;JER", 	200, 0., 6.5, 		200, -5, 5);
	TH2D *hJER_per_eta_all = new TH2D("hJER_per_eta_all", "#DeltaE/E for distance;#DeltaR;JER", 		14,	-6.6,-5.2, 	200, -5, 5);
	TH2D *hEnergy_per_eta_all = new TH2D("hEnergy_per_eta_all", "E_{gen} vs. #eta of leading jet", 				200, 100., 3000.,  	14,-6.6,-5.2);
	
        TH1D *hJER_had_pi = new TH1D("hJER_had_pi", TString::Format("Hadronic (gen and det) Jet Energy Resolution, " + string_gen_radius + " " + string_det_radius ), 200, -5, 5);
        TH1D *hJER_had_e  = new TH1D("hJER_had_e", TString::Format("Hadronic (gen) and em (det) Jet Energy Resolution, " + string_gen_radius + " " + string_det_radius ), 200, -5, 5);
        TH1D *hJER_em_pi = new TH1D("hJER_em_pi", TString::Format("EM (gen) and hadronic (det) Jet Energy Resolution, " + string_gen_radius + " " + string_det_radius ), 200, -5, 5);
        TH1D *hJER_em_e  = new TH1D("hJER_em_e", TString::Format("EM (gen and det) Jet Energy Resolution, " + string_gen_radius + " " + string_det_radius ), 200, -5, 5);
        TH1D *hJER_both_pi = new TH1D("hJER_both_pi", TString::Format("Hybrid (gen) and hadronic (det) Jet Energy Resolution, " + string_gen_radius + " " + string_det_radius ), 200, -5, 5);
        TH1D *hJER_both_e  = new TH1D("hJER_both_e", TString::Format("Hybrid (gen) and em (det) Jet Energy Resolution, " + string_gen_radius + " " + string_det_radius ), 200, -5, 5);

	// RooUnfold.

//	RooUnfoldResponse response (Ebins, Emin, Emax);
	RooUnfoldResponse response 		(hCastorJet_energy, hCastorJet_energy, "response");
	RooUnfoldResponse response_all		(hCastorJet_energy, hCastorJet_energy, "response_all");

	// Jet distance distribution.
	TH1D *hDistance 	= new TH1D("hJetDistance",	"Distance between matched jets;#DeltaR;dN/d#DeltaR", 200, 0., 2.);
	TH1D *hPhiDiff 		= new TH1D("hPhiDiff", 		"Distance in #varphi;#Delta#varphi;dN/d#Delta#varphi", 200, 0., 3.15);
        TH1D *hEtaDiff 		= new TH1D("hEtaDiff", 		"Distance in #eta;#Delta#eta;dN/d#Delta#eta", 200, 0., 0.8);
	TH2D *hEtaPhiDiff 	= new TH2D("hEtaPhiDiff", 	"Distance in #eta and #varphi;#eta;#varphi",200,0.,0.8,200,0.,3.15);
        TH2D *hEtaRDiff 	= new TH2D("hEtaRDiff", 	"Distance in #eta and R;#Delta#eta;#DeltaR",200,0.,0.8,200,0.,3.3);
        TH2D *hPhiRDiff 	= new TH2D("hPhiRDiff", 	"Distance in #eta and #varphi;#Delta#varphi;#DeltaR",200,0.,3.3,200,0.,3.3);

        TH1D *hDistance_all 	= new TH1D("hJetDistance_all", 	"Distance between matched jets;#DeltaR;dN/d#DeltaR", 200, 0., 2.);
        TH1D *hPhiDiff_all 	= new TH1D("hPhiDiff_all", 	"Distance in #varphi;#Delta#varphi;dN/d#Delta#varphi", 200, 0., 3.15);
        TH1D *hEtaDiff_all 	= new TH1D("hEtaDiff_all", 	"Distance in #eta;#Delta#eta;dN/d#Delta#eta", 200, 0., 0.8);
        TH2D *hEtaPhiDiff_all 	= new TH2D("hEtaPhiDiff_all", 	"Distance in #eta and #varphi;#eta;#varphi",200,0.,0.8,200,0.,3.15);
        TH2D *hEtaRDiff_all 	= new TH2D("hEtaRDiff_all", 	"Distance in #eta and R;#Delta#eta;#DeltaR",200,0.,0.8,200,0.,3.3);
        TH2D *hPhiRDiff_all 	= new TH2D("hPhiRDiff_all", 	"Distance in #eta and #varphi;#Delta#varphi;#DeltaR",200,0.,3.3,200,0.,3.3);

	// Pion to electron ratio.
	TH1D *hElectron_energy	= new TH1D("hElectron_energy",	"Energy of electron jets;E_{e} (GeV)",	20, 0., 1000.);
        TH1D *hPion_energy  	= new TH1D("hPion_energy",  	"Energy of Pion jets;E_{#pi} (GeV)",   	20, 0., 1000.);
	TH1D *hPi_e_ratio	= new TH1D("hPi_e_ratio",	"Ratio of pions to electrons;E (GeV)",	20, 0., 1000.);
	
	// Study the jets' energy versus the number of jets.
        TH2D *hNjet_vs_Ejets_gen = new TH2D("hNjet_vs_Ejets_gen", "Number of jets versus E leading jet (gen)",        10, -0.5, 9.5, 50, 100., 3500.);
	TH2D *hNjet_vs_Ejets_det = new TH2D("hNjet_vs_Ejets_det", "Number of jets versus E leading jet (det)", 	10, -0.5, 9.5, 50, 100., 3500.);								
	// Valid gen jets versus Castor jets.
	TH2D *hNumber_of_match_jets = new TH2D("hNumber_of_match_jets", "Castor jet versus Gen jets;N_{Castor};N_{GEN}", 11, -0.5, 10.5, 11, -0.5, 10.5);

	hCastorJet_energy->Sumw2();
	hCastorJet_pt->Sumw2();
	hCastorJet_em->Sumw2();
	hCastorJet_had->Sumw2();
	hCastorJet_fem->Sumw2();
	hCastorJet_fhot->Sumw2();
	hCastorJet_width->Sumw2();
	hCastorJet_depth->Sumw2();
	hCastorJet_sigmaz->Sumw2();
	hCastorJet_ntower->Sumw2();
	hCastorJet_eta->Sumw2();
	hCastorJet_phi->Sumw2();
	hCastorJet_multi->Sumw2();
	
        hCastorJet_energy_response->Sumw2();
        hCastorJet_energy_fakes->Sumw2();
        hCastorJet_energy_misses->Sumw2();
	hCastorJet_Matrix->Sumw2();

	hCastorJet_cf_energy->Sumw2();
        hCastorJet_cf_energy_gen->Sumw2();
	hCastorJet_cf_energy_response->Sumw2();
	hCastorJet_cf_energy_fakes->Sumw2();
	hCastorJet_cf_energy_misses->Sumw2();
	hCastorJet_Matrix->Sumw2();

	hJRE->Sumw2();
	hJER->Sumw2();
	hJER_per_energy->Sumw2();
	hJER_per_distance->Sumw2();
	hMatched->Sumw2();
	hUnmatched->Sumw2();
	
	TIter       next(filelist_); 
	TObjString* fn = 0;
	
	bool isMC = false;
    
	int counter_jer = 0;
	int counter_events = 0;
	int counter_selected_events = 0;
	int counter_files = 0;
	
	std::cout << "start looping over files" << std::endl;
	
	//////////////////////////////
	// create output root file. //
	//////////////////////////////
	
	Char_t filename[200];
	const char* part = currentfile_.Data();
	std::string first(outputname_);
        float etamargin = GenJetContained;         
	
	first += part;

        if( !isData_) { sprintf(filename, date_ +"_STRIPPED_TREE_GEN_" + string_gen_radius + "_DET_" + string_det_radius + "_margin_%f_%i_out_%i_%i_file_%i.root", etamargin, counter_events, totalEvents_, startFile_); }
        else{ sprintf(filename, date_ + "_STRIPPED_TREE_GEN_" + string_det_radius + "_margin_%f_%i_out_%i_%i_file_%i.root", etamargin, counter_events, totalEvents_, startFile_); } 
	TFile* output = new TFile(filename,"RECREATE");	

	////////////////////////////////////////////////////////
	// Prepare the tree with reduced gen jets and events. //
	////////////////////////////////////////////////////////

        std::vector<MyGenJet> good_genJets;// = NULL;
	std::vector<MyCastorJet> *good_CastorJets;

	output->cd();

	TTree * stripped_tree = new TTree("CastorTree","");
	TBranch *b_good_castorJets = stripped_tree->Branch("CastorJets", &good_CastorJets);
        TBranch *b_good_genJets = stripped_tree->Branch("CastorGenJets", &good_genJets);

	cout << "Tree created" << endl;

	
	// start file loop
	while((fn = (TObjString*)next()) && counter_events < totalEvents_) { 
		if( counter_files < startFile_) 		{ counter_files++; continue; }
		if( counter_files > startFile_ + fileSteps) 	{ counter_files++; break; }
		counter_files++;
		currentfile_.Clear();
		currentTFile_->Clear();
		
		currentfile_ = fn->GetString();
		
		std::cout << "opening file " << currentfile_ << " ... " << std::endl;
		
		TStopwatch *timer = new TStopwatch();
		TThread *fileopener;
		fileopener = new TThread("fileopener",(void(*) (void *))&OpenROOTFile,(JetAnalyzer_stripTheTree*) this);
		TThread::SetCancelAsynchronous();
		TThread::SetCancelOn();
		fileopener->Run();
		
		bool fileopen = false;
		bool timeout = true;
		while (timer->RealTime() < 1200) {
			if (fileopener->GetState() != TThread::kRunningState) {
				// file open thread is not running anymore, was it successful?
				// set the static TFile to the instance one
				setCurrentTFile();
				timeout = false;
				break;
			} 
			
			timer->Continue();
			gSystem->Sleep(1000);
		}
		
		// after the time, stop & delete the open thread
		TThread::Kill("fileopener");
		TThread::Delete(fileopener);
		delete fileopener;
		delete timer;
		
		if (timeout) {
			std::cout << "file open thread has timed out, go to next file in the list" << std::endl;
			continue;
		} else {
			// is this instance open?
			if (currentTFile_->IsOpen()) {
				std::cout << "the currentTFile_ is correctly opened" << std::endl;
				fileopen = true;
			}
		}
		
		if (!fileopen) {
			std::cout << "no timeout, but the file could not be opened correctly: going to the next file in the list" << std::endl;
			continue;
		}
		
		
		//////////////////////////////////////////////////
		// Get tree from the files and define all branches
		//////////////////////////////////////////////////
		
		std::cout << "We are working with Det " << det_radius_ << " and Gen " << gen_radius_ << endl;

		// get tree from file
		TTree *tree = new TTree("CastorTree","");
		std::cout << "Extract tree" << endl;
		currentTFile_->GetObject("castortree/CastorTree",tree);
		std::cout << "Tree extracted" << endl;	
	
		
		currentTFile_->cd();
		
		// define objects and branches
		MyEvtId *evtid = NULL;
        MyGenKin *evtkin = NULL;
		MyBeamSpot *BeamSpot = NULL;
		MyHLTrig *HLTrig = NULL;
		MyL1Trig *L1Trig = NULL;
		std::vector<MyVertex> *Vertices = NULL;
		std::vector<MyCastorRecHit> *CastorRecHits = NULL;
		std::vector<MyCastorTower> *CastorTowers = NULL;
		std::vector<MyCastorJet> *CastorJets = NULL;
		std::vector<MyJet> *PFJets = NULL;
		std::vector<MyGenPart> *genParts = NULL;
		std::vector<MyCaloTower> *caloTowers = NULL;
		std::vector<MyGenJet> *genJets = NULL;

		cout << "First batch of stuff set" << endl;

//                std::vector<MyGenJet> *good_genJets = NULL;

		std::vector<MyGenJet> *chargedGenJets = NULL;
		std::vector<MyTrackJet> *trackJets = NULL;
		  // Radius-energy study.
		std::vector<MyCastorJet> *CastorJets_ak3 = NULL;
                std::vector<MyCastorJet> *CastorJets_ak5 = NULL;
                std::vector<MyCastorJet> *CastorJets_ak7 = NULL;
                std::vector<MyGenJet> *genJets_ak3 = NULL;
                std::vector<MyGenJet> *genJets_ak4 = NULL;
                std::vector<MyGenJet> *genJets_ak5 = NULL;
                std::vector<MyGenJet> *genJets_ak7 = NULL;
	
		cout << "Prepare branches" << endl;

		
		TBranch *b_evtid = tree->GetBranch("EvtId");
		cout << "Passed event id" << endl;
        TBranch *b_evtkin = NULL;
        if (!isData_) b_evtkin = tree->GetBranch("GenKin");
		TBranch *b_BeamSpot = tree->GetBranch("beamSpot");
		TBranch *b_HLTrig = tree->GetBranch("HLTrig");
		TBranch *b_L1Trig = tree->GetBranch("L1Trig");
		TBranch *b_vertices = tree->GetBranch("primaryVertex");
		TBranch *b_castorrechits = tree->GetBranch("castorRecHit");
		TBranch *b_castortowers = tree->GetBranch("castorTower");
		TBranch *b_castorjets = tree->GetBranch("castorJet");
	
		TBranch *b_castorjets_ak3 = tree->GetBranch("ak3castorJet");	b_castorjets_ak3->SetAddress(&CastorJets_ak3);
                TBranch *b_castorjets_ak5 = tree->GetBranch("ak5castorJet");    b_castorjets_ak5->SetAddress(&CastorJets_ak5);
                TBranch *b_castorjets_ak7 = tree->GetBranch("ak7castorJet");    b_castorjets_ak7->SetAddress(&CastorJets_ak7);

		if( string_det_radius.Contains("ak3") ){ b_castorjets = tree->GetBranch("ak3castorJet"); }
                if( string_det_radius.Contains("ak4") ){ b_castorjets = tree->GetBranch("ak4castorJet"); }
                if( string_det_radius.Contains("ak5") ){ b_castorjets = tree->GetBranch("ak5castorJet"); }
                if( string_det_radius.Contains("ak7") ){ b_castorjets = tree->GetBranch("ak7castorJet"); }
//		else{ b_castorjets = tree->GetBranch("castorJet"); }
		
//                TBranch *b_castorjets = tree->GetBranch("ak7castorJet");
		TBranch *b_PFJets = tree->GetBranch("pfJet");
		TBranch *b_genParts = NULL;
		if (!isData_) b_genParts = tree->GetBranch("GenPart");
		TBranch *b_caloTowers = tree->GetBranch("caloTower");
		TBranch *b_genJets = NULL;
		TBranch *b_chargedGenJets = NULL;
//		if (!isData_) b_genJets = tree->GetBranch("GenJet");
                if (!isData_ && string_gen_radius.Contains("ak3")) b_genJets = tree->GetBranch("ak3GenJet");
                if (!isData_ && string_gen_radius.Contains("ak4")) b_genJets = tree->GetBranch("ak4GenJet");
                if (!isData_ && string_gen_radius.Contains("ak5")) b_genJets = tree->GetBranch("ak5GenJet");
                if (!isData_ && string_gen_radius.Contains("ak7")) b_genJets = tree->GetBranch("ak7GenJet");		
		if (!isData_) b_chargedGenJets = tree->GetBranch("ChargedGenJet");
		TBranch *b_trackJets = tree->GetBranch("trackJet");

                TBranch *b_genJets_ak3 = NULL;
                TBranch *b_genJets_ak4 = NULL;
                TBranch *b_genJets_ak5 = NULL;
                TBranch *b_genJets_ak7 = NULL;

		cout << "Ready to attack genjets" << endl;

		if(!isData_){
		  b_genJets_ak3 = tree->GetBranch("ak3GenJet");
                  b_genJets_ak4 = tree->GetBranch("ak4GenJet");
                  b_genJets_ak5 = tree->GetBranch("ak5GenJet");
                  b_genJets_ak7 = tree->GetBranch("ak7GenJet");
		}
		
		cout << "Attacked gen jets" << endl;

		b_evtid->SetAddress(&evtid);
        	if (!isData_) b_evtkin->SetAddress(&evtkin);
		b_BeamSpot->SetAddress(&BeamSpot);
		b_HLTrig->SetAddress(&HLTrig);
		b_L1Trig->SetAddress(&L1Trig);
		b_vertices->SetAddress(&Vertices);
		b_castorrechits->SetAddress(&CastorRecHits);
		b_castortowers->SetAddress(&CastorTowers);
		cout << "Before" << endl;
		b_castorjets->SetAddress(&CastorJets);

		cout << "After" << endl;
		b_PFJets->SetAddress(&PFJets);
		if (!isData_) b_genParts->SetAddress(&genParts);
		b_caloTowers->SetAddress(&caloTowers);
		if (!isData_) b_genJets->SetAddress(&genJets);
		if (!isData_) b_chargedGenJets->SetAddress(&chargedGenJets);
		b_trackJets->SetAddress(&trackJets);
		
		
		
	/*	
		if( counter_events == 0){
		  cout << "Let's clone this tree!" << endl;
		  output->cd();
		  //stripped_tree = tree->CloneTree(0); 
		  //cout << "Castorjets" << endl;
		  //stripped_tree->Branch("CastorDetJets",good_CastorJets);
	          //cout << "Genjets" << endl;
		  //stripped_tree->Branch("CastorGenJets",&good_genJets);		  
  
		  cout << "This tree's been cloned!" << endl;
		}		
	*/	
		
		int Nevents = tree->GetEntriesFast();
		std::cout << "file opened, events in this file = " << Nevents << std::endl;
		totalevents += Nevents;
		
		// start event loop
		for (int i=0; i<Nevents && i < totalEvents_;i++) {

			counter_events++;
			if( counter_events%1000== 0){ cout << "\t" << counter_events << "\tpassed" << endl; }
            
			bool passedHadronCuts = false;
			bool passedDetectorCuts = false;
			
			double xix = 10;
			double xiy = 10;
			double xi = 10;
			double xidd = 10e10;
			double ymax = -1;
            
	
			
			/////////////////////////////////////////
			// Do stuff before filters
			/////////////////////////////////////////
						
			
			b_evtid->GetEntry(i);
			b_HLTrig->GetEntry(i);
			b_L1Trig->GetEntry(i);
			b_vertices->GetEntry(i);
			b_caloTowers->GetEntry(i);
			//b_castorrechits->GetEntry(i);
			b_castortowers->GetEntry(i);


			////////////////////////////////////////////////////////////////////
			// Count the number of events with Castor jets on gen. level without looking at detector level.
			////////////////////////////////////////////////////////////////////
                        if(!isData_){

			  // Loop over the generator level jets.
		 	  // Check if they lie in CASTOR, and if they pass a certain threshold.
			  bool  castor_0GeV = false;
			  bool  castor_150GeV = false;
	  		  bool castor_300GeV = false;

			  b_genJets->GetEntry(i);
			  total_events_MC++;

 			  for(int jet = 0; jet < genJets->size(); jet++){
			    MyGenJet currentJet = (*genJets)[jet];
			    if( currentJet.Eta() < -5.2  && currentJet.Eta() > -6.6  ){
			      castor_0GeV = true;
			      double castor_gen_energy = currentJet.Energy();
			      hGenJet_energy_noCuts->Fill( castor_gen_energy );
			      if( castor_gen_energy > 300. ){	
				castor_300GeV = true; 
				break;
			      }
			      else if( castor_gen_energy > 150. ){
				castor_150GeV = true; 
			      }			    
			    }// Jet in Castor.
			  }// Loop over genjets.
		
			  if( castor_300GeV ){
			    castor_300GeV_events_MC++;
			    castor_150GeV_events_MC++;
			    castor_events_MC++;
			  }
			  else if( castor_150GeV ){
			    castor_150GeV_events_MC++;
			    castor_events_MC++;
			  }
			  else if( castor_0GeV ){
			    castor_events_MC++;
			  }			
			} // !isData_

			////////////////////////////////////////////////////////////////////
			// Done.
			////////////////////////////////////////////////////////////////////

			if( Vertices->size() != 1) continue;
			
			// only process a certain run
			//if (evtid->Run != 135521) continue; // go to next event
			
			/////////////////////////////////////////
			// Filter the results
			/////////////////////////////////////////
			
			if (isData_) {
			
				// filter results
				// filter on phys declared bit
				bool physDeclresult = HLTrig->HLTmap["physDeclpath"];
				
				// filter on castor invalid data
				bool castorInvalidDataFilterresult = HLTrig->HLTmap["castorInvalidDataFilterpath"];
				
				// filter out scraping events
				bool noscrapingresult = HLTrig->HLTmap["noscrapingpath"];
				
				bool gooddata = physDeclresult && castorInvalidDataFilterresult && noscrapingresult;
				
				// L1 filter

				bool L1_BX = L1Trig->fTechDecisionBefore[0];
				bool L1_Veto = !L1Trig->fTechDecisionBefore[36] && !L1Trig->fTechDecisionBefore[37] && !L1Trig->fTechDecisionBefore[38] && !L1Trig->fTechDecisionBefore[39];
				
				bool L1_BSC = L1Trig->fTechDecisionBefore[40] || L1Trig->fTechDecisionBefore[41]; // default value for 900 GeV or 7000 GeV 
				
				bool HLT_BSC = true; // default value for 900 GeV or 7000 GeV (do we need an HLT here?)
				
				bool TriggerSelection = false;
                
				TriggerSelection = L1_BX && L1_Veto && L1_BSC && HLT_BSC;
				
				// ask for activity in HF
				bool HF_Activity = true;
				bool HFplus = false;
				bool HFminus = false;
				for (unsigned int itow=0;itow<caloTowers->size();itow++) {
					MyCaloTower mytow = (*caloTowers)[itow];
					if (mytow.hasHF && mytow.Energy() > 4. && mytow.zside == 1 && mytow.Eta() > 3.23 && mytow.Eta() < 4.65) HFplus = true;
					if (mytow.hasHF && mytow.Energy() > 4. && mytow.zside == -1 && mytow.Eta() < -3.23 && mytow.Eta() > -4.65) HFminus = true;
				}
				if (!HFplus || !HFminus) HF_Activity = false; // real HF condition... 
				
				// ask for activity in CASTOR
				bool CASTOR_Activity = false;
				// if at least 1 tower is above noise threshold, set CASTOR activity to true
				if (CastorTowers->size() > 0) CASTOR_Activity = true;
				
				// get vertex info
				// do the oneGoodVertexFilter
				bool wehaveGoodVertex = false;
				for (unsigned int iVert=0;iVert<Vertices->size();iVert++) {
					MyVertex vertex = (*Vertices)[iVert];
					if (vertex.isGoodVertex) wehaveGoodVertex = true;
				}
				
				if (gooddata && TriggerSelection && HF_Activity && CASTOR_Activity && wehaveGoodVertex) passedDetectorCuts = true;
				
			} else {
				
				// MC
				
				// ask for activity in HF
				bool HF_Activity = true; 
				bool HFplus = false;
				bool HFminus = false;
				for (unsigned int itow=0;itow<caloTowers->size();itow++) {
					MyCaloTower mytow = (*caloTowers)[itow];
					if (mytow.hasHF && mytow.Energy() > 4. && mytow.zside == 1 && mytow.Eta() > 3.23 && mytow.Eta() < 4.65) HFplus = true;
					if (mytow.hasHF && mytow.Energy() > 4. && mytow.zside == -1 && mytow.Eta() < -3.23 && mytow.Eta() > -4.65) HFminus = true;
				}
				if (!HFplus || !HFminus) HF_Activity = false; // real HF condition... 
				
				// ask for activity in CASTOR
				bool CASTOR_Activity = false;
				// if at least 1 tower is above noise threshold, set CASTOR activity to true
				if (CastorTowers->size() > 0) CASTOR_Activity = true;
                
				// get vertex info
				// do the oneGoodVertexFilter
				bool wehaveGoodVertex = false;
				for (unsigned int iVert=0;iVert<Vertices->size();iVert++) {
					MyVertex vertex = (*Vertices)[iVert];
					if (vertex.isGoodVertex) wehaveGoodVertex = true;
				}
				
				if (HF_Activity && CASTOR_Activity && wehaveGoodVertex) passedDetectorCuts = true;
                
			}
				
			/////////////////////////////////////////
			// Do stuff after filters
			/////////////////////////////////////////
			
			// get all the remaining branch entries
			b_trackJets->GetEntry(i);
                        if(!isData_){
			  b_genJets->GetEntry(i);
                          b_chargedGenJets->GetEntry(i);
                        }
				
			if (passedDetectorCuts) {
				
				
				/////////////////////////////////////////
				// Start Nvertex == 1 part of the code 
				/////////////////////////////////////////

//				cout << "Event\t" << counter_events << "\tpassed cuts\t" << Vertices->size() << "\tvertices\t" << endl;				
				// only fill the histograms when there's 1 vertex (filter out pile-up)
//				if (Vertices->size() == 1) {				
	
					// get event id stuff
					if( ((i+1) % 10000) == 0) cout << " run " << evtid->Run << " isData = " << evtid->IsData << " lumiblock " << 
						evtid->LumiBlock << " event " << evtid->Evt << endl; 
					if (!evtid->IsData) isMC = true;
								
					good_genJets.clear();
					b_castorjets->GetEntry(i);

					// analyse castor jets
					int NCastorJets = 0;

					vector<double> det_casjet, gen_casjet_energy;
					b_castorjets->GetEntry(i);

					good_CastorJets = CastorJets;		

					// -------------------------
					// AVS - no Central jet required					
					if(!isData_){
					  for(int jet = 0; jet < genJets->size(); jet++){
					    MyGenJet currentJet = (*genJets)[jet];
					    if( currentJet.Eta() < (-5.2 - GenJetContained) && currentJet.Eta() > (-6.6 + GenJetContained) ){
					      good_genJets.push_back( currentJet );
					    }
					  }
					}

// AVS					if( CastorJets->size() > 0 && ( good_genJets.size() > 0 || isData_ )){
					if( CastorJets->size() > 0 ){
					  stripped_tree->Fill();
					  counter_selected_events++;
					}
	
					/**********************
					 * ********************
					 * *******************/
						
					// end of event, print status
					if( ((i+1) % 10000) == 0) std::cout << i+1 <<"events done in file " << it << std::endl;
					//totalevents++;
					
//				} // end if statement for 1 vertex
								
			} // end if passed detector cuts
			
			// combined hadron and detector level code	
		} // end event loop
		
		delete tree;
		currentTFile_->Close();
		it++;
	} // end file loop
	
	std::cout << "file loop has ended with " << counter_selected_events << " / " << counter_events << " selected events" << std::endl;
    
    currentTFile_->Close();
	
	
    // check all histo's for overflow
	
	hf_.checkFlow(hCASTOReflow);
	hf_.checkFlow(hCASTORTowerMulti);
    
    // check distribution of each channel on under or overflow
	for (int icha=0;icha<224;icha++) {
		hf_.checkFlow(hCASTOReflow_channel[icha]);
	}
	
	hf_.checkFlow(hCastorJet_energy);
	hf_.checkFlow(hCastorJet_pt);
	hf_.checkFlow(hCastorJet_em);
	hf_.checkFlow(hCastorJet_had);
	hf_.checkFlow(hCastorJet_fem);
	hf_.checkFlow(hCastorJet_fhot);
	hf_.checkFlow(hCastorJet_width);
	hf_.checkFlow(hCastorJet_depth);
	hf_.checkFlow(hCastorJet_sigmaz);
	hf_.checkFlow(hCastorJet_ntower);
	hf_.checkFlow(hCastorJet_eta);
	hf_.checkFlow(hCastorJet_phi);
	hf_.checkFlow(hCastorJet_multi);
    

	// Create response matrix including misses and fakes.
	hCastorJet_energy_complete_response->Add( hCastorJet_energy_response);
        hCastorJet_energy_complete_response->Add( hCastorJet_energy_misses);
        hCastorJet_energy_complete_response->Add( hCastorJet_energy_fakes);
   
    // write all histo's to file
    
	std::cout << "total number of events = " << totalevents << " from " << it << " file(s)" << endl;

	output->cd();	
  	// Write away useful numbers.
	TTree *tree_numbers = new TTree("useful_numbers", "useful_numbers");
	/*
	tree_numbers->Branch("castor_events", &castor_events, "castor_events/I");
	tree_numbers->Branch("castor_150GeV_events", &castor_150GeV_events, "castor_150GeV_events/I");
	tree_numbers->Branch("castor_300GeV_events", &castor_300GeV_events, "castor_300GeV_events/I");
	*/
	tree_numbers->Branch("total_events", &totalevents, "total_events/I");
	tree_numbers->Branch("total_events_MC", &total_events_MC, "total_events_MC/I");
	tree_numbers->Branch("castor_events_MC", &castor_events_MC, "castor_events_MC/I");
	tree_numbers->Branch("castor_150GeV_events_MC", &castor_150GeV_events_MC, "castor_150GeV_events_MC/I");
	tree_numbers->Branch("castor_300GeV_events_MC", &castor_300GeV_events_MC, "castor_300GeV_events_MC/I");
	tree_numbers->Fill();
	tree_numbers->Write();
	

	hGenJet_energy_noCuts->Write();
		
	//////////////////////////////////////////
	// Save all your histograms in a root file
	//////////////////////////////////////////
	
	// detector level histograms
   /* 
	// eflow histos
	hCASTORTowerMulti->Write();
	hCASTOReflow->Write();
	h2CASTOReflow_grid->Write();
	
	for (int icha=0;icha<224;icha++) {
		hCASTOReflow_channel[icha]->Write();
	}
*
*/
/*
	std::vector<TString> Useless_branches;
        Useless_branches.push_back("EvtId");

        TBranch *brnch;

	for(int branch_n = 0; branch_n < Useless_branches.size(); branch_n++){
                brnch = stripped_tree->GetBranch("EvtId");
		cout << "Got the Bracnhg" << endl;
		stripped_tree->GetListOfBranches()->Remove(brnch);
		cout << "Removed the branch" << endl;
	}
	cout << "To write?" << endl;
*/
	stripped_tree->Write();

	output->Close();
	std::cout << "file " << filename << " created." << std::endl;

        cout << totalevents << "\tevents" << endl;
	cout << counter_jer << "\tJER events" << endl;

}
	
void JetAnalyzer_stripTheTree::AfterLoopCalculations(TString file) {

	///////////////////////////////////////////////////////
	// perform calculations after event by event filling //
	///////////////////////////////////////////////////////
    
    // get the histograms from the file
    // get all the histograms
    HistoRetriever histogetter;
    std::vector<TH1D*> histovector = histogetter.getHistos(file);
    //std::vector<TH2D*> h2Dhistovector = histogetter_.get2DHistos(inputdir+file);
    
    char name [100];
    
    TH1D *hCASTOReflow_channel[224];
	for (int i=0;i<224;i++) {
		sprintf(name,"hCASTOReflow_channel_%d",i+1);
		hCASTOReflow_channel[i] = new TH1D(hf_.getHistoByName(histovector,name)); 
	}
    
    TH1D *hCASTOReflow = new TH1D(hf_.getHistoByName(histovector,"hCASTOReflow"));
	
	std::cout << "starting AfterLoop calculations" << std::endl;
	
	// get mean and error's from all channels and put it in one histo
    TH1D *hCASTOReflow_channels = new TH1D("hCASTOReflow_channels","average energy in used channels",224,1,225);
	for (int icha=0;icha<224;icha++) {
		hCASTOReflow_channels->SetBinContent(icha+1,hCASTOReflow_channel[icha]->GetMean());
		hCASTOReflow_channels->SetBinError(icha+1,hCASTOReflow_channel[icha]->GetMeanError());
	}
  	
	std::cout << "Mean CASTOR energy flow in first 5 modules = " << hCASTOReflow->GetMean() << " +/- " << hCASTOReflow->GetMeanError() << std::endl;
	
	//////////////////////////////////////////////////
	
	// create output root file
	Char_t filename[200];
	sprintf(filename,"AfterLoop_%s",file.Data());
	TFile* output = new TFile(filename,"RECREATE");
	output->cd();
	
	//////////////////////////////////////////
	// Save all your histograms in a root file
	//////////////////////////////////////////
        
  	//hCASTOReflow_channels->Write();
		
	output->Close();
	std::cout << "file " << filename << " created." << std::endl;
	
}

TString JetAnalyzer_stripTheTree::getOutputFile() {
    return LoopOutputFile_;
}

TString JetAnalyzer_stripTheTree::getInputDir() {
	return inputdir_;
}

TString JetAnalyzer_stripTheTree::getCurrentFile() {
	return currentfile_;
}

void JetAnalyzer_stripTheTree::setCurrentTFile() {
	currentTFile_ = currentStaticTFile_;
}

void* JetAnalyzer_stripTheTree::OpenROOTFile(JetAnalyzer_stripTheTree* arg) {
	currentStaticTFile_ = TFile::Open(arg->getInputDir()+arg->getCurrentFile(),"READ");
	return 0;
}

int JetAnalyzer_stripTheTree::posLeadingGenJet(std::vector<MyGenJet> JetVector,double etacut,double minptcut) {
	
	// search for leading jets
	int posLeadingChargedGenJet = -1;

	double tempptchargedgenjet = 0;
	for (unsigned int ijet=0;ijet<JetVector.size();ijet++) {
		if (JetVector[ijet].Pt() > minptcut && fabs(JetVector[ijet].Eta()) < etacut) {
			if (JetVector[ijet].Pt() > tempptchargedgenjet) {
				tempptchargedgenjet = JetVector[ijet].Pt();
				posLeadingChargedGenJet = ijet;
			}
		}
	}
	
	return posLeadingChargedGenJet;
}

int JetAnalyzer_stripTheTree::posLeadingTrackJet(std::vector<MyTrackJet> JetVector,double etacut,double minptcut) {
	
	// search for leading jets
	int posLeadingTrackJetresult = -1;
	double temppttrack = 0;
	for (unsigned int ijet=0;ijet<JetVector.size();ijet++) {
		if (JetVector[ijet].pt_raw > minptcut && fabs(JetVector[ijet].eta_raw) < etacut && JetVector[ijet].pv) {
			if (JetVector[ijet].pt_raw > temppttrack) {
				temppttrack = JetVector[ijet].pt_raw;
				posLeadingTrackJetresult = ijet;
			}
		}
	}
	
	return posLeadingTrackJetresult;
}

