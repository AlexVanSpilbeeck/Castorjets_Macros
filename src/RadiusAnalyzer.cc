////////////////////////////////////////
//////// New CMSSW_4_2_X version ///////
////////////////////////////////////////

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "RadiusAnalyzer.h"
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
#define jetEThreshold 100. 
#define EbinWidth 5.
#define EbinWidth_rel 1.4

#define jet_distance 10
#define jet_distance_string "JER_allPairs__JetSorted_ak5"
#define GenJetRadius "ak7"
#define DetJetRadius "ak7"
#define GenJetContained 0.
#define nEvents 10000
#define PI 3.14159265359

#include "../../../RooUnfold-1.1.1/src/RooUnfold.h"
//#include "../../../RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#endif



TFile *RadiusAnalyzer::currentStaticTFile_ = new TFile();

RadiusAnalyzer::RadiusAnalyzer(TString inputdir, TObjArray* filelist, bool isData, const char* outputname) {
    
	std::cout << "constructing RadiusAnalyzer class..." << std::endl;
	
    inputdir_ = inputdir;
    filelist_ = filelist;
    isData_ = isData;
    outputname_ = outputname;
	
	std::cout << "initializing basic variables..." << std::endl;
    
    // initialize basic variables

    
    LoopOutputFile_ = "";

	currentfile_ = "";
	currentTFile_ = new TFile();
	
	std::cout << "all initialisations done, class constructed!" << std::endl;
    
}

RadiusAnalyzer::~RadiusAnalyzer() { }

void RadiusAnalyzer::Loop() {

#ifdef __CINT__
  gSystem->Load("../../../RooUnfold-1.1.1/libRooUnfold.so");
#endif
	
	
	std::cout << " RadiusAnalyzer Loop function is started " << std::endl;
	
	TString tstring = outputname_;
	std::cout << " TString outputname_ = " << tstring << std::endl;
	TString string_tag = jet_distance_string;
	TString string_det_radius = DetJetRadius;
        TString string_gen_radius = GenJetRadius;
 
	
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
	
        // AVS - number of bins.
	/* We get binwidth from peak JER. */
	double Emin = jetEThreshold + 100., Emax = 3000.;
	int Ebins = static_cast<int> ( (Emax - Emin)/EbinWidth );


	/* If we want variable Ebins. */
	double current_lowE = 200.;
	int count_bincenters = 3;

	vector<double> binEdges;
	binEdges.push_back( -50.);
	binEdges.push_back( 0. );
	binEdges.push_back( 50. );
	binEdges.push_back( jetEThreshold );
	binEdges.push_back( current_lowE );
	while( current_lowE/EbinWidth_rel < Emax ){

	  binEdges.push_back( current_lowE * EbinWidth_rel );

 	  cout << "\tLower edge\t" << count_bincenters << "\t" << current_lowE << "\tBinwidth\t" << current_lowE *0.2 << endl;

	  current_lowE += current_lowE*0.2;
	  count_bincenters++;
	}

	Ebins = binEdges.size();
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
        TH1D *hCastorJet_energy = new TH1D("hCastorJet_energy","CastorJet energy distribution",Ebins-1, Ebins_var);
        TH1D *hGenJet_energy = new TH1D("hGenJet_energy","GenJet energy distribution",Ebins, Emin, Emax);
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
        TH1D *hCastorJet_cf_energy = new TH1D("hCastorJet_cf_energy", "CastorJet energy in with leading central jet",Ebins+10,Emin-10*EbinWidth,Emax);
        TH1D *hCastorJet_cf_energy_gen = new TH1D("hCastorJet_cf_energy_gen", "CastorJet energy in with leading central jet",Ebins+10,Emin-10*EbinWidth,Emax);
	TH2D *hCastorJet_cf_energy_response = new TH2D("hCastorJet_cf_energy_response", "CastorJet energy response matrix for leading central jet",Ebins+10,Emin-10*EbinWidth,Emax,Ebins+10,Emin-10*EbinWidth,Emax);
        TH1D *hCastorJet_cf_energy_fakes = new TH1D("hCastorJet_cf_energy_fakes", "CastorJet energy in with leading central jet - fakes",Ebins+10,Emin-10*EbinWidth,Emax);
        TH1D *hCastorJet_cf_energy_misses = new TH1D("hCastorJet_cf_energy_misses", "CastorJet energy in with leading central jet - misses",Ebins+10,Emin-10*EbinWidth,Emax);

	// All Casto jets response matrix.
//        TH2D *hCastorJet_energy_response = new TH2D("hCastorJet_energy_response","CastorJet energy distribution  - Response",Ebins+10,Emin-10*EbinWidth,Emax,Ebins+10,Emin-10*EbinWidth,Emax);
//        TH1D *hCastorJet_energy_fakes = new TH1D("hCastorJet_energy_fakes","CastorJet energy distribution - Fakes",Ebins+10,Emin-10*EbinWidth,Emax);
//        TH1D *hCastorJet_energy_misses = new TH1D("hCastorJet_energy_misses","CastorJet energy distribution - Misses",Ebins+10,Emin-10*EbinWidth,Emax);

	// Response with variable binning.
        TH2D *hCastorJet_energy_response = new TH2D("hCastorJet_energy_response","CastorJet energy distribution  - Response",Ebins-1, Ebins_var,Ebins-1, Ebins_var);
        TH1D *hCastorJet_energy_fakes = new TH1D("hCastorJet_energy_fakes","CastorJet energy distribution - Fakes", Ebins-1, Ebins_var);
        TH1D *hCastorJet_energy_misses = new TH1D("hCastorJet_energy_misses","CastorJet energy distribution - Misses",Ebins-1, Ebins_var);

//	TH2D *hCastorJet_energy_ratio = new TH2D("hCastorJet_energy_ratio", "Ratio of generator to detector", 1000, 0., 10.,Ebins-1, Ebins_var);
        TH2D *hCastorJet_energy_ratio = new TH2D("hCastorJet_energy_ratio", "Ratio of generator to detector", 1000, 0., 10.,Ebins, Emin, Emax);

	TH2D *hCastorJet_cf_Matrix = new TH2D("hCastorJet_cf_responseMatrix", "Central-forward jet Response Matrix",Ebins+10,Emin-10*EbinWidth,Emax,Ebins+10,Emin-10*EbinWidth,Emax);
//        TH2D *hCastorJet_Matrix = new TH2D("hCastorJet_responseMatrix", "Castor jet Response Matrix",Ebins+10,Emin-10*EbinWidth,Emax,Ebins+10,Emin-10*EbinWidth,Emax);
        TH2D *hCastorJet_Matrix = new TH2D("hCastorJet_responseMatrix", "Castor jet Response Matrix",Ebins-1, Ebins_var,Ebins-1, Ebins_var);

	// Number of trackjets.
	TH2D *hTrackjets_2D_number = new TH2D("hTrackjets_2D_number", "Trackjets vs. Charged Gen Jet", 17,0.,17., 17,0.,17.);
        TH2D *hTrackjets_2D_pt = new TH2D("hTrackjets_2D_pt", "Trackjets vs. Charged Gen Jet", 30,0.,30., 30,0.,30.);

        TH1D *hCastorJet_energy_gen = new TH1D("hCastorJet_energy_gen","CastorJet energy distribution",Ebins-1, Ebins_var);
        TH1D *hCastorJet_pt_gen = new TH1D("hCastorJet_pt_gen","CastorJet pt distribution",500,0,10.);

	// Efficiency control.
	TH1D *hJER = new TH1D("hJER", "Jet Energy Resolution", 200, -5, 5);
        TH2D *hJER_per_energy = new TH2D("hJER_per_energy", "Jet Energy Resolution for fixed energies;E_{gen};JER", 	200.,100.,3000.,  
															200, -5, 5);
        TH2D *hJER_per_distance = new TH2D("hJER_per_distance", "Jet Energy Resolution for distance;#DeltaR;JER", 	200, 0., 6.5, 
															200, -5, 5);
//	TH1D *hMatched = new TH1D("hMatched", "Castor and Gen jets match", Ebins, Emin, Emax);
        TH1D *hMatched = new TH1D("hMatched", "Castor and Gen jets match", 10, -0.5, 9.5);
        TH1D *hUnmatched = new TH1D("hUnmatched", "Castor and Gen jets don't match", Ebins, Emin, Emax);
        TH1D *hJRE = new TH1D("hJRE", "Jet Reconstruction Efficiency", Ebins, Emin, Emax);

        TH1D *hJER_2 = new TH1D("hJER_2", "Jet Energy Resolution", 200, -5, 5);
        TH2D *hJER_per_energy_2 = new TH2D("hJER_per_energy_2", "Jet Energy Resolution for fixed energies;E_{gen};JER",	200.,100.,3000.,
															200, -5, 5);
        TH2D *hJER_per_distance_2 = new TH2D("hJER_per_distance_2", "Jet Energy Resolution for distance;#DeltaR;JER", 	200, 0., 6.5,
															200, -5, 5);
	TH2D *hJER_per_eta = new TH2D("hJER_per_eta", "Jet Energy Resolution for distance;#DeltaR;JER", 		14,-6.6,-5.2,
															200, -5, 5);
	TH2D *hEnergy_per_eta = new TH2D("hEnergy_per_eta", "E_{gen} vs. #eta of leading jet", 200, 100., 3000.,  14,-6.6,-5.2);
	
	// RooUnfold.

//	RooUnfoldResponse response (Ebins, Emin, Emax);
	RooUnfoldResponse response (hCastorJet_energy, hCastorJet_energy, "response");
	RooUnfoldResponse response_onlyMatches(Ebins, Emin, Emax, "response_onlyMatches");
        RooUnfoldResponse response_cf (Ebins, Emin, Emax, "response_cf");
        RooUnfoldResponse response_cf_onlyMatches (Ebins, Emin, Emax, "response_cf_onlyMatches");

	// Jet distance distribution.
	TH1D *hDistance = new TH1D("hJetDistance", "Distance between matched jets;#DeltaR;dN/d#DeltaR", 200, 0., 2.);
	TH1D *hPhiDiff = new TH1D("hPhiDiff", "Distance in #varphi;#Delta#varphi;dN/d#Delta#varphi", 200, 0., 3.15);
        TH1D *hEtaDiff = new TH1D("hEtaDiff", "Distance in #eta;#Delta#eta;dN/d#Delta#eta", 200, 0., 0.8);
	TH2D *hEtaPhiDiff = new TH2D("hEtaPhiDiff", "Distance in #eta and #varphi;#eta;#varphi",200,0.,0.8,200,0.,3.15);
        TH2D *hEtaRDiff = new TH2D("hEtaRDiff", "Distance in #eta and R;#Delta#eta;#DeltaR",200,0.,0.8,200,0.,3.3);
        TH2D *hPhiRDiff = new TH2D("hPhiRDiff", "Distance in #eta and #varphi;#Delta#varphi;#DeltaR",200,0.,3.3,200,0.,3.3);

	// Study the jets' energy versus the number of jets.
        TH2D *hNjet_vs_Ejets_gen = new TH2D("hNjet_vs_Ejets_gen", "Number of jets versus E leading jet (gen)",        10, -0.5, 9.5, 50, 100., 3500.);
	TH2D *hNjet_vs_Ejets_det = new TH2D("hNjet_vs_Ejets_det", "Number of jets versus E leading jet (det)", 	10, -0.5, 9.5, 50, 100., 3500.);								
	// Valid gen jets versus Castor jets.
	TH2D *hNumber_of_match_jets = new TH2D("hNumber_of_match_jets", "Castor jet versus Gen jets;N_{Castor};N_{GEN}", 11, -0.5, 10.5, 11, -0.5, 10.5);

        // Study of the energetic content of a jet versus its radius.
        TH1D *hNjet_3_gen = new TH1D("hNjet_3_gen", "Jets: Gen (ak3)", 11, -0.5, 10.5);
        TH1D *hNjet_5_gen = new TH1D("hNjet_5_gen", "Jets: Gen (ak5)", 11, -0.5, 10.5);
        TH1D *hNjet_7_gen = new TH1D("hNjet_7_gen", "Jets: Gen (ak7)", 11, -0.5, 10.5);

	TH2D *hJet_7to3_gen 	= new TH2D("hJet_7to3_gen",	"Jets: energy of Gen(ak7) vs. Gen(ak3);E_{Gen(ak7)} (GeV);E_{Gen(ak3)}/E_{Gen(ak7)}", 50, 100., 3500., 20, 0., 1.);
          TH1D *hDistance_7to3_gen= new TH1D("hJetDistance_7to3_gen", 	"Distance between matched jets Gen(ak3) - Gen(ak7);#DeltaR;dN/d#DeltaR", 	200, 0., 2.);
          TH1D *hPhiDiff_7to3_gen = new TH1D("hPhiDiff_7to3_gen", 		"Distance in #varphi Gen(ak3) - Gen(ak7);#Delta#varphi;dN/d#Delta#varphi", 	200, 0., 3.15);
          TH1D *hEtaDiff_7to3_gen = new TH1D("hEtaDiff_7to3_gen", 		"Distance in #eta Gen(ak3) - Gen(ak7);#Delta#eta;dN/d#Delta#eta", 		200, 0., 0.8);

        TH2D *hJet_7to5_gen = new TH2D("hJet_7to5_gen", "Jets: energy of Gen(ak7) vs. Gen(ak5);E_{Gen(ak7)} (GeV);E_{Gen(ak5)}/E_{Gen(ak7)}", 50, 100., 3500., 20, 0., 1.);
          TH1D *hDistance_7to5_gen= new TH1D("hJetDistance_7to5_gen",      "Distance between matched jets Gen(ak5) - Gen(ak7);#DeltaR;dN/d#DeltaR",        200, 0., 2.);
          TH1D *hPhiDiff_7to5_gen = new TH1D("hPhiDiff_7to5_gen",          "Distance in #varphi Gen(ak5) - Gen(ak7);#Delta#varphi;dN/d#Delta#varphi",      200, 0., 3.15);
          TH1D *hEtaDiff_7to5_gen = new TH1D("hEtaDiff_7to5_gen",          "Distance in #eta Gen(ak5) - Gen(ak7);#Delta#eta;dN/d#Delta#eta",               200, 0., 0.8);

        TH1D *hNjet_3_det = new TH1D("hNjet_3_det", "Jets: Det (ak3)", 11, -0.5, 10.5);
        TH1D *hNjet_5_det = new TH1D("hNjet_5_det", "Jets: Det (ak5)", 11, -0.5, 10.5);
        TH1D *hNjet_7_det = new TH1D("hNjet_7_det", "Jets: Det (ak7)", 11, -0.5, 10.5);

        TH2D *hJet_7to3_det = new TH2D("hJet_7to3_det", "Jets: energy of Det(ak7) vs. Det(ak3);E_{Det(ak7)} (GeV);E_{Det(ak3)}/E_{Det(ak7)}", 50, 100., 3500., 20, 0., 1.);
          TH1D *hDistance_7to3_det= new TH1D("hJetDistance_7to3_det",       "Distance between matched jets Det(ak3) - Det(ak7);#DeltaR;dN/d#DeltaR",        200, 0., 2.);
          TH1D *hPhiDiff_7to3_det = new TH1D("hPhiDiff_7to3_det",           "Distance in #varphi Det(ak3) - Det(ak7);#Delta#varphi;dN/d#Delta#varphi",      200, 0., 3.15);
          TH1D *hEtaDiff_7to3_det = new TH1D("hEtaDiff_7to3_det",           "Distance in #eta Det(ak3) - Det(ak7);#Delta#eta;dN/d#Delta#eta",               200, 0., 0.8);

        TH2D *hJet_7to5_det = new TH2D("hJet_7to5_det", "Jets: energy of Det(ak7) vs. Det(ak5);E_{Det(ak7)} (GeV);E_{Det(ak5)}/E_{Det(ak7)}", 50, 100., 3500., 20, 0., 1.);
          TH1D *hDistance_7to5_det= new TH1D("hJetDistance_7to5_det",      "Distance between matched jets Det(ak5) - Det(ak7);#DeltaR;dN/d#DeltaR",        200, 0., 2.);
          TH1D *hPhiDiff_7to5_det = new TH1D("hPhiDiff_7to5_det",          "Distance in #varphi Det(ak5) - Det(ak7);#Delta#varphi;dN/d#Delta#varphi",      200, 0., 3.15);
          TH1D *hEtaDiff_7to5_det = new TH1D("hEtaDiff_7to5_det",          "Distance in #eta Det(ak5) - Det(ak7);#Delta#eta;dN/d#Delta#eta",               200, 0., 0.8);


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
	
	std::cout << "start looping over files" << std::endl;
	
	// start file loop
	while((fn = (TObjString*)next()) && counter_events < nEvents) { 
//      while((fn = (TObjString*)next()) ) {
		
		currentfile_.Clear();
		currentTFile_->Clear();
		
		currentfile_ = fn->GetString();
		
		std::cout << "opening file " << currentfile_ << " ... " << std::endl;
		
		TStopwatch *timer = new TStopwatch();
		TThread *fileopener;
		fileopener = new TThread("fileopener",(void(*) (void *))&OpenROOTFile,(RadiusAnalyzer*) this);
		TThread::SetCancelAsynchronous();
		TThread::SetCancelOn();
		fileopener->Run();
		
		bool fileopen = false;
		bool timeout = true;
		while (timer->RealTime() < 30) {
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
		
		std::cout << "We are working with Det " << DetJetRadius << " and Gen " << GenJetRadius << endl;

		// get tree from file
		TTree *tree = new TTree("CastorTree","");
		currentTFile_->GetObject("castortree/CastorTree",tree);
		
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
		std::vector<MyGenJet> *chargedGenJets = NULL;
		std::vector<MyTrackJet> *trackJets = NULL;
		  // Radius-energy study.
		std::vector<MyCastorJet> *CastorJets_ak3 = NULL;
                std::vector<MyCastorJet> *CastorJets_ak5 = NULL;
                std::vector<MyCastorJet> *CastorJets_ak7 = NULL;
                std::vector<MyGenJet> *genJets_ak3 = NULL;
                std::vector<MyGenJet> *genJets_ak5 = NULL;
                std::vector<MyGenJet> *genJets_ak7 = NULL;

		
		TBranch *b_evtid = tree->GetBranch("EvtId");
        TBranch *b_evtkin = NULL;
        if (!isData_) b_evtkin = tree->GetBranch("GenKin");
		TBranch *b_BeamSpot = tree->GetBranch("beamSpot");
		TBranch *b_HLTrig = tree->GetBranch("HLTrig");
		TBranch *b_L1Trig = tree->GetBranch("L1Trig");
		TBranch *b_vertices = tree->GetBranch("primaryVertex");
		TBranch *b_castorrechits = tree->GetBranch("castorRecHit");
		TBranch *b_castortowers = tree->GetBranch("castorTower");
		TBranch *b_castorjets = tree->GetBranch("castorJet");

		TBranch *b_castorjets_ak3 = tree->GetBranch("ak3castorJet");
                TBranch *b_castorjets_ak5 = tree->GetBranch("ak5castorJet");   
                TBranch *b_castorjets_ak7 = tree->GetBranch("ak7castorJet");   

		if( string_det_radius.Contains("ak3") ){ b_castorjets = tree->GetBranch("ak3castorJet"); b_castorjets_ak3 = b_castorjets;}
                if( string_det_radius.Contains("ak5") ){ b_castorjets = tree->GetBranch("ak5castorJet"); b_castorjets_ak5 = b_castorjets;}
                if( string_det_radius.Contains("ak7") ){ b_castorjets = tree->GetBranch("ak7castorJet"); b_castorjets_ak7 = b_castorjets;}
		
		TBranch *b_PFJets = tree->GetBranch("pfJet");
		TBranch *b_genParts = NULL;
		if (!isData_) b_genParts = tree->GetBranch("GenPart");
		TBranch *b_caloTowers = tree->GetBranch("caloTower");
		TBranch *b_genJets = NULL;
		TBranch *b_chargedGenJets = NULL;
//		if (!isData_) b_genJets = tree->GetBranch("GenJet");
                TBranch *b_genJets_ak3 = NULL;
                TBranch *b_genJets_ak5 = NULL;
                TBranch *b_genJets_ak7 = NULL;

                if(!isData_){
                  b_genJets_ak3 = tree->GetBranch("ak3GenJet");
                  b_genJets_ak5 = tree->GetBranch("ak5GenJet");
                  b_genJets_ak7 = tree->GetBranch("ak7GenJet");
                }
		
                if (!isData_ && string_gen_radius.Contains("ak3")){ b_genJets = tree->GetBranch("ak3GenJet"); b_genJets_ak3 = b_genJets;}
                if (!isData_ && string_gen_radius.Contains("ak5")){ b_genJets = tree->GetBranch("ak5GenJet"); b_genJets_ak5 = b_genJets;}
                if (!isData_ && string_gen_radius.Contains("ak7")){ b_genJets = tree->GetBranch("ak7GenJet"); b_genJets_ak7 = b_genJets;}
		if (!isData_) b_chargedGenJets = tree->GetBranch("ChargedGenJet");
		TBranch *b_trackJets = tree->GetBranch("trackJet");
	
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
                b_castorjets_ak3->SetAddress(&CastorJets_ak3);
                b_castorjets_ak5->SetAddress(&CastorJets_ak5);
		b_castorjets_ak7->SetAddress(&CastorJets_ak7);
		cout << "After" << endl;
		b_PFJets->SetAddress(&PFJets);
		if (!isData_) b_genParts->SetAddress(&genParts);
		b_caloTowers->SetAddress(&caloTowers);
		if (!isData_) b_genJets->SetAddress(&genJets);
                if (!isData_) b_genJets_ak3->SetAddress(&genJets_ak3);
                if (!isData_) b_genJets_ak5->SetAddress(&genJets_ak5);
                if (!isData_) b_genJets_ak7->SetAddress(&genJets_ak7);
		if (!isData_) b_chargedGenJets->SetAddress(&chargedGenJets);
		b_trackJets->SetAddress(&trackJets);
		
		int Nevents = tree->GetEntriesFast();
		std::cout << "file opened, events in this file = " << Nevents << std::endl;
		totalevents += Nevents;
		
		// start event loop
		for (int i=0; i<Nevents && i < nEvents;i++) {
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
			// Hadron level code
			/////////////////////////////////////////
			
			if (!isData_) {
				b_genParts->GetEntry(i);
				
				// calculate xi of the event
				
				// sort genParticles in y, from y_min to y_max
				std::vector<MyGenPart> myTempParticles;
				std::vector<MyGenPart> myRapiditySortedParticles;
				// copy only final stable particles with realistic Rapidity in tempvector
				for (unsigned int ipart=0;ipart<genParts->size();ipart++) {
					if ((*genParts)[ipart].status == 1) 
						myTempParticles.push_back((*genParts)[ipart]);
				}
				// do actual sorting
				while (myTempParticles.size() != 0) {
					double min_y = 100000;
					int min_y_pos = -1;
					for (unsigned int ipart = 0;ipart<myTempParticles.size();ipart++) {
						if (myTempParticles[ipart].Rapidity() < min_y) {
							min_y = myTempParticles[ipart].Rapidity();
							min_y_pos = ipart;
						}
					}
					myRapiditySortedParticles.push_back(myTempParticles[min_y_pos]);
					myTempParticles.erase(myTempParticles.begin()+min_y_pos);
				}
				
				// find deltaymax
				double deltaymax = 0;
				int deltaymax_pos = -1;
				for (unsigned int ipart=0;ipart<myRapiditySortedParticles.size()-1;ipart++) {
					double deltay = myRapiditySortedParticles[ipart+1].Rapidity() - myRapiditySortedParticles[ipart].Rapidity();
					if (deltay > deltaymax) {
						deltaymax = deltay;
						deltaymax_pos = ipart;
					}
				}
				ymax = deltaymax;
				
				// calculate Mx2 and My2
				long double XEtot = 0;
				long double XPxtot = 0;
				long double XPytot = 0;
				long double XPztot = 0;
				long double YEtot = 0;
				long double YPxtot = 0;
				long double YPytot = 0;
				long double YPztot = 0;
				
				for (int ipart=0;ipart<=deltaymax_pos;ipart++) {
					XEtot += myRapiditySortedParticles[ipart].E();
					XPxtot += myRapiditySortedParticles[ipart].Px();
					XPytot += myRapiditySortedParticles[ipart].Py();
					XPztot += myRapiditySortedParticles[ipart].Pz();
				}
				long double Mx2 = -1.;
				Mx2 = XEtot*XEtot - XPxtot*XPxtot - XPytot*XPytot - XPztot*XPztot;
				
				for (unsigned int ipart=deltaymax_pos+1;ipart<myRapiditySortedParticles.size();ipart++) {
					YEtot += myRapiditySortedParticles[ipart].E();
					YPxtot += myRapiditySortedParticles[ipart].Px();
					YPytot += myRapiditySortedParticles[ipart].Py();
					YPztot += myRapiditySortedParticles[ipart].Pz();
				}
				long double My2 = YEtot*YEtot - YPxtot*YPxtot - YPytot*YPytot - YPztot*YPztot;
				
				// calculate xix and xiy
				xix = Mx2/(7000*7000);
				xiy = My2/(7000*7000);
				
				// xi of event is max
				xi = std::max(xix,xiy);
				xidd=xix*xiy*7000*7000/(0.938*0.938);
                
                    
                // combine selection
        		// 7 TeV Xi cuts
				if (xix>0.04 || xiy>0.1 || xidd>0.5) passedHadronCuts = true;

				
				// if event passed hadron cut, execute analysis code
				if (passedHadronCuts) {
					
					
					
				} // end if passedHadronCuts
				
			} // end if not data
             
			
			
			/////////////////////////////////////////
			// Do stuff before filters
			/////////////////////////////////////////
						
			b_evtid->GetEntry(i);
			b_HLTrig->GetEntry(i);
			b_L1Trig->GetEntry(i);
			b_vertices->GetEntry(i);
			b_caloTowers->GetEntry(i);
			b_castorrechits->GetEntry(i);
			b_castortowers->GetEntry(i);
			
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
				
				// only fill the histograms when there's 1 vertex (filter out pile-up)
				if (Vertices->size() == 1) {
					
					// get event id stuff
					if( ((i+1) % 10000) == 0) cout << " run " << evtid->Run << " isData = " << evtid->IsData << " lumiblock " << 
						evtid->LumiBlock << " event " << evtid->Evt << endl; 
					if (!evtid->IsData) isMC = true;
					
					// calculate energy flow in CASTOR
					double CASTOReflow = 0;
					
					for (unsigned int j=0;j<CastorRecHits->size();j++) {
						MyCastorRecHit rechit = (*CastorRecHits)[j];
						//if (rechit.bad) std::cout << " found bad rechit: channel " << rechit.cha << std::endl;
						if (!rechit.bad) {
							hCASTOReflow_channel[rechit.cha-1]->Fill(rechit.energy);
							if (rechit.cha <= 80) CASTOReflow += rechit.energy;
							h2CASTOReflow_grid->Fill(rechit.sec,rechit.mod,rechit.energy);
						}
					}
				
					// tower multiplicity code
					hCASTORTowerMulti->Fill(CastorTowers->size());
					
					// fill all minbias eflow histos
					hCASTOReflow->Fill(CASTOReflow);
										
					// ------------------------
					// No central jet required.		

					vector<MyCastorJet> good_castorJets;
					vector<MyGenJet> good_genJets;

					// analyse castor jets
					int NCastorJets = 0;
					vector<double> det_casjet, gen_casjet_energy;
					b_castorjets->GetEntry(i);
					for (unsigned int j=0;j<CastorJets->size();j++) {
						MyCastorJet casjet = (*CastorJets)[j];
						if (casjet.energy > jetEThreshold) {
							NCastorJets++;
							hCastorJet_energy->Fill(casjet.energy);
								det_casjet.push_back(casjet.energy); 
								good_castorJets.push_back( casjet );

							hCastorJet_pt->Fill(casjet.energy*sin(2*atan(exp(5.9))));
							hCastorJet_em->Fill(casjet.eem);
							hCastorJet_had->Fill(casjet.ehad);
							hCastorJet_eta->Fill(casjet.eta);
							hCastorJet_phi->Fill(casjet.phi);
							hCastorJet_fem->Fill(casjet.fem);
							hCastorJet_fhot->Fill(casjet.fhot);
							hCastorJet_width->Fill(casjet.width);
							hCastorJet_depth->Fill(casjet.depth);
							hCastorJet_sigmaz->Fill(casjet.sigmaz);
							hCastorJet_ntower->Fill(casjet.ntower);
						}
					}
					hCastorJet_multi->Fill(NCastorJets);

					// -------------------------
					// AVS - no Central jet required					
					if(!isData_){
					  for(int jet = 0; jet < genJets->size(); jet++){
					    MyGenJet currentJet = (*genJets)[jet];
					    //if( counter_events%1000 == 0) { cout << "\tGenjet\t" << jet << "\twith energy\t" << currentJet.Energy() << "\tand pT\t" << currentJet.Pt() << endl; }
					    hGenJet_energy->Fill( currentJet.Energy() );
					    if( currentJet.Eta() < (-5.2 - GenJetContained) && currentJet.Eta() > (-6.6 + GenJetContained) && currentJet.Energy() > jetEThreshold){
					      gen_casjet_energy.push_back( currentJet.Energy() );
					      hCastorJet_energy_gen->Fill( currentJet.Energy() );
                                              hCastorJet_pt_gen->Fill( currentJet.Pt() );
					      good_genJets.push_back( currentJet );
					    }
					  }
					}

					// Sort Generator jets in Energy to match Detector jet sorting.
					/*
					bool unsorted = true;
					if( gen_casjet_energy.size() == 1 || gen_casjet_energy.size() == 0){ unsorted = false; }
					
					while (unsorted){
					   bool have_sorted = false;
					  for(int element = 0; element < gen_casjet_energy.size()-1; element++){

					    if( gen_casjet_energy[element] < gen_casjet_energy[element+1]){
					      double temp = gen_casjet_energy[element+1];
					      gen_casjet_energy[element+1] = gen_casjet_energy[element];
					      gen_casjet_energy[element] = temp;
					      have_sorted = true;

					      MyGenJet tempjet = good_genJets[element+1];
					      good_genJets[element+1] = good_genJets[element];
					      good_genJets[element] = tempjet;
					    }
					  }
					  if( !have_sorted ){ unsorted = false; };
					  
					}
					*/
/*
					if( counter_events%1 == 0) {
					  for(int element = 0; element < gen_casjet_energy.size(); element++){
					    cout << "Element\t" << element << " has energy " << gen_casjet_energy[element] << endl;
					  }
					}
*/
					/**********************
					 * ********************
					 * *******************/

					// Get the number of jets and the leading jet energy.
				        if( good_genJets.size() > 0) { hNjet_vs_Ejets_gen->Fill( good_genJets.size(), 	(good_genJets[ 0 ]).Energy() ); }
                                        else{ hNjet_vs_Ejets_gen->Fill(0.,0.); }
                                        if( good_castorJets.size() > 0) { hNjet_vs_Ejets_det->Fill( good_castorJets.size(), (good_castorJets[ 0 ]).energy ); }
					else{ hNjet_vs_Ejets_det->Fill(0., 0.); }

					// Match DET and GEN jets. 
					hNumber_of_match_jets->Fill( good_castorJets.size(), good_genJets.size());                             
					int matched_pairs = 0;



					hMatched->Fill( matched_pairs );
					



					/* Study of energy-radius. */
					b_castorjets_ak3->GetEntry(i);
                                          hNjet_3_det->Fill( CastorJets_ak3->size() );

                                        b_castorjets_ak5->GetEntry(i);
                                          hNjet_5_det->Fill( CastorJets_ak5->size() );

                                        b_castorjets_ak7->GetEntry(i);
                                          hNjet_7_det->Fill( CastorJets_ak7->size() );		
	
					for(int detjet_ak7 = 0; detjet_ak7 < CastorJets_ak7->size(); detjet_ak7++){
					  double detjet_ak7_energy = ( (*CastorJets_ak7)[detjet_ak7] ).energy;
					
                                          if( detjet_ak7_energy > jetEThreshold ){
//					    cout << "\tak7 jet\t" << detjet_ak7 << "\tabove threshold" << endl;
                                            double detjet_ak7_phi = ( (*CastorJets_ak7)[detjet_ak7] ).phi;
                                            double detjet_ak7_eta = ( (*CastorJets_ak7)[detjet_ak7] ).eta;					  

					    // Look at ak3 jets. //
					    double distance_7to3 = 1.;
					    double curr_distance = 1.;
					    int ideal_3 = -1;
					    double ideal_phi3 = 0.;
					    double ideal_eta3 = 0.;
					    double ideal_energy3 = 0.;

					    for(int detjet_ak3 = 0; detjet_ak3 < CastorJets_ak3->size(); detjet_ak3++){
					      double detjet_ak3_energy = ( (*CastorJets_ak3)[detjet_ak3] ).energy;

					      if(  detjet_ak3_energy > jetEThreshold ){
//	                                            cout << "\t\tak3 jet\t" << detjet_ak3 << "\tabove threshold" << endl;
					        double detjet_ak3_phi = ( (*CastorJets_ak3)[detjet_ak3] ).phi;
                                                double detjet_ak3_eta = ( (*CastorJets_ak3)[detjet_ak3] ).eta;

						curr_distance = pow( detjet_ak3_phi-detjet_ak7_phi, 2.) + pow(detjet_ak3_eta-detjet_ak7_eta, 2.); // Do not take sqrt yet to save CPU time.
						if( curr_distance < distance_7to3 ){
						  distance_7to3 = curr_distance;	
						  ideal_phi3 = detjet_ak3_phi;
						  ideal_eta3 = detjet_ak3_eta;
						  ideal_energy3 = detjet_ak3_energy;
						}
					      }
					      else{ detjet_ak3 = CastorJets_ak3->size() + 1; }
					    } // Loop over ak3.

				            hDistance_7to3_det->Fill( sqrt(distance_7to3) );	// Now take sqrt on distance.
					    hPhiDiff_7to3_det->Fill( fabs(ideal_phi3 - detjet_ak7_phi) );
                                            hEtaDiff_7to3_det->Fill( fabs(ideal_eta3 - detjet_ak7_eta) );
					    
					    hJet_7to3_det->Fill( detjet_ak7_energy, ideal_energy3/detjet_ak7_energy);
//					    cout << "Fill energy\t" << detjet_ak7_energy << ideal_energy3 << endl;
					  

					  // Look at ak5 jets.
                                            double distance_7to5 = 1.;
                                            curr_distance = 1.;
                                            int ideal_5 = -1;
                                            double ideal_phi5 = 0.;
                                            double ideal_eta5 = 0.;
                                            double ideal_energy5 = 0.;

                                            for(int detjet_ak5 = 0; detjet_ak5 < CastorJets_ak5->size(); detjet_ak5++){
                                              double detjet_ak5_energy = ( (*CastorJets_ak5)[detjet_ak5] ).energy;

                                              if(  detjet_ak5_energy > jetEThreshold ){
                                                double detjet_ak5_phi = ( (*CastorJets_ak5)[detjet_ak5] ).phi;
                                                double detjet_ak5_eta = ( (*CastorJets_ak5)[detjet_ak5] ).eta;

                                                curr_distance = pow( detjet_ak5_phi-detjet_ak7_phi, 2.) + pow(detjet_ak5_eta-detjet_ak7_eta, 2.); // Do not take sqrt yet to save CPU time.
                                                if( curr_distance < distance_7to5 ){
                                                  distance_7to5 = curr_distance;
                                                  ideal_phi5 = detjet_ak5_phi;
                                                  ideal_eta5 = detjet_ak5_eta;
                                                  ideal_energy5 = detjet_ak5_energy;
                                                }
                                              }
                                              else{ detjet_ak5 = CastorJets_ak5->size() + 1; }
					    } // loop over ak5.

                                            hDistance_7to5_det->Fill( sqrt(distance_7to5) );        // Now take sqrt on distance.
                                            hPhiDiff_7to5_det->Fill( fabs(ideal_phi5 - detjet_ak7_phi) );
                                            hEtaDiff_7to5_det->Fill( fabs(ideal_eta5 - detjet_ak7_eta) );

                                            hJet_7to5_det->Fill( detjet_ak7_energy, ideal_energy5/detjet_ak7_energy);
					  
					  }
					  else{ detjet_ak7 = CastorJets_ak7->size() + 1; }
					} // Loop over Det ak7 jets.


//					cout << endl << endl;
						
					// end of event, print status
					if( ((i+1) % 10000) == 0) std::cout << i+1 <<"events done in file " << it << std::endl;
					totalevents++;
					
				} // end if statement for 1 vertex
								
			} // end if passed detector cuts
			
			// combined hadron and detector level code
			
			if (passedHadronCuts && passedDetectorCuts) {
			}
			
		} // end event loop
		
		delete tree;
		currentTFile_->Close();
		it++;
	} // end file loop
	
	std::cout << "file loop has ended" << std::endl;
    
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
    

	// Create proper response matrices including fakes and misses.
	for(int bins_x = 0; bins_x < hCastorJet_energy_response->GetNbinsX(); bins_x++){

	  for(int bins_y = 0; bins_y < hCastorJet_energy_response->GetNbinsY(); bins_y++){
	    double bin_value = hCastorJet_energy_response->GetBinContent( bins_x, bins_y );
	
	    hCastorJet_Matrix->SetBinContent( bins_x, bins_y, bin_value );

	  }
	}

	for(int bins_x = 0; bins_x < hCastorJet_energy_fakes->GetNbinsX(); bins_x++){
	  double miss = hCastorJet_energy_misses->GetBinContent( bins_x );
	  hCastorJet_Matrix->SetBinContent( 1, bins_x, miss);

          double fake = hCastorJet_energy_fakes->GetBinContent( bins_x );
cout << bins_x << "\t" << hCastorJet_Matrix->GetBinCenter( bins_x) << "\tFakes\t" << fake << "\t" << endl;
          hCastorJet_Matrix->SetBinContent( bins_x, 1, fake);
	}

	/*
	// Do the same for central forward stuff.
        for(int bins_x = 0; bins_x < hCastorJet_cf_energy_response->GetNbinsX(); bins_x++){

          for(int bins_y = 0; bins_y < hCastorJet_cf_energy_response->GetNbinsY(); bins_y++){
            double bin_value = hCastorJet_cf_energy_response->GetBinContent( bins_x, bins_y );

            hCastorJet_cf_Matrix->SetBinContent( bins_x, bins_y, bin_value );
          }
        }

        for(int bins_x = 0; bins_x < hCastorJet_cf_energy_fakes->GetNbinsX(); bins_x++){
          double miss = hCastorJet_cf_energy_misses->GetBinContent( bins_x );
          hCastorJet_cf_Matrix->SetBinContent( 10, bins_x, miss);

          double fake = hCastorJet_cf_energy_fakes->GetBinContent( bins_x );
          hCastorJet_cf_Matrix->SetBinContent( bins_x, 10, fake);
        }
	*/
    
    // write all histo's to file
    
	std::cout << "total number of events = " << totalevents << " from " << it << " file(s)" << endl;
	
	
	// create output root file
	Char_t filename[200];
	const char* part = currentfile_.Data();
	std::string first(outputname_);
         
	first += part;
	//strcat(temp,part);
//	sprintf(filename,"output_RadiusAnalyzer_%s",first.c_str());
        if( !isData_) { sprintf(filename,"20141103_Output_RadiusAnalyzer_GEN_" + string_gen_radius + "_DET_" + string_det_radius + "_unsortedE_%i_%s", counter_events, first.c_str()); }
        else{ sprintf(filename,"20141103_Output_RadiusAnalyzer_Data_" + string_det_radius + "_unsortedE_%i_%s", counter_events, first.c_str()); }
	TFile* output = new TFile(filename,"RECREATE");
	output->cd();
		
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

	hJRE->Divide( hMatched, hUnmatched );
	
	hCastorJet_energy->Write();
	hCastorJet_pt->Write();
	hCastorJet_em->Write();
	hCastorJet_had->Write();
	hCastorJet_fem->Write();
	hCastorJet_fhot->Write();
	hCastorJet_width->Write();
	hCastorJet_depth->Write();
	hCastorJet_sigmaz->Write();
	hCastorJet_ntower->Write();
	hCastorJet_eta->Write();
	hCastorJet_phi->Write();
	hCastorJet_multi->Write();
	hGenJet_energy->Write();

        hCastorJet_energy_gen->Write();
	hCastorJet_pt_gen->Write();

        hCastorJet_energy_response->Write();
        hCastorJet_energy_fakes->Write();
        hCastorJet_energy_misses->Write();
	hCastorJet_Matrix->Write();

	hCastorJet_energy_ratio->Write();

        hCastorJet_cf_energy->Write();
        hCastorJet_cf_energy_gen->Write();
        hCastorJet_cf_energy_response->Write();
        hCastorJet_cf_energy_fakes->Write();
        hCastorJet_cf_energy_misses->Write();
        hCastorJet_cf_Matrix->Write();

	hDistance->Write();
	hPhiDiff->Write();
        hEtaDiff->Write();
	hEtaPhiDiff->Write();
        hEtaRDiff->Write();
        hPhiRDiff->Write();
	hNumber_of_match_jets->Write();

	hTrackjets_2D_number->Write();
	hTrackjets_2D_pt->Write();

	hJER->Write();
	hJER_per_energy->Write();
        hJER_per_distance->Write();
        hJER_2->Write();
        hJER_per_energy_2->Write();
        hJER_per_distance_2->Write();

	hJER_per_eta->Write();
	hEnergy_per_eta->Write();

	hMatched->Write();
	hUnmatched->Write();
	hJRE->Write();		
*/
	hNjet_vs_Ejets_gen->Write();
        hNjet_vs_Ejets_det->Write();

	// Study of jet radius-energy
	hNjet_3_gen->Write();
	hNjet_5_gen->Write();
	hNjet_7_gen->Write();
	
	hJet_7to3_gen->Write();
	hDistance_7to3_gen->Write();
	hPhiDiff_7to3_gen->Write(); 
	hEtaDiff_7to3_gen->Write(); 

	hJet_7to5_gen->Write(); 
	hDistance_7to5_gen->Write();
	hPhiDiff_7to5_gen->Write(); 
	hEtaDiff_7to5_gen->Write(); 

	hNjet_3_det->Write(); 
	hNjet_5_det->Write(); 
	hNjet_7_det->Write(); 

	hJet_7to3_det->Write(); 
	hDistance_7to3_det->Write();
	hPhiDiff_7to3_det->Write(); 
	hEtaDiff_7to3_det->Write(); 

	hJet_7to5_det->Write(); 
	hDistance_7to5_det->Write();
	hPhiDiff_7to5_det->Write();
	hEtaDiff_7to5_det->Write(); 


	response.Write();
	response_onlyMatches.Write();
        response_cf.Write();
	response_cf_onlyMatches.Write();

	output->Close();
	std::cout << "file " << filename << " created." << std::endl;
 LoopOutputFile_ = filename;
        cout << totalevents << "\tevents" << endl;
	cout << counter_jer << "\tJER events" << endl;

}
	
void RadiusAnalyzer::AfterLoopCalculations(TString file) {

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
        
  	hCASTOReflow_channels->Write();
		
	output->Close();
	std::cout << "file " << filename << " created." << std::endl;
	
}

TString RadiusAnalyzer::getOutputFile() {
    return LoopOutputFile_;
}

TString RadiusAnalyzer::getInputDir() {
	return inputdir_;
}

TString RadiusAnalyzer::getCurrentFile() {
	return currentfile_;
}

void RadiusAnalyzer::setCurrentTFile() {
	currentTFile_ = currentStaticTFile_;
}

void* RadiusAnalyzer::OpenROOTFile(RadiusAnalyzer* arg) {
	currentStaticTFile_ = TFile::Open(arg->getInputDir()+arg->getCurrentFile(),"READ");
	return 0;
}

int RadiusAnalyzer::posLeadingGenJet(std::vector<MyGenJet> JetVector,double etacut,double minptcut) {
	
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

int RadiusAnalyzer::posLeadingTrackJet(std::vector<MyTrackJet> JetVector,double etacut,double minptcut) {
	
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

