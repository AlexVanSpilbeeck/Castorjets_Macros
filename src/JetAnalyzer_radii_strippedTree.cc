////////////////////////////////////////
//////// New CMSSW_4_2_X version ///////
////////////////////////////////////////

/* Updated by Alex Van Spilbeeck.

 This code looks at:
	-> Gen Jets (-6.1 < eta < -5.7), E above Ethreshold, pT > 1 GeV
	-> Castor Jets, E > 0 GeV, pT > 1 GeV
	-> Gen jets are matched to closest detector level jet in \Delta phi. 
*/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "JetAnalyzer_radii_strippedTree.h"
#include "HistoRetriever.h"



#include "IsolationCut.cpp"
#include "CalibCorrection.h"
#include "CastorSector.h"
#include "PlaceThreshold.h"
#include "../Function_JetType.h"

//STANDARD ROOT INCLUDES

#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
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
#define jetEThreshold_gen 100. 
#define jetEThreshold_det 100. 
#define jetEThreshold 0.
#define EbinWidth 5.
#define EbinWidth_rel 1.4
#define PI 3.14159265359

//#define cut_EI false

#define comments_ false
#define manualUnfold_ true
#define do_calibration_discrete false
//#define do_calibration_function true
#define sector_dependence true

//#define prepare_unfolding true

#include "../../../RooUnfold-1.1.1/src/RooUnfold.h"
//#include "../../../RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#endif





TFile *JetAnalyzer_radii_strippedTree::currentStaticTFile_ = new TFile();

JetAnalyzer_radii_strippedTree::JetAnalyzer_radii_strippedTree(TString inputdir, bool isData, const char* outputname, int totalEvents, TString date, TString filename, TString jettype, double threshold, TString setup, double deltaphimax, double etawidth, TString match) {

    
	std::cout << "constructing JetAnalyzer_radii_strippedTree class..." << std::endl;
	
    inputdir_ = inputdir;
    filename_ = filename;
    isData_ = isData;
    outputname_ = outputname;
    totalEvents_ = totalEvents;
    date_ = date;	
    jettype_ = jettype;
    sectors_ = 0;
    threshold_ = threshold;
    setup_ = setup;
    deltaphimax_ = deltaphimax;
    genetamin_ = -6.6 - etawidth;
    genetamax_ = -5.2 + etawidth;
    match_ = match;

    prepare_unfolding = false;
  
    etaband_ = 0.4;
    if( setup_ == "raw_wide" || setup_ == "unfold" || isData_) etaband_ = 1.4;

    cut_EI = false;
    if( setup_ == "isolated" || setup_ == "calibrated") cut_EI = true;

    do_calibration_function = false;
    if( setup_ == "calibrated" || setup_ == "unfold") do_calibration_function = true;

    if( setup_ == "unfold") prepare_unfolding = true;

    std::cout << "initializing basic variables..." << std::endl;
    
    // initialize basic variables
    LoopOutputFile_ = "";

    if( date == "0" ){       
      LoopOutputFile_ += filename;
      LoopOutputFile_.ReplaceAll( "Stripped_trees/", "");
      LoopOutputFile_.ReplaceAll( "STRIPPED_TREE.root", "");
      LoopOutputFile_ += setup;    
      LoopOutputFile_ += TString::Format("_Emin_%f", threshold_);
      LoopOutputFile_ += TString::Format("_deltaPhiMax_%f", deltaphimax_ );
      LoopOutputFile_ += TString::Format("_etaband_%f", etawidth);
      LoopOutputFile_ += "_" + match_;
      LoopOutputFile_ += ".root"; 
    }
    cout << "Output name\t" << LoopOutputFile_ << "\tfrom\t" << date << endl;
    currentfile_ = "";
    currentTFile_ = new TFile();
	
    std::cout << "all initialisations done, class constructed!" << std::endl;


    if( LoopOutputFile_.Contains( "data" ) ){ fileLabel_ = "data"; }
    else if( LoopOutputFile_.Contains( "displaced_down" )) { fileLabel_ = "displaced_down"; }
    else if( LoopOutputFile_.Contains( "displaced_up" )){ fileLabel_ = "displaced_up"; }
    else if( LoopOutputFile_.Contains( "displaced" )){ fileLabel_ = "displaced"; }
    else if( LoopOutputFile_.Contains( "FullSimulation" )){ fileLabel_ = "FullSimulation"; }
    else if( LoopOutputFile_.Contains( "Pythia84C" )){ fileLabel_ = "Pythia84C"; }
    else{ fileLabel_ = "Pythia6Z2star"; }

    cout << "@@@ Filelabel is\t" << fileLabel_ << endl;
}

JetAnalyzer_radii_strippedTree::~JetAnalyzer_radii_strippedTree() { }

void JetAnalyzer_radii_strippedTree::Loop() {

#ifdef __CINT__
  gSystem->Load("../../../RooUnfold-1.1.1/libRooUnfold.so");
#endif
	
	
	std::cout << " JetAnalyzer_radii_strippedTree Loop function is started " << std::endl;
	
//	TString tstring = outputname_;
//	std::cout << " TString outputname_ = " << tstring << std::endl;


 
/*	
    // reweight the MC in this case
    bool reweightMC = false;
    if (tstring.Contains("Reweighted")) {
        std::cout << "We will reweight the MC now !!!" << std::endl;
        reweightMC = true;
    }
    
*/	//
	// Unfolding can only happen if we have MC.
	//

	bool prepare_unfolding_ = prepare_unfolding;
	if( isData_ ){ prepare_unfolding_ = false; }
	
	using namespace std;
	int it = 0;
	int totalevents = 0;
	
	/////////////////////////////////////
	// Define all histograms
	/////////////////////////////////////


	int Ebins_fix = 100;
	double Emin_fix = -100.;
	double Emax_fix = 3000.;
	
        // AVS - number of bins.
	/* We get binwidth from peak JER. */
//	double Emin = jetEThreshold_det + 0., Emax = 3000.;
//	int Ebins = static_cast<int> ( (Emax - Emin)/EbinWidth );
//	double Emin = 0., Emax = 3000.;
//	int Ebins = 50;
        double Emin, Emax = 2100.;
        int Ebins;

        if( threshold_ == 0. ){ Emin = threshold_; Ebins = 28; }
	if( threshold_ == 150.){ Emin = threshold_; Ebins = 26;}
	if( threshold_ == 300.){ Emin = threshold_; Ebins = 24;}

         

        if( comments_ ){ cout << Ebins << " bins from " << Emin << " GeV to " << Emax << " GeV" << endl; }

	/* If we want variable Ebins. */
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

	  current_lowE += current_lowE*0.2;
	  count_bincenters++;
	}

//	Ebins = binEdges.size();
        float *Ebins_var = new float[Ebins];
//float Ebins_var[Ebins];
	for(int i = 0; i < binEdges.size(); i++){
	  Ebins_var[i] = binEdges[i];
	}

/*
	// We need an extended bin range for our matrix + hits & misses
	// Let us take bins 0 - 50, 50 - 100
	float *Ebins_ext = new float[Ebins+2];
	Ebins_ext[0] = 0.;
	Ebins_ext[1] = 50.;
	for(int i = 0; i < binEdges.size(); i++){
	  Ebins_ext[i+2] = binEdges[i];
	}
        cout << "Variable bins done" << endl;
*/
        /* End variable ebins. */


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
        TH1D *hCastorJet_energy = new TH1D("hCastorJet_energy","CastorJet energy distribution",	Ebins, Emin, Emax);		hCastorJet_energy->Sumw2();
        TH1D *hCastorJet_energy_lead = new TH1D("hCastorJet_energy_lead","CastorJet energy distribution", Ebins, Emin, Emax);	hCastorJet_energy_lead->Sumw2();
        TH1D *hCastorJet_energy_fine = new TH1D("hCastorJet_energy_fine","CastorJet energy distribution", Ebins*10, Emin, Emax);hCastorJet_energy_fine->Sumw2();


        TH1D *hGenJet_energy = new TH1D("hGenJet_energy","GenJet energy distribution",		Ebins, Emin, Emax);
	TH1D *hGenJet_energy_lead = new TH1D("hGenJet_energy_lead","Leading GenJet energy distribution",          Ebins, Emin, Emax);
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

		// All Casto jets response matrix.
		TH2D *hCastorJet_energy_response = new TH2D("hCastorJet_energy_response","CastorJet energy distribution - Response;E_{det};E_{gen};", 	Ebins, Emin, Emax, Ebins, Emin, Emax);	// ,Ebins+10,Emin-10*EbinWidth,Emax,Ebins+10,Emin-10*EbinWidth,Emax);
		TH2D *hCastorJet_energy_response_lead = new TH2D("hCastorJet_energy_response_lead","CastorJet energy distribution - Response;E_{det};E_{gen};", 	Ebins, Emin, Emax, Ebins, Emin, Emax);	// ,Ebins+10,Emin-10*EbinWidth,Emax,Ebin
		TH2D *hCastorJet_energy_response_fine = new TH2D("hCastorJet_energy_response_fine","CastorJet energy distribution - Response;E_{det};E_{gen};", 	Ebins*10, Emin, Emax, Ebins*10, Emin, Emax);

		TH2D *hCastorJet_energy_fakes 	 = new TH2D("hCastorJet_energy_fakes",	 "CastorJet energy distribution - Fakes", 	Ebins, Emin, Emax, Ebins, Emin, Emax);	// ,Ebins+10,Emin-10*EbinWidth,Emax);
		TH2D *hCastorJet_energy_misses	 = new TH2D("hCastorJet_energy_misses",	 "CastorJet energy distribution - Misses", 	Ebins, Emin, Emax, Ebins, Emin, Emax);	// ,Ebins+10,Emin-10*EbinWidth,Emax);

		TH1D *hCastorJet_miss_all   = new TH1D("hCastorJet_miss_all", "Misses", Ebins, Emin, Emax);
                TH1D *hCastorJet_fake_all   = new TH1D("hCastorJet_fake_all", "Fakes", Ebins, Emin, Emax);
                TH1D *hCastorJet_miss_lead   = new TH1D("hCastorJet_miss_lead", "Misses", Ebins, Emin, Emax);
                TH1D *hCastorJet_fake_lead   = new TH1D("hCastorJet_fake_lead", "Fakes", Ebins, Emin, Emax);

		TH1D *hCastorJet_matchedGen_all   = new TH1D("hCastorJet_matchedGen_all", "MatchedGen", Ebins, Emin, Emax);
		
		TH2D *hCastorJet_energy_complete_response = new TH2D("hCastorJet_energy_complete_response","CastorJet energy distribution - Complete Response", Ebins_fix, Emin_fix, Emax_fix, Ebins_fix, Emin_fix, Emax_fix);

		int Bins_response_matrix[7]	= {Ebins*10, 	Ebins*10, 16,  16,  32,  16,  32};
		double Min_response_matrix[7]	= {Emin, 	Emin, 	  0.,  0.,  0.,  0.,  0.};
		double Max_response_matrix[7]	= {Emax, 	Emax, 	  16., 16., 32., 16., 32.};

		THnSparseD *hResponse_leading 	= new THnSparseD("hResponse_leading", 	"Multidimensional Response Matrix;E_{det};E_{gen};N_{match};N_{fake};N_{miss};N_{det};N_{gen};",	7, Bins_response_matrix, Min_response_matrix, Max_response_matrix); hResponse_leading->Sumw2();
//		THnSparseD *hResponse_inclusive = new THnSparseD("hResponse_inclusive", "E_{det};E_{gen};N_{fake};N_{miss}",    4, Bins_response_matrix, Min_response_matrix, Max_response_matrix); hResponse_inclusive->Sumw2();		

		// Response with variable binning.
		TH2D *hCastorJet_energy_ratio = new TH2D("hCastorJet_energy_ratio", "Ratio of generator to detector", 1000, 0., 10.,Ebins, Emin, Emax);

		TH2D *hCastorJet_Matrix = new TH2D("hCastorJet_responseMatrix", "Castor jet Response Matrix",Ebins-1, Ebins_var,Ebins-1, Ebins_var);

		  TH2D *hCastorJet_Matrix_had_det = new TH2D("hCastorJet_responseMatrix_had_det", 	"Castor jet Response Matrix (had. det. jet)"	,Ebins, Emin, Emax, Ebins, Emin, Emax);
		  TH2D *hCastorJet_Matrix_em_det = new TH2D("hCastorJet_responseMatrix_em_det", 	"Castor jet Response Matrix (em. det. jets)"	,Ebins, Emin, Emax, Ebins, Emin, Emax);
		  
		  TH2D *hCastorJet_Matrix_had_det_1sector = new TH2D("hCastorJet_responseMatrix_had_det_1sector", 	"Castor jet Response Matrix (had. det. jet)"	,Ebins, Emin, Emax, Ebins, Emin, Emax);
		  TH2D *hCastorJet_Matrix_em_det_1sector = new TH2D("hCastorJet_responseMatrix_em_det_1sector", 	"Castor jet Response Matrix (em. det. jets)"	,Ebins, Emin, Emax, Ebins, Emin, Emax);

		  TH2D *hCastorJet_Matrix_had_det_nsector = new TH2D("hCastorJet_responseMatrix_had_det_nsector", 	"Castor jet Response Matrix (had. det. jet)"	,Ebins, Emin, Emax, Ebins, Emin, Emax);
		  TH2D *hCastorJet_Matrix_em_det_nsector = new TH2D("hCastorJet_responseMatrix_em_det_nsector", 	"Castor jet Response Matrix (em. det. jets)"	,Ebins, Emin, Emax, Ebins, Emin, Emax);
			  
		TH1D *hCastorJet_energy_gen = new TH1D("hCastorJet_energy_gen","CastorJet energy distribution",Ebins-1, Ebins_var);
		TH1D *hCastorJet_pt_gen = new TH1D("hCastorJet_pt_gen","CastorJet pt distribution",500,0,10.);

		cout << "pre HJER" << endl;

		// Efficiency control.
		TH1D *hJER = new TH1D("hJER", 	"Jet Energy Resolution", 200, -1, 5);
		TH1D *hJER_had = new TH1D("hJER_had", 	"Jet Energy Resolution", 200, -1, 5);
		TH1D *hJER_em = new TH1D("hJER_em", 	"Jet Energy Resolution", 200, -1, 5);

		cout << "hJER" << endl;


		// Delta E/E_gen vs. E_gen
		TH2D *hJER_per_energy = new TH2D("hJER_per_energy", "#DeltaE/E for fixed energies;E_{gen};JER", 			Ebins, Emin, Emax, 200, -5, 5);
		     TH2D *hJER_per_energy_had_det= new TH2D("hJER_per_energy_had_det", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);
		     TH2D *hJER_per_energy_em_det= new TH2D("hJER_per_energy_em_det", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);
		     
		     // 1 Castor sector.
		     TH2D *hJER_per_energy_had_det_1sector = new TH2D("hJER_per_energy_had_det_1sector", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);
		     TH2D *hJER_per_energy_em_det_1sector = new TH2D("hJER_per_energy_em_det_1sector", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);	

		     // 2+ Castor sector
		     TH2D *hJER_per_energy_had_det_nsector = new TH2D("hJER_per_energy_had_det_nsector", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);
		     TH2D *hJER_per_energy_em_det_nsector = new TH2D("hJER_per_energy_em_det_nsector", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);		          
		// Delta E/E_det vs. E_det
		TH2D *hJER_per_eDet = new TH2D("hJER_per_eDet", "#DeltaE/E for fixed energies;E_{det};#DeltaE/E",     Ebins, Emin, Emax, 200, -5, 5);
		     TH2D *hJER_per_eDet_had_det= new TH2D("hJER_per_eDet_had_det", "#DeltaE/E for fixed energies;E_{det};#DeltaE/E",   Ebins, Emin, Emax, 200, -5, 5);
		     TH2D *hJER_per_eDet_em_det= new TH2D("hJER_per_eDet_em_det", "#DeltaE/E for fixed energies;E_{det};#DeltaE/E",  Ebins, Emin, Emax, 200, -5, 5);

		     // 1 Castor sector
		     TH2D *hJER_per_eDet_had_det_1sector = new TH2D("hJER_per_eDet_had_det_1sector", "#DeltaE/E for fixed energies;E_{det};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);
		     TH2D *hJER_per_eDet_em_det_1sector = new TH2D("hJER_per_eDet_em_det_1sector", "#DeltaE/E for fixed energies;E_{det};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);	

		     // 2+ Castor sector
		     TH2D *hJER_per_eDet_had_det_nsector = new TH2D("hJER_per_eDet_had_det_nsector", "#DeltaE/E for fixed energies;E_{det};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);
		     TH2D *hJER_per_eDet_em_det_nsector = new TH2D("hJER_per_eDet_em_det_nsector", "#DeltaE/E for fixed energies;E_{det};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);	

		// Delta E/E_det vs. E_gen
		TH2D *hJER_per_eGen = new TH2D("hJER_per_eGen", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",     Ebins, Emin, Emax, 200, -5, 5);
		     TH2D *hJER_per_eGen_had_det= new TH2D("hJER_per_eGen_had_det", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",   Ebins, Emin, Emax, 200, -5, 5);
		     TH2D *hJER_per_eGen_em_det= new TH2D("hJER_per_eGen_em_det", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",  Ebins, Emin, Emax, 200, -5, 5);

		     // 1 Castor sector
		     TH2D *hJER_per_eGen_had_det_1sector = new TH2D("hJER_per_eGen_had_det_1sector", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",  Ebins, Emin, Emax, 200, -5, 5);
		     TH2D *hJER_per_eGen_em_det_1sector = new TH2D("hJER_per_eGen_em_det_1sector", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);	

		     // 2+ Castor sector
		     TH2D *hJER_per_eGen_had_det_nsector = new TH2D("hJER_per_eGen_had_det_nsector", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E", 	Ebins, Emin, Emax, 200, -5, 5);
		     TH2D *hJER_per_eGen_em_det_nsector = new TH2D("hJER_per_eGen_em_det_nsector", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);	

		
	        cout << "hJER 2D" << endl;

		TH2D *hJER_per_distance = new TH2D("hJER_per_distance", "#DeltaE/E for distance;#DeltaR;#DeltaE/E", 	200, 0., 6.5, 
																200, -5, 5);
		TH2D *hJER_per_eta = new TH2D("hJER_per_eta", "#DeltaE/E for distance;#DeltaR;#DeltaE/E",               14,-6.6,-5.2,
																200, -5, 5);
		TH2D *hEnergy_per_eta = new TH2D("hEnergy_per_eta", "E_{gen} vs. #eta of leading jet", Ebins, Emin, Emax,  14,-6.6,-5.2);

	        cout << "Matching histograms" << endl;

		// Matches.
		TH1D *hMatched = new TH1D("hMatched", "Castor and Gen jets match", 10, -0.5, 9.5);
		TH1D *hUnmatched = new TH1D("hUnmatched", "Castor and Gen jets don't match", Ebins, Emin, Emax);
		TH1D *hJRE = new TH1D("hJRE", "Jet Reconstruction Efficiency", Ebins, Emin, Emax);

		// Towers vs. phi
		TH2D *hTower_phi = new TH2D("hTower_phi", "N_{towers} vs. #Delta#varphi;N_{towers};#Delta#varphi", 6,-0.5,5.5, 50,-3.15,3.15);

		TH2D *hIsolationEnergy = new TH2D("hIsolationEnergy", "Energy in vicinity;Towers;E_{gen};", 5, -0.5, 4.5, 100, -1000., 1800.);

		TH1D *hIsolationEnergy_1D = new TH1D("hIsolationEnergy_1D", "Isolation energy;E_{s}", 10000, 0., 1800.);
		TH1D *hIsolationEnergy_Egen = new TH1D("hIsolationEnergy_Egen", "Isolation energy;E_{s}/E_{gen}", 10000, 0., 100.);
		TH1D *hIsolationEnergy_Edet = new TH1D("hIsolationEnergy_Edet", "Isolation energy;E_{s}/E_{det}", 10000, 0., 100.);

		// RooUnfold.
		cout << "RooUnfold" << endl;

		RooUnfoldResponse response 		(hCastorJet_energy, hCastorJet_energy, "response");
		RooUnfoldResponse response_lead		(hCastorJet_energy, hCastorJet_energy, "response_lead");

		RooUnfoldResponse response_fine 	(hCastorJet_energy_fine, hCastorJet_energy_fine, "response_fine");
		RooUnfoldResponse response_lead_fine	(hCastorJet_energy_fine, hCastorJet_energy_fine, "response_lead_fine");

		RooUnfoldResponse response_match	(hCastorJet_energy, hCastorJet_energy, "response_match");

		// Jet distance distribution.
		TH1D *hDistance 	= new TH1D("hJetDistance",	"Distance between matched jets;#DeltaR;dN/d#DeltaR", 200, 0., 2.);
		TH1D *hPhiDiff 		= new TH1D("hPhiDiff", 		"Distance in #varphi;#Delta#varphi;dN/d#Delta#varphi", 200, 0., 3.15);
		TH1D *hEtaDiff 		= new TH1D("hEtaDiff", 		"Distance in #eta;#Delta#eta;dN/d#Delta#eta", 200, -0.8, 0.8);
		TH2D *hEtaPhiDiff 	= new TH2D("hEtaPhiDiff", 	"Distance in #eta and #varphi;#eta;#varphi",200,0.,0.8,200,0.,3.15);
		TH2D *hEtaRDiff 	= new TH2D("hEtaRDiff", 	"Distance in #eta and R;#Delta#eta;#DeltaR",200,0.,0.8,200,0.,3.3);
		TH2D *hPhiRDiff 	= new TH2D("hPhiRDiff", 	"Distance in #eta and #varphi;#Delta#varphi;#DeltaR",200,0.,3.3,200,0.,3.3);

		// Pion to electron ratio.
		TH1D *hElectron_energy	= new TH1D("hElectron_energy",	"Energy of electron jets;E_{e} (GeV)",	20, 0., 1000.);
		TH1D *hPion_energy  	= new TH1D("hPion_energy",  	"Energy of Pion jets;E_{#pi} (GeV)",   	20, 0., 1000.);
		TH1D *hPi_e_ratio	= new TH1D("hPi_e_ratio",	"Ratio of pions to electrons;E (GeV)",	20, 0., 1000.);
		
		// Study the jets' energy versus the number of jets.
		TH2D *hNjet_vs_Ejets_gen = new TH2D("hNjet_vs_Ejets_gen", "Number of jets versus E leading jet (gen)",        10, -0.5, 9.5, 50, 100., 3500.);
		TH2D *hNjet_vs_Ejets_det = new TH2D("hNjet_vs_Ejets_det", "Number of jets versus E leading jet (det)", 	10, -0.5, 9.5, 50, 100., 3500.);								
		// Valid gen jets versus Castor jets.
		TH2D *hNumber_of_match_jets = new TH2D("hNumber_of_match_jets", "Castor jet versus Gen jets;N_{Castor};N_{GEN}", 11, -0.5, 10.5, 11, -0.5, 10.5);

		// Compare phi on gen and det level.
		TH2D *hPhi_gen_det = new TH2D("hPhi_gen_det", "#varphi_{gen} vs. #varphi_{det};#varphi_{gen};#varphi_{det}", 50, -6.29, 6.29, 50, -6.29, 6.29);
		
		// Needed for correct calibration.
		TH2D *hResponse 	= new TH2D("hResponse", 	"E_{det}/E_{gen};E_{det};#frac{E_{det}}{E_{gen}}", 		Ebins, Emin, Emax, 100, 0., 5.);			hResponse->Sumw2();
		TH2D *hResponse_gen 	= new TH2D("hResponse_gen", 	"E_{det}/E_{gen};E_{gen};#frac{E_{det}}{E_{gen}}", 		Ebins, Emin, Emax, 100, 0., 5.);			hResponse_gen->Sumw2();

		int Bins_response[3]	= {Ebins, 100, 16};
		double Min_response[3]	= {Emin, 0., -1.*PI};
		double Max_response[3]	= {Emax, 5., PI};

		THnSparseD *hResponse_phi 	= new THnSparseD("hResponse_phi", 	"E_{det}/E_{gen};E_{det};#frac{E_{det}}{E_{gen}};#varphi",	3, Bins_response, Min_response, Max_response); hResponse_phi->Sumw2();
		THnSparseD *hResponse_gen_phi 	= new THnSparseD("hResponse_gen_phi", 	"E_{det}/E_{gen};E_{gen};#frac{E_{det}}{E_{gen}};#varphi",    3, Bins_response, Min_response, Max_response); hResponse_gen_phi->Sumw2();
	//        TH3D *hResponse_phi     = new TH3D("hResponse_phi", 	"E_{det}/E_{gen};E_{det};#frac{E_{det}}{E_{gen}};#varphi", 	Ebins, Emin, Emax, 100, 0., 5., 48, -3.15, 3.15);     	hResponse_phi->Sumw2();        
	//        TH3D *hResponse_gen_phi = new TH3D("hResponse_gen_phi",	"E_{det}/E_{gen};E_{gen};#frac{E_{det}}{E_{gen}};#varphi", 	Ebins, Emin, Emax, 100, 0., 5., 48, -3.15, 3.15);//	hResponse_gen_phi->Sumw2();

		// Diego
		
		TH2D *hsumEreco_sumEgen 	= new TH2D("hsumEreco_sumEgen",	";#sum E_{gen};#frac{#sum E_{reco}}{#sum E_{gen}}", 50, 0., 3000., 50, 0., 3.);
		TH2D *hsumEreco_vs_sumEgen	= new TH2D("hsumEreco_vs_sumEgen", ";#sum E_{gen};#sum E_{reco}", 50, 0., 3000., 50, 0., 3000.);
		TH1D *hEreco_Egen_pions		= new TH1D("hEreco_Egen_pions", ";#frac{E_{reco}}{E_{gen}};N_{events}", 50., 0., 3.);


		int Bins_fine[3]	= {Ebins, Ebins*1000, 16};
		double Min_fine[3]	= {Emin, Emin, -1.*PI};
		double Max_fine[3]	= {Emax, Emax, PI};
		
		TH2D *hGen_fine 	= new TH2D("hGen_fine", 	"E_{gen}", Ebins, Emin, Emax, Ebins * 1000, Emin, Emax);			hGen_fine->Sumw2();
		THnSparseD * hGen_fine_phi = new THnSparseD("hGen_fine_phi",     "E_{gen}", 3, Bins_fine, Min_fine, Max_fine);				hGen_fine_phi->Sumw2();
		
		TH2D *hDet_fine 	= new TH2D("hDet_fine", 	"E_{det}", Ebins, Emin, Emax, Ebins * 1000, Emin, Emax);                        hDet_fine->Sumw2();
		THnSparseD * hDet_fine_phi = new THnSparseD("hDet_fine_phi",     "E_{gen}", 3, Bins_fine, Min_fine, Max_fine);                               hGen_fine_phi->Sumw2();


	cout << "TH3";

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
        hCastorJet_energy_response_fine->Sumw2();
        hCastorJet_energy_fakes->Sumw2();
        hCastorJet_energy_misses->Sumw2();
	hCastorJet_Matrix->Sumw2();

	hCastorJet_Matrix->Sumw2();

	hJRE->Sumw2();
	hJER->Sumw2();
	hJER_per_energy->Sumw2();
	hJER_per_distance->Sumw2();
	hMatched->Sumw2();
	hUnmatched->Sumw2();
	
	// -- First vector contains our energy values, second vector our calibration values.
	//vector<double> lowedge;
	//vector<double> muval;
        if( comments_ ){ cout << "// -- Filelabel is\t" << fileLabel_ << endl; }

	if( comments_ ) cout << "// -- Prepare calibration" << endl;
	std::map<int, vector<double> > lowedge;
	std::map<int, vector<double> > muval;

	if( comments_ ) cout << "// -- Prepare calibration II" << endl;
	FillCorrectionVectors( lowedge, muval);

	if( comments_ ) cout << "// -- Prepared calibation" << endl;

	cout << (lowedge[4]).size() << "\tis size of first vector" << endl;
        cout << (muval[4]).size() << "\tis size of first vector" << endl;


	TObjString* fn = 0;
	
	bool isMC = false;
    
	int counter_jer = 0;
	int counter_events = 0;
	int counter_match = 0;
	
	std::cout << "start looping over files\t" << filename_ << std::endl;
	if( isData_ ) cout << "This is data" << endl;

	TString EIstring_ = "off";
	TString calibrationfunc_ = "off";
	TString calibrationdisc_ = "off";
	TString prepareunfold_ = "off";

	if( cut_EI )	{ EIstring_= "on"; }
	if( do_calibration_function ){ calibrationfunc_ = "on"; }
        if( do_calibration_discrete ){ calibrationdisc_ = "on"; }
	if( prepare_unfolding ){ prepareunfold_ = "on"; }

	cout << "\n\n\n\t=====" << endl;
	cout << "\tEIcut turned " << EIstring_ <<  endl;
	cout << "\tCalibration with function turned " << calibrationfunc_ << endl;
	cout << "\tCalibration with values turned " << calibrationdisc_ << endl;
	cout << "\tPreparation for unfolding turned " << prepareunfold_ << endl;
        cout << "\n\n\n\t=====" << endl;	

	TFile * currentfile_ = new TFile( filename_, "Read");

	//////////////////////////////////////////////////
	// Get tree from the files and define all branches
	//////////////////////////////////////////////////
		
	// get tree from file
	TTree *tree;// = new TTree("CastorTree","");
	tree = (TTree*) currentfile_->Get("CastorTree");
		
	// start file loop
	int treesize = tree->GetEntriesFast();
	
	// for(int counter_events = 0; counter_events < totalEvents_ && counter_events < treesize; counter_events++ ) {		
	
 	// define objects and branches

	cout << "Got tree" << endl;
		
	std::vector<MyCastorJet> *CastorJets = NULL;
	TBranch *b_CastorJets = tree->GetBranch("CastorJets");
	b_CastorJets->SetAddress(&CastorJets);

        std::vector<MyGenJet> *CastorGenJets = NULL;
	TBranch *b_CastorGenJets = NULL;	
	if (!isData_) b_CastorGenJets = tree->GetBranch("CastorGenJets");
	if (!isData_) b_CastorGenJets->SetAddress(&CastorGenJets);
		
	int Nevents = tree->GetEntriesFast();
	std::cout << "file opened, events in this file = " << Nevents << std::endl;
	totalevents += Nevents;
	
	cout << "counter_events\t" << counter_events << "\ttotalEvents_\t" << totalEvents_ << "\tNevents\t" << Nevents << "\ttreesize\t" << treesize << endl;
	
	// start event loop

	int absurdFake = 0;
        int emptyCastor = 0;
	for( counter_events = 0; counter_events < totalEvents_ && counter_events < treesize; counter_events++ ) {
	       
	  if( comments_ ){ cout << "******************************************" << endl; }
	  if( counter_events%1000== 0){ cout << "\t" << counter_events << "\tpassed" << endl; }
	
	  /////////////////////////////////////////
	  // Do stuff before filters
	  /////////////////////////////////////////
	
   	  if( !isData_) { b_CastorGenJets->GetEntry( counter_events ); }
	  b_CastorJets->GetEntry( counter_events );	


				
	  /////////////////////////////////////////
	  // Start Nvertex == 1 part of the code 
	  /////////////////////////////////////////	  		

  	  // Match DET and GEN jets. 
  	  if( !isData_ ) { 
	    hNumber_of_match_jets->Fill( CastorJets->size(), CastorGenJets->size());             
	    if( ((*CastorGenJets)[0]).Energy() > threshold_ ){ hGenJet_energy_lead->Fill( ((*CastorGenJets)[0]).Energy() ); }
	    hCastorJet_pt->Fill( ((*CastorGenJets)[0]).Pt() );
	  }
  	  int matched_pairs = 0;

          // -- MATCHING

	  /////////////////////////////////////////////////////////////////////////////////////////////////////
	  // A first look into the jets: prepare the matching by removing too soft jets and wrong-type jets. //
	  /////////////////////////////////////////////////////////////////////////////////////////////////////

	  if( comments_) cout << "\n\n\n\t" << counter_events << endl;

	  for(int det_jet = CastorJets->size()-1; det_jet >= 0; det_jet-- ){
	    MyCastorJet castor_jet = (*CastorJets)[ det_jet ];
	    double det_energy = castor_jet.energy;

	    if( do_calibration_function){
	      double det_phi = castor_jet.phi;
	      int sector = CastorSector( det_phi ) ; 
 	      det_energy = CalibratedDet( det_energy, sector, fileLabel_, threshold_ );
	    }
	    if( det_jet == 0 && det_energy < threshold_ ){ emptyCastor++; }		

	    if( comments_ ){ cout << "Type is\t" << GetJetType( castor_jet ) << "\tfor\t" << counter_events << "\t" << det_jet << "\t" << det_energy << endl; }
/*
	    if(  det_energy < threshold_ ){ 
	      CastorJets->erase( CastorJets->begin() + det_jet ); 
	      if( comments_ ){ cout << "Type is\t" << GetJetType( castor_jet ) << "\tfor\t" << counter_events << "\t" << det_jet << "\t" << det_energy << endl; }
	    }
*/
	  } // Close loop over detector level jets.




	  if( !isData_ ){
	    for( int gen_jet =  CastorGenJets->size()-1; gen_jet >=0; gen_jet--){
	      MyGenJet castor_gen = (*CastorGenJets)[gen_jet];
	      double gen_energy = castor_gen.Energy();
	      hGenJet_energy->Fill( gen_energy );
		
	      //if( gen_energy < threshold_){
	      //  CastorGenJets->erase( CastorGenJets->begin() + gen_jet ); 
	      //}
	    }
          }

	  if( comments_ && !isData_){ cout << "\n\n\n$$$\tEvent\t" << counter_events << "\tdet size\t" << CastorJets->size() << "\tgen size\t" << CastorGenJets->size() << endl;}


  	  while( matched_pairs == 0 && CastorJets->size() > 0 ){
	    if( comments_ ){ cout << "First pair to match" << endl; }	

	    ////////////////////////////////////////////////////
	    // Get information on leading jet (type, energy). //
	    ////////////////////////////////////////////////////

	    if( !prepare_unfolding_ ){
	      if( comments_ ){ cout << "time to not prepare unfolding" << endl; }
	      MyCastorJet castorjet = (*CastorJets)[ 0 ];
	      double castor_sectors = castorjet.ntower;

  	      TString detjettype = GetJetType( castorjet );
	      if( comments_ ){ cout << "Jettype is\t" << detjettype << "\tfor required\t" << jettype_ << "\tsector number is\t" << castor_sectors << endl; }
	      if( jettype_ != "all"){ 			// We wish to differentiate jet types.
 	        if( detjettype == "other"){ break; }	// Not had or em.
	        if( detjettype != jettype_ ){ break; }	// Not our jettype.
	      }
	      if( comments_ ){ cout << "Survived type cut" << endl; }

	      if( sectors_ != 0 ){
	        if( (sectors_ != 1 && castor_sectors == 1) ){ break; }
	        if( (sectors_ == 1 && castor_sectors != 1) ){ break; }
	      }
	      if( comments_ ){ cout << "Survived sectors cut" << endl; } 
	    } // Get type and Nsectors of leading jet.

	    ////////////////////////////////////////////////////////
	    // -- Get the detector level jet energy distribution. //
	    ////////////////////////////////////////////////////////

	    for(int det_jet = 0; det_jet < CastorJets->size(); det_jet++){


	      MyCastorJet castor_det = (*CastorJets)[det_jet];
	      TString detjettype = GetJetType(castor_det);
	      double det_energy = castor_det.energy;	

	      if( do_calibration_function ){
	        double det_phi = castor_det.phi;
	        int sector = CastorSector( det_phi ) ; 
 	        det_energy = CalibratedDet( det_energy, sector, fileLabel_, threshold_ );
	      }
	      if( det_energy > threshold_ ) {
 	       hCastorJet_energy->Fill( det_energy );
	       if( det_jet == 0)  hCastorJet_energy_lead->Fill( det_energy ); 
	      }
	    } // End detector level jet energy distribution.
	
    	    /////////////////////////////////////////////////
    	    //-- Enter all jets in response object.        //
    	    //-- Only needed as preparation for unfolding. //
    	    /////////////////////////////////////////////////	
   			
	    MyCastorJet leading_det;
	    MyGenJet leading_gen;
	    int leading_gen_ = -1;

	    //-------------------------------------//
	    //-------------------------------------//
	    // Will we be unfolding? Yes, we will! //
	    //-------------------------------------//
	    //-------------------------------------//

	    if( prepare_unfolding_ ){				
              if( comments_ ){ cout << "\n\n\t\tUnfolding" << endl; }	

	      double nFake = 0;
 	      double nMiss = 0;
	      double nMatch = 0;
	      double nGen = 0;
	      double nDet = 0;

	      for(int gen_jet = 0; gen_jet < CastorGenJets->size(); gen_jet++){ 
		MyGenJet castor_gen = (*CastorGenJets)[gen_jet];
		if( castor_gen.Energy() > threshold_ ){ nGen++; }
	      }

	      double leading_det_E, leading_gen_E;
	      bool isLeadingPair = true;

   	      for(int det_jet = 0; det_jet < CastorJets->size(); det_jet++){

       	        MyCastorJet castor_det = (*CastorJets)[det_jet];

		// CALIBRATION: Jet energy is calibrated before taking any further action.
                double det_energy = castor_det.energy;
		double det_phi = castor_det.phi; 
	        int sector = CastorSector( det_phi ) ;
		TString jettype = GetJetType( castor_det );

		det_energy = CalibratedDet( det_energy, sector, fileLabel_, threshold_ ); 

		if( comments_ ){ cout << counter_events << "\tDet\t" << det_jet << " of\t" << CastorJets->size() << "\t" << det_energy << "\t" << det_phi << "\t" << GetJetType(castor_det) << endl; }

		// CUT: cut jet if energy is below threshold.
		if( det_energy < threshold_ ){ 
                  if( comments_ ){ cout << "\t\t\t\t\t\t\t\tLow Det E\t" << det_jet << "\t" << det_energy << endl; }
		  continue; 
		} 
		nDet++;
		
    	        // MATCHING: Match closest in phi by looping over all jets and find one with lowest delta phi.
    	        double min_delta_phi = deltaphimax_ ;
		if( sectors_ == 1 ){ min_delta_phi = 0.1; }

    	        int match_gen = -1;		
		double max_energy = 0.;			  

    	        for(int gen_jet =  CastorGenJets->size()-1; gen_jet >=0; gen_jet--){

    	          MyGenJet castor_gen = (*CastorGenJets)[gen_jet];
		  double gen_energy = castor_gen.Energy();
		  double gen_eta = castor_gen.Eta();

		  // Remove and skip gen jet if not hard enough.
		  if( gen_energy < threshold_ || gen_eta < genetamin_ || gen_eta > genetamax_ ){ 
		    if( comments_ ){ cout << "\t\t\t\t\t\t\t\t\t\tgen_energy below threshold\t" << endl; }
		    continue; 
		  }
	
    	          double gen_phi = castor_gen.Phi();
    	          double delta_phi = fabs( det_phi - gen_phi ); if( delta_phi > PI ){ delta_phi = 2.*PI - delta_phi; }

		  if( comments_ ){ cout << "\t\t\t\t\tGen\t" << gen_jet << " of\t" << CastorGenJets->size() << "\t" << gen_energy << "\t" << gen_phi << endl; }
			
		  // Jet lies closest: save Delta Phi and index of the jet.
		  if( match_ == "matchPhi" ){		
    	            if( delta_phi < min_delta_phi ){
		      if( comments_  ){ cout << "\t\t\t\t\t\t\t\t\t\t\t\t\tNew delta phi minimum" << endl; }
    	              min_delta_phi = delta_phi;
    	              match_gen = gen_jet;
    	            } // New minimum.
		    if( comments_ ){ cout << "\t\t\t\t\t\t\t\t\t\t\t" << delta_phi << "\t" << min_delta_phi << endl; }
    	          } //Match E.


		  // Hardest jet close to det. jet: save Delta Phi, energy and index of the jet.
		  if( match_ == "matchE" ){		
    	            if( gen_energy > max_energy && delta_phi < deltaphimax_ ){
		      if( comments_  ){ cout << "\t\t\t\t\t\t\t\t\t\t\t\t\tNew delta phi minimum" << endl; }
    	              max_energy = gen_energy;
    	              match_gen = gen_jet;
    	            } // New minimum.
		    if( comments_ ){ cout << "\t\t\t\t\t\t\t\t\t\t\t" << delta_phi << "\t" << min_delta_phi << endl; }
    	          } //Match E.

		}// Loop over possible match candidates (gen).




    	        // We have a matching gen. jet
    	        if( match_gen != -1 ){
    	          MyGenJet castor_gen = (*CastorGenJets)[match_gen];
    	          double gen_energy = castor_gen.Energy();	

		  if( det_jet == 0){
		    leading_det = (*CastorJets)[0];
		    leading_gen = (*CastorGenJets)[match_gen];
		    leading_gen_ = match_gen;
		    response_lead.Fill( det_energy, gen_energy );
		    response_lead_fine.Fill( det_energy, gen_energy );
		    hCastorJet_energy_response_lead->Fill( det_energy, gen_energy );
		    leading_gen_E = gen_energy; 
		    leading_det_E = det_energy;
		    isLeadingPair = false;
		  }

		  // Feed information to response matrix and histograms.	  

    	          response.Fill( det_energy, gen_energy );
    	          response_match.Fill( det_energy, gen_energy );
		  response_fine.Fill( det_energy, gen_energy );
		  // TH2D!
		  hCastorJet_energy_response->Fill( det_energy, gen_energy );
		  //hGenJet_energy->Fill(gen_energy);	
                  hCastorJet_matchedGen_all->Fill( gen_energy );	
 
		  nMatch++;

	          if( comments_ ){ cout << "\t\t\t\t\tMatch" << endl << endl; }

		  // Because the Gen jets are not sorted as the det jets, it is easiest to delete the matched gen jet.
    	          CastorGenJets->erase( CastorGenJets->begin() + match_gen );
    	        } // Matching gen jet.

		//-----------------------------------------------------------//
    	        // FAKES: We don't have a matching gen. jet. This is a fake. //
		//-----------------------------------------------------------//

    	        else{
    	            response.Fake( det_energy );
		    hCastorJet_fake_all->Fill( det_energy );
		    response_fine.Fake( det_energy ); 

		    if( det_jet == 0 ){ 
		      hCastorJet_fake_lead->Fill( det_energy ); 
		      response_lead.Fake( det_energy ); 
		     response_lead_fine.Fake( det_energy ); 
		    }
		    if( comments_ ){ cout << counter_events << "\t\t\t\t\tFake" << endl << endl; }
		    nFake++;

		    if( det_energy > 1400. ){ cout << "Absurd fake at\t" << counter_events << endl; }
    	        } // FAKES.
    	      }
 
	      //----------------------------------------------------------------------------------//
              // MISSES: We had all det. jets, time for the remaining gen. jets which are misses. //
	      //----------------------------------------------------------------------------------//

    	      if(  CastorGenJets->size() > 0 ){

    	        for(int gen_jet = 0; gen_jet < CastorGenJets->size(); gen_jet++){

    	          MyGenJet castor_gen = (*CastorGenJets)[gen_jet];
    	          double gen_energy = castor_gen.Energy();
		  double gen_eta = castor_gen.Eta();
		  if( gen_energy < threshold_ || gen_eta < genetamin_ || gen_eta > genetamax_ ){ 
		    continue; 
		  }
		  response.Miss(gen_energy ); 
	 	  response_fine.Miss( gen_energy ); 
	 	  //hGenJet_energy->Fill( gen_energy );	
		  hCastorJet_miss_all->Fill( gen_energy );
		  nMiss++;

	          if( comments_ ){ cout << "\t\t\t\t\t\t\t" << gen_jet << " of\t" << CastorGenJets->size() << "\t" << gen_energy << "\t" << "Miss" <<  endl; }
    	        }
    	      } // MISSES.
	      double response_input[7] = {leading_det_E, leading_gen_E, nMatch, nFake, nMiss, nDet, nGen};
	      hResponse_leading->Fill( response_input );

	    } // prepare_unfolding_					
			  
	    //----------------------------//
	    //----------------------------//
 	    // We done 'folding? We done. //
	    //----------------------------//
	    //----------------------------//

	    //------------------------------------------------//
	    // -- Match jets without filling response matrix. //
	    //------------------------------------------------//

	    if( !prepare_unfolding_ && !isData_ ){
	      if( comments_ ){ cout << "Time to look at jets" << endl; }

	      // -- Leading detector level jet.
	      MyCastorJet castor_det = (*CastorJets)[0];
	      double det_energy = castor_det.energy;
	      double phi_det = castor_det.phi;
	      double min_delta_phi = 0.2;
	     
	      if( sectors_ == 1 ){ min_delta_phi = 0.1; }	     

	      // -- Look for closest generator level jet.
	      for(int gen_jet = 0; gen_jet < CastorGenJets->size(); gen_jet++){

		MyGenJet current_gen = (*CastorGenJets)[gen_jet];
		double e_gen = current_gen.Energy();
		if( e_gen < threshold_ ) { continue; }	// Do not continue if E is below threshold.

		double phi_gen = current_gen.Phi();

		double delta_phi = fabs(phi_det - phi_gen);	if( delta_phi > 2.*PI ){ delta_phi = delta_phi - 2.*PI; }
		if( delta_phi < min_delta_phi ){
		  min_delta_phi = delta_phi;
		  leading_gen_ = gen_jet;
		   if( comments_ ){ cout << "Match candidate " << gen_jet << "\t" << current_gen.Eta() << endl; }
		} // Minimum in phi.
	      } // Loop over gen. jets.
	    } // -- No unfolding.

	    ////////////////////////////////////////////////////////////
	    //-- Continue if there has been a match with leading jet. //
	    ////////////////////////////////////////////////////////////

	    if( leading_gen_ == -1 ){
	      if( comments_ ){ cout << "$$$**\tEvent\t" << counter_events << "No match with leading detector level jet\t" << endl; }		
	      break; 
	    }

	    MyCastorJet castor_jet = (*CastorJets)[0];
	    double phi_det_lead = castor_jet.phi;
	    double eta_det_lead = castor_jet.eta;
	    double det_energy_lead = castor_jet.energy;
	    int sectors_lead = castor_jet.ntower ;
	    int sector = CastorSector( phi_det_lead ) ;
	
	    if( !prepare_unfolding_){ leading_gen = (*CastorGenJets)[leading_gen_]; }

	    if( comments_ ){ cout << "Events\t" << counter_events << "\tMatch with\t" << "\tEta\t" << leading_gen.Eta() << endl; cout << "Genjet\t" << leading_gen.Phi() << endl; }
					  
            double eta_gen = leading_gen.Eta(); 		
	    if ( eta_gen < (-5.9 - etaband_/2.) || eta_gen > (-5.9 + etaband_/2.) ){ 
	      if( comments_) {cout << "Non-contained" << endl;}
	      break; }
            double phi_gen = leading_gen.Phi();
            double phidiff = fabs(phi_det_lead-phi_gen);     		if( phidiff > PI ){ phidiff = 2.*PI - phidiff; }
            double etadiff = eta_det_lead-eta_gen;
            double R_diff = sqrt( etadiff*etadiff + phidiff*phidiff );
            double gen_energy = leading_gen.Energy();	
				  
            //////////////////////////////////
	    // Request isolated jet energy. //
	    //////////////////////////////////

	    // Do for isolated and calibrated.	      
            double E_gen_cut = 0.;	    
            if( cut_EI && !prepare_unfolding_){
              int genjet = 0;
	      IsolationCut isolation_energy( phi_det_lead, sectors_lead );
	      
              while (E_gen_cut == 0. && genjet < CastorGenJets->size()){
                if( genjet == leading_gen_ ){ genjet++; continue; }
                MyGenJet genjet_ = (*CastorGenJets)[ genjet ];
					 
                int genpart = 0;
                while( genpart <  (genjet_.JetPart).size() && E_gen_cut == 0.){
                  double phipart = ( (genjet_.JetPart)[genpart]).Phi();
				       					 
                    if( isolation_energy.EnergyContamination( phipart ) ){ 
                      E_gen_cut += ( (genjet_.JetPart)[genpart]).Energy();
		    } // Check if energy is present in Castor sectors.
		    genpart++;
		  } // Loop over Gen. parts.					    					  
		  genjet++;
		} // Loop over Gen. jets.

                hIsolationEnergy->Fill( sectors_lead, E_gen_cut );
                hIsolationEnergy_1D->Fill( E_gen_cut );
                hIsolationEnergy_Egen->Fill( E_gen_cut / gen_energy );
                hIsolationEnergy_Edet->Fill( E_gen_cut / det_energy_lead );				  
              } // Cut on EI.

	      if( E_gen_cut > 0 && cut_EI ){ break; }

	      //////////////////////
	      // Fill histograms. //
	      //////////////////////
					  
	      if( comments_){ cout << "Events\t" << counter_events << "\tActual match\t" << endl; }

	      if( 	do_calibration_discrete ){ 			det_energy_lead = CalibratedDet(lowedge, muval, det_energy_lead, sector);       	}
	      else if( 	do_calibration_function && !sector_dependence){	det_energy_lead = CalibratedDet( det_energy_lead ); 				}
              else if( 	do_calibration_function && sector_dependence){ 	det_energy_lead = CalibratedDet( det_energy_lead, sector, fileLabel_, threshold_ ); 			}					  
				  
	     if( comments_ ){  	cout << "\tFilling response matrix\t" << det_energy_lead << "\t" << gen_energy << "\tsector\t" << sector << endl; }

		hCastorJet_energy_response->Fill( det_energy_lead, gen_energy);
		hCastorJet_energy_response_fine->Fill( det_energy_lead, gen_energy);	

		if( comments_){ cout << "\tIntegrals\t" << hCastorJet_energy_response->Integral() << "\t" << hCastorJet_energy_response_fine->Integral() << endl; }

  		//hGenJet_energy->Fill( gen_energy );
	      
	      	if( comments_ ){ cout << "// -- FILL\t" << counter_events << "\t" << matched_pairs << endl; }

	      	//-- Fill the multidimensional histograms with detector level, generator level energy and the phi value of detector level jet.
	      	//-- Allows for phi dependent calibration later on.

	      	hGen_fine->Fill( gen_energy, gen_energy );
              	double gen_fine[3] = {gen_energy, gen_energy, phi_det_lead};	hGen_fine_phi->Fill( gen_fine );

	      	hDet_fine->Fill( gen_energy, det_energy_lead );
              	double det_fine[3] = {gen_energy, det_energy_lead, phi_det_lead};   hDet_fine_phi->Fill( det_fine );
	      	hResponse->Fill( det_energy_lead, det_energy_lead/gen_energy );
	  
	      	double response_det[3] = {det_energy_lead,det_energy_lead/gen_energy, phi_det_lead};
	      	double response_gen[3] = {gen_energy,	det_energy_lead/gen_energy, phi_det_lead};
	      	hResponse_phi->Fill( response_det );
	      	hResponse_gen->Fill( gen_energy, det_energy_lead/gen_energy );
	      	hResponse_gen_phi->Fill( response_gen );

	      	// -- Information on the matched jet pair.
	      	hDistance->Fill( R_diff );
	      	hPhiDiff -> Fill( phidiff );
	      	hEtaDiff -> Fill( etadiff );
	      	hEtaPhiDiff -> Fill( etadiff, phidiff );
	      	hEtaRDiff -> Fill( etadiff, R_diff );
	      	hPhiRDiff -> Fill( phidiff, R_diff );

	      	hCastorJet_energy_ratio->Fill( gen_energy/det_energy_lead, det_energy_lead);
		
  	        // -- 
	        hTower_phi		->Fill( sectors_lead, phi_det_lead );
					  
	        if( matched_pairs == 0){
		  double JER = -(gen_energy - det_energy_lead)/gen_energy;
                  double JER_eDet = (gen_energy - det_energy_lead)/det_energy_lead;

	          hJER		->Fill( JER );
	          hJER_per_energy	->Fill( gen_energy, JER);
	          hJER_per_eDet	->Fill( det_energy_lead, JER_eDet);
	          hJER_per_eGen       ->Fill( gen_energy, JER_eDet);
	          hJER_per_distance	->Fill( sqrt(pow(phi_det_lead - phi_gen, 2.) + pow(eta_det_lead - eta_gen, 2.)), JER );
	          hJER_per_eta	->Fill( eta_gen, JER);
	          hEnergy_per_eta	->Fill( gen_energy, eta_gen );
	          hEtaDiff		->Fill( eta_gen - eta_det_lead );

	          // Look at jet types and fill in DeltaE/E..
	          if( jettype_ == "had" ){
	            hCastorJet_Matrix_had_det	->Fill(det_energy_lead, gen_energy);
	            hJER_per_energy_had_det  	->Fill( gen_energy, JER);
	            hJER_per_eDet_had_det  	->Fill( det_energy_lead, JER_eDet);
	            hJER_had			->Fill( JER );    
	            if( sectors_lead == 1 ){					      
	              hCastorJet_Matrix_had_det_1sector->Fill(det_energy_lead, gen_energy);
	              hJER_per_energy_had_det_1sector  ->Fill( gen_energy, JER);
	              hJER_per_eDet_had_det_1sector  	 ->Fill( det_energy_lead, JER_eDet); 
	            }
	            else{
	              hCastorJet_Matrix_had_det_nsector->Fill(det_energy_lead, gen_energy);
	              hJER_per_energy_had_det_nsector  ->Fill( gen_energy, JER);
	              hJER_per_eDet_had_det_nsector  	 ->Fill( det_energy_lead, JER_eDet); 					      
	            }
	          }
	          else if( jettype_ == "em" ){
	            hCastorJet_Matrix_em_det	->Fill(det_energy_lead, gen_energy); 
	            hJER_per_energy_em_det  	->Fill( gen_energy, JER);
	            hJER_per_eDet_em_det	->Fill( det_energy_lead, JER_eDet);
	            hJER_em			->Fill( JER );					      

	            if( sectors_lead == 1 ){					      
	              hCastorJet_Matrix_had_det_1sector->Fill(det_energy_lead, gen_energy);
	              hJER_per_energy_had_det_1sector  ->Fill( gen_energy, JER);
	              hJER_per_eDet_had_det_1sector  	 ->Fill( det_energy_lead, JER_eDet); 
	            }
	            else{
	              hCastorJet_Matrix_had_det_nsector->Fill(det_energy_lead, gen_energy);
	              hJER_per_energy_had_det_nsector  ->Fill( gen_energy, JER);
	              hJER_per_eDet_had_det_nsector  	 ->Fill( det_energy_lead, JER_eDet); 					      
	            }
                  }
		} // No matched pairs yet.
	     // } // Cut on EI, if applied already.
	      break; // End loop over det jets.
    	    } // For loop over det jets.
				
            // end of event, print status
    	    if( ((counter_events + 1) % 10000) == 0) std::cout << counter_events+1 <<"events done in file " << std::endl;				
          } // end event loop		
		//delete tree;
//	} // end file loop
	
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
    

	// Create response matrix including misses and fakes.
	hCastorJet_energy_complete_response->Add( hCastorJet_energy_response);
        hCastorJet_energy_complete_response->Add( hCastorJet_energy_misses);
        hCastorJet_energy_complete_response->Add( hCastorJet_energy_fakes);
   
    // write all histo's to file
    
	std::cout << "total number of events = " << totalevents << " from " << it << " file(s)" << endl;
	std::cout << emptyCastor << " events without Castor jets.\t" << endl;
	
	// create output root file
	Char_t filename[200];
	//std::string first(outputname_);
        float etamargin = etaband_;         
	//TString datestring = date;

//	filename.Replace("Histograms/","");
//	filename.Replace("_STRIPPEDTREE.root","");
	if( LoopOutputFile_ == ""){
	  sprintf(filename, date_ + "_JetAnalyzer_etaband_%f_%i_%i_sectors_" + jettype_ + ".root", etamargin, counter_events, sectors_);
	LoopOutputFile_ = TString::Format(filename); 
        }
//	TFile* output = new TFile(filename,"RECREATE");
	TFile* output = new TFile(LoopOutputFile_ , "Recreate");
	cout << "Save as\t" << LoopOutputFile_ << endl;	
	output->cd();
		
	cout << "Absurd fakes:\t" << absurdFake << endl;

	//////////////////////////////////////////
	// Save all your histograms in a root file
	//////////////////////////////////////////
	
	// detector level histograms
    
	// eflow histos
	hCASTORTowerMulti->Write();
	hCASTOReflow->Write();
	h2CASTOReflow_grid->Write();
	
	for (int icha=0;icha<224;icha++) {
		hCASTOReflow_channel[icha]->Write();
	}

	hJRE->Divide( hMatched, hUnmatched );
	
	hCastorJet_energy->Write();
	hCastorJet_energy_lead->Write();
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
        hGenJet_energy_lead->Write();


        hCastorJet_energy_gen->Write();
	hCastorJet_pt_gen->Write();

        hCastorJet_energy_response->Write();
        hCastorJet_energy_response_fine->Write();
	hCastorJet_energy_response_lead->Write();
        hCastorJet_energy_fakes->Write();
        hCastorJet_energy_misses->Write();

	hCastorJet_miss_all->Write();
	hCastorJet_fake_all->Write();
	hCastorJet_fake_lead->Write();

        hCastorJet_matchedGen_all->Write();

        hCastorJet_energy_complete_response->Write();
	hResponse_leading->Write();

	hCastorJet_Matrix->Write();

	hCastorJet_energy_ratio->Write();

	hDistance->Write();
	hPhiDiff->Write();
        hEtaDiff->Write();
	hEtaPhiDiff->Write();
        hEtaRDiff->Write();
        hPhiRDiff->Write();
/*
        hDistance_all->Write();
        hPhiDiff_all->Write();
        hEtaDiff_all->Write();
        hEtaPhiDiff_all->Write();
        hEtaRDiff_all->Write();
        hPhiRDiff_all->Write();
*/
	hTower_phi->Write();

	hElectron_energy	->Write();
	hPion_energy		->Write();
	
	  hPi_e_ratio->Divide( hPion_energy, hElectron_energy );
	hPi_e_ratio		->Write();

	hNumber_of_match_jets	->Write();

	hJER			->Write();
	hJER_per_energy		->Write();
        hJER_per_distance	->Write();
        hJER_per_eta		->Write();
        hEnergy_per_eta		->Write();
/*
        hJER_all		->Write();
        hJER_per_energy_all	->Write();
        hJER_per_distance_all	->Write();
	hJER_per_eta_all	->Write();
	hEnergy_per_eta_all	->Write();
*/

	/* DeltaE/E as function of generator energy. */
	hJER_per_energy->Write();
	
	hJER_per_energy_had_det->Write();
	hJER_per_energy_em_det->Write();

	hJER_per_energy_had_det_1sector->Write();
	hJER_per_energy_em_det_1sector->Write();

	hJER_per_energy_had_det_nsector->Write();
	hJER_per_energy_em_det_nsector->Write();			
	

	/* DeltaE/E as function of detector energy. */

        hJER_per_eDet_had_det->Write();
        hJER_per_eDet_em_det->Write();

	hJER_per_eDet->Write();
	
	hJER_per_eDet_had_det_1sector->Write();
	hJER_per_eDet_em_det_1sector->Write();

	hJER_per_eDet_had_det_nsector->Write();
	hJER_per_eDet_em_det_nsector->Write();	
	
	/* DeltaE/E_det vs. E_gen   */
	hJER_per_eGen->Write();
	

	hCastorJet_Matrix_had_det->Write();
        hCastorJet_Matrix_em_det ->Write();
	
	hCastorJet_Matrix_had_det_1sector->Write();
        hCastorJet_Matrix_em_det_1sector ->Write();
	
	hCastorJet_Matrix_had_det_nsector->Write();
        hCastorJet_Matrix_em_det_nsector ->Write();			

	// Proper calibration.
	hGen_fine	->Write();
        hGen_fine_phi   ->Write();

	hDet_fine	->Write();
        hDet_fine_phi   ->Write();

	hResponse	->Write();
	hResponse_gen	->Write();
        hResponse_phi   ->Write();
        hResponse_gen_phi->Write();


	// Test phi.
	 hPhi_gen_det->Write();
	
	hMatched->Write();
	hUnmatched->Write();
	hJRE->Write();		

	hNjet_vs_Ejets_gen->Write();
        hNjet_vs_Ejets_det->Write();

	hIsolationEnergy->Write();

	hIsolationEnergy_1D->Write();
	hIsolationEnergy_Egen->Write();
	hIsolationEnergy_Edet->Write();

	hsumEreco_sumEgen->Write();
	hsumEreco_vs_sumEgen->Write();

	response.Write();
	response_lead.Write();
        response_fine.Write();
	response_match.Write();
	response_lead_fine.Write();

	output->Close();
	std::cout << "file " << filename << " created." << std::endl;
	std::cout << "LoopOutputFile_ is\t" << LoopOutputFile_ << endl;
	LoopOutputFile_ = filename;
        cout << totalevents << "\tevents" << endl;
	cout << counter_match << "\tmatched events" << endl;

}
	
void JetAnalyzer_radii_strippedTree::AfterLoopCalculations(TString file) {

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

TString JetAnalyzer_radii_strippedTree::getOutputFile() {
    return LoopOutputFile_;
}

TString JetAnalyzer_radii_strippedTree::getInputDir() {
	return inputdir_;
}

TString JetAnalyzer_radii_strippedTree::getCurrentFile() {
	return currentfile_;
}

void JetAnalyzer_radii_strippedTree::setCurrentTFile() {
	currentTFile_ = currentStaticTFile_;
}

void* JetAnalyzer_radii_strippedTree::OpenROOTFile(JetAnalyzer_radii_strippedTree* arg) {
	currentStaticTFile_ = TFile::Open(arg->getInputDir()+arg->getCurrentFile(),"READ");
	return 0;
}

int JetAnalyzer_radii_strippedTree::posLeadingGenJet(std::vector<MyGenJet> JetVector,double etacut,double minptcut) {
	
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

int JetAnalyzer_radii_strippedTree::posLeadingTrackJet(std::vector<MyTrackJet> JetVector,double etacut,double minptcut) {
	
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



