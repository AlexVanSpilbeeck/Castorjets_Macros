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
#define jetEThreshold 200.
#define EbinWidth 5.
#define EbinWidth_rel 1.4
#define PI 3.14159265359

#define cut_EI false

#define comments_ false
#define do_calibration_discrete false
#define do_calibration_function true
#define sector_dependence true

#define prepare_unfolding true

#include "../../../RooUnfold-1.1.1/src/RooUnfold.h"
//#include "../../../RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#endif





TFile *JetAnalyzer_radii_strippedTree::currentStaticTFile_ = new TFile();

JetAnalyzer_radii_strippedTree::JetAnalyzer_radii_strippedTree(TString inputdir, bool isData, const char* outputname, int totalEvents, TString date, TString filename, TString jettype, int sectors, float etaband) {

    
	std::cout << "constructing JetAnalyzer_radii_strippedTree class..." << std::endl;
	
    inputdir_ = inputdir;
    filename_ = filename;
    isData_ = isData;
    outputname_ = outputname;
    totalEvents_ = totalEvents;
    date_ = date;	
    jettype_ = jettype;
    sectors_ = sectors;
    etaband_ = etaband;
	std::cout << "initializing basic variables..." << std::endl;
    
    // initialize basic variables

    
    LoopOutputFile_ = "";

	currentfile_ = "";
	currentTFile_ = new TFile();
	
	std::cout << "all initialisations done, class constructed!" << std::endl;
    
}

JetAnalyzer_radii_strippedTree::~JetAnalyzer_radii_strippedTree() { }

void JetAnalyzer_radii_strippedTree::Loop() {

#ifdef __CINT__
  gSystem->Load("../../../RooUnfold-1.1.1/libRooUnfold.so");
#endif
	
	
	std::cout << " JetAnalyzer_radii_strippedTree Loop function is started " << std::endl;
	
	TString tstring = outputname_;
	std::cout << " TString outputname_ = " << tstring << std::endl;


 
	
    // reweight the MC in this case
    bool reweightMC = false;
    if (tstring.Contains("Reweighted")) {
        std::cout << "We will reweight the MC now !!!" << std::endl;
        reweightMC = true;
    }
    
	//
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
	double Emin = 0., Emax = 3000.;
	int Ebins = 50;


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
        TH1D *hCastorJet_energy = new TH1D("hCastorJet_energy","CastorJet energy distribution",	Ebins, Emin, Emax);
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
        TH2D *hCastorJet_energy_response = new TH2D("hCastorJet_energy_response","CastorJet energy distribution - Response", 	Ebins, Emin, Emax, Ebins, Emin, Emax);	// ,Ebins+10,Emin-10*EbinWidth,Emax,Ebins+10,Emin-10*EbinWidth,Emax);
        TH2D *hCastorJet_energy_fakes 	 = new TH2D("hCastorJet_energy_fakes",	 "CastorJet energy distribution - Fakes", 	Ebins, Emin, Emax, Ebins, Emin, Emax);	// ,Ebins+10,Emin-10*EbinWidth,Emax);
        TH2D *hCastorJet_energy_misses	 = new TH2D("hCastorJet_energy_misses",	 "CastorJet energy distribution - Misses", 	Ebins, Emin, Emax, Ebins, Emin, Emax);	// ,Ebins+10,Emin-10*EbinWidth,Emax);
	
	TH2D *hCastorJet_energy_complete_response = new TH2D("hCastorJet_energy_complete_response","CastorJet energy distribution - Complete Response", Ebins_fix, Emin_fix, Emax_fix, Ebins_fix, Emin_fix, Emax_fix);

	

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
	     TH2D *hJER_per_eGen_had_det_1sector = new TH2D("hJER_per_eGen_had_det_1sector", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);
             TH2D *hJER_per_eGen_em_det_1sector = new TH2D("hJER_per_eGen_em_det_1sector", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);	

	     // 2+ Castor sector
	     TH2D *hJER_per_eGen_had_det_nsector = new TH2D("hJER_per_eGen_had_det_nsector", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);
             TH2D *hJER_per_eGen_em_det_nsector = new TH2D("hJER_per_eGen_em_det_nsector", "#DeltaE/E for fixed energies;E_{gen};#DeltaE/E",   	Ebins, Emin, Emax, 200, -5, 5);	

        TH2D *hJER_per_distance = new TH2D("hJER_per_distance", "#DeltaE/E for distance;#DeltaR;#DeltaE/E", 	200, 0., 6.5, 
															200, -5, 5);
        TH2D *hJER_per_eta = new TH2D("hJER_per_eta", "#DeltaE/E for distance;#DeltaR;#DeltaE/E",               14,-6.6,-5.2,
                                                                                                                        200, -5, 5);
        TH2D *hEnergy_per_eta = new TH2D("hEnergy_per_eta", "E_{gen} vs. #eta of leading jet", Ebins, Emin, Emax,  14,-6.6,-5.2);

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

	RooUnfoldResponse response 		(hCastorJet_energy, hCastorJet_energy, "response");

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
cout << "Half";

	int Bins_fine[3]	= {Ebins, Ebins*1000, 16};
	double Min_fine[3]	= {0., 0., -1.*PI};
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
	for( counter_events = 0; counter_events < totalEvents_ && counter_events < treesize; counter_events++ ) {		
//	for( counter_events = 0; counter_events < treesize; counter_events++ ) {
	       
	  if( counter_events%1000== 0){ cout << "\t" << counter_events << "\tpassed" << endl; }
	
	  /////////////////////////////////////////
	  // Do stuff before filters
	  /////////////////////////////////////////
	
   	  if( !isData_) { b_CastorGenJets->GetEntry( counter_events ); }
	  b_CastorJets->GetEntry( counter_events );			
				
	  /////////////////////////////////////////
	  // Start Nvertex == 1 part of the code 
	  /////////////////////////////////////////
	  			
	  // only fill the histograms when there's 1 vertex (filter out pile-up)

	  // analyse castor jets
	  /*
	  int NCastorJets = 0;
	  vector<double> det_casjet, gen_casjet_energy;
	  for (unsigned int j=0;j<CastorJets->size();j++) {
	    MyCastorJet casjet = (*CastorJets)[j];
	    if (casjet.energy > jetEThreshold_det) {
	      NCastorJets++;
	      hCastorJet_energy->Fill(casjet.energy);
	      det_casjet.push_back(casjet.energy); 
	      CastorJets.push_back( casjet );
  	    
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
          */

	  ///////////////////////////////////
	  // Cut out jets below threshold. //
	  ///////////////////////////////////

	  if( comments_ && jetEThreshold > 0. ){
	    for( int idet = CastorJets->size()-1; idet > -1; idet--){
	      cout << "// --DET\t" <<  (*CastorJets)[idet].energy << endl;
	    }
            for( int igen = CastorGenJets->size()-1; igen > -1; igen--){
              cout << "// --GEN\t" <<  (*CastorGenJets)[igen].Energy() << endl;
            }	   
	    cout << "\n\n";
	  }

	  
	  if( jetEThreshold > 0. ){
	    for( int idet = CastorJets->size()-1; idet > -1; idet--){
	      if( (*CastorJets)[idet].energy < jetEThreshold ){ CastorJets->erase( CastorJets->begin() + idet); }
	    }
	    if( !isData_ ){
              for( int igen = CastorGenJets->size()-1; igen > -1; igen--){
                if( (*CastorGenJets)[igen].Energy() < jetEThreshold ){ CastorGenJets->erase( CastorGenJets->begin() + igen) ; }
              }
	    }
	  }

          if( comments_ && jetEThreshold > 0. ){
            for( int idet = CastorJets->size()-1; idet > -1; idet--){
              cout << "// --DET\t" <<  (*CastorJets)[idet].energy << endl;
            }
	    if( !isData_ ){
              for( int igen = CastorGenJets->size()-1; igen > -1; igen--){
                cout << "// --GEN\t" <<  (*CastorGenJets)[igen].Energy() << endl;
              }
	    }
	    cout << "\n\n\n\n" << endl;
          }


  	  // Match DET and GEN jets. 
  	  if( !isData_ ) { 
	    hNumber_of_match_jets->Fill( CastorJets->size(), CastorGenJets->size());             
	    hGenJet_energy_lead->Fill( ((*CastorGenJets)[0]).Energy() );
	    hCastorJet_pt->Fill( ((*CastorGenJets)[0]).Pt() );
	  }
  	  int matched_pairs = 0;

          // -- MATCHING
          
	  if( comments_ && !isData_){ cout << "\n\n\n$$$\tEvent\t" << counter_events << "\tdet size\t" << CastorJets->size() << "\tgen size\t" << CastorGenJets->size() << endl;}

  	  while( matched_pairs == 0 && CastorJets->size() > 0 ){	
	    if( comments_ ){ cout << "$$$\tEvents\t" << counter_events << endl; }
			  
	      /////////////////////////
	      // Leading Castor Jet. //
	      /////////////////////////
			
	    MyCastorJet castorjet = (*CastorJets)[ 0 ];
	
	    // Analyze jet type.

	    double depth_jet = castorjet.depth;	// depth
	    double fhot_jet = castorjet.fhot;		// fhot
	    double fem_jet = castorjet.fem;
	    double sigmaz_jet = castorjet.sigmaz;
	    double width_jet = castorjet.width;
	    double det_energy = castorjet.energy;				  
	    double det_phi = castorjet.phi;

	    TString detjettype = "other";
	    // Count pions.
	    if( ! (depth_jet > -14450. && det_energy < 175.) ){ // Most likely not a pion.
	      if( ! (depth_jet > -14460. && det_energy > 175.) ){ // Most likely not a pion.
	        if( ! (fem_jet > 0.95) ){ // Most likely no pion.
	          hPion_energy->Fill( det_energy );
	          detjettype = "had";
	        } // Not a pion.
	      } // Not a pion.
	    } // Not a pion.
				
	    // Count electrons.
	    if( ! (fhot_jet < 0.45) ){
	      if( ! (fem_jet < 0.9) ){
	        if( ! ( sigmaz_jet > 30. && det_energy < 75.) ){
	          if( ! ( sigmaz_jet > 40. && det_energy > 75. ) ){
	            if( ! ( width_jet > 11.5) ){
		      if( ! ( depth_jet < -14450. && det_energy < 125.) ){
			if( ! ( depth_jet < -14460. && det_energy > 125.) ){						    
			  hElectron_energy->Fill( det_energy );
			  detjettype = "em";
			}
	              }				
		    }
		  }
	        }
	      }
	    } // Not an electron.

	    double eta_det = castorjet.eta;
	    double phi_det = castorjet.phi;
	    double castor_sectors = castorjet.ntower;
	    double phi_diff_max = 0.2;
	    if( castor_sectors == 1 ){ phi_diff_max = 0.1; }
					  
	    if( comments_ ){ cout << "\t\tnTowers\t" << castor_sectors << "\tphi\t" << phi_det << "\ttype\t" << detjettype << endl; }

	    if( jettype_ != "all"){ 			// We wish to differentiate jet types.
 	      if( detjettype == "other"){ break; }	// Not had or em.
	      if( detjettype != jettype_ ){ break; }	// Not our jettype.
	    }

	    if( sectors_ != 0 ){
	      if( (sectors_ != 1 && castor_sectors == 1) ){ break; }
	      if( (sectors_ == 1 && castor_sectors != 1) ){ break; }
	    }
	
	    if( isData_ ){ 
	      for(int det_jet = 0; det_jet < CastorJets->size(); det_jet++){
		MyCastorJet castor_det = (*CastorJets)[det_jet];
  	        hCastorJet_energy->Fill( castor_det.energy ); 
	      }
	    }
	
    	    ////////////////////////////////////////
    	    // Enter all jets in response matrix. //
    	    ////////////////////////////////////////				

	    MyCastorJet leading_det;
	    MyGenJet leading_gen;
	    int leading_gen_ = -1;

	    if( prepare_unfolding_ ){				
              if( comments_ ){ cout << "\t\tUnfolding" << endl; }	
    	      for(int det_jet = 0; det_jet < CastorJets->size(); det_jet++){
	        //  cout << "$**\tEvent\t" << counter_events << "\tdet\t" << det_jet << "\tgen\t" << CastorGenJets->size() << endl;
       	        MyCastorJet castor_det = (*CastorJets)[det_jet];
    	        double det_phi = castor_det.phi;
    	        double det_energy = castor_det.energy;

    	        //Match closest in phi.
    	        double min_delta_phi = 0.2;
    	        int match_gen = -1;					  

    	        for(int gen_jet = 0; gen_jet < CastorGenJets->size(); gen_jet++){
    	          MyGenJet castor_gen = (*CastorGenJets)[gen_jet];
    	          double gen_phi = castor_gen.Phi();
    	          double delta_phi = fabs( det_phi - gen_phi ); if( delta_phi > PI ){ delta_phi = 2.*PI - delta_phi; }
	
		  if( comments_ ){ 
		    cout << "$$$\tEvent\t" << counter_events << "\t(det, phi)\t" 
			<< det_jet << "\t" << det_phi << "\t" <<  castor_det.energy << "\t(gen, phi, E, pT)\t" 
			<< gen_jet << "\t" << gen_phi << "\t" << castor_gen.Energy() << "\t" 
			<< castor_gen.Pt() << "\t" << castor_gen.Eta() << endl;
		  }
					
    	          if( delta_phi < min_delta_phi ){
    	            min_delta_phi = delta_phi;
    	            match_gen = gen_jet;
    	          } // End of new minimum.
    	        } // Loop over possible match candidates (gen).

    	        // We have a matching gen. jet
    	        if( match_gen != -1 ){
    	          MyGenJet castor_gen = (*CastorGenJets)[match_gen];
    	          double gen_energy = castor_gen.Energy();
    	          double det_energy = castor_det.energy;
		  double det_phi = castor_det.phi;

		  if( det_jet == 0){
		    leading_det = (*CastorJets)[0];
		    leading_gen = (*CastorGenJets)[match_gen];
		    leading_gen_ = match_gen;
		  }

    	          response.Fill( det_energy, gen_energy );
		  hCastorJet_energy->Fill( det_energy);
		  hGenJet_energy->Fill(gen_energy);
	          if( comments_ ){ cout << "$**\tEvent\t" << counter_events << "\tMatch\t" << det_jet << "\t" << det_energy << "\t\t" << match_gen << "\t" << gen_energy << endl; }

		  // Because the Gen jets are not sorted as the det jets, it is easiest to delete the matched gen jet.
    	          CastorGenJets->erase( CastorGenJets->begin() + match_gen );
		  if( comments_ ){ cout << "\t\t\t\t" << match_gen << "\terased" << endl; }
    	        }

    	        // We don't have a matching gen. jet. This is a fake.
    	        else{
    	          response.Fake( det_energy );
		  hCastorJet_energy->Fill( det_energy );
		  if( comments_ ){ cout << "$**\tEvent\t" << counter_events << "\tFake\t" << det_jet << "\t" << det_energy << endl; }
    	        }
    	      }
 
              // -- We had all det. jets, time for the remaining gen. jets.
    	      if(  CastorGenJets->size() > 0 ){
    	        for(int gen_jet = 0; gen_jet < CastorGenJets->size(); gen_jet++){
    	          MyGenJet castor_gen = (*CastorGenJets)[gen_jet];
    	          double gen_energy = castor_gen.Energy();
		  response.Miss( gen_energy );
	 	  hGenJet_energy->Fill( gen_energy );	
	          if( comments_ ){ cout << "$**\tEvent\t" << counter_events << "\tMiss\t" << gen_jet << "\t" << gen_energy << endl; }
    	        }
    	      }
	    } // prepare_unfolding_					
			  

	    ////////////////////////////////////////////////////
	    // -- Match jets without filling response matrix. //
	    ////////////////////////////////////////////////////

	    if( !prepare_unfolding_ && !isData_ ){

	      MyCastorJet castor_det = (*CastorJets)[0];
	      double phi_det = castor_det.phi;
	       double min_delta_phi = 0.2;	     

	      for(int gen_jet = 0; gen_jet < CastorGenJets->size(); gen_jet++){
		MyGenJet current_gen = (*CastorGenJets)[gen_jet];
		double phi_gen = current_gen.Phi();

		double delta_phi = fabs(phi_det - phi_gen);	if( delta_phi > 2.*PI ){ delta_phi = delta_phi - 2.*PI; }
		if( delta_phi < min_delta_phi ){
		  min_delta_phi = delta_phi;
		  leading_gen_ = gen_jet;
		   if( comments_ ){ cout << "Match candidate " << gen_jet << "\t" << current_gen.Eta() << endl; }
		} // Minimum in phi.
	      } // Loop over gen. jets.
	    } // Loop over det. jets.

	    /////////////////////////////////////////
	    // Continue if there has been a match. //
	    /////////////////////////////////////////

	    if( leading_gen_ == -1 ){
	      if( comments_ ){ cout << "$$$**\tEvent\t" << counter_events << "No match\t" << det_energy << endl; } 		
	      break; 
	    }

	    if( !prepare_unfolding_){ leading_gen = (*CastorGenJets)[leading_gen_]; }

            if( comments_ ){ 
	      cout << "Events\t" << counter_events << "\tMatch with\t" << "\tEta\t" << leading_gen.Eta() << endl; 
	      cout << "Genjet\t" << leading_gen.Phi() << endl;
	    }
					  
            double eta_gen = leading_gen.Eta(); 		
	    if ( eta_gen < (-5.9 - etaband_/2.) || eta_gen > (-5.9 + etaband_/2.) ){ 
	      if( comments_) {cout << "Non-contained" << endl;}
	      break; }
            double phi_gen = leading_gen.Phi();
            double phidiff = fabs(phi_det-phi_gen);     		if( phidiff > PI ){ phidiff = 2.*PI - phidiff; }
            double etadiff = eta_det-eta_gen;
            double R_diff = sqrt( etadiff*etadiff + phidiff*phidiff );
            double gen_energy = leading_gen.Energy();	
				  
            //////////////////////////////////
	    // Request isolated jet energy. //
	    //////////////////////////////////
	      
            if( cut_EI && !prepare_unfolding_){
              double E_gen_cut = 0.;
              int genjet = 0;
	      IsolationCut isolation_energy( phi_det, castor_sectors );
	      
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

                hIsolationEnergy->Fill( castor_sectors, E_gen_cut );
                hIsolationEnergy_1D->Fill( E_gen_cut );
                hIsolationEnergy_Egen->Fill( E_gen_cut / gen_energy );
                hIsolationEnergy_Edet->Fill( E_gen_cut / det_energy );

                // If energy around jet is too high, continue.
                if( E_gen_cut > 0.0 * det_energy ){ 
		  if( comments_ ){ cout << "\t\tNon-isolated energy" << endl; }
		  break; }					  
              } // Cut on EI.

	      //////////////////////
	      // Fill histograms. //
	      //////////////////////
					  
	      if( comments_){ cout << "Events\t" << counter_events << "\tActual match\t" << endl; }
	      int sector = CastorSector( det_phi ) ;

	      if( do_calibration_discrete ){ 
		det_energy = CalibratedDet(lowedge, muval, det_energy, sector);
	      }
	      else if( do_calibration_function && !sector_dependence){ det_energy = CalibratedDet( det_energy ); }
              else if( do_calibration_function && !sector_dependence){ det_energy = CalibratedDet( det_energy, sector ); }

					  
	      //response.Fill( det_energy, gen_energy); 
					    
	      hCastorJet_energy_response->Fill( det_energy, gen_energy);

	      if(!prepare_unfolding_){
		hCastorJet_energy->Fill( det_energy );
	        hGenJet_energy->Fill( gen_energy );
//		hCastorJet_pt->Fill( gen_pt ); 
	      }

	      if( comments_ ){ cout << "// -- FILL\t" << counter_events << "\t" << matched_pairs << endl; }

	      hGen_fine->Fill( gen_energy, gen_energy );
              double gen_fine[3] = {gen_energy, gen_energy, phi_det};	hGen_fine_phi->Fill( gen_fine );

	      hDet_fine->Fill( gen_energy, det_energy );
              double det_fine[3] = {gen_energy, det_energy, phi_det};   hDet_fine_phi->Fill( det_fine );

	      hResponse->Fill( det_energy, det_energy/gen_energy );
	  
	      double response_det[3] = {det_energy, det_energy/gen_energy, phi_det};
	      double response_gen[3] = {gen_energy, det_energy/gen_energy, phi_det};
	      hResponse_phi->Fill( response_det );
	      hResponse_gen->Fill( gen_energy, det_energy/gen_energy );
	      hResponse_gen_phi->Fill( response_gen );

	      hDistance->Fill( R_diff );
	      hPhiDiff -> Fill( phidiff );
	      hEtaDiff -> Fill( etadiff );
	      hEtaPhiDiff -> Fill( etadiff, phidiff );
	      hEtaRDiff -> Fill( etadiff, R_diff );
	      hPhiRDiff -> Fill( phidiff, R_diff );

	      hCastorJet_energy_ratio->Fill( gen_energy/det_energy, det_energy);
	      /* Remove det and gen jet from vector. */


	      // matched = true;
	      
	      counter_match++;
		
	      // -- looking into JER information.
//	      double JER = -(gen_energy - det_energy)/gen_energy;
//	      double JER_eDet = (gen_energy - det_energy)/det_energy;					

	      hTower_phi		->Fill( castor_sectors, phi_det );
					  
	      if( matched_pairs == 0){
		double JER = -(gen_energy - det_energy)/gen_energy;
                double JER_eDet = (gen_energy - det_energy)/det_energy;

	        hJER		->Fill( JER );
	        hJER_per_energy	->Fill( gen_energy, JER);
	        hJER_per_eDet	->Fill( det_energy, JER_eDet);
	        hJER_per_eGen       ->Fill( gen_energy, JER_eDet);
	        hJER_per_distance	->Fill( sqrt(pow(phi_det - phi_gen, 2.) + pow(eta_det - eta_gen, 2.)), JER );
	        hJER_per_eta	->Fill( eta_gen, JER);
	        hEnergy_per_eta	->Fill( gen_energy, eta_gen );
	        hEtaDiff		->Fill( eta_gen - eta_det );

	        // Look at jet types and fill in DeltaE/E..
	        if( detjettype == "had" ){
	          hCastorJet_Matrix_had_det	->Fill(det_energy, gen_energy);
	          hJER_per_energy_had_det  	->Fill( gen_energy, JER);
	          hJER_per_eDet_had_det  	->Fill( det_energy, JER_eDet);
	          hJER_had			->Fill( JER );    
	          if( castor_sectors == 1 ){					      
	            hCastorJet_Matrix_had_det_1sector->Fill(det_energy, gen_energy);
	            hJER_per_energy_had_det_1sector  ->Fill( gen_energy, JER);
	            hJER_per_eDet_had_det_1sector  	 ->Fill( det_energy, JER_eDet); 
	          }
	          else{
	            hCastorJet_Matrix_had_det_nsector->Fill(det_energy, gen_energy);
	            hJER_per_energy_had_det_nsector  ->Fill( gen_energy, JER);
	            hJER_per_eDet_had_det_nsector  	 ->Fill( det_energy, JER_eDet); 					      
	          }
	        }
	        else if( detjettype == "em" ){
	          hCastorJet_Matrix_em_det	->Fill(det_energy, gen_energy); 
	          hJER_per_energy_em_det  	->Fill( gen_energy, JER);
	          hJER_per_eDet_em_det	->Fill( det_energy, JER_eDet);
	          hJER_em			->Fill( JER );					      

	          if( castor_sectors == 1 ){					      
	            hCastorJet_Matrix_had_det_1sector->Fill(det_energy, gen_energy);
	            hJER_per_energy_had_det_1sector  ->Fill( gen_energy, JER);
	            hJER_per_eDet_had_det_1sector  	 ->Fill( det_energy, JER_eDet); 
	          }
	          else{
	            hCastorJet_Matrix_had_det_nsector->Fill(det_energy, gen_energy);
	            hJER_per_energy_had_det_nsector  ->Fill( gen_energy, JER);
	            hJER_per_eDet_had_det_nsector  	 ->Fill( det_energy, JER_eDet); 					      
	          }
                }
	      } // No matched pairs yet.
	      matched_pairs++;

	      break; // End loop over gen jets.
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
	
	
	// create output root file
	Char_t filename[200];
	std::string first(outputname_);
        float etamargin = etaband_;         
	//TString datestring = date;

	sprintf(filename, date_ + "_JetAnalyzer_etaband_%f_%i_%i_sectors_" + jettype_ + ".root", etamargin, counter_events, sectors_);
	TFile* output = new TFile(filename,"RECREATE");
	output->cd();
		
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
        hCastorJet_energy_fakes->Write();
        hCastorJet_energy_misses->Write();
        hCastorJet_energy_complete_response->Write();

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

	response.Write();
//	response_all.Write();
//        response_cf.Write();
//	response_cf_onlyMatches.Write();

	output->Close();
	std::cout << "file " << filename << " created." << std::endl;
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

