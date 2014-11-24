////////////////////////////////////////
//////// New CMSSW_4_2_X version ///////
////////////////////////////////////////

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "SystematicsMin.h"
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
#define GenJetRadius "ak5"
#define DetJetRadius "ak5"
#define GenJetContained 0.
#define nEvents 1000000
#define PI 3.14159265359
#define JES_corr 0.79

#include "../../../RooUnfold-1.1.1/src/RooUnfold.h"
//#include "../../../RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#endif



TFile *SystematicsMin::currentStaticTFile_ = new TFile();

SystematicsMin::SystematicsMin(TString inputdir, TObjArray* filelist, bool isData, const char* outputname) {
    
	std::cout << "constructing SystematicsMin class..." << std::endl;
	
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

SystematicsMin::~SystematicsMin() { }

void SystematicsMin::Loop() {

#ifdef __CINT__
  gSystem->Load("../../../RooUnfold-1.1.1/libRooUnfold.so");
#endif
	
	
	std::cout << " SystematicsMin Loop function is started " << std::endl;
	
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
	TH1D *hMatched = new TH1D("hMatched", "Castor and Gen jets match", Ebins, Emin, Emax);
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
		fileopener = new TThread("fileopener",(void(*) (void *))&OpenROOTFile,(SystematicsMin*) this);
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
		
		TBranch *b_evtid = tree->GetBranch("EvtId");
	        TBranch *b_evtkin = NULL;
        	if (!isData_) b_evtkin = tree->GetBranch("GenKin");
		TBranch *b_BeamSpot = tree->GetBranch("beamSpot");
		TBranch *b_HLTrig = tree->GetBranch("HLTrig");
		TBranch *b_L1Trig = tree->GetBranch("L1Trig");
		TBranch *b_vertices = tree->GetBranch("primaryVertex");
		TBranch *b_castorrechits = tree->GetBranch("castorRecHit");
		TBranch *b_castortowers = tree->GetBranch("castorTower");
		TBranch *b_castorjets = NULL; //tree->GetBranch("castorJet");

		if( string_det_radius.Contains("ak3") ) b_castorjets = tree->GetBranch("ak3castorJet");
                if( string_det_radius.Contains("ak5") ) b_castorjets = tree->GetBranch("ak5castorJet");
                if( string_det_radius.Contains("ak7") ) b_castorjets = tree->GetBranch("ak7castorJet");
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
                if (!isData_ && string_gen_radius.Contains("ak5")) b_genJets = tree->GetBranch("ak5GenJet");
                if (!isData_ && string_gen_radius.Contains("ak7")) b_genJets = tree->GetBranch("ak7GenJet");		
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
		cout << "After" << endl;
		b_PFJets->SetAddress(&PFJets);
		if (!isData_) b_genParts->SetAddress(&genParts);
		b_caloTowers->SetAddress(&caloTowers);
		if (!isData_) b_genJets->SetAddress(&genJets);
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
			
//%%%%%%%%%%%%%%%%%%%%%%%

//  vector<fastjet::PseudoJet> particles;
  // an event with three particles:   px    py  pz      E
  //particles.push_back( PseudoJet(   99.0,  0.1,  0, 100.0) ); 
  //particles.push_back( PseudoJet(    4.0, -0.1,  0,   5.0) ); 
  //articles.push_back( PseudoJet(  -99.0,    0,  0,  99.0) );
 

//%%%%%%%%%%%%%%%%%%%%%%%





				
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
						double casjetE = casjet.energy * JES_corr; 
						if (casjetE > jetEThreshold) {
							NCastorJets++;
							hCastorJet_energy->Fill(casjetE);
								det_casjet.push_back(casjetE); 
								good_castorJets.push_back( casjet );
//cout << "\tCastorjet with energy\t" << casjet.energy << endl;
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
//cout << "+" << endl;
					// -------------------------
					// AVS - no Central jet required					
					if(!isData_){
					  for(int jet = 0; jet < genJets->size(); jet++){
					    MyGenJet currentJet = (*genJets)[jet];
					    //if( counter_events%1000 == 0) { cout << "\tGenjet\t" << jet << "\twith energy\t" << currentJet.Energy() << "\tand pT\t" << currentJet.Pt() << endl; }
					    double genjetE = currentJet.Energy() * JES_corr;
					    hGenJet_energy->Fill( genjetE );					    
					    if( currentJet.Eta() < (-5.2 - GenJetContained) && currentJet.Eta() > (-6.6 + GenJetContained) && genjetE > jetEThreshold){
					      gen_casjet_energy.push_back( genjetE );
					      hCastorJet_energy_gen->Fill( genjetE );
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

					// Match DET and GEN jets. 
					hNumber_of_match_jets->Fill( good_castorJets.size(), good_genJets.size());                                       

					int matched_pairs = 0;
					while( good_castorJets.size() > 0 && good_genJets.size() > 0){ // All jets need to be matched, or at least tried to be.

    					  MyCastorJet castorjet = good_castorJets[ 0 ];
    					  double eta_det = castorjet.eta;
    					  double phi_det = castorjet.phi;
    					  double det_energy = castorjet.energy * JES_corr;
			
					  int i_gen = 0;
					  bool matched = false;

  					  for( ; i_gen < good_genJets.size(); i_gen++){

  					    MyGenJet genjet_castor = (good_genJets)[i_gen];
					    double eta_gen = genjet_castor.Eta();
					    double phi_gen = genjet_castor.Phi();
					    double gen_energy = genjet_castor.Energy() * JES_corr;
					    double phidiff = fabs(phi_det-phi_gen);
					      if( phidiff > PI ){ phidiff = 2.*PI - phidiff; }
					    double etadiff = fabs(eta_det-eta_gen);

					    double distance = sqrt( etadiff*etadiff + phidiff*phidiff );

					    /* Fill in fake or match. */

					      response.Fill( det_energy, gen_energy);
					      response_onlyMatches.Fill( det_energy, gen_energy);
					      if( eta_gen < -5.2 && eta_gen > -5.6 && matched_pairs == 0){ hCastorJet_energy_response->Fill( det_energy, gen_energy); }
					      hCastorJet_energy_ratio->Fill( gen_energy/det_energy, det_energy);
					      /* Remove det and gen jet from vector. */
					      good_genJets.erase( good_genJets.begin() + i_gen );
					      good_castorJets.erase( good_castorJets.begin() + 0 );
					      matched = true;

					      hDistance->Fill( distance );
					      hPhiDiff -> Fill( phidiff );
					      hEtaDiff -> Fill( etadiff );					      
					      hEtaPhiDiff -> Fill( etadiff, phidiff );
                                              hEtaRDiff -> Fill( etadiff, distance );
                                              hPhiRDiff -> Fill( phidiff, distance );

					      matched_pairs++;
//					      if( matched_pairs == 1){
                                                double JER = (gen_energy - det_energy)/gen_energy;
                                                hJER_2->Fill( JER );
                                                hJER_per_energy_2->Fill( gen_energy, JER);
                                                hJER_per_distance_2->Fill( sqrt(pow(phi_det - phi_gen, 2.) + pow(eta_det - eta_gen, 2.)), JER );
						hJER_per_eta->Fill( eta_gen, JER);
						hEnergy_per_eta->Fill( gen_energy, eta_gen );
//					      }
					      break; // End loop over gen jets.
					  }// Loop gen jets.
					  if( !matched ){ // Det jet is a fake.
                                            response.Fake( det_energy );
                                            hCastorJet_energy_fakes->Fill( det_energy);
                                            good_castorJets.erase( good_castorJets.begin() +  0 );
					  }
					} // While loop.
					
					/* Fakes */
					for( int i_det = 0; i_det < good_castorJets.size(); i_det++){
					  MyCastorJet castorjet = (good_castorJets)[i_det];
					  double det_energy = castorjet.energy * JES_corr;
					  response.Fake( det_energy );
                                          hCastorJet_energy_fakes->Fill( det_energy);
					}

                                        /* Fakes */
                                        for( int i_gen = 0; i_gen < good_genJets.size(); i_gen++){
					  MyGenJet genjet_castor = (good_genJets)[i_gen];
                                          double gen_energy = genjet_castor.Energy() * JES_corr;
                                          response.Miss( gen_energy );
                                          hCastorJet_energy_misses->Fill( gen_energy);
                                        }


                                        
					//----------------------------
					// End AVS.			
					/*
					// AVS - Jet Energy Resolution & Jet Reconstruction Efficiency
					
					  1. Select only events with both Castorjets and GenJets in the Castor acceptance.
					  2. Couple the most energetic Castorjet with the most energetic GenJet.
					  3. Voilà: Jet Energy Resolution

                                          1. Select only events with both Castorjets and GenJets in the Castor acceptance.
					  2. Select hardest of each.
					  	3. Are they located in the same phi segment of Castor? \Delta\phi < 0.393
					  4. Voilà: Jet Reconstruction Energy.

					*/
/*	
	
					int castor_energetic = -1, gen_energetic = -1;
					int castor_hardest = -1, gen_hardest = -1;

					double castor_energy = 0., gen_energy = 0.;
					double castor_pt = 0., gen_pt = 0.;
					double castor_phi = 0., gen_phi = 0.;
					double gen_phi_2 = 0., gen_eta_2 = 0.;
					double castor_phi_2 = 0., castor_eta_2 = 0.;				
					double gen_eta = 0.;	


                                        if(!isData_){
					  // Find the right GenJets.
                                          for(int jet = 0; jet < genJets->size(); jet++){
                                            MyGenJet currentJet = (*genJets)[jet];
					    double curr_eta = currentJet.Eta();

					    // Right acceptance?
					    if( curr_eta > -5.2 || curr_eta < -6.6) continue;

					    double curr_energy = currentJet.E();
                                            double curr_pt = currentJet.Pt();
					    double curr_phi = currentJet.Phi();

					    if( curr_energy < jetEThreshold ) continue;

					    if( curr_energy > gen_energy ){
					      gen_energy = curr_energy;
					      gen_energetic = jet;
					      gen_phi_2 = curr_phi;
					      gen_eta_2 = curr_eta; 
  					    }

                                            if( curr_pt > gen_pt ){
                                              gen_pt = curr_pt;
					      //gen_energy = curr_energy;
                                              gen_hardest = jet;
					      gen_phi = curr_phi;
					      gen_eta = curr_eta;
                                            }
					  } // Finding the GenJets.
					}
								
                                        // Find the right CastorJets.
                                        for(int jet = 0; jet < CastorJets->size(); jet++){
                                          MyCastorJet currentJet = (*CastorJets)[jet];
                                          double curr_eta = currentJet.eta;
                                          double curr_energy = currentJet.energy;
                                          double curr_pt = curr_energy*sin(2*atan(exp(5.9)));
					  double curr_phi = currentJet.phi;

   					  if( curr_energy < jetEThreshold) continue;

                                          if( curr_energy > castor_energy ){
                                            castor_energy = curr_energy;
                                            castor_energetic = jet;
					    castor_phi_2 = curr_phi;
                                            castor_eta_2 = curr_eta;
                                          }
                                          if( curr_pt > castor_pt ){
                                            castor_pt = curr_pt;
                                            castor_hardest = jet;
					    castor_phi = curr_phi;
                                    
                                          }
                                        } // Finding the Castorjets.



					if( castor_energetic >= 0 && gen_energetic >= 0){ // Found jets!
					  double JER = (gen_energy - castor_energy)/gen_energy;
					  hJER->Fill( JER );
					  hJER_per_energy->Fill( JER, gen_energy);
					  hJER_per_distance->Fill( JER, sqrt(pow(castor_phi - gen_phi, 2.) + pow(castor_eta_2 - gen_eta, 2.) ) );

					  if( fabs( castor_phi - gen_phi ) < 0.393 ) hMatched->Fill( castor_energy );
					  hUnmatched->Fill( castor_energy );
					}
*/ // End JER/JRE separate section.
					// End Jet Energy Resolution & Jet Reconstruction Efficiency

// cout << "Past JRE and JER" << endl;
										
					//---------------------------				
					// AVS - Central + foward jet.
/*					
					if(!isData_){

  					  // Det.
					  // Do we see a leading trackjet?
					  MyTrackJet trackjet_cf_det = (*trackJets)[0];
					    double trackpt = trackjet_cf_det.pt_cal;	
//					    double track_E = trackjet_cf_det.energy;			
					  MyCastorJet castorjet_cf_det = (*CastorJets)[0];
					    double castor_E = castorjet_cf_det.energy;
					    double castorpt = castor_E*sin(2*atan(exp(5.9)));

					  if( castor_E < 500. ){ continue; } // Castor jet weak.
					  if( trackpt < castorpt ){ continue; } // Castor jet leads.
// cout << "\t\t\tLeadin trackjet" << endl;
					  // Gen
					  // Jet 0 needs to be a trackjet.
                                          if( genJets->size() == 0){ continue; }
					  MyGenJet currentJet = (*genJets)[0];	
					    double gen_E = currentJet.E();
					    double gen_pt = currentJet.Pt();
					    double gen_eta = currentJet.Eta();

					  if( fabs( gen_eta ) > 2.5 ){ continue; } // Not a trackjet.
// cout << "\t\t\tActually a trackjet" << endl;
				          //
					  // We have a leading trackjet on both gen and det level.
					  //

					  // We know our castorjets at DET level, we need those at GEN level.
					  vector<double> castorjet_gen_energy;
					  for( int gen = 0; gen < genJets->size(); gen++){
					    MyGenJet currentJet = (*genJets)[gen];
					    if( currentJet.E() > 500. ){
					      castorjet_gen_energy.push_back( currentJet.E() );
					    }
					  } // Loop over genjets to determine Castor( gen )

					  if( castorjet_gen_energy.size() == 0 || CastorJets->size() == 0){ continue; } // No castorjets at gen level
				
					  bool all_gen_jets = false, all_det_jets = false;
					  int det_index = 0, gen_index = 0;
					  while( !all_gen_jets || !all_det_jets ){
					
					    double det_energy = -1.;
					    if( det_index < CastorJets->size() ){ 
					      if( (*CastorJets)[det_index].energy > 500. ){
					      det_energy = (*CastorJets)[det_index].energy; 
					      }
					    }
					    else{ all_det_jets = true; }
					    det_index++;

                                            double gen_energy = -1.;
                                            if( gen_index < castorjet_gen_energy.size() ){ gen_energy = castorjet_gen_energy[gen_index]; }
					    else{ all_gen_jets = true;}
                                            gen_index++;

					    if( i%1000 == 0){ cout << "Gen - size, index, energy\t" << castorjet_gen_energy.size() << "\t" << gen_index << "\t" << gen_energy << endl; }
					    if( i%1000 == 0){ cout << "Det - size, index, energy\t" << CastorJets->size() << "\t" << det_index << "\t" << det_energy << endl; }
		
					    // Fill an actual match. //
					    if( det_index <= CastorJets->size() && gen_index <= castorjet_gen_energy.size() ){
					      response_cf.Fill( det_energy, gen_energy);
					      response_cf_onlyMatches.Fill( det_energy, gen_energy);
					      hCastorJet_cf_energy_response->Fill( det_energy, gen_energy);
					      if( i%1000 == 0){ cout << "\tMatch" << endl; }
					    }
					    // No matching det jet: Miss. //
					    else if( det_index > CastorJets->size() && gen_index <= castorjet_gen_energy.size() ){
					      response_cf.Miss( gen_energy );
                                              hCastorJet_cf_energy_misses->Fill( gen_energy);
					      if( i%1000 == 0){ cout << "\tMiss" << endl; }
					    }
					    // No matching gen jet: Fake. //
					    else if( det_index <= CastorJets->size() && gen_index > castorjet_gen_energy.size() ){
					      response_cf.Fake( det_energy );
                                              hCastorJet_cf_energy_fakes->Fill( det_energy);
					      if( i%1000 == 0){ cout << "\tFake" << endl; }
					    } 
					    if( i%1000 == 0){ cout << endl; }
					  } // While not all gen/det jet have been investigated

					} // if(!isData_)
*/
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
//	sprintf(filename,"output_SystematicsMin_%s",first.c_str());
        if( !isData_) { sprintf(filename,"20141114_Output_SystematicsMin_GEN_" + string_gen_radius + "_DET_" + string_det_radius + "_unsortedE_%i_%s", counter_events, first.c_str()); }
        else{ sprintf(filename,"20141114_Output_SystematicsMin_Data_" + string_det_radius + "_Nevents_%i_%s", counter_events, first.c_str()); }
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
	
void SystematicsMin::AfterLoopCalculations(TString file) {

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

TString SystematicsMin::getOutputFile() {
    return LoopOutputFile_;
}

TString SystematicsMin::getInputDir() {
	return inputdir_;
}

TString SystematicsMin::getCurrentFile() {
	return currentfile_;
}

void SystematicsMin::setCurrentTFile() {
	currentTFile_ = currentStaticTFile_;
}

void* SystematicsMin::OpenROOTFile(SystematicsMin* arg) {
	currentStaticTFile_ = TFile::Open(arg->getInputDir()+arg->getCurrentFile(),"READ");
	return 0;
}

int SystematicsMin::posLeadingGenJet(std::vector<MyGenJet> JetVector,double etacut,double minptcut) {
	
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

int SystematicsMin::posLeadingTrackJet(std::vector<MyTrackJet> JetVector,double etacut,double minptcut) {
	
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
