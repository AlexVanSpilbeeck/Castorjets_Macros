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
#include "../Functions/Function_JetType.h"

//STANDARD ROOT INCLUDES

#include <TROOT.h>
#include <TBranch.h>
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
#define JES_unc 0.150

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

    if (filename.Contains("ak3ak3") || filename.Contains("ak3ak5") || filename.Contains("ak3ak7") ) etaband_ = 0.8;

    if( setup_ == "raw_wide" || setup_ == "unfold" || isData_) etaband_ = 1.4;

    cut_EI = false;
    if( setup_ == "isolated" || setup_ == "calibrated") cut_EI = true;

    do_calibration_function = false;
    if( setup_ == "calibrated" || setup_ == "unfold") do_calibration_function = true;

    if( setup_ == "unfold") prepare_unfolding = true;

    std::cout << "initializing basic variables..." << std::endl;
    
    // initialize basic variables
    LoopOutputFile_ = "";

    if( date == "0" && !isData_){       
      LoopOutputFile_ += filename;
      LoopOutputFile_.ReplaceAll( "Stripped_trees/", "");
      LoopOutputFile_.ReplaceAll( "STRIPPED_TREE.root", "");
      LoopOutputFile_ += setup;    
      LoopOutputFile_ += TString::Format("_Emin_%f", threshold_);
      LoopOutputFile_ += TString::Format("_deltaPhiMax_%f", deltaphimax_ );
      LoopOutputFile_ += TString::Format("_etaband_%f", etawidth);
      LoopOutputFile_ += "_" + jettype_;
      LoopOutputFile_ += "_" + match_;
      LoopOutputFile_ += ".root"; 
    }

    if( date == "0" && isData_ ){
      LoopOutputFile_ += filename;
      LoopOutputFile_.ReplaceAll( "Stripped_trees/", "");
      LoopOutputFile_.ReplaceAll( "STRIPPED_TREE.root", "");
      LoopOutputFile_ += setup;    
      LoopOutputFile_ += TString::Format("_Emin_%f", threshold_);
      LoopOutputFile_ += "_" + jettype_;
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
        int Ebins = 0;

        if( threshold_ == 0. ){ Emin = threshold_; Ebins = 28; }
	if( threshold_ == 150.){ Emin = threshold_; Ebins = 26;}
	if( threshold_ == 300.){ Emin = threshold_; Ebins = 24;}

         

        if( comments_ ){ cout << Ebins << " bins from " << Emin << " GeV to " << Emax << " GeV" << endl; }

	/* If we want variable Ebins. */

	//-------------------------------------------------------------------//
	//-- We want bin width 75 GeV, and bin width 300 GeV above 1200 GeV. //
	//-------------------------------------------------------------------//
        /*
	double current_lowE = 150.;
	vector<double> binEdges;
	binEdges.push_back( 150. );
	double binwidth = 75.;
	Ebins = 1;
        double switchwidth1 = 1050.; // when to increase binwidth.
	double switchwidth2 = 1500.;

	while( current_lowE < 2100. ){
	  if( current_lowE >= switchwidth1 && current_lowE < switchwidth2 &&switchwidth1 < 2100. ){
	    binwidth = 150.;
	  }	
	  if( current_lowE >= switchwidth2 && switchwidth2 < 2100. ){
	    binwidth = 300.;
	  }	
	  current_lowE += binwidth;
	  binEdges.push_back( current_lowE );

          if( current_lowE < 2100.){ Ebins++; }
	  cout << Ebins << "\t" << current_lowE << endl;
	}
//        double *Ebins_var = new double[Ebins];
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
//	  binwidth = .7 * current_lowE;
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
	//-- We want bin width 20% of high bin edge for EI. //
	//-------------------------------------------// 
	double current_lowEI = 0.;
        double current_highEI = 0.;
	vector<double> binEdges_EI;
	binEdges_EI.push_back( 0. );
	binwidth = 0.01;
	int EIbins = 1;

	while( current_highEI < 2100. ){
          // Binwidth = 25% bin_low_edge;
          current_highEI = current_lowEI + binwidth;
	  current_lowEI += binwidth;
	  binEdges_EI.push_back( current_lowEI );
	  binwidth = binwidth * 3.;
	}

	double* EIbins_var = &binEdges_EI[0];

	for(int i = 0; i < binEdges_EI.size(); i++){
	  EIbins_var[i] = binEdges_EI[i];
	}

	EIbins = binEdges_EI.size() - 1;


 	//------------------//
	//-- Bins for eta --//
	//------------------//
	double eta_min = -6.6;
	double eta_max = -5.2; 
	vector<double> binEdges_eta;
	binEdges_eta.push_back( eta_min );
	double binwidth_eta = (6.6-5.2)/20.;
	double current_eta = eta_min;

	while( current_eta < eta_max ){
	  current_eta += binwidth_eta;
	  binEdges_eta.push_back( current_eta );
	}
	double* etabins_var = &binEdges_eta[0];

	for(int i = 0; i < binEdges_eta.size(); i++){
	  etabins_var[i] = binEdges_eta[i];
	}

	int etabins = binEdges_eta.size() - 1;

	//-------------------------------------------//
	//-- We want bin width 20% of high bin edge for EI relative. //
	//-------------------------------------------// 
	double current_lowEI_rel = 0.;
        double current_highEI_rel = 0.;
	vector<double> binEdges_EI_rel;
	binEdges_EI_rel.push_back( 0. );
	binwidth = 0.01;
	int EIbins_rel = 1;

	while( current_highEI_rel < 100. ){
          // Binwidth = 25% bin_low_edge;
          current_highEI_rel = current_lowEI_rel + binwidth;
	  current_lowEI_rel += binwidth;
	  binEdges_EI_rel.push_back( current_lowEI_rel );
	  binwidth = binwidth * 3.;
	}

	double* EIbins_var_rel = &binEdges_EI_rel[0];

	for(int i = 0; i < binEdges_EI_rel.size(); i++){
	  EIbins_var_rel[i] = binEdges_EI_rel[i];
	}

	EIbins_rel = binEdges_EI_rel.size() - 1;


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

	//-- Throw out all that has to do with eflow.
	/*
        cout << "Prepare eflow channels" << endl;

	TH1D *hCASTOReflow_channel[224];
	cout << "It's this bitch";
	for (int i=0;i<224;i++) {
		TString name = TString::Format("hCASTOReflow_channel_%i",i+1 );
		TString title = TString::Format( "CASTOR Energy distribution for channel %i",i+1 );

		cout << "Title\tname\t" << title << "\t" << name << "\t";

		hCASTOReflow_channel[i] = new TH1D(name,title,100,-10,1800);
		hCASTOReflow_channel[i]->Sumw2();
	}
	cout << endl;
	*/

	cout << "hCastorTowerMulti" << endl;
	TH1D *hCASTORTowerMulti = new TH1D("hCASTORTowerMulti","CASTOR Tower Multiplicity (N above threshold) distribution",17,0.,17.);
	hCASTORTowerMulti->Sumw2();

	cout << "Castor tower multi" << endl;
	
	// default energy flow histos - using 5 modules
	TH1D *hCASTOReflow = new TH1D("hCASTOReflow","Total CASTOR energy flow in first 5 modules",252,-30,3750);
	hCASTOReflow->Sumw2();

	cout << "Castor e flow" << endl;
	
	TH2D *h2CASTOReflow_grid = new TH2D("h2CASTOReflow_grid","CASTOR energy weighted module vs sector distribution",16,1,17,14,1,15);
	
	cout << "Halfway" << endl;



	
	// Castor jet histograms
//        TH1D *hCastorJet_energy = new TH1D("hCastorJet_energy","CastorJet energy distribution",	Ebins, Emin, Emax);		hCastorJet_energy->Sumw2();
//        TH1D *hCastorJet_energy_lead = new TH1D("hCastorJet_energy_lead","CastorJet energy distribution", Ebins, Emin, Emax);	hCastorJet_energy_lead->Sumw2();

	//-- Variable bin width;
cout << "sectors" << endl;
	vector<double> sector_vector;
	for(int isector=0; isector<=16; isector++){
	  sector_vector.push_back( isector );
	}
	int sector_vector_size = sector_vector.size() - 1;
	double* sector_array = &sector_vector[0];
	for(int isector=0; isector<=16; isector++){
	  sector_array[isector] = isector;
	} 

        TH1D *hCastorJet_energy = new TH1D("hCastorJet_energy","CastorJet energy distribution",	Ebins, Ebins_var);
		hCastorJet_energy->Sumw2();
        TH1D *hCastorJet_energy_unCalib = new TH1D("hCastorJet_energy_unCalib","CastorJet energy distribution",	Ebins, Ebins_var);
		hCastorJet_energy_unCalib->Sumw2();

        TH2D *hCastorJet_energy_sectors = new TH2D("hCastorJet_energy_sectors","CastorJet energy distribution per sector", Ebins, Ebins_var, 16, sector_array);
		hCastorJet_energy_sectors->Sumw2();

cout << "hist created " << 		hCastorJet_energy_sectors->GetNbinsY() <<  endl;

        TH1D *hCastorJet_energy_JES_up = new TH1D("hCastorJet_energy_JES_up","CastorJet energy distribution plus JES",	Ebins, Ebins_var);
        TH1D *hCastorJet_energy_JES_down = new TH1D("hCastorJet_energy_JES_down","CastorJet energy distribution minus JES",	Ebins, Ebins_var);

        TH1D *hCastorJet_energy_lead = new TH1D("hCastorJet_energy_lead","CastorJet energy distribution", Ebins, Ebins_var);
		hCastorJet_energy_lead->Sumw2();

        TH1D *hCastorJet_energy_fine = new TH1D("hCastorJet_energy_fine","CastorJet energy distribution", Ebins*10, Emin, Emax);hCastorJet_energy_fine->Sumw2();

        cout << "Done the hCastorJet_energy" << endl;

        TH1D *hGenJet_energy = new TH1D("hGenJet_energy","GenJet energy distribution",		Ebins, Ebins_var); //Emin, Emax);
	TH1D *hGenJet_energy_lead = new TH1D("hGenJet_energy_lead","Leading GenJet energy distribution",          Ebins, Ebins_var); // Emin, Emax);
        TH1D *hGenJet_eta = new TH1D("hGenJet_eta","GenJet eta distribution;#eta;;",		20, -7.1, -4.7);
	TH1D *hCastorJet_pt = new TH1D("hCastorJet_pt","CastorJet pt distribution",30,0,30);
	TH1D *hCastorJet_em = new TH1D("hCastorJet_em","CastorJet EM energy distribution",150,0,1500);
	TH1D *hCastorJet_had = new TH1D("hCastorJet_had","CastorJet HAD energy distribution",150,0,1500);
	TH1D *hCastorJet_fem = new TH1D("hCastorJet_fem","CastorJet EM/(EM+HAD) distribution",60,-0.1,1.1);
	TH1D *hCastorJet_fhot = new TH1D("hCastorJet_fhot","CastorJet Fhot distribution",60,-0.1,1.1);
	TH1D *hCastorJet_width = new TH1D("hCastorJet_width","CastorJet width distribution",100,0,1);
	TH1D *hCastorJet_depth = new TH1D("hCastorJet_depth","CastorJet depth distribution",100,-16000,-14000);
	TH1D *hCastorJet_sigmaz = new TH1D("hCastorJet_sigmaz","CastorJet sigmaz distribution",100,0,500);
	TH1D *hCastorJet_ntower = new TH1D("hCastorJet_ntower","CastorJet ntower distribution",16,1,17);
	TH1D *hCastorJet_eta = new TH1D("hCastorJet_eta","CastorJet eta distribution",100,-6.0,-5.8);
	TH1D *hCastorJet_phi = new TH1D("hCastorJet_phi","CastorJet phi distribution",16,-M_PI,+M_PI);
	TH1D *hCastorJet_multi = new TH1D("hCastorJet_multi","CastorJet multiplicity distribution",17,0,17);

	TH2D* hJetTypes = new TH2D("hJetTypes", "Types of matched leading jet;Type_{det};Type_{gen};", 3, 0., 3., 3, 0., 3.);

	// All Casto jets response matrix.
	TH2D *hCastorJet_energy_response = new TH2D("hCastorJet_energy_response","CastorJet energy distribution - Response;E_{det};E_{gen};", 	Ebins, Emin, Emax, Ebins, Emin, Emax);	// ,Ebins+10,Emin-10*EbinWidth,Emax,Ebins+10,Emin-10*EbinWidth,Emax);
	TH2D *hCastorJet_energy_response_lead = new TH2D("hCastorJet_energy_response_lead","CastorJet energy distribution - Response;E_{det};E_{gen};", 	Ebins, Emin, Emax, Ebins, Emin, Emax);	// ,Ebins+10,Emin-10*EbinWidth,Emax,Ebin
	TH2D *hCastorJet_energy_response_fine = new TH2D("hCastorJet_energy_response_fine","CastorJet energy distribution - Response;E_{det};E_{gen};", 	Ebins*10, Emin, Emax, Ebins*10, Emin, Emax);

	TH2D *hCastorJet_energy_fakes 	 = new TH2D("hCastorJet_energy_fakes",	 "CastorJet energy distribution - Fakes", 	Ebins, Emin, Emax, Ebins, Emin, Emax);	// ,Ebins+10,Emin-10*EbinWidth,Emax);
	TH2D *hCastorJet_energy_misses	 = new TH2D("hCastorJet_energy_misses",	 "CastorJet energy distribution - Misses", 	Ebins, Emin, Emax, Ebins, Emin, Emax);	// ,Ebins+10,Emin-10*EbinWidth,Emax);

	TH1D *hCastorJet_miss_all   = new TH1D("hCastorJet_miss_all", "Misses", Ebins, Ebins_var);
        TH1D *hCastorJet_fake_all   = new TH1D("hCastorJet_fake_all", "Fakes", Ebins, Ebins_var);
        TH1D *hCastorJet_miss_lead   = new TH1D("hCastorJet_miss_lead", "Misses", Ebins, Ebins_var);
        TH1D *hCastorJet_fake_lead   = new TH1D("hCastorJet_fake_lead", "Fakes", Ebins, Ebins_var);

	TH1D *hCastorJet_matchedGen_all   = new TH1D("hCastorJet_matchedGen_all", "MatchedGen", Ebins, Emin, Emax);
		
	TH2D *hCastorJet_energy_complete_response = new TH2D("hCastorJet_energy_complete_response","CastorJet energy distribution - Complete Response", Ebins_fix, Emin_fix, Emax_fix, Ebins_fix, Emin_fix, Emax_fix);

	int Bins_response_matrix[7]	= {Ebins*10, 	Ebins*10, 16,  16,  32,  16,  32};
	double Min_response_matrix[7]	= {Emin, 	Emin, 	  0.,  0.,  0.,  0.,  0.};
	double Max_response_matrix[7]	= {Emax, 	Emax, 	  16., 16., 32., 16., 32.};

	THnSparseD *hResponse_leading 	= new THnSparseD("hResponse_leading", 	"Multidimensional Response Matrix;E_{det};E_{gen};N_{match};N_{fake};N_{miss};N_{det};N_{gen};",	7, Bins_response_matrix, Min_response_matrix, Max_response_matrix); hResponse_leading->Sumw2();
//	THnSparseD *hResponse_inclusive = new THnSparseD("hResponse_inclusive", "E_{det};E_{gen};N_{fake};N_{miss}",    4, Bins_response_matrix, Min_response_matrix, Max_response_matrix); hResponse_inclusive->Sumw2();		

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

	// Generator level energy distributions versus eta.
	TH2D* hGenJet_energy_vs_eta  = new TH2D("hGenJet_energy_vs_eta", "E_{gen}  vs. #eta;E_{gen} (GeV);#eta_{gen};", Ebins-1, Ebins_var, 24, -7.1,-4.7);

	// Efficiency control.
	TH1D *hJER = new TH1D("hJER", 	"Jet Energy Resolution", 200, -1, 5);
	TH1D *hJER_had = new TH1D("hJER_had", 	"Jet Energy Resolution", 200, -1, 5);
	TH1D *hJER_em = new TH1D("hJER_em", 	"Jet Energy Resolution", 200, -1, 5);
	TH1D *hJER_other = new TH1D("hJER_other", 	"Jet Energy Resolution", 200, -1, 5);

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

	TH2D *hJER_per_distance = new TH2D("hJER_per_distance", "#DeltaE/E for distance;#DeltaR;#DeltaE/E", 	200, 0., 6.5, 200, -5, 5);
	TH2D *hJER_per_eta = new TH2D("hJER_per_eta", "#DeltaE/E for distance;#DeltaR;#DeltaE/E",               14,-6.6,-5.2, 200, -5, 5);
	TH2D *hEnergy_per_eta = new TH2D("hEnergy_per_eta", "E_{gen} vs. #eta of leading jet", Ebins, Emin, Emax,  14,-6.6,-5.2);


	// Matches.
	TH1D *hMatched = new TH1D("hMatched", "Castor and Gen jets match", 10, -0.5, 9.5);
	TH1D *hUnmatched = new TH1D("hUnmatched", "Castor and Gen jets don't match", Ebins, Emin, Emax);
	TH1D *hJRE = new TH1D("hJRE", "Jet Reconstruction Efficiency", Ebins, Emin, Emax);

		// Towers vs. phi
	TH2D *hTower_phi = new TH2D("hTower_phi", "N_{towers} vs. #Delta#varphi;N_{towers};#Delta#varphi", 6,-0.5,5.5, 50,-3.15,3.15);

	TH2D *hIsolationEnergy = new TH2D("hIsolationEnergy", "Energy in vicinity;Towers;E_{gen};", 6, -0.5, 5.5, 100, -1000., 1800.);
	TH1D *hIsolationEnergy_1D = new TH1D("hIsolationEnergy_1D", "Isolation energy;E_{s}", 10000, 0., 1800.);
	TH1D *hIsolatedEnergy 	= new TH1D("hIsolatedEnergy", "Isolated energy;E_{s}", 10000, 0., 1800.);
	TH1D *hIsolationEnergy_Egen = new TH1D("hIsolationEnergy_Egen", "Isolation energy;E_{s}/E_{gen}", 10000, 0., 100.);
	TH1D *hIsolationEnergy_Edet = new TH1D("hIsolationEnergy_Edet", "Isolation energy;E_{s}/E_{det}", 10000, 0., 100.);

	TH1D *hIsolationEnergy_1D_log = new TH1D("hIsolationEnergy_1D_log", "Isolation energy;E_{s}", EIbins, EIbins_var);
	TH1D *hIsolatedEnergy_log 	= new TH1D("hIsolatedEnergy_log", "Isolated energy;E_{s}", EIbins, EIbins_var);
	TH1D *hIsolationEnergy_Egen_log = new TH1D("hIsolationEnergy_Egen_log", "Isolation energy;E_{s}/E_{gen}",EIbins_rel, EIbins_var_rel);
	TH1D *hIsolationEnergy_Edet_log = new TH1D("hIsolationEnergy_Edet_log", "Isolation energy;E_{s}/E_{det}", EIbins_rel, EIbins_var_rel);


	// CastorJetID distributions.
	TH1D* h_ehad = new TH1D( "ehad",	"Hadronic energy; E_{had.};",  	Ebins, Ebins_var);
	TH1D* h_eem  = new TH1D( "eem", 	"E.M. energy; E_{had.};",   	Ebins, Ebins_var);
	TH1D* h_nTowers	  = new TH1D("nTowers", 	";N_{tower};", 		4, 0.5, 4.5); 
	TH1D* h_sigmaz	  = new TH1D("sigmaz", 	"#sigma_{z};#sigma_{z};", 	100, 0., 700.  );
	TH1D* h_width	  = new TH1D( "width", 	"width;width;", 		100, 0., 0.35 );
	TH1D* h_depth	  = new TH1D( "depth", 	"depth;<z>;", 			100, -16000., -14000. );
	TH1D* h_fhot	  = new TH1D( "fhot", 	"fhot;E_{hot}/E_{tot};",	100, 0., 1. );
	TH1D* h_fem	  = new TH1D( "fem", 	"fem;E_{em}/E_{tot};",  	100, 0., 1. );
	TH1D* h_phi	  = new TH1D( "phi",	"#phi;#phi;", 			16, -3.2, 3.2);


		// RooUnfold.
		cout << "RooUnfold" << endl;

	RooUnfoldResponse response 		(hCastorJet_energy, hCastorJet_energy, "response");
	RooUnfoldResponse response_lead		(hCastorJet_energy, hCastorJet_energy, "response_lead");

	RooUnfoldResponse response_fine 	(hCastorJet_energy_fine, hCastorJet_energy_fine, "response_fine");
	RooUnfoldResponse response_lead_fine	(hCastorJet_energy_fine, hCastorJet_energy_fine, "response_lead_fine");

	RooUnfoldResponse response_match	(hCastorJet_energy, hCastorJet_energy, "response_match");

	cout << "\n\nPreparation of RooUnfold" << endl;
	std::map<int, RooUnfoldResponse* > sector_response;
	cout << "Preparation of RooUnfold" << endl;
	RooUnfoldResponse response_sector1 	(hCastorJet_energy, hCastorJet_energy, "response_sector1");		sector_response[1] = &response_sector1;
	RooUnfoldResponse response_sector2 	(hCastorJet_energy, hCastorJet_energy, "response_sector2");	sector_response[2] = &response_sector2;	
	RooUnfoldResponse response_sector3 	(hCastorJet_energy, hCastorJet_energy, "response_sector3");	sector_response[3] = &response_sector3;	
	RooUnfoldResponse response_sector4 	(hCastorJet_energy, hCastorJet_energy, "response_sector4");	sector_response[4] = &response_sector4;	
	RooUnfoldResponse response_sector5 	(hCastorJet_energy, hCastorJet_energy, "response_sector5");	sector_response[5] = &response_sector5;	
	RooUnfoldResponse response_sector6 	(hCastorJet_energy, hCastorJet_energy, "response_sector6");	sector_response[6] = &response_sector6;	
	RooUnfoldResponse response_sector7 	(hCastorJet_energy, hCastorJet_energy, "response_sector7");	sector_response[7] = &response_sector7;	
	RooUnfoldResponse response_sector8 	(hCastorJet_energy, hCastorJet_energy, "response_sector8");	sector_response[8] = &response_sector8;	
	RooUnfoldResponse response_sector9 	(hCastorJet_energy, hCastorJet_energy, "response_sector9");	sector_response[9] = &response_sector9;	
	RooUnfoldResponse response_sector10 	(hCastorJet_energy, hCastorJet_energy, "response_sector10");	sector_response[10] = &response_sector10;	
	RooUnfoldResponse response_sector11 	(hCastorJet_energy, hCastorJet_energy, "response_sector11");	sector_response[11] = &response_sector11;	
	RooUnfoldResponse response_sector12 	(hCastorJet_energy, hCastorJet_energy, "response_sector12");	sector_response[12] = &response_sector12;	
	RooUnfoldResponse response_sector13 	(hCastorJet_energy, hCastorJet_energy, "response_sector13");	sector_response[13] = &response_sector13;	
	RooUnfoldResponse response_sector14 	(hCastorJet_energy, hCastorJet_energy, "response_sector14");	sector_response[14] = &response_sector14;	
	RooUnfoldResponse response_sector15 	(hCastorJet_energy, hCastorJet_energy, "response_sector15");	sector_response[15] = &response_sector15;	
	RooUnfoldResponse response_sector0 	(hCastorJet_energy, hCastorJet_energy, "response_sector0");	sector_response[0] = &response_sector0;	

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
	TH2D *hResponse 	= new TH2D("hResponse", 	"E_{det}/E_{gen};E_{det};#frac{E_{det}}{E_{gen}}", 	Ebins, Emin, Emax, 100, 0., 5.);			hResponse->Sumw2();
	TH2D *hResponse_gen 	= new TH2D("hResponse_gen", 	"E_{det}/E_{gen};E_{gen};#frac{E_{det}}{E_{gen}}", 	Ebins, Emin, Emax, 100, 0., 5.);			hResponse_gen->Sumw2();
	TH2D *hEgen_pt 		= new TH2D("hEgen_pt",		"E_{gen};p_{T};", 					Ebins, Ebins_var, 100, 0., 25.);			hResponse_gen->Sumw2();

		// Needed for average eta per Egen.
	TH3D *hEgen_Edet_eta 	= new TH3D("hEgen_Edet_eta", 	";E_{gen} [GeV];E_{det} [GeV];#eta", 			Ebins, Ebins_var, Ebins, Ebins_var, etabins, etabins_var);			hEgen_Edet_eta->Sumw2();
	TH2D *hEgen_eta 	= new TH2D("hEgen_eta", 	";E_{gen} [GeV];#eta", 			Ebins, Ebins_var, etabins, etabins_var);			hEgen_eta->Sumw2();
	TH2D *hEta_pt 		= new TH2D("hEta_pt", 		";#eta;p_{T}", 		etabins, etabins_var, 100, 0., 25.);			hEgen_eta->Sumw2();


	int Bins_response[3]	= {Ebins, 100, 16};
	double Min_response[3]	= {Emin, 0., -1.*PI};
	double Max_response[3]	= {Emax, 5., PI};

	THnSparseD *hResponse_phi 	= new THnSparseD("hResponse_phi", 	"E_{det}/E_{gen};E_{det};#frac{E_{det}}{E_{gen}};#varphi",	3, Bins_response, Min_response, Max_response); hResponse_phi->Sumw2();
	THnSparseD *hResponse_gen_phi 	= new THnSparseD("hResponse_gen_phi", 	"E_{det}/E_{gen};E_{gen};#frac{E_{det}}{E_{gen}};#varphi",    3, Bins_response, Min_response, Max_response); hResponse_gen_phi->Sumw2();

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

	// Counting events.
	int castor_events = 0, castor_150GeV_events = 0, castor_300GeV_events = 0;
	int castor_events_MC = 0, castor_150GeV_events_MC = 0, castor_300GeV_events_MC = 0;
	bool castor_0GeV = false, castor_150GeV = false, castor_300GeV = false;
	bool castor_0GeV_MC = false, castor_150GeV_MC = false, castor_300GeV_MC = false;
        int jets_uncalibrated = 0, jets_calibrated = 0;

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
	
	  ////////////////////////////////
	  // Extract this event's jets. //
	  ////////////////////////////////
	
   	  if( !isData_) { b_CastorGenJets->GetEntry( counter_events ); }
	  b_CastorJets->GetEntry( counter_events );	
				
	  /////////////////////////////////////////
	  // Start Nvertex == 1 part of the code //
	  /////////////////////////////////////////






	  //------------------------------------------//  		
	  //------------------------------------------//
  	  // Get some general information on MC jets. //
	  //------------------------------------------//
	  //------------------------------------------//

  	  if( !isData_ ) { 

	    // Count the number of selected events with jets in CASTOR.
	    int ndet = 0, ngen = 0;

	    //-- Loop over detector level jets.
   	    for(int det_jet = CastorJets->size()-1; det_jet >= 0; det_jet-- ){
	      MyCastorJet castor_jet = (*CastorJets)[ det_jet ];

	      double det_energy = castor_jet.energy;           
	      double det_phi = castor_jet.phi;
	      int sector = CastorSector( det_phi ) ; 
 	      det_energy = CalibratedDet( det_energy, sector, fileLabel_, threshold_ );

	      if( det_energy> threshold_ ){ ndet++; }
	    } // Loop over det. jets.

	    //-- Loop over CASTOR (gen) jets.
	    for( int gen_jet =  CastorGenJets->size()-1; gen_jet >=0; gen_jet--){
	      MyGenJet castor_gen = (*CastorGenJets)[gen_jet];
	      double gen_energy = castor_gen.Energy();
  	      castor_0GeV_MC = true;

	      //-- Place various thresholds.
	      if( gen_energy > threshold_ ){
		double gen_eta = castor_gen.Eta();
		hGenJet_eta->Fill( gen_eta );

		if( gen_eta > genetamin_ && gen_eta < genetamax_ ){ 
    	          hGenJet_energy->Fill( gen_energy );
		  hGenJet_energy_vs_eta->Fill( gen_energy, gen_eta );
		  ngen++;
                  hEgen_pt->Fill( gen_energy, castor_gen.Pt() );
		  hEta_pt->Fill( gen_eta, castor_gen.Pt());
		}
		
	    	if( gen_energy > 150. ){
	      	  castor_150GeV_MC = true;
	      	  if( gen_energy > 300. ){
		    castor_300GeV_MC = true;
	          }	   
	        }
	      }
	    }// Loop over gen. jets.

 	    // Close loop over detector level jets.
	    if( castor_0GeV_MC ){ castor_events_MC++; }
	    if( castor_150GeV_MC ){ castor_150GeV_events_MC++; }
	    if( castor_300GeV_MC ){ castor_300GeV_events_MC++; }

	    castor_0GeV_MC = false;
	    castor_150GeV_MC = false;
	    castor_300GeV_MC = false;

	    hNumber_of_match_jets->Fill( ndet, ngen);  	    
                     
	    if( ((*CastorGenJets)[0]).Energy() > threshold_ ){ hGenJet_energy_lead->Fill( ((*CastorGenJets)[0]).Energy() ); }
	    hCastorJet_pt->Fill( ((*CastorGenJets)[0]).Pt() );

	  }
	  
	  //------------------------------------------------------//  		
	  //------------------------------------------------------//
  	  // Get some general information on detector level jets. //
	  //------------------------------------------------------//
	  //------------------------------------------------------//

	  if( comments_) cout << "\n\n\n\t" << counter_events << endl;
	  if( CastorJets->size() > 0 ){ castor_0GeV = true; }

	  for(int det_jet = CastorJets->size()-1; det_jet >= 0; det_jet-- ){
	    MyCastorJet castor_jet = (*CastorJets)[ det_jet ];
	    double det_energy = castor_jet.energy;

	    if( do_calibration_function){
	      double det_phi = castor_jet.phi;
	      int sector = CastorSector( det_phi ) ; 
 	      det_energy = CalibratedDet( det_energy, sector, fileLabel_, threshold_ );
	    }

	    if( det_energy > 150. ){
	      castor_150GeV = true;
	      if( det_energy > 300. ){
		castor_300GeV = true;
	      }
	    }

	    if( det_jet == 0 ){
  	      if( det_energy < threshold_ ){ emptyCastor++; break;}		
	      if( det_energy > threshold_ ){ hCastorJet_eta->Fill( castor_jet.eta ); }
	    }

	    if( comments_ ){ cout << "Type is\t" << GetJetType( castor_jet ) << "\tfor\t" << counter_events << "\t" << det_jet << "\t" << det_energy << endl; }
	  } // Close loop over detector level jets.
	  if( castor_0GeV ){ castor_events++; }
	  if( castor_150GeV ){ castor_150GeV_events++; }
	  if( castor_300GeV ){ castor_300GeV_events++; }

	  castor_0GeV = false;
	  castor_150GeV = false;
	  castor_300GeV = false;



	

	  if( comments_ && !isData_){ cout << "\n\n\n$$$\tEvent\t" << counter_events << "\tdet size\t" << CastorJets->size() << "\tgen size\t" << CastorGenJets->size() << endl;}




	  //---------------------//  		
	  //---------------------//
  	  // Matching procedure. //
	  //---------------------//
	  //---------------------//

  	  int matched_pairs = 0;
		  
  	  while( matched_pairs == 0 && CastorJets->size() > 0 ){ 	// Continue as long as there are Castor jets and the switched matched_pairs has not been turned off.
	    if( comments_ ){ cout << "First pair to match" << endl; }	

	    ////////////////////////////////////////////////////
	    // Get information on leading jet (type, energy). //
	    ////////////////////////////////////////////////////
	
	    if( analyze_type_ ){
	      MyCastorJet castorjet = (*CastorJets)[ 0 ];
  	      TString detjettype = GetJetType( castorjet );	      
	    }

	    //-- Match only the leading jet: determine jet type and number of sectors.
	    //-- (TO DO) The cut on the number of sectors is a relic and needs to be removed.
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
	      double det_phi = castor_det.phi;



 
	      //== Energy

	      int sector = CastorSector( det_phi ) ;
	      if( det_energy > threshold_ ) {

		//== Energy
		hCastorJet_energy_unCalib->Fill( det_energy ); 
		jets_uncalibrated++;
	      }	      

	      if( do_calibration_function && detjettype == "had"){
	        double det_phi = castor_det.phi;
	        int sector = CastorSector( det_phi ) ; 
 	        det_energy = CalibratedDet( det_energy, sector, fileLabel_, threshold_ );
	      }
	      if( det_energy > threshold_ ) {

	        //== CastorJetID plots.
	        h_ehad->Fill( castor_det.ehad );
	        h_eem->Fill( castor_det.eem );
	        h_nTowers	->Fill( castor_det.ntower );
	        h_sigmaz	->Fill( castor_det.sigmaz );
	        h_width	->Fill( castor_det.width );
	        h_depth	->Fill( castor_det.depth );
	        h_fhot	->Fill( castor_det.fhot );
	        h_fem	->Fill( castor_det.fem );
	        h_phi	->Fill( det_phi );


	       jets_calibrated++;
 	       hCastorJet_energy->Fill( det_energy ); // Inclusive jet energy spectrum.
 	       hCastorJet_energy_sectors->Fill( det_energy, sector ); // Inclusive jet energy spectrum.
	       if( det_jet == 0)  hCastorJet_energy_lead->Fill( det_energy ); // Leading jet energy spectrum.
	      }

	      //-- The following fills the detector level jet energy distributions with a shift according to JES.
	      //-- energy - JES
	      if( det_energy * (1. - JES_unc) > threshold_ ){
 	       hCastorJet_energy_JES_down->Fill( det_energy*(1. - JES_unc)  );
	      }
	      //-- energy + JES
	      if( det_energy * (1. - JES_unc)  > threshold_ ){
 	       hCastorJet_energy_JES_up->Fill( det_energy*(1. + JES_unc)  );
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

	      // Count the number of Gen. jets above the energy threshold.
	      for(int gen_jet = 0; gen_jet < CastorGenJets->size(); gen_jet++){ 
		MyGenJet castor_gen = (*CastorGenJets)[gen_jet];
		if( castor_gen.Energy() > threshold_ ){ nGen++; }
	      }

	      double leading_det_E, leading_gen_E;

   	      for(int det_jet = 0; det_jet < CastorJets->size(); det_jet++){

       	        MyCastorJet castor_det = (*CastorJets)[det_jet];

		//-- CALIBRATION: Jet energy is calibrated before taking any further action.
                double det_energy = castor_det.energy;
		double det_phi = castor_det.phi; 
	        int sector = CastorSector( det_phi ) ;
		TString jettype = GetJetType( castor_det );
	
		//Calibration is only needed for non-em jets.
		if( jettype != "em" ){ det_energy = CalibratedDet( det_energy, sector, fileLabel_, threshold_ ); }

		if( comments_ ){ cout << counter_events << "\tDet\t" << det_jet << " of\t" << CastorJets->size() << "\t" << det_energy << "\t" << det_phi << "\t" << GetJetType(castor_det) << endl; }

		// CUT: cut jet if energy is below threshold.
		if( det_energy < threshold_ ){ 
                  if( comments_ ){ cout << "\t\t\t\t\t\t\t\tLow Det E\t" << det_jet << "\t" << det_energy << endl; }
		  continue; 
		} 
		nDet++;
		
    	        // MATCHING
		//-- Match closest in phi by looping over all jets and find one with lowest delta phi.
		//-- Match gen. jet with highest E in interval Delta phi around Castor jet.
    	        double min_delta_phi = deltaphimax_ ;
		if( sectors_ == 1 ){ min_delta_phi = 0.1; }

    	        int match_gen = -1;		
		double max_energy = 0.;			  

    	        for(int gen_jet =  CastorGenJets->size()-1; gen_jet >=0; gen_jet--){

    	          MyGenJet castor_gen = (*CastorGenJets)[gen_jet];
		  double gen_energy = castor_gen.Energy();
		  double gen_eta = castor_gen.Eta();

		  hEgen_eta->Fill( gen_energy, gen_eta );


		  // Remove and skip gen jet if 
		  //-- OR gen. jet not hard enough 
		  //-- OR gen. jet outside of etaband.
		
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
    	          } //Match Ephi


		  // Hardest jet close to det. jet: save Delta Phi, energy and index of the jet.
		  if( match_ == "matchE" ){		
    	            if( gen_energy > max_energy && delta_phi < deltaphimax_ ){
		      if( comments_  ){ cout << "\t\t\t\t\t\t\t\t\t\t\t\t\tNew E maximum" << endl; }
    	              max_energy = gen_energy;
    	              match_gen = gen_jet;
    	            } // New maximum.
		    if( comments_ ){ cout << "\t\t\t\t\t\t\t\t\t\t\t" << delta_phi << "\t" << min_delta_phi << endl; }
    	          } //Match E.

		}// Loop over possible match candidates (gen).




    	        // We have a matching gen. jet
    	        if( match_gen != -1 ){
    	          MyGenJet castor_gen = (*CastorGenJets)[match_gen];
    	          double gen_energy = castor_gen.Energy();	
		  double gen_eta = castor_gen.Eta();

		  //-- If this is a leading detector level jet, store the information in the leading jet histograms as well.
		  if( det_jet == 0){
		    leading_det = (*CastorJets)[0];
		    leading_gen = (*CastorGenJets)[match_gen];
		    leading_gen_ = match_gen;

		    response_lead.Fill( det_energy, gen_energy );
		    response_lead_fine.Fill( det_energy, gen_energy );
		    hCastorJet_energy_response_lead->Fill( det_energy, gen_energy );
		    leading_gen_E = gen_energy; 
		    leading_det_E = det_energy;
		  }

		  // Feed information to response matrix and histograms.	  

    	          response.Fill( det_energy, gen_energy );	// RooUnfoldResponse
    	          response_match.Fill( det_energy, gen_energy );// RooUnfoldResponse (only matches)
		  response_fine.Fill( det_energy, gen_energy ); // RooUnfoldResponse with fine binning for calibration.
		  hCastorJet_energy_response->Fill( det_energy, gen_energy );	// 2D histogram
                  hCastorJet_matchedGen_all->Fill( gen_energy );	
		  (sector_response[sector])->Fill( det_energy, gen_energy );	// Sector dependent response.
 
		  hEgen_Edet_eta->Fill( gen_energy, det_energy, gen_eta);

		  nMatch++;

	          if( comments_ ){ cout << "\t\t\t\t\tMatch" << endl << endl; }

		  // Erase the matched gen. jet from the list of available jets.
    	          CastorGenJets->erase( CastorGenJets->begin() + match_gen ); 
    	        } // Matching gen jet.



		//-----------------------------------------------------------//
    	        // FAKES: We don't have a matching gen. jet. This is a fake. //
		//-----------------------------------------------------------//

    	        else{
    	            response.Fake( det_energy );
		    hCastorJet_fake_all->Fill( det_energy );
		    response_fine.Fake( det_energy ); 
		    (sector_response[sector])->Fake( det_energy );

		    if( det_jet == 0 ){ 
		      hCastorJet_fake_lead->Fill( det_energy ); 
		      response_lead.Fake( det_energy ); 
		     response_lead_fine.Fake( det_energy ); 
		    }
		    if( comments_ ){ cout << counter_events << "\t\t\t\t\tFake" << endl << endl; }
		    nFake++;

		    if( det_energy > 1400. ){ cout << "Absurd fake at\t" << counter_events << endl; }
    	        } // FAKES.
    	      }// Loop over detector level jets.
 

	      //----------------------------------------------------------------------------------//
              // MISSES: We had all det. jets, time for the remaining gen. jets which are misses. //
	      //----------------------------------------------------------------------------------//

    	      if(  CastorGenJets->size() > 0 ){

    	        for(int gen_jet = 0; gen_jet < CastorGenJets->size(); gen_jet++){

    	          MyGenJet castor_gen = (*CastorGenJets)[gen_jet];
    	          double gen_energy = castor_gen.Energy();
		  double gen_eta = castor_gen.Eta();

	          double gen_phi = castor_gen.Phi();
	          int sector = CastorSector( gen_phi ) ; 

		  // Our gen. jets have to lie in our eta band AND must pass an energy cut.
		  if( gen_energy < threshold_ || gen_eta < genetamin_ || gen_eta > genetamax_ ){ 
		    continue; 
		  }
		  response.Miss(gen_energy ); 
	 	  response_fine.Miss( gen_energy ); 
	 	  //hGenJet_energy->Fill( gen_energy );	
		  hCastorJet_miss_all->Fill( gen_energy );
		  (sector_response[sector])->Miss( gen_energy );
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
	    //------------------------------------------------//
	    // -- Match jets without filling response matrix. //
	    //----Needed for calibration ---------------------//
	    //------------------------------------------------//
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

//	      cout << "\t\t\t" << counter_events << "\t" << GetJetType( castor_jet ) << "\t" << GetJetType( leading_gen ) << endl;
              int genjet = 0;
	      IsolationCut isolation_energy( phi_det_lead, sectors_lead );
if( sectors_lead == 5) { cout << "Sectors for leading det.\t" << sectors_lead << endl; }
	      
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
	        if( E_gen_cut == 0.){ 
		  hIsolatedEnergy->Fill( 0. ); 
		  hIsolatedEnergy_log->Fill( 0. ); 
		}

	        if( E_gen_cut > 0. ){ 
		  hIsolationEnergy_1D->Fill( E_gen_cut );	
                  hIsolationEnergy_Egen->Fill( E_gen_cut / gen_energy );
                  hIsolationEnergy_Edet->Fill( E_gen_cut / det_energy_lead );	

		  hIsolationEnergy_1D_log->Fill( E_gen_cut );	
                  hIsolationEnergy_Egen_log->Fill( E_gen_cut / gen_energy );
                  hIsolationEnergy_Edet_log->Fill( E_gen_cut / det_energy_lead );				  
		}
              } // Cut on EI.

	      if( E_gen_cut > 0 && cut_EI ){ break; }

	      //////////////////////
	      // Fill histograms. //
	      //////////////////////

	      // Jet types of the matched jets.
	      int det_id = 1, gen_id = 1;
   	      TString detjettype = GetJetType( castor_jet );
	      if( GetJetType( castor_jet ) == "had" ){ det_id = 0; }
	      if( GetJetType( castor_jet ) == "em" ){ det_id = 2; }
	      if( GetJetType( leading_gen ) == "had" ){ gen_id = 0; }
	      if( GetJetType( leading_gen ) == "em" ){ gen_id = 2; }

	      hJetTypes->Fill( det_id, gen_id );	

	      // If needed, calibrate.					  
	      if( comments_){ cout << "Events\t" << counter_events << "\tActual match\t" << endl; }

	      if( 	do_calibration_discrete ){ 			det_energy_lead = CalibratedDet(lowedge, muval, det_energy_lead, sector);       	}
	      else if( 	do_calibration_function && !sector_dependence){	det_energy_lead = CalibratedDet( det_energy_lead ); 				}
              else if( 	do_calibration_function && sector_dependence){ 	det_energy_lead = CalibratedDet( det_energy_lead, sector, fileLabel_, threshold_ ); 			}					  
				  
	     if( comments_ ){  	cout << "\tFilling response matrix\t" << det_energy_lead << "\t" << gen_energy << "\tsector\t" << sector << endl; }

		hCastorJet_energy_response_lead->Fill( det_energy_lead, gen_energy);
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
	          if( detjettype == "had" ){
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
	          else if( detjettype == "em" ){
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
		  else{
		    hJER_other			->Fill( JER );	
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
	/*
	for (int icha=0;icha<224;icha++) {
		hf_.checkFlow(hCASTOReflow_channel[icha]);
	}
	*/
	
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
	
	/*
	for (int icha=0;icha<224;icha++) {
		hCASTOReflow_channel[icha]->Write();
	}
	*/

//	hJRE->Divide( hMatched, hUnmatched );
	
	hCastorJet_energy->Write();
	hCastorJet_energy_unCalib->Write();
	hCastorJet_energy_sectors->Write();
	hCastorJet_energy_JES_up->Write();
	hCastorJet_energy_JES_down->Write();
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
	hGenJet_eta->Write();

	hGenJet_energy_vs_eta->Write();


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

	hJetTypes->Write();
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
	
cout << "pi_e" << endl;
	  hPi_e_ratio->Divide( hPion_energy, hElectron_energy );
cout << "pi_e II" << endl;
	hPi_e_ratio		->Write();

	hNumber_of_match_jets	->Write();

	hJER			->Write();
	hJER_had		->Write();
	hJER_em			->Write();
	hJER_other		->Write();

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

	// Egen vs. pt.
	hEgen_pt->Write();
	hEta_pt->Write();


	// Test phi.
	 hPhi_gen_det->Write();
	
	hMatched->Write();
	hUnmatched->Write();
	hJRE->Write();		

	hNjet_vs_Ejets_gen->Write();
        hNjet_vs_Ejets_det->Write();

	hIsolationEnergy->Write();
	hIsolatedEnergy->Write();
	hIsolationEnergy_1D->Write();
	hIsolationEnergy_Egen->Write();
	hIsolationEnergy_Edet->Write();

	hIsolatedEnergy_log->Write();
	hIsolationEnergy_1D_log->Write();
	hIsolationEnergy_Egen_log->Write();
	hIsolationEnergy_Edet_log->Write();

	hsumEreco_sumEgen->Write();
	hsumEreco_vs_sumEgen->Write();

	hEgen_Edet_eta->Write();
	hEgen_eta->Write();

	// CastorJetID

	h_ehad->Write();
	h_eem->Write();
	h_nTowers	->Write();
	h_sigmaz	->Write();
	h_width	->Write();
	h_depth	->Write();
	h_fhot	->Write();
	h_fem	->Write();
	h_phi	->Write();


	// Save response objects to file.
	response.Write();
	response_lead.Write();
        response_fine.Write();
	response_match.Write();
	response_lead_fine.Write();

	for(int sector = 0; sector < 16; sector++){
	  (sector_response[sector])->Write();
	}
//	sector_response->Write();


	cout << castor_events << "\t" << castor_150GeV_events << "\t" << castor_300GeV_events << "\t" << endl;


	// Read in tree from the source file and add its total number of events to the ne tree.
	TTree *old_tree;
	int total_events, total_events_nocuts = 0;
	int castor_events_nocuts_branch, castor_events_nocuts = 0;
	int castor_150GeV_events_nocuts_branch, castor_150GeV_events_nocuts = 0;
	old_tree = (TTree*)currentfile_->Get("useful_numbers");

        if( old_tree ){
  	  old_tree->SetBranchAddress("total_events", &total_events );
	  old_tree->SetBranchAddress("castor_events_MC", &castor_events_nocuts_branch );
	  old_tree->SetBranchAddress("castor_150GeV_events_MC", &castor_150GeV_events_nocuts_branch );

	  for( int entry = 0; entry < old_tree->GetEntries(); entry++){
	    old_tree->GetEntry( entry );

	    total_events_nocuts += total_events;
	    castor_events_nocuts += castor_events_nocuts_branch;
	    castor_150GeV_events_nocuts += castor_150GeV_events_nocuts_branch;
          }
        }

	output->cd();

  	// Write away useful numbers.
	TTree *tree_numbers = new TTree("useful_numbers", "useful_numbers");

        if( old_tree ){
  	  tree_numbers->Branch("total_events_nocuts", &total_events_nocuts, "total_events_nocuts/I");
        }

	tree_numbers->Branch("castor_events", &castor_events, "castor_events/I");
	tree_numbers->Branch("castor_150GeV_events", &castor_150GeV_events, "castor_150GeV_events/I");
	tree_numbers->Branch("castor_300GeV_events", &castor_300GeV_events, "castor_300GeV_events/I");
	tree_numbers->Branch("Njets_uncalibrated", &jets_uncalibrated, "jets_uncalibrated/I");
	tree_numbers->Branch("Njets_calibrated", &jets_calibrated, "jets_calibrated/I");

	if( !isData_ ){
	  tree_numbers->Branch("castor_events_MC", &castor_events_MC, "castor_events_MC/I");
	  tree_numbers->Branch("castor_150GeV_events_MC", &castor_150GeV_events_MC, "castor_150GeV_events_MC/I");
	  tree_numbers->Branch("castor_300GeV_events_MC", &castor_300GeV_events_MC, "castor_300GeV_events_MC/I");

	  if (old_tree ){
	    tree_numbers->Branch("castor_events_nocuts", &castor_events_nocuts, "castor_events_nocuts/I");
	    tree_numbers->Branch("castor_150GeV_events_nocuts", &castor_150GeV_events_nocuts, "castor_150GeV_events_nocuts/I");
	  }
	}

	tree_numbers->Fill();
	tree_numbers->Write();
	if( !isData_ && old_tree){
	  // Transport the gen level histogram of all Castor jets before detector level cuts.
	  TH1D* hGenJet_energy_noCuts;
	  hGenJet_energy_noCuts = (TH1D*)currentfile_->Get("hGenJet_energy_noCuts");
	  hGenJet_energy_noCuts->Write();
	}

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




