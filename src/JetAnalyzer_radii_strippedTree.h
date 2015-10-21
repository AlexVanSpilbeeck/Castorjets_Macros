#ifndef JetAnalyzer_radii_strippedTree_h
#define JetAnalyzer_radii_strippedTree_h

#include <TString.h>
#include <TRegexp.h>
#include <TObjArray.h>
#include <TH1.h>
#include <TFile.h>

#include "../src/MyJet.h"
#include "../src/MyGenJet.h"
#include "../src/MyTrackJet.h"

#include "HelperFunctions.h"
//#include "IsolationCut.h"

class JetAnalyzer_radii_strippedTree {
public:
	JetAnalyzer_radii_strippedTree(TString inputdir, bool isData, const char* outputname,int totalEvents, TString date, TString filename, TString jettype, double threshold, TString setup, double deltaphimax);
	virtual ~JetAnalyzer_radii_strippedTree();
    
    // Basic functions
	void Loop();
    void AfterLoopCalculations(TString file);
    
	// getter functions
    TString getOutputFile();
	TString getInputDir();
	TString getCurrentFile();
	
	//setter functions
	void setCurrentTFile();
	
	// static functions
	static void* OpenROOTFile(JetAnalyzer_radii_strippedTree* arg);
    
    // analysis specific helper functions
    int posLeadingGenJet(std::vector<MyGenJet> JetVector,double etacut,double minptcut);
	int posLeadingTrackJet(std::vector<MyTrackJet> JetVector,double etacut,double minptcut);
    
	
private:
    
    // basic variables
    TString inputdir_; 
    bool isData_; 
    const char* outputname_;
    int totalEvents_;
    TString date_;
    TString filename_;
    TString jettype_;
    int sectors_;
    float etaband_;
    double threshold_;
    TString setup_;
    double deltaphimax_;

    bool do_calibration_function;
    bool cut_EI;
    bool prepare_unfolding;
	
	TString currentfile_;
	TFile* currentTFile_;
	static TFile* currentStaticTFile_;
    
    TString LoopOutputFile_;
    
    HelperFunctions hf_;
    
    unsigned int Nruns;
    int nLumiBins;
    
    std::vector<int> runs;
    std::vector<double> LumiBins;

    TString fileLabel_;
	
};

#endif
