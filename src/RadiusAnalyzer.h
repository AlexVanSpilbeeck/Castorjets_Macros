#ifndef RadiusAnalyzer_h
#define RadiusAnalyzer_h

#include <TString.h>
#include <TRegexp.h>
#include <TObjArray.h>
#include <TH1.h>
#include <TFile.h>

#include "../src/MyJet.h"
#include "../src/MyGenJet.h"
#include "../src/MyTrackJet.h"

#include "HelperFunctions.h"

class RadiusAnalyzer {
public:
	RadiusAnalyzer(TString inputdir, TObjArray* filelist, bool isData, const char* outputname);
	virtual ~RadiusAnalyzer();
    
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
	static void* OpenROOTFile(RadiusAnalyzer* arg);
    
    // analysis specific helper functions
    int posLeadingGenJet(std::vector<MyGenJet> JetVector,double etacut,double minptcut);
	int posLeadingTrackJet(std::vector<MyTrackJet> JetVector,double etacut,double minptcut);
    
	
private:
    
    // basic variables
    TString inputdir_; 
    TObjArray* filelist_;
    bool isData_; 
    const char* outputname_;
	
	TString currentfile_;
	TFile* currentTFile_;
	static TFile* currentStaticTFile_;
    
    TString LoopOutputFile_;
    
    HelperFunctions hf_;
    
    unsigned int Nruns;
    int nLumiBins;
    
    std::vector<int> runs;
    std::vector<double> LumiBins;
	
};

#endif
