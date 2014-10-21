#ifndef HadronAnalyzer_h
#define HadronAnalyzer_h

#include <TString.h>
#include <TRegexp.h>
#include <TObjArray.h>
#include <TH1.h>
#include "../src/MyGenJet.h"

class HadronAnalyzer {
public:
	HadronAnalyzer();
	virtual ~HadronAnalyzer();
	void Loop(TString inputdir, TObjArray* filelist, double cmenergy,const char* outputname);
	Double_t getRatioError(TH1D * hMB, TH1D * hQCD);	
	Double_t getRatioError(double a, double b, double errora, double errorb);
	bool isGenDiJet(std::vector<MyGenJet> JetVector, bool backtoback, double usedetacut, double usedptcut);
	bool isGenDiJet(std::vector<MyGenJet> JetVector, bool backtoback, double usedetacut, double lowptcut, double highptcut, bool mean);
	Double_t deltaPhi2(double phi1, double phi2);
	void checkFlow(TH1D *histo);
	int posLeadingGenJet(std::vector<MyGenJet> JetVector,double etacut,double minptcut);
	
private:

	
};

#endif
