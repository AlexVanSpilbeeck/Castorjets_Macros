#ifndef HelperFunctions_h
#define HelperFunctions_h

#include <TString.h>
#include <TRegexp.h>
#include <TObjArray.h>
#include <TH1.h>


class HelperFunctions {
public:
	HelperFunctions();
	virtual ~HelperFunctions();
    
    // statistical
	Double_t getRatioError(TH1D * hMB, TH1D * hQCD);	
	Double_t getRatioError(double a, double b, double errora, double errorb);
	Double_t getMultiError(double a, double b, double errora, double errorb);
	
	// histogram operations
    void checkFlow(TH1D *histo);
    TH1D getHistoByName(std::vector<TH1D*> allhistograms,TString name);
	
    // geometrical
    Double_t deltaPhi2(double phi1, double phi2);
	
	// CASTOR specific
    int rad2Sector(double radians);
    int AddSector(int sector,int add);
    int SubtractSector(int sector,int sub);
    
	
private:
	
	
};

#endif
