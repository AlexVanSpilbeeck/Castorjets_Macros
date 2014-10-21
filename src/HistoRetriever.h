#ifndef HistoRetriever_h
#define HistoRetriever_h

#include <TString.h>
#include <TObjArray.h>
#include <TH1D.h>
#include <TH2D.h>

class HistoRetriever {
	public:
        HistoRetriever();
	virtual ~HistoRetriever();
        std::vector<TH1D*> getHistos(TString inputfile);
	std::vector<TH2D*> get2DHistos(TString inputfile);

	private:
};

#endif
