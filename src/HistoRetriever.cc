#include "HistoRetriever.h"

#include "TObjString.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"
#include "TKey.h"
#include "TFile.h"
#include "TDirectory.h"

#include <iostream>

HistoRetriever::HistoRetriever() {

}

HistoRetriever::~HistoRetriever() { 

}

std::vector<TH1D*> HistoRetriever::getHistos(TString inputfile) {

	std::vector<TH1D*> histoList;
	
	TFile *file = new TFile(inputfile);
	TList *list = file->GetListOfKeys();
    //TDirectory *d = (TDirectory*)((TKey*)file->GetListOfKeys()->At(0))->ReadObj();
    //TList *list = d->GetListOfKeys();
	
	for (int i=0;i<list->GetSize();i++) {		
		TKey *key = (TKey*)list->At(i);
        //std::cout << " key class = " << key->GetClassName() << std::endl;
        
        if (strcmp(key->GetClassName(),"TH1D")==0) {
            TH1D *hist = (TH1D*)key->ReadObj();
            histoList.push_back(hist);
        }
	}	
	
	return histoList;

}

std::vector<TH2D*> HistoRetriever::get2DHistos(TString inputfile) {
    
	std::vector<TH2D*> histoList;
	
	TFile *file = new TFile(inputfile);
	TList *list = file->GetListOfKeys();
	
	for (int i=0;i<list->GetSize();i++) {		
		TKey *key = (TKey*)list->At(i);
        //std::cout << " key class = " << key->GetClassName() << std::endl;
        
        if (strcmp(key->GetClassName(),"TH2D")==0) {
            TH2D *hist = (TH2D*)key->ReadObj();
            histoList.push_back(hist);
        }
	}	
	
	return histoList;
    
}