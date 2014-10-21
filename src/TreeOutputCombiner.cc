
#include "TreeOutputCombiner.h"
#include "HistoRetriever.h"

//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
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
#include "TObjString.h"
#include "TSystem.h"
#include "TList.h"
#include "TKey.h"
#include "TString.h"
#include "TPaveText.h"
#include "TLegend.h"

//STANDARD C++ INCLUDES
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <ctime>
#include<sys/stat.h>
#include<sys/types.h>

TreeOutputCombiner::TreeOutputCombiner() { }

TreeOutputCombiner::~TreeOutputCombiner() { }


void TreeOutputCombiner::Combine(TString inputdir, TObjArray* filelist, double usedWeight) {
	
	
	using namespace std;
	
	/////////////////////////////////////
	// Define all histograms
	/////////////////////////////////////
    
    //std::vector<TH1D*> alladdedhistograms;
    std::vector<TH2D*> alladded2Dhistograms;
    
   	
	TIter       next(filelist); 
	TObjString* fn = 0;
    TIter       next2(filelist); 
	TObjString* fn2 = 0;
	TString currentfile = "";
    TString lastfile = "";
	
	//TH1::AddDirectory(kFALSE);
    
    
    unsigned int maxN1dhistos = 0;
    unsigned int maxN2dhistos = 0;
    
    // start file loop to check for number of histograms in array
	while((fn2 = (TObjString*)next2())) { 
		
		currentfile.Clear();
		currentfile = fn2->GetString();
        
        if (!currentfile.Contains("TreeOutputCombiner") && currentfile.Contains("TreeAnalyzer")) {
			
            HistoRetriever histogetter_;
            std::vector<TH1D*> histovector = histogetter_.getHistos(inputdir+currentfile);
            std::vector<TH2D*> h2Dhistovector = histogetter_.get2DHistos(inputdir+currentfile);
            
            if (maxN1dhistos < histovector.size()) maxN1dhistos = histovector.size();
            if (maxN2dhistos < h2Dhistovector.size()) maxN2dhistos = h2Dhistovector.size();
            
        }
		
	} // end file loop
    
    
    TH1D* alladdedhistograms[maxN1dhistos];
    for (unsigned int ihisto=0;ihisto<maxN1dhistos;ihisto++) {
        alladdedhistograms[ihisto] = new TH1D();
        //alladdedhistograms[ihisto]->Sumw2();
    }
    
	
	std::cout << "start looping over files" << std::endl;
    
    int filenumber = 1;
	
	// start file loop
	while((fn = (TObjString*)next())) { 
		
		currentfile.Clear();
		currentfile = fn->GetString();
        
        if (!currentfile.Contains("TreeOutputCombiner") && currentfile.Contains("TreeAnalyzer")) {
			
            HistoRetriever histogetter_;
            std::vector<TH1D*> histovector = histogetter_.getHistos(inputdir+currentfile);
            std::vector<TH2D*> h2Dhistovector = histogetter_.get2DHistos(inputdir+currentfile);
            
            unsigned int N1dhistos = histovector.size();
            //int N2dhistos = h2Dhistovector.size();
            
            TH1D* histoarray[N1dhistos];
            //TH2D* h2Dhistoarray[N2dhistos];
            
            for (unsigned int ihisto=0;ihisto<N1dhistos;ihisto++) {
                histoarray[ihisto] = new TH1D();
                //histoarray[ihisto]->Sumw2();
                histoarray[ihisto] = (TH1D*)histovector[ihisto]->Clone();
                if (strcmp(histoarray[ihisto]->GetName(),"heflow_leadingjet")==0) std::cout << "Cloned heflow_leadingjet into array with mean = " << histoarray[ihisto]->GetMean() << " and error = " << histoarray[ihisto]->GetMeanError() << std::endl;
            }
            
            
            // put all histograms of the first file in the general vector
            if (filenumber == 1) {
                double pthatweight = 1.;
                if (!currentfile.Contains("LowPtHat")) {
                    pthatweight = usedWeight;
                }
                for (unsigned int ihist=0;ihist<N1dhistos;ihist++) {
                    histoarray[ihist]->Scale(pthatweight);
                    alladdedhistograms[ihist] = histoarray[ihist];
                }
                std::cout << "pushed back " << histovector.size() << " 1D histograms from file " << currentfile << " with weight =  " << pthatweight << std::endl;
                
                for (unsigned int ihist=0;ihist<h2Dhistovector.size();ihist++) {
                    h2Dhistovector[ihist]->Scale(pthatweight);
                    alladded2Dhistograms.push_back(h2Dhistovector[ihist]);
                }
                std::cout << "pushed back " << h2Dhistovector.size() << " 2D histograms from file " << currentfile << " with weight =  " << pthatweight << std::endl;
            }
            
            // add histograms from following files 
            if (filenumber != 1) {
                double pthatweight = 1.;
                if (!currentfile.Contains("LowPtHat")) {
                    pthatweight = usedWeight;
                }
                for (unsigned int ihist=0;ihist<N1dhistos;ihist++) {
                    bool added = false;
                    for (unsigned int jhist=0;jhist<maxN1dhistos;jhist++) {
                        if (strcmp(alladdedhistograms[jhist]->GetName(),histoarray[ihist]->GetName()) == 0) {
                            alladdedhistograms[jhist]->Add(histoarray[ihist],pthatweight);
                            added = true;
                        } 
                    }
                    if (!added) std::cout << " could not add histogram name: " << histoarray[ihist]->GetName() << std::endl;
                }
                std::cout << "added " << histovector.size() << " 1D histograms from file " << currentfile << " with weight =  " << pthatweight << std::endl;
                
                
                for (unsigned int ihist=0;ihist<h2Dhistovector.size();ihist++) {
                    bool added = false;
                    for (unsigned int jhist=0;jhist<alladded2Dhistograms.size();jhist++) {
                        if (strcmp(alladded2Dhistograms[jhist]->GetName(),h2Dhistovector[ihist]->GetName()) == 0) {
                            alladded2Dhistograms[jhist]->Add(h2Dhistovector[ihist],pthatweight);
                            added = true;
                        } else {
                            //std::cout << " error with histogram name: " << alladded2Dhistograms[ihist]->GetName() << std::endl;
                        }
                    }
                    if (!added) std::cout << " could not add histogram name: " << h2Dhistovector[ihist]->GetName() << std::endl;
                }
                std::cout << "added " << h2Dhistovector.size() << " 2D histograms from file " << currentfile << " with weight =  " << pthatweight << std::endl;
            }
            
            
            filenumber++;
            lastfile = currentfile;
        }
		
	} // end file loop
    
    
    
	
	//////////////////////////////////////////////////
	
	// create output root file
	Char_t filename[100];
	const char* part = lastfile.Data();
	sprintf(filename,"output_TreeOutputCombiner_%s",part);
	TFile* output = new TFile(inputdir+filename,"RECREATE");
	output->cd();
	
	//////////////////////////////////////////
	// Save all your histograms in a root file
	//////////////////////////////////////////
    
	for (unsigned int ihist=0;ihist<maxN1dhistos;ihist++) {
        alladdedhistograms[ihist]->Write();
    }
    
    for (unsigned int ihist=0;ihist<alladded2Dhistograms.size();ihist++) {
        alladded2Dhistograms[ihist]->Write();
    }
	
	output->Close();
	std::cout << "file " << inputdir+filename << " created." << std::endl;
	
}

Double_t TreeOutputCombiner::getRatioError(double	a, double b, double errora, double errorb) {
	
	Double_t error = 0;
	Double_t ratio = a/b;
	
	error = ratio*ratio*((errora*errora)/(a*a) + (errorb*errorb)/(b*b));
	error = sqrt(error);
	return error;
	
}
