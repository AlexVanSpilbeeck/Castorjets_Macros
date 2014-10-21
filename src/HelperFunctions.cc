
#include "HelperFunctions.h"
#include "HistoRetriever.h"

//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>

//STANDARD C++ INCLUDES
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include<sys/stat.h>
#include<sys/types.h>

HelperFunctions::HelperFunctions() { }

HelperFunctions::~HelperFunctions() { }

Double_t HelperFunctions::getRatioError(TH1D * hMB, TH1D * hQCD) {
	
	Double_t error = 0;
	Double_t ratio = hQCD->GetMean()/hMB->GetMean();
	Double_t QCDmeanerror = hQCD->GetMeanError();
	Double_t MBmeanerror = hMB->GetMeanError();
	Double_t QCDmean = hQCD->GetMean();
	Double_t MBmean = hMB->GetMean();
	
	error = ratio*ratio*((QCDmeanerror*QCDmeanerror)/(QCDmean*QCDmean) + (MBmeanerror*MBmeanerror)/(MBmean*MBmean));
	error = sqrt(error);
	return error;
	
}

Double_t HelperFunctions::getRatioError(double	a, double b, double errora, double errorb) {
	
	Double_t error = 0;
	Double_t ratio = a/b;
	
	error = ratio*ratio*((errora*errora)/(a*a) + (errorb*errorb)/(b*b));
	error = sqrt(error);
	return error;
	
}

Double_t HelperFunctions::getMultiError(double	a, double b, double errora, double errorb) {
	
	Double_t error = 0;
	
	error = (a*a)*(errorb*errorb) + (errora*errora)*(b*b);
	error = sqrt(error);
	return error;
	
}

TH1D HelperFunctions::getHistoByName(std::vector<TH1D*> allhistograms,TString name) {
    
    TH1D histo;
    for (unsigned int ihist=0;ihist<allhistograms.size();ihist++) {
        if (strcmp(allhistograms[ihist]->GetName(),name) == 0) {
            histo = (*allhistograms[ihist]);
        }
    }
    
    return histo;
    
}

Double_t HelperFunctions::deltaPhi2(double phi1, double phi2) { 
    double result = phi1 - phi2;
    while (result > M_PI) result -= 2*M_PI;
    while (result <= -M_PI) result += 2*M_PI;
    return result;
}

void HelperFunctions::checkFlow(TH1D *histo) {
	if (histo->GetBinContent(histo->GetNbinsX()+1) !=0) std::cout << "!!WARNING!! histogram " << histo->GetName() << " has overflow" << std::endl;
	if (histo->GetBinContent(0) !=0) std::cout << "!!WARNING!! histogram " << histo->GetName() << " has underflow" << std::endl;
}

int HelperFunctions::rad2Sector(double radians) {
	// convert radian angle to CASTOR sector number
	int result = -1;
	double sectorwidthinrad = M_PI/8;
	if (radians > 0.) {
		for (int i=0;i<8;i++) {
			if (radians > i*sectorwidthinrad && radians < (i+1)*sectorwidthinrad) result = i+1;
		}
	} else {
		for (int i=8;i<16;i++) {
			if (radians > -1*(16-i)*sectorwidthinrad && radians < -1*(16-i-1)*sectorwidthinrad) result = i+1;
		}
	}
	return result;
	
} 

int HelperFunctions::AddSector(int sector,int add) {
	// return the next sector: sector+1
	if (sector+add < 17) {
		return sector+add;
	} else {
		return sector+add-16;
	}
}

int HelperFunctions::SubtractSector(int sector,int sub) {
	// return the next sector: sector+1
	if (sector-sub > 0) {
		return sector-sub;
	} else {
		return sector-sub+16;
	}
}
