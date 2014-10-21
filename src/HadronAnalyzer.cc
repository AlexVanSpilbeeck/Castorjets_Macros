
#include "HadronAnalyzer.h"

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

// own classes includes
#include "../src/MyEvtId.h"
#include "../src/MyGenPart.h"
#include "../src/MyGenJet.h"

HadronAnalyzer::HadronAnalyzer() { }

HadronAnalyzer::~HadronAnalyzer() { }


void HadronAnalyzer::Loop(TString inputdir, TObjArray* filelist, double cmenergy, const char* outputname) {
	
	
	TString tstring = outputname;
	std::cout << " TString outputname = " << tstring << std::endl;
	
	using namespace std;
	int it = 0;
	int totalevents = 0;
	
	/////////////////////////////////////
	// Define all histograms
	/////////////////////////////////////
	
	// eflow histograms
	
	// ranges for different CMS energies & data-MC
	double min = 0;
	double max = 0;
	double nBins = 0;
	double ptcut = 10;
	
	if (cmenergy == 900) {
		// ranges for 900 GeV - use this
		// data
		min = 0.015*-2000;
		max = 0.015*60000;
		nBins = 124;
		ptcut = 10;
	} else if (cmenergy == 2760) {
        // ranges for 2760 GeV - use this
		min = 0.015*-2000;
        max = 0.015*184000; // was 100000
        nBins = 279; // was 204
		ptcut = 10;
	} else if (cmenergy == 7000) {
        // ranges for 7 TeV - use this
		min = -30; //0.015*-2000; // to get whole spectrum: -2000
        max = 3750; //0.015*250000; // to get whole spectrum: 250000
		nBins = 252; // was 252
		ptcut = 10;
	}
	
	int NptBins = 8;
	double ptbinning[9] = {0.3,1,2,3,5,7.5,10,15,25};
	
	// hadron level histograms
	TH1D *heflow_hadronlevel;
    TH1D *hEMeflow_hadronlevel;
    TH1D *hHADeflow_hadronlevel;
    TH1D *heflow_dijet_hadronlevel;
	TH1D *heflow_leadingjet_hadronlevel;
    TH1D *hEMeflow_leadingjet_hadronlevel;
    TH1D *hHADeflow_leadingjet_hadronlevel;
	TH1D *heflow_LeadingjetRatio_hadronlevel = new TH1D("heflow_LeadingjetRatio_hadronlevel","Leadingjet Ratio Hadron level",1,0,1);
	
	heflow_hadronlevel = new TH1D("heflow_hadronlevel","Hadron level MB energy flow in CASTOR",nBins,min,max);
    hEMeflow_hadronlevel = new TH1D("hEMeflow_hadronlevel","EM Hadron level MB energy flow in CASTOR",nBins,min,max);
    hHADeflow_hadronlevel = new TH1D("hHADeflow_hadronlevel","HAD Hadron level MB energy flow in CASTOR",nBins,min,max);
	heflow_dijet_hadronlevel = new TH1D("heflow_dijet_hadronlevel","Hadron level dijet energy flow in CASTOR",nBins,min,max);
	heflow_leadingjet_hadronlevel = new TH1D("heflow_leadingjet_hadronlevel","Hadron level leadingjet energy flow in CASTOR",nBins,min,max);
    hEMeflow_leadingjet_hadronlevel = new TH1D("hEMeflow_leadingjet_hadronlevel","EM Hadron level leadingjet energy flow in CASTOR",nBins,min,max);
    hHADeflow_leadingjet_hadronlevel = new TH1D("hHADeflow_leadingjet_hadronlevel","HAD Hadron level leadingjet energy flow in CASTOR",nBins,min,max);
	
	TH1D *hhadron_energy = new TH1D("hhadron_energy","Hadron level energy distribution",nBins,min,max);
	TH1D *hhadron_eta = new TH1D("hhadron_eta","Hadron level eta distribution",100,-15,15);
	TH1D *hhadron_phi = new TH1D("hhadron_phi","Hadron level phi distribution",100,-3.15,3.15);
    
    TH1D *hhadron_eta_energyweighted_leadingjet = new TH1D("hhadron_eta_energyweighted_leadingjet","Energy weighted eta distribution of stable GenParts - leading jet events",100,-15,15);
    TH1D *hhadron_eta_energyweighted_MB = new TH1D("hhadron_eta_energyweighted_MB","Energy weighted eta distribution of stable GenParts - MB events",100,-15,15);
	
    //double HFbins[6] = {3.2,3.54,3.88,4.22,4.56,4.9};
    //TH1D *hhadron_eta_energyweighted_MB = new TH1D("hhadron_eta_energyweighted_MB","Energy weighted eta distribution of stable GenParts - MB events",5,HFbins);
    
	heflow_hadronlevel->Sumw2();
	heflow_dijet_hadronlevel->Sumw2();
	heflow_leadingjet_hadronlevel->Sumw2();
	hhadron_energy->Sumw2();
	hhadron_eta->Sumw2();
	hhadron_phi->Sumw2();
    hhadron_eta_energyweighted_leadingjet->Sumw2();
    hhadron_eta_energyweighted_MB->Sumw2();
	
	char name [100];
	char title [100];
	
	TH1D *heflow_dijetpt_hadronlevel[NptBins];
	TH1D *heflow_leadingjet_chargedgenjets_hadronlevel[NptBins];
	
	for (int i=0;i<NptBins;i++) {
		sprintf(name,"heflow_dijetpt_hadronlevel_%d",i+1);
		sprintf(title,"Hadron level Dijet Energy flow distribution with pt cut %d",i+1);
		heflow_dijetpt_hadronlevel[i] = new TH1D(name,title,nBins,min,max);
		heflow_dijetpt_hadronlevel[i]->Sumw2();
		
		sprintf(name,"heflow_leadingjet_chargedgenjets_hadronlevel_%d",i+1);
		sprintf(title,"ChargedGenJets Hadron level Leadingjet Energy flow distribution with pt cut %d",i+1);
		heflow_leadingjet_chargedgenjets_hadronlevel[i] = new TH1D(name,title,nBins,min,max);
		heflow_leadingjet_chargedgenjets_hadronlevel[i]->Sumw2();
		
	}
	
	
	
	// put sumw2's for some histos
	heflow_LeadingjetRatio_hadronlevel->Sumw2();
	
	TIter       next(filelist); 
	TObjString* fn = 0;
	TString currentfile = "";
	
	
	std::cout << "start looping over files" << std::endl;
	
	// start file loop
	while((fn = (TObjString*)next()) ) { 
		
		currentfile.Clear();
		currentfile = fn->GetString();
		
		std::cout << "opening the file..." << std::endl;
		
		TFile* f = TFile::Open(inputdir+fn->GetString(),"READ");
		if (!f) { std::cout << "Error in HadronAnalyzer: could not open file " << fn->GetString() << std::endl; continue; } 
		// do what ever 
		std::cout << "opened " << fn->GetString() << std::endl;
		
		//////////////////////////////////////////////////
		// Get tree from the files and define all branches
		//////////////////////////////////////////////////
		
		// get tree from file
		TTree *tree = new TTree("CastorTree","");
		f->GetObject("castortree/CastorTree",tree);
		
		// define objects and branches
		std::vector<MyGenPart> *genParts = NULL;
		std::vector<MyGenJet> *genJets = NULL;
		std::vector<MyGenJet> *chargedGenJets = NULL;
		
		TBranch *b_genParts = NULL;
		b_genParts = tree->GetBranch("GenPart");
		TBranch *b_genJets = NULL;
		TBranch *b_chargedGenJets = NULL;
		b_genJets = tree->GetBranch("GenJet");
		b_chargedGenJets = tree->GetBranch("ChargedGenJet");
		
		b_genParts->SetAddress(&genParts);
		b_genJets->SetAddress(&genJets);
		b_chargedGenJets->SetAddress(&chargedGenJets);
		
		int Nevents = tree->GetEntriesFast();
		std::cout << "events in this file = " << Nevents << std::endl;
		//totalevents += Nevents;
		
		// start event loop
		for (int i=0;i<Nevents;i++) {
			
			bool passedHadronCuts = false;
			
			double hadron_eflow = -99999;
			double hadron_eflow_dijet = -99999;
			double hadron_eflow_leadingjet = -99999;
			
			double xix = 10;
			double xiy = 10;
			double xi = 10;
			double xidd = 10e10;
			double ymax = -1;
			
			/////////////////////////////////////////
			// Hadron level code
			/////////////////////////////////////////
			
			b_genParts->GetEntry(i);
			
			
			// look for particles in HF+ & HF- in BSC range
			bool HadronHFplus = false;
			bool HadronHFminus = false;
			for (unsigned int ipart = 0;ipart<genParts->size();ipart++) {
				MyGenPart particle = (*genParts)[ipart];
				if (particle.status == 1) {
					if (fabs(particle.pdgId) != 12 && fabs(particle.pdgId) != 14 && fabs(particle.pdgId) != 16) {
						if (particle.Eta() > 3.9 && particle.Eta() < 4.4 && particle.charge != 0 /*particle.E() > 4.*/) HadronHFplus = true;
						if (particle.Eta() < -3.9 && particle.Eta() > -4.4 && particle.charge != 0 /*particle.E() > 4.*/) HadronHFminus = true;
					}
				}
			}
			/*
			// CASTOR activity
			bool HadronCASTOR_Activity = false;
			for (unsigned int ipart = 0;ipart<genParts->size();ipart++) {
				MyGenPart particle = (*genParts)[ipart];
				if (particle.status == 1) {
					if (fabs(particle.pdgId) != 12 && fabs(particle.pdgId) != 14 && fabs(particle.pdgId) != 16) {
						if (particle.Eta() < -5.2 && particle.Eta() > -6.6 && particle.E() > 3.) HadronCASTOR_Activity = true;
					}
				} 
			}
			*/
            /*
			// vertex condition
			bool Hadronvertexok = false;
			int Ncharged = 0;
			for (unsigned int ipart = 0;ipart<genParts->size();ipart++) {
				MyGenPart particle = (*genParts)[ipart];
				if (particle.status == 1) {
					if (fabs(particle.Eta()) < 2.5 && particle.Pt() > 0.3 && particle.charge != 0) Ncharged++; // Pt cut ?
				}
			}
			if (Ncharged > 3) Hadronvertexok = true;
			*/
            
			// combine selection
			//if (HadronHFplus && HadronHFminus /* && HadronCASTOR_Activity*/) passedHadronCuts = true;
			
			
			
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
				double min_y = 10000;
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
			xix = Mx2/(cmenergy*cmenergy);
			xiy = My2/(cmenergy*cmenergy);
			
			// xi of event is max
			xi = std::max(xix,xiy);
			xidd=xix*xiy*cmenergy*cmenergy/(0.938*0.938);
			
			// combine selection
			if (cmenergy == 900) {
                if (xix>0.1 || xiy>0.4 || xidd>0.5) passedHadronCuts = true;
            } else if (cmenergy == 2760) {
                if (xix>0.07 || xiy>0.2 || xidd>0.5) passedHadronCuts = true;
            } else if (cmenergy == 7000) {
                if (xix>0.04 || xiy>0.1 || xidd>0.5) passedHadronCuts = true;
            }			
			
			// calculate CASTOR hadron level energy flow
			
			double CASTOReflow = 0;
            double CASTOREMeflow = 0;
            double CASTORHADeflow = 0;
			for (unsigned int ipart = 0;ipart<genParts->size();ipart++) {
				MyGenPart particle = (*genParts)[ipart];
				if (particle.status == 1) {
					// control plots
					hhadron_energy->Fill(particle.E());
					hhadron_eta->Fill(particle.Eta());
					hhadron_phi->Fill(particle.Phi());
					// fill castor eflow
					if (particle.Eta() < -5.2 && particle.Eta() > -6.6 && particle.E() > 0.) {
						bool cond1 = false;
						bool cond2 = false;
						bool cond3 = false;
						if (particle.pdgId != 12 && particle.pdgId != 14 && particle.pdgId != 16) cond1 = true;
						if (particle.pdgId != -12 && particle.pdgId != -14 && particle.pdgId != -16) cond2 = true;
						if (particle.pdgId != 13 && particle.pdgId != -13) cond3 = true;
						if (cond1 && cond2 && cond3) {
							CASTOReflow += particle.E();
                            // count EM energy
                            if (particle.pdgId == 11 || particle.pdgId == -11 || particle.pdgId == 22) {
                                CASTOREMeflow += particle.E();
                            } else {
                                CASTORHADeflow += particle.E();
                            }
						}
					}
				}
			}
			
			hadron_eflow = CASTOReflow;
			
			// if event passed hadron cut, execute analysis code
			if (passedHadronCuts) {
				
				heflow_hadronlevel->Fill(CASTOReflow);
                hEMeflow_hadronlevel->Fill(CASTOREMeflow);
                hHADeflow_hadronlevel->Fill(CASTORHADeflow);
                
                // fill MB energy weighted eta histogram for all stable GenParts
                for (unsigned int ipart = 0;ipart<genParts->size();ipart++) {
                    MyGenPart particle = (*genParts)[ipart];
                    if (particle.status == 1) {
                        hhadron_eta_energyweighted_MB->Fill(particle.Eta(),particle.E());
                    }
                }
				
				// do dijet selection 
				b_genJets->GetEntry(i);
				b_chargedGenJets->GetEntry(i);
				
				// check for dijet system
				if (isGenDiJet(*genJets,true,2.5,ptcut)) {
					heflow_dijet_hadronlevel->Fill(CASTOReflow);
					hadron_eflow_dijet = CASTOReflow;
				}
				// check for leading charged genjet
				int posLeadingChargedGenJet = posLeadingGenJet(*chargedGenJets,2.,1.);
				if (posLeadingChargedGenJet != -1) {
					if ((*chargedGenJets)[posLeadingChargedGenJet].Pt() > ptcut) {
						heflow_leadingjet_hadronlevel->Fill(CASTOReflow);
                        hEMeflow_leadingjet_hadronlevel->Fill(CASTOREMeflow);
                        hHADeflow_leadingjet_hadronlevel->Fill(CASTORHADeflow);
						hadron_eflow_leadingjet = CASTOReflow;
                        
                        // fill leadingjet energy weighted eta histogram for all stable GenParts
                        for (unsigned int ipart = 0;ipart<genParts->size();ipart++) {
                            MyGenPart particle = (*genParts)[ipart];
                            if (particle.status == 1) {
                                hhadron_eta_energyweighted_leadingjet->Fill(particle.Eta(),particle.E());
                            }
                        }
					}
                    
					for (int ipt=0;ipt<NptBins;ipt++) {
						if ((*chargedGenJets)[posLeadingChargedGenJet].Pt() > ptbinning[ipt] && (*chargedGenJets)[posLeadingChargedGenJet].Pt() < ptbinning[ipt+1]) 
							heflow_leadingjet_chargedgenjets_hadronlevel[ipt]->Fill(CASTOReflow);
					}
				}
				
				// pt evolution, do manual dijet selection
				for (int ipt=0;ipt<NptBins;ipt++) {
					if ( isGenDiJet(*genJets,true,2.5,ptbinning[ipt],ptbinning[ipt+1],true)) heflow_dijetpt_hadronlevel[ipt]->Fill(CASTOReflow);
				}
				
				// end of event, print status
				if( ((i+1) % 10000) == 0) std::cout << i+1 <<"events done in file " << it << std::endl;
				totalevents++;
				
				
			} // end if passedHadronCuts
			
		} // end event loop
		
		delete tree;
		f->Close();
		it++;
	} // end file loop
	
	
	//////////////////////////////////////////////////
	// energy flow stuff needed after event filling //
	//////////////////////////////////////////////////
	
	std::cout << "starting eflow calculations" << std::endl;
	
	// check all dijet pt histos for under or overflow
	for (int ipt=0;ipt<NptBins;ipt++) {
		checkFlow(heflow_dijetpt_hadronlevel[ipt]);
		checkFlow(heflow_leadingjet_chargedgenjets_hadronlevel[ipt]);
	}
	
	checkFlow(heflow_hadronlevel);
	checkFlow(heflow_dijet_hadronlevel);
	checkFlow(heflow_leadingjet_hadronlevel);
	
	// do eflow ratio stuff	
	double leadingjetratio_hadronlevel = heflow_leadingjet_hadronlevel->GetMean()/heflow_hadronlevel->GetMean();
	double leadingjetratio_error_hadronlevel = getRatioError(heflow_hadronlevel,heflow_leadingjet_hadronlevel);
	heflow_LeadingjetRatio_hadronlevel->SetBinContent(1,leadingjetratio_hadronlevel);
	heflow_LeadingjetRatio_hadronlevel->SetBinError(1,leadingjetratio_error_hadronlevel);
	
	// now print all the leadingjet ratio's and their errors
	std::cout << "Hadron level MinBias energy flow (" << cmenergy << " GeV) = " << heflow_hadronlevel->GetMean() << " and Leadingjet energy flow (" << cmenergy << " GeV, ptcut = " << ptcut << " GeV) = " << heflow_leadingjet_hadronlevel->GetMean() << std::endl;
	std::cout << "Hadron level Leadingjet ratio (default at " << cmenergy << " GeV, ptcut = " << ptcut << " GeV) = " << leadingjetratio_hadronlevel << " +- " << leadingjetratio_error_hadronlevel << std::endl;
	
	// make leadingjet ratio vs pt plot on hadron level
	TH1D *hDijetRatios_vs_pt_hadronlevel = new TH1D("hDijetRatios_vs_pt_hadronlevel","Hadron level Dijet ratios vs dijet ptcut",NptBins,ptbinning);
	TH1D *hLeadingJetRatios_vs_pt_chargedgenjets_hadronlevel = new TH1D("hLeadingJetRatios_vs_pt_chargedgenjets_hadronlevel","Hadron level LeadingJet ratios vs leadingjet ptcut - charged genjets",NptBins,ptbinning);
	for (int ipt=0;ipt<NptBins;ipt++) {
		hDijetRatios_vs_pt_hadronlevel->SetBinContent(ipt+1,heflow_dijetpt_hadronlevel[ipt]->GetMean()/heflow_hadronlevel->GetMean());
		hDijetRatios_vs_pt_hadronlevel->SetBinError(ipt+1,getRatioError(heflow_hadronlevel,heflow_dijetpt_hadronlevel[ipt]));
		
		hLeadingJetRatios_vs_pt_chargedgenjets_hadronlevel->SetBinContent(ipt+1,heflow_leadingjet_chargedgenjets_hadronlevel[ipt]->GetMean()/heflow_hadronlevel->GetMean());
		hLeadingJetRatios_vs_pt_chargedgenjets_hadronlevel->SetBinError(ipt+1,getRatioError(heflow_hadronlevel,heflow_leadingjet_chargedgenjets_hadronlevel[ipt]));
	}
	
	
	//////////////////////////////////////////////////
	
	cout << "total number of events = " << totalevents << " from " << it << " file(s)" << endl;
	
	
	// create output root file
	Char_t filename[200];
	//char* temp = const_cast<char*>(outputname);
	const char* part = currentfile.Data();
	string first(outputname);
	first += part;
	//strcat(temp,part);
	sprintf(filename,"output_HadronAnalyzer_%s",first.c_str());
	TFile* output = new TFile(filename,"RECREATE");
	output->cd();
	
	//////////////////////////////////////////
	// Save all your histograms in a root file
	//////////////////////////////////////////
	
	// eflow histos
	
	// hadron level
	heflow_hadronlevel->Write();
    hEMeflow_hadronlevel->Write();
    hHADeflow_hadronlevel->Write();
	heflow_dijet_hadronlevel->Write();
	heflow_leadingjet_hadronlevel->Write();
    hEMeflow_leadingjet_hadronlevel->Write();
    hHADeflow_leadingjet_hadronlevel->Write();
	heflow_LeadingjetRatio_hadronlevel->Write();
	hhadron_eta->Write();
	hhadron_phi->Write();
	hhadron_energy->Write();
    hhadron_eta_energyweighted_MB->Write();
    hhadron_eta_energyweighted_leadingjet->Write();
	
	for (int ipt=0;ipt<NptBins;ipt++) {
		heflow_dijetpt_hadronlevel[ipt]->Write();
		heflow_leadingjet_chargedgenjets_hadronlevel[ipt]->Write();
	}
	hDijetRatios_vs_pt_hadronlevel->Write();
	hLeadingJetRatios_vs_pt_chargedgenjets_hadronlevel->Write();
	
	output->Close();
	std::cout << "file " << filename << " created." << std::endl;
	
}

Double_t HadronAnalyzer::getRatioError(TH1D * hMB, TH1D * hQCD) {
	
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

Double_t HadronAnalyzer::getRatioError(double	a, double b, double errora, double errorb) {
	
	Double_t error = 0;
	Double_t ratio = a/b;
	
	error = ratio*ratio*((errora*errora)/(a*a) + (errorb*errorb)/(b*b));
	error = sqrt(error);
	return error;
	
}

bool HadronAnalyzer::isGenDiJet(std::vector<MyGenJet> JetVector, bool backtoback, double usedetacut, double usedptcut) {
	
	using namespace std;
	
	//-- central jet selection
	
	bool accept = false;
	
	short int posJet1 = -1;
	short int posJet2 = -1;
	
	double refPt1 = 0;
	double refPt2 = 0;
	
	//-- find the two highest pt jets (corrected pt)
	
	for(vector<MyGenJet>::const_iterator jet = JetVector.begin(); jet < JetVector.end(); ++jet) {
		
		Double_t pt = jet->Pt();
		
		if(pt > refPt1) {
			refPt2 = refPt1;
			posJet2 = posJet1;      
			refPt1 = pt;
			posJet1 = jet - JetVector.begin();
		}
		
		else if(pt > refPt2) {
			refPt2 = pt;
			posJet2 = jet - JetVector.begin();
		}
		
	} 
	
	//-- apply the tight selection to them
	
	if(posJet1 >= 0 && posJet2 >= 0) {
		
		bool accept_jet1 = false;
		bool accept_jet2 = false;
		
		MyGenJet jet1 = JetVector.at(posJet1);
		MyGenJet jet2 = JetVector.at(posJet2);
		
		//-- jet 1 selection
		if(jet1.Pt() > usedptcut && fabs(jet1.Eta()) < usedetacut) accept_jet1 = true;
		
		//-- jet 2 selection
		if(jet2.Pt() > usedptcut && fabs(jet2.Eta()) < usedetacut) accept_jet2 = true;
		
		//-- final selection (back-to-back)
		
		if(accept_jet1 == true && accept_jet2 == true) {
			if (backtoback) {
				double deltaPhi = fabs(deltaPhi2(jet1.Phi(), jet2.Phi()));
				if (fabs(deltaPhi - M_PI) < 1.0) accept = true;
			} else {
				accept = true;
			}
			
		}
		
	} //-- posJet1 >= 0 && posJet2 >= 0
	return accept;
}

bool HadronAnalyzer::isGenDiJet(std::vector<MyGenJet> JetVector, bool backtoback, double usedetacut, double lowptcut, double highptcut, bool mean) {
	
	using namespace std;
	
	//-- central jet selection
	
	bool accept = false;
	
	short int posJet1 = -1;
	short int posJet2 = -1;
	
	double refPt1 = 0;
	double refPt2 = 0;
	
	//-- find the two highest pt jets (corrected pt)
	
	for(vector<MyGenJet>::const_iterator jet = JetVector.begin(); jet < JetVector.end(); ++jet) {
		
		Double_t pt = jet->Pt();
		
		if(pt > refPt1) {
			refPt2 = refPt1;
			posJet2 = posJet1;      
			refPt1 = pt;
			posJet1 = jet - JetVector.begin();
		}
		
		else if(pt > refPt2) {
			refPt2 = pt;
			posJet2 = jet - JetVector.begin();
		}
		
	} 
	
	//-- apply the tight selection to them
	
	if(posJet1 >= 0 && posJet2 >= 0) {
		
		bool accept_jet1 = false;
		bool accept_jet2 = false;
		
		MyGenJet jet1 = JetVector.at(posJet1);
		MyGenJet jet2 = JetVector.at(posJet2);
		
		if (!mean) {
			//-- jet 1 selection
			if(jet1.Pt() > lowptcut && jet1.Pt() < highptcut && fabs(jet1.Eta()) < usedetacut) accept_jet1 = true;
			
			//-- jet 2 selection
			if(jet2.Pt() > lowptcut && jet2.Pt() < highptcut && fabs(jet2.Eta()) < usedetacut) accept_jet2 = true;
			
			//-- final selection (back-to-back)
			
			if(accept_jet1 == true && accept_jet2 == true) {
				if (backtoback) {
					double deltaPhi = fabs(deltaPhi2(jet1.Phi(), jet2.Phi()));
					if (fabs(deltaPhi - M_PI) < 1.0) accept = true;
				} else {
					accept = true;
				}
				
			}
		} else {
			// using mean dijet system pt
			//-- jet 1 selection
			if(fabs(jet1.Eta()) < usedetacut) accept_jet1 = true;
			
			//-- jet 2 selection
			if(fabs(jet2.Eta()) < usedetacut) accept_jet2 = true;
			
			//-- final selection (back-to-back)
			
			if(accept_jet1 == true && accept_jet2 == true) {
				double meanpt = (jet1.Pt() + jet2.Pt())/2;
				bool meanOK = false;
				if (meanpt > lowptcut && meanpt < highptcut) meanOK = true;
				if (backtoback && meanOK) {
					double deltaPhi = fabs(deltaPhi2(jet1.Phi(), jet2.Phi()));
					if (fabs(deltaPhi - M_PI) < 1.0) accept = true;
				} else if (meanOK) {
					accept = true;
				}
				
			}
		}
		
	} //-- posJet1 >= 0 && posJet2 >= 0
	return accept;
}


Double_t HadronAnalyzer::deltaPhi2(double phi1, double phi2) { 
    double result = phi1 - phi2;
    while (result > M_PI) result -= 2*M_PI;
    while (result <= -M_PI) result += 2*M_PI;
    return result;
}

void HadronAnalyzer::checkFlow(TH1D *histo) {
	if (histo->GetBinContent(histo->GetNbinsX()+1) !=0) std::cout << "!!WARNING!! histogram " << histo->GetName() << " has overflow" << std::endl;
	if (histo->GetBinContent(0) !=0) std::cout << "!!WARNING!! histogram " << histo->GetName() << " has underflow" << std::endl;
}


int HadronAnalyzer::posLeadingGenJet(std::vector<MyGenJet> JetVector,double etacut,double minptcut) {
	
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




