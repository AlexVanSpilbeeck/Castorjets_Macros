////////////////////////////////////////
//////// New CMSSW_4_2_X version ///////
////////////////////////////////////////


#include "JetAnalyzer.h"
#include "HistoRetriever.h"

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
#include <TThread.h>
#include <TStopwatch.h>

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

TFile *JetAnalyzer::currentStaticTFile_ = new TFile();

JetAnalyzer::JetAnalyzer(TString inputdir, TObjArray* filelist, bool isData, const char* outputname) {
    
	std::cout << "constructing JetAnalyzer class..." << std::endl;
	
    inputdir_ = inputdir;
    filelist_ = filelist;
    isData_ = isData;
    outputname_ = outputname;
	
	std::cout << "initializing basic variables..." << std::endl;
    
    // initialize basic variables

    
    LoopOutputFile_ = "";
	
	currentfile_ = "";
	currentTFile_ = new TFile();
	
	std::cout << "all initialisations done, class constructed!" << std::endl;
    
}

JetAnalyzer::~JetAnalyzer() { }

void JetAnalyzer::Loop() {
	
	
	std::cout << " JetAnalyzer Loop function is started " << std::endl;
	
	TString tstring = outputname_;
	std::cout << " TString outputname_ = " << tstring << std::endl;
	
    // reweight the MC in this case
    bool reweightMC = false;
    if (tstring.Contains("Reweighted")) {
        std::cout << "We will reweight the MC now !!!" << std::endl;
        reweightMC = true;
    }
    
	
	using namespace std;
	int it = 0;
	int totalevents = 0;
	
	/////////////////////////////////////
	// Define all histograms
	/////////////////////////////////////
	
	// detector level histograms
	
	char name [100];
	char title [100];
		
	TH1D *hCASTORTowerMulti = new TH1D("hCASTORTowerMulti","CASTOR Tower Multiplicity (N above threshold) distribution",17,0,17);
	hCASTORTowerMulti->Sumw2();
	
	// default energy flow histos - using 5 modules
	TH1D *hCASTOReflow = new TH1D("hCASTOReflow","Total CASTOR energy flow in first 5 modules",252,-30,3750);
	hCASTOReflow->Sumw2();
	
	TH2D *h2CASTOReflow_grid = new TH2D("h2CASTOReflow_grid","CASTOR energy weighted module vs sector distribution",16,1,17,14,1,15);
	
	TH1D *hCASTOReflow_channel[224];
	for (int i=0;i<224;i++) {
		sprintf(name,"hCASTOReflow_channel_%d",i+1);
		sprintf(title,"CASTOR Energy distribution for channel %d",i+1);
		hCASTOReflow_channel[i] = new TH1D(name,title,100,-10,1800);
		hCASTOReflow_channel[i]->Sumw2();
	}
	
	// Castor jet histograms
	TH1D *hCastorJet_energy = new TH1D("hCastorJet_energy","CastorJet energy distribution",200,0,2000);
	TH1D *hCastorJet_pt = new TH1D("hCastorJet_pt","CastorJet pt distribution",30,0,30);
	TH1D *hCastorJet_em = new TH1D("hCastorJet_em","CastorJet EM energy distribution",150,0,1500);
	TH1D *hCastorJet_had = new TH1D("hCastorJet_had","CastorJet HAD energy distribution",150,0,1500);
	TH1D *hCastorJet_fem = new TH1D("hCastorJet_fem","CastorJet EM/(EM+HAD) distribution",60,-0.1,1.1);
	TH1D *hCastorJet_fhot = new TH1D("hCastorJet_fhot","CastorJet Fhot distribution",60,-0.1,1.1);
	TH1D *hCastorJet_width = new TH1D("hCastorJet_width","CastorJet width distribution",100,0,1);
	TH1D *hCastorJet_depth = new TH1D("hCastorJet_depth","CastorJet depth distribution",100,-16000,-14000);
	TH1D *hCastorJet_sigmaz = new TH1D("hCastorJet_sigmaz","CastorJet sigmaz distribution",100,0,500);
	TH1D *hCastorJet_ntower = new TH1D("hCastorJet_ntower","CastorJet ntower distribution",16,1,17);
	TH1D *hCastorJet_eta = new TH1D("hCastorJet_eta","CastorJet eta distribution",14,-6.6,-5.2);
	TH1D *hCastorJet_phi = new TH1D("hCastorJet_phi","CastorJet phi distribution",16,-M_PI,+M_PI);
	TH1D *hCastorJet_multi = new TH1D("hCastorJet_multi","CastorJet multiplicity distribution",17,0,17);
	
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
	
	TIter       next(filelist_); 
	TObjString* fn = 0;
	
	bool isMC = false;
    
	
	std::cout << "start looping over files" << std::endl;
	
	// start file loop
	while((fn = (TObjString*)next()) ) { 
		
		currentfile_.Clear();
		currentTFile_->Clear();
		
		currentfile_ = fn->GetString();
		
		std::cout << "opening file " << currentfile_ << " ... " << std::endl;
		
		TStopwatch *timer = new TStopwatch();
		TThread *fileopener;
		fileopener = new TThread("fileopener",(void(*) (void *))&OpenROOTFile,(JetAnalyzer*) this);
		TThread::SetCancelAsynchronous();
		TThread::SetCancelOn();
		fileopener->Run();
		
		bool fileopen = false;
		bool timeout = true;
		while (timer->RealTime() < 30) {
			if (fileopener->GetState() != TThread::kRunningState) {
				// file open thread is not running anymore, was it successful?
				// set the static TFile to the instance one
				setCurrentTFile();
				timeout = false;
				break;
			} 
			
			timer->Continue();
			gSystem->Sleep(1000);
		}
		
		// after the time, stop & delete the open thread
		TThread::Kill("fileopener");
		TThread::Delete(fileopener);
		delete fileopener;
		delete timer;
		
		if (timeout) {
			std::cout << "file open thread has timed out, go to next file in the list" << std::endl;
			continue;
		} else {
			// is this instance open?
			if (currentTFile_->IsOpen()) {
				std::cout << "the currentTFile_ is correctly opened" << std::endl;
				fileopen = true;
			}
		}
		
		if (!fileopen) {
			std::cout << "no timeout, but the file could not be opened correctly: going to the next file in the list" << std::endl;
			continue;
		}
		
		
		//////////////////////////////////////////////////
		// Get tree from the files and define all branches
		//////////////////////////////////////////////////
		
		// get tree from file
		TTree *tree = new TTree("CastorTree","");
		currentTFile_->GetObject("castortree/CastorTree",tree);
		
		// define objects and branches
		MyEvtId *evtid = NULL;
        MyGenKin *evtkin = NULL;
		MyBeamSpot *BeamSpot = NULL;
		MyHLTrig *HLTrig = NULL;
		MyL1Trig *L1Trig = NULL;
		std::vector<MyVertex> *Vertices = NULL;
		std::vector<MyCastorRecHit> *CastorRecHits = NULL;
		std::vector<MyCastorTower> *CastorTowers = NULL;
		std::vector<MyCastorJet> *CastorJets = NULL;
		std::vector<MyJet> *PFJets = NULL;
		std::vector<MyGenPart> *genParts = NULL;
		std::vector<MyCaloTower> *caloTowers = NULL;
		std::vector<MyGenJet> *genJets = NULL;
		std::vector<MyGenJet> *chargedGenJets = NULL;
		std::vector<MyTrackJet> *trackJets = NULL;
		
		TBranch *b_evtid = tree->GetBranch("EvtId");
        TBranch *b_evtkin = NULL;
        if (!isData_) b_evtkin = tree->GetBranch("GenKin");
		TBranch *b_BeamSpot = tree->GetBranch("beamSpot");
		TBranch *b_HLTrig = tree->GetBranch("HLTrig");
		TBranch *b_L1Trig = tree->GetBranch("L1Trig");
		TBranch *b_vertices = tree->GetBranch("primaryVertex");
		TBranch *b_castorrechits = tree->GetBranch("castorRecHit");
		TBranch *b_castortowers = tree->GetBranch("castorTower");
		TBranch *b_castorjets = tree->GetBranch("castorJet");
		TBranch *b_PFJets = tree->GetBranch("pfJet");
		TBranch *b_genParts = NULL;
		if (!isData_) b_genParts = tree->GetBranch("GenPart");
		TBranch *b_caloTowers = tree->GetBranch("caloTower");
		TBranch *b_genJets = NULL;
		TBranch *b_chargedGenJets = NULL;
		if (!isData_) b_genJets = tree->GetBranch("GenJet");
		if (!isData_) b_chargedGenJets = tree->GetBranch("ChargedGenJet");
		TBranch *b_trackJets = tree->GetBranch("trackJet");
		
		b_evtid->SetAddress(&evtid);
        if (!isData_) b_evtkin->SetAddress(&evtkin);
		b_BeamSpot->SetAddress(&BeamSpot);
		b_HLTrig->SetAddress(&HLTrig);
		b_L1Trig->SetAddress(&L1Trig);
		b_vertices->SetAddress(&Vertices);
		b_castorrechits->SetAddress(&CastorRecHits);
		b_castortowers->SetAddress(&CastorTowers);
		b_castorjets->SetAddress(&CastorJets);
		b_PFJets->SetAddress(&PFJets);
		if (!isData_) b_genParts->SetAddress(&genParts);
		b_caloTowers->SetAddress(&caloTowers);
		if (!isData_) b_genJets->SetAddress(&genJets);
		if (!isData_) b_chargedGenJets->SetAddress(&chargedGenJets);
		b_trackJets->SetAddress(&trackJets);
		
		int Nevents = tree->GetEntriesFast();
		std::cout << "file opened, events in this file = " << Nevents << std::endl;
		//totalevents += Nevents;
		
		// start event loop
		for (int i=0;i<Nevents;i++) {
            
			bool passedHadronCuts = false;
			bool passedDetectorCuts = false;
			
			double xix = 10;
			double xiy = 10;
			double xi = 10;
			double xidd = 10e10;
			double ymax = -1;
            
			
			/////////////////////////////////////////
			// Hadron level code
			/////////////////////////////////////////
			
			if (!isData_) {
				b_genParts->GetEntry(i);
				
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
					double min_y = 100000;
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
				xix = Mx2/(7000*7000);
				xiy = My2/(7000*7000);
				
				// xi of event is max
				xi = std::max(xix,xiy);
				xidd=xix*xiy*7000*7000/(0.938*0.938);
                
                    
                // combine selection
        		// 7 TeV Xi cuts
				if (xix>0.04 || xiy>0.1 || xidd>0.5) passedHadronCuts = true;

				
				// if event passed hadron cut, execute analysis code
				if (passedHadronCuts) {
					
					
					
				} // end if passedHadronCuts
				
			} // end if not data
             
			
			
			/////////////////////////////////////////
			// Do stuff before filters
			/////////////////////////////////////////
						
			b_evtid->GetEntry(i);
			b_HLTrig->GetEntry(i);
			b_L1Trig->GetEntry(i);
			b_vertices->GetEntry(i);
			b_caloTowers->GetEntry(i);
			b_castorrechits->GetEntry(i);
			b_castortowers->GetEntry(i);
			
			// only process a certain run
			//if (evtid->Run != 135521) continue; // go to next event
			
			/////////////////////////////////////////
			// Filter the results
			/////////////////////////////////////////
			
			if (isData_) {
				
				// filter results
				// filter on phys declared bit
				bool physDeclresult = HLTrig->HLTmap["physDeclpath"];
				
				// filter on castor invalid data
				bool castorInvalidDataFilterresult = HLTrig->HLTmap["castorInvalidDataFilterpath"];
				
				// filter out scraping events
				bool noscrapingresult = HLTrig->HLTmap["noscrapingpath"];
				
				bool gooddata = physDeclresult && castorInvalidDataFilterresult && noscrapingresult;
				
				// L1 filter

				bool L1_BX = L1Trig->fTechDecisionBefore[0];
				bool L1_Veto = !L1Trig->fTechDecisionBefore[36] && !L1Trig->fTechDecisionBefore[37] && !L1Trig->fTechDecisionBefore[38] && !L1Trig->fTechDecisionBefore[39];
				
				bool L1_BSC = L1Trig->fTechDecisionBefore[40] || L1Trig->fTechDecisionBefore[41]; // default value for 900 GeV or 7000 GeV 
				
				bool HLT_BSC = true; // default value for 900 GeV or 7000 GeV (do we need an HLT here?)
				
				bool TriggerSelection = false;
                
				TriggerSelection = L1_BX && L1_Veto && L1_BSC && HLT_BSC;
				
				// ask for activity in HF
				bool HF_Activity = true;
				bool HFplus = false;
				bool HFminus = false;
				for (unsigned int itow=0;itow<caloTowers->size();itow++) {
					MyCaloTower mytow = (*caloTowers)[itow];
					if (mytow.hasHF && mytow.Energy() > 4. && mytow.zside == 1 && mytow.Eta() > 3.23 && mytow.Eta() < 4.65) HFplus = true;
					if (mytow.hasHF && mytow.Energy() > 4. && mytow.zside == -1 && mytow.Eta() < -3.23 && mytow.Eta() > -4.65) HFminus = true;
				}
				if (!HFplus || !HFminus) HF_Activity = false; // real HF condition... 
				
				// ask for activity in CASTOR
				bool CASTOR_Activity = false;
				// if at least 1 tower is above noise threshold, set CASTOR activity to true
				if (CastorTowers->size() > 0) CASTOR_Activity = true;
				
				// get vertex info
				// do the oneGoodVertexFilter
				bool wehaveGoodVertex = false;
				for (unsigned int iVert=0;iVert<Vertices->size();iVert++) {
					MyVertex vertex = (*Vertices)[iVert];
					if (vertex.isGoodVertex) wehaveGoodVertex = true;
				}
				
				if (gooddata && TriggerSelection && HF_Activity && CASTOR_Activity && wehaveGoodVertex) passedDetectorCuts = true;
				
			} else {
				
				// MC
				
				// ask for activity in HF
				bool HF_Activity = true; 
				bool HFplus = false;
				bool HFminus = false;
				for (unsigned int itow=0;itow<caloTowers->size();itow++) {
					MyCaloTower mytow = (*caloTowers)[itow];
					if (mytow.hasHF && mytow.Energy() > 4. && mytow.zside == 1 && mytow.Eta() > 3.23 && mytow.Eta() < 4.65) HFplus = true;
					if (mytow.hasHF && mytow.Energy() > 4. && mytow.zside == -1 && mytow.Eta() < -3.23 && mytow.Eta() > -4.65) HFminus = true;
				}
				if (!HFplus || !HFminus) HF_Activity = false; // real HF condition... 
				
				// ask for activity in CASTOR
				bool CASTOR_Activity = false;
				// if at least 1 tower is above noise threshold, set CASTOR activity to true
				if (CastorTowers->size() > 0) CASTOR_Activity = true;
                
				// get vertex info
				// do the oneGoodVertexFilter
				bool wehaveGoodVertex = false;
				for (unsigned int iVert=0;iVert<Vertices->size();iVert++) {
					MyVertex vertex = (*Vertices)[iVert];
					if (vertex.isGoodVertex) wehaveGoodVertex = true;
				}
				
				if (HF_Activity && CASTOR_Activity && wehaveGoodVertex) passedDetectorCuts = true;
                
			}
			
				
			/////////////////////////////////////////
			// Do stuff after filters
			/////////////////////////////////////////
			
			// get all the remaining branch entries
			b_trackJets->GetEntry(i);
				
			if (passedDetectorCuts) {
				
				
				/////////////////////////////////////////
				// Start Nvertex == 1 part of the code 
				/////////////////////////////////////////
				
				// only fill the histograms when there's 1 vertex (filter out pile-up)
				if (Vertices->size() == 1) {
					
					// get event id stuff
					if( ((i+1) % 10000) == 0) cout << " run " << evtid->Run << " isData = " << evtid->IsData << " lumiblock " << 
						evtid->LumiBlock << " event " << evtid->Evt << endl; 
					if (!evtid->IsData) isMC = true;
					
					// calculate energy flow in CASTOR
					double CASTOReflow = 0;
					
					for (unsigned int j=0;j<CastorRecHits->size();j++) {
						MyCastorRecHit rechit = (*CastorRecHits)[j];
						//if (rechit.bad) std::cout << " found bad rechit: channel " << rechit.cha << std::endl;
						if (!rechit.bad) {
							hCASTOReflow_channel[rechit.cha-1]->Fill(rechit.energy);
							if (rechit.cha <= 80) CASTOReflow += rechit.energy;
							h2CASTOReflow_grid->Fill(rechit.sec,rechit.mod,rechit.energy);
						}
					}
					
					// tower multiplicity code
					hCASTORTowerMulti->Fill(CastorTowers->size());
					
					// fill all minbias eflow histos
					hCASTOReflow->Fill(CASTOReflow);
										
		
					// analyse castor jets
					int NCastorJets = 0;
					b_castorjets->GetEntry(i);
					for (unsigned int j=0;j<CastorJets->size();j++) {
						MyCastorJet casjet = (*CastorJets)[j];
						if (casjet.energy > 500.) {
							NCastorJets++;
							hCastorJet_energy->Fill(casjet.energy);
							hCastorJet_pt->Fill(casjet.energy*sin(2*atan(exp(5.9))));
							hCastorJet_em->Fill(casjet.eem);
							hCastorJet_had->Fill(casjet.ehad);
							hCastorJet_eta->Fill(casjet.eta);
							hCastorJet_phi->Fill(casjet.phi);
							hCastorJet_fem->Fill(casjet.fem);
							hCastorJet_fhot->Fill(casjet.fhot);
							hCastorJet_width->Fill(casjet.width);
							hCastorJet_depth->Fill(casjet.depth);
							hCastorJet_sigmaz->Fill(casjet.sigmaz);
							hCastorJet_ntower->Fill(casjet.ntower);
						}
					}
					hCastorJet_multi->Fill(NCastorJets);
					
					
					// end of event, print status
					if( ((i+1) % 10000) == 0) std::cout << i+1 <<"events done in file " << it << std::endl;
					totalevents++;
					
				} // end if statement for 1 vertex
								
			} // end if passed detector cuts
			
			// combined hadron and detector level code
			
			if (passedHadronCuts && passedDetectorCuts) {
			}
			
		} // end event loop
		
		delete tree;
		currentTFile_->Close();
		it++;
	} // end file loop
	
	std::cout << "file loop has ended" << std::endl;
    
    currentTFile_->Close();
	
	
    // check all histo's for overflow
	
	hf_.checkFlow(hCASTOReflow);
	hf_.checkFlow(hCASTORTowerMulti);
    
    // check distribution of each channel on under or overflow
	for (int icha=0;icha<224;icha++) {
		hf_.checkFlow(hCASTOReflow_channel[icha]);
	}
	
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
    
    
    // write all histo's to file
    
	std::cout << "total number of events = " << totalevents << " from " << it << " file(s)" << endl;
	
	
	// create output root file
	Char_t filename[200];
	const char* part = currentfile_.Data();
	string first(outputname_);
	first += part;
	//strcat(temp,part);
	sprintf(filename,"output_JetAnalyzer_%s",first.c_str());
	TFile* output = new TFile(filename,"RECREATE");
	output->cd();
		
	//////////////////////////////////////////
	// Save all your histograms in a root file
	//////////////////////////////////////////
	
	// detector level histograms
    
	// eflow histos
	hCASTORTowerMulti->Write();
	hCASTOReflow->Write();
	h2CASTOReflow_grid->Write();
	
	for (int icha=0;icha<224;icha++) {
		hCASTOReflow_channel[icha]->Write();
	}
	
	hCastorJet_energy->Write();
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
		
	output->Close();
	std::cout << "file " << filename << " created." << std::endl;
    LoopOutputFile_ = filename;


}
	
void JetAnalyzer::AfterLoopCalculations(TString file) {

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

TString JetAnalyzer::getOutputFile() {
    return LoopOutputFile_;
}

TString JetAnalyzer::getInputDir() {
	return inputdir_;
}

TString JetAnalyzer::getCurrentFile() {
	return currentfile_;
}

void JetAnalyzer::setCurrentTFile() {
	currentTFile_ = currentStaticTFile_;
}

void* JetAnalyzer::OpenROOTFile(JetAnalyzer* arg) {
	currentStaticTFile_ = TFile::Open(arg->getInputDir()+arg->getCurrentFile(),"READ");
	return 0;
}

int JetAnalyzer::posLeadingGenJet(std::vector<MyGenJet> JetVector,double etacut,double minptcut) {
	
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

int JetAnalyzer::posLeadingTrackJet(std::vector<MyTrackJet> JetVector,double etacut,double minptcut) {
	
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
