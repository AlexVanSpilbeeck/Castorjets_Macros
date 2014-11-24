#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <stdio.h>
#include <cstring>

#include "./src/MainAnalyzer.h"

int main(int argc, char *argv[])
{
	
	//int i;
	//printf("argc = %d\n", argc);
	//for (int i = 0; i<argc; i++) printf("argv[%d] = %s\n", i, argv[i]);
	
	MainAnalyzer* m = new MainAnalyzer();
	
	///////////////////////
	////// DATA ///////////
	///////////////////////
	
	// default data running
	if (strcmp(argv[1],"data") == 0) {
		if (strcmp(argv[2],"7000") == 0) {
			// 7 TeV data
			if (strcmp(argv[3],"JetAnalyzer") == 0) {
				std::cout << "We'll process the 7TeV data tree now with the JetAnalyzer" << std::endl;
				m->makeJetHistos("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/hvanhaev/"
							  //"/MinimumBias/CastorTree_data_MinBias_Commissioning10-May19ReReco-v1_7TeV_53XRECOwithCorrector_v2/99f6c87cf56e3a4a62f56364fdfe2385/",
							  "MinimumBias/CastorTree_data_MinBias_Commissioning10-May19ReReco-v1_7TeV_53XRECOwithCorrector_v3/97203ab154b3c9bdada743fe0880f7f7/",
							  "CastorTree_data_MinBias_Commissioning10-May19ReReco-v1_7TeV_53XRECOwithCorrector_",true,"");
			}
		}	
	}
		
    
	
	///////////////////////////////////////
	//////  MC Default + models ///////////
	///////////////////////////////////////
	
	// default Pythia6 Z2star running
	if (strcmp(argv[1],"Pythia6Z2star")==0) {
		if (strcmp(argv[2],"7000")==0) {
			// 7 TeV Pythia6 Z2star
			if (strcmp(argv[3],"JetAnalyzer") == 0) {
				std::cout << "We'll process the Pythia6 Z2star 7TeV MC tree now with the JetAnalyzer" << std::endl;
				m->makeJetHistos("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/hvanhaev/"
							  "/MinBias_TuneZ2star_HFshowerLibrary_7TeV_pythia6/"
							  "CastorTree_MC_MinBias_TuneZ2star_7TeV_pythia6_LowPU2010_53XRECOwithCorrector_v2/632755d48b0d40f274e8ee11d45c259a/",
							  "CastorTree_MC_7TeV_42X_53XRECOwithCorrector_",false,"Pythia6_Z2star_Default_");
			}
		}
	}
	
	
	
	/////////////////////////////////////////
	//////  MC with different akR ///////////
	/////////////////////////////////////////
	
        if (strcmp(argv[1],"Pythia6Z2star_diffR")==0) {
                if (strcmp(argv[2],"7000")==0) {
                        if (strcmp(argv[3],"JetAnalyzer") == 0) {
                                std::cout << "We'll process the Pythia6 Z2star 7TeV MC tree now with the JetAnalyzer" << std::endl;
                                m->makeJetHistos("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/hvanhaev/"
                                                          "/MinBias_TuneZ2star_HFshowerLibrary_7TeV_pythia6/"
                                                          "CastorTree_MC_MinBias_TuneZ2star_7TeV_pythia6_LowPU2010_53XRECOwithCorrector_v3/d5197a65c5d23484e969be9804f84888/",
                                                          "CastorTree_MC_7TeV_42X_53XRECOwithCorrector_",false,"Pythia6_Z2star_Default_varyR_");
                        }
                }
        }	

        ////////////////////////////////////////////
        //////  My MC with different akR ///////////
        ////////////////////////////////////////////

        if (strcmp(argv[1],"Pythia6Z2star_new")==0) {
                if (strcmp(argv[2],"7000")==0) {
                        if (strcmp(argv[3],"JetAnalyzer") == 0) {
                                std::cout << "We'll process the Pythia6 Z2star 7TeV MC tree now with the JetAnalyzer" << std::endl;
                                m->makeJetHistos("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/avanspil/"
                                                          "/MinBias_TuneZ2star_HFshowerLibrary_7TeV_pythia6/"
                                                          "CastorTree_MC_MinBias_TuneZ2star_7TeV_pythia6_LowPU2010_53XRECOwithCorrector_v3/45ff2e600e4cec01b439799f3c950bc6/",
                                                          "CastorTree_MC_7TeV_42X_53XRECOwithCorrector_",false,"Pythia6_Z2star_Default_varyR_");
                        }
			else if(strcmp(argv[3],"RadiusAnalyzer") == 0) {
                                std::cout << "We'll process the Pythia6 Z2star 7TeV MC tree now with the RadiusAnalyzer" << std::endl;
                                m->makeRadiusHistos("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/avanspil/"
                                                          "/MinBias_TuneZ2star_HFshowerLibrary_7TeV_pythia6/"
                                                          "CastorTree_MC_MinBias_TuneZ2star_7TeV_pythia6_LowPU2010_53XRECOwithCorrector_v3/45ff2e600e4cec01b439799f3c950bc6/",
                                                          "CastorTree_MC_7TeV_42X_53XRECOwithCorrector_",false,"Pythia6_Z2star_Default_varyR_");
                        }
                        else if(strcmp(argv[3],"SystematicsMin") == 0) {
                                std::cout << "We'll process the Pythia6 Z2star 7TeV MC tree now with the minimum systematics" << std::endl;
                                m->makeSysMinHistos("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/avanspil/"
                                                          "/MinBias_TuneZ2star_HFshowerLibrary_7TeV_pythia6/"
                                                          "CastorTree_MC_MinBias_TuneZ2star_7TeV_pythia6_LowPU2010_53XRECOwithCorrector_v3/45ff2e600e4cec01b439799f3c950bc6/",
                                                          "CastorTree_MC_7TeV_42X_53XRECOwithCorrector_",false,"Pythia6_Z2star_Default_varyR_");
                        }
                        else if(strcmp(argv[3],"SystematicsMax") == 0) {
                                std::cout << "We'll process the Pythia6 Z2star 7TeV MC tree now with the maximum systematics" << std::endl;
                                m->makeSysMaxHistos("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/avanspil/"
                                                          "/MinBias_TuneZ2star_HFshowerLibrary_7TeV_pythia6/"
                                                          "CastorTree_MC_MinBias_TuneZ2star_7TeV_pythia6_LowPU2010_53XRECOwithCorrector_v3/45ff2e600e4cec01b439799f3c950bc6/",
                                                          "CastorTree_MC_7TeV_42X_53XRECOwithCorrector_",false,"Pythia6_Z2star_Default_varyR_");
                        }

                }
        }

	
	delete m;
	
	return(0);
}

