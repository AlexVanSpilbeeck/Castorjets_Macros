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

                        if (strcmp(argv[3],"JetAnalyzer_radii") == 0) {
                                std::cout << "We'll process the data with the JetAnalyzer" << std::endl;
                                m->makeJetHistos_radii("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/avanspil/"
							  "/MinimumBias/2015_01_CastorTree_data_MinBias_Commissioning10-May19ReReco-v1_7TeV_53XRECOwithCorrector_v1/032c1e3d0d409c6175040143a103f571/",
                                                          "CastorTree_data_MinBias_Commissioning10-May19ReReco-v1_7TeV_53XRECOwithCorrector_",true,"",
                                                          argv[4],
                                                          argv[5],
                                                          atoi(argv[6]),
                                                          argv[7]);
                        }
                        if (strcmp(argv[3],"JetAnalyzer_radii_strippedTree") == 0) {
                                std::cout << "We'll process the Data tree now with the JetAnalyzer" << std::endl;
                                m->makeJetHistos_radii_strippedTree("/user/avanspil/Castor_Analysis/CMSSW_4_2_10_patch2/src/UACastor/CastorTree/Analysis/LoopRootFiles/",                               
                                                          true,"",
                                                          atoi(argv[4]),// N events
                                                          argv[5],      // label/date
                                                          argv[6],      // Name of inpu file
                                                          argv[7],      // had or em
                                                          //atoi(argv[8]),// Castor sectors. 1 (one only) - 2 (2 or 3) - 0 (any number)
                                                          atof( argv[8]), // Energy threshold.
                                                          argv[9] );      // setup
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
                                std::cout << "We'll process the Pythia6 Z2star 7TeV MC tree now with the JetAnalyzer (new)" << std::endl;
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
/*
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
 */  
                        if (strcmp(argv[3],"JetAnalyzer_radii") == 0) {
                                std::cout << "We'll process the Pythia6 Z2star 7TeV MC tree now with the JetAnalyzer" << std::endl;
                                m->makeJetHistos_radii("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/avanspil/"
                                                          "/MinBias_TuneZ2star_HFshowerLibrary_7TeV_pythia6/"
                                                          "CastorTree_MC_MinBias_TuneZ2star_7TeV_pythia6_LowPU2010_53XRECOwithCorrector_v3/45ff2e600e4cec01b439799f3c950bc6/",
                                                          "CastorTree_MC_7TeV_42X_53XRECOwithCorrector_",false,"Pythia6_Z2star_Default_varyR_", 
							  argv[4], 
							  argv[5], 
							  atoi(argv[6]),
							  argv[7]);
                        }

                }
        }

	// MC without a pt cut on gen level.
//pnfs/iihe/cms/store/user/avanspil/MinBias_TuneZ2star_HFshowerLibrary_7TeV_pythia6/CastorTree_MC_MinBias_TuneZ2star_7TeV_pythia6_LowPU2010_53XRECOwithCorrector_no_gen_pTcut/56a77bcc5e2daa0d1a676919111916e5/CastorTree_MC_7TeV_42X_53XRECOwithCorrector_110_1_ZP8.root
        if (strcmp(argv[1],"Pythia6Z2star_noPtCut")==0) {
                if (strcmp(argv[2],"7000")==0) {	
                        if (strcmp(argv[3],"JetAnalyzer_radii") == 0) {
                                std::cout << "We'll process the Pythia6 Z2star 7TeV MC tree now with the JetAnalyzer" << std::endl;
                                m->makeJetHistos_radii("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/avanspil/"
                                                          "/MinBias_TuneZ2star_HFshowerLibrary_7TeV_pythia6/"
                                                          "CastorTree_MC_MinBias_TuneZ2star_7TeV_pythia6_LowPU2010_53XRECOwithCorrector_no_gen_pTcut/56a77bcc5e2daa0d1a676919111916e5/",
                                                          "CastorTree_MC_7TeV_42X_53XRECOwithCorrector_",false,"Pythia6_Z2star_Default_varyR_",
                                                          argv[4],
                                                          argv[5],
                                                          atoi(argv[6]),
                                                          argv[7]);
                        }
                }
        }	

	//////////////////////////////////////////////////////////////////////////
	// Make a stripped down tree: only good events with only good gen jets. //
	//////////////////////////////////////////////////////////////////////////
        if (strcmp(argv[1],"Pythia6Z2star_noPtCut")==0) {
                if (strcmp(argv[2],"7000")==0) {	
                        if (strcmp(argv[3],"JetAnalyzer_stripTheTree") == 0) {
                                std::cout << "We'll process the Pythia6 Z2star 7TeV MC tree now with the JetAnalyzer" << std::endl;
/*
                                m->makeJetHistos_stripTheTree("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/avanspil/"
                                                          "/MinBias_TuneZ2star_HFshowerLibrary_7TeV_pythia6/"
                                                          "CastorTree_MC_MinBias_TuneZ2star_7TeV_pythia6_LowPU2010_53XRECOwithCorrector_no_gen_pTcut/56a77bcc5e2daa0d1a676919111916e5/",
                                                          "CastorTree_MC_7TeV_42X_53XRECOwithCorrector_",false,"Pythia6_Z2star_Default_varyR_",
*/
                                m->makeJetHistos_stripTheTree("/user/avanspil/public/temp/MinBias_TuneZ2star_HFshowerLibrary_7TeV_pythia6/",			  
							  "CastorTree_MC_7TeV_42X_53XRECOwithCorrector_",
							  false,
                                                          "",
                                                          argv[4],	// GEN radius
                                                          argv[5],	// Det radius
                                                          atoi(argv[6]),// N events
                                                          argv[7],	// Label
							  atoi(argv[8])); // file to start from.
                        }
                }
        }	
        //////////////////////////////////////////////////////////
        // Investigate the tree created with the command above. //
        //////////////////////////////////////////////////////////

        if (strcmp(argv[1],"Pythia6Z2star_noPtCut")==0) {
                if (strcmp(argv[2],"7000")==0) {
                        if (strcmp(argv[3],"JetAnalyzer_radii_strippedTree") == 0) {
			  std::cout << "\t--\tYou have chosen setup:-" << argv[9] << "-" <<  std::endl;

			  if( strcmp(argv[9],"raw") != 0 && strcmp(argv[9],"raw_wide") != 0 && strcmp(argv[9],"isolated") != 0 && strcmp(argv[9],"calibrated") != 0 && strcmp(argv[9],"unfold") != 0    ) {			         std::cout << "Can not compute" << std::endl;
			    std::cout << "Please choose one of the following options\n" <<
					 "\traw_wide\t" << "no isolation, generator jets in all of CASTOR's eta range, no calibration\n" <<
					 "\traw\t\t"	<< "no isolation, generator jets completely contained, no calibration\n" <<
					 "\tisolated\t"	<< "isolated generator jets completely contained by CASTOR\n" <<
					 "\tcalibrated\t" << "isolated generator jets completely contained by CASTOR, calibrated detector level jets\n" <<
					 "\tunfold\t\t"	<< "preparation of RooUnfoldResponse object" << std::endl;
					 
			  }
			  else{
                                std::cout << "We'll process the Pythia6 Z2star 7TeV MC tree now with the JetAnalyzer" << std::endl;
                                m->makeJetHistos_radii_strippedTree("/user/avanspil/Castor_Analysis/CMSSW_4_2_10_patch2/src/UACastor/CastorTree/Analysis/LoopRootFiles/",
//                                m->makeJetHistos_radii_strippedTree("/user/avanspil/public/temp/MinBias_TuneZ2star_HFshowerLibrary_7TeV_pythia6/",                                
                                                          false,"",
                                                          atoi(argv[4]),// N events
                                                          argv[5],	// label/date
							  argv[6],	// Name of inpu file
							  argv[7], 	// had or em
							  // atoi(argv[8]),// Castor sectors. 1 (one only) - 2 (2 or 3) - 0 (any number)
							  atof(argv[8]), // Energy threshold 
							  argv[9]);   // Setup
			  }
                        }
                }
       }
	


	delete m;
	
	return(0);
}

