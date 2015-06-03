SHELL = /bin/bash
CFLAGS = $(shell root-config --cflags)
program_CXX_SRCS := $(wildcard ./src/*.cc)
program_TREE_SOS := $(wildcard ../src/*.so)
program_UNFOLD := $(wildcard ../../../RooUnfold-1.1.1/*.so)
program_FASTJET := $(wilcard /user/avanspil/fastjet-install/lib/*so)
program_CXX_OBJS := ${program_CXX_SRCS:.cc=.o}
program_OBJS := $(program_CXX_OBJS)
LIBDIR = ${PWD}/src
TREEDIR =${PWD}/../src
UNFOLDDIR =${PWD}/../../../RooUnfold-1.1.1

LIBHISTPAINTER := $(wildcard /cvmfs/cms.cern.ch/slc5_amd64_gcc434/cms/cmssw-patch/CMSSW_4_2_10_patch2/external/slc5_amd64_gcc434/lib/libHistPainter.so)

#FASTJETDIR = ${PWD}/Fastjet/fastjet-install
#CFLAGS+=`$(FASTJETDIR)/bin/fastjet-config --cxxflags`
#LDFLAGS+=`$(FASTJETDIR)/bin/fastjet-config --libs` -L$(FASTJETDIR)/lib -lfastjet

#FASTJETDIR = $(PYTHIA8)/Fastjet/fastjet-install/
#FASTJET = `$(FASTJETDIR)/bin/fastjet-config --cxxflags --libs` -L$(FASTJETDIR)/lib -lfastjet



program_INCLUDE_DIRS := ${TREEDIR}
program_LIBRARY_DIRS := ${TREEDIR}
program_UNFOLD_DIRS := ${UNFOLDDIR}
#program_FASTJET_DIRS := ${FASTJERDIR}

#CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir)) $(shell root-config --cflags)
CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir)) $(CFLAGS)
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(shell root-config --cflags --libs --glibs) -lPhysics -lThread -lMinuit -lHtml -lVMC -lEG -lGeom -Wl,-rpath -Wl,${LIBDIR} -Wl,-rpath -Wl,${TREEDIR},-rpath -Wl,${UNFOLDDIR}
#LDFLAGS += $(LDFLAGS) -lPhysics -lThread -lMinuit -lHtml -lVMC -lEG -lGeom -Wl,-rpath -Wl,${LIBDIR} -Wl,-rpath -Wl,${TREEDIR},-rpath -Wl,${UNFOLDDIR}

.PHONY: all clean distclean

all: Run

Plots: JER Systematics Corrections Unfolding

Run: Run.o $(program_OBJS) $(program_TREE_SOS)
	g++ -g Run.o -I IsolationCut.h $(program_OBJS) $(program_TREE_SOS) $(program_UNFOLD) -o Run ${LDFLAGS}

JER: MacroPlot_JER.o $(program_OBJS) $(program_TREE_SOS)
	g++ -g MacroPlot_JER.o -I color.h $(program_OBJS) $(program_TREE_SOS) $(program_UNFOLD) -o JER ${LDFLAGS}

JER_2D: MacroPlot_Matrix_JER.o $(program_OBJS) $(program_TREE_SOS)
	g++ -g MacroPlot_Matrix_JER.o -I color.h $(program_OBJS) $(program_TREE_SOS) $(program_UNFOLD) -o JER_2D ${LDFLAGS}

Response: MacroPlot_Matrix_ResponseMatrix_2by3.o $(program_OBJS) $(program_TREE_SOS)
	g++ -g MacroPlot_Matrix_ResponseMatrix_2by3.o -I color.h $(program_OBJS) $(program_TREE_SOS) $(program_UNFOLD) $(LIBHISTPAINTER) -o Response ${LDFLAGS}

Response_1by3: MacroPlot_Matrix_ResponseMatrix_1by3.o $(program_OBJS) $(program_TREE_SOS)
	g++ -g MacroPlot_Matrix_ResponseMatrix_1by3.o -I color.h $(program_OBJS) $(program_TREE_SOS) $(program_UNFOLD) $(LIBHISTPAINTER) -o Response_1by3 ${LDFLAGS}

Systematics: MacroPlot_CompareSystematics.o $(program_OBJS) $(program_TREE_SOS)
	g++ -g MacroPlot_CompareSystematics.o -I color.h $(program_OBJS) $(program_TREE_SOS) $(program_UNFOLD) -o Systematics ${LDFLAGS}

Corrections:  MacroPlot_CorrectionFactors.o $(program_OBJS) $(program_TREE_SOS)
	g++ -g MacroPlot_CorrectionFactors.o -I color.h $(program_OBJS) $(program_TREE_SOS) $(program_UNFOLD) -o Corrections ${LDFLAGS}

Unfolding:  MacroPlot_Unfold_energy.o $(program_OBJS) $(program_TREE_SOS) $(program_UNFOLD)
	g++ -g MacroPlot_Unfold_energy.o -I color.h $(program_OBJS) $(program_TREE_SOS) $(program_UNFOLD) -o Unfolding ${LDFLAGS}


Distances:  MacroPlot_JetDistances.o $(program_OBJS) $(program_TREE_SOS)
	g++ -g MacroPlot_JetDistances.o -I color.h $(program_OBJS) $(program_TREE_SOS) $(program_UNFOLD) -o Distances ${LDFLAGS}
	
JER_JES:  MacroPlot_JER_mean_spread_numbers.o $(program_OBJS) $(program_TREE_SOS)
	g++ -g MacroPlot_JER_mean_spread_numbers.o -I color.h $(program_OBJS) $(program_TREE_SOS) $(program_UNFOLD) -o JER_JES ${LDFLAGS}	
	
Calibration:  MacroPlot_Apply_Calibration.o $(program_OBJS) $(program_TREE_SOS)
	g++ -g MacroPlot_Apply_Calibration.o -I color.h $(program_OBJS) $(program_TREE_SOS) $(program_UNFOLD) -o Calibration ${LDFLAGS}		

Calibration_phi:  MacroPlot_Apply_Calibration_azimuthal.o $(program_OBJS) $(program_TREE_SOS)
	g++ -g MacroPlot_Apply_Calibration_azimuthal.o -I color.h $(program_OBJS) $(program_TREE_SOS) $(program_UNFOLD) -o Calibration_azimuthal ${LDFLAGS}

Calibration_function:  MacroPlot_Apply_Calibration_function.o $(program_OBJS) $(program_TREE_SOS)
	g++ -g MacroPlot_Apply_Calibration_function.o -I color.h $(program_OBJS) $(program_TREE_SOS) $(program_UNFOLD) -o Calibration_function ${LDFLAGS}

JER_Matrix:  MacroPlot_Matrix_JER_2by2.o $(program_OBJS) $(program_TREE_SOS)
	g++ -g MacroPlot_Matrix_JER_2by2.o -I color.h $(program_OBJS) $(program_TREE_SOS) $(LIBHISTPAINTER) $(program_UNFOLD) -o JER_Matrix ${LDFLAGS} 

ListOfFiles: MacroPlot_ListOfFiles.o $(program_OBJS) $(program_TREE_SOS)
	g++ -g  MacroPlot_ListOfFiles.o $(program_OBJS) $(program_TREE_SOS) $(LIBHISTPAINTER) $(program_UNFOLD) -o ListOfFiles ${LDFLAGS}

Edet:  MacroPlot_E_det.o $(program_OBJS) $(program_TREE_SOS)
	g++ -g MacroPlot_E_det.o -I color.h $(program_OBJS) $(program_TREE_SOS) $(LIBHISTPAINTER) $(program_UNFOLD) -o Edet ${LDFLAGS}

clean:
	@- $(RM) Run *.o 
	@- $(RM) $(program_OBJS) 
	@- $(RM) Calib*.o

distclean: clean
