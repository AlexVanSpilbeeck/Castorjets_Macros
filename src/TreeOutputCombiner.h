/*
 *  TreeOutputCombiner.h
 *  
 *
 *  Created by Hans Van Haevermaet on 26/01/12.
 *  Copyright 2012 Universiteit Antwerpen. All rights reserved.
 *
 */


#ifndef TreeOutputCombiner_h
#define TreeOutputCombiner_h

#include <TString.h>
#include <TRegexp.h>
#include <TObjArray.h>
#include <TH1.h>
#include "../src/MyJet.h"
#include "../src/MyGenJet.h"
#include "../src/MyTrackJet.h"

class TreeOutputCombiner {
public:
	TreeOutputCombiner();
	virtual ~TreeOutputCombiner();
	void Combine(TString inputdir, TObjArray* filelist, double usedWeight);
	Double_t getRatioError(double	a, double b, double errora, double errorb);
	
private:
	
	
};

#endif


