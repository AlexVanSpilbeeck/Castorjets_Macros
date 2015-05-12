// Code by  Tom Cornelis.

#ifndef PLACETHRESHOLD_H
#define PLACETHRESHOLD_H

#include <vector>

bool PlaceThreshold(MyGenJet jet, double Ethreshold){
 return( jet.Energy() > Ethreshold);
}

bool PlaceThreshold(MyCastorJet jet, double Ethreshold){
 return( jet.energy > Ethreshold);
}


#endif
