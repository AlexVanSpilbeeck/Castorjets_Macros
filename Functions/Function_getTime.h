
#ifndef FUNCTION_GETTIME_H
#define FUNCTION_GETTIME_H

#include <ctime>
#include <TString.h>

TString getTime(){

   //-- Prepare the directory to store the plots.
   time_t theTime = time(NULL);
   struct tm *aTime = localtime(&theTime);

   int day = aTime->tm_mday;
   int month = aTime->tm_mon + 1; // Month is 0 - 11, add 1 to get a jan-dec 1-12 concept
   int year = aTime->tm_year + 1900;

   TString datum = TString::Format( "%i", year );
   if( month < 10 ) datum += "0";
   datum += TString::Format( "%i", month );;
   if( day < 10 ) datum += "0";
   datum += TString::Format( "%i", day );;

   return datum;


}

#endif
