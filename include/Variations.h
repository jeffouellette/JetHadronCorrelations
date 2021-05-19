#ifndef __Variations_h__
#define __Variations_h__

#include <MyStyle.h>

#include <TString.h>

#include <vector>
#include <map>


//const int nVar = 7;
//std::vector <TString> variations = {"Nominal", "JetES5PercUpVar", "JetES5PercDownVar", "JetES5PercSmearVar", "JetES2PercUpVar", "JetES2PercDownVar", "JetES2PercSmearVar"};
const int nVar = 6;
std::vector <TString> variations = {"Nominal", "JetES5PercUpVar", "JetES5PercDownVar", "JetES2PercUpVar", "JetES2PercDownVar", "FcalCentVar"};
//const int nVar = 5;
//std::vector <TString> variations = {"Nominal", "JetES5PercUpVar", "JetES5PercDownVar", "JetES2PercUpVar", "JetES2PercDownVar"};
//const int nVar = 1;
//std::vector <TString> variations = {"Nominal"};


std::map <TString, MyStyle> varStyles = {
  {"JetES5PercUpVar",   MyStyle (kPink+5, 3)},
  {"JetES5PercDownVar", MyStyle (kPink+5, 2)},
  {"JetES5PercSmearVar",MyStyle (kGreen+2, 4)},
  {"JetES2PercUpVar",   MyStyle (kAzure+2, 3)},
  {"JetES2PercDownVar", MyStyle (kAzure+2, 2)},
  {"JetES2PercSmearVar",MyStyle (kOrange+7, 4)},
  {"FcalCentVar",       MyStyle (kViolet-5, 5)}
};

std::map <TString, TString> varFullNames = {
  {"JetES5PercUpVar",   "JES 5\% up"},
  {"JetES5PercDownVar", "JES 5\% down"},
  {"JetES5PercSmearVar","JES 5\% smear"},
  {"JetES2PercUpVar",   "JES 2\% up"},
  {"JetES2PercDownVar", "JES 2\% down"},
  {"JetES2PercSmearVar","JES 2\% smear"},
  {"FcalCentVar",       "FCal 0-20\%"}
};


#endif
