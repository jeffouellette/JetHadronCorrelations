#ifndef __Variations_h__
#define __Variations_h__

#include <MyColors.h>
#include <MyStyle.h>

#include <TString.h>

#include <vector>
#include <map>


//const int nVar = 7;
//std::vector <TString> variations = {"Nominal", "JetES5PercUpVar", "JetES5PercDownVar", "JetES5PercSmearVar", "JetES2PercUpVar", "JetES2PercDownVar", "JetES2PercSmearVar"};
//const int nVar = 9;
//std::vector <TString> variations = {"Nominal", "JetES5PercUpVar", "JetES5PercDownVar", "JetES2PercUpVar", "JetES2PercDownVar", "FcalCentVar", "FcalMixCatVar", "pPbFcalMixCatVar", "ppFcalMixCatVar"};
const int nVar = 10;
std::vector <TString> variations = {"Nominal", "JetES5PercUpVar", "JetES5PercDownVar", "FcalCentVar", "NoFcalMixCatVar", "pPbFcalMixCatVar", "ppFcalMixCatVar", "MixCatVar1", "MixCatVar2", "MixCatVar3"};
//const int nVar = 5;
//std::vector <TString> variations = {"Nominal", "JetES5PercUpVar", "JetES5PercDownVar", "JetES2PercUpVar", "JetES2PercDownVar"};
//const int nVar = 1;
//std::vector <TString> variations = {"Nominal"};


std::map <TString, MyStyle> varStyles = {
  {"JetES5PercUpVar",   MyStyle (myLiteGreen, 3)},
  {"JetES5PercDownVar", MyStyle (myLiteGreen, 2)},
  {"JetES5PercSmearVar",MyStyle (kGreen+2, 4)},
  {"JetES2PercUpVar",   MyStyle (kAzure+2, 3)},
  {"JetES2PercDownVar", MyStyle (kAzure+2, 2)},
  {"JetES2PercSmearVar",MyStyle (kOrange+7, 4)},
  {"FcalCentVar",       MyStyle (kViolet-5, 5)},
  {"NoFcalMixCatVar",   MyStyle (myLiteRed, 6)},
  {"pPbFcalMixCatVar",  MyStyle (myCyan, 6)},
  {"ppFcalMixCatVar",   MyStyle (myLitePurple, 6)},
  {"MixCatVar1",        MyStyle (myOrange, 7)},
  {"MixCatVar2",        MyStyle (kAzure+2, 7)},
  {"MixCatVar3",        MyStyle (kGreen+2, 7)},
};

std::map <TString, TString> varFullNames = {
  {"JetES5PercUpVar",   "JES 5\% up"},
  {"JetES5PercDownVar", "JES 5\% down"},
  {"JetES5PercSmearVar","JES 5\% smear"},
  {"JetES2PercUpVar",   "JES 2\% up"},
  {"JetES2PercDownVar", "JES 2\% down"},
  {"JetES2PercSmearVar","JES 2\% smear"},
  {"FcalCentVar",       "FCal 0-20\%"},
  {"NoFcalMixCatVar",   "No FCal Matching"},
  {"pPbFcalMixCatVar",  "FCal Matching (#it{p}+Pb only)"},
  {"ppFcalMixCatVar",   "FCal Matching (#it{pp} only)"},
  {"MixCatVar1",        "Mixing variation 1"},
  {"MixCatVar2",        "Mixing variation 2"},
  {"MixCatVar3",        "Mixing variation 3"},
};


#endif
