#ifndef __Variations_h__
#define __Variations_h__

#include <MyColors.h>
#include <MyStyle.h>

#include <TString.h>

#include <vector>
#include <map>
#include <set>


std::vector <TString> variations = {
  "Nominal",
//  "HITightVar",
//  "HILooseVar",
  "TrkEffVar",
  "FakeRateVar",
  "PrimFitVar",
  "PartSpcVar",
  //"FcalCentVar",
  //"FineFcalCentVar",
  "MixCatVar1",
  "MixCatVar2",
  "MixCatVar3",
  "JESVar0",
  "JESVar1",
  "JESVar2",
  "JESVar3",
  "JESVar4",
  "JESVar5",
  "JESVar6",
  "JESVar7",
  "JESVar8",
  "JESVar9",
  "JESVar10",
  "JESVar10",
  "JESVar11",
  "JESVar12",
  "JESVar13",
  "JESVar14",
  "JESVar15",
  //"JESVar16",
  "JESVar17",
  "JESVar18",
  "JESVar19",
  //"JESVar20",
};
const int nVar = (int)variations.size ();

//const int nVar = 1;
//std::vector <TString> variations = {"Nominal"};


// variations to consider in MC -- skip if not listed here
std::set <TString> mcVariations = {
  "Nominal",
  "JESVar0",
  "JESVar1",
  "JESVar2",
  "JESVar3",
  "JESVar4",
  "JESVar5",
  "JESVar6",
  "JESVar7",
  "JESVar8",
  "JESVar9",
  "JESVar10",
  "JESVar10",
  "JESVar11",
  "JESVar12",
  "JESVar13",
  "JESVar14",
  "JESVar15",
  //"JESVar16",
  "JESVar17",
  "JESVar18",
  "JESVar19",
  //"JESVar20",
};


// variations to consider in data -- skip if not listed here
std::set <TString> dataVariations = {
  "Nominal",
//  "HITightVar",
//  "HILooseVar",
  "TrkEffVar",
  "FakeRateVar",
  "PrimFitVar",
  "PartSpcVar",
  "MixCatVar1",
  "MixCatVar2",
  "MixCatVar3"
};


std::map <TString, MyStyle> varStyles = {
  {"HITightVar",        MyStyle (myLiteRed,   3)},
  {"HILooseVar",        MyStyle (myCyan,      3)},
  {"TrkEffVar",         MyStyle (myLitePurple,3)},
  {"FakeRateVar",       MyStyle (myOrange,    3)},
  {"PrimFitVar",        MyStyle (myBlue,      3)},
  {"PartSpcVar",        MyStyle (myLiteGreen, 3)},
  {"FcalCentVar",       MyStyle (kViolet-5, 5)},
  {"FineFcalCentVar",   MyStyle (kViolet,   5)},
  {"MixCatVar1",        MyStyle (myOrange,        7)},
  {"MixCatVar2",        MyStyle (kAzure+2,        7)},
  {"MixCatVar3",        MyStyle (kGreen+2,        7)},
  {"JESVar0",           MyStyle (manyColors[0],   2)},
  {"JESVar1",           MyStyle (manyColors[1],   2)},
  {"JESVar2",           MyStyle (manyColors[2],   2)},
  {"JESVar3",           MyStyle (manyColors[3],   2)},
  {"JESVar4",           MyStyle (manyColors[4],   2)},
  {"JESVar5",           MyStyle (manyColors[5],   2)},
  {"JESVar6",           MyStyle (manyColors[6],   2)},
  {"JESVar7",           MyStyle (manyColors[7],   2)},
  {"JESVar8",           MyStyle (manyColors[8],   2)},
  {"JESVar9",           MyStyle (manyColors[9],   2)},
  {"JESVar10",          MyStyle (manyColors[10],  2)},
  {"JESVar10",          MyStyle (manyColors[11],  2)},
  {"JESVar11",          MyStyle (manyColors[12],  2)},
  {"JESVar12",          MyStyle (manyColors[13],  2)},
  {"JESVar13",          MyStyle (manyColors[14],  2)},
  {"JESVar14",          MyStyle (manyColors[15],  2)},
  {"JESVar15",          MyStyle (manyColors[16],  2)},
  {"JESVar16",          MyStyle (kWhite,          2)},
  {"JESVar17",          MyStyle (manyColors[17],  2)},
  {"JESVar18",          MyStyle (manyColors[18],  2)},
  {"JESVar19",          MyStyle (manyColors[19],  2)},
  {"JESVar20",          MyStyle (kWhite,          2)},

  //{"NoFcalMixCatVar",   MyStyle (myLiteRed, 6)},
  //{"pPbFcalMixCatVar",  MyStyle (myCyan, 6)},
  //{"ppFcalMixCatVar",   MyStyle (myLitePurple, 6)},
};


std::map <TString, TString> varFullNames = {
  {"HITightVar",        "HITight tracks"},
  {"HILooseVar",        "HILoose tracks"},
  {"TrkEffVar",         "Tracking efficiency"},
  {"FakeRateVar",       "Fake rate"},
  {"PrimFitVar",        "Primary fraction fit"},
  {"PartSpcVar",        "Particle species"},
  {"FcalCentVar",       "FCal-based centrality"},
  {"FineFcalCentVar",   "Fine FCal centrality"},
  {"MixCatVar1",        "Mixing variation 1"},
  {"MixCatVar2",        "Mixing variation 2"},
  {"MixCatVar3",        "Mixing variation 3"},
  {"FcalCentVar",       "FCal 0-20\%"},
  {"NoFcalMixCatVar",   "No FCal Matching"},
  {"pPbFcalMixCatVar",  "FCal Matching (#it{p}+Pb only)"},
  {"ppFcalMixCatVar",   "FCal Matching (#it{pp} only)"},
  {"MixCatVar1",        "Mixing variation 1"},
  {"MixCatVar2",        "Mixing variation 2"},
  {"MixCatVar3",        "Mixing variation 3"},
  {"JESVar0",           "JET_EtaIntercalibration_Modelling"},
  {"JESVar1",           "JET_EtaIntercalibration_TotalStat"},
  {"JESVar2",           "JET_EtaIntercalibration_NonClosure_highE"},
  {"JESVar3",           "JET_EtaIntercalibration_NonClosure_negEta"},
  {"JESVar4",           "JET_EtaIntercalibration_NonClosure_posEta"},
  {"JESVar5",           "JET_PunchThrough_MC16"},
  {"JESVar6",           "JET_EffectiveNP_1"},
  {"JESVar7",           "JET_EffectiveNP_2"},
  {"JESVar8",           "JET_EffectiveNP_3"},
  {"JESVar9",           "JET_EffectiveNP_4"},
  {"JESVar10",          "JET_EffectiveNP_5"},
  {"JESVar11",          "JET_EffectiveNP_6"},
  {"JESVar12",          "JET_EffectiveNP_7"},
  {"JESVar13",          "JET_EffectiveNP_8restTerm"},
  {"JESVar14",          "JET_SingleParticle_HighPt"},
  {"JESVar15",          "JET_RelativeNonClosure_MC16"},
  {"JESVar16",          "JET_EtaIntercalibration_NonClosure_2018data"},
  {"JESVar17",          "Cross-calibration"},
  {"JESVar18",          "Flavor-dep. response"},
  {"JESVar19",          "Flavor fraction"},
  {"JESVar20",          "R-scan"},
};


#endif
