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


  "HITightVar",
  "HILooseVar",
  "TrkEffVar",
  "FakeRateVar",
  "PrimFitVar",
  "JetPrimFracVar",
  "PartSpcVar",

  //"FcalCentVar",
  //"FineFcalCentVar",


  "MixCatVar1",
  "MixCatVar2",
  "MixCatVar3",

  //"NonClosureVar",

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
  "JESVar11",
  "JESVar12",
  "JESVar13",
  "JESVar14",
  "JESVar15",
  //"JESVar16",
  "JESVar17",
  "JESVar18",
  "JESVar19"
  //"JESVar20",


  //"MCTruthJetsTruthParts",     // Truth-level jets and truth-level charged particles
  //"MCRecoJetsTruthParts",       // Reco-level jets and truth-level charged particles
  //"MCRecoJetsTruthMatchedParts",  // Reco-level jets and reco-level, truth-matched charged particles
};
const int nVar = (int)variations.size ();

int GetVariationN (const TString& s) {
  int i = 0;
  while (i < nVar && variations[i] != s) i++;
  return i;
}


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

  "HITightVar",
  "HILooseVar",
  "TrkEffVar",
  "FakeRateVar",
  "PrimFitVar",
  "JetPrimFracVar",
  "PartSpcVar",

  "MixCatVar1",
  "MixCatVar2",
  "MixCatVar3",

  "NonClosureVar",
};


// other variations in MC which are not to be considered as systematics
std::set <TString> otherMCVariations = {
  "MCTruthJetsTruthParts",
  "MCRecoJetsTruthParts",
  "MCRecoJetsTruthMatchedParts",
};

// these variations don't need a background subtraction in pp
std::set <TString> variationsWithNoppBkgd = {
  "MCTruthJetsTruthParts",
  "MCRecoJetsTruthParts",
  "MCRecoJetsTruthMatchedParts",
};

// these variations don't need a background subtraction in p+Pb
std::set <TString> variationsWithNopPbBkgd = {
  "MCTruthJetsTruthParts",
  "MCRecoJetsTruthParts",
  "MCRecoJetsTruthMatchedParts",
};

// these variations don't need an unfold performed on them
std::set <TString> variationsWithNoUnfold = {
  "MCTruthJetsTruthParts",
};

// these variations don't cancel in finding the signal yield
std::set <TString> variationsThatDontCancelInSig = {
  "PartSpcVar",
};

// these variations don't cancel in p+Pb/pp ratio
std::set <TString> variationsThatDontCancelInRatio = {
  "JESVar19",
};

// for each group of variations below take the maximum from that group
std::vector <std::vector <TString>> variationGroups = {
  {"HITightVar", "HILooseVar"},
  {"TrkEffVar"},
  {"FakeRateVar"},
  {"PrimFitVar"},
  {"JetPrimFracVar"},
  {"PartSpcVar"},

  //{"FcalCentVar"},
  //{"FineFcalCentVar"},

  {"MixCatVar1", "MixCatVar3"},
  {"MixCatVar2"},
  {"NonClosureVar"},

  {"JESVar0"},
  {"JESVar1"},
  {"JESVar2"},
  {"JESVar3"},
  {"JESVar4"},
  {"JESVar5"},
  {"JESVar6"},
  {"JESVar7"},
  {"JESVar8"},
  {"JESVar9"},
  {"JESVar10"},
  {"JESVar11"},
  {"JESVar12"},
  {"JESVar13"},
  {"JESVar14"},
  {"JESVar15"},
  //{"JESVar16"},
  {"JESVar17"},
  {"JESVar18"},
  {"JESVar19"},
  //{"JESVar20"},
};


std::map <TString, MyStyle> varStyles = {
  {"HITightVar",        MyStyle (myLiteRed,       4)},
  {"HILooseVar",        MyStyle (myCyan,          4)},
  {"TrkEffVar",         MyStyle (myLitePurple,    4)},
  {"FakeRateVar",       MyStyle (myLiteYellow,    4)},
  {"PrimFitVar",        MyStyle (myLiteBlue,      4)},
  {"JetPrimFracVar",    MyStyle (kRed,            4)},
  {"PartSpcVar",        MyStyle (myLiteGreen,     4)},

  {"FcalCentVar",       MyStyle (kViolet-5,       5)},
  {"FineFcalCentVar",   MyStyle (kViolet,         5)},

  {"MixCatVar1",        MyStyle (myRed,           7)},
  {"MixCatVar2",        MyStyle (myGreen,         7)},
  {"MixCatVar3",        MyStyle (kViolet-3,       7)},
  {"NonClosureVar",     MyStyle (myOrange,        7)},

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
  //{"JESVar16",          MyStyle (kWhite,          2)},
  {"JESVar17",          MyStyle (manyColors[17],  2)},
  {"JESVar18",          MyStyle (manyColors[18],  2)},
  {"JESVar19",          MyStyle (manyColors[19],  2)},
  //{"JESVar20",          MyStyle (kWhite,          2)},

  //{"NoFcalMixCatVar",   MyStyle (myLiteRed, 6)},
  //{"pPbFcalMixCatVar",  MyStyle (myCyan, 6)},
  //{"ppFcalMixCatVar",   MyStyle (myLitePurple, 6)},

  {"Mixing",            MyStyle (myBlue,      1)},
  {"Tracking",          MyStyle (myViolet,    1)},
  {"Jets",              MyStyle (myMaroon,    1)},
};


// for plot labeling -- what is this systematic called?
std::map <TString, TString> varFullNames = {
  {"HITightVar",        "HITight tracks"},
  {"HILooseVar",        "HILoose tracks"},
  {"TrkEffVar",         "Tracking efficiency"},
  {"FakeRateVar",       "Fake rate"},
  {"PrimFitVar",        "Primary fraction fit"},
  {"JetPrimFracVar",    "Jet vs. MB primaries"},
  {"PartSpcVar",        "Particle species"},

  {"FcalCentVar",       "FCal-based centrality"},
  {"FineFcalCentVar",   "Fine FCal centrality"},

  {"MixCatVar1",        "Mixing variation 1"},
  {"MixCatVar2",        "Mixing variation 2"},
  {"MixCatVar3",        "Mixing variation 3"},
  {"NonClosureVar",     "MC non-closure (50\%)"},

  {"FcalCentVar",       "FCal 0-20\%"},
  {"NoFcalMixCatVar",   "No FCal Matching"},
  {"pPbFcalMixCatVar",  "FCal Matching (#it{p}+Pb only)"},
  {"ppFcalMixCatVar",   "FCal Matching (#it{pp} only)"},

  {"MixCatVar1",        "Mixing variation 1"},
  {"MixCatVar2",        "Mixing variation 2"},
  {"MixCatVar3",        "Mixing variation 3"},

  {"JESVar0",           "EtaIntercalibration_Modelling"},
  {"JESVar1",           "EtaIntercalibration_TotalStat"},
  {"JESVar2",           "EtaIntercalibration_NonClosure_highE"},
  {"JESVar3",           "EtaIntercalibration_NonClosure_negEta"},
  {"JESVar4",           "EtaIntercalibration_NonClosure_posEta"},
  {"JESVar5",           "PunchThrough_MC16"},
  {"JESVar6",           "EffectiveNP_1"},
  {"JESVar7",           "EffectiveNP_2"},
  {"JESVar8",           "EffectiveNP_3"},
  {"JESVar9",           "EffectiveNP_4"},
  {"JESVar10",          "EffectiveNP_5"},
  {"JESVar11",          "EffectiveNP_6"},
  {"JESVar12",          "EffectiveNP_7"},
  {"JESVar13",          "EffectiveNP_8restTerm"},
  {"JESVar14",          "SingleParticle_HighPt"},
  {"JESVar15",          "RelativeNonClosure_MC16"},
  //{"JESVar16",          "EtaIntercalibration_NonClosure_2018data"},
  {"JESVar17",          "Cross-calibration"},
  {"JESVar18",          "Flavor-dep. response"},
  {"JESVar19",          "Flavor fraction"},
  //{"JESVar20",          "R-scan"},

  {"Mixing",            "#bf{Total mixing}"},
  {"Tracking",          "#bf{Total tracking}"},
  {"Jets",              "#bf{Total jets}"}
};



bool IsMixingVariation (const TString& s) {
  return s.Contains ("Mix") || s.Contains ("NonClosureVar");
}


bool IsTrackingVariation (const TString& s) {
  return (dataVariations.count (s) > 0 && !IsMixingVariation (s));
}


bool IsJetsVariation (const TString& s) {
  return (s.Contains ("JES"));
}


std::vector <TString> totalVariations = {
  "Mixing",
  "Tracking",
  "Jets"
};


short GetVarN (const TString& s) {
  int iVar = 0;
  while (iVar < nVar && strcmp (variations[iVar], s.Data ()) != 0) iVar++;
  if (iVar == nVar) {
    std::cout << "Cannot find " << s << " result? Please check!" << std::endl;
    return -1;
  }
  return iVar;
}


TString GetVarBkg (const TString& s) {
  if (s == "JESVar0" ||
      s == "JESVar1" ||
      s == "JESVar2" ||
      s == "JESVar3" ||
      s == "JESVar4" ||
      s == "JESVar5" ||
      s == "JESVar6" ||
      s == "JESVar7" ||
      s == "JESVar8" ||
      s == "JESVar9" ||
      s == "JESVar10" ||
      s == "JESVar11" ||
      s == "JESVar12" ||
      s == "JESVar13" ||
      s == "JESVar14" ||
      s == "JESVar15" ||
      s == "JESVar16" ||
      s == "JESVar17" ||
      s == "JESVar18" ||
      s == "JESVar19" ||
      s == "JetPrimFracVar") {
    return "Nominal";
  } else {
    return s; 
  }
}


#endif
