#ifndef __Variations_h__
#define __Variations_h__

#include <MyColors.h>
#include <MyStyle.h>

#include <TString.h>

#include <vector>
#include <map>
#include <set>


const std::vector <TString> variations = {
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

  "MixCatVar4",
  "MixCatVar5",

  "UnfoldingAltWgtsVar",
  "NonClosureVar",
  "RefoldingVar",


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
  //"JESVar14",
  "JESVar15",
  "JESVar16",
  "JESVar17",
  "JESVar19",
  "JESVar18",
  //"JESVar20",


  //"MCTruthJetsTruthParts",     // Truth-level jets and truth-level charged particles
  //"MCRecoJetsTruthParts",       // Reco-level jets and truth-level charged particles
  //"MCRecoJetsTruthMatchedParts",  // Reco-level jets and reco-level, truth-matched charged particles
};
const short nVar = (short)variations.size ();

short GetVarN (const TString& s) {
  short iVar = 0;
  while (iVar < nVar && variations[iVar] != s) iVar++;
  if (iVar == nVar) {
    std::cout << "Cannot find " << s << " result? Please check!" << std::endl;
    return -1;
  }
  return iVar;
}


// variations to consider in MC -- skip if not listed here
const std::set <TString> mcVariations = {
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
  //"JESVar14",
  "JESVar15",
  "JESVar16",
  "JESVar17",
  "JESVar18",
  "JESVar19",
  //"JESVar20",
};


// variations to consider in data -- skip if not listed here
const std::set <TString> dataVariations = {
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
  "MixCatVar4",
  "MixCatVar5",

  "NonClosureVar",
  "UnfoldingAltWgtsVar",
  "RefoldingVar",
};


// other variations in MC which are not to be considered as systematics
const std::set <TString> otherMCVariations = {
  "MCTruthJetsTruthParts",
  "MCRecoJetsTruthParts",
  "MCRecoJetsTruthMatchedParts",
};

// these variations don't need a background subtraction in pp
const std::set <TString> variationsWithNoppBkgd = {
  "MCTruthJetsTruthParts",
  "MCRecoJetsTruthParts",
  "MCRecoJetsTruthMatchedParts",
};

// these variations don't need a background subtraction in p+Pb
const std::set <TString> variationsWithNopPbBkgd = {
  "MCTruthJetsTruthParts",
  "MCRecoJetsTruthParts",
  "MCRecoJetsTruthMatchedParts",
};

// these variations don't need an unfold performed on them
const std::set <TString> variationsWithNoUnfold = {
  "MCTruthJetsTruthParts",
};

// these variations don't cancel in finding the signal yield
const std::set <TString> variationsThatDontCancelInSig = {
  "PartSpcVar",
};

// these variations don't cancel in p+Pb/pp ratio
const std::set <TString> variationsThatDontCancelInRatio = {
  "JESVar19",
};

// thses variations use the jet pT unweighted response matrices
const std::set <TString> variationsWithUnwgtdRespMatrix {
  "UnfoldingAltWgtsVar",
};

const std::set <TString> variationsToSmooth {
  "HITightVar",
  "HILooseVar",
  "JESVar6",
  "JESVar17",
  "JESVar18",

  "JESVar0",
  "JESVar1",
  "JESVar2",
  "JESVar3",
  "JESVar4",
  "JESVar5",
  "JESVar7",
  "JESVar8",
  "JESVar9",
  "JESVar10",
  "JESVar11",
  "JESVar12",
  "JESVar13",
  "JESVar14",
  "JESVar15",
  "JESVar16",
};
const TString hermitePolyFunc = "[0] + [1]*log(x) + [2]*(pow(log(x),2)-1)"; // Hermite polynomials
//std::map <TString, TString> variationSmoothFuncs {
//  {"HITightVar", "ppcf"},
//  {"HILooseVar", "ppcf"},
//  {"JESVar6", hermitePolyFunc},
//  {"JESVar17", hermitePolyFunc},
//  {"JESVar18", hermitePolyFunc},
//};
  
//const TString smoothingFunc = "[0] + [1]*(1-log(x)) + [2]*(0.5*(pow(log(x),2)-4*log(x)+2)) + [3]*(-pow(log(x),3)+9*pow(log(x),2)-18*log(x)+6)/6 + [4]*(pow(log(x),4)-16*pow(log(x),3)+72*pow(log(x),2)-96*log(x)+24)/24"; // Laguerre polynomials
//const TString smoothingFunc = "[0] + [1]*log(x) + [2]*(pow(log(x),2)-1) + [3]*(pow(log(x),3)-log(x)) + [4]*(pow(log(x),4)-6*pow(log(x),2)+3)"; // Hermite polynomials

// for each group of variations below take the maximum from that group
const std::vector <std::vector <TString>> variationGroups = {
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
  {"MixCatVar4"},
  {"MixCatVar5"},
  {"NonClosureVar"},
  {"RefoldingVar"},
  {"UnfoldingAltWgtsVar"},

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
  //{"JESVar14"},
  {"JESVar15"},
  {"JESVar16"},
  {"JESVar17"},
  {"JESVar18"},
  {"JESVar19"},
  //{"JESVar20"},
};


std::map <TString, MyStyle> varStyles = {
  {"HITightVar",          MyStyle (manyColors[0],   4)},
  {"HILooseVar",          MyStyle (manyColors[2],   4)},
  {"TrkEffVar",           MyStyle (manyColors[4],   4)},
  {"FakeRateVar",         MyStyle (manyColors[6],   4)},
  {"PrimFitVar",          MyStyle (manyColors[8],   4)},
  {"JetPrimFracVar",      MyStyle (manyColors[10],  4)},
  {"PartSpcVar",          MyStyle (manyColors[12],  4)},

  {"MixCatVar1",          MyStyle (manyColors[1],   7)},
  {"MixCatVar2",          MyStyle (manyColors[3],   7)},
  {"MixCatVar3",          MyStyle (manyColors[5],   7)},
  {"MixCatVar4",          MyStyle (manyColors[7],   7)},
  {"MixCatVar5",          MyStyle (manyColors[9],   7)},
  {"UnfoldingAltWgtsVar", MyStyle (manyColors[12],  7)},
  {"NonClosureVar",       MyStyle (manyColors[13],  7)},
  {"RefoldingVar",        MyStyle (manyColors[15],  7)},

  {"FcalCentVar",         MyStyle (manyColors[1],   5)},
  {"FineFcalCentVar",     MyStyle (manyColors[3],   5)},

  {"JESVar0",             MyStyle (manyColors[0],   2)},
  {"JESVar1",             MyStyle (manyColors[1],   2)},
  {"JESVar2",             MyStyle (manyColors[2],   2)},
  {"JESVar3",             MyStyle (manyColors[3],   2)},
  {"JESVar4",             MyStyle (manyColors[4],   2)},
  {"JESVar5",             MyStyle (manyColors[5],   2)},
  {"JESVar6",             MyStyle (manyColors[6],   2)},
  {"JESVar7",             MyStyle (manyColors[7],   2)},
  {"JESVar8",             MyStyle (manyColors[8],   2)},
  {"JESVar9",             MyStyle (manyColors[9],   2)},
  {"JESVar10",            MyStyle (manyColors[10],  2)},
  {"JESVar10",            MyStyle (manyColors[11],  2)},
  {"JESVar11",            MyStyle (manyColors[12],  2)},
  {"JESVar12",            MyStyle (manyColors[13],  2)},
  {"JESVar13",            MyStyle (manyColors[14],  2)},
  //{"JESVar14",            MyStyle (kWhite,  2)},
  {"JESVar15",            MyStyle (manyColors[15],  2)},
  {"JESVar16",            MyStyle (manyColors[16],  2)},
  {"JESVar17",            MyStyle (manyColors[17],  2)},
  {"JESVar18",            MyStyle (manyColors[19],  2)},
  {"JESVar19",            MyStyle (manyColors[18],  2)},
  //{"JESVar20",            MyStyle (kWhite,          2)},

  {"Mixing",            MyStyle (myBlue,      1)},
  {"Tracking",          MyStyle (myViolet,    1)},
  {"Jets",              MyStyle (myMaroon,    1)},
};


// for plot labeling -- what is this systematic called?
std::map <TString, TString> varFullNames = {
  {"HITightVar",          "HITight tracks"},
  {"HILooseVar",          "HILoose tracks"},
  {"TrkEffVar",           "Tracking efficiency"},
  {"FakeRateVar",         "Fake rate"},
  {"PrimFitVar",          "Primary fraction fit"},
  {"JetPrimFracVar",      "Jet vs. MB primaries"},
  {"PartSpcVar",          "Particle species"},

  {"FcalCentVar",         "FCal-based centrality"},
  {"FineFcalCentVar",     "Fine FCal centrality"},

  {"NonClosureVar",       "MC non-closure"},
  {"UnfoldingAltWgtsVar", "Unfolding prior"},
  {"RefoldingVar",        "Refolding non-closure"},

  {"FcalCentVar",         "FCal 0-20\%"},

  {"MixCatVar1",          "Mixing variation 1"},
  {"MixCatVar2",          "Mixing variation 2"},
  {"MixCatVar3",          "Mixing variation 3"},
  {"MixCatVar4",          "Mixing variation 4"},
  {"MixCatVar5",          "Mixing variation 5"},

  {"JESVar0",             "EtaIntercalibration_Modelling"},
  {"JESVar1",             "EtaIntercalibration_TotalStat"},
  {"JESVar2",             "EtaIntercalibration_NonClosure_highE"},
  {"JESVar3",             "EtaIntercalibration_NonClosure_negEta"},
  {"JESVar4",             "EtaIntercalibration_NonClosure_posEta"},
  {"JESVar5",             "PunchThrough_MC16"},
  {"JESVar6",             "EffectiveNP_1"},
  {"JESVar7",             "EffectiveNP_2"},
  {"JESVar8",             "EffectiveNP_3"},
  {"JESVar9",             "EffectiveNP_4"},
  {"JESVar10",            "EffectiveNP_5"},
  {"JESVar11",            "EffectiveNP_6"},
  {"JESVar12",            "EffectiveNP_7"},
  {"JESVar13",            "EffectiveNP_8restTerm"},
  //{"JESVar14",            "Cross-calibration"},
  {"JESVar15",            "SingleParticle_HighPt"},
  {"JESVar16",            "BJES_Response"},
  {"JESVar17",            "Cross-calibration"},
  {"JESVar18",            "Flavor-dep. response"},
  {"JESVar19",            "Flavor fraction"},
  //{"JESVar20",            "R-scan"},

  {"Mixing",              "#bf{Total signal}"},
  {"Tracking",            "#bf{Total tracking}"},
  {"Jets",                "#bf{Total jets}"}
};



bool IsMixingVariation (const TString& s) {
  return (s.Contains ("Mix") || s.Contains ("NonClosureVar") || s.Contains ("RefoldingVar") || s.Contains ("UnfoldingAltWgtsVar"));
}


bool IsJetsVariation (const TString& s) {
  return (s.Contains ("JES"));
}


bool IsTrackingVariation (const TString& s) {
  return (!IsJetsVariation (s) && !IsMixingVariation (s));
}


const std::vector <TString> totalVariations = {
  "Tracking",
  "Jets",
  "Mixing",
};
const short nTotVar = (short)totalVariations.size ();


short GetTotVarN (const TString& s) {
  short iTotVar = 0;
  while (iTotVar < nTotVar && totalVariations[iTotVar] != s.Data ()) iTotVar++;
  if (iTotVar == nTotVar) {
    std::cout << "Cannot find " << s << " total variation result? Please check!" << std::endl;
    return -1;
  }
  return iTotVar;
}

TString GetTotVar (const TString& s) {
  if      (IsMixingVariation (s))   return TString ("Mixing");
  else if (IsJetsVariation (s))     return TString ("Jets");
  else if (IsTrackingVariation (s)) return TString ("Tracking");
  else return "???";
}


TString GetVarBkg (const TString& s) {
  if (//s == "JESVar0" ||
      //s == "JESVar1" ||
      //s == "JESVar2" ||
      //s == "JESVar3" ||
      //s == "JESVar4" ||
      //s == "JESVar5" ||
      //s == "JESVar6" ||
      //s == "JESVar7" ||
      //s == "JESVar8" ||
      //s == "JESVar9" ||
      //s == "JESVar10" ||
      //s == "JESVar11" ||
      //s == "JESVar12" ||
      //s == "JESVar13" ||
      //s == "JESVar14" ||
      //s == "JESVar15" ||
      //s == "JESVar16" ||
      //s == "JESVar17" ||
      //s == "JESVar18" ||
      //s == "JESVar19" ||
      s == "JetPrimFracVar" ||
      s == "UnfoldingAltWgtsVar") {
    return "Nominal";
  } else {
    return s; 
  }
}


#endif
