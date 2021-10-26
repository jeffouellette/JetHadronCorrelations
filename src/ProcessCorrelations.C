#ifndef __JetHadronCorrelator_ProcessCorrelations_C__
#define __JetHadronCorrelator_ProcessCorrelations_C__

#include "CentralityDefs.h"
#include "Params.h"
#include "OutTree.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Variations.h"
#include "PrimaryFractionFit.h"

#include <ArrayTemplates.h>
#include <Utilities.h>

#include <TH1D.h>
#include <TH2D.h>

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif

#include <iostream>
#include <string.h>
#include <math.h>

using namespace JetHadronCorrelations;

PrimaryFractionFit bbbf;

short GetJetSpectraNIters (const short iCent) {
  switch (iCent) {
    case -1: return 1;
    case 0:  return 1;
    case 1:  return 1;
    case 2:  return 1;
    case 3:  return 1;
    case 4:  return 1;
    case 5:  return 1;
    default: return 0;
  }
}

short GetTrkSpectraNIters (const short iPtJInt, const short iDir, const short iCent) {
  switch (iPtJInt) {
    case 0: // pTJ > 30 GeV
      //return GetTrkSpectraNIters (1, iDir, iCent); // should default to same nIters for 30 GeV and 60 GeV jets since both use the same unfolded distribution.
    switch (iDir) {
      case 0: // near-side
      switch (iCent) {
        case -1: return 1;
        case 0:  return 3;
        case 1:  return 3;
        case 2:  return 3;
        case 3:  return 2;
        case 4:  return 1;
        case 5:  return 10;
        default: return 0;
      }
      case 1: // perpendicular
      switch (iCent) {
        case -1: return 1;
        case 0:  return 1;
        case 1:  return 1;
        case 2:  return 1;
        case 3:  return 1;
        case 4:  return 1;
        case 5:  return 1;
        default: return 0;
      }
      case 2: // away-side
      switch (iCent) {
        case -1: return 11;
        case 0:  return 9;
        case 1:  return 9;
        case 2:  return 5;
        case 3:  return 4;
        case 4:  return 3;
        case 5:  return 15;
        default: return 0;
      }
      default: return 0;
    }
    case 1: // pTJ > 60 GeV
    switch (iDir) {
      case 0: // near-side
      switch (iCent) {
        case -1: return 2;
        case 0:  return 2;
        case 1:  return 2;
        case 2:  return 2;
        case 3:  return 1;
        case 4:  return 1;
        case 5:  return 2;
        default: return 0;
      }
      case 1: // perpendicular
      switch (iCent) {
        case -1: return 3;
        case 0:  return 1;
        case 1:  return 1;
        case 2:  return 1;
        case 3:  return 1;
        case 4:  return 1;
        case 5:  return 1;
        default: return 0;
      }
      case 2: // away-side
      switch (iCent) {
        case -1: return 3;
        case 0:  return 2;
        case 1:  return 2;
        case 2:  return 2;
        case 3:  return 1;
        case 4:  return 1;
        case 5:  return 3;
        default: return 0;
      }
      default: return 0;
    }
    default: return 0;
  }
}



void ProcessCorrelations (const char* tag, const char* outFileTag, const int nItersMax = 20) {

  TFile* inFile = nullptr;

  const int nItersMin = 1;
  const double* nItersVals = linspace (nItersMin, nItersMax, nItersMax-nItersMin);

  const bool useJetWgts = true;

  TH1D***   h_evt_counts_ref              = Get2DArray <TH1D*> (2, nVar);
  TH1D****  h_jet_counts_ref              = Get3DArray <TH1D*> (2, nPtJBins, nVar);
  //TH1D***   h_evt_counts_ref_bkg          = Get2DArray <TH1D*> (2, nVar);
  TH1D****  h_jet_counts_ref_bkg          = Get3DArray <TH1D*> (2, nPtJBins, nVar);
  TH1D****  h_evt_counts                  = Get3DArray <TH1D*> (2, nZdcCentBins+1, nVar);
  TH1D***** h_jet_counts                  = Get4DArray <TH1D*> (2, nPtJBins, nZdcCentBins+1, nVar);
  //TH1D****  h_evt_counts_bkg              = Get3DArray <TH1D*> (2, nZdcCentBins+1, nVar);
  TH1D***** h_jet_counts_bkg              = Get4DArray <TH1D*> (2, nPtJBins, nZdcCentBins+1, nVar);

  TH1D***   h_jet_pt_ref                  = Get2DArray <TH1D*> (2, nVar);
  TH2D***   h2_jet_pt_cov_ref             = Get2DArray <TH2D*> (2, nVar);
  TH1D****  h_jet_pt                      = Get3DArray <TH1D*> (2, nZdcCentBins+1, nVar);
  TH2D****  h2_jet_pt_cov                 = Get3DArray <TH2D*> (2, nZdcCentBins+1, nVar);

  TH1D***** h_jet_trk_pt_ref              = Get4DArray <TH1D*> (2, nPtJBins, nDir, nVar);
  TH2D***** h2_jet_trk_pt_cov_ref         = Get4DArray <TH2D*> (2, nPtJBins, nDir, nVar);
  TH1D***** h_jet_trk_dphi_ref            = Get4DArray <TH1D*> (2, nPtJBins, nPtChSelections, nVar);
  TH2D***** h2_jet_trk_dphi_cov_ref       = Get4DArray <TH2D*> (2, nPtJBins, nPtChSelections, nVar);

  TH1D***** h_jet_trk_pt_ref_bkg          = Get4DArray <TH1D*> (2, nPtJBins, nDir, nVar);
  TH2D***** h2_jet_trk_pt_cov_ref_bkg     = Get4DArray <TH2D*> (2, nPtJBins, nDir, nVar);
  TH1D***** h_jet_trk_dphi_ref_bkg        = Get4DArray <TH1D*> (2, nPtJBins, nPtChSelections, nVar);
  TH2D***** h2_jet_trk_dphi_cov_ref_bkg   = Get4DArray <TH2D*> (2, nPtJBins, nPtChSelections, nVar);

  TH1D******  h_jet_trk_pt                = Get5DArray <TH1D*> (2, nPtJBins, nDir, nZdcCentBins+1, nVar);
  TH2D******  h2_jet_trk_pt_cov           = Get5DArray <TH2D*> (2, nPtJBins, nDir, nZdcCentBins+1, nVar);
  TH1D******  h_jet_trk_dphi              = Get5DArray <TH1D*> (2, nPtJBins, nPtChSelections, nZdcCentBins+1, nVar);
  TH2D******  h2_jet_trk_dphi_cov         = Get5DArray <TH2D*> (2, nPtJBins, nPtChSelections, nZdcCentBins+1, nVar);

  TH1D******  h_jet_trk_pt_bkg            = Get5DArray <TH1D*> (2, nPtJBins, nDir, nZdcCentBins+1, nVar);
  TH2D******  h2_jet_trk_pt_cov_bkg       = Get5DArray <TH2D*> (2, nPtJBins, nDir, nZdcCentBins+1, nVar);
  TH1D******  h_jet_trk_dphi_bkg          = Get5DArray <TH1D*> (2, nPtJBins, nPtChSelections, nZdcCentBins+1, nVar);
  TH2D******  h2_jet_trk_dphi_cov_bkg     = Get5DArray <TH2D*> (2, nPtJBins, nPtChSelections, nZdcCentBins+1, nVar);

  TH1D*****   h_jet_trk_pt_ref_sig        = Get4DArray <TH1D*> (2, nPtJBins, nDir, nVar);
  TH1D******  h_jet_trk_pt_sig            = Get5DArray <TH1D*> (2, nPtJBins, nDir, nZdcCentBins+1, nVar);
  TH1D*****   h_jet_trk_dphi_ref_sig      = Get4DArray <TH1D*> (2, nPtJBins, nPtChSelections, nVar);
  TH1D******  h_jet_trk_dphi_sig          = Get5DArray <TH1D*> (2, nPtJBins, nPtChSelections, nZdcCentBins+1, nVar);

  TGraph*     g_jet_pt_ref_unfIterUnc     = nullptr;                                // sums of iterations uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_unfIterUnc         = Get1DArray <TGraph*> (nZdcCentBins+1);  // sums of iterations uncertainties as a function of nIter -- data only
  TGraph*     g_jet_pt_ref_unfSumUnc      = nullptr;                                // sums of statistical uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_unfSumUnc          = Get1DArray <TGraph*> (nZdcCentBins+1);  // sums of statistical uncertainties as a function of nIter -- data only
  TGraph*     g_jet_pt_ref_unfTotUnc      = nullptr;                                // sums of total uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_unfTotUnc          = Get1DArray <TGraph*> (nZdcCentBins+1);  // sums of total uncertainties as a function of nIter -- data only

  TGraph*     g_jet_pt_ref_unfIterRelUnc  = nullptr;                                // sums of iterations relative uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_unfIterRelUnc      = Get1DArray <TGraph*> (nZdcCentBins+1);  // sums of iterations relative uncertainties as a function of nIter -- data only
  TGraph*     g_jet_pt_ref_unfSumRelUnc   = nullptr;                                // sums of statistical relative uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_unfSumRelUnc       = Get1DArray <TGraph*> (nZdcCentBins+1);  // sums of statistical relative uncertainties as a function of nIter -- data only
  TGraph*     g_jet_pt_ref_unfTotRelUnc   = nullptr;                                // sums of total relative uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_unfTotRelUnc       = Get1DArray <TGraph*> (nZdcCentBins+1);  // sums of total relative uncertainties as a function of nIter -- data only

  TH1D***     h_jet_pt_ref_unf            = Get2DArray <TH1D*> (2, nVar);
  TH1D****    h_jet_pt_unf                = Get3DArray <TH1D*> (2, nZdcCentBins+1, nVar);

  TH1D****    h_jet_trk_pt_ref_unf_nIters = Get3DArray <TH1D*> (nPtJBins, nDir, nItersMax-nItersMin+2);
  TH1D*****   h_jet_trk_pt_unf_nIters     = Get4DArray <TH1D*> (nPtJBins, nDir, nZdcCentBins+1, nItersMax-nItersMin+2);

  //TH1D*****   h_jet_trk_pt_ref_unf      = Get4DArray <TH1D*> (2, nPtJBins, nDir, nVar);
  //TH1D******  h_jet_trk_pt_unf          = Get5DArray <TH1D*> (2, nPtJBins, nDir, nZdcCentBins+1, nVar);

  //TH1D**    h_jet_trk_pt_ref_sig_bbb    = Get1DArray <TH1D*> (nDir);
  //TH1D***   h_jet_trk_pt_sig_bbb        = Get2DArray <TH1D*> (nDir, nZdcCentBins+1);

  //TF1**    f_jet_trk_pt_ref_sig_bbb     = Get1DArray <TF1*> (nDir);
  //TF1***   f_jet_trk_pt_sig_bbb         = Get2DArray <TF1*> (nDir, nZdcCentBins+1);

  // now the pTJet-integrated histograms (e.g. > 30 GeV and > 60 GeV)
  TH1D*****   h_jetInt_trk_pt_ref                 = Get4DArray <TH1D*> (2, 2, nDir, nVar);
  TH1D******  h_jetInt_trk_pt                     = Get5DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1, nVar);
  TH1D*****   h_jetInt_trk_pt_ref_bkg             = Get4DArray <TH1D*> (2, 2, nDir, nVar);
  TH1D******  h_jetInt_trk_pt_bkg                 = Get5DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1, nVar);
  TH1D*****   h_jetInt_trk_pt_ref_sig             = Get4DArray <TH1D*> (2, 2, nDir, nVar);
  TH1D******  h_jetInt_trk_pt_sig                 = Get5DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1, nVar);

  TGraph***   g_jetInt_trk_pt_ref_unfIterUnc      = Get2DArray <TGraph*> (2, nDir);                           // sums of iterations uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfIterUnc          = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of iterations uncertainties as a function of nIter -- data only
  TGraph***   g_jetInt_trk_pt_ref_unfSumUnc       = Get2DArray <TGraph*> (2, nDir);                           // sums of statistical uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfSumUnc           = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of statistical uncertainties as a function of nIter -- data only
  TGraph***   g_jetInt_trk_pt_ref_unfTotUnc       = Get2DArray <TGraph*> (2, nDir);                           // sums of total uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfTotUnc           = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of total uncertainties as a function of nIter -- data only

  TGraph***   g_jetInt_trk_pt_ref_unfIterRelUnc   = Get2DArray <TGraph*> (2, nDir);                           // sums of iterations relative uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfIterRelUnc       = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of iterations relative uncertainties as a function of nIter -- data only
  TGraph***   g_jetInt_trk_pt_ref_unfSumRelUnc    = Get2DArray <TGraph*> (2, nDir);                           // sums of statistical relative uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfSumRelUnc        = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of statistical relative uncertainties as a function of nIter -- data only
  TGraph***   g_jetInt_trk_pt_ref_unfTotRelUnc    = Get2DArray <TGraph*> (2, nDir);                           // sums of total relative uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfTotRelUnc        = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of total relative uncertainties as a function of nIter -- data only

  TH1D*****   h_jetInt_trk_pt_ref_unf             = Get4DArray <TH1D*> (2, 2, nDir, nVar);
  TH1D******  h_jetInt_trk_pt_unf                 = Get5DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1, nVar);

  TH1D******  h_jetInt_trk_pt_iaa                 = Get5DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1, nVar);
  TH1D******  h_jetInt_trk_pt_iaaNoUnf            = Get5DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1, nVar);

  TH1D*****   h_jetInt_trk_dphi_ref               = Get4DArray <TH1D*> (2, 2, nPtChSelections, nVar);
  TH1D******  h_jetInt_trk_dphi                   = Get5DArray <TH1D*> (2, 2, nPtChSelections, nZdcCentBins+1, nVar);
  TH1D*****   h_jetInt_trk_dphi_ref_bkg           = Get4DArray <TH1D*> (2, 2, nPtChSelections, nVar);
  TH1D******  h_jetInt_trk_dphi_bkg               = Get5DArray <TH1D*> (2, 2, nPtChSelections, nZdcCentBins+1, nVar);
  TH1D*****   h_jetInt_trk_dphi_ref_sig           = Get4DArray <TH1D*> (2, 2, nPtChSelections, nVar);
  TH1D******  h_jetInt_trk_dphi_sig               = Get5DArray <TH1D*> (2, 2, nPtChSelections, nZdcCentBins+1, nVar);


  TH1D******  h_jetInt_trk_dphi_iaa               = Get5DArray <TH1D*> (2, 2, nPtChSelections, nZdcCentBins+1, nVar);

  TH1D***     h_jetInt_trk_pt_ref_unf_closure = Get2DArray <TH1D*> (2, nDir);
  TF1***      f_jetInt_trk_pt_ref_unf_closure = Get2DArray <TF1*> (2, nDir);
  TH1D****    h_jetInt_trk_pt_unf_closure     = Get3DArray <TH1D*> (2, nDir, nZdcCentBins+1);
  TF1****     f_jetInt_trk_pt_unf_closure     = Get3DArray <TF1*> (2, nDir, nZdcCentBins+1);
  TH1D****    h_jetInt_trk_pt_iaa_closure     = Get3DArray <TH1D*> (2, nDir, nZdcCentBins+1);
  TF1****     f_jetInt_trk_pt_iaa_closure     = Get3DArray <TF1*> (2, nDir, nZdcCentBins+1);


  RooUnfoldResponse*    rooUnfResp_jet_pt_ref                     = nullptr;
  RooUnfoldResponse**   rooUnfResp_jet_pt                         = Get1DArray <RooUnfoldResponse*> (nFcalCentBins+1);
  RooUnfoldResponse**   rooUnfResp_jet_trk_pt_ref_sig             = Get1DArray <RooUnfoldResponse*> (nDir);
  RooUnfoldResponse***  rooUnfResp_jet_trk_pt_sig                 = Get2DArray <RooUnfoldResponse*> (nDir, nFcalCentBins+1);


  TGAE****  g_jetInt_trk_pt_ref_syst          = Get3DArray <TGAE*> (2, nDir, nVar);
  TGAE****  g_jetInt_trk_pt_ref_bkg_syst      = Get3DArray <TGAE*> (2, nDir, nVar);
  TGAE****  g_jetInt_trk_pt_ref_sig_syst      = Get3DArray <TGAE*> (2, nDir, nVar);
  TGAE****  g_jetInt_trk_pt_ref_unf_syst      = Get3DArray <TGAE*> (2, nDir, nVar);

  TGAE***** g_jetInt_trk_pt_syst              = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, nVar);
  TGAE***** g_jetInt_trk_pt_bkg_syst          = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, nVar);
  TGAE***** g_jetInt_trk_pt_sig_syst          = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, nVar);
  TGAE***** g_jetInt_trk_pt_unf_syst          = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, nVar);
  TGAE***** g_jetInt_trk_pt_iaa_syst          = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, nVar);

  TGAE****  g_jetInt_trk_dphi_ref_syst        = Get3DArray <TGAE*> (2, nPtChSelections, nVar);
  TGAE****  g_jetInt_trk_dphi_ref_bkg_syst    = Get3DArray <TGAE*> (2, nPtChSelections, nVar);
  TGAE****  g_jetInt_trk_dphi_ref_sig_syst    = Get3DArray <TGAE*> (2, nPtChSelections, nVar);

  TGAE***** g_jetInt_trk_dphi_syst            = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins+1, nVar);
  TGAE***** g_jetInt_trk_dphi_bkg_syst        = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins+1, nVar);
  TGAE***** g_jetInt_trk_dphi_sig_syst        = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins+1, nVar);
  TGAE***** g_jetInt_trk_dphi_iaa_syst        = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins+1, nVar);


  TGAE****  g_jetInt_trk_pt_ref_systTot       = Get3DArray <TGAE*> (2, nDir, 3);
  TGAE****  g_jetInt_trk_pt_ref_bkg_systTot   = Get3DArray <TGAE*> (2, nDir, 3);
  TGAE****  g_jetInt_trk_pt_ref_sig_systTot   = Get3DArray <TGAE*> (2, nDir, 3);
  TGAE****  g_jetInt_trk_pt_ref_unf_systTot   = Get3DArray <TGAE*> (2, nDir, 3);

  TGAE***** g_jetInt_trk_pt_systTot           = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, 3);
  TGAE***** g_jetInt_trk_pt_bkg_systTot       = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, 3);
  TGAE***** g_jetInt_trk_pt_sig_systTot       = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, 3);
  TGAE***** g_jetInt_trk_pt_unf_systTot       = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, 3);
  TGAE***** g_jetInt_trk_pt_iaa_systTot       = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, 3);
 
  TGAE****  g_jetInt_trk_dphi_ref_systTot     = Get3DArray <TGAE*> (2, nPtChSelections, 3);
  TGAE****  g_jetInt_trk_dphi_ref_bkg_systTot = Get3DArray <TGAE*> (2, nPtChSelections, 3);
  TGAE****  g_jetInt_trk_dphi_ref_sig_systTot = Get3DArray <TGAE*> (2, nPtChSelections, 3);

  TGAE***** g_jetInt_trk_dphi_systTot         = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins+1, 3);
  TGAE***** g_jetInt_trk_dphi_bkg_systTot     = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins+1, 3);
  TGAE***** g_jetInt_trk_dphi_sig_systTot     = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins+1, 3);
  TGAE***** g_jetInt_trk_dphi_iaa_systTot     = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins+1, 3);


  TString outFileName = outFileTag;
  outFileName.ReplaceAll (".root", "");
  outFileName = Form ("%s/Results/ProcessCorrelations_%s.root", rootPath.Data (), outFileName.Data ());
  std::cout << "Writing " << outFileName.Data () << std::endl;
  TFile* outFile = new TFile (outFileName.Data (), "recreate");

  for (short iDType = 0; iDType < 2; iDType++) {

    const TString dType = (iDType == 0 ? "data" : "mc");

    for (short iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
        continue;

      {
        TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/%s17_5TeV_hists.root", rootPath.Data (), tag, var.Data (), dType.Data ());
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        h_evt_counts_ref[iDType][iVar]  = (TH1D*) inFile->Get (Form ("h_evt_counts_%s17",   dType.Data ()))->Clone (Form ("h_evt_counts_ref_%s_%s",   dType.Data (), var.Data ()));

        h_jet_pt_ref[iDType][iVar]      = (TH1D*) inFile->Get (Form ("h_jet_pt_%s17",       dType.Data ()))->Clone (Form ("h_jet_pt_ref_%s_%s",       dType.Data (), var.Data ()));
        h2_jet_pt_cov_ref[iDType][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_pt_cov_%s17",  dType.Data ()))->Clone (Form ("h2_jet_pt_cov_ref_%s_%s",  dType.Data (), var.Data ()));

        CalcUncertainties (h_jet_pt_ref[iDType][iVar], h2_jet_pt_cov_ref[iDType][iVar], h_evt_counts_ref[iDType][iVar]);
        h_jet_pt_ref[iDType][iVar]->Scale (h_evt_counts_ref[iDType][iVar]->GetBinContent (2)); // convert distribution to total number of jets by un-scaling 1/N_evt factor

        SaferDelete (&h2_jet_pt_cov_ref[iDType][iVar]);


        for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

          h_jet_counts_ref[iDType][iPtJ][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_%s17", pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_counts_ref_%s_%s_%s", dType.Data (), pTJ.Data (), var.Data ()));

          for (short iDir = 0; iDir < nDir; iDir++) {

            const TString dir = directions[iDir];

            h_jet_trk_pt_ref[iDType][iPtJ][iDir][iVar]      = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_%s_%s17",       dir.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_trk_pt_%s_ref_%s_%s_%s",       dir.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            h2_jet_trk_pt_cov_ref[iDType][iPtJ][iDir][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_cov_%s_%s_%s17",  dir.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_trk_pt_cov_%s_ref_%s_%s_%s",  dir.Data (), dType.Data (), pTJ.Data (), var.Data ()));

            CalcUncertainties (h_jet_trk_pt_ref[iDType][iPtJ][iDir][iVar], h2_jet_trk_pt_cov_ref[iDType][iPtJ][iDir][iVar], h_jet_counts_ref[iDType][iPtJ][iVar]);

            //h_jet_trk_pt_ref[iDType][iPtJ][iDir][iVar]->Rebin (iDir == 1 ? 5 : 2);
            h_jet_trk_pt_ref[iDType][iPtJ][iDir][iVar]->Scale (1., "width");

            SaferDelete (&h2_jet_trk_pt_cov_ref[iDType][iPtJ][iDir][iVar]);

          } // end loop over iDir

          for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

            const TString ptch = pTChSelections[iPtCh].Data ();

            h_jet_trk_dphi_ref[iDType][iPtJ][iPtCh][iVar]       = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_%s_%s17",       ptch.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_trk_dphi_%s_ref_%s_%s_%s",      ptch.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            h2_jet_trk_dphi_cov_ref[iDType][iPtJ][iPtCh][iVar]  = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_cov_%s_%s_%s17",  ptch.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_trk_dphi_cov_%s_ref_%s_%s_%s", ptch.Data (), dType.Data (), pTJ.Data (), var.Data ()));

            CalcUncertainties (h_jet_trk_dphi_ref[iDType][iPtJ][iPtCh][iVar], h2_jet_trk_dphi_cov_ref[iDType][iPtJ][iPtCh][iVar], h_jet_counts_ref[iDType][iPtJ][iVar]);

            h_jet_trk_dphi_ref[iDType][iPtJ][iPtCh][iVar]->Scale (1., "width");

            SaferDelete (&h2_jet_trk_dphi_cov_ref[iDType][iPtJ][iPtCh][iVar]);

          } // end loop over iPtCh

        } // end loop over iPtJ

        inFile->Close ();

      }



      {
        TString inFileName = Form ("%s/Histograms/%s/MixedHists/%s/%s17_5TeV_hists.root", rootPath.Data (), tag, var.Data (), dType.Data ());
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        //h_evt_counts_ref_bkg[iDType][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_mixed_%s17", dType.Data ()))->Clone (Form ("h_evt_counts_ref_bkg_%s_%s", dType.Data (), var.Data ()));


        for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

          h_jet_counts_ref_bkg[iDType][iPtJ][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_mixed_%s17", pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_counts_ref_bkg_%s_%s_%s", dType.Data (), pTJ.Data (), var.Data ()));

          for (short iDir = 0; iDir < nDir; iDir++) {

            const TString dir = directions[iDir];

            h_jet_trk_pt_ref_bkg[iDType][iPtJ][iDir][iVar]      = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_%s_mixed_%s17",       dir.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_trk_pt_%s_ref_bkg_%s_%s_%s",       dir.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            h2_jet_trk_pt_cov_ref_bkg[iDType][iPtJ][iDir][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_cov_%s_%s_mixed_%s17",  dir.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_trk_pt_cov_%s_ref_bkg_%s_%s_%s",  dir.Data (), dType.Data (), pTJ.Data (), var.Data ()));

            CalcUncertainties (h_jet_trk_pt_ref_bkg[iDType][iPtJ][iDir][iVar], h2_jet_trk_pt_cov_ref_bkg[iDType][iPtJ][iDir][iVar], h_jet_counts_ref_bkg[iDType][iPtJ][iVar]);

            //h_jet_trk_pt_ref_bkg[iDType][iPtJ][iDir][iVar]->Rebin (iDir == 1 ? 5 : 2);
            h_jet_trk_pt_ref_bkg[iDType][iPtJ][iDir][iVar]->Scale (1., "width");

            SaferDelete (&h2_jet_trk_pt_cov_ref_bkg[iDType][iPtJ][iDir][iVar]);

          } // end loop over iDir

          for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

            const TString ptch = pTChSelections[iPtCh].Data ();

            h_jet_trk_dphi_ref_bkg[iDType][iPtJ][iPtCh][iVar]       = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_%s_mixed_%s17",       ptch.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_trk_dphi_%s_ref_bkg_%s_%s_%s",      ptch.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            h2_jet_trk_dphi_cov_ref_bkg[iDType][iPtJ][iPtCh][iVar]  = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_cov_%s_%s_mixed_%s17",  ptch.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_trk_dphi_cov_%s_ref_bkg_%s_%s_%s", ptch.Data (), dType.Data (), pTJ.Data (), var.Data ()));

            CalcUncertainties (h_jet_trk_dphi_ref_bkg[iDType][iPtJ][iPtCh][iVar], h2_jet_trk_dphi_cov_ref_bkg[iDType][iPtJ][iPtCh][iVar], h_jet_counts_ref_bkg[iDType][iPtJ][iVar]);

            h_jet_trk_dphi_ref_bkg[iDType][iPtJ][iPtCh][iVar]->Scale (1., "width");

            SaferDelete (&h2_jet_trk_dphi_cov_ref_bkg[iDType][iPtJ][iPtCh][iVar]);

          } // end loop over iPtCh

        } // end loop over iPtJ

        inFile->Close ();

      }



      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

        TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/%s16_5TeV_%s_hists.root", rootPath.Data (), tag, var.Data (), dType.Data (), cent.Data ());
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        h_evt_counts[iDType][iCent][iVar]   = (TH1D*) inFile->Get (Form ("h_evt_counts_%s16",   dType.Data ()))->Clone (Form ("h_evt_counts_pPb_%s_%s_%s",  cent.Data (), dType.Data (), var.Data ()));

        h_jet_pt[iDType][iCent][iVar]       = (TH1D*) inFile->Get (Form ("h_jet_pt_%s16",       dType.Data ()))->Clone (Form ("h_jet_pt_pPb_%s_%s_%s",      cent.Data (), dType.Data (), var.Data ()));
        h2_jet_pt_cov[iDType][iCent][iVar]  = (TH2D*) inFile->Get (Form ("h2_jet_pt_cov_%s16",  dType.Data ()))->Clone (Form ("h2_jet_pt_cov_pPb_%s_%s_%s", cent.Data (), dType.Data (), var.Data ()));

        CalcUncertainties (h_jet_pt[iDType][iCent][iVar], h2_jet_pt_cov[iDType][iCent][iVar], h_evt_counts[iDType][iCent][iVar]);
        h_jet_pt[iDType][iCent][iVar]->Scale (h_evt_counts[iDType][iCent][iVar]->GetBinContent (2)); // convert distribution to total number of jets by un-scaling 1/N_evt factor

        SaferDelete (&h2_jet_pt_cov[iDType][iCent][iVar]);

        for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

          h_jet_counts[iDType][iPtJ][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_%s16",  pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_counts_pPb_%s_%s_%s_%s", cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));

          for (short iDir = 0; iDir < nDir; iDir++) {

            const TString dir = directions[iDir];

            h_jet_trk_pt[iDType][iPtJ][iDir][iCent][iVar]       = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_%s_%s16",       dir.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_trk_pt_%s_pPb_%s_%s_%s_%s",       dir.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            h2_jet_trk_pt_cov[iDType][iPtJ][iDir][iCent][iVar]  = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_cov_%s_%s_%s16",  dir.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_trk_pt_cov_%s_pPb_%s_%s_%s_%s",  dir.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));

            CalcUncertainties (h_jet_trk_pt[iDType][iPtJ][iDir][iCent][iVar], h2_jet_trk_pt_cov[iDType][iPtJ][iDir][iCent][iVar], h_jet_counts[iDType][iPtJ][iCent][iVar]);

            //h_jet_trk_pt[iDType][iPtJ][iDir][iCent][iVar]->Rebin (iDir == 1 ? 5 : 2);
            h_jet_trk_pt[iDType][iPtJ][iDir][iCent][iVar]->Scale (1., "width");

            SaferDelete (&h2_jet_trk_pt_cov[iDType][iPtJ][iDir][iCent][iVar]);

          } // end loop over iDir

          for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

            const TString ptch = pTChSelections[iPtCh].Data ();

            h_jet_trk_dphi[iDType][iPtJ][iPtCh][iCent][iVar]      = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_%s_%s16",       ptch.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_trk_dphi_%s_pPb_%s_%s_%s_%s",      ptch.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            h2_jet_trk_dphi_cov[iDType][iPtJ][iPtCh][iCent][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_cov_%s_%s_%s16",  ptch.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_trk_dphi_cov_%s_pPb_%s_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));

            CalcUncertainties (h_jet_trk_dphi[iDType][iPtJ][iPtCh][iCent][iVar], h2_jet_trk_dphi_cov[iDType][iPtJ][iPtCh][iCent][iVar], h_jet_counts[iDType][iPtJ][iCent][iVar]);

            h_jet_trk_dphi[iDType][iPtJ][iPtCh][iCent][iVar]->Scale (1., "width");

            SaferDelete (&h2_jet_trk_dphi_cov[iDType][iPtJ][iPtCh][iCent][iVar]);

          } // end loop over iPtCh

        } // end loop over iPtJ

        inFile->Close ();

      } // end loop over iCent



      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

        TString inFileName = Form ("%s/Histograms/%s/MixedHists/%s/%s16_5TeV_%s_hists.root", rootPath.Data (), tag, var.Data (), dType.Data (), cent.Data ());
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        //h_evt_counts_bkg[iDType][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_mixed_%s16", dType.Data ()))->Clone (Form ("h_evt_counts_pPb_bkg_%s_%s_%s", cent.Data (), dType.Data (), var.Data ()));


        for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

          h_jet_counts_bkg[iDType][iPtJ][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_mixed_%s16",  pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_counts_pPb_bkg_%s_%s_%s_%s", cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));


          for (short iDir = 0; iDir < nDir; iDir++) {

            const TString dir = directions[iDir];

            h_jet_trk_pt_bkg[iDType][iPtJ][iDir][iCent][iVar]       = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_%s_mixed_%s16",       dir.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_trk_pt_%s_pPb_bkg_%s_%s_%s_%s",       dir.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            h2_jet_trk_pt_cov_bkg[iDType][iPtJ][iDir][iCent][iVar]  = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_cov_%s_%s_mixed_%s16",  dir.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_trk_pt_cov_%s_pPb_bkg_%s_%s_%s_%s",  dir.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));

            CalcUncertainties (h_jet_trk_pt_bkg[iDType][iPtJ][iDir][iCent][iVar], h2_jet_trk_pt_cov_bkg[iDType][iPtJ][iDir][iCent][iVar], h_jet_counts_bkg[iDType][iPtJ][iCent][iVar]);

            //h_jet_trk_pt_bkg[iDType][iPtJ][iDir][iCent][iVar]->Rebin (iDir == 1 ? 5 : 2);
            h_jet_trk_pt_bkg[iDType][iPtJ][iDir][iCent][iVar]->Scale (1., "width");

            SaferDelete (&h2_jet_trk_pt_cov_bkg[iDType][iPtJ][iDir][iCent][iVar]);

          } // end loop over iDir

          for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

            const TString ptch = pTChSelections[iPtCh].Data ();

            h_jet_trk_dphi_bkg[iDType][iPtJ][iPtCh][iCent][iVar]      = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_%s_mixed_%s16",       ptch.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_trk_dphi_%s_pPb_bkg_%s_%s_%s_%s",      ptch.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            h2_jet_trk_dphi_cov_bkg[iDType][iPtJ][iPtCh][iCent][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_cov_%s_%s_mixed_%s16",  ptch.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_trk_dphi_cov_%s_pPb_bkg_%s_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));

            CalcUncertainties (h_jet_trk_dphi_bkg[iDType][iPtJ][iPtCh][iCent][iVar], h2_jet_trk_dphi_cov_bkg[iDType][iPtJ][iPtCh][iCent][iVar], h_jet_counts_bkg[iDType][iPtJ][iCent][iVar]);

            h_jet_trk_dphi_bkg[iDType][iPtJ][iPtCh][iCent][iVar]->Scale (1., "width");

            SaferDelete (&h2_jet_trk_dphi_cov_bkg[iDType][iPtJ][iPtCh][iCent][iVar]);

          } // end loop over iPtCh

        } // end loop over iPtJ

        inFile->Close ();

      } // end loop over iCent

    } // end loop over iVar

  } // end loop over iDType




  {
    TFile* inFile = new TFile (Form ("%s/MakeResponseMatrix/Nominal/allSamples.root", rootPath.Data ()), "read");

    if (useJetWgts) rooUnfResp_jet_pt_ref = (RooUnfoldResponse*) inFile->Get ("rooUnfResp_jet_pt_ref_mc_wgts");
    else            rooUnfResp_jet_pt_ref = (RooUnfoldResponse*) inFile->Get ("rooUnfResp_jet_pt_ref_mc_fullClosure");

    for (short iDir = 0; iDir < nDir; iDir++) {

      const TString dir = directions[iDir];

      if (useJetWgts) rooUnfResp_jet_trk_pt_ref_sig[iDir] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_trk_pt_%s_ref_sig_mc_wgts", dir.Data ()));
      else            rooUnfResp_jet_trk_pt_ref_sig[iDir] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_trk_pt_%s_ref_sig_mc_fullClosure", dir.Data ()));

    } // end loop over iDir

    for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

      const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

      if (useJetWgts) rooUnfResp_jet_pt[iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_pt_%s_mc_wgts", cent));
      else            rooUnfResp_jet_pt[iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_pt_%s_mc_fullClosure", cent));

      for (short iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        if (useJetWgts) rooUnfResp_jet_trk_pt_sig[iDir][iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_trk_pt_%s_%s_sig_mc_wgts", dir.Data (), cent));
        else            rooUnfResp_jet_trk_pt_sig[iDir][iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_trk_pt_%s_%s_sig_mc_fullClosure", dir.Data (), cent));

      } // end loop over iFile

    } // end loop over iDir
  }




  {
    TFile* inFile = new TFile (Form ("%s/aux/MCClosureHists.root", workPath.Data ()), "read");

    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      for (short iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        h_jetInt_trk_pt_ref_unf_closure[iPtJInt][iDir] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_sig_unf_%s_halfClosure", dir.Data (), pTJInt.Data ()));
        f_jetInt_trk_pt_ref_unf_closure[iPtJInt][iDir] = (TF1*)  inFile->Get (Form ("f_jetInt_trk_pt_%s_ref_sig_unf_%s_halfClosure", dir.Data (), pTJInt.Data ()));

        for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

          const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

          h_jetInt_trk_pt_unf_closure[iPtJInt][iDir][iCent] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_%s_sig_unf_%s_halfClosure", dir.Data (), cent, pTJInt.Data ()));
          f_jetInt_trk_pt_unf_closure[iPtJInt][iDir][iCent] = (TF1*)  inFile->Get (Form ("f_jetInt_trk_pt_%s_%s_sig_unf_%s_halfClosure", dir.Data (), cent, pTJInt.Data ()));

          h_jetInt_trk_pt_iaa_closure[iPtJInt][iDir][iCent] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_%s_iaa_unf_%s_halfClosure", dir.Data (), cent, pTJInt.Data ()));
          f_jetInt_trk_pt_iaa_closure[iPtJInt][iDir][iCent] = (TF1*)  inFile->Get (Form ("f_jetInt_trk_pt_%s_%s_iaa_unf_%s_halfClosure", dir.Data (), cent, pTJInt.Data ()));

        } // end loop over iCent

      } // end loop over iDir

    } // end loop over iPtJInt
  }




  ////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //// DERIVE BIN-BY-BIN UNFOLDING FACTORS
  ////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //{
  //  const short iVarNum = GetVarN ("MCTruthJetsTruthParts");
  //  const short iVarDen = GetVarN ("MCRecoJetTruthParts");

  //  TF1* f = nullptr;

  //  const TString funcStr = "[0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)+[4]*pow(log(x),4)";

  //  //bbbf.SetNDeriv (2);
  //  //bbbf.SetDegree (6);

  //  //const int ndf = bbbf.NDF ();

  //  for (short iDir = 0; iDir < nDir; iDir++) {

  //    const TString dir = directions[iDir];

  //    TH1D* h = nullptr;

  //    h = (TH1D*) h_jet_trk_pt_ref[1][iDir][iVarNum]->Clone (Form ("h_jet_trk_pt_%s_ref_sig_bbb", dir.Data ()));
  //    h->Divide (h_jet_trk_pt_ref[1][iDir][iVarDen]);
  //    h_jet_trk_pt_ref_sig_bbb[iDir] = h;

  //    f = new TF1 (Form ("f_jet_trk_pt_%s_ref_sig_bbb", dir.Data ()), funcStr, h->GetXaxis ()->GetBinLowEdge (1), h->GetXaxis ()->GetBinLowEdge (h->GetNbinsX () + 1));
  //    //f = new TF1 (Form ("f_jet_trk_pt_%s_ref_sig_bbb", dir.Data ()), &bbbf, h->GetXaxis ()->GetBinLowEdge (1), h->GetXaxis ()->GetBinLowEdge (h->GetNbinsX () + 1), ndf);
  //    if (iDir == 1) {
  //      double x1, x2;
  //      f->GetRange (x1, x2);
  //      f->SetRange (x1, 20);
  //    }
  //    //f->SetParameter (0, 100);
  //    //f->SetParameter (1, 0);
  //    //f->SetParameter (2, 0);

  //    f_jet_trk_pt_ref_sig_bbb[iDir] = f;
  //    h->Fit (f, "RN0Q");

  //    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

  //      const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

  //      h = (TH1D*) h_jet_trk_pt[1][iDir][iCent][iVarNum]->Clone (Form ("h_jet_trk_pt_%s_pPb_sig_bbb_%s", dir.Data (), cent.Data ()));
  //      h->Divide (h_jet_trk_pt[1][iDir][iCent][iVarDen]);
  //      h_jet_trk_pt_sig_bbb[iDir][iCent] = h;

  //      f = new TF1 (Form ("f_jet_trk_pt_%s_pPb_sig_bbb_%s", dir.Data (), cent.Data ()), funcStr, h->GetXaxis ()->GetBinLowEdge (1), h->GetXaxis ()->GetBinLowEdge (h->GetNbinsX () + 1));
  //      //f = new TF1 (Form ("f_jet_trk_pt_%s_pPb_sig_bbb_%s", dir.Data (), cent.Data ()), &bbbf, h->GetXaxis ()->GetBinLowEdge (1), h->GetXaxis ()->GetBinLowEdge (h->GetNbinsX () + 1), ndf);
  //      if (iDir == 1) {
  //        double x1, x2;
  //        f->GetRange (x1, x2);
  //        f->SetRange (x1, 20);
  //      }
  //      //f->SetParameter (0, 100);
  //      //f->SetParameter (1, 0);
  //      //f->SetParameter (2, 0);

  //      f_jet_trk_pt_sig_bbb[iDir][iCent] = f;
  //      h->Fit (f, "RN0Q");

  //    } // end loop over iCent

  //  } // end loop over iDir
  //}




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CALCULATE SIGNAL YIELDS BY SUBTRACTING THE BACKGROUND COMPONENT
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iDType = 0; iDType < 2; iDType++) {

    const TString dType = (iDType == 0 ? "data" : "mc");

    for (short iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
        continue;

      const bool hasRefBkg = (iDType != 1 && variationsWithNoppBkgd.count (var) == 0);
      //const bool hasRefBkg = (variationsWithNoppBkgd.count (var) == 0);
      const bool hasBkg = (variationsWithNopPbBkgd.count (var) == 0);

      for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          h_jet_trk_pt_ref_sig[iDType][iPtJ][iDir][iVar] = (TH1D*) h_jet_trk_pt_ref[iDType][iPtJ][iDir][iVar]->Clone (Form ("h_jet_trk_pt_%s_ref_sig_%s_%s_%s", dir.Data (), dType.Data (), pTJ.Data (), var.Data ()));
          if (hasRefBkg)
            h_jet_trk_pt_ref_sig[iDType][iPtJ][iDir][iVar]->Add (h_jet_trk_pt_ref_bkg[iDType][iPtJ][iDir][iVar], -1);

          //if (doUnfold)
          //   MultiplyByTF1 (h_jet_trk_pt_ref_sig[iDType][iDir][iVar], f_jet_trk_pt_ref_sig_bbb[iDir]);

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            h_jet_trk_pt_sig[iDType][iPtJ][iDir][iCent][iVar] = (TH1D*) h_jet_trk_pt[iDType][iPtJ][iDir][iCent][iVar]->Clone (Form ("h_jet_trk_pt_%s_pPb_sig_%s_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            if (hasBkg)
              h_jet_trk_pt_sig[iDType][iPtJ][iDir][iCent][iVar]->Add (h_jet_trk_pt_bkg[iDType][iPtJ][iDir][iCent][iVar], -1);

            //if (doUnfold)
            //  MultiplyByTF1 (h_jet_trk_pt_sig[iDType][iDir][iCent][iVar], f_jet_trk_pt_sig_bbb[iDir][iCent]);

            //h_jet_trk_pt_iaa[iDType][iPtJ][iDir][iCent][iVar] = (TH1D*) h_jet_trk_pt_sig[iDType][iPtJ][iDir][iCent][iVar]->Clone (Form ("h_jet_trk_pt_%s_iaa_%s_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            //h_jet_trk_pt_iaa[iDType][iPtJ][iDir][iCent][iVar]->Divide (h_jet_trk_pt_ref_sig[iDType][iPtJ][iDir][iVar]);

          } // end loop over iCent

        } // end loop over iDir


        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          const TString ptch = pTChSelections[iPtCh].Data ();

          h_jet_trk_dphi_ref_sig[iDType][iPtJ][iPtCh][iVar] = (TH1D*) h_jet_trk_dphi_ref[iDType][iPtJ][iPtCh][iVar]->Clone (Form ("h_jet_trk_dphi_%s_ref_sig_%s_%s_%s", ptch.Data (), dType.Data (), pTJ.Data (), var.Data ()));
          if (hasRefBkg)
            h_jet_trk_dphi_ref_sig[iDType][iPtJ][iPtCh][iVar]->Add (h_jet_trk_dphi_ref_bkg[iDType][iPtJ][iPtCh][iVar], -1);

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            h_jet_trk_dphi_sig[iDType][iPtJ][iPtCh][iCent][iVar] = (TH1D*) h_jet_trk_dphi[iDType][iPtJ][iPtCh][iCent][iVar]->Clone (Form ("h_jet_trk_dphi_%s_pPb_sig_%s_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            if (hasBkg)
              h_jet_trk_dphi_sig[iDType][iPtJ][iPtCh][iCent][iVar]->Add (h_jet_trk_dphi_bkg[iDType][iPtJ][iPtCh][iCent][iVar], -1);

            //h_jet_trk_dphi_iaa[iDType][iPtJ][iPtCh][iCent][iVar] = (TH1D*) h_jet_trk_dphi_sig[iDType][iPtJ][iPtCh][iCent][iVar]->Clone (Form ("h_jet_trk_dphi_%s_iaa_%s_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            //h_jet_trk_dphi_iaa[iDType][iPtJ][iPtCh][iCent][iVar]->Divide (h_jet_trk_dphi_ref_sig[iDType][iPtJ][iPtCh][iVar]);

          } // end loop over iCent

        } // end loop over iPtCh

      } // end loop over iPtJ

    } // end loop over iVar

  } // end loop over iDType




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CREATE 2D RESULT HISTOGRAMS, DO 2D BAYES UNFOLD
  // THEN CONVERT BACK TO 1D HISTOGRAMS FOR INTEGRATED JET PT SELECTIONS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iDType = 0; iDType < 2; iDType++) {

    const TString dType = (iDType == 0 ? "data" : "mc");

    for (short iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
        continue;

      const bool doUnfold = (variationsWithNoUnfold.count (var) == 0);

      if (doUnfold) {

        if (iDType == 0 && iVar == 0) {

          g_jet_pt_ref_unfSumUnc      = new TGraph (nItersMax - nItersMin + 1);
          g_jet_pt_ref_unfIterUnc     = new TGraph (nItersMax - nItersMin + 1);
          g_jet_pt_ref_unfTotUnc      = new TGraph (nItersMax - nItersMin + 1);
          g_jet_pt_ref_unfSumRelUnc   = new TGraph (nItersMax - nItersMin + 1);
          g_jet_pt_ref_unfIterRelUnc  = new TGraph (nItersMax - nItersMin + 1);
          g_jet_pt_ref_unfTotRelUnc   = new TGraph (nItersMax - nItersMin + 1);

          TH1D* h_unf = nullptr;
          TH1D* h_unf_prev = (TH1D*) h_jet_pt_ref[iDType][iVar]->Clone ("h_unf_0Iters");

          for (short iIter = 0; iIter < nItersMax-nItersMin+2; iIter++) {

            const short nIters = (short) nItersVals[iIter];

            RooUnfoldBayes* bayesUnf = new RooUnfoldBayes (rooUnfResp_jet_pt_ref, h_jet_pt_ref[iDType][iVar], nIters);
            bayesUnf->SetVerbose (-1);
            h_unf = (TH1D*) bayesUnf->Hreco ()->Clone (Form ("h_unf_%iIters", nIters));
            SaferDelete (&bayesUnf);

            double totErr = 0, totRelErr = 0;
            for (int iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
              if (h_unf->GetXaxis ()->GetBinCenter (iX) < 30 || 300 < h_unf->GetXaxis ()->GetBinCenter (iX)) continue;
              totErr += std::pow (h_unf->GetBinError (iX), 2);
              totRelErr += std::pow (h_unf->GetBinError (iX), 2) / h_unf->GetBinContent (iX);
            }

            double iterErr = 0, iterRelErr = 0;
            for (int iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
              if (h_unf->GetXaxis ()->GetBinCenter (iX) < 30 || 300 < h_unf->GetXaxis ()->GetBinCenter (iX)) continue;
              iterErr += std::fabs (h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX));
              iterRelErr += std::fabs (h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX)) / std::sqrt (h_unf->GetBinContent (iX));
            }

            g_jet_pt_ref_unfSumUnc->SetPoint      (nIters - nItersMin, nIters, std::sqrt (totErr));
            g_jet_pt_ref_unfIterUnc->SetPoint     (nIters - nItersMin, nIters, std::sqrt (iterErr));
            g_jet_pt_ref_unfTotUnc->SetPoint      (nIters - nItersMin, nIters, std::sqrt (totErr + iterErr));
            g_jet_pt_ref_unfSumRelUnc->SetPoint   (nIters - nItersMin, nIters, std::sqrt (totRelErr));
            g_jet_pt_ref_unfIterRelUnc->SetPoint  (nIters - nItersMin, nIters, std::sqrt (iterRelErr));
            g_jet_pt_ref_unfTotRelUnc->SetPoint   (nIters - nItersMin, nIters, std::sqrt (totRelErr + iterRelErr));

            SaferDelete (&h_unf_prev);
            h_unf_prev = h_unf;

          } // end loop over iIter

          SaferDelete (&h_unf);


          {
            TH1D* h = h_jet_pt_ref[iDType][iVar];
            double totErr = 0, totRelErr = 0;
            for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
              if (h->GetXaxis ()->GetBinCenter (iX) < 30 || 300 < h->GetXaxis ()->GetBinCenter (iX)) continue;
              totErr += std::pow (h->GetBinError (iX), 2);
              totRelErr += std::pow (h->GetBinError (iX), 2) / h->GetBinContent (iX);
            }

            g_jet_pt_ref_unfSumUnc->SetPoint    (g_jet_pt_ref_unfSumUnc->GetN (),     0, std::sqrt (totErr));
            g_jet_pt_ref_unfTotUnc->SetPoint    (g_jet_pt_ref_unfTotUnc->GetN (),     0, std::sqrt (totErr));
            g_jet_pt_ref_unfSumRelUnc->SetPoint (g_jet_pt_ref_unfSumRelUnc->GetN (),  0, std::sqrt (totRelErr));
            g_jet_pt_ref_unfTotRelUnc->SetPoint (g_jet_pt_ref_unfTotRelUnc->GetN (),  0, std::sqrt (totRelErr));
          }

        }
  
        RooUnfoldBayes* bayesUnf = new RooUnfoldBayes (rooUnfResp_jet_pt_ref, h_jet_pt_ref[iDType][iVar], GetJetSpectraNIters (-1));
        bayesUnf->SetVerbose (-1);
        h_jet_pt_ref_unf[iDType][iVar] = (TH1D*) bayesUnf->Hreco ()->Clone (Form ("h_jet_pt_ref_unf_%s_%s", dType.Data (), var.Data ()));
        SaferDelete (&bayesUnf);

      }
      else {

        h_jet_pt_ref_unf[iDType][iVar] = (TH1D*) h_jet_pt_ref[iDType][iVar];

      }

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

        if (doUnfold) {

          if (iDType == 0 && iVar == 0) {

            g_jet_pt_unfSumUnc[iCent]     = new TGraph (nItersMax - nItersMin + 1);
            g_jet_pt_unfIterUnc[iCent]    = new TGraph (nItersMax - nItersMin + 1);
            g_jet_pt_unfTotUnc[iCent]     = new TGraph (nItersMax - nItersMin + 1);
            g_jet_pt_unfSumRelUnc[iCent]  = new TGraph (nItersMax - nItersMin + 1);
            g_jet_pt_unfIterRelUnc[iCent] = new TGraph (nItersMax - nItersMin + 1);
            g_jet_pt_unfTotRelUnc[iCent]  = new TGraph (nItersMax - nItersMin + 1);
 
            TH1D* h_unf = nullptr;
            TH1D* h_unf_prev = (TH1D*) h_jet_pt[iDType][iCent][iVar]->Clone ("h_unf_0Iters");

            for (short iIter = 0; iIter < nItersMax-nItersMin+2; iIter++) {
  
              const short nIters = (short) nItersVals[iIter];

              RooUnfoldBayes* bayesUnf = new RooUnfoldBayes (rooUnfResp_jet_pt[iCent], h_jet_pt[iDType][iCent][iVar], nIters);
              bayesUnf->SetVerbose (-1);
              h_unf = (TH1D*) bayesUnf->Hreco ()->Clone (Form ("h_unf_%iIters", nIters));
              SaferDelete (&bayesUnf);

              double totErr = 0, totRelErr = 0;
              for (int iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
                if (h_unf->GetXaxis ()->GetBinCenter (iX) < 30 || 300 < h_unf->GetXaxis ()->GetBinCenter (iX)) continue;
                totErr += std::pow (h_unf->GetBinError (iX), 2);
                totRelErr += std::pow (h_unf->GetBinError (iX), 2) / h_unf->GetBinContent (iX);
              }
  
              double iterErr = 0, iterRelErr = 0;
              for (int iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
                if (h_unf->GetXaxis ()->GetBinCenter (iX) < 30 || 300 < h_unf->GetXaxis ()->GetBinCenter (iX)) continue;
                iterErr += std::fabs (h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX));
                iterRelErr += std::fabs (h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX)) / std::sqrt (h_unf->GetBinContent (iX));
              }
  
              g_jet_pt_unfSumUnc[iCent]->SetPoint     (nIters - nItersMin, nIters, std::sqrt (totErr));
              g_jet_pt_unfIterUnc[iCent]->SetPoint    (nIters - nItersMin, nIters, std::sqrt (iterErr));
              g_jet_pt_unfTotUnc[iCent]->SetPoint     (nIters - nItersMin, nIters, std::sqrt (totErr + iterErr));
              g_jet_pt_unfSumRelUnc[iCent]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (totRelErr));
              g_jet_pt_unfIterRelUnc[iCent]->SetPoint (nIters - nItersMin, nIters, std::sqrt (iterRelErr));
              g_jet_pt_unfTotRelUnc[iCent]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (totRelErr + iterRelErr));

              SaferDelete (&h_unf_prev);
              h_unf_prev = h_unf;

            } // end loop over iIter

            SaferDelete (&h_unf);


            {
              TH1D* h = h_jet_pt[iDType][iCent][iVar];
              double totErr = 0, totRelErr = 0;
              for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
                if (h->GetXaxis ()->GetBinCenter (iX) < 30 || 300 < h->GetXaxis ()->GetBinCenter (iX)) continue;
                totErr += std::pow (h->GetBinError (iX), 2);
                totRelErr += std::pow (h->GetBinError (iX), 2) / h->GetBinContent (iX);
              }

              g_jet_pt_unfSumUnc[iCent]->SetPoint     (g_jet_pt_unfSumUnc[iCent]->GetN (),    0, std::sqrt (totErr));
              g_jet_pt_unfTotUnc[iCent]->SetPoint     (g_jet_pt_unfTotUnc[iCent]->GetN (),    0, std::sqrt (totErr));
              g_jet_pt_unfSumRelUnc[iCent]->SetPoint  (g_jet_pt_unfSumRelUnc[iCent]->GetN (), 0, std::sqrt (totRelErr));
              g_jet_pt_unfTotRelUnc[iCent]->SetPoint  (g_jet_pt_unfTotRelUnc[iCent]->GetN (), 0, std::sqrt (totRelErr));
            }
  
          }

          RooUnfoldBayes* bayesUnf = new RooUnfoldBayes (rooUnfResp_jet_pt[iCent], h_jet_pt[iDType][iCent][iVar], GetJetSpectraNIters (iCent));
          bayesUnf->SetVerbose (-1);
          h_jet_pt_unf[iDType][iCent][iVar] = (TH1D*) bayesUnf->Hreco ()->Clone (Form ("h_jet_pt_unf_%s_%s_%s", dType.Data (), cent.Data (), var.Data ()));
          SaferDelete (&bayesUnf);

        }
        else {
  
          h_jet_pt_unf[iDType][iCent][iVar] = (TH1D*) h_jet_pt[iDType][iCent][iVar];
  
        }

      } // end loop over iCent


      if (doUnfold) {
        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          TH2D* h2 = new TH2D ("h2_temp", "", nPtChBins, pTChBins, nPtJBins, pTJBins);
          h2->Sumw2 ();

          for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

            TH1D* h = h_jet_trk_pt_ref_sig[iDType][iPtJ][iDir][iVar];
            const float nJet = h_jet_pt_ref[iDType][iVar]->GetBinContent (iPtJ+1);

            for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {

              h2->SetBinContent (iPtCh+1, iPtJ+1, nJet * h->GetBinContent (iPtCh+1) * h->GetBinWidth (iPtCh+1));
              h2->SetBinError   (iPtCh+1, iPtJ+1, nJet * h->GetBinError   (iPtCh+1) * h->GetBinWidth (iPtCh+1));

            } // end loop over iPtCh

          } // end loop over iPtJ


          if (iDType == 0 && iVar == 0) {

            for (short iPtJInt : {0, 1}) {

              g_jetInt_trk_pt_ref_unfSumUnc[iPtJInt][iDir]      = new TGraph (nItersMax - nItersMin + 1);
              g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir]     = new TGraph (nItersMax - nItersMin + 1);
              g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir]      = new TGraph (nItersMax - nItersMin + 1);
              g_jetInt_trk_pt_ref_unfSumRelUnc[iPtJInt][iDir]   = new TGraph (nItersMax - nItersMin + 1);
              g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir]  = new TGraph (nItersMax - nItersMin + 1);
              g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir]   = new TGraph (nItersMax - nItersMin + 1);

            } // end loop over iPtJInt

            TH2D* h2_unf = nullptr;
            TH2D* h2_unf_prev = (TH2D*) h2->Clone ("h2_unf_0Iters");
  
            for (short iIter = 0; iIter < nItersMax-nItersMin+2; iIter++) {
  
              const short nIters = (short) nItersVals[iIter];

              RooUnfoldBayes* bayesUnf = new RooUnfoldBayes (rooUnfResp_jet_trk_pt_ref_sig[iDir], h2, nIters);
              bayesUnf->SetVerbose (-1);
              h2_unf = (TH2D*) bayesUnf->Hreco ()->Clone (Form ("h2_unf_%iIters", nIters));
              SaferDelete (&bayesUnf);

              for (short iPtJInt : {0, 1}) {

                const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
                const float minJetPt = (iPtJInt == 0 ? 30. : 60.);

                TH1D* h = new TH1D (Form ("h_jet_trk_pt_%s_ref_unf_data_%s_Nominal_nIters%i", dir.Data (), pTJInt.Data (), nIters), "", nPtChBins, pTChBins);
                h->Sumw2 ();
                h_jet_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter] = h;
  
                double totalJetsUF = 0;
                for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
    
                  if (0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]) < minJetPt || 300 < 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1])) continue;
    
                  totalJetsUF += h_jet_pt_ref_unf[iDType][iVar]->GetBinContent (iPtJ+1);
    
                  for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
    
                    h->SetBinContent (iPtCh+1, h->GetBinContent (iPtCh+1) + h2_unf->GetBinContent (iPtCh+1, iPtJ+1));
                    h->SetBinError (iPtCh+1, h->GetBinError (iPtCh+1) + std::pow (h2_unf->GetBinError (iPtCh+1, iPtJ+1), 2));
    
                  } // end loop over iPtCh
    
                } // end loop over iPtJ

                for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
  
                  h->SetBinContent (iPtCh+1, h->GetBinContent (iPtCh+1) / (totalJetsUF * h->GetBinWidth (iPtCh+1)));
                  h->SetBinError (iPtCh+1, std::sqrt (h->GetBinError (iPtCh+1)) / (totalJetsUF * h->GetBinWidth (iPtCh+1)));
  
                } // end loop over iPtCh

              } // end loop over iPtJInt

              double totErr[2] = {0, 0};
              double totRelErr[2] = {0, 0};
              for (int iX = 1; iX <= h2_unf->GetNbinsX (); iX++) {
                const double ptch = h2_unf->GetXaxis ()->GetBinCenter (iX);
                for (int iY = 1; iY <= h2_unf->GetNbinsY (); iY++) {
                  for (short iPtJInt : {0, 1}) {
                    const double jpt = h2_unf->GetYaxis ()->GetBinCenter (iY);
                    if (ptch < 5 || ptch > (iPtJInt == 0 ? 40 : 75)) continue;
                    totErr[iPtJInt] += ((iPtJInt == 0 ? 30. : 60.) < jpt && jpt < 300 ? std::pow (h2_unf->GetBinError (iX, iY), 2) : 0);
                    totRelErr[iPtJInt] += ((iPtJInt == 0 ? 30. : 60.) < jpt && jpt < 300 && h2_unf->GetBinContent (iX, iY) != 0 ? std::pow (h2_unf->GetBinError (iX, iY), 2) / h2_unf->GetBinContent (iX, iY) : 0);
                  }
                } // end loop over iY
              } // end loop over iX
  
              double iterErr[2] = {0, 0};
              double iterRelErr[2] = {0, 0};
              for (int iX = 1; iX <= h2_unf->GetNbinsX (); iX++) {
                const double ptch = h2_unf->GetXaxis ()->GetBinCenter (iX);
                for (int iY = 1; iY <= h2_unf->GetNbinsY (); iY++) {
                  for (short iPtJInt : {0, 1}) {
                    const double jpt = h2_unf->GetYaxis ()->GetBinCenter (iY);
                    if (ptch < 5 || ptch > (iPtJInt == 0 ? 40 : 75)) continue;
                    iterErr[iPtJInt] += ((iPtJInt == 0 ? 30. : 60.) < jpt && jpt < 300 ? std::fabs (h2_unf->GetBinContent (iX, iY) - h2_unf_prev->GetBinContent (iX, iY)) : 0);
                    iterRelErr[iPtJInt] += ((iPtJInt == 0 ? 30. : 60.) < jpt && jpt < 300 && h2_unf->GetBinContent (iX, iY) != 0 ? std::fabs (h2_unf->GetBinContent (iX, iY) - h2_unf_prev->GetBinContent (iX, iY)) / std::sqrt (h2_unf->GetBinContent (iX, iY)) : 0);
                  }
                } // end loop over iY
              } // end loop over iX
  
 
              for (short iPtJInt : {0, 1}) {

                g_jetInt_trk_pt_ref_unfSumUnc[iPtJInt][iDir]->SetPoint      (nIters - nItersMin, nIters, std::sqrt (totErr[iPtJInt]));
                g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir]->SetPoint     (nIters - nItersMin, nIters, std::sqrt (std::fabs (iterErr[iPtJInt])));
                g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir]->SetPoint      (nIters - nItersMin, nIters, std::sqrt (totErr[iPtJInt] + std::fabs (iterErr[iPtJInt])));
                g_jetInt_trk_pt_ref_unfSumRelUnc[iPtJInt][iDir]->SetPoint   (nIters - nItersMin, nIters, std::sqrt (totRelErr[iPtJInt]));
                g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (std::fabs (iterRelErr[iPtJInt])));
                g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir]->SetPoint   (nIters - nItersMin, nIters, std::sqrt (totRelErr[iPtJInt] + std::fabs (iterRelErr[iPtJInt])));

              } // end loop over iPtJInt

              SaferDelete (&h2_unf_prev);
              h2_unf_prev = h2_unf;

            } // end loop over iIter

            SaferDelete (&h2_unf);


            double totErr[2] = {0, 0};
            double totRelErr[2] = {0, 0};
            for (short iPtJInt : {0, 1}) {

              for (int iX = 1; iX <= h2->GetNbinsX (); iX++) {
                const double ptch = h2->GetXaxis ()->GetBinCenter (iX);
                for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {
                  for (short iPtJInt : {0, 1}) {
                    const double jpt = h2->GetYaxis ()->GetBinCenter (iY);
                    if (ptch < 5 || ptch > (iPtJInt == 0 ? 40 : 75)) continue;
                    totErr[iPtJInt] += ((iPtJInt == 0 ? 30. : 60.) < jpt && jpt < 300 ? std::pow (h2->GetBinError (iX, iY), 2) : 0);
                    totRelErr[iPtJInt] += ((iPtJInt == 0 ? 30. : 60.) < jpt && jpt < 300 && h2->GetBinContent (iX, iY) != 0 ? std::pow (h2->GetBinError (iX, iY), 2) / h2->GetBinContent (iX, iY) : 0);
                  }
                } // end loop over iY
              } // end loop over iX

              g_jetInt_trk_pt_ref_unfSumUnc[iPtJInt][iDir]->SetPoint    (g_jetInt_trk_pt_ref_unfSumUnc[iPtJInt][iDir]->GetN (),     0, std::sqrt (totErr[iPtJInt]));
              g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir]->SetPoint    (g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir]->GetN (),     0, std::sqrt (totErr[iPtJInt]));
              g_jetInt_trk_pt_ref_unfSumRelUnc[iPtJInt][iDir]->SetPoint (g_jetInt_trk_pt_ref_unfSumRelUnc[iPtJInt][iDir]->GetN (),  0, std::sqrt (totRelErr[iPtJInt]));
              g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir]->SetPoint (g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir]->GetN (),  0, std::sqrt (totRelErr[iPtJInt]));

            } // end loop over iPtJInt
  
          }


          for (short iPtJInt : {0, 1}) {

            const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
            const float minJetPt = (iPtJInt == 0 ? 30. : 60.);

            RooUnfoldBayes* bayesUnf = new RooUnfoldBayes (rooUnfResp_jet_trk_pt_ref_sig[iDir], h2, GetTrkSpectraNIters (iPtJInt, iDir, -1));
            bayesUnf->SetVerbose (-1);
            TH2D* h2_unf = (TH2D*) bayesUnf->Hreco ()->Clone ("h2_unf");
            SaferDelete (&bayesUnf);

            TH1D* h = new TH1D (Form ("h_jetInt_trk_pt_%s_ref_unf_%s_%s_%s", dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);
            h->Sumw2 ();
            h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar] = h;

            double totalJetsUF = 0;
            for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

              if (0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]) < minJetPt || 300 < 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1])) continue;

              totalJetsUF += h_jet_pt_ref_unf[iDType][iVar]->GetBinContent (iPtJ+1);

              for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {

                h->SetBinContent (iPtCh+1, h->GetBinContent (iPtCh+1) + h2_unf->GetBinContent (iPtCh+1, iPtJ+1));
                h->SetBinError (iPtCh+1, h->GetBinError (iPtCh+1) + std::pow (h2_unf->GetBinError (iPtCh+1, iPtJ+1), 2));

              } // end loop over iPtCh

            } // end loop over iPtJ

            SaferDelete (&h2_unf);

            for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {

              h->SetBinContent (iPtCh+1, h->GetBinContent (iPtCh+1) / (totalJetsUF * h->GetBinWidth (iPtCh+1)));
              h->SetBinError (iPtCh+1, std::sqrt (h->GetBinError (iPtCh+1)) / (totalJetsUF * h->GetBinWidth (iPtCh+1)));

            } // end loop over iPtCh

          } // end loop over iPtJInt

          SaferDelete (&h2);

        } // end loop over iDir


        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent)); 
          for (short iDir = 0; iDir < nDir; iDir++) {

            const TString dir = directions[iDir];

            TH2D* h2 = new TH2D ("h2_temp", "", nPtChBins, pTChBins, nPtJBins, pTJBins);
            h2->Sumw2 ();

            for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

              TH1D* h = h_jet_trk_pt_sig[iDType][iPtJ][iDir][iCent][iVar];
              const float nJet = h_jet_pt[iDType][iCent][iVar]->GetBinContent (iPtJ+1);

              for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {

                h2->SetBinContent (iPtCh+1, iPtJ+1, nJet * h->GetBinContent (iPtCh+1) * h->GetBinWidth (iPtCh+1));
                h2->SetBinError   (iPtCh+1, iPtJ+1, nJet * h->GetBinError   (iPtCh+1) * h->GetBinWidth (iPtCh+1));

              } // end loop over iPtCh

            } // end loop over iPtJ


            if (iDType == 0 && iVar == 0) {

              for (short iPtJInt : {0, 1}) {

                g_jetInt_trk_pt_unfSumUnc[iPtJInt][iDir][iCent]     = new TGraph (nItersMax - nItersMin + 1);
                g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent]    = new TGraph (nItersMax - nItersMin + 1);
                g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent]     = new TGraph (nItersMax - nItersMin + 1);
                g_jetInt_trk_pt_unfSumRelUnc[iPtJInt][iDir][iCent]  = new TGraph (nItersMax - nItersMin + 1);
                g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent] = new TGraph (nItersMax - nItersMin + 1);
                g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent]  = new TGraph (nItersMax - nItersMin + 1);

              } // end loop over iPtJInt

              TH2D* h2_unf = nullptr;
              TH2D* h2_unf_prev = (TH2D*) h2->Clone ("h2_unf_0Iters");
  
              for (short iIter = 0; iIter < nItersMax-nItersMin+2; iIter++) {

                const short nIters = (short) nItersVals[iIter];

                RooUnfoldBayes* bayesUnf = new RooUnfoldBayes (rooUnfResp_jet_trk_pt_sig[iDir][iCent], h2, nIters);
                bayesUnf->SetVerbose (-1);
                h2_unf = (TH2D*) bayesUnf->Hreco ()->Clone (Form ("h2_unf_%iIters", nIters));
                SaferDelete (&bayesUnf);

                for (short iPtJInt : {0, 1}) {

                  const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
                  const float minJetPt = (iPtJInt == 0 ? 30. : 60.);

                  TH1D* h = new TH1D (Form ("h_jet_trk_pt_%s_pPb_unf_%s_data_%s_Nominal_nIters%i", dir.Data (), cent.Data (), pTJInt.Data (), nIters), "", nPtChBins, pTChBins);
                  h->Sumw2 ();
                  h_jet_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter] = h;

                  double totalJetsUF = 0;
                  for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
    
                    if (0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]) < minJetPt || 300 < 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1])) continue;
    
                    totalJetsUF += h_jet_pt_unf[iDType][iCent][iVar]->GetBinContent (iPtJ+1);
    
                    for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
    
                      h->SetBinContent (iPtCh+1, h->GetBinContent (iPtCh+1) + h2_unf->GetBinContent (iPtCh+1, iPtJ+1));
                      h->SetBinError (iPtCh+1, h->GetBinError (iPtCh+1) + std::pow (h2_unf->GetBinError (iPtCh+1, iPtJ+1), 2));
    
                    } // end loop over iPtCh
    
                  } // end loop over iPtJ

                  for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
  
                    h->SetBinContent (iPtCh+1, h->GetBinContent (iPtCh+1) / (totalJetsUF * h->GetBinWidth (iPtCh+1)));
                    h->SetBinError (iPtCh+1, std::sqrt (h->GetBinError (iPtCh+1)) / (totalJetsUF * h->GetBinWidth (iPtCh+1)));
  
                  } // end loop over iPtCh

                } // end loop over iPtJInt

                double totErr[2] = {0, 0};
                double totRelErr[2] = {0, 0};
                for (int iX = 1; iX <= h2_unf->GetNbinsX (); iX++) {
                  const double ptch = h2_unf->GetXaxis ()->GetBinCenter (iX);
                  for (int iY = 1; iY <= h2_unf->GetNbinsY (); iY++) {
                    for (short iPtJInt : {0, 1}) {
                      const double jpt = h2_unf->GetYaxis ()->GetBinCenter (iY);
                      if (ptch < 5 || ptch > (iPtJInt == 0 ? 40 : 75)) continue;
                      totErr[iPtJInt] += ((iPtJInt == 0 ? 30. : 60.) < jpt && jpt < 300 ? std::pow (h2_unf->GetBinError (iX, iY), 2) : 0);
                      totRelErr[iPtJInt] += ((iPtJInt == 0 ? 30. : 60.) < jpt && jpt < 300 && h2_unf->GetBinContent (iX, iY) != 0 ? std::pow (h2_unf->GetBinError (iX, iY), 2) / h2_unf->GetBinContent (iX, iY) : 0);
                    }
                  } // end loop over iY
                } // end loop over iX
 
                double iterErr[2] = {0, 0};
                double iterRelErr[2] = {0, 0};
                for (int iX = 1; iX <= h2_unf->GetNbinsX (); iX++) {
                  const double ptch = h2_unf->GetXaxis ()->GetBinCenter (iX);
                  for (int iY = 1; iY <= h2_unf->GetNbinsY (); iY++) {
                    for (short iPtJInt : {0, 1}) {
                      const double jpt = h2_unf->GetYaxis ()->GetBinCenter (iY);
                      if (ptch < 5 || ptch > (iPtJInt == 0 ? 40 : 75)) continue;
                      iterErr[iPtJInt] += ((iPtJInt == 0 ? 30. : 60.) < jpt && jpt < 300 ? std::fabs (h2_unf->GetBinContent (iX, iY) - h2_unf_prev->GetBinContent (iX, iY)) : 0);
                      iterRelErr[iPtJInt] += ((iPtJInt == 0 ? 30. : 60.) < jpt && jpt < 300 && h2_unf->GetBinContent (iX, iY) != 0 ? std::fabs (h2_unf->GetBinContent (iX, iY) - h2_unf_prev->GetBinContent (iX, iY)) / std::sqrt (h2_unf->GetBinContent (iX, iY)) : 0);
                    }
                  } // end loop over iY
                } // end loop over iX
  
                for (short iPtJInt : {0, 1}) {

                  g_jetInt_trk_pt_unfSumUnc[iPtJInt][iDir][iCent]->SetPoint     (nIters - nItersMin, nIters, std::sqrt (totErr[iPtJInt]));
                  g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent]->SetPoint    (nIters - nItersMin, nIters, std::sqrt (std::fabs (iterErr[iPtJInt])));
                  g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent]->SetPoint     (nIters - nItersMin, nIters, std::sqrt (totErr[iPtJInt] + std::fabs (iterErr[iPtJInt])));
                  g_jetInt_trk_pt_unfSumRelUnc[iPtJInt][iDir][iCent]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (totRelErr[iPtJInt]));
                  g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent]->SetPoint (nIters - nItersMin, nIters, std::sqrt (std::fabs (iterRelErr[iPtJInt])));
                  g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (totRelErr[iPtJInt] + std::fabs (iterRelErr[iPtJInt])));

                } // end loop over iPtJInt

                SaferDelete (&h2_unf_prev);
                h2_unf_prev = h2_unf;

              } // end loop over iIter

              SaferDelete (&h2_unf);


              double totErr[2] = {0, 0};
              double totRelErr[2] = {0, 0};
              for (short iPtJInt : {0, 1}) {
  
                for (int iX = 1; iX <= h2->GetNbinsX (); iX++) {
                  const double ptch = h2->GetXaxis ()->GetBinCenter (iX);
                  for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {
                    for (short iPtJInt : {0, 1}) {
                      const double jpt = h2->GetYaxis ()->GetBinCenter (iY);
                      if (ptch < 5 || ptch > (iPtJInt == 0 ? 40 : 75)) continue;
                      totErr[iPtJInt] += ((iPtJInt == 0 ? 30. : 60.) < jpt && jpt < 300 ? std::pow (h2->GetBinError (iX, iY), 2) : 0);
                      totRelErr[iPtJInt] += ((iPtJInt == 0 ? 30. : 60.) < jpt && jpt < 300 && h2->GetBinContent (iX, iY) != 0 ? std::pow (h2->GetBinError (iX, iY), 2) / h2->GetBinContent (iX, iY) : 0);
                    }
                  } // end loop over iY
                } // end loop over iX
  
                g_jetInt_trk_pt_unfSumUnc[iPtJInt][iDir][iCent]->SetPoint     (g_jetInt_trk_pt_unfSumUnc[iPtJInt][iDir][iCent]->GetN (),    0, std::sqrt (totErr[iPtJInt]));
                g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent]->SetPoint     (g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent]->GetN (),    0, std::sqrt (totErr[iPtJInt]));
                g_jetInt_trk_pt_unfSumRelUnc[iPtJInt][iDir][iCent]->SetPoint  (g_jetInt_trk_pt_unfSumRelUnc[iPtJInt][iDir][iCent]->GetN (), 0, std::sqrt (totRelErr[iPtJInt]));
                g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent]->SetPoint  (g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent]->GetN (), 0, std::sqrt (totRelErr[iPtJInt]));
  
              } // end loop over iPtJInt
  
            }


            for (short iPtJInt : {0, 1}) {

              const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
              const float minJetPt = (iPtJInt == 0 ? 30. : 60.);

              RooUnfoldBayes* bayesUnf = new RooUnfoldBayes (rooUnfResp_jet_trk_pt_sig[iDir][iCent], h2, GetTrkSpectraNIters (iPtJInt, iDir, iCent));
              bayesUnf->SetVerbose (-1);
              TH2D* h2_unf = (TH2D*) bayesUnf->Hreco ()->Clone ("h2_unf");
              SaferDelete (&bayesUnf);

              TH1D* h = new TH1D (Form ("h_jetInt_trk_pt_%s_pPb_unf_%s_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);
              h->Sumw2 ();
              h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar] = h;

              double totalJetsUF = 0;
              for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

                if (0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]) < minJetPt || 300 < 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1])) continue;

                totalJetsUF += h_jet_pt_unf[iDType][iCent][iVar]->GetBinContent (iPtJ+1);

                for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {

                  h->SetBinContent (iPtCh+1, h->GetBinContent (iPtCh+1) + h2_unf->GetBinContent (iPtCh+1, iPtJ+1));
                  h->SetBinError (iPtCh+1, h->GetBinError (iPtCh+1) + std::pow (h2_unf->GetBinError (iPtCh+1, iPtJ+1), 2));

                } // end loop over iPtCh

              } // end loop over iPtJ

              SaferDelete (&h2_unf);

              for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {

                h->SetBinContent (iPtCh+1, h->GetBinContent (iPtCh+1) / (totalJetsUF * h->GetBinWidth (iPtCh+1)));
                h->SetBinError (iPtCh+1, std::sqrt (h->GetBinError (iPtCh+1)) / (totalJetsUF * h->GetBinWidth (iPtCh+1)));

              } // end loop over iPtCh

            } // end loop over iPtJInt

            SaferDelete (&h2);

          } // end loop over iDir

        } // end loop over iCent

      }
      //else {

      //  for (short iDir = 0; iDir < nDir; iDir++) {

      //    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

      //      h_jet_trk_pt_ref_unf[iDType][iPtJ][iDir][iVar] = h_jet_trk_pt_ref_sig[iDType][iPtJ][iDir][iVar];

      //    } // end loop over iPtJ

      //  } // end loop over iDir

      //  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      //    for (short iDir = 0; iDir < nDir; iDir++) {

      //      for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

      //        h_jet_trk_pt_unf[iDType][iPtJ][iDir][iCent][iVar] = h_jet_trk_pt_sig[iDType][iPtJ][iDir][iCent][iVar];

      //      } // end loop over iPtJ

      //    } // end loop over iDir

      //  } // end loop over iCent

      //}

    } // end loop over iVar

  } // end loop over iDType




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // INTEGRATE OVER JET PT SELECTIONS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iDType = 0; iDType < 2; iDType++) {

    const TString dType = (iDType == 0 ? "data" : "mc");

    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
      const float minJetPt = (iPtJInt == 0 ? 30. : 60.);

      for (short iVar = 0; iVar < nVar; iVar++) {

        const TString var = variations[iVar];

        if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
          continue;

        const bool doUnfold = (variationsWithNoUnfold.count (var) == 0);

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          {
            h_jetInt_trk_pt_ref[iDType][iPtJInt][iDir][iVar]     = new TH1D (Form ("h_jetInt_trk_pt_%s_ref_%s_%s_%s",     dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);
            h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][iDir][iVar] = new TH1D (Form ("h_jetInt_trk_pt_%s_ref_bkg_%s_%s_%s", dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);
            h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][iVar] = new TH1D (Form ("h_jetInt_trk_pt_%s_ref_sig_%s_%s_%s", dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);
            //h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar] = new TH1D (Form ("h_jetInt_trk_pt_%s_ref_unf_%s_%s_%s", dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);

            h_jetInt_trk_pt_ref[iDType][iPtJInt][iDir][iVar]->Sumw2 ();
            h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][iDir][iVar]->Sumw2 ();
            h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][iVar]->Sumw2 ();
            //h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar]->Sumw2 ();

            double totalJets = 0;//, totalJetsUF = 0;
            for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

              if (0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]) < minJetPt || 300 < 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1])) continue;

              const double nJets = h_jet_pt_ref[iDType][iVar]->GetBinContent (iPtJ+1);
              h_jetInt_trk_pt_ref[iDType][iPtJInt][iDir][iVar]->Add (h_jet_trk_pt_ref[iDType][iPtJ][iDir][iVar], nJets);
              h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][iDir][iVar]->Add (h_jet_trk_pt_ref_bkg[iDType][iPtJ][iDir][iVar], nJets);
              h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][iVar]->Add (h_jet_trk_pt_ref_sig[iDType][iPtJ][iDir][iVar], nJets);
              totalJets += nJets;

              //const double nJetsUF = h_jet_pt_ref_unf[iDType][iVar]->GetBinContent (iPtJ+1);
              //h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar]->Add (h_jet_trk_pt_ref_unf[iDType][iPtJ][iDir][iVar], nJetsUF);
              //totalJetsUF += nJetsUF;

            } // end loop over iPtJ

            h_jetInt_trk_pt_ref[iDType][iPtJInt][iDir][iVar]->Scale (1./totalJets);
            h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][iDir][iVar]->Scale (1./totalJets);
            h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][iVar]->Scale (1./totalJets);
            //h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar]->Scale (1./totalJetsUF);

            if (!doUnfold)
              h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar] = (TH1D*) h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][iVar]->Clone (Form ("h_jetInt_trk_pt_%s_ref_unf_%s_%s_%s", dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
          }

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
    
            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            h_jetInt_trk_pt[iDType][iPtJInt][iDir][iCent][iVar]     = new TH1D (Form ("h_jetInt_trk_pt_%s_pPb_%s_%s_%s_%s",     dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);
            h_jetInt_trk_pt_bkg[iDType][iPtJInt][iDir][iCent][iVar] = new TH1D (Form ("h_jetInt_trk_pt_%s_pPb_bkg_%s_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);
            h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][iVar] = new TH1D (Form ("h_jetInt_trk_pt_%s_pPb_sig_%s_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);
            //h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar] = new TH1D (Form ("h_jetInt_trk_pt_%s_pPb_unf_%s_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);

            h_jetInt_trk_pt[iDType][iPtJInt][iDir][iCent][iVar]->Sumw2 ();
            h_jetInt_trk_pt_bkg[iDType][iPtJInt][iDir][iCent][iVar]->Sumw2 ();
            h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][iVar]->Sumw2 ();
            //h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar]->Sumw2 ();

            double totalJets = 0;//, totalJetsUF = 0;
            for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

              if (0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]) < minJetPt || 300 < 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1])) continue;

              const double nJets = h_jet_pt[iDType][iCent][iVar]->GetBinContent (iPtJ+1);
              h_jetInt_trk_pt[iDType][iPtJInt][iDir][iCent][iVar]->Add (h_jet_trk_pt[iDType][iPtJ][iDir][iCent][iVar], nJets);
              h_jetInt_trk_pt_bkg[iDType][iPtJInt][iDir][iCent][iVar]->Add (h_jet_trk_pt_bkg[iDType][iPtJ][iDir][iCent][iVar], nJets);
              h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][iVar]->Add (h_jet_trk_pt_sig[iDType][iPtJ][iDir][iCent][iVar], nJets);
              totalJets += nJets;

              //const double nJetsUF = h_jet_pt_unf[iDType][iCent][iVar]->GetBinContent (iPtJ+1);
              //h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar]->Add (h_jet_trk_pt_unf[iDType][iPtJ][iDir][iCent][iVar], nJetsUF);
              //totalJetsUF += nJetsUF;

            } // end loop over iPtJ

            h_jetInt_trk_pt[iDType][iPtJInt][iDir][iCent][iVar]->Scale (1./totalJets);
            h_jetInt_trk_pt_bkg[iDType][iPtJInt][iDir][iCent][iVar]->Scale (1./totalJets);
            h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][iVar]->Scale (1./totalJets);
            //h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar]->Scale (1./totalJetsUF);

            if (!doUnfold)
              h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar] = (TH1D*) h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][iVar]->Clone (Form ("h_jetInt_trk_pt_%s_pPb_unf_%s_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iDir


        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
  
          const TString ptch = pTChSelections[iPtCh].Data ();

          {
            h_jetInt_trk_dphi_ref[iDType][iPtJInt][iPtCh][iVar]     = new TH1D (Form ("h_jetInt_trk_dphi_%s_ref_%s_%s_%s",     ptch.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nDPhiBins, dPhiBins);
            h_jetInt_trk_dphi_ref_bkg[iDType][iPtJInt][iPtCh][iVar] = new TH1D (Form ("h_jetInt_trk_dphi_%s_ref_bkg_%s_%s_%s", ptch.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nDPhiBins, dPhiBins);
            h_jetInt_trk_dphi_ref_sig[iDType][iPtJInt][iPtCh][iVar] = new TH1D (Form ("h_jetInt_trk_dphi_%s_ref_sig_%s_%s_%s", ptch.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nDPhiBins, dPhiBins);

            h_jetInt_trk_dphi_ref[iDType][iPtJInt][iPtCh][iVar]->Sumw2 ();
            h_jetInt_trk_dphi_ref_bkg[iDType][iPtJInt][iPtCh][iVar]->Sumw2 ();
            h_jetInt_trk_dphi_ref_sig[iDType][iPtJInt][iPtCh][iVar]->Sumw2 ();

            double totalJets = 0;
            for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

              if (0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]) < minJetPt || 300 < 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1])) continue;

              const double nJets = h_jet_pt_ref[iDType][iVar]->GetBinContent (iPtJ+1);
              h_jetInt_trk_dphi_ref[iDType][iPtJInt][iPtCh][iVar]->Add (h_jet_trk_dphi_ref[iDType][iPtJ][iPtCh][iVar], nJets);
              h_jetInt_trk_dphi_ref_bkg[iDType][iPtJInt][iPtCh][iVar]->Add (h_jet_trk_dphi_ref_bkg[iDType][iPtJ][iPtCh][iVar], nJets);
              h_jetInt_trk_dphi_ref_sig[iDType][iPtJInt][iPtCh][iVar]->Add (h_jet_trk_dphi_ref_sig[iDType][iPtJ][iPtCh][iVar], nJets);
              totalJets += nJets;

            } // end loop over iPtJ

            h_jetInt_trk_dphi_ref[iDType][iPtJInt][iPtCh][iVar]->Scale (1./totalJets);
            h_jetInt_trk_dphi_ref_bkg[iDType][iPtJInt][iPtCh][iVar]->Scale (1./totalJets);
          }

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
    
            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            h_jetInt_trk_dphi[iDType][iPtJInt][iPtCh][iCent][iVar]     = new TH1D (Form ("h_jetInt_trk_dphi_%s_pPb_%s_%s_%s_%s",     ptch.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nDPhiBins, dPhiBins);
            h_jetInt_trk_dphi_bkg[iDType][iPtJInt][iPtCh][iCent][iVar] = new TH1D (Form ("h_jetInt_trk_dphi_%s_pPb_bkg_%s_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nDPhiBins, dPhiBins);
            h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent][iVar] = new TH1D (Form ("h_jetInt_trk_dphi_%s_pPb_sig_%s_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nDPhiBins, dPhiBins);

            h_jetInt_trk_dphi[iDType][iPtJInt][iPtCh][iCent][iVar]->Sumw2 ();
            h_jetInt_trk_dphi_bkg[iDType][iPtJInt][iPtCh][iCent][iVar]->Sumw2 ();
            h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent][iVar]->Sumw2 ();

            double totalJets = 0;
            for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

              if (0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]) < minJetPt || 300 < 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1])) continue;

              const double nJets = h_jet_pt[iDType][iCent][iVar]->GetBinContent (iPtJ+1);
              h_jetInt_trk_dphi[iDType][iPtJInt][iPtCh][iCent][iVar]->Add (h_jet_trk_dphi[iDType][iPtJ][iPtCh][iCent][iVar], nJets);
              h_jetInt_trk_dphi_bkg[iDType][iPtJInt][iPtCh][iCent][iVar]->Add (h_jet_trk_dphi_bkg[iDType][iPtJ][iPtCh][iCent][iVar], nJets);
              h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent][iVar]->Add (h_jet_trk_dphi_sig[iDType][iPtJ][iPtCh][iCent][iVar], nJets);
              totalJets += nJets;

            } // end loop over iPtJ

            h_jetInt_trk_dphi[iDType][iPtJInt][iPtCh][iCent][iVar]->Scale (1./totalJets);
            h_jetInt_trk_dphi_bkg[iDType][iPtJInt][iPtCh][iCent][iVar]->Scale (1./totalJets);
            h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent][iVar]->Scale (1./totalJets);

          } // end loop over iCent

        } // end loop over iPtCh

      } // end loop over iVar

    } // end loop over iPtJInt

  } // end loop over iDType




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CALCULATE IAA RATIOS FOR INTEGRATED JET PT SELECTIONS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iDType = 0; iDType < 2; iDType++) {

    const TString dType = (iDType == 0 ? "data" : "mc");

    for (short iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
        continue;

      for (short iPtJInt : {0, 1}) {

        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent][iVar] = (TH1D*) h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar]->Clone (Form ("h_jetInt_trk_pt_%s_iaa_%s_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
            h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent][iVar]->Divide (h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar]);

            h_jetInt_trk_pt_iaaNoUnf[iDType][iPtJInt][iDir][iCent][iVar] = (TH1D*) h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][iVar]->Clone (Form ("h_jetInt_trk_pt_%s_iaaNoUnf_%s_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
            h_jetInt_trk_pt_iaaNoUnf[iDType][iPtJInt][iDir][iCent][iVar]->Divide (h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][iVar]);

          } // end loop over iCent

        } // end loop over iDir


        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          const TString ptch = pTChSelections[iPtCh].Data ();

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            h_jetInt_trk_dphi_iaa[iDType][iPtJInt][iPtCh][iCent][iVar] = (TH1D*) h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent][iVar]->Clone (Form ("h_jetInt_trk_dphi_%s_iaa_%s_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
            h_jetInt_trk_dphi_iaa[iDType][iPtJInt][iPtCh][iCent][iVar]->Divide (h_jetInt_trk_dphi_ref_sig[iDType][iPtJInt][iPtCh][iVar]);

          } // end loop over iCent

        } // end loop over iPtCh

      } // end loop over iPtJInt

    } // end loop over iVar

  } // end loop over iDType




  ////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //// CALCULATE MC CLOSURE AND SMOOTH OVER FLUCTUATIONS
  ////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //{
  //  const short iVarNum = 0;//GetVarN ("Nominal");
  //  const short iVarDen = GetVarN ("MCTruthJetsTruthParts");

  //  TF1* f = nullptr;

  //  const TString funcStr = "[0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)+[4]*pow(log(x),4)";

  //  //bbbf.SetNDeriv (2);
  //  //bbbf.SetDegree (6);

  //  //const int ndf = bbbf.NDF ();

  //  for (short iPtJInt : {0, 1}) {

  //    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

  //    for (short iDir = 0; iDir < nDir; iDir++) {

  //      const TString dir = directions[iDir];

  //      TH1D* h = nullptr;

  //      h = (TH1D*) h_jetInt_trk_pt_ref_unf[1][iPtJInt][iDir][iVarNum]->Clone (Form ("h_jetInt_trk_pt_%s_ref_unf_closure_%s", dir.Data (), pTJInt.Data ()));
  //      h->Divide (h_jetInt_trk_pt_ref_unf[1][iPtJInt][iDir][iVarDen]);
  //      h_jetInt_trk_pt_ref_unf_closure[iPtJInt][iDir] = h;

  //      f = new TF1 (Form ("f_jetInt_trk_pt_%s_%s_ref_unf_closure", dir.Data (), pTJInt.Data ()), funcStr, h->GetXaxis ()->GetBinLowEdge (1), h->GetXaxis ()->GetBinLowEdge (h->GetNbinsX () + 1));
  //      if (iDir == 1) {
  //        double x1, x2;
  //        f->GetRange (x1, x2);
  //        f->SetRange (x1, 20);
  //      }
  //      //f->SetParameter (0, 100);
  //      //f->SetParameter (1, 0);
  //      //f->SetParameter (2, 0);

  //      f_jetInt_trk_pt_ref_unf_closure[iPtJInt][iDir] = f;
  //      h->Fit (f, "RN0Q");

  //      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

  //        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

  //        h = (TH1D*) h_jetInt_trk_pt_unf[1][iPtJInt][iDir][iCent][iVarNum]->Clone (Form ("h_jet_trk_pt_%s_pPb_unf_closure_%s_%s", dir.Data (), cent.Data (), pTJInt.Data ()));
  //        h->Divide (h_jetInt_trk_pt_unf[1][iPtJInt][iDir][iCent][iVarDen]);
  //        h_jetInt_trk_pt_unf_closure[iPtJInt][iDir][iCent] = h;

  //        f = new TF1 (Form ("f_jetInt_trk_pt_%s_pPb_unf_closure_%s_%s", dir.Data (), cent.Data (), pTJInt.Data ()), funcStr, h->GetXaxis ()->GetBinLowEdge (1), h->GetXaxis ()->GetBinLowEdge (h->GetNbinsX () + 1));
  //        if (iDir == 1) {
  //          double x1, x2;
  //          f->GetRange (x1, x2);
  //          f->SetRange (x1, 20);
  //        }
  //        //f->SetParameter (0, 100);
  //        //f->SetParameter (1, 0);
  //        //f->SetParameter (2, 0);

  //        f_jetInt_trk_pt_unf_closure[iPtJInt][iDir][iCent] = f;
  //        h->Fit (f, "RN0Q");


  //        h = (TH1D*) h_jetInt_trk_pt_iaa[1][iPtJInt][iDir][iCent][iVarNum]->Clone (Form ("h_jet_trk_pt_%s_iaa_closure_%s_%s", dir.Data (), cent.Data (), pTJInt.Data ()));
  //        h->Divide (h_jetInt_trk_pt_iaa[1][iPtJInt][iDir][iCent][iVarDen]);
  //        h_jetInt_trk_pt_iaa_closure[iPtJInt][iDir][iCent] = h;

  //        f = new TF1 (Form ("f_jetInt_trk_pt_%s_iaa_closure_%s_%s", dir.Data (), cent.Data (), pTJInt.Data ()), funcStr, h->GetXaxis ()->GetBinLowEdge (1), h->GetXaxis ()->GetBinLowEdge (h->GetNbinsX () + 1));
  //        if (iDir == 1) {
  //          double x1, x2;
  //          f->GetRange (x1, x2);
  //          f->SetRange (x1, 20);
  //        }
  //        //f->SetParameter (0, 100);
  //        //f->SetParameter (1, 0);
  //        //f->SetParameter (2, 0);

  //        f_jetInt_trk_pt_iaa_closure[iPtJInt][iDir][iCent] = f;
  //        h->Fit (f, "RN0Q");

  //      } // end loop over iCent

  //    } // end loop over iDir

  //  } // end loop over iPtJInt
  //}




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CORRECT HISTOGRAMS BY NON-CLOSURE FITS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iDType = 0; iDType < 2; iDType++) {

    for (short iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if ((variationsWithNoUnfold.count (var) > 0) || (iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
        continue;

      const float scaleFactor = (var == "NonClosureVar" ? 1.0 : 0.5); // for non-closure systematic correct by 100% of non-closure; nominal only corrects by 50%

      for (short iPtJInt : {0, 1}) {

        for (short iDir = 0; iDir < nDir; iDir++) {

          DivideByTF1 (h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar], f_jetInt_trk_pt_ref_unf_closure[iPtJInt][iDir], scaleFactor);

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            DivideByTF1 (h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar], f_jetInt_trk_pt_unf_closure[iPtJInt][iDir][iCent], scaleFactor);
            DivideByTF1 (h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent][iVar], f_jetInt_trk_pt_iaa_closure[iPtJInt][iDir][iCent], scaleFactor);

          } // end loop over iCent

        } // end loop over iDir

      } // end loop over iPtJInt

    } // end loop over iVar

  } // end loop over iDType




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CREATE GRAPHS THAT WILL STORE TOTAL SYSTEMATIC UNCERTAINTIES
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iPtJInt : {0, 1}) {

    for (short iDir = 0; iDir < nDir; iDir++) {

      g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][0]     = make_graph (h_jetInt_trk_pt_ref[0][iPtJInt][iDir][0]);
      g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][0] = make_graph (h_jetInt_trk_pt_ref_bkg[0][iPtJInt][iDir][0]);
      g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][0] = make_graph (h_jetInt_trk_pt_ref_sig[0][iPtJInt][iDir][0]);
      g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][0] = make_graph (h_jetInt_trk_pt_ref_unf[0][iPtJInt][iDir][0]);

      ResetTGAEErrors (g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][0]);
      ResetTGAEErrors (g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][0]);
      ResetTGAEErrors (g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][0]);
      ResetTGAEErrors (g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][0]);

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][0]        = make_graph (h_jetInt_trk_pt[0][iPtJInt][iDir][iCent][0]);
        g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][0]    = make_graph (h_jetInt_trk_pt_bkg[0][iPtJInt][iDir][iCent][0]);
        g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][0]    = make_graph (h_jetInt_trk_pt_sig[0][iPtJInt][iDir][iCent][0]);
        g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][0]    = make_graph (h_jetInt_trk_pt_unf[0][iPtJInt][iDir][iCent][0]);
        g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0]    = make_graph (h_jetInt_trk_pt_iaa[0][iPtJInt][iDir][iCent][0]);

        ResetTGAEErrors (g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][0]);
        ResetTGAEErrors (g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][0]);
        ResetTGAEErrors (g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][0]);
        ResetTGAEErrors (g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][0]);
        ResetTGAEErrors (g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0]);

      } // end loop over iCent

    } // end loop over iDir

    for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

      g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][0]     = make_graph (h_jetInt_trk_dphi_ref[0][iPtJInt][iPtCh][0]);
      g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][0] = make_graph (h_jetInt_trk_dphi_ref_bkg[0][iPtJInt][iPtCh][0]);
      g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][0] = make_graph (h_jetInt_trk_dphi_ref_sig[0][iPtJInt][iPtCh][0]);

      ResetTGAEErrors (g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][0]);
      ResetTGAEErrors (g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][0]);
      ResetTGAEErrors (g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][0]);

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][0]      = make_graph (h_jetInt_trk_dphi[0][iPtJInt][iPtCh][iCent][0]);
        g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][0]  = make_graph (h_jetInt_trk_dphi_bkg[0][iPtJInt][iPtCh][iCent][0]);
        g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][0]  = make_graph (h_jetInt_trk_dphi_sig[0][iPtJInt][iPtCh][iCent][0]);
        g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][0]  = make_graph (h_jetInt_trk_dphi_iaa[0][iPtJInt][iPtCh][iCent][0]);

        ResetTGAEErrors (g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][0]);
        ResetTGAEErrors (g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][0]);
        ResetTGAEErrors (g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][0]);
        ResetTGAEErrors (g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][0]);

      } // end loop over iCent

    } // end loop over iPtCh

  } // end loop over iPtJInt




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CREATE GRAPHS THAT WILL STORE SYSTEMATIC UNCERTAINTIES FROM EACH SOURCE SEPARATELY.
  // THEN CALCULATES SYSTEMATIC UNCERTAINTIES FROM EACH SOURCE BY TAKING DIFFERENCES
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iDType = 0; iDType < 2; iDType++) {

    for (short iPtJInt : {0, 1}) {

      for (short iVar = 1; iVar < nVar; iVar++) {

        const TString var = variations[iVar];

        if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) == 0))
          continue;

        for (short iDir = 0; iDir < nDir; iDir++) {

          g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][iVar]     = new TGAE ();
          g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][iVar] = new TGAE ();
          g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar] = new TGAE ();
          g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar] = new TGAE ();

          CalcSystematics (g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][iVar],     h_jetInt_trk_pt_ref[iDType][iPtJInt][iDir][0],      h_jetInt_trk_pt_ref[iDType][iPtJInt][iDir][iVar]);
          CalcSystematics (g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][iVar], h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][iDir][0],  h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][iDir][iVar]);
          CalcSystematics (g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar], h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][0],  h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][iVar]);
          CalcSystematics (g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar], h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][0],  h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar]);

          // allow some uncertainties to not cancel in the bkgd. subtracted (signal) yield by overwriting the current uncertainties
          if (variationsThatDontCancelInSig.count (var) != 0) {
            ResetTGAEErrors (g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar]);
            AddErrorsInQuadrature (g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar], g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][iVar], false, false);
            AddErrorsInQuadrature (g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar], g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][iVar], false, false);
          }

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][iVar]     = new TGAE ();
            g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][iVar] = new TGAE ();
            g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar] = new TGAE ();
            g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar] = new TGAE ();
            g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar] = new TGAE ();

            CalcSystematics (g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][iVar],     h_jetInt_trk_pt[iDType][iPtJInt][iDir][iCent][0],     h_jetInt_trk_pt[iDType][iPtJInt][iDir][iCent][iVar]);
            CalcSystematics (g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][iVar], h_jetInt_trk_pt_bkg[iDType][iPtJInt][iDir][iCent][0], h_jetInt_trk_pt_bkg[iDType][iPtJInt][iDir][iCent][iVar]);
            CalcSystematics (g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar], h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][0], h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][iVar]);
            CalcSystematics (g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar], h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][0], h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar]);
            CalcSystematics (g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar], h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent][0], h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent][iVar]);

            // allow some uncertainties to not cancel in the bkgd. subtracted (signal) yield by overwriting the current uncertainties
            if (variationsThatDontCancelInSig.count (var) != 0) {
              ResetTGAEErrors (g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar]);
              AddErrorsInQuadrature (g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar], g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][iVar], false, false);
              AddErrorsInQuadrature (g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar], g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][iVar], false, false);
            }
            // allow some uncertainties to not cancel in the p+Pb / pp ratio by overwriting the current uncertainties
            if (variationsThatDontCancelInRatio.count (var) != 0) {
              ResetTGAEErrors (g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar]);
              AddRelErrorsInQuadrature (g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar], g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar], false, false);
              AddRelErrorsInQuadrature (g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar], g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar], false, false);
            }

          } // end loop over iCent

        } // end loop over iDir

        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][iVar]      = new TGAE ();
          g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][iVar]  = new TGAE ();
          g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar]  = new TGAE ();

          CalcSystematics (g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][iVar],      h_jetInt_trk_dphi_ref[iDType][iPtJInt][iPtCh][0],     h_jetInt_trk_dphi_ref[iDType][iPtJInt][iPtCh][iVar]);
          CalcSystematics (g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][iVar],  h_jetInt_trk_dphi_ref_bkg[iDType][iPtJInt][iPtCh][0], h_jetInt_trk_dphi_ref_bkg[iDType][iPtJInt][iPtCh][iVar]);
          CalcSystematics (g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar],  h_jetInt_trk_dphi_ref_sig[iDType][iPtJInt][iPtCh][0], h_jetInt_trk_dphi_ref_sig[iDType][iPtJInt][iPtCh][iVar]);

          // allow some uncertainties to not cancel in the bkgd. subtracted (signal) yield by overwriting the current uncertainties
          if (variationsThatDontCancelInSig.count (var) != 0) {
            ResetTGAEErrors (g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar]);
            AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar], g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][iVar], false, false);
            AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar], g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][iVar], false, false);
          }

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][iVar]     = new TGAE ();
            g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][iVar] = new TGAE ();
            g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar] = new TGAE ();
            g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][iVar] = new TGAE ();

            CalcSystematics (g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][iVar],     h_jetInt_trk_dphi[iDType][iPtJInt][iPtCh][iCent][0],      h_jetInt_trk_dphi[iDType][iPtJInt][iPtCh][iCent][iVar]);
            CalcSystematics (g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][iVar], h_jetInt_trk_dphi_bkg[iDType][iPtJInt][iPtCh][iCent][0],  h_jetInt_trk_dphi_bkg[iDType][iPtJInt][iPtCh][iCent][iVar]);
            CalcSystematics (g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar], h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent][0],  h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent][iVar]);
            CalcSystematics (g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][iVar], h_jetInt_trk_dphi_iaa[iDType][iPtJInt][iPtCh][iCent][0],  h_jetInt_trk_dphi_iaa[iDType][iPtJInt][iPtCh][iCent][iVar]);

            // allow some uncertainties to not cancel in the bkgd. subtracted (signal) yield by overwriting the current uncertainties
            if (variationsThatDontCancelInSig.count (var) != 0) {
              ResetTGAEErrors (g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar]);
              AddErrorsInQuadrature (g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar], g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][iVar], false, false);
              AddErrorsInQuadrature (g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar], g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][iVar], false, false);
            }
            // allow some uncertainties to not cancel in the p+Pb / pp ratio by overwriting the current uncertainties
            if (variationsThatDontCancelInRatio.count (var) != 0) {
              ResetTGAEErrors (g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][iVar]);
              AddRelErrorsInQuadrature (g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][iVar], g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar], false, false);
              AddRelErrorsInQuadrature (g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][iVar], g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar], false, false);
            }

          } // end loop over iCent

        } // end loop over iPtCh

      } // end loop over iVar

    } // end loop over iPtJInt

  } // end loop over iDType




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // SYSTEMATIC UNCERTAINTIES DERIVED IN MC MUST HAVE CENTRAL VALUES SET BY CENTRAL VALUES IN DATA
  // THE FINAL UNCERTAINTY IS ASSIGNED TO MATCH THE FRACTIONAL UNCERTAINTY IN MC
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iVar = 1; iVar < nVar; iVar++) {

    const TString var = variations[iVar];

    if (dataVariations.count (var) > 0 || mcVariations.count (var) == 0)
      continue; // skip variations already evaluated in data or that are not evaluated in MC

    for (short iPtJInt : {0, 1}) {

      for (short iDir = 0; iDir < nDir; iDir++) {

        SetCentralValuesKeepRelativeErrors (g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][iVar],      h_jetInt_trk_pt_ref[0][iPtJInt][iDir][0]);
        SetCentralValuesKeepRelativeErrors (g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][iVar],  h_jetInt_trk_pt_ref_bkg[0][iPtJInt][iDir][0]);
        SetCentralValuesKeepRelativeErrors (g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar],  h_jetInt_trk_pt_ref_sig[0][iPtJInt][iDir][0]);
        SetCentralValuesKeepRelativeErrors (g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar],  h_jetInt_trk_pt_ref_unf[0][iPtJInt][iDir][0]);

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          SetCentralValuesKeepRelativeErrors (g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][iVar],     h_jetInt_trk_pt[0][iPtJInt][iDir][iCent][0]);
          SetCentralValuesKeepRelativeErrors (g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][iVar], h_jetInt_trk_pt_bkg[0][iPtJInt][iDir][iCent][0]);
          SetCentralValuesKeepRelativeErrors (g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar], h_jetInt_trk_pt_sig[0][iPtJInt][iDir][iCent][0]);
          SetCentralValuesKeepRelativeErrors (g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar], h_jetInt_trk_pt_unf[0][iPtJInt][iDir][iCent][0]);
          SetCentralValuesKeepRelativeErrors (g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar], h_jetInt_trk_pt_iaa[0][iPtJInt][iDir][iCent][0]);

        } // end loop over iCent

      } // end loop over iDir

      for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

        SetCentralValuesKeepRelativeErrors (g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][iVar],     h_jetInt_trk_dphi_ref[0][iPtJInt][iPtCh][0]);
        SetCentralValuesKeepRelativeErrors (g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][iVar], h_jetInt_trk_dphi_ref_bkg[0][iPtJInt][iPtCh][0]);
        SetCentralValuesKeepRelativeErrors (g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar], h_jetInt_trk_dphi_ref_sig[0][iPtJInt][iPtCh][0]);

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          SetCentralValuesKeepRelativeErrors (g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][iVar],      h_jetInt_trk_dphi[0][iPtJInt][iPtCh][iCent][0]);
          SetCentralValuesKeepRelativeErrors (g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][iVar],  h_jetInt_trk_dphi_bkg[0][iPtJInt][iPtCh][iCent][0]);
          SetCentralValuesKeepRelativeErrors (g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar],  h_jetInt_trk_dphi_sig[0][iPtJInt][iPtCh][iCent][0]);
          SetCentralValuesKeepRelativeErrors (g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][iVar],  h_jetInt_trk_dphi_iaa[0][iPtJInt][iPtCh][iCent][0]);

        } // end loop over iCent

      } // end loop over iPtCh

    } // end loop over iPtJInt
    
  } // end loop over iVar




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // THESE GRAPHS STORE SUMMARY SYSTEMATIC UNCERTAINTIES FOR EACH CATEGORY: TRACKING, JETS, & MIXING.
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iPtJInt : {0, 1}) {

    for (short iTotVar = 0; iTotVar < 3; iTotVar++) {

      for (short iDir = 0; iDir < nDir; iDir++) {

        g_jetInt_trk_pt_ref_systTot[iPtJInt][iDir][iTotVar]     = (TGAE*) g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][0]->Clone ();
        g_jetInt_trk_pt_ref_bkg_systTot[iPtJInt][iDir][iTotVar] = (TGAE*) g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][0]->Clone ();
        g_jetInt_trk_pt_ref_sig_systTot[iPtJInt][iDir][iTotVar] = (TGAE*) g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][0]->Clone ();
        g_jetInt_trk_pt_ref_unf_systTot[iPtJInt][iDir][iTotVar] = (TGAE*) g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][0]->Clone ();

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          g_jetInt_trk_pt_systTot[iPtJInt][iDir][iCent][iTotVar]     = (TGAE*) g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][0]->Clone ();
          g_jetInt_trk_pt_bkg_systTot[iPtJInt][iDir][iCent][iTotVar] = (TGAE*) g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][0]->Clone ();
          g_jetInt_trk_pt_sig_systTot[iPtJInt][iDir][iCent][iTotVar] = (TGAE*) g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][0]->Clone ();
          g_jetInt_trk_pt_unf_systTot[iPtJInt][iDir][iCent][iTotVar] = (TGAE*) g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][0]->Clone ();
          g_jetInt_trk_pt_iaa_systTot[iPtJInt][iDir][iCent][iTotVar] = (TGAE*) g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0]->Clone ();

        } // end loop over iCent

      } // end loop over iDir

      for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

        g_jetInt_trk_dphi_ref_systTot[iPtJInt][iPtCh][iTotVar]     = (TGAE*) g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][0]->Clone ();
        g_jetInt_trk_dphi_ref_bkg_systTot[iPtJInt][iPtCh][iTotVar] = (TGAE*) g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][0]->Clone ();
        g_jetInt_trk_dphi_ref_sig_systTot[iPtJInt][iPtCh][iTotVar] = (TGAE*) g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][0]->Clone ();

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          g_jetInt_trk_dphi_systTot[iPtJInt][iPtCh][iCent][iTotVar]     = (TGAE*) g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][0]->Clone ();
          g_jetInt_trk_dphi_bkg_systTot[iPtJInt][iPtCh][iCent][iTotVar] = (TGAE*) g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][0]->Clone ();
          g_jetInt_trk_dphi_sig_systTot[iPtJInt][iPtCh][iCent][iTotVar] = (TGAE*) g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][0]->Clone ();
          g_jetInt_trk_dphi_iaa_systTot[iPtJInt][iPtCh][iCent][iTotVar] = (TGAE*) g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][0]->Clone ();

        } // end loop over iCent

      } // end loop over iPtCh

    } // end loop over iTotVar

  } // end loop over iPtJInt




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // ADD SYSTEMATIC UNCERTAINTIES FROM EACH CATEGORY INTO ONE OF THREE GRAPHS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iPtJInt : {0, 1}) {

    for (std::vector <TString> vgroup : variationGroups) {

      std::vector <int> iVars = {};
      for (TString s : vgroup) {
        const short iVar = GetVariationN (s);
        if (0 <= iVar && iVar < nVar)
          iVars.push_back (iVar);
      }

      if (iVars.size () == 0)
        continue;

      else if (iVars.size () == 1) {

        const short iVar = iVars[0];
        const TString var = variations[iVar];
        const short iTotVar = (IsJetsVariation (var) ? 2 : (IsTrackingVariation (var) ? 1 : 0));

        for (short iDir = 0; iDir < nDir; iDir++) {

          AddErrorsInQuadrature (g_jetInt_trk_pt_ref_systTot[iPtJInt][iDir][iTotVar],      g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][iVar], false, true);
          AddErrorsInQuadrature (g_jetInt_trk_pt_ref_bkg_systTot[iPtJInt][iDir][iTotVar],  g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][iVar], false, true);
          AddErrorsInQuadrature (g_jetInt_trk_pt_ref_sig_systTot[iPtJInt][iDir][iTotVar],  g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar], false, true);
          AddErrorsInQuadrature (g_jetInt_trk_pt_ref_unf_systTot[iPtJInt][iDir][iTotVar],  g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar], false, true);

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            AddErrorsInQuadrature (g_jetInt_trk_pt_systTot[iPtJInt][iDir][iCent][iTotVar],     g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][iVar], false, true);
            AddErrorsInQuadrature (g_jetInt_trk_pt_bkg_systTot[iPtJInt][iDir][iCent][iTotVar], g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][iVar], false, true);
            AddErrorsInQuadrature (g_jetInt_trk_pt_sig_systTot[iPtJInt][iDir][iCent][iTotVar], g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar], false, true);
            AddErrorsInQuadrature (g_jetInt_trk_pt_unf_systTot[iPtJInt][iDir][iCent][iTotVar], g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar], false, true);
            AddErrorsInQuadrature (g_jetInt_trk_pt_iaa_systTot[iPtJInt][iDir][iCent][iTotVar], g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar], false, true);

          } // end loop over iCent

        } // end loop over iDir

        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_systTot[iPtJInt][iPtCh][iTotVar],     g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][iVar], false, true);
          AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_bkg_systTot[iPtJInt][iPtCh][iTotVar], g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][iVar], false, true);
          AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_sig_systTot[iPtJInt][iPtCh][iTotVar], g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar], false, true);

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            AddErrorsInQuadrature (g_jetInt_trk_dphi_systTot[iPtJInt][iPtCh][iCent][iTotVar],      g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][iVar], false, true);
            AddErrorsInQuadrature (g_jetInt_trk_dphi_bkg_systTot[iPtJInt][iPtCh][iCent][iTotVar],  g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][iVar], false, true);
            AddErrorsInQuadrature (g_jetInt_trk_dphi_sig_systTot[iPtJInt][iPtCh][iCent][iTotVar],  g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar], false, true);
            AddErrorsInQuadrature (g_jetInt_trk_dphi_iaa_systTot[iPtJInt][iPtCh][iCent][iTotVar],  g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][iVar], false, true);

          } // end loop over iCent

        } // end loop over iPtCh

      }

      else {

        const TString var = variations[iVars[0]];
        const short iTotVar = (IsJetsVariation (var) ? 2 : (IsTrackingVariation (var) ? 1 : 0));

        for (short iDir = 0; iDir < nDir; iDir++) {

          AddErrorsInQuadrature (g_jetInt_trk_pt_ref_systTot[iPtJInt][iDir][iTotVar],      g_jetInt_trk_pt_ref_syst[iPtJInt][iDir],     &iVars, true);
          AddErrorsInQuadrature (g_jetInt_trk_pt_ref_bkg_systTot[iPtJInt][iDir][iTotVar],  g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir], &iVars, true);
          AddErrorsInQuadrature (g_jetInt_trk_pt_ref_sig_systTot[iPtJInt][iDir][iTotVar],  g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir], &iVars, true);
          AddErrorsInQuadrature (g_jetInt_trk_pt_ref_unf_systTot[iPtJInt][iDir][iTotVar],  g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir], &iVars, true);

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            AddErrorsInQuadrature (g_jetInt_trk_pt_systTot[iPtJInt][iDir][iCent][iTotVar],     g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent],      &iVars, true);
            AddErrorsInQuadrature (g_jetInt_trk_pt_bkg_systTot[iPtJInt][iDir][iCent][iTotVar], g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent],  &iVars, true);
            AddErrorsInQuadrature (g_jetInt_trk_pt_sig_systTot[iPtJInt][iDir][iCent][iTotVar], g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent],  &iVars, true);
            AddErrorsInQuadrature (g_jetInt_trk_pt_unf_systTot[iPtJInt][iDir][iCent][iTotVar], g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent],  &iVars, true);
            AddErrorsInQuadrature (g_jetInt_trk_pt_iaa_systTot[iPtJInt][iDir][iCent][iTotVar], g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent],  &iVars, true);

          } // end loop over iCent

        } // end loop over iDir

        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_systTot[iPtJInt][iPtCh][iTotVar],     g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh],      &iVars, true);
          AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_bkg_systTot[iPtJInt][iPtCh][iTotVar], g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh],  &iVars, true);
          AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_sig_systTot[iPtJInt][iPtCh][iTotVar], g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh],  &iVars, true);

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            AddErrorsInQuadrature (g_jetInt_trk_dphi_systTot[iPtJInt][iPtCh][iCent][iTotVar],      g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent],     &iVars, true);
            AddErrorsInQuadrature (g_jetInt_trk_dphi_bkg_systTot[iPtJInt][iPtCh][iCent][iTotVar],  g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent], &iVars, true);
            AddErrorsInQuadrature (g_jetInt_trk_dphi_sig_systTot[iPtJInt][iPtCh][iCent][iTotVar],  g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent], &iVars, true);
            AddErrorsInQuadrature (g_jetInt_trk_dphi_iaa_systTot[iPtJInt][iPtCh][iCent][iTotVar],  g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent], &iVars, true);

          } // end loop over iCent

        } // end loop over iPtCh

      }

    } // end loop over iVar

  } // end loop over iPtJInt




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // ADD SYSTEMATIC UNCERTAINTIES FROM ALL SOURCES IN QUADRATURE, STORING RESULTS IN A SINGLE GRAPH
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iPtJInt : {0, 1}) {

    for (short iTotVar = 0; iTotVar < 3; iTotVar++) {

      for (short iDir = 0; iDir < nDir; iDir++) {
    
        AddErrorsInQuadrature (g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][0],      g_jetInt_trk_pt_ref_systTot[iPtJInt][iDir][iTotVar]);
        AddErrorsInQuadrature (g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][0],  g_jetInt_trk_pt_ref_bkg_systTot[iPtJInt][iDir][iTotVar]);
        AddErrorsInQuadrature (g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][0],  g_jetInt_trk_pt_ref_sig_systTot[iPtJInt][iDir][iTotVar]);
        AddErrorsInQuadrature (g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][0],  g_jetInt_trk_pt_ref_unf_systTot[iPtJInt][iDir][iTotVar]);
    
        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
    
          AddErrorsInQuadrature (g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][0],     g_jetInt_trk_pt_systTot[iPtJInt][iDir][iCent][iTotVar]);
          AddErrorsInQuadrature (g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][0], g_jetInt_trk_pt_bkg_systTot[iPtJInt][iDir][iCent][iTotVar]);
          AddErrorsInQuadrature (g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][0], g_jetInt_trk_pt_sig_systTot[iPtJInt][iDir][iCent][iTotVar]);
          AddErrorsInQuadrature (g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][0], g_jetInt_trk_pt_unf_systTot[iPtJInt][iDir][iCent][iTotVar]);
          AddErrorsInQuadrature (g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0], g_jetInt_trk_pt_iaa_systTot[iPtJInt][iDir][iCent][iTotVar]);
    
        } // end loop over iCent
    
      } // end loop over iDir
    
      for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
    
        AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][0],     g_jetInt_trk_dphi_ref_systTot[iPtJInt][iPtCh][iTotVar]);
        AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][0], g_jetInt_trk_dphi_ref_bkg_systTot[iPtJInt][iPtCh][iTotVar]);
        AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][0], g_jetInt_trk_dphi_ref_sig_systTot[iPtJInt][iPtCh][iTotVar]);
    
        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
    
          AddErrorsInQuadrature (g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][0],      g_jetInt_trk_dphi_systTot[iPtJInt][iPtCh][iCent][iTotVar]);
          AddErrorsInQuadrature (g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][0],  g_jetInt_trk_dphi_bkg_systTot[iPtJInt][iPtCh][iCent][iTotVar]);
          AddErrorsInQuadrature (g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][0],  g_jetInt_trk_dphi_sig_systTot[iPtJInt][iPtCh][iCent][iTotVar]);
          AddErrorsInQuadrature (g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][0],  g_jetInt_trk_dphi_iaa_systTot[iPtJInt][iPtCh][iCent][iTotVar]);
    
        } // end loop over iCent
    
      } // end loop over iPtCh

    } // end loop over iTotVar

  } // end loop over iPtJInt




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // FINALLY WRITE OUT EVERYTHING TO A SINGLE ROOT FILE WITH ALL THE RESULTS.
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  {
    outFile->cd ();

    for (short iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if (dataVariations.count (var) == 0 && mcVariations.count (var) == 0)
        continue;

      for (short iPtJInt : {0, 1}) {

        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][iVar]->Write     (Form ("g_jetInt_trk_pt_%s_ref_syst_%s_%s",      dir.Data (), pTJInt.Data (), var.Data ()));
          g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][iVar]->Write (Form ("g_jetInt_trk_pt_%s_ref_bkg_syst_%s_%s",  dir.Data (), pTJInt.Data (), var.Data ()));
          g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar]->Write (Form ("g_jetInt_trk_pt_%s_ref_sig_syst_%s_%s",  dir.Data (), pTJInt.Data (), var.Data ()));
          g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar]->Write (Form ("g_jetInt_trk_pt_%s_ref_unf_syst_%s_%s",  dir.Data (), pTJInt.Data (), var.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][iVar]->Write      (Form ("g_jetInt_trk_pt_%s_syst_%s_%s_%s",     dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][iVar]->Write  (Form ("g_jetInt_trk_pt_%s_bkg_syst_%s_%s_%s", dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar]->Write  (Form ("g_jetInt_trk_pt_%s_sig_syst_%s_%s_%s", dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar]->Write  (Form ("g_jetInt_trk_pt_%s_unf_syst_%s_%s_%s", dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar]->Write  (Form ("g_jetInt_trk_pt_%s_iaa_syst_%s_%s_%s", dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iDir

        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          const TString ptch = pTChSelections[iPtCh].Data ();

          g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][iVar]->Write      (Form ("g_jetInt_trk_dphi_%s_ref_syst_%s_%s",      ptch.Data (), pTJInt.Data (), var.Data ()));
          g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][iVar]->Write  (Form ("g_jetInt_trk_dphi_%s_ref_bkg_syst_%s_%s",  ptch.Data (), pTJInt.Data (), var.Data ()));
          g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar]->Write  (Form ("g_jetInt_trk_dphi_%s_ref_sig_syst_%s_%s",  ptch.Data (), pTJInt.Data (), var.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][iVar]->Write     (Form ("g_jetInt_trk_dphi_%s_syst_%s_%s_%s",     ptch.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][iVar]->Write (Form ("g_jetInt_trk_dphi_%s_bkg_syst_%s_%s_%s", ptch.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar]->Write (Form ("g_jetInt_trk_dphi_%s_sig_syst_%s_%s_%s", ptch.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][iVar]->Write (Form ("g_jetInt_trk_dphi_%s_iaa_syst_%s_%s_%s", ptch.Data (), cent.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iPtCh

      } // end loop over iPtJInt

    } // end loop over iVar


    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      for (short iTotVar = 0; iTotVar < 3; iTotVar++) {

        const TString totVar = totalVariations[iTotVar];

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          g_jetInt_trk_pt_ref_systTot[iPtJInt][iDir][iTotVar]->Write     (Form ("g_jetInt_trk_pt_%s_ref_%s_systTot_%s",      dir.Data (), totVar.Data (), pTJInt.Data ()));
          g_jetInt_trk_pt_ref_bkg_systTot[iPtJInt][iDir][iTotVar]->Write (Form ("g_jetInt_trk_pt_%s_ref_bkg_%s_systTot_%s",  dir.Data (), totVar.Data (), pTJInt.Data ()));
          g_jetInt_trk_pt_ref_sig_systTot[iPtJInt][iDir][iTotVar]->Write (Form ("g_jetInt_trk_pt_%s_ref_sig_%s_systTot_%s",  dir.Data (), totVar.Data (), pTJInt.Data ()));
          g_jetInt_trk_pt_ref_unf_systTot[iPtJInt][iDir][iTotVar]->Write (Form ("g_jetInt_trk_pt_%s_ref_unf_%s_systTot_%s",  dir.Data (), totVar.Data (), pTJInt.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            g_jetInt_trk_pt_systTot[iPtJInt][iDir][iCent][iTotVar]->Write      (Form ("g_jetInt_trk_pt_%s_%s_systTot_%s_%s",     dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            g_jetInt_trk_pt_bkg_systTot[iPtJInt][iDir][iCent][iTotVar]->Write  (Form ("g_jetInt_trk_pt_%s_bkg_%s_systTot_%s_%s", dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            g_jetInt_trk_pt_sig_systTot[iPtJInt][iDir][iCent][iTotVar]->Write  (Form ("g_jetInt_trk_pt_%s_sig_%s_systTot_%s_%s", dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            g_jetInt_trk_pt_unf_systTot[iPtJInt][iDir][iCent][iTotVar]->Write  (Form ("g_jetInt_trk_pt_%s_unf_%s_systTot_%s_%s", dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            g_jetInt_trk_pt_iaa_systTot[iPtJInt][iDir][iCent][iTotVar]->Write  (Form ("g_jetInt_trk_pt_%s_iaa_%s_systTot_%s_%s", dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));

          } // end loop over iCent

        } // end loop over iDir

        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          const TString ptch = pTChSelections[iPtCh].Data ();

          g_jetInt_trk_dphi_ref_systTot[iPtJInt][iPtCh][iTotVar]->Write      (Form ("g_jetInt_trk_dphi_%s_ref_%s_systTot_%s",      ptch.Data (), totVar.Data (), pTJInt.Data ()));
          g_jetInt_trk_dphi_ref_bkg_systTot[iPtJInt][iPtCh][iTotVar]->Write  (Form ("g_jetInt_trk_dphi_%s_ref_bkg_%s_systTot_%s",  ptch.Data (), totVar.Data (), pTJInt.Data ()));
          g_jetInt_trk_dphi_ref_sig_systTot[iPtJInt][iPtCh][iTotVar]->Write  (Form ("g_jetInt_trk_dphi_%s_ref_sig_%s_systTot_%s",  ptch.Data (), totVar.Data (), pTJInt.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            g_jetInt_trk_dphi_systTot[iPtJInt][iPtCh][iCent][iTotVar]->Write     (Form ("g_jetInt_trk_dphi_%s_%s_systTot_%s_%s",      ptch.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            g_jetInt_trk_dphi_bkg_systTot[iPtJInt][iPtCh][iCent][iTotVar]->Write (Form ("g_jetInt_trk_dphi_%s_bkg_%s_systTot_%s_%s",  ptch.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            g_jetInt_trk_dphi_sig_systTot[iPtJInt][iPtCh][iCent][iTotVar]->Write (Form ("g_jetInt_trk_dphi_%s_sig_%s_systTot_%s_%s",  ptch.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            g_jetInt_trk_dphi_iaa_systTot[iPtJInt][iPtCh][iCent][iTotVar]->Write (Form ("g_jetInt_trk_dphi_%s_iaa_%s_systTot_%s_%s",  ptch.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));

          } // end loop over iCent

        } // end loop over iPtCh

      } // end loop over iTotVar

    } // end loop over iPtJInt



    for (short iDType = 0; iDType < 2; iDType++) {

      for (short iVar = 0; iVar < nVar; iVar++) {

        const TString var = variations[iVar];

        if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
          continue;

        h_evt_counts_ref[iDType][iVar]->Write ();

        for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          h_jet_counts_ref[iDType][iPtJ][iVar]->Write ();

        } // end loop over iPtJ

        h_jet_pt_ref[iDType][iVar]->Write ();
        h_jet_pt_ref_unf[iDType][iVar]->Write ();

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          h_evt_counts[iDType][iCent][iVar]->Write ();

          h_jet_pt[iDType][iCent][iVar]->Write ();
          h_jet_pt_unf[iDType][iCent][iVar]->Write ();

          for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

            h_jet_counts[iDType][iPtJ][iCent][iVar]->Write ();

          } // end loop over iPtJ

        } // end loop over iCent

      } // end loop over iVar


      for (short iPtJInt : {0, 1}) {

        for (short iVar = 0; iVar < nVar; iVar++) {

          const TString var = variations[iVar];

          if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
            continue;

          for (short iDir = 0; iDir < nDir; iDir++) {

            h_jetInt_trk_pt_ref[iDType][iPtJInt][iDir][iVar]->Write ();
            h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][iDir][iVar]->Write ();
            h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][iVar]->Write ();
            h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar]->Write ();

            for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

              h_jetInt_trk_pt[iDType][iPtJInt][iDir][iCent][iVar]->Write ();
              h_jetInt_trk_pt_bkg[iDType][iPtJInt][iDir][iCent][iVar]->Write ();
              h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][iVar]->Write ();
              h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar]->Write ();
              h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent][iVar]->Write ();
              h_jetInt_trk_pt_iaaNoUnf[iDType][iPtJInt][iDir][iCent][iVar]->Write ();

            } // end loop over iCent

          } // end loop over iDir

          for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

            h_jetInt_trk_dphi_ref[iDType][iPtJInt][iPtCh][iVar]->Write ();
            h_jetInt_trk_dphi_ref_bkg[iDType][iPtJInt][iPtCh][iVar]->Write ();
            h_jetInt_trk_dphi_ref_sig[iDType][iPtJInt][iPtCh][iVar]->Write ();

            for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

              h_jetInt_trk_dphi[iDType][iPtJInt][iPtCh][iCent][iVar]->Write ();
              h_jetInt_trk_dphi_bkg[iDType][iPtJInt][iPtCh][iCent][iVar]->Write ();
              h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent][iVar]->Write ();
              h_jetInt_trk_dphi_iaa[iDType][iPtJInt][iPtCh][iCent][iVar]->Write ();

            } // end loop over iCent

          } // end loop over iPtCh

        } // end loop over iVar

      } // end loop over iPtJInt

    } // end loop over iDType



    //for (short iPtJInt : {0, 1}) {

    //  for (short iDir = 0; iDir < nDir; iDir++) {

    //    //h_jetInt_trk_pt_ref_sig_bbb[iDir]->Write ();
    //    //f_jet_trk_pt_ref_sig_bbb[iDir]->Write ();

    //    h_jetInt_trk_pt_ref_unf_closure[iPtJInt][iDir]->Write ();
    //    f_jetInt_trk_pt_ref_unf_closure[iPtJInt][iDir]->Write ();

    //    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

    //      //h_jetInt_trk_pt_sig_bbb[iDir][iCent]->Write ();
    //      //f_jet_trk_pt_sig_bbb[iDir][iCent]->Write ();

    //      h_jetInt_trk_pt_unf_closure[iPtJInt][iDir][iCent]->Write ();
    //      f_jetInt_trk_pt_unf_closure[iPtJInt][iDir][iCent]->Write ();

    //      h_jetInt_trk_pt_iaa_closure[iPtJInt][iDir][iCent]->Write ();
    //      f_jetInt_trk_pt_iaa_closure[iPtJInt][iDir][iCent]->Write ();

    //    } // end loop over iCent

    //  } // end loop over iDir

    //} // end loop over iPtJInt



    g_jet_pt_ref_unfSumUnc->Write     ("g_jet_pt_ref_unfSumUnc");
    g_jet_pt_ref_unfIterUnc->Write    ("g_jet_pt_ref_unfIterUnc");
    g_jet_pt_ref_unfTotUnc->Write     ("g_jet_pt_ref_unfTotUnc");
    g_jet_pt_ref_unfSumRelUnc->Write  ("g_jet_pt_ref_unfSumRelUnc");
    g_jet_pt_ref_unfIterRelUnc->Write ("g_jet_pt_ref_unfIterRelUnc");
    g_jet_pt_ref_unfTotRelUnc->Write  ("g_jet_pt_ref_unfTotRelUnc");

    for (short iDir = 0; iDir < nDir; iDir++) {

      const TString dir = directions[iDir];

      for (short iPtJInt : {0, 1}) {

        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

        for (short iIter = 0; iIter < nItersMax-nItersMin+2; iIter++) {
  
          h_jet_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter]->Write ();

        } // end loop over iIter

        g_jetInt_trk_pt_ref_unfSumUnc[iPtJInt][iDir]->Write     (Form ("g_jetInt_trk_pt_ref_unfSumUnc_%s_%s",     dir.Data (), pTJInt.Data ()));
        g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir]->Write    (Form ("g_jetInt_trk_pt_ref_unfIterUnc_%s_%s",    dir.Data (), pTJInt.Data ()));
        g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir]->Write     (Form ("g_jetInt_trk_pt_ref_unfTotUnc_%s_%s",     dir.Data (), pTJInt.Data ()));
        g_jetInt_trk_pt_ref_unfSumRelUnc[iPtJInt][iDir]->Write  (Form ("g_jetInt_trk_pt_ref_unfSumRelUnc_%s_%s",  dir.Data (), pTJInt.Data ()));
        g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir]->Write (Form ("g_jetInt_trk_pt_ref_unfIterRelUnc_%s_%s", dir.Data (), pTJInt.Data ()));
        g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir]->Write  (Form ("g_jetInt_trk_pt_ref_unfTotRelUnc_%s_%s",  dir.Data (), pTJInt.Data ()));

      } // end loop over iPtJInt

    } // end loop over iDir

    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

      g_jet_pt_unfSumUnc[iCent]->Write      (Form ("g_jet_pt_pPb_%s_unfSumUnc",     cent.Data ()));
      g_jet_pt_unfIterUnc[iCent]->Write     (Form ("g_jet_pt_pPb_%s_unfIterUnc",    cent.Data ()));
      g_jet_pt_unfTotUnc[iCent]->Write      (Form ("g_jet_pt_pPb_%s_unfTotUnc",     cent.Data ()));
      g_jet_pt_unfSumRelUnc[iCent]->Write   (Form ("g_jet_pt_pPb_%s_unfSumRelUnc",  cent.Data ()));
      g_jet_pt_unfIterRelUnc[iCent]->Write  (Form ("g_jet_pt_pPb_%s_unfIterRelUnc", cent.Data ()));
      g_jet_pt_unfTotRelUnc[iCent]->Write   (Form ("g_jet_pt_pPb_%s_unfTotRelUnc",  cent.Data ()));

      for (short iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          for (short iIter = 0; iIter < nItersMax-nItersMin+2; iIter++) {
    
            h_jet_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter]->Write ();
  
          } // end loop over iIter

          g_jetInt_trk_pt_unfSumUnc[iPtJInt][iDir][iCent]->Write      (Form ("g_jetInt_trk_pt_%s_unfSumUnc_%s_%s",      dir.Data (), cent.Data (), pTJInt.Data ()));
          g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent]->Write     (Form ("g_jetInt_trk_pt_%s_unfIterUnc_%s_%s",     dir.Data (), cent.Data (), pTJInt.Data ()));
          g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent]->Write      (Form ("g_jetInt_trk_pt_%s_unfTotUnc_%s_%s",      dir.Data (), cent.Data (), pTJInt.Data ()));
          g_jetInt_trk_pt_unfSumRelUnc[iPtJInt][iDir][iCent]->Write   (Form ("g_jetInt_trk_pt_%s_unfSumRelUnc_%s_%s",   dir.Data (), cent.Data (), pTJInt.Data ()));
          g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent]->Write  (Form ("g_jetInt_trk_pt_%s_unfIterRelUnc_%s_%s",  dir.Data (), cent.Data (), pTJInt.Data ()));
          g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent]->Write   (Form ("g_jetInt_trk_pt_%s_unfTotRelUnc_%s_%s",   dir.Data (), cent.Data (), pTJInt.Data ()));

        } // end loop over iPtJInt

      } // end loop over iDir

    } // end loop over iCent


    outFile->Close ();
  }
}


#endif
