#ifndef __JetHadronCorrelator_TriggerStudy_C__
#define __JetHadronCorrelator_TriggerStudy_C__

#include "CentralityDefs.h"
#include "Params.h"
#include "OutTree.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Variations.h"

#include <ArrayTemplates.h>
#include <Utilities.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <string.h>
#include <math.h>

using namespace JetHadronCorrelations;


void TriggerStudy () {

  TString inFilePattern = Form ("%s/Trees/J50/data16_5TeV_iCent*.root", rootPath.Data ());
  //TString inFilePattern = Form ("%s/Trees/J50/data16_5TeV_iCent0.root", rootPath.Data ());
  std::cout << "Reading " << inFilePattern.Data () << std::endl;
  TChain* inTree = new TChain ("pPbTree", "pPbTree");
  inTree->Add (inFilePattern.Data ());


  TString outFileName = Form ("%s/TriggerStudy.root", rootPath.Data ());
  std::cout << "Writing " << outFileName.Data () << std::endl;
  TFile* outFile = new TFile (outFileName.Data (), "recreate");
  

  int jet_n = 0;
  float jet_pt[40] = {};

  inTree->SetBranchAddress ("mbTrig",             &mbTrig);
  inTree->SetBranchAddress ("mbTrigPS",           &mbTrigPS);
  inTree->SetBranchAddress ("jetTrig",            &jetTrig);
  inTree->SetBranchAddress ("jetTrigPS",          &jetTrigPS);
  inTree->SetBranchAddress ("zdc_calibE_Pb",      &zdc_calibE_Pb);
  inTree->SetBranchAddress ("leading_jet",        &leading_jet);
  inTree->SetBranchAddress ("akt4_hi_jet_n",      &jet_n);
  inTree->SetBranchAddress ("akt4_hi_jet_pt",     &jet_pt);


  TEfficiency* e_pPb_j50_eff = new TEfficiency ("e_pPb_j50_eff", ";#it{p}_{T}^{lead} [GeV];Efficiency", 80, 20, 100);


  const long nEvts = inTree->GetEntries ();
  for (long iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      std::cout << "Events " << iEvt / (nEvts / 100) << "\% done...\r" << flush;

    inTree->GetEntry (iEvt);

    if (!mbTrig)
      continue;

    if (leading_jet < 0 || jet_pt[leading_jet] < 20)
      continue;

    e_pPb_j50_eff->Fill (jetTrig, jet_pt[leading_jet]);

  }

  outFile->cd ();

  e_pPb_j50_eff->Write ();

  outFile->Close ();

}

#endif
