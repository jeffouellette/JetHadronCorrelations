#ifndef __JetHadronCorrelator_PeripheralStudy_C__
#define __JetHadronCorrelator_PeripheralStudy_C__

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


void PeripheralStudy () {

  //TString inFileName = Form ("%s/Trees/MinBias/Nominal/data16_5TeV_iCent0.root", rootPath.Data ());
  //std::cout << "Reading " << inFileName.Data () << std::endl;
  //TFile* inFile = new TFile (inFileName.Data (), "read");
  //TTree* inTree = (TTree*) inFile->Get ("pPbTree");

  TString inFilePattern = Form ("%s/Trees/MinBias/Nominal/data16_5TeV_iCent*.root", rootPath.Data ());
  std::cout << "Reading " << inFilePattern.Data () << std::endl;
  TChain* inTree = new TChain ("pPbTree", "pPbTree");
  inTree->Add (Form ("%s/Trees/MinBias/Nominal/data16_5TeV_iCent*.root", rootPath.Data ()));
  

  int jet_n = 0;
  float jet_pt[100] = {};

  inTree->SetBranchAddress ("zdc_calibE_Pb",      &zdc_calibE_Pb);
  inTree->SetBranchAddress ("cluster_sumGap_Pb",  &cluster_sumGap_Pb);
  inTree->SetBranchAddress ("sumGap_Pb",          &sumGap_Pb);
  inTree->SetBranchAddress ("leading_jet",        &leading_jet);
  inTree->SetBranchAddress ("akt4_hi_jet_n",      &jet_n);
  inTree->SetBranchAddress ("akt4_hi_jet_pt",     &jet_pt);


  TH1D* h_sumGaps_Pb_0n0n = new TH1D ("h_sumGaps_Pb_0n0n", "", 50, -0.05, 4.95);
  h_sumGaps_Pb_0n0n->Sumw2 ();
  TH1D* h_sumGaps_Pb_0nXn = new TH1D ("h_sumGaps_Pb_0nXn", "", 50, -0.05, 4.95);
  h_sumGaps_Pb_0nXn->Sumw2 ();

  TH1D* h_cluster_sumGaps_Pb_0n0n = new TH1D ("h_cluster_sumGaps_Pb_0n0n", "", 50, -0.05, 4.95);
  h_cluster_sumGaps_Pb_0n0n->Sumw2 ();
  TH1D* h_cluster_sumGaps_Pb_0nXn = new TH1D ("h_cluster_sumGaps_Pb_0nXn", "", 50, -0.05, 4.95);
  h_cluster_sumGaps_Pb_0nXn->Sumw2 ();

  TH1D* h_sumGaps_Pb_0n0n_30GeVJet = new TH1D ("h_sumGaps_Pb_0n0n_30GeVJet", "", 50, -0.05, 4.95);
  h_sumGaps_Pb_0n0n_30GeVJet->Sumw2 ();
  TH1D* h_sumGaps_Pb_0nXn_30GeVJet = new TH1D ("h_sumGaps_Pb_0nXn_30GeVJet", "", 50, -0.05, 4.95);
  h_sumGaps_Pb_0nXn_30GeVJet->Sumw2 ();

  TH1D* h_cluster_sumGaps_Pb_0n0n_30GeVJet = new TH1D ("h_cluster_sumGaps_Pb_0n0n_30GeVJet", "", 50, -0.05, 4.95);
  h_cluster_sumGaps_Pb_0n0n_30GeVJet->Sumw2 ();
  TH1D* h_cluster_sumGaps_Pb_0nXn_30GeVJet = new TH1D ("h_cluster_sumGaps_Pb_0nXn_30GeVJet", "", 50, -0.05, 4.95);
  h_cluster_sumGaps_Pb_0nXn_30GeVJet->Sumw2 ();


  const long nEvts = inTree->GetEntries ();
  for (long iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      std::cout << "Events " << iEvt / (nEvts / 100) << "\% done...\r" << flush;

    inTree->GetEntry (iEvt);

    if (zdc_calibE_Pb > 18.0971) 
      continue;

    bool is0nXn = (zdc_calibE_Pb > 0);

    bool has30GeVJet = (leading_jet >= 0 && jet_pt[leading_jet] >= 30);

    if (is0nXn) {
      h_sumGaps_Pb_0nXn->Fill (sumGap_Pb);
      h_cluster_sumGaps_Pb_0nXn->Fill (cluster_sumGap_Pb);
      if (has30GeVJet) {
        h_sumGaps_Pb_0nXn_30GeVJet->Fill (sumGap_Pb);
        h_cluster_sumGaps_Pb_0nXn_30GeVJet->Fill (cluster_sumGap_Pb);
      }
    }
    else {
      h_sumGaps_Pb_0n0n->Fill (sumGap_Pb);
      h_cluster_sumGaps_Pb_0n0n->Fill (cluster_sumGap_Pb);
      if (has30GeVJet) {
        h_sumGaps_Pb_0n0n_30GeVJet->Fill (sumGap_Pb);
        h_cluster_sumGaps_Pb_0n0n_30GeVJet->Fill (cluster_sumGap_Pb);
      }
    }

  }


  TString outFileName = Form ("%s/PeripheralStudy/MinBias.root", rootPath.Data ());
  TFile* outFile = new TFile (outFileName.Data (), "recreate");

  h_sumGaps_Pb_0nXn->Write ();
  h_sumGaps_Pb_0n0n->Write ();
  h_cluster_sumGaps_Pb_0nXn->Write ();
  h_cluster_sumGaps_Pb_0n0n->Write ();
  h_sumGaps_Pb_0nXn_30GeVJet->Write ();
  h_sumGaps_Pb_0n0n_30GeVJet->Write ();
  h_cluster_sumGaps_Pb_0nXn_30GeVJet->Write ();
  h_cluster_sumGaps_Pb_0n0n_30GeVJet->Write ();

  outFile->Close ();

}

#endif
