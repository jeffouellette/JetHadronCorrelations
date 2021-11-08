#ifndef __JetHadronCorrelator_JetEtaPhiCheck_C__
#define __JetHadronCorrelator_JetEtaPhiCheck_C__

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


void JetEtaPhiCheck () {

  const TString inFileName = Form ("%s/Trees/J50/data17_5TeV.root", rootPath.Data ());
  std::cout << "Reading " << inFileName.Data () << std::endl;
  TFile* inFile = new TFile (inFileName.Data (), "read");
  TTree* ppTree = (TTree*) inFile->Get ("ppTree");
  

  int jet_n = 0;
  float jet_pt[40] = {};
  float jet_eta[40] = {};
  float jet_phi[40] = {};

  ppTree->SetBranchAddress ("leading_jet",        &leading_jet);
  ppTree->SetBranchAddress ("akt4_hi_jet_n",      &jet_n);
  ppTree->SetBranchAddress ("akt4_hi_jet_pt",     &jet_pt);
  ppTree->SetBranchAddress ("akt4_hi_jet_eta",    &jet_eta);
  ppTree->SetBranchAddress ("akt4_hi_jet_phi",    &jet_phi);


  TH1D* h_jet_phi = new TH1D ("h_jet_phi", ";#phi^{jet};Counts", 64, -M_PI, M_PI);
  h_jet_phi->Sumw2 ();

  double nJetOutsideHEC = 0;
  double nJetInHEC = 0;

  const long nEvts = ppTree->GetEntries ();
  for (long iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      std::cout << "Events " << iEvt / (nEvts / 100) << "\% done...\r" << flush;

    ppTree->GetEntry (iEvt);

    for (short iJ = 0; iJ < jet_n; iJ++) {

      if (jet_pt[iJ] > 60 && 1.1 < jet_eta[iJ] && jet_eta[iJ] < 2.8) {
        h_jet_phi->Fill (jet_phi[iJ]);
        if ((-M_PI <= jet_phi[iJ] && jet_phi[iJ] < -0.5*M_PI+0.4) || (M_PI-0.4 < jet_phi[iJ] && jet_phi[iJ] <= M_PI))
          nJetInHEC += 1.;
        else
          nJetOutsideHEC += 1.;
      }


    } // end loop over iJ

  } // end loop over iEvt


  const double fracIn = nJetInHEC / (nJetInHEC + nJetOutsideHEC);
  const double fracInErr = std::sqrt (fracIn * (1-fracIn) / (nJetInHEC + nJetOutsideHEC));
  const double fracOut = nJetOutsideHEC / (nJetInHEC + nJetOutsideHEC);
  const double fracOutErr = std::sqrt (fracOut * (1-fracOut) / (nJetInHEC + nJetOutsideHEC));


  std::cout << "Jets inside HEC: " << nJetInHEC << " (" << 100 * fracIn << " +/- " << 100 * fracInErr << "%)" << std::endl;
  std::cout << "Jets outside HEC: " << nJetOutsideHEC << " (" << 100 * fracOut << " +/- " << 100 * fracOutErr << "%)" << std::endl;


  TString outFileName = Form ("%s/JetEtaPhiCheck.root", rootPath.Data ());
  TFile* outFile = new TFile (outFileName.Data (), "recreate");

  h_jet_phi->Write ();

  outFile->Close ();

}

#endif
