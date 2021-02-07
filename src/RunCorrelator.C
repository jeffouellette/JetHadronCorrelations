#ifndef __RunCorrelator_C__
#define __RunCorrelator_C__

#include "Params.h"
#include "OutTree.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"

#include <Utilities.h>

#include <TChain.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLorentzVector.h>

#include <iostream>
#include <math.h>

using namespace HadronYieldsAnalysis;



void Correlator (const char* tag, const char* outFileName, TTree* jetsTree, TTree* tracksTree = nullptr) {

  bool doMixing = true;
  if (tracksTree == nullptr) {
    tracksTree = jetsTree;
    doMixing = false;
  }

  SetupDirectories ("Data");

  // get event tagging & jet information from main tree
  jetsTree->Branch ("event_number",  &event_number);
  jetsTree->Branch ("lumi_block",    &lumi_block);
  jetsTree->Branch ("run_number",    &run_number);

  jetsTree->Branch ("event_weight",  &event_weight); 

  jetsTree->SetBranchAddress ("fcal_et_Pb",     &fcal_et_Pb);
  jetsTree->SetBranchAddress ("fcal_et_p",      &fcal_et_p);
  jetsTree->SetBranchAddress ("zdc_calibE_Pb",  &zdc_calibE_Pb);
  jetsTree->SetBranchAddress ("zdc_calibE_p",   &zdc_calibE_p);


  jetsTree->SetBranchAddress ("leading_jet",                &leading_jet);
  jetsTree->SetBranchAddress ("subleading_jet",             &subleading_jet);
  jetsTree->SetBranchAddress ("akt4_hi_jet_n",              &akt4_hi_jet_n);
  jetsTree->SetBranchAddress ("akt4_hi_jet_pt",             &akt4_hi_jet_pt);
  jetsTree->SetBranchAddress ("akt4_hi_jet_eta",            &akt4_hi_jet_eta);
  jetsTree->SetBranchAddress ("akt4_hi_jet_phi",            &akt4_hi_jet_phi);
  jetsTree->SetBranchAddress ("akt4_hi_jet_e",              &akt4_hi_jet_e);
  jetsTree->SetBranchAddress ("leading_jet_phi_transmin",   &leading_jet_phi_transmin);
  jetsTree->SetBranchAddress ("leading_jet_phi_transmax",   &leading_jet_phi_transmax);


  // get event matching & track information from the mixed event tree
  float fcal_et_Pb_matching;
  float fcal_et_p_matching;
  float zdc_calibE_Pb_matching;
  float zdc_calibE_p_matching;
  if (doMixing) {
    tracksTree->SetBranchAddress ("fcal_et_Pb",     &fcal_et_Pb_matching);
    tracksTree->SetBranchAddress ("fcal_et_p",      &fcal_et_p_matching);
    tracksTree->SetBranchAddress ("zdc_calibE_Pb",  &zdc_calibE_Pb_matching);
    tracksTree->SetBranchAddress ("zdc_calibE_p",   &zdc_calibE_p_matching);
  }


  tracksTree->SetBranchAddress ("trk_n",          &trk_n);
  tracksTree->SetBranchAddress ("trk_pt",         &trk_pt);
  tracksTree->SetBranchAddress ("trk_eta",        &trk_eta);
  tracksTree->SetBranchAddress ("trk_phi",        &trk_phi);



  TFile* outFile = new TFile (outFileName, "recreate");


  TH1D* h_jet_counts = new TH1D (Form ("h_jet_counts_%s", tag), "", 3, 0.5, 3.5);
  TH1D* h_ljet_counts = new TH1D (Form ("h_ljet_counts_%s", tag), "", 3, 0.5, 3.5);
  TH1D* h_sljet_counts = new TH1D (Form ("h_sljet_counts_%s", tag), "", 3, 0.5, 3.5);

  TH1D* h_jet_pt = new TH1D (Form ("h_jet_pt_%s", tag), ";#it{p}_{T}^{jet} [GeV];(1/N_{jet}) (dN_{jet}/d#it{p}_{T}) [GeV^{-1}]", nPtJBins, pTJBins);
  TH2D* h2_jet_pt_cov = new TH2D (Form ("h2_jet_pt_cov_%s", tag), ";#it{p}_{T}^{jet} [GeV];#it{p}_{T}^{jet} [GeV];Covariance", nPtJBins, pTJBins, nPtJBins, pTJBins);

  TH2D* h2_jet_eta_phi = new TH2D (Form ("h2_jet_eta_phi_%s", tag), ";#eta^{jet};#phi^{jet};Counts", 40, -3.2, 3.2, 40, -pi, pi);

  TH1D* h_jet_trk_dphi = new TH1D (Form ("h_jet_trk_dphi_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{trk}/d#Delta#phi)", nDPhiBins, 0, pi);
  TH2D* h2_jet_trk_dphi_cov = new TH2D (Form ("h2_jet_trk_dphi_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, 0, pi, nDPhiBins, 0, pi);
  TH1D* h_ljet_trk_dphi = new TH1D (Form ("h_ljet_trk_dphi_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{trk}/d#Delta#phi)", nDPhiBins, 0, pi);
  TH2D* h2_ljet_trk_dphi_cov = new TH2D (Form ("h2_ljet_trk_dphi_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, 0, pi, nDPhiBins, 0, pi);
  TH1D* h_sljet_trk_dphi = new TH1D (Form ("h_sljet_trk_dphi_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{trk}/d#Delta#phi)", nDPhiBins, 0, pi);
  TH2D* h2_sljet_trk_dphi_cov = new TH2D (Form ("h2_sljet_trk_dphi_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, 0, pi, nDPhiBins, 0, pi);

  TH1D* h_jet_trk_pt_ns = new TH1D (Form ("h_jet_trk_pt_ns_%s", tag), ";#it{p}_{T}^{ch} [GeV];(1/N_{jet}) (dN_{trk}/d#it{p}_{T})", nPtChBins, pTChBins);
  TH2D* h2_jet_trk_pt_ns_cov = new TH2D (Form ("h2_jet_trk_pt_ns_cov_%s", tag), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{ch} [GeV];Covariance", nPtChBins, pTChBins, nPtChBins, pTChBins);
  TH1D* h_ljet_trk_pt_ns = new TH1D (Form ("h_ljet_trk_pt_ns_%s", tag), ";#it{p}_{T}^{ch} [GeV];(1/N_{jet}) (dN_{trk}/d#it{p}_{T})", nPtChBins, pTChBins);
  TH2D* h2_ljet_trk_pt_ns_cov = new TH2D (Form ("h2_ljet_trk_pt_ns_cov_%s", tag), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{ch} [GeV];Covariance", nPtChBins, pTChBins, nPtChBins, pTChBins);
  TH1D* h_sljet_trk_pt_ns = new TH1D (Form ("h_sljet_trk_pt_ns_%s", tag), ";#it{p}_{T}^{ch} [GeV];(1/N_{jet}) (dN_{trk}/d#it{p}_{T})", nPtChBins, pTChBins);
  TH2D* h2_sljet_trk_pt_ns_cov = new TH2D (Form ("h2_sljet_trk_pt_ns_cov_%s", tag), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{ch} [GeV];Covariance", nPtChBins, pTChBins, nPtChBins, pTChBins);

  TH1D* h_jet_trk_pt_as = new TH1D (Form ("h_jet_trk_pt_as_%s", tag), ";#it{p}_{T}^{ch} [GeV];(1/N_{jet}) (dN_{trk}/d#it{p}_{T})", nPtChBins, pTChBins);
  TH2D* h2_jet_trk_pt_as_cov = new TH2D (Form ("h2_jet_trk_pt_as_cov_%s", tag), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{ch} [GeV];Covariance", nPtChBins, pTChBins, nPtChBins, pTChBins);
  TH1D* h_ljet_trk_pt_as = new TH1D (Form ("h_ljet_trk_pt_as_%s", tag), ";#it{p}_{T}^{ch} [GeV];(1/N_{jet}) (dN_{trk}/d#it{p}_{T})", nPtChBins, pTChBins);
  TH2D* h2_ljet_trk_pt_as_cov = new TH2D (Form ("h2_ljet_trk_pt_as_cov_%s", tag), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{ch} [GeV];Covariance", nPtChBins, pTChBins, nPtChBins, pTChBins);
  TH1D* h_sljet_trk_pt_as = new TH1D (Form ("h_sljet_trk_pt_as_%s", tag), ";#it{p}_{T}^{ch} [GeV];(1/N_{jet}) (dN_{trk}/d#it{p}_{T})", nPtChBins, pTChBins);
  TH2D* h2_sljet_trk_pt_as_cov = new TH2D (Form ("h2_sljet_trk_pt_as_cov_%s", tag), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{ch} [GeV];Covariance", nPtChBins, pTChBins, nPtChBins, pTChBins);



  // arrays for filling histograms & covariance matrices correctly
  double jet_pt_counts[nPtJBins];
  double jet_trk_dphi_counts[nDPhiBins];
  double ljet_trk_dphi_counts[nDPhiBins];
  double sljet_trk_dphi_counts[nDPhiBins];

  double jet_trk_pt_ns_counts[nPtChBins];
  double ljet_trk_pt_ns_counts[nPtChBins];
  double sljet_trk_pt_ns_counts[nPtChBins];
  double jet_trk_pt_as_counts[nPtChBins];
  double ljet_trk_pt_as_counts[nPtChBins];
  double sljet_trk_pt_as_counts[nPtChBins];

  const int nEvts = (doMixing ? 10. : 1.) * (jetsTree->GetEntries ());
  const int nTrkEvts = tracksTree->GetEntries ();


  int iTrkEvt = 0;
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)                                
      cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;  

    for (int iX = 0; iX < nPtJBins; iX++)
      jet_pt_counts[iX] = 0;
    for (int iX = 0; iX < nDPhiBins; iX++) {
      jet_trk_dphi_counts[iX] = 0;
      ljet_trk_dphi_counts[iX] = 0;
      sljet_trk_dphi_counts[iX] = 0;
    }
    for (int iX = 0; iX < nPtChBins; iX++) {
      jet_trk_pt_ns_counts[iX] = 0;
      ljet_trk_pt_ns_counts[iX] = 0;
      sljet_trk_pt_ns_counts[iX] = 0;
      jet_trk_pt_as_counts[iX] = 0;
      ljet_trk_pt_as_counts[iX] = 0;
      sljet_trk_pt_as_counts[iX] = 0;
    }

    jetsTree->GetEntry (iEvt % jetsTree->GetEntries ());

    const double yboost = GetBoost (run_number);

    const short iCent = (GetZdcCentBin (zdc_calibE_Pb));
    if (IspPb () && iCent != numZdcCentBins-1 && iCent != numZdcCentBins-2) // only look at 0-10% or 10-20% central events
      continue;

    h_jet_counts->Fill (1, akt4_hi_jet_n); // adds to number of jets (i.e. denominator)
    if (leading_jet != -1)
      h_ljet_counts->Fill (1);
    if (subleading_jet != -1)
      h_sljet_counts->Fill (1);


    // do a mixed event
    if (doMixing) {
      const int oldTrkEvt = iTrkEvt;
      bool goodEvent = false;
      do {
        iTrkEvt = (iTrkEvt + 1) % nTrkEvts;
        tracksTree->GetEntry (iTrkEvt);
        goodEvent = (GetZdcCentBin (zdc_calibE_Pb_matching) == iCent); // TODO -- mixing categories!
      }
      while (!goodEvent && iTrkEvt != oldTrkEvt);
    }


    for (int iJet = 0; iJet < akt4_hi_jet_n; iJet++) {
      h2_jet_eta_phi->Fill (akt4_hi_jet_eta[iJet], akt4_hi_jet_phi[iJet]);

      const short iPtJ = GetPtJBin (akt4_hi_jet_pt[iJet]);
      if (0 <= iPtJ && iPtJ < nPtJBins)
        jet_pt_counts[iPtJ]++;
    }

    for (int iTrk = 0; iTrk < trk_n; iTrk++) {

      const short iPtCh = GetPtChBin (trk_pt[iTrk]);
      if (iPtCh < 0 || iPtCh >= nPtChBins)
        continue;

      if (fabs (trk_eta[iTrk] - yboost) > 2.5 - 0.465)
        continue;

      for (int iJet = 0; iJet < akt4_hi_jet_n; iJet++) {
        const float dphi = DeltaPhi (akt4_hi_jet_phi[iJet], trk_phi[iTrk]);
        const short iDPhi = GetDPhiBin (dphi);
        if (iDPhi < 0 || nDPhiBins < iDPhi)
          continue;

        if (trk_pt[iTrk] > 4)
          jet_trk_dphi_counts[iDPhi]++;
        if (dphi < pi/8.)
          jet_trk_pt_ns_counts[iPtCh]++;
        else if (dphi > 7.*pi/8.)
          jet_trk_pt_as_counts[iPtCh]++;
      }

      if (leading_jet != -1) {
        const float dphi = DeltaPhi (akt4_hi_jet_phi[leading_jet], trk_phi[iTrk]);
        const short iDPhi = GetDPhiBin (dphi);
        if (iDPhi < 0 || nDPhiBins < iDPhi)
          continue;

        if (trk_pt[iTrk] > 4)
          ljet_trk_dphi_counts[iDPhi]++;
        if (dphi < pi/8.)
          ljet_trk_pt_ns_counts[iPtCh]++;
        else if (dphi > 7.*pi/8.)
          ljet_trk_pt_as_counts[iPtCh]++;
      }

      if (subleading_jet != -1) {
        const float dphi = DeltaPhi (akt4_hi_jet_phi[subleading_jet], trk_phi[iTrk]);
        const short iDPhi = GetDPhiBin (dphi);
        if (iDPhi < 0 || nDPhiBins < iDPhi)
          continue;

        if (trk_pt[iTrk] > 4)
          sljet_trk_dphi_counts[iDPhi]++;
        if (dphi < pi/8.)
          sljet_trk_pt_ns_counts[iPtCh]++;
        else if (dphi > 7.*pi/8.)
          sljet_trk_pt_as_counts[iPtCh]++;
      }
    }

    for (int iX = 0; iX < nPtJBins; iX++) {
      h_jet_pt->SetBinContent (iX+1, h_jet_pt->GetBinContent (iX+1) + jet_pt_counts[iX]);
      for (int iY = 0; iY < nPtJBins; iY++)
        h2_jet_pt_cov->SetBinContent (iX+1, iY+1, h2_jet_pt_cov->GetBinContent (iX+1, iY+1) + (jet_pt_counts[iX])*(jet_pt_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi->SetBinContent (iX+1, h_jet_trk_dphi->GetBinContent (iX+1) + jet_trk_dphi_counts[iX]);
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_cov->GetBinContent (iX+1, iY+1) + (jet_trk_dphi_counts[iX])*(jet_trk_dphi_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_ljet_trk_dphi->SetBinContent (iX+1, h_ljet_trk_dphi->GetBinContent (iX+1) + ljet_trk_dphi_counts[iX]);
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_ljet_trk_dphi_cov->SetBinContent (iX+1, iY+1, h2_ljet_trk_dphi_cov->GetBinContent (iX+1, iY+1) + (ljet_trk_dphi_counts[iX])*(ljet_trk_dphi_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_sljet_trk_dphi->SetBinContent (iX+1, h_sljet_trk_dphi->GetBinContent (iX+1) + sljet_trk_dphi_counts[iX]);
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_sljet_trk_dphi_cov->SetBinContent (iX+1, iY+1, h2_sljet_trk_dphi_cov->GetBinContent (iX+1, iY+1) + (sljet_trk_dphi_counts[iX])*(sljet_trk_dphi_counts[iY]));
    }

    for (int iX = 0; iX < nPtChBins; iX++) {
      h_jet_trk_pt_ns->SetBinContent (iX+1, h_jet_trk_pt_ns->GetBinContent (iX+1) + jet_trk_pt_ns_counts[iX]);
      for (int iY = 0; iY < nPtChBins; iY++)
        h2_jet_trk_pt_ns_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_pt_ns_cov->GetBinContent (iX+1, iY+1) + (jet_trk_pt_ns_counts[iX])*(jet_trk_pt_ns_counts[iY]));
    }
    for (int iX = 0; iX < nPtChBins; iX++) {
      h_ljet_trk_pt_ns->SetBinContent (iX+1, h_ljet_trk_pt_ns->GetBinContent (iX+1) + ljet_trk_pt_ns_counts[iX]);
      for (int iY = 0; iY < nPtChBins; iY++)
        h2_ljet_trk_pt_ns_cov->SetBinContent (iX+1, iY+1, h2_ljet_trk_pt_ns_cov->GetBinContent (iX+1, iY+1) + (ljet_trk_pt_ns_counts[iX])*(ljet_trk_pt_ns_counts[iY]));
    }
    for (int iX = 0; iX < nPtChBins; iX++) {
      h_sljet_trk_pt_ns->SetBinContent (iX+1, h_sljet_trk_pt_ns->GetBinContent (iX+1) + sljet_trk_pt_ns_counts[iX]);
      for (int iY = 0; iY < nPtChBins; iY++)
        h2_sljet_trk_pt_ns_cov->SetBinContent (iX+1, iY+1, h2_sljet_trk_pt_ns_cov->GetBinContent (iX+1, iY+1) + (sljet_trk_pt_ns_counts[iX])*(sljet_trk_pt_ns_counts[iY]));
    }

    for (int iX = 0; iX < nPtChBins; iX++) {
      h_jet_trk_pt_as->SetBinContent (iX+1, h_jet_trk_pt_as->GetBinContent (iX+1) + jet_trk_pt_as_counts[iX]);
      for (int iY = 0; iY < nPtChBins; iY++)
        h2_jet_trk_pt_as_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_pt_as_cov->GetBinContent (iX+1, iY+1) + (jet_trk_pt_as_counts[iX])*(jet_trk_pt_as_counts[iY]));
    }
    for (int iX = 0; iX < nPtChBins; iX++) {
      h_ljet_trk_pt_as->SetBinContent (iX+1, h_ljet_trk_pt_as->GetBinContent (iX+1) + ljet_trk_pt_as_counts[iX]);
      for (int iY = 0; iY < nPtChBins; iY++)
        h2_ljet_trk_pt_as_cov->SetBinContent (iX+1, iY+1, h2_ljet_trk_pt_as_cov->GetBinContent (iX+1, iY+1) + (ljet_trk_pt_as_counts[iX])*(ljet_trk_pt_as_counts[iY]));
    }
    for (int iX = 0; iX < nPtChBins; iX++) {
      h_sljet_trk_pt_as->SetBinContent (iX+1, h_sljet_trk_pt_as->GetBinContent (iX+1) + sljet_trk_pt_as_counts[iX]);
      for (int iY = 0; iY < nPtChBins; iY++)
        h2_sljet_trk_pt_as_cov->SetBinContent (iX+1, iY+1, h2_sljet_trk_pt_as_cov->GetBinContent (iX+1, iY+1) + (sljet_trk_pt_as_counts[iX])*(sljet_trk_pt_as_counts[iY]));
    }
  }
  cout << "Finished event loop." << endl;


  h_jet_counts->Write ();
  h_ljet_counts->Write ();
  h_sljet_counts->Write ();
  h_jet_pt->Write ();
  h2_jet_pt_cov->Write ();
  h2_jet_eta_phi->Write ();
  h_jet_trk_dphi->Write ();
  h2_jet_trk_dphi_cov->Write ();
  h_ljet_trk_dphi->Write ();
  h2_ljet_trk_dphi_cov->Write ();
  h_sljet_trk_dphi->Write ();
  h2_sljet_trk_dphi_cov->Write ();
  h_jet_trk_pt_ns->Write ();
  h2_jet_trk_pt_ns_cov->Write ();
  h_ljet_trk_pt_ns->Write ();
  h2_ljet_trk_pt_ns_cov->Write ();
  h_sljet_trk_pt_ns->Write ();
  h2_sljet_trk_pt_ns_cov->Write ();
  h_jet_trk_pt_as->Write ();
  h2_jet_trk_pt_as_cov->Write ();
  h_ljet_trk_pt_as->Write ();
  h2_ljet_trk_pt_as_cov->Write ();
  h_sljet_trk_pt_as->Write ();
  h2_sljet_trk_pt_as_cov->Write ();


  outFile->Close ();
}




void RunCorrelator (const char* collSys, const char* tag, const char* outFileName, const char* jetsInFileName, const char* tracksInFileName = nullptr) {

  if      (TString (collSys) == ToTString (CollisionSystem::pPb16s5TeV))  { collisionSystem = CollisionSystem::pPb16s5TeV; }
  else if (TString (collSys) == ToTString (CollisionSystem::pp17))        { collisionSystem = CollisionSystem::pp17;       }
  else {
    std::cout << "In RunCorrelator.C: Invalid or unsupported collision system, exiting." << std::endl;
    return;
  }

  TFile* jetsInFile = new TFile (jetsInFileName, "read");
  TFile* tracksInFile = (tracksInFileName != nullptr ? new TFile (tracksInFileName, "read") : nullptr);

  TTree* jetsInTree = (TTree*) jetsInFile->Get (IspPb () ? "pPbTree" : "ppTree");
  TTree* tracksInTree = (tracksInFileName != nullptr ? (TTree*) tracksInFile->Get ("pPbTree") : nullptr);

  Correlator (tag, outFileName, jetsInTree, tracksInTree);

  jetsInFile->Close ();
  if (tracksInFile)
    tracksInFile->Close ();

  return;
}

#endif
