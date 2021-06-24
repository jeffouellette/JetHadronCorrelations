#ifndef __RunCorrelator_C__
#define __RunCorrelator_C__

#include "Params.h"
#include "CentralityDefs.h"
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

using namespace JetHadronCorrelations;

bool doMixing = false;
bool doMixVar1 = false;
bool doMixVar2 = false;
bool doMixVar3 = false;


void Correlator (const char* tag, const char* outFileName, TTree* jetsTree, TTree* tracksTree = nullptr) {

  // get event tagging & jet information from main tree
  jetsTree->SetBranchAddress ("event_number",  &event_number);
  jetsTree->SetBranchAddress ("lumi_block",    &lumi_block);
  jetsTree->SetBranchAddress ("run_number",    &run_number);

  jetsTree->SetBranchAddress ("event_weight",  &event_weight); 

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
  //jetsTree->SetBranchAddress ("leading_jet_phi_transmin",   &leading_jet_phi_transmin);
  //jetsTree->SetBranchAddress ("leading_jet_phi_transmax",   &leading_jet_phi_transmax);


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

  TH2D* h2_trk_eff = LoadTrackingEfficiency ();
  TH2D* h2_trk_pur = LoadTrackingPurity ();


  TFile* outFile = new TFile (outFileName, "recreate");


  TH1D* h_evt_counts = new TH1D (Form ("h_evt_counts_%s", tag), "", 3, -0.5, 2.5);
  TH1D* h_jet_counts = new TH1D (Form ("h_jet_counts_%s", tag), "", 3, -0.5, 2.5);

  TH1D* h_jet_pt = new TH1D (Form ("h_jet_pt_%s", tag), ";#it{p}_{T}^{jet} [GeV];(1/N_{jet}) (dN_{jet}/d#it{p}_{T}) [GeV^{-1}]", nPtJBins, pTJBins);
  TH2D* h2_jet_pt_cov = new TH2D (Form ("h2_jet_pt_cov_%s", tag), ";#it{p}_{T}^{jet} [GeV];#it{p}_{T}^{jet} [GeV];Covariance", nPtJBins, pTJBins, nPtJBins, pTJBins);

  TH2D* h2_jet_eta_phi = new TH2D (Form ("h2_jet_eta_phi_%s", tag), ";#eta^{jet};#phi^{jet};Counts", 40, -3.2, 3.2, 40, -M_PI, M_PI);

  TH1D* h_jet_trk_dphi_gt0p5_lt1 = new TH1D (Form ("h_jet_trk_dphi_gt0p5_lt1_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
  TH2D* h2_jet_trk_dphi_gt0p5_lt1_cov = new TH2D (Form ("h2_jet_trk_dphi_gt0p5_lt1_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
  TH1D* h_jet_trk_dphi_gt1_lt1p5 = new TH1D (Form ("h_jet_trk_dphi_gt1_lt1p5_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
  TH2D* h2_jet_trk_dphi_gt1_lt1p5_cov = new TH2D (Form ("h2_jet_trk_dphi_gt1_lt1p5_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
  TH1D* h_jet_trk_dphi_gt1p5_lt2 = new TH1D (Form ("h_jet_trk_dphi_gt1p5_lt2_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
  TH2D* h2_jet_trk_dphi_gt1p5_lt2_cov = new TH2D (Form ("h2_jet_trk_dphi_gt1p5_lt2_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
  TH1D* h_jet_trk_dphi_gt2_lt4 = new TH1D (Form ("h_jet_trk_dphi_gt2_lt4_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
  TH2D* h2_jet_trk_dphi_gt2_lt4_cov = new TH2D (Form ("h2_jet_trk_dphi_gt2_lt4_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
  TH1D* h_jet_trk_dphi_gt4_lt6 = new TH1D (Form ("h_jet_trk_dphi_gt4_lt6_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
  TH2D* h2_jet_trk_dphi_gt4_lt6_cov = new TH2D (Form ("h2_jet_trk_dphi_gt4_lt6_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
  TH1D* h_jet_trk_dphi_gt6_lt8 = new TH1D (Form ("h_jet_trk_dphi_gt6_lt8_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
  TH2D* h2_jet_trk_dphi_gt6_lt8_cov = new TH2D (Form ("h2_jet_trk_dphi_gt6_lt8_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
  TH1D* h_jet_trk_dphi_gt8_lt10 = new TH1D (Form ("h_jet_trk_dphi_gt8_lt10_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
  TH2D* h2_jet_trk_dphi_gt8_lt10_cov = new TH2D (Form ("h2_jet_trk_dphi_gt8_lt10_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
  TH1D* h_jet_trk_dphi_gt10_lt15 = new TH1D (Form ("h_jet_trk_dphi_gt10_lt15_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
  TH2D* h2_jet_trk_dphi_gt10_lt15_cov = new TH2D (Form ("h2_jet_trk_dphi_gt10_lt15_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
  TH1D* h_jet_trk_dphi_gt15_lt20 = new TH1D (Form ("h_jet_trk_dphi_gt15_lt20_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
  TH2D* h2_jet_trk_dphi_gt15_lt20_cov = new TH2D (Form ("h2_jet_trk_dphi_gt15_lt20_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
  TH1D* h_jet_trk_dphi_gt20_lt30 = new TH1D (Form ("h_jet_trk_dphi_gt20_lt30_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
  TH2D* h2_jet_trk_dphi_gt20_lt30_cov = new TH2D (Form ("h2_jet_trk_dphi_gt20_lt30_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);

  TH1D* h_jet_trk_pt_ns = new TH1D (Form ("h_jet_trk_pt_ns_%s", tag), ";#it{p}_{T}^{ch} [GeV];(1/N_{jet}) (dN_{ch}/d#it{p}_{T}) [GeV^{-1}}", nPtChBins, pTChBins);
  TH2D* h2_jet_trk_pt_ns_cov = new TH2D (Form ("h2_jet_trk_pt_ns_cov_%s", tag), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{ch} [GeV];Covariance", nPtChBins, pTChBins, nPtChBins, pTChBins);
  TH1D* h_jet_trk_pt_perp = new TH1D (Form ("h_jet_trk_pt_perp_%s", tag), ";#it{p}_{T}^{ch} [GeV];(1/N_{jet}) (dN_{ch}/d#it{p}_{T}) [GeV^{-1}}", nPtChBins, pTChBins);
  TH2D* h2_jet_trk_pt_perp_cov = new TH2D (Form ("h2_jet_trk_pt_perp_cov_%s", tag), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{ch} [GeV];Covariance", nPtChBins, pTChBins, nPtChBins, pTChBins);
  TH1D* h_jet_trk_pt_as = new TH1D (Form ("h_jet_trk_pt_as_%s", tag), ";#it{p}_{T}^{ch} [GeV];(1/N_{jet}) (dN_{ch}/d#it{p}_{T}) [GeV^{-1}}", nPtChBins, pTChBins);
  TH2D* h2_jet_trk_pt_as_cov = new TH2D (Form ("h2_jet_trk_pt_as_cov_%s", tag), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{ch} [GeV];Covariance", nPtChBins, pTChBins, nPtChBins, pTChBins);


  // arrays for filling histograms & covariance matrices correctly
  double jet_pt_counts[nPtJBins];

  double jet_trk_dphi_gt0p5_lt1_counts[nDPhiBins];
  double jet_trk_dphi_gt1_lt1p5_counts[nDPhiBins];
  double jet_trk_dphi_gt1p5_lt2_counts[nDPhiBins];
  double jet_trk_dphi_gt2_lt4_counts[nDPhiBins];
  double jet_trk_dphi_gt4_lt6_counts[nDPhiBins];
  double jet_trk_dphi_gt6_lt8_counts[nDPhiBins];
  double jet_trk_dphi_gt8_lt10_counts[nDPhiBins];
  double jet_trk_dphi_gt10_lt15_counts[nDPhiBins];
  double jet_trk_dphi_gt15_lt20_counts[nDPhiBins];
  double jet_trk_dphi_gt20_lt30_counts[nDPhiBins];

  double jet_trk_pt_ns_counts[nPtChBins];
  double jet_trk_pt_perp_counts[nPtChBins];
  double jet_trk_pt_as_counts[nPtChBins];

  const int nEvts = (doMixing ? 20. : 1.) * (jetsTree->GetEntries ());
  const int nTrkEvts = tracksTree->GetEntries ();


  // setup centrality matching bins for event mixing according to job configuration
  double* mixCentBins = nullptr;
  int nMixCentBins = 0;

  if (Ispp ()) {
    if (doMixVar1)      { mixCentBins = ppMixVar1Bins;  nMixCentBins = nppMixVar1Bins;  }
    else if (doMixVar3) { mixCentBins = ppMixVar3Bins;  nMixCentBins = nppMixVar3Bins;  }
    else                { mixCentBins = ppMixBins;      nMixCentBins = nppMixBins;      }
  }
  else {
    if (doMixVar2)      { mixCentBins = fcalMixVar2Bins;  nMixCentBins = nFcalMixVar2Bins;  }
    else                { mixCentBins = fcalMixBins;      nMixCentBins = nFcalMixBins;      }
  }


  int iTrkEvt = 0;
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)                                
      cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;  

    for (int iX = 0; iX < nPtJBins; iX++)
      jet_pt_counts[iX] = 0;
    for (int iX = 0; iX < nDPhiBins; iX++) {
      jet_trk_dphi_gt0p5_lt1_counts[iX] = 0;
      jet_trk_dphi_gt1_lt1p5_counts[iX] = 0;
      jet_trk_dphi_gt1p5_lt2_counts[iX] = 0;
      jet_trk_dphi_gt2_lt4_counts[iX] = 0;
      jet_trk_dphi_gt4_lt6_counts[iX] = 0;
      jet_trk_dphi_gt6_lt8_counts[iX] = 0;
      jet_trk_dphi_gt8_lt10_counts[iX] = 0;
      jet_trk_dphi_gt10_lt15_counts[iX] = 0;
      jet_trk_dphi_gt15_lt20_counts[iX] = 0;
      jet_trk_dphi_gt20_lt30_counts[iX] = 0;
    }
    for (int iX = 0; iX < nPtChBins; iX++) {
      jet_trk_pt_ns_counts[iX] = 0;
      jet_trk_pt_perp_counts[iX] = 0;
      jet_trk_pt_as_counts[iX] = 0;
    }

    jetsTree->GetEntry (iEvt % jetsTree->GetEntries ());

    if (leading_jet == -1) {
      continue; // require a leading jet (i.e. any jet at all)
    }

    const double yboost = GetBoost (run_number);

    // TODO -- add event level weights!
    //const double ewgt = event_weight;
    const double ewgt = 1.;
    if (ewgt <= 0.) {
      continue;
    }

    h_evt_counts->Fill (0);
    h_evt_counts->Fill (1, ewgt);
    h_evt_counts->Fill (2, ewgt*ewgt);

    double nj = 0;
    double jwgt = 0;

    for (int iJet = 0; iJet < akt4_hi_jet_n; iJet++) {
      if (!MeetsJetPtCut (akt4_hi_jet_pt[iJet]))
        continue;
      const double thisjwgt = GetAkt4JetWeight (akt4_hi_jet_pt[iJet], akt4_hi_jet_eta[iJet], akt4_hi_jet_phi[iJet], 0.4);
      if (thisjwgt <= 0.)
        continue;

      nj += 1;
      jwgt += thisjwgt;
    }

    // skip events with no jets
    if (nj == 0) {
      continue;
    }

    h_jet_counts->Fill (0); // adds to number of jets (i.e. denominator)
    h_jet_counts->Fill (1, ewgt*jwgt);
    h_jet_counts->Fill (2, pow (ewgt*jwgt, 2));


    // do a mixed event
    if (doMixing) {

      // get the "centrality" bin for the current event
      const short iMixCent = GetBin (mixCentBins, nMixCentBins, Ispp () ? fcal_et_p : fcal_et_Pb);
      if (iMixCent < 0 || nMixCentBins <= iMixCent)
        continue;

      const int oldTrkEvt = iTrkEvt;
      bool goodEvent = false;
      do {
        iTrkEvt = (iTrkEvt + 1) % nTrkEvts;
        tracksTree->GetEntry (iTrkEvt);
        //goodEvent = true; // for disabling any mixing categories (must comment out mixing categories lines below!)

        // mixing categories
        goodEvent = goodEvent || GetBin (mixCentBins, nMixCentBins, Ispp () ? fcal_et_p_matching : fcal_et_Pb_matching) == iMixCent;
      }
      while (!goodEvent && iTrkEvt != oldTrkEvt);

      if (!goodEvent) {
        std::cout << "Could not find a good mixed event, so I am skipping this very interesting event!" << std::endl;
        continue;
      }
    } // now we have a good mixed event match (if applicable)


    for (int iJet = 0; iJet < akt4_hi_jet_n; iJet++) {
      h2_jet_eta_phi->Fill (akt4_hi_jet_eta[iJet], akt4_hi_jet_phi[iJet]);

      const double jwgt = GetAkt4JetWeight (akt4_hi_jet_pt[iJet], akt4_hi_jet_eta[iJet], akt4_hi_jet_phi[iJet], 0.4);
      if (jwgt <= 0.)
        continue;

      const short iPtJ = GetPtJBin (akt4_hi_jet_pt[iJet]);
      if (0 <= iPtJ && iPtJ < nPtJBins)
        jet_pt_counts[iPtJ] += jwgt;
    }

    for (int iX = 0; iX < nPtJBins; iX++) {
      h_jet_pt->SetBinContent (iX+1, h_jet_pt->GetBinContent (iX+1) + (ewgt)*(jet_pt_counts[iX]));
      for (int iY = 0; iY < nPtJBins; iY++)
        h2_jet_pt_cov->SetBinContent (iX+1, iY+1, h2_jet_pt_cov->GetBinContent (iX+1, iY+1) + (ewgt)*(jet_pt_counts[iX])*(jet_pt_counts[iY]));
    }


    // initialize all yields to 0
    for (int iX = 0; iX < nDPhiBins; iX++) {
      jet_trk_dphi_gt0p5_lt1_counts[iX] = 0;
      jet_trk_dphi_gt1_lt1p5_counts[iX] = 0;
      jet_trk_dphi_gt1p5_lt2_counts[iX] = 0;
      jet_trk_dphi_gt2_lt4_counts[iX] = 0;
      jet_trk_dphi_gt4_lt6_counts[iX] = 0;
      jet_trk_dphi_gt6_lt8_counts[iX] = 0;
      jet_trk_dphi_gt8_lt10_counts[iX] = 0;
      jet_trk_dphi_gt10_lt15_counts[iX] = 0;
      jet_trk_dphi_gt15_lt20_counts[iX] = 0;
      jet_trk_dphi_gt20_lt30_counts[iX] = 0;
    }
    for (int iX = 0; iX < nPtChBins; iX++) {
      jet_trk_pt_ns_counts[iX] = 0;
      jet_trk_pt_perp_counts[iX] = 0;
      jet_trk_pt_as_counts[iX] = 0;
    }

    // loop over all jets in the event
    for (int iJet = 0; iJet < akt4_hi_jet_n; iJet++) {
      if (!MeetsJetPtCut (akt4_hi_jet_pt[iJet]))
        continue;
      const double thisjwgt = GetAkt4JetWeight (akt4_hi_jet_pt[iJet], akt4_hi_jet_eta[iJet], akt4_hi_jet_phi[iJet], 0.4);
      if (thisjwgt <= 0.)
        continue;

      // correlate charged particles with this jet  
      for (int iTrk = 0; iTrk < trk_n; iTrk++) {
        if (trk_pt[iTrk] < pTChBins[0])
          continue;

        const short iPtCh = GetPtChBin (trk_pt[iTrk]);
        if (iPtCh < 0 || iPtCh >= nPtChBins)
          continue;

        //// TODO -- better rapidity selection by considering mass more carefully
        //const double trk_m = pion_mass; // assume everyone is a pi^+
        //const double trk_en = std::sqrt (std::pow (trk_m, 2) + trk_pt[iTrk] * std::cosh ((double)(trk_eta[iTrk])));
        //const double trk_pz = ((double)trk_pt[iTrk]) * std::sinh ((double)trk_eta[iTrk]);
        //const double trk_y = (trk_en > trk_pz ? 0.5 * std::log ((trk_en + trk_pz)/(trk_en - trk_pz)) : 0.);

        const double trk_y = trk_eta[iTrk];
        
        if (fabs (trk_y - yboost) > 2.5 - 0.465)
          continue;

        const float dphi = DeltaPhi (akt4_hi_jet_phi[iJet], trk_phi[iTrk]);
        const short iDPhi = GetDPhiBin (dphi);
        if (iDPhi < 0 || nDPhiBins < iDPhi)
          continue;

        const float teff = h2_trk_eff->GetBinContent (h2_trk_eff->FindBin (trk_eta[iTrk], trk_pt[iTrk]));
        const float tpur = h2_trk_pur->GetBinContent (h2_trk_pur->FindBin (trk_eta[iTrk], trk_pt[iTrk]));
        const float twgt = (teff > 0. ? tpur / teff : 0.);

        if (0.5 < trk_pt[iTrk] && trk_pt[iTrk] < 1)
          jet_trk_dphi_gt0p5_lt1_counts[iDPhi] += twgt;
        else if (1 < trk_pt[iTrk] && trk_pt[iTrk] < 1.5)
          jet_trk_dphi_gt1_lt1p5_counts[iDPhi] += twgt;
        else if (1.5 < trk_pt[iTrk] && trk_pt[iTrk] < 2)
          jet_trk_dphi_gt1p5_lt2_counts[iDPhi] += twgt;
        else if (2 < trk_pt[iTrk] && trk_pt[iTrk] < 4)
          jet_trk_dphi_gt2_lt4_counts[iDPhi] += twgt;
        else if (4 < trk_pt[iTrk] && trk_pt[iTrk] < 6)
          jet_trk_dphi_gt4_lt6_counts[iDPhi] += twgt;
        else if (6 < trk_pt[iTrk] && trk_pt[iTrk] < 8)
          jet_trk_dphi_gt6_lt8_counts[iDPhi] += twgt;
        else if (8 < trk_pt[iTrk] && trk_pt[iTrk] < 10)
          jet_trk_dphi_gt8_lt10_counts[iDPhi] += twgt;
        else if (10 < trk_pt[iTrk] && trk_pt[iTrk] < 15)
          jet_trk_dphi_gt10_lt15_counts[iDPhi] += twgt;
        else if (15 < trk_pt[iTrk] && trk_pt[iTrk] < 20)
          jet_trk_dphi_gt15_lt20_counts[iDPhi] += twgt;
        else if (20 < trk_pt[iTrk] && trk_pt[iTrk] < 30)
          jet_trk_dphi_gt20_lt30_counts[iDPhi] += twgt;

        if (dphi < M_PI/8.)
          jet_trk_pt_ns_counts[iPtCh] += twgt;
        else if (M_PI/3. < dphi && dphi < 2.*M_PI/3.)
          jet_trk_pt_perp_counts[iPtCh] += twgt;
        else if (dphi > 7.*M_PI/8.)
          jet_trk_pt_as_counts[iPtCh] += twgt;

      } // end loop over tracks
    } // end loop over jets

    if (nj == 0) {
      continue;
    }


    // calculate per-jet hadron yields for that event by dividing out the number of jets
    for (int iX = 0; iX < nDPhiBins; iX++) {
      jet_trk_dphi_gt0p5_lt1_counts[iX] = jet_trk_dphi_gt0p5_lt1_counts[iX] / nj;
      jet_trk_dphi_gt1_lt1p5_counts[iX] = jet_trk_dphi_gt1_lt1p5_counts[iX] / nj;
      jet_trk_dphi_gt1p5_lt2_counts[iX] = jet_trk_dphi_gt1p5_lt2_counts[iX] / nj;
      jet_trk_dphi_gt2_lt4_counts[iX] = jet_trk_dphi_gt2_lt4_counts[iX] / nj;
      jet_trk_dphi_gt4_lt6_counts[iX] = jet_trk_dphi_gt4_lt6_counts[iX] / nj;
      jet_trk_dphi_gt6_lt8_counts[iX] = jet_trk_dphi_gt6_lt8_counts[iX] / nj;
      jet_trk_dphi_gt8_lt10_counts[iX] = jet_trk_dphi_gt8_lt10_counts[iX] / nj;
      jet_trk_dphi_gt10_lt15_counts[iX] = jet_trk_dphi_gt10_lt15_counts[iX] / nj;
      jet_trk_dphi_gt15_lt20_counts[iX] = jet_trk_dphi_gt15_lt20_counts[iX] / nj;
      jet_trk_dphi_gt20_lt30_counts[iX] = jet_trk_dphi_gt20_lt30_counts[iX] / nj;
    }
    for (int iX = 0; iX < nPtChBins; iX++) {
      jet_trk_pt_ns_counts[iX] = jet_trk_pt_ns_counts[iX] / nj;
      jet_trk_pt_perp_counts[iX] = jet_trk_pt_perp_counts[iX] / nj;
      jet_trk_pt_as_counts[iX] = jet_trk_pt_as_counts[iX] / nj;
    }


    // store results in averaging, covariance histograms
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt0p5_lt1->SetBinContent (iX+1, h_jet_trk_dphi_gt0p5_lt1->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt0p5_lt1_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt0p5_lt1_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt0p5_lt1_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt0p5_lt1_counts[iX])*(jet_trk_dphi_gt0p5_lt1_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt1_lt1p5->SetBinContent (iX+1, h_jet_trk_dphi_gt1_lt1p5->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt1_lt1p5_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt1_lt1p5_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt1_lt1p5_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt1_lt1p5_counts[iX])*(jet_trk_dphi_gt1_lt1p5_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt1p5_lt2->SetBinContent (iX+1, h_jet_trk_dphi_gt1p5_lt2->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt1p5_lt2_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt1p5_lt2_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt1p5_lt2_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt1p5_lt2_counts[iX])*(jet_trk_dphi_gt1p5_lt2_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt2_lt4->SetBinContent (iX+1, h_jet_trk_dphi_gt2_lt4->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt2_lt4_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt2_lt4_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt2_lt4_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt2_lt4_counts[iX])*(jet_trk_dphi_gt2_lt4_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt4_lt6->SetBinContent (iX+1, h_jet_trk_dphi_gt4_lt6->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt4_lt6_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt4_lt6_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt4_lt6_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt4_lt6_counts[iX])*(jet_trk_dphi_gt4_lt6_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt6_lt8->SetBinContent (iX+1, h_jet_trk_dphi_gt6_lt8->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt6_lt8_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt6_lt8_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt6_lt8_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt6_lt8_counts[iX])*(jet_trk_dphi_gt6_lt8_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt8_lt10->SetBinContent (iX+1, h_jet_trk_dphi_gt8_lt10->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt8_lt10_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt8_lt10_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt8_lt10_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt8_lt10_counts[iX])*(jet_trk_dphi_gt8_lt10_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt10_lt15->SetBinContent (iX+1, h_jet_trk_dphi_gt10_lt15->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt10_lt15_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt10_lt15_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt10_lt15_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt10_lt15_counts[iX])*(jet_trk_dphi_gt10_lt15_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt15_lt20->SetBinContent (iX+1, h_jet_trk_dphi_gt15_lt20->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt15_lt20_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt15_lt20_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt15_lt20_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt15_lt20_counts[iX])*(jet_trk_dphi_gt15_lt20_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt20_lt30->SetBinContent (iX+1, h_jet_trk_dphi_gt20_lt30->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt20_lt30_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt20_lt30_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt20_lt30_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt20_lt30_counts[iX])*(jet_trk_dphi_gt20_lt30_counts[iY]));
    }


    for (int iX = 0; iX < nPtChBins; iX++) {
      h_jet_trk_pt_ns->SetBinContent (iX+1, h_jet_trk_pt_ns->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_pt_ns_counts[iX]));
      for (int iY = 0; iY < nPtChBins; iY++)
        h2_jet_trk_pt_ns_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_pt_ns_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_pt_ns_counts[iX])*(jet_trk_pt_ns_counts[iY]));
    }
    for (int iX = 0; iX < nPtChBins; iX++) {
      h_jet_trk_pt_perp->SetBinContent (iX+1, h_jet_trk_pt_perp->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_pt_perp_counts[iX]));
      for (int iY = 0; iY < nPtChBins; iY++)
        h2_jet_trk_pt_perp_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_pt_perp_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_pt_perp_counts[iX])*(jet_trk_pt_perp_counts[iY]));
    }
    for (int iX = 0; iX < nPtChBins; iX++) {
      h_jet_trk_pt_as->SetBinContent (iX+1, h_jet_trk_pt_as->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_pt_as_counts[iX]));
      for (int iY = 0; iY < nPtChBins; iY++)
        h2_jet_trk_pt_as_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_pt_as_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_pt_as_counts[iX])*(jet_trk_pt_as_counts[iY]));
    }

  } // end loop over events
  cout << "Finished event loop." << endl;


  SaferDelete (&h2_trk_eff);
  SaferDelete (&h2_trk_pur);


  h_evt_counts->Write ();
  h_jet_counts->Write ();

  h_jet_pt->Write ();
  h2_jet_pt_cov->Write ();
  h2_jet_eta_phi->Write ();

  h_jet_trk_dphi_gt0p5_lt1->Write ();
  h2_jet_trk_dphi_gt0p5_lt1_cov->Write ();
  h_jet_trk_dphi_gt1_lt1p5->Write ();
  h2_jet_trk_dphi_gt1_lt1p5_cov->Write ();
  h_jet_trk_dphi_gt1p5_lt2->Write ();
  h2_jet_trk_dphi_gt1p5_lt2_cov->Write ();
  h_jet_trk_dphi_gt2_lt4->Write ();
  h2_jet_trk_dphi_gt2_lt4_cov->Write ();
  h_jet_trk_dphi_gt4_lt6->Write ();
  h2_jet_trk_dphi_gt4_lt6_cov->Write ();
  h_jet_trk_dphi_gt6_lt8->Write ();
  h2_jet_trk_dphi_gt6_lt8_cov->Write ();
  h_jet_trk_dphi_gt8_lt10->Write ();
  h2_jet_trk_dphi_gt8_lt10_cov->Write ();
  h_jet_trk_dphi_gt10_lt15->Write ();
  h2_jet_trk_dphi_gt10_lt15_cov->Write ();
  h_jet_trk_dphi_gt15_lt20->Write ();
  h2_jet_trk_dphi_gt15_lt20_cov->Write ();
  h_jet_trk_dphi_gt20_lt30->Write ();
  h2_jet_trk_dphi_gt20_lt30_cov->Write ();

  h_jet_trk_pt_ns->Write ();
  h2_jet_trk_pt_ns_cov->Write ();
  h_jet_trk_pt_perp->Write ();
  h2_jet_trk_pt_perp_cov->Write ();
  h_jet_trk_pt_as->Write ();
  h2_jet_trk_pt_as_cov->Write ();


  outFile->Close ();
}




void RunCorrelator (const char* collSys, const char* tag, const char* outFileName, const char* jetsInFileName, const char* tracksInFileName = nullptr) {

  if      (TString (collSys) == ToTString (CollisionSystem::pPb16s5TeV))  { collisionSystem = CollisionSystem::pPb16s5TeV; }
  else if (TString (collSys) == ToTString (CollisionSystem::pp17))        { collisionSystem = CollisionSystem::pp17;       }
  else {
    std::cout << "In RunCorrelator.C: Invalid or unsupported collision system, exiting." << std::endl;
    return;
  }

  if (jet_min_pt == -2)
    jet_min_pt = (double) std::atof (std::getenv ("JET_MIN_PT"));
  if (jet_max_pt == -2)
    jet_max_pt = (double) std::atof (std::getenv ("JET_MAX_PT"));
  doMixing = (tracksInFileName != nullptr);


  doMixVar1 = (std::string (jetsInFileName).find ("MixCatVar1") != std::string::npos || (doMixing && std::string (tracksInFileName).find ("MixCatVar1") != std::string::npos));
  doMixVar2 = (std::string (jetsInFileName).find ("MixCatVar2") != std::string::npos || (doMixing && std::string (tracksInFileName).find ("MixCatVar2") != std::string::npos));
  doMixVar3 = (std::string (jetsInFileName).find ("MixCatVar3") != std::string::npos || (doMixing && std::string (tracksInFileName).find ("MixCatVar3") != std::string::npos));


  TFile* jetsInFile = new TFile (jetsInFileName, "read");
  TFile* tracksInFile = (tracksInFileName != nullptr ? new TFile (tracksInFileName, "read") : nullptr);

  TTree* jetsInTree = (TTree*) jetsInFile->Get (IspPb () ? "pPbTree" : "ppTree");
  TTree* tracksInTree = (tracksInFileName != nullptr ? (TTree*) tracksInFile->Get (IspPb () ? "pPbTree" : "ppTree") : jetsInTree);

  Correlator (TString (tag), outFileName, jetsInTree, tracksInTree);

  jetsInFile->Close ();
  if (tracksInFile)
    tracksInFile->Close ();

  return;
}




int main (int argc, char** argv) {
  assert (argc >= 7);

  jet_min_pt = (double) std::atof (argv[1]);
  jet_max_pt = (double) std::atof (argv[2]);

  if (argc == 8) {
    std::cout << "Main will execute RunCorrelator (const char* collSys, const char* tag, const char* outFileName, const char* jetsInFileName, const char* tracksInFileName)" << std::endl;
    RunCorrelator (argv[3], argv[4], argv[5], argv[6], argv[7]);
  }

  else if (argc == 7) {
    std::cout << "Main will execute RunCorrelator (const char* collSys, const char* tag, const char* outFileName, const char* jetsInFileName)" << std::endl;
    RunCorrelator (argv[3], argv[4], argv[5], argv[6]);
  }

  else {
    std::cout << "Undefined behavior for " << argc << " arguments, exiting" << std::endl;
    return 1;
  }

  return 0;
}

#endif
