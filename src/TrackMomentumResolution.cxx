#ifndef __TrackMomentumResolution_cxx__
#define __TrackMomentumResolution_cxx__

#include "TrackMomentumResolution.h"
#include "Params.h"
#include "CentralityDefs.h"
#include "TreeVariables.h"
#include "Trigger.h"
#include "LocalUtilities.h"

#include <Utilities.h>

#include <TChain.h>
#include <TSystem.h>
#include <TH2D.h>
#include <TF1.h>

#include <iostream>

using namespace std;

namespace JetHadronCorrelations {

bool TrackMomentumResolution (const char* directory,
                                 const int dataSet,
                                 const char* inFileName,
                                 const char* eventWeightsFileName) {
 
  std::cout << "Info: In TrackMomentumResolution.cxx: Entered TrackMomentumResolution routine." << std::endl;

  SetupDirectories ("TrackMomentumResolution");

  if (IsCollisions ()) {
    std::cout << "Error: In TrackMomentumResolution.cxx: Trying to calculate tracking performance in data! Quitting." << std::endl;
    return false;
  }

  if (IsDataOverlay ())
    std::cout << "Info: In TrackMomentumResolution.cxx: Running over data overlay, will check data conditions" << std::endl;
  if (IsHijing ())
    std::cout << "Info: In TrackMomentumResolution.cxx: Running over Hijing sample" << std::endl;

  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  std::cout << "Info: In TrackMomentumResolution.cxx: File Identifier: " << identifier << std::endl;
  std::cout << "Info: In TrackMomentumResolution.cxx: Saving output to " << rootPath << std::endl;


  // creates a file identifier pattern that we will use to identify the directory containing the input files
  TString fileIdentifier;
  if (TString (inFileName) == "") {
    if (IsCollisions ()) {
      fileIdentifier = to_string (dataSet);
    }
    else {
      std::cout << "Error: In TrackMomentumResolution.cxx: Cannot identify this MC file! Quitting." << std::endl;
      return false;
    }
  }
  else
    fileIdentifier = inFileName;


  // opens a TTree as a TChain from all files in a directory matching the file identifier
  TChain* tree = new TChain ("bush", "bush");
  {
    TString pattern = "*.root";
    std::cout << "DataPath = " << dataPath << std::endl;
    auto dir = gSystem->OpenDirectory (dataPath + directory);
    while (const char* f = gSystem->GetDirEntry (dir)) {
      TString file = TString (f);
      if (!file.Contains (fileIdentifier))
        continue;
      std::cout << "Adding " << dataPath + directory + "/" + file + "/*.root" << " to TChain" << std::endl;
      tree->Add (dataPath + directory + "/" + file + "/*.root");
      break;
    }
    std::cout << "Info: In TrackMomentumResolution.cxx: Chain has " << tree->GetListOfFiles ()->GetEntries () << " files, " << tree->GetEntries () << " entries" << std::endl;
  }


  if (!IsHijing ()) {
    assert (crossSectionPicoBarns > 0);
    assert (mcFilterEfficiency > 0);
    assert (mcNumberEvents > 0);
  }


  // variables for filtering MC truth
  const float truth_jet_min_pt = GetJZXR04MinPt (TString (inFileName));
  const float truth_jet_max_pt = GetJZXR04MaxPt (TString (inFileName));
  if (truth_jet_min_pt != 0)
    std::cout << "Checking for leading truth jet with pT > " << truth_jet_min_pt << std::endl;
  if (truth_jet_max_pt != FLT_MAX)
    std::cout << "Checking for leading truth jet with pT < " << truth_jet_max_pt << std::endl;



  tree->SetBranchAddress ("run_number",     &run_number);
  tree->SetBranchAddress ("event_number",   &event_number);
  tree->SetBranchAddress ("lumi_block",     &lumi_block);
  tree->SetBranchAddress ("isOOTPU",        &isOOTPU);
  tree->SetBranchAddress ("BlayerDesyn",    &BlayerDesyn);
  tree->SetBranchAddress ("passes_toroid",  &passes_toroid);


  tree->SetBranchAddress ("actualInteractionsPerCrossing",  &actualInteractionsPerCrossing);
  tree->SetBranchAddress ("averageInteractionsPerCrossing", &averageInteractionsPerCrossing);


  if (!IsCollisions ()) {
    tree->SetBranchAddress ("mcEventWeights", &mcEventWeights);
  }
  if (!IsCollisions () && IsHijing ()) {
    tree->SetBranchAddress ("truth_event_n",      &(truth_event_n));
    tree->SetBranchAddress ("nPart1",             &nPart1);
    tree->SetBranchAddress ("nPart2",             &nPart2);
    tree->SetBranchAddress ("impactParameter",    &impactParameter);
    tree->SetBranchAddress ("nColl",              &nColl);
    tree->SetBranchAddress ("nSpectatorNeutrons", &nSpectatorNeutrons);
    tree->SetBranchAddress ("nSpectatorProtons",  &nSpectatorProtons);
    tree->SetBranchAddress ("eccentricity",       &eccentricity);
    tree->SetBranchAddress ("eventPlaneAngle",    &eventPlaneAngle);
  }


  tree->SetBranchAddress ("nvert",     &nvert);
  tree->SetBranchAddress ("vert_x",    &vert_x);
  tree->SetBranchAddress ("vert_y",    &vert_y);
  tree->SetBranchAddress ("vert_z",    &vert_z);
  tree->SetBranchAddress ("vert_ntrk", &vert_ntrk);
  tree->SetBranchAddress ("vert_type", &vert_type);


  tree->SetBranchAddress ("fcalA_et",       &fcalA_et);
  tree->SetBranchAddress ("fcalC_et",       &fcalC_et);


  if (!Ispp ()) {
    tree->SetBranchAddress ("ZdcCalibEnergy_A",   &ZdcCalibEnergy_A);
    tree->SetBranchAddress ("ZdcCalibEnergy_C",   &ZdcCalibEnergy_C);
    tree->SetBranchAddress ("ZdcRawEnergy_A",     &ZdcRawEnergy_A);
    tree->SetBranchAddress ("ZdcRawEnergy_C",     &ZdcRawEnergy_C);
  }


  if (!IsCollisions ()) {
    tree->SetBranchAddress ("truth_trk_n",        &truth_trk_n);
    tree->SetBranchAddress ("truth_trk_pt",       &truth_trk_pt);
    tree->SetBranchAddress ("truth_trk_eta",      &truth_trk_eta);
    tree->SetBranchAddress ("truth_trk_phi",      &truth_trk_phi);
    tree->SetBranchAddress ("truth_trk_charge",   &truth_trk_charge);
    tree->SetBranchAddress ("truth_trk_pdgid",    &truth_trk_pdgid);
    tree->SetBranchAddress ("truth_trk_barcode",  &truth_trk_barcode);
    tree->SetBranchAddress ("truth_trk_isHadron", &truth_trk_isHadron);
  }

  tree->SetBranchAddress ("ntrk",                   &trk_n);
  tree->SetBranchAddress ("trk_pt",                 &trk_pt);
  tree->SetBranchAddress ("trk_eta",                &trk_eta);
  tree->SetBranchAddress ("trk_phi",                &trk_phi);
  tree->SetBranchAddress ("trk_charge",             &trk_charge);
  tree->SetBranchAddress ("trk_HITight",            &trk_HITight);
  tree->SetBranchAddress ("trk_HILoose",            &trk_HILoose);
  tree->SetBranchAddress ("trk_TightPrimary",       &trk_TightPrimary);
  tree->SetBranchAddress ("trk_d0",                 &trk_d0);
  tree->SetBranchAddress ("trk_d0sig",              &trk_d0sig);
  tree->SetBranchAddress ("trk_z0",                 &trk_z0);
  tree->SetBranchAddress ("trk_z0sig",              &trk_z0sig);
  tree->SetBranchAddress ("trk_theta",              &trk_theta);
  tree->SetBranchAddress ("trk_vz",                 &trk_vz);
  if (!IsCollisions ()) {
    tree->SetBranchAddress ("trk_prob_truth",     &trk_prob_truth);
    tree->SetBranchAddress ("trk_truth_pt",       &trk_truth_pt);
    tree->SetBranchAddress ("trk_truth_eta",      &trk_truth_eta);
    tree->SetBranchAddress ("trk_truth_phi",      &trk_truth_phi);
    tree->SetBranchAddress ("trk_truth_charge",   &trk_truth_charge);
    tree->SetBranchAddress ("trk_truth_type",     &trk_truth_type);
    tree->SetBranchAddress ("trk_truth_orig",     &trk_truth_orig);
    tree->SetBranchAddress ("trk_truth_barcode",  &trk_truth_barcode);
    tree->SetBranchAddress ("trk_truth_pdgid",    &trk_truth_pdgid);
    tree->SetBranchAddress ("trk_truth_vz",       &trk_truth_vz);
    tree->SetBranchAddress ("trk_truth_nIn",      &trk_truth_nIn);
    tree->SetBranchAddress ("trk_truth_isHadron", &trk_truth_isHadron);
  }


  std::cout << "Info : In TrackMomentumResolution.cxx: Saving histograms to " << Form ("%s/%s.root", rootPath.Data (), identifier.Data ()) << std::endl;
  TFile* outFile = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");

  TString sys;
  if (Ispp ())        sys = "pp";
  else if (IspPb ())  sys = "pPb";
  else                sys = "???";

  const int nFinerEtaTrkBins = 40;
  const double* finerEtaTrkBins = linspace (-2.5, 2.5, nFinerEtaTrkBins);

  const double pTchBins[] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 50, 55, 60, 70, 80, 100};
  const int nPtchBins = sizeof (pTchBins) / sizeof (pTchBins[0]) - 1;

  TH1D*** h_tmr = new TH1D**[nPtchBins];

  TH2D** h2_ptreco_pttruth = new TH2D*[nFinerEtaTrkBins];

  const int nrBins = 200;
  double* rBins = new double[nrBins+1];
  for (int i = 0; i <= nrBins; i++) {
    if (i < 0.05*nrBins)
      rBins[i] = -1 + 0.5*i/nrBins;
    else if (i < 0.85*nrBins)
      rBins[i] = -0.5 + (i-0.05*nrBins)/(0.8*nrBins);
    else
      rBins[i] = 0.5 + (i-0.85*nrBins)*3.5/(0.15*nrBins);
  }

  for (int iPtch = 0; iPtch < nPtchBins; iPtch++) {
    h_tmr[iPtch] = new TH1D*[nFinerEtaTrkBins];
    for (int iEta = 0; iEta < nFinerEtaTrkBins; iEta++) {
      h_tmr[iPtch][iEta] = new TH1D (Form ("h_tmr_%s_iPtch%i_iEta%i", sys.Data (), iPtch, iEta), "Track resolution, #it{p}_{T}^{truth} / #it{p}_{T}^{reco} - 1", nrBins, rBins);
      h_tmr[iPtch][iEta]->Sumw2 ();
    }
  }

  for (int iEta = 0; iEta < nFinerEtaTrkBins; iEta++) {
    h2_ptreco_pttruth[iEta] = new TH2D (Form ("h2_ptreco_pttruth_%s_iEta%i", sys.Data (), iEta), ";#it{p}_{T}^{truth} [GeV];#it{p}_{T}^{reco} [GeV];Counts", nPtChBins, pTChBins, nPtChBins, pTChBins);
    h2_ptreco_pttruth[iEta]->Sumw2 ();
  }
  

  const JetRadius r0p4 = JetRadius::R0p4;
  const int nEvts = tree->GetEntries ();

  // Loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      std::cout << "Info: In TrackMomentumResolution.cxx: Event loop " << iEvt / (nEvts / 100) << "\% done...\r" << std::flush;
    tree->GetEntry (iEvt);

    // vertexing cuts, require no pileup vertices and primary vertex with |vz| < 150mm
    {
      bool hasPrimary = false;
      bool hasPileup = false;
      float vz = -999;
      for (int iVert = 0; iVert < nvert; iVert++) {
        if (vert_type[iVert] == 1) {
          hasPrimary = true;
          vz = vert_z[iVert];
        }
        else if (vert_type[iVert] == 3)
          hasPileup = true;
      }
      if (hasPileup || std::fabs (vz) > 150 || !hasPrimary)
        continue;
    }


    // Filter sample based on min/max of pThat range
    if (!IsHijing ()) {
      int iLTJ = -1;
      const int nTJ = GetAktTruthJetN (r0p4);
      for (int iTJ = 0; iTJ < nTJ; iTJ++) {
        if (iLTJ == -1 || GetAktTruthJetPt (iTJ, r0p4) > GetAktTruthJetPt (iLTJ, r0p4))
          iLTJ = iTJ;
      }

      if (iLTJ == -1 || GetAktTruthJetPt (iLTJ, r0p4) < truth_jet_min_pt || GetAktTruthJetPt (iLTJ, r0p4) > truth_jet_max_pt)
        continue;
    }


    //const float eventWeight = 1;
    //const float eventWeight = ((isPbPb && !isHijing) ? h_weights->GetBinContent (h_weights->FindFixBin (doNchWeighting ? ntrk : fcalA_et+fcalC_et)) : 1); // weight is 1 for pp
    //const float eventWeight = (isPbPb ? h_weights->GetBinContent (h_weights->FindFixBin (doNchWeighting ? ntrk : fcalA_et+fcalC_et)) : 1); // weight is 1 for pp
    for (int iTrk = 0; iTrk < trk_n; iTrk++) {
      if (!MeetsTrackCuts (iTrk))
        continue;

      const bool isTruthMatched = (trk_prob_truth[iTrk] > 0.5);

      const bool isFake = !isTruthMatched;
      const bool isSecondary = isTruthMatched && (trk_truth_barcode[iTrk] <= 0 || 200000 <= trk_truth_barcode[iTrk]);
      //const bool isChargedPion = isTruthMatched && std::abs (trk_truth_pdgid[iTrk]) == 211;
      const bool isStrangeBaryon = isTruthMatched && (std::abs (trk_truth_pdgid[iTrk]) == 3112 || std::abs (trk_truth_pdgid[iTrk]) == 3222 || std::abs (trk_truth_pdgid[iTrk]) == 3312 || std::abs (trk_truth_pdgid[iTrk]) == 3334);

      // primary tracks are non-fake, non-secondary tracks. Strange baryons are excluded too.
      const bool isPrimary = !isFake && !isSecondary && !isStrangeBaryon;

      if (!isPrimary)
        continue; // restrict to only primary tracks

      // now cut on the truth-level info
      if (trk_truth_charge[iTrk] == 0 ||
          std::fabs (trk_truth_eta[iTrk]) > 2.5 ||
          !(trk_truth_isHadron[iTrk]))
        continue;

      short iPtch = -1;
      if (pTchBins[0] <= trk_truth_pt[iTrk]) {
        iPtch = 0;
        while (iPtch < nPtchBins && pTchBins[iPtch+1] < trk_truth_pt[iTrk]) iPtch++;
      }

      short iEta = -1;
      if (finerEtaTrkBins[0] <= trk_truth_eta[iTrk]) {
        iEta = 0;
        while (iEta < nFinerEtaTrkBins && finerEtaTrkBins[iEta+1] < trk_truth_eta[iTrk]) iEta++;
      }

      h2_ptreco_pttruth[iEta]->Fill (trk_truth_pt[iTrk], trk_pt[iTrk]);

      if (iPtch >= 0 && iPtch < nPtchBins && iEta >= 0 && iEta < nFinerEtaTrkBins) {
        h_tmr[iPtch][iEta]->Fill (trk_truth_pt[iTrk]/trk_pt[iTrk] - 1.);
      }
    }
  } // end event loop
  std::cout << std::endl << "Info: In TrackMomentumResolution.cxx: Finished event loop." << std::endl;


  SaferDelete (&tree);


  outFile->cd ();

  for (int iPtch = 0; iPtch < nPtchBins; iPtch++) {

    for (int iEta = 0; iEta < nFinerEtaTrkBins; iEta++) {

      h_tmr[iPtch][iEta]->Write ();
      SaferDelete (&(h_tmr[iPtch][iEta]));

    } // end loop over iEta

  } // end loop over iPtch


  for (int iEta = 0; iEta < nFinerEtaTrkBins; iEta++) {

    h2_ptreco_pttruth[iEta]->Write ();
    SaferDelete (&(h2_ptreco_pttruth[iEta]));

  } // end loop over iPtch


  outFile->Write (0, TObject::kOverwrite);
  outFile->Close ();
  SaferDelete (&outFile);

  return true;
}

} // end namespace

#endif
