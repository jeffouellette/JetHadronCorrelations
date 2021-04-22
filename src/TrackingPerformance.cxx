#ifndef __TrackingPerformance_cxx__
#define __TrackingPerformance_cxx__

#include "TrackingPerformance.h"
#include "Params.h"
#include "TreeVariables.h"
#include "Trigger.h"
#include "LocalUtilities.h"

#include <Utilities.h>

#include <TChain.h>
#include <TSystem.h>
#include <TH2D.h>

#include <iostream>

using namespace std;

namespace JetHadronCorrelations {

bool TrackingPerformance (const char* directory,
                          const int dataSet,
                          const char* inFileName,
                          const char* eventWeightsFileName) {
 
  cout << "Info: In TrackingPerformance.cxx: Entered TrackingPerformance routine." << endl;
  cout << "Info: In TrackingPerformance.cxx: Printing systematic onfiguration:";
  cout << "\n\tDoHITightVar (): " << DoHITightVar () << endl;
  cout << "\n\tDoPionsOnlyVar (): " << DoPionsOnlyVar () << endl;

  SetupDirectories ("TrackingPerformance");

  if (IsCollisions ()) {
    cout << "Error: In TrackingPerformance.cxx: Trying to calculate tracking performance in data! Quitting." << endl;
    return false;
  }

  if (IsDataOverlay ())
    cout << "Info: In TrackingPerformance.cxx: Running over data overlay, will check data conditions" << endl;
  if (IsHijing ())
    cout << "Info: In TrackingPerformance.cxx: Running over Hijing sample" << endl;

  // this identifier here corresponds to the name of the output file
  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  std::cout << "Info: In TrackingPerformance.cxx: File Identifier: " << identifier << std::endl;
  std::cout << "Info: In TrackingPerformance.cxx: Saving output to " << rootPath << std::endl;


  // creates a file identifier pattern that we will use to identify the directory containing the input files
  TString fileIdentifier;
  if (TString (inFileName) == "") {
    if (IsCollisions ()) {
      fileIdentifier = to_string (dataSet);
    }
    else {
      std::cout << "Error: In TrackingPerformance.cxx: Cannot identify this MC file! Quitting." << std::endl;
      return false;
    }
  }
  else
    fileIdentifier = inFileName;


  // opens a TTree as a TChain from all files in a directory matching the file identifier
  TChain* tree = new TChain ("bush", "bush");
  {
    TString pattern = "*.root";
    std::cout << "DataPath = " << dataPath;
    auto dir = gSystem->OpenDirectory (dataPath + directory);
    while (const char* f = gSystem->GetDirEntry (dir)) {
      TString file = TString (f);
      if (!file.Contains (fileIdentifier))
        continue;
      std::cout << "Adding " << dataPath + directory + "/" + file + "/*.root" << " to TChain" << std::endl;
      tree->Add (dataPath + directory + "/" + file + "/*.root");
      break;
    }
    std::cout << "Info: In TrackingPerformance.cxx: Chain has " << tree->GetListOfFiles ()->GetEntries () << " files, " << tree->GetEntries () << " entries" << std::endl;
  }

  // TODO: event weights
  //TFile* eventWeightsFile = nullptr;
  //TH1D* h_weights = nullptr;

  //eventWeightsFile = new TFile (eventWeightsFileName, "read");
  //h_weights = (TH1D*) eventWeightsFile->Get (Form ("h_PbPb%s_weights_%s", doNchWeighting ? "Nch" : "FCal", isHijing ? "hijing" : "mc"));
  //cout << "Info: In TrackingPerformance.cxx: Found FCal weighting histogram, " << h_weights->GetName () << endl;


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
  tree->SetBranchAddress ("trk_HItight",            &trk_HItight);
  tree->SetBranchAddress ("trk_HIloose",            &trk_HIloose);
  tree->SetBranchAddress ("trk_TightPrimary",       &trk_TightPrimary);
  tree->SetBranchAddress ("trk_d0",                 &trk_d0);
  tree->SetBranchAddress ("trk_d0sig",              &trk_d0sig);
  tree->SetBranchAddress ("trk_z0",                 &trk_z0);
  tree->SetBranchAddress ("trk_z0sig",              &trk_z0sig);
  tree->SetBranchAddress ("trk_theta",              &trk_theta);
  tree->SetBranchAddress ("trk_vz",                 &trk_vz);
  tree->SetBranchAddress ("trk_nBLayerHits",        &trk_nBLayerHits);
  tree->SetBranchAddress ("trk_nBLayerSharedHits",  &trk_nBLayerSharedHits);
  tree->SetBranchAddress ("trk_nPixelHits",         &trk_nPixelHits);
  tree->SetBranchAddress ("trk_nPixelDeadSensors",  &trk_nPixelDeadSensors);
  tree->SetBranchAddress ("trk_nPixelSharedHits",   &trk_nPixelSharedHits);
  tree->SetBranchAddress ("trk_nSCTHits",           &trk_nSCTHits);
  tree->SetBranchAddress ("trk_nSCTDeadSensors",    &trk_nSCTDeadSensors);
  tree->SetBranchAddress ("trk_nSCTSharedHits",     &trk_nSCTSharedHits);
  tree->SetBranchAddress ("trk_nTRTHits",           &trk_nTRTHits);
  tree->SetBranchAddress ("trk_nTRTSharedHits",     &trk_nTRTSharedHits);
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


  cout << "Info : In TrackingPerformance.cxx: Saving histograms to " << Form ("%s/%s/%s.root", rootPath.Data (), IsHijing () ? "../Hijing/" : "", identifier.Data ()) << endl;
  TFile* outFile = new TFile (Form ("%s/%s/%s.root", rootPath.Data (), IsHijing () ? "../Hijing/" : "", identifier.Data ()), "recreate");

  TString sys;
  if (Ispp ())        sys = "pp";
  else if (IspPb ())  sys = "pPb";
  else                sys = "???";

  const int numFinerEtaTrkBins = 40;
  const double* finerEtaTrkBins = linspace (-2.5, 2.5, numFinerEtaTrkBins);

  const double etaTrkBins[] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};
  const int numEtaTrkBins = sizeof (etaTrkBins) / sizeof (etaTrkBins[0]) - 1;

//const double pTchBins[32] = {0.5, 0.7, 1, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 60, 80, 100};
//const int nPtChBins = 106;
  const double pTchBins[] = {0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 65, 70, 75, 80, 90, 100};
  const int numPtchBins = sizeof (pTchBins) / sizeof (pTchBins[0]) - 1;


  TH1D*  h_truth_matching_prob = new TH1D (Form ("h_truth_matching_prob_%s", sys.Data ()), ";Truth matching prob.;N_{ch}^{rec}", 200, 0, 1);
  h_truth_matching_prob->Sumw2 ();

  TH2D* h2_truth_matched_reco_tracks = new TH2D (Form ("h2_truth_matched_reco_tracks_%s", sys.Data ()), ";#eta;#it{p}_{T} [GeV]", numFinerEtaTrkBins, finerEtaTrkBins, numPtchBins, pTchBins);
  h2_truth_matched_reco_tracks->Sumw2 ();

  TH2D* h2_truth_tracks = new TH2D (Form ("h2_truth_tracks_%s", sys.Data ()), ";#eta;#it{p}_{T} [GeV]", numFinerEtaTrkBins, finerEtaTrkBins, numPtchBins, pTchBins);
  h2_truth_tracks->Sumw2 ();


  TH2D* h2_primary_reco_tracks = new TH2D (Form ("h2_primary_reco_tracks_%s", sys.Data ()), ";#eta;#it{p}_{T} [GeV]", numFinerEtaTrkBins, finerEtaTrkBins, numPtchBins, pTchBins);
  h2_primary_reco_tracks->Sumw2 ();

  TH2D* h2_reco_tracks = new TH2D (Form ("h2_reco_tracks_%s", sys.Data ()), ";#eta;#it{p}_{T} [GeV]", numFinerEtaTrkBins, finerEtaTrkBins, numPtchBins, pTchBins);
  h2_reco_tracks->Sumw2 ();


  TH1D** h_truth_matched_reco_tracks = new TH1D*[numEtaTrkBins];
  TH1D** h_truth_tracks = new TH1D*[numEtaTrkBins];


  TH1D** h_primary_reco_tracks = new TH1D*[numEtaTrkBins];
  TH1D** h_reco_tracks = new TH1D*[numEtaTrkBins];


  for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
    h_truth_matched_reco_tracks[iEta] = new TH1D (Form ("h_truth_matched_reco_tracks_%s_iEta%i", sys.Data (), iEta), ";#it{p}_{T} [GeV]", numPtchBins, pTchBins);
    h_truth_matched_reco_tracks[iEta]->Sumw2 ();

    h_truth_tracks[iEta] = new TH1D (Form ("h_truth_tracks_%s_iEta%i", sys.Data (), iEta), ";#it{p}_{T} [GeV]", numPtchBins, pTchBins);
    h_truth_tracks[iEta]->Sumw2 ();


    h_primary_reco_tracks[iEta] = new TH1D (Form ("h_primary_reco_tracks_%s_iEta%i", sys.Data (), iEta), ";#it{p}_{T} [GeV]", numPtchBins, pTchBins);
    h_primary_reco_tracks[iEta]->Sumw2 ();

    h_reco_tracks[iEta] = new TH1D (Form ("h_reco_tracks_%s_iEta%i", sys.Data (), iEta), ";#it{p}_{T} [GeV]", numPtchBins, pTchBins);
    h_reco_tracks[iEta]->Sumw2 ();
  }


  const int nEvts = tree->GetEntries ();

  // Loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In TrackingPerformance.cxx: Event loop " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    tree->GetEntry (iEvt);

    bool hasPrimary = false;
    bool hasPileup = false;
    float vz = -999;
    for (int iVert = 0; iVert < nvert; iVert++) {
      const bool isPrimary = (vert_type[iVert] == 1);
      hasPrimary = hasPrimary || isPrimary;
      hasPileup = hasPileup || (vert_type[iVert] == 3);
      if (isPrimary)
        vz = vert_z[iVert];
    }
    if (!hasPrimary || hasPileup || fabs (vz) > 150)
    //if (!hasPrimary || fabs (vz) > 150)
      continue;


    const float eventWeight = 1;
    //const float eventWeight = ((isPbPb && !isHijing) ? h_weights->GetBinContent (h_weights->FindFixBin (doNchWeighting ? trk_n : fcalA_et+fcalC_et)) : 1); // weight is 1 for pp

    for (int iTrk = 0; iTrk < trk_n; iTrk++) {
      if (!MeetsTrackCuts (iTrk))
        continue;

      if (trk_pt[iTrk] > 1)
        h_truth_matching_prob->Fill (trk_prob_truth[iTrk]);

      const bool isTruthMatched = (trk_prob_truth[iTrk] > 0.5);

      const bool isFake = !isTruthMatched;
      const bool isSecondary = isTruthMatched && (trk_truth_barcode[iTrk] <= 0 || 200000 <= trk_truth_barcode[iTrk]);

      // primary tracks are non-fake, non-secondary tracks. Strange baryons are excluded too.
      const bool isPrimary = !isFake && !isSecondary && (abs (trk_truth_pdgid[iTrk]) != 3112 && abs (trk_truth_pdgid[iTrk]) != 3222 && abs (trk_truth_pdgid[iTrk]) != 3312 && abs (trk_truth_pdgid[iTrk]) != 3334);

      short iEtaReco = 0;
      while (etaTrkBins[iEtaReco+1] < fabs (trk_eta[iTrk]))
        iEtaReco++;

      if (isPrimary) {
        h2_primary_reco_tracks->Fill (trk_eta[iTrk], trk_pt[iTrk], eventWeight);
        h_primary_reco_tracks[iEtaReco]->Fill (trk_pt[iTrk], eventWeight);

        if (trk_truth_charge[iTrk] != 0 &&
            fabs (trk_truth_eta[iTrk]) <= 2.5 &&
            trk_truth_isHadron[iTrk] &&
            (!DoPionsOnlyVar () || abs (trk_truth_pdgid[iTrk] == 211))) {

          short iEtaTruth = 0;
          while (etaTrkBins[iEtaTruth+1] < fabs (trk_truth_eta[iTrk]))
            iEtaTruth++;

          h2_truth_matched_reco_tracks->Fill (trk_truth_eta[iTrk], trk_truth_pt[iTrk], eventWeight); 
          h_truth_matched_reco_tracks[iEtaTruth]->Fill (trk_truth_pt[iTrk], eventWeight);
          //h2_truth_matched_reco_tracks->Fill (trk_eta[iTrk], trk_pt[iTrk], eventWeight); // do a simple unfolding
          //h_truth_matched_reco_tracks[iEtaReco]->Fill (trk_pt[iTrk], eventWeight);
        }

      }

      h2_reco_tracks->Fill (trk_eta[iTrk], trk_pt[iTrk], eventWeight);
      h_reco_tracks[iEtaReco]->Fill (trk_pt[iTrk], eventWeight);
    }


    for (int iTTrk = 0; iTTrk < truth_trk_n; iTTrk++) {
      if (truth_trk_barcode[iTTrk] <= 0 || 200000 <= truth_trk_barcode[iTTrk] || abs (truth_trk_pdgid[iTTrk]) == 3112 || abs (truth_trk_pdgid[iTTrk]) == 3222 || abs (truth_trk_pdgid[iTTrk]) == 3312 || abs (truth_trk_pdgid[iTTrk]) == 3334)
        continue; // select primary truth particles

      if (truth_trk_charge[iTTrk] != 0 &&
          fabs (truth_trk_eta[iTTrk]) <= 2.5 &&
          truth_trk_isHadron[iTTrk] &&
          (!DoPionsOnlyVar () || abs (truth_trk_pdgid[iTTrk] == 211))) {

        short iEta = 0;
        while (etaTrkBins[iEta+1] < fabs (truth_trk_eta[iTTrk]))
          iEta++;

        //h2_truth_tracks->Fill (truth_trk_eta[iTTrk], truth_trk_pt[iTTrk], 1);
        h2_truth_tracks->Fill (truth_trk_eta[iTTrk], truth_trk_pt[iTTrk], eventWeight);
        h_truth_tracks[iEta]->Fill (truth_trk_pt[iTTrk], eventWeight);
      }
    }

  } // end event loop
  cout << endl << "Info: In TrackingPerformance.cxx: Finished event loop." << endl;


  h_truth_matching_prob->Write ();

  h2_truth_matched_reco_tracks->Write ();
  SaferDelete (&h2_truth_matched_reco_tracks);

  h2_truth_tracks->Write ();
  SaferDelete (&h2_truth_tracks);

  h2_primary_reco_tracks->Write ();
  SaferDelete (&h2_primary_reco_tracks);

  h2_reco_tracks->Write ();
  SaferDelete (&h2_reco_tracks);

  for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
    h_truth_matched_reco_tracks[iEta]->Write ();
    SaferDelete (&h_truth_matched_reco_tracks[iEta]);
    h_truth_tracks[iEta]->Write ();
    SaferDelete (&h_truth_tracks[iEta]);

    h_primary_reco_tracks[iEta]->Write ();
    SaferDelete (&h_primary_reco_tracks[iEta]);
    h_reco_tracks[iEta]->Write ();
    SaferDelete (&h_reco_tracks[iEta]);
  }

  delete [] h_truth_matched_reco_tracks;
  delete [] h_truth_tracks;
  delete [] h_primary_reco_tracks;
  delete [] h_reco_tracks;

  outFile->Write (0, TObject::kOverwrite);
  outFile->Close ();
  SaferDelete (&outFile);

  return true;
}

} // end namespace

#endif
