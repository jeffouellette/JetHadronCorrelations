#ifndef __TrackingPerformance_cxx__
#define __TrackingPerformance_cxx__

#include "TrackingPerformance.h"
#include "Params.h"
#include "CentralityDefs.h"
#include "TreeVariables.h"
#include "Trigger.h"
#include "LocalUtilities.h"

#include <Utilities.h>
#include <ArrayTemplates.h>

#include <TChain.h>
#include <TSystem.h>
#include <TH2D.h>

#include <iostream>

namespace JetHadronCorrelations {

bool TrackingPerformance (const char* directory,
                          const int dataSet,
                          const char* inFileName,
                          const char* eventWeightsFileName) {
 
  std::cout << "Info: In TrackingPerformance.cxx: Entered TrackingPerformance routine." << std::endl;
  std::cout << "Info: In TrackingPerformance.cxx: Printing systematic onfiguration:";

  SetupDirectories ("TrackingPerformance");

  if (IsCollisions ()) {
    std::cout << "Error: In TrackingPerformance.cxx: Trying to calculate tracking performance in data! Quitting." << std::endl;
    return false;
  }

  if (IsDataOverlay ())
    std::cout << "Info: In TrackingPerformance.cxx: Running over data overlay, will check data conditions" << std::endl;
  if (IsHijing ())
    std::cout << "Info: In TrackingPerformance.cxx: Running over Hijing sample" << std::endl;

  // this identifier here corresponds to the name of the output file
  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  std::cout << "Info: In TrackingPerformance.cxx: File Identifier: " << identifier << std::endl;
  std::cout << "Info: In TrackingPerformance.cxx: Saving output to " << rootPath << std::endl;


  // creates a file identifier pattern that we will use to identify the directory containing the input files
  TString fileIdentifier;
  if (TString (inFileName) == "") {
    if (IsCollisions ()) {
      fileIdentifier = std::to_string (dataSet);
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
    std::cout << "Info: In TrackingPerformance.cxx: Chain has " << tree->GetListOfFiles ()->GetEntries () << " files, " << tree->GetEntries () << " entries" << std::endl;
  }



  if (!IsHijing ()) {
    assert (crossSectionPicoBarns > 0);
    assert (mcFilterEfficiency > 0);
    assert (mcNumberEvents > 0);
  }


  TH1D* mcWgtsHist = GetMCFCalWeights ();


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


  tree->SetBranchAddress ("actualInteractionsPerCrossing",  &actualInteractionsPerCrossing);
  tree->SetBranchAddress ("averageInteractionsPerCrossing", &averageInteractionsPerCrossing);


  tree->SetBranchAddress ("mcEventWeights", &mcEventWeights);

  if (IsHijing ()) {
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


  tree->SetBranchAddress ("truth_trk_n",        &truth_trk_n);
  tree->SetBranchAddress ("truth_trk_pt",       &truth_trk_pt);
  tree->SetBranchAddress ("truth_trk_eta",      &truth_trk_eta);
  tree->SetBranchAddress ("truth_trk_phi",      &truth_trk_phi);
  tree->SetBranchAddress ("truth_trk_charge",   &truth_trk_charge);
  tree->SetBranchAddress ("truth_trk_pdgid",    &truth_trk_pdgid);
  tree->SetBranchAddress ("truth_trk_barcode",  &truth_trk_barcode);
  tree->SetBranchAddress ("truth_trk_isHadron", &truth_trk_isHadron);

  tree->SetBranchAddress ("ntrk",                   &trk_n);
  tree->SetBranchAddress ("trk_pt",                 &trk_pt);
  tree->SetBranchAddress ("trk_eta",                &trk_eta);
  tree->SetBranchAddress ("trk_phi",                &trk_phi);
  tree->SetBranchAddress ("trk_charge",             &trk_charge);
  tree->SetBranchAddress ("trk_TightPrimary",       &trk_TightPrimary);
  tree->SetBranchAddress ("trk_HITight",            &trk_HITight);
  tree->SetBranchAddress ("trk_HILoose",            &trk_HILoose);
  tree->SetBranchAddress ("trk_d0",                 &trk_d0);
  tree->SetBranchAddress ("trk_d0sig",              &trk_d0sig);
  tree->SetBranchAddress ("trk_z0",                 &trk_z0);
  tree->SetBranchAddress ("trk_z0sig",              &trk_z0sig);
  tree->SetBranchAddress ("trk_theta",              &trk_theta);
  tree->SetBranchAddress ("trk_vz",                 &trk_vz);

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


  tree->SetBranchAddress ("akt4_hi_jet_n",            &akt4_hi_jet_n);
  tree->SetBranchAddress ("akt4_hi_jet_pt_etajes",    &akt4_hi_jet_pt_etajes);
  tree->SetBranchAddress ("akt4_hi_jet_pt_xcalib",    &akt4_hi_jet_pt_xcalib);
  tree->SetBranchAddress ("akt4_hi_jet_eta_etajes",   &akt4_hi_jet_eta_etajes);
  tree->SetBranchAddress ("akt4_hi_jet_eta_xcalib",   &akt4_hi_jet_eta_xcalib);
  tree->SetBranchAddress ("akt4_hi_jet_phi",          &akt4_hi_jet_phi);
  tree->SetBranchAddress ("akt4_hi_jet_e_etajes",     &akt4_hi_jet_e_etajes);
  tree->SetBranchAddress ("akt4_hi_jet_e_xcalib",     &akt4_hi_jet_e_xcalib);
  tree->SetBranchAddress ("akt4_hi_jet_timing",       &akt4_hi_jet_timing);


  tree->SetBranchAddress ("akt4_truth_jet_n",     &akt4_truth_jet_n);
  tree->SetBranchAddress ("akt4_truth_jet_pt",    &akt4_truth_jet_pt);
  tree->SetBranchAddress ("akt4_truth_jet_eta",   &akt4_truth_jet_eta);
  tree->SetBranchAddress ("akt4_truth_jet_phi",   &akt4_truth_jet_phi);
  tree->SetBranchAddress ("akt4_truth_jet_e",     &akt4_truth_jet_e);


  std::cout << "Info : In TrackingPerformance.cxx: Saving histograms to " << Form ("%s/%s.root", rootPath.Data (), identifier.Data ()) << std::endl;
  TFile* outFile = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");

  TString sys;
  if (Ispp ())        sys = "pp";
  else if (IspPb ())  sys = "pPb";
  else                sys = "???";

  const int nFinerEtaTrkBins = 40;
  const double* finerEtaTrkBins = linspace (-2.5, 2.5, nFinerEtaTrkBins);

  const double pTchBins[] = {0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130}; // new binning, 113 elements or 112 bins
  //const double pTchBins[] = {0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 50, 55, 60, 70, 80, 90, 100}; // old binning
  const int nPtchBins = sizeof (pTchBins) / sizeof (pTchBins[0]) - 1;

                    //pi+, k+,  p+,   e-, mu-, sigma+, sigma-, xi-,  omega-, everyone
  const int PIDs[] = {211, 321, 2212, 11, 13,  3222,   3112,   3312, 3334,   0};
  const int nPIDs = sizeof (PIDs) / sizeof (PIDs[0]);


  // Initialize a bunch of histograms --
  TH1D** h_truth_matching_prob = Get1DArray <TH1D*> (trackWPs.size ());

  TH2D**** h2_truth_matched_primary_tracks = Get3DArray <TH2D*> (nPIDs, trackWPs.size (), nMultBins);
  TH2D***  h2_truth_tracks                 = Get2DArray <TH2D*> (nPIDs, nMultBins);
  TH2D***  h2_truth_tracks_wgt2            = Get2DArray <TH2D*> (nPIDs, nMultBins);

  TH2D*** h2_fake_tracks            = Get2DArray <TH2D*> (trackWPs.size (), nDRBins);
  TH2D*** h2_secondary_tracks       = Get2DArray <TH2D*> (trackWPs.size (), nDRBins);
  TH2D*** h2_strange_tracks         = Get2DArray <TH2D*> (trackWPs.size (), nDRBins);
  TH2D*** h2_primary_tracks         = Get2DArray <TH2D*> (trackWPs.size (), nDRBins);
  TH2D*** h2_reco_tracks            = Get2DArray <TH2D*> (trackWPs.size (), nDRBins);
  TH2D*** h2_reco_tracks_wgt2       = Get2DArray <TH2D*> (trackWPs.size (), nDRBins);


  TH1D***** h_truth_matched_primary_tracks = Get4DArray <TH1D*> (nPIDs, trackWPs.size (), nMultBins, nEtaTrkBins);
  TH1D****  h_truth_tracks                 = Get3DArray <TH1D*> (nPIDs, nMultBins, nEtaTrkBins);
  TH1D****  h_truth_tracks_wgt2            = Get3DArray <TH1D*> (nPIDs, nMultBins, nEtaTrkBins);

  TH1D**** h_fake_tracks            = Get3DArray <TH1D*> (trackWPs.size (), nDRBins, nEtaTrkBins);
  TH1D**** h_secondary_tracks       = Get3DArray <TH1D*> (trackWPs.size (), nDRBins, nEtaTrkBins);
  TH1D**** h_strange_tracks         = Get3DArray <TH1D*> (trackWPs.size (), nDRBins, nEtaTrkBins);
  TH1D**** h_primary_tracks         = Get3DArray <TH1D*> (trackWPs.size (), nDRBins, nEtaTrkBins);
  TH1D**** h_reco_tracks            = Get3DArray <TH1D*> (trackWPs.size (), nDRBins, nEtaTrkBins);
  TH1D**** h_reco_tracks_wgt2       = Get3DArray <TH1D*> (trackWPs.size (), nDRBins, nEtaTrkBins);


  for (int iPID = 0; iPID < nPIDs; iPID++) {

    for (int iMult = 0; iMult < nMultBins; iMult++) {

      h2_truth_tracks[iPID][iMult] = new TH2D (Form ("h2_truth_tracks_%s_PID%i_iMult%i", sys.Data (), PIDs[iPID], iMult), ";#eta;#it{p}_{T} [GeV]", nFinerEtaTrkBins, finerEtaTrkBins, nPtchBins, pTchBins);
      h2_truth_tracks[iPID][iMult]->Sumw2 ();

      h2_truth_tracks_wgt2[iPID][iMult] = new TH2D (Form ("h2_truth_tracks_wgt2_%s_PID%i_iMult%i", sys.Data (), PIDs[iPID], iMult), ";#eta;#it{p}_{T} [GeV]", nFinerEtaTrkBins, finerEtaTrkBins, nPtchBins, pTchBins);
      h2_truth_tracks_wgt2[iPID][iMult]->Sumw2 ();

      for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {

        h_truth_tracks[iPID][iMult][iEta] = new TH1D (Form ("h_truth_tracks_%s_PID%i_iMult%i_iEta%i", sys.Data (), PIDs[iPID], iMult, iEta), ";#it{p}_{T} [GeV]", nPtchBins, pTchBins);
        h_truth_tracks[iPID][iMult][iEta]->Sumw2 ();

        h_truth_tracks_wgt2[iPID][iMult][iEta] = new TH1D (Form ("h_truth_tracks_wgt2_%s_PID%i_iMult%i_iEta%i", sys.Data (), PIDs[iPID], iMult, iEta), ";#it{p}_{T} [GeV]", nPtchBins, pTchBins);
        h_truth_tracks_wgt2[iPID][iMult][iEta]->Sumw2 ();

      } // end loop over iEta

    } // end loop over iMult


    for (int iWP = 0; iWP < (int)trackWPs.size (); iWP++) {

      for (int iMult = 0; iMult < nMultBins; iMult++) {

        h2_truth_matched_primary_tracks[iPID][iWP][iMult] = new TH2D (Form ("h2_truth_matched_primary_tracks_%s_PID%i_%s_iMult%i", sys.Data (), PIDs[iPID], trackWPNames[iWP].c_str (), iMult), ";#eta;#it{p}_{T} [GeV]", nFinerEtaTrkBins, finerEtaTrkBins, nPtchBins, pTchBins);
        h2_truth_matched_primary_tracks[iPID][iWP][iMult]->Sumw2 ();

        for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {

          h_truth_matched_primary_tracks[iPID][iWP][iMult][iEta] = new TH1D (Form ("h_truth_matched_primary_tracks_%s_PID%i_%s_iMult%i_iEta%i", sys.Data (), PIDs[iPID], trackWPNames[iWP].c_str (), iMult, iEta), ";#it{p}_{T} [GeV]", nPtchBins, pTchBins);
          h_truth_matched_primary_tracks[iPID][iWP][iMult][iEta]->Sumw2 ();

        } // end loop over iEta

      } // end loop over iMult

    } // end loop over iWP

  } // end loop over iPID


  for (int iWP = 0; iWP < (int)trackWPs.size (); iWP++) {

    h_truth_matching_prob[iWP] = new TH1D (Form ("h_truth_matching_prob_%s%s_%s", sys.Data (), IsHijing () ? "_Hijing" : "", trackWPNames[iWP].c_str ()), ";Truth matching prob.;N_{ch}^{rec}", 200, 0, 1);
    h_truth_matching_prob[iWP]->Sumw2 ();

    for (int iDR = 0; iDR < nDRBins; iDR++) {

      h2_fake_tracks[iWP][iDR] = new TH2D (Form ("h2_fake_tracks_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR), ";#eta;#it{p}_{T} [GeV]", nFinerEtaTrkBins, finerEtaTrkBins, nPtchBins, pTchBins);
      h2_fake_tracks[iWP][iDR]->Sumw2 ();
      h2_secondary_tracks[iWP][iDR] = new TH2D (Form ("h2_secondary_tracks_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR), ";#eta;#it{p}_{T} [GeV]", nFinerEtaTrkBins, finerEtaTrkBins, nPtchBins, pTchBins);
      h2_secondary_tracks[iWP][iDR]->Sumw2 ();
      h2_strange_tracks[iWP][iDR] = new TH2D (Form ("h2_strange_tracks_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR), ";#eta;#it{p}_{T} [GeV]", nFinerEtaTrkBins, finerEtaTrkBins, nPtchBins, pTchBins);
      h2_strange_tracks[iWP][iDR]->Sumw2 ();
      h2_primary_tracks[iWP][iDR] = new TH2D (Form ("h2_primary_tracks_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR), ";#eta;#it{p}_{T} [GeV]", nFinerEtaTrkBins, finerEtaTrkBins, nPtchBins, pTchBins);
      h2_primary_tracks[iWP][iDR]->Sumw2 ();
      h2_reco_tracks[iWP][iDR] = new TH2D (Form ("h2_reco_tracks_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR), ";#eta;#it{p}_{T} [GeV]", nFinerEtaTrkBins, finerEtaTrkBins, nPtchBins, pTchBins);
      h2_reco_tracks[iWP][iDR]->Sumw2 ();

      h2_reco_tracks_wgt2[iWP][iDR] = new TH2D (Form ("h2_reco_tracks_wgt2_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR), ";#eta;#it{p}_{T} [GeV]", nFinerEtaTrkBins, finerEtaTrkBins, nPtchBins, pTchBins);
      h2_reco_tracks_wgt2[iWP][iDR]->Sumw2 ();

      for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {

        h_fake_tracks[iWP][iDR][iEta] = new TH1D (Form ("h_fake_tracks_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta), ";#it{p}_{T} [GeV]", nPtchBins, pTchBins);
        h_fake_tracks[iWP][iDR][iEta]->Sumw2 ();
        h_secondary_tracks[iWP][iDR][iEta] = new TH1D (Form ("h_secondary_tracks_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta), ";#it{p}_{T} [GeV]", nPtchBins, pTchBins);
        h_secondary_tracks[iWP][iDR][iEta]->Sumw2 ();
        h_strange_tracks[iWP][iDR][iEta] = new TH1D (Form ("h_strange_tracks_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta), ";#it{p}_{T} [GeV]", nPtchBins, pTchBins);
        h_strange_tracks[iWP][iDR][iEta]->Sumw2 ();
        h_primary_tracks[iWP][iDR][iEta] = new TH1D (Form ("h_primary_tracks_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta), ";#it{p}_{T} [GeV]", nPtchBins, pTchBins);
        h_primary_tracks[iWP][iDR][iEta]->Sumw2 ();
        h_reco_tracks[iWP][iDR][iEta] = new TH1D (Form ("h_reco_tracks_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta), ";#it{p}_{T} [GeV]", nPtchBins, pTchBins);
        h_reco_tracks[iWP][iDR][iEta]->Sumw2 ();

        h_reco_tracks_wgt2[iWP][iDR][iEta] = new TH1D (Form ("h_reco_tracks_wgt2_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta), ";#it{p}_{T} [GeV]", nPtchBins, pTchBins);
        h_reco_tracks_wgt2[iWP][iDR][iEta]->Sumw2 ();

      } // end loop over iEta

    } // end loop over iDR

  } // end loop over iWP


  const JetRadius r0p4 = JetRadius::R0p4;
  const int nEvts = tree->GetEntries ();

  // Loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      std::cout << "Info: In TrackingPerformance.cxx: Event loop " << iEvt / (nEvts / 100) << "\% done...\r" << std::flush;
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


    // Get multiplicity bin for event
    int iMult = 0;
    while (iMult < nMultBins-1 && multBins[iMult+1] < trk_n)
      iMult++;
    if (nMultBins < iMult)
      continue;


    const double eventWeight = (UseMCFCalWeights () ? mcWgtsHist->GetBinContent (mcWgtsHist->FindBin (fcalA_et)) : 1);
    //const double eventWeight = (IsHijing () ? 1 : mcEventWeights->at (0) * crossSectionPicoBarns * mcFilterEfficiency * GetJetLuminosity () / mcNumberEvents); // sigma * f * L_int

    // minimum and maximum barcodes, 0 < barcode < 10000 in Pythia and 10000 < barcode < 200000 in Hijing
    const int minBarcode = 0;
    const int maxBarcode = 200000;

    const int njet = GetAktHIJetN (r0p4);

    for (int iTrk = 0; iTrk < trk_n; iTrk++) {

      float mindr = M_PI;
      for (int iJet = 0; iJet < njet; iJet++) {
        if (GetAktHIJetPt (iJet, r0p4) > 15)
          mindr = std::fmin (mindr, DeltaR (trk_eta[iTrk], GetAktHIJetEta (iJet, r0p4), trk_phi[iTrk], GetAktHIJetPhi (iJet, r0p4)));
      }
      int iDR = (mindr > 0.4 ? 0 : (mindr > 0.2 ? 1 : 2));

      short iEtaReco = 0;
      while (etaTrkBins[iEtaReco+1] < std::fabs (trk_eta[iTrk]))
        iEtaReco++;

      // has trigger jet?
      bool hasTrigJet = false;
      for (int iJet = 0; !hasTrigJet && iJet < GetAktHIJetN (r0p4); iJet++) {
        if (!MeetsJetAcceptanceCuts (iJet, r0p4))
          continue; // jet eta/phi & timing cuts

        if (!MeetsJetPtCut (GetAktHIJetPt  (iJet, r0p4)))
          continue; // jet pT cuts

        hasTrigJet = true;
      }
      if (!hasTrigJet)
        continue;


      for (int iWP = 0; iWP < (int)trackWPs.size (); iWP++) {

        if (!MeetsTrackCuts (iTrk, iWP))
          continue;

        if (trk_pt[iTrk] > 1)
          h_truth_matching_prob[iWP]->Fill (trk_prob_truth[iTrk]);

        const bool isTruthMatched = (trk_prob_truth[iTrk] > 0.5);

        const bool isFake = !isTruthMatched;

        const bool isSecondary = isTruthMatched && (trk_truth_barcode[iTrk] <= minBarcode || maxBarcode <= trk_truth_barcode[iTrk]);

        const bool isStrangeBaryon = !isFake && !isSecondary && (std::abs (trk_truth_pdgid[iTrk]) == 3112 || std::abs (trk_truth_pdgid[iTrk]) == 3222 || std::abs (trk_truth_pdgid[iTrk]) == 3312 || std::abs (trk_truth_pdgid[iTrk]) == 3334);

        // primary tracks are non-fake, non-secondary tracks.
        const bool isPrimary = !isFake && !isSecondary;// && !isStrangeBaryon;

        if (!IsDataOverlay ()) {

          if (isFake) {
            h2_fake_tracks[iWP][iDR]->Fill (trk_eta[iTrk], trk_pt[iTrk]);
            h_fake_tracks[iWP][iDR][iEtaReco]->Fill (trk_pt[iTrk]);
            h2_fake_tracks[iWP][nDRBins-1]->Fill (trk_eta[iTrk], trk_pt[iTrk]);
            h_fake_tracks[iWP][nDRBins-1][iEtaReco]->Fill (trk_pt[iTrk]);
          }
          else if (isSecondary) {
            h2_secondary_tracks[iWP][iDR]->Fill (trk_eta[iTrk], trk_pt[iTrk], eventWeight);
            h_secondary_tracks[iWP][iDR][iEtaReco]->Fill (trk_pt[iTrk], eventWeight);
            h2_secondary_tracks[iWP][nDRBins-1]->Fill (trk_eta[iTrk], trk_pt[iTrk], eventWeight);
            h_secondary_tracks[iWP][nDRBins-1][iEtaReco]->Fill (trk_pt[iTrk], eventWeight);
          }
          else if (isStrangeBaryon) {
            h2_strange_tracks[iWP][iDR]->Fill (trk_eta[iTrk], trk_pt[iTrk]);
            h_strange_tracks[iWP][iDR][iEtaReco]->Fill (trk_pt[iTrk]);
            h2_strange_tracks[iWP][nDRBins-1]->Fill (trk_eta[iTrk], trk_pt[iTrk]);
            h_strange_tracks[iWP][nDRBins-1][iEtaReco]->Fill (trk_pt[iTrk]);
          }
          else if (isPrimary) {
            h2_primary_tracks[iWP][iDR]->Fill (trk_eta[iTrk], trk_pt[iTrk]);
            h_primary_tracks[iWP][iDR][iEtaReco]->Fill (trk_pt[iTrk]);
            h2_primary_tracks[iWP][nDRBins-1]->Fill (trk_eta[iTrk], trk_pt[iTrk]);
            h_primary_tracks[iWP][nDRBins-1][iEtaReco]->Fill (trk_pt[iTrk]);
          }

          h2_reco_tracks[iWP][iDR]->Fill (trk_eta[iTrk], trk_pt[iTrk]);
          h_reco_tracks[iWP][iDR][iEtaReco]->Fill (trk_pt[iTrk]);
          h2_reco_tracks[iWP][nDRBins-1]->Fill (trk_eta[iTrk], trk_pt[iTrk]);
          h_reco_tracks[iWP][nDRBins-1][iEtaReco]->Fill (trk_pt[iTrk]);

          h2_reco_tracks_wgt2[iWP][iDR]->Fill (trk_eta[iTrk], trk_pt[iTrk], eventWeight*eventWeight);
          h_reco_tracks_wgt2[iWP][iDR][iEtaReco]->Fill (trk_pt[iTrk], eventWeight*eventWeight);
          h2_reco_tracks_wgt2[iWP][nDRBins-1]->Fill (trk_eta[iTrk], trk_pt[iTrk], eventWeight*eventWeight);
          h_reco_tracks_wgt2[iWP][nDRBins-1][iEtaReco]->Fill (trk_pt[iTrk], eventWeight*eventWeight);

        }


        if (!IsHijing () && (isPrimary || isStrangeBaryon)) {
        //if (isPrimary || isStrangeBaryon) {

          if (trk_truth_charge[iTrk] == 0 || std::fabs (trk_truth_eta[iTrk]) > 2.5)
            continue; // truth-level acceptance cut

          short iEtaTruth = 0;
          while (etaTrkBins[iEtaTruth+1] < std::fabs (trk_truth_eta[iTrk]))
            iEtaTruth++;

          short iPID = 0;
          while (iPID < nPIDs-1 && PIDs[iPID] != trk_truth_pdgid[iTrk])
            iPID++;

          if (0 <= iPID && iPID < nPIDs-1) {
            h2_truth_matched_primary_tracks[iPID][iWP][iMult]->Fill (trk_truth_eta[iTrk], trk_truth_pt[iTrk], eventWeight); 
            h_truth_matched_primary_tracks[iPID][iWP][iMult][iEtaTruth]->Fill (trk_truth_pt[iTrk], eventWeight);
            h2_truth_matched_primary_tracks[iPID][iWP][nMultBins-1]->Fill (trk_truth_eta[iTrk], trk_truth_pt[iTrk], eventWeight); 
            h_truth_matched_primary_tracks[iPID][iWP][nMultBins-1][iEtaTruth]->Fill (trk_truth_pt[iTrk], eventWeight);
          }

          if (isPrimary && trk_truth_isHadron[iTrk]) {
            h2_truth_matched_primary_tracks[nPIDs-1][iWP][iMult]->Fill (trk_truth_eta[iTrk], trk_truth_pt[iTrk], eventWeight); 
            h_truth_matched_primary_tracks[nPIDs-1][iWP][iMult][iEtaTruth]->Fill (trk_truth_pt[iTrk], eventWeight);
            h2_truth_matched_primary_tracks[nPIDs-1][iWP][nMultBins-1]->Fill (trk_truth_eta[iTrk], trk_truth_pt[iTrk], eventWeight); 
            h_truth_matched_primary_tracks[nPIDs-1][iWP][nMultBins-1][iEtaTruth]->Fill (trk_truth_pt[iTrk], eventWeight);
          }

        }

      } // end loop over WPs

    } // end loop over reco tracks


    if (!IsHijing ()) {

      for (int iTTrk = 0; iTTrk < truth_trk_n; iTTrk++) {

        if (truth_trk_charge[iTTrk] == 0 || std::fabs (truth_trk_eta[iTTrk]) > 2.5)
          continue; // truth-level acceptance cut

        const bool isSecondary = (truth_trk_barcode[iTTrk] <= minBarcode || maxBarcode <= truth_trk_barcode[iTTrk]);
        if (isSecondary)
          continue; // don't work with secondaries

        //const bool isStrangeBaryon = !isSecondary && (std::abs (truth_trk_pdgid[iTTrk]) == 3112 || std::abs (truth_trk_pdgid[iTTrk]) == 3222 || std::abs (truth_trk_pdgid[iTTrk]) == 3312 || std::abs (truth_trk_pdgid[iTTrk]) == 3334);

        const bool isPrimary = !isSecondary;// && !isStrangeBaryon;

        short iEta = 0;
        while (etaTrkBins[iEta+1] < std::fabs (truth_trk_eta[iTTrk]))
          iEta++;

        short iPID = 0;
        while (iPID < nPIDs-1 && PIDs[iPID] != truth_trk_pdgid[iTTrk])
          iPID++;

        if (0 <= iPID && iPID < nPIDs-1) {
          h2_truth_tracks[iPID][iMult]->Fill (truth_trk_eta[iTTrk], truth_trk_pt[iTTrk], eventWeight); 
          h_truth_tracks[iPID][iMult][iEta]->Fill (truth_trk_pt[iTTrk], eventWeight);
          h2_truth_tracks[iPID][nMultBins-1]->Fill (truth_trk_eta[iTTrk], truth_trk_pt[iTTrk], eventWeight); 
          h_truth_tracks[iPID][nMultBins-1][iEta]->Fill (truth_trk_pt[iTTrk], eventWeight);

          h2_truth_tracks_wgt2[iPID][iMult]->Fill (truth_trk_eta[iTTrk], truth_trk_pt[iTTrk], eventWeight*eventWeight);
          h_truth_tracks_wgt2[iPID][iMult][iEta]->Fill (truth_trk_pt[iTTrk], eventWeight*eventWeight);
          h2_truth_tracks_wgt2[iPID][nMultBins-1]->Fill (truth_trk_eta[iTTrk], truth_trk_pt[iTTrk], eventWeight*eventWeight);
          h_truth_tracks_wgt2[iPID][nMultBins-1][iEta]->Fill (truth_trk_pt[iTTrk], eventWeight*eventWeight);
        }

        if (isPrimary && truth_trk_isHadron[iTTrk]) {
          h2_truth_tracks[nPIDs-1][iMult]->Fill (truth_trk_eta[iTTrk], truth_trk_pt[iTTrk], eventWeight);
          h_truth_tracks[nPIDs-1][iMult][iEta]->Fill (truth_trk_pt[iTTrk], eventWeight);
          h2_truth_tracks[nPIDs-1][nMultBins-1]->Fill (truth_trk_eta[iTTrk], truth_trk_pt[iTTrk], eventWeight);
          h_truth_tracks[nPIDs-1][nMultBins-1][iEta]->Fill (truth_trk_pt[iTTrk], eventWeight);

          h2_truth_tracks_wgt2[nPIDs-1][nMultBins-1]->Fill (truth_trk_eta[iTTrk], truth_trk_pt[iTTrk], eventWeight*eventWeight);
          h_truth_tracks_wgt2[nPIDs-1][nMultBins-1][iEta]->Fill (truth_trk_pt[iTTrk], eventWeight*eventWeight);
        }

      } // end loop over truth tracks

    }

  } // end event loop
  std::cout << std::endl << "Info: In TrackingPerformance.cxx: Finished event loop." << std::endl;


  SaferDelete (&tree);


  outFile->cd ();


  for (int iPID = 0; iPID < nPIDs; iPID++) {

    for (int iMult = 0; iMult < nMultBins; iMult++) {

      h2_truth_tracks[iPID][iMult]->Write ();
      SaferDelete (&h2_truth_tracks[iPID][iMult]);

      h2_truth_tracks_wgt2[iPID][iMult]->Write ();
      SaferDelete (&h2_truth_tracks_wgt2[iPID][iMult]);

      for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {

        h_truth_tracks[iPID][iMult][iEta]->Write ();
        SaferDelete (&h_truth_tracks[iPID][iMult][iEta]);

        h_truth_tracks_wgt2[iPID][iMult][iEta]->Write ();
        SaferDelete (&h_truth_tracks_wgt2[iPID][iMult][iEta]);

      } // end loop over iEta

    } // end loop over iMult

    for (int iWP = 0; iWP < (int)trackWPs.size (); iWP++) {

      for (int iMult = 0; iMult < nMultBins; iMult++) {

        h2_truth_matched_primary_tracks[iPID][iWP][iMult]->Write ();
        SaferDelete (&h2_truth_matched_primary_tracks[iPID][iWP][iMult]);

        for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {

          h_truth_matched_primary_tracks[iPID][iWP][iMult][iEta]->Write ();
          SaferDelete (&h_truth_matched_primary_tracks[iPID][iWP][iMult][iEta]);

        } // end loop over iEta

      } // end loop over iMult

    } // end loop over iWP

  } // end loop over iPID


  for (int iWP = 0; iWP < (int)trackWPs.size (); iWP++) {

    h_truth_matching_prob[iWP]->Write ();
    SaferDelete (&h_truth_matching_prob[iWP]);

    //for (int iMult = 0; iMult < nMultBins; iMult++) {
    for (int iDR = 0; iDR < nDRBins; iDR++) {

      h2_fake_tracks[iWP][iDR]->Write ();
      SaferDelete (&h2_fake_tracks[iWP][iDR]);
      h2_secondary_tracks[iWP][iDR]->Write ();
      SaferDelete (&h2_secondary_tracks[iWP][iDR]);
      h2_strange_tracks[iWP][iDR]->Write ();
      SaferDelete (&h2_strange_tracks[iWP][iDR]);
      h2_primary_tracks[iWP][iDR]->Write ();
      SaferDelete (&h2_primary_tracks[iWP][iDR]);

      h2_reco_tracks[iWP][iDR]->Write ();
      SaferDelete (&h2_reco_tracks[iWP][iDR]);

      h2_reco_tracks_wgt2[iWP][iDR]->Write ();
      SaferDelete (&h2_reco_tracks_wgt2[iWP][iDR]);

      for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {

        h_fake_tracks[iWP][iDR][iEta]->Write ();
        SaferDelete (&h_fake_tracks[iWP][iDR][iEta]);
        h_secondary_tracks[iWP][iDR][iEta]->Write ();
        SaferDelete (&h_secondary_tracks[iWP][iDR][iEta]);
        h_strange_tracks[iWP][iDR][iEta]->Write ();
        SaferDelete (&h_strange_tracks[iWP][iDR][iEta]);
        h_primary_tracks[iWP][iDR][iEta]->Write ();
        SaferDelete (&h_primary_tracks[iWP][iDR][iEta]);

        h_reco_tracks[iWP][iDR][iEta]->Write ();
        SaferDelete (&h_reco_tracks[iWP][iDR][iEta]);

        h_reco_tracks_wgt2[iWP][iDR][iEta]->Write ();
        SaferDelete (&h_reco_tracks_wgt2[iWP][iDR][iEta]);

      } // end loop over iEta

    } // end loop over iDR

  } // end loop over iWP

  outFile->Write (0, TObject::kOverwrite);
  outFile->Close ();
  SaferDelete (&outFile);

  return true;
}

} // end namespace

#endif
