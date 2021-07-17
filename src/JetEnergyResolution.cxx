#ifndef __JetEnergyResolution_cxx__
#define __JetEnergyResolution_cxx__

#include "JetEnergyResolution.h"
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
#include <TF1.h>

#include <iostream>

using namespace std;

namespace JetHadronCorrelations {

bool JetEnergyResolution (const char* directory,
                          const int dataSet,
                          const char* inFileName,
                          const char* eventWeightsFileName) {
 
  cout << "Info: In JetEnergyResolution.cxx: Entered JetEnergyResolution routine." << endl;
  cout << "Info: In JetEnergyResolution.cxx: Printing systematic onfiguration:";

  SetupDirectories ("JetEnergyResolution");

  if (IsCollisions ()) {
    cout << "Error: In JetEnergyResolution.cxx: Trying to calculate jet performance in data! Quitting." << endl;
    return false;
  }

  if (IsDataOverlay ())
    cout << "Info: In JetEnergyResolution.cxx: Running over data overlay, will check data conditions" << endl;
  if (IsHijing ())
    cout << "Info: In JetEnergyResolution.cxx: Running over Hijing sample" << endl;

  // this identifier here corresponds to the name of the output file
  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  cout << "Info: In JetEnergyResolution.cxx: File Identifier: " << identifier << endl;
  cout << "Info: In JetEnergyResolution.cxx: Saving output to " << rootPath << endl;


  // creates a file identifier pattern that we will use to identify the directory containing the input files
  TString fileIdentifier;
  if (TString (inFileName) == "") {
    if (IsCollisions ()) {
      fileIdentifier = to_string (dataSet);
    }
    else {
      std::cout << "Error: In JetEnergyResolution.cxx: Cannot identify this MC file! Quitting." << std::endl;
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
    std::cout << "Info: In JetEnergyResolution.cxx: Chain has " << tree->GetListOfFiles ()->GetEntries () << " files, " << tree->GetEntries () << " entries" << std::endl;
  }


  if (!IsHijing ()) {
    assert (crossSectionPicoBarns > 0);
    assert (mcFilterEfficiency > 0);
    assert (mcNumberEvents > 0);
  }


  // variables for filtering MC truth
  double truth_jet_min_pt = 0, truth_jet_max_pt = DBL_MAX;
  if (TString (inFileName).Contains ("JZ0")) {
    truth_jet_min_pt = 0;
    truth_jet_max_pt = 20;
  }
  else if (TString (inFileName).Contains ("JZ1")) {
    truth_jet_min_pt = 20;
    truth_jet_max_pt = 60;
  }
  else if (TString (inFileName).Contains ("JZ2")) {
    truth_jet_min_pt = 60;
    truth_jet_max_pt = 160;
  }
  else if (TString (inFileName).Contains ("JZ3")) {
    truth_jet_min_pt = 160;
    truth_jet_max_pt = 400;
  }
  else if (TString (inFileName).Contains ("JZ4")) {
    truth_jet_min_pt = 400;
    truth_jet_max_pt = 800;
  }
  else if (TString (inFileName).Contains ("JZ5")) {
    truth_jet_min_pt = 800;
    truth_jet_max_pt = 1300;
  }
  if (truth_jet_min_pt != 0)
    std::cout << "Checking for leading truth jet with pT > " << truth_jet_min_pt << std::endl;
  if (truth_jet_max_pt != DBL_MAX)
    std::cout << "Checking for leading truth jet with pT < " << truth_jet_max_pt << std::endl;


  std::cout << "Anti-kT R=0.2 jets truth matched to within dR < " << akt2_TruthMatchMaxDR << std::endl;
  std::cout << "Anti-kT R=0.4 jets truth matched to within dR < " << akt4_TruthMatchMaxDR << std::endl;


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


  tree->SetBranchAddress ("akt2_truth_jet_n",     &akt2_truth_jet_n);
  tree->SetBranchAddress ("akt2_truth_jet_pt",    &akt2_truth_jet_pt);
  tree->SetBranchAddress ("akt2_truth_jet_eta",   &akt2_truth_jet_eta);
  tree->SetBranchAddress ("akt2_truth_jet_phi",   &akt2_truth_jet_phi);
  tree->SetBranchAddress ("akt2_truth_jet_e",     &akt2_truth_jet_e);


  tree->SetBranchAddress ("akt2_hi_jet_n",            &akt2_hi_jet_n);
  tree->SetBranchAddress ("akt2_hi_jet_pt_precalib",  &akt2_hi_jet_pt_precalib);
  tree->SetBranchAddress ("akt2_hi_jet_pt_etajes",    &akt2_hi_jet_pt_etajes);
  tree->SetBranchAddress ("akt2_hi_jet_pt_xcalib",    &akt2_hi_jet_pt_xcalib);
  tree->SetBranchAddress ("akt2_hi_jet_eta_precalib", &akt2_hi_jet_eta_precalib);
  tree->SetBranchAddress ("akt2_hi_jet_eta_etajes",   &akt2_hi_jet_eta_etajes);
  tree->SetBranchAddress ("akt2_hi_jet_eta_xcalib",   &akt2_hi_jet_eta_xcalib);
  tree->SetBranchAddress ("akt2_hi_jet_phi",          &akt2_hi_jet_phi);
  tree->SetBranchAddress ("akt2_hi_jet_e_precalib",   &akt2_hi_jet_e_precalib);
  tree->SetBranchAddress ("akt2_hi_jet_e_etajes",     &akt2_hi_jet_e_etajes);
  tree->SetBranchAddress ("akt2_hi_jet_e_xcalib",     &akt2_hi_jet_e_xcalib);
  tree->SetBranchAddress ("akt2_hi_jet_sub_et",       &akt2_hi_jet_sub_et);
  tree->SetBranchAddress ("akt2_hi_jet_sub_e",        &akt2_hi_jet_sub_e);


  tree->SetBranchAddress ("akt4_truth_jet_n",     &akt4_truth_jet_n);
  tree->SetBranchAddress ("akt4_truth_jet_pt",    &akt4_truth_jet_pt);
  tree->SetBranchAddress ("akt4_truth_jet_eta",   &akt4_truth_jet_eta);
  tree->SetBranchAddress ("akt4_truth_jet_phi",   &akt4_truth_jet_phi);
  tree->SetBranchAddress ("akt4_truth_jet_e",     &akt4_truth_jet_e);


  tree->SetBranchAddress ("akt4_hi_jet_n",            &akt4_hi_jet_n);
  tree->SetBranchAddress ("akt4_hi_jet_pt_precalib",  &akt4_hi_jet_pt_precalib);
  tree->SetBranchAddress ("akt4_hi_jet_pt_etajes",    &akt4_hi_jet_pt_etajes);
  tree->SetBranchAddress ("akt4_hi_jet_pt_xcalib",    &akt4_hi_jet_pt_xcalib);
  tree->SetBranchAddress ("akt4_hi_jet_eta_precalib", &akt4_hi_jet_eta_precalib);
  tree->SetBranchAddress ("akt4_hi_jet_eta_etajes",   &akt4_hi_jet_eta_etajes);
  tree->SetBranchAddress ("akt4_hi_jet_eta_xcalib",   &akt4_hi_jet_eta_xcalib);
  tree->SetBranchAddress ("akt4_hi_jet_phi",          &akt4_hi_jet_phi);
  tree->SetBranchAddress ("akt4_hi_jet_e_precalib",   &akt4_hi_jet_e_precalib);
  tree->SetBranchAddress ("akt4_hi_jet_e_etajes",     &akt4_hi_jet_e_etajes);
  tree->SetBranchAddress ("akt4_hi_jet_e_xcalib",     &akt4_hi_jet_e_xcalib);
  tree->SetBranchAddress ("akt4_hi_jet_sub_et",       &akt4_hi_jet_sub_et);
  tree->SetBranchAddress ("akt4_hi_jet_sub_e",        &akt4_hi_jet_sub_e);


  cout << "Info : In JetEnergyResolution.cxx: Saving histograms to " << Form ("%s/%s.root", rootPath.Data (), identifier.Data ()) << endl;
  TFile* outFile = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");

  TString sys;
  if (Ispp ())        sys = "pp";
  else if (IspPb ())  sys = "pPb";
  else                sys = "???";

  const int nFinerEtaBins = 56;
  const double* finerEtaBins = linspace (-2.8, 2.8, nFinerEtaBins);

  const int nRespBins = 240;
  const double* respBins = linspace (0, 2.4, nRespBins);

  const int nEtaRespBins = 80;
  const double* etaRespBins = linspace (-0.2, 0.2, nEtaRespBins);

  const double pTJBins[] = {10, 12, 15, 18, 22, 26, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 360, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1100, 1200, 1300};
  const int nPtJBins = sizeof (pTJBins) / sizeof (pTJBins[0]) - 1;


  // Initialize a bunch of histograms --
  TH1D*** h_r2_jpts     = Get2DArray <TH1D*> (nPtJBins, nFinerEtaBins);
  TH1D*** h_r2_jes      = Get2DArray <TH1D*> (nPtJBins, nFinerEtaBins);
  TH1D*** h_r2_jetacorr = Get2DArray <TH1D*> (nPtJBins, nFinerEtaBins);
  TH1D*** h_r4_jpts     = Get2DArray <TH1D*> (nPtJBins, nFinerEtaBins);
  TH1D*** h_r4_jes      = Get2DArray <TH1D*> (nPtJBins, nFinerEtaBins);
  TH1D*** h_r4_jetacorr = Get2DArray <TH1D*> (nPtJBins, nFinerEtaBins);

  for (int iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

    for (int iEta = 0; iEta < nFinerEtaBins; iEta++) {

      h_r2_jpts[iPtJ][iEta] = new TH1D (Form ("h_r2_jpts_%s_iPtJ%i_iEta%i", sys.Data (), iPtJ, iEta), ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};Counts", nRespBins, respBins);
      h_r2_jpts[iPtJ][iEta]->Sumw2 ();
      h_r2_jes[iPtJ][iEta] = new TH1D (Form ("h_r2_jes_%s_iPtJ%i_iEta%i", sys.Data (), iPtJ, iEta), ";#it{E}_{reco} / #it{E}_{truth};Counts", nRespBins, respBins);
      h_r2_jes[iPtJ][iEta]->Sumw2 ();
      h_r2_jetacorr[iPtJ][iEta] = new TH1D (Form ("h_r2_jetacorr_%s_iPtJ%i_iEta%i", sys.Data (), iPtJ, iEta), ";#eta_{reco} - #eta_{truth};Counts", nEtaRespBins, etaRespBins);
      h_r2_jetacorr[iPtJ][iEta]->Sumw2 ();

      h_r4_jpts[iPtJ][iEta] = new TH1D (Form ("h_r4_jpts_%s_iPtJ%i_iEta%i", sys.Data (), iPtJ, iEta), ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};Counts", nRespBins, respBins);
      h_r4_jpts[iPtJ][iEta]->Sumw2 ();
      h_r4_jes[iPtJ][iEta] = new TH1D (Form ("h_r4_jes_%s_iPtJ%i_iEta%i", sys.Data (), iPtJ, iEta), ";#it{E}_{reco} / #it{E}_{truth};Counts", nRespBins, respBins);
      h_r4_jes[iPtJ][iEta]->Sumw2 ();
      h_r4_jetacorr[iPtJ][iEta] = new TH1D (Form ("h_r4_jetacorr_%s_iPtJ%i_iEta%i", sys.Data (), iPtJ, iEta), ";#eta_{reco} - #eta_{truth};Counts", nEtaRespBins, etaRespBins);
      h_r4_jetacorr[iPtJ][iEta]->Sumw2 ();

    } // end loop over iEta

  } // end loop over iPtJ

  

  const int nEvts = tree->GetEntries ();

  // Loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In JetEnergyResolution.cxx: Event loop " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    tree->GetEntry (iEvt);

    // reject events with pileup vertices or too high z-vertex
    {
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
    }


    // Filter sample based on min/max of pThat range
    if (!IsHijing ()) {
      int iLTJ = -1;
      for (int iTJ = 0; iTJ < akt4_truth_jet_n; iTJ++) {
        if (iLTJ == -1 || akt4_truth_jet_pt[iTJ] > akt4_truth_jet_pt[iLTJ])
          iLTJ = iTJ;
      }

      if (iLTJ == -1 || akt4_truth_jet_pt[iLTJ] < truth_jet_min_pt || akt4_truth_jet_pt[iLTJ] > truth_jet_max_pt)
        continue;
    }


    for (int iJet = 0; iJet < akt2_hi_jet_n; iJet++) {
      if (!MeetsJetAcceptanceCuts (iJet, JetRadius::R0p2))
        continue;

      const float jpt = GetAktHIJetPt (iJet, JetRadius::R0p2);
      const float jeta = GetAktHIJetEta (iJet, JetRadius::R0p2);
      const float jphi = GetAktHIJetPhi (iJet, JetRadius::R0p2);
      const float jen = GetAktHIJetEn (iJet, JetRadius::R0p2);

      const int iTJet = GetAktTruthJetMatch (iJet, JetRadius::R0p2);
      if (iTJet == -1 || DeltaR (jeta, akt2_truth_jet_eta[iTJet], jphi, akt2_truth_jet_phi[iTJet]) > akt2_TruthMatchMaxDR)
        continue;

      const float tjpt = akt2_truth_jet_pt[iTJet];
      const float tjeta = akt2_truth_jet_eta[iTJet];
      const float tjphi = akt2_truth_jet_phi[iTJet];
      const float tjen = akt2_truth_jet_e[iTJet];

      short iPtJ = -1;
      if (pTJBins[0] <= tjpt) {
        iPtJ = 0;
        while (iPtJ < nPtJBins && pTJBins[iPtJ+1] < tjpt) iPtJ++;
      }

      short iEta = -1;
      if (finerEtaBins[0] <= tjeta) {
        iEta = 0;
        while (iEta < nFinerEtaBins && finerEtaBins[iEta+1] < tjeta) iEta++;
      }

      if (iPtJ >= 0 && iPtJ < nPtJBins && iEta >= 0 && iEta < nFinerEtaBins) {
        h_r2_jpts[iPtJ][iEta]->Fill (jpt / tjpt);
        h_r2_jes[iPtJ][iEta]->Fill (jen / tjen);
        h_r2_jetacorr[iPtJ][iEta]->Fill (jeta - tjeta);
      }
    } // end loop over R=0.2 jets


    for (int iJet = 0; iJet < akt4_hi_jet_n; iJet++) {
      if (!MeetsJetAcceptanceCuts (iJet, JetRadius::R0p4))
        continue;

      const float jpt = GetAktHIJetPt (iJet, JetRadius::R0p4);
      const float jeta = GetAktHIJetEta (iJet, JetRadius::R0p4);
      const float jphi = GetAktHIJetPhi (iJet, JetRadius::R0p4);
      const float jen = GetAktHIJetEn (iJet, JetRadius::R0p4);

      const int iTJet = GetAktTruthJetMatch (iJet, JetRadius::R0p4);
      if (iTJet == -1 || DeltaR (jeta, akt4_truth_jet_eta[iTJet], jphi, akt4_truth_jet_phi[iTJet]) > akt4_TruthMatchMaxDR)
        continue;

      const float tjpt = akt4_truth_jet_pt[iTJet];
      const float tjeta = akt4_truth_jet_eta[iTJet];
      const float tjphi = akt4_truth_jet_phi[iTJet];
      const float tjen = akt4_truth_jet_e[iTJet];

      short iPtJ = -1;
      if (pTJBins[0] <= tjpt) {
        iPtJ = 0;
        while (iPtJ < nPtJBins && pTJBins[iPtJ+1] < tjpt) iPtJ++;
      }

      short iEta = -1;
      if (finerEtaBins[0] <= tjeta) {
        iEta = 0;
        while (iEta < nFinerEtaBins && finerEtaBins[iEta+1] < tjeta) iEta++;
      }

      if (iPtJ >= 0 && iPtJ < nPtJBins && iEta >= 0 && iEta < nFinerEtaBins) {
        h_r4_jpts[iPtJ][iEta]->Fill (jpt / tjpt);
        h_r4_jes[iPtJ][iEta]->Fill (jen / tjen);
        h_r4_jetacorr[iPtJ][iEta]->Fill (jeta - tjeta);
      }
    } // end loop over R=0.4 jets
  } // end event loop
  cout << endl << "Info: In JetEnergyResolution.cxx: Finished event loop." << endl;


  SaferDelete (&tree);


  outFile->cd ();

  for (int iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
    for (int iEta = 0; iEta < nFinerEtaBins; iEta++) {
      h_r2_jpts[iPtJ][iEta]->Write ();
      SaferDelete (&(h_r2_jpts[iPtJ][iEta]));
      h_r2_jes[iPtJ][iEta]->Write ();
      SaferDelete (&(h_r2_jes[iPtJ][iEta]));
      h_r2_jetacorr[iPtJ][iEta]->Write ();
      SaferDelete (&(h_r2_jetacorr[iPtJ][iEta]));

      h_r4_jpts[iPtJ][iEta]->Write ();
      SaferDelete (&(h_r4_jpts[iPtJ][iEta]));
      h_r4_jes[iPtJ][iEta]->Write ();
      SaferDelete (&(h_r4_jes[iPtJ][iEta]));
      h_r4_jetacorr[iPtJ][iEta]->Write ();
      SaferDelete (&(h_r4_jetacorr[iPtJ][iEta]));
    } // end loop over iEta
  } // end loop over iPtJ

  outFile->Write (0, TObject::kOverwrite);
  outFile->Close ();
  SaferDelete (&outFile);

  return true;
}

} // end namespace

#endif
