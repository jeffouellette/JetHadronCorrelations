#ifndef __JetPerformance_cxx__
#define __JetPerformance_cxx__

#include "JetPerformance.h"
#include "Params.h"
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

bool JetPerformance (const char* directory,
                     const int dataSet,
                     const char* inFileName,
                     const char* eventWeightsFileName) {
 
  cout << "Info: In JetPerformance.cxx: Entered JetPerformance routine." << endl;
  cout << "Info: In JetPerformance.cxx: Printing systematic onfiguration:";

  SetupDirectories ("JetPerformance");

  if (IsCollisions ()) {
    cout << "Error: In JetPerformance.cxx: Trying to calculate jet performance in data! Quitting." << endl;
    return false;
  }

  if (IsDataOverlay ())
    cout << "Info: In JetPerformance.cxx: Running over data overlay, will check data conditions" << endl;
  if (IsHijing ())
    cout << "Info: In JetPerformance.cxx: Running over Hijing sample" << endl;

  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  cout << "Info: In JetPerformance.cxx: File Identifier: " << identifier << endl;
  cout << "Info: In JetPerformance.cxx: Saving output to " << rootPath << endl;


  // creates a file identifier pattern that we will use to identify the directory containing the input files
  TString fileIdentifier;
  if (TString (inFileName) == "") {
    if (IsCollisions ()) {
      fileIdentifier = to_string (dataSet);
    }
    else {
      std::cout << "Error: In JetPerformance.cxx: Cannot identify this MC file! Quitting." << std::endl;
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
    std::cout << "Info: In JetPerformance.cxx: Chain has " << tree->GetListOfFiles ()->GetEntries () << " files, " << tree->GetEntries () << " entries" << std::endl;
  }

  //TFile* eventWeightsFile = nullptr;
  //TH1D* h_weights = nullptr;

  //eventWeightsFile = new TFile (eventWeightsFileName, "read");
  //h_weights = (TH1D*) eventWeightsFile->Get (Form ("h_PbPb%s_weights_%s", doNchWeighting ? "Nch" : "FCal", isHijing ? "hijing" : "mc"));
  //cout << "Info: In JetPerformance.cxx: Found FCal weighting histogram, " << h_weights->GetName () << endl;

  //First sort jets & tracks into many, smaller TTrees.
  //This is where the sorting based on event information (e.g. centrality, Ntrk, jet pT) will go.
  //Event mixing will take place based on these categories so that total memory usage at any point in time is minimized.

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
    tree->SetBranchAddress ("akt4_truth_jet_n",     &akt4_truth_jet_n);
    tree->SetBranchAddress ("akt4_truth_jet_pt",    &akt4_truth_jet_pt);
    tree->SetBranchAddress ("akt4_truth_jet_eta",   &akt4_truth_jet_eta);
    tree->SetBranchAddress ("akt4_truth_jet_phi",   &akt4_truth_jet_phi);
    tree->SetBranchAddress ("akt4_truth_jet_e",     &akt4_truth_jet_e);
  }

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


  cout << "Info : In JetPerformance.cxx: Saving histograms to " << Form ("%s/%s.root", rootPath.Data (), identifier.Data ()) << endl;
  TFile* outFile = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");

  TString sys;
  if (Ispp ())        sys = "pp";
  else if (IspPb ())  sys = "pPb";
  else                sys = "???";

  const int numFinerEtaBins = 90;
  const double* finerEtaBins = linspace (-4.5, 4.5, numFinerEtaBins);

  const double enJBins[] = {10, 12, 15, 18, 22, 26, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 360, 400, 450, 500};
  const int numEnJBins = sizeof (enJBins) / sizeof (enJBins[0]) - 1;

  TH1D*** h_jpts = new TH1D**[numEnJBins];
  TH1D*** h_jes = new TH1D**[numEnJBins];
  TH1D*** h_jetacorr = new TH1D**[numEnJBins];

  for (int iEnJ = 0; iEnJ < numEnJBins; iEnJ++) {
    h_jpts[iEnJ] = new TH1D*[numFinerEtaBins];
    h_jes[iEnJ] = new TH1D*[numFinerEtaBins];
    h_jetacorr[iEnJ] = new TH1D*[numFinerEtaBins];
    for (int iEta = 0; iEta < numFinerEtaBins; iEta++) {
      h_jpts[iEnJ][iEta] = new TH1D (Form ("h_jpts_%s_iEnJ%i_iEta%i", sys.Data (), iEnJ, iEta), "#it{p}_{T}^{reco} / #it{p}_{T}^{truth}", 140, 0.3, 1.7);
      h_jpts[iEnJ][iEta]->Sumw2 ();
      h_jes[iEnJ][iEta] = new TH1D (Form ("h_jes_%s_iEnJ%i_iEta%i", sys.Data (), iEnJ, iEta), "#it{E}_{reco} / #it{E}_{truth}", 140, 0.3, 1.7);
      h_jes[iEnJ][iEta]->Sumw2 ();
      h_jetacorr[iEnJ][iEta] = new TH1D (Form ("h_jetacorr_%s_iEnJ%i_iEta%i", sys.Data (), iEnJ, iEta), "#eta_{reco} - #eta_{truth}", 80, -0.2, 0.2);
      h_jetacorr[iEnJ][iEta]->Sumw2 ();
    } // end loop over iEta
  } // end loop over iEnJ

  

  const int nEvts = tree->GetEntries ();

  // Loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In JetPerformance.cxx: Event loop " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
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

    for (int iJet = 0; iJet < akt4_hi_jet_n; iJet++) {
      if (!MeetsJetAcceptanceCuts (iJet))
        continue;

      const float jpt = akt4_hi_jet_pt_xcalib[iJet];
      const float jeta = akt4_hi_jet_eta_xcalib[iJet];
      const float jphi = akt4_hi_jet_phi[iJet];
      const float jen = akt4_hi_jet_e_xcalib[iJet];

      int iTJet = -1;
      for (int jTJet = 0; jTJet < akt4_truth_jet_n; jTJet++) {
        if (iTJet == -1 || DeltaR (jeta, akt4_truth_jet_eta[jTJet], jphi, akt4_truth_jet_phi[jTJet]) < DeltaR (jeta, akt4_truth_jet_eta[iTJet], jphi, akt4_truth_jet_phi[iTJet]))
          iTJet = jTJet;
      }

      if (iTJet == -1 || DeltaR (jeta, akt4_truth_jet_eta[iTJet], jphi, akt4_truth_jet_phi[iTJet]) > 0.2)
        continue;

      const float tjpt = akt4_truth_jet_pt[iTJet];
      const float tjeta = akt4_truth_jet_eta[iTJet];
      //const float tjphi = akt4_truth_jet_phi[iTJet];
      const float tjen = akt4_truth_jet_e[iTJet];

      short iEnJ = -1;
      if (enJBins[0] <= tjpt) {
        iEnJ = 0;
        while (iEnJ < numEnJBins && enJBins[iEnJ+1] < tjpt) iEnJ++;
      }

      short iEta = -1;
      if (finerEtaBins[0] <= tjeta) {
        iEta = 0;
        while (iEta < numFinerEtaBins && finerEtaBins[iEta+1] < tjeta) iEta++;
      }

      if (iEnJ >= 0 && iEnJ < numEnJBins && iEta >= 0 && iEta < numFinerEtaBins) {
        h_jpts[iEnJ][iEta]->Fill (jpt / tjpt);
        h_jes[iEnJ][iEta]->Fill (jen / tjen);
        h_jetacorr[iEnJ][iEta]->Fill (jeta - tjeta);
      }
    }
  } // end event loop
  cout << endl << "Info: In JetPerformance.cxx: Finished event loop." << endl;


  SaferDelete (&tree);


  outFile->cd ();

  for (int iEnJ = 0; iEnJ < numEnJBins; iEnJ++) {
    for (int iEta = 0; iEta < numFinerEtaBins; iEta++) {
      h_jpts[iEnJ][iEta]->Write ();
      SaferDelete (&(h_jpts[iEnJ][iEta]));
      h_jes[iEnJ][iEta]->Write ();
      SaferDelete (&(h_jes[iEnJ][iEta]));
      h_jetacorr[iEnJ][iEta]->Write ();
      SaferDelete (&(h_jetacorr[iEnJ][iEta]));
    } // end loop over iEta
  } // end loop over iEnJ

  outFile->Write (0, TObject::kOverwrite);
  outFile->Close ();
  SaferDelete (&outFile);

  return true;
}

} // end namespace

#endif
