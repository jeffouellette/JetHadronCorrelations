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
 
  std::cout << "Info: In JetEnergyResolution.cxx: Entered JetEnergyResolution routine." << std::endl;

  SetupDirectories ("JetEnergyResolution");

  if (IsCollisions ()) {
    std::cout << "Error: In JetEnergyResolution.cxx: Trying to calculate jet performance in data! Quitting." << std::endl;
    return false;
  }

  if (IsDataOverlay ())
    std::cout << "Info: In JetEnergyResolution.cxx: Running over data overlay, will check data conditions" << std::endl;
  if (IsHijing ()) {
    std::cout << "Info: In JetEnergyResolution.cxx: Not configured to process Hijing samples, exiting gracefully!" << std::endl;
    return true;
  }

  // this identifier here corresponds to the name of the output file
  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  std::cout << "Info: In JetEnergyResolution.cxx: File Identifier: " << identifier << std::endl;
  std::cout << "Info: In JetEnergyResolution.cxx: Saving output to " << rootPath << std::endl;


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
  const float truth_jet_min_pt = GetJZXR04MinPt (TString (inFileName));
  const float truth_jet_max_pt = GetJZXR04MaxPt (TString (inFileName));
  if (truth_jet_min_pt != 0)
    std::cout << "Checking for leading truth jet with pT > " << truth_jet_min_pt << std::endl;
  if (truth_jet_max_pt != FLT_MAX)
    std::cout << "Checking for leading truth jet with pT < " << truth_jet_max_pt << std::endl;


  std::cout << std::endl << "Jet acceptance config.:" << std::endl;
  std::cout << "Anti-kT R=0.2:" << std::endl;
  std::cout << "  --> jets truth matched to within dR < " << akt2_TruthMatchMaxDR << std::endl;
  std::cout << "  --> HI jets must be isolated to dR > " << akt2_hi_IsoMinDR << std::endl;
  std::cout << "  --> truth jets must be isolated to dR > " << akt2_truth_IsoMinDR << std::endl;
  std::cout << "  --> HI jets must have minimum pT in isolation calc. " << akt2_hi_IsoMinPt << std::endl;
  std::cout << "  --> truth jets must have minimum pT in isolation calc. " << akt2_truth_IsoMinPt << std::endl;
  std::cout << "Anti-kT R=0.4:" << std::endl;
  std::cout << "  -->jets truth matched to within dR < " << akt4_TruthMatchMaxDR << std::endl;
  std::cout << "  -->HI jets must be isolated to dR > " << akt4_hi_IsoMinDR << std::endl;
  std::cout << "  -->truth jets must be isolated to dR > " << akt4_truth_IsoMinDR << std::endl;
  std::cout << "  -->HI jets must have minimum pT in isolation calc. " << akt4_hi_IsoMinPt << std::endl;
  std::cout << "  -->truth jets must have minimum pT in isolation calc. " << akt4_truth_IsoMinPt << std::endl << std::endl;


  tree->SetBranchAddress ("run_number",     &run_number);
  tree->SetBranchAddress ("event_number",   &event_number);
  tree->SetBranchAddress ("lumi_block",     &lumi_block);


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
  tree->SetBranchAddress ("akt2_hi_jet_LooseBad",     &akt2_hi_jet_LooseBad);


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
  tree->SetBranchAddress ("akt4_hi_jet_LooseBad",     &akt4_hi_jet_LooseBad);


  std::cout << "Info : In JetEnergyResolution.cxx: Saving histograms to " << Form ("%s/%s.root", rootPath.Data (), identifier.Data ()) << std::endl;
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

  const double pTJBins[] = {8, 10, 12, 16, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 190, 220, 250, 280, 310, 340, 370, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1100, 1200, 1300};
  const int nPtJBins = sizeof (pTJBins) / sizeof (pTJBins[0]) - 1;

  const std::vector <JetRadius> radii = {JetRadius::R0p2, JetRadius::R0p4};
  const int nRadii = radii.size ();


  // Initialize a bunch of histograms --
  TH1D**** h_jpts     = Get3DArray <TH1D*> (nRadii, nPtJBins, nFinerEtaBins);
  TH1D**** h_jes      = Get3DArray <TH1D*> (nRadii, nPtJBins, nFinerEtaBins);
  TH1D**** h_jetacorr = Get3DArray <TH1D*> (nRadii, nPtJBins, nFinerEtaBins);

  TH2D** h2_jet_ptreco_pttruth  = Get1DArray <TH2D*> (nRadii);
  TH2D** h2_jet_ereco_etruth    = Get1DArray <TH2D*> (nRadii);


  for (int iR = 0; iR < nRadii; iR++) {

    const int r = (radii[iR] == JetRadius::R0p4 ? 4 : 2);

    for (int iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

      for (int iEta = 0; iEta < nFinerEtaBins; iEta++) {

        h_jpts[iR][iPtJ][iEta] = new TH1D (Form ("h_r%i_jpts_%s_iPtJ%i_iEta%i", r, sys.Data (), iPtJ, iEta), ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};Counts", nRespBins, respBins);
        h_jpts[iR][iPtJ][iEta]->Sumw2 ();
        h_jes[iR][iPtJ][iEta] = new TH1D (Form ("h_r%i_jes_%s_iPtJ%i_iEta%i", r, sys.Data (), iPtJ, iEta), ";#it{E}_{reco} / #it{E}_{truth};Counts", nRespBins, respBins);
        h_jes[iR][iPtJ][iEta]->Sumw2 ();
        h_jetacorr[iR][iPtJ][iEta] = new TH1D (Form ("h_r%i_jetacorr_%s_iPtJ%i_iEta%i", r, sys.Data (), iPtJ, iEta), ";#eta_{reco} - #eta_{truth};Counts", nEtaRespBins, etaRespBins);
        h_jetacorr[iR][iPtJ][iEta]->Sumw2 ();

      } // end loop over iEta

    } // end loop over iPtJ

    h2_jet_ptreco_pttruth[iR] = new TH2D (Form ("h2_r%i_jet_ptreco_pttruth_%s", r, sys.Data ()), ";#it{p}_{T}^{truth} [GeV];#it{p}_{T}^{reco} [GeV];Counts", nPtJBins, pTJBins, nPtJBins, pTJBins);
    h2_jet_ptreco_pttruth[iR]->Sumw2 ();
    h2_jet_ereco_etruth[iR] = new TH2D (Form ("h2_r%i_jet_ereco_etruth_%s", r, sys.Data ()), ";#it{E}_{truth} [GeV];#it{E}_{reco} [GeV];Counts", nPtJBins, pTJBins, nPtJBins, pTJBins);
    h2_jet_ereco_etruth[iR]->Sumw2 ();

  } // end loop over iR

  

  const JetRadius r0p4 = JetRadius::R0p4;
  const int nEvts = tree->GetEntries ();

  // Loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      std::cout << "Info: In JetEnergyResolution.cxx: Event loop " << iEvt / (nEvts / 100) << "\% done...\r" << std::flush;
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


    //const double eventWeight = 1;
    // JZ weighting scheme: sigma * eff * L_int / N_evt
    const double eventWeight = (IsHijing () ? 1 : mcEventWeights->at (0) * crossSectionPicoBarns * mcFilterEfficiency * GetJetLuminosity () / mcNumberEvents);


    // Fill histograms for all jet radii
    for (int iR = 0; iR < nRadii; iR++) {

      const JetRadius radius = radii[iR];

      const float jet_TruthMatchMaxDR = GetAktTruthMatchMaxDR (radius); // maximum truth match delta R cut


      // NEW VERSION -- loop over truth jets
      const int nTJets = GetAktTruthJetN (radius);
      for (int iTJet = 0; iTJet < nTJets; iTJet++) {

        const float tjpt = GetAktTruthJetPt (iTJet, radius);
        const float tjeta = GetAktTruthJetEta (iTJet, radius);
        const float tjphi = GetAktTruthJetPhi (iTJet, radius);
        const float tjen = GetAktTruthJetEn (iTJet, radius);

        if (GetAktTruthJetIso (iTJet, radius) < GetAktTruthIsoMinDR (radius))
          continue; // truth jet isolation cut

        const int iJet = GetAktHIJetMatch (iTJet, radius);
        if (iJet == -1)
          continue; // in case of no reco. jet match

        const float jpt = GetAktHIJetPt (iJet, radius);
        const float jeta = GetAktHIJetEta (iJet, radius);
        const float jphi = GetAktHIJetPhi (iJet, radius);
        const float jen = GetAktHIJetEn (iJet, radius);

        if (DeltaR (jeta, tjeta, jphi, tjphi) > jet_TruthMatchMaxDR)
          continue; // reco.-truth jet matching cut

        if (!MeetsJetAcceptanceCuts (iJet, radius))
          continue; // reco. jet acceptance cut

        if (GetAktHIJetIso (iJet, radius) < GetAktHIIsoMinDR (radius))
          continue; // reco. jet isolation cut

        if (jpt / tjpt < 0.4 && GetAktHIJetPt (iJet, radius, -1, 2) < 15)
          continue; // cut out very low pT jets (before calibration), except in cases with a correspondingly low pT truth match.

        h2_jet_ptreco_pttruth[iR]->Fill (tjpt, jpt, eventWeight);
        h2_jet_ereco_etruth[iR]->Fill (tjen, jen, eventWeight);

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
          h_jpts[iR][iPtJ][iEta]->Fill (jpt / tjpt, eventWeight);
          h_jes[iR][iPtJ][iEta]->Fill (jen / tjen, eventWeight);
          h_jetacorr[iR][iPtJ][iEta]->Fill (jeta - tjeta, eventWeight);
        }

      } // end loop over jets

      // OLD VERSION -- loop over reco. jets
      /*
      const int jn = GetAktHIJetN (radius);
      for (int iJet = 0; iJet < jn; iJet++) {
        if (!MeetsJetAcceptanceCuts (iJet, radius))
          continue;

        const float jpt = GetAktHIJetPt (iJet, radius);
        const float jeta = GetAktHIJetEta (iJet, radius);
        const float jphi = GetAktHIJetPhi (iJet, radius);
        const float jen = GetAktHIJetEn (iJet, radius);

        if (GetAktHIJetIso (iJet, radius) < GetAktHIIsoMinDR (radius))
          continue; // reco. jet isolation cut

        const int iTJet = GetAktTruthJetMatch (iJet, radius);
        if (iTJet == -1)
          continue; // in case of no truth jet match

        const float tjpt = GetAktTruthJetPt (iTJet, radius);
        const float tjeta = GetAktTruthJetEta (iTJet, radius);
        const float tjphi = GetAktTruthJetPhi (iTJet, radius);
        const float tjen = GetAktTruthJetEn (iTJet, radius);

        if (DeltaR (jeta, tjeta, jphi, tjphi) > jet_TruthMatchMaxDR)
          continue; // reco.-truth jet matching cut

        if (GetAktTruthJetIso (iTJet, radius) < GetAktTruthIsoMinDR (radius))
          continue; // truth jet isolation cut

        if (jpt / tjpt < 0.4 && GetAktHIJetPt (iJet, radius, -1, 2) < 15) // cuts on pre-calibrated jet pT
          continue; // cut out very low pT jets (before calibration), except in cases with a correspondingly low pT truth match.

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
          h_jpts[iR][iPtJ][iEta]->Fill (jpt / tjpt);
          h_jes[iR][iPtJ][iEta]->Fill (jen / tjen);
          h_jetacorr[iR][iPtJ][iEta]->Fill (jeta - tjeta);
        }

      } // end loop over jets
      */

    } // end loop over jet radii

  } // end event loop

  std::cout << std::endl << "Info: In JetEnergyResolution.cxx: Finished event loop." << std::endl;


  SaferDelete (&tree);


  outFile->cd ();

  for (int iR = 0; iR < nRadii; iR++) {

    for (int iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

      for (int iEta = 0; iEta < nFinerEtaBins; iEta++) {

        h_jpts[iR][iPtJ][iEta]->Write ();
        SaferDelete (&(h_jpts[iR][iPtJ][iEta]));
        h_jes[iR][iPtJ][iEta]->Write ();
        SaferDelete (&(h_jes[iR][iPtJ][iEta]));
        h_jetacorr[iR][iPtJ][iEta]->Write ();
        SaferDelete (&(h_jetacorr[iR][iPtJ][iEta]));

      } // end loop over iEta

    } // end loop over iPtJ

    h2_jet_ptreco_pttruth[iR]->Write ();
    SaferDelete (&(h2_jet_ptreco_pttruth[iR]));
    h2_jet_ereco_etruth[iR]->Write ();
    SaferDelete (&(h2_jet_ereco_etruth[iR]));

  } // end loop over iR

  outFile->Write (0, TObject::kOverwrite);
  outFile->Close ();
  SaferDelete (&outFile);

  return true;
}

} // end namespace

#endif
