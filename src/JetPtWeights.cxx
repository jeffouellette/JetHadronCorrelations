#ifndef __JetPtWeights_cxx__
#define __JetPtWeights_cxx__

#include "JetPtWeights.h"
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

bool JetPtWeights (const char* directory,
                   const int dataSet,
                   const char* inFileName,
                   const char* eventWeightsFileName) {
 
  std::cout << "Info: In JetPtWeights.cxx: Entered JetPtWeights routine." << std::endl;

  SetupDirectories ("JetPtWeights");

  if (IsDataOverlay ())
    std::cout << "Info: In JetPtWeights.cxx: Running over data overlay, will check data conditions" << std::endl;
  if (IsHijing ()) {
    std::cout << "Info: In JetPtWeights.cxx: Not configured to process Hijing samples, exiting gracefully!" << std::endl;
    return true;
  }

  // this identifier here corresponds to the name of the output file
  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  std::cout << "Info: In JetPtWeights.cxx: File Identifier: " << identifier << std::endl;
  std::cout << "Info: In JetPtWeights.cxx: Saving output to " << rootPath << std::endl;


  // creates a file identifier pattern that we will use to identify the directory containing the input files
  TString fileIdentifier;
  if (TString (inFileName) == "") {
    if (IsCollisions ()) {
      fileIdentifier = to_string (dataSet);
    }
    else {
      std::cout << "Error: In JetPtWeights.cxx: Cannot identify this MC file! Quitting." << std::endl;
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

      if (IsCollisions ()) {
        if (Ispp ()) {
          if (UseMinBiasTriggers () && !file.Contains ("MinBias"))
            continue;
          if (!UseMinBiasTriggers () && !file.Contains ("Main"))
            continue;
        }
        if (IspPb ()) {
          if (!file.Contains ("Main"))
            continue;
        }
      }

      if (!file.Contains (fileIdentifier))
        continue;
      std::cout << "Adding " << dataPath + directory + "/" + file + "/*.root" << " to TChain" << std::endl;
      tree->Add (dataPath + directory + "/" + file + "/*.root");
      break;
    }
    std::cout << "Chain has " << tree->GetListOfFiles ()->GetEntries () << " files, " << tree->GetEntries () << " entries" << std::endl;
  }

  if (!IsCollisions () && !IsHijing ()) {
    assert (crossSectionPicoBarns > 0);
    assert (mcFilterEfficiency > 0);
    assert (mcNumberEvents > 0);
  }



  // Identify appropriate trigger arrays (e.g. use pp triggers instead of p+Pb triggers)
  if (IspPb16 ()) {
    minbias_trig_n = minbias_trig_n_pPb16;
    minbias_trig_name = minbias_trig_name_pPb16;
    jet_trig_n = jet_trig_n_pPb16s5TeV;
    jet_trig_name = jet_trig_name_pPb16s5TeV;
  }
  else if (Ispp17 ()) {
    minbias_trig_n = minbias_trig_n_pp17;
    minbias_trig_name = minbias_trig_name_pp17;
    jet_trig_n = jet_trig_n_pp17;
    jet_trig_name = jet_trig_name_pp17;
  }
  else {
    std::cout << "Error: In RunCorrelator.cxx: Invalid or unsupported collision system, exiting." << std::endl;
    return false;
  }


  // pointers to trigger objects
  Trigger* trigger = nullptr;

  // assign branches for triggers
  if (IsCollisions ()) {
    if (UseJ50Triggers ()) {
      trigger = new Trigger (jet_trig_name[0]);
      std::cout << "Info: In RunCorrelator.cxx: Looking for " << jet_trig_name[0] << " trigger" << std::endl;
      tree->SetBranchAddress ((jet_trig_name[0]+"_decision").c_str (), &(trigger->trigDecision));
      tree->SetBranchAddress ((jet_trig_name[0]+"_prescale").c_str (), &(trigger->trigPrescale));
    }

    else if (UseMinBiasTriggers ()) {
      trigger = new Trigger (minbias_trig_name[0]);
      std::cout << "Info: In RunCorrelator.cxx: Looking for " << minbias_trig_name[0] << " trigger" << std::endl;
      tree->SetBranchAddress ((minbias_trig_name[0]+"_decision").c_str (), &(trigger->trigDecision));
      tree->SetBranchAddress ((minbias_trig_name[0]+"_prescale").c_str (), &(trigger->trigPrescale));
    }
  }


  // variables for filtering MC truth
  const float truth_jet_min_pt = GetJZXR04MinPt (TString (inFileName));
  const float truth_jet_max_pt = GetJZXR04MaxPt (TString (inFileName));
  if (truth_jet_min_pt != 0)
    std::cout << "Checking for leading truth jet with pT > " << truth_jet_min_pt << std::endl;
  if (truth_jet_max_pt != FLT_MAX)
    std::cout << "Checking for leading truth jet with pT < " << truth_jet_max_pt << std::endl;


  // weight for this JZ slice
  const float jzScaleFactor = GetJZScaleFactor (TString (inFileName)) * crossSectionPicoBarns * mcFilterEfficiency * GetJetLuminosity () / mcNumberEvents; // sigma * f * L_int


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

  if (IsCollisions () && !Ispp ()) {
    tree->SetBranchAddress ("ZdcCalibEnergy_A", &ZdcCalibEnergy_A);
    tree->SetBranchAddress ("ZdcCalibEnergy_C", &ZdcCalibEnergy_C);
  }


  // get truth jets from MC
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
  tree->SetBranchAddress ("akt4_hi_jet_LooseBad",   &akt4_hi_jet_LooseBad);


  const short nJESVar = GetNJESVar ();
  if (!IsCollisions () && nJESVar != -1) {
    std::cout << "Info: In JetPtWeights.cxx: Branching JES variation " << nJESVar << std::endl;
    tree->SetBranchAddress (Form ("akt4_hi_jet_pt_sys_JES_%i", nJESVar), akt4_hi_jet_pt_sys_JES_ALL[nJESVar]);
  }


  // setup centrality bins (only relevant for p+Pb)
  double* centBins = (DoFcalCentVar () ? fcalCentBins : (DoFineFcalCentVar () ? fineFcalCentBins : (!IsCollisions () ? fcalCentBins : zdcCentBins)));
  const short nCentBins = (DoFcalCentVar () ? nFcalCentBins : (DoFineFcalCentVar () ? nFineFcalCentBins : (!IsCollisions () ? nFcalCentBins : nZdcCentBins)));

  std::cout << "Centrality bin cuts: ";
  for (short iCent = 0; iCent < nCentBins; iCent++)
    std::cout << centBins[iCent] << ", ";
  std::cout << centBins[nCentBins] << std::endl;


  // for jet energy resolution cuts
  TF1* f_jer = LoadJetEnergyResFunction ();


  //// for corrected FCal Et values in data overlay
  //std::map <const unsigned int, float>* m_overlay_fcalet = (IsDataOverlay () ? GetOverlayFCalMap () : nullptr);


  TString sys;
  if (Ispp ())        sys = "pp";
  else if (IspPb ())  sys = "pPb";
  else                sys = "???";

  const double finePtJBins[] = {10, 11, 12, 13, 14, 15, 17.5, 20, 22.5, 25, 27.5, 30, 33, 36, 40, 45, 50, 55, 60, 65, 70, 75, 82.5, 90, 100, 110, 120, 130, 145, 160, 180, 200, 220, 240, 260, 280, 300, 325, 350, 375, 400};
  //const double finePtJBins[] = {7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5, 21, 21.5, 22, 22.5, 23, 23.5, 24, 24.5, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 47.5, 50, 52.5, 55, 57.5, 60, 62.5, 65, 67.5, 70, 72.5, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 240, 260, 280, 300, 325, 350, 375, 400};
  const short nFinePtJBins = sizeof (finePtJBins) / sizeof (finePtJBins[0]) - 1;


  const double etaJBins[] = {-2.8, -2.1, -1.2, -0.8, -0.3, 0, 0.3, 0.8, 1.2, 2.1, 2.8};
  const short nEtaJBins = sizeof (etaJBins) / sizeof (etaJBins[0]) - 1;


  // create output histograms in output files
  const short nFiles = (Ispp () ? 1 : nCentBins);
  TFile** outFiles = new TFile*[nFiles];

  TH1D**  h_evt_counts    = Get1DArray <TH1D*> (nFiles);

  TH1D**  h_jet_pt        = Get1DArray <TH1D*> (nFiles);
  TH2D**  h2_jet_pt_cov   = Get1DArray <TH2D*> (nFiles);

  TH2D*** h2_jet_eta_phi  = Get2DArray <TH2D*> (nFiles, nFinePtJBins);

  TH2D**  h2_jet_pt_eta_jer_frac_num = Get1DArray <TH2D*> (nFiles);
  TH2D**  h2_jet_pt_eta_jer_frac_den = Get1DArray <TH2D*> (nFiles);


  for (short iFile = 0; iFile < nFiles; iFile++) {

    TString outFileName = Form ("%s/%s%s%s.root", rootPath.Data (), identifier.Data (), nFiles > 1 ? Form ("_iCent%i", iFile) : "", IsCollisions () ? (UseMinBiasTriggers () ? "_mb" : "_j50") : "");
    std::cout << "Info : In JetPtWeights.cxx: Saving histograms to " << outFileName.Data () << std::endl;
    TFile* outFile = new TFile (outFileName.Data (), "recreate");
    outFiles[iFile] = outFile;

    h_evt_counts[iFile] = new TH1D ("h_evt_counts",  "", 3, -0.5, 2.5);

    h_jet_pt[iFile] = new TH1D ("h_jet_pt", ";#it{p}_{T}^{jet} [GeV];(1/N_{jet}) (dN_{jet}/d#it{p}_{T}) [GeV^{-1}]", nFinePtJBins, finePtJBins);
    h2_jet_pt_cov[iFile] = new TH2D ("h2_jet_pt_cov", ";#it{p}_{T}^{jet} [GeV];#it{p}_{T}^{jet} [GeV];Covariance", nFinePtJBins, finePtJBins, nFinePtJBins, finePtJBins);

    h2_jet_pt_eta_jer_frac_num[iFile] = new TH2D ("h2_jet_pt_eta_jer_frac_num", ";#it{p}_{T}^{jet} [GeV];#it{#eta}^{jet};N_{jet}", nFinePtJBins, finePtJBins, nEtaJBins, etaJBins);
    h2_jet_pt_eta_jer_frac_num[iFile]->Sumw2 ();

    h2_jet_pt_eta_jer_frac_den[iFile] = new TH2D ("h2_jet_pt_eta_jer_frac_den", ";#it{p}_{T}^{jet} [GeV];#it{#eta}^{jet};N_{jet}", nFinePtJBins, finePtJBins, nEtaJBins, etaJBins);
    h2_jet_pt_eta_jer_frac_den[iFile]->Sumw2 ();

    for (short iPtJ = 0; iPtJ < nFinePtJBins; iPtJ++) {

      const TString pTJ = Form ("%g-%gGeVJets", finePtJBins[iPtJ], finePtJBins[iPtJ+1]);

      h2_jet_eta_phi[iFile][iPtJ] = new TH2D (Form ("h2_jet_eta_phi_%s", pTJ.Data ()), ";#eta^{jet};#phi^{jet};Counts", 64, -3.2, 3.2, 64, -M_PI, M_PI);
    } // end loop over iPtJ

  } // end loop over iFile

  // arrays for filling histograms & covariance matrices correctly
  double* jet_pt_counts         = Get1DArray <double> (nFinePtJBins);


  const JetRadius r0p4 = JetRadius::R0p4;


  // Load histograms with flavour-related uncertainties. (Only actually used if nJESVar = 18 or 19.)
  TH2D* h2_flavourFracUnc = GetFlavorFractionUnc (r0p4);
  TH2D* h2_flavourRespUnc = GetFlavorResponseUnc (r0p4);


  const int nEvts = tree->GetEntries ();


  // Loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      std::cout << "Info: In JetPtWeights.cxx: Event loop " << iEvt / (nEvts / 100) << "\% done...\r" << std::flush;

    tree->GetEntry (iEvt);


    float prescale = 1;


    if (IsCollisions ()) {
      const short jn = UseTruthJets () ? GetAktTruthJetN (r0p4) :  GetAktHIJetN (r0p4);
      float ljpt = 0;
      for (short iJet = 0; iJet < jn; iJet++) {

        if (!MeetsJetAcceptanceCuts (iJet, r0p4, nJESVar))
          continue; // jet eta/phi & timing cuts

        ljpt = std::fmax (ljpt, GetAktHIJetPt  (iJet, r0p4, nJESVar));

      } // end loop over iJet

      if (ljpt < 60 && UseJ50Triggers ())
        continue;
      else if (ljpt >= 60 && UseMinBiasTriggers ())
        continue;


      // triggering cut, require appropriate jet trigger to have fired
      if (!trigger->trigDecision)
        continue;
      event_weight = 1; // optionally set to trigger prescale, but not really necessary
      prescale = trigger->trigPrescale;
    }
    else {
      prescale = (IsHijing () ? 1 : jzScaleFactor * mcEventWeights->at (0));
    }


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


    {
      bool isWellTimed = true;
      if (IsCollisions () && Ispp () && UseMinBiasTriggers ()) {
        for (short iJ = 0; iJ < GetAktHIJetN (r0p4); iJ++) {
          if (!MeetsJetAcceptanceCuts (iJ, r0p4, nJESVar))
            continue; // jet eta/phi & timing cuts
          if (GetAktHIJetPt (iJ, r0p4) < 15)
            continue; // minimum pT cut
          if (GetAktHIJetTiming (iJ, r0p4) > 10)
            isWellTimed = false;
        }
      }
      if (!isWellTimed)
        continue; // skip events with a mistimed jet in pp
    }


    // Filter sample based on min/max of pThat range
    if (!IsCollisions () && !IsHijing ()) {
      int iLTJ = -1;
      const int nTJ = GetAktTruthJetN (r0p4);
      for (int iTJ = 0; iTJ < nTJ; iTJ++) {
        if (iLTJ == -1 || GetAktTruthJetPt (iTJ, r0p4) > GetAktTruthJetPt (iLTJ, r0p4))
          iLTJ = iTJ;
      }

      if (iLTJ == -1 || GetAktTruthJetPt (iLTJ, r0p4) < truth_jet_min_pt || GetAktTruthJetPt (iLTJ, r0p4) > truth_jet_max_pt)
        continue;
    }


    // p+Pb only -- split events by centrality
    short iCent = 0;
    if (IspPb ()) {
      const float centVar = ((!IsCollisions () || DoFcalCentVar ()  || DoFineFcalCentVar ()) ? fcalA_et : (ZdcCalibEnergy_A * 1e3));
      iCent = GetBin (centBins, nCentBins, centVar);
      if (iCent < 0 || nCentBins <= iCent)
        continue;
    }
    const short iFile = iCent;


    // event weights -- these are not 1 in MC due to JZ weighting (they are 1 in data)
    const double ewgt = prescale;
    if (ewgt <= 0.)
      continue;


    // find truth-reco jet matches uniquely
    const std::vector <short> recoJetMatches = GetAktRecoJetMatches (r0p4, nJESVar);
    const std::vector <short> truthJetMatches = GetAktTruthJetMatches (recoJetMatches, r0p4);


    // initialize all histogramming bins to 0 for this event
    for (short iPtJ = 0; iPtJ < nFinePtJBins; iPtJ++)
      jet_pt_counts[iPtJ] = 0;

   
    short nj = 0;

    const short jn = UseTruthJets () ? GetAktTruthJetN (r0p4) :  GetAktHIJetN (r0p4);
    for (short iJet = 0; iJet < jn; iJet++) {

      float rjpt  = (UseTruthJets () ? GetAktTruthJetPt  (iJet, r0p4) : GetAktHIJetPt  (iJet, r0p4, nJESVar, IsCollisions () || (IsDataOverlay () && truthJetMatches[iJet] == -1) ? 0 : 1));
      const float rjeta = (UseTruthJets () ? GetAktTruthJetEta (iJet, r0p4) : GetAktHIJetEta (iJet, r0p4, nJESVar));
      const float rjphi = (UseTruthJets () ? GetAktTruthJetPhi (iJet, r0p4) : GetAktHIJetPhi (iJet, r0p4, nJESVar));

      if (nJESVar == 18) // Flavour-dependent response unc.
        rjpt *= 1 + (h2_flavourRespUnc->GetBinContent (h2_flavourRespUnc->GetXaxis ()->FindBin (rjpt), h2_flavourRespUnc->GetYaxis ()->FindBin (IspPb () ? rjeta : std::fabs (rjeta))));
      else if (nJESVar == 19) // Flavour fraction unc.
        rjpt *= 1 + (h2_flavourFracUnc->GetBinContent (h2_flavourFracUnc->GetXaxis ()->FindBin (rjpt), h2_flavourFracUnc->GetYaxis ()->FindBin (IspPb () ? rjeta : std::fabs (rjeta))));

      if (!MeetsJetAcceptanceCuts (iJet, r0p4, nJESVar))
        continue; // jet eta/phi & timing cuts
      if (!MeetsJetPtCut (rjpt))
        continue; // jet pT cuts, for trigger relevance purposes

      const float thisjwgt = GetAktJetWeight (rjpt, rjeta, rjphi, r0p4);
      if (thisjwgt <= 0.)
        continue; // sanity check

      // in MC, check that the jet is truth-matched. Otherwise problems can emerge, particularly in the overlay.
      if (!IsCollisions () && !UseTruthJets ()) {
        //const short iTJet = GetAktTruthJetMatch (iJet, r0p4);
        const short iTJet = truthJetMatches[iJet];
        if (iTJet == -1)
          continue; // in case of no truth jet match

        h2_jet_pt_eta_jer_frac_den[iFile]->Fill (rjpt, rjeta, ewgt * thisjwgt);
        const float tjpt = GetAktTruthJetPt (iTJet, r0p4);
        if (std::fabs (rjpt/tjpt - 1) > nJERSigma*0.01*f_jer->Eval (tjpt))
          continue; // cut on jets reconstructed well outside the JER
        h2_jet_pt_eta_jer_frac_num[iFile]->Fill (rjpt, rjeta, ewgt * thisjwgt);

        //if (rjpt < truth_jet_min_pt || rjpt > truth_jet_max_pt)
        //  continue; // cut on jets outside truth range. Ensures that only one sample fills each pT region.
      } 

      // find which pTJ bin this jet belongs in
      short iPtJ = 0;
      if (rjpt < finePtJBins[0])
        iPtJ = -1;
      else {
        while (iPtJ < nFinePtJBins) {
          if (rjpt < finePtJBins[iPtJ+1])
            break;
          else
            iPtJ++;
        }
      }


      if (0 <= iPtJ && iPtJ < nFinePtJBins) {
        jet_pt_counts[iPtJ] += thisjwgt;
        h2_jet_eta_phi[iFile][iPtJ]->Fill (rjeta, rjphi);
        nj++; // add to total number of "trigger" jets
      }
    }

    // skip events with no jets (otherwise the "trigger" didn't fire). This just saves time with mixing really.
    if (nj == 0)
      continue;

    for (short iPtJX = 0; iPtJX < nFinePtJBins; iPtJX++) {
      h_jet_pt[iCent]->SetBinContent (iPtJX+1, h_jet_pt[iCent]->GetBinContent (iPtJX+1) + (ewgt)*(jet_pt_counts[iPtJX]));
      for (short iPtJY = 0; iPtJY < nFinePtJBins; iPtJY++)
        h2_jet_pt_cov[iCent]->SetBinContent (iPtJX+1, iPtJY+1, h2_jet_pt_cov[iCent]->GetBinContent (iPtJX+1, iPtJY+1) + (ewgt)*(jet_pt_counts[iPtJX])*(jet_pt_counts[iPtJY]));
    }


    // the weights for the per-trigger jet yield are "w_evt"
    h_evt_counts[iFile]->SetBinContent (1, h_evt_counts[iFile]->GetBinContent (1) + 1);
    h_evt_counts[iFile]->SetBinContent (2, h_evt_counts[iFile]->GetBinContent (2) + ewgt);
    h_evt_counts[iFile]->SetBinContent (3, h_evt_counts[iFile]->GetBinContent (3) + ewgt*ewgt);

  } // end event loop

  std::cout << std::endl << "Info: In JetPtWeights.cxx: Finished event loop." << std::endl;


  SaferDelete (&tree);

  SaferDelete (&f_jer);

  //SaferDelete (&m_overlay_fcalet);


  for (short iFile = 0; iFile < nFiles; iFile++) {

    outFiles[iFile]->cd ();

    h_evt_counts[iFile]->Write ();

    h_jet_pt[iFile]->Write ();
    h2_jet_pt_cov[iFile]->Write ();

    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

      h2_jet_eta_phi[iFile][iPtJ]->Write ();

    } // end loop over iPtJ

    h2_jet_pt_eta_jer_frac_num[iFile]->Write ();
    h2_jet_pt_eta_jer_frac_den[iFile]->Write ();


    outFiles[iFile]->Close ();
  }

  SaferDelete (&trigger);

  return true;
}

} // end namespace

#endif
