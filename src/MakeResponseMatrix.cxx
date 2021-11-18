#ifndef __JetHadronCorrelations_MakeResponseMatrix_cxx__
#define __JetHadronCorrelations_MakeResponseMatrix_cxx__

#include "Params.h"
#include "CentralityDefs.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Variations.h"
#include "Process.h"

#include <ArrayTemplates.h>
#include <Utilities.h>

#include <TChain.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TString.h>

//#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
//#endif

#include <iostream>
#include <string.h>
#include <math.h>


using namespace JetHadronCorrelations;



bool MakeResponseMatrix (const char* directory,
                         const char* inFileName,
                         const char* outFileName) {

  std::cout << "Info: In MakeResponseMatrix.cxx: Entered MakeResponseMatrix routine." << std::endl;

  SetupDirectories ("MakeResponseMatrix");


  const int nFiles = (Ispp () ? 1 : nFcalCentBins+1);

  //const double pTJBins[] = {20, 30, 45, 60, 80, 100, 130, 160, 200, 240, 320};
  //const short nPtJBins = sizeof (pTJBins) / sizeof (pTJBins[0]) - 1;

  TH1D***    h_jet_pt_wgts                  = Get2DArray <TH1D*> (nFiles, 2);
  TH1D***    h_jet_pt_fullClosure           = Get2DArray <TH1D*> (nFiles, 2);
  TH1D***    h_jet_pt_halfClosure           = Get2DArray <TH1D*> (nFiles, 2);
  TH2D****   h2_jet_trk_pt_sig_wgts         = Get3DArray <TH2D*> (nDir, nFiles, 2);
  TH2D****   h2_jet_trk_pt_sig_fullClosure  = Get3DArray <TH2D*> (nDir, nFiles, 2);
  TH2D****   h2_jet_trk_pt_sig_halfClosure  = Get3DArray <TH2D*> (nDir, nFiles, 2);


  RooUnfoldResponse**   rooUnfResp_jet_pt_wgts                = Get1DArray <RooUnfoldResponse*> (nFiles);
  RooUnfoldResponse**   rooUnfResp_jet_pt_altwgts             = Get1DArray <RooUnfoldResponse*> (nFiles);
  RooUnfoldResponse**   rooUnfResp_jet_pt_fullClosure         = Get1DArray <RooUnfoldResponse*> (nFiles);
  RooUnfoldResponse**   rooUnfResp_jet_pt_halfClosure         = Get1DArray <RooUnfoldResponse*> (nFiles);
  RooUnfoldResponse***  rooUnfResp_jet_trk_pt_sig_wgts        = Get2DArray <RooUnfoldResponse*> (nDir, nFiles);
  RooUnfoldResponse***  rooUnfResp_jet_trk_pt_sig_altwgts     = Get2DArray <RooUnfoldResponse*> (nDir, nFiles);
  RooUnfoldResponse***  rooUnfResp_jet_trk_pt_sig_fullClosure = Get2DArray <RooUnfoldResponse*> (nDir, nFiles);
  RooUnfoldResponse***  rooUnfResp_jet_trk_pt_sig_halfClosure = Get2DArray <RooUnfoldResponse*> (nDir, nFiles);


  const TString identifier = GetIdentifier (0, directory, inFileName);
  std::cout << "Info: In MakeResponseMatrix.cxx: File Identifier: " << identifier << std::endl;
  std::cout << "Info: In MakeResponseMatrix.cxx: Saving output to " << rootPath << std::endl;

  // Load files for output
  //TFile* outFile = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");
  TFile* outFile = new TFile (Form ("%s/%s", rootPath.Data (), outFileName), "recreate");


  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CREATE SIGNAL YIELDS 2D HISTOGRAM WITH JET PT on x-axis and PARTICLE PT on y-axis
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iVar = 0; iVar < 2; iVar++) {

    const TString var = (iVar == 0 ? "Nominal" : "MCTruthJetsTruthParts");

    for (int iFile = 0; iFile < nFiles; iFile++) {

      const char* cent = (Ispp () ? "ref" : (iFile == nZdcCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iFile)));

      h_jet_pt_wgts[iFile][iVar] = new TH1D (Form ("h_jet_pt_wgts_%s_mc_%s",  cent, var.Data ()), ";#it{p}_{T}^{jet} [GeV]", nPtJBins, pTJBins);
      h_jet_pt_wgts[iFile][iVar]->Sumw2 ();
      h_jet_pt_fullClosure[iFile][iVar] = new TH1D (Form ("h_jet_pt_fullClosure_%s_mc_%s",  cent, var.Data ()), ";#it{p}_{T}^{jet} [GeV]", nPtJBins, pTJBins);
      h_jet_pt_fullClosure[iFile][iVar]->Sumw2 ();
      h_jet_pt_halfClosure[iFile][iVar] = new TH1D (Form ("h_jet_pt_halfClosure_%s_mc_%s",  cent, var.Data ()), ";#it{p}_{T}^{jet} [GeV]", nPtJBins, pTJBins);
      h_jet_pt_halfClosure[iFile][iVar]->Sumw2 ();

      for (short iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];
  
        h2_jet_trk_pt_sig_wgts[iDir][iFile][iVar]  = new TH2D (Form ("h2_jet_trk_pt_wgts_%s_%s_sig_mc_%s",  dir.Data (), cent, var.Data ()), ";#it{p}_{T}^{jet} [GeV];#it{p}_{T}^{ch} [GeV]", nPtChBins, pTChBins, nPtJBins, pTJBins);
        h2_jet_trk_pt_sig_wgts[iDir][iFile][iVar]->Sumw2 ();
        h2_jet_trk_pt_sig_fullClosure[iDir][iFile][iVar]  = new TH2D (Form ("h2_jet_trk_pt_fullClosure_%s_%s_sig_mc_%s",  dir.Data (), cent, var.Data ()), ";#it{p}_{T}^{jet} [GeV];#it{p}_{T}^{ch} [GeV]", nPtChBins, pTChBins, nPtJBins, pTJBins);
        h2_jet_trk_pt_sig_fullClosure[iDir][iFile][iVar]->Sumw2 ();
        h2_jet_trk_pt_sig_halfClosure[iDir][iFile][iVar]  = new TH2D (Form ("h2_jet_trk_pt_halfClosure_%s_%s_sig_mc_%s",  dir.Data (), cent, var.Data ()), ";#it{p}_{T}^{jet} [GeV];#it{p}_{T}^{ch} [GeV]", nPtChBins, pTChBins, nPtJBins, pTJBins);
        h2_jet_trk_pt_sig_halfClosure[iDir][iFile][iVar]->Sumw2 ();

      } // end loop over iFile

    } // end loop over iDir

  } // end loop over iVar




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // BUILD RESPONSE MATRICES
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (int iFile = 0; iFile < nFiles; iFile++) {

    const char* cent = (Ispp () ? "ref" : (iFile == nZdcCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iFile)));

    rooUnfResp_jet_pt_wgts[iFile]         = new RooUnfoldResponse (h_jet_pt_fullClosure[iFile][0], h_jet_pt_fullClosure[iFile][1], Form ("rooUnfResp_jet_pt_%s_mc_wgts", cent), "");
    rooUnfResp_jet_pt_altwgts[iFile]      = new RooUnfoldResponse (h_jet_pt_fullClosure[iFile][0], h_jet_pt_fullClosure[iFile][1], Form ("rooUnfResp_jet_pt_%s_mc_altwgts", cent), "");
    rooUnfResp_jet_pt_fullClosure[iFile]  = new RooUnfoldResponse (h_jet_pt_fullClosure[iFile][0], h_jet_pt_fullClosure[iFile][1], Form ("rooUnfResp_jet_pt_%s_mc_fullClosure", cent), "");
    rooUnfResp_jet_pt_halfClosure[iFile]  = new RooUnfoldResponse (h_jet_pt_halfClosure[iFile][0], h_jet_pt_halfClosure[iFile][1], Form ("rooUnfResp_jet_pt_%s_mc_halfClosure", cent), "");

    for (short iDir = 0; iDir < nDir; iDir++) {

      const TString dir = directions[iDir];

      rooUnfResp_jet_trk_pt_sig_wgts[iDir][iFile]         = new RooUnfoldResponse (h2_jet_trk_pt_sig_fullClosure[iDir][iFile][0], h2_jet_trk_pt_sig_fullClosure[iDir][iFile][1], Form ("rooUnfResp_jet_trk_pt_%s_%s_sig_mc_wgts", dir.Data (), cent), "");
      rooUnfResp_jet_trk_pt_sig_altwgts[iDir][iFile]      = new RooUnfoldResponse (h2_jet_trk_pt_sig_fullClosure[iDir][iFile][0], h2_jet_trk_pt_sig_fullClosure[iDir][iFile][1], Form ("rooUnfResp_jet_trk_pt_%s_%s_sig_mc_altwgts", dir.Data (), cent), "");
      rooUnfResp_jet_trk_pt_sig_fullClosure[iDir][iFile]  = new RooUnfoldResponse (h2_jet_trk_pt_sig_fullClosure[iDir][iFile][0], h2_jet_trk_pt_sig_fullClosure[iDir][iFile][1], Form ("rooUnfResp_jet_trk_pt_%s_%s_sig_mc_fullClosure", dir.Data (), cent), "");
      rooUnfResp_jet_trk_pt_sig_halfClosure[iDir][iFile]  = new RooUnfoldResponse (h2_jet_trk_pt_sig_halfClosure[iDir][iFile][0], h2_jet_trk_pt_sig_halfClosure[iDir][iFile][1], Form ("rooUnfResp_jet_trk_pt_%s_%s_sig_mc_halfClosure", dir.Data (), cent), "");

    } // end loop over iFile

  } // end loop over iDir




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // NOW LOOP OVER TTREE AND FILL RESPONSE MATRICES
  //////////////////////////////////////////////////////////////////////////////////////////////////// 

  // creates a file identifier pattern that we will use to identify the directory containing the input files
  TString fileIdentifier;
  if (TString (inFileName) == "") {
    if (IsCollisions ()) {
      fileIdentifier = to_string (0);
    }
    else {
      std::cout << "Error: In JetHadronSkimmer.C: Cannot identify this MC file! Quitting." << std::endl;
      return false;
    }
  }
  else
    fileIdentifier = inFileName;


  // opens a TTree as a TChain from all files in a directory matching the file identifierK
  TChain* tree = new TChain ("bush", "bush");
  {
    TString pattern = "*.root";
    std::cout << "DataPath = " << dataPath << std::endl;
    if (TString (inFileName).EndsWith (".root")) 
      tree->Add (dataPath + directory + "/" + inFileName);

    else {
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
    }
    std::cout << "Chain has " << tree->GetListOfFiles ()->GetEntries () << " files, " << tree->GetEntries () << " entries" << std::endl;
  }


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
  //tree->SetBranchAddress ("vert_x",    &vert_x);
  //tree->SetBranchAddress ("vert_y",    &vert_y);
  tree->SetBranchAddress ("vert_z",    &vert_z);
  //tree->SetBranchAddress ("vert_ntrk", &vert_ntrk);
  tree->SetBranchAddress ("vert_type", &vert_type);


  tree->SetBranchAddress ("fcalA_et",       &fcalA_et);
  tree->SetBranchAddress ("fcalC_et",       &fcalC_et);

  tree->SetBranchAddress ("ntrk",                   &trk_n);
  tree->SetBranchAddress ("trk_pt",                 &trk_pt);
  tree->SetBranchAddress ("trk_eta",                &trk_eta);
  tree->SetBranchAddress ("trk_phi",                &trk_phi);
  tree->SetBranchAddress ("trk_charge",             &trk_charge);
  tree->SetBranchAddress ("trk_TightPrimary",       &trk_TightPrimary);
  tree->SetBranchAddress ("trk_HITight",            &trk_HITight);
  tree->SetBranchAddress ("trk_HILoose",            &trk_HILoose);
  tree->SetBranchAddress ("trk_prob_truth",         &trk_prob_truth);
  tree->SetBranchAddress ("trk_truth_pt",           &trk_truth_pt);
  tree->SetBranchAddress ("trk_truth_eta",          &trk_truth_eta);
  tree->SetBranchAddress ("trk_truth_phi",          &trk_truth_phi);
  tree->SetBranchAddress ("trk_truth_charge",       &trk_truth_charge);
  tree->SetBranchAddress ("trk_truth_type",         &trk_truth_type);
  tree->SetBranchAddress ("trk_truth_orig",         &trk_truth_orig);
  tree->SetBranchAddress ("trk_truth_barcode",      &trk_truth_barcode);
  tree->SetBranchAddress ("trk_truth_pdgid",        &trk_truth_pdgid);
  tree->SetBranchAddress ("trk_truth_vz",           &trk_truth_vz);
  tree->SetBranchAddress ("trk_truth_nIn",          &trk_truth_nIn);
  tree->SetBranchAddress ("trk_truth_isHadron",     &trk_truth_isHadron);

  tree->SetBranchAddress ("truth_trk_n",            &truth_trk_n);
  tree->SetBranchAddress ("truth_trk_pt",           &truth_trk_pt);
  tree->SetBranchAddress ("truth_trk_eta",          &truth_trk_eta);
  tree->SetBranchAddress ("truth_trk_phi",          &truth_trk_phi);
  tree->SetBranchAddress ("truth_trk_charge",       &truth_trk_charge);
  tree->SetBranchAddress ("truth_trk_pdgid",        &truth_trk_pdgid);
  tree->SetBranchAddress ("truth_trk_barcode",      &truth_trk_barcode);
  tree->SetBranchAddress ("truth_trk_isHadron",     &truth_trk_isHadron);

  tree->SetBranchAddress ("akt4_hi_jet_n",          &akt4_hi_jet_n);
  tree->SetBranchAddress ("akt4_hi_jet_pt_etajes",  &akt4_hi_jet_pt_etajes);
  tree->SetBranchAddress ("akt4_hi_jet_pt_xcalib",  &akt4_hi_jet_pt_xcalib);
  tree->SetBranchAddress ("akt4_hi_jet_eta_etajes", &akt4_hi_jet_eta_etajes);
  tree->SetBranchAddress ("akt4_hi_jet_eta_xcalib", &akt4_hi_jet_eta_xcalib);
  tree->SetBranchAddress ("akt4_hi_jet_phi",        &akt4_hi_jet_phi);
  tree->SetBranchAddress ("akt4_hi_jet_e_etajes",   &akt4_hi_jet_e_etajes);
  tree->SetBranchAddress ("akt4_hi_jet_e_xcalib",   &akt4_hi_jet_e_xcalib);
  tree->SetBranchAddress ("akt4_hi_jet_timing",     &akt4_hi_jet_timing);

  tree->SetBranchAddress ("akt4_truth_jet_n",       &akt4_truth_jet_n);
  tree->SetBranchAddress ("akt4_truth_jet_pt",      &akt4_truth_jet_pt);
  tree->SetBranchAddress ("akt4_truth_jet_eta",     &akt4_truth_jet_eta);
  tree->SetBranchAddress ("akt4_truth_jet_phi",     &akt4_truth_jet_phi);
  tree->SetBranchAddress ("akt4_truth_jet_e",       &akt4_truth_jet_e);


  const int nJESVar = GetNJESVar ();
  if (nJESVar != -1) {
    std::cout << "Info: In RunCorrelator.cxx: Branching JES variation " << nJESVar << std::endl;
    tree->SetBranchAddress (Form ("akt4_hi_jet_pt_sys_JES_%i", nJESVar), akt4_hi_jet_pt_sys_JES_ALL[nJESVar]);
  }


  const int nTrkWPVar = (DoHITightVar () ? 1 : (DoHILooseVar () ? 2 : 0));


  // for jet energy resolution cuts
  TF1* f_jer = LoadJetEnergyResFunction ();


  // for tracking corections
  TH2D** h2_trk_eff = LoadTrackingEfficiency ();
  TGAE** g_trk_pur = LoadTrackingPurity (!DoJetPrimFracVar ());
  TGAE** gf_trk_pur = LoadTrackingPurityFuncs (!DoJetPrimFracVar ());


  // for MC reweighting to data
  TF1** f_jet_wgts = LoadJetPtWgtFuncs ();
  TH1D** h_jet_wgts = LoadJetPtWgtHists ();


  // for corrected FCal Et values in data overlay
  std::map <const unsigned int, float>* m_overlay_fcalet = (IsDataOverlay () ? GetOverlayFCalMap () : nullptr);


  const long nEvts = tree->GetEntries ();

  // Which jet radius to use. We use R=0.4. Could also plausibly use R=0.2.
  const JetRadius r0p4 = JetRadius::R0p4;


  // minimum and maximum barcodes, 0 < barcode < 10000 in Pythia and 10000 < barcode < 200000 in Hijing
  const int minBarcode = 0;
  const int maxBarcode = 200000;


  const float truth_jet_min_pt = GetJZXR04MinPt (TString (inFileName));
  const float truth_jet_max_pt = GetJZXR04MaxPt (TString (inFileName));
  if (truth_jet_min_pt != 0)
    std::cout << "Checking for leading truth jet with pT > " << truth_jet_min_pt << std::endl;
  if (truth_jet_max_pt != FLT_MAX)
    std::cout << "Checking for leading truth jet with pT < " << truth_jet_max_pt << std::endl;



  ////////////////////////////////////////////////////////////////////////////////////////////////////  
  // Main loop over events
  ////////////////////////////////////////////////////////////////////////////////////////////////////  
  for (long iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)                                
      std::cout << iEvt / (nEvts / 100) << "\% done...\r" << std::flush;

    tree->GetEntry (iEvt % tree->GetEntries ());


    ////////////////////////////////////////////////////////////////////////////////////////////////////  
    // Event selection -- vertexing, truth pT filter, etc.
    ////////////////////////////////////////////////////////////////////////////////////////////////////  
    // vertexing cuts, require no pileup vertices and primary vertex with |vz| < 150mm
    float vz = -999;
    {
      bool hasPrimary = false;
      bool hasPileup = false;
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


    // MC only -- filter events in sample based on min/max of pThat range
    // also sets the appropriate JZ weight
    if (!IsCollisions ()) {
      int iLTJ = -1;
      const int nTJ = GetAktTruthJetN (r0p4);
      for (int iTJ = 0; iTJ < nTJ; iTJ++) {
        if (iLTJ == -1 || GetAktTruthJetPt (iTJ, r0p4) > GetAktTruthJetPt (iLTJ, r0p4))
          iLTJ = iTJ;
      }

      if (iLTJ == -1 || GetAktTruthJetPt (iLTJ, r0p4) < truth_jet_min_pt || GetAktTruthJetPt (iLTJ, r0p4) > truth_jet_max_pt)
        continue;

      event_weight = mcEventWeights->at (0) * crossSectionPicoBarns * mcFilterEfficiency * GetJetLuminosity () / mcNumberEvents; // sigma * f * L_int
    }


    // MC only -- get pure overlay A-side FCal Sum ET value to avoid truth event bias
    if (IspPb () && IsDataOverlay ())
      fcalA_et = m_overlay_fcalet->at (event_number);


    // p+Pb only -- split events by centrality
    short iCent = 0;
    if (IspPb ()) {
      iCent = GetBin (fcalCentBins, nFcalCentBins, fcalA_et);
      if (iCent < 0 || nFcalCentBins <= iCent)
        continue;
    }
    const short iFile = iCent;


    // event weights -- these are not 1 in MC due to JZ weighting (they are 1 in data)
    const double ewgt = event_weight;
    if (ewgt <= 0.)
      continue;


    ////////////////////////////////////////////////////////////////////////////////////////////////////  
    // Fill response matrices
    ////////////////////////////////////////////////////////////////////////////////////////////////////  
    // System boost in p+Pb or pp -- requires a run number to be specified to check for pPb vs Pbp
    const float yboost = GetBoost (run_number);


    //// Get multiplicity bin for event, must be done after mixing (important)!
    //short iMult = 0;
    //while (iMult < nMultBins-1 && multBins[iMult+1] < trk_n)
    //  iMult++;
    //if (nMultBins < iMult)
    //  continue;
    const short iMult = nMultBins-1;


    // loop over all jets again in the event but now correlate reco-matched truth tracks and fill the 2D response matrix
    for (int iTJet = 0; iTJet < GetAktTruthJetN (r0p4); iTJet++) {

      if (!MeetsTruthJetAcceptanceCuts (iTJet, r0p4))
        continue; // if jet is outside of acceptance just ignore it
     
 
      const float tjpt  = GetAktTruthJetPt  (iTJet, r0p4);
      const float tjeta = GetAktTruthJetEta (iTJet, r0p4);
      const float tjphi = GetAktTruthJetPhi (iTJet, r0p4);


      const int iRJet = GetAktHIJetMatch (iTJet, r0p4, nJESVar);
      const float rjpt  = (iRJet < 0 ? 0. : GetAktHIJetPt  (iRJet, r0p4, nJESVar));


      const float jwgt = GetAktJetWeight (tjpt, tjeta, tjphi, r0p4);
      if (jwgt <= 0)
        continue; // sanity check


      if (std::fabs (rjpt/tjpt - 1) > 3*0.01*f_jer->Eval (tjpt))
        continue; // cut on jets reconstructed well outside the JER -- these are bad matches in overlay


      bool isReconstructed = (iRJet >= 0);
      isReconstructed &= MeetsJetAcceptanceCuts (iRJet, r0p4, nJESVar); // reco jet must be accepted too
      //isReconstructed &= std::fabs (rjpt/tjpt - 1) < 3*0.01*f_jer->Eval (tjpt); // cut on jets reconstructed well outside the JER
      isReconstructed &= (pTJBins[0] < rjpt && rjpt < pTJBins[nPtJBins]);

      if (isReconstructed) {
        rooUnfResp_jet_pt_wgts[iFile]->Fill (rjpt, tjpt, ewgt*(iRJet < 0 ? 0. : f_jet_wgts[iFile]->Eval (rjpt)));
        rooUnfResp_jet_pt_altwgts[iFile]->Fill (rjpt, tjpt, ewgt*(iRJet < 0 ? 0. : h_jet_wgts[iFile]->GetBinContent (h_jet_wgts[iFile]->FindBin (rjpt))));
        rooUnfResp_jet_pt_fullClosure[iFile]->Fill (rjpt, tjpt, ewgt);
        if (iEvt % 2 == 0)
          rooUnfResp_jet_pt_halfClosure[iFile]->Fill (rjpt, tjpt, ewgt);

        // fill centrality-integrated response too
        if (IspPb ()) {
          rooUnfResp_jet_pt_wgts[nFiles-1]->Fill (rjpt, tjpt, ewgt*(iRJet < 0 ? 0. : f_jet_wgts[nFiles-1]->Eval (rjpt)));
          rooUnfResp_jet_pt_altwgts[nFiles-1]->Fill (rjpt, tjpt, ewgt*(iRJet < 0 ? 0. : h_jet_wgts[nFiles-1]->GetBinContent (h_jet_wgts[nFiles-1]->FindBin (rjpt))));
          rooUnfResp_jet_pt_fullClosure[nFiles-1]->Fill (rjpt, tjpt, ewgt);
          if (iEvt % 2 == 0)
            rooUnfResp_jet_pt_halfClosure[nFiles-1]->Fill (rjpt, tjpt, ewgt);
        }
      }


      // correlate charged particles with this jet  
      for (int iTrk = 0; iTrk < trk_n; iTrk++) {

        const bool isTruthMatched = (trk_prob_truth[iTrk] > 0.5);

        const bool isFake = !isTruthMatched;

        const bool isSecondary = isTruthMatched && (trk_truth_barcode[iTrk] <= minBarcode || maxBarcode <= trk_truth_barcode[iTrk]);

        const bool isPrimary = !isFake && !isSecondary;

        if (!isPrimary)
          continue;

        if (!MeetsTrackCuts (iTrk, nTrkWPVar))
          continue; // cut on bad quality tracks

        const float trk_y = trk_eta[iTrk];
        if (std::fabs (trk_y - yboost) > 2.035) // 2.035 = 2.5 - 0.465
          continue; // rapidity acceptance cut

        short iDir = -1;
        const float dphi = DeltaPhi (tjphi, trk_phi[iTrk]);
        if (dphi < M_PI/8.)
          iDir = 0;
        else if (M_PI/3. < dphi && dphi < 2.*M_PI/3.)
          iDir = 1;
        else if (dphi > 7.*M_PI/8.)
          iDir = 2;
        if (iDir == -1)
          continue;
        
        if (isReconstructed) {
          rooUnfResp_jet_trk_pt_sig_wgts[iDir][iFile]->Fill (trk_pt[iTrk], rjpt, trk_truth_pt[iTrk], tjpt, ewgt*(iRJet < 0 ? 0. : f_jet_wgts[iFile]->Eval (rjpt)));
          rooUnfResp_jet_trk_pt_sig_altwgts[iDir][iFile]->Fill (trk_pt[iTrk], rjpt, trk_truth_pt[iTrk], tjpt, ewgt*(iRJet < 0 ? 0. : h_jet_wgts[iFile]->GetBinContent (h_jet_wgts[iFile]->FindBin (rjpt))));
          rooUnfResp_jet_trk_pt_sig_fullClosure[iDir][iFile]->Fill (trk_pt[iTrk], rjpt, trk_truth_pt[iTrk], tjpt, ewgt);
          if (iEvt % 2 == 0)
            rooUnfResp_jet_trk_pt_sig_halfClosure[iDir][iFile]->Fill (trk_pt[iTrk], rjpt, trk_truth_pt[iTrk], tjpt, ewgt);
  
          // fill centrality-integrated response too
          if (IspPb ()) {
            rooUnfResp_jet_trk_pt_sig_wgts[iDir][nFiles-1]->Fill (trk_pt[iTrk], rjpt, trk_truth_pt[iTrk], tjpt, ewgt*(iRJet < 0 ? 0. : f_jet_wgts[nFiles-1]->Eval (rjpt)));
            rooUnfResp_jet_trk_pt_sig_altwgts[iDir][nFiles-1]->Fill (trk_pt[iTrk], rjpt, trk_truth_pt[iTrk], tjpt, ewgt*(iRJet < 0 ? 0. : h_jet_wgts[nFiles-1]->GetBinContent (h_jet_wgts[nFiles-1]->FindBin (rjpt))));
            rooUnfResp_jet_trk_pt_sig_fullClosure[iDir][nFiles-1]->Fill (trk_pt[iTrk], rjpt, trk_truth_pt[iTrk], tjpt, ewgt);
            if (iEvt % 2 == 0)
              rooUnfResp_jet_trk_pt_sig_halfClosure[iDir][nFiles-1]->Fill (trk_pt[iTrk], rjpt, trk_truth_pt[iTrk], tjpt, ewgt);
          }
        }

      } // end loop over tracks


      // fill truth jet pT spectrum
      if (isReconstructed) {
        h_jet_pt_wgts[iFile][1]->Fill (tjpt, ewgt*(iRJet < 0 ? 0. : f_jet_wgts[iFile]->Eval (rjpt)));
        h_jet_pt_fullClosure[iFile][1]->Fill (tjpt, ewgt);
        if (iEvt % 2 == 1)
          h_jet_pt_halfClosure[iFile][1]->Fill (tjpt, ewgt);

        if (IspPb ()) {
          h_jet_pt_wgts[nFiles-1][1]->Fill (tjpt, ewgt*(iRJet < 0 ? 0. : f_jet_wgts[nFiles-1]->Eval (rjpt)));
          h_jet_pt_fullClosure[nFiles-1][1]->Fill (tjpt, ewgt);
          if (iEvt % 2 == 1)
            h_jet_pt_halfClosure[nFiles-1][1]->Fill (tjpt, ewgt);
        }
      }


      // correlate truth charged particles with this jet  
      for (int iTTrk = 0; iTTrk < truth_trk_n; iTTrk++) {

        if (truth_trk_charge[iTTrk] == 0)
          continue; // cut on neutrals

        if (!truth_trk_isHadron[iTTrk])
          continue; // cut on non-hadrons

        const bool isSecondary = (truth_trk_barcode[iTTrk] <= minBarcode || maxBarcode <= truth_trk_barcode[iTTrk]);

        const bool isPrimary = !isSecondary;

        if (!isPrimary)
          continue;

        const float truth_trk_y = truth_trk_eta[iTTrk];
        if (std::fabs (truth_trk_y - yboost) > 2.035) // 2.035 = 2.5 - 0.465
          continue; // rapidity acceptance cut

        short iDir = -1;
        const float dphi = DeltaPhi (tjphi, truth_trk_phi[iTTrk]);
        if (dphi < M_PI/8.)
          iDir = 0;
        else if (M_PI/3. < dphi && dphi < 2.*M_PI/3.)
          iDir = 1;
        else if (dphi > 7.*M_PI/8.)
          iDir = 2;
        if (iDir == -1)
          continue;

        // fill truth jet FF plots
        h2_jet_trk_pt_sig_wgts[iDir][iFile][1]->Fill (truth_trk_pt[iTTrk], tjpt, ewgt*(iRJet < 0 ? 0. : f_jet_wgts[iFile]->Eval (rjpt)));
        h2_jet_trk_pt_sig_fullClosure[iDir][iFile][1]->Fill (truth_trk_pt[iTTrk], tjpt, ewgt);
        if (iEvt % 2 == 1)
          h2_jet_trk_pt_sig_halfClosure[iDir][iFile][1]->Fill (truth_trk_pt[iTTrk], tjpt, ewgt);

        if (IspPb ()) {
          h2_jet_trk_pt_sig_wgts[iDir][nFiles-1][1]->Fill (truth_trk_pt[iTTrk], tjpt, ewgt*(iRJet < 0 ? 0. : f_jet_wgts[nFiles-1]->Eval (rjpt)));
          h2_jet_trk_pt_sig_fullClosure[iDir][nFiles-1][1]->Fill (truth_trk_pt[iTTrk], tjpt, ewgt);
          if (iEvt % 2 == 1)
            h2_jet_trk_pt_sig_halfClosure[iDir][nFiles-1][1]->Fill (truth_trk_pt[iTTrk], tjpt, ewgt);
        }

      } // end loop over truth tracks

    } // end loop over truth jets




    ////////////////////////////////////////////////////////////////////////////////////////////////////  
    // Fill histograms at reco. level for unfolding tests
    ////////////////////////////////////////////////////////////////////////////////////////////////////  
    for (int iRJet = 0; iRJet < GetAktHIJetN (r0p4); iRJet++) {

      if (!MeetsJetAcceptanceCuts (iRJet, r0p4, nJESVar))
        continue; // jet eta/phi & timing cuts


      const float rjpt  = GetAktHIJetPt  (iRJet, r0p4, nJESVar);
      const float rjeta = GetAktHIJetEta (iRJet, r0p4, nJESVar);
      const float rjphi = GetAktHIJetPhi (iRJet, r0p4, nJESVar);


      if (rjpt < pTJBins[0] || pTJBins[nPtJBins] < rjpt)
        continue;


      const float jwgt = GetAktJetWeight (rjpt, rjeta, rjphi, r0p4);
      if (jwgt <= 0)
        continue; // sanity check


      // In MC, check that the jet is truth-matched. Otherwise problems can emerge, particularly in the overlay.
      const int iTJet = GetAktTruthJetMatch (iRJet, r0p4);
      if (iTJet == -1)
        continue; // in case of no truth jet match
      const float tjpt = GetAktTruthJetPt (iTJet, r0p4);
      if (std::fabs (rjpt/tjpt - 1) > 3*0.01*f_jer->Eval (tjpt))
        continue; // cut on jets reconstructed well outside the JER -- these are bad matches in overlay

      if (tjpt < pTJBins[0] || pTJBins[nPtJBins] < tjpt)
        continue;


      // fill reco jet spectrum
      h_jet_pt_wgts[iFile][0]->Fill (rjpt, ewgt*(iRJet < 0 ? 0. : f_jet_wgts[iFile]->Eval (rjpt)));
      h_jet_pt_fullClosure[iFile][0]->Fill (rjpt, ewgt);
      if (iEvt % 2 == 1)
        h_jet_pt_halfClosure[iFile][0]->Fill (rjpt, ewgt);
      if (IspPb ()) {
        h_jet_pt_wgts[nFiles-1][0]->Fill (rjpt, ewgt*(iRJet < 0 ? 0. : f_jet_wgts[nFiles-1]->Eval (rjpt)));
        h_jet_pt_fullClosure[nFiles-1][0]->Fill (rjpt, ewgt);
        if (iEvt % 2 == 1)
          h_jet_pt_halfClosure[nFiles-1][0]->Fill (rjpt, ewgt);
      }


      // correlate charged particles with this jet  
      for (int iTrk = 0; iTrk < trk_n; iTrk++) {

        //const bool isTruthMatched = (trk_prob_truth[iTrk] > 0.5);

        //const bool isFake = !isTruthMatched;

        //const bool isSecondary = isTruthMatched && (trk_truth_barcode[iTrk] <= minBarcode || maxBarcode <= trk_truth_barcode[iTrk]);

        //const bool isPrimary = !isFake && !isSecondary;

        //if (!isPrimary)
        //  continue;

        if (!MeetsTrackCuts (iTrk, nTrkWPVar))
          continue; // cut on bad quality tracks

        if (trk_pt[iTrk] < pTChBins[0])
          continue; // histogram pT cut on tracks

        const float trk_y = trk_eta[iTrk];
        if (std::fabs (trk_y - yboost) > 2.035) // 2.035 = 2.5 - 0.465
          continue; // rapidity acceptance cut

        short iDir = -1;
        const float dphi = DeltaPhi (rjphi, trk_phi[iTrk]);
        if (dphi < M_PI/8.)
          iDir = 0;
        else if (M_PI/3. < dphi && dphi < 2.*M_PI/3.)
          iDir = 1;
        else if (dphi > 7.*M_PI/8.)
          iDir = 2;
        if (iDir == -1)
          continue;

        short iEta = 0;
        while (iEta < nEtaTrkBins && etaTrkBins[iEta+1] < std::fabs (trk_eta[iTrk]))
          iEta++;
        if (iEta == nEtaTrkBins)
          continue;

        const float teff = h2_trk_eff[iMult]->GetBinContent (h2_trk_eff[iMult]->FindBin (trk_eta[iTrk], trk_pt[iTrk]));
        const float tpur = (DoPrimFitVar () ? g_trk_pur[iEta]->Eval (trk_pt[iTrk]) : gf_trk_pur[iEta]->Eval (trk_pt[iTrk]));
        //const float tpur = 1;
        const float twgt = (teff > 0. ? tpur / teff : 0.);

        // fill reco jet FF plots
        h2_jet_trk_pt_sig_wgts[iDir][iFile][0]->Fill (trk_pt[iTrk], rjpt, ewgt*twgt*(iRJet < 0 ? 0. : f_jet_wgts[iFile]->Eval (rjpt)));
        h2_jet_trk_pt_sig_fullClosure[iDir][iFile][0]->Fill (trk_pt[iTrk], rjpt, ewgt*twgt);
        if (iEvt % 2 == 1)
          h2_jet_trk_pt_sig_halfClosure[iDir][iFile][0]->Fill (trk_pt[iTrk], rjpt, ewgt*twgt);

        if (IspPb ()) {
          h2_jet_trk_pt_sig_wgts[iDir][nFiles-1][0]->Fill (trk_pt[iTrk], rjpt, ewgt*twgt*(iRJet < 0 ? 0. : f_jet_wgts[nFiles-1]->Eval (rjpt)));
          h2_jet_trk_pt_sig_fullClosure[iDir][nFiles-1][0]->Fill (trk_pt[iTrk], rjpt, ewgt*twgt);
          if (iEvt % 2 == 1)
            h2_jet_trk_pt_sig_halfClosure[iDir][nFiles-1][0]->Fill (trk_pt[iTrk], rjpt, ewgt*twgt);
        }

      } // end loop over tracks

    } // end loop over reco jets


  } // end loop over events
  std::cout << "Finished event loop." << std::endl;



  SaferDelete (&f_jer);

  for (int iMult = 0; iMult < nMultBins; iMult++) {
    SaferDelete (&h2_trk_eff[iMult]);
  }
  Delete1DArray (h2_trk_eff, nMultBins);

  for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {
    SaferDelete (&g_trk_pur[iEta]);
    SaferDelete (&gf_trk_pur[iEta]);
  }
  Delete1DArray (g_trk_pur, nEtaTrkBins);
  Delete1DArray (gf_trk_pur, nEtaTrkBins);

  for (int iFile = 0; iFile < nFiles; iFile++) {
    SaferDelete (&f_jet_wgts[iFile]);
  }
  Delete1DArray (f_jet_wgts, nFiles);


  SaferDelete (&m_overlay_fcalet);
 



  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // FINALLY WRITE OUT EVERYTHING TO A SINGLE ROOT FILE WITH ALL THE RESULTS.
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  {
    outFile->cd ();

    for (short iFile = 0; iFile < nFiles; iFile++) {

      rooUnfResp_jet_pt_wgts[iFile]->Write ();
      rooUnfResp_jet_pt_altwgts[iFile]->Write ();
      rooUnfResp_jet_pt_fullClosure[iFile]->Write ();
      rooUnfResp_jet_pt_halfClosure[iFile]->Write ();

      h_jet_pt_wgts[iFile][0]->Write ();
      h_jet_pt_wgts[iFile][1]->Write ();
      h_jet_pt_fullClosure[iFile][0]->Write ();
      h_jet_pt_fullClosure[iFile][1]->Write ();
      h_jet_pt_halfClosure[iFile][0]->Write ();
      h_jet_pt_halfClosure[iFile][1]->Write ();

      for (short iDir = 0; iDir < nDir; iDir++) {
  
        rooUnfResp_jet_trk_pt_sig_wgts[iDir][iFile]->Write ();
        rooUnfResp_jet_trk_pt_sig_altwgts[iDir][iFile]->Write ();
        rooUnfResp_jet_trk_pt_sig_fullClosure[iDir][iFile]->Write ();
        rooUnfResp_jet_trk_pt_sig_halfClosure[iDir][iFile]->Write ();

        h2_jet_trk_pt_sig_wgts[iDir][iFile][0]->Write ();
        h2_jet_trk_pt_sig_wgts[iDir][iFile][1]->Write ();
        h2_jet_trk_pt_sig_fullClosure[iDir][iFile][0]->Write ();
        h2_jet_trk_pt_sig_fullClosure[iDir][iFile][1]->Write ();
        h2_jet_trk_pt_sig_halfClosure[iDir][iFile][0]->Write ();
        h2_jet_trk_pt_sig_halfClosure[iDir][iFile][1]->Write ();

      } // end loop over iDir

    } // end loop over iFile

    outFile->Close ();
  }


  return true;
}


int main (int argc, char** argv) {

  int argn                = 1;
  const char* subdir      = argv[argn++];

  TString collSys = TString (argv[argn++]);
  if (!SetCollisionSystem (collSys)) {
    std::cout << "In RunCorrelator.cxx: Invalid collision system, exiting." << std::endl;
    return 1;
  }

  TString dType = TString (argv[argn++]);
  if (!SetDataType (dType)) {
    std::cout << "In RunCorrelator.cxx: Invalid data type, exiting." << std::endl;
    return 2;
  }

  TString sFlag = TString (argv[argn++]);
  if (!SetSystFlag (sFlag)) {
    std::cout << "In RunCorrelator.cxx: Invalid systematic flag, exiting." << std::endl;
    return 3;
  }

  const char* inFileName            = (argc > argn && argv[argn] ? argv[argn++] : "");
  const char* outFileName           = (argc > argn && argv[argn] ? argv[argn++] : "");
  if (!IsCollisions ())
    GetMCWeights (inFileName);

  std::cout << "Configuration set to";
  std::cout << "\n\tsubdir = " << subdir;
  std::cout << "\n\tCollisionSystem = " << ToTString (collisionSystem);
  std::cout << "\n\tDataType = " << ToTString (dataType);
  std::cout << "\n\tTriggerType = " << ToTString (triggerType);
  std::cout << "\n\tSystFlag = " << ToTString (systFlag);
  std::cout << "\n\tinFileName = " << inFileName;
  if (crossSectionPicoBarns != 0.)
    std::cout << "\n\t  --> deduced crossSectionPicoBarns = " << crossSectionPicoBarns;
  if (mcFilterEfficiency != 0.)
    std::cout << "\n\t  --> deduced mcFilterEfficiency = " << mcFilterEfficiency;
  if (mcNumberEvents != 0.)
    std::cout << "\n\t  --> deduced mcNumberEvents = " << mcNumberEvents;
  std::cout << std::endl;


  std::cout << "Running MakeResponseMatrix algorithm..." << std::endl;
  bool success = MakeResponseMatrix (subdir, inFileName, outFileName);

  delete [] zdcCentBins;
  zdcCentBins = nullptr;
  delete [] fcalCentBins;
  fcalCentBins = nullptr;
  delete [] fcalMixBins;
  fcalMixBins = nullptr;
  delete [] fineFcalCentBins;
  fineFcalCentBins = nullptr;
  delete [] ppMixBins;
  ppMixBins = nullptr;

  if (success) {
    std::cout << "Finished algorithm!" << std::endl;
    return 0;
  }
  else {
    std::cout << "Algorithm failed!" << std::endl;
    return 1;
  }
}


#endif
