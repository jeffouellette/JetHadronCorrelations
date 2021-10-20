#ifndef __RunCorrelator_C__
#define __RunCorrelator_C__

#include "Params.h"
#include "CentralityDefs.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Trigger.h"
#include "Process.h"

#include <Utilities.h>
#include <ArrayTemplates.h>

#include <TChain.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TPRegexp.h>
#include <TRandom3.h>

#include <iostream>
#include <math.h>
#include <map>
#include <chrono>


using namespace JetHadronCorrelations;


bool doMixing = false;
bool sameTreeMixing = false;

float truth_jet_min_pt = 0;
float truth_jet_max_pt = FLT_MAX;



// Main function for performing jet-hadron correlations. (Wrapper function available below.)
bool Correlator (const char* tag, const char* outFilePattern, TTree* jetsTree, TTree* tracksTree = nullptr) {

  auto start = std::chrono::steady_clock::now (); // time at beginning of correlator function.


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
  Trigger* jetTrigger = nullptr;
  Trigger* trackTreeTrigger = nullptr;

  // assign branches for triggers
  if (IsCollisions ()) {
    if (UseJ50Triggers ()) {
      jetTrigger = new Trigger (jet_trig_name[0]);
      std::cout << "Info: In RunCorrelator.cxx: Looking for " << jet_trig_name[0] << " trigger" << std::endl;
      jetsTree->SetBranchAddress ((jet_trig_name[0]+"_decision").c_str (), &(jetTrigger->trigDecision));
      jetsTree->SetBranchAddress ((jet_trig_name[0]+"_prescale").c_str (), &(jetTrigger->trigPrescale));
    }
    else if (UseJ100Triggers ()) {
      jetTrigger = new Trigger (jet_trig_name[1]);
      std::cout << "Info: In RunCorrelator.cxx: Looking for " << jet_trig_name[1] << " trigger" << std::endl;
      jetsTree->SetBranchAddress ((jet_trig_name[1]+"_decision").c_str (), &(jetTrigger->trigDecision));
      jetsTree->SetBranchAddress ((jet_trig_name[1]+"_prescale").c_str (), &(jetTrigger->trigPrescale));
    }
    else if (!UseJetTriggers ()) {
      jetTrigger = new Trigger (minbias_trig_name[0]);
      std::cout << "Info: In RunCorrelator.cxx: Looking for " << minbias_trig_name[0] << " trigger" << std::endl;
      jetsTree->SetBranchAddress ((minbias_trig_name[0]+"_decision").c_str (), &(jetTrigger->trigDecision));
      jetsTree->SetBranchAddress ((minbias_trig_name[0]+"_prescale").c_str (), &(jetTrigger->trigPrescale));
    }
    else {
      std::cout << "Error: In RunCorrelator.cxx: Invalid trigger scheme? Please debug! Exiting." << std::endl;
      return false;
    }
  }

  if (doMixing) {
    trackTreeTrigger = new Trigger (minbias_trig_name[0]);
    tracksTree->SetBranchAddress ((minbias_trig_name[0]+"_decision").c_str (), &(trackTreeTrigger->trigDecision));
    tracksTree->SetBranchAddress ((minbias_trig_name[0]+"_prescale").c_str (), &(trackTreeTrigger->trigPrescale));
  }


  // setup branch addresses
  jetsTree->SetBranchAddress ("event_number",  &event_number);
  jetsTree->SetBranchAddress ("lumi_block",    &lumi_block);
  jetsTree->SetBranchAddress ("run_number",    &run_number);


  jetsTree->SetBranchAddress ("actualInteractionsPerCrossing",  &actualInteractionsPerCrossing);
  jetsTree->SetBranchAddress ("averageInteractionsPerCrossing", &averageInteractionsPerCrossing);


  if (!IsCollisions ())
    jetsTree->SetBranchAddress ("mcEventWeights", &mcEventWeights);

  jetsTree->SetBranchAddress ("nvert",     &nvert);
  //jetsTree->SetBranchAddress ("vert_x",    &vert_x);
  //jetsTree->SetBranchAddress ("vert_y",    &vert_y);
  jetsTree->SetBranchAddress ("vert_z",    &vert_z);
  //jetsTree->SetBranchAddress ("vert_ntrk", &vert_ntrk);
  jetsTree->SetBranchAddress ("vert_type", &vert_type);

  jetsTree->SetBranchAddress ("fcalA_et",         &fcalA_et);
  jetsTree->SetBranchAddress ("fcalC_et",         &fcalC_et);
  jetsTree->SetBranchAddress ("fcalA_et_Cos2",    &fcalA_et_Cos2);
  jetsTree->SetBranchAddress ("fcalC_et_Cos2",    &fcalC_et_Cos2);
  jetsTree->SetBranchAddress ("fcalA_et_Sin2",    &fcalA_et_Sin2);
  jetsTree->SetBranchAddress ("fcalC_et_Sin2",    &fcalC_et_Sin2);

  if (IsCollisions () && !Ispp ()) {
    jetsTree->SetBranchAddress ("ZdcCalibEnergy_A", &ZdcCalibEnergy_A);
    jetsTree->SetBranchAddress ("ZdcCalibEnergy_C", &ZdcCalibEnergy_C);
  }

  // get event matching & track information from the mixed event tree
  if (doMixing) {
    tracksTree->SetBranchAddress ("nvert",     &nvert_matching);
    //tracksTree->SetBranchAddress ("vert_x",    &vert_x_matching);
    //tracksTree->SetBranchAddress ("vert_y",    &vert_y_matching);
    tracksTree->SetBranchAddress ("vert_z",    &vert_z_matching);
    //tracksTree->SetBranchAddress ("vert_ntrk", &vert_ntrk_matching);
    tracksTree->SetBranchAddress ("vert_type", &vert_type_matching);

    tracksTree->SetBranchAddress ("fcalA_et",         &fcalA_et_matching);
    tracksTree->SetBranchAddress ("fcalC_et",         &fcalC_et_matching);
    tracksTree->SetBranchAddress ("fcalA_et_Cos2",    &fcalA_et_Cos2_matching);
    tracksTree->SetBranchAddress ("fcalC_et_Cos2",    &fcalC_et_Cos2_matching);
    tracksTree->SetBranchAddress ("fcalA_et_Sin2",    &fcalA_et_Sin2_matching);
    tracksTree->SetBranchAddress ("fcalC_et_Sin2",    &fcalC_et_Sin2_matching);

    if (IsCollisions () && !Ispp ()) {
      tracksTree->SetBranchAddress ("ZdcCalibEnergy_A", &ZdcCalibEnergy_A_matching);
      tracksTree->SetBranchAddress ("ZdcCalibEnergy_C", &ZdcCalibEnergy_C_matching);
    }
  }

  if (UseTruthParticles ()) {
    if (IsCollisions ()) {
      std::cout << "Error: In RunCorrelator.cxx: Configured to use truth particles but running over data? Please investigate. Exiting gracefully." << std::endl;
      return true;
    }
    tracksTree->SetBranchAddress ("truth_trk_n",            &trk_n);
    tracksTree->SetBranchAddress ("truth_trk_pt",           &trk_pt);
    tracksTree->SetBranchAddress ("truth_trk_eta",          &trk_eta);
    tracksTree->SetBranchAddress ("truth_trk_phi",          &trk_phi);
    tracksTree->SetBranchAddress ("truth_trk_charge",       &trk_charge);
    tracksTree->SetBranchAddress ("truth_trk_barcode",      &trk_truth_barcode);
    tracksTree->SetBranchAddress ("truth_trk_isHadron",     &trk_truth_isHadron);
  }
  else {
    tracksTree->SetBranchAddress ("ntrk",                   &trk_n);
    tracksTree->SetBranchAddress ("trk_pt",                 &trk_pt);
    tracksTree->SetBranchAddress ("trk_eta",                &trk_eta);
    tracksTree->SetBranchAddress ("trk_phi",                &trk_phi);
    tracksTree->SetBranchAddress ("trk_charge",             &trk_charge);
    tracksTree->SetBranchAddress ("trk_TightPrimary",       &trk_TightPrimary);
    tracksTree->SetBranchAddress ("trk_HITight",            &trk_HITight);
    tracksTree->SetBranchAddress ("trk_HILoose",            &trk_HILoose);
    if (!IsCollisions () && !doMixing)
      tracksTree->SetBranchAddress ("trk_prob_truth",       &trk_prob_truth);

    //// get truth particles from MC (i.e. jetsTree, not tracksTree!)
    //if (!IsCollisions ()) {
    //  jetsTree->SetBranchAddress ("truth_trk_n",            &truth_trk_n);
    //  jetsTree->SetBranchAddress ("truth_trk_pt",           &truth_trk_pt);
    //  jetsTree->SetBranchAddress ("truth_trk_eta",          &truth_trk_eta);
    //  jetsTree->SetBranchAddress ("truth_trk_phi",          &truth_trk_phi);
    //  jetsTree->SetBranchAddress ("truth_trk_charge",       &truth_trk_charge);
    //  jetsTree->SetBranchAddress ("truth_trk_pdgid",        &truth_trk_pdgid);
    //  jetsTree->SetBranchAddress ("truth_trk_barcode",      &truth_trk_barcode);
    //  jetsTree->SetBranchAddress ("truth_trk_isHadron",     &truth_trk_isHadron);
    //}
  }

  // get truth jets from MC
  if (!IsCollisions ()) {
    jetsTree->SetBranchAddress ("akt4_truth_jet_n",     &akt4_truth_jet_n);
    jetsTree->SetBranchAddress ("akt4_truth_jet_pt",    &akt4_truth_jet_pt);
    jetsTree->SetBranchAddress ("akt4_truth_jet_eta",   &akt4_truth_jet_eta);
    jetsTree->SetBranchAddress ("akt4_truth_jet_phi",   &akt4_truth_jet_phi);
    jetsTree->SetBranchAddress ("akt4_truth_jet_e",     &akt4_truth_jet_e);
  }

  jetsTree->SetBranchAddress ("akt4_hi_jet_n",            &akt4_hi_jet_n);
  jetsTree->SetBranchAddress ("akt4_hi_jet_pt_etajes",    &akt4_hi_jet_pt_etajes);
  jetsTree->SetBranchAddress ("akt4_hi_jet_pt_xcalib",    &akt4_hi_jet_pt_xcalib);
  jetsTree->SetBranchAddress ("akt4_hi_jet_eta_etajes",   &akt4_hi_jet_eta_etajes);
  jetsTree->SetBranchAddress ("akt4_hi_jet_eta_xcalib",   &akt4_hi_jet_eta_xcalib);
  jetsTree->SetBranchAddress ("akt4_hi_jet_phi",          &akt4_hi_jet_phi);
  jetsTree->SetBranchAddress ("akt4_hi_jet_e_etajes",     &akt4_hi_jet_e_etajes);
  jetsTree->SetBranchAddress ("akt4_hi_jet_e_xcalib",     &akt4_hi_jet_e_xcalib);
  jetsTree->SetBranchAddress ("akt4_hi_jet_timing",       &akt4_hi_jet_timing);


  const short nJESVar = GetNJESVar ();
  if (!IsCollisions () && nJESVar != -1) {
    std::cout << "Info: In RunCorrelator.cxx: Branching JES variation " << nJESVar << std::endl;
    jetsTree->SetBranchAddress (Form ("akt4_hi_jet_pt_sys_JES_%i", nJESVar), akt4_hi_jet_pt_sys_JES_ALL[nJESVar]);
  }


  const short nTrkWPVar = (DoHITightVar () ? 1 : (DoHILooseVar () ? 2 : 0));


  // setup centrality bins (only relevant for p+Pb)
  double* centBins = (DoFcalCentVar () ? fcalCentBins : (DoFineFcalCentVar () ? fineFcalCentBins : (!IsCollisions () ? fcalCentBins : zdcCentBins)));
  const short nCentBins = (DoFcalCentVar () ? nFcalCentBins : (DoFineFcalCentVar () ? nFineFcalCentBins : (!IsCollisions () ? nFcalCentBins : nZdcCentBins)));

  std::cout << "Centrality bin cuts: ";
  for (short iCent = 0; iCent < nCentBins; iCent++)
    std::cout << centBins[iCent] << ", ";
  std::cout << centBins[nCentBins] << std::endl;


  // for jet energy resolution cuts
  TF1* f_jer = LoadJetEnergyResFunction ();


  // for tracking corections
  TH2D** h2_trk_eff = LoadTrackingEfficiency ();
  TH1D** h_trk_pur = LoadTrackingPurity ();
  TF1** f_trk_pur = LoadTrackingPurityFuncs ();


  // for corrected FCal Et values in data overlay
  std::map <const unsigned int, float>* m_overlay_fcalet = (IsDataOverlay () ? GetOverlayFCalMap () : nullptr);


  // create output histograms in output files
  const short nFiles = (Ispp () ? 1 : nCentBins);
  TFile** outFiles = new TFile*[nFiles];

  TH1D**  h_evt_counts    = Get1DArray <TH1D*> (nFiles);
  TH1D*** h_jet_counts    = Get2DArray <TH1D*> (nFiles, nPtJBins);

  TH1D**  h_jet_pt        = Get1DArray <TH1D*> (nFiles);
  TH2D**  h2_jet_pt_cov   = Get1DArray <TH2D*> (nFiles);

  TH2D*** h2_jet_eta_phi  = Get2DArray <TH2D*> (nFiles, nPtJBins);

  TH1D**** h_jet_trk_dphi       = Get3DArray <TH1D*> (nFiles, nPtJBins, nPtChSelections);
  TH2D**** h2_jet_trk_dphi_cov  = Get3DArray <TH2D*> (nFiles, nPtJBins, nPtChSelections);

  TH1D**** h_jet_trk_pt       = Get3DArray <TH1D*> (nFiles, nPtJBins, nDir);
  TH2D**** h2_jet_trk_pt_cov  = Get3DArray <TH2D*> (nFiles, nPtJBins, nDir);


  for (short iFile = 0; iFile < nFiles; iFile++) {

    TString outFileName = outFilePattern;
    outFileName.ReplaceAll ("iCent_", Form ("iCent%i_", iFile));
    TFile* outFile = new TFile (outFileName.Data (), "recreate");
    outFiles[iFile] = outFile;

    h_evt_counts[iFile] = new TH1D (Form ("h_evt_counts_%s", tag),  "", 3, -0.5, 2.5);

    h_jet_pt[iFile] = new TH1D (Form ("h_jet_pt_%s", tag), ";#it{p}_{T}^{jet} [GeV];(1/N_{jet}) (dN_{jet}/d#it{p}_{T}) [GeV^{-1}]", nPtJBins, pTJBins);
    h2_jet_pt_cov[iFile] = new TH2D (Form ("h2_jet_pt_cov_%s", tag), ";#it{p}_{T}^{jet} [GeV];#it{p}_{T}^{jet} [GeV];Covariance", nPtJBins, pTJBins, nPtJBins, pTJBins);

    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

      const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

      h_jet_counts[iFile][iPtJ] = new TH1D (Form ("h_jet_counts_%s_%s", pTJ.Data (), tag),  "", 3, -0.5, 2.5);

      h2_jet_eta_phi[iFile][iPtJ] = new TH2D (Form ("h2_jet_eta_phi_%s_%s", pTJ.Data (), tag), ";#eta^{jet};#phi^{jet};Counts", 64, -3.2, 3.2, 64, -M_PI, M_PI);

      for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

        const TString pTCh = pTChSelections[iPtCh];

        h_jet_trk_dphi[iFile][iPtJ][iPtCh] = new TH1D (Form ("h_jet_trk_dphi_%s_%s_%s", pTCh.Data (), pTJ.Data (), tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
        h_jet_trk_dphi[iFile][iPtJ][iPtCh]->Sumw2 ();

        h2_jet_trk_dphi_cov[iFile][iPtJ][iPtCh] = new TH2D (Form ("h2_jet_trk_dphi_cov_%s_%s_%s", pTCh.Data (), pTJ.Data (), tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
        h2_jet_trk_dphi_cov[iFile][iPtJ][iPtCh]->Sumw2 ();

      } // end loop over iPtCh

      for (short iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        h_jet_trk_pt[iFile][iPtJ][iDir] = new TH1D (Form ("h_jet_trk_pt_%s_%s_%s", dir.Data (), pTJ.Data (), tag), ";#it{p}_{T}^{ch} [GeV];(1/N_{jet}) (dN_{ch}/d#it{p}_{T}) [GeV^{-1}]", nPtChBins, pTChBins);
        h_jet_trk_pt[iFile][iPtJ][iDir]->Sumw2 ();

        h2_jet_trk_pt_cov[iFile][iPtJ][iDir] = new TH2D (Form ("h2_jet_trk_pt_cov_%s_%s_%s", dir.Data (), pTJ.Data (), tag), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{ch} [GeV];Covariance", nPtChBins, pTChBins, nPtChBins, pTChBins);
        h2_jet_trk_pt_cov[iFile][iPtJ][iDir]->Sumw2 ();

      } // end loop over iDir

    } // end loop over iPtJ

  } // end loop over iFile


  // arrays for filling histograms & covariance matrices correctly
  double* jet_pt_counts         = Get1DArray <double> (nPtJBins);

  double*** jet_trk_dphi_counts = Get3DArray <double> (nPtJBins, nPtChSelections, nDPhiBins);
  double*** jet_trk_pt_counts   = Get3DArray <double> (nPtJBins, nDir, nPtChBins);

  //const long nEvts = (doMixing ? 1. : 1.) * (jetsTree->GetEntries ()); // for debugging... TODO remove me!
  const long nEvts = (doMixing ? 20. : 1.) * (jetsTree->GetEntries ());
  const long nTrkEvts = tracksTree->GetEntries ();

  std::cout << "Correlator will process " << nEvts << " events" << std::endl;


  // setup centrality matching bins for event mixing according to job configuration
  double* mixCentBins = nullptr;
  short nMixCentBins = 0;

  if (Ispp ()) {
    if      (DoMixCatVar1 ()) { mixCentBins = ppMixVar1Bins;    nMixCentBins = nppMixVar1Bins;    }
    else if (DoMixCatVar3 ()) { mixCentBins = ppMixVar3Bins;    nMixCentBins = nppMixVar3Bins;    }
    else                      { mixCentBins = ppMixBins;        nMixCentBins = nppMixBins;        }
  }
  else {
    if      (DoMixCatVar2 ()) { mixCentBins = fcalMixVar2Bins;  nMixCentBins = nFcalMixVar2Bins;  }
    //else if (DoMixCatVar6 ()) { mixCentBins = fcalMixVar6Bins;  nMixCentBins = nFcalMixVar6Bins;  }
    else                      { mixCentBins = fcalMixBins;      nMixCentBins = nFcalMixBins;      }
  }


  // counter for mixed events (declared outside of loop scope to keep its value from one loop to the next)
  long iTrkEvt = 0;

  // Which jet radius to use. We use R=0.4. Could also plausibly use R=0.2.
  const JetRadius r0p4 = JetRadius::R0p4;



  ////////////////////////////////////////////////////////////////////////////////////////////////////  
  // Main loop over events
  ////////////////////////////////////////////////////////////////////////////////////////////////////  
  for (long iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)                                
      std::cout << iEvt / (nEvts / 100) << "\% done...\r" << std::flush;

    jetsTree->GetEntry (iEvt % jetsTree->GetEntries ());


    float prescale = 1;


    // triggering cut, require appropriate jet trigger to have fired
    if (IsCollisions ()) {
      if (!jetTrigger->trigDecision)
        continue;
      event_weight = 1; // optionally set to trigger prescale, but not really necessary
      prescale = jetTrigger->trigPrescale;
    }


    // vertexing cuts, require no pileup vertices and primary vertex with |vz| < 150mm
    float vz = -999;
    {
      bool hasPrimary = false;
      bool hasPileup = false;
      for (short iVert = 0; iVert < nvert; iVert++) {
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


    // MC only -- filter events in sample based on min/max of pThat range
    // also sets the appropriate JZ weight
    if (!IsCollisions ()) {
      short iLTJ = -1;
      const short nTJ = GetAktTruthJetN (r0p4);
      for (short iTJ = 0; iTJ < nTJ; iTJ++) {
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
      const float centVar = ((!IsCollisions () || DoFcalCentVar ()  || DoFineFcalCentVar ()) ? fcalA_et : (ZdcCalibEnergy_A * 1e3));
      iCent = GetBin (centBins, nCentBins, centVar);
      if (iCent < 0 || nCentBins <= iCent)
        continue;
    }
    const short iFile = iCent;


    // calculate event plane angle based on FCal q-vector moments
    QnVector Q2 = (IspPb () ? GetPbQ2Vec () : GetProtonQ2Vec ());
    const float psi2 = std::atan2 (Q2.second, Q2.first)/2.;


    // event weights -- these are not 1 in MC due to JZ weighting (they are 1 in data)
    const double ewgt = event_weight * prescale;
    if (ewgt <= 0.)
      continue;


    // initialize all histogramming bins to 0 for this event
    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++)
      jet_pt_counts[iPtJ] = 0;

   
    short nj = 0;

    const short jn = UseTruthJets () ? GetAktTruthJetN (r0p4) :  GetAktHIJetN (r0p4);
    for (short iJet = 0; iJet < jn; iJet++) {

      const float jpt  = (UseTruthJets () ? GetAktTruthJetPt  (iJet, r0p4) : GetAktHIJetPt  (iJet, r0p4, nJESVar));
      const float jeta = (UseTruthJets () ? GetAktTruthJetEta (iJet, r0p4) : GetAktHIJetEta (iJet, r0p4, nJESVar));
      const float jphi = (UseTruthJets () ? GetAktTruthJetPhi (iJet, r0p4) : GetAktHIJetPhi (iJet, r0p4, nJESVar));

      if (!MeetsJetAcceptanceCuts (iJet, r0p4, nJESVar))
        continue; // jet eta/phi & timing cuts
      if (!MeetsJetPtCut (jpt))
        continue; // jet pT cuts, for trigger relevance purposes

      const float thisjwgt = GetAktJetWeight (jpt, jeta, jphi, r0p4);
      if (thisjwgt <= 0.)
        continue; // sanity check

      // in MC, check that the jet is truth-matched. Otherwise problems can emerge, particularly in the overlay.
      if (!IsCollisions () && !UseTruthJets ()) {
        const short iTJet = GetAktTruthJetMatch (iJet, r0p4);
        if (iTJet == -1)
          continue; // in case of no truth jet match
        const float tjpt = GetAktTruthJetPt (iTJet, r0p4);
        if (std::fabs (jpt/tjpt - 1) > 3*0.01*f_jer->Eval (tjpt))
          continue; // cut on jets reconstructed well outside the JER
      }

      const short iPtJ = GetPtJBin (jpt);
      if (0 <= iPtJ && iPtJ < nPtJBins) {
        jet_pt_counts[iPtJ] += thisjwgt;
        h2_jet_eta_phi[iFile][iPtJ]->Fill (jeta, jphi);
        nj++; // add to total number of "trigger" jets
      }
    }

    // skip events with no jets (otherwise the "trigger" didn't fire). This just saves time with mixing really.
    if (nj == 0)
      continue;

    for (short iPtJX = 0; iPtJX < nPtJBins; iPtJX++) {
      h_jet_pt[iCent]->SetBinContent (iPtJX+1, h_jet_pt[iCent]->GetBinContent (iPtJX+1) + (ewgt)*(jet_pt_counts[iPtJX]));
      for (short iPtJY = 0; iPtJY < nPtJBins; iPtJY++)
        h2_jet_pt_cov[iCent]->SetBinContent (iPtJX+1, iPtJY+1, h2_jet_pt_cov[iCent]->GetBinContent (iPtJX+1, iPtJY+1) + (ewgt)*(jet_pt_counts[iPtJX])*(jet_pt_counts[iPtJY]));
    }


    // the weights for the per-trigger jet yield are "w_evt"
    h_evt_counts[iFile]->SetBinContent (1, h_evt_counts[iFile]->GetBinContent (1) + 1);
    h_evt_counts[iFile]->SetBinContent (2, h_evt_counts[iFile]->GetBinContent (2) + ewgt);
    h_evt_counts[iFile]->SetBinContent (3, h_evt_counts[iFile]->GetBinContent (3) + ewgt*ewgt);

    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
      const float jwgt = jet_pt_counts[iPtJ];

      if (jwgt <= 0)
        continue; // IMPORTANT! Skip entries where there were no jets. These shouldn't contribute to the per-jet track yield at all.

      // the weights for the per-jet track yield are "w_evt * sum (w_jet)"
      h_jet_counts[iFile][iPtJ]->SetBinContent (1, h_jet_counts[iFile][iPtJ]->GetBinContent (1) + 1);
      h_jet_counts[iFile][iPtJ]->SetBinContent (2, h_jet_counts[iFile][iPtJ]->GetBinContent (2) + ewgt*jwgt);
      h_jet_counts[iFile][iPtJ]->SetBinContent (3, h_jet_counts[iFile][iPtJ]->GetBinContent (3) + std::pow (ewgt*jwgt, 2));

    } // end loop over iPtJ


    // System boost in p+Pb or pp -- requires a run number to be specified to check for pPb vs Pbp
    const float yboost = GetBoost (run_number);


    // do mixed event procedure: get a good event to mix with here
    if (doMixing) {

      // get the "centrality" bin for the current event
      const short iMixCent = GetBin (mixCentBins, nMixCentBins, fcalA_et + (Ispp () ? fcalC_et : 0));
      if (iMixCent < 0 || nMixCentBins <= iMixCent)
        continue; // sanity check (not too peripheral or too central)

      const long oldTrkEvt = iTrkEvt;
      bool goodEvent = false;
      do {
        iTrkEvt = (iTrkEvt + 1) % nTrkEvts; // make sure to wrap around

        // make sure we are not mixing an event with itself
        if (sameTreeMixing && iTrkEvt == iEvt % jetsTree->GetEntries ())
          continue;

        tracksTree->GetEntry (iTrkEvt);

        // triggering cut, require appropriate track selection trigger to have fired
        if (!trackTreeTrigger->trigDecision)
          continue;

        // vertexing cuts, require no pileup vertices and primary vertex with |vz| < 150mm
        float vz_matching = -999;
        {
          bool hasPrimary = false;
          bool hasPileup = false;
          for (short iVert = 0; iVert < nvert_matching; iVert++) {
            if (vert_type_matching[iVert] == 1) {
              hasPrimary = true;
              vz_matching = vert_z_matching[iVert];
            }
            else if (vert_type_matching[iVert] == 3)
              hasPileup = true;
          }
          if (hasPileup || std::fabs (vz_matching) > 150 || !hasPrimary)
            continue;
        }

        // further mixing categories -- Zdc centrality in p+Pb data and FCal centrality in MC
        if (IsCollisions () && IspPb () && GetBin (zdcCentBins, nZdcCentBins, ZdcCalibEnergy_A_matching*1e3) != GetBin (zdcCentBins, nZdcCentBins, ZdcCalibEnergy_A*1e3))
          continue; // require the same ZDC centrality if doing the p+Pb mixed event
        //else if (!IsCollisions () && IspPb () && GetBin (fcalCentBins, nFcalCentBins, fcalA_et_matching) != iCent) 
        //  continue; // require the same FCal centrality in MC
        if (GetBin (mixCentBins, nMixCentBins, fcalA_et_matching + (Ispp () ? fcalC_et_matching : 0)) != iMixCent)
          continue; // then match FCal centrality with jet ("trigger") event

        if (DoMixCatVar4 () || DoMixCatVar5 ()) {
          // calculate Q2 vector for matched event
          QnVector Q2mix = (IspPb () ? GetPbQ2Vec (true) : GetProtonQ2Vec (true));
          if (DeltaPhi (psi2, std::atan2 (Q2mix.second, Q2mix.first)/2.) > M_PI / 4)
            continue; // cut on psi2 matching for mixed event -- exploratory for now!
        }

        // if we've made it this far without proceeding to the next loop then its a good event!
        goodEvent = true;
      }
      while (!goodEvent && iTrkEvt != oldTrkEvt);

      if (!goodEvent) {
        std::cout << "Sorry, I could not find a good mixed event, so I will be skipping this very interesting trigger event!" << std::endl;
        continue;
      }
    } // now we have a good mixed event match (if applicable)


    // initialize all yields to 0
    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
      for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
        for (short iDPhi = 0; iDPhi < nDPhiBins; iDPhi++) {
          jet_trk_dphi_counts[iPtJ][iPtCh][iDPhi] = 0;
        }
      } // end loop over iPtCh
      for (short iDir = 0; iDir < nDir; iDir++) {
        for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
          jet_trk_pt_counts[iPtJ][iDir][iPtCh] = 0;
        }
      } // end loop over iDir
    } // end loop over iPtJ


    //// Get multiplicity bin for event, must be done after mixing (important)!
    //short iMult = 0;
    //while (iMult < nMultBins-1 && multBins[iMult+1] < trk_n)
    //  iMult++;
    //if (nMultBins < iMult)
    //  continue;
    const short iMult = nMultBins-1;


    // loop over all jets again in the event but now correlate tracks
    for (short iJet = 0; iJet < jn; iJet++) {

      const float jpt  = (UseTruthJets () ? GetAktTruthJetPt  (iJet, r0p4) : GetAktHIJetPt  (iJet, r0p4, nJESVar));
      const float jeta = (UseTruthJets () ? GetAktTruthJetEta (iJet, r0p4) : GetAktHIJetEta (iJet, r0p4, nJESVar));
      const float jphi = (UseTruthJets () ? GetAktTruthJetPhi (iJet, r0p4) : GetAktHIJetPhi (iJet, r0p4, nJESVar));

      if (!MeetsJetAcceptanceCuts (iJet, r0p4, nJESVar))
        continue; // jet eta/phi & timing cuts
      if (!MeetsJetPtCut (jpt))
        continue; // jet pT cuts, for trigger relevance purposes

      const float thisjwgt = GetAktJetWeight (jpt, jeta, jphi, r0p4);
      if (thisjwgt <= 0.)
        continue; // sanity check

      // In MC, check that the jet is truth-matched. Otherwise problems can emerge, particularly in the overlay.
      if (!IsCollisions () && !UseTruthJets ()) {
        const short iTJet = GetAktTruthJetMatch (iJet, r0p4);
        if (iTJet == -1)
          continue; // in case of no truth jet match
        const float tjpt = GetAktTruthJetPt (iTJet, r0p4);
        if (std::fabs (jpt/tjpt - 1) > 3*0.01*f_jer->Eval (tjpt))
          continue; // cut on jets reconstructed well outside the JER
      }

      const short iPtJ = GetPtJBin (jpt);
      if (iPtJ < 0 || nPtJBins <= iPtJ) 
        continue; // sanity check

      // correlate charged particles with this jet  
      for (short iTrk = 0; iTrk < trk_n; iTrk++) {

        if (!MeetsTrackCuts (iTrk, nTrkWPVar))
          continue; // cut on bad quality tracks

        if (trk_pt[iTrk] < pTChBins[0])
          continue; // histogram pT cut on tracks

        //// better rapidity selection by considering mass more carefully
        //const double trk_m = pion_mass; // assume everyone is a pi^+
        //const double trk_en = std::sqrt (std::pow (trk_m, 2) + trk_pt[iTrk] * std::cosh ((double)(trk_eta[iTrk])));
        //const double trk_pz = ((double)trk_pt[iTrk]) * std::sinh ((double)trk_eta[iTrk]);
        //const double trk_y = (trk_en > trk_pz ? 0.5 * std::log ((trk_en + trk_pz)/(trk_en - trk_pz)) : 0.);

        const float trk_y = trk_eta[iTrk];
        if (std::fabs (trk_y - yboost) > 2.035) // 2.035 = 2.5 - 0.465
          continue; // rapidity acceptance cut

        const float dphi = DeltaPhi (jphi, trk_phi[iTrk]);
        const short iDPhi = GetDPhiBin (dphi);
        if (iDPhi < 0 || nDPhiBins < iDPhi)
          continue; // sanity check to avoid under/over-flow

        short iEta = 0;
        while (iEta < nEtaTrkBins && etaTrkBins[iEta+1] < std::fabs (trk_eta[iTrk]))
          iEta++;
        if (iEta == nEtaTrkBins)
          continue;


        const float teff = h2_trk_eff[iMult]->GetBinContent (h2_trk_eff[iMult]->FindBin (trk_eta[iTrk], trk_pt[iTrk]));
        const float tpur = (DoPrimFitVar () ? h_trk_pur[iEta]->GetBinContent (h_trk_pur[iEta]->FindBin (trk_pt[iTrk])) : f_trk_pur[iEta]->Eval (trk_pt[iTrk]));
        const float twgt = thisjwgt * (UseTruthParticles () ? 1 : (teff > 0. ? tpur / teff : 0.));

        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
          if (pTChStrCuts[iPtCh].first < trk_pt[iTrk] && trk_pt[iTrk] < pTChStrCuts[iPtCh].second) {
            jet_trk_dphi_counts[iPtJ][iPtCh][iDPhi] += twgt;
            break;
          }
        }

        const short iPtCh = GetPtChBin (trk_pt[iTrk]);
        if (iPtCh < 0 || iPtCh >= nPtChBins)
          continue; // sanity check to avoid under/over-flow

        short iDir = -1;
        if (dphi < M_PI/8.)
          iDir = 0;
          //jet_trk_pt_ns_counts[iPtCh]           += twgt;
        else if (M_PI/3. < dphi && dphi < 2.*M_PI/3.)
          iDir = 1;
          //jet_trk_pt_perp_counts[iPtCh]         += twgt;
        else if (dphi > 7.*M_PI/8.)
          iDir = 2;
          //jet_trk_pt_as_counts[iPtCh]           += twgt;
        if (iDir != -1)
          jet_trk_pt_counts[iPtJ][iDir][iPtCh]  += twgt;

      } // end loop over tracks

    } // end loop over jets


    // evaluate the per-jet hadron yields for that event by dividing out the number of jets
    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
      const float jwgt = jet_pt_counts[iPtJ];

      if (jwgt == 0)
        continue; // IMPORTANT! Skip entries where there were no jets. These shouldn't contribute to the per-jet track yield at all.

      for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
        for (short iDPhi = 0; iDPhi < nDPhiBins; iDPhi++) {
          jet_trk_dphi_counts[iPtJ][iPtCh][iDPhi] = jet_trk_dphi_counts[iPtJ][iPtCh][iDPhi] / jwgt;
        }
      } // end loop over iPtCh
      for (short iDir = 0; iDir < nDir; iDir++) {
        for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
          jet_trk_pt_counts[iPtJ][iDir][iPtCh] = jet_trk_pt_counts[iPtJ][iDir][iPtCh] / jwgt;
        }
      } // end loop over iDir
    } // end loop over iPtJ


    // store results in averaging & covariance histograms
    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
      const float jwgt = jet_pt_counts[iPtJ];

      if (jwgt == 0)
        continue; // IMPORTANT! Skip entries where there were no jets. These shouldn't contribute to the per-jet track yield at all.

      for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
        TH1D* h = h_jet_trk_dphi[iFile][iPtJ][iPtCh];
        TH2D* h2 = h2_jet_trk_dphi_cov[iFile][iPtJ][iPtCh];
        double* counts = jet_trk_dphi_counts[iPtJ][iPtCh];
        for (short iDPhiX = 0; iDPhiX < nDPhiBins; iDPhiX++) {
          h->SetBinContent (iDPhiX+1, h->GetBinContent (iDPhiX+1) + (ewgt)*(jwgt)*(counts[iDPhiX]));
          for (short iDPhiY = 0; iDPhiY < nDPhiBins; iDPhiY++)
            h2->SetBinContent (iDPhiX+1, iDPhiY+1, h2->GetBinContent (iDPhiX+1, iDPhiY+1) + (ewgt)*(jwgt)*(counts[iDPhiX])*(counts[iDPhiY]));
        }
      } // end loop over iPtCh

      for (short iDir = 0; iDir < nDir; iDir++) {
        TH1D* h = h_jet_trk_pt[iFile][iPtJ][iDir];
        TH2D* h2 = h2_jet_trk_pt_cov[iFile][iPtJ][iDir];
        double* counts = jet_trk_pt_counts[iPtJ][iDir];
        for (short iPtChX = 0; iPtChX < nPtChBins; iPtChX++) {
          h->SetBinContent (iPtChX+1, h->GetBinContent (iPtChX+1) + (ewgt)*(jwgt)*(counts[iPtChX]));
          for (short iPtChY = 0; iPtChY < nPtChBins; iPtChY++)
            h2->SetBinContent (iPtChX+1, iPtChY+1, h2->GetBinContent (iPtChX+1, iPtChY+1) + (ewgt)*(jwgt)*(counts[iPtChX])*(counts[iPtChY]));
        }
      } // end loop over iDir

    } // end loop over iPtJ

  } // end loop over events
  std::cout << "Finished event loop." << std::endl;


  SaferDelete (&f_jer);

  SaferDelete (&h2_trk_eff);

  for (short iEta = 0; iEta < nEtaTrkBins; iEta++) {
    SaferDelete (&f_trk_pur[iEta]);
    SaferDelete (&h_trk_pur[iEta]);
  }
  Delete1DArray (f_trk_pur, nEtaTrkBins);
  Delete1DArray (h_trk_pur, nEtaTrkBins);

  //SaferDelete (&h_probs);

  SaferDelete (&m_overlay_fcalet);


  for (short iFile = 0; iFile < nFiles; iFile++) {

    outFiles[iFile]->cd ();

    h_evt_counts[iFile]->Write ();

    h_jet_pt[iFile]->Write ();
    h2_jet_pt_cov[iFile]->Write ();

    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

      h_jet_counts[iFile][iPtJ]->Write ();

      h2_jet_eta_phi[iFile][iPtJ]->Write ();

      for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
        h_jet_trk_dphi[iFile][iPtJ][iPtCh]->Write ();
        h2_jet_trk_dphi_cov[iFile][iPtJ][iPtCh]->Write ();
      }
      for (short iDir = 0; iDir < nDir; iDir++) {
        h_jet_trk_pt[iFile][iPtJ][iDir]->Write ();
        h2_jet_trk_pt_cov[iFile][iPtJ][iDir]->Write ();
      }
    }

    outFiles[iFile]->Close ();
  }

  SaferDelete (&jetTrigger);
  SaferDelete (&trackTreeTrigger);


  auto end = std::chrono::steady_clock::now (); // time at end of correlator function.
  auto nseconds = std::chrono::duration_cast <chrono::seconds> (end - start).count ();
  std::cout << "Correlator completed in " << nseconds << " seconds" << std::endl;

  return true;
}




bool RunCorrelator (const char* directory,
                    const char* tag,
                    const char* jetsInFileName,
                    const char* outFilePattern,
                    const char* mixdirectory = nullptr,
                    const char* tracksInFileName = nullptr) {

  if (!Is5TeV () || (!Ispp17 () && !IspPb16 ())) {
    std::cout << "In RunCorrelator.cxx: Invalid or unsupported collision system, exiting." << std::endl;
    return false;
  }

  if (jet_min_pt == -2)
    jet_min_pt = (float) std::atof (std::getenv ("JET_MIN_PT"));
  if (jet_max_pt == -2)
    jet_max_pt = (float) std::atof (std::getenv ("JET_MAX_PT"));
  doMixing = (tracksInFileName != nullptr && strcmp (tracksInFileName, "") != 0);

  // check mixing configuration makes sense
  if (!doMixing && (DoMixCatVar1 () || DoMixCatVar2 () || DoMixCatVar3 () || DoMixCatVar4 () || DoMixCatVar5 () || DoMixCatVar6 ())) {
    std::cout << "In RunCorrelator.cxx: Configured to run mixing variations but no mixing configuration specified? Please investigate. Exiting." << std::endl;
    return false;
  }

  // histograms for some mixing variations are symlinks to central values to avoid file copy issues, make sure not to overwrite them by exiting gracefully
  if ((DoMixCatVar2 () || DoMixCatVar4 () || DoMixCatVar6 ()) && Ispp ()) {
    std::cout << "In RunCorrelator.cxx: MixCatVar2 and MixCatVar4 do nothing in pp, exiting gracefully." << std::endl;
    return true;
  }
  else if ((DoMixCatVar1 () || DoMixCatVar3 () || DoMixCatVar5 ()) && IspPb ()) {
    std::cout << "In RunCorrelator.cxx: MixCatVar1, MixCatVar3, and MixCatVar5 do nothing in p+Pb, exiting gracefully." << std::endl;
    return true;
  }


  // keeps track of whether the jets and tracks are taken from the same tree. This is important for making sure we don't mix with the same event.
  sameTreeMixing = doMixing && strcmp (tracksInFileName, jetsInFileName) == 0;
  if (sameTreeMixing) {
    std::cout << "In RunCorrelator.cxx: Will be mixing events from the same tree -- will check that we are not mixing events." << std::endl;
  }

  // cannot run mixing at truth-level (not yet defined)
  if (doMixing && (UseTruthJets () || UseTruthParticles ())) {
    std::cout << "In RunCorrelator.cxx: Configured to run mixing but at MC truth level. This behavior is not defined, please investigate. Exiting gracefully." << std::endl;
    return true;
  }


  //////////////////////////////////////////////////////////////////////////////////////////
  // Get TChain object for input jets tree
  //////////////////////////////////////////////////////////////////////////////////////////
  TChain* jetsInTree = new TChain ("bush", "input jets trees");

  {
    // opens a bunch of TTrees as a TChain from all files in a directory matching the file identifier
    std::cout << "DataPath = " << dataPath << std::endl;
    auto dir = gSystem->OpenDirectory (dataPath + directory);
    while (const char* f = gSystem->GetDirEntry (dir)) {
      TString file = TString (f);

      if (!file.Contains (jetsInFileName))
        continue;
      if (IsCollisions ()) {
        if ((UseMinBiasTriggers () && !file.Contains ("MinBias")) || (UseJetTriggers () && !file.Contains ("Main")))
          continue;
      }

      std::cout << "Adding " << dataPath + directory + "/" + file + "/*.root" << " to input jets trees TChain" << std::endl;
      jetsInTree->Add (dataPath + directory + "/" + file + "/*.root");
      break;
    }
    if (jetsInTree->GetEntries () == 0) {
      std::cout << "Info: In RunCorrelator.cxx: Jet input chain has " << jetsInTree->GetEntries () << " entries, exiting gracefully." << std::endl;
      return true;
    }
    if (jetsInTree->GetListOfFiles ()->GetEntries () == 0) {
      std::cout << "Info: In RunCorrelator.cxx: Jet input chain has " << jetsInTree->GetListOfFiles () ->GetEntries () << " files, exiting gracefully." << std::endl;
      return true;
    }
    std::cout << "Info: In RunCorrelator.cxx: Jet input chain has " << jetsInTree->GetListOfFiles ()->GetEntries () << " files, " << jetsInTree->GetEntries () << " entries" << std::endl;
  }



  //////////////////////////////////////////////////////////////////////////////////////////
  // Get TChain object for input tracks tree
  //////////////////////////////////////////////////////////////////////////////////////////
  TChain* tracksInTree = (doMixing ? new TChain ("bush", "input tracks trees") : nullptr);

  if (tracksInTree) {
    // opens a bunch of TTrees as a TChain from all files in a directory matching the file identifier
    std::cout << "DataPath = " << dataPath << std::endl;
    auto dir = gSystem->OpenDirectory (dataPath + mixdirectory);

    TString patternString = tracksInFileName;
    patternString.ReplaceAll (".", "\\."); // save '.' as normal character
    patternString.ReplaceAll ("*",".*"); // wildcard '.' of any length '*'
    std::cout << patternString << std::endl;

    TPRegexp pattern = TPRegexp (patternString);

    while (const char* f = gSystem->GetDirEntry (dir)) {
      TString file = TString (f);

      if (IsCollisions ()) {
        if (!file.Contains ("MinBias"))
          continue;
      }
      if (!pattern.MatchB (file))
        continue;
      //if (!file.Contains (tracksInFileName))
      //  continue;

      std::cout << "Adding " << dataPath + mixdirectory + "/" + file + "/*.root" << " to input tracks trees TChain" << std::endl;
      tracksInTree->Add (dataPath + mixdirectory + "/" + file + "/*.root");
    }
    if (tracksInTree->GetEntries () == 0) {
      std::cout << "Info: In RunCorrelator.cxx: Tracks input chain has " << tracksInTree->GetEntries () << " entries, exiting gracefully." << std::endl;
      return true;
    }
    if (tracksInTree->GetListOfFiles ()->GetEntries () == 0) {
      std::cout << "Info: In RunCorrelator.cxx: Tracks input chain has " << tracksInTree->GetListOfFiles () ->GetEntries () << " files, exiting gracefully." << std::endl;
      return true;
    }
    std::cout << "Info: In RunCorrelator.cxx: Tracks input chain has " << tracksInTree->GetListOfFiles ()->GetEntries () << " files, " << tracksInTree->GetEntries () << " entries" << std::endl;
  }
  else
    tracksInTree = jetsInTree;


  truth_jet_min_pt = GetJZXR04MinPt (TString (jetsInFileName));
  truth_jet_max_pt = GetJZXR04MaxPt (TString (jetsInFileName));
  if (truth_jet_min_pt != 0)
    std::cout << "Checking for leading truth jet with pT > " << truth_jet_min_pt << std::endl;
  if (truth_jet_max_pt != FLT_MAX)
    std::cout << "Checking for leading truth jet with pT < " << truth_jet_max_pt << std::endl;

  return Correlator (TString (tag), outFilePattern, jetsInTree, tracksInTree);

}



int main (int argc, char** argv) {

  int argn                = 1;
  const char* subdir      = argv[argn++];
  const char* tag         = argv[argn++];
  //jet_min_pt = (double) std::atof (argv[argn++]);
  //jet_max_pt = (double) std::atof (argv[argn++]);

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

  TString tType = TString (argv[argn++]);
  if (!SetTriggerType (tType)) {
    std::cout << "In RunCorrelator.cxx: Invalid trigger type, exiting." << std::endl;
    return 3;
  }

  //jet_min_pt = (!IsCollisions () || !UseJetTriggers ()) ? pTJBins[0] : (UseJ50Triggers () ? 60 : 130); // if we want to use j100 for pp
  jet_min_pt = (!IsCollisions () || !UseJetTriggers ()) ? pTJBins[0] : 60;
  //jet_max_pt = (!IsCollisions () || (!Ispp () && UseJ50Triggers ()) || (Ispp () && UseJ100Triggers ())) ? pTJBins[nPtJBins] : (UseMinBiasTriggers () ? 60 : 130); // if we want to use j100 for pp
  jet_max_pt = (!IsCollisions () || UseJ50Triggers ()) ? pTJBins[nPtJBins] : 60;

  TString sFlag = TString (argv[argn++]);
  if (!SetSystFlag (sFlag)) {
    std::cout << "In RunCorrelator.cxx: Invalid systematic flag, exiting." << std::endl;
    return 4;
  }

  const char* outFilePattern        = (argc > argn && argv[argn] ? argv[argn++] : "");

  const char* inFileName            = (argc > argn && argv[argn] ? argv[argn++] : "");
  if (!IsCollisions ()) {
    GetMCWeights (inFileName);
  }

  const char* mixdir                = (argc > argn && argv[argn] ? argv[argn++] : "");
  const char* assocInFileName       = (argc > argn && argv[argn] ? argv[argn++] : "");

  std::cout << "Configuration set to";
  std::cout << "\n\tsubdir = " << subdir;
  std::cout << "\n\tCorrelatorTag = " << tag;
  std::cout << "\n\tjet_min_pt = " << jet_min_pt;
  std::cout << "\n\tjet_max_pt = " << jet_max_pt;
  std::cout << "\n\tCollisionSystem = " << ToTString (collisionSystem);
  std::cout << "\n\tDataType = " << ToTString (dataType);
  std::cout << "\n\tTriggerType = " << ToTString (triggerType);
  std::cout << "\n\tSystFlag = " << ToTString (systFlag);
  std::cout << "\n\toutFilePattern = " << outFilePattern;
  std::cout << "\n\tinFileName = " << inFileName;
  if (crossSectionPicoBarns != 0.)
    std::cout << "\n\t  --> deduced crossSectionPicoBarns = " << crossSectionPicoBarns;
  if (mcFilterEfficiency != 0.)
    std::cout << "\n\t  --> deduced mcFilterEfficiency = " << mcFilterEfficiency;
  if (mcNumberEvents != 0.)
    std::cout << "\n\t  --> deduced mcNumberEvents = " << mcNumberEvents;
  std::cout << "\n\tmixingDir = " << mixdir;
  std::cout << "\n\tassocInFileName = " << assocInFileName;
  std::cout << std::endl;


  std::cout << "Running Correlator algorithm..." << std::endl;
  bool success = RunCorrelator (subdir, tag, inFileName, outFilePattern, mixdir, assocInFileName);

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
