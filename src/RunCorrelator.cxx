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


using namespace JetHadronCorrelations;


bool doMixing = false;
bool sameTreeMixing = false;

float truth_jet_min_pt = 0;
float truth_jet_max_pt = FLT_MAX;



// Main function for performing jet-hadron correlations. (Wrapper function available below.)
bool Correlator (const char* tag, const char* outFilePattern, TTree* jetsTree, TTree* tracksTree = nullptr) {

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
    if (UseJet50GeVTriggers ()) {
      jetTrigger = new Trigger (jet_trig_name[0]);
      std::cout << "Info: In RunCorrelator.cxx: Looking for " << jet_trig_name[0] << " trigger" << std::endl;
      jetsTree->SetBranchAddress ((jet_trig_name[0]+"_decision").c_str (), &(jetTrigger->trigDecision));
      jetsTree->SetBranchAddress ((jet_trig_name[0]+"_prescale").c_str (), &(jetTrigger->trigPrescale));
    }
    else if (UseJet100GeVTriggers ()) {
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

    if (doMixing) {
      trackTreeTrigger = new Trigger (minbias_trig_name[0]);
      tracksTree->SetBranchAddress ((minbias_trig_name[0]+"_decision").c_str (), &(trackTreeTrigger->trigDecision));
      tracksTree->SetBranchAddress ((minbias_trig_name[0]+"_prescale").c_str (), &(trackTreeTrigger->trigPrescale));
    }
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
  }
  //tracksTree->SetBranchAddress ("trk_d0",                 &trk_d0);
  //tracksTree->SetBranchAddress ("trk_d0sig",              &trk_d0sig);
  //tracksTree->SetBranchAddress ("trk_z0",                 &trk_z0);
  //tracksTree->SetBranchAddress ("trk_z0sig",              &trk_z0sig);
  //tracksTree->SetBranchAddress ("trk_theta",              &trk_theta);
  //tracksTree->SetBranchAddress ("trk_vz",                 &trk_vz);
  //tracksTree->SetBranchAddress ("trk_nBLayerHits",        &trk_nBLayerHits);
  //tracksTree->SetBranchAddress ("trk_nBLayerSharedHits",  &trk_nBLayerSharedHits);
  //tracksTree->SetBranchAddress ("trk_nPixelHits",         &trk_nPixelHits);
  //tracksTree->SetBranchAddress ("trk_nPixelDeadSensors",  &trk_nPixelDeadSensors);
  //tracksTree->SetBranchAddress ("trk_nPixelSharedHits",   &trk_nPixelSharedHits);
  //tracksTree->SetBranchAddress ("trk_nSCTHits",           &trk_nSCTHits);
  //tracksTree->SetBranchAddress ("trk_nSCTDeadSensors",    &trk_nSCTDeadSensors);
  //tracksTree->SetBranchAddress ("trk_nSCTSharedHits",     &trk_nSCTSharedHits);
  //tracksTree->SetBranchAddress ("trk_nTRTHits",           &trk_nTRTHits);
  //tracksTree->SetBranchAddress ("trk_nTRTSharedHits",     &trk_nTRTSharedHits);
  //if (!IsCollisions ()) {
  //  tracksTree->SetBranchAddress ("trk_prob_truth",     &trk_prob_truth);
  //  tracksTree->SetBranchAddress ("trk_truth_pt",       &trk_truth_pt);
  //  tracksTree->SetBranchAddress ("trk_truth_eta",      &trk_truth_eta);
  //  tracksTree->SetBranchAddress ("trk_truth_phi",      &trk_truth_phi);
  //  tracksTree->SetBranchAddress ("trk_truth_charge",   &trk_truth_charge);
  //  tracksTree->SetBranchAddress ("trk_truth_type",     &trk_truth_type);
  //  tracksTree->SetBranchAddress ("trk_truth_orig",     &trk_truth_orig);
  //  tracksTree->SetBranchAddress ("trk_truth_barcode",  &trk_truth_barcode);
  //  tracksTree->SetBranchAddress ("trk_truth_pdgid",    &trk_truth_pdgid);
  //  tracksTree->SetBranchAddress ("trk_truth_vz",       &trk_truth_vz);
  //  tracksTree->SetBranchAddress ("trk_truth_nIn",      &trk_truth_nIn);
  //  tracksTree->SetBranchAddress ("trk_truth_isHadron", &trk_truth_isHadron);
  //}


  if (!IsCollisions ()) {
    jetsTree->SetBranchAddress ("akt4_truth_jet_n",     &akt4_truth_jet_n);
    jetsTree->SetBranchAddress ("akt4_truth_jet_pt",    &akt4_truth_jet_pt);
    jetsTree->SetBranchAddress ("akt4_truth_jet_eta",   &akt4_truth_jet_eta);
    jetsTree->SetBranchAddress ("akt4_truth_jet_phi",   &akt4_truth_jet_phi);
    jetsTree->SetBranchAddress ("akt4_truth_jet_e",     &akt4_truth_jet_e);
  }

  jetsTree->SetBranchAddress ("akt4_hi_jet_n",            &akt4_hi_jet_n);
  //jetsTree->SetBranchAddress ("akt4_hi_jet_pt_precalib",  &akt4_hi_jet_pt_precalib);
  jetsTree->SetBranchAddress ("akt4_hi_jet_pt_etajes",    &akt4_hi_jet_pt_etajes);
  jetsTree->SetBranchAddress ("akt4_hi_jet_pt_xcalib",    &akt4_hi_jet_pt_xcalib);
  //jetsTree->SetBranchAddress ("akt4_hi_jet_eta_precalib", &akt4_hi_jet_eta_precalib);
  jetsTree->SetBranchAddress ("akt4_hi_jet_eta_etajes",   &akt4_hi_jet_eta_etajes);
  jetsTree->SetBranchAddress ("akt4_hi_jet_eta_xcalib",   &akt4_hi_jet_eta_xcalib);
  jetsTree->SetBranchAddress ("akt4_hi_jet_phi",          &akt4_hi_jet_phi);
  //jetsTree->SetBranchAddress ("akt4_hi_jet_e_precalib",   &akt4_hi_jet_e_precalib);
  jetsTree->SetBranchAddress ("akt4_hi_jet_e_etajes",     &akt4_hi_jet_e_etajes);
  jetsTree->SetBranchAddress ("akt4_hi_jet_e_xcalib",     &akt4_hi_jet_e_xcalib);
  //jetsTree->SetBranchAddress ("akt4_hi_jet_sub_et",       &akt4_hi_jet_sub_et);
  //jetsTree->SetBranchAddress ("akt4_hi_jet_sub_e",        &akt4_hi_jet_sub_e);
  jetsTree->SetBranchAddress ("akt4_hi_jet_timing",       &akt4_hi_jet_timing);


  const int nJESVar = GetNJESVar ();
  if (!IsCollisions () && nJESVar != -1) {
    std::cout << "Info: In RunCorrelator.cxx: Branching JES variation " << nJESVar << std::endl;
    jetsTree->SetBranchAddress (Form ("akt4_hi_jet_pt_sys_JES_%i", nJESVar), akt4_hi_jet_pt_sys_JES_ALL[nJESVar]);
  }


  const int nTrkWPVar = (DoHITightVar () ? 1 : (DoHILooseVar () ? 2 : 0));


  // setup centrality bins (only relevant for p+Pb)
  double* centBins = (DoFcalCentVar () ? fcalCentBins : (DoFineFcalCentVar () ? fineFcalCentBins : (!IsCollisions () ? fcalCentBins : zdcCentBins)));
  const int nCentBins = (DoFcalCentVar () ? nFcalCentBins : (DoFineFcalCentVar () ? nFineFcalCentBins : (!IsCollisions () ? nFcalCentBins : nZdcCentBins)));


  TH2D* h2_trk_eff = LoadTrackingEfficiency ();
  TH1D** h_trk_pur = LoadTrackingPurity ();

  TF1** f_trk_pur = LoadTrackingPurityFuncs ();


  // create output histograms in output files
  const int nFiles = (Ispp () ? 1 : nCentBins);
  TFile** outFiles = new TFile*[nFiles];

  TH1D** h_evt_counts = new TH1D*[nFiles];
  TH1D** h_jet_counts = new TH1D*[nFiles];

  TH1D** h_jet_pt = new TH1D*[nFiles];
  TH2D** h2_jet_pt_cov = new TH2D*[nFiles];

  TH2D** h2_jet_eta_phi = new TH2D*[nFiles];

  TH1D** h_jet_trk_dphi_gt0p5_lt1 = new TH1D*[nFiles];
  TH2D** h2_jet_trk_dphi_gt0p5_lt1_cov = new TH2D*[nFiles];
  TH1D** h_jet_trk_dphi_gt1_lt1p5 = new TH1D*[nFiles];
  TH2D** h2_jet_trk_dphi_gt1_lt1p5_cov = new TH2D*[nFiles];
  TH1D** h_jet_trk_dphi_gt1p5_lt2 = new TH1D*[nFiles];
  TH2D** h2_jet_trk_dphi_gt1p5_lt2_cov = new TH2D*[nFiles];
  TH1D** h_jet_trk_dphi_gt2_lt4 = new TH1D*[nFiles];
  TH2D** h2_jet_trk_dphi_gt2_lt4_cov = new TH2D*[nFiles];
  TH1D** h_jet_trk_dphi_gt4_lt6 = new TH1D*[nFiles];
  TH2D** h2_jet_trk_dphi_gt4_lt6_cov = new TH2D*[nFiles];
  TH1D** h_jet_trk_dphi_gt6_lt8 = new TH1D*[nFiles];
  TH2D** h2_jet_trk_dphi_gt6_lt8_cov = new TH2D*[nFiles];
  TH1D** h_jet_trk_dphi_gt8_lt10 = new TH1D*[nFiles];
  TH2D** h2_jet_trk_dphi_gt8_lt10_cov = new TH2D*[nFiles];
  TH1D** h_jet_trk_dphi_gt10_lt15 = new TH1D*[nFiles];
  TH2D** h2_jet_trk_dphi_gt10_lt15_cov = new TH2D*[nFiles];
  TH1D** h_jet_trk_dphi_gt15_lt20 = new TH1D*[nFiles];
  TH2D** h2_jet_trk_dphi_gt15_lt20_cov = new TH2D*[nFiles];
  TH1D** h_jet_trk_dphi_gt20_lt30 = new TH1D*[nFiles];
  TH2D** h2_jet_trk_dphi_gt20_lt30_cov = new TH2D*[nFiles];

  TH1D** h_jet_trk_pt_ns = new TH1D*[nFiles];
  TH2D** h2_jet_trk_pt_ns_cov = new TH2D*[nFiles];
  TH1D** h_jet_trk_pt_perp = new TH1D*[nFiles];
  TH2D** h2_jet_trk_pt_perp_cov = new TH2D*[nFiles];
  TH1D** h_jet_trk_pt_as = new TH1D*[nFiles];
  TH2D** h2_jet_trk_pt_as_cov = new TH2D*[nFiles];


  for (int iFile = 0; iFile < nFiles; iFile++) {

    TString outFileName = outFilePattern;
    outFileName.ReplaceAll ("iCent_", Form ("iCent%i_", iFile));
    TFile* outFile = new TFile (outFileName.Data (), "recreate");
    outFiles[iFile] = outFile;

    h_evt_counts[iFile] = new TH1D (Form ("h_evt_counts_%s", tag), "", 3, -0.5, 2.5);
    h_jet_counts[iFile] = new TH1D (Form ("h_jet_counts_%s", tag), "", 3, -0.5, 2.5);

    h_jet_pt[iFile] = new TH1D (Form ("h_jet_pt_%s", tag), ";#it{p}_{T}^{jet} [GeV];(1/N_{jet}) (dN_{jet}/d#it{p}_{T}) [GeV^{-1}]", nPtJBins, pTJBins);
    h2_jet_pt_cov[iFile] = new TH2D (Form ("h2_jet_pt_cov_%s", tag), ";#it{p}_{T}^{jet} [GeV];#it{p}_{T}^{jet} [GeV];Covariance", nPtJBins, pTJBins, nPtJBins, pTJBins);

    h2_jet_eta_phi[iFile] = new TH2D (Form ("h2_jet_eta_phi_%s", tag), ";#eta^{jet};#phi^{jet};Counts", 64, -3.2, 3.2, 64, -M_PI, M_PI);

    h_jet_trk_dphi_gt0p5_lt1[iFile] = new TH1D (Form ("h_jet_trk_dphi_gt0p5_lt1_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
    h2_jet_trk_dphi_gt0p5_lt1_cov[iFile] = new TH2D (Form ("h2_jet_trk_dphi_gt0p5_lt1_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
    h_jet_trk_dphi_gt1_lt1p5[iFile] = new TH1D (Form ("h_jet_trk_dphi_gt1_lt1p5_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
    h2_jet_trk_dphi_gt1_lt1p5_cov[iFile] = new TH2D (Form ("h2_jet_trk_dphi_gt1_lt1p5_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
    h_jet_trk_dphi_gt1p5_lt2[iFile] = new TH1D (Form ("h_jet_trk_dphi_gt1p5_lt2_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
    h2_jet_trk_dphi_gt1p5_lt2_cov[iFile] = new TH2D (Form ("h2_jet_trk_dphi_gt1p5_lt2_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
    h_jet_trk_dphi_gt2_lt4[iFile] = new TH1D (Form ("h_jet_trk_dphi_gt2_lt4_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
    h2_jet_trk_dphi_gt2_lt4_cov[iFile] = new TH2D (Form ("h2_jet_trk_dphi_gt2_lt4_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
    h_jet_trk_dphi_gt4_lt6[iFile] = new TH1D (Form ("h_jet_trk_dphi_gt4_lt6_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
    h2_jet_trk_dphi_gt4_lt6_cov[iFile] = new TH2D (Form ("h2_jet_trk_dphi_gt4_lt6_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
    h_jet_trk_dphi_gt6_lt8[iFile] = new TH1D (Form ("h_jet_trk_dphi_gt6_lt8_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
    h2_jet_trk_dphi_gt6_lt8_cov[iFile] = new TH2D (Form ("h2_jet_trk_dphi_gt6_lt8_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
    h_jet_trk_dphi_gt8_lt10[iFile] = new TH1D (Form ("h_jet_trk_dphi_gt8_lt10_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
    h2_jet_trk_dphi_gt8_lt10_cov[iFile] = new TH2D (Form ("h2_jet_trk_dphi_gt8_lt10_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
    h_jet_trk_dphi_gt10_lt15[iFile] = new TH1D (Form ("h_jet_trk_dphi_gt10_lt15_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
    h2_jet_trk_dphi_gt10_lt15_cov[iFile] = new TH2D (Form ("h2_jet_trk_dphi_gt10_lt15_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
    h_jet_trk_dphi_gt15_lt20[iFile] = new TH1D (Form ("h_jet_trk_dphi_gt15_lt20_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
    h2_jet_trk_dphi_gt15_lt20_cov[iFile] = new TH2D (Form ("h2_jet_trk_dphi_gt15_lt20_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);
    h_jet_trk_dphi_gt20_lt30[iFile] = new TH1D (Form ("h_jet_trk_dphi_gt20_lt30_%s", tag), ";#Delta#phi;(1/N_{jet}) (dN_{ch}/d#Delta#phi)", nDPhiBins, dPhiBins);
    h2_jet_trk_dphi_gt20_lt30_cov[iFile] = new TH2D (Form ("h2_jet_trk_dphi_gt20_lt30_cov_%s", tag), ";#Delta#phi;#Delta#phi;Covariance", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);

    h_jet_trk_pt_ns[iFile] = new TH1D (Form ("h_jet_trk_pt_ns_%s", tag), ";#it{p}_{T}^{ch} [GeV];(1/N_{jet}) (dN_{ch}/d#it{p}_{T}) [GeV^{-1}]", nPtChBins, pTChBins);
    h2_jet_trk_pt_ns_cov[iFile] = new TH2D (Form ("h2_jet_trk_pt_ns_cov_%s", tag), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{ch} [GeV];Covariance", nPtChBins, pTChBins, nPtChBins, pTChBins);
    h_jet_trk_pt_perp[iFile] = new TH1D (Form ("h_jet_trk_pt_perp_%s", tag), ";#it{p}_{T}^{ch} [GeV];(1/N_{jet}) (dN_{ch}/d#it{p}_{T}) [GeV^{-1}]", nPtChBins, pTChBins);
    h2_jet_trk_pt_perp_cov[iFile] = new TH2D (Form ("h2_jet_trk_pt_perp_cov_%s", tag), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{ch} [GeV];Covariance", nPtChBins, pTChBins, nPtChBins, pTChBins);
    h_jet_trk_pt_as[iFile] = new TH1D (Form ("h_jet_trk_pt_as_%s", tag), ";#it{p}_{T}^{ch} [GeV];(1/N_{jet}) (dN_{ch}/d#it{p}_{T}) [GeV^{-1}]", nPtChBins, pTChBins);
    h2_jet_trk_pt_as_cov[iFile] = new TH2D (Form ("h2_jet_trk_pt_as_cov_%s", tag), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{ch} [GeV];Covariance", nPtChBins, pTChBins, nPtChBins, pTChBins);

  }


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

  const long nEvts = (doMixing ? 20. : 1.) * (jetsTree->GetEntries ());
  const long nTrkEvts = tracksTree->GetEntries ();


  // setup centrality matching bins for event mixing according to job configuration
  double* mixCentBins = nullptr;
  int nMixCentBins = 0;

  if (Ispp ()) {
    if      (DoMixCatVar1 ()) { mixCentBins = ppMixVar1Bins;    nMixCentBins = nppMixVar1Bins;    }
    else if (DoMixCatVar3 ()) { mixCentBins = ppMixVar3Bins;    nMixCentBins = nppMixVar3Bins;    }
    else                      { mixCentBins = ppMixBins;        nMixCentBins = nppMixBins;        }
  }
  else {
    if (DoMixCatVar2 ())      { mixCentBins = fcalMixVar2Bins;  nMixCentBins = nFcalMixVar2Bins;  }
    else                      { mixCentBins = fcalMixBins;      nMixCentBins = nFcalMixBins;      }
  }


  // for smearing truth jets if specified
  TRandom3* rndm = new TRandom3 ();
  TF1* f_jer = LoadJetEnergyResFunction ();


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


    // triggering cut, require appropriate jet trigger to have fired
    if (IsCollisions ()) {
      if (!jetTrigger->trigDecision)
        continue;
      else
        event_weight = 1; // optionally set to trigger prescale, but not really necessary
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


    // MC only -- filter sample based on min/max of pThat range
    // also set appropriate JZ weight
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


    int iCent = 0;
    // p+Pb only -- split events by centrality
    if (IspPb ()) {
      const double centVar = ((!IsCollisions () || DoFcalCentVar ()  || DoFineFcalCentVar ()) ? fcalA_et : (ZdcCalibEnergy_A * 1e3));
      iCent = GetBin (centBins, nCentBins, centVar);
      if (iCent < 0 || iCent > nCentBins-1)
        continue;
    }
    const int iFile = iCent;


    // calculate event plane angle based on FCal q-vector moments
    QnVector Q2 = (IspPb () ? GetPbQ2Vec () : GetProtonQ2Vec ());
    const float psi2 = std::atan2 (Q2.second, Q2.first)/2.;


    // event weights -- these are not 1 in MC due to JZ weighting (they are 1 in data)
    const double ewgt = event_weight;
    if (ewgt <= 0.)
      continue;


    // initialize all histogramming bins to 0 for this event
    for (int iX = 0; iX < nPtJBins; iX++)
      jet_pt_counts[iX] = 0;
    for (int iX = 0; iX < nDPhiBins; iX++) {
      jet_trk_dphi_gt0p5_lt1_counts[iX] = 0;
      jet_trk_dphi_gt1_lt1p5_counts[iX] = 0;
      jet_trk_dphi_gt1p5_lt2_counts[iX] = 0;
      jet_trk_dphi_gt2_lt4_counts[iX]   = 0;
      jet_trk_dphi_gt4_lt6_counts[iX]   = 0;
      jet_trk_dphi_gt6_lt8_counts[iX]   = 0;
      jet_trk_dphi_gt8_lt10_counts[iX]  = 0;
      jet_trk_dphi_gt10_lt15_counts[iX] = 0;
      jet_trk_dphi_gt15_lt20_counts[iX] = 0;
      jet_trk_dphi_gt20_lt30_counts[iX] = 0;
    }
    for (int iX = 0; iX < nPtChBins; iX++) {
      jet_trk_pt_ns_counts[iX]          = 0;
      jet_trk_pt_perp_counts[iX]        = 0;
      jet_trk_pt_as_counts[iX]          = 0;
    }

   
    // possibly smear truth jet pTs by JER
    if (DoMCTruthJESSmear ()) {
      const int jn = GetAktTruthJetN (r0p4);
      for (int iJet = 0; iJet < jn; iJet++) {
        const float sf = rndm->Gaus (1, 0.01*f_jer->Eval (akt4_truth_jet_pt[iJet]));
        akt4_truth_jet_pt[iJet] = sf * akt4_truth_jet_pt[iJet];
        akt4_truth_jet_e[iJet]  = sf * akt4_truth_jet_e[iJet];
      }
    }


    float nj = 0;
    float jwgt = 0;

    const int jn = UseTruthJets () ? GetAktTruthJetN (r0p4) :  GetAktHIJetN (r0p4);
    for (int iJet = 0; iJet < jn; iJet++) {

      const float jpt  = (UseTruthJets () ? GetAktTruthJetPt  (iJet, r0p4) : GetAktHIJetPt  (iJet, r0p4, nJESVar));
      const float jeta = (UseTruthJets () ? GetAktTruthJetEta (iJet, r0p4) : GetAktHIJetEta (iJet, r0p4, nJESVar));
      const float jphi = (UseTruthJets () ? GetAktTruthJetPhi (iJet, r0p4) : GetAktHIJetPhi (iJet, r0p4, nJESVar));

      if (!MeetsJetAcceptanceCuts (iJet, r0p4, nJESVar))
        continue; // jet eta/phi & timing cuts

      const float thisjwgt = GetAktJetWeight (jpt, jeta, jphi, r0p4);
      if (thisjwgt <= 0.)
        continue; // sanity check

      const short iPtJ = GetPtJBin (jpt);
      if (0 <= iPtJ && iPtJ < nPtJBins)
        jet_pt_counts[iPtJ] += jwgt;

      if (!MeetsJetPtCut (jpt))
        continue; // jet pT cuts

      h2_jet_eta_phi[iFile]->Fill (jeta, jphi);

      nj += 1;
      jwgt += thisjwgt;
    }

    for (int iX = 0; iX < nPtJBins; iX++) {
      h_jet_pt[iCent]->SetBinContent (iX+1, h_jet_pt[iCent]->GetBinContent (iX+1) + (ewgt)*(jet_pt_counts[iX]));
      for (int iY = 0; iY < nPtJBins; iY++)
        h2_jet_pt_cov[iCent]->SetBinContent (iX+1, iY+1, h2_jet_pt_cov[iCent]->GetBinContent (iX+1, iY+1) + (ewgt)*(jet_pt_counts[iX])*(jet_pt_counts[iY]));
    }


    // skip events with no jets (otherwise will have divide-by-zero errors later on)
    if (nj == 0)
      continue;

    const float yboost = GetBoost (run_number);

    h_evt_counts[iFile]->Fill (0);
    h_evt_counts[iFile]->Fill (1, ewgt);
    h_evt_counts[iFile]->Fill (2, ewgt*ewgt);

    h_jet_counts[iFile]->Fill (0); // adds to number of jets (i.e. denominator)
    h_jet_counts[iFile]->Fill (1, ewgt*jwgt);
    h_jet_counts[iFile]->Fill (2, std::pow (ewgt*jwgt, 2));


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
        if (IsCollisions () && !trackTreeTrigger->trigDecision)
          continue;

        // vertexing cuts, require no pileup vertices and primary vertex with |vz| < 150mm
        {
          bool hasPrimary = false;
          bool hasPileup = false;
          float vz = -999;
          for (int iVert = 0; iVert < nvert_matching; iVert++) {
            if (vert_type_matching[iVert] == 1) {
              hasPrimary = true;
              vz = vert_z_matching[iVert];
            }
            else if (vert_type_matching[iVert] == 3)
              hasPileup = true;
          }
          if (hasPileup || std::fabs (vz) > 150 || !hasPrimary)
            continue;
        }

        // further mixing categories -- Zdc centrality in p+Pb data and FCal centrality
        if (IsCollisions () && IspPb () && GetBin (zdcCentBins, nZdcCentBins, ZdcCalibEnergy_A_matching*1e3) != GetBin (zdcCentBins, nZdcCentBins, ZdcCalibEnergy_A*1e3))
          continue; // require the same ZDC centrality if doing the p+Pb mixed event
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
    for (int iX = 0; iX < nDPhiBins; iX++) {
      jet_trk_dphi_gt0p5_lt1_counts[iX] = 0;
      jet_trk_dphi_gt1_lt1p5_counts[iX] = 0;
      jet_trk_dphi_gt1p5_lt2_counts[iX] = 0;
      jet_trk_dphi_gt2_lt4_counts[iX]   = 0;
      jet_trk_dphi_gt4_lt6_counts[iX]   = 0;
      jet_trk_dphi_gt6_lt8_counts[iX]   = 0;
      jet_trk_dphi_gt8_lt10_counts[iX]  = 0;
      jet_trk_dphi_gt10_lt15_counts[iX] = 0;
      jet_trk_dphi_gt15_lt20_counts[iX] = 0;
      jet_trk_dphi_gt20_lt30_counts[iX] = 0;
    }
    for (int iX = 0; iX < nPtChBins; iX++) {
      jet_trk_pt_ns_counts[iX]          = 0;
      jet_trk_pt_perp_counts[iX]        = 0;
      jet_trk_pt_as_counts[iX]          = 0;
    }


    // loop over all jets again in the event but now correlate tracks
    for (int iJet = 0; iJet < jn; iJet++) {

      const float jpt  = (UseTruthJets () ? GetAktTruthJetPt  (iJet, r0p4) : GetAktHIJetPt  (iJet, r0p4, nJESVar));
      const float jeta = (UseTruthJets () ? GetAktTruthJetEta (iJet, r0p4) : GetAktHIJetEta (iJet, r0p4, nJESVar));
      const float jphi = (UseTruthJets () ? GetAktTruthJetPhi (iJet, r0p4) : GetAktHIJetPhi (iJet, r0p4, nJESVar));

      if (!MeetsJetAcceptanceCuts (iJet, r0p4, nJESVar))
        continue; // jet eta/phi & timing cuts
      if (!MeetsJetPtCut (jpt))
        continue; // jet pT cuts

      const float thisjwgt = GetAktJetWeight (jpt, jeta, jphi, r0p4);
      if (thisjwgt <= 0.)
        continue; // sanity check

      // correlate charged particles with this jet  
      for (int iTrk = 0; iTrk < trk_n; iTrk++) {

        if (!MeetsTrackCuts (iTrk, nTrkWPVar))
          continue; // cut on bad quality tracks

        if (trk_pt[iTrk] < pTChBins[0])
          continue; // histogram pT cut on tracks

        const short iPtCh = GetPtChBin (trk_pt[iTrk]);
        if (iPtCh < 0 || iPtCh >= nPtChBins)
          continue; // sanity check to avoid under/over-flow

        //// TODO -- better rapidity selection by considering mass more carefully
        //const double trk_m = pion_mass; // assume everyone is a pi^+
        //const double trk_en = std::sqrt (std::pow (trk_m, 2) + trk_pt[iTrk] * std::cosh ((double)(trk_eta[iTrk])));
        //const double trk_pz = ((double)trk_pt[iTrk]) * std::sinh ((double)trk_eta[iTrk]);
        //const double trk_y = (trk_en > trk_pz ? 0.5 * std::log ((trk_en + trk_pz)/(trk_en - trk_pz)) : 0.);

        const float trk_y = trk_eta[iTrk];
        if (std::fabs (trk_y - yboost) > 2.5 - 0.465)
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


        const float teff = h2_trk_eff->GetBinContent (h2_trk_eff->FindBin (trk_eta[iTrk], trk_pt[iTrk]));
        const float tpur = (DoPrimFitVar () ? h_trk_pur[iEta]->GetBinContent (h_trk_pur[iEta]->FindBin (trk_pt[iTrk])) : f_trk_pur[iEta]->Eval (trk_pt[iTrk]));
        const float twgt = thisjwgt * (UseTruthParticles () ? 1 : (teff > 0. ? tpur / teff : 0.));

        if (0.5 < trk_pt[iTrk] && trk_pt[iTrk] < 1)
          jet_trk_dphi_gt0p5_lt1_counts[iDPhi]  += twgt;
        else if (1 < trk_pt[iTrk] && trk_pt[iTrk] < 1.5)
          jet_trk_dphi_gt1_lt1p5_counts[iDPhi]  += twgt;
        else if (1.5 < trk_pt[iTrk] && trk_pt[iTrk] < 2)
          jet_trk_dphi_gt1p5_lt2_counts[iDPhi]  += twgt;
        else if (2 < trk_pt[iTrk] && trk_pt[iTrk] < 4)
          jet_trk_dphi_gt2_lt4_counts[iDPhi]    += twgt;
        else if (4 < trk_pt[iTrk] && trk_pt[iTrk] < 6)
          jet_trk_dphi_gt4_lt6_counts[iDPhi]    += twgt;
        else if (6 < trk_pt[iTrk] && trk_pt[iTrk] < 8)
          jet_trk_dphi_gt6_lt8_counts[iDPhi]    += twgt;
        else if (8 < trk_pt[iTrk] && trk_pt[iTrk] < 10)
          jet_trk_dphi_gt8_lt10_counts[iDPhi]   += twgt;
        else if (10 < trk_pt[iTrk] && trk_pt[iTrk] < 15)
          jet_trk_dphi_gt10_lt15_counts[iDPhi]  += twgt;
        else if (15 < trk_pt[iTrk] && trk_pt[iTrk] < 20)
          jet_trk_dphi_gt15_lt20_counts[iDPhi]  += twgt;
        else if (20 < trk_pt[iTrk] && trk_pt[iTrk] < 30)
          jet_trk_dphi_gt20_lt30_counts[iDPhi]  += twgt;

        if (dphi < M_PI/8.)
          jet_trk_pt_ns_counts[iPtCh]           += twgt;
        else if (M_PI/3. < dphi && dphi < 2.*M_PI/3.)
          jet_trk_pt_perp_counts[iPtCh]         += twgt;
        else if (dphi > 7.*M_PI/8.)
          jet_trk_pt_as_counts[iPtCh]           += twgt;

      } // end loop over tracks
    } // end loop over jets


    // calculate per-jet hadron yields for that event by dividing out the number of jets
    for (int iX = 0; iX < nDPhiBins; iX++) {
      jet_trk_dphi_gt0p5_lt1_counts[iX] = jet_trk_dphi_gt0p5_lt1_counts[iX] / jwgt;
      jet_trk_dphi_gt1_lt1p5_counts[iX] = jet_trk_dphi_gt1_lt1p5_counts[iX] / jwgt;
      jet_trk_dphi_gt1p5_lt2_counts[iX] = jet_trk_dphi_gt1p5_lt2_counts[iX] / jwgt;
      jet_trk_dphi_gt2_lt4_counts[iX]   = jet_trk_dphi_gt2_lt4_counts[iX]   / jwgt;
      jet_trk_dphi_gt4_lt6_counts[iX]   = jet_trk_dphi_gt4_lt6_counts[iX]   / jwgt;
      jet_trk_dphi_gt6_lt8_counts[iX]   = jet_trk_dphi_gt6_lt8_counts[iX]   / jwgt;
      jet_trk_dphi_gt8_lt10_counts[iX]  = jet_trk_dphi_gt8_lt10_counts[iX]  / jwgt;
      jet_trk_dphi_gt10_lt15_counts[iX] = jet_trk_dphi_gt10_lt15_counts[iX] / jwgt;
      jet_trk_dphi_gt15_lt20_counts[iX] = jet_trk_dphi_gt15_lt20_counts[iX] / jwgt;
      jet_trk_dphi_gt20_lt30_counts[iX] = jet_trk_dphi_gt20_lt30_counts[iX] / jwgt;
    }
    for (int iX = 0; iX < nPtChBins; iX++) {
      jet_trk_pt_ns_counts[iX]          = jet_trk_pt_ns_counts[iX]          / jwgt;
      jet_trk_pt_perp_counts[iX]        = jet_trk_pt_perp_counts[iX]        / jwgt;
      jet_trk_pt_as_counts[iX]          = jet_trk_pt_as_counts[iX]          / jwgt;
    }


    // store results in averaging & covariance histograms
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt0p5_lt1[iFile]->SetBinContent (iX+1, h_jet_trk_dphi_gt0p5_lt1[iFile]->GetBinContent (iX+1) + (ewgt)*(jet_trk_dphi_gt0p5_lt1_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt0p5_lt1_cov[iFile]->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt0p5_lt1_cov[iFile]->GetBinContent (iX+1, iY+1) + (ewgt)*(jet_trk_dphi_gt0p5_lt1_counts[iX])*(jet_trk_dphi_gt0p5_lt1_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt1_lt1p5[iFile]->SetBinContent (iX+1, h_jet_trk_dphi_gt1_lt1p5[iFile]->GetBinContent (iX+1) + (ewgt)*(jet_trk_dphi_gt1_lt1p5_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt1_lt1p5_cov[iFile]->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt1_lt1p5_cov[iFile]->GetBinContent (iX+1, iY+1) + (ewgt)*(jet_trk_dphi_gt1_lt1p5_counts[iX])*(jet_trk_dphi_gt1_lt1p5_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt1p5_lt2[iFile]->SetBinContent (iX+1, h_jet_trk_dphi_gt1p5_lt2[iFile]->GetBinContent (iX+1) + (ewgt)*(jet_trk_dphi_gt1p5_lt2_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt1p5_lt2_cov[iFile]->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt1p5_lt2_cov[iFile]->GetBinContent (iX+1, iY+1) + (ewgt)*(jet_trk_dphi_gt1p5_lt2_counts[iX])*(jet_trk_dphi_gt1p5_lt2_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt2_lt4[iFile]->SetBinContent (iX+1, h_jet_trk_dphi_gt2_lt4[iFile]->GetBinContent (iX+1) + (ewgt)*(jet_trk_dphi_gt2_lt4_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt2_lt4_cov[iFile]->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt2_lt4_cov[iFile]->GetBinContent (iX+1, iY+1) + (ewgt)*(jet_trk_dphi_gt2_lt4_counts[iX])*(jet_trk_dphi_gt2_lt4_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt4_lt6[iFile]->SetBinContent (iX+1, h_jet_trk_dphi_gt4_lt6[iFile]->GetBinContent (iX+1) + (ewgt)*(jet_trk_dphi_gt4_lt6_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt4_lt6_cov[iFile]->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt4_lt6_cov[iFile]->GetBinContent (iX+1, iY+1) + (ewgt)*(jet_trk_dphi_gt4_lt6_counts[iX])*(jet_trk_dphi_gt4_lt6_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt6_lt8[iFile]->SetBinContent (iX+1, h_jet_trk_dphi_gt6_lt8[iFile]->GetBinContent (iX+1) + (ewgt)*(jet_trk_dphi_gt6_lt8_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt6_lt8_cov[iFile]->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt6_lt8_cov[iFile]->GetBinContent (iX+1, iY+1) + (ewgt)*(jet_trk_dphi_gt6_lt8_counts[iX])*(jet_trk_dphi_gt6_lt8_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt8_lt10[iFile]->SetBinContent (iX+1, h_jet_trk_dphi_gt8_lt10[iFile]->GetBinContent (iX+1) + (ewgt)*(jet_trk_dphi_gt8_lt10_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt8_lt10_cov[iFile]->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt8_lt10_cov[iFile]->GetBinContent (iX+1, iY+1) + (ewgt)*(jet_trk_dphi_gt8_lt10_counts[iX])*(jet_trk_dphi_gt8_lt10_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt10_lt15[iFile]->SetBinContent (iX+1, h_jet_trk_dphi_gt10_lt15[iFile]->GetBinContent (iX+1) + (ewgt)*(jet_trk_dphi_gt10_lt15_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt10_lt15_cov[iFile]->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt10_lt15_cov[iFile]->GetBinContent (iX+1, iY+1) + (ewgt)*(jet_trk_dphi_gt10_lt15_counts[iX])*(jet_trk_dphi_gt10_lt15_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt15_lt20[iFile]->SetBinContent (iX+1, h_jet_trk_dphi_gt15_lt20[iFile]->GetBinContent (iX+1) + (ewgt)*(jet_trk_dphi_gt15_lt20_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt15_lt20_cov[iFile]->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt15_lt20_cov[iFile]->GetBinContent (iX+1, iY+1) + (ewgt)*(jet_trk_dphi_gt15_lt20_counts[iX])*(jet_trk_dphi_gt15_lt20_counts[iY]));
    }
    for (int iX = 0; iX < nDPhiBins; iX++) {
      h_jet_trk_dphi_gt20_lt30[iFile]->SetBinContent (iX+1, h_jet_trk_dphi_gt20_lt30[iFile]->GetBinContent (iX+1) + (ewgt)*(jet_trk_dphi_gt20_lt30_counts[iX]));
      for (int iY = 0; iY < nDPhiBins; iY++)
        h2_jet_trk_dphi_gt20_lt30_cov[iFile]->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt20_lt30_cov[iFile]->GetBinContent (iX+1, iY+1) + (ewgt)*(jet_trk_dphi_gt20_lt30_counts[iX])*(jet_trk_dphi_gt20_lt30_counts[iY]));
    }


    for (int iX = 0; iX < nPtChBins; iX++) {
      h_jet_trk_pt_ns[iFile]->SetBinContent (iX+1, h_jet_trk_pt_ns[iFile]->GetBinContent (iX+1) + (ewgt)*(jet_trk_pt_ns_counts[iX]));
      for (int iY = 0; iY < nPtChBins; iY++)
        h2_jet_trk_pt_ns_cov[iFile]->SetBinContent (iX+1, iY+1, h2_jet_trk_pt_ns_cov[iFile]->GetBinContent (iX+1, iY+1) + (ewgt)*(jet_trk_pt_ns_counts[iX])*(jet_trk_pt_ns_counts[iY]));
    }
    for (int iX = 0; iX < nPtChBins; iX++) {
      h_jet_trk_pt_perp[iFile]->SetBinContent (iX+1, h_jet_trk_pt_perp[iFile]->GetBinContent (iX+1) + (ewgt)*(jet_trk_pt_perp_counts[iX]));
      for (int iY = 0; iY < nPtChBins; iY++)
        h2_jet_trk_pt_perp_cov[iFile]->SetBinContent (iX+1, iY+1, h2_jet_trk_pt_perp_cov[iFile]->GetBinContent (iX+1, iY+1) + (ewgt)*(jet_trk_pt_perp_counts[iX])*(jet_trk_pt_perp_counts[iY]));
    }
    for (int iX = 0; iX < nPtChBins; iX++) {
      h_jet_trk_pt_as[iFile]->SetBinContent (iX+1, h_jet_trk_pt_as[iFile]->GetBinContent (iX+1) + (ewgt)*(jet_trk_pt_as_counts[iX]));
      for (int iY = 0; iY < nPtChBins; iY++)
        h2_jet_trk_pt_as_cov[iFile]->SetBinContent (iX+1, iY+1, h2_jet_trk_pt_as_cov[iFile]->GetBinContent (iX+1, iY+1) + (ewgt)*(jet_trk_pt_as_counts[iX])*(jet_trk_pt_as_counts[iY]));
    }

  } // end loop over events
  std::cout << "Finished event loop." << std::endl;


  SaferDelete (&h2_trk_eff);

  for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {
    SaferDelete (&f_trk_pur[iEta]);
    SaferDelete (&h_trk_pur[iEta]);
  }
  Delete1DArray (f_trk_pur, nEtaTrkBins);
  Delete1DArray (h_trk_pur, nEtaTrkBins);


  for (int iFile = 0; iFile < nFiles; iFile++) {

    outFiles[iFile]->cd ();

    h_evt_counts[iFile]->Write ();
    h_jet_counts[iFile]->Write ();

    h_jet_pt[iFile]->Write ();
    h2_jet_pt_cov[iFile]->Write ();
    h2_jet_eta_phi[iFile]->Write ();

    h_jet_trk_dphi_gt0p5_lt1[iFile]->Write ();
    h2_jet_trk_dphi_gt0p5_lt1_cov[iFile]->Write ();
    h_jet_trk_dphi_gt1_lt1p5[iFile]->Write ();
    h2_jet_trk_dphi_gt1_lt1p5_cov[iFile]->Write ();
    h_jet_trk_dphi_gt1p5_lt2[iFile]->Write ();
    h2_jet_trk_dphi_gt1p5_lt2_cov[iFile]->Write ();
    h_jet_trk_dphi_gt2_lt4[iFile]->Write ();
    h2_jet_trk_dphi_gt2_lt4_cov[iFile]->Write ();
    h_jet_trk_dphi_gt4_lt6[iFile]->Write ();
    h2_jet_trk_dphi_gt4_lt6_cov[iFile]->Write ();
    h_jet_trk_dphi_gt6_lt8[iFile]->Write ();
    h2_jet_trk_dphi_gt6_lt8_cov[iFile]->Write ();
    h_jet_trk_dphi_gt8_lt10[iFile]->Write ();
    h2_jet_trk_dphi_gt8_lt10_cov[iFile]->Write ();
    h_jet_trk_dphi_gt10_lt15[iFile]->Write ();
    h2_jet_trk_dphi_gt10_lt15_cov[iFile]->Write ();
    h_jet_trk_dphi_gt15_lt20[iFile]->Write ();
    h2_jet_trk_dphi_gt15_lt20_cov[iFile]->Write ();
    h_jet_trk_dphi_gt20_lt30[iFile]->Write ();
    h2_jet_trk_dphi_gt20_lt30_cov[iFile]->Write ();

    h_jet_trk_pt_ns[iFile]->Write ();
    h2_jet_trk_pt_ns_cov[iFile]->Write ();
    h_jet_trk_pt_perp[iFile]->Write ();
    h2_jet_trk_pt_perp_cov[iFile]->Write ();
    h_jet_trk_pt_as[iFile]->Write ();
    h2_jet_trk_pt_as_cov[iFile]->Write ();


    outFiles[iFile]->Close ();
  }

  SaferDelete (&jetTrigger);
  SaferDelete (&trackTreeTrigger);

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
  if (!doMixing && (DoMixCatVar1 () || DoMixCatVar2 () || DoMixCatVar3 () || DoMixCatVar4 () || DoMixCatVar5 ())) {
    std::cout << "In RunCorrelator.cxx: Configured to run mixing variations but no mixing configuration specified? Please investigate. Exiting." << std::endl;
    return false;
  }

  // histograms for some mixing variations are symlinks to central values to avoid file copy issues, make sure not to overwrite them by exiting gracefully
  if ((DoMixCatVar2 () || DoMixCatVar4 ()) && Ispp ()) {
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

      if (IsCollisions ()) {
        if ((UseMinBiasTriggers () && !file.Contains ("MinBias")) || (UseJetTriggers () && !file.Contains ("Main")))
          continue;
      }
      if (!file.Contains (jetsInFileName))
        continue;

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
  jet_min_pt = (double) std::atof (argv[argn++]);
  jet_max_pt = (double) std::atof (argv[argn++]);

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
