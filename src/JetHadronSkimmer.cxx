#ifndef __JetHadronSkimmer_cxx__
#define __JetHadronSkimmer_cxx__

#include "JetHadronSkimmer.h"
#include "Params.h"
#include "TreeVariables.h"
#include "OutTree.h"
#include "LocalUtilities.h"
#include "Trigger.h"

#include <Utilities.h>
#include <AtlasUtils.h>

#include <TChain.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TLorentzVector.h>

#include <iostream>
#include <random>

using namespace std;

namespace JetHadronCorrelations {


bool JetHadronSkimmer (const char* directory,
                       const int dataSet,
                       const char* inFileName) {
 
  std::cout << "Info: In JetHadronSkimmer.cxx: Entered JetHadronSkimmer routine." << std::endl;
  std::cout << "Info: In JetHadronSkimmer.cxx: Printing systematic onfiguration:";
  std::cout << std::endl;

  if (IsHijing ())
    std::cout << "Info: In JetHadronSkimmer.cxx: File detected as Hijing, will not check for data quality whatsoever" << std::endl;
  else if (IsDataOverlay ())
    std::cout << "Info: In JetHadronSkimmer.cxx: File detected as data overlay, will check data quality (i.e. pile-up & detector defects)" << std::endl;

  if (!IsCollisions () && IsHijing ())
    SetupDirectories ("MC/Hijing", false);
  else if (!IsCollisions ())
    SetupDirectories ("MC");
  else {
    if (UseJet50GeVTriggers ())
      SetupDirectories ("50GeVJetsData");
    else if (UseJet100GeVTriggers ())
      SetupDirectories ("100GeVJetsData");
    else if (UseMinBiasTriggers ())
      SetupDirectories ("MinBiasData");
  }
    

  // this identifier here corresponds to the name of the output file
  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  std::cout << "Info: In JetHadronSkimmer.cxx: File Identifier: " << identifier << std::endl;
  std::cout << "Info: In JetHadronSkimmer.cxx: Saving output to " << rootPath << std::endl;


  // creates a file identifier pattern that we will use to identify the directory containing the input files
  TString fileIdentifier;
  if (TString (inFileName) == "") {
    if (IsCollisions ()) {
      fileIdentifier = to_string (dataSet);
    }
    else {
      std::cout << "Error: In JetHadronSkimmer.cxx: Cannot identify this MC file! Quitting." << std::endl;
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
    std::cout << "Info: In JetHadronSkimmer.cxx: Chain has " << tree->GetListOfFiles ()->GetEntries () << " files, " << tree->GetEntries () << " entries" << std::endl;
  }


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
  tree->SetBranchAddress ("fcalA_et_Cos2",  &fcalA_et_Cos2);
  tree->SetBranchAddress ("fcalC_et_Cos2",  &fcalC_et_Cos2);
  tree->SetBranchAddress ("fcalA_et_Sin2",  &fcalA_et_Sin2);
  tree->SetBranchAddress ("fcalC_et_Sin2",  &fcalC_et_Sin2);
  tree->SetBranchAddress ("fcalA_et_Cos3",  &fcalA_et_Cos3);
  tree->SetBranchAddress ("fcalC_et_Cos3",  &fcalC_et_Cos3);
  tree->SetBranchAddress ("fcalA_et_Sin3",  &fcalA_et_Sin3);
  tree->SetBranchAddress ("fcalC_et_Sin3",  &fcalC_et_Sin3);
  tree->SetBranchAddress ("fcalA_et_Cos4",  &fcalA_et_Cos4);
  tree->SetBranchAddress ("fcalC_et_Cos4",  &fcalC_et_Cos4);
  tree->SetBranchAddress ("fcalA_et_Sin4",  &fcalA_et_Sin4);
  tree->SetBranchAddress ("fcalC_et_Sin4",  &fcalC_et_Sin4);


  if (!Ispp ()) {
    tree->SetBranchAddress ("ZdcCalibEnergy_A",   &ZdcCalibEnergy_A);
    tree->SetBranchAddress ("ZdcCalibEnergy_C",   &ZdcCalibEnergy_C);
    tree->SetBranchAddress ("ZdcRawEnergy_A",     &ZdcRawEnergy_A);
    tree->SetBranchAddress ("ZdcRawEnergy_C",     &ZdcRawEnergy_C);
  }


  if (IsPbPb ()) {
    tree->SetBranchAddress ("cluster_sumGap_A",  &cluster_sumGap_A);
    tree->SetBranchAddress ("cluster_sumGap_C",  &cluster_sumGap_C);
    tree->SetBranchAddress ("cluster_edgeGap_A", &cluster_edgeGap_A);
    tree->SetBranchAddress ("cluster_edgeGap_C", &cluster_edgeGap_C);
    tree->SetBranchAddress ("sumGap_A",          &sumGap_A);
    tree->SetBranchAddress ("sumGap_C",          &sumGap_C);
    tree->SetBranchAddress ("edgeGap_A",         &edgeGap_A);
    tree->SetBranchAddress ("edgeGap_C",         &edgeGap_C);
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
  


  TString outTreeName = "";
  if (!Is5TeV ()) {
    std::cout << "Error: In JetHadronSkimmer.cxx::JetHadronSkimmer (const char*, const int, const char*): Unsupported beam collision energy, quitting." << std::endl;
    return false;
  }
  else if (IspPb16 ()) {
    outTreeName = "pPbTree";
    minbias_trig_n = minbias_trig_n_pPb16;
    minbias_trig_name = minbias_trig_name_pPb16;
    jet_trig_n = jet_trig_n_pPb16s5TeV;
    jet_trig_name = jet_trig_name_pPb16s5TeV;
  }
  else if (Ispp17 ()) {
    outTreeName = "ppTree";
    minbias_trig_n = minbias_trig_n_pp17;
    minbias_trig_name = minbias_trig_name_pp17;
    jet_trig_n = jet_trig_n_pp17;
    jet_trig_name = jet_trig_name_pp17;
  }
  else if (IsPbPb18 ()) {
    outTreeName = "PbPbTree";
    minbias_trig_n = minbias_trig_n_PbPb18;
    minbias_trig_name = minbias_trig_name_PbPb18;
    jet_trig_n = jet_trig_n_PbPb18;
    jet_trig_name = jet_trig_name_PbPb18;
  }
  else {
    std::cout << "Error: In JetHadronSkimmer.cxx::JetHadronSkimmer (const char*, const int, const char*): Beam configuration not recognized, quitting." << std::endl;
    return false;
  }


  Trigger* jetTriggers[jet_trig_n] = {};
  Trigger* minbiasTriggers[minbias_trig_n] = {};

  if (IsCollisions ()) {
    for (int iTrig = 0; iTrig < jet_trig_n; iTrig++) {
      jetTriggers[iTrig] = new Trigger (jet_trig_name[iTrig]);
      tree->SetBranchAddress ((jet_trig_name[iTrig]+"_decision").c_str (), &(jetTriggers[iTrig]->trigDecision));
      tree->SetBranchAddress ((jet_trig_name[iTrig]+"_prescale").c_str (), &(jetTriggers[iTrig]->trigPrescale));
    }
    for (int iTrig = 0; iTrig < minbias_trig_n; iTrig++) {
      minbiasTriggers[iTrig] = new Trigger (minbias_trig_name[iTrig]);
      tree->SetBranchAddress ((minbias_trig_name[iTrig]+"_decision").c_str (), &(minbiasTriggers[iTrig]->trigDecision));
      tree->SetBranchAddress ((minbias_trig_name[iTrig]+"_prescale").c_str (), &(minbiasTriggers[iTrig]->trigPrescale));
    }
  }


  centBins = (DoFcalCentVar () ? fcalCentBins : zdcCentBins);
  numCentBins = (DoFcalCentVar () ? numFcalCentBins : numZdcCentBins);


  // Load files for output
  const int numFileBins = (Ispp () ? 1 : numCentBins);
  TFile* outFiles[numFileBins];
  OutTree* outTrees[numFileBins];
  for (int iCent = 0; iCent < numFileBins; iCent++) {
    TString fName = (numFileBins == 1 ? Form ("%s/%s.root", rootPath.Data (), identifier.Data ()) : Form ("%s/%s_iCent%i.root", rootPath.Data (), identifier.Data (), iCent));
    outFiles[iCent] = new TFile (fName.Data (), "recreate");
    outFiles[iCent]->Delete (Form ("%s;*", outTreeName.Data ()));

    outTrees[iCent] = new OutTree (outTreeName.Data (), outFiles[iCent]);
    outTrees[iCent]->SetBranchEventInfo ();
    outTrees[iCent]->SetBranchJets ();
    outTrees[iCent]->SetBranchTracks ();
    outTrees[iCent]->SetBranches ();
  }


  const int nEvts = tree->GetEntries ();

  std::random_device rndm;
  std::mt19937 mt19937_gen (rndm ());
  std::normal_distribution <double> jetES2PercSmearDist (1.0, 0.02);
  std::normal_distribution <double> jetES5PercSmearDist (1.0, 0.05);

  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      std::cout << "Info: In JetHadronSkimmer.cxx: Events " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    tree->GetEntry (iEvt);

    if (IsCollisions ()) {
      event_weight = -1;
      if (UseJet50GeVTriggers () && jetTriggers[0]->trigDecision)
        event_weight = jetTriggers[0]->trigPrescale;
      else if (UseJet100GeVTriggers () && jetTriggers[1]->trigDecision)
        event_weight = jetTriggers[1]->trigPrescale;
      else if (UseMinBiasTriggers () && minbiasTriggers[0]->trigDecision)
        event_weight = minbiasTriggers[0]->trigPrescale;
    }
    else {
      event_weight = crossSectionPicoBarns * mcFilterEfficiency * GetJetLuminosity (); // sigma * f * L_int
      //event_weight = mcEventWeights->at (0) * mcFilterEfficiency * GetJetLuminosity (); // sigma * f * L_int
      if (IsHijing ()) {
        assert (truth_event_n > 0);
        ip = impactParameter[0];
        eventPlane = eventPlaneAngle[0];
      }
    }

    if (!IsHijing ()) {
      bool hasPrimary = false;
      bool hasPileup = false;
      vz = -999;
      for (int iVert = 0; iVert < nvert; iVert++) {
        if (vert_type[iVert] == 1) {
          hasPrimary = true;
          vz = vert_z[iVert];
        }
        else if (vert_type[iVert] == 3)
          hasPileup = true;
      }
      if (!hasPrimary || (!DoWithPileupVar () && hasPileup) || fabs (vz) > 150)
        continue;
    }
    else {
      if (nvert > 0) vz = vert_z[0];
      else vz = 0;
    }

    short iCent = -1;
    if (IspPb ()) {
      fcal_et_Pb = IsPeriodA () ? fcalA_et : fcalC_et; // Pb-going side is the A side in period A (all 5.02 TeV runs and 8.16 TeV runs before or including 313435)
      fcal_et_p = IsPeriodA () ? fcalC_et : fcalA_et;
      q2x_Pb = IsPeriodA () ? fcalA_et_Cos2 : fcalC_et_Cos2;
      q2y_Pb = IsPeriodA () ? fcalA_et_Sin2 : fcalC_et_Sin2;
      q2x_p = IsPeriodA () ? fcalC_et_Cos2 : fcalA_et_Cos2;
      q2y_p = IsPeriodA () ? fcalC_et_Sin2 : fcalA_et_Sin2;
      zdc_calibE_Pb = IsPeriodA () ? ZdcCalibEnergy_A : ZdcCalibEnergy_C;
      zdc_calibE_p = IsPeriodA () ? ZdcCalibEnergy_C : ZdcCalibEnergy_A;
      zdc_calibE_Pb *= 1e3;
      zdc_calibE_p *= 1e3;
      iCent = GetCentBin (DoFcalCentVar () ? fcal_et_Pb : zdc_calibE_Pb);
    }
    else if (Ispp ()) {
      fcal_et_Pb = -999;
      fcal_et_p = fcalA_et + fcalC_et;
      q2x_Pb = -999;
      q2y_Pb = -999;
      q2x_p = fcalA_et_Cos2 + fcalC_et_Cos2;
      q2y_p = fcalA_et_Sin2 + fcalC_et_Sin2;
      zdc_calibE_Pb = -999;
      zdc_calibE_p = 0;
      iCent = 0;
    }

    if (iCent < 0 || iCent > numCentBins-1)
      continue;

    if (event_weight != -1) { // trigger requirement

      out_akt4_hi_jet_n = 0;
      for (int iJ = 0; iJ < akt4_hi_jet_n; iJ++) {
        if (!MeetsJetAcceptanceCuts (iJ))
          continue;

        // Crude systematic -- smear jet energy scale by 2% or 5% (ad-hoc)
        double jpt = akt4_hi_jet_pt_xcalib[iJ];
        double jen = akt4_hi_jet_e_xcalib[iJ]; 

        if (DoJetES5PercUpVar ())    { jpt *= 1.05;  jen *= 1.05; }
        if (DoJetES5PercDownVar ())  { jpt *= 0.95;  jen *= 0.95; }
        if (DoJetES5PercSmearVar ()) { 
          double sf = jetES5PercSmearDist (mt19937_gen);
          jpt *= sf;
          jen *= sf;
        }
        if (DoJetES2PercUpVar ())    { jpt *= 1.02;  jen *= 1.02; }
        if (DoJetES2PercDownVar ())  { jpt *= 0.98;  jen *= 0.98; }
        if (DoJetES2PercSmearVar ()) { 
          double sf = jetES2PercSmearDist (mt19937_gen);
          jpt *= sf;
          jen *= sf;
        }

        out_akt4_hi_jet_pt[out_akt4_hi_jet_n] = jpt;
        out_akt4_hi_jet_eta[out_akt4_hi_jet_n] = akt4_hi_jet_eta_xcalib[iJ];
        out_akt4_hi_jet_phi[out_akt4_hi_jet_n] = akt4_hi_jet_phi[iJ];
        out_akt4_hi_jet_e[out_akt4_hi_jet_n] = jen;
        out_akt4_hi_jet_n++;
      }

      int lJ = -1, slJ = -1;
      for (int iJ = 0; iJ < out_akt4_hi_jet_n; iJ++) {
        if (lJ == -1 || out_akt4_hi_jet_pt[lJ] < out_akt4_hi_jet_pt[iJ])
          lJ = iJ;
      } // finds leading jet
      for (int iJ = 0; iJ < out_akt4_hi_jet_n; iJ++) {
        if (iJ == lJ)
          continue;
        if (slJ == -1 || out_akt4_hi_jet_pt[slJ] < out_akt4_hi_jet_pt[iJ])
          slJ = iJ;
      } // finds subleading jet

      leading_jet = lJ;
      subleading_jet = slJ;

      if (UseJetTriggers () && leading_jet == -1)
        continue; // skip jet triggered events with no leading jet


      out_trk_n = 0;
      for (int iTrk = 0; iTrk < trk_n; iTrk++) {
        if (!MeetsTrackCuts (iTrk))
          continue;

        out_trk_pt[out_trk_n] = trk_pt[iTrk];
        out_trk_eta[out_trk_n] = trk_eta[iTrk];
        out_trk_phi[out_trk_n] = trk_phi[iTrk];
        out_trk_charge[out_trk_n] = trk_charge[iTrk];
        out_trk_d0[out_trk_n] = trk_d0[iTrk];
        out_trk_z0[out_trk_n] = trk_z0[iTrk];
        if (!IsCollisions ())
          out_trk_truth_matched[out_trk_n] = (trk_prob_truth[iTrk] > 0.5);
        out_trk_n++;
      }

      float sumptcw = 0, sumptccw = 0;
      for (int iTrk = 0; iTrk < out_trk_n; iTrk++) {
        float dphi = DeltaPhi (out_trk_phi[iTrk], out_akt4_hi_jet_phi[leading_jet], true);
        if (pi/3. < fabs (dphi) && fabs (dphi) < 2.*pi/3.) {
          if (dphi < 0) // then trackphi is ccw of zphi
            sumptccw += out_trk_pt[iTrk];
          else
            sumptcw += out_trk_pt[iTrk];
        }
      } // end loop over tracks

      if (sumptccw > sumptcw) {
        leading_jet_phi_transmax = out_akt4_hi_jet_phi[leading_jet] + pi/2.;
        leading_jet_phi_transmin = out_akt4_hi_jet_phi[leading_jet] - pi/2.;
      }
      else {
        leading_jet_phi_transmax = out_akt4_hi_jet_phi[leading_jet] - pi/2.;
        leading_jet_phi_transmin = out_akt4_hi_jet_phi[leading_jet] + pi/2.;
      }

      outTrees[iCent]->Fill ();

    } // end jet trigger scope

  } // end event selection
  std::cout << std::endl << "Info: In JetHadronSkimmer.cxx: Finished processing events." << std::endl;

  if (IspPb ()) {
    for (int iCent = 0; iCent < numFileBins; iCent++) {
      outFiles[iCent]->Write (0, TObject::kOverwrite);
      outFiles[iCent]->Close ();
      SaferDelete (&(outFiles[iCent]));
    }
  }
  else if (Ispp ()) {
    outFiles[0]->Write ();//0, TObject::kOverwrite);
    outFiles[0]->Close ();
    SaferDelete (&(outFiles[0]));
  }


  return true;
}

} // end namespace

#endif
