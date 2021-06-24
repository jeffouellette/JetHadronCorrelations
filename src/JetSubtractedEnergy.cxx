#ifndef __JetSubtractedEnergy_cxx__
#define __JetSubtractedEnergy_cxx__

#include "JetSubtractedEnergy.h"
#include "Params.h"
#include "CentralityDefs.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Trigger.h"

#include <Utilities.h>
#include <AtlasUtils.h>

#include <TH2D.h>
#include <TH3D.h>
#include <TChain.h>
#include <TSystem.h>

#include <iostream>

using namespace std;

namespace JetHadronCorrelations {


bool JetSubtractedEnergy (const char* directory,
                          const int dataSet,
                          const char* inFileName) {

  cout << "Info: In JetSubtractedEnergy.cxx: Entered JetSubtractedEnergy routine." << endl;

  SetupDirectories ("JetSubtractedEnergy");


  // this identifier here corresponds to the name of the output file
  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  std::cout << "Info: In JetSubtractedEnergy.cxx: File Identifier: " << identifier << std::endl;
  std::cout << "Info: In JetSubtractedEnergy.cxx: Saving output to " << rootPath << std::endl;


  // creates a file identifier pattern that we will use to identify the directory containing the input files
  TString fileIdentifier;
  if (TString (inFileName) == "") {
    if (IsCollisions ()) {
      fileIdentifier = to_string (dataSet);
    }
    else {
      std::cout << "Error: In JetHadronSkimmer.C: Cannot identify this MC file! Quitting." << std::endl;
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
    std::cout << "Chain has " << tree->GetListOfFiles ()->GetEntries () << " files, " << tree->GetEntries () << " entries" << std::endl;
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


  if (!Ispp ()) {
    tree->SetBranchAddress ("ZdcCalibEnergy_A",   &ZdcCalibEnergy_A);
    tree->SetBranchAddress ("ZdcCalibEnergy_C",   &ZdcCalibEnergy_C);
    tree->SetBranchAddress ("ZdcRawEnergy_A",     &ZdcRawEnergy_A);
    tree->SetBranchAddress ("ZdcRawEnergy_C",     &ZdcRawEnergy_C);
  }


  tree->SetBranchAddress ("cluster_sumGap_A",  &cluster_sumGap_A);
  tree->SetBranchAddress ("cluster_sumGap_C",  &cluster_sumGap_C);
  tree->SetBranchAddress ("cluster_edgeGap_A", &cluster_edgeGap_A);
  tree->SetBranchAddress ("cluster_edgeGap_C", &cluster_edgeGap_C);
  tree->SetBranchAddress ("sumGap_A",          &sumGap_A);
  tree->SetBranchAddress ("sumGap_C",          &sumGap_C);
  tree->SetBranchAddress ("edgeGap_A",         &edgeGap_A);
  tree->SetBranchAddress ("edgeGap_C",         &edgeGap_C);


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


  if (!Is5TeV ()) {
    std::cout << "Error: In JetSubtractedEnergy.cxx::JetSubtractedEnergy (const char*, const int, const char*): Unsupported beam collision energy, quitting." << std::endl;
    return false;
  }
  else if (IspPb16 ()) {
    jet_trig_n = jet_trig_n_pPb16s5TeV;
    jet_trig_name = jet_trig_name_pPb16s5TeV;
  }
  else if (Ispp17 ()) {
    jet_trig_n = jet_trig_n_pp17;
    jet_trig_name = jet_trig_name_pp17;
  }
  else if (IsPbPb18 ()) {
    jet_trig_n = jet_trig_n_PbPb18;
    jet_trig_name = jet_trig_name_PbPb18;
  }
  else {
    std::cout << "Error: In JetSubtractedEnergy.cxx::JetSubtractedEnergy (const char*, const int, const char*): Beam configuration not recognized, quitting." << std::endl;
    return false;
  }


  Trigger* jetTriggers[jet_trig_n];
  Trigger* zdcL1Triggers[zdc_L1_trig_n];

  if (IsCollisions ()) {
    for (int iTrig = 0; iTrig < jet_trig_n; iTrig++) {
      jetTriggers[iTrig] = new Trigger (jet_trig_name[iTrig]);
      tree->SetBranchAddress ((jet_trig_name[iTrig]+"_decision").c_str (), &(jetTriggers[iTrig]->trigDecision));
      tree->SetBranchAddress ((jet_trig_name[iTrig]+"_prescale").c_str (), &(jetTriggers[iTrig]->trigPrescale));
    }
    if (!Ispp ()) {
      for (int iTrig = 0; iTrig < zdc_L1_trig_n; iTrig++) {
        zdcL1Triggers[iTrig] = new Trigger (zdc_L1_trig_name[iTrig]);
        tree->SetBranchAddress ((zdc_L1_trig_name[iTrig]+"_decision").c_str (), &(zdcL1Triggers[iTrig]->trigDecision));
        tree->SetBranchAddress ((zdc_L1_trig_name[iTrig]+"_prescale").c_str (), &(zdcL1Triggers[iTrig]->trigPrescale));
      }
    }
  }


  // Load files for output
  TFile* outFile = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");

  TH2D* h2_zdc_calibE_jet_subEt = new TH2D ("h2_zdc_calibE_jet_subEt", "", 128, 0, 160, 150, 0, 15);
  TH2D* h2_zdc_calibE_jet_subE  = new TH2D ("h2_zdc_calibE_jet_subE", "", 128, 0, 160, 150, 0, 100);
  TH2D* h2_Pb_fcal_et_jet_subEt = new TH2D ("h2_Pb_fcal_et_jet_subEt", "", 125, -30, 220, 150, 0, 15);
  TH2D* h2_Pb_fcal_et_jet_subE  = new TH2D ("h2_Pb_fcal_et_jet_subE", "", 125, -30, 220, 150, 0, 100);

  const int nCentBins = (DoFcalCentVar () ? nFcalCentBins : (DoFineFcalCentVar () ? nFineFcalCentBins : nZdcCentBins));
  const int nFileBins = (Ispp () ? 1 : nCentBins);
  TH1D** h_jet_subEt = new TH1D* [nFileBins];
  TH1D** h_jet_subE = new TH1D* [nFileBins];
  for (int iFile = 0; iFile < nFileBins; iFile++) {
    h_jet_subEt[iFile] = new TH1D (Form ("h_jet_subEt%s", nFileBins == 1 ? "" : Form ("_iCent%i", iFile)), "", 150, 0, 15);
    h_jet_subE[iFile] = new TH1D (Form ("h_jet_subE%s", nFileBins == 1 ? "" : Form ("_iCent%i", iFile)), "", 150, 0, 100);
  }


  const int nEvts = tree->GetEntries ();

  // First loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In JetSubtractedEnergy.cxx: Events " << iEvt / (nEvts / 100) << "\% done...\r" << flush;

    tree->GetEntry (iEvt);

    if (IsPbPb18 () && (IsCollisions () || IsDataOverlay ()) && BlayerDesyn)
      continue; // check for B layer desynchronization
    if (IsPbPb () && (IsCollisions () || IsDataOverlay ()) && isOOTPU)
      continue; // check for out-of-time pile-up

    bool hasPrimaryVert = false;
    bool hasPileup = false;
    float vz = -999;
    for (int iVert = 0; iVert < nvert; iVert++) {
      if (vert_type[iVert] == 1) {
        hasPrimaryVert = true;
        vz = vert_z[iVert];
      }
      if (vert_type[iVert] == 3 && vert_ntrk[iVert] > 6)
        hasPileup = true;
    }
    if (!hasPrimaryVert || hasPileup || fabs (vz) > 150)
      continue;

    float fcal_et_Pb = 0, fcal_et_p = 0;
    float zdc_calibE_Pb = 0, zdc_calibE_p = 0;
    //float edgeGap_Pb = 0, edgeGap_p = 0;
    //bool zdc_Pb_decision = false, zdc_p_decision = false;

    short iCent = -1;

    if (IspPb ()) {
      fcal_et_Pb = IsPeriodA () ? fcalA_et : fcalC_et; // Pb-going side is the A side in period A (runs < 313435)
      //fcal_et_p = IsPeriodA () ? fcalC_et : fcalA_et;

      zdc_calibE_Pb = (IsPeriodA () ? ZdcCalibEnergy_A : ZdcCalibEnergy_C) * 1e3;
      //zdc_calibE_p = (IsPeriodA () ? ZdcCalibEnergy_C : ZdcCalibEnergy_A) * 1e3;

      //edgeGap_Pb = IsPeriodA () ? edgeGap_A : edgeGap_C;
      //edgeGap_p = IsPeriodA () ? edgeGap_C : edgeGap_A;

      //zdc_Pb_decision = IsPeriodA () ? zdcL1Triggers[0]->trigDecision : zdcL1Triggers[1]->trigDecision;
      //zdc_p_decision = IsPeriodA () ? zdcL1Triggers[1]->trigDecision : zdcL1Triggers[0]->trigDecision;

      if (DoFcalCentVar ())
        iCent = GetBin (fcalCentBins, nFcalCentBins, fcal_et_Pb);
      else if (DoFineFcalCentVar ())
        iCent = GetBin (fineFcalCentBins, nFineFcalCentBins, fcal_et_Pb);
      else
        iCent = GetBin (zdcCentBins, nZdcCentBins, zdc_calibE_Pb);
    }
    else if (Ispp ()) {
      //fcal_et_p = fcalA_et + fcalC_et;
      fcal_et_Pb = -999;

      //zdc_calibE_p = ZdcCalibEnergy_A + ZdcCalibEnergy_C;
      zdc_calibE_Pb = -999;

      //edgeGap_p = -999;
      //edgeGap_Pb = -999;

      iCent = 0;
    }

    if (iCent < 0 || iCent > nCentBins-1)
      continue;


    if (jetTriggers[0]->trigDecision) {

      for (int iJ = 0; iJ < akt4_hi_jet_n; iJ++) {
        if (!MeetsJetAcceptanceCuts (iJ))
          continue;

        if (akt4_hi_jet_pt_xcalib[iJ] < 60)
          continue;

        if (!Ispp ()) {
          h2_zdc_calibE_jet_subEt->Fill (zdc_calibE_Pb, akt4_hi_jet_sub_et[iJ]);
          h2_zdc_calibE_jet_subE->Fill (zdc_calibE_Pb, akt4_hi_jet_sub_e[iJ]);
          h2_Pb_fcal_et_jet_subEt->Fill (fcal_et_Pb, akt4_hi_jet_sub_et[iJ]);
          h2_Pb_fcal_et_jet_subE->Fill (fcal_et_Pb, akt4_hi_jet_sub_e[iJ]);
        }

        h_jet_subEt[iCent]->Fill (akt4_hi_jet_sub_et[iJ]);
        h_jet_subE[iCent]->Fill (akt4_hi_jet_sub_e[iJ]);
      }
    }
  } // end event loop
  cout << endl << "Info: In JetSubtractedEnergy.cxx: Finished processing events." << endl;


  if (IsCollisions ()) {
    for (int iTrig = 0; iTrig < jet_trig_n; iTrig++)
      SaferDelete (&jetTriggers[iTrig]);
    if (!Ispp ()) {
      for (int iTrig = 0; iTrig < zdc_L1_trig_n; iTrig++)
        SaferDelete (&zdcL1Triggers[iTrig]);
    }
  }
  SaferDelete (&tree);


  outFile->cd ();

  if (!Ispp ()) {
    h2_zdc_calibE_jet_subEt->Write ();
    h2_zdc_calibE_jet_subE->Write ();
    h2_Pb_fcal_et_jet_subEt->Write ();
    h2_Pb_fcal_et_jet_subE->Write ();
  }
  SaferDelete (&h2_zdc_calibE_jet_subEt);
  SaferDelete (&h2_zdc_calibE_jet_subE);
  SaferDelete (&h2_Pb_fcal_et_jet_subEt);
  SaferDelete (&h2_Pb_fcal_et_jet_subE);

  for (int iFile = 0; iFile < nFileBins; iFile++) {
    h_jet_subEt[iFile]->Write ();
    SaferDelete (&h_jet_subEt[iFile]);
    h_jet_subE[iFile]->Write ();
    SaferDelete (&h_jet_subE[iFile]);
  }
  delete [] h_jet_subEt;
  delete [] h_jet_subE;

  outFile->Close ();
  SaferDelete (&outFile);

  return true;
}

} // end namespace

#endif
