#ifndef __CentralityAnalysis_cxx__
#define __CentralityAnalysis_cxx__

#include "CentralityAnalysis.h"
#include "Params.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Trigger.h"

#include <Utilities.h>
#include <AtlasUtils.h>

#include <TH2D.h>
#include <TChain.h>
#include <TSystem.h>

#include <iostream>

using namespace std;

namespace JetHadronCorrelations {


bool CentralityAnalysis (const char* directory,
                         const int dataSet,
                         const char* inFileName) {

  cout << "Info: In CentralityAnalysis.cxx: Entered CentralityAnalysis routine." << endl;

  SetupDirectories ("CentralityAnalysis");


  // this identifier here corresponds to the name of the output file
  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  std::cout << "Info: In CentralityAnalysis.cxx: File Identifier: " << identifier << std::endl;
  std::cout << "Info: In CentralityAnalysis.cxx: Saving output to " << rootPath << std::endl;


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


  TString outTreeName = "";
  if (!Is5TeV ()) {
    std::cout << "Error: In CentralityAnalysis.cxx::CentralityAnalysis (const char*, const int, const char*): Unsupported beam collision energy, quitting." << std::endl;
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
    std::cout << "Error: In CentralityAnalysis.cxx::CentralityAnalysis (const char*, const int, const char*): Beam configuration not recognized, quitting." << std::endl;
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


  // Load files for output
  TFile* outFile = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");
  TH2D* h2_mb_Pb_fcal_et_zdc_calibE = new TH2D (Form ("h2_mb_Pb_fcal_et_zdc_calibE_run%i", dataSet), "", 80, -100, 300, 1400, 0, 175);
  TH2D* h2_jet_Pb_fcal_et_zdc_calibE = new TH2D (Form ("h2_jet_Pb_fcal_et_zdc_calibE_run%i", dataSet), "", 80, -100, 300, 1400, 0, 175);

  TH1D* h_mb_Pb_fcal_et = new TH1D (Form ("h_mb_Pb_fcal_et_run%i", dataSet), "", 250, -30, 220);
  h_mb_Pb_fcal_et->Sumw2 ();
  TH1D* h_mb_p_fcal_et = new TH1D (Form ("h_mb_p_fcal_et_run%i", dataSet), "", 250, -30, 220);
  h_mb_p_fcal_et->Sumw2 ();

  TH1D* h_mb_Pb_zdc_calibE = new TH1D (Form ("h_mb_Pb_zdc_calibE_run%i", dataSet), "", 1400, 0, 175);
  h_mb_Pb_zdc_calibE->Sumw2 ();
  TH1D* h_mb_Pb_zdc_calibE_cut = new TH1D (Form ("h_mb_Pb_zdc_calibE_cut_run%i", dataSet), "", 1400, 0, 175);
  h_mb_Pb_zdc_calibE_cut->Sumw2 ();
  TH1D* h_mb_p_zdc_calibE = new TH1D (Form ("h_mb_p_zdc_calibE_run%i", dataSet), "", 1400, 0, 175);
  h_mb_p_zdc_calibE->Sumw2 ();

  TH1D* h_jet_Pb_fcal_et = new TH1D (Form ("h_jet_Pb_fcal_et_run%i", dataSet), "", 250, -30, 220);
  h_jet_Pb_fcal_et->Sumw2 ();
  TH1D* h_jet_p_fcal_et = new TH1D (Form ("h_jet_p_fcal_et_run%i", dataSet), "", 250, -30, 220);
  h_jet_p_fcal_et->Sumw2 ();

  TH1D* h_jet_Pb_zdc_calibE = new TH1D (Form ("h_jet_Pb_zdc_calibE_run%i", dataSet), "", 1400, 0, 175);
  h_jet_Pb_zdc_calibE->Sumw2 ();
  TH1D* h_jet_Pb_zdc_calibE_cut = new TH1D (Form ("h_jet_Pb_zdc_calibE_cut_run%i", dataSet), "", 1400, 0, 175);
  h_jet_Pb_zdc_calibE_cut->Sumw2 ();
  TH1D* h_jet_p_zdc_calibE = new TH1D (Form ("h_jet_p_zdc_calibE_run%i", dataSet), "", 1400, 0, 175);
  h_jet_p_zdc_calibE->Sumw2 ();


  const int nEvts = tree->GetEntries ();

  // First loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In CentralityAnalysis.cxx: Events " << iEvt / (nEvts / 100) << "\% done...\r" << flush;

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
      if (vert_type[iVert] == 3)
        hasPileup = true;
    }
    if (!hasPrimaryVert || hasPileup || fabs (vz) > 150)
      continue;

    float fcal_et_Pb = 0, fcal_et_p = 0, zdc_calibE_Pb = 0, zdc_calibE_p = 0;
    if (IspPb ()) {
      fcal_et_Pb = IsPeriodA () ? fcalA_et : fcalC_et; // Pb-going side is the A side in period A (runs < 313435)
      fcal_et_p = IsPeriodA () ? fcalC_et : fcalA_et;
      zdc_calibE_Pb = IsPeriodA () ? ZdcCalibEnergy_A : ZdcCalibEnergy_C;
      zdc_calibE_p = IsPeriodA () ? ZdcCalibEnergy_C : ZdcCalibEnergy_A;
    }
    else if (Ispp ()) {
      fcal_et_Pb = -999;
      fcal_et_p = fcalA_et + fcalC_et;
      zdc_calibE_Pb = -999;
      zdc_calibE_p = ZdcCalibEnergy_A + ZdcCalibEnergy_C;
    }
    zdc_calibE_Pb *= 1e3;
    zdc_calibE_p *= 1e3;

    if (minbiasTriggers[0]->trigDecision) {
      h2_mb_Pb_fcal_et_zdc_calibE->Fill (fcal_et_Pb, zdc_calibE_Pb);

      h_mb_Pb_fcal_et->Fill (fcal_et_Pb);
      h_mb_p_fcal_et->Fill (fcal_et_p);

      h_mb_Pb_zdc_calibE->Fill (zdc_calibE_Pb);
      if (zdc_calibE_Pb > 10)
        h_mb_Pb_zdc_calibE_cut->Fill (zdc_calibE_Pb);
      h_mb_p_zdc_calibE->Fill (zdc_calibE_p);
    }

    if (jetTriggers[0]->trigDecision) {
      h2_jet_Pb_fcal_et_zdc_calibE->Fill (fcal_et_Pb, zdc_calibE_Pb);

      h_jet_Pb_fcal_et->Fill (fcal_et_Pb);
      h_jet_p_fcal_et->Fill (fcal_et_p);

      h_jet_Pb_zdc_calibE->Fill (zdc_calibE_Pb);
      if (zdc_calibE_Pb > 10)
        h_jet_Pb_zdc_calibE_cut->Fill (zdc_calibE_Pb);
      h_jet_p_zdc_calibE->Fill (zdc_calibE_p);
    }
  } // end event loop
  cout << endl << "Info: In CentralityAnalysis.cxx: Finished processing events." << endl;

  h2_mb_Pb_fcal_et_zdc_calibE->Write ();
  h_mb_Pb_fcal_et->Write ();
  h_mb_p_fcal_et->Write ();
  h_mb_Pb_zdc_calibE->Write ();
  h_mb_Pb_zdc_calibE_cut->Write ();
  h_mb_p_zdc_calibE->Write ();

  h2_jet_Pb_fcal_et_zdc_calibE->Write ();
  h_jet_Pb_fcal_et->Write ();
  h_jet_p_fcal_et->Write ();
  h_jet_Pb_zdc_calibE->Write ();
  h_jet_Pb_zdc_calibE_cut->Write ();
  h_jet_p_zdc_calibE->Write ();

  outFile->Close ();
  SaferDelete (&outFile);

  return true;
}

} // end namespace

#endif
