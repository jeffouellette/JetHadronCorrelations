#ifndef __TriggerPerformance_cxx__
#define __TriggerPerformance_cxx__

#include "TriggerPerformance.h"
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

bool TriggerPerformance (const char* directory,
                         const int dataSet,
                         const char* inFileName) {
 
  std::cout << "Info: In TriggerPerformance.cxx: Entered TriggerPerformance routine." << std::endl;

  SetupDirectories ("TriggerPerformance");

  if (!IsCollisions ()) {
    std::cout << "Info: In TriggerPerformance.cxx: Not configured to evaluate trigger performance in MC, exiting gracefully." << std::endl;
    return true;
  }


  // this identifier here corresponds to the name of the output file
  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  std::cout << "Info: In TriggerPerformance.cxx: File Identifier: " << identifier << std::endl;
  std::cout << "Info: In TriggerPerformance.cxx: Saving output to " << rootPath << std::endl;


  // creates a file identifier pattern that we will use to identify the directory containing the input files
  TString fileIdentifier;
  if (TString (inFileName) == "")
    fileIdentifier = to_string (dataSet);
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
    std::cout << "Info: In TriggerPerformance.cxx: Chain has " << tree->GetListOfFiles ()->GetEntries () << " files, " << tree->GetEntries () << " entries" << std::endl;
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
    std::cout << "Error: In Correlator.C: Invalid or unsupported collision system, exiting." << std::endl;
    return false;
  }


  Trigger** triggers = new Trigger*[jet_trig_n];


  tree->SetBranchAddress ("run_number",     &run_number);
  tree->SetBranchAddress ("event_number",   &event_number);
  tree->SetBranchAddress ("lumi_block",     &lumi_block);
  tree->SetBranchAddress ("isOOTPU",        &isOOTPU);
  tree->SetBranchAddress ("BlayerDesyn",    &BlayerDesyn);
  tree->SetBranchAddress ("passes_toroid",  &passes_toroid);


  tree->SetBranchAddress ("actualInteractionsPerCrossing",  &actualInteractionsPerCrossing);
  tree->SetBranchAddress ("averageInteractionsPerCrossing", &averageInteractionsPerCrossing);


  tree->SetBranchAddress ("nvert",     &nvert);
  tree->SetBranchAddress ("vert_x",    &vert_x);
  tree->SetBranchAddress ("vert_y",    &vert_y);
  tree->SetBranchAddress ("vert_z",    &vert_z);
  tree->SetBranchAddress ("vert_ntrk", &vert_ntrk);
  tree->SetBranchAddress ("vert_type", &vert_type);


  tree->SetBranchAddress ("fcalA_et",       &fcalA_et);
  tree->SetBranchAddress ("fcalC_et",       &fcalC_et);


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


  std::cout << "Info : In TriggerPerformance.cxx: Saving histograms to " << Form ("%s/%s.root", rootPath.Data (), identifier.Data ()) << std::endl;
  TFile* outFile = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");

  TString sys;
  if (Ispp ())        sys = "pp";
  else if (IspPb ())  sys = "pPb";
  else                sys = "???";


  // Initialize a bunch of histograms --
  TH1D** h_jet_trig = Get1DArray <TH1D*> (nTrig); // trigger events
  TH1D** h_jet_all  = Get1DArray <TH1D*> (nTrig); // all events

  for (int iTrig = 0; iTrig < nTrig; iTrig++) {

    h_jet_trig[iTrig] = new TH1D (Form ("h_jet_trig_%s_%s", sys.Data (), triggers[iTrig].name.c_str ()), ";#it{p}_{T}^{leading} [GeV];Counts", nRespBins, respBins);
    h_jet_trig[iTrig]->Sumw2 ();
    h_jet_all[iTrig] = new TH1D (Form ("h_jet_all_%s_%s", sys.Data (), triggers[iTrig].name.c_str ()), ";#it{p}_{T}^{leading} [GeV];Counts", nRespBins, respBins);
    h_jet_all[iTrig]->Sumw2 ();

  } // end loop over iTrig
  

  const JetRadius r0p4 = JetRadius::R0p4;
  const int nEvts = tree->GetEntries ();

  // Loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      std::cout << "Info: In TriggerPerformance.cxx: Event loop " << iEvt / (nEvts / 100) << "\% done...\r" << std::flush;
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


  

    const int jn = GetAktHIJetN (r0p4);
    for (int iJet = 0; iJet < jn; iJet++) {

      if (!MeetsJetAcceptanceCuts (iJet, r0p4))
        continue;

      const float jpt = GetAktHIJetPt (iJet, r0p4);
      const float jeta = GetAktHIJetEta (iJet, r0p4);
      const float jphi = GetAktHIJetPhi (iJet, r0p4);
      const float jen = GetAktHIJetEn (iJet, r0p4);


    } // end loop over jets

  } // end event loop

  std::cout << std::endl << "Info: In TriggerPerformance.cxx: Finished event loop." << std::endl;


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

  } // end loop over iR

  outFile->Write (0, TObject::kOverwrite);
  outFile->Close ();
  SaferDelete (&outFile);

  return true;
}

} // end namespace

#endif
