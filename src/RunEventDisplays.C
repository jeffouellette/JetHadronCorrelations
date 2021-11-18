#ifndef __RunCorrelator_C__
#define __RunCorrelator_C__

#include "Params.h"
#include "OutTree.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"

#include <Utilities.h>

#include <TChain.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLorentzVector.h>

#include <iostream>
#include <math.h>

using namespace JetHadronCorrelations;

bool doMixing = false;


void Correlator (const char* tag, const char* outFileName, TTree* jetsTree, TTree* tracksTree = nullptr) {

  // get event tagging & jet information from main tree
  jetsTree->SetBranchAddress ("event_number",  &event_number);
  jetsTree->SetBranchAddress ("lumi_block",    &lumi_block);
  jetsTree->SetBranchAddress ("run_number",    &run_number);

  jetsTree->SetBranchAddress ("event_weight",  &event_weight); 

  jetsTree->SetBranchAddress ("fcal_et_Pb",     &fcal_et_Pb);
  jetsTree->SetBranchAddress ("fcal_et_p",      &fcal_et_p);
  jetsTree->SetBranchAddress ("zdc_calibE_Pb",  &zdc_calibE_Pb);
  jetsTree->SetBranchAddress ("zdc_calibE_p",   &zdc_calibE_p);


  jetsTree->SetBranchAddress ("leading_jet",                &leading_jet);
  jetsTree->SetBranchAddress ("subleading_jet",             &subleading_jet);
  jetsTree->SetBranchAddress ("akt4_hi_jet_n",              &akt4_hi_jet_n);
  jetsTree->SetBranchAddress ("akt4_hi_jet_pt",             &akt4_hi_jet_pt);
  jetsTree->SetBranchAddress ("akt4_hi_jet_eta",            &akt4_hi_jet_eta);
  jetsTree->SetBranchAddress ("akt4_hi_jet_phi",            &akt4_hi_jet_phi);
  jetsTree->SetBranchAddress ("akt4_hi_jet_e",              &akt4_hi_jet_e);


  jetsTree->SetBranchAddress ("trk_n",          &trk_n);
  jetsTree->SetBranchAddress ("trk_pt",         &trk_pt);
  jetsTree->SetBranchAddress ("trk_eta",        &trk_eta);
  jetsTree->SetBranchAddress ("trk_phi",        &trk_phi);


  TCanvas* canvas = new TCanvas ("canvas", "", 800, 800);

  TEllipse* circ = new TEllipse ();
  circ->SetLineWidth (2);
  circ->SetLineColor (kGrey+3);
  

  const int nEvts = jetsTree->GetEntries ();


  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)                                
      cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;  


    jetsTree->GetEntry (iEvt % jetsTree->GetEntries ());

    if (leading_jet == -1) {
      continue; // require a leading jet (i.e. any jet at all)
    }

    const double yboost = GetBoost (run_number);


    for (int iJet = 0; iJet < akt4_hi_jet_n; iJet++) {
      if (!MeetsJetPtCut (akt4_hi_jet_pt[iJet]))
        continue;


      for (int iTrk = 0; iTrk < trk_n; iTrk++) {
        if (trk_pt[iTrk] < pTChBins[0])
          continue;

      }

      for (int iX = 0; iX < nDPhiBins; iX++) {
        h_jet_trk_dphi_gt0p5_lt1->SetBinContent (iX+1, h_jet_trk_dphi_gt0p5_lt1->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt0p5_lt1_counts[iX]));
        for (int iY = 0; iY < nDPhiBins; iY++)
          h2_jet_trk_dphi_gt0p5_lt1_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt0p5_lt1_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt0p5_lt1_counts[iX])*(jet_trk_dphi_gt0p5_lt1_counts[iY]));
      }
      for (int iX = 0; iX < nDPhiBins; iX++) {
        h_jet_trk_dphi_gt1_lt1p5->SetBinContent (iX+1, h_jet_trk_dphi_gt1_lt1p5->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt1_lt1p5_counts[iX]));
        for (int iY = 0; iY < nDPhiBins; iY++)
          h2_jet_trk_dphi_gt1_lt1p5_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt1_lt1p5_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt1_lt1p5_counts[iX])*(jet_trk_dphi_gt1_lt1p5_counts[iY]));
      }
      for (int iX = 0; iX < nDPhiBins; iX++) {
        h_jet_trk_dphi_gt1p5_lt2->SetBinContent (iX+1, h_jet_trk_dphi_gt1p5_lt2->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt1p5_lt2_counts[iX]));
        for (int iY = 0; iY < nDPhiBins; iY++)
          h2_jet_trk_dphi_gt1p5_lt2_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt1p5_lt2_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt1p5_lt2_counts[iX])*(jet_trk_dphi_gt1p5_lt2_counts[iY]));
      }
      for (int iX = 0; iX < nDPhiBins; iX++) {
        h_jet_trk_dphi_gt2_lt4->SetBinContent (iX+1, h_jet_trk_dphi_gt2_lt4->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt2_lt4_counts[iX]));
        for (int iY = 0; iY < nDPhiBins; iY++)
          h2_jet_trk_dphi_gt2_lt4_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt2_lt4_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt2_lt4_counts[iX])*(jet_trk_dphi_gt2_lt4_counts[iY]));
      }
      for (int iX = 0; iX < nDPhiBins; iX++) {
        h_jet_trk_dphi_gt4_lt6->SetBinContent (iX+1, h_jet_trk_dphi_gt4_lt6->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt4_lt6_counts[iX]));
        for (int iY = 0; iY < nDPhiBins; iY++)
          h2_jet_trk_dphi_gt4_lt6_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt4_lt6_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt4_lt6_counts[iX])*(jet_trk_dphi_gt4_lt6_counts[iY]));
      }
      for (int iX = 0; iX < nDPhiBins; iX++) {
        h_jet_trk_dphi_gt6_lt8->SetBinContent (iX+1, h_jet_trk_dphi_gt6_lt8->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt6_lt8_counts[iX]));
        for (int iY = 0; iY < nDPhiBins; iY++)
          h2_jet_trk_dphi_gt6_lt8_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt6_lt8_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt6_lt8_counts[iX])*(jet_trk_dphi_gt6_lt8_counts[iY]));
      }
      for (int iX = 0; iX < nDPhiBins; iX++) {
        h_jet_trk_dphi_gt8_lt10->SetBinContent (iX+1, h_jet_trk_dphi_gt8_lt10->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt8_lt10_counts[iX]));
        for (int iY = 0; iY < nDPhiBins; iY++)
          h2_jet_trk_dphi_gt8_lt10_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt8_lt10_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt8_lt10_counts[iX])*(jet_trk_dphi_gt8_lt10_counts[iY]));
      }
      for (int iX = 0; iX < nDPhiBins; iX++) {
        h_jet_trk_dphi_gt10_lt15->SetBinContent (iX+1, h_jet_trk_dphi_gt10_lt15->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt10_lt15_counts[iX]));
        for (int iY = 0; iY < nDPhiBins; iY++)
          h2_jet_trk_dphi_gt10_lt15_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt10_lt15_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt10_lt15_counts[iX])*(jet_trk_dphi_gt10_lt15_counts[iY]));
      }
      for (int iX = 0; iX < nDPhiBins; iX++) {
        h_jet_trk_dphi_gt15_lt20->SetBinContent (iX+1, h_jet_trk_dphi_gt15_lt20->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt15_lt20_counts[iX]));
        for (int iY = 0; iY < nDPhiBins; iY++)
          h2_jet_trk_dphi_gt15_lt20_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt15_lt20_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt15_lt20_counts[iX])*(jet_trk_dphi_gt15_lt20_counts[iY]));
      }
      for (int iX = 0; iX < nDPhiBins; iX++) {
        h_jet_trk_dphi_gt20_lt30->SetBinContent (iX+1, h_jet_trk_dphi_gt20_lt30->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_dphi_gt20_lt30_counts[iX]));
        for (int iY = 0; iY < nDPhiBins; iY++)
          h2_jet_trk_dphi_gt20_lt30_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_dphi_gt20_lt30_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_dphi_gt20_lt30_counts[iX])*(jet_trk_dphi_gt20_lt30_counts[iY]));
      }

      for (int iX = 0; iX < nPtChBins; iX++) {
        h_jet_trk_pt_ns->SetBinContent (iX+1, h_jet_trk_pt_ns->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_pt_ns_counts[iX]));
        for (int iY = 0; iY < nPtChBins; iY++)
          h2_jet_trk_pt_ns_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_pt_ns_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_pt_ns_counts[iX])*(jet_trk_pt_ns_counts[iY]));
      }
      for (int iX = 0; iX < nPtChBins; iX++) {
        h_jet_trk_pt_as->SetBinContent (iX+1, h_jet_trk_pt_as->GetBinContent (iX+1) + (ewgt*jwgt)*(jet_trk_pt_as_counts[iX]));
        for (int iY = 0; iY < nPtChBins; iY++)
          h2_jet_trk_pt_as_cov->SetBinContent (iX+1, iY+1, h2_jet_trk_pt_as_cov->GetBinContent (iX+1, iY+1) + (ewgt*jwgt)*(jet_trk_pt_as_counts[iX])*(jet_trk_pt_as_counts[iY]));
      }
    }
  }
  cout << "Finished event loop." << endl;


  SaferDelete (&h2_trk_eff);
  SaferDelete (&h2_trk_pur);


  h_evt_counts->Write ();
  h_jet_counts->Write ();

  h_jet_pt->Write ();
  h2_jet_pt_cov->Write ();
  h2_jet_eta_phi->Write ();

  h_jet_trk_dphi_gt0p5_lt1->Write ();
  h2_jet_trk_dphi_gt0p5_lt1_cov->Write ();
  h_jet_trk_dphi_gt1_lt1p5->Write ();
  h2_jet_trk_dphi_gt1_lt1p5_cov->Write ();
  h_jet_trk_dphi_gt1p5_lt2->Write ();
  h2_jet_trk_dphi_gt1p5_lt2_cov->Write ();
  h_jet_trk_dphi_gt2_lt4->Write ();
  h2_jet_trk_dphi_gt2_lt4_cov->Write ();
  h_jet_trk_dphi_gt4_lt6->Write ();
  h2_jet_trk_dphi_gt4_lt6_cov->Write ();
  h_jet_trk_dphi_gt6_lt8->Write ();
  h2_jet_trk_dphi_gt6_lt8_cov->Write ();
  h_jet_trk_dphi_gt8_lt10->Write ();
  h2_jet_trk_dphi_gt8_lt10_cov->Write ();
  h_jet_trk_dphi_gt10_lt15->Write ();
  h2_jet_trk_dphi_gt10_lt15_cov->Write ();
  h_jet_trk_dphi_gt15_lt20->Write ();
  h2_jet_trk_dphi_gt15_lt20_cov->Write ();
  h_jet_trk_dphi_gt20_lt30->Write ();
  h2_jet_trk_dphi_gt20_lt30_cov->Write ();

  h_jet_trk_pt_ns->Write ();
  h2_jet_trk_pt_ns_cov->Write ();
  h_jet_trk_pt_as->Write ();
  h2_jet_trk_pt_as_cov->Write ();


  outFile->Close ();
}




void RunCorrelator (const char* collSys, const char* tag, const char* outFileName, const char* jetsInFileName, const char* tracksInFileName = nullptr) {

  if      (TString (collSys) == ToTString (CollisionSystem::pPb16s5TeV))  { collisionSystem = CollisionSystem::pPb16s5TeV; }
  else if (TString (collSys) == ToTString (CollisionSystem::pp17))        { collisionSystem = CollisionSystem::pp17;       }
  else {
    std::cout << "In RunCorrelator.C: Invalid or unsupported collision system, exiting." << std::endl;
    return;
  }

  if (jet_min_pt == -2)
    jet_min_pt = (double) std::atof (std::getenv ("JET_MIN_PT"));
  if (jet_max_pt == -2)
    jet_max_pt = (double) std::atof (std::getenv ("JET_MAX_PT"));
  doMixing = (tracksInFileName != nullptr);

  TFile* jetsInFile = new TFile (jetsInFileName, "read");
  TFile* tracksInFile = (tracksInFileName != nullptr ? new TFile (tracksInFileName, "read") : nullptr);

  TTree* jetsInTree = (TTree*) jetsInFile->Get (IspPb () ? "pPbTree" : "ppTree");
  TTree* tracksInTree = (tracksInFileName != nullptr ? (TTree*) tracksInFile->Get (IspPb () ? "pPbTree" : "ppTree") : jetsInTree);

  Correlator (TString (tag), outFileName, jetsInTree, tracksInTree);

  jetsInFile->Close ();
  if (tracksInFile)
    tracksInFile->Close ();

  return;
}




int main (int argc, char** argv) {
  assert (argc >= 7);

  jet_min_pt = (double) std::atof (argv[1]);
  jet_max_pt = (double) std::atof (argv[2]);

  if (argc == 8) {
    std::cout << "Main will execute RunCorrelator (const char* collSys, const char* tag, const char* outFileName, const char* jetsInFileName, const char* tracksInFileName)" << std::endl;
    RunCorrelator (argv[3], argv[4], argv[5], argv[6], argv[7]);
  }

  else if (argc == 7) {
    std::cout << "Main will execute RunCorrelator (const char* collSys, const char* tag, const char* outFileName, const char* jetsInFileName)" << std::endl;
    RunCorrelator (argv[3], argv[4], argv[5], argv[6]);
  }

  else {
    std::cout << "Undefined behavior for " << argc << " arguments, exiting" << std::endl;
    return 1;
  }

  return 0;
}

#endif
