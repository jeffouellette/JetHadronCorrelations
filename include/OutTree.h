#ifndef __OutTree_h__
#define __OutTree_h__

#include "TreeVariables.h"

#include <TTree.h>

#include <vector>

using namespace std;

namespace JetHadronCorrelations {

// Here we define additional analysis-level variables unavailable from the AOD output,
// plus "out_" versions of arrays with slimmed entries meeting basic selections.
float event_weight = 0;

float ip = 0;
float eventPlane = 0;

float zdcEnergy = 0;
float vz = 0;

float fcal_et_Pb = 0;
float fcal_et_p = 0;
float q2x_A = 0;
float q2y_A = 0;
float q2x_C = 0;
float q2y_C = 0;
float q2x_Pb = 0;
float q2y_Pb = 0;
float q2x_p = 0;
float q2y_p = 0;
float zdc_calibE_Pb = 0;
float zdc_calibE_p = 0;

float cluster_sumGap_Pb = 0;
float cluster_sumGap_p = 0;
float cluster_edgeGap_Pb = 0;
float cluster_edgeGap_p = 0;
float sumGap_Pb = 0;
float sumGap_p = 0;
float edgeGap_Pb = 0;
float edgeGap_p = 0;

int   leading_jet = 0;
int   subleading_jet = 0;

float leading_jet_phi_transmin = 0;
float leading_jet_phi_transmax = 0;

int   out_truth_trk_n = 0;
float out_truth_trk_pt[max_truth_trk_n];
float out_truth_trk_eta[max_truth_trk_n];
float out_truth_trk_phi[max_truth_trk_n];
int   out_truth_trk_charge[max_truth_trk_n];

int   out_trk_n = 0;
float out_trk_pt[max_trk_n];
float out_trk_eta[max_trk_n];
float out_trk_phi[max_trk_n];
float out_trk_charge[max_trk_n];
float out_trk_d0[max_trk_n];
float out_trk_z0[max_trk_n];
bool  out_trk_TightPrimary[max_trk_n];
bool  out_trk_HITight[max_trk_n];
bool  out_trk_HILoose[max_trk_n];
bool  out_trk_truth_matched[max_trk_n];

int   out_akt4_hi_jet_n;
float out_akt4_hi_jet_pt[max_akt4_hi_jet_n];
float out_akt4_hi_jet_eta[max_akt4_hi_jet_n];
float out_akt4_hi_jet_phi[max_akt4_hi_jet_n];
float out_akt4_hi_jet_e[max_akt4_hi_jet_n];

int   out_akt4_truth_jet_n;
float out_akt4_truth_jet_pt[max_akt4_truth_jet_n];
float out_akt4_truth_jet_eta[max_akt4_truth_jet_n];
float out_akt4_truth_jet_phi[max_akt4_truth_jet_n];
float out_akt4_truth_jet_e[max_akt4_truth_jet_n];



struct OutTree {
  private:
  bool branchEventInfo = false;
  bool branchTruthTracks = false;
  bool branchTracks = false;
  bool branchTruthJets = false;
  bool branchJets = false;

  public:
  TTree* tree = nullptr;

  OutTree () {
    tree = nullptr;
  }
  OutTree (const char* name, TFile* file) {
    tree = new TTree (name, name);
    tree->SetDirectory (file);
  }
  ~OutTree () {
    if (tree != nullptr)
      SaferDelete (&tree);
  }

  void SetBranchEventInfo     (const bool _branchEventInfo = true)    { branchEventInfo = _branchEventInfo; }
  void SetBranchTruthTracks   (const bool _branchTruthTracks = true)  { branchTruthTracks = _branchTruthTracks; }
  void SetBranchTracks        (const bool _branchTracks = true)       { branchTracks = _branchTracks; }
  void SetBranchTruthJets     (const bool _branchTruthJets = true)    { branchTruthJets = _branchTruthJets; }
  void SetBranchJets          (const bool _branchJets = true)         { branchJets = _branchJets; }

  void SetBranches () {
    if (branchEventInfo) {
      tree->Branch ("event_number",  &event_number,  "event_number/i");
      tree->Branch ("lumi_block",    &lumi_block,    "lumi_block/i");
      tree->Branch ("run_number",    &run_number,    "run_number/i");

      tree->Branch ("event_weight",  &event_weight,  "event_weight/F");

      if (IsHijing ()) {
        tree->Branch ("impactParameter",  &ip,          "impactParameter/F");
        tree->Branch ("eventPlane",       &eventPlane,  "eventPlane/F");
      }

      tree->Branch ("vz",            &vz,            "vz/F");

      if (!Ispp ()) {
        tree->Branch ("fcal_et_Pb",    &fcal_et_Pb,    "fcal_et_Pb/F");
        tree->Branch ("fcal_et_p",     &fcal_et_p,     "fcal_et_p/F");
        tree->Branch ("q2x_Pb",        &q2x_Pb,        "q2x_Pb/F");
        tree->Branch ("q2y_Pb",        &q2y_Pb,        "q2y_Pb/F");
        tree->Branch ("q2x_p",         &q2x_p,         "q2x_p/F");
        tree->Branch ("q2y_p",         &q2y_p,         "q2y_p/F");
        tree->Branch ("zdc_calibE_Pb", &zdc_calibE_Pb, "zdc_calibE_Pb/F");
        tree->Branch ("zdc_calibE_p",  &zdc_calibE_p,  "zdc_calibE_p/F");

        tree->Branch ("cluster_sumGap_Pb",  &cluster_sumGap_Pb,   "cluster_sumGap_Pb/F");
        tree->Branch ("cluster_sumGap_p",   &cluster_sumGap_p,    "cluster_sumGap_p/F");
        tree->Branch ("cluster_edgeGap_Pb", &cluster_edgeGap_Pb,  "cluster_edgeGap_Pb/F");
        tree->Branch ("cluster_edgeGap_p",  &cluster_edgeGap_p,   "cluster_edgeGap_p/F");
        tree->Branch ("sumGap_Pb",          &sumGap_Pb,           "sumGap_Pb/F");
        tree->Branch ("sumGap_p",           &sumGap_p,            "sumGap_p/F");
        tree->Branch ("edgeGap_Pb",         &edgeGap_Pb,          "edgeGap_Pb/F");
        tree->Branch ("edgeGap_p",          &edgeGap_p,           "edgeGap_p/F");
      }
      else {
        tree->Branch ("fcal_et_A",  &fcalA_et,  "fcal_et_A/F");
        tree->Branch ("fcal_et_C",  &fcalC_et,  "fcal_et_C/F");
        tree->Branch ("q2x_A",      &q2x_A,     "q2x_A/F");
        tree->Branch ("q2y_A",      &q2y_A,     "q2y_A/F");
        tree->Branch ("q2x_C",      &q2x_C,     "q2x_C/F");
        tree->Branch ("q2y_C",      &q2y_C,     "q2y_C/F");

        tree->Branch ("cluster_sumGap_A",   &cluster_sumGap_A,  "cluster_sumGap_A/F");
        tree->Branch ("cluster_sumGap_C",   &cluster_sumGap_C,  "cluster_sumGap_C/F");
        tree->Branch ("cluster_edgeGap_A",  &cluster_edgeGap_A, "cluster_edgeGap_A/F");
        tree->Branch ("cluster_edgeGap_C",  &cluster_edgeGap_C, "cluster_edgeGap_C/F");
        tree->Branch ("sumGap_A",           &sumGap_A,          "sumGap_A/F");
        tree->Branch ("sumGap_C",           &sumGap_C,          "sumGap_C/F");
        tree->Branch ("edgeGap_A",          &edgeGap_A,         "edgeGap_A/F");
        tree->Branch ("edgeGap_C",          &edgeGap_C,         "edgeGap_C/F");
      }
    }

    if (!IsCollisions () && branchTruthTracks) {
      tree->Branch ("truth_trk_n",          &out_truth_trk_n,       "truth_trk_n/I");
      //tree->Branch ("truth_trk_pt",         &out_truth_trk_pt,      "truth_trk_pt[truth_trk_n]/F");
      //tree->Branch ("truth_trk_eta",        &out_truth_trk_eta,     "truth_trk_eta[truth_trk_n]/F");
      //tree->Branch ("truth_trk_phi",        &out_truth_trk_phi,     "truth_trk_phi[truth_trk_n]/F");
      //tree->Branch ("truth_trk_charge",     &out_truth_trk_charge,  "truth_trk_charge[truth_trk_n]/F");
    }
    if (branchTracks) {
      tree->Branch ("trk_n",                &out_trk_n,             "trk_n/I");
      //tree->Branch ("trk_pt",               &out_trk_pt,            "trk_pt[trk_n]/F");
      //tree->Branch ("trk_eta",              &out_trk_eta,           "trk_eta[trk_n]/F");
      //tree->Branch ("trk_phi",              &out_trk_phi,           "trk_phi[trk_n]/F");
      //tree->Branch ("trk_charge",           &out_trk_charge,        "trk_charge[trk_n]/F");
      //tree->Branch ("trk_d0",               &out_trk_d0,            "trk_d0[trk_n]/F");
      //tree->Branch ("trk_z0",               &out_trk_z0,            "trk_z0[trk_n]/F");
      //tree->Branch ("trk_TightPrimary",     &out_trk_TightPrimary,  "trk_TightPrimary[trk_n]/O");
      //tree->Branch ("trk_HITight",          &out_trk_HITight,       "trk_HITight[trk_n]/O");
      //tree->Branch ("trk_HILoose",          &out_trk_HILoose,       "trk_HILoose[trk_n]/O");
      //if (!IsCollisions ())
      //  tree->Branch ("trk_truth_matched",  &out_trk_truth_matched, "trk_truth_matched[trk_n]/O");
    }

    if (!IsCollisions () && branchTruthJets) {
      tree->Branch ("akt4_truth_jet_n",    &out_akt4_truth_jet_n,   "akt4_truth_jet_n/I");
      tree->Branch ("akt4_truth_jet_pt",   &out_akt4_truth_jet_pt,  "akt4_truth_jet_pt[akt4_truth_jet_n]/F");
      tree->Branch ("akt4_truth_jet_eta",  &out_akt4_truth_jet_eta, "akt4_truth_jet_eta[akt4_truth_jet_n]/F");
      tree->Branch ("akt4_truth_jet_phi",  &out_akt4_truth_jet_phi, "akt4_truth_jet_phi[akt4_truth_jet_n]/F");
      tree->Branch ("akt4_truth_jet_e",    &out_akt4_truth_jet_e,   "akt4_truth_jet_e[akt4_truth_jet_n]/F");
    }
    if (branchJets) {
      tree->Branch ("leading_jet",              &leading_jet,               "leading_jet/I");
      tree->Branch ("subleading_jet",           &subleading_jet,            "subleading_jet/I");
      tree->Branch ("akt4_hi_jet_n",            &out_akt4_hi_jet_n,         "akt4_hi_jet_n/I");
      tree->Branch ("akt4_hi_jet_pt",           &out_akt4_hi_jet_pt,        "akt4_hi_jet_pt[akt4_hi_jet_n]/F");
      tree->Branch ("akt4_hi_jet_eta",          &out_akt4_hi_jet_eta,       "akt4_hi_jet_eta[akt4_hi_jet_n]/F");
      tree->Branch ("akt4_hi_jet_phi",          &out_akt4_hi_jet_phi,       "akt4_hi_jet_phi[akt4_hi_jet_n]/F");
      tree->Branch ("akt4_hi_jet_e",            &out_akt4_hi_jet_e,         "akt4_hi_jet_e[akt4_hi_jet_n]/F");
      tree->Branch ("leading_jet_phi_transmin", &leading_jet_phi_transmin,  "leading_jet_phi_transmin/F");
      tree->Branch ("leading_jet_phi_transmax", &leading_jet_phi_transmax,  "leading_jet_phi_transmax/F");
    }
    return;
  }

  void Fill () {
    tree->Fill ();
    return;
  }
};

}

#endif
