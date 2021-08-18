#ifndef __TreeVariables_h__
#define __TreeVariables_h__

#include <TTree.h>
#include <TH1.h>

#include <vector>


float event_weight = 0;


// event info 
unsigned int run_number;
unsigned int lumi_block;
unsigned int event_number;
bool passes_toroid; 
bool isOOTPU;
bool BlayerDesyn;

// collision info
float actualInteractionsPerCrossing;
float averageInteractionsPerCrossing;

std::vector<float>* mcEventWeights;

// truth event info
const static int max_truth_event_n = 5;
int truth_event_n;
int nPart1[max_truth_event_n];
int nPart2[max_truth_event_n];
float impactParameter[max_truth_event_n];
int nColl[max_truth_event_n];
float sigmaInelasticNN[max_truth_event_n];
int nSpectatorNeutrons[max_truth_event_n];
int nSpectatorProtons[max_truth_event_n];
float eccentricity[max_truth_event_n];
float eventPlaneAngle[max_truth_event_n];

// Minimum bias trigger names
int minbias_trig_n;
const std::string* minbias_trig_name;

const static int minbias_trig_n_PbPb18 = 3;
const std::string minbias_trig_name_PbPb18 [minbias_trig_n_PbPb18] = {
  "HLT_noalg_cc_L1TE600.0ETA49",
  "HLT_mb_sptrk_L1ZDC_A_C_VTE50",
  "HLT_noalg_pc_L1TE50_VTE600.0ETA49"
};

const static int minbias_trig_n_pp17 = 1;
const std::string minbias_trig_name_pp17 [minbias_trig_n_pp17] = {
  "HLT_mb_sptrk"
};

const static int minbias_trig_n_pPb16 = 2;
const std::string minbias_trig_name_pPb16 [minbias_trig_n_pPb16] = {
  "HLT_mb_sptrk_L1MBTS_1",
  "HLT_mb_sptrk"
};

const static int minbias_trig_n_PbPb15 = 2;
const std::string minbias_trig_name_PbPb15 [minbias_trig_n_PbPb15] = {
  "HLT_noalg_mb_L1TE50",
  "HLT_mb_sptrk_ion_L1ZDC_A_C_VTE50"
};

const static int max_minbias_trig_n = std::max (minbias_trig_n_PbPb18, std::max (minbias_trig_n_pp17, std::max (minbias_trig_n_pPb16, minbias_trig_n_PbPb15)));


// Electron trigger names
int electron_trig_n;
const std::string* electron_trig_name;

const static int electron_trig_n_PbPb18 = 1;
const std::string electron_trig_name_PbPb18 [electron_trig_n_PbPb18] = {
  "HLT_e15_lhloose_ion_L1EM12"
};

const static int electron_trig_n_pp17 = 1;
const std::string electron_trig_name_pp17 [electron_trig_n_pp17] = {
  "HLT_e15_lhloose_L1EM12"
};

const static int electron_trig_n_pPb16 = 1;
const std::string electron_trig_name_pPb16 [electron_trig_n_pPb16] = {
  "HLT_e15_lhloose_nod0"
};

const static int electron_trig_n_PbPb15 = 1;
const std::string electron_trig_name_PbPb15 [electron_trig_n_PbPb15] = {
  "HLT_e15_loose_ion_L1EM12"
};

const static int max_electron_trig_n = std::max (electron_trig_n_PbPb18, std::max (electron_trig_n_pp17, std::max (electron_trig_n_pPb16, electron_trig_n_PbPb15)));


// Muon trigger names
int muon_trig_n;
const std::string* muon_trig_name;

const static int muon_trig_n_PbPb18 = 1;
const std::string muon_trig_name_PbPb18 [muon_trig_n_PbPb18] = {
  "HLT_mu14"
};

const static int muon_trig_n_pp17 = 1;
const std::string muon_trig_name_pp17 [muon_trig_n_pp17] = {
  "HLT_mu14"
};

const static int muon_trig_n_pPb16 = 1;
const std::string muon_trig_name_pPb16 [muon_trig_n_pPb16] = {
  "HLT_mu15"
};

const static int muon_trig_n_PbPb15 = 1;
const std::string muon_trig_name_PbPb15 [muon_trig_n_PbPb15] = {
  "HLT_mu14"
};

const static int max_muon_trig_n = std::max (muon_trig_n_PbPb18, std::max (muon_trig_n_pp17, std::max (muon_trig_n_pPb16, muon_trig_n_PbPb15)));


// Photon trigger names
int photon_trig_n;
const std::string* photon_trig_name;

const static int photon_trig_n_PbPb18 = 2;
const std::string photon_trig_name_PbPb18 [photon_trig_n_PbPb18] = {
  "HLT_g30_loose_ion",
  "HLT_g50_loose_ion"
};

const static int photon_trig_n_pp17 = 2;
const std::string photon_trig_name_pp17 [photon_trig_n_pp17] = {
  //"HLT_g15_loose_L1EM7",
  //"HLT_g20_loose_L1EM15",
  //"HLT_g25_loose_L1EM15",
  "HLT_g30_loose_L1EM15",
  "HLT_g35_loose_L1EM15"
};

const static int photon_trig_n_pPb16 = 2;
const std::string photon_trig_name_pPb16 [photon_trig_n_pPb16] = {
  //"HLT_g10_loose",
  //"HLT_g15_loose",
  //"HLT_g20_loose",
  //"HLT_g25_loose",
  "HLT_g30_loose",
  "HLT_g35_loose"
};

const static int photon_trig_n_PbPb15 = 1;
const std::string photon_trig_name_PbPb15 [photon_trig_n_PbPb15] = {
  "HLT_g20_loose_ion"
};

const static int max_photon_trig_n = std::max (photon_trig_n_PbPb18, std::max (photon_trig_n_pp17, std::max (photon_trig_n_pPb16, photon_trig_n_PbPb15)));


// Jet trigger names
int jet_trig_n;
const std::string* jet_trig_name;

const static int jet_trig_n_PbPb18 = 1;
const std::string jet_trig_name_PbPb18 [jet_trig_n_PbPb18] = {
  "HLT_j100_ion_L1J20"
};

const static int jet_trig_n_pp17 = 3;
const std::string jet_trig_name_pp17 [jet_trig_n_pp17] = {
  //"HLT_j110"
  "HLT_j50_L1J15", // to match 5 TeV p+Pb trigger
  "HLT_j75_L1J20",
  "HLT_j100_L1J20"
};

const static int jet_trig_n_pPb16s5TeV = 1;
const std::string jet_trig_name_pPb16s5TeV [jet_trig_n_pPb16s5TeV] = {
  "HLT_j50_ion_L1J10"
};

const static int jet_trig_n_pPb16 = 1;
const std::string jet_trig_name_pPb16 [jet_trig_n_pPb16] = {
  "HLT_j100_ion_L1J20"
};

const static int jet_trig_n_Pbp16 = 1;
const std::string jet_trig_name_Pbp16 [jet_trig_n_Pbp16] = {
  "HLT_j100_L1J20"
};

const static int jet_trig_n_PbPb15 = 1;
const std::string jet_trig_name_PbPb15 [jet_trig_n_PbPb15] = {
  "HLT_j75_ion_L1TE50"
};


// ZDC L1 trigger names
const static int zdc_L1_trig_n = 2;
const std::string zdc_L1_trig_name [zdc_L1_trig_n] = {
  "L1_ZDC_A",
  "L1_ZDC_C"
};


// Minimum bias trigger info
bool* minbias_trig_decision;
float* minbias_trig_prescale;

// Electron trigger info
bool* electron_trig_decision;
float* electron_trig_prescale;

// Muon trigger info
bool* muon_trig_decision;
float* muon_trig_prescale;

// Photon trigger info
bool* photon_trig_decision;
float* photon_trig_prescale;

// Jet trigger info
bool* jet_trig_decision;
float* jet_trig_prescale;

// ZDC L1 trigger info
bool* zdc_L1_trig_decision;
bool* zdc_L1_trig_tbp;
bool* zdc_L1_trig_tap;
bool* zdc_L1_trig_tav;
float* zdc_L1_trig_prescale;

// vertex info
const static int max_nvert = 30;
int nvert;
float vert_x[max_nvert];
float vert_y[max_nvert];
float vert_z[max_nvert];
int vert_ntrk[max_nvert];
int vert_type[max_nvert];
float vert_sumpt[max_nvert];

int nvert_matching;
float vert_x_matching[max_nvert];
float vert_y_matching[max_nvert];
float vert_z_matching[max_nvert];
int vert_ntrk_matching[max_nvert];
int vert_type_matching[max_nvert];
float vert_sumpt_matching[max_nvert];

// FCal Et info
float fcalA_et;
float fcalC_et;
float fcalA_et_Cos2;
float fcalC_et_Cos2;
float fcalA_et_Sin2;
float fcalC_et_Sin2;
float fcalA_et_Cos3;
float fcalC_et_Cos3;
float fcalA_et_Sin3;
float fcalC_et_Sin3;
float fcalA_et_Cos4;
float fcalC_et_Cos4;
float fcalA_et_Sin4;
float fcalC_et_Sin4;

float fcalA_et_matching;
float fcalC_et_matching;
float fcalA_et_Cos2_matching;
float fcalC_et_Cos2_matching;
float fcalA_et_Sin2_matching;
float fcalC_et_Sin2_matching;
float fcalA_et_Cos3_matching;
float fcalC_et_Cos3_matching;
float fcalA_et_Sin3_matching;
float fcalC_et_Sin3_matching;
float fcalA_et_Cos4_matching;
float fcalC_et_Cos4_matching;
float fcalA_et_Sin4_matching;
float fcalC_et_Sin4_matching;

// ZDC energies
float  ZdcRawEnergy_A;
float  ZdcRawEnergy_C;
float  ZdcCalibEnergy_A;
float  ZdcCalibEnergy_C;

float  ZdcRawEnergy_A_matching;
float  ZdcRawEnergy_C_matching;
float  ZdcCalibEnergy_A_matching;
float  ZdcCalibEnergy_C_matching;

// Sum of gaps and edge gaps for UPC background
float cluster_sumGap_A;
float cluster_sumGap_C;
float cluster_edgeGap_A;
float cluster_edgeGap_C;
float sumGap_A;
float sumGap_C;
float edgeGap_A;
float edgeGap_C;

//// Photon info
//const static int max_photon_n = 1000;
//int photon_n;
//float photon_pt_precalib[max_photon_n];
//float photon_pt[max_photon_n];
//float photon_eta[max_photon_n];
//float photon_etaBE[max_photon_n];
//float photon_phi[max_photon_n];
//bool photon_matched[max_photon_trig_n][max_photon_n];
//bool photon_tight[max_photon_n];
//bool photon_medium[max_photon_n];
//bool photon_loose[max_photon_n];
//unsigned int photon_isem[max_photon_n];
//int photon_convFlag[max_photon_n];
//float photon_Rconv[max_photon_n];
//float photon_etcone20[max_photon_n];
//float photon_etcone30[max_photon_n];
//float photon_etcone40[max_photon_n];
//float photon_topoetcone20[max_photon_n];
//float photon_topoetcone30[max_photon_n];
//float photon_topoetcone40[max_photon_n];
//float photon_Rhad1[max_photon_n];
//float photon_Rhad[max_photon_n];
//float photon_e277[max_photon_n];
//float photon_Reta[max_photon_n];
//float photon_Rphi[max_photon_n];
//float photon_weta1[max_photon_n];
//float photon_weta2[max_photon_n];
//float photon_wtots1[max_photon_n];
//float photon_f1[max_photon_n];
//float photon_f3[max_photon_n];
//float photon_fracs1[max_photon_n];
//float photon_DeltaE[max_photon_n];
//float photon_Eratio[max_photon_n];
//float photon_pt_sys[max_photon_n];
//float photon_eta_sys[max_photon_n];
//float photon_phi_sys[max_photon_n];

//// Truth photon info
//const static int max_truth_photon_n = 1000;
//int truth_photon_n;
//float truth_photon_pt[max_truth_photon_n];
//float truth_photon_eta[max_truth_photon_n];
//float truth_photon_phi[max_truth_photon_n];
//int truth_photon_barcode[max_truth_photon_n];

//// Electrons info 
//const static int max_electron_n = 1000;
//int electron_n;
//float electron_pt_precalib[max_electron_n];
//float electron_pt[max_electron_n];
//float electron_eta[max_electron_n];
//float electron_etaBE[max_electron_n];
//float electron_phi[max_electron_n];
//int electron_charge[max_electron_n];
//bool electron_lhtight[max_electron_n];
//bool electron_lhmedium[max_electron_n];
//bool electron_lhloose[max_electron_n];
//bool electron_lhmediuhi[max_electron_n];
//bool electron_lhloose_hi[max_electron_n];
//bool electron_matched[max_electron_trig_n][max_electron_n];
//float electron_etcone20[max_electron_n];
//float electron_etcone30[max_electron_n];
//float electron_etcone40[max_electron_n];
//float electron_topoetcone20[max_electron_n];
//float electron_topoetcone30[max_electron_n];
//float electron_topoetcone40[max_electron_n];
//int electron_ntrk[max_electron_n];
//float electron_id_track_pt[max_electron_n];
//float electron_id_track_eta[max_electron_n];
//float electron_id_track_phi[max_electron_n];
//float electron_id_track_charge[max_electron_n];
//float electron_id_track_d0[max_electron_n];
//float electron_id_track_d0sig[max_electron_n];
//float electron_id_track_z0[max_electron_n];
//float electron_id_track_z0sig[max_electron_n];
//float electron_id_track_theta[max_electron_n];
//float electron_id_track_vz[max_electron_n];
//bool electron_id_track_tightprimary[max_electron_n];
//bool electron_id_track_hiloose[max_electron_n];
//bool electron_id_track_hitight[max_electron_n];
//float electron_Rhad1[max_electron_n];
//float electron_Rhad[max_electron_n];
//float electron_e277[max_electron_n];
//float electron_Reta[max_electron_n];
//float electron_Rphi[max_electron_n];
//float electron_weta1[max_electron_n];
//float electron_weta2[max_electron_n];
//float electron_wtots1[max_electron_n];
//float electron_f1[max_electron_n];
//float electron_f3[max_electron_n];
//float electron_fracs1[max_electron_n];
//float electron_DeltaE[max_electron_n];
//float electron_Eratio[max_electron_n];
//float electron_pt_sys[max_electron_n];
//float electron_eta_sys[max_electron_n];
//float electron_phi_sys[max_electron_n];

//// Truth electrons info
//const static int max_truth_electron_n = 1000;
//int truth_electron_n;
//float truth_electron_pt[max_truth_electron_n];
//float truth_electron_eta[max_truth_electron_n];
//float truth_electron_phi[max_truth_electron_n];
//int truth_electron_charge[max_truth_electron_n];
//int truth_electron_barcode[max_truth_electron_n];

//// EGamma calorimeter cluster info, see https://ucatlas.github.io/RootCoreDocumentation/2.4.28/dc/d4b/CaloCluster__v1_8h_source.html
//const static int max_cluster_n = 250;
//int cluster_n;
//float cluster_pt[max_cluster_n];
//float cluster_et[max_cluster_n];
//float cluster_eta[max_cluster_n];
//float cluster_phi[max_cluster_n];
//float cluster_energyBE[max_cluster_n];
//float cluster_etaBE[max_cluster_n];
//float cluster_phiBE[max_cluster_n];
//float cluster_calE[max_cluster_n];
//float cluster_calEta[max_cluster_n];
//float cluster_calPhi[max_cluster_n];
//int cluster_size[max_cluster_n];
//int cluster_status[max_cluster_n];

//// Truth muons info
//const static int max_truth_muon_n = 1000;
//int truth_muon_n;
//float truth_muon_pt[max_truth_muon_n];
//float truth_muon_eta[max_truth_muon_n];
//float truth_muon_phi[max_truth_muon_n];
//int truth_muon_charge[max_truth_muon_n];
//int truth_muon_barcode[max_truth_muon_n];

//// Muons info
//const static int max_muon_n = 40;
//int muon_n;
//float muon_pt_precalib[max_muon_n];
//float muon_pt[max_muon_n];
//float muon_ms_pt_precalib[max_muon_n];
//float muon_ms_pt[max_muon_n];
//float muon_eta[max_muon_n];
//float muon_phi[max_muon_n];
//int muon_charge[max_muon_n];
//bool muon_tight[max_muon_n];
//bool muon_medium[max_muon_n];
//bool muon_loose[max_muon_n];
//bool muon_matched[max_muon_trig_n][max_muon_n];
//float muon_etcone20[max_muon_n];
//float muon_etcone30[max_muon_n];
//float muon_etcone40[max_muon_n];
//float muon_topoetcone20[max_muon_n];
//float muon_topoetcone30[max_muon_n];
//float muon_topoetcone40[max_muon_n];
//float muon_id_track_pt[max_muon_n];
//float muon_id_track_eta[max_muon_n];
//float muon_id_track_phi[max_muon_n];
//float muon_id_track_charge[max_muon_n];
//float muon_id_track_d0[max_muon_n];
//float muon_id_track_d0sig[max_muon_n];
//float muon_id_track_z0[max_muon_n];
//float muon_id_track_z0sig[max_muon_n];
//float muon_id_track_theta[max_muon_n];
//float muon_id_track_vz[max_muon_n];
//bool muon_id_track_tightprimary[max_muon_n];
//bool muon_id_track_hiloose[max_muon_n];
//bool muon_id_track_hitight[max_muon_n];
//float muon_ms_track_pt[max_muon_n];
//float muon_ms_track_eta[max_muon_n];
//float muon_ms_track_phi[max_muon_n];
//float muon_ms_track_charge[max_muon_n];
//float muon_ms_track_d0[max_muon_n];
//float muon_ms_track_d0sig[max_muon_n];
//float muon_ms_track_z0[max_muon_n];
//float muon_ms_track_z0sig[max_muon_n];
//float muon_ms_track_theta[max_muon_n];
//float muon_ms_track_vz[max_muon_n];
//std::vector <std::vector <double>> muon_pt_sys;
//std::vector <std::vector <double>> muon_eta_sys;
//std::vector <std::vector <double>> muon_phi_sys;

// Truth tracks info
const static int max_truth_trk_n = 10000;
int truth_trk_n;
float truth_trk_pt[max_truth_trk_n];
float truth_trk_eta[max_truth_trk_n];
float truth_trk_phi[max_truth_trk_n];
float truth_trk_charge[max_truth_trk_n];
int truth_trk_pdgid[max_truth_trk_n];
int truth_trk_barcode[max_truth_trk_n];
bool truth_trk_isHadron[max_truth_trk_n];

// Tracks info
const static int max_trk_n = 10000;
int trk_n;
float trk_pt[max_trk_n];
float trk_eta[max_trk_n];
float trk_phi[max_trk_n];
float trk_charge[max_trk_n];
bool trk_Loose[max_trk_n]; //!
bool trk_LoosePrimary[max_trk_n]; //!
bool trk_TightPrimary[max_trk_n]; //!
bool trk_HITight[max_trk_n]; //!
bool trk_HILoose[max_trk_n]; //!
float trk_d0[max_trk_n];
float trk_d0sig[max_trk_n];
float trk_z0[max_trk_n];
float trk_z0sig[max_trk_n];
float trk_theta[max_trk_n];
float trk_vz[max_trk_n];
char trk_nBLayerHits[max_trk_n]; //!
char trk_nBLayerSharedHits[max_trk_n]; //!
char trk_nPixelHits[max_trk_n]; //!
char trk_nPixelHoles[max_trk_n]; //!
char trk_nPixelSharedHits[max_trk_n]; //!
char trk_nPixelDeadSensors[max_trk_n]; //!
char trk_nSCTHits[max_trk_n]; //!
char trk_nSCTHoles[max_trk_n]; //!
char trk_nSCTDoubleHoles[max_trk_n]; //!
char trk_nSCTSharedHits[max_trk_n]; //!
char trk_nSCTDeadSensors[max_trk_n]; //!
float trk_pixeldEdx[max_trk_n]; //!
float trk_prob_truth[max_trk_n];
float trk_truth_pt[max_trk_n];
float trk_truth_eta[max_trk_n];
float trk_truth_phi[max_trk_n];
float trk_truth_charge[max_trk_n];
int trk_truth_type[max_trk_n];
int trk_truth_orig[max_trk_n];
int trk_truth_barcode[max_trk_n];
int trk_truth_pdgid[max_trk_n];
float trk_truth_vz[max_trk_n];
int trk_truth_nIn[max_trk_n];
bool trk_truth_isHadron[max_trk_n];

// list of track working points
const std::vector <bool*> trackWPs = {trk_TightPrimary, trk_HITight, trk_HILoose};
const std::vector <std::string> trackWPStrs = {"TightPrimary", "HITight", "HILoose"};

// Anti-kT R=0.2 Truth Jets info
const static int max_akt2_truth_jet_n = 40;
int akt2_truth_jet_n;
float akt2_truth_jet_pt[max_akt2_truth_jet_n];
float akt2_truth_jet_eta[max_akt2_truth_jet_n];
float akt2_truth_jet_phi[max_akt2_truth_jet_n];
float akt2_truth_jet_e[max_akt2_truth_jet_n];

// Anti-kT R=0.4 Truth Jets info
const static int max_akt4_truth_jet_n = 40;
int akt4_truth_jet_n;
float akt4_truth_jet_pt[max_akt4_truth_jet_n];
float akt4_truth_jet_eta[max_akt4_truth_jet_n];
float akt4_truth_jet_phi[max_akt4_truth_jet_n];
float akt4_truth_jet_e[max_akt4_truth_jet_n];

// Anti-kT R=0.4 EMTopo Jets info
const static int max_akt4_emtopo_jet_n = 40;
int akt4_emtopo_jet_n;
float akt4_emtopo_jet_pt[max_akt4_emtopo_jet_n];
float akt4_emtopo_jet_eta[max_akt4_emtopo_jet_n];
float akt4_emtopo_jet_phi[max_akt4_emtopo_jet_n];
float akt4_emtopo_jet_e[max_akt4_emtopo_jet_n];
bool akt4_emtopo_jet_LooseBad[max_akt4_emtopo_jet_n];

const static int nJESSys = 21;
const static int nJERSys = 10;

// Anti-kT R=0.2 HI Jets info
const static int max_akt2_hi_jet_n = 40;
int akt2_hi_jet_n;
float akt2_hi_jet_phi[max_akt2_hi_jet_n];
float akt2_hi_jet_pt_precalib[max_akt2_hi_jet_n];
float akt2_hi_jet_eta_precalib[max_akt2_hi_jet_n];
float akt2_hi_jet_e_precalib[max_akt2_hi_jet_n];
float akt2_hi_jet_pt_etajes[max_akt2_hi_jet_n];
float akt2_hi_jet_eta_etajes[max_akt2_hi_jet_n];
float akt2_hi_jet_e_etajes[max_akt2_hi_jet_n];
float akt2_hi_jet_pt_xcalib[max_akt2_hi_jet_n];
float akt2_hi_jet_eta_xcalib[max_akt2_hi_jet_n];
float akt2_hi_jet_e_xcalib[max_akt2_hi_jet_n];
float akt2_hi_jet_sub_et[max_akt2_hi_jet_n];
float akt2_hi_jet_sub_e[max_akt2_hi_jet_n];
bool akt2_hi_jet_LooseBad[max_akt2_hi_jet_n];
float akt2_hi_jet_pt_sys_JES_ALL[nJESSys][max_akt2_hi_jet_n];
float akt2_hi_jet_pt_sys_JER_ALL[nJERSys][max_akt2_hi_jet_n];
float akt2_hi_jet_timing[max_akt2_hi_jet_n];

// Anti-kT R=0.4 HI Jets info
const static int max_akt4_hi_jet_n = 40;
int akt4_hi_jet_n;
float akt4_hi_jet_phi[max_akt4_hi_jet_n];
float akt4_hi_jet_pt_precalib[max_akt4_hi_jet_n];
float akt4_hi_jet_eta_precalib[max_akt4_hi_jet_n];
float akt4_hi_jet_e_precalib[max_akt4_hi_jet_n];
float akt4_hi_jet_pt_etajes[max_akt4_hi_jet_n];
float akt4_hi_jet_eta_etajes[max_akt4_hi_jet_n];
float akt4_hi_jet_e_etajes[max_akt4_hi_jet_n];
float akt4_hi_jet_pt_xcalib[max_akt4_hi_jet_n];
float akt4_hi_jet_eta_xcalib[max_akt4_hi_jet_n];
float akt4_hi_jet_e_xcalib[max_akt4_hi_jet_n];
float akt4_hi_jet_sub_et[max_akt4_hi_jet_n];
float akt4_hi_jet_sub_e[max_akt4_hi_jet_n];
bool akt4_hi_jet_LooseBad[max_akt4_hi_jet_n];
float akt4_hi_jet_pt_sys_JES_ALL[nJESSys][max_akt4_hi_jet_n];
float akt4_hi_jet_pt_sys_JER_ALL[nJERSys][max_akt4_hi_jet_n];
float akt4_hi_jet_timing[max_akt4_hi_jet_n];


#endif
