#ifndef __JetHadronCorrelatorPlotPtCh_C__
#define __JetHadronCorrelatorPlotPtCh_C__

#include "Params.h"
#include "OutTree.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"

#include <ArrayTemplates.h>
#include <Utilities.h>
#include <MyStyle.h>

#include <TColor.h>
#include <TLine.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>

#include <iostream>
#include <math.h>

using namespace JetHadronCorrelations;

TColor* tcolor = new TColor ();
const Color_t myBlue = (Color_t) tcolor->GetColor (45, 64, 245);
const Color_t myPurple = (Color_t) tcolor->GetColor (130,  10, 130);
const Color_t myRed = (Color_t) tcolor->GetColor (255,  12,  73);

const Color_t systColors[10] = {kRed+1, kAzure-2, kGreen+2, kViolet-3, kMagenta, kCyan+1, kOrange-3, kGreen-7, kBlue+1, kPink-5};

TLine* l = new TLine ();
TLatex* tl = new TLatex ();

const int nVar = 4;
std::vector <TString> variations = {"JetES5PercUpVar", "JetES5PercDownVar", "JetES2PercUpVar", "JetES2PercDownVar"};

std::map <TString, MyStyle> varStyles = {
  {"JetES5PercUpVar",   MyStyle (kPink+5, 3)},
  {"JetES5PercDownVar", MyStyle (kPink+5, 2)},
  {"JetES2PercUpVar",   MyStyle (kAzure+2, 3)},
  {"JetES2PercDownVar", MyStyle (kAzure+2, 2)}
};
std::map <TString, TString> varFullNames = {
  {"JetES5PercUpVar",   "JES 5\% up"},
  {"JetES5PercDownVar", "JES 5\% down"},
  {"JetES2PercUpVar",   "JES 2\% up"},
  {"JetES2PercDownVar", "JES 2\% down"}
};


void PlotPtCh (const char* tag, const char* inFileTag) {

  SetupDirectories ("Data");

  TFile* inFile = nullptr;

  TH1D** h_evt_counts = new TH1D*[2];
  TH1D** h_jet_counts = new TH1D*[2];
  TH1D** h_ljet_counts = new TH1D*[2];
  TH1D** h_sljet_counts = new TH1D*[2];
  TH1D** h_evt_counts_bkg = new TH1D*[2];
  TH1D** h_jet_counts_bkg = new TH1D*[2];
  TH1D** h_ljet_counts_bkg = new TH1D*[2];
  TH1D** h_sljet_counts_bkg = new TH1D*[2];

  TH1D** h_jet_trk_pt_ns = new TH1D*[2];
  TH1D** h_ljet_trk_pt_ns = new TH1D*[2];
  TH1D** h_sljet_trk_pt_ns = new TH1D*[2];
  TH1D** h_jet_trk_pt_as = new TH1D*[2];
  TH1D** h_ljet_trk_pt_as = new TH1D*[2];
  TH1D** h_sljet_trk_pt_as = new TH1D*[2];
  TH1D** h_jet_trk_pt_ns_bkg = new TH1D*[2];
  TH1D** h_ljet_trk_pt_ns_bkg = new TH1D*[2];
  TH1D** h_sljet_trk_pt_ns_bkg = new TH1D*[2];
  TH1D** h_jet_trk_pt_as_bkg = new TH1D*[2];
  TH1D** h_ljet_trk_pt_as_bkg = new TH1D*[2];
  TH1D** h_sljet_trk_pt_as_bkg = new TH1D*[2];

  TH1D** h_jet_trk_pt_ns_sig = new TH1D*[2];
  TH1D** h_ljet_trk_pt_ns_sig = new TH1D*[2];
  TH1D** h_sljet_trk_pt_ns_sig = new TH1D*[2];
  TH1D** h_jet_trk_pt_as_sig = new TH1D*[2];
  TH1D** h_ljet_trk_pt_as_sig = new TH1D*[2];
  TH1D** h_sljet_trk_pt_as_sig = new TH1D*[2];

  TH1D* h_jet_trk_pt_ns_iaa = nullptr;
  TH1D* h_ljet_trk_pt_ns_iaa = nullptr;
  TH1D* h_sljet_trk_pt_ns_iaa = nullptr;
  TH1D* h_jet_trk_pt_as_iaa = nullptr;
  TH1D* h_ljet_trk_pt_as_iaa = nullptr;
  TH1D* h_sljet_trk_pt_as_iaa = nullptr;

  TGAE** g_jet_trk_pt_ns_syst = new TGAE*[2];
  TGAE** g_ljet_trk_pt_ns_syst = new TGAE*[2];
  TGAE** g_sljet_trk_pt_ns_syst = new TGAE*[2];
  TGAE** g_jet_trk_pt_as_syst = new TGAE*[2];
  TGAE** g_ljet_trk_pt_as_syst = new TGAE*[2];
  TGAE** g_sljet_trk_pt_as_syst = new TGAE*[2];
  TGAE** g_jet_trk_pt_ns_bkg_syst = new TGAE*[2];
  TGAE** g_ljet_trk_pt_ns_bkg_syst = new TGAE*[2];
  TGAE** g_sljet_trk_pt_ns_bkg_syst = new TGAE*[2];
  TGAE** g_jet_trk_pt_as_bkg_syst = new TGAE*[2];
  TGAE** g_ljet_trk_pt_as_bkg_syst = new TGAE*[2];
  TGAE** g_sljet_trk_pt_as_bkg_syst = new TGAE*[2];

  TGAE** g_jet_trk_pt_ns_sig_syst = new TGAE*[2];
  TGAE** g_ljet_trk_pt_ns_sig_syst = new TGAE*[2];
  TGAE** g_sljet_trk_pt_ns_sig_syst = new TGAE*[2];
  TGAE** g_jet_trk_pt_as_sig_syst = new TGAE*[2];
  TGAE** g_ljet_trk_pt_as_sig_syst = new TGAE*[2];
  TGAE** g_sljet_trk_pt_as_sig_syst = new TGAE*[2];

  TGAE* g_jet_trk_pt_ns_iaa_syst = nullptr;
  TGAE* g_ljet_trk_pt_ns_iaa_syst = nullptr;
  TGAE* g_sljet_trk_pt_ns_iaa_syst = nullptr;
  TGAE* g_jet_trk_pt_as_iaa_syst = nullptr;
  TGAE* g_ljet_trk_pt_as_iaa_syst = nullptr;
  TGAE* g_sljet_trk_pt_as_iaa_syst = nullptr;


  TH1D*** h_jet_trk_pt_ns_syst = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_ljet_trk_pt_ns_syst = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_sljet_trk_pt_ns_syst = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_jet_trk_pt_as_syst = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_ljet_trk_pt_as_syst = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_sljet_trk_pt_as_syst = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_jet_trk_pt_ns_bkg_syst = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_ljet_trk_pt_ns_bkg_syst = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_sljet_trk_pt_ns_bkg_syst = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_jet_trk_pt_as_bkg_syst = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_ljet_trk_pt_as_bkg_syst = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_sljet_trk_pt_as_bkg_syst = Get2DArray <TH1D*> (2, nVar);

  TH1D*** h_jet_trk_pt_ns_sig_syst = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_ljet_trk_pt_ns_sig_syst = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_sljet_trk_pt_ns_sig_syst = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_jet_trk_pt_as_sig_syst = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_ljet_trk_pt_as_sig_syst = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_sljet_trk_pt_as_sig_syst = Get2DArray <TH1D*> (2, nVar);

  TH1D** h_jet_trk_pt_ns_iaa_syst = Get1DArray <TH1D*> (nVar);
  TH1D** h_ljet_trk_pt_ns_iaa_syst = Get1DArray <TH1D*> (nVar);
  TH1D** h_sljet_trk_pt_ns_iaa_syst = Get1DArray <TH1D*> (nVar);
  TH1D** h_jet_trk_pt_as_iaa_syst = Get1DArray <TH1D*> (nVar);
  TH1D** h_ljet_trk_pt_as_iaa_syst = Get1DArray <TH1D*> (nVar);
  TH1D** h_sljet_trk_pt_as_iaa_syst = Get1DArray <TH1D*> (nVar);


  {
    std::cout << Form ("Reading ./rootFiles/Results/PlotPtCh_%s.root", inFileTag) << std::endl;
    inFile = new TFile (Form ("./rootFiles/Results/PlotPtCh_%s.root", inFileTag), "read");

    h_evt_counts[0] = (TH1D*) inFile->Get ("h_evt_counts_ref");
    h_jet_counts[0] = (TH1D*) inFile->Get ("h_jet_counts_ref");
    h_ljet_counts[0] = (TH1D*) inFile->Get ("h_ljet_counts_ref");
    h_sljet_counts[0] = (TH1D*) inFile->Get ("h_sljet_counts_ref");
    h_evt_counts[1] = (TH1D*) inFile->Get ("h_evt_counts");
    h_jet_counts[1] = (TH1D*) inFile->Get ("h_jet_counts");
    h_ljet_counts[1] = (TH1D*) inFile->Get ("h_ljet_counts");
    h_sljet_counts[1] = (TH1D*) inFile->Get ("h_sljet_counts");
    h_evt_counts_bkg[1] = (TH1D*) inFile->Get ("h_evt_counts_bkg");
    h_jet_counts_bkg[1] = (TH1D*) inFile->Get ("h_jet_counts_bkg");
    h_ljet_counts_bkg[1] = (TH1D*) inFile->Get ("h_ljet_counts_bkg");
    h_sljet_counts_bkg[1] = (TH1D*) inFile->Get ("h_sljet_counts_bkg");

    h_jet_trk_pt_ns[0] = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_ref");
    h_ljet_trk_pt_ns[0] = (TH1D*) inFile->Get ("h_ljet_trk_pt_ns_ref");
    h_sljet_trk_pt_ns[0] = (TH1D*) inFile->Get ("h_sljet_trk_pt_ns_ref");
    h_jet_trk_pt_as[0] = (TH1D*) inFile->Get ("h_jet_trk_pt_as_ref");
    h_ljet_trk_pt_as[0] = (TH1D*) inFile->Get ("h_ljet_trk_pt_as_ref");
    h_sljet_trk_pt_as[0] = (TH1D*) inFile->Get ("h_sljet_trk_pt_as_ref");
    h_jet_trk_pt_ns[1] = (TH1D*) inFile->Get ("h_jet_trk_pt_ns");
    h_ljet_trk_pt_ns[1] = (TH1D*) inFile->Get ("h_ljet_trk_pt_ns");
    h_sljet_trk_pt_ns[1] = (TH1D*) inFile->Get ("h_sljet_trk_pt_ns");
    h_jet_trk_pt_as[1] = (TH1D*) inFile->Get ("h_jet_trk_pt_as");
    h_ljet_trk_pt_as[1] = (TH1D*) inFile->Get ("h_ljet_trk_pt_as");
    h_sljet_trk_pt_as[1] = (TH1D*) inFile->Get ("h_sljet_trk_pt_as");
    h_jet_trk_pt_ns_bkg[1] = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_bkg");
    h_ljet_trk_pt_ns_bkg[1] = (TH1D*) inFile->Get ("h_ljet_trk_pt_ns_bkg");
    h_sljet_trk_pt_ns_bkg[1] = (TH1D*) inFile->Get ("h_sljet_trk_pt_ns_bkg");
    h_jet_trk_pt_as_bkg[1] = (TH1D*) inFile->Get ("h_jet_trk_pt_as_bkg");
    h_ljet_trk_pt_as_bkg[1] = (TH1D*) inFile->Get ("h_ljet_trk_pt_as_bkg");
    h_sljet_trk_pt_as_bkg[1] = (TH1D*) inFile->Get ("h_sljet_trk_pt_as_bkg");

    h_jet_trk_pt_ns_sig[0] = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_ref_sig");
    h_ljet_trk_pt_ns_sig[0] = (TH1D*) inFile->Get ("h_ljet_trk_pt_ns_ref_sig");
    h_sljet_trk_pt_ns_sig[0] = (TH1D*) inFile->Get ("h_sljet_trk_pt_ns_ref_sig");
    h_jet_trk_pt_as_sig[0] = (TH1D*) inFile->Get ("h_jet_trk_pt_as_ref_sig");
    h_ljet_trk_pt_as_sig[0] = (TH1D*) inFile->Get ("h_ljet_trk_pt_as_ref_sig");
    h_sljet_trk_pt_as_sig[0] = (TH1D*) inFile->Get ("h_sljet_trk_pt_as_ref_sig");
    h_jet_trk_pt_ns_sig[1] = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_sig");
    h_ljet_trk_pt_ns_sig[1] = (TH1D*) inFile->Get ("h_ljet_trk_pt_ns_sig");
    h_sljet_trk_pt_ns_sig[1] = (TH1D*) inFile->Get ("h_sljet_trk_pt_ns_sig");
    h_jet_trk_pt_as_sig[1] = (TH1D*) inFile->Get ("h_jet_trk_pt_as_sig");
    h_ljet_trk_pt_as_sig[1] = (TH1D*) inFile->Get ("h_ljet_trk_pt_as_sig");
    h_sljet_trk_pt_as_sig[1] = (TH1D*) inFile->Get ("h_sljet_trk_pt_as_sig");

    h_jet_trk_pt_ns_iaa = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_iaa");
    h_ljet_trk_pt_ns_iaa = (TH1D*) inFile->Get ("h_ljet_trk_pt_ns_iaa");
    h_sljet_trk_pt_ns_iaa = (TH1D*) inFile->Get ("h_sljet_trk_pt_ns_iaa");
    h_jet_trk_pt_as_iaa = (TH1D*) inFile->Get ("h_jet_trk_pt_as_iaa");
    h_ljet_trk_pt_as_iaa = (TH1D*) inFile->Get ("h_ljet_trk_pt_as_iaa");
    h_sljet_trk_pt_as_iaa = (TH1D*) inFile->Get ("h_sljet_trk_pt_as_iaa");

    g_jet_trk_pt_ns_syst[0] = (TGAE*) inFile->Get ("g_jet_trk_pt_ns_ref_syst");
    g_ljet_trk_pt_ns_syst[0] = (TGAE*) inFile->Get ("g_ljet_trk_pt_ns_ref_syst");
    g_sljet_trk_pt_ns_syst[0] = (TGAE*) inFile->Get ("g_sljet_trk_pt_ns_ref_syst");
    g_jet_trk_pt_as_syst[0] = (TGAE*) inFile->Get ("g_jet_trk_pt_as_ref_syst");
    g_ljet_trk_pt_as_syst[0] = (TGAE*) inFile->Get ("g_ljet_trk_pt_as_ref_syst");
    g_sljet_trk_pt_as_syst[0] = (TGAE*) inFile->Get ("g_sljet_trk_pt_as_ref_syst");
    g_jet_trk_pt_ns_syst[1] = (TGAE*) inFile->Get ("g_jet_trk_pt_ns_syst");
    g_ljet_trk_pt_ns_syst[1] = (TGAE*) inFile->Get ("g_ljet_trk_pt_ns_syst");
    g_sljet_trk_pt_ns_syst[1] = (TGAE*) inFile->Get ("g_sljet_trk_pt_ns_syst");
    g_jet_trk_pt_as_syst[1] = (TGAE*) inFile->Get ("g_jet_trk_pt_as_syst");
    g_ljet_trk_pt_as_syst[1] = (TGAE*) inFile->Get ("g_ljet_trk_pt_as_syst");
    g_sljet_trk_pt_as_syst[1] = (TGAE*) inFile->Get ("g_sljet_trk_pt_as_syst");
    g_jet_trk_pt_ns_bkg_syst[1] = (TGAE*) inFile->Get ("g_jet_trk_pt_ns_bkg_syst");
    g_ljet_trk_pt_ns_bkg_syst[1] = (TGAE*) inFile->Get ("g_ljet_trk_pt_ns_bkg_syst");
    g_sljet_trk_pt_ns_bkg_syst[1] = (TGAE*) inFile->Get ("g_sljet_trk_pt_ns_bkg_syst");
    g_jet_trk_pt_as_bkg_syst[1] = (TGAE*) inFile->Get ("g_jet_trk_pt_as_bkg_syst");
    g_ljet_trk_pt_as_bkg_syst[1] = (TGAE*) inFile->Get ("g_ljet_trk_pt_as_bkg_syst");
    g_sljet_trk_pt_as_bkg_syst[1] = (TGAE*) inFile->Get ("g_sljet_trk_pt_as_bkg_syst");

    g_jet_trk_pt_ns_sig_syst[0] = (TGAE*) inFile->Get ("g_jet_trk_pt_ns_ref_sig_syst");
    g_ljet_trk_pt_ns_sig_syst[0] = (TGAE*) inFile->Get ("g_ljet_trk_pt_ns_ref_sig_syst");
    g_sljet_trk_pt_ns_sig_syst[0] = (TGAE*) inFile->Get ("g_sljet_trk_pt_ns_ref_sig_syst");
    g_jet_trk_pt_as_sig_syst[0] = (TGAE*) inFile->Get ("g_jet_trk_pt_as_ref_sig_syst");
    g_ljet_trk_pt_as_sig_syst[0] = (TGAE*) inFile->Get ("g_ljet_trk_pt_as_ref_sig_syst");
    g_sljet_trk_pt_as_sig_syst[0] = (TGAE*) inFile->Get ("g_sljet_trk_pt_as_ref_sig_syst");
    g_jet_trk_pt_ns_sig_syst[1] = (TGAE*) inFile->Get ("g_jet_trk_pt_ns_sig_syst");
    g_ljet_trk_pt_ns_sig_syst[1] = (TGAE*) inFile->Get ("g_ljet_trk_pt_ns_sig_syst");
    g_sljet_trk_pt_ns_sig_syst[1] = (TGAE*) inFile->Get ("g_sljet_trk_pt_ns_sig_syst");
    g_jet_trk_pt_as_sig_syst[1] = (TGAE*) inFile->Get ("g_jet_trk_pt_as_sig_syst");
    g_ljet_trk_pt_as_sig_syst[1] = (TGAE*) inFile->Get ("g_ljet_trk_pt_as_sig_syst");
    g_sljet_trk_pt_as_sig_syst[1] = (TGAE*) inFile->Get ("g_sljet_trk_pt_as_sig_syst");

    g_jet_trk_pt_ns_iaa_syst = (TGAE*) inFile->Get ("g_jet_trk_pt_ns_iaa_syst");
    g_ljet_trk_pt_ns_iaa_syst = (TGAE*) inFile->Get ("g_ljet_trk_pt_ns_iaa_syst");
    g_sljet_trk_pt_ns_iaa_syst = (TGAE*) inFile->Get ("g_sljet_trk_pt_ns_iaa_syst");
    g_jet_trk_pt_as_iaa_syst = (TGAE*) inFile->Get ("g_jet_trk_pt_as_iaa_syst");
    g_ljet_trk_pt_as_iaa_syst = (TGAE*) inFile->Get ("g_ljet_trk_pt_as_iaa_syst");
    g_sljet_trk_pt_as_iaa_syst = (TGAE*) inFile->Get ("g_sljet_trk_pt_as_iaa_syst");

    for (int iVar = 0; iVar < nVar; iVar++) {
      h_jet_trk_pt_ns_syst[0][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_ref_%s", variations[iVar].Data ()));
      h_ljet_trk_pt_ns_syst[0][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_pt_ns_ref_%s", variations[iVar].Data ()));
      h_sljet_trk_pt_ns_syst[0][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_pt_ns_ref_%s", variations[iVar].Data ()));
      h_jet_trk_pt_as_syst[0][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_ref_%s", variations[iVar].Data ()));
      h_ljet_trk_pt_as_syst[0][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_pt_as_ref_%s", variations[iVar].Data ()));
      h_sljet_trk_pt_as_syst[0][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_pt_as_ref_%s", variations[iVar].Data ()));
      h_jet_trk_pt_ns_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_%s", variations[iVar].Data ()));
      h_ljet_trk_pt_ns_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_pt_ns_%s", variations[iVar].Data ()));
      h_sljet_trk_pt_ns_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_pt_ns_%s", variations[iVar].Data ()));
      h_jet_trk_pt_as_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_%s", variations[iVar].Data ()));
      h_ljet_trk_pt_as_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_pt_as_%s", variations[iVar].Data ()));
      h_sljet_trk_pt_as_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_pt_as_%s", variations[iVar].Data ()));
      h_jet_trk_pt_ns_bkg_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_bkg_%s", variations[iVar].Data ()));
      h_ljet_trk_pt_ns_bkg_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_pt_ns_bkg_%s", variations[iVar].Data ()));
      h_sljet_trk_pt_ns_bkg_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_pt_ns_bkg_%s", variations[iVar].Data ()));
      h_jet_trk_pt_as_bkg_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_bkg_%s", variations[iVar].Data ()));
      h_ljet_trk_pt_as_bkg_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_pt_as_bkg_%s", variations[iVar].Data ()));
      h_sljet_trk_pt_as_bkg_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_pt_as_bkg_%s", variations[iVar].Data ()));

      h_jet_trk_pt_ns_sig_syst[0][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_ref_sig_%s", variations[iVar].Data ()));
      h_ljet_trk_pt_ns_sig_syst[0][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_pt_ns_ref_sig_%s", variations[iVar].Data ()));
      h_sljet_trk_pt_ns_sig_syst[0][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_pt_ns_ref_sig_%s", variations[iVar].Data ()));
      h_jet_trk_pt_as_sig_syst[0][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_ref_sig_%s", variations[iVar].Data ()));
      h_ljet_trk_pt_as_sig_syst[0][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_pt_as_ref_sig_%s", variations[iVar].Data ()));
      h_sljet_trk_pt_as_sig_syst[0][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_pt_as_ref_sig_%s", variations[iVar].Data ()));
      h_jet_trk_pt_ns_sig_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_sig_%s", variations[iVar].Data ()));
      h_ljet_trk_pt_ns_sig_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_pt_ns_sig_%s", variations[iVar].Data ()));
      h_sljet_trk_pt_ns_sig_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_pt_ns_sig_%s", variations[iVar].Data ()));
      h_jet_trk_pt_as_sig_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_sig_%s", variations[iVar].Data ()));
      h_ljet_trk_pt_as_sig_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_pt_as_sig_%s", variations[iVar].Data ()));
      h_sljet_trk_pt_as_sig_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_pt_as_sig_%s", variations[iVar].Data ()));

      h_jet_trk_pt_ns_iaa_syst[iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_iaa_%s", variations[iVar].Data ()));
      h_ljet_trk_pt_ns_iaa_syst[iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_pt_ns_iaa_%s", variations[iVar].Data ()));
      h_sljet_trk_pt_ns_iaa_syst[iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_pt_ns_iaa_%s", variations[iVar].Data ()));
      h_jet_trk_pt_as_iaa_syst[iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_iaa_%s", variations[iVar].Data ()));
      h_ljet_trk_pt_as_iaa_syst[iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_pt_as_iaa_%s", variations[iVar].Data ()));
      h_sljet_trk_pt_as_iaa_syst[iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_pt_as_iaa_%s", variations[iVar].Data ()));
    }
  }



  l->SetLineStyle (7);
  l->SetLineWidth (1);
  l->SetLineColor (kBlack);



  {
    const char* canvasName = "c_jet_trk_pt_as";
    TCanvas* c = new TCanvas (canvasName, "", 1200, 1120);
    c->cd ();
    const double fuPad = 480./1120.;
    const double fdPad = 320./1120.;
    const double fcPad = 1.0 - fuPad - fdPad;
    const double fxPad = 0.42;
    TPad* ulPad = new TPad (Form ("%s_ulPad", canvasName), "", 0.0, 1.0-fuPad, fxPad+0.1, 1.0);
    TPad* clPad = new TPad (Form ("%s_clPad", canvasName), "", 0.0, fdPad, fxPad+0.1, 1.0-fuPad);
    TPad* dlPad = new TPad (Form ("%s_dlPad", canvasName), "", 0.0, 0.0, fxPad+0.1, fdPad);
    TPad* urPad = new TPad (Form ("%s_urPad", canvasName), "", fxPad+0.1, 1.0-fuPad, 1.0, 1.0);
    TPad* crPad = new TPad (Form ("%s_crPad", canvasName), "", fxPad+0.1, fdPad, 1.0, 1.0-fuPad);
    TPad* drPad = new TPad (Form ("%s_drPad", canvasName), "", fxPad+0.1, 0.0, 1.0, fdPad);

    ulPad->SetTopMargin (0.14);
    urPad->SetTopMargin (0.14);
    ulPad->SetBottomMargin (0);
    urPad->SetBottomMargin (0);
    clPad->SetTopMargin (0);
    crPad->SetTopMargin (0);
    clPad->SetBottomMargin (0);
    crPad->SetBottomMargin (0);
    dlPad->SetTopMargin (0);
    drPad->SetTopMargin (0);
    dlPad->SetBottomMargin (0.25);
    drPad->SetBottomMargin (0.25);

    ulPad->SetLeftMargin (0.1 / 0.52);
    clPad->SetLeftMargin (0.1 / 0.52);
    dlPad->SetLeftMargin (0.1 / 0.52);
    ulPad->SetRightMargin (0);
    clPad->SetRightMargin (0);
    dlPad->SetRightMargin (0);
    urPad->SetLeftMargin (0);
    crPad->SetLeftMargin (0);
    drPad->SetLeftMargin (0);
    urPad->SetRightMargin (0.03 / 0.48);
    crPad->SetRightMargin (0.03 / 0.48);
    drPad->SetRightMargin (0.03 / 0.48);
    ulPad->Draw ();
    clPad->Draw ();
    dlPad->Draw ();
    urPad->Draw ();
    crPad->Draw ();
    drPad->Draw ();

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    ulPad->cd (); 
    ulPad->SetLogx ();
    ulPad->SetLogy ();

    float ymin = 8e-6;
    float ymax = 3e1;

    h = (TH1D*) h_jet_trk_pt_ns[1]->Clone ("h");
    h->Reset ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fuPad);
    h->GetXaxis ()->SetLabelSize (0.028/fuPad);
    h->GetXaxis ()->SetTitleOffset (2.1*fuPad);
    h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]");
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fuPad);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    g = g_jet_trk_pt_ns_syst[0];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_jet_trk_pt_ns[0];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myBlue);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    g = g_jet_trk_pt_ns_syst[1];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myRed);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_jet_trk_pt_ns[1];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myRed);
    g->SetLineColor (myRed);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    g = g_jet_trk_pt_ns_bkg_syst[1];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (kBlack);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_jet_trk_pt_ns_bkg[1];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    //myText (0.24, 0.78, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
    myText (0.24, 0.12, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{0-20%}", 0.020/fuPad);
    myText (0.24, 0.06, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.020/fuPad);

    tl->SetTextAlign (22);
    tl->SetTextFont (42);
    tl->SetTextSize (0.022/fuPad);
    tl->DrawLatexNDC (0.1/0.52 + 0.5*(1.-0.1/0.52), 0.94, "#Delta#phi < #pi/8 (near-side)");


    urPad->cd (); 
    urPad->SetLogx ();
    urPad->SetLogy ();

    ymin = 8e-6;
    ymax = 3e1;

    h = (TH1D*) h_jet_trk_pt_as[1]->Clone ("h");
    h->Reset ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fuPad);
    h->GetXaxis ()->SetLabelSize (0.028/fuPad);
    h->GetXaxis ()->SetTitleOffset (2.1*fuPad);
    h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]");
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fuPad);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    g = g_jet_trk_pt_as_syst[0];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_jet_trk_pt_as[0];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myBlue);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    g = g_jet_trk_pt_as_syst[1];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myRed);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_jet_trk_pt_as[1];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myRed);
    g->SetLineColor (myRed);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    g = g_jet_trk_pt_as_bkg_syst[1];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (kBlack);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_jet_trk_pt_as_bkg[1];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myText (0.58, 0.78, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
    myMarkerText (0.10, 0.18, myRed, kFullCircle, "#it{p}+Pb jet-tagged events", 0.8, 0.020/fuPad, true);
    myMarkerText (0.10, 0.12, myBlue, kFullCircle, "#it{pp} jet-tagged events", 0.8, 0.020/fuPad, true);
    myMarkerText (0.10, 0.06, kBlack, kOpenCircle, "#it{p}+Pb mixed events", 0.8, 0.020/fuPad);

    tl->DrawLatexNDC (0.5*(1-0.03/0.48), 0.94, "#Delta#phi > 7#pi/8 (away-side)");


    clPad->cd (); 
    clPad->SetLogx ();
    clPad->SetLogy ();

    ymin = 8e-6;
    ymax = 3e1;

    h = (TH1D*) h_jet_trk_pt_ns_sig[1]->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    g = g_jet_trk_pt_ns_sig_syst[0];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_jet_trk_pt_ns_sig[0];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    g = g_jet_trk_pt_ns_sig_syst[1];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myRed);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_jet_trk_pt_ns_sig[1];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (myRed);
    g->SetLineWidth (2);
    g->SetMarkerColor (myRed);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myText (0.24, 0.26, kBlack, "Jet-hadron correlations", 0.020/fcPad);
    myText (0.24, 0.18, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV", strcmp (tag, "30GeVJets") == 0 ? 30 : 60), 0.020/fcPad);
    myText (0.24, 0.10, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fcPad);


    crPad->cd (); 
    crPad->SetLogx ();
    crPad->SetLogy ();

    ymin = 8e-6;
    ymax = 3e1;

    h = (TH1D*) h_jet_trk_pt_as_sig[1]->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    g = g_jet_trk_pt_as_sig_syst[0];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_jet_trk_pt_as_sig[0];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    g = g_jet_trk_pt_as_sig_syst[1];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myRed);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_jet_trk_pt_as_sig[1];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (myRed);
    g->SetLineWidth (2);
    g->SetMarkerColor (myRed);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);


    dlPad->cd (); 
    dlPad->SetLogx ();

    ymin = 0.53;
    ymax = (strcmp (tag, "30GeVJets") == 0 ? 1.45 : 1.17);

    h = (TH1D*) h_jet_trk_pt_ns_iaa->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}A} = #it{p}+Pb / #it{pp}");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    g = g_jet_trk_pt_ns_iaa_syst;
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_jet_trk_pt_ns_iaa;
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);


    drPad->cd (); 
    drPad->SetLogx ();

    ymin = 0.53;
    ymax = (strcmp (tag, "30GeVJets") == 0 ? 1.45 : 1.17);

    h = (TH1D*) h_jet_trk_pt_as_iaa->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}A} = #it{p}+Pb / #it{pp}");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    g = g_jet_trk_pt_as_iaa_syst;
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_jet_trk_pt_as_iaa;
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    c->SaveAs (Form ("Plots/JetTagged_HadronYields_Central_comparison_PtCh_%s.pdf", tag)); 
  }



  {
    const char* canvasName = "c_ljet_trk_pt";
    TCanvas* c = new TCanvas (canvasName, "", 1200, 1120);
    c->cd ();
    const double fuPad = 480./1120.;
    const double fdPad = 320./1120.;
    const double fcPad = 1.0 - fuPad - fdPad;
    const double fxPad = 0.42;
    TPad* ulPad = new TPad (Form ("%s_ulPad", canvasName), "", 0.0, 1.0-fuPad, fxPad+0.1, 1.0);
    TPad* clPad = new TPad (Form ("%s_clPad", canvasName), "", 0.0, fdPad, fxPad+0.1, 1.0-fuPad);
    TPad* dlPad = new TPad (Form ("%s_dlPad", canvasName), "", 0.0, 0.0, fxPad+0.1, fdPad);
    TPad* urPad = new TPad (Form ("%s_urPad", canvasName), "", fxPad+0.1, 1.0-fuPad, 1.0, 1.0);
    TPad* crPad = new TPad (Form ("%s_crPad", canvasName), "", fxPad+0.1, fdPad, 1.0, 1.0-fuPad);
    TPad* drPad = new TPad (Form ("%s_drPad", canvasName), "", fxPad+0.1, 0.0, 1.0, fdPad);

    ulPad->SetTopMargin (0.14);
    urPad->SetTopMargin (0.14);
    ulPad->SetBottomMargin (0);
    urPad->SetBottomMargin (0);
    clPad->SetTopMargin (0);
    crPad->SetTopMargin (0);
    clPad->SetBottomMargin (0);
    crPad->SetBottomMargin (0);
    dlPad->SetTopMargin (0);
    drPad->SetTopMargin (0);
    dlPad->SetBottomMargin (0.25);
    drPad->SetBottomMargin (0.25);

    ulPad->SetLeftMargin (0.1 / 0.52);
    clPad->SetLeftMargin (0.1 / 0.52);
    dlPad->SetLeftMargin (0.1 / 0.52);
    ulPad->SetRightMargin (0);
    clPad->SetRightMargin (0);
    dlPad->SetRightMargin (0);
    urPad->SetLeftMargin (0);
    crPad->SetLeftMargin (0);
    drPad->SetLeftMargin (0);
    urPad->SetRightMargin (0.03 / 0.48);
    crPad->SetRightMargin (0.03 / 0.48);
    drPad->SetRightMargin (0.03 / 0.48);
    ulPad->Draw ();
    clPad->Draw ();
    dlPad->Draw ();
    urPad->Draw ();
    crPad->Draw ();
    drPad->Draw ();

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    ulPad->cd (); 
    ulPad->SetLogx ();
    ulPad->SetLogy ();

    float ymin = 8e-6;
    float ymax = 3e1;

    h = (TH1D*) h_ljet_trk_pt_ns[1]->Clone ("h");
    h->Reset ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fuPad);
    h->GetXaxis ()->SetLabelSize (0.028/fuPad);
    h->GetXaxis ()->SetTitleOffset (2.1*fuPad);
    h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]");
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fuPad);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    g = g_ljet_trk_pt_ns_syst[0];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_ljet_trk_pt_ns[0];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myBlue);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    g = g_ljet_trk_pt_ns_syst[1];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myRed);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_ljet_trk_pt_ns[1];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myRed);
    g->SetLineColor (myRed);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    g = g_ljet_trk_pt_ns_bkg_syst[1];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (kBlack);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_ljet_trk_pt_ns_bkg[1];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    //myText (0.24, 0.78, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
    myText (0.24, 0.12, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{0-20%}", 0.020/fuPad);
    myText (0.24, 0.06, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.020/fuPad);

    tl->SetTextAlign (22);
    tl->SetTextFont (42);
    tl->SetTextSize (0.022/fuPad);
    tl->DrawLatexNDC (0.1/0.52 + 0.5*(1.-0.1/0.52), 0.94, "#Delta#phi < #pi/8 (near-side)");


    urPad->cd (); 
    urPad->SetLogx ();
    urPad->SetLogy ();

    ymin = 8e-6;
    ymax = 3e1;

    h = (TH1D*) h_ljet_trk_pt_as[1]->Clone ("h");
    h->Reset ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fuPad);
    h->GetXaxis ()->SetLabelSize (0.028/fuPad);
    h->GetXaxis ()->SetTitleOffset (2.1*fuPad);
    h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]");
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fuPad);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    g = g_ljet_trk_pt_as_syst[0];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_ljet_trk_pt_as[0];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myBlue);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    g = g_ljet_trk_pt_as_syst[1];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myRed);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_ljet_trk_pt_as[1];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myRed);
    g->SetLineColor (myRed);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    g = g_ljet_trk_pt_as_bkg_syst[1];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (kBlack);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_ljet_trk_pt_as_bkg[1];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myText (0.58, 0.78, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
    myMarkerText (0.10, 0.18, myRed, kFullCircle, "#it{p}+Pb jet-tagged events", 0.8, 0.020/fuPad, true);
    myMarkerText (0.10, 0.12, myBlue, kFullCircle, "#it{pp} jet-tagged events", 0.8, 0.020/fuPad, true);
    myMarkerText (0.10, 0.06, kBlack, kOpenCircle, "#it{p}+Pb mixed events", 0.8, 0.020/fuPad);

    tl->DrawLatexNDC (0.5*(1-0.03/0.48), 0.94, "#Delta#phi > 7#pi/8 (away-side)");


    clPad->cd (); 
    clPad->SetLogx ();
    clPad->SetLogy (); 

    ymin = 8e-6;
    ymax = 3e1;

    h = (TH1D*) h_ljet_trk_pt_ns_sig[1]->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    g = g_ljet_trk_pt_ns_sig_syst[0];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_ljet_trk_pt_ns_sig[0];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    g = g_ljet_trk_pt_ns_sig_syst[1];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myRed);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_ljet_trk_pt_ns_sig[1];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (myRed);
    g->SetLineWidth (2);
    g->SetMarkerColor (myRed);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myText (0.24, 0.26, kBlack, "Leading jet-hadron correlations", 0.020/fcPad);
    myText (0.24, 0.18, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV", strcmp (tag, "30GeVJets") == 0 ? 30 : 60), 0.020/fcPad);
    myText (0.24, 0.10, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fcPad);


    crPad->cd (); 
    crPad->SetLogx ();
    crPad->SetLogy (); 

    ymin = 8e-6;
    ymax = 3e1;

    h = (TH1D*) h_ljet_trk_pt_as_sig[1]->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    g = g_ljet_trk_pt_as_sig_syst[0];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_ljet_trk_pt_as_sig[0];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    g = g_ljet_trk_pt_as_sig_syst[1];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myRed);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_ljet_trk_pt_as_sig[1];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (myRed);
    g->SetLineWidth (2);
    g->SetMarkerColor (myRed);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);


    dlPad->cd (); 
    dlPad->SetLogx ();

    ymin = 0.53;
    ymax = (strcmp (tag, "30GeVJets") == 0 ? 1.45 : 1.17);

    h = (TH1D*) h_ljet_trk_pt_ns_iaa->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}A} = #it{p}+Pb / #it{pp}");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    g = g_ljet_trk_pt_ns_iaa_syst;
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_ljet_trk_pt_ns_iaa;
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);


    drPad->cd (); 
    drPad->SetLogx ();

    ymin = 0.53;
    ymax = (strcmp (tag, "30GeVJets") == 0 ? 1.45 : 1.17);

    h = (TH1D*) h_ljet_trk_pt_as_iaa->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}A} = #it{p}+Pb / #it{pp}");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    g = g_ljet_trk_pt_as_iaa_syst;
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_ljet_trk_pt_as_iaa;
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    c->SaveAs (Form ("Plots/LeadingJetTagged_HadronYields_Central_comparison_PtCh_%s.pdf", tag)); 
  }



  {
    const char* canvasName = "c_jet_trk_pt_ns_syst";
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    c->cd (); 
    c->SetLogx ();
    //c->SetLogy ();

    float ymin = -20;
    float ymax = 20;


    for (int iSys : {0, 1}) {
      c->Clear ();

      h = (TH1D*) h_jet_trk_pt_ns[iSys]->Clone ("h");
      h->Reset ();
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      //h->GetXaxis ()->SetTitleSize (0.028);
      //h->GetXaxis ()->SetLabelSize (0.028);
      h->GetYaxis ()->SetTitle ("#delta N_{ch} / N_{ch} [%]");
      //h->GetYaxis ()->SetTitleSize (0.028);
      //h->GetYaxis ()->SetLabelSize (0.028);

      h->SetLineWidth (1);
      h->SetLineStyle (2);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      for (int iVar = 0; iVar < nVar; iVar++) {
        h = (TH1D*) h_jet_trk_pt_ns_syst[iSys][iVar]->Clone ("htemp");
        SaveRelativeErrors (h, h_jet_trk_pt_ns[iSys], true);
        for (int iX = 1; iX <= h->GetNbinsX (); iX++) h->SetBinContent (iX, 100*h->GetBinContent (iX) - 100);
        g = make_graph (h);
        ResetXErrors (g);
        ResetTGAEErrors (g);

        g->SetLineColor (varStyles[variations[iVar]].first);
        g->SetLineStyle (varStyles[variations[iVar]].second);
        g->SetLineWidth (3);
        ((TGAE*) g->Clone ())->Draw ("L");
        SaferDelete (&g);
        SaferDelete (&h);
      }

      myText (0.22, 0.87, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
      if (iSys == 0) {
        myText (0.22, 0.828, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
      }
      else {
        myText (0.22, 0.828, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, 0-20%", 0.032);
      }
      myText (0.22, 0.786, kBlack, "#Delta#phi < #pi/8 (near-side)", 0.032);
      for (int iVar = 0; iVar < nVar; iVar++)
        myLineColorText (0.25, 0.744-iVar*0.040, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 1.0, 0.028);
      

      c->SaveAs (Form ("Plots/Systematics/TotalJetTaggedYield_%s_nearside_ptch_%s_syst.pdf", iSys == 0 ? "pp" : "pPb", tag));
    }
  }


  {
    const char* canvasName = "c_jet_trk_pt_as_syst";
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    c->cd (); 
    c->SetLogx ();
    //c->SetLogy ();

    float ymin = -20;
    float ymax = 20;


    for (int iSys : {0, 1}) {
      c->Clear ();

      h = (TH1D*) h_jet_trk_pt_as[iSys]->Clone ("h");
      h->Reset ();
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      //h->GetXaxis ()->SetTitleSize (0.028);
      //h->GetXaxis ()->SetLabelSize (0.028);
      h->GetYaxis ()->SetTitle ("#delta N_{ch} / N_{ch} [%]");
      //h->GetYaxis ()->SetTitleSize (0.028);
      //h->GetYaxis ()->SetLabelSize (0.028);

      h->SetLineWidth (1);
      h->SetLineStyle (2);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      for (int iVar = 0; iVar < nVar; iVar++) {
        h = (TH1D*) h_jet_trk_pt_as_syst[iSys][iVar]->Clone ("htemp");
        SaveRelativeErrors (h, h_jet_trk_pt_as[iSys], true);
        for (int iX = 1; iX <= h->GetNbinsX (); iX++) h->SetBinContent (iX, 100*h->GetBinContent (iX) - 100);
        g = make_graph (h);
        ResetXErrors (g);
        ResetTGAEErrors (g);

        g->SetLineColor (varStyles[variations[iVar]].first);
        g->SetLineStyle (varStyles[variations[iVar]].second);
        g->SetLineWidth (3);
        ((TGAE*) g->Clone ())->Draw ("L");
        SaferDelete (&g);
        SaferDelete (&h);
      }

      myText (0.22, 0.87, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
      if (iSys == 0) {
        myText (0.22, 0.828, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
      }
      else {
        myText (0.22, 0.828, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, 0-20%", 0.032);
      }
      myText (0.22, 0.786, kBlack, "#Delta#phi > 7#pi/8 (away-side)", 0.032);
      for (int iVar = 0; iVar < nVar; iVar++)
        myLineColorText (0.25, 0.744-iVar*0.040, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 1.0, 0.028);
      

      c->SaveAs (Form ("Plots/Systematics/TotalJetTaggedYield_%s_awayside_ptch_%s_syst.pdf", iSys == 0 ? "pp" : "pPb", tag));
    }
  }


  {
    const char* canvasName = "c_jet_trk_pt_ns_sig_syst";
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    c->cd (); 
    c->SetLogx ();
    //c->SetLogy ();

    float ymin = -20;
    float ymax = 20;


    for (int iSys : {0, 1}) {
      c->Clear ();

      h = (TH1D*) h_jet_trk_pt_ns[iSys]->Clone ("h");
      h->Reset ();
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      //h->GetXaxis ()->SetTitleSize (0.028);
      //h->GetXaxis ()->SetLabelSize (0.028);
      h->GetYaxis ()->SetTitle ("#delta N_{ch} / N_{ch} [%]");
      //h->GetYaxis ()->SetTitleSize (0.028);
      //h->GetYaxis ()->SetLabelSize (0.028);

      h->SetLineWidth (1);
      h->SetLineStyle (2);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      for (int iVar = 0; iVar < nVar; iVar++) {
        h = (TH1D*) h_jet_trk_pt_ns_sig_syst[iSys][iVar]->Clone ("htemp");
        SaveRelativeErrors (h, h_jet_trk_pt_ns_sig[iSys], true);
        for (int iX = 1; iX <= h->GetNbinsX (); iX++) h->SetBinContent (iX, 100*h->GetBinContent (iX) - 100);
        g = make_graph (h);
        ResetXErrors (g);
        ResetTGAEErrors (g);

        g->SetLineColor (varStyles[variations[iVar]].first);
        g->SetLineStyle (varStyles[variations[iVar]].second);
        g->SetLineWidth (3);
        ((TGAE*) g->Clone ())->Draw ("L");
        SaferDelete (&g);
        SaferDelete (&h);
      }

      myText (0.22, 0.87, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
      if (iSys == 0) {
        myText (0.22, 0.828, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
      }
      else {
        myText (0.22, 0.828, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, 0-20%", 0.032);
      }
      myText (0.22, 0.768, kBlack, "#Delta#phi < #pi/8 (near-side)", 0.032);
      for (int iVar = 0; iVar < nVar; iVar++)
        myLineColorText (0.25, 0.744-iVar*0.040, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 1.0, 0.028);
      

      c->SaveAs (Form ("Plots/Systematics/SignalJetTaggedYield_%s_nearside_ptch_%s_syst.pdf", iSys == 0 ? "pp" : "pPb", tag));
    }
  }


  {
    const char* canvasName = "c_jet_trk_pt_as_sig_syst";
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    c->cd (); 
    c->SetLogx ();
    //c->SetLogy ();

    float ymin = -20;
    float ymax = 20;


    for (int iSys : {0, 1}) {
      c->Clear ();

      h = (TH1D*) h_jet_trk_pt_as[iSys]->Clone ("h");
      h->Reset ();
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      //h->GetXaxis ()->SetTitleSize (0.028);
      //h->GetXaxis ()->SetLabelSize (0.028);
      h->GetYaxis ()->SetTitle ("#delta N_{ch} / N_{ch} [%]");
      //h->GetYaxis ()->SetTitleSize (0.028);
      //h->GetYaxis ()->SetLabelSize (0.028);

      h->SetLineWidth (1);
      h->SetLineStyle (2);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      for (int iVar = 0; iVar < nVar; iVar++) {
        h = (TH1D*) h_jet_trk_pt_as_sig_syst[iSys][iVar]->Clone ("htemp");
        SaveRelativeErrors (h, h_jet_trk_pt_as_sig[iSys], true);
        for (int iX = 1; iX <= h->GetNbinsX (); iX++) h->SetBinContent (iX, 100*h->GetBinContent (iX) - 100);
        g = make_graph (h);
        ResetXErrors (g);
        ResetTGAEErrors (g);

        g->SetLineColor (varStyles[variations[iVar]].first);
        g->SetLineStyle (varStyles[variations[iVar]].second);
        g->SetLineWidth (3);
        ((TGAE*) g->Clone ())->Draw ("L");
        SaferDelete (&g);
        SaferDelete (&h);
      }

      myText (0.22, 0.87, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
      if (iSys == 0) {
        myText (0.22, 0.828, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
      }
      else {
        myText (0.22, 0.828, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, 0-20%", 0.032);
      }
      myText (0.22, 0.786, kBlack, "#Delta#phi > 7#pi/8 (away-side)", 0.032);
      for (int iVar = 0; iVar < nVar; iVar++)
        myLineColorText (0.25, 0.744-iVar*0.040, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 1.0, 0.028);
      

      c->SaveAs (Form ("Plots/Systematics/SignalJetTaggedYield_%s_awayside_ptch_%s_syst.pdf", iSys == 0 ? "pp" : "pPb", tag));
    }
  }


  {
    const char* canvasName = "c_jet_trk_pt_ns_iaa_syst";
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    c->cd (); 
    c->SetLogx ();
    //c->SetLogy ();

    float ymin = -20;
    float ymax = 20;


    h = (TH1D*) h_jet_trk_pt_ns_iaa->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    //h->GetXaxis ()->SetTitleSize (0.028);
    //h->GetXaxis ()->SetLabelSize (0.028);
    h->GetYaxis ()->SetTitle ("#delta I_{pA} / I_{pA} [%]");
    //h->GetYaxis ()->SetTitleSize (0.028);
    //h->GetYaxis ()->SetLabelSize (0.028);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    for (int iVar = 0; iVar < nVar; iVar++) {
      h = (TH1D*) h_jet_trk_pt_ns_iaa_syst[iVar]->Clone ("htemp");
      SaveRelativeErrors (h, h_jet_trk_pt_ns_iaa, true);
      for (int iX = 1; iX <= h->GetNbinsX (); iX++) h->SetBinContent (iX, 100*h->GetBinContent (iX) - 100);
      g = make_graph (h);
      ResetXErrors (g);
      ResetTGAEErrors (g);

      g->SetLineColor (varStyles[variations[iVar]].first);
      g->SetLineStyle (varStyles[variations[iVar]].second);
      g->SetLineWidth (3);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
      SaferDelete (&h);
    }

    myText (0.22, 0.87, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.828, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    myText (0.22, 0.786, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, 0-20%", 0.032);
    myText (0.22, 0.744, kBlack, "#Delta#phi < #pi/8 (near-side)", 0.032);
    for (int iVar = 0; iVar < nVar; iVar++)
      myLineColorText (0.25, 0.702-iVar*0.040, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 1.0, 0.028);
    

    c->SaveAs (Form ("Plots/Systematics/JetTagged_IpA_nearside_ptch_%s_syst.pdf", tag));
  }


  {
    const char* canvasName = "c_jet_trk_pt_as_iaa_syst";
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    c->cd (); 
    c->SetLogx ();
    //c->SetLogy ();

    float ymin = -20;
    float ymax = 20;


    h = (TH1D*) h_jet_trk_pt_as_iaa->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    //h->GetXaxis ()->SetTitleSize (0.028);
    //h->GetXaxis ()->SetLabelSize (0.028);
    h->GetYaxis ()->SetTitle ("#delta I_{pA} / I_{pA} [%]");
    //h->GetYaxis ()->SetTitleSize (0.028);
    //h->GetYaxis ()->SetLabelSize (0.028);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    for (int iVar = 0; iVar < nVar; iVar++) {
      h = (TH1D*) h_jet_trk_pt_as_iaa_syst[iVar]->Clone ("htemp");
      SaveRelativeErrors (h, h_jet_trk_pt_as_iaa, true);
      for (int iX = 1; iX <= h->GetNbinsX (); iX++) h->SetBinContent (iX, 100*h->GetBinContent (iX) - 100);
      g = make_graph (h);
      ResetXErrors (g);
      ResetTGAEErrors (g);

      g->SetLineColor (varStyles[variations[iVar]].first);
      g->SetLineStyle (varStyles[variations[iVar]].second);
      g->SetLineWidth (3);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
      SaferDelete (&h);
    }

    myText (0.22, 0.87, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.828, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    myText (0.22, 0.786, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, 0-20%", 0.032);
    myText (0.22, 0.744, kBlack, "#Delta#phi > 7#pi/8 (away-side)", 0.032);
    for (int iVar = 0; iVar < nVar; iVar++)
      myLineColorText (0.25, 0.702-iVar*0.040, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 1.0, 0.028);
    

    c->SaveAs (Form ("Plots/Systematics/JetTagged_IpA_awayside_ptch_%s_syst.pdf", tag));
  }

}


#endif
