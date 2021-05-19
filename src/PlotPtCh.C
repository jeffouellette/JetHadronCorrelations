#ifndef __JetHadronCorrelatorPlotPtCh_C__
#define __JetHadronCorrelatorPlotPtCh_C__

#include "Params.h"
#include "OutTree.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Variations.h"

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
const Color_t myGreen = (Color_t) tcolor->GetColor ( 54, 167,  80);
const Color_t myOrange = (Color_t) tcolor->GetColor (255,  68,   0);

const Color_t systColors[10] = {kRed+1, kAzure-2, kGreen+2, kViolet-3, kMagenta, kCyan+1, kOrange-3, kGreen-7, kBlue+1, kPink-5};

TLine* l = new TLine ();
TLatex* tl = new TLatex ();


void DoOffset (TH1D* h, const char* tag, const int iCent) {

  ifstream offsetsFile;
  offsetsFile.open (Form ("%s/aux/%s_IAAOffsets.dat", workPath.Data (), tag));

  string inTag, pTChSel, centStr, offpercerr;
  double offset = 0, offerr = 0;

  while (!offsetsFile.eof ()) {
    offsetsFile >> inTag >> pTChSel >> centStr >> offset >> offerr >> offpercerr;

    if (TString (centStr.c_str ()) == TString (Form ("%i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]))) {

      for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
        if (pTChStrCuts[TString (pTChSel.c_str ())].first <= h->GetBinCenter (iX) && h->GetBinCenter (iX) <= pTChStrCuts[TString (pTChSel.c_str ())].second)
          h->SetBinContent (iX, h->GetBinContent (iX)-offset);
      }
    }
  }

  return;
}



void DoOffset (TGAE* g, const char* tag, const int iCent) {

  ifstream offsetsFile;
  offsetsFile.open (Form ("%s/aux/%s_IAAOffsets.dat", workPath.Data (), tag));

  string inTag, pTChSel, centStr, offpercerr;
  double offset = 0, offerr = 0;

  while (!offsetsFile.eof ()) {
    offsetsFile >> inTag >> pTChSel >> centStr >> offset >> offerr >> offpercerr;

    if (TString (centStr.c_str ()) == Form ("%i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent])) {

      double x, y;
      for (int i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (pTChStrCuts[TString (pTChSel.c_str ())].first <= x && x <= pTChStrCuts[TString (pTChSel.c_str ())].second)
          g->SetPoint (i, x, y-offset);
      }
    }
  }

  return;
}


void PlotPtCh (const char* tag, const char* inFileTag) {

  TFile* inFile = nullptr;

  TH1D*  h_evt_counts_ref = nullptr;
  TH1D*  h_jet_counts_ref = nullptr;
  TH1D*  h_evt_counts_ref_bkg = nullptr;
  TH1D*  h_jet_counts_ref_bkg = nullptr;
  TH1D** h_evt_counts = new TH1D*[numZdcCentBins];
  TH1D** h_jet_counts = new TH1D*[numZdcCentBins];
  TH1D** h_evt_counts_bkg = new TH1D*[numZdcCentBins];
  TH1D** h_jet_counts_bkg = new TH1D*[numZdcCentBins];

  TH1D*  h_jet_trk_pt_ns_ref = nullptr;
  TH1D*  h_jet_trk_pt_as_ref = nullptr;
  TH1D*  h_jet_trk_pt_ns_ref_bkg = nullptr;
  TH1D*  h_jet_trk_pt_as_ref_bkg = nullptr;
  TH1D** h_jet_trk_pt_ns = new TH1D*[numZdcCentBins];
  TH1D** h_jet_trk_pt_as = new TH1D*[numZdcCentBins];
  TH1D** h_jet_trk_pt_ns_bkg = new TH1D*[numZdcCentBins];
  TH1D** h_jet_trk_pt_as_bkg = new TH1D*[numZdcCentBins];

  TH1D*  h_jet_trk_pt_ns_ref_sig = nullptr;
  TH1D*  h_jet_trk_pt_as_ref_sig = nullptr;
  TH1D** h_jet_trk_pt_ns_sig = new TH1D*[numZdcCentBins];
  TH1D** h_jet_trk_pt_as_sig = new TH1D*[numZdcCentBins];

  TH1D** h_jet_trk_pt_ns_iaa = new TH1D*[numZdcCentBins];
  TH1D** h_jet_trk_pt_as_iaa = new TH1D*[numZdcCentBins];

  TGAE*  g_jet_trk_pt_ns_ref_syst = nullptr;
  TGAE*  g_jet_trk_pt_as_ref_syst = nullptr;
  TGAE*  g_jet_trk_pt_ns_ref_bkg_syst = nullptr;
  TGAE*  g_jet_trk_pt_as_ref_bkg_syst = nullptr;
  TGAE** g_jet_trk_pt_ns_syst = new TGAE*[numZdcCentBins];
  TGAE** g_jet_trk_pt_as_syst = new TGAE*[numZdcCentBins];
  TGAE** g_jet_trk_pt_ns_bkg_syst = new TGAE*[numZdcCentBins];
  TGAE** g_jet_trk_pt_as_bkg_syst = new TGAE*[numZdcCentBins];

  TGAE*  g_jet_trk_pt_ns_ref_sig_syst = nullptr;
  TGAE*  g_jet_trk_pt_as_ref_sig_syst = nullptr;
  TGAE** g_jet_trk_pt_ns_sig_syst = new TGAE*[numZdcCentBins];
  TGAE** g_jet_trk_pt_as_sig_syst = new TGAE*[numZdcCentBins];

  TGAE** g_jet_trk_pt_ns_iaa_syst = new TGAE*[numZdcCentBins];
  TGAE** g_jet_trk_pt_as_iaa_syst = new TGAE*[numZdcCentBins];


  TH1D**  h_jet_trk_pt_ns_ref_syst = Get1DArray <TH1D*> (nVar);
  TH1D**  h_jet_trk_pt_as_ref_syst = Get1DArray <TH1D*> (nVar);
  TH1D**  h_jet_trk_pt_ns_ref_bkg_syst = Get1DArray <TH1D*> (nVar);
  TH1D**  h_jet_trk_pt_as_ref_bkg_syst = Get1DArray <TH1D*> (nVar);
  TH1D*** h_jet_trk_pt_ns_syst = Get2DArray <TH1D*> (numZdcCentBins, nVar);
  TH1D*** h_jet_trk_pt_as_syst = Get2DArray <TH1D*> (numZdcCentBins, nVar);
  TH1D*** h_jet_trk_pt_ns_bkg_syst = Get2DArray <TH1D*> (numZdcCentBins, nVar);
  TH1D*** h_jet_trk_pt_as_bkg_syst = Get2DArray <TH1D*> (numZdcCentBins, nVar);

  TH1D**  h_jet_trk_pt_ns_ref_sig_syst = Get1DArray <TH1D*> (nVar);
  TH1D**  h_jet_trk_pt_as_ref_sig_syst = Get1DArray <TH1D*> (nVar);
  TH1D*** h_jet_trk_pt_ns_sig_syst = Get2DArray <TH1D*> (numZdcCentBins, nVar);
  TH1D*** h_jet_trk_pt_as_sig_syst = Get2DArray <TH1D*> (numZdcCentBins, nVar);

  TH1D*** h_jet_trk_pt_ns_iaa_syst = Get2DArray <TH1D*> (numZdcCentBins, nVar);
  TH1D*** h_jet_trk_pt_as_iaa_syst = Get2DArray <TH1D*> (numZdcCentBins, nVar);


  {
    TString inFileName = inFileTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/PlotPtCh_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    h_evt_counts_ref = (TH1D*) inFile->Get ("h_evt_counts_ref_Nominal");
    h_jet_counts_ref = (TH1D*) inFile->Get ("h_jet_counts_ref_Nominal");

    h_evt_counts_ref_bkg = (TH1D*) inFile->Get ("h_evt_counts_ref_bkg_Nominal");
    h_jet_counts_ref_bkg = (TH1D*) inFile->Get ("h_jet_counts_ref_bkg_Nominal");

    h_jet_trk_pt_ns_ref = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_ref_Nominal");
    h_jet_trk_pt_as_ref = (TH1D*) inFile->Get ("h_jet_trk_pt_as_ref_Nominal");

    h_jet_trk_pt_ns_ref_bkg = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_ref_bkg_Nominal");
    h_jet_trk_pt_as_ref_bkg = (TH1D*) inFile->Get ("h_jet_trk_pt_as_ref_bkg_Nominal");

    h_jet_trk_pt_ns_ref_sig = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_ref_sig_Nominal");
    h_jet_trk_pt_as_ref_sig = (TH1D*) inFile->Get ("h_jet_trk_pt_as_ref_sig_Nominal");

    g_jet_trk_pt_ns_ref_syst = (TGAE*) inFile->Get ("g_jet_trk_pt_ns_ref_syst");
    g_jet_trk_pt_as_ref_syst = (TGAE*) inFile->Get ("g_jet_trk_pt_as_ref_syst");

    g_jet_trk_pt_ns_ref_bkg_syst = (TGAE*) inFile->Get ("g_jet_trk_pt_ns_ref_bkg_syst");
    g_jet_trk_pt_as_ref_bkg_syst = (TGAE*) inFile->Get ("g_jet_trk_pt_as_ref_bkg_syst");

    g_jet_trk_pt_ns_ref_sig_syst = (TGAE*) inFile->Get ("g_jet_trk_pt_ns_ref_sig_syst");
    g_jet_trk_pt_as_ref_sig_syst = (TGAE*) inFile->Get ("g_jet_trk_pt_as_ref_sig_syst");

    for (int iCent = 0; iCent < numZdcCentBins; iCent++) {

      h_evt_counts[iCent] = (TH1D*) inFile->Get (Form ("h_evt_counts_pPb_iCent%i_Nominal", iCent));
      h_jet_counts[iCent] = (TH1D*) inFile->Get (Form ("h_jet_counts_pPb_iCent%i_Nominal", iCent));

      h_evt_counts_bkg[iCent] = (TH1D*) inFile->Get (Form ("h_evt_counts_pPb_bkg_iCent%i_Nominal", iCent));
      h_jet_counts_bkg[iCent] = (TH1D*) inFile->Get (Form ("h_jet_counts_pPb_bkg_iCent%i_Nominal", iCent));

      h_jet_trk_pt_ns[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_pPb_iCent%i_Nominal", iCent));
      h_jet_trk_pt_as[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_pPb_iCent%i_Nominal", iCent));
      h_jet_trk_pt_ns_bkg[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_pPb_bkg_iCent%i_Nominal", iCent));
      h_jet_trk_pt_as_bkg[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_pPb_bkg_iCent%i_Nominal", iCent));

      h_jet_trk_pt_ns_sig[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_pPb_sig_iCent%i_Nominal", iCent));
      h_jet_trk_pt_as_sig[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_pPb_sig_iCent%i_Nominal", iCent));

      h_jet_trk_pt_ns_iaa[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_iaa_iCent%i_Nominal", iCent));
      h_jet_trk_pt_as_iaa[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_iaa_iCent%i_Nominal", iCent));

      g_jet_trk_pt_ns_syst[iCent] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_ns_syst_iCent%i", iCent));
      g_jet_trk_pt_as_syst[iCent] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_as_syst_iCent%i", iCent));

      g_jet_trk_pt_ns_bkg_syst[iCent] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_ns_bkg_syst_iCent%i", iCent));
      g_jet_trk_pt_as_bkg_syst[iCent] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_as_bkg_syst_iCent%i", iCent));

      g_jet_trk_pt_ns_sig_syst[iCent] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_ns_sig_syst_iCent%i", iCent));
      g_jet_trk_pt_as_sig_syst[iCent] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_as_sig_syst_iCent%i", iCent));

      g_jet_trk_pt_ns_iaa_syst[iCent] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_ns_iaa_syst_iCent%i", iCent));
      g_jet_trk_pt_as_iaa_syst[iCent] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_as_iaa_syst_iCent%i", iCent));
    }


    for (int iVar = 1; iVar < nVar; iVar++) {

      h_jet_trk_pt_ns_ref_syst[iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_ref_%s", variations[iVar].Data ()));
      h_jet_trk_pt_as_ref_syst[iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_ref_%s", variations[iVar].Data ()));

      h_jet_trk_pt_ns_ref_sig_syst[iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_ref_sig_%s", variations[iVar].Data ()));
      h_jet_trk_pt_as_ref_sig_syst[iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_ref_sig_%s", variations[iVar].Data ()));

      for (int iCent = 0; iCent < numZdcCentBins; iCent++) {

        h_jet_trk_pt_ns_syst[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_pPb_iCent%i_%s", iCent, variations[iVar].Data ()));
        h_jet_trk_pt_as_syst[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_pPb_iCent%i_%s", iCent, variations[iVar].Data ()));

        h_jet_trk_pt_ns_bkg_syst[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_pPb_bkg_iCent%i_%s", iCent, variations[iVar].Data ()));
        h_jet_trk_pt_as_bkg_syst[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_pPb_bkg_iCent%i_%s", iCent, variations[iVar].Data ()));

        h_jet_trk_pt_ns_sig_syst[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_pPb_sig_iCent%i_%s", iCent, variations[iVar].Data ()));
        h_jet_trk_pt_as_sig_syst[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_pPb_sig_iCent%i_%s", iCent, variations[iVar].Data ()));

        h_jet_trk_pt_ns_iaa_syst[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_iaa_iCent%i_%s", iCent, variations[iVar].Data ()));
        h_jet_trk_pt_as_iaa_syst[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_iaa_iCent%i_%s", iCent, variations[iVar].Data ()));

      }
    }
  }



  l->SetLineStyle (7);
  l->SetLineWidth (1);
  l->SetLineColor (kBlack);



  for (int iCent = 0; iCent < numZdcCentBins; iCent++) {
    const char* canvasName = Form ("c_jet_trk_pt_iCent%i", iCent);
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
    h = (TH1D*) h_jet_trk_pt_ns_ref->Clone ("h");
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

    g = g_jet_trk_pt_ns_ref_syst;
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    h = h_jet_trk_pt_ns_ref;
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

    g = g_jet_trk_pt_ns_ref_bkg_syst;
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myPurple);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    h = h_jet_trk_pt_ns_ref_bkg;
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myPurple);
    g->SetLineColor (myPurple);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    g = g_jet_trk_pt_ns_syst[iCent];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myRed);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    h = h_jet_trk_pt_ns[iCent];
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

    g = g_jet_trk_pt_ns_bkg_syst[iCent];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myGreen);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    h = h_jet_trk_pt_ns_bkg[iCent];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myGreen);
    g->SetLineColor (myGreen);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    //myText (0.24, 0.78, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
    myText (0.24, 0.12, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.020/fuPad);
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

    h = (TH1D*) h_jet_trk_pt_as_ref->Clone ("h");
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

    g = g_jet_trk_pt_as_ref_syst;
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    h = h_jet_trk_pt_as_ref;
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

    g = g_jet_trk_pt_as_ref_bkg_syst;
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myPurple);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    h = h_jet_trk_pt_as_ref_bkg;
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myPurple);
    g->SetLineColor (myPurple);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    g = g_jet_trk_pt_as_syst[iCent];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myRed);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    h = h_jet_trk_pt_as[iCent];
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

    g = g_jet_trk_pt_as_bkg_syst[iCent];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myGreen);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    h = h_jet_trk_pt_as_bkg[iCent];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myGreen);
    g->SetLineColor (myGreen);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myText (0.58, 0.78, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
    myBoxText2 (0.10, 0.24, myRed, kFullCircle, "#it{p}+Pb jet-tagged events", 0.8, 0.020/fuPad, true);
    myBoxText2 (0.10, 0.18, myBlue, kFullCircle, "#it{pp} jet-tagged events", 0.8, 0.020/fuPad, true);
    myBoxText2 (0.10, 0.12, myGreen, kOpenCircle, "#it{p}+Pb bkgd.", 0.8, 0.020/fuPad);
    myBoxText2 (0.10, 0.06, myPurple, kOpenCircle, "#it{pp} bkgd.", 0.8, 0.020/fuPad);

    tl->DrawLatexNDC (0.5*(1-0.03/0.48), 0.94, "#Delta#phi > 7#pi/8 (away-side)");


    clPad->cd (); 
    clPad->SetLogx ();
    clPad->SetLogy ();

    ymin = 8e-6;
    ymax = 3e1;

    h = (TH1D*) h_jet_trk_pt_ns_ref_sig->Clone ("h");
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

    g = g_jet_trk_pt_ns_ref_sig_syst;
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    h = h_jet_trk_pt_ns_ref_sig;
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

    g = g_jet_trk_pt_ns_sig_syst[iCent];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myRed);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    h = h_jet_trk_pt_ns_sig[iCent];
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
    myText (0.24, 0.18, kBlack, Form ("#it{p}_{T}^{jet} > %s", GetJetPtStr (tag).Data ()), 0.020/fcPad);
    myText (0.24, 0.10, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fcPad);


    crPad->cd (); 
    crPad->SetLogx ();
    crPad->SetLogy ();

    ymin = 8e-6;
    ymax = 3e1;

    h = (TH1D*) h_jet_trk_pt_as_ref_sig->Clone ("h");
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

    g = g_jet_trk_pt_as_ref_sig_syst;
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    h = h_jet_trk_pt_as_ref_sig;
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

    g = g_jet_trk_pt_as_sig_syst[iCent];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myRed);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    h = h_jet_trk_pt_as_sig[iCent];
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
    ymax = (strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.45 : 1.17);

    h = (TH1D*) h_jet_trk_pt_ns_iaa[iCent]->Clone ("h");
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

    g = g_jet_trk_pt_ns_iaa_syst[iCent];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myRed);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    h = h_jet_trk_pt_ns_iaa[iCent];
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


    g = (TGAE*) g_jet_trk_pt_ns_iaa_syst[iCent]->Clone ();
    DoOffset (g, tag, iCent);
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = (TH1D*) h_jet_trk_pt_ns_iaa[iCent]->Clone ("htemp");
    DoOffset (h, tag, iCent);
    g = make_graph (h);
    SaferDelete (&h);
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
    ymax = (strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.45 : 1.17);

    h = (TH1D*) h_jet_trk_pt_as_iaa[iCent]->Clone ("h");
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

    g = g_jet_trk_pt_as_iaa_syst[iCent];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myRed);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    h = h_jet_trk_pt_as_iaa[iCent];
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

    g = (TGAE*) g_jet_trk_pt_as_iaa_syst[iCent]->Clone ("g");
    DoOffset (g, tag, iCent);
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = (TH1D*) h_jet_trk_pt_as_iaa[iCent]->Clone ("htemp");
    DoOffset (h, tag, iCent);
    g = make_graph (h);
    SaferDelete (&h);
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

    c->SaveAs (Form ("%s/Plots/PtCh/JetTagged_HadronYields_%i-%i%%_comparison_PtCh_%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], tag)); 
  }



  for (int iCent = 0; iCent < numZdcCentBins; iCent++) {
    const char* canvasName = Form ("c_jet_trk_pt_FcalVsZdc_iCent%i", iCent);
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

    h = (TH1D*) h_jet_trk_pt_ns_ref->Clone ("h");
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

    h = h_jet_trk_pt_ns_ref;
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_pt_ns_ref_bkg;
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_pt_ns[iCent];
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

    h = h_jet_trk_pt_ns_bkg[iCent];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myGreen);
    g->SetLineColor (myGreen);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    int iVar = 0;
    while (iVar < nVar && strcmp (variations[iVar], "FcalCentVar") != 0) iVar++;
    if (iVar == nVar)
      std::cout << "Cannot find FCal centrality binned result??? Please check!" << std::endl;

    h = h_jet_trk_pt_ns_syst[iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullSquare);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myBlue);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_pt_ns_bkg_syst[iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullSquare);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myOrange);
    g->SetLineColor (myOrange);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);



    //myText (0.24, 0.78, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
    myText (0.24, 0.12, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.020/fuPad);
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

    h = (TH1D*) h_jet_trk_pt_as_ref->Clone ("h");
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

    h = h_jet_trk_pt_as_ref;
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");

    h = h_jet_trk_pt_as_ref_bkg;
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");

    h = h_jet_trk_pt_as[iCent];
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

    h = h_jet_trk_pt_as_bkg[iCent];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myGreen);
    g->SetLineColor (myGreen);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_pt_as_syst[iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullSquare);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myBlue);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_pt_as_bkg_syst[iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullSquare);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myOrange);
    g->SetLineColor (myOrange);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myText (0.58, 0.78, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
    myMarkerText (0.10, 0.36, myRed, kFullCircle, Form ("#it{p}+Pb #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.8, 0.020/fuPad, true);
    myMarkerText (0.10, 0.30, myGreen, kFullCircle, Form ("#it{p}+Pb bkgd., #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.8, 0.020/fuPad, true);
    myMarkerText (0.10, 0.24, myBlue, kFullSquare, Form ("#it{p}+Pb #bf{FCal %i-%i%%}", fcalCentPercs[iCent+1], fcalCentPercs[iCent]), 0.8, 0.020/fuPad, true);
    myMarkerText (0.10, 0.18, myOrange, kFullSquare, Form ("#it{p}+Pb bkgd., #bf{FCal %i-%i%%}", fcalCentPercs[iCent+1], fcalCentPercs[iCent]), 0.8, 0.020/fuPad, true);
    myMarkerText (0.10, 0.12, kBlack, kFullCircle, "#it{pp} jet-tagged events", 0.8, 0.020/fuPad);
    myMarkerText (0.10, 0.06, kBlack, kOpenCircle, "#it{pp} bkgd.", 0.8, 0.020/fuPad);

    tl->DrawLatexNDC (0.5*(1-0.03/0.48), 0.94, "#Delta#phi > 7#pi/8 (away-side)");


    clPad->cd (); 
    clPad->SetLogx ();
    clPad->SetLogy ();

    ymin = 8e-6;
    ymax = 3e1;

    h = (TH1D*) h_jet_trk_pt_ns_ref_sig->Clone ("h");
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

    h = h_jet_trk_pt_ns_ref_sig;
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_pt_ns_sig[iCent];
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

    h = h_jet_trk_pt_ns_sig_syst[iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullSquare);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myText (0.24, 0.26, kBlack, "Jet-hadron correlations", 0.020/fcPad);
    myText (0.24, 0.18, kBlack, Form ("#it{p}_{T}^{jet} > %s", GetJetPtStr (tag).Data ()), 0.020/fcPad);
    myText (0.24, 0.10, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fcPad);


    crPad->cd (); 
    crPad->SetLogx ();
    crPad->SetLogy ();

    ymin = 8e-6;
    ymax = 3e1;

    h = (TH1D*) h_jet_trk_pt_as_ref_sig->Clone ("h");
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

    h = h_jet_trk_pt_as_ref_sig;
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_pt_as_sig[iCent];
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

    h = h_jet_trk_pt_as_sig_syst[iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullSquare);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);


    dlPad->cd (); 
    dlPad->SetLogx ();

    ymin = 0.53;
    ymax = (strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.45 : 1.17);

    h = (TH1D*) h_jet_trk_pt_ns_iaa[iCent]->Clone ("h");
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

    h = h_jet_trk_pt_ns_iaa[iCent];
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

    h = h_jet_trk_pt_ns_iaa_syst[iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullSquare);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);


    drPad->cd (); 
    drPad->SetLogx ();

    ymin = 0.53;
    ymax = (strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.45 : 1.17);

    h = (TH1D*) h_jet_trk_pt_as_iaa[iCent]->Clone ("h");
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

    h = h_jet_trk_pt_as_iaa[iCent];
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

    h = h_jet_trk_pt_as_iaa_syst[iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullSquare);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    c->SaveAs (Form ("%s/Plots/PtCh/JetTagged_HadronYields_%i-%i%%_FCalvsZDC_PtCh_%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], tag)); 
  }



  /*
  for (int iCent = 0; iCent < numZdcCentBins; iCent++) {
    const char* canvasName = Form ("c_jet_trk_pt_ns_iCent%i_syst", iCent);
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    c->cd (); 
    c->SetLogx ();
    //c->SetLogy ();

    float ymin = -20;
    float ymax = 20;

    c->Clear ();

    h = (TH1D*) h_jet_trk_pt_ns[iCent]->Clone ("h");
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

    for (int iVar = 1; iVar < nVar; iVar++) {
      h = (TH1D*) h_jet_trk_pt_ns_syst[iCent][iVar]->Clone ("htemp");
      SaveRelativeErrors (h, h_jet_trk_pt_ns[iCent], true);
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

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    //if (iSys == 0) {
    //  myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    //}
    //else {
      myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
    //}
    myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
    for (int iVar = 1; iVar < nVar; iVar++)
      myLineColorText (0.25, 0.77-iVar*0.040, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 1.0, 0.028);
    

    //c->SaveAs (Form ("%s/Plots/Systematics/TotalJetTaggedYield_%s_nearside_ptch_%s_syst.pdf", workPath.Data (), iSys == 0 ? "pp" : Form ("pPb_%i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), tag));
    c->SaveAs (Form ("%s/Plots/Systematics/TotalJetTaggedYield_%s_nearside_ptch_%s_syst.pdf", workPath.Data (), Form ("pPb_%i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), tag));
    
  }


  for (int iCent = 0; iCent < numZdcCentBins; iCent++) {
    const char* canvasName = Form ("c_jet_trk_pt_as_iCent%i_syst", iCent);
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    c->cd (); 
    c->SetLogx ();
    //c->SetLogy ();

    float ymin = -20;
    float ymax = 20;

    c->Clear ();

    h = (TH1D*) h_jet_trk_pt_as[iCent]->Clone ("h");
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

    for (int iVar = 1; iVar < nVar; iVar++) {
      h = (TH1D*) h_jet_trk_pt_as_syst[iCent][iVar]->Clone ("htemp");
      SaveRelativeErrors (h, h_jet_trk_pt_as[iCent], true);
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

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    //if (iSys == 0) {
    //  myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    //}
    //else {
      myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
    //}
    myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);
    for (int iVar = 1; iVar < nVar; iVar++)
      myLineColorText (0.25, 0.77-iVar*0.040, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 1.0, 0.028);
    

    //c->SaveAs (Form ("%s/Plots/Systematics/TotalJetTaggedYield_%s_awayside_ptch_%s_syst.pdf", workPath.Data (), iSys == 0 ? "pp" : Form ("pPb_%i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), tag));
    c->SaveAs (Form ("%s/Plots/Systematics/TotalJetTaggedYield_%s_awayside_ptch_%s_syst.pdf", workPath.Data (), Form ("pPb_%i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), tag));
  }


  for (int iCent = 0; iCent < numZdcCentBins; iCent++) {
    const char* canvasName = Form ("c_jet_trk_pt_ns_sig_iCent%i_syst", iCent);
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    c->cd (); 
    c->SetLogx ();
    //c->SetLogy ();

    float ymin = -20;
    float ymax = 20;

    c->Clear ();

    h = (TH1D*) h_jet_trk_pt_ns_sig[iCent]->Clone ("h");
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

    for (int iVar = 1; iVar < nVar; iVar++) {
      h = (TH1D*) h_jet_trk_pt_ns_sig_syst[iCent][iVar]->Clone ("htemp");
      SaveRelativeErrors (h, h_jet_trk_pt_ns_sig[iCent], true);
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

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    //if (iSys == 0) {
    //  myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    //}
    //else {
      myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
    //}
    myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
    for (int iVar = 1; iVar < nVar; iVar++)
      myLineColorText (0.25, 0.77-iVar*0.040, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 1.0, 0.028);
    

    //c->SaveAs (Form ("%s/Plots/Systematics/SignalJetTaggedYield_%s_nearside_ptch_%s_syst.pdf", workPath.Data (), iSys == 0 ? "pp" : Form ("pPb_%i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), tag));
    c->SaveAs (Form ("%s/Plots/Systematics/SignalJetTaggedYield_%s_nearside_ptch_%s_syst.pdf", workPath.Data (), Form ("pPb_%i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), tag));
  }


  for (int iCent = 0; iCent < numZdcCentBins; iCent++) {
    const char* canvasName = Form ("c_jet_trk_pt_as_sig_iCent%i_syst", iCent);
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    c->cd (); 
    c->SetLogx ();
    //c->SetLogy ();

    float ymin = -20;
    float ymax = 20;

    c->Clear ();

    h = (TH1D*) h_jet_trk_pt_as_sig[iCent]->Clone ("h");
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

    for (int iVar = 1; iVar < nVar; iVar++) {
      h = (TH1D*) h_jet_trk_pt_as_sig_syst[iCent][iVar]->Clone ("htemp");
      SaveRelativeErrors (h, h_jet_trk_pt_as_sig[iCent], true);
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

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    //if (iSys == 0) {
    //  myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    //}
    //else {
      myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
    //}
    myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);
    for (int iVar = 1; iVar < nVar; iVar++)
      myLineColorText (0.25, 0.77-iVar*0.040, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 1.0, 0.028);
    

    //c->SaveAs (Form ("%s/Plots/Systematics/SignalJetTaggedYield_%s_awayside_ptch_%s_syst.pdf", workPath.Data (), iSys == 0 ? "pp" : Form ("pPb_%i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), tag));
    c->SaveAs (Form ("%s/Plots/Systematics/SignalJetTaggedYield_%s_awayside_ptch_%s_syst.pdf", workPath.Data (), Form ("pPb_%i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), tag));
  }


  for (int iCent = 0; iCent < numZdcCentBins; iCent++) {
    const char* canvasName = Form ("c_jet_trk_pt_ns_iaa_iCent%i_syst", iCent);
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    c->cd (); 
    c->SetLogx ();
    //c->SetLogy ();

    float ymin = -20;
    float ymax = 20;


    h = (TH1D*) h_jet_trk_pt_ns_iaa[iCent]->Clone ("h");
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

    for (int iVar = 1; iVar < nVar; iVar++) {
      h = (TH1D*) h_jet_trk_pt_ns_iaa_syst[iCent][iVar]->Clone ("htemp");
      SaveRelativeErrors (h, h_jet_trk_pt_ns_iaa[iCent], true);
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

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    myText (0.22, 0.81, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
    myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
    for (int iVar = 1; iVar < nVar; iVar++)
      myLineColorText (0.25, 0.73-iVar*0.040, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 1.0, 0.028);
    

    c->SaveAs (Form ("%s/Plots/Systematics/JetTagged_IpA_%i-%i%%_nearside_ptch_%s_syst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], tag));
  }


  for (int iCent = 0; iCent < numZdcCentBins; iCent++) {
    const char* canvasName = Form ("c_jet_trk_pt_as_iaa_iCent%i_syst", iCent);
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    c->cd (); 
    c->SetLogx ();
    //c->SetLogy ();

    float ymin = -20;
    float ymax = 20;


    h = (TH1D*) h_jet_trk_pt_as_iaa[iCent]->Clone ("h");
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

    for (int iVar = 1; iVar < nVar; iVar++) {
      h = (TH1D*) h_jet_trk_pt_as_iaa_syst[iCent][iVar]->Clone ("htemp");
      SaveRelativeErrors (h, h_jet_trk_pt_as_iaa[iCent], true);
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

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    myText (0.22, 0.81, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
    myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);
    for (int iVar = 1; iVar < nVar; iVar++)
      myLineColorText (0.25, 0.73-iVar*0.040, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 1.0, 0.028);
    

    c->SaveAs (Form ("%s/Plots/Systematics/JetTagged_IpA_%i-%i%%_awayside_ptch_%s_syst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], tag));
  }
  */

}


#endif
