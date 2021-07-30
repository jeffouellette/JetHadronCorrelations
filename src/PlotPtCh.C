#ifndef __JetHadronCorrelatorPlotPtCh_C__
#define __JetHadronCorrelatorPlotPtCh_C__

#include <iostream>
#include <math.h>

#include <TColor.h>
#include <TLine.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>

#include <ArrayTemplates.h>
#include <Utilities.h>
#include <MyStyle.h>
#include <MyColors.h>

#include "CentralityDefs.h"
#include "Params.h"
#include "OutTree.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Variations.h"


using namespace JetHadronCorrelations;


TLine* l = new TLine ();
TLatex* tl = new TLatex ();

const std::vector <TString> directions = {"ns", "perp", "as"};
const int nDir = (int) directions.size ();


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

  TH1D**  h_evt_counts_ref      = Get1DArray <TH1D*> (2);
  TH1D**  h_jet_counts_ref      = Get1DArray <TH1D*> (2);
  TH1D**  h_evt_counts_ref_bkg  = Get1DArray <TH1D*> (2);
  TH1D**  h_jet_counts_ref_bkg  = Get1DArray <TH1D*> (2);
  TH1D*** h_evt_counts          = Get2DArray <TH1D*> (2, nZdcCentBins);
  TH1D*** h_jet_counts          = Get2DArray <TH1D*> (2, nZdcCentBins);
  TH1D*** h_evt_counts_bkg      = Get2DArray <TH1D*> (2, nZdcCentBins);
  TH1D*** h_jet_counts_bkg      = Get2DArray <TH1D*> (2, nZdcCentBins);

  TH1D***  h_jet_trk_pt_ref     = Get2DArray <TH1D*> (2, nDir);
  TH1D***  h_jet_trk_pt_ref_bkg = Get2DArray <TH1D*> (2, nDir);
  TH1D**** h_jet_trk_pt         = Get3DArray <TH1D*> (2, nDir, nZdcCentBins);
  TH1D**** h_jet_trk_pt_bkg     = Get3DArray <TH1D*> (2, nDir, nZdcCentBins);

  TH1D***  h_jet_trk_pt_ref_sig = Get2DArray <TH1D*> (2, nDir);
  TH1D**** h_jet_trk_pt_sig     = Get3DArray <TH1D*> (2, nDir, nZdcCentBins);

  TH1D**** h_jet_trk_pt_iaa     = Get3DArray <TH1D*> (2, nDir, nZdcCentBins);

  TGAE**  g_jet_trk_pt_ref_syst     = Get1DArray <TGAE*> (nDir);
  TGAE**  g_jet_trk_pt_ref_bkg_syst = Get1DArray <TGAE*> (nDir);
  TGAE*** g_jet_trk_pt_syst         = Get2DArray <TGAE*> (nDir, nZdcCentBins);
  TGAE*** g_jet_trk_pt_bkg_syst     = Get2DArray <TGAE*> (nDir, nZdcCentBins);

  TGAE**  g_jet_trk_pt_ref_sig_syst = Get1DArray <TGAE*> (nDir);
  TGAE*** g_jet_trk_pt_sig_syst     = Get2DArray <TGAE*> (nDir, nZdcCentBins);

  TGAE*** g_jet_trk_pt_iaa_syst     = Get2DArray <TGAE*> (nDir, nZdcCentBins);


  TH1D***  h_jet_trk_pt_ref_syst      = Get2DArray <TH1D*> (nDir, nVar);
  TH1D***  h_jet_trk_pt_ref_bkg_syst  = Get2DArray <TH1D*> (nDir, nVar);
  TH1D**** h_jet_trk_pt_syst          = Get3DArray <TH1D*> (nDir, nZdcCentBins, nVar);
  TH1D**** h_jet_trk_pt_bkg_syst      = Get3DArray <TH1D*> (nDir, nZdcCentBins, nVar);

  TH1D***  h_jet_trk_pt_ref_sig_syst  = Get2DArray <TH1D*> (nDir, nVar);
  TH1D**** h_jet_trk_pt_sig_syst      = Get3DArray <TH1D*> (nDir, nZdcCentBins, nVar);

  TH1D**** h_jet_trk_pt_iaa_syst      = Get3DArray <TH1D*> (nDir, nZdcCentBins, nVar);


  {
    TString inFileName = inFileTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/PlotPtCh_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    for (int iDType = 0; iDType < 2; iDType++) {

      const TString dType = (iDType == 0 ? "data" : "mc");

      h_evt_counts_ref[iDType] = (TH1D*) inFile->Get (Form ("h_evt_counts_ref_%s_Nominal", dType.Data ()));
      h_jet_counts_ref[iDType] = (TH1D*) inFile->Get (Form ("h_jet_counts_ref_%s_Nominal", dType.Data ()));

      h_evt_counts_ref_bkg[iDType] = (TH1D*) inFile->Get (Form ("h_evt_counts_ref_bkg_%s_Nominal", dType.Data ()));
      h_jet_counts_ref_bkg[iDType] = (TH1D*) inFile->Get (Form ("h_jet_counts_ref_bkg_%s_Nominal", dType.Data ()));

      for (int iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        h_jet_trk_pt_ref[iDType][iDir] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_ref_%s_Nominal", dir.Data (), dType.Data ()));

        h_jet_trk_pt_ref_bkg[iDType][iDir] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_ref_bkg_%s_Nominal", dir.Data (), dType.Data ()));

        h_jet_trk_pt_ref_sig[iDType][iDir] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_ref_sig_%s_Nominal", dir.Data (), dType.Data ()));

      } // end loop over iDir

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        h_evt_counts[iDType][iCent] = (TH1D*) inFile->Get (Form ("h_evt_counts_pPb_iCent%i_%s_Nominal", iCent, dType.Data ()));
        h_jet_counts[iDType][iCent] = (TH1D*) inFile->Get (Form ("h_jet_counts_pPb_iCent%i_%s_Nominal", iCent, dType.Data ()));

        h_evt_counts_bkg[iDType][iCent] = (TH1D*) inFile->Get (Form ("h_evt_counts_pPb_bkg_iCent%i_%s_Nominal", iCent, dType.Data ()));
        h_jet_counts_bkg[iDType][iCent] = (TH1D*) inFile->Get (Form ("h_jet_counts_pPb_bkg_iCent%i_%s_Nominal", iCent, dType.Data ()));

        for (int iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          h_jet_trk_pt[iDType][iDir][iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_pPb_iCent%i_%s_Nominal", dir.Data (), iCent, dType.Data ()));

          h_jet_trk_pt_bkg[iDType][iDir][iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_pPb_bkg_iCent%i_%s_Nominal", dir.Data (), iCent, dType.Data ()));

          h_jet_trk_pt_sig[iDType][iDir][iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_pPb_sig_iCent%i_%s_Nominal", dir.Data (), iCent, dType.Data ()));

          h_jet_trk_pt_iaa[iDType][iDir][iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_iaa_iCent%i_%s_Nominal", dir.Data (), iCent, dType.Data ()));

        } // end loop over iDir

      } // end loop over iCent

    } // end loop over iDType



    for (int iDir = 0; iDir < nDir; iDir++) {

      const TString dir = directions[iDir];

      g_jet_trk_pt_ref_syst[iDir] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_ref_syst", dir.Data ()));

      g_jet_trk_pt_ref_bkg_syst[iDir] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_ref_bkg_syst", dir.Data ()));

      g_jet_trk_pt_ref_sig_syst[iDir] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_ref_sig_syst", dir.Data ()));

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        g_jet_trk_pt_syst[iDir][iCent] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_syst_iCent%i", dir.Data (), iCent));

        g_jet_trk_pt_bkg_syst[iDir][iCent] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_bkg_syst_iCent%i", dir.Data (), iCent));

        g_jet_trk_pt_sig_syst[iDir][iCent] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_sig_syst_iCent%i", dir.Data (), iCent));

        g_jet_trk_pt_iaa_syst[iDir][iCent] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_iaa_syst_iCent%i", dir.Data (), iCent));

      } // end loop over iCent


      for (int iVar = 1; iVar < nVar; iVar++) {

        const TString dType = (dataVariations.count (variations[iVar]) > 0 ? "data" : "mc");

        h_jet_trk_pt_ref_syst[iDir][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_ref_%s_%s", dir.Data (), dType.Data (), variations[iVar].Data ()));

        h_jet_trk_pt_ref_bkg_syst[iDir][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_ref_bkg_%s_%s", dir.Data (), dType.Data (), variations[iVar].Data ()));

        h_jet_trk_pt_ref_sig_syst[iDir][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_ref_sig_%s_%s", dir.Data (), dType.Data (), variations[iVar].Data ()));

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          h_jet_trk_pt_syst[iDir][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_pPb_iCent%i_%s_%s", dir.Data (), iCent, dType.Data (), variations[iVar].Data ()));

          h_jet_trk_pt_bkg_syst[iDir][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_pPb_bkg_iCent%i_%s_%s", dir.Data (), iCent, dType.Data (), variations[iVar].Data ()));

          h_jet_trk_pt_sig_syst[iDir][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_pPb_sig_iCent%i_%s_%s", dir.Data (), iCent, dType.Data (), variations[iVar].Data ()));

          h_jet_trk_pt_iaa_syst[iDir][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_iaa_iCent%i_%s_%s", dir.Data (), iCent, dType.Data (), variations[iVar].Data ()));

        } // end loop over iCent

      } // end loop over iVar

    } // end loop over iDir

  }



  l->SetLineStyle (7);
  l->SetLineWidth (1);
  l->SetLineColor (kBlack);



  for (int iDType = 0; iDType < 2; iDType++) {

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

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
      h = (TH1D*) h_jet_trk_pt_ref[iDType][0]->Clone ("h");
      h->Reset ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]");
      h->GetYaxis ()->SetTitleSize (0.028/fuPad);
      h->GetYaxis ()->SetLabelSize (0.028/fuPad);
      h->GetYaxis ()->SetTitleOffset (3.0*fuPad);

      h->SetLineWidth (1);
      h->SetLineStyle (2);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      myDrawSyst (g_jet_trk_pt_ref_syst[0], myBlue);
      h = h_jet_trk_pt_ref[iDType][0];
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myBlue, kFullCircle, 0.8);
      SaferDelete (&g);

      myDrawSyst (g_jet_trk_pt_ref_bkg_syst[0], myPurple);
      h = h_jet_trk_pt_ref_bkg[iDType][0];
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myPurple, kOpenCircle, 0.8);
      SaferDelete (&g);

      myDrawSyst (g_jet_trk_pt_syst[0][iCent], myRed);
      h = h_jet_trk_pt[iDType][0][iCent];
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myRed, kFullCircle, 0.8);
      SaferDelete (&g);

      myDrawSyst (g_jet_trk_pt_bkg_syst[0][iCent], myGreen);
      h = h_jet_trk_pt_bkg[iDType][0][iCent];
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myGreen, kOpenCircle, 0.8);
      SaferDelete (&g);

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

      h = (TH1D*) h_jet_trk_pt_ref[iDType][2]->Clone ("h");
      h->Reset ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]");
      h->GetYaxis ()->SetTitleSize (0.028/fuPad);
      h->GetYaxis ()->SetLabelSize (0.028/fuPad);
      h->GetYaxis ()->SetTitleOffset (3.0*fuPad);

      h->SetLineWidth (1);
      h->SetLineStyle (2);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      myDrawSyst (g_jet_trk_pt_ref_syst[2], myBlue);
      h = h_jet_trk_pt_ref[iDType][2];
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myBlue, kFullCircle, 0.8);
      SaferDelete (&g);

      myDrawSyst (g_jet_trk_pt_ref_bkg_syst[2], myPurple);
      h = h_jet_trk_pt_ref_bkg[iDType][2];
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myPurple, kOpenCircle, 0.8);
      SaferDelete (&g);

      myDrawSyst (g_jet_trk_pt_syst[2][iCent], myRed);
      h = h_jet_trk_pt[iDType][2][iCent];
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myRed, kFullCircle, 0.8);
      SaferDelete (&g);

      myDrawSyst (g_jet_trk_pt_bkg_syst[2][iCent], myGreen);
      h = h_jet_trk_pt_bkg[iDType][2][iCent];
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myGreen, kOpenCircle, 0.8);
      SaferDelete (&g);

      if (iDType == 0) myText (0.58, 0.78, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
      else {
        myText (0.42, 0.78, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.022/fuPad);
        myText (0.42, 0.72, kBlack, Form ("Pythia8 JZ%i (+ #it{p}+Pb Overlay)", TString (tag) == TString ("30GeVJets") ? 1 : 2), 0.022/fuPad);
      }
      myBoxText2 (0.10, 0.24, myRed, kFullCircle, "#it{p}+Pb total", 0.8, 0.020/fuPad, true);
      myBoxText2 (0.10, 0.18, myBlue, kFullCircle, "#it{pp} total", 0.8, 0.020/fuPad, true);
      myBoxText2 (0.10, 0.12, myGreen, kOpenCircle, "#it{p}+Pb bkgd.", 0.8, 0.020/fuPad);
      myBoxText2 (0.10, 0.06, myPurple, kOpenCircle, "#it{pp} bkgd.", 0.8, 0.020/fuPad);

      tl->DrawLatexNDC (0.5*(1-0.03/0.48), 0.94, "#Delta#phi > 7#pi/8 (away-side)");


      clPad->cd (); 
      clPad->SetLogx ();
      clPad->SetLogy ();

      ymin = 8e-6;
      ymax = 3e1;

      h = (TH1D*) h_jet_trk_pt_ref_sig[iDType][0]->Clone ("h");
      h->Reset ();
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
      h->GetYaxis ()->SetTitleSize (0.028/fdPad);
      h->GetYaxis ()->SetLabelSize (0.028/fdPad);
      h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
      h->GetYaxis ()->CenterTitle ();

      h->SetLineWidth (1);
      h->SetLineStyle (2);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      myDrawSyst (g_jet_trk_pt_ref_sig_syst[0], myBlue);
      h = h_jet_trk_pt_ref_sig[iDType][0];
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myBlue, kFullCircle, 0.8);
      SaferDelete (&g);

      myDrawSyst (g_jet_trk_pt_sig_syst[0][iCent], myRed);
      h = h_jet_trk_pt_sig[iDType][0][iCent];
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myRed, kFullCircle, 0.8);
      SaferDelete (&g);

      myText (0.24, 0.26, kBlack, "Jet-hadron correlations", 0.020/fcPad);
      myText (0.24, 0.18, kBlack, Form ("#it{p}_{T}^{jet} > %s", GetJetPtStr (tag).Data ()), 0.020/fcPad);
      myText (0.24, 0.10, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fcPad);


      crPad->cd (); 
      crPad->SetLogx ();
      crPad->SetLogy ();

      ymin = 8e-6;
      ymax = 3e1;

      h = (TH1D*) h_jet_trk_pt_ref_sig[iDType][2]->Clone ("h");
      h->Reset ();
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
      h->GetYaxis ()->SetTitleSize (0.028/fdPad);
      h->GetYaxis ()->SetLabelSize (0.028/fdPad);
      h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
      h->GetYaxis ()->CenterTitle ();

      h->SetLineWidth (1);
      h->SetLineStyle (2);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      myDrawSyst (g_jet_trk_pt_ref_sig_syst[2], myBlue);
      h = h_jet_trk_pt_ref_sig[iDType][2];
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myBlue, kFullCircle, 0.8);
      SaferDelete (&g);

      myDrawSyst (g_jet_trk_pt_sig_syst[2][iCent], myRed);
      h = h_jet_trk_pt_sig[iDType][2][iCent];
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myRed, kFullCircle, 0.8);
      SaferDelete (&g);


      dlPad->cd (); 
      dlPad->SetLogx ();

      ymin = 0.53;
      ymax = (strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.45 : 1.17);

      h = (TH1D*) h_jet_trk_pt_iaa[iDType][0][iCent]->Clone ("h");
      h->Reset ();
      for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      h->GetXaxis ()->SetTitleSize (0.028/fdPad);
      h->GetXaxis ()->SetLabelSize (0.028/fdPad);
      h->GetXaxis ()->SetTitleOffset (3.7*fdPad);
      h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
      h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}Pb} = #it{p}+Pb / #it{pp}");
      h->GetYaxis ()->SetTitleSize (0.028/fdPad);
      h->GetYaxis ()->SetLabelSize (0.028/fdPad);
      h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
      h->GetYaxis ()->CenterTitle ();

      h->SetLineWidth (1);
      h->SetLineStyle (2);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      myDrawSyst (g_jet_trk_pt_iaa_syst[0][iCent], myRed);
      h = h_jet_trk_pt_iaa[iDType][0][iCent];
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myRed, kFullCircle, 0.8);
      SaferDelete (&g);


      drPad->cd (); 
      drPad->SetLogx ();

      ymin = 0.53;
      ymax = (strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.45 : 1.17);

      h = (TH1D*) h_jet_trk_pt_iaa[iDType][2][iCent]->Clone ("h");
      h->Reset ();
      for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      h->GetXaxis ()->SetTitleSize (0.028/fdPad);
      h->GetXaxis ()->SetLabelSize (0.028/fdPad);
      h->GetXaxis ()->SetTitleOffset (3.7*fdPad);
      h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
      h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}Pb} = #it{p}+Pb / #it{pp}");
      h->GetYaxis ()->SetTitleSize (0.028/fdPad);
      h->GetYaxis ()->SetLabelSize (0.028/fdPad);
      h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
      h->GetYaxis ()->CenterTitle ();

      h->SetLineWidth (1);
      h->SetLineStyle (2);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      myDrawSyst (g_jet_trk_pt_iaa_syst[2][iCent], myRed);
      h = h_jet_trk_pt_iaa[iDType][2][iCent];
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myRed, kFullCircle, 0.8);
      SaferDelete (&g);

      c->SaveAs (Form ("%s/Plots/PtCh/JetTagged_HadronYields_%i-%i%%_comparison_PtCh_%s%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], tag, iDType == 1 ? "_mc" : "")); 

    } // end loop over iCent

  } // end loop over iDType


/*
  for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
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

    int iVar = 0;
    while (iVar < nVar && strcmp (variations[iVar], "FcalCentVar") != 0) iVar++;
    if (iVar == nVar) {
      std::cout << "Cannot find FCal centrality binned result??? Please check!" << std::endl;
      return;
    }

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    ulPad->cd (); 
    ulPad->SetLogx ();
    ulPad->SetLogy ();

    float ymin = 8e-6;
    float ymax = 3e1;

    h = (TH1D*) h_jet_trk_pt_ref[0]->Clone ("h");
    h->Reset ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]");
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fuPad);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_trk_pt_ref[0][0];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_ref_bkg[0][0];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kOpenCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt[0][0][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_bkg[0][0][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myGreen, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_syst[0][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_bkg_syst[0][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myOrange, kFullSquare, 0.8);
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

    h = (TH1D*) h_jet_trk_pt_ref[0][2]->Clone ("h");
    h->Reset ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]");
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fuPad);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_trk_pt_ref[0][2];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_ref_bkg[0][2];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kOpenCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt[0][2][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_bkg[0][2][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myGreen, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_syst[2][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_bkg_syst[2][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myOrange, kFullSquare, 0.8);
    SaferDelete (&g);

    myText (0.58, 0.78, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
    myLineText2 (0.10, 0.36, myRed, kFullCircle, Form ("#it{p}+Pb #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.8, 0.020/fuPad, true);
    myLineText2 (0.10, 0.30, myGreen, kFullCircle, Form ("#it{p}+Pb bkgd., #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.8, 0.020/fuPad, true);
    myLineText2 (0.10, 0.24, myBlue, kFullSquare, Form ("#it{p}+Pb #bf{FCal %i-%i%%}", fcalCentPercs[iCent+1], fcalCentPercs[iCent]), 0.8, 0.020/fuPad, true);
    myLineText2 (0.10, 0.18, myOrange, kFullSquare, Form ("#it{p}+Pb bkgd., #bf{FCal %i-%i%%}", fcalCentPercs[iCent+1], fcalCentPercs[iCent]), 0.8, 0.020/fuPad, true);
    myLineText2 (0.10, 0.12, kBlack, kFullCircle, "#it{pp} total", 0.8, 0.020/fuPad);
    myLineText2 (0.10, 0.06, kBlack, kOpenCircle, "#it{pp} bkgd.", 0.8, 0.020/fuPad);

    tl->DrawLatexNDC (0.5*(1-0.03/0.48), 0.94, "#Delta#phi > 7#pi/8 (away-side)");


    clPad->cd (); 
    clPad->SetLogx ();
    clPad->SetLogy ();

    ymin = 8e-6;
    ymax = 3e1;

    h = (TH1D*) h_jet_trk_pt_ref_sig[0][0]->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_trk_pt_ref_sig[0][0];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_sig[0][0][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_sig_syst[0][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);

    myText (0.24, 0.26, kBlack, "Jet-hadron correlations", 0.020/fcPad);
    myText (0.24, 0.18, kBlack, Form ("#it{p}_{T}^{jet} > %s", GetJetPtStr (tag).Data ()), 0.020/fcPad);
    myText (0.24, 0.10, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fcPad);


    crPad->cd (); 
    crPad->SetLogx ();
    crPad->SetLogy ();

    ymin = 8e-6;
    ymax = 3e1;

    h = (TH1D*) h_jet_trk_pt_ref_sig[0][2]->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_trk_pt_ref_sig[0][2];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_sig[0][2][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_sig_syst[2][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);


    dlPad->cd (); 
    dlPad->SetLogx ();

    ymin = 0.53;
    ymax = (strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.45 : 1.17);

    h = (TH1D*) h_jet_trk_pt_iaa[0][0][iCent]->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.7*fdPad);
    h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}Pb} = #it{p}+Pb / #it{pp}");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_trk_pt_iaa[0][0][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_iaa_syst[0][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);


    drPad->cd (); 
    drPad->SetLogx ();

    ymin = 0.53;
    ymax = (strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.45 : 1.17);

    h = (TH1D*) h_jet_trk_pt_iaa[0][2][iCent]->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.7*fdPad);
    h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}Pb} = #it{p}+Pb / #it{pp}");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_trk_pt_iaa[0][2][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_iaa_syst[2][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);

    c->SaveAs (Form ("%s/Plots/PtCh/JetTagged_HadronYields_%i-%i%%_FCalvsZDC_PtCh_%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], tag)); 
  }



  for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
    const char* canvasName = Form ("c_jet_trk_pt_NoFcalMixCatVar_iCent%i", iCent);
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

    int iVar = 0;
    while (iVar < nVar && strcmp (variations[iVar], "NoFcalMixCatVar") != 0) iVar++;
    if (iVar == nVar) {
      std::cout << "Cannot find FCal matching result??? Please check!" << std::endl;
      return;
    }

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    ulPad->cd (); 
    ulPad->SetLogx ();
    ulPad->SetLogy ();

    float ymin = 8e-6;
    float ymax = 3e1;

    h = (TH1D*) h_jet_trk_pt_ref[0][0]->Clone ("h");
    h->Reset ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]");
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fuPad);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_trk_pt_ref[0][0];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kOpenCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_ref_bkg[0][0];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myGreen, kOpenCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_ref_bkg_syst[0][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myOrange, kOpenSquare, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt[0][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_bkg[0][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_bkg_syst[0][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
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

    h = (TH1D*) h_jet_trk_pt_ref[0][2]->Clone ("h");
    h->Reset ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]");
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fuPad);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_trk_pt_ref[0][2];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kOpenCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_ref_bkg[0][2];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myGreen, kOpenCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_ref_bkg_syst[2][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myOrange, kOpenSquare, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt[0][2][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_bkg[0][2][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_bkg_syst[2][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);

    myText (0.58, 0.78, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
    myLineText2 (0.10, 0.36, kBlack, kFullCircle, "#it{p}+Pb total", 0.8, 0.020/fuPad, true);
    myLineText2 (0.10, 0.30, myRed, kFullCircle, "#it{p}+Pb bkgd., #bf{w/} FCal matching", 0.8, 0.020/fuPad, true);
    myLineText2 (0.10, 0.24, myBlue, kFullSquare, "#it{p}+Pb bkgd., #bf{w/out} FCal matching", 0.8, 0.020/fuPad, true);
    myLineText2 (0.10, 0.18, kBlack, kOpenCircle, "#it{pp} total", 0.8, 0.020/fuPad);
    myLineText2 (0.10, 0.12, myGreen, kOpenCircle, "#it{pp} bkgd., #bf{w/} FCal matching", 0.8, 0.020/fuPad);
    myLineText2 (0.10, 0.06, myOrange, kOpenSquare, "#it{pp} bkgd., #bf{w/out} FCal matching", 0.8, 0.020/fuPad);

    tl->DrawLatexNDC (0.5*(1-0.03/0.48), 0.94, "#Delta#phi > 7#pi/8 (away-side)");


    clPad->cd (); 
    clPad->SetLogx ();
    clPad->SetLogy ();

    ymin = 8e-6;
    ymax = 3e1;

    h = (TH1D*) h_jet_trk_pt_ref_sig[0][0]->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_trk_pt_ref_sig[0][0];
    g = make_graph (h);
    ResetXErrors (g); myDraw (g, myGreen, kOpenCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_ref_sig_syst[0][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myOrange, kOpenSquare, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_sig[0][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_sig_syst[0][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);

    myText (0.24, 0.26, kBlack, "Jet-hadron correlations", 0.020/fcPad);
    myText (0.24, 0.18, kBlack, Form ("#it{p}_{T}^{jet} > %s", GetJetPtStr (tag).Data ()), 0.020/fcPad);
    myText (0.24, 0.10, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fcPad);


    crPad->cd (); 
    crPad->SetLogx ();
    crPad->SetLogy ();

    ymin = 8e-6;
    ymax = 3e1;

    h = (TH1D*) h_jet_trk_pt_ref_sig[0][2]->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_trk_pt_ref_sig[0][2];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myGreen, kOpenCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_ref_sig_syst[2][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myOrange, kOpenSquare, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_sig[0][2][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_sig_syst[2][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);

    myLineText2 (0.10, 0.34, myRed, kFullCircle, "#it{p}+Pb #bf{w/} FCal matching", 0.8, 0.020/fcPad, true);
    myLineText2 (0.10, 0.26, myBlue, kFullSquare, "#it{p}+Pb #bf{w/out} FCal matching", 0.8, 0.020/fcPad, true);
    myLineText2 (0.10, 0.18, myGreen, kOpenCircle, "#it{pp} #bf{w/} FCal matching", 0.8, 0.020/fcPad);
    myLineText2 (0.10, 0.10, myOrange, kOpenSquare, "#it{pp} #bf{w/out} FCal matching", 0.8, 0.020/fcPad);


    dlPad->cd (); 
    dlPad->SetLogx ();

    ymin = 0.53;
    ymax = (strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.45 : 1.17);

    h = (TH1D*) h_jet_trk_pt_iaa[0][0][iCent]->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.7*fdPad);
    h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}Pb} = #it{p}+Pb / #it{pp}");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_trk_pt_iaa[0][0][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_iaa_syst[0][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);


    drPad->cd (); 
    drPad->SetLogx ();

    ymin = 0.53;
    ymax = (strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.45 : 1.17);

    h = (TH1D*) h_jet_trk_pt_iaa[0][2][iCent]->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.7*fdPad);
    h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}Pb} = #it{p}+Pb / #it{pp}");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_trk_pt_iaa[0][2][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jet_trk_pt_iaa_syst[2][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);

    myLineText2 (0.10, 0.48, myRed, kFullCircle, "#bf{w/} FCal matching", 0.8, 0.020/fdPad, true);
    myLineText2 (0.10, 0.38, myBlue, kFullSquare, "#bf{w/out} FCal matching", 0.8, 0.020/fdPad, true);

    c->SaveAs (Form ("%s/Plots/PtCh/JetTagged_HadronYields_%i-%i%%_FCalMixCatVar_PtCh_%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], tag)); 
  }*/


 
  for (int iDir : {0, 2}) {

    const char* canvasName = Form ("c_jet_trk_pt_%s_sig2bkg", directions[iDir].Data ());
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    c->cd (); 
    c->SetLogx ();
    c->SetLogy ();

    float ymin = 1e-1;
    float ymax = 1e6;

    c->Clear ();

    h = (TH1D*) h_jet_trk_pt_ref_sig[0][iDir]->Clone ("h");
    h->Divide (h_jet_trk_pt_ref_bkg[0][iDir]);
    g = make_graph (h);
    ResetXErrors (g);

    g->GetXaxis ()->SetMoreLogLabels ();
    g->GetYaxis ()->SetRangeUser (ymin, ymax);
    g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    //g->GetXaxis ()->SetTitleSize (0.028);
    //g->GetXaxis ()->SetLabelSize (0.028);
    g->GetYaxis ()->SetTitle ("Sig. / Bkgd.");
    //g->GetYaxis ()->SetTitleSize (0.028);
    //g->GetYaxis ()->SetLabelSize (0.028);

    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);

    ((TGAE*) g->Clone ())->Draw ("AP");
    SaferDelete (&g);
    SaferDelete (&h);

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      h = (TH1D*) h_jet_trk_pt_sig[0][iDir][iCent]->Clone ("htemp");
      h->Divide (h_jet_trk_pt_bkg[0][iDir][iCent]);
      g = make_graph (h);
      SaferDelete (&h);
      ResetXErrors (g);
      myDraw (g, colors[iCent], kFullCircle, 0.8);
      SaferDelete (&g);
    }

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    myText (0.22, 0.81, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    if (directions[iDir] == "ns")
      myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
    else if (directions[iDir] == "as")
      myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);
    myLineText2 (0.25, 0.72, kBlack, kFullCircle, "#bf{#it{pp}}", 0.8, 0.032, true);
    for (int iCent = 0; iCent < nZdcCentBins; iCent++)
      myLineText2 (0.25, 0.68-iCent*0.04, colors[iCent], kFullCircle, Form ("#bf{#it{p}+Pb, %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.8, 0.032, true);

    c->SaveAs (Form ("%s/Plots/PtCh/SigToBkgd_%s_%s_ptch.pdf", workPath.Data (), tag, directions[iDir] == "ns" ? "nearside" : "awayside"));
  }




  for (int iDType : {0, 1}) {

    for (int iDir : {0, 2}) {

      {
        const char* canvasName = Form ("c_jet_trk_pt_%s_pp_%s_syst", directions[iDir].Data (), iDType == 0 ? "data" : "mc");
        TCanvas* c = new TCanvas (canvasName, "", 800, 800);

        TH1D* h = nullptr, *h_tot = nullptr;
        TGAE* g = nullptr;

        c->cd (); 
        c->SetLogx ();

        float ymin = (iDType == 0 ? -20 : -10);
        float ymax = (iDType == 0 ?  20 :  10);

        h = (TH1D*) h_jet_trk_pt_ref[iDType][iDir]->Clone ("h");
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

        h_tot = (TH1D*) h_jet_trk_pt_ref[iDType][iDir]->Clone ("h_total_unc");
        h_tot->Reset ();

        for (int iVar = 1; iVar < nVar; iVar++) {
          if ((iDType == 0  && dataVariations.count (variations[iVar]) == 0) || (iDType == 1 && mcVariations.count (variations[iVar]) == 0))
            continue;

          h = (TH1D*) h_jet_trk_pt_ref_syst[iDir][iVar]->Clone ("htemp");
          SaveRelativeErrors (h, h_jet_trk_pt_ref[iDType][iDir], true);
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            h->SetBinContent (iX, 100*h->GetBinContent (iX) - 100);
            h_tot->SetBinContent (iX, std::sqrt (std::pow (h_tot->GetBinContent (iX), 2) + std::pow (h->GetBinContent (iX), 2)));
          }
          g = make_graph (h);
          SaferDelete (&h);
          ResetXErrors (g);
          ResetTGAEErrors (g);

          g->SetLineColor (varStyles[variations[iVar]].first);
          g->SetLineStyle (varStyles[variations[iVar]].second);
          g->SetLineWidth (3);
          ((TGAE*) g->Clone ())->Draw ("L");
          SaferDelete (&g);
        }

        g = make_graph (h_tot);
        SaferDelete (&h_tot);
        ResetXErrors (g);
        ResetTGAEErrors (g);

        g->SetLineColor (kBlack);
        g->SetLineStyle (1);
        g->SetLineWidth (3);
        ((TGAE*) g->Clone ())->Draw ("L");
        SaferDelete (&g);

        myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
        myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
        if (directions[iDir] == "ns")
          myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
        else if (directions[iDir] == "as")
          myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);

        int count = 0;
        for (int iVar = 1; iVar < nVar; iVar++) {
          if ((iDType == 0  && dataVariations.count (variations[iVar]) > 0) || (iDType == 1 && mcVariations.count (variations[iVar]) > 0)) {
            myLineColorText (0.25+0.40*(count>=12), 0.45-(count%12)*0.020, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 0.5, 0.014);
            count++;
          }
        }

        c->SaveAs (Form ("%s/Plots/Systematics/TotalJetTaggedYield_pp_%s_ptch_%s_%ssyst.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : "awayside", tag, iDType == 0 ? "":"mcBased_"));
      }


      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        const char* canvasName = Form ("c_jet_trk_pt_%s_%s_pPb_iCent%i_syst", directions[iDir].Data (), iDType == 0 ? "data" : "mc", iCent);
        TCanvas* c = new TCanvas (canvasName, "", 800, 800);

        TH1D* h = nullptr, *h_tot = nullptr;
        TGAE* g = nullptr;

        c->cd (); 
        c->SetLogx ();

        float ymin = (iDType == 0 ? -20 : -10);
        float ymax = (iDType == 0 ?  20 :  10);

        h = (TH1D*) h_jet_trk_pt[iDType][iDir][iCent]->Clone ("h");
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

        h_tot = (TH1D*) h_jet_trk_pt[iDType][iDir][iCent]->Clone ("h_total_unc");
        h_tot->Reset ();

        for (int iVar = 1; iVar < nVar; iVar++) {
          if ((iDType == 0  && dataVariations.count (variations[iVar]) == 0) || (iDType == 1 && mcVariations.count (variations[iVar]) == 0))
            continue;

          h = (TH1D*) h_jet_trk_pt_syst[iDir][iCent][iVar]->Clone ("htemp");
          SaveRelativeErrors (h, h_jet_trk_pt[iDType][iDir][iCent], true);
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            h->SetBinContent (iX, 100*h->GetBinContent (iX) - 100);
            h_tot->SetBinContent (iX, std::sqrt (std::pow (h_tot->GetBinContent (iX), 2) + std::pow (h->GetBinContent (iX), 2)));
          }
          g = make_graph (h);
          SaferDelete (&h);
          ResetXErrors (g);
          ResetTGAEErrors (g);

          g->SetLineColor (varStyles[variations[iVar]].first);
          g->SetLineStyle (varStyles[variations[iVar]].second);
          g->SetLineWidth (3);
          ((TGAE*) g->Clone ())->Draw ("L");
          SaferDelete (&g);
        }

        g = make_graph (h_tot);
        SaferDelete (&h_tot);
        ResetXErrors (g);
        ResetTGAEErrors (g);

        g->SetLineColor (kBlack);
        g->SetLineStyle (1);
        g->SetLineWidth (3);
        ((TGAE*) g->Clone ())->Draw ("L");
        SaferDelete (&g);

        myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
        myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
        if (directions[iDir] == "ns")
          myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
        else if (directions[iDir] == "as")
          myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);

        int count = 0;
        for (int iVar = 1; iVar < nVar; iVar++) {
          if ((iDType == 0  && dataVariations.count (variations[iVar]) > 0) || (iDType == 1 && mcVariations.count (variations[iVar]) > 0)) {
            myLineColorText (0.25+0.40*(count>=12), 0.45-(count%12)*0.020, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 0.5, 0.014);
            count++;
          }
        }
        
        c->SaveAs (Form ("%s/Plots/Systematics/TotalJetTaggedYield_pPb_%i-%i%%_%s_ptch_%s_%ssyst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", tag, iDType == 0 ? "":"mcBased_"));
      } // end loop over iCent



      {
        const char* canvasName = Form ("c_jet_trk_pt_%s_sig_%s_pp_syst", directions[iDir].Data (), iDType == 0 ? "data" : "mc");
        TCanvas* c = new TCanvas (canvasName, "", 800, 800);

        TH1D* h = nullptr, *h_tot = nullptr;
        TGAE* g = nullptr;

        c->cd (); 
        c->SetLogx ();

        float ymin = (iDType == 0 ? -20 : -10);
        float ymax = (iDType == 0 ?  20 :  10);

        h = (TH1D*) h_jet_trk_pt_ref_sig[iDType][iDir]->Clone ("h");
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

        h_tot = (TH1D*) h_jet_trk_pt_ref_sig[iDType][iDir]->Clone ("h_total_unc");
        h_tot->Reset ();

        for (int iVar = 1; iVar < nVar; iVar++) {
          if ((iDType == 0  && dataVariations.count (variations[iVar]) == 0) || (iDType == 1 && mcVariations.count (variations[iVar]) == 0))
            continue;

          h = (TH1D*) h_jet_trk_pt_ref_sig_syst[iDir][iVar]->Clone ("htemp");
          SaveRelativeErrors (h, h_jet_trk_pt_ref_sig[iDType][iDir], true);
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            h->SetBinContent (iX, 100*h->GetBinContent (iX) - 100);
            h_tot->SetBinContent (iX, std::sqrt (std::pow (h_tot->GetBinContent (iX), 2) + std::pow (h->GetBinContent (iX), 2)));
          }
          g = make_graph (h);
          SaferDelete (&h);
          ResetXErrors (g);
          ResetTGAEErrors (g);

          g->SetLineColor (varStyles[variations[iVar]].first);
          g->SetLineStyle (varStyles[variations[iVar]].second);
          g->SetLineWidth (3);
          ((TGAE*) g->Clone ())->Draw ("L");
          SaferDelete (&g);
        }

        g = make_graph (h_tot);
        SaferDelete (&h_tot);
        ResetXErrors (g);
        ResetTGAEErrors (g);

        g->SetLineColor (kBlack);
        g->SetLineStyle (1);
        g->SetLineWidth (3);
        ((TGAE*) g->Clone ())->Draw ("L");
        SaferDelete (&g);

        myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
        myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
        if (directions[iDir] == "ns")
          myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
        else if (directions[iDir] == "as")
          myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);

        int count = 0;
        for (int iVar = 1; iVar < nVar; iVar++) {
          if ((iDType == 0  && dataVariations.count (variations[iVar]) > 0) || (iDType == 1 && mcVariations.count (variations[iVar]) > 0)) {
            myLineColorText (0.25+0.40*(count>=12), 0.45-(count%12)*0.020, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 0.5, 0.014);
            count++;
          }
        }

        c->SaveAs (Form ("%s/Plots/Systematics/SignalJetTaggedYield_pp_%s_ptch_%s_%ssyst.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : "awayside", tag, iDType == 0 ? "":"mcBased_"));
      }


      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        const char* canvasName = Form ("c_jet_trk_pt_%s_sig_%s_pPb_iCent%i_syst", directions[iDir].Data (), iDType == 0 ? "data" : "mc", iCent);
        TCanvas* c = new TCanvas (canvasName, "", 800, 800);

        TH1D* h = nullptr, *h_tot = nullptr;
        TGAE* g = nullptr;

        c->cd (); 
        c->SetLogx ();
        //c->SetLogy ();

        float ymin = (iDType == 0 ? -20 : -10);
        float ymax = (iDType == 0 ?  20 :  10);

        h = (TH1D*) h_jet_trk_pt_sig[iDType][iDir][iCent]->Clone ("h");
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

        h_tot = (TH1D*) h_jet_trk_pt_sig[iDType][iDir][iCent]->Clone ("h_total_unc");
        h_tot->Reset ();

        for (int iVar = 1; iVar < nVar; iVar++) {
          if ((iDType == 0  && dataVariations.count (variations[iVar]) == 0) || (iDType == 1 && mcVariations.count (variations[iVar]) == 0))
            continue;

          h = (TH1D*) h_jet_trk_pt_sig_syst[iDir][iCent][iVar]->Clone ("htemp");
          SaveRelativeErrors (h, h_jet_trk_pt_sig[iDType][iDir][iCent], true);
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            h->SetBinContent (iX, 100*h->GetBinContent (iX) - 100);
            h_tot->SetBinContent (iX, std::sqrt (std::pow (h_tot->GetBinContent (iX), 2) + std::pow (h->GetBinContent (iX), 2)));
          }
          g = make_graph (h);
          SaferDelete (&h);
          ResetXErrors (g);
          ResetTGAEErrors (g);

          g->SetLineColor (varStyles[variations[iVar]].first);
          g->SetLineStyle (varStyles[variations[iVar]].second);
          g->SetLineWidth (3);
          ((TGAE*) g->Clone ())->Draw ("L");
          SaferDelete (&g);
        }

        g = make_graph (h_tot);
        SaferDelete (&h_tot);
        ResetXErrors (g);
        ResetTGAEErrors (g);

        g->SetLineColor (kBlack);
        g->SetLineStyle (1);
        g->SetLineWidth (3);
        ((TGAE*) g->Clone ())->Draw ("L");
        SaferDelete (&g);

        myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
        myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
        if (directions[iDir] == "ns")
          myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
        else if (directions[iDir] == "as")
          myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);

        int count = 0;
        for (int iVar = 1; iVar < nVar; iVar++) {
          if ((iDType == 0  && dataVariations.count (variations[iVar]) > 0) || (iDType == 1 && mcVariations.count (variations[iVar]) > 0)) {
            myLineColorText (0.25+0.40*(count>=12), 0.45-(count%12)*0.020, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 0.5, 0.014);
            count++;
          }
        }

        c->SaveAs (Form ("%s/Plots/Systematics/SignalJetTaggedYield_pPb_%i-%i%%_%s_ptch_%s_%ssyst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", tag, iDType == 0 ? "":"mcBased_"));
      } // end loop over iCent



      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        const char* canvasName = Form ("c_jet_trk_pt_%s_iaa_%s_pPb_iCent%i_syst", directions[iDir].Data (), iDType == 0 ? "data" : "mc", iCent);
        TCanvas* c = new TCanvas (canvasName, "", 800, 800);

        TH1D* h = nullptr, *h_tot = nullptr;
        TGAE* g = nullptr;

        c->cd (); 
        c->SetLogx ();
        //c->SetLogy ();

        float ymin = (iDType == 0 ? -20 : -10);
        float ymax = (iDType == 0 ?  20 :  10);

        h = (TH1D*) h_jet_trk_pt_iaa[iDType][iDir][iCent]->Clone ("h");
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

        h_tot = (TH1D*) h_jet_trk_pt_iaa[iDType][iDir][iCent]->Clone ("h_total_unc");
        h_tot->Reset ();

        for (int iVar = 1; iVar < nVar; iVar++) {
          if ((iDType == 0  && dataVariations.count (variations[iVar]) == 0) || (iDType == 1 && mcVariations.count (variations[iVar]) == 0))
            continue;

          h = (TH1D*) h_jet_trk_pt_iaa_syst[iDir][iCent][iVar]->Clone ("htemp");
          SaveRelativeErrors (h, h_jet_trk_pt_iaa[iDType][iDir][iCent], true);
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            h->SetBinContent (iX, 100*h->GetBinContent (iX) - 100);
            h_tot->SetBinContent (iX, std::sqrt (std::pow (h_tot->GetBinContent (iX), 2) + std::pow (h->GetBinContent (iX), 2)));
          }
          g = make_graph (h);
          SaferDelete (&h);
          ResetXErrors (g);
          ResetTGAEErrors (g);

          g->SetLineColor (varStyles[variations[iVar]].first);
          g->SetLineStyle (varStyles[variations[iVar]].second);
          g->SetLineWidth (3);
          ((TGAE*) g->Clone ())->Draw ("L");
          SaferDelete (&g);
        }

        g = make_graph (h_tot);
        SaferDelete (&h_tot);
        ResetXErrors (g);
        ResetTGAEErrors (g);

        g->SetLineColor (kBlack);
        g->SetLineStyle (1);
        g->SetLineWidth (3);
        ((TGAE*) g->Clone ())->Draw ("L");
        SaferDelete (&g);

        myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
        myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
        myText (0.22, 0.81, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
        if (directions[iDir] == "ns")
          myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
        else if (directions[iDir] == "as")
          myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);

        int count = 0;
        for (int iVar = 1; iVar < nVar; iVar++) {
          if ((iDType == 0  && dataVariations.count (variations[iVar]) > 0) || (iDType == 1 && mcVariations.count (variations[iVar]) > 0)) {
            myLineColorText (0.25+0.40*(count>=12), 0.45-(count%12)*0.020, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 0.5, 0.014);
            count++;
          }
        }

        c->SaveAs (Form ("%s/Plots/Systematics/JetTagged_IpPb_%i-%i%%_%s_ptch_%s_%ssyst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", tag, iDType == 0 ? "":"mcBased_"));
      } // end loop over iCent

    } // end loop over iDir

  } // end loop over iDType

}


#endif
