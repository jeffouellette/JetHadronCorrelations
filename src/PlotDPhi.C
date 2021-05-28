#ifndef __JetHadronCorrelatorPlotDPhi_C__
#define __JetHadronCorrelatorPlotDPhi_C__

#include <iostream>
#include <vector>
#include <map>
#include <math.h>

#include <TColor.h>
#include <TLine.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TString.h>
#include <TLorentzVector.h>

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


void PlotDPhi (const char* tag, const char* inFileTag) {

  TFile* inFile = nullptr;

  TH1D*  h_evt_counts_ref = nullptr;
  TH1D*  h_jet_counts_ref = nullptr;
  TH1D*  h_evt_counts_ref_bkg = nullptr;
  TH1D*  h_jet_counts_ref_bkg = nullptr;
  TH1D** h_evt_counts = new TH1D*[numZdcCentBins];
  TH1D** h_jet_counts = new TH1D*[numZdcCentBins];
  TH1D** h_evt_counts_bkg = new TH1D*[numZdcCentBins];
  TH1D** h_jet_counts_bkg = new TH1D*[numZdcCentBins];

  TH1D**  h_jet_trk_dphi_ref = Get1DArray <TH1D*> (nPtChSelections);
  TH2D**  h2_jet_trk_dphi_cov_ref = Get1DArray <TH2D*> (nPtChSelections);
  TH1D**  h_jet_trk_dphi_ref_bkg = Get1DArray <TH1D*> (nPtChSelections);
  TH2D**  h2_jet_trk_dphi_cov_ref_bkg = Get1DArray <TH2D*> (nPtChSelections);
  TH1D*** h_jet_trk_dphi = Get2DArray <TH1D*> (numZdcCentBins, nPtChSelections);
  TH2D*** h2_jet_trk_dphi_cov = Get2DArray <TH2D*> (numZdcCentBins, nPtChSelections);
  TH1D*** h_jet_trk_dphi_bkg = Get2DArray <TH1D*> (numZdcCentBins, nPtChSelections);
  TH2D*** h2_jet_trk_dphi_cov_bkg = Get2DArray <TH2D*> (numZdcCentBins, nPtChSelections);

  TH1D**  h_jet_trk_dphi_ref_sig = Get1DArray <TH1D*> (nPtChSelections);
  TH1D*** h_jet_trk_dphi_sig = Get2DArray <TH1D*> (numZdcCentBins, nPtChSelections);

  TH1D*** h_jet_trk_dphi_iaa = Get2DArray <TH1D*> (numZdcCentBins, nPtChSelections);

  TGAE**  g_jet_trk_dphi_ref_syst = Get1DArray <TGAE*> (nPtChSelections);
  TGAE**  g_jet_trk_dphi_ref_bkg_syst = Get1DArray <TGAE*> (nPtChSelections);
  TGAE*** g_jet_trk_dphi_syst = Get2DArray <TGAE*> (numZdcCentBins, nPtChSelections);
  TGAE*** g_jet_trk_dphi_bkg_syst = Get2DArray <TGAE*> (numZdcCentBins, nPtChSelections);

  TGAE**  g_jet_trk_dphi_ref_sig_syst = Get1DArray <TGAE*> (nPtChSelections);
  TGAE*** g_jet_trk_dphi_sig_syst = Get2DArray <TGAE*> (numZdcCentBins, nPtChSelections);

  TGAE*** g_jet_trk_dphi_iaa_syst = Get2DArray <TGAE*> (numZdcCentBins, nPtChSelections);

  TH1D***  h_jet_trk_dphi_ref_syst = Get2DArray <TH1D*> (nPtChSelections, nVar);
  TH1D***  h_jet_trk_dphi_ref_bkg_syst = Get2DArray <TH1D*> (nPtChSelections, nVar);
  TH1D**** h_jet_trk_dphi_syst = Get3DArray <TH1D*> (numZdcCentBins, nPtChSelections, nVar);
  TH1D**** h_jet_trk_dphi_bkg_syst = Get3DArray <TH1D*> (numZdcCentBins, nPtChSelections, nVar);

  TH1D***  h_jet_trk_dphi_ref_sig_syst = Get2DArray <TH1D*> (nPtChSelections, nVar);
  TH1D**** h_jet_trk_dphi_sig_syst = Get3DArray <TH1D*> (numZdcCentBins, nPtChSelections, nVar);

  TH1D**** h_jet_trk_dphi_iaa_syst = Get3DArray <TH1D*> (numZdcCentBins, nPtChSelections, nVar);



  {
    TString inFileName = inFileTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/PlotDPhi_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    h_evt_counts_ref = (TH1D*) inFile->Get ("h_evt_counts_ref_Nominal");
    h_jet_counts_ref = (TH1D*) inFile->Get ("h_jet_counts_ref_Nominal");
    h_evt_counts_ref_bkg = (TH1D*) inFile->Get ("h_evt_counts_ref_bkg_Nominal");
    h_jet_counts_ref_bkg = (TH1D*) inFile->Get ("h_jet_counts_ref_bkg_Nominal");

    for (int iCent = 0; iCent < numZdcCentBins; iCent++) { 
      h_evt_counts[iCent] = (TH1D*) inFile->Get (Form ("h_evt_counts_pPb_iCent%i_Nominal", iCent));
      h_jet_counts[iCent] = (TH1D*) inFile->Get (Form ("h_jet_counts_pPb_iCent%i_Nominal", iCent));
      h_evt_counts_bkg[iCent] = (TH1D*) inFile->Get (Form ("h_evt_counts_pPb_bkg_iCent%i_Nominal", iCent));
      h_jet_counts_bkg[iCent] = (TH1D*) inFile->Get (Form ("h_jet_counts_pPb_bkg_iCent%i_Nominal", iCent));
    }

    for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

      g_jet_trk_dphi_ref_syst[iPtCh] = (TGAE*) inFile->Get (Form ("g_jet_trk_dphi_%s_ref_syst", pTChSelections[iPtCh].Data ()));
      g_jet_trk_dphi_ref_bkg_syst[iPtCh] = (TGAE*) inFile->Get (Form ("g_jet_trk_dphi_%s_ref_bkg_syst", pTChSelections[iPtCh].Data ()));
      g_jet_trk_dphi_ref_sig_syst[iPtCh] = (TGAE*) inFile->Get (Form ("g_jet_trk_dphi_%s_ref_sig_syst", pTChSelections[iPtCh].Data ()));

      h_jet_trk_dphi_ref[iPtCh] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_pp_Nominal", pTChSelections[iPtCh].Data ()));
      h_jet_trk_dphi_ref_bkg[iPtCh] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_ref_bkg_Nominal", pTChSelections[iPtCh].Data ()));
      h_jet_trk_dphi_ref_sig[iPtCh] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_ref_sig_Nominal", pTChSelections[iPtCh].Data ()));

      for (int iCent = 0; iCent < numZdcCentBins; iCent++) { 

        g_jet_trk_dphi_syst[iCent][iPtCh] = (TGAE*) inFile->Get (Form ("g_jet_trk_dphi_%s_syst_iCent%i", pTChSelections[iPtCh].Data (), iCent));
        g_jet_trk_dphi_bkg_syst[iCent][iPtCh] = (TGAE*) inFile->Get (Form ("g_jet_trk_dphi_%s_bkg_syst_iCent%i", pTChSelections[iPtCh].Data (), iCent));
        g_jet_trk_dphi_sig_syst[iCent][iPtCh] = (TGAE*) inFile->Get (Form ("g_jet_trk_dphi_%s_sig_syst_iCent%i", pTChSelections[iPtCh].Data (), iCent));
        g_jet_trk_dphi_iaa_syst[iCent][iPtCh] = (TGAE*) inFile->Get (Form ("g_jet_trk_dphi_%s_iaa_syst_iCent%i", pTChSelections[iPtCh].Data (), iCent));

        h_jet_trk_dphi[iCent][iPtCh] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_pPb_iCent%i_Nominal", pTChSelections[iPtCh].Data (), iCent));
        h_jet_trk_dphi_bkg[iCent][iPtCh] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_pPb_bkg_iCent%i_Nominal", pTChSelections[iPtCh].Data (), iCent));
        h_jet_trk_dphi_sig[iCent][iPtCh] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_pPb_sig_iCent%i_Nominal", pTChSelections[iPtCh].Data (), iCent));
        h_jet_trk_dphi_iaa[iCent][iPtCh] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_iaa_iCent%i_Nominal", pTChSelections[iPtCh].Data (), iCent));

      }

      for (int iVar = 1; iVar < nVar; iVar++) {

        h_jet_trk_dphi_ref_syst[iPtCh][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_pp_%s", pTChSelections[iPtCh].Data (), variations[iVar].Data ()));
        h_jet_trk_dphi_ref_bkg_syst[iPtCh][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_ref_bkg_%s", pTChSelections[iPtCh].Data (), variations[iVar].Data ()));
        h_jet_trk_dphi_ref_sig_syst[iPtCh][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_ref_sig_%s", pTChSelections[iPtCh].Data (), variations[iVar].Data ()));

        for (int iCent = 0; iCent < numZdcCentBins; iCent++) {

          h_jet_trk_dphi_syst[iCent][iPtCh][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_pPb_iCent%i_%s", pTChSelections[iPtCh].Data (), iCent, variations[iVar].Data ()));
          h_jet_trk_dphi_bkg_syst[iCent][iPtCh][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_pPb_bkg_iCent%i_%s", pTChSelections[iPtCh].Data (), iCent, variations[iVar].Data ()));
          h_jet_trk_dphi_sig_syst[iCent][iPtCh][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_pPb_sig_iCent%i_%s", pTChSelections[iPtCh].Data (), iCent, variations[iVar].Data ()));
          h_jet_trk_dphi_iaa_syst[iCent][iPtCh][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_iaa_iCent%i_%s", pTChSelections[iPtCh].Data (), iCent, variations[iVar].Data ()));

        }
      }
    }
  }



  l->SetLineStyle (7);
  l->SetLineWidth (1);
  l->SetLineColor (kBlack);



  for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
    for (int iCent = 0; iCent < numZdcCentBins; iCent++) {
      const char* canvasName = Form ("c_jet_trk_dphi_%s_iCent%i", pTChSelections[iPtCh].Data (), iCent);
      TCanvas* c = new TCanvas (canvasName, "", 800, 1120);
      c->cd ();
      const double fuPad = 480./1120.;
      const double fdPad = 320./1120.;
      const double fcPad = 1.0 - fuPad - fdPad;
      TPad* uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 1.0-fuPad, 1.0, 1.0);
      TPad* cPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, fdPad, 1.0, 1.0-fuPad);
      TPad* dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, fdPad);

      uPad->SetBottomMargin (0);
      cPad->SetTopMargin (0);
      cPad->SetBottomMargin (0);
      dPad->SetTopMargin (0);
      dPad->SetBottomMargin (0.25);
      uPad->Draw ();
      cPad->Draw ();
      dPad->Draw ();

      TH1D* h = nullptr; 
      TGAE* g = nullptr;

      uPad->cd (); 

      float ymin = -4;
      float ymax = 33;

      h = (TH1D*) h_jet_trk_dphi_ref[iPtCh]->Clone ("h");
      h->Reset ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
      h->GetXaxis ()->SetTitleSize (0.028/fuPad);
      h->GetXaxis ()->SetLabelSize (0.028/fuPad);
      h->GetXaxis ()->SetTitleOffset (2.1*fuPad);
      h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#Delta#phi)");
      h->GetYaxis ()->SetTitleSize (0.028/fuPad);
      h->GetYaxis ()->SetLabelSize (0.028/fuPad);
      h->GetYaxis ()->SetTitleOffset (2.1*fuPad);

      h->SetLineWidth (1);
      h->SetLineStyle (2);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      TBox* shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
      shadedBox->SetFillColorAlpha (kGray, 0.3);
      shadedBox->Draw ();
      l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

      shadedBox = new TBox (0, ymin, M_PI/8., ymax);
      shadedBox->SetFillColorAlpha (kGray, 0.3);
      shadedBox->Draw ();
      l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

      g = g_jet_trk_dphi_ref_syst[iPtCh];
      g->SetMarkerStyle (0);
      g->SetMarkerSize (0.);
      g->SetLineColor (myBlue);
      g->SetLineWidth (1);
      ((TGAE*) g->Clone ())->Draw ("5");
      h = h_jet_trk_dphi_ref[iPtCh];
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

      g = g_jet_trk_dphi_ref_bkg_syst[iPtCh];
      g->SetMarkerStyle (0);
      g->SetMarkerSize (0.);
      g->SetLineColor (myPurple);
      g->SetLineWidth (1);
      ((TGAE*) g->Clone ())->Draw ("5");
      h = h_jet_trk_dphi_ref_bkg[iPtCh];
      g = make_graph (h);
      ResetXErrors (g);
      g->SetMarkerStyle (kOpenCircle);
      g->SetMarkerSize (0.8);
      g->SetMarkerColor (myPurple);
      g->SetLineColor (myPurple);
      g->SetLineWidth (2);
      ((TGAE*) g->Clone ())->Draw ("p");
      SaferDelete (&g);

      g = g_jet_trk_dphi_syst[iCent][iPtCh];
      g->SetMarkerStyle (0);
      g->SetMarkerSize (0.);
      g->SetLineColor (myRed);
      g->SetLineWidth (1);
      ((TGAE*) g->Clone ())->Draw ("5");
      h = h_jet_trk_dphi[iCent][iPtCh];
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

      g = g_jet_trk_dphi_bkg_syst[iCent][iPtCh];
      g->SetMarkerStyle (0);
      g->SetMarkerSize (0.);
      g->SetLineColor (myGreen);
      g->SetLineWidth (1);
      ((TGAE*) g->Clone ())->Draw ("5");
      h = h_jet_trk_dphi_bkg[iCent][iPtCh];
      g = make_graph (h);
      ResetXErrors (g);
      g->SetMarkerStyle (kOpenCircle);
      g->SetMarkerSize (0.8);
      g->SetMarkerColor (myGreen);
      g->SetLineColor (myGreen);
      g->SetLineWidth (2);
      ((TGAE*) g->Clone ())->Draw ("p");
      SaferDelete (&g);

      myText (0.30, 0.83, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
      myText (0.30, 0.77, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.020/fuPad);
      myText (0.30, 0.71, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.020/fuPad);
      myText (0.30, 0.65, kBlack, "Jet-hadron correlations", 0.020/fuPad);
      myText (0.30, 0.59, kBlack, Form ("#it{p}_{T}^{jet} > %s, %s", GetJetPtStr (tag).Data (), pTChStrs[pTChSelections[iPtCh]].Data ()), 0.020/fuPad);
      myText (0.30, 0.53, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fuPad);


      cPad->cd (); 

      ymin = -4;
      ymax = 28;

      h = (TH1D*) h_jet_trk_dphi_ref_sig[iPtCh]->Clone ("h");
      h->Reset ();
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
      h->GetXaxis ()->SetTitleSize (0.028/fdPad);
      h->GetXaxis ()->SetLabelSize (0.028/fdPad);
      h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
      //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
      h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
      h->GetYaxis ()->SetTitleSize (0.028/fdPad);
      h->GetYaxis ()->SetLabelSize (0.028/fdPad);
      h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
      h->GetYaxis ()->CenterTitle ();

      h->SetLineWidth (1);
      h->SetLineStyle (2);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
      shadedBox->SetFillColorAlpha (kGray, 0.3);
      shadedBox->Draw ();
      l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

      shadedBox = new TBox (0, ymin, M_PI/8., ymax);
      shadedBox->SetFillColorAlpha (kGray, 0.3);
      shadedBox->Draw ();
      l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

      g = g_jet_trk_dphi_ref_sig_syst[iPtCh];
      g->SetMarkerStyle (0);
      g->SetMarkerSize (0.);
      g->SetLineColor (myBlue);
      g->SetLineWidth (1);
      ((TGAE*) g->Clone ())->Draw ("5");
      h = h_jet_trk_dphi_ref_sig[iPtCh];
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

      g = g_jet_trk_dphi_sig_syst[iCent][iPtCh];
      g->SetMarkerStyle (0);
      g->SetMarkerSize (0.);
      g->SetLineColor (myRed);
      g->SetLineWidth (1);
      ((TGAE*) g->Clone ())->Draw ("5");
      h = h_jet_trk_dphi_sig[iCent][iPtCh];
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

      myBoxText2 (0.36, 0.84, myRed, kFullCircle, "#it{p}+Pb jet-tagged events", 0.8, 0.020/fcPad, true);
      myBoxText2 (0.36, 0.75, myBlue, kFullCircle, "#it{pp} jet-tagged events", 0.8, 0.020/fcPad, true);
      myBoxText2 (0.36, 0.66, myGreen, kOpenCircle, "#it{p}+Pb bkgd.", 0.8, 0.020/fcPad);
      myBoxText2 (0.36, 0.57, myPurple, kOpenCircle, "#it{pp} bkgd.", 0.8, 0.020/fcPad);


      dPad->cd (); 

      ymin = strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 0.8 : 0.83;
      ymax = strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.4 : 1.17;

      h = (TH1D*) h_jet_trk_dphi_iaa[iCent][iPtCh]->Clone ("h");
      h->Reset ();
      for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
      h->GetXaxis ()->SetTitleSize (0.028/fdPad);
      h->GetXaxis ()->SetLabelSize (0.028/fdPad);
      h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
      //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
      h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}A} = #it{p}+Pb / #it{pp}");
      h->GetYaxis ()->SetTitleSize (0.028/fdPad);
      h->GetYaxis ()->SetLabelSize (0.028/fdPad);
      h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
      h->GetYaxis ()->CenterTitle ();

      h->SetLineWidth (1);
      h->SetLineStyle (2);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
      shadedBox->SetFillColorAlpha (kGray, 0.3);
      shadedBox->Draw ();
      l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

      shadedBox = new TBox (0, ymin, M_PI/8., ymax);
      shadedBox->SetFillColorAlpha (kGray, 0.3);
      shadedBox->Draw ();
      l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

      g = g_jet_trk_dphi_iaa_syst[iCent][iPtCh];
      g->SetMarkerStyle (0);
      g->SetMarkerSize (0.);
      g->SetLineColor (myRed);
      g->SetLineWidth (1);
      ((TGAE*) g->Clone ())->Draw ("5");
      h = h_jet_trk_dphi_iaa[iCent][iPtCh];
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

      c->SaveAs (Form ("%s/Plots/DPhi/JetTagged_HadronYields_%i-%i%%_comparison_dphi_%s_%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], pTChSelections[iPtCh].Data (), tag));
    }
  }



  for (int iCent = 0; iCent < numZdcCentBins; iCent++) {
    const char* canvasName = Form ("c_jet_trk_dphi_FcalvsZdc_iCent%i", iCent);
    TCanvas* c = new TCanvas (canvasName, "", 800, 1120);
    c->cd ();
    const double fuPad = 480./1120.;
    const double fdPad = 320./1120.;
    const double fcPad = 1.0 - fuPad - fdPad;
    TPad* uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 1.0-fuPad, 1.0, 1.0);
    TPad* cPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, fdPad, 1.0, 1.0-fuPad);
    TPad* dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, fdPad);

    uPad->SetBottomMargin (0);
    cPad->SetTopMargin (0);
    cPad->SetBottomMargin (0);
    dPad->SetTopMargin (0);
    dPad->SetBottomMargin (0.25);
    uPad->Draw ();
    cPad->Draw ();
    dPad->Draw ();

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    uPad->cd (); 

    float ymin = -4;
    float ymax = 33;

    h = (TH1D*) h_jet_trk_dphi_ref[3]->Clone ("h");
    h->Reset ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
    h->GetXaxis ()->SetTitleSize (0.028/fuPad);
    h->GetXaxis ()->SetLabelSize (0.028/fuPad);
    h->GetXaxis ()->SetTitleOffset (2.1*fuPad);
    h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#Delta#phi)");
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fuPad);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    TBox* shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
    shadedBox->SetFillColorAlpha (kGray, 0.3);
    shadedBox->Draw ();
    l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

    shadedBox = new TBox (0, ymin, M_PI/8., ymax);
    shadedBox->SetFillColorAlpha (kGray, 0.3);
    shadedBox->Draw ();
    l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

    h = h_jet_trk_dphi_ref[3];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_dphi_ref_bkg[3];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_dphi[iCent][3];
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

    h = h_jet_trk_dphi_bkg[iCent][3];
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

    h = h_jet_trk_dphi_syst[iCent][3][iVar];
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

    h = h_jet_trk_dphi_bkg_syst[iCent][3][iVar];
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

    myText (0.30, 0.83, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
    myText (0.30, 0.77, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.020/fuPad);
    myText (0.30, 0.71, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.020/fuPad);
    myText (0.30, 0.65, kBlack, "Jet-hadron correlations", 0.020/fuPad);
    myText (0.30, 0.59, kBlack, Form ("#it{p}_{T}^{jet} > %s, #it{p}_{T}^{ch} > 2 GeV", GetJetPtStr (tag).Data ()), 0.020/fuPad);
    myText (0.30, 0.53, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fuPad);


    cPad->cd (); 

    ymin = -4;
    ymax = 28;

    h = (TH1D*) h_jet_trk_dphi_ref_sig[3]->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
    shadedBox->SetFillColorAlpha (kGray, 0.3);
    shadedBox->Draw ();
    l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

    shadedBox = new TBox (0, ymin, M_PI/8., ymax);
    shadedBox->SetFillColorAlpha (kGray, 0.3);
    shadedBox->Draw ();
    l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

    h = h_jet_trk_dphi_ref_sig[3];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_dphi_sig[iCent][3];
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

    h = h_jet_trk_dphi_sig_syst[iCent][3][iVar];
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

    myMarkerText (0.36, 0.84, myRed, kFullCircle, Form ("#it{p}+Pb #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.8, 0.020/fcPad, true);
    myMarkerText (0.36, 0.75, myGreen, kFullCircle, Form ("#it{p}+Pb bkgd., #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.8, 0.020/fcPad, true);
    myMarkerText (0.36, 0.66, myBlue, kFullSquare, Form ("#it{p}+Pb #bf{FCal %i-%i%%}", fcalCentPercs[iCent+1], fcalCentPercs[iCent]), 0.8, 0.020/fcPad, true);
    myMarkerText (0.36, 0.57, myOrange, kFullSquare, Form ("#it{p}+Pb bkgd., #bf{FCal %i-%i%%}", fcalCentPercs[iCent+1], fcalCentPercs[iCent]), 0.8, 0.020/fcPad, true);
    myMarkerText (0.36, 0.48, kBlack, kFullCircle, "#it{pp} jet-tagged events", 0.8, 0.020/fcPad);
    myMarkerText (0.36, 0.39, kBlack, kOpenCircle, "#it{pp} bkgd.", 0.8, 0.020/fcPad);


    dPad->cd (); 

    ymin = 0.5;
    ymax = 1.4;

    h = (TH1D*) h_jet_trk_dphi_iaa[iCent][3]->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}A} = #it{p}+Pb / #it{pp}");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
    shadedBox->SetFillColorAlpha (kGray, 0.3);
    shadedBox->Draw ();
    l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

    shadedBox = new TBox (0, ymin, M_PI/8., ymax);
    shadedBox->SetFillColorAlpha (kGray, 0.3);
    shadedBox->Draw ();
    l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

    h = h_jet_trk_dphi_iaa[iCent][3];
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

    h = h_jet_trk_dphi_iaa_syst[iCent][3][iVar];
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

    c->SaveAs (Form ("%s/Plots/DPhi/JetTagged_HadronYields_%i-%i%%_FCalvsZDC_dphi_gt2_lt4_%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], tag));
  }



  for (int iCent = 0; iCent < numZdcCentBins; iCent++) {
    const char* canvasName = Form ("c_jet_trk_dphi_gt0p5_lt1_FcalvsZdc_iCent%i", iCent);
    TCanvas* c = new TCanvas (canvasName, "", 800, 1120);
    c->cd ();
    const double fuPad = 480./1120.;
    const double fdPad = 320./1120.;
    const double fcPad = 1.0 - fuPad - fdPad;
    TPad* uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 1.0-fuPad, 1.0, 1.0);
    TPad* cPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, fdPad, 1.0, 1.0-fuPad);
    TPad* dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, fdPad);

    uPad->SetBottomMargin (0);
    cPad->SetTopMargin (0);
    cPad->SetBottomMargin (0);
    dPad->SetTopMargin (0);
    dPad->SetBottomMargin (0.25);
    uPad->Draw ();
    cPad->Draw ();
    dPad->Draw ();

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    uPad->cd (); 

    float ymin = -4;
    float ymax = 33;

    h = (TH1D*) h_jet_trk_dphi_ref[0]->Clone ("h");
    h->Reset ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
    h->GetXaxis ()->SetTitleSize (0.028/fuPad);
    h->GetXaxis ()->SetLabelSize (0.028/fuPad);
    h->GetXaxis ()->SetTitleOffset (2.1*fuPad);
    h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#Delta#phi)");
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fuPad);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    TBox* shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
    shadedBox->SetFillColorAlpha (kGray, 0.3);
    shadedBox->Draw ();
    l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

    shadedBox = new TBox (0, ymin, M_PI/8., ymax);
    shadedBox->SetFillColorAlpha (kGray, 0.3);
    shadedBox->Draw ();
    l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

    h = h_jet_trk_dphi_ref[0];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_dphi_ref_bkg[0];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_dphi[iCent][0];
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

    h = h_jet_trk_dphi_bkg[iCent][0];
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

    h = h_jet_trk_dphi_syst[iCent][0][iVar];
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

    h = h_jet_trk_dphi_bkg_syst[iCent][0][iVar];
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

    myText (0.30, 0.83, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
    myText (0.30, 0.77, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.020/fuPad);
    myText (0.30, 0.71, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.020/fuPad);
    myText (0.30, 0.65, kBlack, "Jet-hadron correlations", 0.020/fuPad);
    myText (0.30, 0.59, kBlack, Form ("#it{p}_{T}^{jet} > %s, 0.5 < #it{p}_{T}^{ch} < 1 GeV", GetJetPtStr (tag).Data ()), 0.020/fuPad);
    myText (0.30, 0.53, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fuPad);


    cPad->cd (); 

    ymin = -4;
    ymax = 28;

    h = (TH1D*) h_jet_trk_dphi_ref_sig[0]->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
    shadedBox->SetFillColorAlpha (kGray, 0.3);
    shadedBox->Draw ();
    l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

    shadedBox = new TBox (0, ymin, M_PI/8., ymax);
    shadedBox->SetFillColorAlpha (kGray, 0.3);
    shadedBox->Draw ();
    l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

    h = h_jet_trk_dphi_ref_sig[0];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_dphi_sig[iCent][0];
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

    h = h_jet_trk_dphi_sig_syst[iCent][0][iVar];
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

    myMarkerText (0.36, 0.84, myRed, kFullCircle, Form ("#it{p}+Pb #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.8, 0.020/fcPad, true);
    myMarkerText (0.36, 0.75, myGreen, kFullCircle, Form ("#it{p}+Pb bkgd., #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.8, 0.020/fcPad, true);
    myMarkerText (0.36, 0.66, myBlue, kFullSquare, Form ("#it{p}+Pb #bf{FCal %i-%i%%}", fcalCentPercs[iCent+1], fcalCentPercs[iCent]), 0.8, 0.020/fcPad, true);
    myMarkerText (0.36, 0.57, myOrange, kFullSquare, Form ("#it{p}+Pb bkgd., #bf{FCal %i-%i%%}", fcalCentPercs[iCent+1], fcalCentPercs[iCent]), 0.8, 0.020/fcPad, true);
    myMarkerText (0.36, 0.48, kBlack, kFullCircle, "#it{pp} jet-tagged events", 0.8, 0.020/fcPad);
    myMarkerText (0.36, 0.39, kBlack, kOpenCircle, "#it{pp} bkgd.", 0.8, 0.020/fcPad);


    dPad->cd (); 

    ymin = 0.1;
    ymax = 1.2;

    h = (TH1D*) h_jet_trk_dphi_iaa[iCent][0]->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}A} = #it{p}+Pb / #it{pp}");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
    shadedBox->SetFillColorAlpha (kGray, 0.3);
    shadedBox->Draw ();
    l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

    shadedBox = new TBox (0, ymin, M_PI/8., ymax);
    shadedBox->SetFillColorAlpha (kGray, 0.3);
    shadedBox->Draw ();
    l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

    h = h_jet_trk_dphi_iaa[iCent][0];
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

    h = h_jet_trk_dphi_iaa_syst[iCent][0][iVar];
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

    c->SaveAs (Form ("%s/Plots/DPhi/JetTagged_HadronYields_%i-%i%%_FCalvsZDC_dphi_gt0p5_lt1_%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], tag));
  }



  /*
  for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
    for (int iCent = 0; iCent < numZdcCentBins; iCent++) {
      const char* canvasName = Form ("c_jet_trk_dphi_%s_pPb_iCent%i_syst", pTChSelections[iPtCh].Data (), iCent);
      TCanvas* c = new TCanvas (canvasName, "", 800, 800);

      TH1D* h = nullptr; 
      TGAE* g = nullptr;

      c->cd (); 

      float ymin = -20;
      float ymax = 20;

      c->Clear ();

      h = (TH1D*) h_jet_trk_dphi[iCent][iPtCh]->Clone ("h");
      h->Reset ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
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
        h = (TH1D*) h_jet_trk_dphi_syst[iCent][iPtCh][iVar]->Clone ("htemp");
        SaveRelativeErrors (h, h_jet_trk_dphi[iCent][iPtCh], true);
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
      myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
      myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, %s", GetJetPtStr (tag).Data (), pTChStrs[pTChSelections[iPtCh]].Data ()), 0.032);
      for (int iVar = 1; iVar < nVar; iVar++)
        myLineColorText (0.25, 0.77-iVar*0.040, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 1.0, 0.028);
      

      //c->SaveAs (Form ("%s/Plots/Systematics/TotalJetTaggedYield_%s_dphi_%s_syst.pdf", workPath.Data (), iSys == 0 ? "pp" : Form ("pPb_%i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), tag));
      c->SaveAs (Form ("%s/Plots/Systematics/TotalJetTaggedYield_%s_dphi_%s_%s_syst.pdf", workPath.Data (), Form ("pPb_%i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), pTChSelections[iPtCh].Data (), tag));
    }
  }



  for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
    for (int iCent = 0; iCent < numZdcCentBins; iCent++) {
      const char* canvasName = Form ("c_jet_trk_dphi_%s_pPb_sig_iCent%i_syst", pTChSelections[iPtCh].Data (), iCent);
      TCanvas* c = new TCanvas (canvasName, "", 800, 800);

      TH1D* h = nullptr; 
      TGAE* g = nullptr;

      c->cd (); 

      float ymin = -20;
      float ymax = 20;

      c->Clear ();

      h = (TH1D*) h_jet_trk_dphi_sig[iCent][iPtCh]->Clone ("h");
      h->Reset ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
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
        h = (TH1D*) h_jet_trk_dphi_sig_syst[iCent][iPtCh][iVar]->Clone ("htemp");
        SaveRelativeErrors (h, h_jet_trk_dphi_sig[iCent][iPtCh], true);
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
      myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
      myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, %s", GetJetPtStr (tag).Data (), pTChStrs[pTChSelections[iPtCh]].Data ()), 0.032);
      for (int iVar = 1; iVar < nVar; iVar++)
        myLineColorText (0.25, 0.77-iVar*0.040, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 1.0, 0.028);
      

      //c->SaveAs (Form ("%s/Plots/Systematics/SignalJetTaggedYield_%s_dphi_%s_syst.pdf", workPath.Data (), iSys == 0 ? "pp" : Form ("pPb_%i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), tag));
      c->SaveAs (Form ("%s/Plots/Systematics/SignalJetTaggedYield_%s_dphi_%s_%s_syst.pdf", workPath.Data (), Form ("pPb_%i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), pTChSelections[iPtCh].Data (), tag));
    }
  }



  for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
    for (int iCent = 0; iCent < numZdcCentBins; iCent++) {
      const char* canvasName = Form ("c_jet_trk_dphi_%s_iaa_iCent%i_syst", pTChSelections[iPtCh].Data (), iCent);
      TCanvas* c = new TCanvas (canvasName, "", 800, 800);

      TH1D* h = nullptr; 
      TGAE* g = nullptr;

      c->cd (); 

      float ymin = -20;
      float ymax = 20;


      h = (TH1D*) h_jet_trk_dphi_iaa[iCent][iPtCh]->Clone ("h");
      h->Reset ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
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
        h = (TH1D*) h_jet_trk_dphi_iaa_syst[iCent][iPtCh][iVar]->Clone ("htemp");
        SaveRelativeErrors (h, h_jet_trk_dphi_iaa[iCent][iPtCh], true);
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
      myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, %s", GetJetPtStr (tag).Data (), pTChStrs[pTChSelections[iPtCh]].Data ()), 0.032);
      for (int iVar = 1; iVar < nVar; iVar++)
        myLineColorText (0.25, 0.73-iVar*0.040, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 1.0, 0.028);
      

      c->SaveAs (Form ("%s/Plots/Systematics/JetTagged_IpA_%i-%i%%_dphi_%s_%s_syst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], pTChSelections[iPtCh].Data (), tag));
    }
  }
  */

}


#endif
