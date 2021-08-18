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


const bool makeTotalSystPlots = false;
const bool makeBkgdSystPlots  = false;
const bool makeSigSystPlots   = false;
const bool makeIpPbSystPlots  = false;
const bool makeMCClosurePlots = true;
const bool makeUnfoldingPlots = true;

const float maxDataSyst = 20; // maximum y-axis for data-driven systematics
const float maxMCSyst = 10; // maximum y-axis for MC-driven systematics


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

  TH1D**  h_jet_trk_pt_ref_bbb  = Get1DArray <TH1D*> (nDir);
  TH1D*** h_jet_trk_pt_bbb      = Get2DArray <TH1D*> (nDir, nZdcCentBins);
  TF1**  f_jet_trk_pt_ref_bbb   = Get1DArray <TF1*> (nDir);
  TF1*** f_jet_trk_pt_bbb       = Get2DArray <TF1*> (nDir, nZdcCentBins);


  TGAE***  g_jet_trk_pt_ref_syst     = Get2DArray <TGAE*> (nDir, nVar);
  TGAE***  g_jet_trk_pt_ref_bkg_syst = Get2DArray <TGAE*> (nDir, nVar);
  TGAE**** g_jet_trk_pt_syst         = Get3DArray <TGAE*> (nDir, nZdcCentBins, nVar);
  TGAE**** g_jet_trk_pt_bkg_syst     = Get3DArray <TGAE*> (nDir, nZdcCentBins, nVar);

  TGAE***  g_jet_trk_pt_ref_sig_syst = Get2DArray <TGAE*> (nDir, nVar);
  TGAE**** g_jet_trk_pt_sig_syst     = Get3DArray <TGAE*> (nDir, nZdcCentBins, nVar);

  TGAE**** g_jet_trk_pt_iaa_syst     = Get3DArray <TGAE*> (nDir, nZdcCentBins, nVar);


  TGAE***  g_jet_trk_pt_ref_systTot     = Get2DArray <TGAE*> (nDir, 3);
  TGAE***  g_jet_trk_pt_ref_bkg_systTot = Get2DArray <TGAE*> (nDir, 3);
  TGAE**** g_jet_trk_pt_systTot         = Get3DArray <TGAE*> (nDir, nZdcCentBins, 3);
  TGAE**** g_jet_trk_pt_bkg_systTot     = Get3DArray <TGAE*> (nDir, nZdcCentBins, 3);

  TGAE***  g_jet_trk_pt_ref_sig_systTot = Get2DArray <TGAE*> (nDir, 3);
  TGAE**** g_jet_trk_pt_sig_systTot     = Get3DArray <TGAE*> (nDir, nZdcCentBins, 3);

  TGAE**** g_jet_trk_pt_iaa_systTot     = Get3DArray <TGAE*> (nDir, nZdcCentBins, 3);


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
    inFileName = Form ("%s/Results/ProcessCorrelations_%s.root", rootPath.Data (), inFileName.Data ());
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

      h_jet_trk_pt_ref_bbb[iDir] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_ref_bbb", dir.Data ()));
      f_jet_trk_pt_ref_bbb[iDir] = (TF1*)  inFile->Get (Form ("f_jet_trk_pt_%s_ref_bbb", dir.Data ()));

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        h_jet_trk_pt_bbb[iDir][iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_pPb_bbb_iCent%i", dir.Data (), iCent));
        f_jet_trk_pt_bbb[iDir][iCent] = (TF1*)  inFile->Get (Form ("f_jet_trk_pt_%s_pPb_bbb_iCent%i", dir.Data (), iCent));

      } // end loop over iCent

    } // end loop over iDir



    for (int iDir = 0; iDir < nDir; iDir++) {

      const TString dir = directions[iDir];

      for (int iVar = 0; iVar < nVar; iVar++) {

        const TString var = variations[iVar];

        g_jet_trk_pt_ref_syst[iDir][iVar]     = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_ref_syst_%s",      dir.Data (), var.Data ()));
        g_jet_trk_pt_ref_bkg_syst[iDir][iVar] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_ref_bkg_syst_%s",  dir.Data (), var.Data ()));
        g_jet_trk_pt_ref_sig_syst[iDir][iVar] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_ref_sig_syst_%s",  dir.Data (), var.Data ()));

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          g_jet_trk_pt_syst[iDir][iCent][iVar]      = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_syst_iCent%i_%s",      dir.Data (), iCent, var.Data ()));
          g_jet_trk_pt_bkg_syst[iDir][iCent][iVar]  = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_bkg_syst_iCent%i_%s",  dir.Data (), iCent, var.Data ()));
          g_jet_trk_pt_sig_syst[iDir][iCent][iVar]  = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_sig_syst_iCent%i_%s",  dir.Data (), iCent, var.Data ()));
          g_jet_trk_pt_iaa_syst[iDir][iCent][iVar]  = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_iaa_syst_iCent%i_%s",  dir.Data (), iCent, var.Data ()));

        } // end loop over iCent

      } // end loop over iVar

      for (int iTotVar = 0; iTotVar < 3; iTotVar++) {

        const TString totVar = totalVariations[iTotVar];

        g_jet_trk_pt_ref_systTot[iDir][iTotVar]     = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_ref_%s_systTot",      dir.Data (), totVar.Data ()));
        g_jet_trk_pt_ref_bkg_systTot[iDir][iTotVar] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_ref_bkg_%s_systTot",  dir.Data (), totVar.Data ()));
        g_jet_trk_pt_ref_sig_systTot[iDir][iTotVar] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_ref_sig_%s_systTot",  dir.Data (), totVar.Data ()));

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          g_jet_trk_pt_systTot[iDir][iCent][iTotVar]      = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_%s_systTot_iCent%i",      dir.Data (), totVar.Data (), iCent));
          g_jet_trk_pt_bkg_systTot[iDir][iCent][iTotVar]  = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_bkg_%s_systTot_iCent%i",  dir.Data (), totVar.Data (), iCent));
          g_jet_trk_pt_sig_systTot[iDir][iCent][iTotVar]  = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_sig_%s_systTot_iCent%i",  dir.Data (), totVar.Data (), iCent));
          g_jet_trk_pt_iaa_systTot[iDir][iCent][iTotVar]  = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_iaa_%s_systTot_iCent%i",  dir.Data (), totVar.Data (), iCent));

        } // end loop over iCent

      } // end loop over iTotVar


      for (int iVar = 1; iVar < nVar; iVar++) {

        const TString var = variations[iVar];
        const TString dType = (dataVariations.count (var) > 0 ? "data" : "mc");

        h_jet_trk_pt_ref_syst[iDir][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_ref_%s_%s", dir.Data (), dType.Data (), var.Data ()));
        h_jet_trk_pt_ref_bkg_syst[iDir][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_ref_bkg_%s_%s", dir.Data (), dType.Data (), var.Data ()));
        h_jet_trk_pt_ref_sig_syst[iDir][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_ref_sig_%s_%s", dir.Data (), dType.Data (), var.Data ()));

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          h_jet_trk_pt_syst[iDir][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_pPb_iCent%i_%s_%s", dir.Data (), iCent, dType.Data (), var.Data ()));
          h_jet_trk_pt_bkg_syst[iDir][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_pPb_bkg_iCent%i_%s_%s", dir.Data (), iCent, dType.Data (), var.Data ()));
          h_jet_trk_pt_sig_syst[iDir][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_pPb_sig_iCent%i_%s_%s", dir.Data (), iCent, dType.Data (), var.Data ()));
          h_jet_trk_pt_iaa_syst[iDir][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_iaa_iCent%i_%s_%s", dir.Data (), iCent, dType.Data (), var.Data ()));

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

      g = (TGAE*) g_jet_trk_pt_ref_syst[0][0]->Clone ();
      h = h_jet_trk_pt_ref[iDType][0];
      SetCentralValuesKeepRelativeErrors (g, h);
      myDrawSyst (g, myBlue);
      SaferDelete (&g);
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myBlue, kFullCircle, 0.8);
      SaferDelete (&g);

      g = (TGAE*) g_jet_trk_pt_ref_bkg_syst[0][0]->Clone ();
      h = h_jet_trk_pt_ref_bkg[iDType][0];
      SetCentralValuesKeepRelativeErrors (g, h);
      myDrawSyst (g, myPurple);
      SaferDelete (&g);
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myPurple, kOpenCircle, 0.8);
      SaferDelete (&g);

      g = (TGAE*) g_jet_trk_pt_syst[0][iCent][0]->Clone ();
      h = h_jet_trk_pt[iDType][0][iCent];
      SetCentralValuesKeepRelativeErrors (g, h);
      myDrawSyst (g, myRed);
      SaferDelete (&g);
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myRed, kFullCircle, 0.8);
      SaferDelete (&g);

      g = (TGAE*) g_jet_trk_pt_bkg_syst[0][iCent][0]->Clone ();
      h = h_jet_trk_pt_bkg[iDType][0][iCent];
      SetCentralValuesKeepRelativeErrors (g, h);
      myDrawSyst (g, myGreen);
      SaferDelete (&g);
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myGreen, kOpenCircle, 0.8);
      SaferDelete (&g);

      myText (0.24, 0.12, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{%s %i-%i%%}", iDType == 0 ? "ZDC" : "FCal", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.020/fuPad);
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

      g = (TGAE*) g_jet_trk_pt_ref_syst[2][0]->Clone ();
      h = h_jet_trk_pt_ref[iDType][2];
      SetCentralValuesKeepRelativeErrors (g, h);
      myDrawSyst (g, myBlue);
      SaferDelete (&g);
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myBlue, kFullCircle, 0.8);
      SaferDelete (&g);

      g = (TGAE*) g_jet_trk_pt_ref_bkg_syst[2][0]->Clone ();
      h = h_jet_trk_pt_ref_bkg[iDType][2];
      SetCentralValuesKeepRelativeErrors (g, h);
      myDrawSyst (g, myPurple);
      SaferDelete (&g);
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myPurple, kOpenCircle, 0.8);
      SaferDelete (&g);

      g = (TGAE*) g_jet_trk_pt_syst[2][iCent][0]->Clone ();
      h = h_jet_trk_pt[iDType][2][iCent];
      SetCentralValuesKeepRelativeErrors (g, h);
      myDrawSyst (g, myRed);
      SaferDelete (&g);
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myRed, kFullCircle, 0.8);
      SaferDelete (&g);

      g = (TGAE*) g_jet_trk_pt_bkg_syst[2][iCent][0]->Clone ();
      h = h_jet_trk_pt_bkg[iDType][2][iCent];
      SetCentralValuesKeepRelativeErrors (g, h);
      myDrawSyst (g, myGreen);
      SaferDelete (&g);
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

      g = (TGAE*) g_jet_trk_pt_ref_sig_syst[0][0]->Clone ();
      h = h_jet_trk_pt_ref_sig[iDType][0];
      SetCentralValuesKeepRelativeErrors (g, h);
      myDrawSyst (g, myBlue);
      SaferDelete (&g);
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myBlue, kFullCircle, 0.8);
      SaferDelete (&g);

      g = (TGAE*) g_jet_trk_pt_sig_syst[0][iCent][0]->Clone ();
      h = h_jet_trk_pt_sig[iDType][0][iCent];
      SetCentralValuesKeepRelativeErrors (g, h);
      myDrawSyst (g, myRed);
      SaferDelete (&g);
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

      g = (TGAE*) g_jet_trk_pt_ref_sig_syst[2][0]->Clone ();
      h = h_jet_trk_pt_ref_sig[iDType][2];
      SetCentralValuesKeepRelativeErrors (g, h);
      myDrawSyst (g, myBlue);
      SaferDelete (&g);
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myBlue, kFullCircle, 0.8);
      SaferDelete (&g);

      g = (TGAE*) g_jet_trk_pt_sig_syst[2][iCent][0]->Clone ();
      h = h_jet_trk_pt_sig[iDType][2][iCent];
      SetCentralValuesKeepRelativeErrors (g, h);
      myDrawSyst (g, myRed);
      SaferDelete (&g);
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myRed, kFullCircle, 0.8);
      SaferDelete (&g);


      dlPad->cd (); 
      dlPad->SetLogx ();

      ymin = 0.53;
      ymax = 1.47;//(strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.45 : 1.17);

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

      g = (TGAE*) g_jet_trk_pt_iaa_syst[0][iCent][0]->Clone ();
      h = h_jet_trk_pt_iaa[iDType][0][iCent];
      SetCentralValuesKeepRelativeErrors (g, h);
      myDrawSyst (g, myRed);
      SaferDelete (&g);
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myRed, kFullCircle, 0.8);
      SaferDelete (&g);


      drPad->cd (); 
      drPad->SetLogx ();

      ymin = 0.53;
      ymax = 1.47;//(strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.45 : 1.17);

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

      g = (TGAE*) g_jet_trk_pt_iaa_syst[2][iCent][0]->Clone ();
      h = h_jet_trk_pt_iaa[iDType][2][iCent];
      SetCentralValuesKeepRelativeErrors (g, h);
      myDrawSyst (g, myRed);
      SaferDelete (&g);
      g = make_graph (h);
      ResetXErrors (g);
      myDraw (g, myRed, kFullCircle, 0.8);
      SaferDelete (&g);

      c->SaveAs (Form ("%s/Plots/PtCh/JetTagged_HadronYields_%i-%iperc_comparison_PtCh_%s%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], tag, iDType == 1 ? "_mc" : "")); 

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

    c->SaveAs (Form ("%s/Plots/PtCh/JetTagged_HadronYields_%i-%iperc_FCalvsZDC_PtCh_%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], tag)); 
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

    c->SaveAs (Form ("%s/Plots/PtCh/JetTagged_HadronYields_%i-%iperc_FCalMixCatVar_PtCh_%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], tag)); 
  }*/




  {
    for (int iDir : {0, 2}) {
  
      const char* canvasName = Form ("c_jet_trk_pt_signal_and_ratio_%s", directions[iDir].Data ());
      TCanvas* c = new TCanvas (canvasName, "", 800, 1000);
      c->cd ();
  
      const double fPad = 600./1000.;
      TPad* uPad = new TPad (Form ("%s_uPad", canvasName), "", 0, 1-fPad, 1, 1.0);
      TPad* dPad = new TPad (Form ("%s_dPad", canvasName), "", 0, 0.0, 1, 1-fPad);
  
      uPad->SetTopMargin (0.04);
      uPad->SetBottomMargin (0);
      dPad->SetTopMargin (0);
      dPad->SetBottomMargin (0.25);
  
      uPad->SetLeftMargin (0.12);
      dPad->SetLeftMargin (0.12);
      uPad->SetRightMargin (0.03);
      dPad->SetRightMargin (0.03);
  
      uPad->Draw ();
      dPad->Draw ();

      TH1D* h = nullptr;
      TGAE* g = nullptr;

      uPad->cd (); 
      uPad->SetLogx ();
      uPad->SetLogy ();

      float ymin = 2e-9;
      float ymax = 3e4;
      h = new TH1D ("htemp", ";#it{p}_{T}^{ch} [GeV];(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]", 1, pTChBins[0], pTChBins[nPtChBins]);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetYaxis ()->SetTitleSize (0.028/fPad);
      h->GetYaxis ()->SetLabelSize (0.028/fPad);
      h->GetYaxis ()->SetTitleOffset (2.0*fPad);

      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);


      h = h_jet_trk_pt_ref_sig[0][iDir];
      g = (TGAE*) g_jet_trk_pt_ref_sig_syst[iDir][0]->Clone ();
      SetCentralValuesKeepRelativeErrors (g, h);
      ScaleGraph (g, nullptr, std::pow (10, 3));
      myDrawSystFill (g, lauraSystColors[0], 1, 1001);
      SaferDelete (&g);

      g = make_graph (h);
      ScaleGraph (g, nullptr, std::pow (10, 3));
      ResetXErrors (g);
      myDraw (g, lauraColors[0], kFullCircle, 1.4, 1, 2, "P", false);
      SaferDelete (&g);


      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        h = (TH1D*) h_jet_trk_pt_ref_sig[0][iDir]->Clone ("htemp");
        h->Scale (std::pow (10, 2-iCent));
        myDrawHist (h, kBlack);
        SaferDelete (&h);

        g = (TGAE*) g_jet_trk_pt_sig_syst[iDir][iCent][0]->Clone ();
        h = h_jet_trk_pt_sig[0][iDir][iCent];
        SetCentralValuesKeepRelativeErrors (g, h);
        ScaleGraph (g, nullptr, std::pow (10, 2-iCent));
        myDrawSystFill (g, lauraSystColors[iCent+1], 1.0, 1001);
        SaferDelete (&g);
  
        g = make_graph (h);
        ResetXErrors (g);
        ScaleGraph (g, nullptr, std::pow (10, 2-iCent));
        myDraw (g, lauraColors[iCent+1], kFullCircle, 1.4, 1, 2, "P", false);
        SaferDelete (&g);

      } // end loop over iCent

      myText (0.56, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.038);
      myText (0.56, 0.84, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} %s", GetJetPtStr (tag).Data (), directions[iDir] == "ns" ? "< #pi/8" : "> 7#pi/8"), 0.038);

      myText (0.18, 0.29, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.038);
      myText (0.18, 0.24, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.038);

      mySimpleMarkerAndBoxAndLineText (0.27, 0.16, 1.4, 1001, lauraSystColors[0], 1.0, lauraColors[0], kFullCircle, 1.6, "#it{pp} (#times10^{3})", 0.022/fPad);
      for (int iCent = 0; iCent < nZdcCentBins; iCent++) { 
        mySimpleMarkerAndBoxAndLineText (0.27 + (iCent >= 2 ? 0.4 : 0), 0.16-((iCent+1)%3)*0.04, 1.4, 1001, lauraSystColors[iCent+1], 1.0, lauraColors[iCent+1], kFullCircle, 1.6, Form ("ZDC %i-%i%% (#times10^{%i})", zdcCentPercs[iCent+1], zdcCentPercs[iCent], 2-iCent), 0.022/fPad);
      }
      mySimpleMarkerAndBoxAndLineText (0.27, 0.04, 1.4, 1001, kWhite, 0.0, kBlack, kDot, 0.0, "#it{pp} (#it{scaled to} #it{p}+Pb)", 0.022/fPad);


      dPad->cd ();
      dPad->SetLogx();
  
      ymin = 0.15;
      ymax = 1.85;

      h = new TH1D ("htemp", ";#it{p}_{T}^{ch} [GeV];#it{I}_{pPb} + Const.", 1, pTChBins[0], pTChBins[nPtChBins]);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetXaxis ()->SetTitleSize (0.028/(1-fPad));
      h->GetXaxis ()->SetLabelSize (0.028/(1-fPad));
      h->GetYaxis ()->SetTitleSize (0.028/(1-fPad));
      h->GetYaxis ()->SetLabelSize (0.028/(1-fPad));
      h->GetYaxis ()->SetTitleOffset (2.0*(1-fPad));
      h->GetYaxis ()->CenterTitle ();

      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      double x, y;
      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        const int maxx = 30;//(strcmp (tag, "30GeVJets") == 0 ? (iCent < 3 ? 15 : 30) : 40);
        const double offset = 0.25*(2-iCent);

        g = (TGAE*) g_jet_trk_pt_iaa_syst[iDir][iCent][0]->Clone ();
        h = h_jet_trk_pt_iaa[0][iDir][iCent];
        SetCentralValuesKeepRelativeErrors (g, h);
        for (int i = 0; i< g->GetN (); i++) {
          g->GetPoint (i, x, y);
          g->SetPoint (i, x, y+offset);
        }
        for (int i = 0; i< g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x > maxx) {
            g->RemovePoint (i);
            i--;
          }
        }
        myDrawSystFill (g, lauraSystColors[iCent+1], 1.0, 1001);
        SaferDelete (&g);
      
        g = make_graph (h);
        ResetXErrors (g);
        for (int i = 0; i< g->GetN (); i++) {
          g->GetPoint (i, x, y);
          g->SetPoint (i, x, y+offset);
        }
        for (int i = 0; i< g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x > maxx) {
            g->RemovePoint (i);
            i--;
          }
        }
        myDraw (g, lauraColors[iCent+1], kFullCircle, 1.4, 1, 2, "P", false);
        SaferDelete (&g);

        l->SetLineWidth (3);
        l->SetLineColor (lauraColors[iCent+1]);
        l->DrawLine (pTChBins[0], 1+offset, pTChBins[nPtChBins], 1+offset);

        tl->SetTextAlign (offset > 0 ? 31 : 33);
        tl->SetTextFont (43);
        tl->SetTextSize (20);
        tl->SetTextColor (lauraColors[iCent+1]);
        tl->DrawLatex (pTChBins[0] * std::exp (0.97 * std::log (pTChBins[nPtChBins]/pTChBins[0])), 1+offset+(offset > 0 ? 0.015:-0.015), Form ("#bf{#it{%s%g}}", offset == 0 ? "" : (offset > 0 ? "+ " : "#minus "), std::fabs (offset)));

      } // end loop over iCent

      c->SaveAs (Form ("%s/Plots/PtCh/SignalJetTaggedYield_%s_%s.pdf", workPath.Data (), tag, directions[iDir] == "ns" ? "nearside" : "awayside"));
  
    } // end loop over iDir

  }


  if (makeUnfoldingPlots) {

    const char* canvasName = "c_jet_trk_pt_bbb";
    TCanvas* c = new TCanvas (canvasName, "", 800, 1000);
    c->Divide (2, 3);

    {
      c->cd (6);

      gPad->SetLogx ();

      TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch, truth} [GeV];Truth-level / Reco.-level", 1, pTChBins[0], pTChBins[nPtChBins]);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (0.5, 1.5);
      h->GetYaxis ()->CenterTitle ();
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      for (int iDir : {0, 2}) {
        myDraw (h_jet_trk_pt_ref_bbb[iDir], directions[iDir] == "ns" ? myRed : myBlue, kOpenCircle, 1.0);
        myDraw (f_jet_trk_pt_ref_bbb[iDir], directions[iDir] == "ns" ? myRed : myBlue, 2, 2);
      }

      myText (0.2, 0.84, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.06);
      myText (0.2, 0.76, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.06);
    }

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      c->cd (nZdcCentBins-iCent);

      gPad->SetLogx ();

      TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch, truth} [GeV];Truth-level / Reco.-level", 1, pTChBins[0], pTChBins[nPtChBins]);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (0.5, 1.5);
      h->GetYaxis ()->CenterTitle ();
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      for (int iDir : {0, 2}) {
        myDraw (h_jet_trk_pt_bbb[iDir][iCent], directions[iDir] == "ns" ? myRed : myBlue, kOpenCircle, 1.0);
        myDraw (f_jet_trk_pt_bbb[iDir][iCent], directions[iDir] == "ns" ? myRed : myBlue, 2, 2);
      }

      myText (0.2, 0.84, kBlack, Form ("FCal %i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.06);
      if (iCent == nZdcCentBins-1) {
        myText (0.2, 0.75, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.06);
        myText (0.2, 0.66, kBlack, Form ("#it{p}_{T}^{jet} > %s", GetJetPtStr (tag).Data ()), 0.06);
      }
      else if (iCent == nZdcCentBins-2) {
        myLineText2 (0.6, 0.84, myRed, kOpenCircle, "Near side", 1.2, 0.06);
        myLineText2 (0.6, 0.75, myBlue, kOpenCircle, "Away side", 1.2, 0.06);
      }

    } // end loop over iCent

    c->SaveAs (Form ("%s/Plots/PtCh/BBBUnfoldingFactors_%s.pdf", workPath.Data (), tag));

  }




  if (makeMCClosurePlots) {
    int iMCTruthLevel = 0;
    while (iMCTruthLevel < nVar && strcmp (variations[iMCTruthLevel], "MCTruthLevel") != 0) iMCTruthLevel++;
    //while (iMCTruthLevel < nVar && strcmp (variations[iMCTruthLevel], "MCTruthJESSmear") != 0) iMCTruthLevel++;
    if (iMCTruthLevel == nVar) {
      std::cout << "Cannot find MC truth-level result? Please check!" << std::endl;
      return;
    }

    for (int iDir : {0, 2}) {
  
      const char* canvasName = Form ("c_jet_trk_pt_mcClosure_%s", directions[iDir].Data ());
      TCanvas* c = new TCanvas (canvasName, "", 800, 1000);
      c->cd ();
  
      const double fPad = 600./1000.;
      TPad* uPad = new TPad (Form ("%s_uPad", canvasName), "", 0, 1-fPad, 1, 1.0);
      TPad* dPad = new TPad (Form ("%s_dPad", canvasName), "", 0, 0.0, 1, 1-fPad);
  
      uPad->SetTopMargin (0.04);
      uPad->SetBottomMargin (0);
      dPad->SetTopMargin (0);
      dPad->SetBottomMargin (0.25);
  
      uPad->SetLeftMargin (0.12);
      dPad->SetLeftMargin (0.12);
      uPad->SetRightMargin (0.03);
      dPad->SetRightMargin (0.03);
  
      uPad->Draw ();
      dPad->Draw ();

      TH1D* h = nullptr;
      TGAE* g = nullptr;

      uPad->cd (); 
      uPad->SetLogx ();
      uPad->SetLogy ();

      float ymin = 8e-9;
      float ymax = 3e4;
      h = new TH1D ("htemp", ";#it{p}_{T}^{ch} [GeV];(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]", 1, pTChBins[0], pTChBins[nPtChBins]);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetYaxis ()->SetTitleSize (0.028/fPad);
      h->GetYaxis ()->SetLabelSize (0.028/fPad);
      h->GetYaxis ()->SetTitleOffset (2.0*fPad);

      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);


      h = h_jet_trk_pt_ref[1][iDir];
      g = (TGAE*) g_jet_trk_pt_ref_syst[iDir][0]->Clone ();
      SetCentralValuesKeepRelativeErrors (g, h);
      ScaleGraph (g, nullptr, std::pow (10, 3));
      //myDrawSyst (g, kBlack);
      myDrawSystFill (g, lauraSystColors[0], 1, 1001);
      SaferDelete (&g);

      g = make_graph (h);
      ScaleGraph (g, nullptr, std::pow (10, 3));
      ResetXErrors (g);
      //myDraw (g, kBlack, kFullCircle, 0.8);
      myDraw (g, lauraColors[0], kFullCircle, 1.4, 1, 2, "P", false);
      SaferDelete (&g);

      h = (TH1D*) h_jet_trk_pt_ref_syst[iDir][iMCTruthLevel]->Clone ("htemp");
      h->Scale (std::pow (10, 3));
      //myDrawHist (h, kBlack);
      myDrawHist (h, lauraColors[0]);
      SaferDelete (&h);

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        h = h_jet_trk_pt_sig[1][iDir][iCent];
        g = (TGAE*) g_jet_trk_pt_sig_syst[iDir][iCent][0]->Clone ();
        SetCentralValuesKeepRelativeErrors (g, h);
        ScaleGraph (g, nullptr, std::pow (10, 2-iCent));
        //myDrawSyst (g, colors[iCent]);
        myDrawSystFill (g, lauraSystColors[iCent+1], 1.0, 1001);
        SaferDelete (&g);

        g = make_graph (h);
        ScaleGraph (g, nullptr, std::pow (10, 2-iCent));
        ResetXErrors (g);
        //myDraw (g, colors[iCent], kFullCircle, 0.8);
        myDraw (g, lauraColors[iCent+1], kFullCircle, 1.4, 1, 2, "P", false);
        SaferDelete (&g);

        h = (TH1D*) h_jet_trk_pt_syst[iDir][iCent][iMCTruthLevel]->Clone ("htemp");
        h->Scale (std::pow (10, 2-iCent));
        //myDrawHist (h, colors[iCent]);
        myDrawHist (h, lauraColors[iCent+1]);
        SaferDelete (&h);

      } // end loop over iCent

      myText (0.52, 0.89, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.024/fPad);
      myText (0.62, 0.12, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.024/fPad);
      myText (0.62, 0.07, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.024/fPad);
      myText (0.52, 0.85, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} %s", GetJetPtStr (tag).Data (), directions[iDir] == "ns" ? "< #pi/8" : "> 7#pi/8"), 0.024/fPad);

      myText (0.170, 0.36, kBlack, "Truth", 0.022/fPad);
      myText (0.245, 0.36, kBlack, "Reco.", 0.022/fPad);

      mySimpleMarkerAndBoxAndLineText (0.24, 0.32, 1.4, 1001, kWhite, 0.0, lauraColors[0], kDot, 0.0, "", 0.022/fPad);
      mySimpleMarkerAndBoxAndLineText (0.32, 0.32, 1.4, 1001, lauraSystColors[0], 1.0, lauraColors[0], kFullCircle, 1.6, "#it{pp} (#times10^{3})", 0.022/fPad);
      for (int iCent = 0; iCent < nZdcCentBins; iCent++) { 
        mySimpleMarkerAndBoxAndLineText (0.24, 0.27-iCent*0.05, 1.4, 1001, kWhite, 0.0, lauraColors[iCent+1], kDot, 0.0, "", 0.022/fPad);
        mySimpleMarkerAndBoxAndLineText (0.32, 0.27-iCent*0.05, 1.4, 1001, lauraSystColors[iCent+1], 1.0, lauraColors[iCent+1], kFullCircle, 1.6, Form ("FCal %i-%i%% (#times10^{%i})", zdcCentPercs[iCent+1], zdcCentPercs[iCent], 2-iCent), 0.022/fPad);
      }


      dPad->cd ();
      dPad->SetLogx();
  
      ymin = 0.0;
      ymax = 2.0;

      h = new TH1D ("htemp", ";#it{p}_{T}^{ch} [GeV];Reco. / Truth", 1, pTChBins[0], pTChBins[nPtChBins]);
      h->SetBinContent (1, 1);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetXaxis ()->SetTitleSize (0.028/(1-fPad));
      h->GetXaxis ()->SetLabelSize (0.028/(1-fPad));
      h->GetYaxis ()->SetTitleSize (0.028/(1-fPad));
      h->GetYaxis ()->SetLabelSize (0.028/(1-fPad));
      h->GetYaxis ()->SetTitleOffset (2.0*(1-fPad));
      h->GetYaxis ()->CenterTitle ();

      h->SetLineWidth (1);
      h->SetLineStyle (2);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      double x, y;
      {
        const int maxx = 30;//(strcmp (tag, "30GeVJets") == 0 ? (iCent < 3 ? 15 : 30) : 40);
        const double offset = 0.25*3;

        h = h_jet_trk_pt_ref[1][iDir];
        g = (TGAE*) g_jet_trk_pt_ref_syst[iDir][0]->Clone ();
        SetCentralValuesKeepRelativeErrors (g, h);
        ScaleGraph (g, h_jet_trk_pt_ref_syst[iDir][iMCTruthLevel]);
        for (int i = 0; i< g->GetN (); i++) {
          g->GetPoint (i, x, y);
          g->SetPoint (i, x, y+offset);
        }
        for (int i = 0; i< g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x > maxx) {
            g->RemovePoint (i);
            i--;
          }
        }
        myDrawSystFill (g, lauraSystColors[0], 1.0, 1001);
        SaferDelete (&g);
      
        g = make_graph (h);
        ScaleGraph (g, h_jet_trk_pt_ref_syst[iDir][iMCTruthLevel]);
        ResetXErrors (g);
        for (int i = 0; i< g->GetN (); i++) {
          g->GetPoint (i, x, y);
          g->SetPoint (i, x, y+offset);
        }
        for (int i = 0; i< g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x > maxx) {
            g->RemovePoint (i);
            i--;
          }
        }
        myDraw (g, lauraColors[0], kFullCircle, 1.4, 1, 2, "P", false);
        SaferDelete (&g);

        l->SetLineWidth (3);
        l->SetLineColor (lauraColors[0]);
        l->DrawLine (pTChBins[0], 1+offset, pTChBins[nPtChBins], 1+offset);

        tl->SetTextAlign (offset > 0 ? 31 : 33);
        tl->SetTextFont (43);
        tl->SetTextSize (20);
        tl->SetTextColor (lauraColors[0]);
        tl->DrawLatex (pTChBins[0] * std::exp (0.97 * std::log (pTChBins[nPtChBins]/pTChBins[0])), 1+offset+(offset > 0 ? 0.015:-0.015), Form ("#bf{#it{%s%g}}", offset == 0 ? "" : (offset > 0 ? "+ " : "#minus "), std::fabs (offset)));

      }


      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        const int maxx = 30;//(strcmp (tag, "30GeVJets") == 0 ? (iCent < 3 ? 15 : 30) : 40);
        const double offset = 0.25*(2-iCent);

        h = h_jet_trk_pt_sig[1][iDir][iCent];
        g = (TGAE*) g_jet_trk_pt_sig_syst[iDir][iCent][0]->Clone ();
        SetCentralValuesKeepRelativeErrors (g, h);
        ScaleGraph (g, h_jet_trk_pt_syst[iDir][iCent][iMCTruthLevel]);
        for (int i = 0; i< g->GetN (); i++) {
          g->GetPoint (i, x, y);
          g->SetPoint (i, x, y+offset);
        }
        for (int i = 0; i< g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x > maxx) {
            g->RemovePoint (i);
            i--;
          }
        }
        myDrawSystFill (g, lauraSystColors[iCent+1], 1.0, 1001);
        SaferDelete (&g);
      
        g = make_graph (h);
        ScaleGraph (g, h_jet_trk_pt_syst[iDir][iCent][iMCTruthLevel]);
        ResetXErrors (g);
        for (int i = 0; i< g->GetN (); i++) {
          g->GetPoint (i, x, y);
          g->SetPoint (i, x, y+offset);
        }
        for (int i = 0; i< g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x > maxx) {
            g->RemovePoint (i);
            i--;
          }
        }
        myDraw (g, lauraColors[iCent+1], kFullCircle, 1.4, 1, 2, "P", false);
        SaferDelete (&g);

        l->SetLineWidth (3);
        l->SetLineColor (lauraColors[iCent+1]);
        l->DrawLine (pTChBins[0], 1+offset, pTChBins[nPtChBins], 1+offset);

        tl->SetTextAlign (offset > 0 ? 31 : 33);
        tl->SetTextFont (43);
        tl->SetTextSize (20);
        tl->SetTextColor (lauraColors[iCent+1]);
        tl->DrawLatex (pTChBins[0] * std::exp (0.97 * std::log (pTChBins[nPtChBins]/pTChBins[0])), 1+offset+(offset > 0 ? 0.015:-0.015), Form ("#bf{#it{%s%g}}", offset == 0 ? "" : (offset > 0 ? "+ " : "#minus "), std::fabs (offset)));

      } // end loop over iCent

      c->SaveAs (Form ("%s/Plots/PtCh/MCClosure_%s_%s.pdf", workPath.Data (), tag, directions[iDir] == "ns" ? "nearside" : "awayside"));
  
    } // end loop over iDir



    const char* canvasName = "c_jet_trk_pt_closureRatios";
    TCanvas* c = new TCanvas (canvasName, "", 800, 1000);
    c->Divide (2, 3);

    TH1D* h = nullptr;
    TGAE* g = nullptr;

    {
      c->cd (6);

      gPad->SetLogx ();

      h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];Closure (Reco. / Truth)", 1, pTChBins[0], pTChBins[nPtChBins]);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (0.7, 1.3);
      h->GetYaxis ()->CenterTitle ();
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      for (int iDir : {0, 2}) {
        const Color_t col = (directions[iDir] == "ns" ? myLitePurple : (directions[iDir] == "perp" ? kBlack : myLiteBlue));

        h = h_jet_trk_pt_ref[1][iDir];
        g = (TGAE*) g_jet_trk_pt_ref_syst[iDir][0]->Clone ();
        SetCentralValuesKeepRelativeErrors (g, h);
        ScaleGraph (g, h_jet_trk_pt_ref_sig_syst[iDir][iMCTruthLevel]);
        myDrawSyst (g, col, 1, 1, 0.3, "3");
        SaferDelete (&g);
      
        g = make_graph (h);
        ScaleGraph (g, h_jet_trk_pt_ref_sig_syst[iDir][iMCTruthLevel]);
        ResetXErrors (g);
        myDraw (g, col, kOpenCircle, 0.8, 1, 2, "P", false);
        SaferDelete (&g);
      }

      myText (0.2, 0.84, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.06);
      myText (0.2, 0.76, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.06);
    }

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      c->cd (nZdcCentBins-iCent);

      gPad->SetLogx ();

      h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];Closure (Reco. / Truth)", 1, pTChBins[0], pTChBins[nPtChBins]);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (0.7, 1.3);
      h->GetYaxis ()->CenterTitle ();
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      for (int iDir : {0, 2}) {
        const Color_t col = (directions[iDir] == "ns" ? myLitePurple : (directions[iDir] == "perp" ? kBlack : myLiteBlue));

        h = h_jet_trk_pt_sig[1][iDir][iCent];
        g = (TGAE*) g_jet_trk_pt_sig_syst[iDir][iCent][0]->Clone ();
        SetCentralValuesKeepRelativeErrors (g, h);
        ScaleGraph (g, h_jet_trk_pt_sig_syst[iDir][iCent][iMCTruthLevel]);
        myDrawSyst (g, col, 1, 1, 0.3, "3");
        SaferDelete (&g);
      
        g = make_graph (h);
        ScaleGraph (g, h_jet_trk_pt_sig_syst[iDir][iCent][iMCTruthLevel]);
        ResetXErrors (g);
        myDraw (g, col, kOpenCircle, 0.8, 1, 2, "P", false);
        SaferDelete (&g);
      }

      myText (0.2, 0.84, kBlack, Form ("FCal %i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.06);
      if (iCent == nZdcCentBins-1) {
        myText (0.2, 0.75, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.06);
        myText (0.2, 0.66, kBlack, Form ("#it{p}_{T}^{jet} > %s", GetJetPtStr (tag).Data ()), 0.06);
      }
      else if (iCent == nZdcCentBins-2) {
        myLineText2 (0.6, 0.84, myLitePurple, kOpenCircle, "Near side", 1.2, 0.06);
        myLineText2 (0.6, 0.75, myLiteBlue,   kOpenCircle, "Away side", 1.2, 0.06);
      }

    } // end loop over iCent

    c->SaveAs (Form ("%s/Plots/PtCh/MCClosureRatios_%s.pdf", workPath.Data (), tag));

  }



 
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
    myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} %s", GetJetPtStr (tag).Data (), directions[iDir] == "ns" ? "< #pi/8" : "> 7#pi/8"), 0.032);
    myLineText2 (0.25, 0.72, kBlack, kFullCircle, "#bf{#it{pp}}", 0.8, 0.032, true);
    for (int iCent = 0; iCent < nZdcCentBins; iCent++)
      myLineText2 (0.25, 0.68-iCent*0.04, colors[iCent], kFullCircle, Form ("#bf{#it{p}+Pb, %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.8, 0.032, true);

    c->SaveAs (Form ("%s/Plots/PtCh/SigToBkgd_%s_%s_ptch.pdf", workPath.Data (), tag, directions[iDir] == "ns" ? "nearside" : "awayside"));
  }




  for (int iDType : {0, 1}) {

    const int maxNSys = (iDType == 0 ? 6 : 10);

    for (int iDir : {0, 2}) {

      if (makeTotalSystPlots) {
        {
          const char* canvasName = Form ("c_jet_trk_pt_%s_pp_%s_syst", directions[iDir].Data (), iDType == 0 ? "data" : "mc");
          TCanvas* c = new TCanvas (canvasName, "", 800, 800);

          TH1D* h = nullptr;
          TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

          c->cd (); 
          c->SetLogx ();

          const float ymin = (iDType == 0 ? -maxDataSyst : -maxMCSyst);
          const float ymax = (iDType == 0 ?  maxDataSyst :  maxMCSyst);

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

          std::vector <int> totVars (0);
          if (iDType == 0) {
            totVars.push_back (0);
            totVars.push_back (1);
          }
          else 
            totVars.push_back (2);
          for (int iTotVar : totVars) {
            const TString totVar = totalVariations[iTotVar];
            g = ConvertErrorsToCentralValues (g_jet_trk_pt_ref_systTot[iDir][iTotVar], true, 100);
            gup = (TGAE*) g->Clone ();
            gdown = (TGAE*) g->Clone ();
            FlipTGAE (gdown);
         
            myDraw (gup, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDraw (gdown, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDrawFill (gup, gdown, varStyles[totVar].first, 0.1);

            SaferDelete (&g);
            SaferDelete (&gup);
            SaferDelete (&gdown);
          }

          for (int iVar = 1; iVar < nVar; iVar++) {
            const TString var = variations[iVar];
            if ((iDType == 0  && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) == 0))
              continue;

            g = (TGAE*) g_jet_trk_pt_ref_syst[iDir][iVar]->Clone ("gtemp");
            SaveRelativeErrors (g, g, g_jet_trk_pt_ref_syst[iDir][0], 100);
            ResetXErrors (g);
            ResetTGAEErrors (g);
            myDraw (g, varStyles[var].first, kDot, 0, varStyles[var].second, 4, "L");
            SaferDelete (&g);
          }

          myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
          myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
          myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} %s", GetJetPtStr (tag).Data (), directions[iDir] == "ns" ? "< #pi/8" : "> 7#pi/8"), 0.032);

          int count = 0;
          for (int iVar = 1; iVar < nVar; iVar++) {
            const TString var = variations[iVar];
            if ((iDType == 0  && dataVariations.count (var) > 0) || (iDType == 1 && mcVariations.count (var) > 0)) {
              myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[var].first, varStyles[var].second, varFullNames[var], 0.5, 0.014);
              count++;
            }
          }
          for (int iTotVar : totVars) {
            const TString totVar = totalVariations[iTotVar];
            myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 0.5, 0.014);
            count++;
          }

          c->SaveAs (Form ("%s/Plots/Systematics/PtCh/TotalJetTaggedYield_pp_%s_ptch_%s_%ssyst.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : "awayside", tag, iDType == 0 ? "":"mcBased_"));
        }

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          const char* canvasName = Form ("c_jet_trk_pt_%s_%s_pPb_iCent%i_syst", directions[iDir].Data (), iDType == 0 ? "data" : "mc", iCent);
          TCanvas* c = new TCanvas (canvasName, "", 800, 800);

          TH1D* h = nullptr;
          TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

          c->cd (); 
          c->SetLogx ();

          const float ymin = (iDType == 0 ? -maxDataSyst : -maxMCSyst);
          const float ymax = (iDType == 0 ?  maxDataSyst :  maxMCSyst);

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

          std::vector <int> totVars (0);
          if (iDType == 0) {
            totVars.push_back (0);
            totVars.push_back (1);
          }
          else 
            totVars.push_back (2);
          for (int iTotVar : totVars) {
            const TString totVar = totalVariations[iTotVar];
            g = ConvertErrorsToCentralValues (g_jet_trk_pt_systTot[iDir][iCent][iTotVar], true, 100);
            gup = (TGAE*) g->Clone ();
            gdown = (TGAE*) g->Clone ();
            FlipTGAE (gdown);
         
            myDraw (gup, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDraw (gdown, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDrawFill (gup, gdown, varStyles[totVar].first, 0.1);

            SaferDelete (&g);
            SaferDelete (&gup);
            SaferDelete (&gdown);
          }

          for (int iVar = 1; iVar < nVar; iVar++) {
            const TString var = variations[iVar];
            if ((iDType == 0  && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) == 0))
              continue;

            g = (TGAE*) g_jet_trk_pt_syst[iDir][iCent][iVar]->Clone ("gtemp");
            SaveRelativeErrors (g, g, g_jet_trk_pt_syst[iDir][iCent][0], 100);
            ResetXErrors (g);
            ResetTGAEErrors (g);
            myDraw (g, varStyles[var].first, kDot, 0, varStyles[var].second, 4, "L");
            SaferDelete (&g);
          }

          myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
          myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{%s %i-%i%%}", iDType == 0 ? "ZDC" : "FCal", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
          myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} %s", GetJetPtStr (tag).Data (), directions[iDir] == "ns" ? "< #pi/8" : "> 7#pi/8"), 0.032);

          int count = 0;
          for (int iVar = 1; iVar < nVar; iVar++) {
            const TString var = variations[iVar];
            if ((iDType == 0  && dataVariations.count (var) > 0) || (iDType == 1 && mcVariations.count (var) > 0)) {
              myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[var].first, varStyles[var].second, varFullNames[var], 0.5, 0.014);
              count++;
            }
          }
          for (int iTotVar : totVars) {
            const TString totVar = totalVariations[iTotVar];
            myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 0.5, 0.014);
            count++;
          }
          
          c->SaveAs (Form ("%s/Plots/Systematics/PtCh/TotalJetTaggedYield_pPb_%i-%iperc_%s_ptch_%s_%ssyst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", tag, iDType == 0 ? "":"mcBased_"));
        } // end loop over iCent
      }



      if (makeBkgdSystPlots && iDType != 1) {
        {
          const char* canvasName = Form ("c_jet_trk_pt_%s_pp_%s_bkg_syst", directions[iDir].Data (), iDType == 0 ? "data" : "mc");
          TCanvas* c = new TCanvas (canvasName, "", 800, 800);

          TH1D* h = nullptr;
          TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

          c->cd (); 
          c->SetLogx ();

          const float ymin = (iDType == 0 ? -maxDataSyst : -maxMCSyst);
          const float ymax = (iDType == 0 ?  maxDataSyst :  maxMCSyst);

          h = (TH1D*) h_jet_trk_pt_ref_bkg[iDType][iDir]->Clone ("h");
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

          std::vector <int> totVars (0);
          if (iDType == 0) {
            totVars.push_back (0);
            totVars.push_back (1);
          }
          else 
            totVars.push_back (2);
          for (int iTotVar : totVars) {
            const TString totVar = totalVariations[iTotVar];
            g = ConvertErrorsToCentralValues (g_jet_trk_pt_ref_bkg_systTot[iDir][iTotVar], true, 100);
            gup = (TGAE*) g->Clone ();
            gdown = (TGAE*) g->Clone ();
            FlipTGAE (gdown);
         
            myDraw (gup, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDraw (gdown, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDrawFill (gup, gdown, varStyles[totVar].first, 0.1);

            SaferDelete (&g);
            SaferDelete (&gup);
            SaferDelete (&gdown);
          }

          for (int iVar = 1; iVar < nVar; iVar++) {
            const TString var = variations[iVar];
            if ((iDType == 0  && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) == 0))
              continue;

            g = (TGAE*) g_jet_trk_pt_ref_bkg_syst[iDir][iVar]->Clone ("gtemp");
            SaveRelativeErrors (g, g, g_jet_trk_pt_ref_bkg_syst[iDir][0], 100);
            ResetXErrors (g);
            ResetTGAEErrors (g);
            myDraw (g, varStyles[var].first, kDot, 0, varStyles[var].second, 4, "L");
            SaferDelete (&g);
          }

          myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
          myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
          myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} %s", GetJetPtStr (tag).Data (), directions[iDir] == "ns" ? "< #pi/8" : "> 7#pi/8"), 0.032);

          int count = 0;
          for (int iVar = 1; iVar < nVar; iVar++) {
            const TString var = variations[iVar];
            if ((iDType == 0  && dataVariations.count (var) > 0) || (iDType == 1 && mcVariations.count (var) > 0)) {
              myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[var].first, varStyles[var].second, varFullNames[var], 0.5, 0.014);
              count++;
            }
          }
          for (int iTotVar : totVars) {
            const TString totVar = totalVariations[iTotVar];
            myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 0.5, 0.014);
            count++;
          }

          c->SaveAs (Form ("%s/Plots/Systematics/PtCh/BkgdJetTaggedYield_pp_%s_ptch_%s_%ssyst.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : "awayside", tag, iDType == 0 ? "":"mcBased_"));
        }

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          const char* canvasName = Form ("c_jet_trk_pt_%s_%s_pPb_iCent%i_bkgd_syst", directions[iDir].Data (), iDType == 0 ? "data" : "mc", iCent);
          TCanvas* c = new TCanvas (canvasName, "", 800, 800);

          TH1D* h = nullptr;
          TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

          c->cd (); 
          c->SetLogx ();

          const float ymin = (iDType == 0 ? -maxDataSyst : -maxMCSyst);
          const float ymax = (iDType == 0 ?  maxDataSyst :  maxMCSyst);

          h = (TH1D*) h_jet_trk_pt_bkg[iDType][iDir][iCent]->Clone ("h");
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

          std::vector <int> totVars (0);
          if (iDType == 0) {
            totVars.push_back (0);
            totVars.push_back (1);
          }
          else 
            totVars.push_back (2);
          for (int iTotVar : totVars) {
            const TString totVar = totalVariations[iTotVar];
            g = ConvertErrorsToCentralValues (g_jet_trk_pt_bkg_systTot[iDir][iCent][iTotVar], true, 100);
            gup = (TGAE*) g->Clone ();
            gdown = (TGAE*) g->Clone ();
            FlipTGAE (gdown);
         
            myDraw (gup, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDraw (gdown, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDrawFill (gup, gdown, varStyles[totVar].first, 0.1);

            SaferDelete (&g);
            SaferDelete (&gup);
            SaferDelete (&gdown);
          }

          for (int iVar = 1; iVar < nVar; iVar++) {
            const TString var = variations[iVar];
            if ((iDType == 0  && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) == 0))
              continue;

            g = (TGAE*) g_jet_trk_pt_bkg_syst[iDir][iCent][iVar]->Clone ("gtemp");
            SaveRelativeErrors (g, g, g_jet_trk_pt_bkg_syst[iDir][iCent][0], 100);
            ResetXErrors (g);
            ResetTGAEErrors (g);
            myDraw (g, varStyles[var].first, kDot, 0, varStyles[var].second, 4, "L");
            SaferDelete (&g);
          }

          myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
          myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{%s %i-%i%%}", iDType == 0 ? "ZDC" : "FCal", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
          myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} %s", GetJetPtStr (tag).Data (), directions[iDir] == "ns" ? "< #pi/8" : "> 7#pi/8"), 0.032);

          int count = 0;
          for (int iVar = 1; iVar < nVar; iVar++) {
            const TString var = variations[iVar];
            if ((iDType == 0  && dataVariations.count (var) > 0) || (iDType == 1 && mcVariations.count (var) > 0)) {
              myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[var].first, varStyles[var].second, varFullNames[var], 0.5, 0.014);
              count++;
            }
          }
          for (int iTotVar : totVars) {
            const TString totVar = totalVariations[iTotVar];
            myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 0.5, 0.014);
            count++;
          }
          
          c->SaveAs (Form ("%s/Plots/Systematics/PtCh/BkgdJetTaggedYield_pPb_%i-%iperc_%s_ptch_%s_%ssyst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", tag, iDType == 0 ? "":"mcBased_"));
        } // end loop over iCent
      }



      if (makeSigSystPlots) {
        {
          const char* canvasName = Form ("c_jet_trk_pt_%s_sig_%s_pp_syst", directions[iDir].Data (), iDType == 0 ? "data" : "mc");
          TCanvas* c = new TCanvas (canvasName, "", 800, 800);

          TH1D* h = nullptr;
          TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

          c->cd (); 
          c->SetLogx ();

          const float ymin = (iDType == 0 ? -maxDataSyst : -maxMCSyst);
          const float ymax = (iDType == 0 ?  maxDataSyst :  maxMCSyst);

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

          std::vector <int> totVars (0);
          if (iDType == 0) {
            totVars.push_back (0);
            totVars.push_back (1);
          }
          else 
            totVars.push_back (2);
          for (int iTotVar : totVars) {
            const TString totVar = totalVariations[iTotVar];
            g = ConvertErrorsToCentralValues (g_jet_trk_pt_ref_sig_systTot[iDir][iTotVar], true, 100);
            gup = (TGAE*) g->Clone ();
            gdown = (TGAE*) g->Clone ();
            FlipTGAE (gdown);
         
            myDraw (gup, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDraw (gdown, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDrawFill (gup, gdown, varStyles[totVar].first, 0.1);

            SaferDelete (&g);
            SaferDelete (&gup);
            SaferDelete (&gdown);
          }

          for (int iVar = 1; iVar < nVar; iVar++) {
            const TString var = variations[iVar];
            if ((iDType == 0  && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) == 0))
              continue;

            g = (TGAE*) g_jet_trk_pt_ref_sig_syst[iDir][iVar]->Clone ("gtemp");
            SaveRelativeErrors (g, g, g_jet_trk_pt_ref_sig_syst[iDir][0], 100);
            ResetXErrors (g);
            ResetTGAEErrors (g);
            myDraw (g, varStyles[var].first, kDot, 0, varStyles[var].second, 4, "L");
            SaferDelete (&g);
          }

          myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
          myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
          if (directions[iDir] == "ns")
            myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
          else if (directions[iDir] == "as")
            myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);

          int count = 0;
          for (int iVar = 1; iVar < nVar; iVar++) {
            const TString var = variations[iVar];
            if ((iDType == 0  && dataVariations.count (var) > 0) || (iDType == 1 && mcVariations.count (var) > 0)) {
              myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[var].first, varStyles[var].second, varFullNames[var], 0.5, 0.014);
              count++;
            }
          }
          for (int iTotVar : totVars) {
            const TString totVar = totalVariations[iTotVar];
            myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 0.5, 0.014);
            count++;
          }

          c->SaveAs (Form ("%s/Plots/Systematics/PtCh/SignalJetTaggedYield_pp_%s_ptch_%s_%ssyst.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : "awayside", tag, iDType == 0 ? "":"mcBased_"));
        }


        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          const char* canvasName = Form ("c_jet_trk_pt_%s_sig_%s_pPb_iCent%i_syst", directions[iDir].Data (), iDType == 0 ? "data" : "mc", iCent);
          TCanvas* c = new TCanvas (canvasName, "", 800, 800);

          TH1D* h = nullptr;
          TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

          c->cd (); 
          c->SetLogx ();
          //c->SetLogy ();

          const float ymin = (iDType == 0 ? -maxDataSyst : -maxMCSyst);
          const float ymax = (iDType == 0 ?  maxDataSyst :  maxMCSyst);

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

          std::vector <int> totVars (0);
          if (iDType == 0) {
            totVars.push_back (0);
            totVars.push_back (1);
          }
          else 
            totVars.push_back (2);
          for (int iTotVar : totVars) {
            const TString totVar = totalVariations[iTotVar];
            g = ConvertErrorsToCentralValues (g_jet_trk_pt_sig_systTot[iDir][iCent][iTotVar], true, 100);
            gup = (TGAE*) g->Clone ();
            gdown = (TGAE*) g->Clone ();
            FlipTGAE (gdown);
         
            myDraw (gup, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDraw (gdown, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDrawFill (gup, gdown, varStyles[totVar].first, 0.1);

            SaferDelete (&g);
            SaferDelete (&gup);
            SaferDelete (&gdown);
          }

          for (int iVar = 1; iVar < nVar; iVar++) {
            const TString var = variations[iVar];
            if ((iDType == 0  && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) == 0))
              continue;

            g = (TGAE*) g_jet_trk_pt_sig_syst[iDir][iCent][iVar]->Clone ("gtemp");
            SaveRelativeErrors (g, g, g_jet_trk_pt_sig_syst[iDir][iCent][0], 100);
            ResetXErrors (g);
            ResetTGAEErrors (g);
            myDraw (g, varStyles[var].first, kDot, 0, varStyles[var].second, 4, "L");
            SaferDelete (&g);
          }

          myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
          myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{%s %i-%i%%}", iDType == 0 ? "ZDC" : "FCal", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
          if (directions[iDir] == "ns")
            myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
          else if (directions[iDir] == "as")
            myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);

          int count = 0;
          for (int iVar = 1; iVar < nVar; iVar++) {
            const TString var = variations[iVar];
            if ((iDType == 0  && dataVariations.count (var) > 0) || (iDType == 1 && mcVariations.count (var) > 0)) {
              myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[var].first, varStyles[var].second, varFullNames[var], 0.5, 0.014);
              count++;
            }
          }
          for (int iTotVar : totVars) {
            const TString totVar = totalVariations[iTotVar];
            myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 0.5, 0.014);
            count++;
          }

          c->SaveAs (Form ("%s/Plots/Systematics/PtCh/SignalJetTaggedYield_pPb_%i-%iperc_%s_ptch_%s_%ssyst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", tag, iDType == 0 ? "":"mcBased_"));
        } // end loop over iCent
      }



      if (makeIpPbSystPlots) {
        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          const char* canvasName = Form ("c_jet_trk_pt_%s_iaa_%s_pPb_iCent%i_syst", directions[iDir].Data (), iDType == 0 ? "data" : "mc", iCent);
          TCanvas* c = new TCanvas (canvasName, "", 800, 800);

          TH1D* h = nullptr;
          TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

          c->cd (); 
          c->SetLogx ();
          //c->SetLogy ();

          const float ymin = (iDType == 0 ? -maxDataSyst : -maxMCSyst);
          const float ymax = (iDType == 0 ?  maxDataSyst :  maxMCSyst);

          h = (TH1D*) h_jet_trk_pt_iaa[iDType][iDir][iCent]->Clone ("h");
          h->Reset ();
          h->GetXaxis ()->SetMoreLogLabels ();
          h->GetYaxis ()->SetRangeUser (ymin, ymax);
          h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
          //h->GetXaxis ()->SetTitleSize (0.028);
          //h->GetXaxis ()->SetLabelSize (0.028);
          h->GetYaxis ()->SetTitle ("#delta I_{pPb} / I_{pPb} [%]");
          //h->GetYaxis ()->SetTitleSize (0.028);
          //h->GetYaxis ()->SetLabelSize (0.028);

          h->SetLineWidth (1);
          h->SetLineStyle (2);
          h->DrawCopy ("hist ][");
          SaferDelete (&h);

          std::vector <int> totVars (0);
          if (iDType == 0) {
            totVars.push_back (0);
            totVars.push_back (1);
          }
          else 
            totVars.push_back (2);
          for (int iTotVar : totVars) {
            const TString totVar = totalVariations[iTotVar];
            g = ConvertErrorsToCentralValues (g_jet_trk_pt_iaa_systTot[iDir][iCent][iTotVar], true, 100);
            gup = (TGAE*) g->Clone ();
            gdown = (TGAE*) g->Clone ();
            FlipTGAE (gdown);
         
            myDraw (gup, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDraw (gdown, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDrawFill (gup, gdown, varStyles[totVar].first, 0.1);

            SaferDelete (&g);
            SaferDelete (&gup);
            SaferDelete (&gdown);
          }

          for (int iVar = 1; iVar < nVar; iVar++) {
            const TString var = variations[iVar];
            if ((iDType == 0  && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) == 0))
              continue;

            g = (TGAE*) g_jet_trk_pt_iaa_syst[iDir][iCent][iVar]->Clone ("gtemp");
            SaveRelativeErrors (g, g, g_jet_trk_pt_iaa_syst[iDir][iCent][0], 100);
            ResetXErrors (g);
            ResetTGAEErrors (g);
            myDraw (g, varStyles[var].first, kDot, 0, varStyles[var].second, 4, "L");
            SaferDelete (&g);
          }

          myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
          myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
          myText (0.22, 0.81, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{%s %i-%i%%}", iDType == 0 ? "ZDC" : "FCal", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
          if (directions[iDir] == "ns")
            myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
          else if (directions[iDir] == "as")
            myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);

          int count = 0;
          for (int iVar = 1; iVar < nVar; iVar++) {
            const TString var = variations[iVar];
            if ((iDType == 0  && dataVariations.count (var) > 0) || (iDType == 1 && mcVariations.count (var) > 0)) {
              myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[var].first, varStyles[var].second, varFullNames[var], 0.5, 0.014);
              count++;
            }
          }
          for (int iTotVar : totVars) {
            const TString totVar = totalVariations[iTotVar];
            myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 0.5, 0.014);
            count++;
          }

          c->SaveAs (Form ("%s/Plots/Systematics/PtCh/JetTagged_IpPb_%i-%iperc_%s_ptch_%s_%ssyst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", tag, iDType == 0 ? "":"mcBased_"));
        } // end loop over iCent
      }

    } // end loop over iDir

  } // end loop over iDType



  {
    const int maxNSys = 6;

    for (int iDir : {0, 2}) {

      if (makeSigSystPlots) {
        {
          const char* canvasName = Form ("c_jet_trk_pt_%s_sig_data_pp_systTot", directions[iDir].Data ());
          TCanvas* c = new TCanvas (canvasName, "", 800, 800);

          TH1D* h = nullptr;
          TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

          c->cd (); 
          c->SetLogx ();

          const float ymin = -maxDataSyst;
          const float ymax =  maxDataSyst;

          h = (TH1D*) h_jet_trk_pt_ref_sig[0][iDir]->Clone ("h");
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

          for (int iTotVar : {0, 1, 2}) {
            const TString totVar = totalVariations[iTotVar];
            g = ConvertErrorsToCentralValues (g_jet_trk_pt_ref_sig_systTot[iDir][iTotVar], true, 100);
            gup = (TGAE*) g->Clone ();
            gdown = (TGAE*) g->Clone ();
            FlipTGAE (gdown);
         
            myDraw (gup, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDraw (gdown, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDrawFill (gup, gdown, varStyles[totVar].first, 0.1);

            SaferDelete (&g);
            SaferDelete (&gup);
            SaferDelete (&gdown);
          }

          g = ConvertErrorsToCentralValues (g_jet_trk_pt_ref_sig_syst[iDir][0], true, 100);
          myDraw (g, kBlack, kDot, 0, 1, 2, "L");
          FlipTGAE (g);
          myDraw (g, kBlack, kDot, 0, 1, 2, "L");
          SaferDelete (&g);

          myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
          myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
          if (directions[iDir] == "ns")
            myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
          else if (directions[iDir] == "as")
            myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);

          int count = 0;
          for (int iTotVar : {0, 1, 2}) {
            const TString totVar = totalVariations[iTotVar];
            myLineColorText (0.35+0.35*(count/maxNSys), 0.35-(count%maxNSys)*0.040, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 1.0, 0.028);
            count++;
          }
          myLineColorText (0.35+0.35*(count/maxNSys), 0.35-(count%maxNSys)*0.040, kBlack, 1, "#bf{Total unc.}", 1.0, 0.028);

          c->SaveAs (Form ("%s/Plots/Systematics/PtCh/SignalJetTaggedYield_pp_%s_ptch_%s_systTot.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : "awayside", tag));
        }


        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          const char* canvasName = Form ("c_jet_trk_pt_%s_sig_data_pPb_iCent%i_systTot", directions[iDir].Data (), iCent);
          TCanvas* c = new TCanvas (canvasName, "", 800, 800);

          TH1D* h = nullptr;
          TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

          c->cd (); 
          c->SetLogx ();
          //c->SetLogy ();

          const float ymin = -maxDataSyst;
          const float ymax =  maxDataSyst;

          h = (TH1D*) h_jet_trk_pt_sig[0][iDir][iCent]->Clone ("h");
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

          for (int iTotVar : {0, 1, 2}) {
            const TString totVar = totalVariations[iTotVar];
            g = ConvertErrorsToCentralValues (g_jet_trk_pt_sig_systTot[iDir][iCent][iTotVar], true, 100);
            gup = (TGAE*) g->Clone ();
            gdown = (TGAE*) g->Clone ();
            FlipTGAE (gdown);
         
            myDraw (gup, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDraw (gdown, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDrawFill (gup, gdown, varStyles[totVar].first, 0.1);

            SaferDelete (&g);
            SaferDelete (&gup);
            SaferDelete (&gdown);
          }

          g = ConvertErrorsToCentralValues (g_jet_trk_pt_sig_syst[iDir][iCent][0], true, 100);
          myDraw (g, kBlack, kDot, 0, 1, 2, "L");
          FlipTGAE (g);
          myDraw (g, kBlack, kDot, 0, 1, 2, "L");
          SaferDelete (&g);

          myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
          myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
          if (directions[iDir] == "ns")
            myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
          else if (directions[iDir] == "as")
            myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);

          int count = 0;
          for (int iTotVar : {0, 1, 2}) {
            const TString totVar = totalVariations[iTotVar];
            myLineColorText (0.35+0.35*(count/maxNSys), 0.35-(count%maxNSys)*0.040, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 1.0, 0.028);
            count++;
          }
          myLineColorText (0.35+0.35*(count/maxNSys), 0.35-(count%maxNSys)*0.040, kBlack, 1, "#bf{Total unc.}", 1.0, 0.028);

          c->SaveAs (Form ("%s/Plots/Systematics/PtCh/SignalJetTaggedYield_pPb_%i-%iperc_%s_ptch_%s_systTot.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", tag));
        } // end loop over iCent
      }



      if (makeIpPbSystPlots) {
        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          const char* canvasName = Form ("c_jet_trk_pt_%s_iaa_data_pPb_iCent%i_systTot", directions[iDir].Data (), iCent);
          TCanvas* c = new TCanvas (canvasName, "", 800, 800);

          TH1D* h = nullptr;
          TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

          c->cd (); 
          c->SetLogx ();
          //c->SetLogy ();

          const float ymin = -maxDataSyst;
          const float ymax =  maxDataSyst;

          h = (TH1D*) h_jet_trk_pt_iaa[0][iDir][iCent]->Clone ("h");
          h->Reset ();
          h->GetXaxis ()->SetMoreLogLabels ();
          h->GetYaxis ()->SetRangeUser (ymin, ymax);
          h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
          //h->GetXaxis ()->SetTitleSize (0.028);
          //h->GetXaxis ()->SetLabelSize (0.028);
          h->GetYaxis ()->SetTitle ("#delta I_{pPb} / I_{pPb} [%]");
          //h->GetYaxis ()->SetTitleSize (0.028);
          //h->GetYaxis ()->SetLabelSize (0.028);

          h->SetLineWidth (1);
          h->SetLineStyle (2);
          h->DrawCopy ("hist ][");
          SaferDelete (&h);

          for (int iTotVar : {0, 1, 2}) {
            const TString totVar = totalVariations[iTotVar];
            g = ConvertErrorsToCentralValues (g_jet_trk_pt_iaa_systTot[iDir][iCent][iTotVar], true, 100);
            gup = (TGAE*) g->Clone ();
            gdown = (TGAE*) g->Clone ();
            FlipTGAE (gdown);
         
            myDraw (gup, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDraw (gdown, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
            myDrawFill (gup, gdown, varStyles[totVar].first, 0.1);

            SaferDelete (&g);
            SaferDelete (&gup);
            SaferDelete (&gdown);
          }

          g = ConvertErrorsToCentralValues (g_jet_trk_pt_iaa_syst[iDir][iCent][0], true, 100);
          myDraw (g, kBlack, kDot, 0, 1, 2, "L");
          FlipTGAE (g);
          myDraw (g, kBlack, kDot, 0, 1, 2, "L");
          SaferDelete (&g);

          myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
          myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
          myText (0.22, 0.81, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}",  zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
          if (directions[iDir] == "ns")
            myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
          else if (directions[iDir] == "as")
            myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);

          int count = 0;
          for (int iTotVar : {0, 1, 2}) {
            const TString totVar = totalVariations[iTotVar];
            myLineColorText (0.35+0.35*(count/maxNSys), 0.35-(count%maxNSys)*0.040, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 1.0, 0.028);
            count++;
          }
          myLineColorText (0.35+0.35*(count/maxNSys), 0.35-(count%maxNSys)*0.040, kBlack, 1, "#bf{Total unc.}", 1.0, 0.028);

          c->SaveAs (Form ("%s/Plots/Systematics/PtCh/JetTagged_IpPb_%i-%iperc_%s_ptch_%s_systTot.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", tag));
        } // end loop over iCent
      }

    } // end loop over iDir

  }

}


#endif
