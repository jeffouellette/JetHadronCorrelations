#ifndef __JetHadronCorrelatorPlotNIters_C__
#define __JetHadronCorrelatorPlotNIters_C__

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


void PlotNIters (const char* inFileTag) {

  TFile* inFile = nullptr;

  const int nItersMin = 1;
  const int nItersMax = 20;
  const double* nItersVals = linspace (nItersMin, nItersMax, nItersMax-nItersMin);

  TH1D****  h_jetInt_trk_pt_ref_sig   = Get3DArray <TH1D*> (2, 2, nDir);
  TH1D***** h_jetInt_trk_pt_sig       = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);

  TH1D****  h_jetInt_trk_pt_ref_unf   = Get3DArray <TH1D*> (2, 2, nDir);
  TH1D***** h_jetInt_trk_pt_unf       = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);

  TH1D****  h_jet_trk_pt_ref_unf_nIters = Get3DArray <TH1D*> (nPtJBins, nDir, nItersMax-nItersMin+2);
  TH1D***** h_jet_trk_pt_unf_nIters     = Get4DArray <TH1D*> (nPtJBins, nDir, nZdcCentBins+1, nItersMax-nItersMin+2);

  TH1D***** h_jetInt_trk_pt_iaa       = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);
  TH1D***** h_jetInt_trk_pt_iaaNoUnf  = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);


  TGraph*     g_jet_pt_ref_unfIterUnc     = nullptr;                                // sums of iterations uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_unfIterUnc         = Get1DArray <TGraph*> (nZdcCentBins+1);  // sums of iterations uncertainties as a function of nIter -- data only
  TGraph*     g_jet_pt_ref_unfSumUnc      = nullptr;                                // sums of statistical uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_unfSumUnc          = Get1DArray <TGraph*> (nZdcCentBins+1);  // sums of statistical uncertainties as a function of nIter -- data only
  TGraph*     g_jet_pt_ref_unfTotUnc      = nullptr;                                // sums of total uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_unfTotUnc          = Get1DArray <TGraph*> (nZdcCentBins+1);  // sums of total uncertainties as a function of nIter -- data only

  TGraph*     g_jet_pt_ref_unfIterRelUnc  = nullptr;                                // sums of iterations relative uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_unfIterRelUnc      = Get1DArray <TGraph*> (nZdcCentBins+1);  // sums of iterations relative uncertainties as a function of nIter -- data only
  TGraph*     g_jet_pt_ref_unfSumRelUnc   = nullptr;                                // sums of statistical relative uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_unfSumRelUnc       = Get1DArray <TGraph*> (nZdcCentBins+1);  // sums of statistical relative uncertainties as a function of nIter -- data only
  TGraph*     g_jet_pt_ref_unfTotRelUnc   = nullptr;                                // sums of total relative uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_unfTotRelUnc       = Get1DArray <TGraph*> (nZdcCentBins+1);  // sums of total relative uncertainties as a function of nIter -- data only


  TGraph***   g_jetInt_trk_pt_ref_unfIterUnc    = Get2DArray <TGraph*> (2, nDir);                           // sums of iterations uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfIterUnc        = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of iterations uncertainties as a function of nIter -- data only
  TGraph***   g_jetInt_trk_pt_ref_unfSumUnc     = Get2DArray <TGraph*> (2, nDir);                           // sums of statistical uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfSumUnc         = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of statistical uncertainties as a function of nIter -- data only
  TGraph***   g_jetInt_trk_pt_ref_unfTotUnc     = Get2DArray <TGraph*> (2, nDir);                           // sums of total uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfTotUnc         = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of total uncertainties as a function of nIter -- data only

  TGraph***   g_jetInt_trk_pt_ref_unfIterRelUnc = Get2DArray <TGraph*> (2, nDir);                           // sums of iterations relative uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfIterRelUnc     = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of iterations relative uncertainties as a function of nIter -- data only
  TGraph***   g_jetInt_trk_pt_ref_unfSumRelUnc  = Get2DArray <TGraph*> (2, nDir);                           // sums of statistical relative uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfSumRelUnc      = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of statistical relative uncertainties as a function of nIter -- data only
  TGraph***   g_jetInt_trk_pt_ref_unfTotRelUnc  = Get2DArray <TGraph*> (2, nDir);                           // sums of total relative uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfTotRelUnc      = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of total relative uncertainties as a function of nIter -- data only


  {
    TString inFileName = inFileTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/ProcessCorrelations_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    for (int iDType = 0; iDType < 2; iDType++) {

      const TString dType = (iDType == 0 ? "data" : "mc");

      for (short iPtJInt : {0, 1}) {

        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

        for (int iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_sig_%s_%s_Nominal",  dir.Data (), dType.Data (), pTJInt.Data ()));
          h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_unf_%s_%s_Nominal",  dir.Data (), dType.Data (), pTJInt.Data ()));

        } // end loop over iDir


        for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

          for (int iDir = 0; iDir < nDir; iDir++) {

            const TString dir = directions[iDir];

            h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent]       = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_sig_%s_%s_%s_Nominal", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));
            h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent]       = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_unf_%s_%s_%s_Nominal", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));
            h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent]       = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_iaa_%s_%s_%s_Nominal",     dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));
            h_jetInt_trk_pt_iaaNoUnf[iDType][iPtJInt][iDir][iCent]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_iaaNoUnf_%s_%s_%s_Nominal",  dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));

          } // end loop over iDir

        } // end loop over iCent

      } // end loop over iPtJInt

    } // end loop over iDType


    g_jet_pt_ref_unfSumUnc      = (TGraph*) inFile->Get ("g_jet_pt_ref_unfSumUnc");
    g_jet_pt_ref_unfIterUnc     = (TGraph*) inFile->Get ("g_jet_pt_ref_unfIterUnc");
    g_jet_pt_ref_unfTotUnc      = (TGraph*) inFile->Get ("g_jet_pt_ref_unfTotUnc");
    g_jet_pt_ref_unfSumRelUnc   = (TGraph*) inFile->Get ("g_jet_pt_ref_unfSumRelUnc");
    g_jet_pt_ref_unfIterRelUnc  = (TGraph*) inFile->Get ("g_jet_pt_ref_unfIterRelUnc");
    g_jet_pt_ref_unfTotRelUnc   = (TGraph*) inFile->Get ("g_jet_pt_ref_unfTotRelUnc");

    for (short iDir = 0; iDir < nDir; iDir++) {

      const TString dir = directions[iDir];

      for (short iPtJInt : {0, 1}) {

        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

        for (short iIter = 0; iIter < nItersMax-nItersMin+2; iIter++) {
  
          const short nIters = (short) nItersVals[iIter];

          h_jet_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_ref_unf_data_%s_Nominal_nIters%i", dir.Data (), pTJInt.Data (), nIters));

        } // end loop over iIter

        g_jetInt_trk_pt_ref_unfSumUnc[iPtJInt][iDir]      = (TGraph*) inFile->Get (Form ("g_jetInt_trk_pt_ref_unfSumUnc_%s_%s", dir.Data (), pTJInt.Data ()));
        g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir]     = (TGraph*) inFile->Get (Form ("g_jetInt_trk_pt_ref_unfIterUnc_%s_%s", dir.Data (), pTJInt.Data ()));
        g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir]      = (TGraph*) inFile->Get (Form ("g_jetInt_trk_pt_ref_unfTotUnc_%s_%s", dir.Data (), pTJInt.Data ()));
        g_jetInt_trk_pt_ref_unfSumRelUnc[iPtJInt][iDir]   = (TGraph*) inFile->Get (Form ("g_jetInt_trk_pt_ref_unfSumRelUnc_%s_%s", dir.Data (), pTJInt.Data ()));
        g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir]  = (TGraph*) inFile->Get (Form ("g_jetInt_trk_pt_ref_unfIterRelUnc_%s_%s", dir.Data (), pTJInt.Data ()));
        g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir]   = (TGraph*) inFile->Get (Form ("g_jetInt_trk_pt_ref_unfTotRelUnc_%s_%s", dir.Data (), pTJInt.Data ()));

      } // end loop over iPtJInt

    } // end loop over iDir

    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

      g_jet_pt_unfSumUnc[iCent]     = (TGraph*) inFile->Get (Form ("g_jet_pt_pPb_%s_unfSumUnc", cent.Data ()));
      g_jet_pt_unfIterUnc[iCent]    = (TGraph*) inFile->Get (Form ("g_jet_pt_pPb_%s_unfIterUnc", cent.Data ()));
      g_jet_pt_unfTotUnc[iCent]     = (TGraph*) inFile->Get (Form ("g_jet_pt_pPb_%s_unfTotUnc", cent.Data ()));
      g_jet_pt_unfSumRelUnc[iCent]  = (TGraph*) inFile->Get (Form ("g_jet_pt_pPb_%s_unfSumRelUnc", cent.Data ()));
      g_jet_pt_unfIterRelUnc[iCent] = (TGraph*) inFile->Get (Form ("g_jet_pt_pPb_%s_unfIterRelUnc", cent.Data ()));
      g_jet_pt_unfTotRelUnc[iCent]  = (TGraph*) inFile->Get (Form ("g_jet_pt_pPb_%s_unfTotRelUnc", cent.Data ()));

      for (short iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          for (short iIter = 0; iIter < nItersMax-nItersMin+2; iIter++) {
    
            const short nIters = (short) nItersVals[iIter];
  
            h_jet_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_pPb_unf_%s_data_%s_Nominal_nIters%i", dir.Data (), cent.Data (), pTJInt.Data (), nIters));
  
          } // end loop over iIter

          g_jetInt_trk_pt_unfSumUnc[iPtJInt][iDir][iCent]     = (TGraph*) inFile->Get (Form ("g_jetInt_trk_pt_%s_unfSumUnc_%s_%s", dir.Data (), cent.Data (), pTJInt.Data ()));
          g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent]    = (TGraph*) inFile->Get (Form ("g_jetInt_trk_pt_%s_unfIterUnc_%s_%s", dir.Data (), cent.Data (), pTJInt.Data ()));
          g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent]     = (TGraph*) inFile->Get (Form ("g_jetInt_trk_pt_%s_unfTotUnc_%s_%s", dir.Data (), cent.Data (), pTJInt.Data ()));
          g_jetInt_trk_pt_unfSumRelUnc[iPtJInt][iDir][iCent]  = (TGraph*) inFile->Get (Form ("g_jetInt_trk_pt_%s_unfSumRelUnc_%s_%s", dir.Data (), cent.Data (), pTJInt.Data ()));
          g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent] = (TGraph*) inFile->Get (Form ("g_jetInt_trk_pt_%s_unfIterRelUnc_%s_%s", dir.Data (), cent.Data (), pTJInt.Data ()));
          g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent]  = (TGraph*) inFile->Get (Form ("g_jetInt_trk_pt_%s_unfTotRelUnc_%s_%s", dir.Data (), cent.Data (), pTJInt.Data ()));

        } // end loop over iPtJInt

      } // end loop over iDir

    } // end loop over iCent

  }



  l->SetLineStyle (7);
  l->SetLineWidth (1);
  l->SetLineColor (kBlack);



  {
    const char* canvasName = "c_jet_pt_unfUncs";
    TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
    c->Divide (4, 2);

    TGraph* g = nullptr;

    TLine* l = new TLine ();
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kBlack);

    double ymax = 0, ymin = 0;
    {
      c->cd (7);

      gPad->SetLogy ();

      double x, y, xmin;
      ymax = 0, ymin = DBL_MAX;
      g = g_jet_pt_ref_unfTotUnc;
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) {
          xmin = x; ymin = y;
        }
        if (ymax > 0 && y > 2. * ymax) continue;
        ymax = std::fmax (y, ymax);
      }
      ymax *= 5;
      
      g = g_jet_pt_ref_unfSumUnc;
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      g = g_jet_pt_ref_unfIterUnc;
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      ymin *= 0.2;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{J}^{2}}", 1, 0, 20);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      g = g_jet_pt_ref_unfTotUnc;
      g->SetMarkerColor (kBlack);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_ref_unfSumUnc;
      g->SetMarkerColor (kBlue);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_ref_unfIterUnc;
      g->SetMarkerColor (kRed);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
      l->DrawLine (xmin, ymin, xmin, ymax/5);
    }


    for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {
      c->cd (nZdcCentBins+1-iCent);

      gPad->SetLogy ();

      double x, y, xmin;
      ymax = 0, ymin = DBL_MAX;
      g = g_jet_pt_unfTotUnc[iCent];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) {
          xmin = x; ymin = y;
        }
        if (ymax > 0 && y > 2. * ymax) continue;
        ymax = std::fmax (y, ymax);
      }
      ymax *= 5;

      g = g_jet_pt_unfSumUnc[iCent];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      g = g_jet_pt_unfIterUnc[iCent];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      ymin *= 0.2;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{J}^{2}}", 1, 0, 20);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      g = g_jet_pt_unfTotUnc[iCent];
      g->SetMarkerColor (kBlack);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_unfSumUnc[iCent];
      g->SetMarkerColor (kBlue);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_unfIterUnc[iCent];
      g->SetMarkerColor (kRed);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      if (iCent < nZdcCentBins)
        myText (0.2, 0.865, kBlack, Form ("#bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
      else
        myText (0.2, 0.865, kBlack, "#bf{All centralities}", 0.05);
      l->DrawLine (xmin, ymin, xmin, ymax/5);

    } // end loop over iCent

    c->cd (8);
    myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
    myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.66, kBlack, "#it{p}_{T}^{jet} [GeV] #in (30, 300)", 0.065);
    myLineText2 (0.15, 0.56, kBlack,        kOpenCircle, "Stat. unc. #oplus Iter. diff.", 1.2, 0.06);
    myLineText2 (0.15, 0.48, kBlue,         kOpenCircle, "Stat. unc. only", 1.2, 0.06);
    myLineText2 (0.15, 0.40, kRed,          kOpenCircle, "Iter. diff. only, i #rightarrow i+1", 1.2, 0.06);

    c->SaveAs (Form ("%s/Plots/PtCh/UnfUncs_Summary_JetSpectra.pdf", workPath.Data ()));
    SaferDelete (&l);
  }




  {
    const char* canvasName = "c_jet_pt_unfRelUncs";
    TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
    c->Divide (4, 2);

    TGraph* g = nullptr;

    TLine* l = new TLine ();
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kBlack);

    double ymax = 0, ymin = 0;
    {
      c->cd (7);

      //gPad->SetLogy ();

      double x, y, xmin;
      ymax = 0, ymin = DBL_MAX;
      g = g_jet_pt_ref_unfTotRelUnc;
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) {
          xmin = x; ymin = y;
        }
        if (ymax > 0 && y > 2. * ymax) continue;
        ymax = std::fmax (y, ymax);
      }
      ymax *= 1.2;

      ymin = 0;  
      //g = g_jet_pt_ref_unfSumRelUnc;
      //for (short i = 0; i < g->GetN (); i++) {
      //  g->GetPoint (i, x, y);
      //  if (x != 0 && y < ymin && y > 0) ymin = y;
      //}
      //g = g_jet_pt_ref_unfIterRelUnc;
      //for (short i = 0; i < g->GetN (); i++) {
      //  g->GetPoint (i, x, y);
      //  if (x != 0 && y < ymin && y > 0) ymin = y;
      //}
      //ymin *= 0.2;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{J}^{2} / N_{J}}", 1, 0, 20);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      g = g_jet_pt_ref_unfTotRelUnc;
      g->SetMarkerColor (kBlack);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_ref_unfSumRelUnc;
      g->SetMarkerColor (kBlue);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_ref_unfIterRelUnc;
      g->SetMarkerColor (kRed);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
      l->DrawLine (xmin, ymin, xmin, ymax/1.2);
    }


    for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {
      c->cd (nZdcCentBins+1-iCent);

      //gPad->SetLogy ();

      double x, y, xmin;
      ymax = 0, ymin = DBL_MAX;
      g = g_jet_pt_unfTotRelUnc[iCent];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) {
          xmin = x; ymin = y;
        }
        if (ymax > 0 && y > 2. * ymax) continue;
        ymax = std::fmax (y, ymax);
      }
      ymax *= 1.2;

      ymin = 0;  
      //g = g_jet_pt_unfSumRelUnc[iCent];
      //for (short i = 0; i < g->GetN (); i++) {
      //  g->GetPoint (i, x, y);
      //  if (x != 0 && y < ymin && y > 0) ymin = y;
      //}
      //g = g_jet_pt_unfIterRelUnc[iCent];
      //for (short i = 0; i < g->GetN (); i++) {
      //  g->GetPoint (i, x, y);
      //  if (x != 0 && y < ymin && y > 0) ymin = y;
      //}
      //ymin *= 0.2;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{J}^{2} / N_{J}}", 1, 0, 20);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      g = g_jet_pt_unfTotRelUnc[iCent];
      g->SetMarkerColor (kBlack);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_unfSumRelUnc[iCent];
      g->SetMarkerColor (kBlue);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_unfIterRelUnc[iCent];
      g->SetMarkerColor (kRed);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      if (iCent < nZdcCentBins)
        myText (0.2, 0.865, kBlack, Form ("#bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
      else
        myText (0.2, 0.865, kBlack, "#bf{All centralities}", 0.05);
      l->DrawLine (xmin, ymin, xmin, ymax/1.2);

    } // end loop over iCent

    c->cd (8);
    myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
    myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.66, kBlack, "#it{p}_{T}^{jet} [GeV] #in (30, 300)", 0.065);
    myLineText2 (0.15, 0.56, kBlack,        kOpenCircle, "Stat. unc. #oplus Iter. diff.", 1.2, 0.06);
    myLineText2 (0.15, 0.48, kBlue,         kOpenCircle, "Stat. unc. only", 1.2, 0.06);
    myLineText2 (0.15, 0.40, kRed,          kOpenCircle, "Iter. diff. only, i #rightarrow i+1", 1.2, 0.06);

    c->SaveAs (Form ("%s/Plots/PtCh/UnfRelUncs_Summary_JetSpectra.pdf", workPath.Data ()));
    SaferDelete (&l);
  }




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (int iDir : {0, 1, 2}) {

      const char* canvasName = Form ("c_jetInt_trk_pt_unfUncs_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
      c->Divide (4, 2);

      TGraph* g = nullptr;

      TLine* l = new TLine ();
      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kBlack);

      double ymax = 0, ymin = 0;
      {
        c->cd (7);

        gPad->SetLogy ();

        double x, y, xmin;
        ymax = 0, ymin = DBL_MAX;
        g = g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) {
            xmin = x; ymin = y;
          }
          if (ymax > 0 && y > 2. * ymax) continue;
          ymax = std::fmax (y, ymax);
        }
        ymax *= 5;

        g = g_jetInt_trk_pt_ref_unfSumUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        g = g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        ymin *= 0.2;

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{Y}^{2}}", 1, 0, 20);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        g = g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir];
        g->SetMarkerColor (kBlack);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_ref_unfSumUnc[iPtJInt][iDir];
        g->SetMarkerColor (kBlue);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir];
        g->SetMarkerColor (kRed);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
        l->DrawLine (xmin, ymin, xmin, ymax/5);
      }


      for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        gPad->SetLogy ();

        double x, y, xmin;
        ymax = 0, ymin = DBL_MAX;
        g = g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) {
            xmin = x; ymin = y;
          }
          if (ymax > 0 && y > 2. * ymax) continue;
          ymax = std::fmax (y, ymax);
        }
        ymax *= 5;

        g = g_jetInt_trk_pt_unfSumUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        g = g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        ymin *= 0.2;

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{Y}^{2}}", 1, 0, 20);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        g = g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kBlack);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_unfSumUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kBlue);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kRed);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, kBlack, Form ("#bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, kBlack, "#bf{All centralities}", 0.05);
        l->DrawLine (xmin, ymin, xmin, ymax/5);

      } // end loop over iCent

      c->cd (8);
      myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
      myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, 300)", iPtJInt == 0 ? 30 : 60), 0.065);
      myText (0.1, 0.57, kBlack, Form ("#it{p}_{T}^{ch} [GeV] #in (5, %i)", iPtJInt == 0 ? 40 : 75), 0.065);
      myText (0.1, 0.48, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.065);
      myLineText2 (0.15, 0.40, kBlack,        kOpenCircle, "Stat. unc. #oplus Iter. diff.", 1.2, 0.06);
      myLineText2 (0.15, 0.31, kBlue,         kOpenCircle, "Stat. unc. only", 1.2, 0.06);
      myLineText2 (0.15, 0.22, kRed,          kOpenCircle, "Iter. diff. only, i #rightarrow i+1", 1.2, 0.06);

      c->SaveAs (Form ("%s/Plots/PtCh/UnfUncs_Summary_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
      SaferDelete (&l);

    } // end loop over iDir

  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (int iDir : {0, 1, 2}) {

      const char* canvasName = Form ("c_jetInt_trk_pt_unfRelUncs_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
      c->Divide (4, 2);

      TGraph* g = nullptr;

      TLine* l = new TLine ();
      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kBlack);

      double ymax = 0, ymin = 0;
      {
        c->cd (7);

        //gPad->SetLogy ();

        double x, y, xmin;
        ymax = 0, ymin = DBL_MAX;
        g = g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) {
            xmin = x; ymin = y;
          }
          if (ymax > 0 && y > 2. * ymax) continue;
          ymax = std::fmax (y, ymax);
        }
        ymax *= 1.2;

        ymin = 0;  
        //g = g_jetInt_trk_pt_ref_unfSumRelUnc[iPtJInt][iDir];
        //for (short i = 0; i < g->GetN (); i++) {
        //  g->GetPoint (i, x, y);
        //  if (x != 0 && y < ymin && y > 0) ymin = y;
        //}
        //g = g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir];
        //for (short i = 0; i < g->GetN (); i++) {
        //  g->GetPoint (i, x, y);
        //  if (x != 0 && y < ymin && y > 0) ymin = y;
        //}
        //ymin *= 0.2;

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{Y}^{2} / Y}", 1, 0, 20);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        g = g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir];
        g->SetMarkerColor (kBlack);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_ref_unfSumRelUnc[iPtJInt][iDir];
        g->SetMarkerColor (kBlue);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir];
        g->SetMarkerColor (kRed);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
        l->DrawLine (xmin, ymin, xmin, ymax/1.2);
      }


      for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        //gPad->SetLogy ();

        double x, y, xmin;
        ymax = 0, ymin = DBL_MAX;
        g = g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) {
            xmin = x; ymin = y;
          }
          if (ymax > 0 && y > 2. * ymax) continue;
          ymax = std::fmax (y, ymax);
        }
        ymax *= 1.2;

        ymin = 0;  
        //g = g_jetInt_trk_pt_unfSumRelUnc[iPtJInt][iDir][iCent];
        //for (short i = 0; i < g->GetN (); i++) {
        //  g->GetPoint (i, x, y);
        //  if (x != 0 && y < ymin && y > 0) ymin = y;
        //}
        //g = g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent];
        //for (short i = 0; i < g->GetN (); i++) {
        //  g->GetPoint (i, x, y);
        //  if (x != 0 && y < ymin && y > 0) ymin = y;
        //}
        //ymin *= 0.2;

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{Y}^{2} / Y}", 1, 0, 20);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        g = g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kBlack);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_unfSumRelUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kBlue);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kRed);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, kBlack, Form ("#bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, kBlack, "#bf{All centralities}", 0.05);
        l->DrawLine (xmin, ymin, xmin, ymax/1.2);

      } // end loop over iCent

      c->cd (8);
      myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
      myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, 300)", iPtJInt == 0 ? 30 : 60), 0.065);
      myText (0.1, 0.57, kBlack, Form ("#it{p}_{T}^{ch} [GeV] #in (5, %i)", iPtJInt == 0 ? 40 : 75), 0.065);
      myText (0.1, 0.48, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.065);
      myLineText2 (0.15, 0.40, kBlack,        kOpenCircle, "Stat. unc. #oplus Iter. diff.", 1.2, 0.06);
      myLineText2 (0.15, 0.31, kBlue,         kOpenCircle, "Stat. unc. only", 1.2, 0.06);
      myLineText2 (0.15, 0.22, kRed,          kOpenCircle, "Iter. diff. only, i #rightarrow i+1", 1.2, 0.06);

      c->SaveAs (Form ("%s/Plots/PtCh/UnfRelUncs_Summary_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
      SaferDelete (&l);

    } // end loop over iDir

  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (int iDir : {0, 1, 2}) {

      const char* canvasName = Form ("c_jetInt_trk_pt_sigVsUnf_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
      c->Divide (4, 2);

      TGAE* g = nullptr;

      {
        c->cd (7);

        gPad->SetLogx ();

        TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];Unfold / No unfold", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 4 : 2)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser (0.74, 1.4);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);


        for (short iIter : {0, 1, 2, 3, 4, 5}) {
          h = h_jet_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter];
          g = make_graph (h);
          ScaleGraph (g, h_jetInt_trk_pt_ref_sig[0][iPtJInt][iDir]);
          ResetXErrors (g);
          if (iDir == 1)
            TrimGraph (g, 0, 10);
          myDraw (g, colorfulColors[0], 24+iIter, 1.0, 1, 2, "P", false);
          SaferDelete (&g);
        }

        myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
      }


      for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        gPad->SetLogx ();

        TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];Unfold / No unfold", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 4 : 2)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser (0.74, 1.4);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        for (short iIter : {0, 1, 2, 3, 4, 5}) {
          h = h_jet_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter];
          g = make_graph (h);
          ScaleGraph (g, h_jetInt_trk_pt_sig[0][iPtJInt][iDir][iCent]);
          ResetXErrors (g);
          if (iDir == 1)
            TrimGraph (g, 0, 10);
          myDraw (g, colorfulColors[nZdcCentBins-iCent+1], 24+iIter, 1.0, 1, 2, "P", false);
          SaferDelete (&g);
        }

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, colorfulColors[nZdcCentBins-iCent+1], Form ("#bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, colorfulColors[nZdcCentBins-iCent+1], "#bf{All centralities}", 0.05);

      } // end loop over iCent

      c->cd (8);
      myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
      myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, 300)", iPtJInt == 0 ? 30 : 60), 0.065);
      myText (0.1, 0.57, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.065);

      c->SaveAs (Form ("%s/Plots/PtCh/UnfComp_Summary_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));

    } // end loop over iDir

  } // end loop over iPtJInt


}


#endif
