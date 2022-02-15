#ifndef __JetHadronCorrelatorPlotJets_C__
#define __JetHadronCorrelatorPlotJets_C__

#include <iostream>
#include <math.h>

#include <TColor.h>
#include <TLine.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
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
TLatex* tl = new TLatex ();


TString GetSamp (const short iDType, const short iSamp) {
  if (iDType == 0)
    return "AllTrigs";
  switch (iSamp) {
    case 0: return "JZ0";
    case 1: return "JZ1";
    case 2: return "JZ2";
    case 3: return "JZ3";
    case 4: return "JZ123";
    case 5: return "JZ0123";
  }
  return "???";
}



void PlotJets (const char* tag, const char* inFileTag) {

  const TString var = variations[0];

  const short nSamps = 6; // JZ0, 1, 2, 3, JZ1-3, JZ0-3

  TFile* inFile = nullptr;

  TH1D***  h_evt_counts_ref     = Get2DArray <TH1D*> (2, nSamps);
  TH1D**** h_evt_counts         = Get3DArray <TH1D*> (2, nZdcCentBins+1, nSamps);

  TH1D*** h_jet_pt_ref          = Get2DArray <TH1D*> (3, nSamps);

  TH1D**** h_jet_pt             = Get3DArray <TH1D*> (3, nZdcCentBins+1, nSamps);

  TH1D**** h_jet_pt_ratio       = Get3DArray <TH1D*> (2, nZdcCentBins+1, nSamps);

  TH1D***   h_jet_pt_datamc_ratio_ref  = Get2DArray <TH1D*> (2, nSamps);
  TH1D****  h_jet_pt_datamc_ratio      = Get3DArray <TH1D*> (2, nZdcCentBins+1, nSamps);

  TF1***    f_jet_pt_datamc_ratio_ref  = Get2DArray <TF1*> (2, nSamps);
  TF1****   f_jet_pt_datamc_ratio      = Get3DArray <TF1*> (2, nZdcCentBins+1, nSamps);

  //TH2D****  h2_jet_eta_phi_ref   = Get3DArray <TH2D*> (2, nPtJBins, nSamps);
  //TH2D***** h2_jet_eta_phi       = Get4DArray <TH2D*> (2, nPtJBins, nZdcCentBins+1, nSamps);


  {
    TString inFileName = inFileTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/ProcessJets_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    for (short iDType = 0; iDType < 2; iDType++) {

      const TString dType = (iDType == 0 ? "data" : "mc");

      for (short iSamp = 0; iSamp < nSamps; iSamp++) {

        if (iDType == 0 && iSamp > 0)
          continue;

        const TString samp = GetSamp (iDType, iSamp);

        h_evt_counts_ref[iDType][iSamp] = (TH1D*) inFile->Get (Form ("h_evt_counts_ref_%s_%s",    dType.Data (), samp.Data ()));

        h_jet_pt_ref[iDType][iSamp]     = (TH1D*) inFile->Get (Form ("h_jet_pt_ref_%s_%s",        dType.Data (), samp.Data ()));

        //for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        //  const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

        //  h2_jet_eta_phi_ref[iDType][iPtJ][iSamp]  = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_ref_%s_%s", dType.Data (), samp.Data ()));

        //} // end loop over iPtJ
  

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));
  
          h_evt_counts[iDType][iCent][iSamp]    = (TH1D*) inFile->Get (Form ("h_evt_counts_pPb_%s_%s_%s",   cent.Data (), dType.Data (), samp.Data ()));

          h_jet_pt[iDType][iCent][iSamp]        = (TH1D*) inFile->Get (Form ("h_jet_pt_pPb_%s_%s_%s",       cent.Data (), dType.Data (), samp.Data ()));
          h_jet_pt_ratio[iDType][iCent][iSamp]  = (TH1D*) inFile->Get (Form ("h_jet_pt_ratio_%s_%s_%s",     cent.Data (), dType.Data (), samp.Data ()));

          //for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          //  const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

          //  h2_jet_eta_phi[iDType][iPtJ][iCent][iSamp]  = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_pPb_%s_%s_%s_%s", cent.Data (), dType.Data (), pTJ.Data (), samp.Data ()));

          //} // end loop over iPtJ
  
        } // end loop over iCent

      } // end loop over iSamp

    } // end loop over iDType


    for (short iSamp = 0; iSamp < nSamps; iSamp++) {

      const TString samp = GetSamp (1, iSamp);

      h_jet_pt_ref[2][iSamp]              = (TH1D*) inFile->Get (Form ("h_jet_pt_ref_mcScaled_%s",            samp.Data ()));

      h_jet_pt_datamc_ratio_ref[0][iSamp] = (TH1D*) inFile->Get (Form ("h_jet_pt_datamc_ratio_ref_%s",        samp.Data ()));
      f_jet_pt_datamc_ratio_ref[0][iSamp] = (TF1*)  inFile->Get (Form ("f_jet_pt_datamc_ratio_ref_%s",        samp.Data ()));

      h_jet_pt_datamc_ratio_ref[1][iSamp] = (TH1D*) inFile->Get (Form ("h_jet_pt_datamcScaled_ratio_ref_%s",  samp.Data ()));
      f_jet_pt_datamc_ratio_ref[1][iSamp] = (TF1*)  inFile->Get (Form ("f_jet_pt_datamcScaled_ratio_ref_%s",  samp.Data ()));

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

        h_jet_pt[2][iCent][iSamp]               = (TH1D*) inFile->Get (Form ("h_jet_pt_pPb_%s_mcScaled_%s",       cent.Data (), samp.Data ()));

        h_jet_pt_datamc_ratio[0][iCent][iSamp]  = (TH1D*) inFile->Get (Form ("h_jet_pt_datamc_ratio_%s_%s",       cent.Data (), samp.Data ()));
        f_jet_pt_datamc_ratio[0][iCent][iSamp]  = (TF1*)  inFile->Get (Form ("f_jet_pt_datamc_ratio_%s_%s",       cent.Data (), samp.Data ()));

        h_jet_pt_datamc_ratio[1][iCent][iSamp]  = (TH1D*) inFile->Get (Form ("h_jet_pt_datamcScaled_ratio_%s_%s", cent.Data (), samp.Data ()));
        f_jet_pt_datamc_ratio[1][iCent][iSamp]  = (TF1*)  inFile->Get (Form ("f_jet_pt_datamcScaled_ratio_%s_%s", cent.Data (), samp.Data ()));

      } // end loop over iCent

    } // end loop over iSamp

  }


  //for (float trigpt : {30., 60.}) {
//  float maxpt = (trigpt == 30. ? 60. : 300.);
  //  float* njet = new float[nZdcCentBins+2];

  //  std::cout << "---------------" << std::endl << "JETS IN DATA > " << trigpt << " GeV" << std::endl << "---------------" << std::endl;

  //  njet[0] = 0;
  //  for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
  //    if (0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) >= trigpt && 0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) < maxpt)
  //      njet[0] += h_jet_counts_ref[0][iPtJ][0]->GetBinContent (1);
  //  } // end loop over iPtJ
  //  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
  //    njet[iCent+1] = 0;
  //    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
  //      if (0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) >= trigpt && 0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) < maxpt)
  //        njet[iCent+1] += h_jet_counts[0][iPtJ][iCent][0]->GetBinContent (1);
  //    } // end loop over iPtJ
  //  } // end loop over iCent

  //  std::cout << "Number of pp jets: " << njet[0] << std::endl;
  //  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++)
  //    std::cout << "Number of p+Pb " << zdcCentPercs[iCent+1] << "\%-" << zdcCentPercs[iCent] << "\% jets: " << njet[iCent+1] << std::endl;

  //  std::cout << "Formatted for latex:" << std::endl;
  //  std::cout << (int) njet[0];
  //  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++)
  //    std::cout << " & " << (int) njet[iCent+1];
  //  std::cout << std::endl << std::endl;


  //  float integral = 0;
  //  for (short iX = h_jet_pt_ref[0][0]->FindBin (trigpt); iX <= h_jet_pt_ref[0][0]->GetNbinsX (); iX++)
  //    integral += h_jet_pt_ref[0][0]->GetBinContent (iX) * h_jet_pt_ref[0][0]->GetBinWidth (iX);
  //  std::cout << "Average number of jets per trigger jet in pp: " << integral << std::endl;
  //  integral = 0;
  //  for (short iX = h_jet_pt[0][nZdcCentBins][0]->FindBin (trigpt); iX <= h_jet_pt[0][nZdcCentBins][0]->GetNbinsX (); iX++)
  //    integral += h_jet_pt[0][nZdcCentBins][0]->GetBinContent (iX) * h_jet_pt[0][nZdcCentBins][0]->GetBinWidth (iX);
  //  std::cout << "Average number of jets per trigger jet in p+Pb: " << integral << std::endl;
  //  std::cout << std::endl << std::endl;


  //  std::cout << "---------------" << std::endl << "JETS IN MC > " << trigpt << " GeV" << std::endl << "---------------" << std::endl;

  //  njet[0] = 0;
  //  for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
  //    if (0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) >= trigpt && 0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) < maxpt)
  //      njet[0] += h_jet_counts_ref[1][iPtJ][0]->GetBinContent (1);
  //  } // end loop over iPtJ
  //  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
  //    njet[iCent+1] = 0;
  //    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
  //      if (0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) >= trigpt && 0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) < maxpt)
  //        njet[iCent+1] += h_jet_counts[1][iPtJ][iCent][0]->GetBinContent (1);
  //    } // end loop over iPtJ
  //  } // end loop over iCent

  //  std::cout << "Number of pp jets: " << njet[0] << std::endl;
  //  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++)
  //    std::cout << "Number of p+Pb " << zdcCentPercs[iCent+1] << "\%-" << zdcCentPercs[iCent] << "\% jets: " << njet[iCent+1] << std::endl;

  //  std::cout << "Formatted for latex:" << std::endl;
  //  std::cout << (int) njet[0];
  //  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++)
  //    std::cout << " & " << (int) njet[iCent+1];
  //  std::cout << std::endl;


  //  integral = 0;
  //  for (short iX = h_jet_pt_ref[1][0]->FindBin (trigpt); iX <= h_jet_pt_ref[1][0]->GetNbinsX (); iX++)
  //    integral += h_jet_pt_ref[1][0]->GetBinContent (iX) * h_jet_pt_ref[1][0]->GetBinWidth (iX);
  //  std::cout << "Average number of jets per trigger jet in pp: " << integral << std::endl;
  //  integral = 0;
  //  for (short iX = h_jet_pt[1][nZdcCentBins][0]->FindBin (trigpt); iX <= h_jet_pt[1][nZdcCentBins][0]->GetNbinsX (); iX++)
  //    integral += h_jet_pt[1][nZdcCentBins][0]->GetBinContent (iX) * h_jet_pt[1][nZdcCentBins][0]->GetBinWidth (iX);
  //  std::cout << "Average number of jets per trigger jet in p+Pb: " << integral << std::endl;
  //  std::cout << std::endl << std::endl;

  //  delete[] njet;
  //}



  l->SetLineStyle (7);
  l->SetLineWidth (1);
  l->SetLineColor (kBlack);



  for (short iDType = 0; iDType < 2; iDType++) {

    const char* canvasName = Form ("c_jet_pt");

    const short iSamp = (iDType == 1 ? nSamps-1 : 0);

    TCanvas* c = new TCanvas (canvasName, "", 800, 1000);
    c->cd ();

    c->SetTopMargin (0.04);
    c->SetBottomMargin (0.12);
    c->SetLeftMargin (0.12);
    c->SetRightMargin (0.03);

    TH1D* h = nullptr;
    TGAE* g = nullptr;

    double ymin=2e-13;
    double ymax=1e3;

    c->SetLogx();
    c->SetLogy ();

    const double maxx = 400;
    h = new TH1D ("htemp", ";#it{p}_{T}^{jet} [GeV];(1 / N_{jet}) (dN_{jet} / d#it{p}_{T}^{jet}) [GeV^{-1}]", 1, pTJBins[0], maxx);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitleSize (0.036);
    h->GetXaxis ()->SetLabelSize (0.036);
    h->GetXaxis ()->SetTitleOffset (1.5);
    h->GetYaxis ()->SetTitleSize (0.036);
    h->GetYaxis ()->SetLabelSize (0.036);
    h->GetYaxis ()->SetTitleOffset (1.5);

    h->SetLineWidth (0);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    {
      h = (TH1D*) h_jet_pt_ref[2*iDType][iSamp]->Clone ("htemp");
      h->Scale (std::pow (10, 3));
      myDrawHist (h, kBlack, 1, 2);
      SaferDelete (&h);

      h = h_jet_pt_ref[2*iDType][iSamp];
   
      g = make_graph (h);
      RecenterGraph (g);
      ScaleGraph (g, nullptr, std::pow (10, 3));
      myDraw (g, colorfulColors[0], kFullCircle, 1.4, 1, 3, "P", false);
      SaferDelete (&g);
    }

    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      h = (TH1D*) h_jet_pt_ref[2*iDType][iSamp]->Clone ("htemp");
      h->Scale (std::pow (10, 2-iCent));
      myDrawHist (h, kBlack, 1, 2);
      SaferDelete (&h);

      h = h_jet_pt[2*iDType][iCent][iSamp];

      g = make_graph (h);
      ScaleGraph (g, nullptr, std::pow (10, 2-iCent));
      RecenterGraph (g);
      myDraw (g, colorfulColors[iCent+1], kFullCircle, 1.4, 1, 3, "P", false);
      SaferDelete (&g);

    } // end loop over iCent

    myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.034);
    myText (0.61, 0.89, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.034);
    myText (0.61, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.034);

    mySimpleMarkerAndBoxAndLineText (0.27, 0.255, 1.4, 1001, colorfulSystColors[0], 1.0, colorfulColors[0], kFullCircle, 1.6, "#it{pp} (#times10^{3})", 0.028);

    for (short iCent = 0; iCent < nZdcCentBins; iCent++) {
      mySimpleMarkerAndBoxAndLineText (0.27 + (iCent >= 2 ? 0.34 : 0), 0.255-((iCent+1)%3)*0.035, 1.4, 1001, colorfulSystColors[iCent+1], 1.0, colorfulColors[iCent+1], kFullCircle, 1.6, Form ("%s %i-%i%% (#times10^{%i})", iDType == 0 ? "ZDC" : "FCal", zdcCentPercs[iCent+1], zdcCentPercs[iCent], 2-iCent), 0.028);
    } // end loop over iCent
    mySimpleMarkerAndBoxAndLineText (0.61, 0.15, 1.4, 1001, colorfulSystColors[nZdcCentBins+1], 1.0, colorfulColors[nZdcCentBins+1], kFullCircle, 1.6, Form ("All cent. (#times10^{%i})", 2-nZdcCentBins), 0.028);
    mySimpleMarkerAndBoxAndLineText (0.27, 0.15, 1.4, 1001, kWhite, 0.0, kBlack, kDot, 0.0, "#it{pp} (#it{scaled})", 0.028);

    c->RedrawAxis();

    c->SaveAs (Form ("%s/Plots/JetDistributions/JetPtSpectrum%s.pdf", workPath.Data (), iDType == 1 ? "_mc" : ""));

  }



  for (short iDType = 0; iDType < 2; iDType++) {
    const char* canvasName = Form ("c_jet_trk_pt_jet_ratio_%s", iDType == 1 ? "mc" : "data");
    TCanvas* c = new TCanvas (canvasName, "", 1200, 800);
    c->Divide (3, 2);

    for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {
      c->cd (nZdcCentBins+1-iCent);

      gPad->SetLogx ();

      TH1D* h = new TH1D ("h", ";#it{p}_{T}^{jet} [GeV];#it{p}+Pb / #it{pp}", 1, pTJBins[0], pTJBins[nPtJBins]);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (0.5, 2.0);
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      myDraw (h_jet_pt_ratio[iDType][iCent][0], colorfulColors[nZdcCentBins-iCent], kFullCircle, 1.0, 1, 2, false);

      if (iCent < nZdcCentBins)
        myText (0.2, 0.865, colorfulColors[nZdcCentBins-iCent], Form ("#bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
      else
        myText (0.2, 0.865, colorfulColors[0], "#bf{All centralities}", 0.05);
      if (iCent == nZdcCentBins) {
        myText (0.2, 0.80, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.05);
        myText (0.2, 0.74, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.05);
      }

    } // end loop over iCent

    c->cd ();
    myText (0.065, 0.971, kBlack, "#bf{#it{ATLAS}} Internal", 0.027);

    c->SaveAs (Form ("%s/Plots/JetDistributions/JetPtSpectrum_RatioSummary%s.pdf", workPath.Data (), iDType == 1 ? "_mc" : ""));
  } // end loop over iDType




  {
    const char* canvasName = "c_jet_trk_pt_jet_datamc_ratio";
    TCanvas* c = new TCanvas (canvasName, "", 1300, 700);
    c->Divide (4, 2);

    {
      c->cd (7);

      gPad->SetLogx ();

      TH1D* h = new TH1D ("h", ";#it{p}_{T}^{jet} [GeV];Data / MC", 1, pTJBins[0], pTJBins[nPtJBins]);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (0.00, 1.8);
      //h->GetYaxis ()->SetRangeUser (0.10, 0.4);
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      myDraw (h_jet_pt_datamc_ratio_ref[0][0], colorfulColors[0], kOpenCircle, 1.0, 1, 2, false);
      myDraw (h_jet_pt_datamc_ratio_ref[1][0], colorfulColors[0], kFullCircle, 1.0, 1, 2, false);
      myDraw (f_jet_pt_datamc_ratio_ref[1][0], colorfulColors[0], 1, 2);

      myText (0.24, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
    }

    for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {
      c->cd (nZdcCentBins+1-iCent);

      gPad->SetLogx ();

      TH1D* h = new TH1D ("h", ";#it{p}_{T}^{jet} [GeV];Data / MC", 1, pTJBins[0], pTJBins[nPtJBins]);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (0.00, 1.8);
      //h->GetYaxis ()->SetRangeUser (0.10, 0.4);
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      myDraw (h_jet_pt_datamc_ratio[0][iCent][0], colorfulColors[iCent+1], kOpenCircle, 1.0, 1, 2, false);
      myDraw (h_jet_pt_datamc_ratio[1][iCent][0], colorfulColors[iCent+1], kFullCircle, 1.0, 1, 2, false);
      myDraw (f_jet_pt_datamc_ratio[1][iCent][0], colorfulColors[iCent+1], 1, 2);

      if (iCent < nZdcCentBins)
        myText (0.24, 0.865, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
      else
        myText (0.24, 0.865, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.05);

    } // end loop over iCent

    c->cd (8);
    myText (0.1, 0.84, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.07);
    myText (0.1, 0.75, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.66, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.07);
    myText (0.1, 0.57, kBlack, "Pythia8 JZ0-3", 0.07);

    c->SaveAs (Form ("%s/Plots/JetDistributions/JetPtSpectrum_DataMC_RatioSummary.pdf", workPath.Data ()));
  }




  /*{
    const char* canvasName = "c_jet_eta_phi_ref";

    TCanvas* c = new TCanvas (canvasName, "", 880, 800);

    c->SetRightMargin (0.15);
    c->SetLeftMargin (0.15);
    c->SetTopMargin (0.04);
    c->SetBottomMargin (0.15);

    TH2D* h2 = h2_jet_eta_phi_ref[0][0];

    h2->Draw ("colz");

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    //myText (0.22, 0.81, kBlack, TString (tag).Contains ("30GeVJets") ? "#it{p}_{T}^{jet} > 30 GeV" : "#it{p}_{T}^{jet} > 60 GeV", 0.032);

    c->SaveAs (Form ("%s/Plots/JetDistributions/JetEtaPhiSpectrum_pp.pdf", workPath.Data ()));
  }



  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
    const char* canvasName = Form ("c_jet_eta_phi_iCent%i", iCent);

    TCanvas* c = new TCanvas (canvasName, "", 880, 800);
    c->SetRightMargin (0.15);
    c->SetLeftMargin (0.15);
    c->SetTopMargin (0.04);
    c->SetBottomMargin (0.15);

    TH2D* h2 = h2_jet_eta_phi[0][iCent][0];

    h2->Draw ("colz");

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    if (iCent < nZdcCentBins)
      myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
    else
      myText (0.22, 0.85, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, All cent.", 0.032);
    //myText (0.22, 0.81, kBlack, TString (tag).Contains ("30GeVJets") ? "#it{p}_{T}^{jet} > 30 GeV" : "#it{p}_{T}^{jet} > 60 GeV", 0.032);

    c->SaveAs (Form ("%s/Plots/JetDistributions/JetEtaPhiSpectrum_pPb_iCent%i_%s.pdf", workPath.Data (), iCent));
  }*/


}


#endif
