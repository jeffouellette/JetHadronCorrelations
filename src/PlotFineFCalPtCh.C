#ifndef __JetHadronCorrelatorPlotFineFCalPtCh_C__
#define __JetHadronCorrelatorPlotFineFCalPtCh_C__

#include <vector>
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

std::vector <int> plotCents = {4, (int) std::floor (0.5*(numFineFcalCentBins+1)), numFineFcalCentBins-1};


TLine* l = new TLine ();
TLatex* tl = new TLatex ();


void PlotFineFCalPtCh (const char* tag, const char* inFileTag, const char* nomFileTag) {

  TFile* inFile = nullptr;

  TH1D*  h_jet_trk_pt_ns_ref        = nullptr;
  TH1D*  h_jet_trk_pt_perp_ref      = nullptr;
  TH1D*  h_jet_trk_pt_as_ref        = nullptr;

  TH1D*  h_jet_trk_pt_ns_ref_bkg    = nullptr;
  TH1D*  h_jet_trk_pt_perp_ref_bkg  = nullptr;
  TH1D*  h_jet_trk_pt_as_ref_bkg    = nullptr;

  TH1D*  h_jet_trk_pt_ns_ref_sig    = nullptr;
  TH1D*  h_jet_trk_pt_perp_ref_sig  = nullptr;
  TH1D*  h_jet_trk_pt_as_ref_sig    = nullptr;

  TH1D** h_jet_trk_pt_ns_nom        = Get1DArray <TH1D*> (numZdcCentBins);
  TH1D** h_jet_trk_pt_perp_nom      = Get1DArray <TH1D*> (numZdcCentBins);
  TH1D** h_jet_trk_pt_as_nom        = Get1DArray <TH1D*> (numZdcCentBins);

  TH1D** h_jet_trk_pt_ns_bkg_nom    = Get1DArray <TH1D*> (numZdcCentBins);
  TH1D** h_jet_trk_pt_perp_bkg_nom  = Get1DArray <TH1D*> (numZdcCentBins);
  TH1D** h_jet_trk_pt_as_bkg_nom    = Get1DArray <TH1D*> (numZdcCentBins);

  TH1D** h_jet_trk_pt_ns_sig_nom    = Get1DArray <TH1D*> (numZdcCentBins);
  TH1D** h_jet_trk_pt_perp_sig_nom  = Get1DArray <TH1D*> (numZdcCentBins);
  TH1D** h_jet_trk_pt_as_sig_nom    = Get1DArray <TH1D*> (numZdcCentBins);

  TH1D** h_jet_trk_pt_ns_iaa_nom    = Get1DArray <TH1D*> (numZdcCentBins);
  TH1D** h_jet_trk_pt_perp_iaa_nom  = Get1DArray <TH1D*> (numZdcCentBins);
  TH1D** h_jet_trk_pt_as_iaa_nom    = Get1DArray <TH1D*> (numZdcCentBins);

  TH1D** h_jet_trk_pt_ns          = Get1DArray <TH1D*> (numFineFcalCentBins+1);
  TH1D** h_jet_trk_pt_perp        = Get1DArray <TH1D*> (numFineFcalCentBins+1);
  TH1D** h_jet_trk_pt_as          = Get1DArray <TH1D*> (numFineFcalCentBins+1);

  TH1D** h_jet_trk_pt_ns_bkg      = Get1DArray <TH1D*> (numFineFcalCentBins+1);
  TH1D** h_jet_trk_pt_perp_bkg    = Get1DArray <TH1D*> (numFineFcalCentBins+1);
  TH1D** h_jet_trk_pt_as_bkg      = Get1DArray <TH1D*> (numFineFcalCentBins+1);

  TH1D** h_jet_trk_pt_ns_sig      = Get1DArray <TH1D*> (numFineFcalCentBins+1);
  TH1D** h_jet_trk_pt_perp_sig    = Get1DArray <TH1D*> (numFineFcalCentBins+1);
  TH1D** h_jet_trk_pt_as_sig      = Get1DArray <TH1D*> (numFineFcalCentBins+1);

  TH1D** h_jet_trk_pt_ns_iaa      = Get1DArray <TH1D*> (numFineFcalCentBins+1);
  TH1D** h_jet_trk_pt_perp_iaa    = Get1DArray <TH1D*> (numFineFcalCentBins+1);
  TH1D** h_jet_trk_pt_as_iaa      = Get1DArray <TH1D*> (numFineFcalCentBins+1);

  TH1D*  h_wgts_rebinned = nullptr;


  {
    TString inFileName = inFileTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/PlotFineFCalPtCh_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    for (int iCent = 0; iCent < numFineFcalCentBins; iCent++) {
      h_jet_trk_pt_ns[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_pPb_FineFcalCent%i", iCent));
      h_jet_trk_pt_perp[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_pPb_FineFcalCent%i", iCent));
      h_jet_trk_pt_as[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_pPb_FineFcalCent%i", iCent));
      h_jet_trk_pt_ns_bkg[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_pPb_bkg_FineFcalCent%i", iCent));
      h_jet_trk_pt_perp_bkg[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_pPb_bkg_FineFcalCent%i", iCent));
      h_jet_trk_pt_as_bkg[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_pPb_bkg_FineFcalCent%i", iCent));

      h_jet_trk_pt_ns_sig[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_pPb_sig_FineFcalCent%i", iCent));
      h_jet_trk_pt_perp_sig[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_pPb_sig_FineFcalCent%i", iCent));
      h_jet_trk_pt_as_sig[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_pPb_sig_FineFcalCent%i", iCent));

      h_jet_trk_pt_ns_iaa[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_iaa_FineFcalCent%i", iCent));
      h_jet_trk_pt_perp_iaa[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_iaa_FineFcalCent%i", iCent));
      h_jet_trk_pt_as_iaa[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_iaa_FineFcalCent%i", iCent));
    }

    h_jet_trk_pt_ns_iaa[numFineFcalCentBins] = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_iaa_FineFcalComb");
    h_jet_trk_pt_perp_iaa[numFineFcalCentBins] = (TH1D*) inFile->Get ("h_jet_trk_pt_perp_iaa_FineFcalComb");
    h_jet_trk_pt_as_iaa[numFineFcalCentBins] = (TH1D*) inFile->Get ("h_jet_trk_pt_as_iaa_FineFcalComb");

    h_wgts_rebinned = (TH1D*) inFile->Get ("h_wgts_rebinned");
  }
  


  {
    TString inFileName = nomFileTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/PlotPtCh_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    h_jet_trk_pt_ns_ref = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_ref_Nominal");
    h_jet_trk_pt_perp_ref = (TH1D*) inFile->Get ("h_jet_trk_pt_perp_ref_Nominal");
    h_jet_trk_pt_as_ref = (TH1D*) inFile->Get ("h_jet_trk_pt_as_ref_Nominal");

    h_jet_trk_pt_ns_ref_bkg = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_ref_bkg_Nominal");
    h_jet_trk_pt_perp_ref_bkg = (TH1D*) inFile->Get ("h_jet_trk_pt_perp_ref_bkg_Nominal");
    h_jet_trk_pt_as_ref_bkg = (TH1D*) inFile->Get ("h_jet_trk_pt_as_ref_bkg_Nominal");

    h_jet_trk_pt_ns_ref_sig = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_ref_sig_Nominal");
    h_jet_trk_pt_perp_ref_sig = (TH1D*) inFile->Get ("h_jet_trk_pt_perp_ref_sig_Nominal");
    h_jet_trk_pt_as_ref_sig = (TH1D*) inFile->Get ("h_jet_trk_pt_as_ref_sig_Nominal");

    for (int iCent = 0; iCent < numZdcCentBins; iCent++) {
      h_jet_trk_pt_ns_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_pPb_iCent%i_Nominal", iCent));
      h_jet_trk_pt_perp_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_pPb_iCent%i_Nominal", iCent));
      h_jet_trk_pt_as_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_pPb_iCent%i_Nominal", iCent));

      h_jet_trk_pt_ns_bkg_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_pPb_bkg_iCent%i_Nominal", iCent));
      h_jet_trk_pt_perp_bkg_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_pPb_bkg_iCent%i_Nominal", iCent));
      h_jet_trk_pt_as_bkg_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_pPb_bkg_iCent%i_Nominal", iCent));

      h_jet_trk_pt_ns_sig_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_pPb_sig_iCent%i_Nominal", iCent));
      h_jet_trk_pt_perp_sig_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_pPb_sig_iCent%i_Nominal", iCent));
      h_jet_trk_pt_as_sig_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_pPb_sig_iCent%i_Nominal", iCent));

      h_jet_trk_pt_ns_iaa_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_iaa_iCent%i_Nominal", iCent));
      h_jet_trk_pt_perp_iaa_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_iaa_iCent%i_Nominal", iCent));
      h_jet_trk_pt_as_iaa_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_iaa_iCent%i_Nominal", iCent));
    }
  }



  int iZdcCent = 0;
  while (iZdcCent < numZdcCentBins && zdcCentPercs[iZdcCent] > 20) iZdcCent++;
  int iZdcPeriph = 0;
  while (iZdcPeriph < numZdcCentBins && zdcCentPercs[iZdcPeriph] > 80) iZdcPeriph++;



  {
    const char* canvasName = "c_jet_trk_pt_ns_iaa";
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    c->cd (); 
    c->SetLogx ();
    //c->SetLogy ();

    float ymin = 0.2;
    float ymax = 1.5;


    h = (TH1D*) h_jet_trk_pt_ns_iaa[1]->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    //h->GetXaxis ()->SetTitleSize (0.028);
    //h->GetXaxis ()->SetLabelSize (0.028);
    h->GetYaxis ()->SetTitle ("I_{pPb} (#it{p}_{T}^{ch})");
    //h->GetYaxis ()->SetTitleSize (0.028);
    //h->GetYaxis ()->SetLabelSize (0.028);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    for (int iCent = 4; iCent < numFineFcalCentBins; iCent+=5) {
      h = (TH1D*) h_jet_trk_pt_ns_iaa[iCent]->Clone ("htemp");
      g = make_graph (h);

      deltaize (g, 0.895+((iCent-4)/5)*0.01, true);
      
      ResetXErrors (g);
      g->SetLineColor (manyColors[(iCent-4)/5]);
      g->SetLineWidth (2);
      g->SetMarkerColor (manyColors[(iCent-4)/5]);
      g->SetMarkerStyle (kFullCircle);
      g->SetMarkerSize (0.8);
      ((TGAE*) g->Clone ())->Draw ("p");
      g->SetMarkerStyle (kOpenCircle);
      g->SetMarkerColor (kBlack);
      g->SetLineWidth (0);
      ((TGAE*) g->Clone ())->Draw ("p");
      SaferDelete (&g);
    }

    h = (TH1D*) h_jet_trk_pt_ns_iaa[numFineFcalCentBins]->Clone ("htemp");
    g = make_graph (h);

    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = (TH1D*) h_jet_trk_pt_ns_iaa_nom[iZdcCent]->Clone ("htemp");
    g = make_graph (h);

    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = (TH1D*) h_jet_trk_pt_ns_iaa_nom[iZdcPeriph]->Clone ("htemp");
    g = make_graph (h);

    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kOpenTriangleUp);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    myText (0.22, 0.81, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
    for (int iCent = 4; iCent < numFineFcalCentBins; iCent+=5)
      myText (iCent < 0.5*numFineFcalCentBins ? 0.58 : 0.74, 0.53 - 0.036*0.2*(iCent < 0.5*numFineFcalCentBins ? iCent : iCent - 0.5*numFineFcalCentBins), manyColors[(iCent-4)/5], Form ("%i-%i%%", fineFcalCentPercs[iCent+1], fineFcalCentPercs[iCent]), 0.03);
    myMarkerText (0.24, 0.54-8*0.036, kBlack, kFullCircle, "FCal 0-100% (wgt'd avg.)", 1.2, 0.03, true);
    myMarkerText (0.24, 0.54-9*0.036, kBlack, kOpenSquare, Form ("Zdc %i-%i%%", zdcCentPercs[iZdcCent+1], zdcCentPercs[iZdcCent]), 1.2, 0.03, true);
    myMarkerText (0.24, 0.54-10*0.036, kBlack, kOpenTriangleUp, Form ("Zdc %i-%i%%", zdcCentPercs[iZdcPeriph+1], zdcCentPercs[iZdcPeriph]), 1.2, 0.03, true);
    

    c->SaveAs (Form ("%s/Plots/CentralityBiasStudy/JetTagged_IpA_nearside_ptch_%s_FineFcal.pdf", workPath.Data (), tag));
  }



  {
    const char* canvasName = "c_jet_trk_pt_perp_iaa";
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    c->cd (); 
    c->SetLogx ();
    //c->SetLogy ();

    float ymin = -0.5;
    float ymax = 2.0;


    h = new TH1D ("h", "", 1, 0.5, 2);
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    //h->GetXaxis ()->SetTitleSize (0.028);
    //h->GetXaxis ()->SetLabelSize (0.028);
    h->GetYaxis ()->SetTitle ("I_{pPb} (#it{p}_{T}^{ch})");
    //h->GetYaxis ()->SetTitleSize (0.028);
    //h->GetYaxis ()->SetLabelSize (0.028);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    for (int iCent = 4; iCent < numFineFcalCentBins; iCent+=5) {
      h = (TH1D*) h_jet_trk_pt_perp_iaa[iCent]->Clone ("htemp");
      g = make_graph (h);

      deltaize (g, 0.895+((iCent-4)/5)*0.01, true);
      
      ResetXErrors (g);
      g->SetLineColor (manyColors[(iCent-4)/5]);
      g->SetLineWidth (2);
      g->SetMarkerColor (manyColors[(iCent-4)/5]);
      g->SetMarkerStyle (kFullCircle);
      g->SetMarkerSize (0.8);
      ((TGAE*) g->Clone ())->Draw ("p");
      g->SetMarkerStyle (kOpenCircle);
      g->SetMarkerColor (kBlack);
      g->SetLineWidth (0);
      ((TGAE*) g->Clone ())->Draw ("p");
      SaferDelete (&g);
    }

    h = (TH1D*) h_jet_trk_pt_perp_iaa[numFineFcalCentBins]->Clone ("htemp");
    g = make_graph (h);

    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = (TH1D*) h_jet_trk_pt_perp_iaa_nom[iZdcCent]->Clone ("htemp");
    g = make_graph (h);

    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = (TH1D*) h_jet_trk_pt_perp_iaa_nom[iZdcPeriph]->Clone ("htemp");
    g = make_graph (h);

    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kOpenTriangleUp);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    myText (0.22, 0.81, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #pi/3 < |#Delta#phi_{ch,jet}| < 2#pi/3", GetJetPtStr (tag).Data ()), 0.032);
    for (int iCent = 4; iCent < numFineFcalCentBins; iCent+=5)
      myText (iCent < 0.5*numFineFcalCentBins ? 0.58 : 0.74, 0.53 - 0.036*0.2*(iCent < 0.5*numFineFcalCentBins ? iCent : iCent - 0.5*numFineFcalCentBins), manyColors[(iCent-4)/5], Form ("%i-%i%%", fineFcalCentPercs[iCent+1], fineFcalCentPercs[iCent]), 0.03);
    myMarkerText (0.24, 0.54-8*0.036, kBlack, kFullCircle, "FCal 0-100% (wgt'd avg.)", 1.2, 0.03, true);
    myMarkerText (0.24, 0.54-9*0.036, kBlack, kOpenSquare, Form ("Zdc %i-%i%%", zdcCentPercs[iZdcCent+1], zdcCentPercs[iZdcCent]), 1.2, 0.03, true);
    myMarkerText (0.24, 0.54-10*0.036, kBlack, kOpenTriangleUp, Form ("Zdc %i-%i%%", zdcCentPercs[iZdcPeriph+1], zdcCentPercs[iZdcPeriph]), 1.2, 0.03, true);
    

    c->SaveAs (Form ("%s/Plots/CentralityBiasStudy/JetTagged_IpA_perpendicular_ptch_%s_FineFcal.pdf", workPath.Data (), tag));
  }



  {
    const char* canvasName = "c_jet_trk_pt_as_iaa";
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    c->cd (); 
    c->SetLogx ();
    //c->SetLogy ();

    float ymin = 0.2;
    float ymax = 1.5;


    h = (TH1D*) h_jet_trk_pt_as_iaa[1]->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    //h->GetXaxis ()->SetTitleSize (0.028);
    //h->GetXaxis ()->SetLabelSize (0.028);
    h->GetYaxis ()->SetTitle ("I_{pPb} (#it{p}_{T}^{ch})");
    //h->GetYaxis ()->SetTitleSize (0.028);
    //h->GetYaxis ()->SetLabelSize (0.028);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    for (int iCent = 4; iCent < numFineFcalCentBins; iCent+=5) {
      h = (TH1D*) h_jet_trk_pt_as_iaa[iCent]->Clone ("htemp");
      g = make_graph (h);

      deltaize (g, 0.895+((iCent-4)/5)*0.01, true);
      
      ResetXErrors (g);
      g->SetLineColor (manyColors[(iCent-4)/5]);
      g->SetLineWidth (2);
      g->SetMarkerColor (manyColors[(iCent-4)/5]);
      g->SetMarkerStyle (kFullCircle);
      g->SetMarkerSize (0.8);
      ((TGAE*) g->Clone ())->Draw ("p");
      g->SetMarkerStyle (kOpenCircle);
      g->SetMarkerColor (kBlack);
      g->SetLineWidth (0);
      ((TGAE*) g->Clone ())->Draw ("p");
      SaferDelete (&g);
    }

    h = (TH1D*) h_jet_trk_pt_as_iaa[numFineFcalCentBins]->Clone ("htemp");
    g = make_graph (h);

    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);


    h = (TH1D*) h_jet_trk_pt_as_iaa_nom[iZdcCent]->Clone ("htemp");
    g = make_graph (h);

    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);


    h = (TH1D*) h_jet_trk_pt_as_iaa_nom[iZdcPeriph]->Clone ("htemp");
    g = make_graph (h);

    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kOpenTriangleUp);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);


    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    myText (0.22, 0.81, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);
    for (int iCent = 4; iCent < numFineFcalCentBins; iCent+=5)
      myText (iCent < 0.5*numFineFcalCentBins ? 0.58 : 0.74, 0.53 - 0.036*0.2*(iCent < 0.5*numFineFcalCentBins ? iCent : iCent - 0.5*numFineFcalCentBins), manyColors[(iCent-4)/5], Form ("%i-%i%%", fineFcalCentPercs[iCent+1], fineFcalCentPercs[iCent]), 0.03);
    myMarkerText (0.24, 0.54-8*0.036, kBlack, kFullCircle, "FCal 0-100% (wgt'd avg.)", 1.2, 0.03, true);
    myMarkerText (0.24, 0.54-9*0.036, kBlack, kOpenSquare, Form ("Zdc %i-%i%%", zdcCentPercs[iZdcCent+1], zdcCentPercs[iZdcCent]), 1.2, 0.03, true);
    myMarkerText (0.24, 0.54-10*0.036, kBlack, kOpenTriangleUp, Form ("Zdc %i-%i%%", zdcCentPercs[iZdcPeriph+1], zdcCentPercs[iZdcPeriph]), 1.2, 0.03, true);
    

    c->SaveAs (Form ("%s/Plots/CentralityBiasStudy/JetTagged_IpA_awayside_ptch_%s_FineFcal.pdf", workPath.Data (), tag));
  }



  {
    const char* canvasName = "c_jet_trk_pt_ns_iaa_summary";
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    c->cd (); 
    c->SetLogx ();
    //c->SetLogy ();

    float ymin = 0.2;
    float ymax = 1.5;


    h = (TH1D*) h_jet_trk_pt_ns_iaa[1]->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    //h->GetXaxis ()->SetTitleSize (0.028);
    //h->GetXaxis ()->SetLabelSize (0.028);
    h->GetYaxis ()->SetTitle ("I_{pPb} (#it{p}_{T}^{ch})");
    //h->GetYaxis ()->SetTitleSize (0.028);
    //h->GetYaxis ()->SetLabelSize (0.028);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    for (int iCent : plotCents) {
      h = (TH1D*) h_jet_trk_pt_ns_iaa[iCent]->Clone ("htemp");
      g = make_graph (h);

      deltaize (g, 0.895+(iCent-4)*0.002, true);
      
      ResetXErrors (g);
      g->SetLineColor (manyColors[(iCent-4)/5]);
      g->SetLineWidth (2);
      g->SetMarkerColor (manyColors[(iCent-4)/5]);
      g->SetMarkerStyle (kFullCircle);
      g->SetMarkerSize (0.8);
      ((TGAE*) g->Clone ())->Draw ("p");
      g->SetMarkerStyle (kOpenCircle);
      g->SetMarkerColor (kBlack);
      g->SetLineWidth (0);
      ((TGAE*) g->Clone ())->Draw ("p");
      SaferDelete (&g);
    }

    h = (TH1D*) h_jet_trk_pt_ns_iaa[numFineFcalCentBins]->Clone ("htemp");
    g = make_graph (h);

    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);


    h = (TH1D*) h_jet_trk_pt_ns_iaa_nom[iZdcCent]->Clone ("htemp");
    g = make_graph (h);

    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);


    h = (TH1D*) h_jet_trk_pt_ns_iaa_nom[iZdcPeriph]->Clone ("htemp");
    g = make_graph (h);

    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kOpenTriangleUp);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    myText (0.22, 0.81, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
    int i = 0;
    for (int iCent : plotCents)
      myText (0.58, 0.53-(i++)*0.036, manyColors[(iCent-4)/5], Form ("%i-%i%%", fineFcalCentPercs[iCent+1], fineFcalCentPercs[iCent]), 0.03);
    myMarkerText (0.58, 0.53-(i++)*0.036, kBlack, kFullCircle, "FCal 0-100% (wgt'd avg.)", 1.2, 0.03, true);
    myMarkerText (0.58, 0.53-(i++)*0.036, kBlack, kOpenSquare, Form ("Zdc %i-%i%%", zdcCentPercs[iZdcCent+1], zdcCentPercs[iZdcCent]), 1.2, 0.03, true);
    myMarkerText (0.58, 0.53-(i++)*0.036, kBlack, kOpenTriangleUp, Form ("Zdc %i-%i%%", zdcCentPercs[iZdcPeriph+1], zdcCentPercs[iZdcPeriph]), 1.2, 0.03, true);
    

    c->SaveAs (Form ("%s/Plots/CentralityBiasStudy/JetTagged_IpA_nearside_ptch_%s_FineFcal_Summary.pdf", workPath.Data (), tag));
  }



  {
    const char* canvasName = "c_jet_trk_pt_as_iaa_summary";
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    c->cd (); 
    c->SetLogx ();
    //c->SetLogy ();

    float ymin = 0.2;
    float ymax = 1.5;


    h = (TH1D*) h_jet_trk_pt_as_iaa[1]->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    //h->GetXaxis ()->SetTitleSize (0.028);
    //h->GetXaxis ()->SetLabelSize (0.028);
    h->GetYaxis ()->SetTitle ("I_{pPb} (#it{p}_{T}^{ch})");
    //h->GetYaxis ()->SetTitleSize (0.028);
    //h->GetYaxis ()->SetLabelSize (0.028);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    for (int iCent : plotCents) {
      h = (TH1D*) h_jet_trk_pt_as_iaa[iCent]->Clone ("htemp");
      g = make_graph (h);

      deltaize (g, 0.895+(iCent-4)*0.002, true);
      
      ResetXErrors (g);
      g->SetLineColor (manyColors[(iCent-4)/5]);
      g->SetLineWidth (2);
      g->SetMarkerColor (manyColors[(iCent-4)/5]);
      g->SetMarkerStyle (kFullCircle);
      g->SetMarkerSize (0.8);
      ((TGAE*) g->Clone ())->Draw ("p");
      g->SetMarkerStyle (kOpenCircle);
      g->SetMarkerColor (kBlack);
      g->SetLineWidth (0);
      ((TGAE*) g->Clone ())->Draw ("p");
      SaferDelete (&g);
    }

    h = (TH1D*) h_jet_trk_pt_as_iaa[numFineFcalCentBins]->Clone ("htemp");
    g = make_graph (h);

    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);


    h = (TH1D*) h_jet_trk_pt_as_iaa_nom[iZdcCent]->Clone ("htemp");
    g = make_graph (h);

    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);


    h = (TH1D*) h_jet_trk_pt_as_iaa_nom[iZdcPeriph]->Clone ("htemp");
    g = make_graph (h);

    ResetXErrors (g);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (kOpenTriangleUp);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);


    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    myText (0.22, 0.81, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);
    int i = 0;
    for (int iCent : plotCents)
      myText (0.58, 0.53-(i++)*0.036, manyColors[(iCent-4)/5], Form ("%i-%i%%", fineFcalCentPercs[iCent+1], fineFcalCentPercs[iCent]), 0.03);
    myMarkerText (0.58, 0.53-(i++)*0.036, kBlack, kFullCircle, "FCal 0-100% (wgt'd avg.)", 1.2, 0.03, true);
    myMarkerText (0.58, 0.53-(i++)*0.036, kBlack, kOpenSquare, Form ("Zdc %i-%i%%", zdcCentPercs[iZdcCent+1], zdcCentPercs[iZdcCent]), 1.2, 0.03, true);
    myMarkerText (0.58, 0.53-(i++)*0.036, kBlack, kOpenTriangleUp, Form ("Zdc %i-%i%%", zdcCentPercs[iZdcPeriph+1], zdcCentPercs[iZdcPeriph]), 1.2, 0.03, true);
    

    c->SaveAs (Form ("%s/Plots/CentralityBiasStudy/JetTagged_IpA_awayside_ptch_%s_FineFcal_Summary.pdf", workPath.Data (), tag));
  }



  {
    const char* canvasName = "c_weights";
    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    TH1D* h = nullptr; 

    c->cd (); 

    h = (TH1D*) h_wgts_rebinned->Clone ("h");
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetXaxis ()->SetTitle ("FCal-based Centrality [%]");
    //h->GetXaxis ()->SetTitleSize (0.028);
    //h->GetXaxis ()->SetLabelSize (0.028);
    h->GetYaxis ()->SetTitle ("Weight");
    //h->GetYaxis ()->SetTitleSize (0.028);
    //h->GetYaxis ()->SetLabelSize (0.028);

    h->SetLineWidth (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    c->SaveAs (Form ("%s/Plots/CentralityBiasStudy/Weights.pdf", workPath.Data ()));
  }
}


#endif
