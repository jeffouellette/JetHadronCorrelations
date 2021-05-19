#ifndef __JetHadronCorrelatorPlotFineFCalPtCh_C__
#define __JetHadronCorrelatorPlotFineFCalPtCh_C__

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

#include <vector>
#include <iostream>
#include <math.h>

using namespace JetHadronCorrelations;

TColor* tcolor = new TColor ();
const Color_t myBlue = (Color_t) tcolor->GetColor (45, 64, 245);
const Color_t myPurple = (Color_t) tcolor->GetColor (130,  10, 130);
const Color_t myRed = (Color_t) tcolor->GetColor (255,  12,  73);
const Color_t myGreen = (Color_t) tcolor->GetColor ( 54, 167,  80);
const Color_t myOrange = (Color_t) tcolor->GetColor (255,  68,   0);

const Color_t colors[] = {kRed+1, kPink-3, kPink+9, kPink+6, kMagenta-3, kViolet-3, kViolet+7, kBlue+1, kAzure-2, kAzure+7, kAzure+6, kCyan+1, kTeal, kTeal-5, kGreen-3, kSpring+7, kSpring+5, kYellow-3, kYellow, kOrange-2, kOrange-3, kOrange+6, kOrange+10};

std::vector <int> plotCents = {4, (int) std::floor (0.5*(numFineFcalCentBins+3)), numFineFcalCentBins-1};

TLine* l = new TLine ();
TLatex* tl = new TLatex ();


void PlotFineFCalPtCh (const char* tag, const char* inFileTag, const char* nomFileTag) {

  TFile* inFile = nullptr;

  TH1D*  h_jet_trk_pt_ns_ref      = nullptr;
  TH1D*  h_jet_trk_pt_as_ref      = nullptr;

  TH1D*  h_jet_trk_pt_ns_ref_bkg  = nullptr;
  TH1D*  h_jet_trk_pt_as_ref_bkg  = nullptr;

  TH1D*  h_jet_trk_pt_ns_ref_sig  = nullptr;
  TH1D*  h_jet_trk_pt_as_ref_sig  = nullptr;

  TH1D** h_jet_trk_pt_ns_nom      = Get1DArray <TH1D*> (numZdcCentBins);
  TH1D** h_jet_trk_pt_as_nom      = Get1DArray <TH1D*> (numZdcCentBins);

  TH1D** h_jet_trk_pt_ns_bkg_nom  = Get1DArray <TH1D*> (numZdcCentBins);
  TH1D** h_jet_trk_pt_as_bkg_nom  = Get1DArray <TH1D*> (numZdcCentBins);

  TH1D** h_jet_trk_pt_ns_sig_nom  = Get1DArray <TH1D*> (numZdcCentBins);
  TH1D** h_jet_trk_pt_as_sig_nom  = Get1DArray <TH1D*> (numZdcCentBins);

  TH1D** h_jet_trk_pt_ns_iaa_nom  = Get1DArray <TH1D*> (numZdcCentBins);
  TH1D** h_jet_trk_pt_as_iaa_nom  = Get1DArray <TH1D*> (numZdcCentBins);

  TH1D** h_jet_trk_pt_ns          = Get1DArray <TH1D*> (numFineFcalCentBins+1);
  TH1D** h_jet_trk_pt_as          = Get1DArray <TH1D*> (numFineFcalCentBins+1);

  TH1D** h_jet_trk_pt_ns_bkg      = Get1DArray <TH1D*> (numFineFcalCentBins+1);
  TH1D** h_jet_trk_pt_as_bkg      = Get1DArray <TH1D*> (numFineFcalCentBins+1);

  TH1D** h_jet_trk_pt_ns_sig      = Get1DArray <TH1D*> (numFineFcalCentBins+1);
  TH1D** h_jet_trk_pt_as_sig      = Get1DArray <TH1D*> (numFineFcalCentBins+1);

  TH1D** h_jet_trk_pt_ns_iaa      = Get1DArray <TH1D*> (numFineFcalCentBins+1);
  TH1D** h_jet_trk_pt_as_iaa      = Get1DArray <TH1D*> (numFineFcalCentBins+1);


  {
    TString inFileName = inFileTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/PlotFineFCalPtCh_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    for (int iCent = 0; iCent < numFineFcalCentBins; iCent++) {
      h_jet_trk_pt_ns[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_pPb_FineFcalCent%i", iCent));
      h_jet_trk_pt_as[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_pPb_FineFcalCent%i", iCent));
      h_jet_trk_pt_ns_bkg[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_pPb_bkg_FineFcalCent%i", iCent));
      h_jet_trk_pt_as_bkg[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_pPb_bkg_FineFcalCent%i", iCent));

      h_jet_trk_pt_ns_sig[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_pPb_sig_FineFcalCent%i", iCent));
      h_jet_trk_pt_as_sig[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_pPb_sig_FineFcalCent%i", iCent));

      h_jet_trk_pt_ns_iaa[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_iaa_FineFcalCent%i", iCent));
      h_jet_trk_pt_as_iaa[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_iaa_FineFcalCent%i", iCent));
    }

    h_jet_trk_pt_ns_iaa[numFineFcalCentBins] = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_iaa_FineFcalComb");
    h_jet_trk_pt_as_iaa[numFineFcalCentBins] = (TH1D*) inFile->Get ("h_jet_trk_pt_as_iaa_FineFcalComb");
  }
  


  {
    TString inFileName = nomFileTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/PlotPtCh_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    h_jet_trk_pt_ns_ref = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_ref_Nominal");
    h_jet_trk_pt_as_ref = (TH1D*) inFile->Get ("h_jet_trk_pt_as_ref_Nominal");

    h_jet_trk_pt_ns_ref_bkg = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_ref_bkg_Nominal");
    h_jet_trk_pt_as_ref_bkg = (TH1D*) inFile->Get ("h_jet_trk_pt_as_ref_bkg_Nominal");

    h_jet_trk_pt_ns_ref_sig = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_ref_sig_Nominal");
    h_jet_trk_pt_as_ref_sig = (TH1D*) inFile->Get ("h_jet_trk_pt_as_ref_sig_Nominal");

    for (int iCent = 0; iCent < numZdcCentBins; iCent++) {
      h_jet_trk_pt_ns_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_pPb_iCent%i_Nominal", iCent));
      h_jet_trk_pt_as_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_pPb_iCent%i_Nominal", iCent));

      h_jet_trk_pt_ns_bkg_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_pPb_bkg_iCent%i_Nominal", iCent));
      h_jet_trk_pt_as_bkg_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_pPb_bkg_iCent%i_Nominal", iCent));

      h_jet_trk_pt_ns_sig_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_pPb_sig_iCent%i_Nominal", iCent));
      h_jet_trk_pt_as_sig_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_pPb_sig_iCent%i_Nominal", iCent));

      h_jet_trk_pt_ns_iaa_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_iaa_iCent%i_Nominal", iCent));
      h_jet_trk_pt_as_iaa_nom[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_iaa_iCent%i_Nominal", iCent));
    }
  }



  int iZdcCent = 0;
  while (iZdcCent < numZdcCentBins && zdcCentPercs[iZdcCent] > 20) iZdcCent++;



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
      g->SetLineColor (colors[(iCent-4)/5]);
      g->SetLineWidth (2);
      g->SetMarkerColor (colors[(iCent-4)/5]);
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
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    myText (0.22, 0.81, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
    for (int iCent = 4; iCent < std::floor (0.5*numFineFcalCentBins)+2; iCent+=5)
      myText (0.58, 0.53-((iCent-4)/5)*0.036, colors[(iCent-4)/5], Form ("%i-%i%%", fineFcalCentPercs[iCent+1], fineFcalCentPercs[iCent]), 0.03);
    for (int iCent = std::floor (0.5*numFineFcalCentBins)+2; iCent < numFineFcalCentBins; iCent+=5)
      myText (0.74, 0.53-((iCent-std::floor (0.5*numFineFcalCentBins)-2)/5)*0.036, colors[(iCent-4)/5], Form ("%i-%i%%", fineFcalCentPercs[iCent+1], fineFcalCentPercs[iCent]), 0.03);
    myMarkerText (0.34, 0.53-8*0.036, kBlack, kFullCircle, Form ("FCal %i-%i%% (wgt'd)", fcalCentPercs[numFcalCentBins], fcalCentPercs[numFcalCentBins-1]), 1.2, 0.03, true);
    myMarkerText (0.34, 0.53-9*0.036, kBlack, kOpenCircle, Form ("Zdc %i-%i%%", zdcCentPercs[numZdcCentBins], zdcCentPercs[numZdcCentBins-1]), 1.2, 0.03, true);
    

    c->SaveAs (Form ("%s/Plots/CentralityBiasStudy/JetTagged_IpA_nearside_ptch_%s_FineFcal.pdf", workPath.Data (), tag));
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
      g->SetLineColor (colors[(iCent-4)/5]);
      g->SetLineWidth (2);
      g->SetMarkerColor (colors[(iCent-4)/5]);
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
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    myText (0.22, 0.81, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);
    for (int iCent = 4; iCent < std::floor (0.5*numFineFcalCentBins)+2; iCent+=5)
      myText (0.58, 0.53-((iCent-4)/5)*0.036, colors[(iCent-4)/5], Form ("%i-%i%%", fineFcalCentPercs[iCent+1], fineFcalCentPercs[iCent]), 0.03);
    for (int iCent = std::floor (0.5*numFineFcalCentBins)+2; iCent < numFineFcalCentBins; iCent+=5)
      myText (0.74, 0.53-((iCent-std::floor (0.5*numFineFcalCentBins)-2)/5)*0.036, colors[(iCent-4)/5], Form ("%i-%i%%", fineFcalCentPercs[iCent+1], fineFcalCentPercs[iCent]), 0.03);
    myMarkerText (0.34, 0.53-8*0.036, kBlack, kFullCircle, Form ("FCal %i-%i%% (wgt'd)", fcalCentPercs[numFcalCentBins], fcalCentPercs[numFcalCentBins-1]), 1.2, 0.03, true);
    myMarkerText (0.34, 0.53-9*0.036, kBlack, kOpenCircle, Form ("Zdc %i-%i%%", zdcCentPercs[numZdcCentBins], zdcCentPercs[numZdcCentBins-1]), 1.2, 0.03, true);
    

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

      deltaize (g, 0.895+iCent*0.01, true);
      
      ResetXErrors (g);
      g->SetLineColor (colors[(iCent-4)/5]);
      g->SetLineWidth (2);
      g->SetMarkerColor (colors[(iCent-4)/5]);
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
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    myText (0.22, 0.81, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} < #pi/8 (near-side)", GetJetPtStr (tag).Data ()), 0.032);
    int i = 0;
    for (int iCent : plotCents)
      myText (0.58, 0.53-(i++)*0.036, colors[(iCent-4)/5], Form ("%i-%i%%", fineFcalCentPercs[iCent+1], fineFcalCentPercs[iCent]), 0.03);
    myMarkerText (0.58, 0.53-(i++)*0.036, kBlack, kFullCircle, Form ("FCal %i-%i%% (wgt'd)", fcalCentPercs[numFcalCentBins], fcalCentPercs[numFcalCentBins-1]), 1.2, 0.03, true);
    myMarkerText (0.58, 0.53-(i++)*0.036, kBlack, kOpenCircle, Form ("Zdc %i-%i%%", zdcCentPercs[numZdcCentBins], zdcCentPercs[numZdcCentBins-1]), 1.2, 0.03, true);
    

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

      deltaize (g, 0.895+iCent*0.01, true);
      
      ResetXErrors (g);
      g->SetLineColor (colors[(iCent-4)/5]);
      g->SetLineWidth (2);
      g->SetMarkerColor (colors[(iCent-4)/5]);
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
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (1.2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    myText (0.22, 0.81, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %s, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", GetJetPtStr (tag).Data ()), 0.032);
    int i = 0;
    for (int iCent : plotCents)
      myText (0.58, 0.53-(i++)*0.036, colors[(iCent-4)/5], Form ("%i-%i%%", fineFcalCentPercs[iCent+1], fineFcalCentPercs[iCent]), 0.03);
    myMarkerText (0.58, 0.53-(i++)*0.036, kBlack, kFullCircle, Form ("FCal %i-%i%% (wgt'd)", fcalCentPercs[numFcalCentBins], fcalCentPercs[numFcalCentBins-1]), 1.2, 0.03, true);
    myMarkerText (0.58, 0.53-(i++)*0.036, kBlack, kOpenCircle, Form ("Zdc %i-%i%%", zdcCentPercs[numZdcCentBins], zdcCentPercs[numZdcCentBins-1]), 1.2, 0.03, true);
    

    c->SaveAs (Form ("%s/Plots/CentralityBiasStudy/JetTagged_IpA_awayside_ptch_%s_FineFcal_Summary.pdf", workPath.Data (), tag));
  }
}


#endif
