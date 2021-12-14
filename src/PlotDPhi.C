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


void PlotDPhi (const char* inFileTag) {

  TFile* inFile = nullptr;

  TH1D****  h_jetInt_trk_dphi_ref             = Get3DArray <TH1D*> (2, 2, nPtChSelections);
  TH1D****  h_jetInt_trk_dphi_ref_bkg         = Get3DArray <TH1D*> (2, 2, nPtChSelections);
  TH1D***** h_jetInt_trk_dphi                 = Get4DArray <TH1D*> (2, 2, nPtChSelections, nZdcCentBins);
  TH1D***** h_jetInt_trk_dphi_bkg             = Get4DArray <TH1D*> (2, 2, nPtChSelections, nZdcCentBins);

  TH1D****  h_jetInt_trk_dphi_ref_sig         = Get3DArray <TH1D*> (2, 2, nPtChSelections);
  TH1D***** h_jetInt_trk_dphi_sig             = Get4DArray <TH1D*> (2, 2, nPtChSelections, nZdcCentBins);

  TH1D***** h_jetInt_trk_dphi_iaa             = Get4DArray <TH1D*> (2, 2, nPtChSelections, nZdcCentBins);


  TGAE****  g_jetInt_trk_dphi_ref_syst        = Get3DArray <TGAE*> (2, nPtChSelections, nVar);
  TGAE****  g_jetInt_trk_dphi_ref_bkg_syst    = Get3DArray <TGAE*> (2, nPtChSelections, nVar);
  TGAE***** g_jetInt_trk_dphi_syst            = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins, nVar);
  TGAE***** g_jetInt_trk_dphi_bkg_syst        = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins, nVar);

  TGAE****  g_jetInt_trk_dphi_ref_sig_syst    = Get3DArray <TGAE*> (2, nPtChSelections, nVar);
  TGAE***** g_jetInt_trk_dphi_sig_syst        = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins, nVar);

  TGAE***** g_jetInt_trk_dphi_iaa_syst        = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins, nVar);


  //TGAE****  g_jetInt_trk_dphi_ref_systTot     = Get3DArray <TGAE*> (2, nPtChSelections, 3);
  //TGAE****  g_jetInt_trk_dphi_ref_bkg_systTot = Get3DArray <TGAE*> (2, nPtChSelections, 3);
  //TGAE***** g_jetInt_trk_dphi_systTot         = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins, 3);
  //TGAE***** g_jetInt_trk_dphi_bkg_systTot     = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins, 3);

  TGAE****  g_jetInt_trk_dphi_ref_sig_systTot = Get3DArray <TGAE*> (2, nPtChSelections, 3);
  TGAE***** g_jetInt_trk_dphi_sig_systTot     = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins, 3);

  TGAE***** g_jetInt_trk_dphi_iaa_systTot     = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins, 3);


  //TH1D****  h_jetInt_trk_dphi_ref_syst        = Get3DArray <TH1D*> (2, nPtChSelections, nVar);
  //TH1D****  h_jetInt_trk_dphi_ref_bkg_syst    = Get3DArray <TH1D*> (2, nPtChSelections, nVar);
  //TH1D***** h_jetInt_trk_dphi_syst            = Get4DArray <TH1D*> (2, nPtChSelections, nZdcCentBins, nVar);
  //TH1D***** h_jetInt_trk_dphi_bkg_syst        = Get4DArray <TH1D*> (2, nPtChSelections, nZdcCentBins, nVar);

  TH1D****  h_jetInt_trk_dphi_ref_sig_syst    = Get3DArray <TH1D*> (2, nPtChSelections, nVar);
  TH1D***** h_jetInt_trk_dphi_sig_syst        = Get4DArray <TH1D*> (2, nPtChSelections, nZdcCentBins, nVar);

  TH1D***** h_jetInt_trk_dphi_iaa_syst        = Get4DArray <TH1D*> (2, nPtChSelections, nZdcCentBins, nVar);


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

        for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
  
          const TString ptch = pTChSelections[iPtCh];
  
          h_jetInt_trk_dphi_ref[iDType][iPtJInt][iPtCh]     = (TH1D*) inFile->Get (Form ("h_jetInt_trk_dphi_%s_ref_%s_%s_Nominal",      ptch.Data (), dType.Data (), pTJInt.Data ()));
          h_jetInt_trk_dphi_ref_bkg[iDType][iPtJInt][iPtCh] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_dphi_%s_ref_bkg_%s_%s_Nominal",  ptch.Data (), dType.Data (), pTJInt.Data ()));
          h_jetInt_trk_dphi_ref_sig[iDType][iPtJInt][iPtCh] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_dphi_%s_ref_sig_%s_%s_Nominal",  ptch.Data (), dType.Data (), pTJInt.Data ()));
  
        } // end loop over iPtCh
  
        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
  
          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));
  
          for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
  
            const TString ptch = pTChSelections[iPtCh];
  
            h_jetInt_trk_dphi[iDType][iPtJInt][iPtCh][iCent]      = (TH1D*) inFile->Get (Form ("h_jetInt_trk_dphi_%s_pPb_%s_%s_%s_Nominal",      ptch.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));
            h_jetInt_trk_dphi_bkg[iDType][iPtJInt][iPtCh][iCent]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_dphi_%s_pPb_bkg_%s_%s_%s_Nominal",  ptch.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));
            h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_dphi_%s_pPb_sig_%s_%s_%s_Nominal",  ptch.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));
  
          } // end loop over iPtCh
  
        } // end loop over iCent

      } // end loop over iPtJInt

    } // end loop over iDType



    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

         const TString ptch = pTChSelections[iPtCh];

        for (int iVar = 0; iVar < nVar; iVar++) {

          const TString var = variations[iVar];

          //g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][iVar]     = (TGAE*) inFile->Get (Form ("g_jetInt_trk_dphi_%s_ref_syst_%s_%s",      ptch.Data (), pTJInt.Data (), var.Data ()));
          //g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][iVar] = (TGAE*) inFile->Get (Form ("g_jetInt_trk_dphi_%s_ref_bkg_syst_%s_%s",  ptch.Data (), pTJInt.Data (), var.Data ()));
          g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar] = (TGAE*) inFile->Get (Form ("g_jetInt_trk_dphi_%s_ref_sig_syst_%s_%s",  ptch.Data (), pTJInt.Data (), var.Data ()));

          for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            //g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][iVar]      = (TGAE*) inFile->Get (Form ("g_jetInt_trk_dphi_%s_syst_%s_%s_%s",      ptch.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            //g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][iVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_dphi_%s_bkg_syst_%s_%s_%s",  ptch.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_dphi_%s_sig_syst_%s_%s_%s",  ptch.Data (), cent.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iVar

        for (int iTotVar = 0; iTotVar < 3; iTotVar++) {

          const TString totVar = totalVariations[iTotVar];

          //g_jetInt_trk_dphi_ref_systTot[iPtJInt][iPtCh][iTotVar]     = (TGAE*) inFile->Get (Form ("g_jetInt_trk_dphi_%s_ref_%s_systTot_%s",      ptch.Data (), totVar.Data (), pTJInt.Data ()));
          //g_jetInt_trk_dphi_ref_bkg_systTot[iPtJInt][iPtCh][iTotVar] = (TGAE*) inFile->Get (Form ("g_jetInt_trk_dphi_%s_ref_bkg_%s_systTot_%s",  ptch.Data (), totVar.Data (), pTJInt.Data ()));
          g_jetInt_trk_dphi_ref_sig_systTot[iPtJInt][iPtCh][iTotVar] = (TGAE*) inFile->Get (Form ("g_jetInt_trk_dphi_%s_ref_sig_%s_systTot_%s",  ptch.Data (), totVar.Data (), pTJInt.Data ()));

          for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            //g_jetInt_trk_dphi_systTot[iPtJInt][iPtCh][iCent][iTotVar]      = (TGAE*) inFile->Get (Form ("g_jetInt_trk_dphi_%s_%s_systTot_%s_%s",      ptch.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            //g_jetInt_trk_dphi_bkg_systTot[iPtJInt][iPtCh][iCent][iTotVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_dphi_%s_bkg_%s_systTot_%s_%s",  ptch.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            g_jetInt_trk_dphi_sig_systTot[iPtJInt][iPtCh][iCent][iTotVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_dphi_%s_sig_%s_systTot_%s_%s",  ptch.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));

          } // end loop over iCent

        } // end loop over iTotVar


        for (int iVar = 1; iVar < nVar; iVar++) {

          const TString var = variations[iVar];
          const TString dType = (dataVariations.count (var) > 0 ? "data" : "mc");

          //h_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][iVar]      = (TH1D*) inFile->Get (Form ("h_jetInt_trk_dphi_%s_ref_%s_%s_%s",     ptch.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
          //h_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][iVar]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_dphi_%s_ref_bkg_%s_%s_%s", ptch.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
          h_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_dphi_%s_ref_sig_%s_%s_%s", ptch.Data (), dType.Data (), pTJInt.Data (), var.Data ()));

          for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            //h_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][iVar]     = (TH1D*) inFile->Get (Form ("h_jetInt_trk_dphi_%s_pPb_%s_%s_%s_%s",     ptch.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
            //h_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_dphi_%s_pPb_bkg_%s_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
            h_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_dphi_%s_pPb_sig_%s_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iVar

      } // end loop over iPtCh

    } // end loop over iPtJInt

  }




  {
    TString inFileName = inFileTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/ProcessUnfolding_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    for (short iDType = 0; iDType < 2; iDType++) {

      const TString dType = (iDType == 0 ? "data" : "mc");

      for (short iPtJInt : {0, 1}) {

        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

          for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
  
            const TString ptch = pTChSelections[iPtCh];
  
            h_jetInt_trk_dphi_iaa[iDType][iPtJInt][iPtCh][iCent]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_dphi_%s_iaa_%s_%s_%s_Nominal", ptch.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));
  
          } // end loop over iPtCh

        } // end loop over iCent

      } // end loop over iPtJInt

    } // end loop over iDType



    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

         const TString ptch = pTChSelections[iPtCh];

        for (int iVar = 0; iVar < nVar; iVar++) {

          const TString var = variations[iVar];

          for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][iVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_dphi_%s_iaa_syst_%s_%s_%s",  ptch.Data (), cent.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iVar

        for (int iTotVar = 0; iTotVar < 3; iTotVar++) {

          const TString totVar = totalVariations[iTotVar];

          for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            g_jetInt_trk_dphi_iaa_systTot[iPtJInt][iPtCh][iCent][iTotVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_dphi_%s_iaa_%s_systTot_%s_%s",  ptch.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));

          } // end loop over iCent

        } // end loop over iTotVar


        for (int iVar = 1; iVar < nVar; iVar++) {

          const TString var = variations[iVar];
          const TString dType = (dataVariations.count (var) > 0 ? "data" : "mc");

          for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            h_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_dphi_%s_iaa_%s_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iVar

      } // end loop over iPtCh

    } // end loop over iPtJInt

  }



  l->SetLineStyle (7);
  l->SetLineWidth (1);
  l->SetLineColor (kBlack);



  for (int iDType = 0; iDType < 2; iDType++) {

    const TString dType = (iDType == 0 ? "data" : "mc");

    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

        const TString ptch = pTChSelections[iPtCh];

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          const char* canvasName = Form ("c_jet_trk_dphi_%s_iCent%i_%s_%s", dType.Data (), iCent, pTJInt.Data (), ptch.Data ());

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

          h = new TH1D ("htemp", ";#Delta#phi_{ch,jet};(1/N_{jet}) (dN_{ch} / d#Delta#phi)", 1, 0, M_PI);
          h->SetBinContent (1, 0);
          h->GetYaxis ()->SetRangeUser (ymin, ymax);
          h->GetXaxis ()->SetTitleSize (0.028/fuPad);
          h->GetXaxis ()->SetLabelSize (0.028/fuPad);
          h->GetXaxis ()->SetTitleOffset (2.1*fuPad);
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

          //g = (TGAE*) g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][0]->Clone ();
          h = h_jetInt_trk_dphi_ref[iDType][iPtJInt][iPtCh];
          //SetCentralValuesKeepRelativeErrors (g, h);
          //myDrawSyst (g, myBlue);
          //SaferDelete (&g);
          g = make_graph (h);
          ResetXErrors (g);
          myDraw (g, myBlue, kFullCircle, 1.2);
          SaferDelete (&g);

          //g = (TGAE*) g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][0]->Clone ();
          h = h_jetInt_trk_dphi_ref_bkg[iDType][iPtJInt][iPtCh];
          //SetCentralValuesKeepRelativeErrors (g, h);
          //myDrawSyst (g, myPurple);
          //SaferDelete (&g);
          g = make_graph (h);
          ResetXErrors (g);
          myDraw (g, myPurple, kOpenCircle, 1.2);
          SaferDelete (&g);

          //g = (TGAE*) g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][0]->Clone ();
          h = h_jetInt_trk_dphi[iDType][iPtJInt][iPtCh][iCent];
          //SetCentralValuesKeepRelativeErrors (g, h);
          //myDrawSyst (g, myRed);
          //SaferDelete (&g);
          g = make_graph (h);
          ResetXErrors (g);
          myDraw (g, myRed, kFullCircle, 1.2);
          SaferDelete (&g);

          //g = (TGAE*) g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][0]->Clone ();
          h = h_jetInt_trk_dphi_bkg[iDType][iPtJInt][iPtCh][iCent];
          //SetCentralValuesKeepRelativeErrors (g, h);
          //myDrawSyst (g, myGreen);
          //SaferDelete (&g);
          g = make_graph (h);
          ResetXErrors (g);
          myDraw (g, myGreen, kOpenCircle, 1.2);
          SaferDelete (&g);

          myText (0.30, 0.83, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
          myText (0.30, 0.77, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.020/fuPad);
          myText (0.30, 0.71, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.020/fuPad);
          myText (0.30, 0.65, kBlack, "Jet-hadron correlations", 0.020/fuPad);
          myText (0.30, 0.59, kBlack, Form ("#it{p}_{T}^{jet} > %s, %s", pTJInt.Data (), pTChStrs[ptch].Data ()), 0.020/fuPad);
          myText (0.30, 0.53, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fuPad);


          cPad->cd (); 

          ymin = -4;
          ymax = 28;

          h = new TH1D ("htemp", ";#Delta#phi_{ch,jet};(Sig.+Bkg.) - Bkg.", 1, 0, M_PI);
          h->SetBinContent (1, 0);
          h->GetXaxis ()->SetMoreLogLabels ();
          h->GetYaxis ()->SetRangeUser (ymin, ymax);
          h->GetXaxis ()->SetTitleSize (0.028/fdPad);
          h->GetXaxis ()->SetLabelSize (0.028/fdPad);
          h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
          //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
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

          //g = (TGAE*) g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][0]->Clone ();
          h = h_jetInt_trk_dphi_ref_sig[iDType][iPtJInt][iPtCh];
          //SetCentralValuesKeepRelativeErrors (g, h);
          //myDrawSyst (g, myBlue);
          //SaferDelete (&g);
          g = make_graph (h);
          ResetXErrors (g);
          myDraw (g, myBlue, kFullCircle, 1.2);
          SaferDelete (&g);

          //g = (TGAE*) g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][0]->Clone ();
          h = h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent];
          //SetCentralValuesKeepRelativeErrors (g, h);
          //myDrawSyst (g, myRed);
          //SaferDelete (&g);
          g = make_graph (h);
          ResetXErrors (g);
          myDraw (g, myRed, kFullCircle, 1.2);
          SaferDelete (&g);

          myBoxText2 (0.36, 0.84, myRed, kFullCircle, "#it{p}+Pb total", 0.8, 0.020/fcPad, true);
          myBoxText2 (0.36, 0.75, myBlue, kFullCircle, "#it{pp} total", 0.8, 0.020/fcPad, true);
          myBoxText2 (0.36, 0.66, myGreen, kOpenCircle, "#it{p}+Pb bkgd.", 0.8, 0.020/fcPad);
          myBoxText2 (0.36, 0.57, myPurple, kOpenCircle, "#it{pp} bkgd.", 0.8, 0.020/fcPad);


          dPad->cd (); 

          ymin = 0.5;
          ymax = 1.5;
          //ymin = strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 0.8 : 0.83;
          //ymax = strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.4 : 1.17;

          h = new TH1D ("htemp", ";#Delta#phi_{ch,jet};#it{I}_{#it{p}Pb} = #it{p}+Pb / #it{pp}", 1, 0, M_PI);
          h->SetBinContent (1, 1);
          h->GetXaxis ()->SetMoreLogLabels ();
          h->GetYaxis ()->SetRangeUser (ymin, ymax);
          h->GetXaxis ()->SetTitleSize (0.028/fdPad);
          h->GetXaxis ()->SetLabelSize (0.028/fdPad);
          h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
          //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
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

          //g = (TGAE*) g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][0]->Clone ();
          h = h_jetInt_trk_dphi_iaa[iDType][iPtJInt][iPtCh][iCent];
          //SetCentralValuesKeepRelativeErrors (g, h);
          //myDrawSyst (g, myRed);
          //SaferDelete (&g);
          g = make_graph (h);
          ResetXErrors (g);
          myDraw (g, myRed, kFullCircle, 1.2);
          SaferDelete (&g);

          c->SaveAs (Form ("%s/Plots/DPhi/JetTagged_HadronYields_%i-%iperc_comparison_dphi_%s_%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], ptch.Data (), pTJInt.Data ()));

        } // end loop over iCent

      } // end loop over iPtCh

    } // end loop over iPtJInt

  } // end loop over iDType



  //for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
  //  const char* canvasName = Form ("c_jet_trk_dphi_gt2_lt4_FCalvsZdc_iCent%i", iCent);
  //  TCanvas* c = new TCanvas (canvasName, "", 800, 1120);
  //  c->cd ();
  //  const double fuPad = 480./1120.;
  //  const double fdPad = 320./1120.;
  //  const double fcPad = 1.0 - fuPad - fdPad;
  //  TPad* uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 1.0-fuPad, 1.0, 1.0);
  //  TPad* cPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, fdPad, 1.0, 1.0-fuPad);
  //  TPad* dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, fdPad);

  //  uPad->SetBottomMargin (0);
  //  cPad->SetTopMargin (0);
  //  cPad->SetBottomMargin (0);
  //  dPad->SetTopMargin (0);
  //  dPad->SetBottomMargin (0.25);
  //  uPad->Draw ();
  //  cPad->Draw ();
  //  dPad->Draw ();

  //  const int iVar = GetVarN ("FcalCentVar");

  //  TH1D* h = nullptr; 
  //  TGAE* g = nullptr;

  //  uPad->cd (); 

  //  float ymin = -4;
  //  float ymax = 33;

  //  h = (TH1D*) h_jetInt_trk_dphi_ref[3]->Clone ("h");
  //  h->Reset ();
  //  h->GetYaxis ()->SetRangeUser (ymin, ymax);
  //  h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
  //  h->GetXaxis ()->SetTitleSize (0.028/fuPad);
  //  h->GetXaxis ()->SetLabelSize (0.028/fuPad);
  //  h->GetXaxis ()->SetTitleOffset (2.1*fuPad);
  //  h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#Delta#phi)");
  //  h->GetYaxis ()->SetTitleSize (0.028/fuPad);
  //  h->GetYaxis ()->SetLabelSize (0.028/fuPad);
  //  h->GetYaxis ()->SetTitleOffset (2.1*fuPad);

  //  h->SetLineWidth (1);
  //  h->SetLineStyle (2);
  //  h->DrawCopy ("hist ][");
  //  SaferDelete (&h);

  //  TBox* shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

  //  shadedBox = new TBox (0, ymin, M_PI/8., ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

  //  h = h_jetInt_trk_dphi_ref[3];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, kBlack, kFullCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_ref_bkg[3];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, kBlack, kOpenCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi[iCent][3];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myRed, kFullCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_bkg[iCent][3];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myGreen, kOpenCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_syst[iCent][3][iVar];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myBlue, kFullSquare, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_bkg_syst[iCent][3][iVar];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myOrange, kOpenCircle, 0.8);
  //  SaferDelete (&g);

  //  myText (0.30, 0.83, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
  //  myText (0.30, 0.77, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.020/fuPad);
  //  myText (0.30, 0.71, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.020/fuPad);
  //  myText (0.30, 0.65, kBlack, "Jet-hadron correlations", 0.020/fuPad);
  //  myText (0.30, 0.59, kBlack, Form ("#it{p}_{T}^{jet} > %s, 2 < #it{p}_{T}^{ch} < 4 GeV", GetJetPtStr (tag).Data ()), 0.020/fuPad);
  //  myText (0.30, 0.53, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fuPad);


  //  cPad->cd (); 

  //  ymin = -4;
  //  ymax = 28;

  //  h = (TH1D*) h_jetInt_trk_dphi_ref_sig[3]->Clone ("h");
  //  h->Reset ();
  //  h->GetXaxis ()->SetMoreLogLabels ();
  //  h->GetYaxis ()->SetRangeUser (ymin, ymax);
  //  h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
  //  h->GetXaxis ()->SetTitleSize (0.028/fdPad);
  //  h->GetXaxis ()->SetLabelSize (0.028/fdPad);
  //  h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
  //  //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
  //  h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
  //  h->GetYaxis ()->SetTitleSize (0.028/fdPad);
  //  h->GetYaxis ()->SetLabelSize (0.028/fdPad);
  //  h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
  //  h->GetYaxis ()->CenterTitle ();

  //  h->SetLineWidth (1);
  //  h->SetLineStyle (2);
  //  h->DrawCopy ("hist ][");
  //  SaferDelete (&h);

  //  shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

  //  shadedBox = new TBox (0, ymin, M_PI/8., ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

  //  h = h_jetInt_trk_dphi_ref_sig[3];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, kBlack, kOpenCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_sig[iCent][3];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myRed, kFullCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_sig_syst[iCent][3][iVar];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myBlue, kFullSquare, 0.8);
  //  SaferDelete (&g);

  //  myLineText2 (0.36, 0.84, myRed, kFullCircle, Form ("#it{p}+Pb #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.8, 0.020/fcPad, true);
  //  myLineText2 (0.36, 0.75, myGreen, kOpenCircle, Form ("#it{p}+Pb bkgd., #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.8, 0.020/fcPad, true);
  //  myLineText2 (0.36, 0.66, myBlue, kFullSquare, Form ("#it{p}+Pb #bf{FCal %i-%i%%}", fcalCentPercs[iCent+1], fcalCentPercs[iCent]), 0.8, 0.020/fcPad, true);
  //  myLineText2 (0.36, 0.57, myOrange, kOpenSquare, Form ("#it{p}+Pb bkgd., #bf{FCal %i-%i%%}", fcalCentPercs[iCent+1], fcalCentPercs[iCent]), 0.8, 0.020/fcPad, true);
  //  myLineText2 (0.36, 0.48, kBlack, kFullCircle, "#it{pp} total", 0.8, 0.020/fcPad);
  //  myLineText2 (0.36, 0.39, kBlack, kOpenCircle, "#it{pp} bkgd.", 0.8, 0.020/fcPad);


  //  dPad->cd (); 

  //  ymin = 0.5;
  //  ymax = 1.4;

  //  h = (TH1D*) h_jetInt_trk_dphi_iaa[iCent][3]->Clone ("h");
  //  h->Reset ();
  //  for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
  //  h->GetXaxis ()->SetMoreLogLabels ();
  //  h->GetYaxis ()->SetRangeUser (ymin, ymax);
  //  h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
  //  h->GetXaxis ()->SetTitleSize (0.028/fdPad);
  //  h->GetXaxis ()->SetLabelSize (0.028/fdPad);
  //  h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
  //  //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
  //  h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}Pb} = #it{p}+Pb / #it{pp}");
  //  h->GetYaxis ()->SetTitleSize (0.028/fdPad);
  //  h->GetYaxis ()->SetLabelSize (0.028/fdPad);
  //  h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
  //  h->GetYaxis ()->CenterTitle ();

  //  h->SetLineWidth (1);
  //  h->SetLineStyle (2);
  //  h->DrawCopy ("hist ][");
  //  SaferDelete (&h);

  //  shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

  //  shadedBox = new TBox (0, ymin, M_PI/8., ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

  //  h = h_jetInt_trk_dphi_iaa[iCent][3];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myRed, kFullCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_iaa_syst[iCent][3][iVar];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myBlue, kFullSquare, 0.8);
  //  SaferDelete (&g);

  //  myLineText2 (0.36, 0.48, myRed, kFullCircle, Form ("#bf{Zdc %i-%i%%}", fcalCentPercs[iCent+1], fcalCentPercs[iCent]), 0.8, 0.020/fdPad, true);
  //  myLineText2 (0.36, 0.38, myBlue, kFullSquare, Form ("#bf{FCal %i-%i%%}", fcalCentPercs[iCent+1], fcalCentPercs[iCent]), 0.8, 0.020/fdPad, true);

  //  c->SaveAs (Form ("%s/Plots/DPhi/JetTagged_HadronYields_%i-%iperc_FCalvsZDC_dphi_gt2_lt4_%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], tag));
  //}



  //for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
  //  const char* canvasName = Form ("c_jet_trk_dphi_gt0p5_lt1_FCalvsZdc_iCent%i", iCent);
  //  TCanvas* c = new TCanvas (canvasName, "", 800, 1120);
  //  c->cd ();
  //  const double fuPad = 480./1120.;
  //  const double fdPad = 320./1120.;
  //  const double fcPad = 1.0 - fuPad - fdPad;
  //  TPad* uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 1.0-fuPad, 1.0, 1.0);
  //  TPad* cPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, fdPad, 1.0, 1.0-fuPad);
  //  TPad* dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, fdPad);

  //  uPad->SetBottomMargin (0);
  //  cPad->SetTopMargin (0);
  //  cPad->SetBottomMargin (0);
  //  dPad->SetTopMargin (0);
  //  dPad->SetBottomMargin (0.25);
  //  uPad->Draw ();
  //  cPad->Draw ();
  //  dPad->Draw ();

  //  const int iVar = GetVarN ("FcalCentVar");

  //  TH1D* h = nullptr; 
  //  TGAE* g = nullptr;

  //  uPad->cd (); 

  //  float ymin = -4;
  //  float ymax = 33;

  //  h = (TH1D*) h_jetInt_trk_dphi_ref[0]->Clone ("h");
  //  h->Reset ();
  //  h->GetYaxis ()->SetRangeUser (ymin, ymax);
  //  h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
  //  h->GetXaxis ()->SetTitleSize (0.028/fuPad);
  //  h->GetXaxis ()->SetLabelSize (0.028/fuPad);
  //  h->GetXaxis ()->SetTitleOffset (2.1*fuPad);
  //  h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#Delta#phi)");
  //  h->GetYaxis ()->SetTitleSize (0.028/fuPad);
  //  h->GetYaxis ()->SetLabelSize (0.028/fuPad);
  //  h->GetYaxis ()->SetTitleOffset (2.1*fuPad);

  //  h->SetLineWidth (1);
  //  h->SetLineStyle (2);
  //  h->DrawCopy ("hist ][");
  //  SaferDelete (&h);

  //  TBox* shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

  //  shadedBox = new TBox (0, ymin, M_PI/8., ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

  //  h = h_jetInt_trk_dphi_ref[0];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, kBlack, kFullCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_ref_bkg[0];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, kBlack, kOpenCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi[iCent][0];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myRed, kFullCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_bkg[iCent][0];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myGreen, kOpenCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_syst[iCent][0][iVar];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myBlue, kFullSquare, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_bkg_syst[iCent][0][iVar];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myOrange, kOpenSquare, 0.8);
  //  SaferDelete (&g);

  //  myText (0.30, 0.83, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
  //  myText (0.30, 0.77, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.020/fuPad);
  //  myText (0.30, 0.71, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.020/fuPad);
  //  myText (0.30, 0.65, kBlack, "Jet-hadron correlations", 0.020/fuPad);
  //  myText (0.30, 0.59, kBlack, Form ("#it{p}_{T}^{jet} > %s, 0.5 < #it{p}_{T}^{ch} < 1 GeV", GetJetPtStr (tag).Data ()), 0.020/fuPad);
  //  myText (0.30, 0.53, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fuPad);


  //  cPad->cd (); 

  //  ymin = -4;
  //  ymax = 28;

  //  h = (TH1D*) h_jetInt_trk_dphi_ref_sig[0]->Clone ("h");
  //  h->Reset ();
  //  h->GetXaxis ()->SetMoreLogLabels ();
  //  h->GetYaxis ()->SetRangeUser (ymin, ymax);
  //  h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
  //  h->GetXaxis ()->SetTitleSize (0.028/fdPad);
  //  h->GetXaxis ()->SetLabelSize (0.028/fdPad);
  //  h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
  //  //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
  //  h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
  //  h->GetYaxis ()->SetTitleSize (0.028/fdPad);
  //  h->GetYaxis ()->SetLabelSize (0.028/fdPad);
  //  h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
  //  h->GetYaxis ()->CenterTitle ();

  //  h->SetLineWidth (1);
  //  h->SetLineStyle (2);
  //  h->DrawCopy ("hist ][");
  //  SaferDelete (&h);

  //  shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

  //  shadedBox = new TBox (0, ymin, M_PI/8., ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

  //  h = h_jetInt_trk_dphi_ref_sig[0];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, kBlack, kFullCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_sig[iCent][0];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myRed, kFullCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_sig_syst[iCent][0][iVar];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myBlue, kFullSquare, 0.8);
  //  SaferDelete (&g);

  //  myLineText2 (0.36, 0.84, myRed, kFullCircle, Form ("#it{p}+Pb #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.8, 0.020/fcPad, true);
  //  myLineText2 (0.36, 0.75, myGreen, kOpenCircle, Form ("#it{p}+Pb bkgd., #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.8, 0.020/fcPad, true);
  //  myLineText2 (0.36, 0.66, myBlue, kFullSquare, Form ("#it{p}+Pb #bf{FCal %i-%i%%}", fcalCentPercs[iCent+1], fcalCentPercs[iCent]), 0.8, 0.020/fcPad, true);
  //  myLineText2 (0.36, 0.57, myOrange, kOpenSquare, Form ("#it{p}+Pb bkgd., #bf{FCal %i-%i%%}", fcalCentPercs[iCent+1], fcalCentPercs[iCent]), 0.8, 0.020/fcPad, true);
  //  myLineText2 (0.36, 0.48, kBlack, kFullCircle, "#it{pp} total", 0.8, 0.020/fcPad);
  //  myLineText2 (0.36, 0.39, kBlack, kOpenCircle, "#it{pp} bkgd.", 0.8, 0.020/fcPad);


  //  dPad->cd (); 

  //  ymin = 0.5;
  //  ymax = 1.4;

  //  h = (TH1D*) h_jetInt_trk_dphi_iaa[iCent][0]->Clone ("h");
  //  h->Reset ();
  //  for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
  //  h->GetXaxis ()->SetMoreLogLabels ();
  //  h->GetYaxis ()->SetRangeUser (ymin, ymax);
  //  h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
  //  h->GetXaxis ()->SetTitleSize (0.028/fdPad);
  //  h->GetXaxis ()->SetLabelSize (0.028/fdPad);
  //  h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
  //  //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
  //  h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}Pb} = #it{p}+Pb / #it{pp}");
  //  h->GetYaxis ()->SetTitleSize (0.028/fdPad);
  //  h->GetYaxis ()->SetLabelSize (0.028/fdPad);
  //  h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
  //  h->GetYaxis ()->CenterTitle ();

  //  h->SetLineWidth (1);
  //  h->SetLineStyle (2);
  //  h->DrawCopy ("hist ][");
  //  SaferDelete (&h);

  //  shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

  //  shadedBox = new TBox (0, ymin, M_PI/8., ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

  //  h = h_jetInt_trk_dphi_iaa[iCent][0];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myRed, kFullCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_iaa_syst[iCent][0][iVar];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myBlue, kFullSquare, 0.8);
  //  SaferDelete (&g);

  //  myLineText2 (0.36, 0.48, myRed, kFullCircle, Form ("#bf{Zdc %i-%i%%}", fcalCentPercs[iCent+1], fcalCentPercs[iCent]), 0.8, 0.020/fdPad, true);
  //  myLineText2 (0.36, 0.38, myBlue, kFullSquare, Form ("#bf{FCal %i-%i%%}", fcalCentPercs[iCent+1], fcalCentPercs[iCent]), 0.8, 0.020/fdPad, true);

  //  c->SaveAs (Form ("%s/Plots/DPhi/JetTagged_HadronYields_%i-%iperc_FCalvsZDC_dphi_gt0p5_lt1_%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], tag));
  //}



  //for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
  //  const char* canvasName = Form ("c_jet_trk_dphi_gt0p5_lt1_NoFcalMixCatVar_iCent%i", iCent);
  //  TCanvas* c = new TCanvas (canvasName, "", 800, 1120);
  //  c->cd ();
  //  const double fuPad = 480./1120.;
  //  const double fdPad = 320./1120.;
  //  const double fcPad = 1.0 - fuPad - fdPad;
  //  TPad* uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 1.0-fuPad, 1.0, 1.0);
  //  TPad* cPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, fdPad, 1.0, 1.0-fuPad);
  //  TPad* dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, fdPad);

  //  uPad->SetBottomMargin (0);
  //  cPad->SetTopMargin (0);
  //  cPad->SetBottomMargin (0);
  //  dPad->SetTopMargin (0);
  //  dPad->SetBottomMargin (0.25);
  //  uPad->Draw ();
  //  cPad->Draw ();
  //  dPad->Draw ();

  //  const int iVar = GetVarN ("NoFcalMixCatVar");

  //  TH1D* h = nullptr; 
  //  TGAE* g = nullptr;

  //  uPad->cd (); 

  //  float ymin = -4;
  //  float ymax = 33;

  //  h = (TH1D*) h_jetInt_trk_dphi_ref[0]->Clone ("h");
  //  h->Reset ();
  //  h->GetYaxis ()->SetRangeUser (ymin, ymax);
  //  h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
  //  h->GetXaxis ()->SetTitleSize (0.028/fuPad);
  //  h->GetXaxis ()->SetLabelSize (0.028/fuPad);
  //  h->GetXaxis ()->SetTitleOffset (2.1*fuPad);
  //  h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#Delta#phi)");
  //  h->GetYaxis ()->SetTitleSize (0.028/fuPad);
  //  h->GetYaxis ()->SetLabelSize (0.028/fuPad);
  //  h->GetYaxis ()->SetTitleOffset (2.1*fuPad);

  //  h->SetLineWidth (1);
  //  h->SetLineStyle (2);
  //  h->DrawCopy ("hist ][");
  //  SaferDelete (&h);

  //  TBox* shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

  //  shadedBox = new TBox (0, ymin, M_PI/8., ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

  //  h = h_jetInt_trk_dphi_ref[0];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, kBlack, kOpenCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_ref_bkg[0];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myGreen, kOpenCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_ref_bkg_syst[0][iVar];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myOrange, kOpenSquare, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi[iCent][0];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, kBlack, kFullCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_bkg[iCent][0];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myRed, kFullCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_bkg_syst[iCent][0][iVar];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myBlue, kFullSquare, 0.8);
  //  SaferDelete (&g);

  //  myText (0.30, 0.83, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
  //  myText (0.30, 0.77, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.020/fuPad);
  //  myText (0.30, 0.71, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.020/fuPad);
  //  myText (0.30, 0.65, kBlack, "Jet-hadron correlations", 0.020/fuPad);
  //  myText (0.30, 0.59, kBlack, Form ("#it{p}_{T}^{jet} > %s, 0.5 < #it{p}_{T}^{ch} < 1 GeV", GetJetPtStr (tag).Data ()), 0.020/fuPad);
  //  myText (0.30, 0.53, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fuPad);


  //  cPad->cd (); 

  //  ymin = -4;
  //  ymax = 28;

  //  h = (TH1D*) h_jetInt_trk_dphi_ref_sig[0]->Clone ("h");
  //  h->Reset ();
  //  h->GetXaxis ()->SetMoreLogLabels ();
  //  h->GetYaxis ()->SetRangeUser (ymin, ymax);
  //  h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
  //  h->GetXaxis ()->SetTitleSize (0.028/fdPad);
  //  h->GetXaxis ()->SetLabelSize (0.028/fdPad);
  //  h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
  //  //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
  //  h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
  //  h->GetYaxis ()->SetTitleSize (0.028/fdPad);
  //  h->GetYaxis ()->SetLabelSize (0.028/fdPad);
  //  h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
  //  h->GetYaxis ()->CenterTitle ();

  //  h->SetLineWidth (1);
  //  h->SetLineStyle (2);
  //  h->DrawCopy ("hist ][");
  //  SaferDelete (&h);

  //  shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

  //  shadedBox = new TBox (0, ymin, M_PI/8., ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

  //  h = h_jetInt_trk_dphi_ref_sig[0];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myGreen, kOpenCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_ref_sig_syst[0][iVar];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myOrange, kOpenSquare, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_sig[iCent][0];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myRed, kFullCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_sig_syst[iCent][0][iVar];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myBlue, kFullSquare, 0.8);
  //  SaferDelete (&g);

  //  myLineText2 (0.36, 0.84, kBlack, kFullCircle, "#it{p}+Pb total", 0.8, 0.020/fcPad, true);
  //  myLineText2 (0.36, 0.75, myRed, kFullCircle, "#it{p}+Pb #bf{w/} FCal matching", 0.8, 0.020/fcPad, true);
  //  myLineText2 (0.36, 0.66, myBlue, kFullSquare, "#it{p}+Pb #bf{w/out} FCal matching", 0.8, 0.020/fcPad, true);
  //  myLineText2 (0.36, 0.57, kBlack, kOpenCircle, "#it{pp} total", 0.8, 0.020/fcPad);
  //  myLineText2 (0.36, 0.48, myGreen, kOpenCircle, "#it{pp} #bf{w/} FCal matching", 0.8, 0.020/fcPad);
  //  myLineText2 (0.36, 0.39, myOrange, kOpenSquare, "#it{pp} #bf{w/out} FCal matching", 0.8, 0.020/fcPad);


  //  dPad->cd (); 

  //  ymin = 0.1;
  //  ymax = 1.2;

  //  h = (TH1D*) h_jetInt_trk_dphi_iaa[iCent][0]->Clone ("h");
  //  h->Reset ();
  //  for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
  //  h->GetXaxis ()->SetMoreLogLabels ();
  //  h->GetYaxis ()->SetRangeUser (ymin, ymax);
  //  h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
  //  h->GetXaxis ()->SetTitleSize (0.028/fdPad);
  //  h->GetXaxis ()->SetLabelSize (0.028/fdPad);
  //  h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
  //  //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
  //  h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}Pb} = #it{p}+Pb / #it{pp}");
  //  h->GetYaxis ()->SetTitleSize (0.028/fdPad);
  //  h->GetYaxis ()->SetLabelSize (0.028/fdPad);
  //  h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
  //  h->GetYaxis ()->CenterTitle ();

  //  h->SetLineWidth (1);
  //  h->SetLineStyle (2);
  //  h->DrawCopy ("hist ][");
  //  SaferDelete (&h);

  //  shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

  //  shadedBox = new TBox (0, ymin, M_PI/8., ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

  //  h = h_jetInt_trk_dphi_iaa[iCent][0];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myRed, kFullCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_iaa_syst[iCent][0][iVar];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myBlue, kFullSquare, 0.8);
  //  SaferDelete (&g);

  //  myLineText2 (0.36, 0.48, myRed, kFullCircle, "#bf{w/} FCal matching", 0.8, 0.020/fdPad, true);
  //  myLineText2 (0.36, 0.38, myBlue, kFullSquare, "#bf{w/out} FCal matching", 0.8, 0.020/fdPad, true);

  //  c->SaveAs (Form ("%s/Plots/DPhi/JetTagged_HadronYields_%i-%iperc_FCalMixCatVar_dphi_gt0p5_lt1_%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], tag));
  //}



  //for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
  //  const char* canvasName = Form ("c_jet_trk_dphi_gt2_lt4_NoFcalMixCatVar_iCent%i", iCent);
  //  TCanvas* c = new TCanvas (canvasName, "", 800, 1120);
  //  c->cd ();
  //  const double fuPad = 480./1120.;
  //  const double fdPad = 320./1120.;
  //  const double fcPad = 1.0 - fuPad - fdPad;
  //  TPad* uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 1.0-fuPad, 1.0, 1.0);
  //  TPad* cPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, fdPad, 1.0, 1.0-fuPad);
  //  TPad* dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, fdPad);

  //  uPad->SetBottomMargin (0);
  //  cPad->SetTopMargin (0);
  //  cPad->SetBottomMargin (0);
  //  dPad->SetTopMargin (0);
  //  dPad->SetBottomMargin (0.25);
  //  uPad->Draw ();
  //  cPad->Draw ();
  //  dPad->Draw ();

  //  const int iVar = GetVarN ("NoFcalMixCatVar");

  //  TH1D* h = nullptr; 
  //  TGAE* g = nullptr;

  //  uPad->cd (); 

  //  float ymin = -4;
  //  float ymax = 33;

  //  h = (TH1D*) h_jetInt_trk_dphi_ref[3]->Clone ("h");
  //  h->Reset ();
  //  h->GetYaxis ()->SetRangeUser (ymin, ymax);
  //  h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
  //  h->GetXaxis ()->SetTitleSize (0.028/fuPad);
  //  h->GetXaxis ()->SetLabelSize (0.028/fuPad);
  //  h->GetXaxis ()->SetTitleOffset (2.1*fuPad);
  //  h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#Delta#phi)");
  //  h->GetYaxis ()->SetTitleSize (0.028/fuPad);
  //  h->GetYaxis ()->SetLabelSize (0.028/fuPad);
  //  h->GetYaxis ()->SetTitleOffset (2.1*fuPad);

  //  h->SetLineWidth (1);
  //  h->SetLineStyle (2);
  //  h->DrawCopy ("hist ][");
  //  SaferDelete (&h);

  //  TBox* shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

  //  shadedBox = new TBox (0, ymin, M_PI/8., ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

  //  h = h_jetInt_trk_dphi_ref[3];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, kBlack, kOpenCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_ref_bkg[3];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myGreen, kOpenCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_ref_bkg_syst[3][iVar];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myOrange, kOpenSquare, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi[iCent][3];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, kBlack, kFullCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_bkg[iCent][3];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myRed, kFullCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_bkg_syst[iCent][3][iVar];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myBlue, kFullSquare, 0.8);
  //  SaferDelete (&g);

  //  myText (0.30, 0.83, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
  //  myText (0.30, 0.77, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.020/fuPad);
  //  myText (0.30, 0.71, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.020/fuPad);
  //  myText (0.30, 0.65, kBlack, "Jet-hadron correlations", 0.020/fuPad);
  //  myText (0.30, 0.59, kBlack, Form ("#it{p}_{T}^{jet} > %s, 2 < #it{p}_{T}^{ch} < 4 GeV", GetJetPtStr (tag).Data ()), 0.020/fuPad);
  //  myText (0.30, 0.53, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fuPad);


  //  cPad->cd (); 

  //  ymin = -4;
  //  ymax = 28;

  //  h = (TH1D*) h_jetInt_trk_dphi_ref_sig[3]->Clone ("h");
  //  h->Reset ();
  //  h->GetXaxis ()->SetMoreLogLabels ();
  //  h->GetYaxis ()->SetRangeUser (ymin, ymax);
  //  h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
  //  h->GetXaxis ()->SetTitleSize (0.028/fdPad);
  //  h->GetXaxis ()->SetLabelSize (0.028/fdPad);
  //  h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
  //  //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
  //  h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
  //  h->GetYaxis ()->SetTitleSize (0.028/fdPad);
  //  h->GetYaxis ()->SetLabelSize (0.028/fdPad);
  //  h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
  //  h->GetYaxis ()->CenterTitle ();

  //  h->SetLineWidth (1);
  //  h->SetLineStyle (2);
  //  h->DrawCopy ("hist ][");
  //  SaferDelete (&h);

  //  shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

  //  shadedBox = new TBox (0, ymin, M_PI/8., ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

  //  h = h_jetInt_trk_dphi_ref_sig[3];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myGreen, kOpenCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_ref_sig_syst[3][iVar];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myOrange, kOpenSquare, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_sig[iCent][3];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myRed, kFullCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_sig_syst[iCent][3][iVar];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myBlue, kFullSquare, 0.8);
  //  SaferDelete (&g);

  //  myLineText2 (0.36, 0.84, kBlack, kFullCircle, "#it{p}+Pb total", 0.8, 0.020/fcPad, true);
  //  myLineText2 (0.36, 0.75, myRed, kFullCircle, "#it{p}+Pb #bf{w/} FCal matching", 0.8, 0.020/fcPad, true);
  //  myLineText2 (0.36, 0.66, myBlue, kFullSquare, "#it{p}+Pb #bf{w/out} FCal matching", 0.8, 0.020/fcPad, true);
  //  myLineText2 (0.36, 0.57, kBlack, kOpenCircle, "#it{pp} total", 0.8, 0.020/fcPad);
  //  myLineText2 (0.36, 0.48, myGreen, kOpenCircle, "#it{pp} #bf{w/} FCal matching", 0.8, 0.020/fcPad);
  //  myLineText2 (0.36, 0.39, myOrange, kOpenSquare, "#it{pp} #bf{w/out} FCal matching", 0.8, 0.020/fcPad);


  //  dPad->cd (); 

  //  ymin = 0.5;
  //  ymax = 1.4;

  //  h = (TH1D*) h_jetInt_trk_dphi_iaa[iCent][3]->Clone ("h");
  //  h->Reset ();
  //  for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
  //  h->GetXaxis ()->SetMoreLogLabels ();
  //  h->GetYaxis ()->SetRangeUser (ymin, ymax);
  //  h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
  //  h->GetXaxis ()->SetTitleSize (0.028/fdPad);
  //  h->GetXaxis ()->SetLabelSize (0.028/fdPad);
  //  h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
  //  //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
  //  h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}Pb} = #it{p}+Pb / #it{pp}");
  //  h->GetYaxis ()->SetTitleSize (0.028/fdPad);
  //  h->GetYaxis ()->SetLabelSize (0.028/fdPad);
  //  h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
  //  h->GetYaxis ()->CenterTitle ();

  //  h->SetLineWidth (1);
  //  h->SetLineStyle (2);
  //  h->DrawCopy ("hist ][");
  //  SaferDelete (&h);

  //  shadedBox = new TBox (7.*M_PI/8., ymin, M_PI, ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (7.*M_PI/8., ymin, 7.*M_PI/8., ymax);

  //  shadedBox = new TBox (0, ymin, M_PI/8., ymax);
  //  shadedBox->SetFillColorAlpha (kGray, 0.3);
  //  shadedBox->Draw ();
  //  l->DrawLine (M_PI/8., ymin, M_PI/8., ymax);

  //  h = h_jetInt_trk_dphi_iaa[iCent][3];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myRed, kFullCircle, 0.8);
  //  SaferDelete (&g);

  //  h = h_jetInt_trk_dphi_iaa_syst[iCent][3][iVar];
  //  g = make_graph (h);
  //  ResetXErrors (g);
  //  myDraw (g, myBlue, kFullSquare, 0.8);
  //  SaferDelete (&g);

  //  myLineText2 (0.36, 0.48, myRed, kFullCircle, "#bf{w/} FCal matching", 0.8, 0.020/fdPad, true);
  //  myLineText2 (0.36, 0.38, myBlue, kFullSquare, "#bf{w/out} FCal matching", 0.8, 0.020/fdPad, true);

  //  c->SaveAs (Form ("%s/Plots/DPhi/JetTagged_HadronYields_%i-%iperc_FCalMixCatVar_dphi_gt2_lt4_%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], tag));
  //}


/*
  for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
    {
      const char* canvasName = Form ("c_jet_trk_dphi_%s_pp_syst", pTChSelections[iPtCh].Data ());
      TCanvas* c = new TCanvas (canvasName, "", 800, 800);

      TH1D* h = nullptr; 
      TGAE* g = nullptr;

      c->cd (); 

      float ymin = -20;
      float ymax = 20;

      c->Clear ();

      h = (TH1D*) h_jetInt_trk_dphi_ref[iPtCh]->Clone ("h");
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
        h = (TH1D*) h_jetInt_trk_dphi_ref_syst[iPtCh][iVar]->Clone ("htemp");
        SaveRelativeErrors (h, h_jetInt_trk_dphi_ref[iPtCh], true);
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
      myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, %s", GetJetPtStr (tag).Data (), pTChStrs[pTChSelections[iPtCh]].Data ()), 0.032);
      for (int iVar = 1; iVar < nVar; iVar++)
        myLineColorText (0.25, 0.77-iVar*0.040, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 1.0, 0.028);

      c->SaveAs (Form ("%s/Plots/Systematics/TotalJetTaggedYield_pp_dphi_%s_%s_syst.pdf", workPath.Data (), pTChSelections[iPtCh].Data (), tag));
    }
    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      const char* canvasName = Form ("c_jet_trk_dphi_%s_pPb_iCent%i_syst", pTChSelections[iPtCh].Data (), iCent);
      TCanvas* c = new TCanvas (canvasName, "", 800, 800);

      TH1D* h = nullptr; 
      TGAE* g = nullptr;

      c->cd (); 

      float ymin = -20;
      float ymax = 20;

      c->Clear ();

      h = (TH1D*) h_jetInt_trk_dphi[iCent][iPtCh]->Clone ("h");
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
        h = (TH1D*) h_jetInt_trk_dphi_syst[iCent][iPtCh][iVar]->Clone ("htemp");
        SaveRelativeErrors (h, h_jetInt_trk_dphi[iCent][iPtCh], true);
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

      c->SaveAs (Form ("%s/Plots/Systematics/TotalJetTaggedYield_%s_dphi_%s_%s_syst.pdf", workPath.Data (), Form ("pPb_%i-%iperc", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), pTChSelections[iPtCh].Data (), tag));
    }
  }



  for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
    {
      const char* canvasName = Form ("c_jet_trk_dphi_%s_pp_sig_syst", pTChSelections[iPtCh].Data ());
      TCanvas* c = new TCanvas (canvasName, "", 800, 800);

      TH1D* h = nullptr; 
      TGAE* g = nullptr;

      c->cd (); 

      float ymin = -20;
      float ymax = 20;

      c->Clear ();

      h = (TH1D*) h_jetInt_trk_dphi_ref_sig[iPtCh]->Clone ("h");
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
        h = (TH1D*) h_jetInt_trk_dphi_ref_sig_syst[iPtCh][iVar]->Clone ("htemp");
        SaveRelativeErrors (h, h_jetInt_trk_dphi_ref_sig[iPtCh], true);
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
      myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %s, %s", GetJetPtStr (tag).Data (), pTChStrs[pTChSelections[iPtCh]].Data ()), 0.032);
      for (int iVar = 1; iVar < nVar; iVar++)
        myLineColorText (0.25, 0.77-iVar*0.040, varStyles[variations[iVar]].first, varStyles[variations[iVar]].second, varFullNames[variations[iVar]], 1.0, 0.028);
      
      c->SaveAs (Form ("%s/Plots/Systematics/SignalJetTaggedYield_pp_dphi_%s_%s_syst.pdf", workPath.Data (), pTChSelections[iPtCh].Data (), tag));
    }
    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      const char* canvasName = Form ("c_jet_trk_dphi_%s_pPb_sig_iCent%i_syst", pTChSelections[iPtCh].Data (), iCent);
      TCanvas* c = new TCanvas (canvasName, "", 800, 800);

      TH1D* h = nullptr; 
      TGAE* g = nullptr;

      c->cd (); 

      float ymin = -20;
      float ymax = 20;

      c->Clear ();

      h = (TH1D*) h_jetInt_trk_dphi_sig[iCent][iPtCh]->Clone ("h");
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
        h = (TH1D*) h_jetInt_trk_dphi_sig_syst[iCent][iPtCh][iVar]->Clone ("htemp");
        SaveRelativeErrors (h, h_jetInt_trk_dphi_sig[iCent][iPtCh], true);
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
      
      c->SaveAs (Form ("%s/Plots/Systematics/SignalJetTaggedYield_%s_dphi_%s_%s_syst.pdf", workPath.Data (), Form ("pPb_%i-%iperc", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), pTChSelections[iPtCh].Data (), tag));
    }
  }



  for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      const char* canvasName = Form ("c_jet_trk_dphi_%s_iaa_iCent%i_syst", pTChSelections[iPtCh].Data (), iCent);
      TCanvas* c = new TCanvas (canvasName, "", 800, 800);

      TH1D* h = nullptr; 
      TGAE* g = nullptr;

      c->cd (); 

      float ymin = -20;
      float ymax = 20;


      h = (TH1D*) h_jetInt_trk_dphi_iaa[iCent][iPtCh]->Clone ("h");
      h->Reset ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetXaxis ()->SetTitle ("#Delta#phi_{ch,jet}");
      //h->GetXaxis ()->SetTitleSize (0.028);
      //h->GetXaxis ()->SetLabelSize (0.028);
      h->GetYaxis ()->SetTitle ("#delta I_{pPb} / I_{pPb} [%]");
      //h->GetYaxis ()->SetTitleSize (0.028);
      //h->GetYaxis ()->SetLabelSize (0.028);

      h->SetLineWidth (1);
      h->SetLineStyle (2);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      for (int iVar = 1; iVar < nVar; iVar++) {
        h = (TH1D*) h_jetInt_trk_dphi_iaa_syst[iCent][iPtCh][iVar]->Clone ("htemp");
        SaveRelativeErrors (h, h_jetInt_trk_dphi_iaa[iCent][iPtCh], true);
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

      c->SaveAs (Form ("%s/Plots/Systematics/JetTagged_IpPb_%i-%iperc_dphi_%s_%s_syst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], pTChSelections[iPtCh].Data (), tag));
    }
  }
*/

}


#endif
