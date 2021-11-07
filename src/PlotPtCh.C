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


//const bool makeTotalSystPlots = false;
const bool makeBkgdSystPlots  = false;
const bool makeSigSystPlots   = true;
const bool makeIpPbSystPlots  = true;
const bool makeMCClosurePlots = true;
const bool makeUnfoldingPlots = true;

const float maxDataSyst = 20; // maximum y-axis for data-driven systematics
const float maxMCSyst = 10; // maximum y-axis for MC-driven systematics


//void DoOffset (TH1D* h, const char* tag, const int iCent) {
//
//  ifstream offsetsFile;
//  offsetsFile.open (Form ("%s/aux/%s_IAAOffsets.dat", workPath.Data (), tag));
//
//  string inTag, pTChSel, centStr, offpercerr;
//  double offset = 0, offerr = 0;
//
//  while (!offsetsFile.eof ()) {
//    offsetsFile >> inTag >> pTChSel >> centStr >> offset >> offerr >> offpercerr;
//
//    if (TString (centStr.c_str ()) == TString (Form ("%i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]))) {
//
//      for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
//        if (pTChStrCuts[TString (pTChSel.c_str ())].first <= h->GetBinCenter (iX) && h->GetBinCenter (iX) <= pTChStrCuts[TString (pTChSel.c_str ())].second)
//          h->SetBinContent (iX, h->GetBinContent (iX)-offset);
//      }
//    }
//  }
//
//  return;
//}
//
//
//
//void DoOffset (TGAE* g, const char* tag, const int iCent) {
//
//  ifstream offsetsFile;
//  offsetsFile.open (Form ("%s/aux/%s_IAAOffsets.dat", workPath.Data (), tag));
//
//  string inTag, pTChSel, centStr, offpercerr;
//  double offset = 0, offerr = 0;
//
//  while (!offsetsFile.eof ()) {
//    offsetsFile >> inTag >> pTChSel >> centStr >> offset >> offerr >> offpercerr;
//
//    if (TString (centStr.c_str ()) == Form ("%i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent])) {
//
//      double x, y;
//      for (int i = 0; i < g->GetN (); i++) {
//        g->GetPoint (i, x, y);
//        if (pTChStrCuts[TString (pTChSel.c_str ())].first <= x && x <= pTChStrCuts[TString (pTChSel.c_str ())].second)
//          g->SetPoint (i, x, y-offset);
//      }
//    }
//  }
//
//  return;
//}


void PlotPtCh (const char* inFileTag) {

  TLine* l = new TLine ();
  TLatex* tl = new TLatex ();

  TFile* inFile = nullptr;

  TH1D****  h_jetInt_trk_pt_ref       = Get3DArray <TH1D*> (2, 2, nDir);
  TH1D****  h_jetInt_trk_pt_ref_bkg   = Get3DArray <TH1D*> (2, 2, nDir);
  TH1D***** h_jetInt_trk_pt           = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);
  TH1D***** h_jetInt_trk_pt_bkg       = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);

  TH1D****  h_jetInt_trk_pt_ref_sig   = Get3DArray <TH1D*> (2, 2, nDir);
  TH1D***** h_jetInt_trk_pt_sig       = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);

  TH1D****  h_jetInt_trk_pt_ref_unf   = Get3DArray <TH1D*> (2, 2, nDir);
  TH1D***** h_jetInt_trk_pt_unf       = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);

  TH1D***** h_jetInt_trk_pt_iaa       = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);
  TH1D***** h_jetInt_trk_pt_iaaNoUnf  = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);


  //TGAE****  g_jetInt_trk_pt_ref_syst     = Get3DArray <TGAE*> (2, nDir, nVar);
  //TGAE****  g_jetInt_trk_pt_ref_bkg_syst = Get3DArray <TGAE*> (2, nDir, nVar);
  //TGAE***** g_jetInt_trk_pt_syst         = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, nVar);
  //TGAE***** g_jetInt_trk_pt_bkg_syst     = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, nVar);

  TGAE****  g_jetInt_trk_pt_ref_sig_syst = Get3DArray <TGAE*> (2, nDir, nVar);
  TGAE***** g_jetInt_trk_pt_sig_syst     = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, nVar);

  TGAE****  g_jetInt_trk_pt_ref_unf_syst = Get3DArray <TGAE*> (2, nDir, nVar);
  TGAE***** g_jetInt_trk_pt_unf_syst     = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, nVar);

  TGAE***** g_jetInt_trk_pt_iaa_syst     = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, nVar);


  //TGAE****  g_jetInt_trk_pt_ref_systTot     = Get3DArray <TGAE*> (2, nDir, 3);
  //TGAE****  g_jetInt_trk_pt_ref_bkg_systTot = Get3DArray <TGAE*> (2, nDir, 3);
  //TGAE***** g_jetInt_trk_pt_systTot         = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, 3);
  //TGAE***** g_jetInt_trk_pt_bkg_systTot     = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, 3);

  TGAE****  g_jetInt_trk_pt_ref_sig_systTot = Get3DArray <TGAE*> (2, nDir, 3);
  TGAE***** g_jetInt_trk_pt_sig_systTot     = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, 3);

  TGAE****  g_jetInt_trk_pt_ref_unf_systTot = Get3DArray <TGAE*> (2, nDir, 3);
  TGAE***** g_jetInt_trk_pt_unf_systTot     = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, 3);

  TGAE***** g_jetInt_trk_pt_iaa_systTot     = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, 3);


  //TH1D****  h_jetInt_trk_pt_ref_syst      = Get3DArray <TH1D*> (2, nDir, nVar);
  //TH1D****  h_jetInt_trk_pt_ref_bkg_syst  = Get3DArray <TH1D*> (2, nDir, nVar);
  //TH1D***** h_jetInt_trk_pt_syst          = Get4DArray <TH1D*> (2, nDir, nZdcCentBins+1, nVar);
  //TH1D***** h_jetInt_trk_pt_bkg_syst      = Get4DArray <TH1D*> (2, nDir, nZdcCentBins+1, nVar);

  TH1D****  h_jetInt_trk_pt_ref_sig_syst  = Get3DArray <TH1D*> (2, nDir, nVar);
  TH1D***** h_jetInt_trk_pt_sig_syst      = Get4DArray <TH1D*> (2, nDir, nZdcCentBins+1, nVar);

  TH1D****  h_jetInt_trk_pt_ref_unf_syst  = Get3DArray <TH1D*> (2, nDir, nVar);
  TH1D***** h_jetInt_trk_pt_unf_syst      = Get4DArray <TH1D*> (2, nDir, nZdcCentBins+1, nVar);

  TH1D***** h_jetInt_trk_pt_iaa_syst      = Get4DArray <TH1D*> (2, nDir, nZdcCentBins+1, nVar);


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

          h_jetInt_trk_pt_ref[iDType][iPtJInt][iDir]      = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_%s_%s_Nominal",      dir.Data (), dType.Data (), pTJInt.Data ()));
          h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][iDir]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_bkg_%s_%s_Nominal",  dir.Data (), dType.Data (), pTJInt.Data ()));
          h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_sig_%s_%s_Nominal",  dir.Data (), dType.Data (), pTJInt.Data ()));

        } // end loop over iDir


        for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

          for (int iDir = 0; iDir < nDir; iDir++) {

            const TString dir = directions[iDir];

            h_jetInt_trk_pt[iDType][iPtJInt][iDir][iCent]           = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_%s_%s_%s_Nominal",     dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));
            h_jetInt_trk_pt_bkg[iDType][iPtJInt][iDir][iCent]       = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_bkg_%s_%s_%s_Nominal", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));
            h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent]       = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_sig_%s_%s_%s_Nominal", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));

          } // end loop over iDir

        } // end loop over iCent

      } // end loop over iPtJInt

    } // end loop over iDType



    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      for (int iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        for (int iVar = 0; iVar < nVar; iVar++) {

          const TString var = variations[iVar];

          //g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][iVar]     = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_ref_syst_%s_%s",      dir.Data (), pTJInt.Data (), var.Data ()));
          //g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][iVar] = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_ref_bkg_syst_%s_%s",  dir.Data (), pTJInt.Data (), var.Data ()));
          g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar] = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_ref_sig_syst_%s_%s",  dir.Data (), pTJInt.Data (), var.Data ()));

          for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            //g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][iVar]      = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_syst_%s_%s_%s",     dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            //g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][iVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_bkg_syst_%s_%s_%s", dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_sig_syst_%s_%s_%s", dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iVar

        for (int iTotVar = 0; iTotVar < 3; iTotVar++) {

          const TString totVar = totalVariations[iTotVar];

          //g_jetInt_trk_pt_ref_systTot[iPtJInt][iDir][iTotVar]     = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_ref_%s_systTot_%s",      dir.Data (), totVar.Data (), pTJInt.Data ()));
          //g_jetInt_trk_pt_ref_bkg_systTot[iPtJInt][iDir][iTotVar] = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_ref_bkg_%s_systTot_%s",  dir.Data (), totVar.Data (), pTJInt.Data ()));
          g_jetInt_trk_pt_ref_sig_systTot[iPtJInt][iDir][iTotVar] = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_ref_sig_%s_systTot_%s",  dir.Data (), totVar.Data (), pTJInt.Data ()));

          for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            //g_jetInt_trk_pt_systTot[iPtJInt][iDir][iCent][iTotVar]      = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_%s_systTot_%s_%s",     dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            //g_jetInt_trk_pt_bkg_systTot[iPtJInt][iDir][iCent][iTotVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_bkg_%s_systTot_%s_%s", dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            g_jetInt_trk_pt_sig_systTot[iPtJInt][iDir][iCent][iTotVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_sig_%s_systTot_%s_%s", dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));

          } // end loop over iCent

        } // end loop over iTotVar


        for (int iVar = 1; iVar < nVar; iVar++) {

          const TString var = variations[iVar];
          const TString dType = (dataVariations.count (var) > 0 ? "data" : "mc");

          //h_jetInt_trk_pt_ref_syst[iPtJInt][iDir][iVar]     = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_%s_%s_%s",     dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
          //h_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][iVar] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_bkg_%s_%s_%s", dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
          h_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_sig_%s_%s_%s", dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()));

          for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            //h_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][iVar]      = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_%s_%s_%s_%s",      dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
            //h_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][iVar]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_bkg_%s_%s_%s_%s",  dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
            h_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_sig_%s_%s_%s_%s",  dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iVar

      } // end loop over iDir

    } // end loop over iPtJInt

  }




  {
    TString inFileName = inFileTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/ProcessUnfolding_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    for (int iDType = 0; iDType < 2; iDType++) {

      const TString dType = (iDType == 0 ? "data" : "mc");

      for (short iPtJInt : {0, 1}) {

        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

        for (int iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_unf_%s_%s_Nominal",  dir.Data (), dType.Data (), pTJInt.Data ()));

        } // end loop over iDir


        for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

          for (int iDir = 0; iDir < nDir; iDir++) {

            const TString dir = directions[iDir];

            h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent]       = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_unf_%s_%s_%s_Nominal", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));
            h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent]       = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_iaa_%s_%s_%s_Nominal",     dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));
            h_jetInt_trk_pt_iaaNoUnf[iDType][iPtJInt][iDir][iCent]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_iaaNoUnf_%s_%s_%s_Nominal",  dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));

          } // end loop over iDir

        } // end loop over iCent

      } // end loop over iPtJInt

    } // end loop over iDType



    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      for (int iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        for (int iVar = 0; iVar < nVar; iVar++) {

          const TString var = variations[iVar];

          g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar] = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_ref_unf_syst_%s_%s",  dir.Data (), pTJInt.Data (), var.Data ()));

          for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_unf_syst_%s_%s_%s", dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_iaa_syst_%s_%s_%s", dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iVar

        for (int iTotVar = 0; iTotVar < 3; iTotVar++) {

          const TString totVar = totalVariations[iTotVar];

          g_jetInt_trk_pt_ref_unf_systTot[iPtJInt][iDir][iTotVar] = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_ref_unf_%s_systTot_%s",  dir.Data (), totVar.Data (), pTJInt.Data ()));

          for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            g_jetInt_trk_pt_unf_systTot[iPtJInt][iDir][iCent][iTotVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_unf_%s_systTot_%s_%s", dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            g_jetInt_trk_pt_iaa_systTot[iPtJInt][iDir][iCent][iTotVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_iaa_%s_systTot_%s_%s", dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));

          } // end loop over iCent

        } // end loop over iTotVar


        for (int iVar = 1; iVar < nVar; iVar++) {

          const TString var = variations[iVar];
          const TString dType = (dataVariations.count (var) > 0 ? "data" : "mc");

          h_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_unf_%s_%s_%s", dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()));

          for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            h_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_unf_%s_%s_%s_%s",  dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
            h_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_iaa_%s_%s_%s_%s",      dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iVar

      } // end loop over iDir

    } // end loop over iPtJInt

  }



  l->SetLineStyle (7);
  l->SetLineWidth (1);
  l->SetLineColor (kBlack);



  for (int iDType = 0; iDType < 2; iDType++) {

    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        const char* canvasName = Form ("c_jetInt_trk_pt_iCent%i_%s", iCent, pTJInt.Data ());

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
        h = new TH1D ("h", ";;(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->GetYaxis ()->SetTitleSize (0.028/fuPad);
        h->GetYaxis ()->SetLabelSize (0.028/fuPad);
        h->GetYaxis ()->SetTitleOffset (3.0*fuPad);

        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        //g = (TGAE*) g_jetInt_trk_pt_ref_syst[iPtJInt][0][0]->Clone ();
        h = h_jetInt_trk_pt_ref[iDType][iPtJInt][0];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myBlue);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myBlue, kFullCircle, 0.8);
        SaferDelete (&g);

        //g = (TGAE*) g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][0][0]->Clone ();
        h = h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][0];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myPurple);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myPurple, kOpenCircle, 0.8);
        SaferDelete (&g);

        //g = (TGAE*) g_jetInt_trk_pt_syst[iPtJInt][0][iCent][0]->Clone ();
        h = h_jetInt_trk_pt[iDType][iPtJInt][0][iCent];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myRed);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myRed, kFullCircle, 0.8);
        SaferDelete (&g);

        //g = (TGAE*) g_jetInt_trk_pt_bkg_syst[iPtJInt][0][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_bkg[iDType][iPtJInt][0][iCent];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myGreen);
        //SaferDelete (&g);
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

        h = new TH1D ("h", ";;(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->GetYaxis ()->SetTitleSize (0.028/fuPad);
        h->GetYaxis ()->SetLabelSize (0.028/fuPad);
        h->GetYaxis ()->SetTitleOffset (3.0*fuPad);

        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        //g = (TGAE*) g_jetInt_trk_pt_ref_syst[iPtJInt][2][0]->Clone ();
        h = h_jetInt_trk_pt_ref[iDType][iPtJInt][2];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myBlue);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myBlue, kFullCircle, 0.8);
        SaferDelete (&g);

        //g = (TGAE*) g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][2][0]->Clone ();
        h = h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][2];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myPurple);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myPurple, kOpenCircle, 0.8);
        SaferDelete (&g);

        //g = (TGAE*) g_jetInt_trk_pt_syst[iPtJInt][2][iCent][0]->Clone ();
        h = h_jetInt_trk_pt[iDType][iPtJInt][2][iCent];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myRed);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myRed, kFullCircle, 0.8);
        SaferDelete (&g);

        //g = (TGAE*) g_jetInt_trk_pt_bkg_syst[iPtJInt][2][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_bkg[iDType][iPtJInt][2][iCent];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myGreen);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myGreen, kOpenCircle, 0.8);
        SaferDelete (&g);

        if (iDType == 0) myText (0.58, 0.78, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
        else {
          myText (0.42, 0.78, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.022/fuPad);
          myText (0.42, 0.72, kBlack, "Pythia8 (+ #it{p}+Pb Overlay)", 0.022/fuPad);
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

        h = new TH1D ("h", ";;(Sig.+Bkg.) - Bkg.", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->GetYaxis ()->SetTitleSize (0.028/fdPad);
        h->GetYaxis ()->SetLabelSize (0.028/fdPad);
        h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
        h->GetYaxis ()->CenterTitle ();

        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        //g = (TGAE*) g_jetInt_trk_pt_ref_sig_syst[iPtJInt][0][0]->Clone ();
        h = h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][0];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myBlue);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myBlue, kFullCircle, 0.8);
        SaferDelete (&g);

        //g = (TGAE*) g_jetInt_trk_pt_sig_syst[iPtJInt][0][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_sig[iDType][iPtJInt][0][iCent];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myRed);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myRed, kFullCircle, 0.8);
        SaferDelete (&g);


        //g = (TGAE*) g_jetInt_trk_pt_ref_unf_syst[iPtJInt][0][0]->Clone ();
        h = h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][0];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myBlue);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myBlue, kOpenCircle, 0.8);
        SaferDelete (&g);

        //g = (TGAE*) g_jetInt_trk_pt_unf_syst[iPtJInt][0][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_unf[iDType][iPtJInt][0][iCent];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myRed);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myRed, kOpenCircle, 0.8);
        SaferDelete (&g);

        myText (0.24, 0.26, kBlack, "Jet-hadron correlations", 0.020/fcPad);
        myText (0.24, 0.18, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV", iPtJInt == 0 ? 30 : 60), 0.020/fcPad);
        myText (0.24, 0.10, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fcPad);


        crPad->cd (); 
        crPad->SetLogx ();
        crPad->SetLogy ();

        ymin = 8e-6;
        ymax = 3e1;

        h = new TH1D ("h", ";;(Sig.+Bkg.) - Bkg.", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->GetYaxis ()->SetTitleSize (0.028/fdPad);
        h->GetYaxis ()->SetLabelSize (0.028/fdPad);
        h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
        h->GetYaxis ()->CenterTitle ();

        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        //g = (TGAE*) g_jetInt_trk_pt_ref_sig_syst[iPtJInt][2][0]->Clone ();
        h = h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][2];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myBlue);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myBlue, kFullCircle, 0.8);
        SaferDelete (&g);

        //g = (TGAE*) g_jetInt_trk_pt_sig_syst[iPtJInt][2][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_sig[iDType][iPtJInt][2][iCent];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myRed);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myRed, kFullCircle, 0.8);
        SaferDelete (&g);


        //g = (TGAE*) g_jetInt_trk_pt_ref_unf_syst[iPtJInt][2][0]->Clone ();
        h = h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][2];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myBlue);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myBlue, kOpenCircle, 0.8);
        SaferDelete (&g);

        //g = (TGAE*) g_jetInt_trk_pt_unf_syst[iPtJInt][2][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_unf[iDType][iPtJInt][2][iCent];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myRed);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myRed, kOpenCircle, 0.8);
        SaferDelete (&g);


        dlPad->cd (); 
        dlPad->SetLogx ();

        ymin = 0.53;
        ymax = 1.47;//(strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.45 : 1.17);

        h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#it{I}_{#it{p}Pb} = #it{p}+Pb / #it{pp}", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
        h->SetBinContent (1, 1);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->GetXaxis ()->SetTitleSize (0.028/fdPad);
        h->GetXaxis ()->SetLabelSize (0.028/fdPad);
        h->GetXaxis ()->SetTitleOffset (3.7*fdPad);
        h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
        h->GetYaxis ()->SetTitleSize (0.028/fdPad);
        h->GetYaxis ()->SetLabelSize (0.028/fdPad);
        h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
        h->GetYaxis ()->CenterTitle ();

        h->SetLineWidth (1);
        h->SetLineStyle (2);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        //g = (TGAE*) g_jetInt_trk_pt_iaa_syst[iPtJInt][0][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_iaa[iDType][iPtJInt][0][iCent];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myRed);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myRed, kFullCircle, 0.8);
        SaferDelete (&g);


        drPad->cd (); 
        drPad->SetLogx ();

        ymin = 0.53;
        ymax = 1.47;//(strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.45 : 1.17);

        h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#it{I}_{#it{p}Pb} = #it{p}+Pb / #it{pp}", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
        h->SetBinContent (1, 1);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->GetXaxis ()->SetTitleSize (0.028/fdPad);
        h->GetXaxis ()->SetLabelSize (0.028/fdPad);
        h->GetXaxis ()->SetTitleOffset (3.7*fdPad);
        h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
        h->GetYaxis ()->SetTitleSize (0.028/fdPad);
        h->GetYaxis ()->SetLabelSize (0.028/fdPad);
        h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
        h->GetYaxis ()->CenterTitle ();

        h->SetLineWidth (1);
        h->SetLineStyle (2);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        //g = (TGAE*) g_jetInt_trk_pt_iaa_syst[iPtJInt][2][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_iaa[iDType][iPtJInt][2][iCent];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myRed);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myRed, kFullCircle, 0.8);
        SaferDelete (&g);

        c->SaveAs (Form ("%s/Plots/PtCh/JetTagged_HadronYields_%i-%iperc_comparison_PtCh_%iGeVJets%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], iPtJInt == 0 ? 30 : 60, iDType == 1 ? "_mc" : "")); 

      } // end loop over iCent

    } // end loop over iPtJInt

  } // end loop over iDType


/*
  for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
    const char* canvasName = Form ("c_jetInt_trk_pt_FcalVsZdc_iCent%i_%s", iCent, pTJInt.Data ());
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

    const int iVar = GetVarN ("FcalCentVar");

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    ulPad->cd (); 
    ulPad->SetLogx ();
    ulPad->SetLogy ();

    float ymin = 8e-6;
    float ymax = 3e1;

    h = new TH1D ("h", ";;(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fuPad);

    h->SetLineWidth (0);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jetInt_trk_pt_ref[0][iPtJInt][0];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_ref_bkg[0][iPtJInt][0];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kOpenCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt[0][iPtJInt][0][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_bkg[0][iPtJInt][0][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myGreen, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_syst[0][iPtJInt][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_bkg_syst[0][iPtJInt][iCent][iVar];
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

    h = new TH1D ("h", ";;(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fuPad);

    h->SetLineWidth (0);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jetInt_trk_pt_ref[0][iPtJInt][2];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_ref_bkg[0][iPtJInt][2];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kOpenCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt[0][iPtJInt][2][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_bkg[0][iPtJInt][2][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myGreen, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_syst[2][iPtJInt][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_bkg_syst[2][iPtJInt][iCent][iVar];
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

    h = new TH1D ("h", ";;(Sig.+Bkg.) - Bkg.", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (0);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jetInt_trk_pt_ref_sig[0][iPtJInt][0];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_sig[0][iPtJInt][0][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_sig_syst[0][iPtJInt][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);

    myText (0.24, 0.26, kBlack, "Jet-hadron correlations", 0.020/fcPad);
    myText (0.24, 0.18, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV", iPtJInt == 0 ? 30 : 60), 0.020/fcPad);
    myText (0.24, 0.10, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fcPad);


    crPad->cd (); 
    crPad->SetLogx ();
    crPad->SetLogy ();

    ymin = 8e-6;
    ymax = 3e1;

    h = new TH1D ("h", ";;(Sig.+Bkg.) - Bkg.", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (0);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jetInt_trk_pt_ref_sig[0][iPtJInt][2];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_sig[0][iPtJInt][2][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_sig_syst[2][iPtJInt][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);


    dlPad->cd (); 
    dlPad->SetLogx ();

    ymin = 0.53;
    ymax = (strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.45 : 1.17);

    h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#it{I}_{#it{p}Pb} = #it{p}+Pb / #it{pp}", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
    h->SetBinContent (1, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.7*fdPad);
    h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jetInt_trk_pt_iaa[0][iPtJInt][0][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_iaa_syst[0][iPtJInt][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);


    drPad->cd (); 
    drPad->SetLogx ();

    ymin = 0.53;
    ymax = (strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.45 : 1.17);

    h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#it{I}_{#it{p}Pb} = #it{p}+Pb / #it{pp}", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
    h->SetBinContent (1, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.7*fdPad);
    h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jetInt_trk_pt_iaa[0][iPtJInt][2][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_iaa_syst[2][iPtJInt][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);

    c->SaveAs (Form ("%s/Plots/PtCh/JetTagged_HadronYields_%i-%iperc_FCalvsZDC_PtCh_%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], tag)); 
  }



  for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
    const char* canvasName = Form ("c_jetInt_trk_pt_NoFcalMixCatVar_iCent%i_%s", iCent, pTJInt.Data ());
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

    const int iVar = GetVarN ("NoFcalMixCatVar");

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    ulPad->cd (); 
    ulPad->SetLogx ();
    ulPad->SetLogy ();

    float ymin = 8e-6;
    float ymax = 3e1;

    h = new TH1D ("h", ";;(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fuPad);

    h->SetLineWidth (0);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jetInt_trk_pt_ref[0][iPtJInt][0];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kOpenCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_ref_bkg[0][iPtJInt][0];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myGreen, kOpenCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_ref_bkg_syst[0][iPtJInt][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myOrange, kOpenSquare, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt[0][iPtJInt][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_bkg[0][iPtJInt][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_bkg_syst[0][iPtJInt][iCent][iVar];
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

    h = new TH1D ("h", ";;(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]");
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fuPad);

    h->SetLineWidth (0);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jetInt_trk_pt_ref[0][iPtJInt][2];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kOpenCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_ref_bkg[0][iPtJInt][2];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myGreen, kOpenCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_ref_bkg_syst[2][iPtJInt][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myOrange, kOpenSquare, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt[0][iPtJInt][2][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, kBlack, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_bkg[0][iPtJInt][2][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_bkg_syst[2][iPtJInt][iCent][iVar];
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

    h = new TH1D ("h", ";;(Sig.+Bkg.) - Bkg.", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (0);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jetInt_trk_pt_ref_sig[0][iPtJInt][0];
    g = make_graph (h);
    ResetXErrors (g); myDraw (g, myGreen, kOpenCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_ref_sig_syst[0][iPtJInt][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myOrange, kOpenSquare, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_sig[0][iPtJInt][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_sig_syst[0][iPtJInt][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);

    myText (0.24, 0.26, kBlack, "Jet-hadron correlations", 0.020/fcPad);
    myText (0.24, 0.18, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV", iPtJInt == 0 ? 30 : 60), 0.020/fcPad);
    myText (0.24, 0.10, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fcPad);


    crPad->cd (); 
    crPad->SetLogx ();
    crPad->SetLogy ();

    ymin = 8e-6;
    ymax = 3e1;

    h = new TH1D ("h", ";;(Sig.+Bkg.) - Bkg.", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (0);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jetInt_trk_pt_ref_sig[0][iPtJInt][2];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myGreen, kOpenCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_ref_sig_syst[2][iPtJInt][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myOrange, kOpenSquare, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_sig[0][iPtJInt][2][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_sig_syst[2][iPtJInt][iCent][iVar];
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

    h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#it{I}_{#it{p}Pb} = #it{p}+Pb / #it{pp}", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
    h->SetBinContent (1, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.7*fdPad);
    h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (0);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jetInt_trk_pt_iaa[0][iPtJInt][0][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_iaa_syst[0][iPtJInt][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);


    drPad->cd (); 
    drPad->SetLogx ();

    ymin = 0.53;
    ymax = (strcmp (tag, "30GeVJets") == 0 || strcmp (tag, "15GeVJets") == 0 ? 1.45 : 1.17);

    h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#it{I}_{#it{p}Pb} = #it{p}+Pb / #it{pp}", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
    h->SetBinContent (1, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.7*fdPad);
    h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (3.0*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (0);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jetInt_trk_pt_iaa[0][iPtJInt][2][iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    h = h_jetInt_trk_pt_iaa_syst[2][iPtJInt][iCent][iVar];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullSquare, 0.8);
    SaferDelete (&g);

    myLineText2 (0.10, 0.48, myRed, kFullCircle, "#bf{w/} FCal matching", 0.8, 0.020/fdPad, true);
    myLineText2 (0.10, 0.38, myBlue, kFullSquare, "#bf{w/out} FCal matching", 0.8, 0.020/fdPad, true);

    c->SaveAs (Form ("%s/Plots/PtCh/JetTagged_HadronYields_%i-%iperc_FCalMixCatVar_PtCh_%iGeVJets.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], iPtJInt == 0 ? 30 : 60)); 
  }*/




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (int iDir : {0, 1, 2}) {
  
      const char* canvasName = Form ("c_jetInt_trk_pt_signalYields_%s_%s", directions[iDir].Data (), pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 800, 1000);
      c->cd ();
  
      c->SetTopMargin (0.04);
      c->SetBottomMargin (0.12);
      c->SetLeftMargin (0.12);
      c->SetRightMargin (0.03);

      TH1D* h = nullptr;
      TGAE* g = nullptr;

      c->SetLogx ();
      c->SetLogy ();

      float ymin = 2e-9;
      float ymax = 3e4;
      h = new TH1D ("htemp", ";#it{p}_{T}^{ch} [GeV];(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->GetXaxis ()->SetTitleSize (0.036);
      h->GetXaxis ()->SetLabelSize (0.036);
      h->GetYaxis ()->SetTitleSize (0.036);
      h->GetYaxis ()->SetLabelSize (0.036);
      h->GetYaxis ()->SetTitleOffset (1.5);

      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);


      h = h_jetInt_trk_pt_ref_unf[0][iPtJInt][iDir];
      g = (TGAE*) g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][0]->Clone ();
      SetCentralValuesKeepRelativeErrors (g, h);
      ScaleGraph (g, nullptr, std::pow (10, 3));
      myDrawSystFill (g, colorfulSystColors[0], 1, 1001);
      SaferDelete (&g);

      g = make_graph (h);
      ScaleGraph (g, nullptr, std::pow (10, 3));
      ResetXErrors (g);
      myDraw (g, colorfulColors[0], kFullCircle, 1.4, 1, 2, "P", false);
      SaferDelete (&g);


      for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        h = (TH1D*) h_jetInt_trk_pt_ref_unf[0][iPtJInt][iDir]->Clone ("htemp");
        h->Scale (std::pow (10, 2-iCent));
        myDrawHist (h, kBlack, 1, 3);
        SaferDelete (&h);

        g = (TGAE*) g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_unf[0][iPtJInt][iDir][iCent];
        SetCentralValuesKeepRelativeErrors (g, h);
        ScaleGraph (g, nullptr, std::pow (10, 2-iCent));
        myDrawSystFill (g, colorfulSystColors[iCent+1], 1.0, 1001);
        SaferDelete (&g);
  
        g = make_graph (h);
        ResetXErrors (g);
        ScaleGraph (g, nullptr, std::pow (10, 2-iCent));
        myDraw (g, colorfulColors[iCent+1], kFullCircle, 1.4, 1, 2, "P", false);
        SaferDelete (&g);

      } // end loop over iCent

      myText (0.56, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
      myText (0.56, 0.86, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} %s", iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.032);

      myText (0.20, 0.36, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.030);
      myText (0.20, 0.32, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.030);

      mySimpleMarkerAndBoxAndLineText (0.27, 0.27, 1.4, 1001, colorfulSystColors[0], 1.0, colorfulColors[0], kFullCircle, 1.6, "#it{pp} (#times10^{3})", 0.030);
      for (int iCent = 0; iCent < nZdcCentBins; iCent++) { 
        mySimpleMarkerAndBoxAndLineText (0.27 + (iCent >= 2 ? 0.4 : 0), 0.27-((iCent+1)%3)*0.040, 1.4, 1001, colorfulSystColors[iCent+1], 1.0, colorfulColors[iCent+1], kFullCircle, 1.6, Form ("ZDC %i-%i%% (#times10^{%i})", zdcCentPercs[iCent+1], zdcCentPercs[iCent], 2-iCent), 0.030);
      }
      mySimpleMarkerAndBoxAndLineText (0.67, 0.15, 1.4, 1001, colorfulSystColors[nZdcCentBins+1], 1.0, colorfulColors[nZdcCentBins+1], kFullCircle, 1.6, Form ("All cent. (#times10^{%i})", 2-nZdcCentBins), 0.030);
      mySimpleMarkerAndBoxAndLineText (0.27, 0.15, 1.4, 1001, kWhite, 0.0, kBlack, kDot, 0.0, "#it{pp} (#it{scaled to} #it{p}+Pb)", 0.030);

      c->RedrawAxis();

      c->SaveAs (Form ("%s/Plots/PtCh/SignalJetTaggedYield_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
  
    } // end loop over iDir

  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (int iDir : {0, 1, 2}) {
      const char* canvasName = Form ("c_jetInt_trk_pt_IpPb_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1200, 800);
      c->Divide (3, 2);

      for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        gPad->SetLogx ();

        TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#it{I}_{#it{p}Pb}", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser (0.74, 1.4);
        //h->GetYaxis ()->SetRangeUser (0.0, 3);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        double x, y;

        TGAE* g = (TGAE*) g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_iaa[0][iPtJInt][iDir][iCent];
        //h = h_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iMCTruthJetsTruthParts];
        SetCentralValuesKeepRelativeErrors (g, h);
        //if (strcmp (tag, "30GeVJets") == 0)
        //  TrimGraph (g, 0, 30);
        if (iDir == 1)
          TrimGraph (g, 0, 10);
        myDrawSystFill (g, colorfulSystColors[nZdcCentBins-iCent], 0.6, 1001);
        SaferDelete (&g);

        g = make_graph (h);
        ResetXErrors (g);
        //if (strcmp (tag, "30GeVJets") == 0)
        //  TrimGraph (g, 0, 30);
        if (iDir == 1)
          TrimGraph (g, 0, 10);
        myDraw (g, colorfulColors[nZdcCentBins-iCent], kFullCircle, 1.0, 1, 2, "P", false);
        SaferDelete (&g);


        h = h_jetInt_trk_pt_iaaNoUnf[0][iPtJInt][iDir][iCent];
        g = make_graph (h);
        ResetXErrors (g);
        if (iDir == 1)
          TrimGraph (g, 0, 10);
        myDraw (g, colorfulColors[nZdcCentBins-iCent], kOpenCircle, 1.0, 1, 2, "P", false);
        SaferDelete (&g);


        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, colorfulColors[nZdcCentBins-iCent], Form ("#bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, colorfulColors[0], "#bf{All centralities}", 0.05);

      } // end loop over iCent

      c->cd ();
      myText (0.065, 0.971, kBlack, "#bf{#it{ATLAS}} Internal", 0.027);
  
      c->cd (1);
      myText (0.2, 0.80, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.05);
      myText (0.2, 0.74, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.05);
      c->cd (2);
      myText (0.2, 0.80, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV", iPtJInt == 0 ? 30 : 60), 0.05);
      myText (0.2, 0.74, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.05);
      c->cd (3);
      myLineText2 (0.26, 0.80, kBlack, kFullCircle, "Unfolded", 1.2, 0.05);
      myLineText2 (0.26, 0.74, kBlack, kOpenCircle, "No unfold", 1.2, 0.05);

      c->SaveAs (Form ("%s/Plots/PtCh/IpPb_Summary_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
    } // end loop over iDir
  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (int iDir : {0, 1, 2}) {

      const char* canvasName = Form ("c_jetInt_trk_pt_%s_sig2bkg_%s", directions[iDir].Data (), pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 800, 800);

      TH1D* h = nullptr; 
      TGAE* g = nullptr;

      c->cd (); 
      c->SetLogx ();
      c->SetLogy ();

      float ymin = 1e-1;
      float ymax = 1e6;

      c->Clear ();

      h = (TH1D*) h_jetInt_trk_pt_ref_sig[0][iPtJInt][iDir]->Clone ("h");
      h->Divide (h_jetInt_trk_pt_ref_bkg[0][iPtJInt][iDir]);
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
      g->SetMarkerStyle (kOpenCircle);
      g->SetMarkerSize (1.4);

      ((TGAE*) g->Clone ())->Draw ("AP");
      SaferDelete (&g);
      SaferDelete (&h);

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
        h = (TH1D*) h_jetInt_trk_pt_sig[0][iPtJInt][iDir][iCent]->Clone ("htemp");
        h->Divide (h_jetInt_trk_pt_bkg[0][iPtJInt][iDir][iCent]);
        g = make_graph (h);
        SaferDelete (&h);
        ResetXErrors (g);
        myDraw (g, colors[iCent], kOpenCircle, 1.4);
        SaferDelete (&g);
      }

      myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
      myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
      myText (0.22, 0.81, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
      myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} %s", iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.032);
      myLineText2 (0.25, 0.72, kBlack, kOpenCircle, "#bf{#it{pp}}", 1.4, 0.032, true);
      for (int iCent = 0; iCent < nZdcCentBins; iCent++)
        myLineText2 (0.25, 0.68-iCent*0.04, colors[iCent], kOpenCircle, Form ("#bf{#it{p}+Pb, %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 1.4, 0.032, true);

      c->SaveAs (Form ("%s/Plots/PtCh/SigToBkgd_%iGeVJets_%s_ptch.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));

    } // end loop over iDir

  } // end loop over iPtJInt




  for (int iDType : {0, 1}) {

    const int maxNSys = (iDType == 0 ? 6 : 10);

    for (int iDir : {0, 2}) {

      //if (makeTotalSystPlots) {
      //  for (short iPtJInt : {0, 1}) {

      //    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      //    const char* canvasName = Form ("c_jetInt_trk_pt_%s_pp_%s_syst_%s", directions[iDir].Data (), iDType == 0 ? "data" : "mc", pTJInt.Data ());
      //    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

      //    TH1D* h = nullptr;
      //    TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

      //    c->cd (); 
      //    c->SetLogx ();

      //    const float ymin = (iDType == 0 ? -maxDataSyst : -maxMCSyst);
      //    const float ymax = (iDType == 0 ?  maxDataSyst :  maxMCSyst);

      //    h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#delta N_{ch} / N_{ch} [%]", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
      //    h->GetXaxis ()->SetMoreLogLabels ();
      //    h->GetYaxis ()->SetRangeUser (ymin, ymax);

      //    h->SetLineWidth (1);
      //    h->SetLineStyle (2);
      //    h->DrawCopy ("hist ][");
      //    SaferDelete (&h);

      //    std::vector <int> totVars (0);
      //    if (iDType == 0) {
      //      totVars.push_back (0);
      //      totVars.push_back (1);
      //    }
      //    else 
      //      totVars.push_back (2);
      //    for (int iTotVar : totVars) {
      //      const TString totVar = totalVariations[iTotVar];
      //      g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_ref_systTot[iPtJInt][iDir][iTotVar], true, 100);
      //      gup = (TGAE*) g->Clone ();
      //      gdown = (TGAE*) g->Clone ();
      //      FlipTGAE (gdown);
      //   
      //      myDraw (gup, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
      //      myDraw (gdown, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
      //      myDrawFill (gup, gdown, varStyles[totVar].first, 0.1);

      //      SaferDelete (&g);
      //      SaferDelete (&gup);
      //      SaferDelete (&gdown);
      //    }

      //    for (int iVar = 1; iVar < nVar; iVar++) {
      //      const TString var = variations[iVar];
      //      if ((iDType == 0  && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) == 0))
      //        continue;

      //      g = (TGAE*) g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][iVar]->Clone ("gtemp");
      //      SaveRelativeErrors (g, g, g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][0], 100);
      //      ResetXErrors (g);
      //      ResetTGAEErrors (g);
      //      myDraw (g, varStyles[var].first, kDot, 0, varStyles[var].second, 4, "L");
      //      SaferDelete (&g);
      //    }

      //    myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
      //    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
      //    myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} %s", iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "< #pi/8" : "> 7#pi/8"), 0.032);

      //    int count = 0;
      //    for (int iVar = 1; iVar < nVar; iVar++) {
      //      const TString var = variations[iVar];
      //      if ((iDType == 0  && dataVariations.count (var) > 0) || (iDType == 1 && mcVariations.count (var) > 0)) {
      //        myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[var].first, varStyles[var].second, varFullNames[var], 0.5, 0.014);
      //        count++;
      //      }
      //    }
      //    for (int iTotVar : totVars) {
      //      const TString totVar = totalVariations[iTotVar];
      //      myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 0.5, 0.014);
      //      count++;
      //    }

      //    c->SaveAs (Form ("%s/Plots/Systematics/PtCh/TotalJetTaggedYield_pp_%s_ptch_%iGeVJets_%ssyst.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60, iDType == 0 ? "":"mcBased_"));

      //  } // end loop over iPtJInt


      //  for (short iPtJInt : {0, 1}) {

      //    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      //    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

      //      const char* canvasName = Form ("c_jetInt_trk_pt_%s_%s_pPb_iCent%i_syst_%s", directions[iDir].Data (), iDType == 0 ? "data" : "mc", iCent, pTJInt.Data ());
      //      TCanvas* c = new TCanvas (canvasName, "", 800, 800);

      //      TH1D* h = nullptr;
      //      TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

      //      c->cd (); 
      //      c->SetLogx ();

      //      const float ymin = (iDType == 0 ? -maxDataSyst : -maxMCSyst);
      //      const float ymax = (iDType == 0 ?  maxDataSyst :  maxMCSyst);

      //      h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#delta N_{ch} / N_{ch} [%]", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
      //      h->GetXaxis ()->SetMoreLogLabels ();
      //      h->GetYaxis ()->SetRangeUser (ymin, ymax);

      //      h->SetLineWidth (1);
      //      h->SetLineStyle (2);
      //      h->DrawCopy ("hist ][");
      //      SaferDelete (&h);

      //      std::vector <int> totVars (0);
      //      if (iDType == 0) {
      //        totVars.push_back (0);
      //        totVars.push_back (1);
      //      }
      //      else 
      //        totVars.push_back (2);
      //      for (int iTotVar : totVars) {
      //        const TString totVar = totalVariations[iTotVar];
      //        g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_systTot[iPtJInt][iDir][iCent][iTotVar], true, 100);
      //        gup = (TGAE*) g->Clone ();
      //        gdown = (TGAE*) g->Clone ();
      //        FlipTGAE (gdown);
      //     
      //        myDraw (gup, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
      //        myDraw (gdown, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
      //        myDrawFill (gup, gdown, varStyles[totVar].first, 0.1);

      //        SaferDelete (&g);
      //        SaferDelete (&gup);
      //        SaferDelete (&gdown);
      //      }

      //      for (int iVar = 1; iVar < nVar; iVar++) {
      //        const TString var = variations[iVar];
      //        if ((iDType == 0  && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) == 0))
      //          continue;

      //        g = (TGAE*) g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][iVar]->Clone ("gtemp");
      //        SaveRelativeErrors (g, g, g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][0], 100);
      //        ResetXErrors (g);
      //        ResetTGAEErrors (g);
      //        myDraw (g, varStyles[var].first, kDot, 0, varStyles[var].second, 4, "L");
      //        SaferDelete (&g);
      //      }

      //      myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
      //      myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{%s %i-%i%%}", iDType == 0 ? "ZDC" : "FCal", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
      //      myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} %s", iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "< #pi/8" : "> 7#pi/8"), 0.032);

      //      int count = 0;
      //      for (int iVar = 1; iVar < nVar; iVar++) {
      //        const TString var = variations[iVar];
      //        if ((iDType == 0  && dataVariations.count (var) > 0) || (iDType == 1 && mcVariations.count (var) > 0)) {
      //          myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[var].first, varStyles[var].second, varFullNames[var], 0.5, 0.014);
      //          count++;
      //        }
      //      }
      //      for (int iTotVar : totVars) {
      //        const TString totVar = totalVariations[iTotVar];
      //        myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 0.5, 0.014);
      //        count++;
      //      }
      //      
      //      c->SaveAs (Form ("%s/Plots/Systematics/PtCh/TotalJetTaggedYield_pPb_%i-%iperc_%s_ptch_%iGeVJets_%ssyst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60, iDType == 0 ? "":"mcBased_"));
      //    } // end loop over iCent

      //  } // end loop over iPtJInt
      //}



      //if (makeBkgdSystPlots && iDType != 1) {
      //  for (short iPtJInt : {0, 1}) {

      //    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      //    const char* canvasName = Form ("c_jetInt_trk_pt_%s_pp_%s_bkg_syst_%s", directions[iDir].Data (), iDType == 0 ? "data" : "mc", pTJInt.Data ());
      //    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

      //    TH1D* h = nullptr;
      //    TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

      //    c->cd (); 
      //    c->SetLogx ();

      //    const float ymin = (iDType == 0 ? -maxDataSyst : -maxMCSyst);
      //    const float ymax = (iDType == 0 ?  maxDataSyst :  maxMCSyst);

      //    h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#delta N_{ch} / N_{ch} [%]", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
      //    h->GetXaxis ()->SetMoreLogLabels ();
      //    h->GetYaxis ()->SetRangeUser (ymin, ymax);

      //    h->SetLineWidth (1);
      //    h->SetLineStyle (2);
      //    h->DrawCopy ("hist ][");
      //    SaferDelete (&h);

      //    std::vector <int> totVars (0);
      //    if (iDType == 0) {
      //      totVars.push_back (0);
      //      totVars.push_back (1);
      //    }
      //    else 
      //      totVars.push_back (2);
      //    for (int iTotVar : totVars) {
      //      const TString totVar = totalVariations[iTotVar];
      //      g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_ref_bkg_systTot[iPtJInt][iDir][iTotVar], true, 100);
      //      gup = (TGAE*) g->Clone ();
      //      gdown = (TGAE*) g->Clone ();
      //      FlipTGAE (gdown);
      //   
      //      myDraw (gup, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
      //      myDraw (gdown, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
      //      myDrawFill (gup, gdown, varStyles[totVar].first, 0.1);

      //      SaferDelete (&g);
      //      SaferDelete (&gup);
      //      SaferDelete (&gdown);
      //    }

      //    for (int iVar = 1; iVar < nVar; iVar++) {
      //      const TString var = variations[iVar];
      //      if ((iDType == 0  && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) == 0))
      //        continue;

      //      g = (TGAE*) g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][iVar]->Clone ("gtemp");
      //      SaveRelativeErrors (g, g, g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][0], 100);
      //      ResetXErrors (g);
      //      ResetTGAEErrors (g);
      //      myDraw (g, varStyles[var].first, kDot, 0, varStyles[var].second, 4, "L");
      //      SaferDelete (&g);
      //    }

      //    myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
      //    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
      //    myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} %s", iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "< #pi/8" : "> 7#pi/8"), 0.032);

      //    int count = 0;
      //    for (int iVar = 1; iVar < nVar; iVar++) {
      //      const TString var = variations[iVar];
      //      if ((iDType == 0  && dataVariations.count (var) > 0) || (iDType == 1 && mcVariations.count (var) > 0)) {
      //        myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[var].first, varStyles[var].second, varFullNames[var], 0.5, 0.014);
      //        count++;
      //      }
      //    }
      //    for (int iTotVar : totVars) {
      //      const TString totVar = totalVariations[iTotVar];
      //      myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 0.5, 0.014);
      //      count++;
      //    }

      //    c->SaveAs (Form ("%s/Plots/Systematics/PtCh/BkgdJetTaggedYield_pp_%s_ptch_%iGeVJets_%ssyst.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60, iDType == 0 ? "":"mcBased_"));

      //  } // end loop over iPtJInt


      //  for (short iPtJInt : {0, 1}) {

      //    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      //    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

      //      const char* canvasName = Form ("c_jetInt_trk_pt_%s_%s_pPb_iCent%i_bkgd_syst_%s", directions[iDir].Data (), iDType == 0 ? "data" : "mc", iCent, pTJInt.Data ());
      //      TCanvas* c = new TCanvas (canvasName, "", 800, 800);

      //      TH1D* h = nullptr;
      //      TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

      //      c->cd (); 
      //      c->SetLogx ();

      //      const float ymin = (iDType == 0 ? -maxDataSyst : -maxMCSyst);
      //      const float ymax = (iDType == 0 ?  maxDataSyst :  maxMCSyst);

      //      h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#delta N_{ch} / N_{ch} [%]", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
      //      h->GetXaxis ()->SetMoreLogLabels ();
      //      h->GetYaxis ()->SetRangeUser (ymin, ymax);

      //      h->SetLineWidth (1);
      //      h->SetLineStyle (2);
      //      h->DrawCopy ("hist ][");
      //      SaferDelete (&h);

      //      std::vector <int> totVars (0);
      //      if (iDType == 0) {
      //        totVars.push_back (0);
      //        totVars.push_back (1);
      //      }
      //      else 
      //        totVars.push_back (2);
      //      for (int iTotVar : totVars) {
      //        const TString totVar = totalVariations[iTotVar];
      //        g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_bkg_systTot[iDir][iPtJInt][iCent][iTotVar], true, 100);
      //        gup = (TGAE*) g->Clone ();
      //        gdown = (TGAE*) g->Clone ();
      //        FlipTGAE (gdown);
      //     
      //        myDraw (gup, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
      //        myDraw (gdown, varStyles[totVar].first, kDot, 0, varStyles[totVar].second, 2, "L");
      //        myDrawFill (gup, gdown, varStyles[totVar].first, 0.1);

      //        SaferDelete (&g);
      //        SaferDelete (&gup);
      //        SaferDelete (&gdown);
      //      }

      //      for (int iVar = 1; iVar < nVar; iVar++) {
      //        const TString var = variations[iVar];
      //        if ((iDType == 0  && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) == 0))
      //          continue;

      //        g = (TGAE*) g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][iVar]->Clone ("gtemp");
      //        SaveRelativeErrors (g, g, g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][0], 100);
      //        ResetXErrors (g);
      //        ResetTGAEErrors (g);
      //        myDraw (g, varStyles[var].first, kDot, 0, varStyles[var].second, 4, "L");
      //        SaferDelete (&g);
      //      }

      //      myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
      //      myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{%s %i-%i%%}", iDType == 0 ? "ZDC" : "FCal", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
      //      myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} %s", iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "< #pi/8" : "> 7#pi/8"), 0.032);

      //      int count = 0;
      //      for (int iVar = 1; iVar < nVar; iVar++) {
      //        const TString var = variations[iVar];
      //        if ((iDType == 0  && dataVariations.count (var) > 0) || (iDType == 1 && mcVariations.count (var) > 0)) {
      //          myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[var].first, varStyles[var].second, varFullNames[var], 0.5, 0.014);
      //          count++;
      //        }
      //      }
      //      for (int iTotVar : totVars) {
      //        const TString totVar = totalVariations[iTotVar];
      //        myLineColorText (0.35+0.35*(count/maxNSys), 0.41-(count%maxNSys)*0.020, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 0.5, 0.014);
      //        count++;
      //      }
      //      
      //      c->SaveAs (Form ("%s/Plots/Systematics/PtCh/BkgdJetTaggedYield_pPb_%i-%iperc_%s_ptch_%iGeVJets_%ssyst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60, iDType == 0 ? "":"mcBased_"));

      //    } // end loop over iCent

      //  } // end loop over iPtJInt
      //}



      if (makeSigSystPlots) {
        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          const char* canvasName = Form ("c_jetInt_trk_pt_%s_sig_%s_pp_syst_%s", directions[iDir].Data (), iDType == 0 ? "data" : "mc", pTJInt.Data ());
          TCanvas* c = new TCanvas (canvasName, "", 800, 800);

          TH1D* h = nullptr;
          TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

          c->cd (); 
          c->SetLogx ();

          const float ymin = (iDType == 0 ? -maxDataSyst : -maxMCSyst);
          const float ymax = (iDType == 0 ?  maxDataSyst :  maxMCSyst);

          h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#delta N_{ch} / N_{ch} [%]", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
          h->GetXaxis ()->SetMoreLogLabels ();
          h->GetYaxis ()->SetRangeUser (ymin, ymax);

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
            g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_ref_sig_systTot[iPtJInt][iDir][iTotVar], true, 100);
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

            g = (TGAE*) g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar]->Clone ("gtemp");
            SaveRelativeErrors (g, g, g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][0], 100);
            ResetXErrors (g);
            ResetTGAEErrors (g);
            myDraw (g, varStyles[var].first, kDot, 0, varStyles[var].second, 4, "L");
            SaferDelete (&g);
          }

          myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
          myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
          if (directions[iDir] == "ns")
            myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} < #pi/8 (near-side)", iPtJInt == 0 ? 30 : 60), 0.032);
          else if (directions[iDir] == "as")
            myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", iPtJInt == 0 ? 30 : 60), 0.032);

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

          c->SaveAs (Form ("%s/Plots/Systematics/PtCh/SignalJetTaggedYield_pp_%s_ptch_%iGeVJets_%ssyst.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60, iDType == 0 ? "":"mcBased_"));

        } // end loop over iPtJInt


        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

            const char* canvasName = Form ("c_jetInt_trk_pt_%s_sig_%s_pPb_iCent%i_syst_%s", directions[iDir].Data (), iDType == 0 ? "data" : "mc", iCent, pTJInt.Data ());
            TCanvas* c = new TCanvas (canvasName, "", 800, 800);

            TH1D* h = nullptr;
            TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

            c->cd (); 
            c->SetLogx ();
            //c->SetLogy ();

            const float ymin = (iDType == 0 ? -maxDataSyst : -maxMCSyst);
            const float ymax = (iDType == 0 ?  maxDataSyst :  maxMCSyst);

            h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#delta N_{ch} / N_{ch} [%]", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
            h->GetXaxis ()->SetMoreLogLabels ();
            h->GetYaxis ()->SetRangeUser (ymin, ymax);

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
              g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_sig_systTot[iPtJInt][iDir][iCent][iTotVar], true, 100);
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

              g = (TGAE*) g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar]->Clone ("gtemp");
              SaveRelativeErrors (g, g, g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][0], 100);
              ResetXErrors (g);
              ResetTGAEErrors (g);
              myDraw (g, varStyles[var].first, kDot, 0, varStyles[var].second, 4, "L");
              SaferDelete (&g);
            }

            myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
            myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{%s %i-%i%%}", iDType == 0 ? "ZDC" : "FCal", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
            if (directions[iDir] == "ns")
              myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} < #pi/8 (near-side)", iPtJInt == 0 ? 30 : 60), 0.032);
            else if (directions[iDir] == "as")
              myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", iPtJInt == 0 ? 30 : 60), 0.032);

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

            c->SaveAs (Form ("%s/Plots/Systematics/PtCh/SignalJetTaggedYield_pPb_%i-%iperc_%s_ptch_%iGeVJets_%ssyst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60, iDType == 0 ? "":"mcBased_"));
          } // end loop over iCent

        } // end loop over iPtJInt
      }



      if (makeIpPbSystPlots) {
        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
  
            const char* canvasName = Form ("c_jetInt_trk_pt_%s_iaa_%s_pPb_iCent%i_syst_%s", directions[iDir].Data (), iDType == 0 ? "data" : "mc", iCent, pTJInt.Data ());
            TCanvas* c = new TCanvas (canvasName, "", 800, 800);
  
            TH1D* h = nullptr;
            TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;
  
            c->cd (); 
            c->SetLogx ();
            //c->SetLogy ();
  
            const float ymin = (iDType == 0 ? -maxDataSyst : -maxMCSyst);
            const float ymax = (iDType == 0 ?  maxDataSyst :  maxMCSyst);
  
            h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#delta I_{pPb} / I_{pPb} [%]", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
            h->GetXaxis ()->SetMoreLogLabels ();
            h->GetYaxis ()->SetRangeUser (ymin, ymax);
  
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
              g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_iaa_systTot[iPtJInt][iDir][iCent][iTotVar], true, 100);
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
  
              g = (TGAE*) g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar]->Clone ("gtemp");
              SaveRelativeErrors (g, g, g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0], 100);
              ResetXErrors (g);
              ResetTGAEErrors (g);
              myDraw (g, varStyles[var].first, kDot, 0, varStyles[var].second, 4, "L");
              SaferDelete (&g);
            }
  
            myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.032);
            myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
            myText (0.22, 0.81, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{%s %i-%i%%}", iDType == 0 ? "ZDC" : "FCal", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
            if (directions[iDir] == "ns")
              myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} < #pi/8 (near-side)", iPtJInt == 0 ? 30 : 60), 0.032);
            else if (directions[iDir] == "as")
              myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", iPtJInt == 0 ? 30 : 60), 0.032);
  
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
  
            c->SaveAs (Form ("%s/Plots/Systematics/PtCh/JetTagged_IpPb_%i-%iperc_%s_ptch_%iGeVJets_%ssyst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60, iDType == 0 ? "":"mcBased_"));
          } // end loop over iCent

        } // end loop over iPtJInt

      }

    } // end loop over iDir

  } // end loop over iDType



  {
    const int maxNSys = 6;

    for (int iDir : {0, 2}) {

      if (makeSigSystPlots) {
        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          const char* canvasName = Form ("c_jetInt_trk_pt_%s_sig_data_pp_systTot_%s", directions[iDir].Data (), pTJInt.Data ());
          TCanvas* c = new TCanvas (canvasName, "", 800, 800);

          TH1D* h = nullptr;
          TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

          c->cd (); 
          c->SetLogx ();

          const float ymin = -maxDataSyst;
          const float ymax =  maxDataSyst;

          h = (TH1D*) h_jetInt_trk_pt_ref_sig[0][iPtJInt][iDir]->Clone ("h");
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
            g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_ref_sig_systTot[iPtJInt][iDir][iTotVar], true, 100);
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

          g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][0], true, 100);
          myDraw (g, kBlack, kDot, 0, 1, 2, "L");
          FlipTGAE (g);
          myDraw (g, kBlack, kDot, 0, 1, 2, "L");
          SaferDelete (&g);

          myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
          myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
          if (directions[iDir] == "ns")
            myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} < #pi/8 (near-side)", iPtJInt == 0 ? 30 : 60), 0.032);
          else if (directions[iDir] == "as")
            myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", iPtJInt == 0 ? 30 : 60), 0.032);

          int count = 0;
          for (int iTotVar : {0, 1, 2}) {
            const TString totVar = totalVariations[iTotVar];
            myLineColorText (0.35+0.35*(count/maxNSys), 0.35-(count%maxNSys)*0.040, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 1.0, 0.028);
            count++;
          }
          myLineColorText (0.35+0.35*(count/maxNSys), 0.35-(count%maxNSys)*0.040, kBlack, 1, "#bf{Total unc.}", 1.0, 0.028);

          c->SaveAs (Form ("%s/Plots/Systematics/PtCh/SignalJetTaggedYield_pp_%s_ptch_%iGeVJets_systTot.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60));

        } // end loop over iPtJInt


        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

            const char* canvasName = Form ("c_jetInt_trk_pt_%s_sig_data_pPb_iCent%i_systTot_%s", directions[iDir].Data (), iCent, pTJInt.Data ());
            TCanvas* c = new TCanvas (canvasName, "", 800, 800);

            TH1D* h = nullptr;
            TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

            c->cd (); 
            c->SetLogx ();
            //c->SetLogy ();

            const float ymin = -maxDataSyst;
            const float ymax =  maxDataSyst;

            h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#delta N_{ch} / N_{ch} [%]", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
            h->GetXaxis ()->SetMoreLogLabels ();
            h->GetYaxis ()->SetRangeUser (ymin, ymax);

            h->SetLineWidth (1);
            h->SetLineStyle (2);
            h->DrawCopy ("hist ][");
            SaferDelete (&h);

            for (int iTotVar : {0, 1, 2}) {
              const TString totVar = totalVariations[iTotVar];
              g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_sig_systTot[iPtJInt][iDir][iCent][iTotVar], true, 100);
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

            g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][0], true, 100);
            myDraw (g, kBlack, kDot, 0, 1, 2, "L");
            FlipTGAE (g);
            myDraw (g, kBlack, kDot, 0, 1, 2, "L");
            SaferDelete (&g);

            myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
            myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
            if (directions[iDir] == "ns")
              myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} < #pi/8 (near-side)", iPtJInt == 0 ? 30 : 60), 0.032);
            else if (directions[iDir] == "as")
              myText (0.22, 0.81, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", iPtJInt == 0 ? 30 : 60), 0.032);

            int count = 0;
            for (int iTotVar : {0, 1, 2}) {
              const TString totVar = totalVariations[iTotVar];
              myLineColorText (0.35+0.35*(count/maxNSys), 0.35-(count%maxNSys)*0.040, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 1.0, 0.028);
              count++;
            }
            myLineColorText (0.35+0.35*(count/maxNSys), 0.35-(count%maxNSys)*0.040, kBlack, 1, "#bf{Total unc.}", 1.0, 0.028);

            c->SaveAs (Form ("%s/Plots/Systematics/PtCh/SignalJetTaggedYield_pPb_%i-%iperc_%s_ptch_%iGeVJets_systTot.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60));
          } // end loop over iCent

        } // end loop over iPtJInt
      }



      if (makeIpPbSystPlots) {
        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

            const char* canvasName = Form ("c_jetInt_trk_pt_%s_iaa_data_pPb_iCent%i_systTot_%s", directions[iDir].Data (), iCent, pTJInt.Data ());
            TCanvas* c = new TCanvas (canvasName, "", 800, 800);

            TH1D* h = nullptr;
            TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

            c->cd (); 
            c->SetLogx ();
            //c->SetLogy ();

            const float ymin = -maxDataSyst;
            const float ymax =  maxDataSyst;

            h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#delta I_{pPb} / I_{pPb} [%]", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
            h->GetXaxis ()->SetMoreLogLabels ();
            h->GetYaxis ()->SetRangeUser (ymin, ymax);

            h->SetLineWidth (1);
            h->SetLineStyle (2);
            h->DrawCopy ("hist ][");
            SaferDelete (&h);

            for (int iTotVar : {0, 1, 2}) {
              const TString totVar = totalVariations[iTotVar];
              g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_iaa_systTot[iPtJInt][iDir][iCent][iTotVar], true, 100);
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

            g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0], true, 100);
            myDraw (g, kBlack, kDot, 0, 1, 2, "L");
            FlipTGAE (g);
            myDraw (g, kBlack, kDot, 0, 1, 2, "L");
            SaferDelete (&g);

            myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
            myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
            myText (0.22, 0.81, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}",  zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
            if (directions[iDir] == "ns")
              myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} < #pi/8 (near-side)", iPtJInt == 0 ? 30 : 60), 0.032);
            else if (directions[iDir] == "as")
              myText (0.22, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} > 7#pi/8 (away-side)", iPtJInt == 0 ? 30 : 60), 0.032);

            int count = 0;
            for (int iTotVar : {0, 1, 2}) {
              const TString totVar = totalVariations[iTotVar];
              myLineColorText (0.35+0.35*(count/maxNSys), 0.35-(count%maxNSys)*0.040, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 1.0, 0.028);
              count++;
            }
            myLineColorText (0.35+0.35*(count/maxNSys), 0.35-(count%maxNSys)*0.040, kBlack, 1, "#bf{Total unc.}", 1.0, 0.028);

            c->SaveAs (Form ("%s/Plots/Systematics/PtCh/JetTagged_IpPb_%i-%iperc_%s_ptch_%iGeVJets_systTot.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60));
          } // end loop over iCent

        } // end loop over iPtJInt

      }

    } // end loop over iDir

  }

}


#endif
