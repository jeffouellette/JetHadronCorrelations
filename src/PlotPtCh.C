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
#include <LinAlg.h>
#include <MyStyle.h>
#include <MyColors.h>

#include "CentralityDefs.h"
#include "Params.h"
#include "OutTree.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Variations.h"
#include "GetComparisons.h"


using namespace JetHadronCorrelations;


//const bool makeTotalSystPlots   = false;
//const bool makeBkgdSystPlots    = false;
const bool makeSigSystPlots     = true;
const bool makeIpPbSystPlots    = true;
const bool makeCovariancePlots  = false;
const bool DrawNoUnfold         = true;
const bool DrawConstFit         = false;


const float maxDataSyst = 20; // maximum y-axis for data-driven systematics
const float maxMCSyst = 20; // maximum y-axis for MC-driven systematics


//void DoOffset (TH1D* h, const char* tag, const short iCent) {
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
//void DoOffset (TGAE* g, const char* tag, const short iCent) {
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


double GetDiagonality (const TH2* h) {

  assert (h->GetNbinsX () == h->GetNbinsY ()); // h should be square

  const int d = h->GetNbinsX ();

  TMatrixD a = GetTMatrixD (h);

  TMatrixD j (d, 1);
  TMatrixD jt (1, d);
  TMatrixD r (d, 1);
  TMatrixD rt (1, d);
  TMatrixD r2 (d, 1);
  TMatrixD r2t (1, d);
  for (int i = 0; i < d; i++) {
    j[i][0] = 1; // all ones
    jt[0][i] = 1;
    r[i][0] = i+1;
    rt[0][i] = i+1;
    r2[i][0] = (i+1)*(i+1);
    r2t[0][i] = (i+1)*(i+1);
  }

  const double n = (jt * a * j)[0][0];
  const double sumx = (rt * a * j)[0][0];
  const double sumy = (jt * a * r)[0][0];
  const double sumx2 = (r2t * a * j)[0][0];
  const double sumy2 = (jt * a * r2)[0][0];
  const double sumxy = (rt * a * r)[0][0];

  return (n*sumxy - sumx*sumy) / (std::sqrt (n*sumx2 - sumx*sumx) * std::sqrt (n*sumy2 - sumy*sumy));
}


void PlotPtCh (const char* rawTag, const char* unfoldTag) {

  TLine* l = new TLine ();
  TLatex* tl = new TLatex ();

  TFile* inFile = nullptr;

  TH1D***   h_jet_pt_ref_unf          = Get2DArray <TH1D*> (2, 2);
  TH1D****  h_jet_pt_unf              = Get3DArray <TH1D*> (2, 2, nZdcCentBins+1);

  TH1D****  h_jetInt_trk_pt_ref       = Get3DArray <TH1D*> (2, 2, nDir);
  TH1D****  h_jetInt_trk_pt_ref_bkg   = Get3DArray <TH1D*> (2, 2, nDir);
  TH1D***** h_jetInt_trk_pt           = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);
  TH1D***** h_jetInt_trk_pt_bkg       = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);

  TH1D****  h_jetInt_trk_pt_ref_sig   = Get3DArray <TH1D*> (2, 2, nDir);
  TH1D***** h_jetInt_trk_pt_sig       = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);

  TH1D****  h_jetInt_trk_pt_ref_unf   = Get3DArray <TH1D*> (2, 2, nDir);
  TH1D***** h_jetInt_trk_pt_unf       = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);

  TH2D****    h2_jetInt_trk_pt_cov_ref_unf  = Get3DArray <TH2D*> (2, 2, nDir);
  TH2D*****   h2_jetInt_trk_pt_cov_unf      = Get4DArray <TH2D*> (2, 2, nDir, nZdcCentBins+1);

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


  TGAE** g_R_DpT_pPb_ATLAS = GetATLASJetFF ();
  TGAE** g_R_DpT_pPb_CMS = GetCMSJetFF ();
  TGAE** g_RpPb_CMS = GetCMSRpPb ();
  TGAE***** g_angantyr_iaa = GetPythiaAngantyrIpPb ();
  TH1D*** h_ampt_iaa = GetAMPTIpPb ();
  TGAE*** g_Zh_iaa = GetZHIAA ();


  TGAE**** g_jh_iaa_GAM = GetJhIAAGAM ();
  TGAE**** g_jh_pPb_unf_GAM = GetJhpPbYieldGAM ();
  TGAE***  g_jh_ref_unf_GAM = GetJhRefYieldGAM ();


  {
    TString inFileName = rawTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/ProcessCorrelations_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    for (short iDType = 0; iDType < 2; iDType++) {

      const TString dType = (iDType == 0 ? "data" : "mc");

      for (short iPtJInt : {0, 1}) {

        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          h_jetInt_trk_pt_ref[iDType][iPtJInt][iDir]      = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_%s_%s_Nominal",      dir.Data (), dType.Data (), pTJInt.Data ()));
          h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][iDir]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_bkg_%s_%s_Nominal",  dir.Data (), dType.Data (), pTJInt.Data ()));
          h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_sig_%s_%s_Nominal",  dir.Data (), dType.Data (), pTJInt.Data ()));

        } // end loop over iDir


        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

          for (short iDir = 0; iDir < nDir; iDir++) {

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

      for (short iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        for (short iVar = 0; iVar < nVar; iVar++) {

          const TString var = variations[iVar];

          //g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][iVar]     = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_ref_syst_%s_%s",      dir.Data (), pTJInt.Data (), var.Data ()));
          //g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][iVar] = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_ref_bkg_syst_%s_%s",  dir.Data (), pTJInt.Data (), var.Data ()));
          g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar] = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_ref_sig_syst_%s_%s",  dir.Data (), pTJInt.Data (), var.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            //g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][iVar]      = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_syst_%s_%s_%s",     dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            //g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][iVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_bkg_syst_%s_%s_%s", dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_sig_syst_%s_%s_%s", dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iVar

        for (short iTotVar = 0; iTotVar < 3; iTotVar++) {

          const TString totVar = totalVariations[iTotVar];

          //g_jetInt_trk_pt_ref_systTot[iPtJInt][iDir][iTotVar]     = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_ref_%s_systTot_%s",      dir.Data (), totVar.Data (), pTJInt.Data ()));
          //g_jetInt_trk_pt_ref_bkg_systTot[iPtJInt][iDir][iTotVar] = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_ref_bkg_%s_systTot_%s",  dir.Data (), totVar.Data (), pTJInt.Data ()));
          g_jetInt_trk_pt_ref_sig_systTot[iPtJInt][iDir][iTotVar] = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_ref_sig_%s_systTot_%s",  dir.Data (), totVar.Data (), pTJInt.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            //g_jetInt_trk_pt_systTot[iPtJInt][iDir][iCent][iTotVar]      = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_%s_systTot_%s_%s",     dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            //g_jetInt_trk_pt_bkg_systTot[iPtJInt][iDir][iCent][iTotVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_bkg_%s_systTot_%s_%s", dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            g_jetInt_trk_pt_sig_systTot[iPtJInt][iDir][iCent][iTotVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_sig_%s_systTot_%s_%s", dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));

          } // end loop over iCent

        } // end loop over iTotVar


        for (short iVar = 1; iVar < nVar; iVar++) {

          const TString var = variations[iVar];
          const TString dType = (dataVariations.count (var) > 0 ? "data" : "mc");

          //h_jetInt_trk_pt_ref_syst[iPtJInt][iDir][iVar]     = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_%s_%s_%s",     dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
          //h_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][iVar] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_bkg_%s_%s_%s", dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
          h_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_sig_%s_%s_%s", dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

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
    TString inFileName = unfoldTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/ProcessUnfolding_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    for (short iDType = 0; iDType < 2; iDType++) {

      const TString dType = (iDType == 0 ? "data" : "mc");

      for (short iPtJInt : {0, 1}) {

        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

        h_jet_pt_ref_unf[iDType][iPtJInt] = (TH1D*) inFile->Get (Form ("h_jet_pt_ref_unf_%s_%s_Nominal", dType.Data (), pTJInt.Data ()));

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir]      = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_unf_%s_%s_Nominal",  dir.Data (), dType.Data (), pTJInt.Data ()));

          h2_jetInt_trk_pt_cov_ref_unf[iDType][iPtJInt][iDir] = (TH2D*) inFile->Get (Form ("h2_jetInt_trk_pt_cov_%s_ref_unf_%s_%s",  dir.Data (), dType.Data (), pTJInt.Data ()));

        } // end loop over iDir


        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

          h_jet_pt_unf[iDType][iPtJInt][iCent] = (TH1D*) inFile->Get (Form ("h_jet_pt_unf_%s_%s_%s_Nominal", dType.Data (), cent.Data (), pTJInt.Data ()));

          for (short iDir = 0; iDir < nDir; iDir++) {

            const TString dir = directions[iDir];

            h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent]       = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_unf_%s_%s_%s_Nominal", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));
            h2_jetInt_trk_pt_cov_unf[iDType][iPtJInt][iDir][iCent]  = (TH2D*) inFile->Get (Form ("h2_jetInt_trk_pt_cov_%s_pPb_unf_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));

            h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent]       = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_iaa_%s_%s_%s_Nominal",     dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));

            h_jetInt_trk_pt_iaaNoUnf[iDType][iPtJInt][iDir][iCent]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_iaaNoUnf_%s_%s_%s_Nominal",  dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));

          } // end loop over iDir

        } // end loop over iCent

      } // end loop over iPtJInt

    } // end loop over iDType



    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      for (short iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        for (short iVar = 0; iVar < nVar; iVar++) {

          const TString var = variations[iVar];

          g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar] = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_ref_unf_syst_%s_%s",  dir.Data (), pTJInt.Data (), var.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_unf_syst_%s_%s_%s", dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_iaa_syst_%s_%s_%s", dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iVar

        for (short iTotVar = 0; iTotVar < 3; iTotVar++) {

          const TString totVar = totalVariations[iTotVar];

          g_jetInt_trk_pt_ref_unf_systTot[iPtJInt][iDir][iTotVar] = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_ref_unf_%s_systTot_%s",  dir.Data (), totVar.Data (), pTJInt.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            g_jetInt_trk_pt_unf_systTot[iPtJInt][iDir][iCent][iTotVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_unf_%s_systTot_%s_%s", dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            g_jetInt_trk_pt_iaa_systTot[iPtJInt][iDir][iCent][iTotVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_pt_%s_iaa_%s_systTot_%s_%s", dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));

          } // end loop over iCent

        } // end loop over iTotVar


        for (short iVar = 1; iVar < nVar; iVar++) {

          const TString var = variations[iVar];
          const TString dType = (dataVariations.count (var) > 0 ? "data" : "mc");

          h_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_unf_%s_%s_%s", dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            h_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_unf_%s_%s_%s_%s",  dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
            h_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_iaa_%s_%s_%s_%s",      dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iVar

      } // end loop over iDir

    } // end loop over iPtJInt

  }



  // print out mean jet pT values
  for (short iDType = 0; iDType < 2; iDType++) {

    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30 GeV" : "60 GeV");

      std::cout << "---------------" << std::endl << "JETS IN " << (iDType == 0 ? "DATA" : "MC") << " > " << pTJInt << std::endl << "---------------" << std::endl;

      TH1D* h = h_jet_pt_ref_unf[iDType][iPtJInt];
      h->GetXaxis ()->SetRange (3+2*iPtJInt, h->GetNbinsX () - 2);
      std::cout << "Average jet pT (pp): " << h->GetMean () << " +/- " << h->GetMeanError () << std::endl;
      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        h = h_jet_pt_unf[iDType][iPtJInt][iCent];
        h->GetXaxis ()->SetRange (3+2*iPtJInt, h->GetNbinsX () - 2);
        std::cout << "Average jet pT (p+Pb, " << zdcCentPercs[iCent+1] << "-" << zdcCentPercs[iCent] << "%): " << h->GetMean () << " +/- " << h->GetMeanError () << std::endl;
      } // end loop over iCent

      std::cout << "Formatted for latex:" << std::endl;
      h = h_jet_pt_ref_unf[iDType][iPtJInt];
      std::cout << Form ("$%.2f \\pm %.2f$", h->GetMean (), h->GetMeanError ());
      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        h = h_jet_pt_unf[iDType][iPtJInt][iCent];
        std::cout << " & " << Form ("$%.2f \\pm %.2f$", h->GetMean (), h->GetMeanError ());
      }
      std::cout << std::endl << std::endl;

    } // end loop over iPtJInt

  } // end loop over iDType



  l->SetLineStyle (7);
  l->SetLineWidth (1);
  l->SetLineColor (kBlack);



  for (short iDType = 0; iDType < 2; iDType++) {

    const TString dType = (iDType == 0 ? "data" : "mc");

    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      const bool hasRefBkg = (iDType != 1 && variationsWithNoppBkgd.count ("Nominal") == 0);
      const bool hasBkg = (variationsWithNopPbBkgd.count ("Nominal") == 0);

      for (short iCent = 0; iCent < nZdcCentBins; iCent++) {

        const char* canvasName = Form ("c_jetInt_trk_pt_%s_iCent%i_%s", dType.Data (), iCent, pTJInt.Data ());

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
        myDraw (g, myBlue, kFullCircle, 1.2);
        SaferDelete (&g);

        if (hasRefBkg) {
          //g = (TGAE*) g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][0][0]->Clone ();
          h = h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][0];
          //SetCentralValuesKeepRelativeErrors (g, h);
          //myDrawSyst (g, myPurple);
          //SaferDelete (&g);
          g = make_graph (h);
          ResetXErrors (g);
          myDraw (g, myPurple, kOpenCircle, 1.2);
          SaferDelete (&g);
        }

        //g = (TGAE*) g_jetInt_trk_pt_syst[iPtJInt][0][iCent][0]->Clone ();
        h = h_jetInt_trk_pt[iDType][iPtJInt][0][iCent];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myRed);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myRed, kFullCircle, 1.2);
        SaferDelete (&g);

        if (hasBkg) {
          //g = (TGAE*) g_jetInt_trk_pt_bkg_syst[iPtJInt][0][iCent][0]->Clone ();
          h = h_jetInt_trk_pt_bkg[iDType][iPtJInt][0][iCent];
          //SetCentralValuesKeepRelativeErrors (g, h);
          //myDrawSyst (g, myGreen);
          //SaferDelete (&g);
          g = make_graph (h);
          ResetXErrors (g);
          myDraw (g, myGreen, kOpenCircle, 1.2);
          SaferDelete (&g);
        }

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
        myDraw (g, myBlue, kFullCircle, 1.2);
        SaferDelete (&g);

        if (hasRefBkg) {
          //g = (TGAE*) g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][2][0]->Clone ();
          h = h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][2];
          //SetCentralValuesKeepRelativeErrors (g, h);
          //myDrawSyst (g, myPurple);
          //SaferDelete (&g);
          g = make_graph (h);
          ResetXErrors (g);
          myDraw (g, myPurple, kOpenCircle, 1.2);
          SaferDelete (&g);
        }

        //g = (TGAE*) g_jetInt_trk_pt_syst[iPtJInt][2][iCent][0]->Clone ();
        h = h_jetInt_trk_pt[iDType][iPtJInt][2][iCent];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myRed);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myRed, kFullCircle, 1.2);
        SaferDelete (&g);

        if (hasBkg) {
          //g = (TGAE*) g_jetInt_trk_pt_bkg_syst[iPtJInt][2][iCent][0]->Clone ();
          h = h_jetInt_trk_pt_bkg[iDType][iPtJInt][2][iCent];
          //SetCentralValuesKeepRelativeErrors (g, h);
          //myDrawSyst (g, myGreen);
          //SaferDelete (&g);
          g = make_graph (h);
          ResetXErrors (g);
          myDraw (g, myGreen, kOpenCircle, 1.2);
          SaferDelete (&g);
        }

        if (iDType == 0) myText (0.58, 0.78, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
        else {
          myText (0.42, 0.78, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.022/fuPad);
          myText (0.42, 0.72, kBlack, "Pythia8 (+ #it{p}+Pb Overlay)", 0.022/fuPad);
        }
        myBoxText2 (0.10, 0.24, myRed, kFullCircle, "#it{p}+Pb total", 1.2, 0.020/fuPad, true);
        myBoxText2 (0.10, 0.18, myBlue, kFullCircle, "#it{pp} total", 1.2, 0.020/fuPad, true);
        myBoxText2 (0.10, 0.12, myGreen, kOpenCircle, "#it{p}+Pb bkgd.", 1.2, 0.020/fuPad);
        myBoxText2 (0.10, 0.06, myPurple, kOpenCircle, "#it{pp} bkgd.", 1.2, 0.020/fuPad);

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
        myDraw (g, myBlue, kFullCircle, 1.2);
        SaferDelete (&g);

        //g = (TGAE*) g_jetInt_trk_pt_sig_syst[iPtJInt][0][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_sig[iDType][iPtJInt][0][iCent];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myRed);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myRed, kFullCircle, 1.2);
        SaferDelete (&g);


        ////g = (TGAE*) g_jetInt_trk_pt_ref_unf_syst[iPtJInt][0][0]->Clone ();
        //h = h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][0];
        ////SetCentralValuesKeepRelativeErrors (g, h);
        ////myDrawSyst (g, myBlue);
        ////SaferDelete (&g);
        //g = make_graph (h);
        //ResetXErrors (g);
        //myDraw (g, myBlue, kOpenCircle, 1.2);
        //SaferDelete (&g);

        ////g = (TGAE*) g_jetInt_trk_pt_unf_syst[iPtJInt][0][iCent][0]->Clone ();
        //h = h_jetInt_trk_pt_unf[iDType][iPtJInt][0][iCent];
        ////SetCentralValuesKeepRelativeErrors (g, h);
        ////myDrawSyst (g, myRed);
        ////SaferDelete (&g);
        //g = make_graph (h);
        //ResetXErrors (g);
        //myDraw (g, myRed, kOpenCircle, 1.2);
        //SaferDelete (&g);

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
        myDraw (g, myBlue, kFullCircle, 1.2);
        SaferDelete (&g);

        //g = (TGAE*) g_jetInt_trk_pt_sig_syst[iPtJInt][2][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_sig[iDType][iPtJInt][2][iCent];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myRed);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myRed, kFullCircle, 1.2);
        SaferDelete (&g);


        ////g = (TGAE*) g_jetInt_trk_pt_ref_unf_syst[iPtJInt][2][0]->Clone ();
        //h = h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][2];
        ////SetCentralValuesKeepRelativeErrors (g, h);
        ////myDrawSyst (g, myBlue);
        ////SaferDelete (&g);
        //g = make_graph (h);
        //ResetXErrors (g);
        //myDraw (g, myBlue, kOpenCircle, 1.2);
        //SaferDelete (&g);

        ////g = (TGAE*) g_jetInt_trk_pt_unf_syst[iPtJInt][2][iCent][0]->Clone ();
        //h = h_jetInt_trk_pt_unf[iDType][iPtJInt][2][iCent];
        ////SetCentralValuesKeepRelativeErrors (g, h);
        ////myDrawSyst (g, myRed);
        ////SaferDelete (&g);
        //g = make_graph (h);
        //ResetXErrors (g);
        //myDraw (g, myRed, kOpenCircle, 1.2);
        //SaferDelete (&g);


        dlPad->cd (); 
        dlPad->SetLogx ();

        ymin = 0.53;
        ymax = 1.47;

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
        h = h_jetInt_trk_pt_iaaNoUnf[iDType][iPtJInt][0][iCent];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myRed);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myRed, kFullCircle, 1.2);
        SaferDelete (&g);


        drPad->cd (); 
        drPad->SetLogx ();

        ymin = 0.53;
        ymax = 1.47;

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
        h = h_jetInt_trk_pt_iaaNoUnf[iDType][iPtJInt][2][iCent];
        //SetCentralValuesKeepRelativeErrors (g, h);
        //myDrawSyst (g, myRed);
        //SaferDelete (&g);
        g = make_graph (h);
        ResetXErrors (g);
        myDraw (g, myRed, kFullCircle, 1.2);
        SaferDelete (&g);

        c->SaveAs (Form ("%s/Plots/PtCh/JetTagged_HadronYields_%i-%iperc_comparison_PtCh_%iGeVJets%s.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], iPtJInt == 0 ? 30 : 60, iDType == 1 ? "_mc" : "")); 

      } // end loop over iCent

    } // end loop over iPtJInt

  } // end loop over iDType


/*
  for (short iCent = 0; iCent < nZdcCentBins; iCent++) {
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

    const short iVar = GetVarN ("FcalCentVar");

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
    ymax = 1.45;

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
    ymax = 1.45;

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



  for (short iCent = 0; iCent < nZdcCentBins; iCent++) {
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

    const short iVar = GetVarN ("NoFcalMixCatVar");

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
    ymax = 1.45;

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
    ymax = 1.45;

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




  if (makeCovariancePlots) {
    gStyle->SetPalette (kBlackBody);
    TColor::InvertPalette();

    const float zmin = 1e-4;
    const float zmax = 1;

    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      for (short iDir : {0, 1, 2}) {
    
        const char* canvasName = Form ("c_jetInt_trk_pt_cov_ref_%s_%s", directions[iDir].Data (), pTJInt.Data ());
        TCanvas* c = new TCanvas (canvasName, "", 960, 800);

        const double lMargin = 0.15;
        const double rMargin = 0.20;
        const double bMargin = 0.15;
        const double tMargin = 0.04;

        c->SetLeftMargin (lMargin);
        c->SetRightMargin (rMargin);
        c->SetBottomMargin (bMargin);
        c->SetTopMargin (tMargin);

        c->cd ();
    
        c->SetLogx ();
        c->SetLogy ();
        c->SetLogz ();

        TH2D* h2 = (TH2D*) h2_jetInt_trk_pt_cov_ref_unf[0][iPtJInt][iDir]->Clone ("htemp");
        //TH1D* h = h_jetInt_trk_pt_ref_unf[0][iPtJInt][iDir];
        for (int ix = 1; ix <= h2->GetNbinsX (); ix++) {
          for (int iy = 1; iy <= h2->GetNbinsY (); iy++) {
            h2->SetBinContent (ix, iy, h2->GetBinContent (ix, iy) / (std::sqrt (h2->GetBinContent (ix, ix)) * std::sqrt (h2->GetBinContent (iy, iy))));
            //if (h->GetBinContent (ix) * h->GetBinContent (iy) > 0)
            //  h2->SetBinContent (ix, iy, h2->GetBinContent (ix, iy) / (h->GetBinContent (ix) * h->GetBinContent (iy)));
            //const double c = h2->GetBinContent (ix, iy);
            //h2->SetBinContent (ix, iy, ((c >= 0) - (c < 0)) * std::sqrt (std::fabs (c)));
          }
        }
        const double r = GetDiagonality (h2);

        TAxis* xax = h2->GetXaxis ();
        TAxis* yax = h2->GetYaxis ();
        TAxis* zax = h2->GetZaxis ();

        xax->SetTitle ("#it{p}_{T}^{ch, 1} [GeV]");
        yax->SetTitle ("#it{p}_{T}^{ch, 2} [GeV]");
        //zax->SetTitle ("cov (dY_{1}, dY_{2}) / dY_{1} dY_{2}");
        //zax->SetTitle ("cov [dY(#it{p}_{T}^{ch, 1}, dY(#it{p}_{T}^{ch, 2})]");
        zax->SetTitle ("#rho [dY(#it{p}_{T}^{ch, 1}), dY(#it{p}_{T}^{ch, 2})]");

        xax->SetRangeUser (pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
        yax->SetRangeUser (pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
        zax->SetRangeUser (zmin, zmax);

        xax->SetMoreLogLabels ();
        yax->SetMoreLogLabels ();

        xax->SetTitleFont (43);
        xax->SetTitleSize (32);
        yax->SetTitleFont (43);
        yax->SetTitleSize (32);
        zax->SetTitleFont (43);
        zax->SetTitleSize (32);
        xax->SetLabelFont (43);
        xax->SetLabelSize (32);
        yax->SetLabelFont (43);
        yax->SetLabelSize (32);
        zax->SetLabelFont (43);
        zax->SetLabelSize (32);

        xax->SetTitleOffset (1.5);
        yax->SetTitleOffset (1.5);
        zax->SetTitleOffset (1.75);

        h2->DrawCopy ("colz");
        SaferDelete (&h2);

        TBox* b = new TBox (0.7, 33, iDir == 0 ? 13 : (iDir == 1 ? 26 : 15), 100);
        b->SetFillColorAlpha (kWhite, 1);
        b->SetLineColor (kBlack);
        b->SetLineWidth (1);
        b->SetLineStyle (1);
        b->Draw ("l");
        myText (0.22, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
        myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.036);
        myText (0.22, 0.80, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} %s", iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.036);
        myText (0.22, 0.75, kBlack, Form ("Diagonality = %.2f", r), 0.036);

        c->SaveAs (Form ("%s/Plots/PtCh/Covariance_SignalJetTaggedYield_pp_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
        SaferDelete (&b);
    
      } // end loop over iDir

    } // end loop over iPtJInt



    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

      for (short iPtJInt : {0, 1}) {

        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

        for (short iDir : {0, 1, 2}) {
      
          const char* canvasName = Form ("c_jetInt_trk_pt_cov_pPb_%s_%s_%s", directions[iDir].Data (), cent.Data (), pTJInt.Data ());
          TCanvas* c = new TCanvas (canvasName, "", 960, 800);

          const double lMargin = 0.15;
          const double rMargin = 0.20;
          const double bMargin = 0.15;
          const double tMargin = 0.04;

          c->SetLeftMargin (lMargin);
          c->SetRightMargin (rMargin);
          c->SetBottomMargin (bMargin);
          c->SetTopMargin (tMargin);

          c->cd ();
      
          c->SetLogx ();
          c->SetLogy ();
          c->SetLogz ();

          TH2D* h2 = (TH2D*) h2_jetInt_trk_pt_cov_unf[0][iPtJInt][iDir][iCent]->Clone ("htemp");
          //TH1D* h = h_jetInt_trk_pt_unf[0][iPtJInt][iDir][iCent];
          for (int ix = 1; ix <= h2->GetNbinsX (); ix++) {
            for (int iy = 1; iy <= h2->GetNbinsY (); iy++) {
              h2->SetBinContent (ix, iy, h2->GetBinContent (ix, iy) / (std::sqrt (h2->GetBinContent (ix, ix)) * std::sqrt (h2->GetBinContent (iy, iy))));
              //if (h->GetBinContent (ix) * h->GetBinContent (iy) > 0)
              //  h2->SetBinContent (ix, iy, h2->GetBinContent (ix, iy) / (h->GetBinContent (ix) * h->GetBinContent (iy)));
              //const double c = h2->GetBinContent (ix, iy);
              //h2->SetBinContent (ix, iy, ((c >= 0) - (c < 0)) * std::sqrt (std::fabs (c)));
            }
          }
          const double r = GetDiagonality (h2);

          TAxis* xax = h2->GetXaxis ();
          TAxis* yax = h2->GetYaxis ();
          TAxis* zax = h2->GetZaxis ();

          xax->SetTitle ("#it{p}_{T}^{ch, 1} [GeV]");
          yax->SetTitle ("#it{p}_{T}^{ch, 2} [GeV]");
          //zax->SetTitle ("cov (dY_{1}, dY_{2}) / dY_{1} dY_{2}");
          zax->SetTitle ("#rho [dY(#it{p}_{T}^{ch, 1}, dY(#it{p}_{T}^{ch, 2})]");

          xax->SetRangeUser (pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
          yax->SetRangeUser (pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
          zax->SetRangeUser (zmin, zmax);

          xax->SetMoreLogLabels ();
          yax->SetMoreLogLabels ();

          xax->SetTitleFont (43);
          xax->SetTitleSize (32);
          yax->SetTitleFont (43);
          yax->SetTitleSize (32);
          zax->SetTitleFont (43);
          zax->SetTitleSize (32);
          xax->SetLabelFont (43);
          xax->SetLabelSize (32);
          yax->SetLabelFont (43);
          yax->SetLabelSize (32);
          zax->SetLabelFont (43);
          zax->SetLabelSize (32);

          xax->SetTitleOffset (1.5);
          yax->SetTitleOffset (1.5);
          zax->SetTitleOffset (1.75);

          h2->DrawCopy ("colz");
          SaferDelete (&h2);

          TBox* b = new TBox (0.7, 24, iDir == 0 ? 13 : (iDir == 1 ? 26 : 15), 100);
          b->SetFillColorAlpha (kWhite, 1);
          b->SetLineColor (kBlack);
          b->SetLineWidth (1);
          b->SetLineStyle (1);
          b->Draw ("l");
          myText (0.22, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
          myText (0.22, 0.85, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.036);
          myText (0.22, 0.80, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} %s", iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.036);
          if (iCent == nZdcCentBins)
            myText (0.22, 0.75, kBlack, "ZDC 0-100%", 0.036);
          else
            myText (0.22, 0.75, kBlack, Form ("ZDC %i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.036);
          myText (0.22, 0.70, kBlack, Form ("Diagonality = %.2f", r), 0.036);

          c->SaveAs (Form ("%s/Plots/PtCh/Covariance_SignalJetTaggedYield_pPb_%s_%iGeVJets_%s.pdf", workPath.Data (), cent.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
          SaferDelete (&b);
      
        } // end loop over iDir

      } // end loop over iPtJInt

    } // end loop over iCent
  }




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (short iDir : {0, 1, 2}) {
  
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
      //h = new TH1D ("htemp", ";#it{p}_{T}^{ch} [GeV];(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]", 1, 0, pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
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
      myDrawSystFill (g, colorfulSystColors[0], 0.8, 1001);
      SaferDelete (&g);

      g = make_graph (h);
      ScaleGraph (g, nullptr, std::pow (10, 3));
      ResetXErrors (g);
      myDraw (g, colorfulColors[0], kFullCircle, 1.4, 1, 2, "P", false);
      SaferDelete (&g);


      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        h = (TH1D*) h_jetInt_trk_pt_ref_unf[0][iPtJInt][iDir]->Clone ("htemp");
        h->Scale (std::pow (10, 2-iCent));
        //myDrawHist (h, kBlack, 1, 3);
        myDraw (h, kBlack, kDot, 0, 1, 3);
        SaferDelete (&h);

        g = (TGAE*) g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_unf[0][iPtJInt][iDir][iCent];
        SetCentralValuesKeepRelativeErrors (g, h);
        ScaleGraph (g, nullptr, std::pow (10, 2-iCent));
        myDrawSystFill (g, colorfulSystColors[iCent+1], 0.8, 1001);
        SaferDelete (&g);
  
        g = make_graph (h);
        ResetXErrors (g);
        ScaleGraph (g, nullptr, std::pow (10, 2-iCent));
        myDraw (g, colorfulColors[iCent+1], kFullCircle, 1.4, 1, 2, "P", false);
        SaferDelete (&g);

      } // end loop over iCent

      myText (0.56, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
      myText (0.56, 0.86, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} %s", iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.032);

      myText (0.18, 0.385, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.030);
      myText (0.18, 0.350, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.030);
      myText (0.18, 0.315, kBlack, "ZDC centralities", 0.030);

      mySimpleMarkerAndBoxAndLineText (0.25, 0.275, 1.4, 1001, colorfulSystColors[0], 0.8, colorfulColors[0], kFullCircle, 1.6, "#it{pp} (#times10^{3})", 0.030);
      for (short iCent = 0; iCent < nZdcCentBins; iCent++) { 
        mySimpleMarkerAndBoxAndLineText (0.25 + (iCent >= 2 ? 0.35 : 0), 0.275-((iCent+1)%3)*0.035, 1.4, 1001, colorfulSystColors[iCent+1], 0.8, colorfulColors[iCent+1], kFullCircle, 1.6, Form ("%i-%i%% (#times10^{%i})", zdcCentPercs[iCent+1], zdcCentPercs[iCent], 2-iCent), 0.030);
      }
      mySimpleMarkerAndBoxAndLineText (0.60, 0.17, 1.4, 1001, colorfulSystColors[nZdcCentBins+1], 0.8, colorfulColors[nZdcCentBins+1], kFullCircle, 1.6, Form ("0-100%% (#times10^{%i})", 2-nZdcCentBins), 0.030);
      mySimpleMarkerAndBoxAndLineText (0.25, 0.17, 1.4, 1001, kWhite, 0.0, kBlack, kDot, 0.0, "#it{pp} (#it{scaled})", 0.030);

      c->RedrawAxis();

      c->SaveAs (Form ("%s/Plots/PtCh/SignalJetTaggedYield_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
  
    } // end loop over iDir

  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (short iDir : {0, 1, 2}) {
      const char* canvasName = Form ("c_jetInt_trk_pt_IpPb_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1200, 800);
      c->Divide (3, 2);

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
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
        if (iDir == 1)
          TrimGraph (g, 0, 10);
        RecenterGraph (g);
        //TGAE* gstat = (TGAE*) g->Clone ("gstat");
        myDrawSystFill (g, colorfulSystColors[nZdcCentBins-iCent], 0.6, 1001);
        SaferDelete (&g);

        g = make_graph (h);
        if (iDir == 1)
          TrimGraph (g, 0, 10);
        RecenterGraph (g);
        ResetXErrors (g);
        myDraw (g, colorfulColors[nZdcCentBins-iCent], kFullCircle, 1.0, 1, 2, "P", false);
        SaferDelete (&g);

        TF1* cfit = nullptr;
        if (DrawConstFit) {
          TF1* cfit = new TF1 ("cfit", "[0]", 6, iPtJInt == 0 ? 30 : 60);
          cfit->SetParameter (0, 1);
          h_jetInt_trk_pt_iaa[0][iPtJInt][iDir][iCent]->Fit (cfit, "RN0Q");
          cfit->SetLineColor (colorfulColors[nZdcCentBins-iCent]);
          cfit->SetLineStyle (2);
          cfit->SetLineWidth (2);
          ((TF1*) cfit->Clone ())->Draw ("same");
          myText (0.56, 0.22, colorfulColors[nZdcCentBins-iCent], Form ("#chi^{2}/dof = %.2f / %i", cfit->GetChisquare (), cfit->GetNDF ()), 0.05); 
          //SaferDelete (&cfit);
        }
        //SaferDelete (&gstat);


        if (DrawNoUnfold) {
          h = h_jetInt_trk_pt_iaaNoUnf[0][iPtJInt][iDir][iCent];
          g = make_graph (h);
          ResetXErrors (g);
          if (iDir == 1)
            TrimGraph (g, 0, 10);
          RecenterGraph (g);
          myDraw (g, colorfulColors[nZdcCentBins-iCent], kOpenCircle, 1.0, 1, 2, "P", false);
          SaferDelete (&g);
        }


        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, colorfulColors[nZdcCentBins-iCent], Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, colorfulColors[0], "#bf{#it{p}+Pb, 0-100%}", 0.05);

      } // end loop over iCent

      c->cd ();
      myText (0.065, 0.971, kBlack, "#bf{#it{ATLAS}} Internal", 0.027);
  
      c->cd (1);
      myText (0.2, 0.80, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.05);
      myText (0.2, 0.74, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.05);
      c->cd (2);
      myText (0.2, 0.80, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV", iPtJInt == 0 ? 30 : 60), 0.05);
      myText (0.2, 0.74, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.05);

      if (DrawNoUnfold) {
        c->cd (3);
        myLineText2 (0.26, 0.80, kBlack, kFullCircle, "Unfolded", 1.2, 0.05);
        myLineText2 (0.26, 0.74, kBlack, kOpenCircle, "No unfold", 1.2, 0.05);
      }

      c->SaveAs (Form ("%s/Plots/PtCh/IpPb_Summary_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
    } // end loop over iDir
  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {
 
    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    //for (short iDir : {0, 1, 2}) {
    for (short iDir : {0, 2}) {

      const char* canvasName = Form ("c_jetInt_trk_pt_yield_GAM_comp_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
      c->Divide (4, 2);

      TGAE* g = nullptr;
      TH1D* h = nullptr;

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kBlack);

      {
        c->cd (7);

        gPad->SetLogx();

        TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];Y(#it{p}_{T}^{ch}) Now / at GAM", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser ((iPtJInt == 0 ? 0.61 : 0.74), (iPtJInt == 0 ? 1.6 : 1.45));
        //h->GetYaxis ()->SetRangeUser (0.0, 3);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        TGAE* g = (TGAE*) g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][0]->Clone ();
        h = h_jetInt_trk_pt_ref_unf[0][iPtJInt][iDir];
        SetCentralValuesKeepRelativeErrors (g, h);
        ScaleByGraph (g, g_jh_ref_unf_GAM[iPtJInt][iDir]);
        if (iDir == 1)
          TrimGraph (g, 0, 10);
        RecenterGraph (g);
        myDrawSystFill (g, colorfulSystColors[nZdcCentBins+1], 0.6, 1001);
        SaferDelete (&g);

        g = make_graph (h);
        ScaleByGraph (g, g_jh_ref_unf_GAM[iPtJInt][iDir]);
        if (iDir == 1)
          TrimGraph (g, 0, 10);
        RecenterGraph (g);
        ResetXErrors (g);
        myDraw (g, colorfulColors[nZdcCentBins+1], kFullCircle, 1.0, 1, 2, "P", false);
        SaferDelete (&g);

        myText (0.2, 0.865, colorfulColors[nZdcCentBins+1], "#bf{#it{pp}}", 0.05);
      }


      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        gPad->SetLogx();

        TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];Y(#it{p}_{T}^{ch}) Now / at GAM", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser ((iPtJInt == 0 ? 0.61 : 0.74), (iPtJInt == 0 ? 1.6 : 1.45));
        //h->GetYaxis ()->SetRangeUser (0.0, 3);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        TGAE* g = (TGAE*) g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_unf[0][iPtJInt][iDir][iCent];
        SetCentralValuesKeepRelativeErrors (g, h);
        ScaleByGraph (g, g_jh_pPb_unf_GAM[iPtJInt][iDir][iCent]);
        if (iDir == 1)
          TrimGraph (g, 0, 10);
        RecenterGraph (g);
        myDrawSystFill (g, colorfulSystColors[nZdcCentBins-iCent], 0.6, 1001);
        SaferDelete (&g);

        g = make_graph (h);
        ScaleByGraph (g, g_jh_pPb_unf_GAM[iPtJInt][iDir][iCent]);
        if (iDir == 1)
          TrimGraph (g, 0, 10);
        RecenterGraph (g);
        ResetXErrors (g);
        myDraw (g, colorfulColors[nZdcCentBins-iCent], kFullCircle, 1.0, 1, 2, "P", false);
        SaferDelete (&g);

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, colorfulColors[nZdcCentBins-iCent], Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, colorfulColors[0], "#bf{#it{p}+Pb, 0-100%}", 0.05);

      } // end loop over iCent

      c->cd (8);
      myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
      myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV", iPtJInt == 0 ? 30 : 60), 0.07);
      myText (0.1, 0.57, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.07);

      c->SaveAs (Form ("%s/Plots/PtCh/Yield_GAM_Comparison_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));

    } // end loop over iDir

  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (short iDir : {0, 2}) {
      const char* canvasName = Form ("c_jetInt_trk_pt_IpPb_GAM_Ratio_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1200, 800);
      c->Divide (3, 2);

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        gPad->SetLogx ();

        TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#it{I}_{#it{p}Pb} Now / at GAM", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
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
        ScaleByGraph (g, g_jh_iaa_GAM[iPtJInt][iDir][iCent]);
        if (iDir == 1)
          TrimGraph (g, 0, 10);
        RecenterGraph (g);
        //TGAE* gstat = (TGAE*) g->Clone ("gstat");
        myDrawSystFill (g, colorfulSystColors[nZdcCentBins-iCent], 0.6, 1001);
        SaferDelete (&g);

        g = make_graph (h);
        ScaleByGraph (g, g_jh_iaa_GAM[iPtJInt][iDir][iCent]);
        if (iDir == 1)
          TrimGraph (g, 0, 10);
        RecenterGraph (g);
        ResetXErrors (g);
        myDraw (g, colorfulColors[nZdcCentBins-iCent], kFullCircle, 1.0, 1, 2, "P", false);
        SaferDelete (&g);

        //g = g_jh_iaa_GAM[iPtJInt][iDir][iCent];
        //ResetXErrors (g);
        //if (iDir == 1)
        //  TrimGraph (g, 0, 10);
        //RecenterGraph (g);
        //myDraw (g, colorfulColors[nZdcCentBins-iCent], kOpenCircle, 1.0, 1, 2, "P", false);
        //SaferDelete (&g);

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, colorfulColors[nZdcCentBins-iCent], Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, colorfulColors[0], "#bf{#it{p}+Pb, 0-100%}", 0.05);

      } // end loop over iCent

      c->cd ();
      myText (0.065, 0.971, kBlack, "#bf{#it{ATLAS}} Internal", 0.027);
  
      c->cd (1);
      myText (0.2, 0.80, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.05);
      myText (0.2, 0.74, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.05);
      c->cd (2);
      myText (0.2, 0.80, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV", iPtJInt == 0 ? 30 : 60), 0.05);
      myText (0.2, 0.74, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.05);

      //if (DrawNoUnfold) {
      //  c->cd (3);
      //  myLineText2 (0.26, 0.80, kBlack, kFullCircle, "Current", 1.2, 0.05);
      //  myLineText2 (0.26, 0.74, kBlack, kOpenCircle, "GAM", 1.2, 0.05);
      //}

      c->SaveAs (Form ("%s/Plots/PtCh/IpPb_GAM_Ratio_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
    } // end loop over iDir
  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (short iDir : {0, 2}) {
      const char* canvasName = Form ("c_jetInt_trk_pt_IpPb_GAM_Comparison_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1200, 800);
      c->Divide (3, 2);

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
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
        if (iDir == 1)
          TrimGraph (g, 0, 10);
        RecenterGraph (g);
        //TGAE* gstat = (TGAE*) g->Clone ("gstat");
        myDrawSystFill (g, colorfulSystColors[nZdcCentBins-iCent], 0.6, 1001);
        SaferDelete (&g);

        g = make_graph (h);
        if (iDir == 1)
          TrimGraph (g, 0, 10);
        RecenterGraph (g);
        ResetXErrors (g);
        myDraw (g, colorfulColors[nZdcCentBins-iCent], kFullCircle, 1.0, 1, 2, "P", false);
        SaferDelete (&g);

        g = g_jh_iaa_GAM[iPtJInt][iDir][iCent];
        ResetXErrors (g);
        if (iDir == 1)
          TrimGraph (g, 0, 10);
        RecenterGraph (g);
        myDraw (g, colorfulColors[nZdcCentBins-iCent], kOpenCircle, 1.0, 1, 2, "P", false);
        SaferDelete (&g);

        h = h_jetInt_trk_pt_iaaNoUnf[0][iPtJInt][iDir][iCent];
        g = make_graph (h);
        ResetXErrors (g);
        if (iDir == 1)
          TrimGraph (g, 0, 10);
        RecenterGraph (g);
        myDraw (g, colorfulColors[nZdcCentBins-iCent], kOpenTriangleUp, 1.0, 1, 2, "P", false);
        SaferDelete (&g);

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, colorfulColors[nZdcCentBins-iCent], Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, colorfulColors[0], "#bf{#it{p}+Pb, 0-100%}", 0.05);

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
      myLineText2 (0.26, 0.80, kBlack, kFullCircle, "Current", 1.2, 0.05);
      myLineText2 (0.26, 0.74, kBlack, kOpenCircle, "GAM", 1.2, 0.05);
      myLineText2 (0.26, 0.68, kBlack, kOpenTriangleUp, "Current, no unfold", 1.2, 0.05);

      c->SaveAs (Form ("%s/Plots/PtCh/IpPb_GAM_Comparison_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
    } // end loop over iDir
  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {
  //{
    //const short iPtJInt = 1;

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (short iDir : {0, 1, 2}) {
      const char* canvasName = Form ("c_jetInt_trk_pt_IpPb_fewCent_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1200, 400);
      c->Divide (3, 1);

      int iCanvas = 1;
      for (short iCent : {0, 2, 4}) {
        c->cd (iCanvas++);

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
        if (iDir == 1)
          TrimGraph (g, 0, 10);
        RecenterGraph (g);
        myDrawSystFill (g, colorfulSystColors[nZdcCentBins-iCent], 0.6, 1001);
        SaferDelete (&g);

        g = make_graph (h);
        if (iDir == 1)
          TrimGraph (g, 0, 10);
        RecenterGraph (g);
        ResetXErrors (g);
        myDraw (g, colorfulColors[nZdcCentBins-iCent], kFullCircle, 1.0, 1, 2, "P", false);
        SaferDelete (&g);


        //h = h_jetInt_trk_pt_iaaNoUnf[0][iPtJInt][iDir][iCent];
        //g = make_graph (h);
        //ResetXErrors (g);
        //if (iDir == 1)
        //  TrimGraph (g, 0, 10);
        //RecenterGraph (g);
        //myDraw (g, colorfulColors[nZdcCentBins-iCent], kOpenCircle, 1.0, 1, 2, "P", false);
        //SaferDelete (&g);


        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, colorfulColors[nZdcCentBins-iCent], Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, colorfulColors[0], "#bf{#it{p}+Pb, 0-100%}", 0.05);

      } // end loop over iCent

      c->cd ();
      myText (0.065, 0.950, kBlack, "#bf{#it{ATLAS}} Internal", 0.054);
  
      c->cd (1);
      myText (0.2, 0.80, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.05);
      myText (0.2, 0.74, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.05);
      c->cd (2);
      myText (0.2, 0.80, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV", iPtJInt == 0 ? 30 : 60), 0.05);
      myText (0.2, 0.74, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.05);
      //c->cd (3);
      //myLineText2 (0.26, 0.80, kBlack, kFullCircle, "Unfolded", 1.2, 0.05);
      //myLineText2 (0.26, 0.74, kBlack, kOpenCircle, "No unfold", 1.2, 0.05);

      c->SaveAs (Form ("%s/Plots/PtCh/IpPb_FewCent_Summary_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
    } // end loop over iDir
  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    {
      const short iDir = 0; // FF should only be compared to near-side.
      const char* canvasName = Form ("c_jetInt_trk_pt_IpPb_ATLAS_FF_Comp_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1200, 800);
      c->Divide (3, 2);

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        gPad->SetLogx ();

        TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#it{I}_{#it{p}Pb}", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser ((iPtJInt == 0 ? 0.61 : 0.74), (iPtJInt == 0 ? 1.6 : 1.45));
        //h->GetYaxis ()->SetRangeUser (0.0, 3);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        TGAE* g = nullptr;

        if (iPtJInt == 0) // draw 45-60GeV plot
          g = (TGAE*) g_R_DpT_pPb_ATLAS[0]->Clone ("gtemp");
        else // draw 60-80GeV plot
          g = (TGAE*) g_R_DpT_pPb_ATLAS[1]->Clone ("gtemp");
        myDrawSyst (g, kBlack);
        RecenterGraph (g);
        ResetTGAEErrors (g);
        ResetXErrors (g);
        myDraw (g, kBlack, kFullSquare, 1.4, 1, 2, "P", false);
        SaferDelete (&g);

        g = (TGAE*) g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_iaa[0][iPtJInt][iDir][iCent];
        //h = h_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iMCTruthJetsTruthParts];
        SetCentralValuesKeepRelativeErrors (g, h);
        RecenterGraph (g);
        myDrawSystFill (g, colorfulSystColors[nZdcCentBins-iCent], 0.6, 1001);
        SaferDelete (&g);

        g = make_graph (h);
        RecenterGraph (g);
        ResetXErrors (g);
        myDraw (g, colorfulColors[nZdcCentBins-iCent], kFullCircle, 1.0, 1, 2, "P", false);
        SaferDelete (&g);

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, colorfulColors[nZdcCentBins-iCent], Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, colorfulColors[0], "#bf{#it{p}+Pb, 0-100%}", 0.05);

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
      myLineText2 (0.26, 0.80, kBlack, kFullCircle, "This result", 1.2, 0.05);
      myLineText2 (0.26, 0.73, kBlack, kFullSquare, Form ("ATLAS #it{R}_{D(#it{p}_{T})}; #it{p}_{T}^{jet} = %i-%i GeV", iPtJInt == 0 ? 45 : 60, iPtJInt == 0 ? 60 : 80), 1.2, 0.05);

      c->SaveAs (Form ("%s/Plots/PtCh/IpPb_ATLAS_FF_Comp_Summary_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
    } // end iDir=0 scope
  } // end loop over iPtJInt



  { 

    const short iPtJInt = 1;
    const TString pTJInt = "60GeV";

    {
      const short iDir = 0; // FF should only be compared to near-side.
      const char* canvasName = Form ("c_jetInt_trk_pt_IpPb_CMS_FF_Comp_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1200, 800);
      c->Divide (3, 2);

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        gPad->SetLogx ();

        TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#it{I}_{#it{p}Pb}", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser ((iPtJInt == 0 ? 0.61 : 0.74), (iPtJInt == 0 ? 1.6 : 1.45));
        //h->GetYaxis ()->SetRangeUser (0.0, 3);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        TGAE* g = nullptr;

        g = (TGAE*) g_R_DpT_pPb_CMS[0]->Clone ("gtemp");
        myDrawSyst (g, kBlack);
        SaferDelete (&g);

        g = (TGAE*) g_R_DpT_pPb_CMS[1]->Clone ("gtemp");
        RecenterGraph (g);
        ResetXErrors (g);
        myDraw (g, kBlack, kFullSquare, 1.4, 1, 2, "P", false);
        SaferDelete (&g);

        g = (TGAE*) g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_iaa[0][iPtJInt][iDir][iCent];
        //h = h_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iMCTruthJetsTruthParts];
        SetCentralValuesKeepRelativeErrors (g, h);
        RecenterGraph (g);
        myDrawSystFill (g, colorfulSystColors[nZdcCentBins-iCent], 0.6, 1001);
        SaferDelete (&g);

        g = make_graph (h);
        RecenterGraph (g);
        ResetXErrors (g);
        myDraw (g, colorfulColors[nZdcCentBins-iCent], kFullCircle, 1.0, 1, 2, "P", false);
        SaferDelete (&g);

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, colorfulColors[nZdcCentBins-iCent], Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, colorfulColors[0], "#bf{#it{p}+Pb, 0-100%}", 0.05);

      } // end loop over iCent

      c->cd ();
      myText (0.065, 0.971, kBlack, "#bf{#it{ATLAS}} Internal", 0.027);
  
      c->cd (1);
      myText (0.2, 0.80, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.05);
      myText (0.2, 0.74, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.05);
      c->cd (2);
      myText (0.2, 0.80, kBlack, "#it{p}_{T}^{jet} > 60 GeV", 0.05);
      myText (0.2, 0.74, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.05);
      c->cd (3);
      myLineText2 (0.26, 0.80, kBlack, kFullCircle, "This result", 1.2, 0.05);
      myLineText2 (0.26, 0.73, kBlack, kFullSquare, "CMS #it{R}_{D(#it{p}_{T})}; #it{p}_{T}^{jet} = 60-80 GeV", 1.2, 0.05);

      c->SaveAs (Form ("%s/Plots/PtCh/IpPb_CMS_FF_Comp_Summary_60GeVJets_%s.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
    } // end iDir=0 scope
  } // end loop over iPtJInt



  { 

    const short iPtJInt = 1;
    const TString pTJInt = "60GeV";

    {
      const short iDir = 0; // FF should only be compared to near-side.
      const char* canvasName = Form ("c_jetInt_trk_pt_IpPb_CMS_RpPb_Comp_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1200, 800);
      c->Divide (3, 2);

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        gPad->SetLogx ();

        TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#it{I}_{#it{p}Pb}", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser ((iPtJInt == 0 ? 0.61 : 0.74), (iPtJInt == 0 ? 1.6 : 1.45));
        //h->GetYaxis ()->SetRangeUser (0.0, 3);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        TGAE* g = nullptr;

        g = (TGAE*) g_RpPb_CMS[0]->Clone ("gtemp");
        myDrawSyst (g, kBlack);
        SaferDelete (&g);

        g = (TGAE*) g_RpPb_CMS[1]->Clone ("gtemp");
        RecenterGraph (g);
        ResetXErrors (g);
        myDraw (g, kBlack, kFullSquare, 1.4, 1, 2, "P", false);
        SaferDelete (&g);

        g = (TGAE*) g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_iaa[0][iPtJInt][iDir][iCent];
        //h = h_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iMCTruthJetsTruthParts];
        SetCentralValuesKeepRelativeErrors (g, h);
        RecenterGraph (g);
        myDrawSystFill (g, colorfulSystColors[nZdcCentBins-iCent], 0.6, 1001);
        SaferDelete (&g);

        g = make_graph (h);
        RecenterGraph (g);
        ResetXErrors (g);
        myDraw (g, colorfulColors[nZdcCentBins-iCent], kFullCircle, 1.0, 1, 2, "P", false);
        SaferDelete (&g);

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, colorfulColors[nZdcCentBins-iCent], Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, colorfulColors[0], "#bf{#it{p}+Pb, 0-100%}", 0.05);

      } // end loop over iCent

      c->cd ();
      myText (0.065, 0.971, kBlack, "#bf{#it{ATLAS}} Internal", 0.027);
  
      c->cd (1);
      myText (0.2, 0.80, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.05);
      myText (0.2, 0.74, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.05);
      c->cd (2);
      myText (0.2, 0.80, kBlack, "#it{p}_{T}^{jet} > 60 GeV", 0.05);
      myText (0.2, 0.74, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.05);
      c->cd (3);
      myLineText2 (0.26, 0.80, kBlack, kFullCircle, "This result", 1.2, 0.05);
      myLineText2 (0.26, 0.73, kBlack, kFullSquare, "CMS h^{#pm} #it{R}_{#it{p}Pb}", 1.2, 0.05);

      c->SaveAs (Form ("%s/Plots/PtCh/IpPb_CMS_RpPb_Comp_Summary_60GeVJets_%s.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
    } // end iDir=0 scope
  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (short iDir : {0, 2}) {
      const char* canvasName = Form ("c_jetInt_trk_pt_IpPb_Angantyr_Comp_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1200, 800);
      c->Divide (3, 2);

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        //gPad->SetLogx ();

        TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#it{I}_{#it{p}Pb}", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
        //h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser ((iPtJInt == 0 ? 0.61 : 0.74), (iPtJInt == 0 ? 1.6 : 1.45));
        //h->GetYaxis ()->SetRangeUser (0.0, 3);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        TGAE* g = nullptr;

        g = (TGAE*) g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_iaa[0][iPtJInt][iDir][iCent];
        //h = h_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iMCTruthJetsTruthParts];
        SetCentralValuesKeepRelativeErrors (g, h);
        TrimGraph (g, 4, 90);
        RecenterGraph (g);
        myDrawSystFill (g, colorfulSystColors[nZdcCentBins-iCent], 0.6, 1001);
        SaferDelete (&g);

        g = make_graph (h);
        TrimGraph (g, 4, 90);
        //RecenterGraph (g);
        ResetXErrors (g);
        myDraw (g, colorfulColors[nZdcCentBins-iCent], kFullCircle, 1.0, 1, 2, "P", false);
        SaferDelete (&g);


        g = (TGAE*) g_angantyr_iaa[iPtJInt][iDir][iCent][0]->Clone ("gtemp");
        TrimGraph (g, 4, iPtJInt == 0 ? 30 : 70);
        TGAE* gup = new TGAE ();
        TGAE* gdown = new TGAE ();
        MakeGupAndGdown (g, gup, gdown);
        myDrawFill (gup, gdown, myLitePurple, 0.7);
        ResetTGAEErrors (g);
        ResetXErrors (g);
        myDraw (g, myViolet, kDot, 0, 2, 2, "L");
        SaferDelete (&g);
        SaferDelete (&gup);
        SaferDelete (&gdown);


        g = (TGAE*) g_angantyr_iaa[iPtJInt][iDir][iCent][1]->Clone ("gtemp");
        TrimGraph (g, 4, iPtJInt == 0 ? 30 : 70);
        gup = new TGAE ();
        gdown = new TGAE ();
        MakeGupAndGdown (g, gup, gdown);
        myDrawFill (gup, gdown, myLiteBlue, 0.5);
        ResetTGAEErrors (g);
        ResetXErrors (g);
        myDraw (g, myBlue, kDot, 0, 1, 2, "L");
        SaferDelete (&g);
        SaferDelete (&gup);
        SaferDelete (&gdown);


        //g = (TGAE*) g_angantyr_iaa[iPtJInt][iDir][iCent][2]->Clone ("gtemp");
        //TrimGraph (g, 4, iPtJInt == 0 ? 30 : 70);
        //gup = new TGAE ();
        //gdown = new TGAE ();
        //MakeGupAndGdown (g, gup, gdown);
        //myDrawFill (gup, gdown, myCyan, 0.7);
        //ResetTGAEErrors (g);
        //ResetXErrors (g);
        //myDraw (g, myCyan, kDot, 0, 3, 2, "L");
        //SaferDelete (&g);
        //SaferDelete (&gup);
        //SaferDelete (&gdown);


        //g = (TGAE*) g_angantyr_iaa[iPtJInt][iDir][iCent][3]->Clone ("gtemp");
        //TrimGraph (g, 4, iPtJInt == 0 ? 30 : 70);
        //gup = new TGAE ();
        //gdown = new TGAE ();
        //MakeGupAndGdown (g, gup, gdown);
        //myDrawFill (gup, gdown, myGreen, 0.5);
        //ResetTGAEErrors (g);
        //ResetXErrors (g);
        //myDraw (g, myGreen, kDot, 0, 4, 2, "L");
        //SaferDelete (&g);
        //SaferDelete (&gup);
        //SaferDelete (&gdown);

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, colorfulColors[nZdcCentBins-iCent], Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, colorfulColors[0], "#bf{#it{p}+Pb, 0-100%}", 0.05);

      } // end loop over iCent

      c->cd ();
      myText (0.065, 0.971, kBlack, "#bf{#it{ATLAS}} Internal", 0.027);
  
      c->cd (1);
      myText (0.2, 0.80, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.05);
      myText (0.2, 0.74, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.05);
      myText (0.2, 0.68, kBlack, "Pythia 8.306", 0.05);
      c->cd (2);
      myText (0.2, 0.80, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV", iPtJInt == 0 ? 30 : 60), 0.05);
      myText (0.2, 0.74, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.05);
      c->cd (3);
      myLineText2 (0.26, 0.80, kBlack, kFullCircle, "Data", 1.2, 0.05);
      myLineText (0.26, 0.73, myPurple, 2, "#scale[0.8]{#bf{ANGANTYR}}, w/ EPPS16 (NLO)", 1.5, 0.05);
      myLineText (0.26, 0.66, myBlue, 1, "#scale[0.8]{#bf{ANGANTYR}}, no nPDF", 1.5, 0.05);
      //myLineText (0.26, 0.59, myCyan, 3, "#scale[0.8]{#bf{ANGANTYR}}, no rescatter, w/ nPDF", 1.5, 0.05);
      //myLineText (0.26, 0.52, myGreen, 4, "#scale[0.8]{#bf{ANGANTYR}}, no rescatter, no nPDF", 1.5, 0.05);

      c->SaveAs (Form ("%s/Plots/PtCh/IpPb_Angantyr_Comp_Summary_%sJets_%s.pdf", workPath.Data (), pTJInt.Data (), directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
    } // end iDir=0 scope
  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (short iDir : {0, 2}) {
      const char* canvasName = Form ("c_jetInt_trk_pt_IpPb_Angantyr_0-20perc_Comp_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 800, 800);
      c->cd ();

      const short iCent = 4;
      gPad->SetLogx ();

      //TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#it{I}_{#it{p}Pb}", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
      TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#it{I}_{#it{p}Pb}", 1, 4, iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (0.87, 1.35);
      //h->GetYaxis ()->SetRangeUser (0.0, 3);
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      TGAE* g = nullptr;

      g = (TGAE*) g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0]->Clone ();
      h = h_jetInt_trk_pt_iaa[0][iPtJInt][iDir][iCent];
      //h = h_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iMCTruthJetsTruthParts];
      SetCentralValuesKeepRelativeErrors (g, h);
      RecenterGraph (g);
      //myDrawSystFill (g, colorfulSystColors[0], 0.6, 1001);
      myDrawSyst (g, kBlack, 1, 2);
      SaferDelete (&g);

      g = make_graph (h);
      RecenterGraph (g);
      ResetXErrors (g);
      //myDraw (g, colorfulColors[0], kFullCircle, 1.5, 1, 3, "P", false);
      myDraw (g, kBlack, 53, 1.8, 1, 3, "P", false);
      SaferDelete (&g);


      g = (TGAE*) g_angantyr_iaa[iPtJInt][iDir][iCent][0]->Clone ("gtemp");
      RecenterGraph (g);
      TrimGraph (g, 4, iPtJInt == 0 ? 30 : 70);
      TGAE* gup = new TGAE ();
      TGAE* gdown = new TGAE ();
      MakeGupAndGdown (g, gup, gdown);
      myDrawFill (gup, gdown, myLitePurple, 0.7);
      ResetTGAEErrors (g);
      ResetXErrors (g);
      myDraw (g, myViolet, kDot, 0, 2, 2, "L");
      SaferDelete (&g);
      SaferDelete (&gup);
      SaferDelete (&gdown);


      g = (TGAE*) g_angantyr_iaa[iPtJInt][iDir][iCent][1]->Clone ("gtemp");
      RecenterGraph (g);
      TrimGraph (g, 4, iPtJInt == 0 ? 30 : 70);
      gup = new TGAE ();
      gdown = new TGAE ();
      MakeGupAndGdown (g, gup, gdown);
      myDrawFill (gup, gdown, myLiteBlue, 0.5);
      ResetTGAEErrors (g);
      ResetXErrors (g);
      myDraw (g, myBlue, kDot, 0, 1, 2, "L");
      SaferDelete (&g);
      SaferDelete (&gup);
      SaferDelete (&gdown);


      //g = (TGAE*) g_angantyr_iaa[iPtJInt][iDir][iCent][2]->Clone ("gtemp");
      //TrimGraph (g, 4, iPtJInt == 0 ? 30 : 70);
      //gup = new TGAE ();
      //gdown = new TGAE ();
      //MakeGupAndGdown (g, gup, gdown);
      //myDrawFill (gup, gdown, myCyan, 0.7);
      //ResetTGAEErrors (g);
      //ResetXErrors (g);
      //myDraw (g, myCyan, kDot, 0, 3, 2, "L");
      //SaferDelete (&g);
      //SaferDelete (&gup);
      //SaferDelete (&gdown);


      //g = (TGAE*) g_angantyr_iaa[iPtJInt][iDir][iCent][3]->Clone ("gtemp");
      //TrimGraph (g, 4, iPtJInt == 0 ? 30 : 70);
      //gup = new TGAE ();
      //gdown = new TGAE ();
      //MakeGupAndGdown (g, gup, gdown);
      //myDrawFill (gup, gdown, myGreen, 0.5);
      //ResetTGAEErrors (g);
      //ResetXErrors (g);
      //myDraw (g, myGreen, kDot, 0, 4, 2, "L");
      //SaferDelete (&g);
      //SaferDelete (&gup);
      //SaferDelete (&gdown);

      myText (0.24, 0.890, kBlack, "#bf{#it{ATLAS}} Internal", 0.034);
      myText (0.24, 0.850, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, ZDC 0-20%", 0.034);
      myText (0.24, 0.810, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.034);
      myText (0.24, 0.770, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} %s", iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.034);
      mySimpleMarkerAndBoxAndLineText (0.32, 0.725, 1.5, 1001, kBlack, 0.0, kBlack, 53, 1.8, "Data", 0.034);

      mySimpleMarkerAndBoxAndLineText (0.32, 0.685, 1.5, 1001, myLitePurple, 0.7, myViolet, kDot, 0.0, "#scale[0.8]{#bf{ANGANTYR}}, w/ EPPS16 (NLO)", 0.034, 2);
      mySimpleMarkerAndBoxAndLineText (0.32, 0.645, 1.5, 1001, myLiteBlue, 0.5, myBlue, kDot, 0.0, "#scale[0.8]{#bf{ANGANTYR}}, no nPDF", 0.034, 1);
      //mySimpleMarkerAndBoxAndLineText (0.32, 0.615, 1.5, 1001, myCyan, 0.7, myCyan, kDot, 0.0, "#scale[0.8]{#bf{ANGANTYR}}, no rescatter, w/ EPPS16 (NLO)", 0.034, 3);
      //mySimpleMarkerAndBoxAndLineText (0.32, 0.580, 1.5, 1001, myGreen, 0.5, myGreen, kDot, 0.0, "#scale[0.8]{#bf{ANGANTYR}}, no rescatter, no nPDF", 0.034, 4);

      c->SaveAs (Form ("%s/Plots/PtCh/IpPb_Angantyr_Comp_0-20perc_%sJets_%s.pdf", workPath.Data (), pTJInt.Data (), directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
    } // end iDir=0 scope
  } // end loop over iPtJInt




  {
    for (short iDir : {0, 2}) {
      const char* canvasName = Form ("c_jetInt_trk_pt_IpPb_Angantyr_0-20perc_Comp_iDir%i_allPtJet", iDir);
      TCanvas* c = new TCanvas (canvasName, "", 800, 800);
      c->cd ();

      const short iCent = 4;
      gPad->SetLogx ();

      //TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#it{I}_{#it{p}Pb}", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
      TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#it{I}_{#it{p}Pb}", 1, 3, iDir == 1 ? 10 : pTChBins[nPtChBins-1]);//pTChBins[nPtChBins]);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (0.67, 1.37);
      //h->GetYaxis ()->SetRangeUser (0.0, 3);
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      TGAE* g = nullptr;

      gStyle->SetEndErrorSize (5);
      for (int iPtJInt : {0, 1}) {

        g = (TGAE*) g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0]->Clone ();
        h = h_jetInt_trk_pt_iaa[0][iPtJInt][iDir][iCent];
        //h = h_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iMCTruthJetsTruthParts];
        SetCentralValuesKeepRelativeErrors (g, h);
        TrimGraph (g, 3, iPtJInt == 0 ? 60 : 90);
        RecenterGraph (g);
        //myDrawSystFill (g, colorfulSystColors[0], 0.6, 1001);
        if (iPtJInt == 1)
          //myDrawSyst (g, kBlack, 1, 2);//, 0., "[]");
          myDrawSyst (g, kGray+1, 1, 0, 0.5);
        else
          myDrawSyst (g, myLiteYellow, 1, 0, 0.35);
        SaferDelete (&g);

        g = make_graph (h);
        TrimGraph (g, 3, iPtJInt == 0 ? 60 : 90);
        RecenterGraph (g);
        deltaize (g, iPtJInt == 0 ? 0.988 : 1.012, true);
        ResetXErrors (g);
        //myDraw (g, colorfulColors[0], kFullCircle, 1.5, 1, 3, "P", false);
        if (iPtJInt == 1)
          myDraw (g, kBlack, 53, 1.8, 1, 3, "P", false);
        else
          myDraw (g, kBlack, 21, 1.8, 1, 3, "P", true);
        SaferDelete (&g);
      }


      for (int iPtJInt : {0, 1}) {
        g = (TGAE*) g_angantyr_iaa[iPtJInt][iDir][iCent][0]->Clone ("gtemp");
        RecenterGraph (g);
        TrimGraph (g, 4, iPtJInt == 0 ? 30 : 70);
        TGAE* gup = new TGAE ();
        TGAE* gdown = new TGAE ();
        MakeGupAndGdown (g, gup, gdown);
        myDrawFill (gup, gdown, iPtJInt == 1 ? myLitePurple : myLiteRed, 0.7);
        ResetTGAEErrors (g);
        ResetXErrors (g);
        myDraw (g, iPtJInt == 1 ? myViolet : myRed, kDot, 0, 2, 2, "L");
        SaferDelete (&g);
        SaferDelete (&gup);
        SaferDelete (&gdown);


        if (iPtJInt == 1) {
          g = (TGAE*) g_angantyr_iaa[iPtJInt][iDir][iCent][1]->Clone ("gtemp");
          RecenterGraph (g);
          TrimGraph (g, 4, iPtJInt == 0 ? 30 : 70);
          gup = new TGAE ();
          gdown = new TGAE ();
          MakeGupAndGdown (g, gup, gdown);
          myDrawFill (gup, gdown, iPtJInt == 1 ? myLiteBlue : myLiteYellow, 0.5);
          ResetTGAEErrors (g);
          ResetXErrors (g);
          myDraw (g, myBlue, kDot, 0, 1, 2, "L");
          SaferDelete (&g);
          SaferDelete (&gup);
          SaferDelete (&gdown);
        }

      }

      myText (0.24, 0.890, kBlack, "#bf{#it{ATLAS}} Internal", 0.034);
      myText (0.24, 0.850, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, ZDC 0-20%", 0.032);
      myText (0.24, 0.810, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
      myText (0.24, 0.770, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.032);

      mySimpleMarkerAndBoxAndLineText (0.29, 0.360, 1.5, 1001, kGray+1,       0.50, kBlack,  53, 1.8, "Data, #it{p}_{T}^{jet} > 60 GeV", 0.028);
      mySimpleMarkerAndBoxAndLineText (0.29, 0.320, 1.5, 1001, myLiteYellow,  0.35, kBlack,  21, 1.8, "Data, #it{p}_{T}^{jet} > 30 GeV", 0.028);

      mySimpleMarkerAndBoxAndLineText (0.29, 0.280, 1.5, 1001, myLitePurple, 0.7, myViolet, kDot, 0.0, "#scale[0.8]{#bf{ANGANTYR}}, EPPS16 (NLO), #it{p}_{T}^{jet} > 60 GeV", 0.028, 2);
      mySimpleMarkerAndBoxAndLineText (0.29, 0.240, 1.5, 1001, myLiteRed,    0.7, myRed,    kDot, 0.0, "#scale[0.8]{#bf{ANGANTYR}}, EPPS16 (NLO), #it{p}_{T}^{jet} > 30 GeV", 0.028, 2);
      mySimpleMarkerAndBoxAndLineText (0.29, 0.200, 1.5, 1001, myLiteBlue,   0.5, myBlue,   kDot, 0.0, "#scale[0.8]{#bf{ANGANTYR}}, no nPDF, #it{p}_{T}^{jet} > 60 GeV",      0.028, 1);

      c->SaveAs (Form ("%s/Plots/PtCh/IpPb_Angantyr_Comp_0-20perc_AllPtJets_%s.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
    } // end iDir=0 scope
    gStyle->SetEndErrorSize (0);
  } // end loop over iPtJInt




  { 
    const short iPtJInt = 1;
    const TString pTJInt = "60GeV";

    for (short iDir : {0, 2}) {
      const char* canvasName = Form ("c_jetInt_trk_pt_IpPb_AMPT_Comp_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 800, 800);
      c->cd ();

      const short iCent = 4;
      gPad->SetLogx ();

      //TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#it{I}_{#it{p}Pb}", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
      TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#it{I}_{#it{p}Pb}", 1, 4, iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (0.87, 1.35);
      //h->GetYaxis ()->SetRangeUser (0.0, 3);
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      TGAE* g = nullptr;

      g = (TGAE*) g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0]->Clone ();
      h = h_jetInt_trk_pt_iaa[0][iPtJInt][iDir][iCent];
      //h = h_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iMCTruthJetsTruthParts];
      SetCentralValuesKeepRelativeErrors (g, h);
      RecenterGraph (g);
      //myDrawSystFill (g, colorfulSystColors[0], 0.6, 1001);
      myDrawSyst (g, kBlack, 1, 2);
      SaferDelete (&g);

      g = make_graph (h);
      RecenterGraph (g);
      ResetXErrors (g);
      //myDraw (g, colorfulColors[0], kFullCircle, 1.5, 1, 3, "P", false);
      myDraw (g, kBlack, 53, 1.8, 1, 3, "P", false);
      SaferDelete (&g);


      g = make_graph (h_ampt_iaa[iDir][0]);
      TrimGraph (g, 5, 100);
      RecenterGraph (g);
      ResetTGAEErrors (g);
      myDraw (g, kMagenta, 72, 1.8, 1, 3, "LP", false);
      SaferDelete (&g);


      g = make_graph (h_ampt_iaa[iDir][1]);
      TrimGraph (g, 5, 100);
      RecenterGraph (g);
      ResetTGAEErrors (g);
      myDraw (g, myViolet, kFullSquare, 1.8, 1, 3, "LP", false);
      SaferDelete (&g);

      myText (0.24, 0.890, kBlack, "#bf{#it{ATLAS}} Internal", 0.034);
      myText (0.24, 0.850, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, ZDC 0-20%", 0.034);
      myText (0.24, 0.810, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.034);
      myText (0.24, 0.770, kBlack, Form ("#it{p}_{T}^{jet} > 60 GeV, #Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.034);
      //mySimpleMarkerAndBoxAndLineText (0.32, 0.725, 1.5, 1001, colorfulSystColors[0], 1.0, colorfulColors[0], kFullCircle, 1.6, "Data", 0.034);
      mySimpleMarkerAndBoxAndLineText (0.32, 0.725, 1.5, 1001, kBlack, 0.0, kBlack, 53, 1.8, "Data", 0.034);
      mySimpleMarkerAndBoxAndLineText (0.32, 0.685, 1.5, 1001, kWhite, 0.0, myViolet, kFullSquare, 1.8, "AMPT, w/ FS interactions", 0.034);
      mySimpleMarkerAndBoxAndLineText (0.32, 0.645, 1.5, 1001, kWhite, 0.0, kMagenta, 72, 1.8, "AMPT, no FS interactions", 0.034);

      c->SaveAs (Form ("%s/Plots/PtCh/IpPb_AMPT_Comp_0-20perc_60GeVJets_%s.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
    } // end iDir=0 scope
  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) { 
    const TString pTJInt = iPtJInt == 0 ? "30GeV" : "60GeV";

    {
      const short iDir = 2;
      const char* canvasName = Form ("c_jetInt_trk_pt_IpPb_Zh_Comp_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 800, 800);
      c->cd ();

      const short iCent = 4;
      gPad->SetLogx ();
      gPad->SetLogy ();

      TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];Ratio to #it{pp}", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (iPtJInt == 0 ? 0.1 : 0.2, iPtJInt == 0 ? 3.5 : 4.5);
      h->GetYaxis ()->SetMoreLogLabels ();
      //h->GetYaxis ()->SetRangeUser (0.0, 3);
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      TGAE* g = nullptr;

      const Color_t jhCol = colorfulColors[5];
      const Color_t jhSystCol = colorfulSystColors[5];

      g = (TGAE*) g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0]->Clone ();
      h = h_jetInt_trk_pt_iaa[0][iPtJInt][iDir][iCent];
      SetCentralValuesKeepRelativeErrors (g, h);
      RecenterGraph (g);
      myDrawSystFill (g, jhSystCol, 0.6, 1001);
      SaferDelete (&g);

      g = make_graph (h);
      TrimGraph (g, 0.5, 90);
      RecenterGraph (g);
      ResetXErrors (g);
      myDraw (g, jhCol, 53, 1.8, 1, 3, "PL", false);
      SaferDelete (&g);


      const Color_t zhCol = colorfulColors[1];
      const Color_t zhSystCol = colorfulSystColors[1];

      g = (TGAE*) g_Zh_iaa[iPtJInt][1]->Clone ();
      myDrawSystFill (g, zhSystCol, 0.5, 1001);
      SaferDelete (&g);

      g = (TGAE*) g_Zh_iaa[iPtJInt][0]->Clone ();
      RecenterGraph (g);
      ResetXErrors (g);
      myDraw (g, zhCol, kFullSquare, 1.8, 1, 3, "PL", false);
      //ResetTGAEErrors (g);
      //myDraw (g, kBlack, 54, 1.8, 1, 3, "P", false);
      //SaferDelete (&g);


      myText (0.52, 0.890, kBlack, "#bf{#it{ATLAS}} Internal", 0.034);
      myText (0.52, 0.845, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.034);
      myText (0.52, 0.800, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.034);
      myText (0.52, 0.755, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.034);
      mySimpleMarkerAndBoxAndLineText (0.28, 0.355, 1.5, 1001, jhSystCol, 0.6, jhCol, 53, 1.8, "0-20\% #it{p}+Pb jet-#it{h}", 0.034);
      myText (0.28, 0.310, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{h,jet} %s", iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.034);
      mySimpleMarkerAndBoxAndLineText (0.28, 0.255, 1.5, 1001, zhSystCol, 0.5, zhCol, kFullSquare, 1.8, "0-10\% Pb+Pb #it{Z}-#it{h}", 0.034);
      myText (0.28, 0.210, kBlack, Form ("#it{p}_{T}^{Z} %s, #Delta#phi_{hZ} > 3#pi/4", iPtJInt == 0 ? "= 30-60 GeV" : "> 60 GeV"), 0.034);
      c->SaveAs (Form ("%s/Plots/PtCh/IpPb_Zh_Comp_0-20perc_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));
    } // end iDir=0 scope
  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (short iDir : {0, 1, 2}) {

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

      for (short iCent = 0; iCent < nZdcCentBins; iCent++) {
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
      for (short iCent = 0; iCent < nZdcCentBins; iCent++)
        myLineText2 (0.25, 0.68-iCent*0.04, colors[iCent], kOpenCircle, Form ("#bf{#it{p}+Pb, %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 1.4, 0.032, true);

      c->SaveAs (Form ("%s/Plots/PtCh/SigToBkgd_%iGeVJets_%s_ptch.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));

    } // end loop over iDir

  } // end loop over iPtJInt



  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (short iDir : {0, 1, 2}) {

      const char* canvasName = Form ("c_jetInt_trk_pt_%s_sig2bkg2ncoll_%s", directions[iDir].Data (), pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 800, 800);

      TH1D* h = nullptr; 
      TGAE* g = nullptr;

      c->cd (); 
      c->SetLogx ();
      c->SetLogy ();

      float ymin = 1e-1;
      float ymax = 1e4;

      c->Clear ();

      h = (TH1D*) h_jetInt_trk_pt_ref_sig[0][iPtJInt][iDir]->Clone ("h");
      h->Divide (h_jetInt_trk_pt_ref_bkg[0][iPtJInt][iDir]);
      g = make_graph (h);
      ResetXErrors (g);

      g->GetXaxis ()->SetRangeUser (0.5, 12);
      g->GetXaxis ()->SetMoreLogLabels ();
      g->GetYaxis ()->SetRangeUser (ymin, ymax);
      g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      //g->GetXaxis ()->SetTitleSize (0.028);
      //g->GetXaxis ()->SetLabelSize (0.028);
      g->GetYaxis ()->SetTitle ("Sig. / Bkgd. #times #sqrt{#LTN_{part}#GT}");
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

      for (short iCent = 0; iCent < nZdcCentBins; iCent++) {
        h = (TH1D*) h_jetInt_trk_pt_sig[0][iPtJInt][iDir][iCent]->Clone ("htemp");
        h->Divide (h_jetInt_trk_pt_bkg[0][iPtJInt][iDir][iCent]);
        g = make_graph (h);
        ScaleGraph (g, nullptr, std::sqrt (zdcNcollValues[iCent]+1));
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
      for (short iCent = 0; iCent < nZdcCentBins; iCent++)
        myLineText2 (0.25, 0.68-iCent*0.04, colors[iCent], kOpenCircle, Form ("#bf{#it{p}+Pb, %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 1.4, 0.032, true);

      c->SaveAs (Form ("%s/Plots/PtCh/SigToBkgd_OverNcoll_%iGeVJets_%s_ptch.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));

    } // end loop over iDir

  } // end loop over iPtJInt




  for (short iSysType : {0, 1, 2}) {

    const short iDType = iSysType % 2;
    const TString sysType = (iSysType == 0 ? "tracking" : (iSysType == 1 ? "jets" : "mixing"));

    const short maxNSys = (iSysType == 1 ? 10 : 6);
    const float tsize = (iSysType == 1 ? 0.024 : 0.030);
    const float x0 = (iSysType == 1 ? 0.24 : 0.35);
    const float dx = (iSysType == 1 ? 0.41 : 0.31);
    const float y0 = (iSysType == 1 ? 0.44 : 0.41);
    const float dy = (iSysType == 1 ? 0.028 : 0.04);


    for (short iDir : {0, 2}) {

      //if (makeTotalSystPlots) {
      //  for (short iPtJInt : {0, 1}) {

      //    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      //    const char* canvasName = Form ("c_jetInt_trk_pt_%s_pp_%s_syst_%s", directions[iDir].Data (), sysType.Data (), pTJInt.Data ());
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

      //    {
      //      const TString totVar = totalVariations[iSysType];
      //      g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_ref_systTot[iPtJInt][iDir][iSysType], true, 100);
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

      //    for (short iVar = 1; iVar < nVar; iVar++) {
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

      //    short count = 0;
      //    for (short iVar = 1; iVar < nVar; iVar++) {
      //      const TString var = variations[iVar];
      //      if ((iSysType == 0  && !IsTrackingVariation (var)) || (iSysType == 1 && !IsJetsVariation (var)) || (iSysType == 2 && !IsMixingVariation (var)))
      //        continue;
      //      myLineColorText (x0+dx*(count/maxNSys), y0-(count%maxNSys)*dy, varStyles[var].first, varStyles[var].second, varFullNames[var], 0.5, tsize);
      //      count++;
      //    }
      //    {
      //      const TString totVar = totalVariations[iSysType];
      //      myLineColorText (x0+dx*(count/maxNSys), y0-(count%maxNSys)*dy, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 0.5, tsize);
      //      count++;
      //    }

      //    c->SaveAs (Form ("%s/Plots/Systematics/PtCh/TotalJetTaggedYield_pp_%s_ptch_%iGeVJets_%s_syst.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60, sysType.Data ()));

      //  } // end loop over iPtJInt


      //  for (short iPtJInt : {0, 1}) {

      //    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      //    for (short iCent = 0; iCent < nZdcCentBins; iCent++) {

      //      const char* canvasName = Form ("c_jetInt_trk_pt_%s_%s_pPb_iCent%i_syst_%s", directions[iDir].Data (), sysType.Data (), iCent, pTJInt.Data ());
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

      //      {
      //        const TString totVar = totalVariations[iSysType];
      //        g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_systTot[iPtJInt][iDir][iCent][iSysType], true, 100);
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

      //      for (short iVar = 1; iVar < nVar; iVar++) {
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

      //      short count = 0;
      //      for (short iVar = 1; iVar < nVar; iVar++) {
      //        const TString var = variations[iVar];
      //        if ((iSysType == 0  && !IsTrackingVariation (var)) || (iSysType == 1 && !IsJetsVariation (var)) || (iSysType == 2 && !IsMixingVariation (var)))
      //          continue;
      //        myLineColorText (x0+dx*(count/maxNSys), y0-(count%maxNSys)*dy, varStyles[var].first, varStyles[var].second, varFullNames[var], 0.5, tsize);
      //        count++;
      //      }
      //        const TString totVar = totalVariations[iSysType];
      //        myLineColorText (x0+dx*(count/maxNSys), y0-(count%maxNSys)*dy, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 0.5, tsize);
      //        count++;
      //      }
      //      
      //      c->SaveAs (Form ("%s/Plots/Systematics/PtCh/TotalJetTaggedYield_pPb_%i-%iperc_%s_ptch_%iGeVJets_%s_syst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60, sysType.Data ()));
      //    } // end loop over iCent

      //  } // end loop over iPtJInt
      //}



      //if (makeBkgdSystPlots && iDType != 1) {
      //  for (short iPtJInt : {0, 1}) {

      //    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      //    const char* canvasName = Form ("c_jetInt_trk_pt_%s_pp_%s_bkg_syst_%s", directions[iDir].Data (), sysType.Data (), pTJInt.Data ());
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

      //    {
      //      const TString totVar = totalVariations[iSysType];
      //      g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_ref_bkg_systTot[iPtJInt][iDir][iSysType], true, 100);
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

      //    for (short iVar = 1; iVar < nVar; iVar++) {
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

      //    short count = 0;
      //    for (short iVar = 1; iVar < nVar; iVar++) {
      //      const TString var = variations[iVar];
      //      if ((iSysType == 0  && !IsTrackingVariation (var)) || (iSysType == 1 && !IsJetsVariation (var)) || (iSysType == 2 && !IsMixingVariation (var)))
      //        continue;
      //      myLineColorText (x0+dx*(count/maxNSys), y0-(count%maxNSys)*dy, varStyles[var].first, varStyles[var].second, varFullNames[var], 0.5, tsize);
      //      count++;
      //    }
      //    {
      //      const TString totVar = totalVariations[iSysType];
      //      myLineColorText (x0+dx*(count/maxNSys), y0-(count%maxNSys)*dy, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 0.5, tsize);
      //      count++;
      //    }

      //    c->SaveAs (Form ("%s/Plots/Systematics/PtCh/BkgdJetTaggedYield_pp_%s_ptch_%iGeVJets_%s_syst.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60, sysType.Data ()));

      //  } // end loop over iPtJInt


      //  for (short iPtJInt : {0, 1}) {

      //    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      //    for (short iCent = 0; iCent < nZdcCentBins; iCent++) {

      //      const char* canvasName = Form ("c_jetInt_trk_pt_%s_%s_pPb_iCent%i_bkgd_syst_%s", directions[iDir].Data (), sysType.Data (), iCent, pTJInt.Data ());
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

      //      {
      //        const TString totVar = totalVariations[iSysType];
      //        g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_bkg_systTot[iDir][iPtJInt][iCent][iSysType], true, 100);
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

      //      for (short iVar = 1; iVar < nVar; iVar++) {
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

      //      short count = 0;
      //      for (short iVar = 1; iVar < nVar; iVar++) {
      //        const TString var = variations[iVar];
      //        if ((iSysType == 0  && !IsTrackingVariation (var)) || (iSysType == 1 && !IsJetsVariation (var)) || (iSysType == 2 && !IsMixingVariation (var)))
      //          continue;
      //        myLineColorText (x0+dx*(count/maxNSys), y0-(count%maxNSys)*dy, varStyles[var].first, varStyles[var].second, varFullNames[var], 0.5, tsize);
      //        count++;
      //      }
      //      {
      //        const TString totVar = totalVariations[iSysType];
      //        myLineColorText (x0+dx*(count/maxNSys), y0-(count%maxNSys)*dy, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 0.5, tsize);
      //        count++;
      //      }
      //      
      //      c->SaveAs (Form ("%s/Plots/Systematics/PtCh/BkgdJetTaggedYield_pPb_%i-%iperc_%s_ptch_%iGeVJets_%s_syst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60, sysType.Data ()));

      //    } // end loop over iCent

      //  } // end loop over iPtJInt
      //}



      if (makeSigSystPlots) {
        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          const char* canvasName = Form ("c_jetInt_trk_pt_%s_unf_%s_pp_syst_%s", directions[iDir].Data (), sysType.Data (), pTJInt.Data ());
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

          {
            const TString totVar = totalVariations[iSysType];
            g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_ref_unf_systTot[iPtJInt][iDir][iSysType], true, 100);
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

          for (short iVar = 1; iVar < nVar; iVar++) {
            const TString var = variations[iVar];
            if ((iSysType == 0  && !IsTrackingVariation (var)) || (iSysType == 1 && !IsJetsVariation (var)) || (iSysType == 2 && !IsMixingVariation (var)))
              continue;

            g = (TGAE*) g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar]->Clone ("gtemp");
            SaveRelativeErrors (g, g, g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][0], 100);
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

          short count = 0;
          for (short iVar = 1; iVar < nVar; iVar++) {
            const TString var = variations[iVar];
            if ((iSysType == 0  && !IsTrackingVariation (var)) || (iSysType == 1 && !IsJetsVariation (var)) || (iSysType == 2 && !IsMixingVariation (var)))
              continue;
            myLineColorText (x0+dx*(count/maxNSys), y0-(count%maxNSys)*dy, varStyles[var].first, varStyles[var].second, varFullNames[var], 0.5, tsize);
            count++;
          }
          {
            const TString totVar = totalVariations[iSysType];
            myLineColorText (x0+dx*(count/maxNSys), y0-(count%maxNSys)*dy, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 0.5, tsize);
            count++;
          }

          c->SaveAs (Form ("%s/Plots/Systematics/PtCh/SignalJetTaggedYield_pp_%s_ptch_%iGeVJets_%s_syst.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60, sysType.Data ()));

        } // end loop over iPtJInt



        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          for (short iCent = 0; iCent < nZdcCentBins; iCent++) {

            const char* canvasName = Form ("c_jetInt_trk_pt_%s_unf_%s_pPb_iCent%i_syst_%s", directions[iDir].Data (), sysType.Data (), iCent, pTJInt.Data ());
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

            {
              const TString totVar = totalVariations[iSysType];
              g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_unf_systTot[iPtJInt][iDir][iCent][iSysType], true, 100);
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

            for (short iVar = 1; iVar < nVar; iVar++) {
              const TString var = variations[iVar];
              if ((iSysType == 0  && !IsTrackingVariation (var)) || (iSysType == 1 && !IsJetsVariation (var)) || (iSysType == 2 && !IsMixingVariation (var)))
                continue;

              g = (TGAE*) g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar]->Clone ("gtemp");
              SaveRelativeErrors (g, g, g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][0], 100);
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

            short count = 0;
            for (short iVar = 1; iVar < nVar; iVar++) {
              const TString var = variations[iVar];
              if ((iSysType == 0  && !IsTrackingVariation (var)) || (iSysType == 1 && !IsJetsVariation (var)) || (iSysType == 2 && !IsMixingVariation (var)))
                continue;
              myLineColorText (x0+dx*(count/maxNSys), y0-(count%maxNSys)*dy, varStyles[var].first, varStyles[var].second, varFullNames[var], 0.5, tsize);
              count++;
            }
            {
              const TString totVar = totalVariations[iSysType];
              myLineColorText (x0+dx*(count/maxNSys), y0-(count%maxNSys)*dy, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 0.5, tsize);
              count++;
            }

            c->SaveAs (Form ("%s/Plots/Systematics/PtCh/SignalJetTaggedYield_pPb_%i-%iperc_%s_ptch_%iGeVJets_%s_syst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60, sysType.Data ()));
          } // end loop over iCent

        } // end loop over iPtJInt



        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          const char* canvasName = Form ("c_jetInt_trk_pt_%s_unf_%s_pp_syst_%s", directions[iDir].Data (), sysType.Data (), pTJInt.Data ());
          TCanvas* c = new TCanvas (canvasName, "", 2800, 1400);
          c->Divide (4, 2);

          std::vector <short> totVars (0);
          if (iDType == 0) {
            totVars.push_back (0);
            totVars.push_back (1);
          }
          else 
            totVars.push_back (2);

          {
            TH1D* h = nullptr;
            TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

            c->cd (7); 
            gPad->SetLogx ();

            const float ymin = (iDType == 0 ? -maxDataSyst : -maxMCSyst);
            const float ymax = (iDType == 0 ?  maxDataSyst :  maxMCSyst);

            h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#delta N_{ch} / N_{ch} [%]", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
            h->GetXaxis ()->SetMoreLogLabels ();
            h->GetYaxis ()->SetRangeUser (ymin, ymax);

            h->SetLineWidth (1);
            h->SetLineStyle (2);
            h->DrawCopy ("hist ][");
            SaferDelete (&h);

            {
              const TString totVar = totalVariations[iSysType];
              g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_ref_unf_systTot[iPtJInt][iDir][iSysType], true, 100);
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

            for (short iVar = 1; iVar < nVar; iVar++) {
              const TString var = variations[iVar];
              if ((iSysType == 0  && !IsTrackingVariation (var)) || (iSysType == 1 && !IsJetsVariation (var)) || (iSysType == 2 && !IsMixingVariation (var)))
                continue;

              g = (TGAE*) g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar]->Clone ("gtemp");
              SaveRelativeErrors (g, g, g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][0], 100);
              ResetXErrors (g);
              ResetTGAEErrors (g);
              myDraw (g, varStyles[var].first, kDot, 0, varStyles[var].second, 2, "L");
              SaferDelete (&g);
            }

            myText (0.25, 0.86, kBlack, "#bf{#it{pp}}", 0.06);
          }


          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            TH1D* h = nullptr;
            TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

            c->cd (nZdcCentBins+1-iCent);
            gPad->SetLogx ();

            const float ymin = (iDType == 0 ? -maxDataSyst : -maxMCSyst);
            const float ymax = (iDType == 0 ?  maxDataSyst :  maxMCSyst);

            h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#delta N_{ch} / N_{ch} [%]", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
            h->GetXaxis ()->SetMoreLogLabels ();
            h->GetYaxis ()->SetRangeUser (ymin, ymax);

            h->SetLineWidth (1);
            h->SetLineStyle (2);
            h->DrawCopy ("hist ][");
            SaferDelete (&h);

            {
              const TString totVar = totalVariations[iSysType];
              g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_unf_systTot[iPtJInt][iDir][iCent][iSysType], true, 100);
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

            for (short iVar = 1; iVar < nVar; iVar++) {
              const TString var = variations[iVar];
              if ((iSysType == 0  && !IsTrackingVariation (var)) || (iSysType == 1 && !IsJetsVariation (var)) || (iSysType == 2 && !IsMixingVariation (var)))
                continue;

              g = (TGAE*) g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar]->Clone ("gtemp");
              SaveRelativeErrors (g, g, g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][0], 100);
              ResetXErrors (g);
              ResetTGAEErrors (g);
              myDraw (g, varStyles[var].first, kDot, 0, varStyles[var].second, 2, "L");
              SaferDelete (&g);
            }

            if (iCent < nZdcCentBins)
              myText (0.25, 0.86, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.06);
            else
              myText (0.25, 0.86, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.06);

          } // end loop over iCent

          c->cd (8);
          myText (0.00, 0.95, kBlack, "#bf{#it{ATLAS}} Internal", 0.06);
          myText (0.00, 0.87, kBlack, "#it{pp} & #it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.06);
          myText (0.00, 0.79, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, #Delta#phi_{ch,jet} %s", iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "< #pi/8" : "> 7#pi/8"), 0.06);

          short count = 0;
          for (short iVar = 1; iVar < nVar; iVar++) {
            const TString var = variations[iVar];
            if ((iSysType == 0  && !IsTrackingVariation (var)) || (iSysType == 1 && !IsJetsVariation (var)) || (iSysType == 2 && !IsMixingVariation (var)))
              continue;
            myLineColorText (0.00+(iDType == 1 ? 0.60 : 0.50)*(count/maxNSys), 0.72-(count%maxNSys)*0.060, varStyles[var].first, varStyles[var].second, varFullNames[var], 2.0, 0.04);
            count++;
          }
          {
            const TString totVar = totalVariations[iSysType];
            myLineColorText (0.00+(iDType == 1 ? 0.60 : 0.50)*(count/maxNSys), 0.72-(count%maxNSys)*0.060, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 2.0, 0.04);
            count++;
          }

          c->SaveAs (Form ("%s/Plots/Systematics/PtCh/SignalJetTaggedYield_%s_ptch_%iGeVJets_%s_syst.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60, sysType.Data ()));

        } // end loop over iPtJInt
      }



      if (makeIpPbSystPlots) {
        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          for (short iCent = 0; iCent < nZdcCentBins; iCent++) {
  
            const char* canvasName = Form ("c_jetInt_trk_pt_%s_iaa_%s_pPb_iCent%i_syst_%s", directions[iDir].Data (), sysType.Data (), iCent, pTJInt.Data ());
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
  
            {
              const TString totVar = totalVariations[iSysType];
              g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_iaa_systTot[iPtJInt][iDir][iCent][iSysType], true, 100);
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
  
            for (short iVar = 1; iVar < nVar; iVar++) {
              const TString var = variations[iVar];
              if ((iSysType == 0  && !IsTrackingVariation (var)) || (iSysType == 1 && !IsJetsVariation (var)) || (iSysType == 2 && !IsMixingVariation (var)))
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
  
            short count = 0;
            for (short iVar = 1; iVar < nVar; iVar++) {
              const TString var = variations[iVar];
              if ((iSysType == 0  && !IsTrackingVariation (var)) || (iSysType == 1 && !IsJetsVariation (var)) || (iSysType == 2 && !IsMixingVariation (var)))
                continue;
              myLineColorText (x0+dx*(count/maxNSys), y0-(count%maxNSys)*dy, varStyles[var].first, varStyles[var].second, varFullNames[var], 0.5, tsize);
              count++;
            }
            {
              const TString totVar = totalVariations[iSysType];
              myLineColorText (x0+dx*(count/maxNSys), y0-(count%maxNSys)*dy, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 0.5, tsize);
              count++;
            }
  
            c->SaveAs (Form ("%s/Plots/Systematics/PtCh/JetTagged_IpPb_%i-%iperc_%s_ptch_%iGeVJets_%s_syst.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60, sysType.Data ()));
          } // end loop over iCent

        } // end loop over iPtJInt



        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          const char* canvasName = Form ("c_jetInt_trk_pt_%s_iaa_%s_pPb_syst_%s", directions[iDir].Data (), sysType.Data (), pTJInt.Data ());
          TCanvas* c = new TCanvas (canvasName, "", 2800, 1400);
          c->Divide (4, 2);

          std::vector <short> totVars (0);
          if (iDType == 0) {
            totVars.push_back (0);
            totVars.push_back (1);
          }
          else 
            totVars.push_back (2);

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
  
            TH1D* h = nullptr;
            TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;
  
            c->cd (nZdcCentBins-iCent+1); 
            gPad->SetLogx ();
  
            const float ymin = (iDType == 0 ? -maxDataSyst : -maxMCSyst);
            const float ymax = (iDType == 0 ?  maxDataSyst :  maxMCSyst);
  
            h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#delta I_{pPb} / I_{pPb} [%]", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
            h->GetXaxis ()->SetMoreLogLabels ();
            h->GetYaxis ()->SetRangeUser (ymin, ymax);
  
            h->SetLineWidth (1);
            h->SetLineStyle (2);
            h->DrawCopy ("hist ][");
            SaferDelete (&h);
  
            {
              const TString totVar = totalVariations[iSysType];
              g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_iaa_systTot[iPtJInt][iDir][iCent][iSysType], true, 100);
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
  
            for (short iVar = 1; iVar < nVar; iVar++) {
              const TString var = variations[iVar];
              if ((iSysType == 0  && !IsTrackingVariation (var)) || (iSysType == 1 && !IsJetsVariation (var)) || (iSysType == 2 && !IsMixingVariation (var)))
                continue;
  
              g = (TGAE*) g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar]->Clone ("gtemp");
              SaveRelativeErrors (g, g, g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0], 100);
              ResetXErrors (g);
              ResetTGAEErrors (g);
              myDraw (g, varStyles[var].first, kDot, 0, varStyles[var].second, 2, "L");
              SaferDelete (&g);
            }
  
            if (iCent < nZdcCentBins)
              myText (0.25, 0.86, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.06);
            else
              myText (0.25, 0.86, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.06);

          } // end loop over iCent

          c->cd (7);
          myText (0.00, 0.95, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
          myText (0.00, 0.87, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV", iPtJInt == 0 ? 30 : 60), 0.07);
          myText (0.00, 0.79, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8 (near-side)" : "> 7#pi/8 (away-side)"), 0.07);

          short count = 0;
          for (short iVar = 1; iVar < nVar; iVar++) {
            const TString var = variations[iVar];
            if ((iSysType == 0  && !IsTrackingVariation (var)) || (iSysType == 1 && !IsJetsVariation (var)) || (iSysType == 2 && !IsMixingVariation (var)))
              continue;
            float xcoord = 0.05;
            float ycoord = 0.72-count*0.07;
            if (count >= 10) {
              c->cd (8);
              ycoord += 0.85;
            }
            myLineColorText (xcoord, ycoord, varStyles[var].first, varStyles[var].second, varFullNames[var], 2.0, 0.06);
            count++;
          }
          {
            const TString totVar = totalVariations[iSysType];
            float xcoord = 0.05;
            float ycoord = 0.72-count*0.07;
            if (count >= 10) {
              c->cd (8);
              ycoord += 0.85;
            }
            myLineColorText (xcoord, ycoord, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 2.0, 0.06);
            count++;
          }

          c->SaveAs (Form ("%s/Plots/Systematics/PtCh/JetTagged_IpPb_%s_ptch_%iGeVJets_%s_syst.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60, sysType.Data ()));

        } // end loop over iPtJInt

      }

    } // end loop over iDir

  } // end loop over iDType



  {
    const short maxNSys = 6;

    for (short iDir : {0, 2}) {

      if (makeSigSystPlots) {
        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          const char* canvasName = Form ("c_jetInt_trk_pt_%s_unf_data_pp_systTot_%s", directions[iDir].Data (), pTJInt.Data ());
          TCanvas* c = new TCanvas (canvasName, "", 800, 800);

          TH1D* h = nullptr;
          TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

          c->cd ();
          gPad->SetLogx ();
          //gPad->SetLogy ();

          const float ymin = -maxDataSyst;
          const float ymax =  maxDataSyst;

          h = (TH1D*) h_jetInt_trk_pt_ref_unf[0][iPtJInt][iDir]->Clone ("h");
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

          for (short iTotVar : {0, 1, 2}) {
            const TString totVar = totalVariations[iTotVar];
            g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_ref_unf_systTot[iPtJInt][iDir][iTotVar], true, 100);
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

          g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][0], true, 100);
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

          short count = 0;
          for (short iTotVar : {0, 1, 2}) {
            const TString totVar = totalVariations[iTotVar];
            myLineColorText (0.35+0.35*(count/maxNSys), 0.35-(count%maxNSys)*0.040, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 1.0, 0.028);
            count++;
          }
          myLineColorText (0.35+0.35*(count/maxNSys), 0.35-(count%maxNSys)*0.040, kBlack, 1, "#bf{Total syst. unc.}", 1.0, 0.028);

          c->SaveAs (Form ("%s/Plots/Systematics/PtCh/SignalJetTaggedYield_pp_%s_ptch_%iGeVJets_systTot.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60));

        } // end loop over iPtJInt


        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          for (short iCent = 0; iCent < nZdcCentBins; iCent++) {

            const char* canvasName = Form ("c_jetInt_trk_pt_%s_unf_data_pPb_iCent%i_systTot_%s", directions[iDir].Data (), iCent, pTJInt.Data ());
            TCanvas* c = new TCanvas (canvasName, "", 800, 800);

            TH1D* h = nullptr;
            TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

            c->cd ();
            gPad->SetLogx ();
            //gPad->SetLogy ();

            const float ymin = -maxDataSyst;
            const float ymax =  maxDataSyst;

            h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#delta N_{ch} / N_{ch} [%]", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
            h->GetXaxis ()->SetMoreLogLabels ();
            h->GetYaxis ()->SetRangeUser (ymin, ymax);

            h->SetLineWidth (1);
            h->SetLineStyle (2);
            h->DrawCopy ("hist ][");
            SaferDelete (&h);

            for (short iTotVar : {0, 1, 2}) {
              const TString totVar = totalVariations[iTotVar];
              g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_unf_systTot[iPtJInt][iDir][iCent][iTotVar], true, 100);
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

            g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][0], true, 100);
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

            short count = 0;
            for (short iTotVar : {0, 1, 2}) {
              const TString totVar = totalVariations[iTotVar];
              myLineColorText (0.35+0.35*(count/maxNSys), 0.35-(count%maxNSys)*0.040, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 1.0, 0.030);
              count++;
            }
            myLineColorText (0.35+0.35*(count/maxNSys), 0.35-(count%maxNSys)*0.040, kBlack, 1, "#bf{Total syst. unc.}", 1.0, 0.030);

            c->SaveAs (Form ("%s/Plots/Systematics/PtCh/SignalJetTaggedYield_pPb_%i-%iperc_%s_ptch_%iGeVJets_systTot.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60));
          } // end loop over iCent

        } // end loop over iPtJInt


        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          const char* canvasName = Form ("c_jetInt_trk_pt_%s_unf_data_systTot_%s", directions[iDir].Data (), pTJInt.Data ());
          TCanvas* c = new TCanvas (canvasName, "", 2800, 1400);
          c->Divide (4, 2);
          

          {
            TH1D* h = nullptr;
            TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

            c->cd (7);
            gPad->SetLogx ();
            //gPad->SetLogy ();

            const float ymin = -maxDataSyst;
            const float ymax =  maxDataSyst;

            h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#delta N_{ch} / N_{ch} [%]", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins - (iPtJInt == 0 ? 3 : 1)]);
            h->GetXaxis ()->SetMoreLogLabels ();
            h->GetYaxis ()->SetRangeUser (ymin, ymax);

            h->SetLineWidth (1);
            h->SetLineStyle (2);
            h->DrawCopy ("hist ][");
            SaferDelete (&h);

            {
              gup = make_graph (h_jetInt_trk_pt_ref_unf[0][iPtJInt][iDir]);
              g = ConvertErrorsToCentralValues (gup, true, 100);
              SaferDelete (&gup);

              gup = (TGAE*) g->Clone ();
              gdown = (TGAE*) g->Clone ();
              FlipTGAE (gdown);

              myDraw (gup, kGray+1, kDot, 0, 2, 2, "L");
              myDraw (gdown, kGray+1, kDot, 0, 2, 2, "L");

              SaferDelete (&g);
              SaferDelete (&gup);
              SaferDelete (&gdown);
            }

            for (short iTotVar : {0, 1, 2}) {
              const TString totVar = totalVariations[iTotVar];
              g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_ref_unf_systTot[iPtJInt][iDir][iTotVar], true, 100);
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

            g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][0], true, 100);
            myDraw (g, kBlack, kDot, 0, 1, 2, "L");
            FlipTGAE (g);
            myDraw (g, kBlack, kDot, 0, 1, 2, "L");
            SaferDelete (&g);

            myText (0.25, 0.86, kBlack, "#bf{#it{pp}}", 0.06);
          }

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            TH1D* h = nullptr;
            TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

            c->cd (nZdcCentBins+1-iCent);
            gPad->SetLogx ();
            //gPad->SetLogy ();

            const float ymin = -maxDataSyst;
            const float ymax =  maxDataSyst;

            h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#delta N_{ch} / N_{ch} [%]", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins - (iPtJInt == 0 ? 3 : 1)]);
            h->GetXaxis ()->SetMoreLogLabels ();
            h->GetYaxis ()->SetRangeUser (ymin, ymax);

            h->SetLineWidth (1);
            h->SetLineStyle (2);
            h->DrawCopy ("hist ][");
            SaferDelete (&h);

            {
              gup = make_graph (h_jetInt_trk_pt_unf[0][iPtJInt][iDir][iCent]);
              g = ConvertErrorsToCentralValues (gup, true, 100);
              SaferDelete (&gup);

              gup = (TGAE*) g->Clone ();
              gdown = (TGAE*) g->Clone ();
              FlipTGAE (gdown);

              myDraw (gup, kGray+1, kDot, 0, 2, 2, "L");
              myDraw (gdown, kGray+1, kDot, 0, 2, 2, "L");

              SaferDelete (&g);
              SaferDelete (&gup);
              SaferDelete (&gdown);
            }

            for (short iTotVar : {0, 1, 2}) {
              const TString totVar = totalVariations[iTotVar];
              g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_unf_systTot[iPtJInt][iDir][iCent][iTotVar], true, 100);
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

            g = ConvertErrorsToCentralValues (g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][0], true, 100);
            myDraw (g, kBlack, kDot, 0, 1, 2, "L");
            FlipTGAE (g);
            myDraw (g, kBlack, kDot, 0, 1, 2, "L");
            SaferDelete (&g);

            if (iCent < nZdcCentBins)
              myText (0.25, 0.86, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.06);
            else
              myText (0.25, 0.86, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.06);

          } // end loop over iCent

          c->cd (8);
          myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
          myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
          myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.07);
          myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV", iPtJInt == 0 ? 30 : 60), 0.07);
          myText (0.1, 0.57, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8 (near-side)" : "> 7#pi/8 (away-side)"), 0.07);

          short count = 0;
          myLineColorText (0.16, 0.49, kGray+1, 2, "#bf{Total stat. unc.}", 2.0, 0.07);
          count++;
          for (short iTotVar : {0, 1, 2}) {
            const TString totVar = totalVariations[iTotVar];
            myLineColorText (0.16, 0.49-count*0.08, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 1.0, 0.07);
            count++;
          }
          myLineColorText (0.16, 0.49-count*0.080, kBlack, 1, "#bf{Total syst. unc.}", 1.0, 0.07);

          c->SaveAs (Form ("%s/Plots/Systematics/PtCh/SignalJetTaggedYield_%s_ptch_%iGeVJets_systTot.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60));

        } // end loop over iPtJInt
      }



      if (makeIpPbSystPlots) {
        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          for (short iCent = 0; iCent < nZdcCentBins; iCent++) {

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

            for (short iTotVar : {0, 1, 2}) {
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

            short count = 0;
            for (short iTotVar : {0, 1, 2}) {
              const TString totVar = totalVariations[iTotVar];
              myLineColorText (0.35+0.35*(count/maxNSys), 0.35-(count%maxNSys)*0.040, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 1.0, 0.030);
              count++;
            }
            myLineColorText (0.35+0.35*(count/maxNSys), 0.35-(count%maxNSys)*0.040, kBlack, 1, "#bf{Total syst. unc.}", 1.0, 0.028);

            c->SaveAs (Form ("%s/Plots/Systematics/PtCh/JetTagged_IpPb_%i-%iperc_%s_ptch_%iGeVJets_systTot.pdf", workPath.Data (), zdcCentPercs[iCent+1], zdcCentPercs[iCent], directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60));
          } // end loop over iCent

        } // end loop over iPtJInt



        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          const char* canvasName = Form ("c_jetInt_trk_pt_%s_iaa_data_systTot_%s", directions[iDir].Data (), pTJInt.Data ());
          TCanvas* c = new TCanvas (canvasName, "", 2800, 1400);
          c->Divide (4, 2);

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            TH1D* h = nullptr;
            TGAE* g = nullptr, *gup = nullptr, *gdown = nullptr;

            c->cd (nZdcCentBins-iCent+1); 
            gPad->SetLogx ();

            const float ymin = -maxDataSyst;
            const float ymax =  maxDataSyst;

            h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];#delta I_{pPb} / I_{pPb} [%]", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
            h->GetXaxis ()->SetMoreLogLabels ();
            h->GetYaxis ()->SetRangeUser (ymin, ymax);

            h->SetLineWidth (1);
            h->SetLineStyle (2);
            h->DrawCopy ("hist ][");
            SaferDelete (&h);

            {
              gup = make_graph (h_jetInt_trk_pt_iaa[0][iPtJInt][iDir][iCent]);
              g = ConvertErrorsToCentralValues (gup, true, 100);
              SaferDelete (&gup);

              gup = (TGAE*) g->Clone ();
              gdown = (TGAE*) g->Clone ();
              FlipTGAE (gdown);

              myDraw (gup, kGray+1, kDot, 0, 2, 2, "L");
              myDraw (gdown, kGray+1, kDot, 0, 2, 2, "L");

              SaferDelete (&g);
              SaferDelete (&gup);
              SaferDelete (&gdown);
            }

            for (short iTotVar : {0, 1, 2}) {
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

            if (iCent < nZdcCentBins)
              myText (0.25, 0.86, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.06);
            else
              myText (0.25, 0.86, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.06);

          } // end loop over iCent

          c->cd (7);
          myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
          myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
          myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.07);
          myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV", iPtJInt == 0 ? 30 : 60), 0.07);
          myText (0.1, 0.57, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8 (near-side)" : "> 7#pi/8 (away-side)"), 0.07);

          short count = 0;
          myLineColorText (0.16, 0.49, kGray+1, 2, "#bf{Total stat. unc.}", 2.0, 0.07);
          count++;
          for (short iTotVar : {0, 1, 2}) {
            const TString totVar = totalVariations[iTotVar];
            myLineColorText (0.16, 0.49-count*0.08, varStyles[totVar].first, varStyles[totVar].second, varFullNames[totVar], 1.0, 0.07);
            count++;
          }
          myLineColorText (0.16, 0.49-count*0.080, kBlack, 1, "#bf{Total syst. unc.}", 1.0, 0.07);

          c->SaveAs (Form ("%s/Plots/Systematics/PtCh/JetTagged_IpPb_%s_ptch_%iGeVJets_systTot.pdf", workPath.Data (), directions[iDir] == "ns" ? "nearside" : "awayside", iPtJInt == 0 ? 30 : 60));

        } // end loop over iPtJInt

      }

    } // end loop over iDir

  }

}


#endif
