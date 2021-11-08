#ifndef __JetHadronCorrelator_NjetNcoll_C__
#define __JetHadronCorrelator_NjetNcoll_C__

#include <iostream>
#include <math.h>

#include <TColor.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TChain.h>
#include <TTree.h>

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


// NOTE: DEFAULT AREA UNITS ARE µb THROUGHOUT!!!

void NjetNcoll () {

  TFile* outFile = new TFile (Form ("%s/NjetNcoll.root", rootPath.Data ()), "recreate");

  float zdcE = 0, fcalE = 0;
  int jet_n = 0;
  float* jet_pt = new float[40];
  float* jet_eta = new float[40];
  bool jetTrig = false;

  const float sigmapp = 67600; // in µb

  const int nCentBins = nZdcCentBins;
  //const int nCentBins = 7;

  //const double fracCent[7] = {0.3, 0.2, 0.1, 0.1, 0.1, 0.05, 0.05};
  const double fracCent[5] = {0.2, 0.2, 0.2, 0.2, 0.2};

  double** jet_counts = Get2DArray <double> (2, nCentBins+1);
  double n_mb = 56759796. / 25.1334; // total number of MB events per 1/µb (units: µb)

  const float yboost = -0.465;

  for (short iPtJInt : {0, 1}) {
    for (short iCent = 0; iCent < nCentBins+1; iCent++) {
      jet_counts[iPtJInt][iCent] = 0.;
    }
  }


  TChain* pPbTree = new TChain ("pPbTree", "pPbTree");
  for (short iPtJInt : {0, 1}) {

    const float minJetPt = (iPtJInt == 0 ? 30 : 60);

    {
      const TString inFilePattern = Form ("%s/Trees/%s/data16_5TeV_iCent*.root", rootPath.Data (), iPtJInt == 0 ? "MinBias" : "J50");
      std::cout << "Reading " << inFilePattern.Data () << std::endl;
      pPbTree->Reset ();
      pPbTree->Add (Form ("%s/Trees/%s/data16_5TeV_iCent0.root", rootPath.Data (), iPtJInt == 0 ? "MinBias" : "J50"));
      pPbTree->Add (Form ("%s/Trees/%s/data16_5TeV_iCent1.root", rootPath.Data (), iPtJInt == 0 ? "MinBias" : "J50"));
      pPbTree->Add (Form ("%s/Trees/%s/data16_5TeV_iCent2.root", rootPath.Data (), iPtJInt == 0 ? "MinBias" : "J50"));
      pPbTree->Add (Form ("%s/Trees/%s/data16_5TeV_iCent3.root", rootPath.Data (), iPtJInt == 0 ? "MinBias" : "J50"));
      pPbTree->Add (Form ("%s/Trees/%s/data16_5TeV_iCent4.root", rootPath.Data (), iPtJInt == 0 ? "MinBias" : "J50"));

      pPbTree->SetBranchAddress ("fcal_et_Pb",         &fcalE);
      pPbTree->SetBranchAddress ("zdc_calibE_Pb",      &zdcE);
      pPbTree->SetBranchAddress ("akt4_hi_jet_n",      &jet_n);
      pPbTree->SetBranchAddress ("akt4_hi_jet_pt",     jet_pt);
      pPbTree->SetBranchAddress ("akt4_hi_jet_eta",    jet_eta);

      if (iPtJInt == 1)
        pPbTree->SetBranchAddress ("jetTrig",          &jetTrig);
      else
        jetTrig = true;

      const long nEvts = pPbTree->GetEntries ();
      for (long iEvt = 0; iEvt < nEvts; iEvt++) {
        if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
          std::cout << "Events " << iEvt / (nEvts / 100) << "\% done...\r" << flush;

        pPbTree->GetEntry (iEvt);

        if (zdcE == 0)
          continue;

        if (iPtJInt == 1 && !jetTrig)
          continue;

        short iCent = 0;
        if (zdcE < 21.3314)
          iCent = 0;
        else if (zdcE < 40.1944)
          iCent = 1;
        else if (zdcE < 54.5756)
          iCent = 2;
        else if (zdcE < 67.1273)
          iCent = 3;
        else
          iCent = 4;

        //short iCent = 0;
        //if (fcalE < 2.55)
        //  continue;
        //else if (fcalE < 13.40)
        //  iCent = 0; // 60-90%
        //else if (fcalE < 24.07)
        //  iCent = 1; // 40-60%
        //else if (fcalE < 31.03)
        //  iCent = 2; // 30-40%
        //else if (fcalE < 40.01)
        //  iCent = 3; // 20-30%
        //else if (fcalE < 53.69)
        //  iCent = 4; // 10-20%
        //else if (fcalE < 65.95)
        //  iCent = 5; // 5-10%
        //else
        //  iCent = 6; // 0-5%

        for (short iJet = 0; iJet < jet_n; iJet++) {
          if (jet_pt[iJet] > minJetPt)// && std::fabs (jet_eta[iJet] - yboost) < 2.0)
            jet_counts[iPtJInt][iCent] = jet_counts[iPtJInt][iCent] + 1;
        }
      }
      std::cout << std::endl;
    }


    {
      const TString inFileName = Form ("rootFiles/Trees/%s/data17_5TeV.root", iPtJInt == 0 ? "MinBias" : "J50");
      std::cout << "Reading " << inFileName.Data () << std::endl;
      TFile* inFile = new TFile (inFileName.Data (), "read");
      TTree* ppTree = (TTree*) inFile->Get ("ppTree");

      ppTree->SetBranchAddress ("akt4_hi_jet_n",      &jet_n);
      ppTree->SetBranchAddress ("akt4_hi_jet_pt",     jet_pt);
      ppTree->SetBranchAddress ("akt4_hi_jet_eta",    jet_eta);

      jetTrig = true;

      const long nEvts = ppTree->GetEntries ();
      for (long iEvt = 0; iEvt < nEvts; iEvt++) {
        if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
          std::cout << "Events " << iEvt / (nEvts / 100) << "\% done...\r" << flush;

        ppTree->GetEntry (iEvt);

        if (iPtJInt == 1 && !jetTrig)
          continue;

        for (short iJet = 0; iJet < jet_n; iJet++) {
          if (jet_pt[iJet] > minJetPt)// && std::fabs (jet_eta[iJet]) < 2.0)
            jet_counts[iPtJInt][nCentBins] = jet_counts[iPtJInt][nCentBins] + 1;
        }
      }
      inFile->Close ();
      std::cout << std::endl;
    }
  }


  TGAE** g_njet_ncoll = new TGAE*[2];
  TGAE** g_njet_ncoll_syst = new TGAE*[2];

  TGAE** g_njetoverncoll_ncoll = new TGAE*[2];
  TGAE** g_njetoverncoll_ncoll_syst = new TGAE*[2];

  for (short iPtJInt : {0, 1}) {

    g_njet_ncoll[iPtJInt] = new TGAE ();
    g_njet_ncoll_syst[iPtJInt] = new TGAE ();

    g_njetoverncoll_ncoll[iPtJInt] = new TGAE ();
    g_njetoverncoll_ncoll_syst[iPtJInt] = new TGAE ();

    double njetlumipp = 1;

    {
      const double njet = jet_counts[iPtJInt][nCentBins];
      const double ncoll = 1;
      const double sigma_ncoll = 0;
      //const double lumi = (iPtJInt == 0 ? 2.65513e3 : 3.57487e6) * 0.244; // multiply by Poisson probability of 1 with lambda=2.2 (average pileup rate for run 340718).
      const double lumi = (iPtJInt == 0 ? 2.65513e3 : 3.57487e6);

      njetlumipp = njet / lumi;

      g_njet_ncoll[iPtJInt]->SetPoint (0, ncoll, 1);
      g_njet_ncoll[iPtJInt]->SetPointEXhigh (0, 0);
      g_njet_ncoll[iPtJInt]->SetPointEXlow (0, 0);
      g_njet_ncoll[iPtJInt]->SetPointEYhigh (0, 0);
      g_njet_ncoll[iPtJInt]->SetPointEYlow (0, 0);

      g_njet_ncoll_syst[iPtJInt]->SetPoint (0, ncoll, 1);
      g_njet_ncoll_syst[iPtJInt]->SetPointEXhigh (0, 0);
      g_njet_ncoll_syst[iPtJInt]->SetPointEXlow (0, 0);
      g_njet_ncoll_syst[iPtJInt]->SetPointEYhigh (0, 0);
      g_njet_ncoll_syst[iPtJInt]->SetPointEYlow (0, 0);

      g_njetoverncoll_ncoll[iPtJInt]->SetPoint (0, ncoll, 1);//njet / (lumi * ncoll));
      g_njetoverncoll_ncoll[iPtJInt]->SetPointEXhigh (0, 0);
      g_njetoverncoll_ncoll[iPtJInt]->SetPointEXlow (0, 0);
      g_njetoverncoll_ncoll[iPtJInt]->SetPointEYhigh (0, 0);
      g_njetoverncoll_ncoll[iPtJInt]->SetPointEYlow (0, 0);

      g_njetoverncoll_ncoll_syst[iPtJInt]->SetPoint (0, ncoll, 1);
      g_njetoverncoll_ncoll_syst[iPtJInt]->SetPointEXhigh (0, 0);
      g_njetoverncoll_ncoll_syst[iPtJInt]->SetPointEXlow (0, 0);
      g_njetoverncoll_ncoll_syst[iPtJInt]->SetPointEYhigh (0, 0);
      g_njetoverncoll_ncoll_syst[iPtJInt]->SetPointEYlow (0, 0);

    }

    for (int iCent = 0; iCent < nCentBins; iCent++) {

      const double njet = jet_counts[iPtJInt][iCent];
      const double ncoll = zdcNcollValues[iCent];
      //const double ncoll = fcalNcollValues[iCent];
      const double taa = ncoll / sigmapp;
      const double sigma_ncoll = zdcNcollErrors[iCent];
      //const double sigma_ncoll = fcalNcollErrors[iCent];
      const double sigma_taa = sigma_ncoll / sigmapp;
      const double lumi = (iPtJInt == 0 ? 25.1334 : 3.55545e2);

      g_njet_ncoll[iPtJInt]->SetPoint (iCent+1, ncoll, (njet * ncoll) / (njetlumipp * lumi * taa * n_mb * fracCent[iCent]));
      g_njet_ncoll[iPtJInt]->SetPointEXhigh (iCent+1, 0);
      g_njet_ncoll[iPtJInt]->SetPointEXlow (iCent+1, 0);
      g_njet_ncoll[iPtJInt]->SetPointEYhigh (iCent+1, (std::sqrt (njet) * ncoll) / (njetlumipp * lumi * taa * n_mb * fracCent[iCent]));
      g_njet_ncoll[iPtJInt]->SetPointEYlow (iCent+1, (std::sqrt (njet) * ncoll) / (njetlumipp * lumi * taa * n_mb * fracCent[iCent]));

      g_njet_ncoll_syst[iPtJInt]->SetPoint (iCent+1, ncoll, (njet * ncoll) / (njetlumipp * lumi * taa * n_mb * fracCent[iCent]));
      g_njet_ncoll_syst[iPtJInt]->SetPointEXhigh (iCent+1, sigma_ncoll);
      g_njet_ncoll_syst[iPtJInt]->SetPointEXlow (iCent+1, sigma_ncoll);
      g_njet_ncoll_syst[iPtJInt]->SetPointEYhigh (iCent+1, 0.4);
      g_njet_ncoll_syst[iPtJInt]->SetPointEYlow (iCent+1, 0.4);

      g_njetoverncoll_ncoll[iPtJInt]->SetPoint (iCent+1, ncoll, njet / (njetlumipp * lumi * taa * n_mb * fracCent[iCent]));
      g_njetoverncoll_ncoll[iPtJInt]->SetPointEXhigh (iCent+1, 0);
      g_njetoverncoll_ncoll[iPtJInt]->SetPointEXlow (iCent+1, 0);
      g_njetoverncoll_ncoll[iPtJInt]->SetPointEYhigh (iCent+1, std::sqrt (njet) / (njetlumipp * lumi * taa * n_mb * fracCent[iCent]));
      g_njetoverncoll_ncoll[iPtJInt]->SetPointEYlow (iCent+1, std::sqrt (njet) / (njetlumipp * lumi * taa * n_mb * fracCent[iCent]));

      g_njetoverncoll_ncoll_syst[iPtJInt]->SetPoint (iCent, ncoll, njet / (njetlumipp * lumi * taa * n_mb * fracCent[iCent]));
      g_njetoverncoll_ncoll_syst[iPtJInt]->SetPointEXhigh (iCent, sigma_ncoll);
      g_njetoverncoll_ncoll_syst[iPtJInt]->SetPointEXlow (iCent, sigma_ncoll);
      g_njetoverncoll_ncoll_syst[iPtJInt]->SetPointEYhigh (iCent, std::fabs (njet * sigma_taa / (njetlumipp * lumi * taa * taa * n_mb * fracCent[iCent])));
      g_njetoverncoll_ncoll_syst[iPtJInt]->SetPointEYlow (iCent, std::fabs (njet * sigma_taa / (njetlumipp * lumi * taa * taa * n_mb * fracCent[iCent])));

    } // end loop over iCent

  }


  outFile->cd ();
  for (short iPtJInt : {0, 1}) {
    g_njet_ncoll[iPtJInt]->Write (Form ("g_njet_ncoll_%s", iPtJInt == 0 ? "30GeVJets" : "60GeVJets"));
    g_njet_ncoll_syst[iPtJInt]->Write (Form ("g_njet_ncoll_syst_%s", iPtJInt == 0 ? "30GeVJets" : "60GeVJets"));
    g_njetoverncoll_ncoll[iPtJInt]->Write (Form ("g_njetoverncoll_ncoll_%s", iPtJInt == 0 ? "30GeVJets" : "60GeVJets"));
    g_njetoverncoll_ncoll_syst[iPtJInt]->Write (Form ("g_njetoverncoll_ncoll_syst_%s", iPtJInt == 0 ? "30GeVJets" : "60GeVJets"));
  }
  outFile->Close ();

  Delete2DArray (&jet_counts, 2, nCentBins+1);
}


#endif
