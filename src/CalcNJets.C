#ifndef __JetHadronCorrelator_CalcNJets_C__
#define __JetHadronCorrelator_CalcNJets_C__

#include "CentralityDefs.h"
#include "Params.h"
#include "OutTree.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Variations.h"

#include <ArrayTemplates.h>
#include <Utilities.h>

#include <TH1D.h>
#include <TH2D.h>

#include <iostream>
#include <string.h>
#include <math.h>

using namespace JetHadronCorrelations;


void CalcNJets () {

  TFile* inFile = nullptr;

  TH1D**    h_evt_counts_ref              = Get1DArray <TH1D*> (3);
  TH1D***   h_jet_counts_ref              = Get2DArray <TH1D*> (3, nPtJBins);
  TH1D***   h_evt_counts                  = Get2DArray <TH1D*> (3, nZdcCentBins+1);
  TH1D****  h_jet_counts                  = Get3DArray <TH1D*> (3, nPtJBins, nZdcCentBins+1);

  TH1D**    h_jet_pt_ref                  = Get1DArray <TH1D*> (3);
  TH2D**    h2_jet_pt_cov_ref             = Get1DArray <TH2D*> (3);
  TH1D***   h_jet_pt                      = Get2DArray <TH1D*> (3, nZdcCentBins+1);
  TH2D***   h2_jet_pt_cov                 = Get2DArray <TH2D*> (3, nZdcCentBins+1);

  TDirectory* dir = gDirectory;

  for (short iDType = 0; iDType < 3; iDType++) {

    const TString dType = (iDType < 2 ? "data" : "mc");
    const TString tType = (iDType == 0 ? "mb" : (iDType == 1 ? "j50" : "JZ0123"));

    TH2D* h2_cov = nullptr;

    {
      TString inFileName = Form ("%s/JetPtWeights/Nominal/%s17_5TeV_%s.root", rootPath.Data (), dType.Data (), tType.Data ());
      std::cout << "Reading " << inFileName.Data () << std::endl;
      inFile = new TFile (inFileName.Data (), "read");
      dir->cd ();

      h_evt_counts_ref[iDType]  = (TH1D*) inFile->Get ("h_evt_counts")->Clone   (Form ("h_evt_counts_ref_%s",   tType.Data ()));

      h_jet_pt_ref[iDType]      = (TH1D*) inFile->Get ("h_jet_pt")->Clone       (Form ("h_jet_pt_ref_%s",       tType.Data ()));
      h2_cov                    = (TH2D*) inFile->Get ("h2_jet_pt_cov")->Clone  (Form ("h2_jet_pt_cov_ref_%s",  tType.Data ()));

      CalcUncertainties (h_jet_pt_ref[iDType], h2_cov, h_evt_counts_ref[iDType]);

      SaferDelete (&h2_cov);

      h_jet_pt_ref[iDType]->Scale (h_evt_counts_ref[iDType]->GetBinContent (2)); // convert distribution to total number of jets by un-scaling 1/N_evt factor

      SaferDelete (&h_evt_counts_ref[iDType]);


      for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

        h_jet_counts_ref[iDType][iPtJ] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s", pTJ.Data ()))->Clone (Form ("h_jet_counts_ref_%s_%s", tType.Data (), pTJ.Data ()));

      } // end loop over iPtJ

      inFile->Close ();
      SaferDelete (&inFile);
    }


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

      TString inFileName = Form ("%s/JetPtWeights/Nominal/%s16_5TeV_%s_%s.root", rootPath.Data (), dType.Data (), tType.Data (), cent.Data ());
      std::cout << "Reading " << inFileName.Data () << std::endl;
      inFile = new TFile (inFileName.Data (), "read");
      dir->cd ();

      h_evt_counts[iDType][iCent] = (TH1D*) inFile->Get ("h_evt_counts")->Clone   (Form ("h_evt_counts_pPb_%s_%s",  cent.Data (), tType.Data ()));

      h_jet_pt[iDType][iCent]     = (TH1D*) inFile->Get ("h_jet_pt")->Clone       (Form ("h_jet_pt_pPb_%s_%s",      cent.Data (), tType.Data ()));
      h2_cov                      = (TH2D*) inFile->Get ("h2_jet_pt_cov")->Clone  (Form ("h2_jet_pt_cov_pPb_%s_%s", cent.Data (), tType.Data ()));

      CalcUncertainties (h_jet_pt[iDType][iCent], h2_cov, h_evt_counts[iDType][iCent]);

      SaferDelete (&h2_cov);

      h_jet_pt[iDType][iCent]->Scale (h_evt_counts[iDType][iCent]->GetBinContent (2)); // convert distribution to total number of jets by un-scaling 1/N_evt factor

      SaferDelete (&h_evt_counts[iDType][iCent]);

      for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

        h_jet_counts[iDType][iPtJ][iCent] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s",  pTJ.Data ()))->Clone (Form ("h_jet_counts_pPb_%s_%s_%s", cent.Data (), tType.Data (), pTJ.Data ()));

      } // end loop over iPtJ

      inFile->Close ();
      SaferDelete (&inFile);
    } // end loop over iCent

  } // end loop over iDType




  float* njet = new float[nZdcCentBins+2];
  {
    const float trigpt = 30;
    const float maxpt = 60;

    std::cout << "---------------" << std::endl << "JETS MEETING MB > " << trigpt << " GeV" << std::endl << "---------------" << std::endl;

    njet[0] = 0;
    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
      if (0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) >= trigpt && 0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) < maxpt)
        njet[0] += h_jet_counts_ref[0][iPtJ]->GetBinContent (1);
    } // end loop over iPtJ
    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
      njet[iCent+1] = 0;
      for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
        if (0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) >= trigpt && 0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) < maxpt)
          njet[iCent+1] += h_jet_counts[0][iPtJ][iCent]->GetBinContent (1);
      } // end loop over iPtJ
    } // end loop over iCent

    std::cout << "Number of pp jets: " << njet[0] << std::endl;
    for (short iCent = 0; iCent < nZdcCentBins; iCent++)
      std::cout << "Number of p+Pb " << zdcCentPercs[iCent+1] << "\%-" << zdcCentPercs[iCent] << "\% jets: " << njet[iCent+1] << std::endl;
    std::cout << "Number of p+Pb 0-100\% jets: " << njet[nZdcCentBins] << std::endl;

    std::cout << "Formatted for latex:" << std::endl;
    std::cout << (int) njet[0];
    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++)
      std::cout << " & " << (int) njet[iCent+1];
    std::cout << std::endl << std::endl;
  }



  {
    const float trigpt = 60;
    const float maxpt = 300;

    std::cout << "---------------" << std::endl << "JETS MEETING J50" << std::endl << "---------------" << std::endl;

    njet[0] = 0;
    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
      if (0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) >= trigpt && 0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) < maxpt)
        njet[0] += h_jet_counts_ref[1][iPtJ]->GetBinContent (1);
    } // end loop over iPtJ
    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
      njet[iCent+1] = 0;
      for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
        if (0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) >= trigpt && 0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) < maxpt)
          njet[iCent+1] += h_jet_counts[1][iPtJ][iCent]->GetBinContent (1);
      } // end loop over iPtJ
    } // end loop over iCent

    std::cout << "Number of pp jets: " << njet[0] << std::endl;
    for (short iCent = 0; iCent < nZdcCentBins; iCent++)
      std::cout << "Number of p+Pb " << zdcCentPercs[iCent+1] << "\%-" << zdcCentPercs[iCent] << "\% jets: " << njet[iCent+1] << std::endl;
    std::cout << "Number of p+Pb 0-100\% jets: " << njet[nZdcCentBins] << std::endl;

    std::cout << "Formatted for latex:" << std::endl;
    std::cout << (int) njet[0];
    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++)
      std::cout << " & " << (int) njet[iCent+1];
    std::cout << std::endl << std::endl;
  }


  {
    const float trigpt = 60;
    const float maxpt = 300;

    std::cout << "---------------" << std::endl << "JETS IN MC > " << trigpt << " GeV" << std::endl << "---------------" << std::endl;
        
    njet[0] = 0;
    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
      if (0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) >= trigpt && 0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) < maxpt)
        njet[0] += h_jet_counts_ref[2][iPtJ]->GetBinContent (1);
    } // end loop over iPtJ
    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
      njet[iCent+1] = 0;
      for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
        if (0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) >= trigpt && 0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) < maxpt)
          njet[iCent+1] += h_jet_counts[2][iPtJ][iCent]->GetBinContent (1);
      } // end loop over iPtJ
    } // end loop over iCent
    
    std::cout << "Number of pp jets: " << njet[0] << std::endl;
    for (short iCent = 0; iCent < nZdcCentBins; iCent++)
      std::cout << "Number of p+Pb " << zdcCentPercs[iCent+1] << "\%-" << zdcCentPercs[iCent] << "\% jets: " << njet[iCent+1] << std::endl;
    std::cout << "Number of p+Pb 0-100\% jets: " << njet[nZdcCentBins] << std::endl;
    
    std::cout << "Formatted for latex:" << std::endl;
    std::cout << (int) njet[0];
    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++)
      std::cout << " & " << (int) njet[iCent+1];
    std::cout << std::endl;
    

  }
  delete[] njet;


}


#endif
