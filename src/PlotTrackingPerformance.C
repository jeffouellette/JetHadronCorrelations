#ifndef __PlotTrackingPerformance_C__
#define __PlotTrackingPerformance_C__

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLine.h>
#include <TCanvas.h>

#include <ArrayTemplates.h>
#include <Utilities.h>
#include <MyStyle.h>
#include <MyColors.h>

#include "CentralityDefs.h"
#include "TreeVariables.h"
#include "Params.h"
#include "LocalUtilities.h"
#include "PrimaryFractionFit.h"

using namespace JetHadronCorrelations;

// Bin edge definitions: eta, pTch, pp/pPb.
const int nFinerEtaTrkBins = 40;
const double* finerEtaTrkBins = linspace (-2.5, 2.5, nFinerEtaTrkBins);

//const double pTchBins[] = {0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 50, 55, 60, 70, 80, 100}; // old binning
const double pTchBins[] = {0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130}; // new binning, 113 elements or 112 bins
const int nPtchBins = sizeof (pTchBins) / sizeof (pTchBins[0]) - 1;

                  //pi+, k+,  p+,   e-, mu-, sigma+, sigma-, xi-,  omega-, everyone
const int PIDs[] = {211, 321, 2212, 11, 13,  3222,   3112,   3312, 3334,   0};
const std::string partNames[] = {"#pi^{+}", "K^{+}", "p^{+}", "e^{-}", "#mu^{-}", "#Sigma^{+}", "#Sigma^{-}", "#Xi^{-}", "#Omega^{-}", "All h^{#pm}"};
const int nPIDs = sizeof (PIDs) / sizeof (PIDs[0]);

const vector <int> systems = {0, 1}; // 0 = pp, 1 = pPb
const int nSystems = systems.size ();


PrimaryFractionFit pff;


void InvertRate (TH1* h) {
  for (int iX = 1; iX <= h->GetNbinsX (); iX++)
    for (int iY = 1; iY <= h->GetNbinsY (); iY++)
      for (int iZ = 1; iZ <= h->GetNbinsZ (); iZ++)
        h->SetBinContent (iX, iY, iZ, 1.-h->GetBinContent (iX, iY, iZ));
  return;
}



TF1* DoPurityFit (TH1D* h, const int degree, const int nderiv) {
  std::string fname = std::string (h->GetName ());
  fname.replace (fname.begin (), fname.begin () + 1, "f");

  pff.SetNDeriv (nderiv);
  pff.SetDegree (degree);

  const int ndf = pff.NDF ();

  TF1* fit = new TF1 (fname.c_str (), &pff, 0.4, 150, ndf);

  double mean = 0, den = 0;
  for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
    mean += h->GetBinContent (ix);// * h->GetBinWidth (ix);
    den += 1;//h->GetBinWidth (ix);

    if (h->GetBinContent (ix) > 0 && h->GetBinError (ix) / h->GetBinContent (ix) < 0.005)
      h->SetBinError (ix, 0.005*h->GetBinError (ix));
  }
  mean = mean / den;

  fit->SetParameter (0, h->GetBinCenter (h->GetNbinsX ()));
  fit->SetParameter (1, mean);

  for (int i = 2; i < ndf; i++)
    fit->SetParameter (i, std::pow (0.1, i-nderiv));

  h->Fit (fit, "RN0Q");

  return fit;
}



//TF1* DoPurityFit (TH1D* h, double x0, const int nN, const int nD) {
//  std::string fname = std::string (h->GetName ());
//  fname.replace (fname.begin (), fname.begin () + 1, "f");
//
//  std::string fitParams[nN+1];
//  fitParams[0] = "[0]";
//
//  for (int n = nD+1; n <= nN; n++) {
//    fitParams[n] = "[" + std::to_string (n-nD) + "]";
//  }
//
//  for (int n = nD; n >= 1; n--) {
//
//    if (n == 1)
//      fitParams[n] = "(-(";
//    else
//      fitParams[n] = "(-1/" + std::to_string ((int) factorial (n)) + "*(";
//
//    bool needsPlus = false;
//    for (int np = n+1; np <= nN; np++) {
//      if (needsPlus)
//        fitParams[n] = fitParams[n] + "+";
//      fitParams[n] = fitParams[n] + std::to_string ((int) (factorial (np) / factorial (np-n))) + "*(log(" + std::to_string ((int) x0) + ")^" + std::to_string (np-n) + ")*";
//      fitParams[n] = fitParams[n] + "(" + fitParams[np] + ")";
//      needsPlus = true;
//    }
//
//    fitParams[n] += "))";
//  }
//
//  //for (int n = 0; n <= nN; n++) {
//  //  std::cout << "  --> param " << n << " is " << fitParams[n] << std::endl;
//  //}
//
//  std::string fitStr = fitParams[0];
//  for (int n = 1; n <= nN; n++) {
//    fitStr = fitStr + "+" + fitParams[n] + "*(log(x))^" + std::to_string (n);
//  }
//
//  //std::cout << "  --> fit str = " << fitStr << std::endl;
//
//  TF1* fit = new TF1 (fname.c_str (), fitStr.c_str (), 0.5, 100);
//
//  h->Fit (fit, "RN0Q");
//
//  return fit;
//}



void PlotTrackingPerformance () {

  TH1D*** h_truth_matching_prob = Get2DArray <TH1D*> (nSystems, 3);

  TH2D***** h2_truth_matched_primary_tracks = Get4DArray <TH2D*> (nSystems, nPIDs, trackWPs.size (), nMultBins);
  TH2D****  h2_truth_tracks                 = Get3DArray <TH2D*> (nSystems, nPIDs, nMultBins);
  TH2D****  h2_truth_tracks_wgt2            = Get3DArray <TH2D*> (nSystems, nPIDs, nMultBins);

  TH2D***** h2_efficiency                   = Get4DArray <TH2D*> (nSystems, nPIDs, trackWPs.size (), nMultBins);

  TH2D**** h2_fake_tracks               = Get3DArray <TH2D*> (nSystems, trackWPs.size (), nDRBins);
  TH2D**** h2_secondary_tracks          = Get3DArray <TH2D*> (nSystems, trackWPs.size (), nDRBins);
  TH2D**** h2_strange_tracks            = Get3DArray <TH2D*> (nSystems, trackWPs.size (), nDRBins);
  TH2D**** h2_primary_tracks            = Get3DArray <TH2D*> (nSystems, trackWPs.size (), nDRBins);
  TH2D**** h2_primary_tracks_fakes_p100 = Get3DArray <TH2D*> (nSystems, trackWPs.size (), nDRBins);
  TH2D**** h2_reco_tracks               = Get3DArray <TH2D*> (nSystems, trackWPs.size (), nDRBins);
  TH2D**** h2_reco_tracks_wgt2          = Get3DArray <TH2D*> (nSystems, trackWPs.size (), nDRBins);

  TH2D**** h2_fake_rate                 = Get3DArray <TH2D*> (nSystems, trackWPs.size (), nDRBins);
  TH2D**** h2_secondary_rate            = Get3DArray <TH2D*> (nSystems, trackWPs.size (), nDRBins);
  TH2D**** h2_strange_rate              = Get3DArray <TH2D*> (nSystems, trackWPs.size (), nDRBins);
  TH2D**** h2_primary_rate              = Get3DArray <TH2D*> (nSystems, trackWPs.size (), nDRBins);
  TH2D**** h2_primary_rate_fakes_p100   = Get3DArray <TH2D*> (nSystems, trackWPs.size (), nDRBins);


  TH1D****** h_truth_matched_primary_tracks  = Get5DArray <TH1D*> (nSystems, nPIDs, trackWPs.size (), nMultBins, nEtaTrkBins);
  TH1D*****  h_truth_tracks                  = Get4DArray <TH1D*> (nSystems, nPIDs, nMultBins, nEtaTrkBins);
  TH1D*****  h_truth_tracks_wgt2             = Get4DArray <TH1D*> (nSystems, nPIDs, nMultBins, nEtaTrkBins);

  TH1D****** h_efficiency                    = Get5DArray <TH1D*> (nSystems, nPIDs, trackWPs.size (), nMultBins, nEtaTrkBins);

  TH1D***** h_fake_tracks               = Get4DArray <TH1D*> (nSystems, trackWPs.size (), nDRBins, nEtaTrkBins);
  TH1D***** h_secondary_tracks          = Get4DArray <TH1D*> (nSystems, trackWPs.size (), nDRBins, nEtaTrkBins);
  TH1D***** h_strange_tracks            = Get4DArray <TH1D*> (nSystems, trackWPs.size (), nDRBins, nEtaTrkBins);
  TH1D***** h_primary_tracks            = Get4DArray <TH1D*> (nSystems, trackWPs.size (), nDRBins, nEtaTrkBins);
  TH1D***** h_primary_tracks_fakes_p100 = Get4DArray <TH1D*> (nSystems, trackWPs.size (), nDRBins, nEtaTrkBins);
  TH1D***** h_reco_tracks               = Get4DArray <TH1D*> (nSystems, trackWPs.size (), nDRBins, nEtaTrkBins);
  TH1D***** h_reco_tracks_wgt2          = Get4DArray <TH1D*> (nSystems, trackWPs.size (), nDRBins, nEtaTrkBins);

  TH1D***** h_fake_rate                 = Get4DArray <TH1D*> (nSystems, trackWPs.size (), nDRBins, nEtaTrkBins);
  TH1D***** h_secondary_rate            = Get4DArray <TH1D*> (nSystems, trackWPs.size (), nDRBins, nEtaTrkBins);
  TH1D***** h_strange_rate              = Get4DArray <TH1D*> (nSystems, trackWPs.size (), nDRBins, nEtaTrkBins);
  TH1D***** h_primary_rate              = Get4DArray <TH1D*> (nSystems, trackWPs.size (), nDRBins, nEtaTrkBins);
  TH1D***** h_primary_rate_fakes_p100   = Get4DArray <TH1D*> (nSystems, trackWPs.size (), nDRBins, nEtaTrkBins);

  TF1***** f_primary_rate               = Get4DArray <TF1*> (nSystems, trackWPs.size (), nDRBins, nEtaTrkBins);
  TF1***** f_primary_rate_fakes_p100    = Get4DArray <TF1*> (nSystems, trackWPs.size (), nDRBins, nEtaTrkBins);


  TFile* inFile = new TFile (Form ("%s/TrackingPerformance/Nominal/outFile.root", rootPath.Data ()), "read");


  for (int iSys : systems) {

    const TString sys = (iSys == 0 ? "pp" : "pPb");

    for (int iPID = 0; iPID < nPIDs; iPID++) {

      for (int iMult = 0; iMult < nMultBins; iMult++) {

        h2_truth_tracks[iSys][iPID][iMult] = (TH2D*) inFile->Get (Form ("h2_truth_tracks_%s_PID%i_iMult%i", sys.Data (), PIDs[iPID], iMult));
        h2_truth_tracks_wgt2[iSys][iPID][iMult] = (TH2D*) inFile->Get (Form ("h2_truth_tracks_wgt2_%s_PID%i_iMult%i", sys.Data (), PIDs[iPID], iMult));

        if (iMult == 0 || iMult == 1 || iMult == 2) {
          h2_truth_tracks[iSys][iPID][iMult]->RebinY (4);
          h2_truth_tracks_wgt2[iSys][iPID][iMult]->RebinY (4);
        }
        //else if (iMult == 1) {
        //  h2_truth_tracks[iSys][iPID][iMult]->RebinY (5);
        //  h2_truth_tracks_wgt2[iSys][iPID][iMult]->RebinY (5);
        //}

        for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {

          h_truth_tracks[iSys][iPID][iMult][iEta] = (TH1D*) inFile->Get (Form ("h_truth_tracks_%s_PID%i_iMult%i_iEta%i", sys.Data (), PIDs[iPID], iMult, iEta));
          h_truth_tracks_wgt2[iSys][iPID][iMult][iEta] = (TH1D*) inFile->Get (Form ("h_truth_tracks_wgt2_%s_PID%i_iMult%i_iEta%i", sys.Data (), PIDs[iPID], iMult, iEta));

          if (iMult == 0 || iMult == 1 || iMult == 2) {
            h_truth_tracks[iSys][iPID][iMult][iEta]->Rebin (4);
            h_truth_tracks_wgt2[iSys][iPID][iMult][iEta]->Rebin (4);
          }
          //else if (iMult == 1) {
          //  h_truth_tracks[iSys][iPID][iMult][iEta]->Rebin (2);
          //  h_truth_tracks_wgt2[iSys][iPID][iMult][iEta]->Rebin (2);
          //}

        } // end loop over iEta

      } // end loop over iMult


      for (int iWP = 0; iWP < trackWPs.size (); iWP++) {

        for (int iMult = 0; iMult < nMultBins; iMult++) {

          h2_truth_matched_primary_tracks[iSys][iPID][iWP][iMult] = (TH2D*) inFile->Get (Form ("h2_truth_matched_primary_tracks_%s_PID%i_%s_iMult%i", sys.Data (), PIDs[iPID], trackWPNames[iWP].c_str (), iMult));

          if (iMult == 0 || iMult == 1 || iMult == 2)
            h2_truth_matched_primary_tracks[iSys][iPID][iWP][iMult]->RebinY (4);
          //else if (iMult == 1)
          //  h2_truth_matched_primary_tracks[iSys][iPID][iWP][iMult]->RebinY (2);

          for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {

            h_truth_matched_primary_tracks[iSys][iPID][iWP][iMult][iEta] = (TH1D*) inFile->Get (Form ("h_truth_matched_primary_tracks_%s_PID%i_%s_iMult%i_iEta%i", sys.Data (), PIDs[iPID], trackWPNames[iWP].c_str (), iMult, iEta));

            if (iMult == 0 || iMult == 1 || iMult == 2)
              h_truth_matched_primary_tracks[iSys][iPID][iWP][iMult][iEta]->Rebin (4);
            //else if (iMult == 1)
            //  h_truth_matched_primary_tracks[iSys][iPID][iWP][iMult][iEta]->Rebin (2);

          } // end loop over iEta

        } // end loop over iMult

      } // end loop over iWP

    } // end loop over iPID


    for (int iWP = 0; iWP < trackWPs.size (); iWP++) {

      h_truth_matching_prob[iSys][iWP] = (TH1D*) inFile->Get (Form ("h_truth_matching_prob_%s_%s", sys.Data (), trackWPNames[iWP].c_str ()));

      //for (int iMult = 0; iMult < nMultBins; iMult++) {
      for (int iDR = 0; iDR < nDRBins; iDR++) {

        h2_fake_tracks[iSys][iWP][iDR] = (TH2D*) inFile->Get (Form ("h2_fake_tracks_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR));
        h2_secondary_tracks[iSys][iWP][iDR] = (TH2D*) inFile->Get (Form ("h2_secondary_tracks_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR));
        h2_strange_tracks[iSys][iWP][iDR] = (TH2D*) inFile->Get (Form ("h2_strange_tracks_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR));
        //h2_primary_tracks[iSys][iWP][iDR] = (TH2D*) inFile->Get (Form ("h2_primary_tracks_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR));

        h2_reco_tracks[iSys][iWP][iDR] = (TH2D*) inFile->Get (Form ("h2_reco_tracks_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR));
        h2_reco_tracks_wgt2[iSys][iWP][iDR] = (TH2D*) inFile->Get (Form ("h2_reco_tracks_wgt2_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR));

        h2_primary_tracks[iSys][iWP][iDR] = (TH2D*) h2_reco_tracks[iSys][iWP][iDR]->Clone (Form ("h2_primary_tracks_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR));
        h2_primary_tracks[iSys][iWP][iDR]->Add (h2_fake_tracks[iSys][iWP][iDR], -1);
        h2_primary_tracks[iSys][iWP][iDR]->Add (h2_secondary_tracks[iSys][iWP][iDR], -1);

        h2_primary_tracks_fakes_p100[iSys][iWP][iDR] = (TH2D*) h2_reco_tracks[iSys][iWP][iDR]->Clone (Form ("h2_primary_tracks_fakes_p100_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR));
        h2_primary_tracks_fakes_p100[iSys][iWP][iDR]->Add (h2_fake_tracks[iSys][iWP][iDR], -2); // 100% increase of fake rate for systematics
        h2_primary_tracks_fakes_p100[iSys][iWP][iDR]->Add (h2_secondary_tracks[iSys][iWP][iDR], -1);

        for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {

          h_fake_tracks[iSys][iWP][iDR][iEta] = (TH1D*) inFile->Get (Form ("h_fake_tracks_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta));
          h_secondary_tracks[iSys][iWP][iDR][iEta] = (TH1D*) inFile->Get (Form ("h_secondary_tracks_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta));
          h_strange_tracks[iSys][iWP][iDR][iEta] = (TH1D*) inFile->Get (Form ("h_strange_tracks_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta));
          //h_primary_tracks[iSys][iWP][iDR][iEta] = (TH1D*) inFile->Get (Form ("h_primary_tracks_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta));
          h_reco_tracks[iSys][iWP][iDR][iEta] = (TH1D*) inFile->Get (Form ("h_reco_tracks_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta));
          h_reco_tracks_wgt2[iSys][iWP][iDR][iEta] = (TH1D*) inFile->Get (Form ("h_reco_tracks_wgt2_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta));

          h_primary_tracks[iSys][iWP][iDR][iEta] = (TH1D*) h_reco_tracks[iSys][iWP][iDR][iEta]->Clone (Form ("h_primary_tracks_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta));
          h_primary_tracks[iSys][iWP][iDR][iEta]->Add (h_fake_tracks[iSys][iWP][iDR][iEta], -1);
          h_primary_tracks[iSys][iWP][iDR][iEta]->Add (h_secondary_tracks[iSys][iWP][iDR][iEta], -1);

          h_primary_tracks_fakes_p100[iSys][iWP][iDR][iEta] = (TH1D*) h_reco_tracks[iSys][iWP][iDR][iEta]->Clone (Form ("h_primary_tracks_fakes_p100_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta));
          h_primary_tracks_fakes_p100[iSys][iWP][iDR][iEta]->Add (h_fake_tracks[iSys][iWP][iDR][iEta], -2); // 100% increase of fake rate for systematics
          h_primary_tracks_fakes_p100[iSys][iWP][iDR][iEta]->Add (h_secondary_tracks[iSys][iWP][iDR][iEta], -1);

          if (iSys == 0) {
            h_fake_tracks[iSys][iWP][iDR][iEta]->Rebin (2);
            h_secondary_tracks[iSys][iWP][iDR][iEta]->Rebin (2);
            h_strange_tracks[iSys][iWP][iDR][iEta]->Rebin (2);
            h_primary_tracks[iSys][iWP][iDR][iEta]->Rebin (2);
            h_primary_tracks_fakes_p100[iSys][iWP][iDR][iEta]->Rebin (2);
            h_reco_tracks[iSys][iWP][iDR][iEta]->Rebin (2);
            h_reco_tracks_wgt2[iSys][iWP][iDR][iEta]->Rebin (2);
          }

        } // end loop over iEta

      } // end loop over iDR

    } // end loop over iWP

  } // end loop over iSys


  {
    //double coarserPtchBins[] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.5, 5, 5.5, 6, 7, 8, 10, 14, 20, 40, 100};
    //const double pTchBins[] = {0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 50, 55, 60, 70, 80, 100};
    double coarserPtchBins[] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 4, 4.5, 5, 6, 8, 12, 150};
    int nCoarserPtchBins = sizeof (coarserPtchBins) / sizeof (coarserPtchBins[0]) - 1;

    for (int iWP = 0; iWP < trackWPs.size (); iWP++) {

      //for (int iMult = 0; iMult < nMultBins; iMult++) {
      for (int iDR = 0; iDR < nDRBins; iDR++) {

        for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {

          RebinSomeBins (&(h_fake_tracks[1][iWP][iDR][iEta]), nCoarserPtchBins, (double*) coarserPtchBins);
          RebinSomeBins (&(h_secondary_tracks[1][iWP][iDR][iEta]), nCoarserPtchBins, (double*) coarserPtchBins);
          RebinSomeBins (&(h_strange_tracks[1][iWP][iDR][iEta]), nCoarserPtchBins, (double*) coarserPtchBins);
          RebinSomeBins (&(h_primary_tracks[1][iWP][iDR][iEta]), nCoarserPtchBins, (double*) coarserPtchBins);
          RebinSomeBins (&(h_primary_tracks_fakes_p100[1][iWP][iDR][iEta]), nCoarserPtchBins, (double*) coarserPtchBins);
          RebinSomeBins (&(h_reco_tracks[1][iWP][iDR][iEta]), nCoarserPtchBins, (double*) coarserPtchBins);
          RebinSomeBins (&(h_reco_tracks_wgt2[1][iWP][iDR][iEta]), nCoarserPtchBins, (double*) coarserPtchBins);

        } // end loop over iEta

      } // end loop over iDR

    } // end loop over iWP
  }


  for (int iSys : systems) {

    const TString sys = (iSys == 0 ? "pp" : "pPb");

    for (int iPID = 0; iPID < nPIDs; iPID++) {

      for (int iWP = 0; iWP < trackWPs.size (); iWP++) {

        for (int iMult = 0; iMult < nMultBins; iMult++) {

          h2_efficiency[iSys][iPID][iWP][iMult] = (TH2D*) h2_truth_matched_primary_tracks[iSys][iPID][iWP][iMult]->Clone (Form ("h2_efficiency_%s_PID%i_%s_iMult%i", sys.Data (), PIDs[iPID], trackWPNames[iWP].c_str (), iMult));
          BinomialDivide (h2_efficiency[iSys][iPID][iWP][iMult], h2_truth_matched_primary_tracks[iSys][iPID][iWP][iMult], h2_truth_tracks[iSys][iPID][iMult], h2_truth_tracks_wgt2[iSys][iPID][iMult]);

          for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {

            h_efficiency[iSys][iPID][iWP][iMult][iEta] = (TH1D*) h_truth_matched_primary_tracks[iSys][iPID][iWP][iMult][iEta]->Clone (Form ("h_efficiency_%s_PID%i_%s_iMult%i_iEta%i", sys.Data (), PIDs[iPID], trackWPNames[iWP].c_str (), iMult, iEta));
            BinomialDivide (h_efficiency[iSys][iPID][iWP][iMult][iEta], h_truth_matched_primary_tracks[iSys][iPID][iWP][iMult][iEta], h_truth_tracks[iSys][iPID][iMult][iEta], h_truth_tracks_wgt2[iSys][iPID][iMult][iEta]);

          } // end loop over iEta

        } // end loop over iMult

      } // end loop over iWP

    } // end loop over iPID


    for (int iWP = 0; iWP < trackWPs.size (); iWP++) {

      //for (int iMult = 0; iMult < nMultBins; iMult++) {
      for (int iDR = 0; iDR < nDRBins; iDR++) {

        h2_fake_rate[iSys][iWP][iDR] = (TH2D*) h2_fake_tracks[iSys][iWP][iDR]->Clone (Form ("h2_fake_rate_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR));
        BinomialDivide (h2_fake_rate[iSys][iWP][iDR], h2_fake_tracks[iSys][iWP][iDR], h2_reco_tracks[iSys][iWP][iDR], h2_reco_tracks_wgt2[iSys][iWP][iDR]);
        InvertRate (h2_fake_rate[iSys][iWP][iDR]);
        h2_secondary_rate[iSys][iWP][iDR] = (TH2D*) h2_secondary_tracks[iSys][iWP][iDR]->Clone (Form ("h2_secondary_rate_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR));
        BinomialDivide (h2_secondary_rate[iSys][iWP][iDR], h2_secondary_tracks[iSys][iWP][iDR], h2_reco_tracks[iSys][iWP][iDR], h2_reco_tracks_wgt2[iSys][iWP][iDR]);
        InvertRate (h2_secondary_rate[iSys][iWP][iDR]);
        h2_strange_rate[iSys][iWP][iDR] = (TH2D*) h2_strange_tracks[iSys][iWP][iDR]->Clone (Form ("h2_strange_rate_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR));
        BinomialDivide (h2_strange_rate[iSys][iWP][iDR], h2_strange_tracks[iSys][iWP][iDR], h2_reco_tracks[iSys][iWP][iDR], h2_reco_tracks_wgt2[iSys][iWP][iDR]);
        InvertRate (h2_strange_rate[iSys][iWP][iDR]);
        h2_primary_rate[iSys][iWP][iDR] = (TH2D*) h2_primary_tracks[iSys][iWP][iDR]->Clone (Form ("h2_primary_rate_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR));
        BinomialDivide (h2_primary_rate[iSys][iWP][iDR], h2_primary_tracks[iSys][iWP][iDR], h2_reco_tracks[iSys][iWP][iDR], h2_reco_tracks_wgt2[iSys][iWP][iDR]);
        h2_primary_rate_fakes_p100[iSys][iWP][iDR] = (TH2D*) h2_primary_tracks_fakes_p100[iSys][iWP][iDR]->Clone (Form ("h2_primary_rate_fakes_p100_%s_%s_iDR%i", sys.Data (), trackWPNames[iWP].c_str (), iDR));
        BinomialDivide (h2_primary_rate_fakes_p100[iSys][iWP][iDR], h2_primary_tracks_fakes_p100[iSys][iWP][iDR], h2_reco_tracks[iSys][iWP][iDR], h2_reco_tracks_wgt2[iSys][iWP][iDR]);

        for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {

          h_fake_rate[iSys][iWP][iDR][iEta] = (TH1D*) h_fake_tracks[iSys][iWP][iDR][iEta]->Clone (Form ("h_fake_rate_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta));
          BinomialDivide (h_fake_rate[iSys][iWP][iDR][iEta], h_fake_tracks[iSys][iWP][iDR][iEta], h_reco_tracks[iSys][iWP][iDR][iEta], h_reco_tracks_wgt2[iSys][iWP][iDR][iEta]);
          InvertRate (h_fake_rate[iSys][iWP][iDR][iEta]);
          h_secondary_rate[iSys][iWP][iDR][iEta] = (TH1D*) h_secondary_tracks[iSys][iWP][iDR][iEta]->Clone (Form ("h_secondary_rate_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta));
          BinomialDivide (h_secondary_rate[iSys][iWP][iDR][iEta], h_secondary_tracks[iSys][iWP][iDR][iEta], h_reco_tracks[iSys][iWP][iDR][iEta], h_reco_tracks_wgt2[iSys][iWP][iDR][iEta]);
          InvertRate (h_secondary_rate[iSys][iWP][iDR][iEta]);
          h_strange_rate[iSys][iWP][iDR][iEta] = (TH1D*) h_strange_tracks[iSys][iWP][iDR][iEta]->Clone (Form ("h_strange_rate_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta));
          BinomialDivide (h_strange_rate[iSys][iWP][iDR][iEta], h_strange_tracks[iSys][iWP][iDR][iEta], h_reco_tracks[iSys][iWP][iDR][iEta], h_reco_tracks_wgt2[iSys][iWP][iDR][iEta]);
          InvertRate (h_strange_rate[iSys][iWP][iDR][iEta]);
          h_primary_rate[iSys][iWP][iDR][iEta] = (TH1D*) h_primary_tracks[iSys][iWP][iDR][iEta]->Clone (Form ("h_primary_rate_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta));
          BinomialDivide (h_primary_rate[iSys][iWP][iDR][iEta], h_primary_tracks[iSys][iWP][iDR][iEta], h_reco_tracks[iSys][iWP][iDR][iEta], h_reco_tracks_wgt2[iSys][iWP][iDR][iEta]);
          h_primary_rate_fakes_p100[iSys][iWP][iDR][iEta] = (TH1D*) h_primary_tracks_fakes_p100[iSys][iWP][iDR][iEta]->Clone (Form ("h_primary_rate_fakes_p100_%s_%s_iDR%i_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iDR, iEta));
          BinomialDivide (h_primary_rate_fakes_p100[iSys][iWP][iDR][iEta], h_primary_tracks_fakes_p100[iSys][iWP][iDR][iEta], h_reco_tracks[iSys][iWP][iDR][iEta], h_reco_tracks_wgt2[iSys][iWP][iDR][iEta]);

        } // end loop over iEta

      } // end loop over iDR

    } // end loop over iWP

  } // end loop over iSys



  for (int iSys : systems) {

    const TString sys = (iSys == 0 ? "pp" : "pPb");

    for (int iWP = 0; iWP < trackWPs.size (); iWP++) {

      //for (int iMult = 0; iMult < nMultBins; iMult++) {
      //for (int iMult = nMultBins-1; iMult < nMultBins; iMult++) {
      for (int iDR = 0; iDR < nDRBins; iDR++) {

        for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {

          f_primary_rate[iSys][iWP][iDR][iEta] = DoPurityFit (h_primary_rate[iSys][iWP][iDR][iEta], 8, 2);
          f_primary_rate_fakes_p100[iSys][iWP][iDR][iEta] = DoPurityFit (h_primary_rate_fakes_p100[iSys][iWP][iDR][iEta], 8, 2);

        } // end loop over iEta

      } // end loop over iDR

    } // end loop over iWP

  } // end loop over iSys




  TLatex* tl = new TLatex ();
  TLine* l = new TLine ();

  const short iWP = 0;



  for (int iSys : systems) {

    const int iMult = nMultBins-1;

    const char* cname = Form ("c_truth_matching_prob_%s", iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (cname, "", 800, 800);

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;
    c->SetLeftMargin (lMargin);
    c->SetRightMargin (rMargin);
    c->SetBottomMargin (bMargin);
    c->SetTopMargin (tMargin);

    c->SetLogy ();

    TH1D* htemp = new TH1D ("htemp", "", 1, 0, 1);
    htemp->SetBinContent (1, 0);

    TAxis* xax = htemp->GetXaxis ();
    TAxis* yax = htemp->GetYaxis ();

    xax->SetTitle ("Truth matching probability");
    xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
    xax->SetLabelFont (43);
    xax->SetLabelSize (32);

    yax->SetTitle ("Counts of hadrons");
    yax->SetLabelFont (43);
    yax->SetLabelSize (32);
    yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
    const double ymin = 0.5;
    const double ymax = 1e10;
    yax->SetRangeUser (ymin, ymax);

    htemp->SetLineWidth (1);
    htemp->SetLineStyle (2);
    
    htemp->DrawCopy ("hist");
    SaferDelete (&htemp);

    TH1D* h_pp = (TH1D*) h_truth_matching_prob[0][iWP]->Clone ("h_pp");
    h_pp->Rebin (10);
    myDraw (h_pp, kAzure+2, kOpenSquare, 1.4);

    TH1D* h_pPb = (TH1D*) h_truth_matching_prob[1][iWP]->Clone ("h_pPb");
    h_pPb->Rebin (10);
    myDraw (h_pPb, kRed+1, kFullCircle, 1.4);

    tl->SetTextFont (43);
    tl->SetTextSize (32);
    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.88, "#bf{#it{ATLAS}} Simulation Internal");
    tl->SetTextSize (24);
    myLineText2 (0.26, 0.84, kAzure+2, kOpenSquare, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV", 1.4, 0.032, true);
    myLineText2 (0.26, 0.80, kRed+1, kFullCircle, "Pythia8 + #it{p}+Pb overlay, #sqrt{s_{NN}} = 5.02 TeV", 1.4, 0.032, true);
    tl->DrawLatexNDC (0.22, 0.75, Form ("%s tracks", trackWPStrs[iWP].c_str ()));

    TLine* l0p3 = new TLine (0.3, 0.5, 0.3, 1.4e5);                                
    TLine* l0p5 = new TLine (0.5, 0.5, 0.5, 1.4e5);                                
    TLine* l0p6 = new TLine (0.6, 0.5, 0.6, 1.4e5);                                
                                                                                   
    l0p3->SetLineWidth (1);                                                        
    l0p3->SetLineStyle (2);                                                        
    l0p5->SetLineWidth (1);                                                        
    l0p5->SetLineStyle (2);                                                        
    l0p6->SetLineWidth (1);                                                        
    l0p6->SetLineStyle (2);                                                        
                                                                                   
    l0p3->Draw ("same");                                                           
    l0p5->Draw ("same");                                                           
    l0p6->Draw ("same");    

    const float i_pp_0p3 = 100 * h_pp->Integral (h_pp->FindFixBin (0.3), h_pp->FindFixBin (1)) / h_pp->Integral (h_pp->FindFixBin (0), h_pp->FindFixBin (1));
    const float i_pp_0p5 = 100 * h_pp->Integral (h_pp->FindFixBin (0.5), h_pp->FindFixBin (1)) / h_pp->Integral (h_pp->FindFixBin (0), h_pp->FindFixBin (1));
    const float i_pp_0p6 = 100 * h_pp->Integral (h_pp->FindFixBin (0.6), h_pp->FindFixBin (1)) / h_pp->Integral (h_pp->FindFixBin (0), h_pp->FindFixBin (1));
    const float i_pPb_0p3 = 100 * h_pPb->Integral (h_pPb->FindFixBin (0.3), h_pPb->FindFixBin (1)) / h_pPb->Integral (h_pPb->FindFixBin (0), h_pPb->FindFixBin (1));
    const float i_pPb_0p5 = 100 * h_pPb->Integral (h_pPb->FindFixBin (0.5), h_pPb->FindFixBin (1)) / h_pPb->Integral (h_pPb->FindFixBin (0), h_pPb->FindFixBin (1));
    const float i_pPb_0p6 = 100 * h_pPb->Integral (h_pPb->FindFixBin (0.6), h_pPb->FindFixBin (1)) / h_pPb->Integral (h_pPb->FindFixBin (0), h_pPb->FindFixBin (1));

    myText (0.35, 0.680, kBlack, "f_{ch}(> 0.3)", 0.026);
    myText (0.50, 0.680, kBlack, "f_{ch}(> 0.5)", 0.026);
    myText (0.62, 0.680, kBlack, "f_{ch}(> 0.6)", 0.026);

    myText (0.35, 0.645, kRed+1, Form ("%.3f%%", i_pPb_0p3), 0.026);
    myText (0.50, 0.645, kRed+1, Form ("%.3f%%", i_pPb_0p5), 0.026);
    myText (0.62, 0.645, kRed+1, Form ("%.3f%%", i_pPb_0p6), 0.026);
    myText (0.35, 0.610, kAzure+2, Form ("%.3f%%", i_pp_0p3), 0.026);
    myText (0.50, 0.610, kAzure+2, Form ("%.3f%%", i_pp_0p5), 0.026);
    myText (0.62, 0.610, kAzure+2, Form ("%.3f%%", i_pp_0p6), 0.026);

    c->SaveAs (Form ("%s/Plots/TrackingPerformance/TruthMatchingProb.pdf", workPath.Data ()));
  }



  //for (int iSys : systems) {
  {

    //const char* cname = Form ("c_eff_sum_vs_mult_%s", iSys == 0 ? "pp" : "pPb");
    const char* cname = Form ("c_eff_sum_vs_mult");

    TCanvas* c = new TCanvas (cname, "", 800*(nMultBins-1), 800);
    c->Divide (3, 1);

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    for (int iMult = 0; iMult < (nMultBins-1); iMult++) {

      c->cd (iMult+1);

      gPad->SetLeftMargin (lMargin);
      gPad->SetRightMargin (rMargin);
      gPad->SetBottomMargin (bMargin);
      gPad->SetTopMargin (tMargin);

      gPad->SetLogx ();

      TH1D* htemp = new TH1D ("htemp", "", 1, pTchBins[0], pTchBins[nPtchBins]);
      htemp->SetBinContent (1, 1);
  
      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetLabelSize (0);

      yax->SetTitle ("Track Reco. Efficiency");
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      const double ymin = 0.5;
      const double ymax = 1.06;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (1);
      htemp->SetLineStyle (2);
  
      htemp->DrawCopy ("hist");
      SaferDelete (&htemp);

      tl->SetTextFont (43);
      tl->SetTextSize (32);
      tl->SetTextAlign (21);
      
      const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
      //tl->DrawLatex (0.5,  yoff, "0.5");
      tl->DrawLatex (0.7,  yoff, "0.7");
      tl->DrawLatex (1,  yoff, "1");
      tl->DrawLatex (2,  yoff, "2");
      tl->DrawLatex (3,  yoff, "3");
      tl->DrawLatex (4,  yoff, "4");
      tl->DrawLatex (5,  yoff, "5");
      tl->DrawLatex (6,  yoff, "6");
      tl->DrawLatex (7,  yoff, "7");
      tl->DrawLatex (10, yoff, "10");
      tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");
      //tl->DrawLatex (80, yoff, "80");
      tl->DrawLatex (100, yoff, "100");

      for (int iSys : systems)
        for (int iEta = 0; iEta < nEtaTrkBins; iEta++)
          myDraw (h_efficiency[iSys][nPIDs-1][iWP][iMult][iEta], colors[iEta], iSys == 0 ? kOpenCircle : kOpenSquare, 1.0);

      tl->SetTextColor (kBlack);
      tl->SetTextAlign (11);
      tl->SetTextSize (28);
      tl->DrawLatexNDC (0.22, 0.22, Form ("N_{ch}^{rec} = %i-%i", (int) std::ceil (multBins[iMult]), (int) std::floor (multBins[iMult+1])));

    } // end loop over iMult

    c->cd (1);
    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.88, "#bf{#it{ATLAS}} Simulation Internal");
    tl->SetTextSize (24);
    //tl->DrawLatexNDC (0.22, 0.84, iSys == 0 ? "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV" : "Pythia8 + #it{p}+Pb overlay, #sqrt{s_{NN}} = 5.02 TeV");
    //tl->DrawLatexNDC (0.22, 0.80, Form ("%s tracks", trackWPStrs[iWP].c_str ()));
    tl->DrawLatexNDC (0.22, 0.84, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");
    tl->DrawLatexNDC (0.22, 0.80, "Pythia8 + #it{p}+Pb overlay, #sqrt{s_{NN}} = 5.02 TeV");
    tl->DrawLatexNDC (0.22, 0.76, Form ("%s tracks", trackWPStrs[iWP].c_str ()));

    c->cd (3);
    tl->SetTextSize (20);
    tl->DrawLatexNDC (0.59, 0.376, "#it{pp}");
    tl->DrawLatexNDC (0.675, 0.376, "#it{p}+Pb");
    for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {
      myLineText2 (0.74, 0.34-iEta*0.036, colors[iEta], kOpenSquare, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.026, true);
      myLineText2 (0.64, 0.34-iEta*0.036, colors[iEta], kOpenCircle, "", 1.2, 0.026, true);
    }

    //c->SaveAs (Form ("%s/Plots/TrackingPerformance/EfficiencySummary_vsMultiplicity_%s.pdf", workPath.Data (), iSys == 0 ? "pp" : "pPb"));
    c->SaveAs (Form ("%s/Plots/TrackingPerformance/EfficiencySummary_vsMultiplicity.pdf", workPath.Data ()));
  } // end loop over iSys



  for (int iSys : systems) {

    const int iMult = nMultBins-1;

    const char* cname = Form ("c_eff_sum_%s", iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (cname, "", 800, 800);

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c->SetLeftMargin (lMargin);
    c->SetRightMargin (rMargin);
    c->SetBottomMargin (bMargin);
    c->SetTopMargin (tMargin);

    c->SetLogx ();

    TH1D* htemp = new TH1D ("htemp", "", 1, pTchBins[0], pTchBins[nPtchBins]);
    htemp->SetBinContent (1, 1);
  
    TAxis* xax = htemp->GetXaxis ();
    TAxis* yax = htemp->GetYaxis ();

    xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
    xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
    xax->SetLabelSize (0);

    yax->SetTitle ("Track Reco. Efficiency");
    yax->SetLabelFont (43);
    yax->SetLabelSize (32);
    yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
    const double ymin = 0.5;
    const double ymax = 1.06;
    yax->SetRangeUser (ymin, ymax);

    htemp->SetLineWidth (1);
    htemp->SetLineStyle (2);
  
    htemp->DrawCopy ("hist");
    SaferDelete (&htemp);

    tl->SetTextFont (43);
    tl->SetTextSize (32);
    tl->SetTextAlign (21);
    
    const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
    //tl->DrawLatex (0.5,  yoff, "0.5");
    tl->DrawLatex (0.7,  yoff, "0.7");
    tl->DrawLatex (1,  yoff, "1");
    tl->DrawLatex (2,  yoff, "2");
    tl->DrawLatex (3,  yoff, "3");
    tl->DrawLatex (4,  yoff, "4");
    tl->DrawLatex (5,  yoff, "5");
    tl->DrawLatex (6,  yoff, "6");
    tl->DrawLatex (7,  yoff, "7");
    tl->DrawLatex (10, yoff, "10");
    tl->DrawLatex (20, yoff, "20");
    tl->DrawLatex (30, yoff, "30");
    tl->DrawLatex (40, yoff, "40");
    tl->DrawLatex (60, yoff, "60");
    //tl->DrawLatex (80, yoff, "80");
    tl->DrawLatex (100, yoff, "100");

    for (int iEta = 0; iEta < nEtaTrkBins; iEta++)
      myDraw (h_efficiency[iSys][nPIDs-1][iWP][iMult][iEta], colors[iEta], kOpenCircle, 0.8);

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.88, "#bf{#it{ATLAS}} Simulation Internal");
    tl->SetTextSize (24);
    tl->DrawLatexNDC (0.22, 0.84, iSys == 0 ? "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV" : "Pythia8 + #it{p}+Pb overlay, #sqrt{s_{NN}} = 5.02 TeV");
    tl->DrawLatexNDC (0.22, 0.80, Form ("%s tracks", trackWPStrs[iWP].c_str ()));
 
    for (int iEta = 0; iEta < nEtaTrkBins; iEta++)
      myLineText2 (0.74, 0.34-iEta*0.036, colors[iEta], kOpenCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.0, 0.026, true);

    c->SaveAs (Form ("%s/Plots/TrackingPerformance/EfficiencySummary_%s.pdf", workPath.Data (), iSys == 0 ? "pp" : "pPb"));
  } // end loop over iSys



  for (int iSys : systems) {

    const char* cname = Form ("c_pur_sum_vs_dr_%s", iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (cname, "", 800*(nDRBins-1), 800);
    c->Divide (3, 1);

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    for (int iDR = 0; iDR < (nDRBins-1); iDR++) {

      c->cd (iDR+1);

      gPad->SetLeftMargin (lMargin);
      gPad->SetRightMargin (rMargin);
      gPad->SetBottomMargin (bMargin);
      gPad->SetTopMargin (tMargin);

      gPad->SetLogx ();

      {
        TH1D* htemp = new TH1D ("htemp", "", 1, pTchBins[0], pTchBins[nPtchBins]);
        htemp->SetBinContent (1, 1);
  
        TAxis* xax = htemp->GetXaxis ();
        TAxis* yax = htemp->GetYaxis ();

        xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
        xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
        xax->SetLabelSize (0);

        yax->SetTitle ("Primary Fraction");
        yax->SetLabelFont (43);
        yax->SetLabelSize (32);
        yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
        const double ymin = 0.64;
        const double ymax = 1.12;
        yax->SetRangeUser (ymin, ymax);

        htemp->SetLineWidth (1);
        htemp->SetLineStyle (2);

        htemp->DrawCopy ("hist");
        SaferDelete (&htemp);

        tl->SetTextFont (43);
        tl->SetTextSize (32);
        tl->SetTextAlign (21);
        
        const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
        //tl->DrawLatex (0.5,  yoff, "0.5");
        tl->DrawLatex (0.7,  yoff, "0.7");
        tl->DrawLatex (1,  yoff, "1");
        tl->DrawLatex (2,  yoff, "2");
        tl->DrawLatex (3,  yoff, "3");
        tl->DrawLatex (4,  yoff, "4");
        tl->DrawLatex (5,  yoff, "5");
        tl->DrawLatex (6,  yoff, "6");
        tl->DrawLatex (7,  yoff, "7");
        tl->DrawLatex (10, yoff, "10");
        tl->DrawLatex (20, yoff, "20");
        tl->DrawLatex (30, yoff, "30");
        tl->DrawLatex (40, yoff, "40");
        tl->DrawLatex (60, yoff, "60");
        //tl->DrawLatex (80, yoff, "80");
        tl->DrawLatex (100, yoff, "100");
      }

      for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {
        myDraw (h_primary_rate[iSys][iWP][iDR][iEta], colors[iEta], kOpenCircle, 0.6);
        TF1* f = f_primary_rate[iSys][iWP][iDR][iEta];
        if (f) {
          f->SetLineColor (colors[iEta]);
          f->SetLineWidth (1);
          f->Draw ("same");
        }
      }

      tl->SetTextColor (kBlack);
      tl->SetTextAlign (11);

      tl->SetTextSize (28);
      tl->DrawLatexNDC (0.22, 0.88, "#bf{#it{ATLAS}} Simulation Internal");
      tl->SetTextSize (24);
      tl->DrawLatexNDC (0.22, 0.84, iSys == 0 ? "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV" : "Hijing #it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV");
      tl->DrawLatexNDC (0.22, 0.80, Form ("%s tracks", trackWPStrs[iWP].c_str ()));

      tl->SetTextSize (28);
      tl->DrawLatexNDC (0.22, 0.22, Form ("%g < #Delta R_{ch, jet} < %g", drBins[nDRBins-iDR-2], drBins[nDRBins-iDR-1]));

      tl->SetTextSize (18);
      tl->SetTextAlign (21);
      tl->DrawLatexNDC (0.705, 0.37, "Primaries");
      for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {
        myLineText2 (0.74, 0.34-iEta*0.036, colors[iEta], kOpenCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.026, true);
      }

    } // end loop over iDR

    c->SaveAs (Form ("%s/Plots/TrackingPerformance/PuritySummary_vs_JetDR_%s.pdf", workPath.Data (), iSys == 0 ? "pp" : "pPb"));
  } // end loop over iSys



  for (int iSys : systems) {

    const int iDR = nDRBins-1;
    const char* cname = Form ("c_pur_sum_%s", iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (cname, "", 800, 800);

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c->SetLeftMargin (lMargin);
    c->SetRightMargin (rMargin);
    c->SetBottomMargin (bMargin);
    c->SetTopMargin (tMargin);

    c->SetLogx ();

    TH1D* htemp = new TH1D ("htemp", "", 1, pTchBins[0], pTchBins[nPtchBins]);
    htemp->SetBinContent (1, 1);
  
    TAxis* xax = htemp->GetXaxis ();
    TAxis* yax = htemp->GetYaxis ();

    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
    xax->SetLabelSize (0);

    yax->SetTitle ("Primary Fraction");
    yax->SetLabelFont (43);
    yax->SetLabelSize (32);
    yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
    const double ymin = 0.82;
    const double ymax = 1.06;
    yax->SetRangeUser (ymin, ymax);

    htemp->SetLineWidth (1);
    htemp->SetLineStyle (2);

    htemp->DrawCopy ("hist");
    SaferDelete (&htemp);

    tl->SetTextFont (43);
    tl->SetTextSize (32);
    tl->SetTextAlign (21);
    
    const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
    //tl->DrawLatex (0.5,  yoff, "0.5");
    tl->DrawLatex (0.7,  yoff, "0.7");
    tl->DrawLatex (1,  yoff, "1");
    tl->DrawLatex (2,  yoff, "2");
    tl->DrawLatex (3,  yoff, "3");
    tl->DrawLatex (4,  yoff, "4");
    tl->DrawLatex (5,  yoff, "5");
    tl->DrawLatex (6,  yoff, "6");
    tl->DrawLatex (7,  yoff, "7");
    tl->DrawLatex (10, yoff, "10");
    tl->DrawLatex (20, yoff, "20");
    tl->DrawLatex (30, yoff, "30");
    tl->DrawLatex (40, yoff, "40");
    tl->DrawLatex (60, yoff, "60");
    //tl->DrawLatex (80, yoff, "80");
    tl->DrawLatex (100, yoff, "100");

    for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {
      myDraw (h_primary_rate[iSys][iWP][iDR][iEta], colors[iEta], kOpenCircle, 0.8);
      TF1* f = f_primary_rate[iSys][iWP][iDR][iEta];
      if (f != nullptr) {
        f->SetLineColor (colors[iEta]);
        f->SetLineWidth (1);
        f->Draw ("same");
      }
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.89, "#bf{#it{ATLAS}} Simulation Internal");
    tl->SetTextSize (24);
    tl->DrawLatexNDC (0.22, 0.85, iSys == 0 ? "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV" : "Hijing #it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV");
    tl->DrawLatexNDC (0.22, 0.81, Form ("%s tracks", trackWPStrs[iWP].c_str ()));

    tl->SetTextSize (18);
    tl->SetTextAlign (21);
    tl->DrawLatexNDC (0.705, 0.37, "Primaries");
    for (int iEta = 0; iEta < nEtaTrkBins; iEta++)
      myLineText2 (0.74, 0.34-iEta*0.036, colors[iEta], kOpenCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.026, true);

    c->SaveAs (Form ("%s/Plots/TrackingPerformance/PuritySummary_%s.pdf", workPath.Data (), iSys == 0 ? "pp" : "pPb"));
  } // end loop over iSys



  for (int iSys : systems) {

    const int iEta = 0;
    const int iMult = nMultBins-1;

    const char* cname = Form ("c_eff_com_%s", iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (cname, "", 800, 800);

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c->SetLeftMargin (lMargin);
    c->SetRightMargin (rMargin);
    c->SetBottomMargin (bMargin);
    c->SetTopMargin (tMargin);

    c->SetLogx ();

    TH1D* htemp = new TH1D ("htemp", "", 1, pTchBins[0], pTchBins[nPtchBins]);
    htemp->SetBinContent (1, 1);
    
    TAxis* xax = htemp->GetXaxis ();
    TAxis* yax = htemp->GetYaxis ();
  
    xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
    xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
    xax->SetLabelSize (0);
  
    yax->SetTitle ("Track Reco. Efficiency");
    yax->SetLabelFont (43);
    yax->SetLabelSize (32);
    yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
    const double ymin = 0.0;
    const double ymax = 1.25;
    yax->SetRangeUser (ymin, ymax);
  
    htemp->SetLineWidth (1);
    htemp->SetLineStyle (2);
  
    htemp->DrawCopy ("hist");
    SaferDelete (&htemp);
  
    tl->SetTextFont (43);
    tl->SetTextSize (32);
    tl->SetTextAlign (21);
    
    const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
    //tl->DrawLatex (0.5,  yoff, "0.5");
    tl->DrawLatex (0.7,  yoff, "0.7");
    tl->DrawLatex (1,  yoff, "1");
    tl->DrawLatex (2,  yoff, "2");
    tl->DrawLatex (3,  yoff, "3");
    tl->DrawLatex (4,  yoff, "4");
    tl->DrawLatex (5,  yoff, "5");
    tl->DrawLatex (6,  yoff, "6");
    tl->DrawLatex (7,  yoff, "7");
    tl->DrawLatex (10, yoff, "10");
    tl->DrawLatex (20, yoff, "20");
    tl->DrawLatex (30, yoff, "30");
    tl->DrawLatex (40, yoff, "40");
    tl->DrawLatex (60, yoff, "60");
    //tl->DrawLatex (80, yoff, "80");
    tl->DrawLatex (100, yoff, "100");

    for (int iPID = 0; iPID < nPIDs; iPID++)
      myDraw (h_efficiency[iSys][iPID][iWP][iMult][iEta], systColors[iPID], kOpenCircle, 0.8);
    myDraw (h_efficiency[iSys][nPIDs-1][iWP][iMult][iEta], kBlack, kFullCircle, 0.8);

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.89, "#bf{#it{ATLAS}} Simulation Internal");
    tl->SetTextSize (24);
    tl->DrawLatexNDC (0.22, 0.85, iSys == 0 ? "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV" : "Pythia8 + #it{p}+Pb overlay, #sqrt{s_{NN}} = 5.02 TeV");
    tl->DrawLatexNDC (0.22, 0.81, Form ("%s tracks", trackWPStrs[iWP].c_str ()));

    for (int iPID = 0; iPID < nPIDs-1; iPID++)
      myLineText2 (0.25+(iPID>=nPIDs/2 ? 0.13 : 0), 0.42-(iPID%(nPIDs/2))*0.036, systColors[iPID], kOpenCircle, partNames[iPID].c_str (), 1.2, 0.026, true);
    myLineText2 (0.38, 0.42-((nPIDs-1)%(nPIDs/2))*0.036, kBlack, kFullCircle, partNames[nPIDs-1].c_str (), 1.2, 0.026, true);

    c->SaveAs (Form ("%s/Plots/TrackingPerformance/EfficiencyComposition_%s.pdf", workPath.Data (), iSys == 0 ? "pp" : "pPb"));
  } // end loop over iSys


  for (int iSys : systems) {

    const int iDR = nDRBins-1;

    const char* cname = Form ("c_pur_comp_%s", iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (cname, "", 800, 800);

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c->SetLeftMargin (lMargin);
    c->SetRightMargin (rMargin);
    c->SetBottomMargin (bMargin);
    c->SetTopMargin (tMargin);

    c->SetLogx ();

    TH1D* htemp = new TH1D ("htemp", "", 1, pTchBins[0], pTchBins[nPtchBins]);
    htemp->SetBinContent (1, 1);

    TAxis* xax = htemp->GetXaxis ();
    TAxis* yax = htemp->GetYaxis ();

    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
    xax->SetLabelSize (0);

    yax->SetTitle ("Fraction of tracks");
    yax->SetLabelFont (43);
    yax->SetLabelSize (32);
    yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
    const double ymin = 0.82;
    const double ymax = 1.06;
    yax->SetRangeUser (ymin, ymax);

    htemp->SetLineWidth (1);
    htemp->SetLineStyle (2);

    htemp->DrawCopy ("hist");
    SaferDelete (&htemp);

    tl->SetTextFont (43);
    tl->SetTextSize (32);
    tl->SetTextAlign (21);
    
    const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
    //tl->DrawLatex (0.5,  yoff, "0.5");
    tl->DrawLatex (0.7,  yoff, "0.7");
    tl->DrawLatex (1,  yoff, "1");
    tl->DrawLatex (2,  yoff, "2");
    tl->DrawLatex (3,  yoff, "3");
    tl->DrawLatex (4,  yoff, "4");
    tl->DrawLatex (5,  yoff, "5");
    tl->DrawLatex (6,  yoff, "6");
    tl->DrawLatex (7,  yoff, "7");
    tl->DrawLatex (10, yoff, "10");
    tl->DrawLatex (20, yoff, "20");
    tl->DrawLatex (30, yoff, "30");
    tl->DrawLatex (40, yoff, "40");
    tl->DrawLatex (60, yoff, "60");
    //tl->DrawLatex (80, yoff, "80");
    tl->DrawLatex (100, yoff, "100");

    for (int iEta : {0, 3})
      myDraw (h_fake_rate[iSys][iWP][iDR][iEta], colors[iEta], kOpenCrossX, 1.0);

    for (int iEta : {0, 3})
      myDraw (h_secondary_rate[iSys][iWP][iDR][iEta], colors[iEta], kOpenCross, 1.0);

    for (int iEta : {0, 3})
      myDraw (h_strange_rate[iSys][iWP][iDR][iEta], colors[iEta], kOpenSquare, 1.0);

    for (int iEta : {0, 3})
      myDraw (h_primary_rate[iSys][iWP][iDR][iEta], colors[iEta], kOpenCircle, 1.0);

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.89, "#bf{#it{ATLAS}} Simulation Internal");
    tl->SetTextSize (24);
    tl->DrawLatexNDC (0.22, 0.85, iSys == 0 ? "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV" : "Hijing #it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV");
    tl->DrawLatexNDC (0.22, 0.81, Form ("%s tracks", trackWPStrs[iWP].c_str ()));

    tl->SetTextSize (18);
    tl->SetTextAlign (21);
    tl->DrawLatexNDC (0.345, 0.262, "Fakes");
    tl->DrawLatexNDC (0.465, 0.262, "Secondaries");
    tl->DrawLatexNDC (0.585, 0.262, "#Sigma^{#pm}, #Xi^{-}, #Omega^{-}");
    tl->DrawLatexNDC (0.705, 0.262, "Primaries");

    for (int iEta : {0, 3}) {
      myLineText2 (0.38, 0.232-(iEta/3)*0.036, colors[iEta], kOpenCrossX, "", 0.8, 0.026, true);
      myLineText2 (0.50, 0.232-(iEta/3)*0.036, colors[iEta], kOpenCross, "", 0.8, 0.026, true);
      myLineText2 (0.62, 0.232-(iEta/3)*0.036, colors[iEta], kOpenSquare, "", 0.8, 0.026, true);
      myLineText2 (0.74, 0.232-(iEta/3)*0.036, colors[iEta], kOpenCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.026, true);
    }

    c->SaveAs (Form ("%s/Plots/TrackingPerformance/PurityComponents_%s.pdf", workPath.Data (), iSys == 0 ? "pp" : "pPb"));
  } // end loop over iSys



  for (int iSys : systems) {

    const int iMult = nMultBins-1;
    const int iPID = nPIDs-1;

    const char* cname = Form ("c_eff_all_%s", iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (cname, "", 880, 800);

    const double lMargin = 0.15;
    const double rMargin = 0.15;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c->SetLeftMargin (lMargin);
    c->SetRightMargin (rMargin);
    c->SetBottomMargin (bMargin);
    c->SetTopMargin (tMargin);

    c->SetLogy ();

    TH2D* h2 = (TH2D*) h2_efficiency[iSys][iPID][iWP][iMult]->Clone ("temp");

    TAxis* xax = h2->GetXaxis ();
    TAxis* yax = h2->GetYaxis ();
    TAxis* zax = h2->GetZaxis ();

    xax->SetTitle ("#eta");
    yax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    zax->SetTitle ("Track Reco. Efficiency");

    const double zmin = 0.5;
    const double zmax = 1;

    h2->SetLineWidth (0);

    yax->SetMoreLogLabels ();
    zax->SetRangeUser (zmin, zmax);

    xax->SetTitleFont (43);
    xax->SetTitleSize (32);
    yax->SetTitleFont (43);
    yax->SetTitleSize (32);
    zax->SetTitleFont (43);
    zax->SetTitleSize (32);
    zax->SetTitleOffset (1.1*zax->GetTitleOffset ());
    xax->SetLabelFont (43);
    xax->SetLabelSize (32);
    yax->SetLabelFont (43);
    yax->SetLabelSize (32);
    zax->SetLabelFont (43);
    zax->SetLabelSize (24);

    h2->DrawCopy ("colz");
    SaferDelete (&h2);

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (12);

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    tl->SetTextSize (24);
    tl->DrawLatexNDC (0.22, 0.850, iSys == 0 ? "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV" : "Pythia8 + #it{p}+Pb Overlay, #sqrt{s_{NN}} = 5.02 TeV");
    tl->DrawLatexNDC (0.22, 0.810, Form ("%s tracks", trackWPStrs[iWP].c_str ()));

    c->SaveAs (Form ("%s/Plots/TrackingPerformance/EfficiencyMap_%s.pdf", workPath.Data (), iSys == 0 ? "pp" : "pPb"));
  } // end loop over iSys



  for (int iSys : systems) {

    const int iDR = nDRBins-1;

    const char* cname = Form ("c_pur_all_%s", iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (cname, "", 880, 800);

    const double lMargin = 0.15;
    const double rMargin = 0.15;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c->SetLeftMargin (lMargin);
    c->SetRightMargin (rMargin);
    c->SetBottomMargin (bMargin);
    c->SetTopMargin (tMargin);

    c->SetLogy ();

    TH2D* h2 = (TH2D*) h2_primary_rate[iSys][iWP][iDR]->Clone ("temp");

    TAxis* xax = h2->GetXaxis ();
    TAxis* yax = h2->GetYaxis ();
    TAxis* zax = h2->GetZaxis ();

    xax->SetTitle ("#eta");
    yax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    zax->SetTitle ("Primary Fraction");

    const double zmin = 0.5;
    const double zmax = 1;

    h2->SetLineWidth (0);

    yax->SetMoreLogLabels ();
    zax->SetRangeUser (zmin, zmax);

    xax->SetTitleFont (43);
    xax->SetTitleSize (32);
    yax->SetTitleFont (43);
    yax->SetTitleSize (32);
    zax->SetTitleFont (43);
    zax->SetTitleSize (32);
    zax->SetTitleOffset (1.1*zax->GetTitleOffset ());
    xax->SetLabelFont (43);
    xax->SetLabelSize (32);
    yax->SetLabelFont (43);
    yax->SetLabelSize (32);
    zax->SetLabelFont (43);
    zax->SetLabelSize (24);

    h2->DrawCopy ("colz");
    SaferDelete (&h2);

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (12);

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    tl->SetTextSize (24);
    tl->DrawLatexNDC (0.22, 0.850, iSys == 0 ? "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV" : "Hijing #it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV");
    tl->DrawLatexNDC (0.22, 0.810, Form ("%s tracks", trackWPStrs[iWP].c_str ()));

    c->SaveAs (Form ("%s/Plots/TrackingPerformance/PurityMap_%s.pdf", workPath.Data (), iSys == 0 ? "pp" : "pPb"));
  } // end loop over iSys




  TFile* outFile = new TFile (Form ("%s/aux/TrackingPerformance.root", workPath.Data ()), "recreate");


  for (int iSys : systems) {

    for (int iWP = 0; iWP < trackWPs.size (); iWP++) {

      for (int iMult = 0; iMult < nMultBins; iMult++) {

        for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {

          if (f_primary_rate[iSys][iWP][iMult][iEta])
            f_primary_rate[iSys][iWP][iMult][iEta]->Write ();

          if (f_primary_rate_fakes_p100[iSys][iWP][iMult][iEta])
            f_primary_rate_fakes_p100[iSys][iWP][iMult][iEta]->Write ();

          if (h_primary_rate[iSys][iWP][iMult][iEta])
            h_primary_rate[iSys][iWP][iMult][iEta]->Write ();

        } // end loop over iEta

  
        for (int iPID = 0; iPID < nPIDs; iPID++) {

          if (h2_efficiency[iSys][iPID][iWP][iMult])
            h2_efficiency[iSys][iPID][iWP][iMult]->Write ();

        } // end loop over iPID

      } // end loop over iMult

    } // end loop over iWP

  } // end loop over iSys


  outFile->Close ();


}

#endif
