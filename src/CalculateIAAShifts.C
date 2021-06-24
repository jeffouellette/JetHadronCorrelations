#ifndef __JetHadronCorrelatorCalculateIAAShifts_C__
#define __JetHadronCorrelatorCalculateIAAShifts_C__

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
#include <TString.h>
#include <TLorentzVector.h>

#include <iostream>
#include <vector>
#include <map>
#include <math.h>

using namespace JetHadronCorrelations;


TLine* l = new TLine ();


void CalculateIAAShifts (const char* tag, const char* inFileTag) {

  TFile* inFile = nullptr;

  TH1D*** h_jet_trk_dphi_iaa = Get2DArray <TH1D*> (nZdcCentBins, nPtChSelections);


  {
    TString inFileName = inFileTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/PlotDPhi_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) { 

        h_jet_trk_dphi_iaa[iCent][iPtCh] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_iaa_iCent%i_Nominal", pTChSelections[iPtCh].Data (), iCent));

      }
    }
  }



  {
   
    ofstream outFile; 
    outFile.open (Form ("%s/aux/%s_IAAOffsets.dat", workPath.Data (), tag));

    for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) { 

        TF1* f = new TF1 ("f", "1+[0]", 0.3*M_PI, 0.7*M_PI);

        h_jet_trk_dphi_iaa[iCent][iPtCh]->Fit (f, "RN0Q");

        const double offset = f->GetParameter (0);
        const double offerr = f->GetParError (0);

        outFile << std::setw (12) << tag << std::setw (12) << pTChSelections[iPtCh].Data () << std::setw (12) << Form ("%i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]) << std::setw (12) << offset << std::setw (12) << offerr << std::setw (12) << (offerr / offset) * 100 << "%" << "\n";
      }
    }
    outFile.close ();
  }


}


#endif
