#include <ArrayTemplates.h>
#include <Utilities.h>

void DerivepPbFlavourUnc () {

  gStyle->SetPalette (kBlackBody);
  TColor::InvertPalette();

  TFile* inFile = new TFile ("FlavorJESUncertainty_R0p4.root", "read");

  const double etaBins[] = {-4.4, -3.6, -2.8, -2.1, -1.2, -0.8, -0.3, 0, 0.3, 0.8, 1.2, 2.1, 2.8, 3.6, 4.4};
  const int nEtaBins = sizeof (etaBins) / sizeof (etaBins[0]) - 1;

  TH2D* h2FgP = (TH2D*) inFile->Get ("hFgP")->Clone ("h2FgP");
  TH2D* h2FgH = (TH2D*) inFile->Get ("hFgH")->Clone ("h2FgH");
  TH2D* h2RgP = (TH2D*) inFile->Get ("hRgP")->Clone ("h2RgP");
  TH2D* h2RqP = (TH2D*) inFile->Get ("hRqP")->Clone ("h2RqP");
  TH2D* h2RgH = (TH2D*) inFile->Get ("hRgH")->Clone ("h2RgH");
  TH2D* h2RqH = (TH2D*) inFile->Get ("hRqH")->Clone ("h2RqH");
  TH2D* h2RP = (TH2D*) inFile->Get ("hRP")->Clone ("h2RP");

  const int nPtBins = h2FgP->GetNbinsX ();
  double* pTBins = new double[nPtBins+1];
  for (int iPt = 0; iPt < nPtBins; iPt++) {
    pTBins[iPt] = h2FgP->GetXaxis ()->GetBinLowEdge (iPt+1);
  }
  pTBins[nPtBins] = h2FgP->GetXaxis ()->GetBinLowEdge (nPtBins+1);


  const float boost = -0.465;


  TFile* outFile = new TFile ("FlavorJESUncertainty_R0p4_pPb.root", "recreate");

  TH2D* h2dR  = new TH2D ("h2dR", ";#it{p}_{T} [GeV];Lab #eta", nPtBins, pTBins, 8800, -4.4, 4.4);
  TH2D* h2dF  = new TH2D ("h2dF", ";#it{p}_{T} [GeV];Lab #eta", nPtBins, pTBins, 8800, -4.4, 4.4);

  TH2D* h2dR_nb  = new TH2D ("h2dR_nb", ";#it{p}_{T} [GeV];Lab #eta", nPtBins, pTBins, 8800, -4.4, 4.4);
  TH2D* h2dF_nb  = new TH2D ("h2dF_nb", ";#it{p}_{T} [GeV];Lab #eta", nPtBins, pTBins, 8800, -4.4, 4.4);

  for (int iPt = 0; iPt < nPtBins; iPt++) {

    const double pT = 0.5 * (pTBins[iPt] + pTBins[iPt+1]);

    for (int iEta = 0; iEta < 8800; iEta++) {

      double etalo = -4.4 + 0.001*iEta;
      double etahi = -4.4 + 0.001*(iEta+1);

      const double eta = 0.5 * (etalo + etahi);
      const double etastar = eta + boost;

      if (etastar < -4.4 || etastar > 4.4) {// || pT * std::fmin (std::cosh (etalo), std::cosh (etahi)) > 5020) {
        h2dR->SetBinContent (iPt+1, iEta+1, 0);
        h2dF->SetBinContent (iPt+1, iEta+1, 0);
        continue;
      }

      // response is *always* evaluated at lab eta
      const double RgP = h2RgP->GetBinContent (h2RgP->GetXaxis ()->FindBin (pT), h2RgP->GetYaxis ()->FindBin (std::fabs (eta)));
      const double RqP = h2RqP->GetBinContent (h2RqP->GetXaxis ()->FindBin (pT), h2RqP->GetYaxis ()->FindBin (std::fabs (eta)));
      const double RgH = h2RgH->GetBinContent (h2RgH->GetXaxis ()->FindBin (pT), h2RgH->GetYaxis ()->FindBin (std::fabs (eta)));
      const double RqH = h2RqH->GetBinContent (h2RqH->GetXaxis ()->FindBin (pT), h2RqH->GetYaxis ()->FindBin (std::fabs (eta)));

      // no-boost flavour fraction is evaluated at lab eta
      const double FgP_nb = h2FgP->GetBinContent (h2FgP->GetXaxis ()->FindBin (pT), h2FgP->GetYaxis ()->FindBin (std::fabs (eta)));
      const double FgH_nb = h2FgH->GetBinContent (h2FgH->GetXaxis ()->FindBin (pT), h2FgH->GetYaxis ()->FindBin (std::fabs (eta)));

      // flavour fraction is evaluated at lab eta + boost
      const double FgP = h2FgP->GetBinContent (h2FgP->GetXaxis ()->FindBin (pT), h2FgP->GetYaxis ()->FindBin (std::fabs (etastar)));
      const double FgH = h2FgH->GetBinContent (h2FgH->GetXaxis ()->FindBin (pT), h2FgH->GetYaxis ()->FindBin (std::fabs (etastar)));


      // flavour response uncertainty is defined as FgP * (RgP - RgH)
      double dR = FgP * (RgH - RgP);
      double dR_nb = FgP_nb * (RgH - RgP);
      h2dR->SetBinContent (iPt+1, iEta+1, dR);
      h2dR_nb->SetBinContent (iPt+1, iEta+1, dR_nb);


      // flavour fraction uncertainty is defined as dFg * dRP / aRP (when defineable)
      double dRP = std::fabs (RgP - RqP);
      double dFg = FgH - FgP;
      double dFg_nb = FgH_nb - FgP_nb;
      double aRP = FgP * RgP + (1.0-FgP) * RqP;
      double aRP_nb = FgP_nb * RgP + (1.0-FgP_nb) * RqP;

      double dF = 0;
      double dF_nb = 0;
      if (aRP != 0)
        dF = dRP * dFg / aRP;
      if (aRP_nb != 0)
        dF_nb = dRP * dFg_nb / aRP_nb;

      //if (iEta == 1000 && pT > 800) {
      //  std::cout << "at eta = " << eta << ", pT = " << pT << ", we have " << std::endl;
      //  std::cout << "  RgP = " << RgP << std::endl;
      //  std::cout << "  RqP = " << RqP << std::endl;
      //  std::cout << "  RgH = " << RgH << std::endl;
      //  std::cout << "  RqH = " << RqH << std::endl;
      //  std::cout << "  FgP = " << FgP << std::endl;
      //  std::cout << "  FgH = " << FgH << std::endl;
      //  std::cout << "  dR  = " << dR  << std::endl;
      //  std::cout << "  dRP = " << dRP << std::endl;
      //  std::cout << "  dFg = " << dFg << std::endl;
      //  std::cout << "  aRP = " << aRP << std::endl;
      //  std::cout << "  dF  = " << dF  << std::endl;
      //}
 
      h2dF->SetBinContent (iPt+1, iEta+1, dF);
      h2dF_nb->SetBinContent (iPt+1, iEta+1, dF_nb);

    } // end loop over iEta

  } // end loop over iPt



  // calculate average differences in same bins as nominal flavour uncertainties
  TH2D* h2dR_avg  = new TH2D ("h2dR_avg", ";#it{p}_{T} [GeV];Lab #eta", nPtBins, pTBins, nEtaBins, etaBins);
  TH2D* h2dF_avg  = new TH2D ("h2dF_avg", ";#it{p}_{T} [GeV];Lab #eta", nPtBins, pTBins, nEtaBins, etaBins);
  TH2D* h2dR_avg_nb  = new TH2D ("h2dR_avg_nb", ";#it{p}_{T} [GeV];Lab #eta", nPtBins, pTBins, nEtaBins, etaBins);
  TH2D* h2dF_avg_nb  = new TH2D ("h2dF_avg_nb", ";#it{p}_{T} [GeV];Lab #eta", nPtBins, pTBins, nEtaBins, etaBins);


  for (int iPt = 0; iPt < nPtBins; iPt++) {

    for (int iEta = 0; iEta < nEtaBins; iEta++) {

      double dR_tot = 0;
      double dF_tot = 0;
      double dR_tot_nb = 0;
      double dF_tot_nb = 0;
      int counts = 0;

      const double etalo = etaBins[iEta];
      const double etahi = etaBins[iEta+1];

      for (int iEtaP = 0; iEtaP < 8800; iEtaP++) {

        double etaplo = -4.4 + 0.001*iEtaP;
        double etaphi = -4.4 + 0.001*(iEtaP+1);
        const double etap = 0.5 * (etaplo + etaphi);

        if (etalo < etap && etap < etahi) {

          dR_tot += h2dR->GetBinContent (iPt+1, iEtaP+1);
          dF_tot += h2dF->GetBinContent (iPt+1, iEtaP+1);
          dR_tot_nb += h2dR_nb->GetBinContent (iPt+1, iEtaP+1);
          dF_tot_nb += h2dF_nb->GetBinContent (iPt+1, iEtaP+1);

          counts++;

        }

      } // end loop over iEtaP

      h2dR_avg->SetBinContent (iPt+1, iEta+1, dR_tot / counts);
      h2dF_avg->SetBinContent (iPt+1, iEta+1, dF_tot / counts);
      h2dR_avg_nb->SetBinContent (iPt+1, iEta+1, dR_tot_nb / counts);
      h2dF_avg_nb->SetBinContent (iPt+1, iEta+1, dF_tot_nb / counts);

    } // end loop over iEta

  } // end loop over iPt


  


  TCanvas* c1 = new TCanvas ("c1", "", 1000, 800);
  c1->SetRightMargin (0.20);
  c1->SetLogx ();
  TH2D* h_ratio1 = (TH2D*) h2dR->Clone ("h_ratio1");
  h_ratio1->Add (h2dR_nb, -1);
  h_ratio1->GetZaxis ()->SetRangeUser (-0.015, 0.015);
  h_ratio1->Draw ("colz");
  myText (0.25, 0.85, kBlack, "Flavour response uncertainty", 0.04);
  myText (0.25, 0.80, kBlack, "#deltaJES_{Boost} - #deltaJES_{No boost}", 0.04);
  c1->SaveAs ("dR_ratio.pdf");


  TCanvas* c2 = new TCanvas ("c2", "", 1000, 800);
  c2->SetRightMargin (0.20);
  c2->SetLogx ();
  TH2D* h_ratio2 = (TH2D*) h2dF->Clone ("h_ratio2");
  h_ratio2->Add (h2dF_nb, -1);
  h_ratio2->GetZaxis ()->SetRangeUser (-0.005, 0.005);
  h_ratio2->Draw ("colz");
  myText (0.25, 0.85, kBlack, "Flavour fraction uncertainty", 0.04);
  myText (0.25, 0.80, kBlack, "#deltaJES_{Boost} - #deltaJES_{No boost}", 0.04);
  c2->SaveAs ("dF_ratio.pdf");


  outFile->cd ();

  h2dR->Write ("term1");
  h2dF->Write ("term2");


  outFile->Close ();

}
