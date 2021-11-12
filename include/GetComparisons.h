#ifndef __GetATLASJetFF_h__
#define __GetATLASJetFF_h__

#include <TGraphAsymmErrors.h>

typedef TGraphAsymmErrors TGAE;


namespace JetHadronCorrelations {


TGAE** GetATLASJetFF () {

  TGAE** g_R_DpT_pPb = new TGAE*[3];

  g_R_DpT_pPb[0] = new TGAE ();
  g_R_DpT_pPb[1] = new TGAE ();
  g_R_DpT_pPb[2] = new TGAE ();

  TGAE* g = nullptr;


  {
    g = g_R_DpT_pPb[0];

    double xvals[] = {1.28575, 2.03642, 3.23818, 5.14907, 8.15397, 12.91355, 20.53669, 32.51752};
    double yvals[] = {1.03748, 1.01022, 0.99523, 0.98296, 0.98705, 0.97479, 0.93663, 0.96525};
    double yups[] = {1.15605, 1.12879, 1.10971, 1.10562, 1.11789, 1.12879, 1.11789, 1.16695};
    double ydowns[] = {0.91346, 0.88620, 0.87257, 0.85486, 0.84532, 0.81942, 0.74583, 0.76218};
    double xdivs[] = {0.99750, 1.57329, 2.49149, 3.96148, 6.29876, 9.97482, 15.79639, 25.21797, 39.61451};
    const short n = sizeof (xvals) / sizeof (xvals[0]) - 1;

    for (short i = 0; i < n; i++) {

      g->SetPoint (i, xvals[i], yvals[i]);
      g->SetPointEXlow (i, xvals[i] - xdivs[i]);
      g->SetPointEXhigh (i, xdivs[i+1] - xvals[i]);
      g->SetPointEYlow (i, yvals[i] - ydowns[i]);
      g->SetPointEYhigh (i, yups[i] - yvals[i]);

    }
  }


  {
    g = g_R_DpT_pPb[1];

    double xvals[] = {1.28892, 2.04600, 3.24817, 5.13620, 8.15338, 12.94303, 20.46714, 32.49547};
    double yvals[] = {1.09048, 1.09293, 1.07085, 1.05559, 1.05123, 1.04550, 1.02206, 0.98499};
    double yups[] = {1.16953, 1.17198, 1.13627, 1.12646, 1.13300, 1.10547, 1.10792, 1.09539};
    double ydowns[] = {1.01416, 1.01661, 1.00543, 0.98608, 0.97763, 0.98145, 0.94847, 0.88141};
    double xdivs[] = {1.00148, 1.58349, 2.49387, 3.95874, 6.28407, 9.97534, 15.77255, 25.03753, 39.58768};
    const short n = sizeof (xvals) / sizeof (xvals[0]) - 1;

    for (short i = 0; i < n; i++) {

      g->SetPoint (i, xvals[i], yvals[i]);
      g->SetPointEXlow (i, xvals[i] - xdivs[i]);
      g->SetPointEXhigh (i, xdivs[i+1] - xvals[i]);
      g->SetPointEYlow (i, yvals[i] - ydowns[i]);
      g->SetPointEYhigh (i, yups[i] - yvals[i]);

    }
  }


  {
    g = g_R_DpT_pPb[2];

    double xvals[] = {1.29084, 2.05353, 3.25204, 5.17329, 8.19260, 12.97403, 20.63981, 32.54119, 51.75679};
    double yvals[] = {1.07108, 1.05213, 1.04683, 1.04289, 1.03759, 1.03366, 1.01334, 0.98347, 1.04233};
    double yups[] = {1.10930, 1.09581, 1.08915, 1.08931, 1.08265, 1.08963, 1.08433, 1.06265, 1.18431};
    double ydowns[] = {1.00844, 1.03285, 1.00314, 0.99511, 0.98708, 0.97632, 0.93689, 0.89337, 0.89216};
    double xdivs[] = {0.99410, 1.56719, 2.51562, 3.94799, 6.30869, 9.94559, 15.89255, 25.05443, 39.67646, 63.11563};
    const short n = sizeof (xvals) / sizeof (xvals[0]) - 1;

    for (short i = 0; i < n; i++) {

      g->SetPoint (i, xvals[i], yvals[i]);
      g->SetPointEXlow (i, xvals[i] - xdivs[i]);
      g->SetPointEXhigh (i, xdivs[i+1] - xvals[i]);
      g->SetPointEYlow (i, yvals[i] - ydowns[i]);
      g->SetPointEYhigh (i, yups[i] - yvals[i]);

    }
  }

  return g_R_DpT_pPb;
}



//TGAE*** GetPythiaAngantyr () {
//
//  TFile* inFile = new TFile (Form ("%s/rootFiles/finalHists.root", std::getenv ("PYTHIA_ANGANTYR_STUDY_PATH")), "read");
//
//  TGAE*** g_trk_pt = Get2DArray <TGAE*> (nDir, 2, 5); // direction, spectra vs. ratio, pPb centrality or pp. Note only 4 pPb centralities: 60-100, 40-60, 20-40, & 0-20%.
//
//
//  for (short iDir : {0, 1, 2}) {
//
//    const TString dir = directions[iDir];
//
//    g_trk_pt[iDir][0][nZdcCentBins+2] = make_graph ((TH1D*) inFile->Get (Form ("h_trk_pt_%s_pp", dir)));
//
//    for (short iCent = 0; iCent < 4; iCent++) {
//
//      g_trk_pt[iDir][0][iCent] = make_graph ((TH1D*) inFile->Get (Form ("h_trk_pt_%s_cent%i_pPb", dir, cent)));
//      g_trk_pt[iDir][1][iCent] = make_graph ((TH1D*) inFile->Get (Form ("h_trk_pt_%s_cent%i_pPb", dir, cent)));
//    }
//
//  }
//}


} // end namespace

#endif
