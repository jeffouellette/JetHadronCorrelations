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



/**
 * Returns CMS preliminary jet FF for jets between 60-80GeV, https://cds.cern.ch/record/2030077/files/HIN-15-004-pas.pdf
 */
TGAE** GetCMSJetFF () {

  TGAE** g_R_DpT_pPb = new TGAE*[2];

  g_R_DpT_pPb[0] = new TGAE ();
  g_R_DpT_pPb[1] = new TGAE ();

  TGAE* g = nullptr;


  {
    double xvals[]  = {0.56679, 0.69944, 0.90020, 1.20835, 1.61348, 2.15443, 2.87676, 3.84127, 5.15619, 6.88493, 9.19328, 12.27555, 16.39124, 21.88682, 29.37897, 38.81871, 52.38152, 69.94374};
    double yvals[]  = {1.05867, 1.07439, 1.06505, 1.04418, 1.03295, 1.01594, 1.02016, 1.01473, 1.01317, 1.01354, 0.97914, 1.03164,  0.99725,  0.97251,  1.03467,  0.99447,  0.96202,  1.30613};
    double yups[]   = {1.15909, 1.17095, 1.12878, 1.10597, 1.09090, 1.07580, 1.07424, 1.07074, 1.06917, 1.07147, 1.03515, 1.09344,  1.07064,  1.05362,  1.12736,  1.14898,  1.21113,  1.87775};
    double ydowns[] = {0.94474, 0.97011, 0.99360, 0.97659, 0.96729, 0.95414, 0.96030, 0.95294, 0.95138, 0.95368, 0.91735, 0.96212,  0.92966,  0.89333,  0.93424,  0.85544,  0.75538,  0.64567};
    double ystats[] = {1.08957, 1.07439, 1.06698, 1.04804, 1.03295, 1.01594, 1.02016, 1.01473, 1.01317, 1.01354, 0.97914, 1.03164,  0.99918,  0.97251,  1.03467,  1.03503,  1.04891,  1.83912};
    double xdivs[]  = {0.49961, 0.62963, 0.77291, 1.02663, 1.37806, 1.83045, 2.46998, 3.28081, 4.35782, 5.88039, 7.85194, 10.53977, 14.07349, 18.69345, 25.09247, 33.32967, 44.73886, 59.73869, 79.76759};

    const short n = sizeof (xvals) / sizeof (xvals[0]) - 1;

    g = g_R_DpT_pPb[0];
    for (short i = 0; i < n; i++) {

      g->SetPoint (i, xvals[i], yvals[i]);
      g->SetPointEXlow (i, xvals[i] - xdivs[i]);
      g->SetPointEXhigh (i, xdivs[i+1] - xvals[i]);
      g->SetPointEYlow (i, yvals[i] - ydowns[i]);
      g->SetPointEYhigh (i, yups[i] - yvals[i]);

    }

    g = g_R_DpT_pPb[1];
    for (short i = 0; i < n; i++) {

      g->SetPoint (i, xvals[i], yvals[i]);
      g->SetPointEXlow (i, xvals[i] - xdivs[i]);
      g->SetPointEXhigh (i, xdivs[i+1] - xvals[i]);
      g->SetPointEYlow (i, std::fabs (ystats[i] - yvals[i]));
      g->SetPointEYhigh (i, std::fabs (ystats[i] - yvals[i]));

    }
  }

  return g_R_DpT_pPb;
}


/**
 * Returns CMS charged hadron RpPb 
 */
TGAE** GetCMSRpPb () {

  TGAE** g_RpPb = new TGAE*[2];

  g_RpPb[0] = new TGAE ();
  g_RpPb[1] = new TGAE ();

  TGAE* g = nullptr;


  {
    double xvals[]  = {0.55162, 0.65385, 0.75189, 0.85026, 0.95390, 1.05144, 1.15624, 1.30751, 1.50845, 1.70430, 1.90137, 2.10260, 2.29928, 2.81171, 3.61571, 4.42152, 5.21174, 5.99325, 6.81534, 8.45523, 10.79645, 13.20910, 16.78089, 21.66607, 26.39137, 32.04344, 38.53031, 44.85591, 54.46242, 67.20433, 80.02889, 94.99302, 112.02820};
    double yvals[]  = {0.60881, 0.63654, 0.68462, 0.70000, 0.74615, 0.76650, 0.79727, 0.82115, 0.87611, 0.90000, 0.94038, 0.95962, 0.97692, 1.01538, 1.05769, 1.07692, 1.08654, 1.08077, 1.07611, 1.08462, 1.05739,  1.07626,  1.09517,  1.13635,  1.14859,  1.13079,  1.14080,  1.19422,  1.19645,  1.19422,  1.19200,  1.18977,  1.17753};
    double yups[]   = {0.64557, 0.67450, 0.72347, 0.74350, 0.79247, 0.81473, 0.84700, 0.88150, 0.94048, 0.97832, 1.02395, 1.05511, 1.07626, 1.12745, 1.18532, 1.20201, 1.20313, 1.18198, 1.17419, 1.17753, 1.15082,  1.16863,  1.19867,  1.24430,  1.25988,  1.23985,  1.25098,  1.30885,  1.31108,  1.30885,  1.30662,  1.30440,  1.29216};
    double ydowns[] = {0.57212, 0.59882, 0.64111, 0.65447, 0.69899, 0.71679, 0.74573, 0.76131, 0.81027, 0.82140, 0.85924, 0.86147, 0.87927, 0.90376, 0.93047, 0.95050, 0.97053, 0.97721, 0.97943, 0.99947, 0.97943,  0.99501,  0.99947,  1.03953,  1.05288,  1.03508,  1.04398,  1.09295,  1.09517,  1.09072,  1.08850,  1.08850,  1.07514};
    double ystats[] = {0.60773, 0.63666, 0.68340, 0.69899, 0.74573, 0.76798, 0.79692, 0.82140, 0.87705, 0.89931, 0.94160, 0.95718, 0.97721, 1.01505, 1.05734, 1.07737, 1.08627, 1.07959, 1.07514, 1.08405, 1.05734,  1.09295,  1.13079,  1.13524,  1.14859,  1.12856,  1.14192,  1.19311,  1.19533,  1.21314,  1.21314,  1.22650,  1.22983}; 
    double xdivs[]  = {0.50049, 0.59987, 0.70061, 0.80255, 0.90165, 1.00321, 1.10542, 1.20630, 1.39980, 1.59829, 1.80147, 1.99790, 2.20147, 2.39460, 3.20368, 3.99175, 4.78434, 5.62410, 6.42159, 7.21455, 9.65217,  11.98766, 14.41444, 19.22247, 24.02857, 28.89291, 35.19437, 41.50580, 47.85332, 60.99000, 73.81273, 86.76854, 103.99706, 120.67953};


    const short n = sizeof (xvals) / sizeof (xvals[0]) - 1;

    g = g_RpPb[0];
    for (short i = 0; i < n; i++) {

      g->SetPoint (i, xvals[i], yvals[i]);
      g->SetPointEXlow (i, xvals[i] - xdivs[i]);
      g->SetPointEXhigh (i, xdivs[i+1] - xvals[i]);
      g->SetPointEYlow (i, yvals[i] - ydowns[i]);
      g->SetPointEYhigh (i, yups[i] - yvals[i]);

    }

    g = g_RpPb[1];
    for (short i = 0; i < n; i++) {

      g->SetPoint (i, xvals[i], yvals[i]);
      g->SetPointEXlow (i, xvals[i] - xdivs[i]);
      g->SetPointEXhigh (i, xdivs[i+1] - xvals[i]);
      g->SetPointEYlow (i, std::fabs (ystats[i] - yvals[i]));
      g->SetPointEYhigh (i, std::fabs (ystats[i] - yvals[i]));

    }
  }

  return g_RpPb;
}



TGAE**** GetPythiaAngantyrIpPb () {

  TFile* inFile = new TFile (Form ("%s/rootFiles/finalHists.root", std::getenv ("PYTHIA_ANGANTYR_STUDY_PATH")), "read");


  TGAE**** g_trk_pt_ratio = Get3DArray <TGAE*> (nDir, 5, 4);
  for (int iDir = 0; iDir < nDir; iDir++) {
    const TString dir = directions[iDir];

    for (int iCent = 0; iCent < 5; iCent++) {
      int iConfig = 0;
      for (TString config : {"_allowRescatter_withNPDF", "_allowRescatter", "_withNPDF", ""}) {
        g_trk_pt_ratio[iDir][iCent][iConfig] = make_graph ((TH1D*) inFile->Get (Form ("h_trk_pt_%s_cent%i_ratio%s", dir.Data (), iCent, config.Data ())));
        if (g_trk_pt_ratio[iDir][iCent][iConfig] == nullptr)
          std::cout << "Can't find graph for " << Form ("h_trk_pt_%s_cent%i_ratio%s", dir.Data (), iCent, config.Data ()) << "??? Please check!" << std::endl;
        iConfig++;
      }
    }
  }

  //inFile->Close ();

  return g_trk_pt_ratio;
}



} // end namespace

#endif
