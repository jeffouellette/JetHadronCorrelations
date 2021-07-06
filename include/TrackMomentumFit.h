#ifndef __TrackMomentumFit_h__
#define __TrackMomentumFit_h__

#include <math.h>

#include <TF1.h>

typedef long double ldouble;

/**
 * Functor that returns a double-sided Crystal ball function with a double Gaussian core (as opposed to the usual single Gaussian). There are 9 free parameters in this functional form.
 */
class TrackMomentumFit {
  public:

    double operator() (double* x, double* p) {
      // param 0: tn
      // param 1: tp
      // param 2: r0
      // param 3: A1
      // param 4: A2
      // param 5: sigma1
      // param 6: sigma2
      // param 7: cn
      // param 8: cp

      //assert (p[0] < p[2]);
      //assert (p[1] > p[2]);

      if (x[0] < p[0]) {
        const long double g_tn = p[3] * std::exp (-0.5 * std::pow ((ldouble)(p[0]-p[2])/p[5], 2)) + p[4] * std::exp (-0.5 * std::pow ((ldouble)(p[0]-p[2])/p[6], 2));
        const long double gp_tn = (p[2]-p[0]) * (p[3] * std::exp (-0.5 * std::pow ((ldouble) (p[0]-p[2])/p[5], 2)) / std::pow (p[5], 2) + p[4] * std::exp (-0.5 * std::pow ((ldouble) (p[0]-p[2])/p[6], 2)) / std::pow (p[6], 2));
        return g_tn * std::exp ((ldouble) (gp_tn * p[0]) / (g_tn * p[7]) * (std::pow ((ldouble) std::fabs (x[0] / p[0]), p[7]) - 1.));
      }
      if (x[0] > p[1]) {
        const long double g_tp = p[3] * std::exp (-0.5 * std::pow ((ldouble)(p[1]-p[2])/p[5], 2)) + p[4] * std::exp (-0.5 * std::pow ((ldouble) (p[1]-p[2])/p[6], 2));
        const long double gp_tp = (p[2]-p[1]) * (p[3] * std::exp (-0.5 * std::pow ((ldouble) (p[1]-p[2])/p[5], 2)) / std::pow (p[5], 2) + p[4] * std::exp (-0.5 * std::pow ((ldouble) (p[1]-p[2])/p[6], 2)) / std::pow (p[6], 2));
        return g_tp * std::exp ((ldouble) (gp_tp * p[1]) / (g_tp * p[8]) * (std::pow ((ldouble) std::fabs (x[0] / p[1]), p[8]) - 1.));
      }
      return p[3] * std::exp ((ldouble) (-0.5 * std::pow ((ldouble) (x[0]-p[2])/p[5], 2))) + p[4] * std::exp ((ldouble) (-0.5 * std::pow ((ldouble) (x[0]-p[2])/p[6], 2)));
    }

    void InitParams (TF1* fit, const double integral = 1e5, const double rms = 0.05) {
      fit->SetParameter (0, -1.5*rms);
      fit->SetParameter (1, 1.5*rms);
      fit->SetParameter (2, 0);
      fit->SetParameter (3, 0.5 * integral / (0.68 * rms * std::sqrt (2*M_PI)));
      fit->SetParameter (4, 0.5 * integral / (0.68 * rms * std::sqrt (2*M_PI)));
      fit->SetParameter (5, rms);
      fit->SetParameter (6, rms);
      fit->SetParameter (7, 1.1);
      fit->SetParameter (8, 1.1);

      //fit->SetParLimits (0, -0.5, 0);
      //fit->SetParLimits (1, 0, 0.5);
      //fit->SetParLimits (2, -0.5, 0.5);
      //fit->SetParLimits (3, 0., 1e16);
      //fit->SetParLimits (4, 0., 1e16);
      fit->SetParLimits (5, 0., 2.);
      fit->SetParLimits (6, 0., 2.);
      //fit->SetParLimits (7, -4, 4); 
      //fit->SetParLimits (8, -4, 4);
    }

    void CopyParams (TF1* f1, TF1* f2) {
      for (int iParam = 0; iParam < 9; iParam++) {
        f2->SetParameter (iParam, f1->GetParameter (iParam));
        f2->SetParError (iParam, f1->GetParError (iParam));
      }
    }
};


/**
 * Same as the above, but returns the logarithm instead. Implemented differently here to exploit some nice logarithm properties.
 */
class LogTrackMomentumFit {
  public:
    double operator() (double* x, double* p) {
      // param 0: tn
      // param 1: tp
      // param 2: r0
      // param 3: A1
      // param 4: A2
      // param 5: sigma1
      // param 6: sigma2
      // param 7: cn
      // param 8: cp

      //assert (p[0] < p[2]);
      //assert (p[1] > p[2]);

      if (x[0] < p[0]) {
        const long double g_tn = p[3] * std::exp (-0.5 * std::pow ((ldouble)(p[0]-p[2])/p[5], 2)) + p[4] * std::exp (-0.5 * std::pow ((ldouble)(p[0]-p[2])/p[6], 2));
        const long double gp_tn = (p[2]-p[0]) * (p[3] * std::exp (-0.5 * std::pow ((ldouble) (p[0]-p[2])/p[5], 2)) / std::pow (p[5], 2) + p[4] * std::exp (-0.5 * std::pow ((ldouble) (p[0]-p[2])/p[6], 2)) / std::pow (p[6], 2));
        return std::log (g_tn) + (gp_tn * p[0]) / (g_tn * p[7]) * (std::pow ((ldouble) std::fabs (x[0] / p[0]), p[7]) - 1);
      }
      if (x[0] > p[1]) {
        const long double g_tp = p[3] * std::exp (-0.5 * std::pow ((ldouble)(p[1]-p[2])/p[5], 2)) + p[4] * std::exp (-0.5 * std::pow ((ldouble) (p[1]-p[2])/p[6], 2));
        const long double gp_tp = (p[2]-p[1]) * (p[3] * std::exp (-0.5 * std::pow ((ldouble) (p[1]-p[2])/p[5], 2)) / std::pow (p[5], 2) + p[4] * std::exp (-0.5 * std::pow ((ldouble) (p[1]-p[2])/p[6], 2)) / std::pow (p[6], 2));
        return std::log (g_tp) + (gp_tp * p[1]) / (g_tp * p[8]) * (std::pow ((ldouble) std::fabs (x[0] / p[1]), p[8]) - 1);
      }
      return std::log (p[3] * std::exp ((ldouble) (-0.5 * std::pow ((ldouble) (x[0]-p[2])/p[5], 2))) + p[4] * std::exp ((ldouble) (-0.5 * std::pow ((ldouble) (x[0]-p[2])/p[6], 2))));
    }

    void InitParams (TF1* fit, const double integral = 1e5, const double rms = 0.05) {
      fit->SetParameter (0, -1.5*rms);
      fit->SetParameter (1, 1.5*rms);
      fit->SetParameter (2, 0);
      fit->SetParameter (3, 0.5 * integral / (0.68 * rms * std::sqrt (2*M_PI)));
      fit->SetParameter (4, 0.5 * integral / (0.68 * rms * std::sqrt (2*M_PI)));
      fit->SetParameter (5, rms);
      fit->SetParameter (6, rms);
      fit->SetParameter (7, 0.5);
      fit->SetParameter (8, 0.5);

      //fit->SetParLimits (0, -0.5, 0);
      //fit->SetParLimits (1, 0, 0.5);
      fit->SetParLimits (2, -0.1, 0.1);
      //fit->SetParLimits (3, 0., 1e16);
      //fit->SetParLimits (4, 0., 1e16);
      fit->SetParLimits (5, 0., 1.);
      fit->SetParLimits (6, 0., 1.);
      //fit->SetParLimits (7, -4, 4); 
      //fit->SetParLimits (8, -4, 4);
    }

    void CopyParams (TF1* f1, TF1* f2) {
      for (int iParam = 0; iParam < 9; iParam++) {
        f2->SetParameter (iParam, f1->GetParameter (iParam));
        f2->SetParError (iParam, f1->GetParError (iParam));
      }
    }
};

#endif
