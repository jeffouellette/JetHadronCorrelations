#ifndef __TrackMomentumFit_h__
#define __TrackMomentumFit_h__

#include <math.h>

#include <TF1.h>

typedef long double ldouble;

/**
 * Functor that returns a double-sided Crystal ball function with a double Gaussian core (as opposed to the usual single Gaussian). There are 9 free parameters in this functional form.
 */
class TrackMomentumFit {
  protected:
    bool doDoubleGaussian = false;

  public:

    void DoDoubleGaussian (bool doDG) {
      doDoubleGaussian = doDG;
    }
    bool DoDoubleGaussian () {
      return doDoubleGaussian;
    }

    int ndf () {
      return DoDoubleGaussian () ? 9 : 7;
    }

    virtual double Eval (double x, double* p) {
      // single gaussian mode
      // param 0: tn
      // param 1: tp
      // param 2: r0
      // param 3: A1
      // param 4: sigma1
      // param 5: cn
      // param 6: cp

      // double gaussian mode
      // param 0: tn
      // param 1: tp
      // param 2: r0
      // param 3: A1
      // param 4: A2
      // param 5: sigma1
      // param 6: sigma2
      // param 7: cn
      // param 8: cp

      const ldouble tn = p[0];
      const ldouble tp = p[1];
      const ldouble r0 = p[2];
      const ldouble a1 = p[3];
      const ldouble a2 = DoDoubleGaussian () ? p[4] : 0; // default to 0
      const ldouble s1 = DoDoubleGaussian () ? p[5] : p[4];
      const ldouble s2 = DoDoubleGaussian () ? p[6] : 1; // default to 1 to avoid divide by 0 errors
      const ldouble cn = DoDoubleGaussian () ? p[7] : p[5];
      const ldouble cp = DoDoubleGaussian () ? p[8] : p[6];

      if (x < tn) {
        const ldouble g_tn = a1 * std::exp (-0.5 * std::pow ((tn-r0)/s1, 2)) + a2 * std::exp (-0.5 * std::pow ((tn-r0)/s2, 2));
        const ldouble gp_tn = (r0-tn) * (a1 * std::exp (-0.5 * std::pow ((tn-r0)/s1, 2)) / std::pow (s1, 2) + a2 * std::exp (-0.5 * std::pow ((tn-r0)/s2, 2)) / std::pow (s2, 2));
        return  g_tn * std::exp ((gp_tn * tn) / (g_tn * cn) * (std::pow (std::fabs (x / tn), cn) - 1));
      }
      if (x > tp) {
        const ldouble g_tp = a1 * std::exp (-0.5 * std::pow ((tp-r0)/s1, 2)) + a2 * std::exp (-0.5 * std::pow ((tp-r0)/s2, 2));
        const ldouble gp_tp = (r0-tp) * (a1 * std::exp (-0.5 * std::pow ((tp-r0)/s1, 2)) / std::pow (s1, 2) + a2 * std::exp (-0.5 * std::pow ((tp-r0)/s2, 2)) / std::pow (s2, 2));
        return g_tp * std::exp ((gp_tp * tp) / (g_tp * cp) * (std::pow (std::fabs (x / tp), cp) - 1));
      }
      return a1 * std::exp ((-0.5 * std::pow ((x-r0)/s1, 2))) + a2 * std::exp ((-0.5 * std::pow ((x-r0)/s2, 2)));
    }

    void InitParams (TF1* fit, const double integral = 1e5, const double rms = 0.05) {
      if (!DoDoubleGaussian ()) {
        fit->SetParameter (0, -1.5*rms);
        fit->SetParameter (1, 1.5*rms);
        fit->SetParameter (2, 0);
        fit->SetParameter (3, integral / (0.68 * rms * std::sqrt (2*M_PI)));
        fit->SetParameter (4, rms);
        fit->SetParameter (5, 1);
        fit->SetParameter (6, 1);

        fit->SetParLimits (0, -0.5, 0);
        fit->SetParLimits (1, 0, 0.5);
        fit->SetParLimits (2, -0.1, 0.1);
        //fit->SetParLimits (3, 0., 1e16);
        fit->SetParLimits (4, 0., 1.);
        //fit->SetParLimits (5, -4, 4); 
        //fit->SetParLimits (6, -4, 4);
      }
      else {
        fit->SetParameter (0, -1.5*rms);
        fit->SetParameter (1, 1.5*rms);
        fit->SetParameter (2, 0);
        fit->SetParameter (3, integral / (0.68 * rms * std::sqrt (2*M_PI)));
        fit->SetParameter (4, 0);
        fit->SetParameter (5, rms);
        fit->SetParameter (6, 0);
        fit->SetParameter (7, 0.5);
        fit->SetParameter (8, 0.5);

        fit->SetParLimits (0, -0.5, 0);
        fit->SetParLimits (1, 0, 0.5);
        fit->SetParLimits (2, -0.1, 0.1);
        //fit->SetParLimits (3, 0., 1e16);
        //fit->SetParLimits (4, 0., 1e16);
        fit->SetParLimits (5, 0., 1.);
        fit->SetParLimits (6, 0., 1.);
        //fit->SetParLimits (7, -4, 4); 
        //fit->SetParLimits (8, -4, 4);
      }
    }

    void CopyParams (TF1* f1, TF1* f2) {
      for (int iParam = 0; iParam < ndf (); iParam++) {
        f2->SetParameter (iParam, f1->GetParameter (iParam));
        f2->SetParError (iParam, f1->GetParError (iParam));
      }
    }

    virtual double operator () (double* x, double* p) {
      return Eval (x[0], p);
    }
};


/**
 * Same as the above, but returns the logarithm instead.
 */
class LogTrackMomentumFit : public TrackMomentumFit {
  public:
    double operator() (double* x, double* p) {
      return std::log (Eval (x[0], p));
    }
};

#endif
