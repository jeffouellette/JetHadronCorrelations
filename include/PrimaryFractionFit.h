#ifndef __PrimaryFractionFit_h__
#define __PrimaryFractionFit_h__

#include <math.h>

#include <TF1.h>

typedef long double ldouble;


long factorial (long i) {
  assert (i >= 0);
  if (i == 0 || i == 1)
    return i;
  return i * factorial (i-1);
}


/**
 * Functor that returns a piecewise n^th order polynomial in log(x) that transitions to a constant above some fitted value.
 */
class PrimaryFractionFit {

  protected:
    int degree = 6; // maximum polynomial degree
    int nderiv = 1; // number of continuous derivatives at transition to constant. nderiv=1 means only continuous (first derivative can be noncontinuous)

  public:

    void SetDegree (int _degree) { degree = _degree; }
    int  Degree    ()            { return degree;    }

    void SetNDeriv (int _nderiv) { nderiv = _nderiv; }
    int  NDeriv    ()            { return nderiv;    }

    int  NDF       ()            { return degree - nderiv + 1; }

    double polyLogN (double x, double* p) {

      ldouble polyCoeffs[degree+1];
      polyCoeffs[0] = p[1];
      for (int n = nderiv+1; n <= degree; n++) {
        polyCoeffs[n] = p[n-nderiv+1];
      }

      for (int n = nderiv; n >= 1; n--) {
        polyCoeffs[n] = 0;
        const long den = factorial (n);
        for (int np = n+1; np <= degree; np++) {
          polyCoeffs[n] -= (ldouble) ((ldouble) factorial (np) / ((ldouble) factorial (np-n) * den)) * polyCoeffs[np] * std::pow (std::log (p[0]), np-n);
        }
      }

      ldouble val = polyCoeffs[0];
      for (int i = 1; i < degree; i++) {
        val += polyCoeffs[i] * std::pow (std::log (x), i);
      }

      return (double) val;
    }

    double Eval (double x, double* p) {
      // param 0: transition to constant
      // param 1-NDF: parameters of polynomial
      return polyLogN (std::fmin (x, p[0]), p);
    }

    double operator() (double* x, double* p) {
      return Eval (x[0], p);
    }
};


#endif
