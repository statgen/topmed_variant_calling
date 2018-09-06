#ifndef __MATHDERIV_H__
#define __MATHDERIV_H__

#include "MathVector.h"

// Evaluates the derivative of function func() at x, using h as an initial guess
// stepsize. An estimate of the error in the derivative is stored in err.

double dfunction(double (* func)(double), double x, double h, double & err);

// Same as above, but without error estimate
//

double dfunction(double (* func)(double), double x, double h);

#endif


