#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "/home/bumquist/Tools/ibex/OUT/include/ibex/ibex.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ginac/ginac.h>
#include "math.h"

using namespace GiNaC;
using namespace ibex;
using namespace std;


class polynomial {
public:
    polynomial(string polyexp,string * variables,unsigned nbvar);
    ~polynomial();
    IntervalVector eval_bernstein(const IntervalVector& box);

private:
    unsigned nbvar;
    unsigned * maxdegs;
    unsigned * occurences;
    ex *coefs;
    symbol *affine_variables;

    double * eval_affine_coef(const IntervalVector& box);
    double compute_bernstein_coef(unsigned * bern_coefs,unsigned* degrees,double * coef_eval,unsigned nbcoef=0,unsigned curvar=0);

};

inline pair<double,double> affine_coefs(Interval a);
inline unsigned factorial(unsigned n);
inline double binomial(unsigned a,unsigned b);

#endif
