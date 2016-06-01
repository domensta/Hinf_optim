#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "/home/kiwi/Tools/ibex-lib-ibex-2.1.18/OUT/include/ibex/ibex.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ginac/ginac.h>
#include "math.h"
#include <vector>

using namespace GiNaC;
using namespace ibex;
using namespace std;


class polynomial {
public:
    polynomial(string polyexp,string * aff_variables,unsigned nbvar_aff,string * variables,unsigned nbvar);
    ~polynomial();
    IntervalVector eval_bernstein(const IntervalVector& box);

private:
    unsigned nbvar;
    unsigned nbvar_aff;
    unsigned * maxdegs;
    unsigned * occurences;
    vector<Function*> coefs;
//    symbol *affine_variables;

    IntervalVector eval_affine_coef(const IntervalVector& box);
    IntervalVector compute_bernstein_coef(unsigned * bern_coefs,unsigned* degrees,IntervalVector coef_eval,unsigned nbcoef=0,unsigned curvar=0);

};

inline pair<double,double> affine_coefs(Interval a);
inline unsigned factorial(unsigned n);
inline double binomial(unsigned a,unsigned b);

#endif
