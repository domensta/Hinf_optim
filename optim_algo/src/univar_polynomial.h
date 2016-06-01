#ifndef UNIVAR_POLYNOMIAL_H
#define UNIVAR_POLYNOMIAL_H

#include "/home/kiwi/Tools/ibex-lib-ibex-2.1.18/OUT/include/ibex/ibex.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ginac/ginac.h>
#include "math.h"
#include <vector>
#include "icomplex.h"
#include <Eigen/Eigenvalues>

using namespace GiNaC;
using namespace ibex;
using namespace std;
using namespace Eigen;


class univar_polynomial {
public:
    univar_polynomial(string polyexp,string pol_var,string * variables,unsigned nbvar);
    univar_polynomial(vector<Function *> coefs_func);
    ~univar_polynomial();

    Interval eval(const IntervalVector& box,const Interval& polbox);
    icomplex eval(const IntervalVector& box,const icomplex& polbox);
    int get_roots(const IntervalVector &box,pair<icomplex,Interval> * rvec,int s);
    MatrixXd buildCompanionMatrix(Vector midcoef);
    void normalize(Vector * vec);

private:
    unsigned nbvar;
    unsigned order;
    vector<Function*> coefs;
//    symbol *affine_variables;

};



#endif
