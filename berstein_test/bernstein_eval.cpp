#include "/home/bumquist/Tools/ibex/OUT/include/ibex/ibex.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ginac/ginac.h>
#include <fstream>
#include <string>
#include <vector>
#include "math.h"

using namespace ibex;
using namespace std;

unsigned factorial(unsigned n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

// put an interval [a,b] into the form x*[0,1]+y, y = a and x = b-a. return (x,y)
pair<double,double> toAffine(Interval a) {
    pair<double,double> p;
    p.first = a.ub()-a.lb();
    p.second = a.lb();
    return p;
}

double binomial(unsigned a,unsigned b) {
    if(b>a) {
        cout<<"error: first argument expect to be greater than the second to compute binomial"<<endl;
        return -1;
    }
    if(a==b||b==0)
        return 1;
    else
        return double(factorial(a))/double(factorial(a-b)*factorial(b));
}



double compute_bernstein_coef(double * coef,unsigned nbvar,unsigned * bern_coefs,unsigned* degrees,unsigned maxdeg,unsigned nbcoef=0,unsigned curvar=0) {
    double res(0);
    if(curvar<=nbvar-1) {
        for(unsigned i=0;i<=bern_coefs[curvar];i++) {
            degrees[curvar] = i;
            res+=compute_bernstein_coef(coef,nbvar,bern_coefs,degrees,maxdeg,nbcoef+pow((maxdeg+1),(nbvar-1-curvar))*(maxdeg-i),curvar+1);
        }
    }
    else {
//        cout<<"nbcoef: "<<nbcoef<<endl<<"coef[nbcoef]: "<<coef[nbcoef]<<endl;
        res = coef[nbcoef];
        for(unsigned j=0;j<nbvar;j++){
//            cout<<" bern_coefs[j]: "<<bern_coefs[j]<<" degrees[j]: "<<degrees[j]<<", binomial(bern_coefs[j],degrees[j]): "<<binomial(bern_coefs[j],degrees[j])<<endl;
//            cout<<" maxdeg: "<<maxdeg<<" degrees[j]: "<<degrees[j]<<", binomial(bern_coefs[j],degrees[j]): "<<binomial(maxdeg,degrees[j])<<endl;
//            cout<<"binomial(bern_coefs[j],degrees[j])/binomial(maxdeg,degrees[j]): "<<binomial(bern_coefs[j],degrees[j])/binomial(maxdeg,degrees[j])<<endl;
            res*=binomial(bern_coefs[j],degrees[j])/binomial(maxdeg,degrees[j]);
        }
    }
    return res;
}


int main() {
    double coef[9] = {0,0,1,0,-1,-1,1,5,7};
    unsigned bern_coefs[2] = {2,2};
    unsigned nbvar(2);
    unsigned maxdeg(2);
    unsigned degrees[2] = {0,0};
    for(unsigned i=0;i<=2;i++){
        for(unsigned j=0;j<=2;j++){
            bern_coefs[0] = i;
            bern_coefs[1]=j;
            cout<<"bernstein("<<i<<","<<j<<"): "<<compute_bernstein_coef(coef,nbvar,bern_coefs,degrees,maxdeg)<<endl;
            cout<<"****************"<<endl;
        }
    }

    return 0;

}
