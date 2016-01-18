#include "polynomial.h"


polynomial::polynomial(string polyexp, string *variables,unsigned nbvar):nbvar(nbvar),maxdeg(0) {

    affine_variables = new symbol[nbvar*2];
    symbol var[nbvar];
    for(unsigned i=0;i<nbvar;i++) {
        var[i] = symbol(variables[i]);
        affine_variables[i*2] = symbol("a_"+ variables[i]);
        affine_variables[i*2+1] = symbol("b_"+ variables[i]);
    }
    symtab tab;
    for(unsigned i=0;i<nbvar;i++)
        tab[variables[i]] = affine_variables[i*2]*var[i]+affine_variables[i*2+1];
    GiNaC::parser reader(tab);
    ex affinexp = reader(polyexp);
    affinexp = affinexp.expand();
    cout<<"affine exp: "<<affinexp<<endl;
    maxdegs = new unsigned[nbvar];
    for(unsigned i=0;i<nbvar;i++) {
        maxdegs[i] = affinexp.degree(affine_variables[i*2]);
        maxdeg = maxdeg>maxdegs[i]?maxdeg:maxdegs[i];
    }
    coefs = new ex[unsigned(std::pow(maxdeg+1,nbvar))];
    ex extmp;
    unsigned nbit;
    unsigned degree;
    for(unsigned i=0;i<unsigned(std::pow(maxdeg+1,nbvar));i++) {
        nbit = i;
        degree = maxdeg-i/unsigned(std::pow(maxdeg+1,nbvar-1));
        extmp = affinexp.expand().coeff(var[0],degree);
//        cout<<"var: 0, "<<"degree: "<<degree<<endl;
        nbit%=unsigned(std::pow(maxdeg+1,nbvar-1));
//        cout<<"extmp: "<<affinexp.coeff(var[0],degree)<<endl;
        for(unsigned j=1;j<nbvar;j++) {
            degree = maxdeg-unsigned(nbit/std::pow(maxdeg+1,nbvar-1-j));
            extmp = extmp.expand().coeff(var[j],degree);
//            cout<<"var: "<<j<<", "<<"degree: "<<degree<<endl;
            nbit%=unsigned(std::pow(maxdeg+1,nbvar-1-j));


        }
        coefs[i] = extmp;
//        cout<<"coef at "<<i<<" :"<<extmp<<endl;
    }
    cout<<"constructor end"<<endl;

}

polynomial::~polynomial() {delete coefs;delete maxdegs;delete affine_variables;}

double * polynomial::eval_affine_coef(const IntervalVector& box) {
    unsigned nbcoef = unsigned(std::pow(maxdeg+1,nbvar));
    pair<double,double> tmp;
    double * res = new double[nbcoef];
    exmap m;
    for(unsigned i=0;i<nbvar;i++) {
        tmp = affine_coefs(box[i]);
//        cout<<"tmp: "<<tmp.first<<","<<tmp.second<<endl;
        m[affine_variables[i*2]] = tmp.first;
        m[affine_variables[i*2+1]] = tmp.second;
    }
    for(unsigned i=0;i<nbcoef;i++) {
        res[i] = ex_to<numeric>(coefs[i].subs(m)).to_double();
//        cout<<"coef "<<i<<": "<<res[i];
    }
    return res;
}

IntervalVector polynomial::eval_bernstein(const IntervalVector& box) {
    if(box.size()!=nbvar)
        return IntervalVector(1,Interval::ALL_REALS);
    unsigned nbcoef = unsigned(std::pow(maxdeg+1,nbvar));
    double * ev_coef = eval_affine_coef(box);
    unsigned bern_coefs[nbvar],tmp_degrees[nbvar];
    unsigned nbit;
    unsigned degree;
    double res;
    double min(POS_INFINITY),max(NEG_INFINITY);
    for(unsigned i=0;i<nbcoef;i++) {
        nbit = i;
        for(unsigned j=0;j<nbvar;j++) {
            degree = maxdeg-unsigned(nbit/std::pow(maxdeg+1,nbvar-1-j));
            bern_coefs[j] = degree;
            nbit%=unsigned(std::pow(maxdeg+1,nbvar-1-j));
        }
        res = compute_bernstein_coef(bern_coefs,tmp_degrees,ev_coef);
        if(min>res) min=res;
        if(max<res) max=res;
    }
    delete ev_coef;
    return IntervalVector(1,Interval(min,max));

}

double polynomial::compute_bernstein_coef(unsigned * bern_coefs,unsigned* degrees,double * coef_eval,unsigned nbcoef,unsigned curvar) {
    double res(0);
    if(curvar<=nbvar-1) {
        for(unsigned i=0;i<=bern_coefs[curvar];i++) {
            degrees[curvar] = i;
            res+=compute_bernstein_coef(bern_coefs,degrees,coef_eval,nbcoef+std::pow((maxdeg+1),(nbvar-1-curvar))*(maxdeg-i),curvar+1);

        }
    }
    else {
        res = coef_eval[nbcoef];
        for(unsigned j=0;j<nbvar;j++){
            res*=binomial(bern_coefs[j],degrees[j])/binomial(maxdeg,degrees[j]);
        }
    }
    return res;
}

pair<double,double> affine_coefs(Interval a) {
    pair<double,double> p;
    p.first = a.ub()-a.lb();
    p.second = a.lb();
    return p;
}

unsigned factorial(unsigned n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
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

