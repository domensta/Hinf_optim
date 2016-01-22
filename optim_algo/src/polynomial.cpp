#include "polynomial.h"


polynomial::polynomial(string polyexp, string *aff_variables,unsigned nbvar_aff, string *variables,unsigned nbvar):nbvar(nbvar),nbvar_aff(nbvar_aff) {
    symbol * affine_variables = new symbol[nbvar_aff*2];
    symbol var[nbvar+nbvar_aff];
    for(unsigned i=0;i<nbvar_aff;i++) {
        var[i] = symbol(aff_variables[i]);
        affine_variables[i*2] = symbol("x("+to_string(i*2+1)+")");
        affine_variables[i*2+1] = symbol("x("+to_string(i*2+2)+")");
    }
    for(unsigned i=0;i<nbvar;i++)
        var[nbvar_aff+i] = symbol("x("+to_string(nbvar_aff*2+i+1)+")");
//    cout<<"variables ok"<<endl;
    symtab tab;
    for(unsigned i=0;i<nbvar_aff;i++)
        tab[aff_variables[i]] = affine_variables[i*2]*var[i]+affine_variables[i*2+1];
    for(unsigned i=0;i<nbvar;i++)
        tab[variables[i]] = var[nbvar_aff+i];
//    cout<<"tab ok"<<endl;
    GiNaC::parser reader(tab);
    ex affinexp = reader(polyexp);
//    cout<<"affinexp: "<<affinexp<<endl;
    affinexp = affinexp.expand();
    maxdegs = new unsigned[nbvar_aff];
    for(unsigned i=0;i<nbvar_aff;i++) {
        maxdegs[i] = affinexp.degree(var[i]);
        cout<<"affine exp polynomial in "<<var[i]<<"? "<<is_polynomial(affinexp,var[i])<<endl;
        cout<<"maxdeg of "<<var[i]<<": "<<maxdegs[i]<<endl;
    }
    occurences = new unsigned[nbvar_aff];
    occurences[nbvar_aff-1]=1;
    for(int i=nbvar_aff-2;i>=0;i--) {
        occurences[i] = occurences[i+1]*(maxdegs[i+1]+1);
        cout<<"occurences "<<i<<": "<<occurences[i]<<endl;
    }
    unsigned nbcoef = occurences[0]*(maxdegs[0]+1);
    cout<<"nbcoef: "<<nbcoef<<endl;
//    coefs = new Function[nbcoef];
    ex extmp;
    unsigned nbit;
    unsigned degree;
    for(unsigned i=0;i<nbcoef;i++) {
//        cout<<"current coef: "<<i<<endl;
        nbit = i;
        degree = maxdegs[0]-i/occurences[0];
        extmp = affinexp.expand().coeff(var[0],degree);
//        cout<<"var: 0, "<<"degree: "<<degree<<endl;
        nbit%=occurences[0];
//        cout<<"extmp: "<<affinexp.coeff(var[0],degree)<<endl;
        for(unsigned j=1;j<nbvar_aff;j++) {
            degree = maxdegs[j]-nbit/occurences[j];
            extmp = extmp.expand().coeff(var[j],degree);
            cout<<"var: "<<j<<", "<<"degree: "<<degree<<endl;
            nbit%=occurences[j];


        }
        cout<<"expression: "<<extmp<<endl;
        stringstream ss;
        ss<<extmp;
        string xvar = "x["+to_string(2*nbvar_aff+nbvar)+"]";
        string sexp = ss.str();
        replace(sexp.begin(),sexp.end(),'E','e');
        char * xvarch = new char[xvar.length()+1];
        strcpy(xvarch,xvar.c_str());
        char *fexp = new char[sexp.length()+1];
        strcpy(fexp,sexp.c_str());
//        cout<<"nbcoef: "<<i<<", first argument: "<<xvarch<<" second: "<<fexp<<endl;
        Function * tmp = new Function(xvarch,fexp);
        coefs.push_back(tmp);
//        cout<<"evaluation of coef "<<i<<": "<<coefs.back()->eval_vector(IntervalVector(nbvar*2,Interval(1)))<<endl;
//        cout<<"coef at "<<i<<": "<<extmp<<endl;
    }
    cout<<"end of constructor"<<endl;
}

polynomial::~polynomial() {
//    delete coefs;
    delete maxdegs;
//    delete affine_variables;
    while(!coefs.empty()) {
        delete coefs.back();
        coefs.pop_back();
    }
}

IntervalVector polynomial::eval_affine_coef(const IntervalVector& box) {
//    cout<<"affine evaluation..."<<endl;
    unsigned nbcoef = occurences[0]*(maxdegs[0]+1);
    pair<double,double> tmp;
//    double * res = new double[nbcoef];
//    exmap m;
    IntervalVector m(nbvar_aff*2+nbvar);
    IntervalVector res(nbcoef);
    for(unsigned i=0;i<nbvar_aff;i++) {
        tmp = affine_coefs(box[i]);
//        cout<<"tmp: "<<tmp.first<<","<<tmp.second<<endl;
//        m[affine_variables[i*2]] = tmp.first;
//        m[affine_variables[i*2+1]] = tmp.second;
        m[i*2] = Interval(tmp.first);
        m[i*2+1] = Interval(tmp.second);
    }
    for(unsigned i=0;i<nbvar;i++) {
        m[2*nbvar_aff+i] = box[nbvar_aff+i];
    }
//    cout<<"affine vector: "<<m<<endl;
    for(unsigned i=0;i<nbcoef;i++) {
        res[i] = (coefs.at(i)->eval_vector(m))[0];
//        cout<<"exp coef "<<i<<": "<<*(coefs.at(i))<<endl<<"eval: "<<res[i]<<endl;
    }
//    cout<<"passed"<<endl;
    return res;
}

IntervalVector polynomial::eval_bernstein(const IntervalVector& box) {
    if(box.size()!=(nbvar+nbvar_aff))
        return IntervalVector(1,Interval::ALL_REALS);
    unsigned nbcoef = occurences[0]*(maxdegs[0]+1);
//    cout<<"affine evaluation:"<<endl;
    IntervalVector ev_coef= eval_affine_coef(box);
//    cout<<"ev_coef: "<<ev_coef<<endl;
//    double * ev_coef = eval_affine_coef(box);
    unsigned bern_coefs[nbvar_aff],tmp_degrees[nbvar_aff];
    unsigned nbit;
    unsigned degree;
    IntervalVector res(1);
    double min(POS_INFINITY),max(NEG_INFINITY);
    for(unsigned i=0;i<nbcoef;i++) {
        nbit = i;
        for(unsigned j=0;j<nbvar_aff;j++) {
            degree = maxdegs[j]-nbit/occurences[j];
            bern_coefs[j] = degree;
            nbit%=occurences[j];
//            cout<<"berncoef at "<<j<<": "<<degree<<endl;
        }
        res = compute_bernstein_coef(bern_coefs,tmp_degrees,ev_coef);
//        cout<<"res: "<<res<<endl;
        if(min>res[0].lb()) min=res[0].lb();
        if(max<res[0].ub()) max=res[0].ub();
    }
//    cout<<"min: "<<min<<", max: "<<max<<endl;
    return IntervalVector(1,Interval(min,max));

}

IntervalVector polynomial::compute_bernstein_coef(unsigned * bern_coefs,unsigned* degrees,IntervalVector coef_eval,unsigned nbcoef,unsigned curvar) {
    IntervalVector res(1,Interval(0));
    if(curvar<=nbvar_aff-1) {
        for(unsigned i=0;i<=bern_coefs[curvar];i++) {
            degrees[curvar] = i;
//            cout<<"degree at "<<curvar<<": "<<i<<endl;
            res+=compute_bernstein_coef(bern_coefs,degrees,coef_eval,nbcoef+occurences[curvar]*(maxdegs[curvar]-i),curvar+1);
//            cout<<"res+compute_bern: "<<res<<endl;
        }
    }
    else {
        res[0] = coef_eval[nbcoef];
//        cout<<"res: "<<coef_eval[nbcoef]<<endl;
        for(unsigned j=0;j<nbvar_aff;j++){
//            cout<<"bern_coefs[j]: "<<bern_coefs[j]<<", degrees[j]: "<<degrees[j]<<", maxdegs[j]"<<maxdegs[j]<<endl;

            res*=binomial(bern_coefs[j],degrees[j])/binomial(maxdegs[j],degrees[j]);
        }
//        cout<<"res*binom: "<<res<<endl;
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
        cout<<"error: first argument expect to be greater than the second to compute binomial"<<endl<<": "<<a<<"<"<<b<<endl;
        return -1;
    }
    if(a==b||b==0)
        return 1;
    else
        return double(factorial(a))/double(factorial(a-b)*factorial(b));
}

