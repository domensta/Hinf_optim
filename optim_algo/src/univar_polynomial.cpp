#include "univar_polynomial.h"


univar_polynomial::univar_polynomial(string polyexp, string pol_var, string *variables,unsigned nbvar):nbvar(nbvar){
    symbol var[nbvar+1];
    for(unsigned i=0;i<nbvar;i++)
        var[i] = symbol("x("+to_string(i+1)+")");
    var[nbvar] = symbol(pol_var);
//    cout<<"variables ok"<<endl;
    symtab tab;

    for(unsigned i=0;i<nbvar;i++)
        tab[variables[i]] = var[i];
    tab[pol_var] = var[nbvar];
//    cout<<"tab ok"<<endl;
    GiNaC::parser reader(tab);
    ex pexp = reader(polyexp);
//    cout<<"affinexp: "<<pexp<<endl;
    pexp = pexp.expand();
    order = pexp.degree(var[nbvar]);

    ex extmp;
    for(unsigned i=0;i<=order;i++) {

        extmp = pexp.expand().coeff(var[nbvar],i);
//        cout<<"degree: "<<order-i<<endl;
//        cout<<"expression: "<<extmp<<endl;
        stringstream ss;
        ex nd = numer_denom(extmp.expand());
        ss<<nd[0]<<"/("<<nd[1]<<")";
        string xvar = "x["+to_string(nbvar)+"]";
        string sexp = ss.str();
//        cout<<"xvar: "<<xvar<<" sexp: "<<sexp<<endl;
        replace(sexp.begin(),sexp.end(),'E','e');
        char * xvarch = new char[xvar.length()+1];
        strcpy(xvarch,xvar.c_str());
        char *fexp = new char[sexp.length()+1];
        strcpy(fexp,sexp.c_str());

//        cout<<"first argument: "<<xvarch<<" second: "<<fexp<<endl;
        Function * tmp = new Function(xvarch,fexp);
        coefs.push_back(tmp);
//        cout<<"evaluation of coef "<<i<<": "<<coefs.back()->eval_vector(IntervalVector(nbvar*2,Interval(1)))<<endl;
//        cout<<"coef at "<<i<<": "<<extmp<<endl;
    }
//    cout<<"end of constructor"<<endl;
}

univar_polynomial::univar_polynomial(vector<Function *> coefs_func) {
    order = coefs_func.size()-1;
    nbvar = coefs_func.at(0)->nb_var();
    coefs = coefs_func;
}

univar_polynomial::~univar_polynomial() {
    while(!coefs.empty()) {
        delete coefs.back();
        coefs.pop_back();
    }
}

Interval univar_polynomial::eval(const IntervalVector& box,const Interval& polbox) {
    if(box.size() != nbvar)
        return Interval::EMPTY_SET;
    IntervalMatrix pol_coefs(1,order),power(order,1);
    for(unsigned i=0;i<=order;i++) {
        pol_coefs[0][i] = coefs.at(i)->eval_vector(box)[0];
        power[i][0] = ibex::pow(polbox,int(i));
        
    }
    return (pol_coefs*power)[0][0];
}

icomplex univar_polynomial::eval(const IntervalVector& box,const icomplex& polbox) {
    if(box.size() != nbvar)
        return icomplex(Interval::EMPTY_SET,Interval::EMPTY_SET);
    icomplex res(Interval(0),Interval(0));
    icomplex power(Interval(1),Interval(0));
    for(unsigned i=0;i<=order;i++) {
        res = res+ power*coefs.at(i)->eval_vector(box)[0];
        power = power*polbox;
    }
    return res;
}

int univar_polynomial::get_roots(const IntervalVector &box,pair<icomplex,Interval> * rvec,int rvecsize) {

    if(box.size() != nbvar || rvecsize != order)
        return 0;
    IntervalVector coefs_eval(order+1);
    for(unsigned i=0;i<=order;i++){
        coefs_eval[i] = coefs.at(i)->eval_vector(box)[0];
//        cout<<"coef eval "<<i<<" : "<<coefs_eval[i]<<endl;
        if(coefs_eval[i].is_empty())
            return 0;
    }
    Vector coefs_eval_mid = coefs_eval.mid();
//    cout<<"midpoint: "<<coefs_eval_mid<<endl;
    normalize(&coefs_eval_mid);
//    cout<<"normalized vector: "<<coefs_eval_mid<< endl;
    MatrixXd compmat = buildCompanionMatrix(coefs_eval_mid);
//    cout<<"companion matrix: "<<compmat<<endl;
    VectorXcd eivals = compmat.eigenvalues();
//    cout<<"eigen values: "<<eivals<<endl;
    icomplex eigval[order];
    for(unsigned i=0;i<order;i++){
        eigval[i] = icomplex(Interval(eivals[i].real()),Interval(eivals[i].imag()));
    }
    for(unsigned i=0;i<order;i++) {
        icomplex alphanu_den(Interval(1),Interval(0));
        for(unsigned j=0;j<order;j++) {
            if(j!=i)
                alphanu_den = alphanu_den*(eigval[i]-eigval[j]);
        }
//        cout<<"alphanu den: "<<alphanu_den<<endl;
        icomplex pol_eval = eval(box,eigval[i]);
//        cout<<"pol_eval at: "<<eigval[i]<<" : "<<pol_eval<<endl;
        icomplex alphanu = pol_eval/alphanu_den;
//        cout<<"alphanu: "<<alphanu<<endl;
        icomplex rnu = alphanu/coefs_eval[order]*order/2;
        rvec[i].first = eigval[i]-rnu;
//        cout<<"center: "<<rvec[i].first <<endl;
        rvec[i].second = rnu.abs();
//        cout<<"radius: "<<rvec[i].second<<endl;

    }
    return 1;

}

MatrixXd univar_polynomial::buildCompanionMatrix(Vector midcoef) {
//    cout<<"midcoef: "<<midcoef<<endl;
    MatrixXd compmat(order,order);
    for(unsigned line=0;line<order;line++) {
        for(unsigned col=0;col<order;col++) {
//            cout<<"line: "<<line<<" col: "<<col<<" : ";
            if(line==col+1){
//                cout<<1<<endl;
                compmat(line,col) = 1;
            }
            else if(col == order-1){
//                cout<<-midcoef[line]<<endl;
                compmat(line,col)=-midcoef[line];
            }

            else{
//                cout<<0<<endl;
                compmat(line,col)=0;
            }
        }

    }
    return compmat;
}

void univar_polynomial::normalize(Vector* vec) {
    *vec = 1/((*vec)[order])*(*vec);
//    cout<<"vec: "<<*vec<<endl;
}




