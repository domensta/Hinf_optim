

#include "/home/bumquist/Tools/ibex/OUT/include/ibex/ibex.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "src/ibex_SetIntervalReg.h"
#include "src/polynomial.h"
#include "src/vibes.h"
#include <cstdlib>
#include <ctime>
#include <Eigen/Eigenvalues>
#include "src/icomplex.h"
#include "src/univar_polynomial.h"



using namespace std;
const double NEWTON_CEIL=5e8;
const double GAUSS_SEIDEL_RATIO=1e-04;
const bool HC4_INCREMENTAL = true;
const double PROPAG_RATIO = 0.01;
const double FIXPOINT_RATIO = 0.2;
const int SAMPLE_SIZE = 1;
const double EQ_EPS= 1.e-8;



using namespace ibex;
using namespace std;
using namespace Eigen;


class bernfunc {
public:
    vector< pair<polynomial *,polynomial *> > modulus;
    bernfunc(vector< pair<polynomial *,polynomial *> > modulus);
    IntervalVector eval_bernstein( const IntervalVector& box);
};
bernfunc::bernfunc(vector< pair<polynomial *,polynomial *> > modulus):modulus(modulus)
{}
IntervalVector bernfunc::eval_bernstein(const IntervalVector &box) {
    Interval res(0);
    for(unsigned i=0;i<modulus.size();i++)
       res += ibex::sqrt(modulus.at(i).first->eval_bernstein(box)[0])/ibex::sqrt(modulus.at(i).second->eval_bernstein(box)[0]);
    return IntervalVector(1,res);
}


class heap_elem {
public:
//    SetNodeReg * node;
    IntervalVector box;
    Interval eval;

//    heap_elem(SetNodeReg * node, IntervalVector box, Interval eval);
    heap_elem(IntervalVector box, Interval eval);
    heap_elem(const heap_elem& he);
};
//heap_elem::heap_elem(SetNodeReg * node, IntervalVector box, Interval eval): node(node),box(box),eval(eval)
//{}
heap_elem::heap_elem(IntervalVector box, Interval eval):box(box),eval(eval)
{}
heap_elem::heap_elem(const heap_elem& he):box(he.box),eval(he.eval)
{}

class list_elem{
public:
    IntervalVector box;
    IntervalVector fmax;
    Heap<heap_elem> heap;
//    list_elem(IntervalVector box,SetIntervalReg tree,IntervalVector fmax);
    list_elem(IntervalVector box,Heap<heap_elem> heap,IntervalVector fmax);
    bool computable;
//    ~list_elem();

};
//list_elem::list_elem(IntervalVector box,SetIntervalReg tree,IntervalVector fmax):box(box),
//    tree(tree),fmax(fmax)
//{}
list_elem::list_elem(IntervalVector box,Heap<heap_elem> heap,IntervalVector fmax):box(box),
    heap(heap),fmax(fmax),computable(true)
{}
//list_elem::~list_elem() {heap.flush();}

class cost_uplo: public CostFunc<list_elem> {
public:
    virtual double cost(const list_elem& elem) const;
};

double cost_uplo::cost(const list_elem& elem) const {
    return elem.fmax[0].lb();
}

class cost_box: public CostFunc<list_elem> {
public:
    virtual double cost(const list_elem& elem) const;
};

double cost_box::cost(const list_elem& elem) const {
    return elem.box.volume();
}

class optiw {
public:

    double lower_ub;
    double preclw;
    Function *f;
    vector< bernfunc > bernwset;
    vector< bernfunc > bernwfree;
    vector< bernfunc> berncst;
    bool entire;
    optiw(double lower_ub, double preclw,Function* f,vector< bernfunc> frac_wset,vector< bernfunc > frac_wfree,vector< bernfunc> frac_cst,bool entire);

};
optiw::optiw(double lower_ub, double preclw,Function* f,vector< bernfunc> frac_wset,vector< bernfunc> frac_wfree,vector< bernfunc> frac_cst,bool entire): lower_ub(lower_ub),preclw(preclw),f(f),bernwset(frac_wset),bernwfree(frac_wfree),berncst(frac_cst),entire(entire)
{}



class costf1 : public CostFunc<heap_elem> {
public:
    virtual double cost(const heap_elem& elem) const;

};
double costf1::cost(const heap_elem& elem) const {
    return -elem.eval.lb();
}
class costf2 : public CostFunc<heap_elem> {
public:
    virtual double cost(const heap_elem& elem) const;

};
double costf2::cost(const heap_elem& elem) const {
    return -elem.eval.ub();
}



pair<polynomial*,polynomial*> get_pol(Function *f,string * aff_var,unsigned nbvar_aff,string * var,unsigned nbvar) {
    symtab tab;
    for(unsigned i=0;i<nbvar_aff;i++)
        tab[aff_var[i]] = symbol(aff_var[i]);
    for(unsigned i=0;i<nbvar;i++)
        tab[var[i]] = symbol(var[i]);
    stringstream expstring;
    expstring<<f->expr();
//    cout<<"expstring: "<<expstring.str()<<endl;
    GiNaC::parser reader(tab);
    ex twz1ex = reader(expstring);
    ex twz1exnum = numer(twz1ex);
    stringstream ss;
    ss<<twz1exnum;
//    cout<<"numerator expression: "<<ss.str()<<endl;
    polynomial *num = new polynomial(ss.str(),aff_var,nbvar_aff,var,nbvar);
//    cout<<"numerator created"<<endl;
    ex twz1exden = denom(twz1ex);
//    cout<<"get denominator"<<endl;
    stringstream ssden;
    ssden<<twz1exden;
//    cout<<"denominator expression: "<<ssden.str()<<endl;
    polynomial *den = new polynomial(ssden.str(),aff_var,nbvar_aff,var,nbvar);
//    cout<<"denominator created"<<endl;
    pair<polynomial*,polynomial*> frac;
    frac.first=num;
    frac.second=den;
    return frac;

}

pair<polynomial*,polynomial*> real_part_pol(string sexp,string *var,unsigned nbvar) {
    symtab tab;
    for(unsigned i=0;i<nbvar;i++){
        tab[var[i]] = realsymbol(var[i]);
        cout<<"var "<<i<<" "<<tab[var[i]]<<endl;
    }
    realsymbol w("w");
    ex s = I*w;
//    symbol s;
    tab["s"] = s;
    GiNaC::parser reader(tab);

    string allvar[nbvar+1];
    for(unsigned i=0;i<nbvar;i++)
        allvar[i] = var[i];
    allvar[nbvar] = "w";
    ex tmp = reader(sexp);
    cout<<"ex: "<<tmp<<endl;
    for(unsigned i=0;i<tmp.degree(s);i++) {
        cout<<"degree "<<i<<": "<<tmp.coeff(s,i)<<endl;
    }
//    cout<<"numer: "<<numer(tmp.expand())<<endl;
    ex exnum = numer(imag_part(tmp.expand()));
    cout<<"real part: "<<exnum<<endl;
    ex exden = denom(tmp.expand());
    cout<<"denom: "<<exden<<endl;
    stringstream ssnum,ssden;
    ssnum<<exnum;
    ssden<<exden;
    polynomial * num = new polynomial(ssnum.str(),allvar,nbvar+1,allvar,0);
    polynomial * den = new polynomial(ssden.str(),allvar,nbvar+1,allvar,0);
    pair<polynomial*,polynomial*> frac;
    frac.first=num;
    frac.second=den;
    return frac;

}

Interval cutoff_bound(vector<Function> pol_coef, const IntervalVector& box) {
    Interval res(0);
    for(unsigned i=0;i<pol_coef.size()-1;i++){
        res = ibex::max(res,pol_coef.at(i).eval_vector(box)[0]);
    }
    res = 1+res/pol_coef.back().eval_vector(box)[0];
    return res;
}

double best_upper_bound_forall_w(optiw * str,list_elem * elem,bool midp,Vector midpt) {

//    costf2 b;
//    Heap<heap_elem> heap(b); // want to increase the uplo fast to reject w
    stack<heap_elem*> save_stack; // save element

    //list.push(*str.Inilw);
    IntervalVector res(1,Interval::EMPTY_SET),resmid(1,Interval::EMPTY_SET);
    heap_elem *elemtmp;
    double upper_ub = NEG_INFINITY;
    double lower_ub = elem->fmax[0].lb();
    LargestFirst lf(0,0.5);

    Variable kp,ki,kd,tf,u;
    NumConstraint infub(kp,ki,kd,tf,u,(*(str->f))(kp,ki,kd,tf,u)<=1000000);
    CtcFwdBwd ctc1(infub);


    IntervalVector box(elem->box.size()+1);
    for(unsigned i=0;i<elem->box.size();i++){
        if(midp)
            box[i] = midpt[i];
        else
            box[i] = elem->box[i];
    }
    double nbit(0),nbitmax;
    if(midp)
        nbitmax = 1000000000;
    else
        nbitmax = 10;

    while(!elem->heap.empty()) {
//        cout<<"lower_ub: "<<lower_ub<<endl;

//    cout<<"nb box to treat: "<<count<<endl;
//    for(unsigned i=0;i<nbitmax;i++){
//        if(heap.empty()) {

//        }
//        cout<<"lower ub: "<<lower_ub<<endl;
        elemtmp = elem->heap.pop();
//        cout<<"box: "<<elemtmp->box<<" nb it: "<<nbit<<endl;
//        cout<<"box: "<<elemtmp->box<<endl;
        box[elem->box.size()] = elemtmp->box[0];

//        box[elem->box.size()] = elemtmp->box.mid()[0];


//            cout<<"first res mid: "<<resmid<<" for box: "<<box<<endl;
        if(midp){
            //********box eval********
            IntervalVector cres = str->f->eval_vector(box);
            IntervalVector ctcbox(box);
//            if (cres[0].lb()<str->lower_ub)
//            {
//                ctc1.contract(box);
//            }
            if(!ctcbox.is_empty())
            {
                //            cout<<"res for box "<<box<<" : "<<res<<endl;
//                IntervalVector round_box(elem->box.size()+1);
//                round_box[0] = ibex::pow(Interval(10),elemtmp->box[0]);
//                for(unsigned i=1;i<=elem->box.size();i++) {
//                    round_box[i] = box[i-1];
//                }
//                res[0] =str->bernwfree.at(0).eval_bernstein(round_box)[0]*str->berncst.at(0).eval_bernstein(IntervalVector(1,round_box[0]))[0];
//                for(unsigned i=1;i<str->bernwfree.size();i++)
//                    res[0] = max(res[0],str->bernwfree.at(i).eval_bernstein(round_box)[0]*str->berncst.at(i).eval_bernstein(IntervalVector(1,round_box[0]))[0]);
//                ////            cout<<"bernstein res for box "<<round_box<<" : "<<res<<endl;
//                //            res.operator &=(str->f->eval_vector(box));

//                if(res.is_disjoint(cres)) {
//                    cout<<"ERROR: disjoint classic and bern eval"<<endl;
//                    cout<<"cres: "<<cres<<" , res: "<<res<<endl;
//                }
//                else
//                    res.operator &=(cres);
                res = cres;

                //****** mid w eval *******
                box[elem->box.size()] = elemtmp->box.mid()[0];
                resmid = str->f->eval_vector(box);
            }
            else{
                cout<<"w failure midpt"<<endl;
                resmid[0]=Interval(1000);
            }
//            if(!resmid[0].is_subset(res[0]))
//            resmid = res.lb();
//            cout<<"eval ok"<<endl;
//            cout<<"res: "<<res<<", resmid: "<<resmid<<endl;
//            cout<<"freq: "<<elemtmp->box<<" res: "<<res<<", resmid: "<<resmid<<endl;

        }
        else
        {
//            cout<<"freq: "<<elemtmp->box<<endl;
            //********box eval********
            res = str->f->eval_vector(box);
            IntervalVector ctcbox(box);
//            cout<<"res: "<<res<<endl;

            if (res[0].lb()<str->lower_ub)
            {
//                cout<<"ctcbox: "<<ctcbox;
                ctc1.contract(ctcbox);
//                cout<<" after contraction: "<<ctcbox<<endl;
            }
            if(ctcbox!=box)
                cout<<"contraction hit! box: "<<box<<" ctcbox: "<<ctcbox<<endl;
            if(!ctcbox.is_empty())
            {
                //****** mid w eval *******
                box[elem->box.size()] = elemtmp->box.mid()[0];
                IntervalVector cresmid = str->f->eval_vector(box);
////                //            cout<<"resmid for box "<<box<<" : "<<cresmid[0]<<endl;
//                box[elem->box.size()] = ibex::pow(Interval(10),box[elem->box.size()]);
//                resmid[0] =str->bernwset.at(0).eval_bernstein(box)[0]*str->berncst.at(0).eval_bernstein(IntervalVector(1,box[elem->box.size()]))[0];
//                for(unsigned i=1;i<str->bernwset.size();i++)
//                    resmid[0] = max(resmid[0],str->bernwset.at(i).eval_bernstein(box)[0]*str->berncst.at(i).eval_bernstein(IntervalVector(1,box[elem->box.size()]))[0]);
//                //           cout<<"bern resmid for box "<<box<<": "<<resmid<<endl;

//                if(resmid.is_disjoint(cresmid)) {
//                    cout<<"ERROR: disjoint classic and bern eval"<<endl;
//                    cout<<"cresmid: "<<cresmid<<" , resmid: "<<resmid<<endl;
//                }
//                else
//                    resmid.operator &=(cresmid);
                resmid=cresmid;
            }
            else{
                cout<<"w failure"<<endl;
                resmid[0]=Interval(10000);
            }

        }
//        cout<<"resmid: "<<resmid<<endl;
        if(resmid[0].lb()>str->lower_ub) {
//            cout<<"stop slave because w failure"<<endl;
            if(!midp){
                elem->fmax = IntervalVector(1,Interval::EMPTY_SET); // box K does not minimize the maximum
                delete elemtmp;
                elem->heap.flush();
                return 0;
            }
            else{
                delete elemtmp;
                elem->heap.flush();
//                cout<<"midpoint failure"<<endl;
                return POS_INFINITY;
            }
        }

        if(res[0].ub()<lower_ub) { // this box w does not contains the maximum for the box K
            delete elemtmp;
            continue;
        }

        //elem->fmax = res[0].lb()>elem->fmax[0].lb() ? IntervalVector(1,Interval(res[0].lb(),elem->fmax[0].ub())): elem->fmax;  // increase upper_lb
        elemtmp->eval = res[0];

        lower_ub = resmid[0].lb()>lower_ub? resmid[0].lb(): lower_ub;
        if(elemtmp->box[0].diam() > str->preclw && nbit<nbitmax) {
            nbit++;
//                cout<<"box: "<<elemtmp->box<<", flat? "<<elemtmp->box.is_flat()<<endl;
            pair<IntervalVector,IntervalVector> boxes= lf.bisect(elemtmp->box);
            elem->heap.push(new heap_elem(boxes.first,elemtmp->eval));
            elem->heap.push(new heap_elem(boxes.second,elemtmp->eval));
            delete elemtmp;
        }
        else {
                save_stack.push(elemtmp);
            }

    }

    while(!save_stack.empty()) {
        elem->heap.push(save_stack.top());
        save_stack.pop();
    }
    upper_ub = elem->heap.top()->eval.ub();

    if(!midp)
        elem->fmax = IntervalVector(1,Interval(lower_ub,upper_ub));
    else{
        elem->heap.flush();

    }

    return upper_ub;
}

void visit_leaves(SetIntervalReg& iset) {

    IntervalVector drawb(2);
    vector<SetNodeReg*> leaves = iset.getLeaf();
    double invol(0),maybevol(0);
    while(!leaves.empty()) {
        IntervalVector box = iset.findNodeBox(leaves.back());
        drawb[0] = box[0];
        drawb[1] = box[1];
        if(leaves.back()->status == __IBEX_IN__)
            vibes::drawBox(drawb,"[blue]");
        else if(leaves.back()->status == __IBEX_UNK__)
            vibes::drawBox(drawb,"[yellow]");
        else
            vibes::drawBox(drawb,"[red]");
        leaves.pop_back();
    }
}

Vector randpt(IntervalVector box, int size) {
    Vector pt(size);
    for(unsigned i=0;i<size;i++) {
        int rnd = std::rand();
        pt[i] = double(rnd)/double(RAND_MAX)*box[i].diam()+box[i].lb();
    }
    return pt;
}

int main() {

    std::srand(std::time(0));

    IntervalVector Iniprob(2,Interval(-10,10));
    IntervalVector Iniprob2(1,Interval(-0.1,0.1));
    // precision:
    double prec(0.1);
    double preclwmin = 1.e-4;// dynamic initialization in loop
    double stop_criterion(1); // stop if distance between uplo and loup lower than stop_criterion

    // init boxes

    IntervalVector IniboxK(4,Interval(-10,10));
    IniboxK[3] = Interval(1,5);
    IntervalVector Inilw(1,Interval(-3,3));


    //function definition
    Variable kp,ki,kd,tf,w,u;

    // pid, tf set to 1
//    Function Twz1("kp","ki","kd","u","(sqrt(((-0.1464*exp(ln(10)*4*u)+33.75*exp(ln(10)*6*u))^2+(-33.7524*exp(ln(10)*5*u)+0.14399999999999999999*exp(ln(10)*3*u))^2)*1/(((0.00366)*exp(ln(10)*2*u)*kp+0.0036*kd*exp(ln(10)*2*u)-0.1464*exp(ln(10)*4*u)+33.75*exp(ln(10)*6*u)-0.0036*ki-exp(ln(10)*4*u)*kp-kd*exp(ln(10)*4*u)+1.00006*exp(ln(10)*2*u)*ki)^2+(-0.00366*exp(ln(10)*u)*ki+(5.9999999999999999995*10^(-5))*kd*exp(ln(10)*3*u)-33.7524*exp(ln(10)*5*u)+1.00006*exp(ln(10)*3*u)*kp+0.14399999999999999999*exp(ln(10)*3*u)+exp(ln(10)*3*u)*ki-0.0036*exp(ln(10)*u)*kp)^2)) + sqrt((((0.036)*exp(ln(10)*u)-10.0006*exp(ln(10)*3*u))^2+(10*exp(ln(10)*4*u)-0.0366*exp(ln(10)*2*u))^2)*1/(((0.00366)*exp(ln(10)*2*u)*kp+0.0036*kd*exp(ln(10)*2*u)-0.1464*exp(ln(10)*4*u)+33.75*exp(ln(10)*6*u)-0.0036*ki-exp(ln(10)*4*u)*kp-kd*exp(ln(10)*4*u)+1.00006*exp(ln(10)*2*u)*ki)^2+(-0.00366*exp(ln(10)*u)*ki+(5.9999999999999999995*10^(-5))*kd*exp(ln(10)*3*u)-33.7524*exp(ln(10)*5*u)+1.00006*exp(ln(10)*3*u)*kp+0.14399999999999999999*exp(ln(10)*3*u)+exp(ln(10)*3*u)*ki-0.0036*exp(ln(10)*u)*kp)^2)))*sqrt(((-0.01248+exp(ln(10)*2*u))^2+0.024964*exp(ln(10)*2*u))*1/((0.001248-3.162*exp(ln(10)*2*u))^2+0.0078960996*exp(ln(10)*2*u)))");
//    Function Twz2("kp","ki","kd","u","(sqrt(((-33.7524*exp(ln(10)*4*u)*ki-0.1464*exp(ln(10)*4*u)*kp+33.75*kd*exp(ln(10)*6*u)-0.14399999999999999999*kd*exp(ln(10)*4*u)+0.14399999999999999999*exp(ln(10)*2*u)*ki+33.75*exp(ln(10)*6*u)*kp)^2+((0.14399999999999999999)*exp(ln(10)*3*u)*kp-33.75*exp(ln(10)*5*u)*ki-33.7524*exp(ln(10)*5*u)*kp+0.1464*exp(ln(10)*3*u)*ki-0.0023999999999999999998*kd*exp(ln(10)*5*u))^2)*1/(((0.00366)*exp(ln(10)*2*u)*kp+0.0036*kd*exp(ln(10)*2*u)-0.1464*exp(ln(10)*4*u)+33.75*exp(ln(10)*6*u)-0.0036*ki-exp(ln(10)*4*u)*kp-kd*exp(ln(10)*4*u)+1.00006*exp(ln(10)*2*u)*ki)^2+(-0.00366*exp(ln(10)*u)*ki+(5.9999999999999999995*10^(-5))*kd*exp(ln(10)*3*u)-33.7524*exp(ln(10)*5*u)+1.00006*exp(ln(10)*3*u)*kp+0.14399999999999999999*exp(ln(10)*3*u)+exp(ln(10)*3*u)*ki-0.0036*exp(ln(10)*u)*kp)^2))+sqrt(((-0.0366*exp(ln(10)*2*u)*kp-0.036*kd*exp(ln(10)*2*u)+0.036*ki+10*exp(ln(10)*4*u)*kp+10*kd*exp(ln(10)*4*u)-10.0006*exp(ln(10)*2*u)*ki)^2+((0.0366)*exp(ln(10)*u)*ki-(5.9999999999999999995*10^(-4))*kd*exp(ln(10)*3*u)-10.0006*exp(ln(10)*3*u)*kp-10*exp(ln(10)*3*u)*ki+0.036*exp(ln(10)*u)*kp)^2)*1/(((0.00366)*exp(ln(10)*2*u)*kp+0.0036*kd*exp(ln(10)*2*u)-0.1464*exp(ln(10)*4*u)+33.75*exp(ln(10)*6*u)-0.0036*ki-exp(ln(10)*4*u)*kp-kd*exp(ln(10)*4*u)+1.00006*exp(ln(10)*2*u)*ki)^2+(-0.00366*exp(ln(10)*u)*ki+(5.9999999999999999995*10^(-5))*kd*exp(ln(10)*3*u)-33.7524*exp(ln(10)*5*u)+1.00006*exp(ln(10)*3*u)*kp+0.14399999999999999999*exp(ln(10)*3*u)+exp(ln(10)*3*u)*ki-0.0036*exp(ln(10)*u)*kp)^2)))*sqrt(1/((7896.0996)*exp(ln(10)*2*u)+(-3948+exp(ln(10)*2*u))^2)*((78.960996)*exp(ln(10)*2*u)+(-39.48+exp(ln(10)*2*u))^2))");
//    Function Tw1z1w("kp","ki","kd","w","1/((-(0.0036)*w*kp+(0.14399999999999999999)*w^3+w^3*ki-(33.7524)*w^5+(1.00006)*w^3*kp-(0.00366)*w*ki+(5.9999999999999999995e-5)*kd*w^3)^2+(kd*w^4-(1.00006)*w^2*ki+w^4*kp+(0.0036)*ki+(0.1464)*w^4-(33.75)*w^6-(0.00366)*w^2*kp-(0.0036)*kd*w^2)^2)*((-(0.1464)*w^4+(33.75)*w^6)^2+((0.14399999999999999999)*w^3-(33.7524)*w^5)^2)");
//    Function Tw2z1w("kp","ki","kd","w","((-(0.0366)*w^2+10*w^4)^2+((0.036)*w-(10.0006)*w^3)^2)*1/((-(0.0036)*w*kp+(0.14399999999999999999)*w^3+w^3*ki-(33.7524)*w^5+(1.00006)*w^3*kp-(0.00366)*w*ki+(5.9999999999999999995e-5)*kd*w^3)^2+(kd*w^4-(1.00006)*w^2*ki+w^4*kp+(0.0036)*ki+(0.1464)*w^4-(33.75)*w^6-(0.00366)*w^2*kp-(0.0036)*kd*w^2)^2)");
//    Function Tw1z2w("kp","ki","kd","w","(((33.75)*kd*w^6+(33.75)*w^6*kp+(0.14399999999999999999)*w^2*ki-(0.14399999999999999999)*kd*w^4-(0.1464)*w^4*kp-(33.7524)*w^4*ki)^2+(-(0.0023999999999999999998)*kd*w^5-(33.7524)*w^5*kp+(0.1464)*w^3*ki+(0.14399999999999999999)*w^3*kp-(33.75)*w^5*ki)^2)*1/((-(0.0036)*w*kp+w^3*ki-(33.7524)*w^5+(1.00006)*w^3*kp+(0.14399999999999999999)*w^3+(5.9999999999999999995e-5)*kd*w^3-(0.00366)*w*ki)^2+((1.00006)*w^2*ki-kd*w^4+(33.75)*w^6-(0.1464)*w^4-w^4*kp-(0.0036)*ki+(0.00366)*w^2*kp+(0.0036)*kd*w^2)^2)");
//    Function Tw2z2w("kp","ki","kd","w","1/(((0.14399999999999999999)*w^3+w^3*ki-(0.0036)*w*kp-(0.00366)*w*ki+(5.9999999999999999995e-5)*kd*w^3-(33.7524)*w^5+(1.00006)*w^3*kp)^2+(w^4*kp-(1.00006)*w^2*ki+kd*w^4+(0.0036)*ki-(0.00366)*w^2*kp-(0.0036)*kd*w^2+(0.1464)*w^4-(33.75)*w^6)^2)*((10*w^4*kp-(10.0006)*w^2*ki+10*kd*w^4+(0.036)*ki-(0.0366)*w^2*kp-(0.036)*kd*w^2)^2+(10*w^3*ki-(0.036)*w*kp-(0.0366)*w*ki+(5.9999999999999999995e-4)*kd*w^3+(10.0006)*w^3*kp)^2)");

    //pid tf free
    Function Tw1z1w("kp","ki","kd","tf","w","1/((kd*w^4-(33.75)*w^6*tf+w^4*tf*kp-(5.9999999999999999995e-5)*w^2*kp+(0.14399999999999999999)*w^4*tf-(0.0036)*w^2*tf*kp-w^2*ki+(0.0036)*ki-(5.9999999999999999995e-5)*w^2*ki*tf+(0.0023999999999999999998)*w^4-(0.0036)*kd*w^2)^2+(w^3*ki*tf+(0.14399999999999999999)*w^3-(0.0036)*w*ki*tf+w^3*kp-(0.0023999999999999999998)*w^5*tf-(5.9999999999999999995e-5)*w*ki-(33.75)*w^5-(0.0036)*w*kp+(5.9999999999999999995e-5)*w^3*tf*kp+(5.9999999999999999995e-5)*kd*w^3)^2)*(((33.75)*w^6*tf-(0.14399999999999999999)*w^4*tf-(0.0023999999999999999998)*w^4)^2+((0.14399999999999999999)*w^3-(0.0023999999999999999998)*w^5*tf-(33.75)*w^5)^2)");
    Function Tw2z1w("kp","ki","kd","tf","w","((10*w^4*tf-(5.9999999999999999995e-4)*w^2-(0.036)*w^2*tf)^2+(-(0.036)*w+10*w^3+(5.9999999999999999995e-4)*w^3*tf)^2)*1/((kd*w^4-(33.75)*w^6*tf+w^4*tf*kp-(5.9999999999999999995e-5)*w^2*kp+(0.14399999999999999999)*w^4*tf-(0.0036)*w^2*tf*kp-w^2*ki+(0.0036)*ki-(5.9999999999999999995e-5)*w^2*ki*tf+(0.0023999999999999999998)*w^4-(0.0036)*kd*w^2)^2+(w^3*ki*tf+(0.14399999999999999999)*w^3-(0.0036)*w*ki*tf+w^3*kp-(0.0023999999999999999998)*w^5*tf-(5.9999999999999999995e-5)*w*ki-(33.75)*w^5-(0.0036)*w*kp+(5.9999999999999999995e-5)*w^3*tf*kp+(5.9999999999999999995e-5)*kd*w^3)^2)");
    Function Tw1z2w("kp","ki","kd","tf","w","1/((kd*w^4-(33.75)*w^6*tf+w^4*tf*kp-(5.9999999999999999995e-5)*w^2*kp+(0.14399999999999999999)*w^4*tf-(0.0036)*w^2*tf*kp-w^2*ki+(0.0036)*ki-(5.9999999999999999995e-5)*w^2*ki*tf+(0.0023999999999999999998)*w^4-(0.0036)*kd*w^2)^2+(w^3*ki*tf+(0.14399999999999999999)*w^3-(0.0036)*w*ki*tf+w^3*kp-(0.0023999999999999999998)*w^5*tf-(5.9999999999999999995e-5)*w*ki-(33.75)*w^5-(0.0036)*w*kp+(5.9999999999999999995e-5)*w^3*tf*kp+(5.9999999999999999995e-5)*kd*w^3)^2)*((-(0.14399999999999999999)*kd*w^4-(0.14399999999999999999)*w^4*tf*kp+(33.75)*kd*w^6+(0.14399999999999999999)*w^2*ki-(0.0023999999999999999998)*w^4*kp+(33.75)*w^6*tf*kp-(33.75)*w^4*ki-(0.0023999999999999999998)*w^4*ki*tf)^2+(-(0.0023999999999999999998)*kd*w^5+(0.14399999999999999999)*w^3*ki*tf+(0.0023999999999999999998)*w^3*ki-(0.0023999999999999999998)*w^5*tf*kp+(0.14399999999999999999)*w^3*kp-(33.75)*w^5*ki*tf-(33.75)*w^5*kp)^2)");
    Function Tw2z2w("kp","ki","kd","tf","w","1/((kd*w^4-(33.75)*w^6*tf+w^4*tf*kp-(5.9999999999999999995e-5)*w^2*kp+(0.14399999999999999999)*w^4*tf-(0.0036)*w^2*tf*kp-w^2*ki+(0.0036)*ki-(5.9999999999999999995e-5)*w^2*ki*tf+(0.0023999999999999999998)*w^4-(0.0036)*kd*w^2)^2+(w^3*ki*tf+(0.14399999999999999999)*w^3-(0.0036)*w*ki*tf+w^3*kp-(0.0023999999999999999998)*w^5*tf-(5.9999999999999999995e-5)*w*ki-(33.75)*w^5-(0.0036)*w*kp+(5.9999999999999999995e-5)*w^3*tf*kp+(5.9999999999999999995e-5)*kd*w^3)^2)*((10*w^3*ki*tf-(0.036)*w*ki*tf+10*w^3*kp-(5.9999999999999999995e-4)*w*ki-(0.036)*w*kp+(5.9999999999999999995e-4)*w^3*tf*kp+(5.9999999999999999995e-4)*kd*w^3)^2+(10*kd*w^4+10*w^4*tf*kp-(5.9999999999999999995e-4)*w^2*kp-(0.036)*w^2*tf*kp-10*w^2*ki+(0.036)*ki-(5.9999999999999999995e-4)*w^2*ki*tf-(0.036)*kd*w^2)^2)");
    cout<<"w ok"<<endl;
    Function Tw1z1u("kp","ki","kd","tf","u","sqrt(1/((kd*exp(ln(10)*4*u)-33.75*exp(ln(10)*6*u)*tf+exp(ln(10)*4*u)*tf*kp-(5.9999999999999999995*10^(-5))*exp(ln(10)*2*u)*kp+0.14399999999999999999*exp(ln(10)*4*u)*tf-0.0036*exp(ln(10)*2*u)*tf*kp-exp(ln(10)*2*u)*ki+0.0036*ki-(5.9999999999999999995*10^(-5))*exp(ln(10)*2*u)*ki*tf+0.0023999999999999999998*exp(ln(10)*4*u)-0.0036*kd*exp(ln(10)*2*u))^2+(exp(ln(10)*3*u)*ki*tf+0.14399999999999999999*exp(ln(10)*3*u)-0.0036*exp(ln(10)*u)*ki*tf+exp(ln(10)*3*u)*kp-0.0023999999999999999998*exp(ln(10)*5*u)*tf-(5.9999999999999999995*10^(-5))*exp(ln(10)*u)*ki-33.75*exp(ln(10)*5*u)-0.0036*exp(ln(10)*u)*kp+(5.9999999999999999995*10^(-5))*exp(ln(10)*3*u)*tf*kp+(5.9999999999999999995*10^(-5))*kd*exp(ln(10)*3*u))^2)*(((33.75)*exp(ln(10)*6*u)*tf-0.14399999999999999999*exp(ln(10)*4*u)*tf-0.0023999999999999999998*exp(ln(10)*4*u))^2+((0.14399999999999999999)*exp(ln(10)*3*u)-0.0023999999999999999998*exp(ln(10)*5*u)*tf-33.75*exp(ln(10)*5*u))^2))");
    cout<<"1";
    Function Tw2z1u("kp","ki","kd","tf","u","sqrt(((10*exp(ln(10)*4*u)*tf-(5.9999999999999999995*10^(-4))*exp(ln(10)*2*u)-0.036*exp(ln(10)*2*u)*tf)^2+(-0.036*exp(ln(10)*u)+10*exp(ln(10)*3*u)+(5.9999999999999999995*10^(-4))*exp(ln(10)*3*u)*tf)^2)*1/((kd*exp(ln(10)*4*u)-33.75*exp(ln(10)*6*u)*tf+exp(ln(10)*4*u)*tf*kp-(5.9999999999999999995*10^(-5))*exp(ln(10)*2*u)*kp+0.14399999999999999999*exp(ln(10)*4*u)*tf-0.0036*exp(ln(10)*2*u)*tf*kp-exp(ln(10)*2*u)*ki+0.0036*ki-(5.9999999999999999995*10^(-5))*exp(ln(10)*2*u)*ki*tf+0.0023999999999999999998*exp(ln(10)*4*u)-0.0036*kd*exp(ln(10)*2*u))^2+(exp(ln(10)*3*u)*ki*tf+0.14399999999999999999*exp(ln(10)*3*u)-0.0036*exp(ln(10)*u)*ki*tf+exp(ln(10)*3*u)*kp-0.0023999999999999999998*exp(ln(10)*5*u)*tf-(5.9999999999999999995*10^(-5))*exp(ln(10)*u)*ki-33.75*exp(ln(10)*5*u)-0.0036*exp(ln(10)*u)*kp+(5.9999999999999999995*10^(-5))*exp(ln(10)*3*u)*tf*kp+(5.9999999999999999995*10^(-5))*kd*exp(ln(10)*3*u))^2))");
    cout<<"2";
    Function Tw1z2u("kp","ki","kd","tf","u","sqrt(1/((kd*exp(ln(10)*4*u)-33.75*exp(ln(10)*6*u)*tf+exp(ln(10)*4*u)*tf*kp-(5.9999999999999999995*10^(-5))*exp(ln(10)*2*u)*kp+0.14399999999999999999*exp(ln(10)*4*u)*tf-0.0036*exp(ln(10)*2*u)*tf*kp-exp(ln(10)*2*u)*ki+0.0036*ki-(5.9999999999999999995*10^(-5))*exp(ln(10)*2*u)*ki*tf+0.0023999999999999999998*exp(ln(10)*4*u)-0.0036*kd*exp(ln(10)*2*u))^2+(exp(ln(10)*3*u)*ki*tf+0.14399999999999999999*exp(ln(10)*3*u)-0.0036*exp(ln(10)*u)*ki*tf+exp(ln(10)*3*u)*kp-0.0023999999999999999998*exp(ln(10)*5*u)*tf-(5.9999999999999999995*10^(-5))*exp(ln(10)*u)*ki-33.75*exp(ln(10)*5*u)-0.0036*exp(ln(10)*u)*kp+(5.9999999999999999995*10^(-5))*exp(ln(10)*3*u)*tf*kp+(5.9999999999999999995*10^(-5))*kd*exp(ln(10)*3*u))^2)*((-0.14399999999999999999*kd*exp(ln(10)*4*u)-0.14399999999999999999*exp(ln(10)*4*u)*tf*kp+33.75*kd*exp(ln(10)*6*u)+0.14399999999999999999*exp(ln(10)*2*u)*ki-0.0023999999999999999998*exp(ln(10)*4*u)*kp+33.75*exp(ln(10)*6*u)*tf*kp-33.75*exp(ln(10)*4*u)*ki-0.0023999999999999999998*exp(ln(10)*4*u)*ki*tf)^2+(-0.0023999999999999999998*kd*exp(ln(10)*5*u)+0.14399999999999999999*exp(ln(10)*3*u)*ki*tf+0.0023999999999999999998*exp(ln(10)*3*u)*ki-0.0023999999999999999998*exp(ln(10)*5*u)*tf*kp+0.14399999999999999999*exp(ln(10)*3*u)*kp-33.75*exp(ln(10)*5*u)*ki*tf-33.75*exp(ln(10)*5*u)*kp)^2))");
    cout<<"3";
    Function Tw2z2u("kp","ki","kd","tf","u","sqrt(1/((kd*exp(ln(10)*4*u)-33.75*exp(ln(10)*6*u)*tf+exp(ln(10)*4*u)*tf*kp-(5.9999999999999999995*10^(-5))*exp(ln(10)*2*u)*kp+0.14399999999999999999*exp(ln(10)*4*u)*tf-0.0036*exp(ln(10)*2*u)*tf*kp-exp(ln(10)*2*u)*ki+0.0036*ki-(5.9999999999999999995*10^(-5))*exp(ln(10)*2*u)*ki*tf+0.0023999999999999999998*exp(ln(10)*4*u)-0.0036*kd*exp(ln(10)*2*u))^2+(exp(ln(10)*3*u)*ki*tf+0.14399999999999999999*exp(ln(10)*3*u)-0.0036*exp(ln(10)*u)*ki*tf+exp(ln(10)*3*u)*kp-0.0023999999999999999998*exp(ln(10)*5*u)*tf-(5.9999999999999999995*10^(-5))*exp(ln(10)*u)*ki-33.75*exp(ln(10)*5*u)-0.0036*exp(ln(10)*u)*kp+(5.9999999999999999995*10^(-5))*exp(ln(10)*3*u)*tf*kp+(5.9999999999999999995*10^(-5))*kd*exp(ln(10)*3*u))^2)*((10*exp(ln(10)*3*u)*ki*tf-0.036*exp(ln(10)*u)*ki*tf+10*exp(ln(10)*3*u)*kp-(5.9999999999999999995*10^(-4))*exp(ln(10)*u)*ki-0.036*exp(ln(10)*u)*kp+(5.9999999999999999995*10^(-4))*exp(ln(10)*3*u)*tf*kp+(5.9999999999999999995*10^(-4))*kd*exp(ln(10)*3*u))^2+(10*kd*exp(ln(10)*4*u)+10*exp(ln(10)*4*u)*tf*kp-(5.9999999999999999995*10^(-4))*exp(ln(10)*2*u)*kp-0.036*exp(ln(10)*2*u)*tf*kp-10*exp(ln(10)*2*u)*ki+0.036*ki-(5.9999999999999999995*10^(-4))*exp(ln(10)*2*u)*ki*tf-0.036*kd*exp(ln(10)*2*u))^2))");
    cout<<"u ok"<<endl;
    Function w1("kp","ki","kd","w","1/((0.0078960996)*w^2+(0.001248-(3.162)*w^2)^2)*((0.024964)*w^2+(-0.01248+w^2)^2)");
    Function w2("kp","ki","kd","w","((78.960996)*w^2+(-39.48+w^2)^2)*1/((7896.0996)*w^2+(-3948+w^2)^2)");
    Function w1usqrt("kp","ki","kd","u","sqrt(1/((0.0078960996)*exp(ln(10)*2*u)+(0.001248-3.162*exp(ln(10)*2*u))^2)*((-0.01248+exp(ln(10)*2*u))^2+0.024964*exp(ln(10)*2*u)))");
    Function w2usqrt("kp","ki","kd","u","sqrt(1/((-3948+exp(ln(10)*2*u))^2+7896.0996*exp(ln(10)*2*u))*((-39.48+exp(ln(10)*2*u))^2+78.960996*exp(ln(10)*2*u)))");
    cout<<"cst ok"<<endl;
    Function Twz1(kp,ki,kd,tf,u,(Tw1z1u(kp,ki,kd,tf,u)+Tw2z1u(kp,ki,kd,tf,u))*w1usqrt(kp,ki,kd,u));
    Function Twz2(kp,ki,kd,tf,u,(Tw1z2u(kp,ki,kd,tf,u)+Tw2z2u(kp,ki,kd,tf,u))*w2usqrt(kp,ki,kd,u));

    Function Max12(kp,ki,kd,tf,u,ibex::max(Twz1(kp,ki,kd,tf,u),Twz2(kp,ki,kd,tf,u)));


    vector<bernfunc> bernwset;
    vector<bernfunc> bernwfree;
    vector<bernfunc> berncst;

    vector< pair<polynomial *,polynomial *> > frac_wfree;
    vector< pair<polynomial *,polynomial *> > frac_wset;
    vector< pair<polynomial *,polynomial *> > frac_cst;

    string aff_var[4];aff_var[0] = "kp";aff_var[1] = "ki";aff_var[2] = "kd";aff_var[3] = "tf";
    string var[4];var[0]="w";
    frac_wset.push_back(get_pol(&Tw1z1w,aff_var,4,var,1));
    frac_wset.push_back(get_pol(&Tw2z1w,aff_var,4,var,1));
    bernwset.push_back(bernfunc(frac_wset));
    frac_wset.clear();
    frac_wset.push_back(get_pol(&Tw1z2w,aff_var,4,var,1));
    frac_wset.push_back(get_pol(&Tw2z2w,aff_var,4,var,1));
    bernwset.push_back(bernfunc(frac_wset));
    aff_var[0] = "w";
    var[0] = "kp";var[1] = "ki";var[2] = "kd";var[3] = "tf";
    frac_wfree.push_back(get_pol(&Tw1z1w,aff_var,1,var,4));
    frac_wfree.push_back(get_pol(&Tw2z1w,aff_var,1,var,4));
    bernwfree.push_back(bernfunc(frac_wfree));
    frac_wfree.clear();
    frac_wfree.push_back(get_pol(&Tw1z2w,aff_var,1,var,4));
    frac_wfree.push_back(get_pol(&Tw2z2w,aff_var,1,var,4));
    bernwfree.push_back(bernfunc(frac_wfree));
    frac_cst.push_back(get_pol(&w1,aff_var,1,var,0));
    berncst.push_back(bernfunc(frac_cst));
    frac_cst.clear();
    frac_cst.push_back(get_pol(&w2,aff_var,1,var,0));
    berncst.push_back(bernfunc(frac_cst));
    cout<<"func ok"<<endl;


    //pid tf =1 routh function
//    Function rS1("kp","ki","kd","0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki");
//    Function rS2("kp","ki","kd","1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*((0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(33.7524)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))");
//    Function rS3("kp","ki","kd","(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*1/((0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(33.7524)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*(1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*((0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(33.7524)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(((0.0036)*kp+(0.00366)*ki)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*1/(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))");
//    Function rS4("kp","ki","kd","-1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*((0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(33.7524)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*((0.0036)*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*((0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(33.7524)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*ki-(((0.0036)*kp+(0.00366)*ki)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*1/((0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(33.7524)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*((0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(33.7524)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(((0.0036)*kp+(0.00366)*ki)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*1/(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)))*1/(1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*((0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(33.7524)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(((0.0036)*kp+(0.00366)*ki)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*1/(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))");
//    Function rS5("kp","ki","kd","(0.0036)*ki");
//    Function rGS1("kp","ki","kd","0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki");
//    Function rGS2("kp","ki","kd","((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)");
//    Function rGS3("kp","ki","kd","(1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)");
//    Function rGS4("kp","ki","kd","-1/(1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*(1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))*((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+(0.0036)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*ki*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)");
//    Function rGS5("kp","ki","kd","(0.0036)*ki");
//    Function rT1("kp","ki","kd","0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki");
//    Function rT2("kp","ki","kd","((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)");
//    Function rT3("kp","ki","kd","1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))");
//    Function rT4("kp","ki","kd","-((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))+(0.0036)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*ki)*1/(((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))");
//    Function rT5("kp","ki","kd", "(0.0036)*ki");

    //pid tf free routh function

        Function rS1("kp","ki","kd","tf","(33.75)*tf");
        Function rS2("kp","ki","kd","tf","33.75+(0.0023999999999999999998)*tf");
        Function rS3("kp","ki","kd","tf","(0.08099999999999999999+(5.760000000000053272e-6)*tf+(3.4559999999999999995e-4)*tf^2+(3.7499999999999999997e-4)*tf^2*kp+(3.7499999999999999997e-4)*tf*kd+(33.75)*kd-(33.75)*tf^2*ki)*1/(33.75+(0.0023999999999999999998)*tf)");
        Function rS4("kp","ki","kd","tf","(0.011663999999999999998+(8.294400000000076711e-7)*tf-(0.08099999999999999999)*tf*ki+(4.976639999999999999e-5)*tf^2+(2.2499999999999999996e-8)*tf^3*kp^2+(33.75)*tf*kd*ki+(0.012656249999999999999)*kp+(2.2499999999999999996e-8)*tf*kd^2+(1.0799999999999999998e-4)*tf^2*kp+(3.7499999999999999997e-4)*tf^2*kp^2+(0.0023999999999999999998)*tf*kp*kd-(5.2919965439999999993e-4)*tf*kd+(0.75937985999999999985)*kd-(0.0016499999999999999999)*tf^3*kp*ki+(4.499999999999999999e-8)*tf^2*kp*kd+(6.371996543999999999e-4)*tf^3*ki-(33.75)*tf^2*kp*ki-(1139.0625)*ki+(33.75)*kp*kd-(0.0016499999999999999999)*tf^2*kd*ki+(9.000000000000532727e-7)*tf*kp+(0.0020249999999999999998)*kd^2-(33.75)*tf^3*ki^2-(0.75937985999999999985)*tf^2*ki)*1/(0.08099999999999999999+(5.760000000000053272e-6)*tf+(3.4559999999999999995e-4)*tf^2+(3.7499999999999999997e-4)*tf^2*kp+(3.7499999999999999997e-4)*tf*kd+(33.75)*kd-(33.75)*tf^2*ki)");
        Function rS5("kp","ki","kd","tf","1/(0.011663999999999999998+(8.294400000000076711e-7)*tf-(0.08099999999999999999)*tf*ki+(4.976639999999999999e-5)*tf^2+(2.2499999999999999996e-8)*tf^3*kp^2+(33.75)*tf*kd*ki+(0.012656249999999999999)*kp+(2.2499999999999999996e-8)*tf*kd^2+(1.0799999999999999998e-4)*tf^2*kp+(3.7499999999999999997e-4)*tf^2*kp^2+(0.0023999999999999999998)*tf*kp*kd-(5.2919965439999999993e-4)*tf*kd+(0.75937985999999999985)*kd-(0.0016499999999999999999)*tf^3*kp*ki+(4.499999999999999999e-8)*tf^2*kp*kd+(6.371996543999999999e-4)*tf^3*ki-(33.75)*tf^2*kp*ki-(1139.0625)*ki+(33.75)*kp*kd-(0.0016499999999999999999)*tf^2*kd*ki+(9.000000000000532727e-7)*tf*kp+(0.0020249999999999999998)*kd^2-(33.75)*tf^3*ki^2-(0.75937985999999999985)*tf^2*ki)*((2.0229749999999999995e-8)*tf*kd^3-(0.20503124999999999999)*tf^3*kp*ki^2+(6.5611337471502335975e-6)*tf^5*ki^2-(2.6729999999999820196e-4)*tf*kp*ki-(0.14934374999999999999)*tf^4*ki^3+(2.4187500000001638245e-6)*tf^3*kd*ki^2+(1.0322634401280047733e-4)*tf*ki+(7.7760000000003164365e-12)*tf^3*kp^2-(0.034183269269999999998)*tf*kd*ki-(1.7495996759999999996e-8)*tf^5*kp^2*ki-(0.055687499999999999994)*tf^2*kp^2*ki+(0.06834374999999999999)*kp^2*kd+(0.09226817243847071999)*tf^4*ki^2+(5.831999999999999998e-13)*tf^4*kp^2*kd-(4.860000000000194026e-6)*tf^5*ki^3+(0.23034374999999999997)*tf*kp*kd*ki-(5.7736716019199999992e-5)*tf*kd^2+(1139.0625)*kp*kd*ki-(2.4603716283749999998e-4)*tf^2*kd^2*ki+(9.956249999999999998e-11)*tf^3*kp^3+(4.9207533716249999993e-4)*tf*kp*kd^2-(3.4991993519999999992e-8)*tf^4*kp*kd*ki-(4.572285014015999999e-9)*tf^3*kp*kd+(1.0935012959999999997e-7)*tf^2*kp^2-(2.3619586003200064722e-6)*tf*kp*kd+(5.831999999999999998e-13)*tf^3*kp*kd^2+(2.5628906249999999996e-5)*kp^2-(6.046624766361599999e-6)*tf^4*ki+(2.01553920000000932e-7)*tf*kd+(5.1175574496049766397e-5)*tf^3*kd*ki-(5.4674999999999999998)*tf*ki^2+(0.0014171759999999999997)*kd-(2.4603763803749999997e-4)*tf^4*kp^2*ki+(1139.0625)*tf*kd*ki^2-(6.62808821175746329e-26)*tf^3*kp+(2.4603793672499999997e-4)*tf^2*kp^2*kd+(7.5937499999999999986e-7)*tf^2*kp^3+(2.0306684707200064723e-6)*tf^3*kp*ki-(4.2998169599999999986e-10)*tf^5*ki+(6.429784199035023358e-5)*tf^4*kp*ki-(1.7495996759999999996e-8)*tf^3*kd^2*ki-(25.629100649999999996)*tf^2*ki^2-(1.8794531249999999998)*kp*ki-(6.4297841990350233585e-5)*tf^2*kp*kd-(38443.359375)*ki^2+(4.5722850140159999985e-9)*tf^4*kd*ki-(6.541874999999999999e-6)*tf^3*kp^2*ki+(25.628742225)*kd*ki+(9.719999999999999998e-6)*tf*kp^2*kd+(3.6450000000001078766e-9)*tf*kp^2+(0.16199999999999999998)*tf^2*kd*ki^2+(1.8662280560639067924e-8)*tf^3*ki-(2.375999999999835458e-10)*tf^5*kp*ki^2+(1.2138750000000247298e-5)*tf^2*kp*kd*ki+(4.1006249999999999996e-6)*kp*kd^2-(1139.0625)*tf^2*kp*ki^2+(6.0466247663615999987e-6)*tf^2*kd+(0.011948744456999999998)*tf^2*kp*ki+(0.72581023133999999987)*ki-(0.016607521408499999998)*kp*kd+(4.5722850140159999985e-9)*tf^5*kp*ki-(4.572285014015999999e-9)*tf^2*kd^2+(2.7337532399999999995e-9)*tf^4*kp^3-(0.18453326606872991997)*tf^2*kd*ki-(1.5534375163356999834e-23)*tf*kp-(4.920748008749999999e-4)*tf^3*kp*kd*ki+(2.5697256479999999993e-8)*tf^3*kp^2*kd+(0.09226465298999999998)*kd^2+(0.029524488335999999997)*tf^3*ki^2-(2.375999999999835458e-10)*tf^4*kd*ki^2-(1.0097419586828951109e-28)*tf^4*kp^2+(4.2998169599999999986e-10)*tf^3*kd+(1.9439999999999999993e-13)*tf^2*kd^3+(1.9439999999999999993e-13)*tf^5*kp^3-(1139.0625)*tf^3*ki^3+(2.4603749999999999997e-4)*kd^3-(1.21612500000000417355e-5)*tf^4*kp*ki^2+(0.0016796179906540093438)*tf^2*ki+(4.3193253239999999988e-8)*tf^2*kp*kd^2+(4.1006250000004754175e-6)*tf*kd^2*ki)*1/(33.75+(0.0023999999999999999998)*tf)");
        Function rS6("kp","ki","kd","tf","1/(0.08099999999999999999+(5.760000000000053272e-6)*tf+(3.4559999999999999995e-4)*tf^2+(3.7499999999999999997e-4)*tf^2*kp+(3.7499999999999999997e-4)*tf*kd+(33.75)*kd-(33.75)*tf^2*ki)*1/((2.0229749999999999995e-8)*tf*kd^3-(0.20503124999999999999)*tf^3*kp*ki^2+(6.5611337471502335975e-6)*tf^5*ki^2-(2.6729999999999820196e-4)*tf*kp*ki-(0.14934374999999999999)*tf^4*ki^3+(2.4187500000001638245e-6)*tf^3*kd*ki^2+(1.0322634401280047733e-4)*tf*ki+(7.7760000000003164365e-12)*tf^3*kp^2-(0.034183269269999999998)*tf*kd*ki-(1.7495996759999999996e-8)*tf^5*kp^2*ki-(0.055687499999999999994)*tf^2*kp^2*ki+(0.06834374999999999999)*kp^2*kd+(0.09226817243847071999)*tf^4*ki^2+(5.831999999999999998e-13)*tf^4*kp^2*kd-(4.860000000000194026e-6)*tf^5*ki^3+(0.23034374999999999997)*tf*kp*kd*ki-(5.7736716019199999992e-5)*tf*kd^2+(1139.0625)*kp*kd*ki-(2.4603716283749999998e-4)*tf^2*kd^2*ki+(9.956249999999999998e-11)*tf^3*kp^3+(4.9207533716249999993e-4)*tf*kp*kd^2-(3.4991993519999999992e-8)*tf^4*kp*kd*ki-(4.572285014015999999e-9)*tf^3*kp*kd+(1.0935012959999999997e-7)*tf^2*kp^2-(2.3619586003200064722e-6)*tf*kp*kd+(5.831999999999999998e-13)*tf^3*kp*kd^2+(2.5628906249999999996e-5)*kp^2-(6.046624766361599999e-6)*tf^4*ki+(2.01553920000000932e-7)*tf*kd+(5.1175574496049766397e-5)*tf^3*kd*ki-(5.4674999999999999998)*tf*ki^2+(0.0014171759999999999997)*kd-(2.4603763803749999997e-4)*tf^4*kp^2*ki+(1139.0625)*tf*kd*ki^2-(6.62808821175746329e-26)*tf^3*kp+(2.4603793672499999997e-4)*tf^2*kp^2*kd+(7.5937499999999999986e-7)*tf^2*kp^3+(2.0306684707200064723e-6)*tf^3*kp*ki-(4.2998169599999999986e-10)*tf^5*ki+(6.429784199035023358e-5)*tf^4*kp*ki-(1.7495996759999999996e-8)*tf^3*kd^2*ki-(25.629100649999999996)*tf^2*ki^2-(1.8794531249999999998)*kp*ki-(6.4297841990350233585e-5)*tf^2*kp*kd-(38443.359375)*ki^2+(4.5722850140159999985e-9)*tf^4*kd*ki-(6.541874999999999999e-6)*tf^3*kp^2*ki+(25.628742225)*kd*ki+(9.719999999999999998e-6)*tf*kp^2*kd+(3.6450000000001078766e-9)*tf*kp^2+(0.16199999999999999998)*tf^2*kd*ki^2+(1.8662280560639067924e-8)*tf^3*ki-(2.375999999999835458e-10)*tf^5*kp*ki^2+(1.2138750000000247298e-5)*tf^2*kp*kd*ki+(4.1006249999999999996e-6)*kp*kd^2-(1139.0625)*tf^2*kp*ki^2+(6.0466247663615999987e-6)*tf^2*kd+(0.011948744456999999998)*tf^2*kp*ki+(0.72581023133999999987)*ki-(0.016607521408499999998)*kp*kd+(4.5722850140159999985e-9)*tf^5*kp*ki-(4.572285014015999999e-9)*tf^2*kd^2+(2.7337532399999999995e-9)*tf^4*kp^3-(0.18453326606872991997)*tf^2*kd*ki-(1.5534375163356999834e-23)*tf*kp-(4.920748008749999999e-4)*tf^3*kp*kd*ki+(2.5697256479999999993e-8)*tf^3*kp^2*kd+(0.09226465298999999998)*kd^2+(0.029524488335999999997)*tf^3*ki^2-(2.375999999999835458e-10)*tf^4*kd*ki^2-(1.0097419586828951109e-28)*tf^4*kp^2+(4.2998169599999999986e-10)*tf^3*kd+(1.9439999999999999993e-13)*tf^2*kd^3+(1.9439999999999999993e-13)*tf^5*kp^3-(1139.0625)*tf^3*ki^3+(2.4603749999999999997e-4)*kd^3-(1.21612500000000417355e-5)*tf^4*kp*ki^2+(0.0016796179906540093438)*tf^2*ki+(4.3193253239999999988e-8)*tf^2*kp*kd^2+(4.1006250000004754175e-6)*tf*kd^2*ki)*(-(3.5307346347016243853e-6)*tf^3*kp*ki^2+(1.6386098812199999994e-13)*tf^4*kp^3*kd^2+(1.3440937499999999996e-16)*tf^5*kp^5-(5.466753942567071844e-10)*tf^8*kp*ki^3+(1.298370758537532027e-5)*tf^4*kd^2*ki^2+(1.3223982101929431184e-10)*tf^7*kp*ki^2+(0.0083037868457343749955)*tf*kp^2*kd^2*ki+(1.0497599999999999995e-18)*tf^6*kp^4*kd+(1.3839609374999999997e-4)*kp^2*kd^2*ki-(1.4366763417624766939e-9)*tf^5*ki^2+(6.470787664715795222e-11)*tf^2*kp*kd^3*ki+(3.749845247119888749e-8)*tf*kp*ki+(9.041874636937499998e-7)*tf^5*kp^3*ki^2-(0.44841490589794288335)*tf^4*ki^3-(2.2553437500000988575e-4)*tf^2*kp^2*kd*ki^2+(1.5746399999999999994e-18)*tf^6*kp^2*kd^2*ki+(0.008303765624999999999)*kp^3*kd^2-(1.1325926768653244006e-10)*tf^5*kp^2*kd*ki^2+(0.005535843749999842173)*tf^3*ki^4-(7.522977585416896509e-12)*tf^6*kp*ki+(1.0497599999999999995e-18)*tf^5*kp*kd^3*ki+(0.61240238269312499995)*kd^2*ki^2+(2.4774342718464095138e-16)*tf^5*kp^4-(7.0568715504919306673e-4)*tf^3*kd*ki^2-(3.5263873843200217417e-9)*tf*ki+(4.587400246292034065e-7)*tf^3*kp^2*kd*ki+(5.1081788334381465576e-15)*tf^7*kd*ki^2-(6.3772894876835312977e-9)*tf^4*kp^2*kd^2*ki-(9.3971249790667045395e-4)*tf^4*kp^2*ki^2-(1.0518095374818172365e-6)*tf^2*kd^3*ki^2-(3.8654753428727007862e-29)*tf^3*kp^2+(2.6243999999999999988e-19)*tf^8*kp^4*ki-(4.5917971728082037458e-8)*tf*kd*ki+(7.2559411199999999964e-19)*tf^7*kp^2*kd*ki-(4.7239195625999999982e-14)*tf^5*kd^3*ki^2-(6.8904427139675999987e-6)*tf*kp*kd^3-(7.687605808274727626e-6)*tf^2*kp^2*kd^2+(5.3496602689535999975e-16)*tf^6*kd*ki+(7.2559411199999999964e-19)*tf^5*kp^2*kd^2+(1.0251562499999999997e-12)*tf^4*kp^5+(0.041518849584937499992)*tf^2*kp*kd^2*ki^2+(4.6589167968868544213e-12)*tf^5*kp^2*ki-(3.8518598887744717061e-34)*tf^6*kp*kd*ki+(2.8231495780328956776e-7)*tf^2*kp^2*ki+(2.3066015624999999997)*tf^4*kp*ki^4-(4.842753242718599999e-6)*kp^2*kd+(3.3978252046660257948e-5)*tf^4*ki^2+(5.0929762499999999982e-14)*tf^3*kd^4*ki-(8.370195784171656298e-9)*tf^5*kp*kd^2*ki^2-(0.0011957422499999999998)*kd^3*ki-(7.1833954467028037244e-11)*tf^4*kp^2*kd+(2.3097795277355433247e-4)*tf^5*ki^3+(2.5095813905277746412e-6)*tf^4*kd^2*ki^3+(5.0929762499999999982e-14)*tf^2*kp*kd^4-(8.484815639892795983e-6)*tf^7*kp*ki^3-(9.71848124999885371e-7)*tf^3*kp*ki^3-(4.982667521585796203e-7)*tf*kp*kd*ki+(0.06576559207498689218)*kp*kd*ki+(1.899646502215679999e-14)*tf^4*kd^3*ki+(5.1081788334381465576e-15)*tf^8*kp*ki^2+(2.457383047951117093e-8)*tf^2*kd^2*ki-(6.5610000000000667786e-4)*tf^4*kd*ki^4+(6.244432598829374998e-9)*tf^2*kp^2*kd^3-(4.3548628751999999984e-14)*tf^7*kp^4*ki+(1.8895711842719999995e-12)*tf^4*kp^4-(5.108178833438146558e-15)*tf^5*kd^2*ki+(9.0699290873858264864e-15)*tf^3*kp^3+(1.4764131566860499998e-7)*tf^2*kp^3*kd+(1.14791623332020509345e-8)*tf*kp*kd^2+(2.5798925946470399988e-14)*tf^5*kp^3*kd-(5.4667539425670718437e-10)*tf^7*kd*ki^3+(8.713828124999600799e-5)*tf^4*kp^2*ki^3-(1.5923571180560192876e-7)*tf^4*kp^3*ki-(1.6989204578945295054e-6)*tf^4*kp*kd*ki+(5.296095665529092643e-6)*tf^3*kp*kd^2*ki+(6.2279700421479110973)*tf^3*kd*ki^3+(5.015307988160764567e-13)*tf^3*kp*kd+(8.740421797958495536e-10)*tf^5*kd^2*ki^2-(5.296080452307382455e-6)*tf^6*kd*ki^3-(3.221208033873737277e-31)*tf^2*kp^2+(9.050960307198425593e-4)*tf^2*kd^2*ki^2-(5.978718335878137252e-5)*tf^4*kp^3*kd*ki+(1.5746399999999999994e-18)*tf^5*kp^3*kd^2-(0.0064008193359370579856)*tf*kd*ki^3+(8.815968460800054354e-11)*tf*kp*kd+(2.9893627969931249997e-5)*tf^2*kp^3*kd^2-(5.3496602689535999975e-16)*tf^8*ki^2-(6.6430184140117709985e-10)*tf^6*kp^4*ki+(2.4186470399999999987e-19)*tf^7*kp^4-(1.14280879147617626355e-11)*tf^3*kp*kd^2+(0.033631727418750618515)*tf^4*kp*kd*ki^2-(3.985807499999952097e-5)*tf*ki^3+(1.6423003983599921879e-6)*tf^3*kp*kd^2*ki^2+(3.4763492772550099038e-8)*tf^6*ki^3+(0.008303762834296874999)*tf^4*kp^3*ki^2-(0.0020178126553904999998)*kp^2*kd^2-(0.011210690411812641075)*tf^6*kp*ki^3+(0.28025172501626953122)*tf^2*kp^2*ki^2-(3.276146064085421545e-7)*tf^5*kp^2*ki^2+(1.08058176821577553936e-10)*tf^6*kp*kd*ki^2-(0.033631298873275313246)*tf^2*kp*kd^2*ki+(2.5798925946470399987e-14)*tf^5*kp*kd^2*ki-(3.0091910341667586039e-10)*tf^4*ki+(77.84780273437499999)*tf^2*ki^4+(4.139416050647039998e-14)*tf^4*kp^2*kd^2+(2.6243999999999999988e-19)*tf^7*kp^5-(7.971612375599999996e-14)*tf^6*kp^3*kd*ki-(4.9158285469705520987e-9)*tf^4*kd^3*ki^2-(7.522977585416896509e-12)*tf^7*ki^2+(0.0017299511718751325873)*kp*kd*ki^2+(2.9893556249999999996e-5)*tf*kd^4*ki+(9.2956042968745690447e-7)*tf*kp*kd*ki^2+(3.9835171700435981843e-9)*tf^3*kd*ki+(5.380949814381601922e-5)*tf*ki^2-(0.033215112675703124996)*tf^4*kp*kd*ki^3+(1.899646502215679999e-14)*tf^3*kp*kd^3-(5.2474981040007645666e-13)*tf^6*ki^2+(3.8955937499998705142e-4)*tf^5*kp*ki^4+(1.3286040518588228994e-9)*tf^3*kp*kd^3*ki-(0.016607560250531249998)*tf^3*kp^2*kd*ki^2+(6.040587622685443461e-11)*tf^8*ki^3+(3.690566873999999999e-15)*tf^6*kp^5+(2.9930776233986544344e-8)*tf^4*kp^2*ki+(2.9893556249999999996e-5)*kp*kd^4-(1.4171758687799999995e-13)*tf^7*kp^2*kd*ki^2+(2.4186470399999999987e-19)*tf^5*kd^3*ki+(0.0086414761191152418275)*tf*kd*ki^2-(6.890444751158099999e-6)*tf^2*kd^3*ki-(5.108178833438146558e-15)*tf^5*kp^2*kd+(2.9893514710656842997e-5)*tf^5*kd^2*ki^3+(3.5367890625000175664e-4)*tf^2*kd^2*ki^3-(6.7260489106851562485e-4)*kp^2*ki-(4.8977807661268991972e-11)*tf^6*kd*ki^2-(3.7498525937559024343e-8)*tf^2*kp^2*kd+(6.377303337407999998e-11)*tf^2*kp^3+(2.1930918608583355385e-10)*tf^4*kp^3*kd+(1.0339708485319860019e-10)*tf^3*kp*ki+(8.494587772564136487e-7)*tf^7*ki^3-(2.139864107581439999e-14)*tf^5*ki+(1.0497599999999999995e-18)*tf^4*kp^2*kd^3-(1.1428097589349922633e-11)*tf^4*kd^2*ki-(5.3496602689535999975e-16)*tf^7*kp*ki+(3.4012304621567999987e-15)*tf^6*kp^4-(0.033215080599703125002)*tf^5*kd*ki^4+(2.3042949609374999995e-5)*kp^3*kd+(7.473389062499999998e-9)*kp^3-(2.6156857451537109368e-5)*tf^2*kp^3*ki+(0.0029575246466160710385)*tf^2*kp^2*kd*ki-(1.8338433372776557145e-11)*tf^5*kp^3*ki-(1.9486170000001059961e-9)*tf*kp*ki^2+(2.6243999999999999993e-11)*tf^3*kp^4*kd-(0.44840908331279999985)*kd^2*ki+(8.494640599825929355e-7)*tf^6*kp*ki^2+(4.982259374999999999e-7)*kp^2*kd^3+(5.978711865091368599e-5)*tf^6*kp*kd*ki^3+(8.7404272567131829223e-10)*tf^5*kp^2*kd*ki+(0.008303774259093749997)*tf^5*kp^2*ki^3-(9.795535991390589888e-10)*tf^4*kp*ki-(2.3066026004644101174e-6)*tf^4*kp^2*kd*ki^2+(2.2592160410025676885e-11)*tf^3*kp^3*kd+(1.09240658747999999954e-13)*tf^4*kp*kd^3*ki-(2.3168529009843749994e-11)*tf^5*kp^4*ki-(8.8163437499999552927e-4)*tf^3*kp*kd*ki^3+(8.4945388150931045225e-7)*tf^3*kd^2*ki+(0.009039859112205835467)*tf^2*ki^2+(1.5647985437399999995e-13)*tf^3*kp^2*kd^3-(0.9445523433918749999)*tf^2*kd*ki^3+(1.7577399595874399995e-4)*kp*ki+(6.6430236764812499984e-10)*tf^4*kp^4*kd+(4.251526571159999999e-9)*tf^7*kp*kd*ki^3+(1.59432300000003637155e-12)*tf*kp^3-(1.5595234559999999993e-14)*tf^6*kd^2*ki^2+(3.5263936534531276787e-9)*tf^2*kp*kd+(5.3496602689535999975e-16)*tf^5*kp*kd-(1.256280555314910781e-34)*tf^6*kp^3+(5.9049000000002357413e-7)*tf^8*ki^5+(0.03363108151220090518)*tf^5*kd*ki^3+(1.2732440624999999997e-6)*tf*kp^3*kd^2-(3.1139714120923145256)*tf*kd^2*ki^2-(0.0112104104303523613375)*tf^7*ki^4-(0.016607522946234375)*kp*kd^3*ki-(5.535835999818172366e-7)*tf*kp*kd^3*ki+(0.25223202254913681237)*ki^2+(2.7900652499999999992e-9)*tf^2*kd^4*ki+(2.4186470399999999987e-19)*tf^4*kp*kd^3+(5.1081788334381465576e-15)*tf^7*kp^2*ki+(2.560327734374999999e-10)*tf^2*kp^4+(7.0585795215360231896e-13)*tf^4*kd*ki-(4.7239195625999999982e-14)*tf^8*kp^3*ki^2+(3.5263977391531434278e-10)*tf^3*kp^2*ki+(2.2143401243999999988e-14)*tf^5*kp^2*kd^2*ki-(0.002543733571165023312)*tf^5*kp*ki^3-(0.007963650272475359998)*kd*ki-(0.0025949267578129209255)*tf^2*kp*ki^3-(7.687650033810499649e-6)*tf^6*kp^2*ki^2+(0.011210227082819999997)*kp*kd^3-(7.971837968963533823e-7)*tf^8*ki^4+(2.6243999999999999988e-19)*tf^3*kp*kd^4-(1.9929059274316562998e-9)*tf^6*kp^2*kd*ki^2+(0.0083037754664999999985)*tf^7*ki^5+(4.059630703322175465e-7)*tf^2*kp^2*kd^2*ki-(2.1774314375999999993e-6)*tf^3*kp^3*kd*ki+(2.3066015624999999997)*tf*kd^2*ki^3-(1.0331206917799749596e-9)*tf*kp^2*kd+(3.0242109375000492122e-4)*tf^6*ki^5-(4.529823797634901151e-27)*tf*kp^2+(2.1257632855799999995e-9)*tf^6*kd^2*ki^3-(3.1190469119999999986e-14)*tf^7*kp*kd*ki^2+(0.89682452064787713413)*tf^2*kd*ki^2-(1.6989136050623459597e-6)*tf^5*kd*ki^2+(0.011210115480209999997)*tf*kd^3*ki-(1.0805776531532163071e-10)*tf^3*kp^2*kd^2-(8.968076499165012251e-5)*tf^4*kp*kd^2*ki^2-(105.098270387974690475)*tf^2*ki^3+(0.0072436778961572203366)*tf^3*kp*kd*ki^2+(0.008303786891296874996)*tf^6*kp*ki^4+(1.0497599999999999995e-18)*tf^7*kp^3*kd*ki-(5.9787094293210852836e-5)*tf^3*kd^3*ki^2-(2.0061231952643058266e-11)*tf^3*ki-(4.0828273033358282404e-4)*tf^6*ki^4-(2.1930837147173763068e-10)*tf^6*kp^3*ki-(1.2081178147747334921e-10)*tf^5*kp*kd*ki+(5.1599597584970948617e-8)*tf^5*kp*ki^2-(1.4171758687799999995e-13)*tf^6*kp*kd^2*ki^2-(1.544533220977059356e-4)*tf^2*kp*kd*ki+(1.9909125681188399994e-4)*kp*kd^2+(1.4024137558978175707e-6)*tf^6*kp^2*ki^3-(0.064597960673630362705)*tf^2*kp*ki^2+(2.3066015624999999997)*kp*kd^2*ki^2+(4.9207500000000847784e-4)*tf*kp*kd^2*ki^2+(1.4614613097694478995e-9)*tf^7*kp^3*ki^2-(0.18683472656250899036)*ki^3-(2.0482618157100434786e-6)*tf^6*kd*ki^4-(5.108178833438146558e-15)*tf^4*kp*kd^2+(2.9893557439171583152e-5)*tf^6*kp^3*ki^2+(1.0079160332942940909e-5)*tf^5*kp*kd*ki^2-(4.6132031249999999995)*tf^3*kd*ki^4-(8.246402229540165526e-32)*tf^5*kp^2+(9.336372358045689739e-7)*tf^2*kp*ki-(77.84780273437499999)*kd*ki^3-(1.8452810791406249996e-7)*tf^4*kp^4*ki-(1.6529940863999999996e-5)*ki+(4.1324852159999999985e-7)*kp*kd+(2.230739999999773888e-11)*tf^7*kd*ki^4+(6.2280382714411020727)*tf^2*kp*kd*ki^2-(3.2042265706255059434e-13)*tf^5*kp*ki+(7.2559411199999999964e-19)*tf^6*kp*kd^2*ki+(2.3066015624999999997)*tf^5*ki^5+(4.969957500000029125e-14)*tf^3*kp^4-(1.2194004097843199995e-14)*tf^7*kp^3*ki-(1.5595234559999999993e-14)*tf^8*kp^2*ki^2-(5.580130130943531298e-9)*tf^5*kp^3*kd*ki-(0.00469994441802286499)*tf*kp*kd^2*ki-(3.1140285696480325486)*kp*kd^2*ki+(7.522977585416896509e-12)*tf^4*kp*kd+(1.3604921848627199992e-13)*tf^4*kp^3+(1.8452812499999999996e-7)*tf^2*kp^4*kd+(2.9893603940256842993e-5)*tf^7*kp^2*ki^3+(6.9892914183125004017e-4)*tf^4*kd*ki^3+(6.200146312199999998e-14)*tf^5*kp^4*kd-(0.016607528174531249998)*tf^2*kp^3*kd*ki-(4.6943946716737762124e-9)*tf^3*kp^3*ki-(1.659088704932107864e-7)*tf*kp^2*kd^2-(0.016607522946234374998)*tf*kd^3*ki^2+(7.522977585416896509e-12)*tf^5*kd*ki-(3.1139955161391131985)*tf^5*ki^4-(3.3978097495031845542e-5)*tf^2*kd*ki+(9.04188184290028645e-7)*tf^7*kp*ki^4-(3.273667036577279999e-10)*tf^2*kp*kd^3+(2.6243999999999999988e-19)*tf^4*kd^4*ki+(7.65984382987727923e-10)*tf^4*kp*kd^2*ki-(7.8214271250607683846e-8)*tf^3*kp*kd*ki+(0.041518828364203124995)*tf^3*kd^2*ki^3-(5.3331149817731264207e-12)*tf^3*kp^2*kd-(3.273667141553279999e-10)*tf^3*kd^3*ki+(7.958923724081594567e-7)*tf^3*ki^2-(3.114006587866471274)*tf^4*kp*ki^3-(7.6598400297427199974e-10)*tf^7*kp^2*ki^2-(6.6905075810994386885e-8)*tf^4*kd*ki^2-(1.37440037159004029905e-33)*tf^4*kp^2-(2.304287079279960937e-5)*tf*kp^2*kd*ki+(2.4186470399999999987e-19)*tf^8*kp^3*ki-(0.2802517057945898437)*kp^2*kd*ki+(2.1257632855799999995e-9)*tf^8*kp^2*ki^3+(2.7900652499999999992e-9)*tf*kp*kd^4+(1.9928973065512499992e-5)*tf^3*kp^2*ki^2+(7.183394479243987723e-11)*tf^6*kp^2*ki+(105.09814582785880195)*kd*ki^2-(4.6132031249999999995)*tf^2*kp*kd*ki^3-(1.9929048535744535967e-6)*tf^5*kp*kd*ki^3-(5.3915431735295999985e-15)*tf^6*kp^2*kd*ki+(9.6745881600003042145e-18)*tf^5*kp^3-(1.4348904342795044232e-7)*tf*kp^2*ki-(0.0074737681570904975972)*tf^3*ki^3+(4.118669582068124999e-9)*tf^3*kp^3*kd^2-(8.968065173237571551e-5)*tf^3*kp^2*kd^2*ki+(1.537525726583923339e-5)*tf^4*kp^2*kd*ki-(4.4638915873173038986e-5)*tf^4*kp*ki^2-(0.033630798156569443665)*tf^3*kd^2*ki^2+(0.33214989918937499992)*tf^4*ki^4-(1.4105574613812510714e-7)*tf^2*ki+(4.7289945234375262124e-9)*tf*kp^3*kd+(8.4945464136791133854e-7)*tf^2*kp*kd^2-(9.134142187499673382e-6)*kp*ki^2+(4.7470896537058799988e-4)*tf*kd^2*ki+(5.978715900108749999e-5)*tf*kp^2*kd^3+(7.2559411199999999964e-19)*tf^6*kp^3*kd+(2.230739999999773888e-11)*tf^8*kp*ki^4)");
        Function rS7("kp","ki","kd","tf","(0.0036)*ki");
        //    Function rGS1("kp","ki","kd","0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki");
    //    Function rGS2("kp","ki","kd","((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)");
    //    Function rGS3("kp","ki","kd","(1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)");
    //    Function rGS4("kp","ki","kd","-1/(1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*(1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))*((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+(0.0036)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*ki*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)");
    //    Function rGS5("kp","ki","kd","(0.0036)*ki");
    //    Function rT1("kp","ki","kd","0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki");
    //    Function rT2("kp","ki","kd","((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)");
    //    Function rT3("kp","ki","kd","1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))");
    //    Function rT4("kp","ki","kd","-((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))+(0.0036)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*ki)*1/(((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))");
    //    Function rT5("kp","ki","kd", "(0.0036)*ki");
    cout<<"routh ok"<<endl;
    vector<Function*> routh_func;
    vector<Ctc*> array_ctc;
    vector<Sep*> sep_array;
    NumConstraint *c1= new NumConstraint(kp,ki,kd,tf,rS1(kp,ki,kd,tf)>=0);
    array_ctc.push_back(new CtcFwdBwd(NumConstraint(rS1,LEQ)));
    routh_func.push_back(&rS1);
    sep_array.push_back(new SepFwdBwd(rS1,GEQ));
    NumConstraint *c2= new NumConstraint(kp,ki,kd,tf,rS2(kp,ki,kd,tf)>=0);
    array_ctc.push_back(new CtcFwdBwd(NumConstraint(rS2,LEQ)));
    routh_func.push_back(&rS2);
    sep_array.push_back(new SepFwdBwd(rS2,GEQ));
    NumConstraint *c3= new NumConstraint(kp,ki,kd,tf,rS3(kp,ki,kd,tf)>=0);
    array_ctc.push_back(new CtcFwdBwd(NumConstraint(rS3,LEQ)));
    routh_func.push_back(&rS3);
    sep_array.push_back(new SepFwdBwd(rS3,GEQ));
    NumConstraint *c4= new NumConstraint(kp,ki,kd,tf,rS4(kp,ki,kd,tf)>=0);
    array_ctc.push_back(new CtcFwdBwd(NumConstraint(rS4,LEQ)));
    routh_func.push_back(&rS4);
    sep_array.push_back(new SepFwdBwd(rS4,GEQ));
    NumConstraint *c5= new NumConstraint(kp,ki,kd,tf,rS5(kp,ki,kd,tf)>=0);
    array_ctc.push_back(new CtcFwdBwd(NumConstraint(rS5,LEQ)));
    sep_array.push_back(new SepFwdBwd(rS5,GEQ));
    routh_func.push_back(&rS5);
    NumConstraint *c6= new NumConstraint(kp,ki,kd,tf,rS6(kp,ki,kd,tf)>=0);
    array_ctc.push_back(new CtcFwdBwd(NumConstraint(rS6,LEQ)));
    sep_array.push_back(new SepFwdBwd(rS6,GEQ));
    routh_func.push_back(&rS6);
    NumConstraint *c7= new NumConstraint(kp,ki,kd,tf,rS7(kp,ki,kd,tf)>=0);
    array_ctc.push_back(new CtcFwdBwd(NumConstraint(rS7,LEQ)));
    sep_array.push_back(new SepFwdBwd(rS7,GEQ));
    routh_func.push_back(&rS7);
//    NumConstraint *c8= new NumConstraint(kp,ki,kd,rT3(kp,ki,kd)>=0);
//    array_ctc.push_back(new CtcFwdBwd(*c8));
//    NumConstraint *c9= new NumConstraint(kp,ki,kd,rT4(kp,ki,kd)>=0);
//    array_ctc.push_back(new CtcFwdBwd(*c9));
//    NumConstraint *c10= new NumConstraint(kp,ki,kd,rT5(kp,ki,kd)>=0);
//    array_ctc.push_back(new CtcFwdBwd(*c10));
//    NumConstraint *c11= new NumConstraint(kp,ki,kd,rGS1(kp,ki,kd)>=0);
//    array_ctc.push_back(new CtcFwdBwd(*c11));
//    NumConstraint *c12= new NumConstraint(kp,ki,kd,rGS2(kp,ki,kd)>=0);
//    array_ctc.push_back(new CtcFwdBwd(*c12));
//    NumConstraint *c13= new NumConstraint(kp,ki,kd,rGS3(kp,ki,kd)>=0);
//    array_ctc.push_back(new CtcFwdBwd(*c13));
//    NumConstraint *c14= new NumConstraint(kp,ki,kd,rGS4(kp,ki,kd)>=0);
//    array_ctc.push_back(new CtcFwdBwd(*c14));
//    NumConstraint *c15= new NumConstraint(kp,ki,kd,rGS5(kp,ki,kd)>=0);
//    array_ctc.push_back(new CtcFwdBwd(*c15));
    CtcUnion ctc_routh(array_ctc);

    SystemFactory fac;
    fac.add_var(kp);
    fac.add_var(ki);
    fac.add_var(kd);
    fac.add_var(tf);
    fac.add_ctr(*c1);
    fac.add_ctr(*c2);
    fac.add_ctr(*c3);
    fac.add_ctr(*c4);
    fac.add_ctr(*c5);
    fac.add_ctr(*c6);
    fac.add_ctr(*c7);

    NormalizedSystem sys(fac) ;

    Ctc* hc4=new CtcHC4(sys,PROPAG_RATIO,HC4_INCREMENTAL);
    Ctc* hc44cid=new CtcHC4(sys.ctrs,0.1,true);
    Ctc* acidhc4=new CtcAcid(sys,*hc44cid,EQ_EPS);
    Ctc* ctc=new CtcCompo(*hc4, *acidhc4);

    vector<LinearRelaxXTaylor::corner_point> cpoints;
    cpoints.push_back(LinearRelaxXTaylor::RANDOM);
    cpoints.push_back(LinearRelaxXTaylor::RANDOM_INV);
    LinearRelax* lr = new LinearRelaxXTaylor(sys,cpoints);
    Ctc* cxn_poly = new CtcPolytopeHull(*lr, CtcPolytopeHull::ALL_BOX);
    Ctc* hc44xn = new CtcHC4(sys.ctrs,PROPAG_RATIO,false);
    Ctc* cxn_compo = new CtcCompo(*cxn_poly, *hc44xn);
    Ctc* cxn = new CtcFixPoint(*cxn_compo, FIXPOINT_RATIO);

    Ctc* ctc_final = new CtcCompo(*ctc, *cxn);


//    var[0] = "kp";var[1] = "ki";var[2] = "kd";var[3] = "T";
//    string cpol = "(0.0001662448486980308*ki + 0.011839198205153703*ki*s + 0.0001662448486980308*kp*s + 0.006649793947921231*T*s^4 + 0.4735679282061481*T*s^5 + 18.421216097060995*T*s^6 + 111.6266663114158*T*s^7 + 3950.3495770125814*T*s^8 + 88.8858291332026*T*s^9 + 1*T*s^10 + 0.0001662448486980308*kd*s^2 + 0.011839198205153703*kd*s^3 + 0.4677458906512658*kd*s^4 + 3.3044004884136675*kd*s^5 + 117.04734585855945*kd*s^6 + 2.6336538673212293*kd*s^7 + 0.02962962962962963*kd*s^8 + 0.4677458906512658*ki*s^2 + 3.3044004884136675*ki*s^3 + 117.04734585855945*ki*s^4 + 2.6336538673212293*ki*s^5 + 0.02962962962962963*ki*s^6 + 0.011839198205153703*kp*s^2 + 0.4677458906512658*kp*s^3 + 3.3044004884136675*kp*s^4 + 117.04734585855945*kp*s^5 + 2.6336538673212293*kp*s^6 + 0.02962962962962963*kp*s^7 + 0.006649793947921231*s^3 + 0.4735679282061481*s^4 + 18.421216097060995*s^5 + 111.6266663114158*s^6 + 3950.3495770125814*s^7 + 88.8858291332026*s^8 + 1*s^9 + 0.0001662448486980308*T*ki*s + 0.011839198205153703*T*ki*s^2 + 0.4677458906512658*T*ki*s^3 + 3.3044004884136675*T*ki*s^4 + 117.04734585855945*T*ki*s^5 + 2.6336538673212293*T*ki*s^6 + 0.02962962962962963*T*ki*s^7 + 0.0001662448486980308*T*kp*s^2 + 0.011839198205153703*T*kp*s^3 + 0.4677458906512658*T*kp*s^4 + 3.3044004884136675*T*kp*s^5 + 117.04734585855945*T*kp*s^6 + 2.6336538673212293*T*kp*s^7 + 0.02962962962962963*T*kp*s^8)/(1*T)";
//    pair<polynomial*,polynomial*> cpolfrac = real_part_pol(cpol,var,4);
    vector<Function> pol_coef;
//    cout<<"create pol coef "<<endl;
//    pol_coef.push_back(Function("kp","ki","kd","T","(1.662448486980308e-4)*ki/T"));
//    cout<<"test1"<<endl;
//    pol_coef.push_back(Function("kp","ki","kd","T", "1/T*((0.011839198205153703)*ki+(1.662448486980308e-4)*ki*T+(1.662448486980308e-4)*kp)"));
//    cout<<"test2"<<endl;
//    pol_coef.push_back(Function("kp","ki","kd","T","((0.4677458906512658)*ki+(1.662448486980308e-4)*kp*T+(0.011839198205153703)*ki*T+(1.662448486980308e-4)*kd+(0.011839198205153703)*kp)/T"));
//    cout<<"test3"<<endl;
//    pol_coef.push_back(Function("kp","ki","kd","T","(0.006649793947921231+(3.3044004884136675)*ki+(0.011839198205153703)*kp*T+(0.4677458906512658)*ki*T+(0.011839198205153703)*kd+(0.4677458906512658)*kp)/T"));
//    cout<<"test4"<<endl;
//    pol_coef.push_back(Function("kp","ki","kd","T"," (0.4735679282061481+(117.04734585855945)*ki+(0.4677458906512658)*kp*T+(3.3044004884136675)*ki*T+(0.4677458906512658)*kd+(3.3044004884136675)*kp+(0.006649793947921231)*T)/T"));
//    cout<<"test5"<<endl;
//    pol_coef.push_back(Function("kp","ki","kd","T","(18.421216097060995+(2.6336538673212293)*ki+(3.3044004884136675)*kp*T+(117.04734585855945)*ki*T+(3.3044004884136675)*kd+(117.04734585855945)*kp+(0.4735679282061481)*T)/T"));
//    cout<<"test6"<<endl;
//    pol_coef.push_back(Function("kp","ki","kd","T","1/T*(111.6266663114158+(0.02962962962962963)*ki+(117.04734585855945)*kp*T+(2.6336538673212293)*ki*T+(117.04734585855945)*kd+(2.6336538673212293)*kp+(18.421216097060995)*T)"));
//    cout<<"test7"<<endl;
//    pol_coef.push_back(Function("kp","ki","kd","T","1/T*(3950.3495770125814+(2.6336538673212293)*kp*T+(0.02962962962962963)*ki*T+(2.6336538673212293)*kd+(0.02962962962962963)*kp+(111.6266663114158)*T)"));
//    cout<<"test8"<<endl;
//    pol_coef.push_back(Function("kp","ki","kd","T","(88.8858291332026+(0.02962962962962963)*kp*T+(0.02962962962962963)*kd+(3950.3495770125814)*T)/T"));

//    cout<<"test9"<<endl;
//    pol_coef.push_back(Function("kp","ki","kd","T","(1+(88.8858291332026)*T)/T"));
//    cout<<"pol coef function ok"<<endl;
//    Function real_part("kp","ki","kd","T","w", "(-(0.02962962962962963)*ki*w^6+(1.662448486980308e-4)*ki+(117.04734585855945)*ki*w^4+(3950.3495770125814)*T*w^8-(0.011839198205153703)*kp*w^2-(2.6336538673212293)*ki*T*w^6+(0.4735679282061481)*w^4-(117.04734585855945)*kd*w^6-(117.04734585855945)*T*kp*w^6+(3.3044004884136675)*ki*T*w^4-(111.6266663114158)*w^6+(0.4677458906512658)*kd*w^4+(0.4677458906512658)*T*kp*w^4+(0.006649793947921231)*T*w^4-(1.662448486980308e-4)*kd*w^2-(1.662448486980308e-4)*T*kp*w^2-(0.011839198205153703)*ki*T*w^2-(18.421216097060995)*T*w^6+(0.02962962962962963)*kd*w^8-(2.6336538673212293)*kp*w^6+(0.02962962962962963)*T*kp*w^8-T*w^10+(3.3044004884136675)*kp*w^4+(88.8858291332026)*w^8-(0.4677458906512658)*ki*w^2)/T");

//    Function imag_part("kp","ki","kd","T","w", "(0.4735679282061481)*T*w^5-(2.6336538673212293)*w^7*kd-(0.4677458906512658)*T*w^3*ki-(111.6266663114158)*T*w^7+(2.6336538673212293)*w^5*ki-(0.011839198205153703)*T*kp*w^3+(117.04734585855945)*kp*w^5+(0.011839198205153703)*w*ki+(3.3044004884136675)*w^5*kd+w^9-(0.02962962962962963)*kp*w^7-(0.006649793947921231)*w^3-(2.6336538673212293)*T*kp*w^7-(0.02962962962962963)*T*w^7*ki+(18.421216097060995)*w^5-(0.011839198205153703)*w^3*kd-(0.4677458906512658)*kp*w^3-(3950.3495770125814)*w^7+(3.3044004884136675)*T*kp*w^5+(1.662448486980308e-4)*T*w*ki+(88.8858291332026)*T*w^9+(117.04734585855945)*T*w^5*ki-(3.3044004884136675)*w^3*ki+(1.662448486980308e-4)*kp*w");
//    NumConstraint wzero(real_part,EQ);
//    CtcFwdBwd fw(wzero);
//    BitSet mask(BitSet::all(5));
//    mask.remove(4);

////    9.39156 ; 0.0358978 ; 3.1466 ; 2.01083
////    [-8.72881, -8.1775] ; [6.975, 7.25065] ; [9.08469, 9.49658] ; [1.1195, 1.16699]
////    [-10, -9.92527] ; [9.42821, 9.49672] ; [9.54245, 9.59852] ; [1.00253, 1.00912]
/// -9.8,-9.58,8.3,0.9
/// -9.95,9.45,9.57,1.006
    Interval Kp(-10, -9.92527),Ki(9.42821, 9.49672),Kd(9.54245, 9.59852),T(1.00253, 1.00912);
    IntervalVector stabox(4),staboxw(5);
    stabox[0] =Kp;stabox[1] =Ki;stabox[2] =Kd;stabox[3] =T;
    vector<Function *> pt_coefs;
//    staboxw[0] =Kp;staboxw[1] =Ki;staboxw[2] =Kd;staboxw[3] =T;

//    Interval co = cutoff_bound(pol_coef,stabox);
//    cout<<"cutoff freq: "<<co<<endl;
//    IntervalVector w_dom(1,Interval(0,1));
//    CtcExist fex(fw,mask,w_dom,0.01);
//    vibes::beginDrawing();
//    vibes::newFigure("real part pol");
//    double wmin(0),wmax(1),freqeps(0.001);
//    Interval freq(wmin,wmin+freqeps);
//    IntervalVector cdraw(2);
//    while(freq.ub()<wmax) {
//        cdraw[0] = freq;
//        staboxw[4] = freq;
//        cdraw[1] = real_part.eval_vector(staboxw)[0];
//        vibes::drawBox(cdraw,"red[red]");
//        cdraw[1] = imag_part.eval_vector(staboxw)[0];
//        vibes::drawBox(cdraw,"blue[blue]");
//        freq+=Interval(freqeps);
//    }
//    vibes::endDrawing();
//    return 0;

    cout<<endl<<"*******************"<<endl;
    pt_coefs.push_back(new Function("kp","ki","kd","T","(10537019111845495613241032739783*ki)/(63382530011411470074835160268800000*T)"));
    pt_coefs.push_back(new Function("kp","ki","kd","T","(6123437479178601467006569798917525860280205615571544779328*ki + 85984702646059087287357082558915403661686951915714248704*kp + 85984702646059087287357082558915403661686951915714248704*T*ki)/(517217245042237489695942121663229349471324453736965734400000*T)"));
    pt_coefs.push_back(new Function("kp","ki","kd","T","(85984702646059087287357082558915403661686951915714248704*kd + 241926240942475369556027638570784249909947691164891995634641*ki + 6123437479178601467006569798917525860280205615571544779328*kp + 6123437479178601467006569798917525860280205615571544779328*T*ki + 85984702646059087287357082558915403661686951915714248704*T*kp)/(517217245042237489695942121663229349471324453736965734400000*T)"));
    pt_coefs.push_back(new Function("kp","ki","kd","T","(6123437479178601467006569798917525860280205615571544779328*kd + 1709092917133541268429605152874087053884471459553034400822336*ki + 241926240942475369556027638570784249909947691164891995634641*kp + 241926240942475369556027638570784249909947691164891995634641*T*ki + 6123437479178601467006569798917525860280205615571544779328*T*kp + 3439388105842363491494283302356616146467478076628569948160)/(517217245042237489695942121663229349471324453736965734400000*T)"));
    pt_coefs.push_back(new Function("kp","ki","kd","T","(3439388105842363491494283302356616146467478076628569948160*T + 241926240942475369556027638570784249909947691164891995634641*kd + 60538905764470062617174426996447380495905857364817563757838336*ki + 1709092917133541268429605152874087053884471459553034400822336*kp + 1709092917133541268429605152874087053884471459553034400822336*T*ki + 241926240942475369556027638570784249909947691164891995634641*T*kp + 244937499167144058674511452919335327040734690096633144672256)/(517217245042237489695942121663229349471324453736965734400000*T)"));
    pt_coefs.push_back(new Function("kp","ki","kd","T","(244937499167144058674511452919335327040734690096633144672256*T + 1709092917133541268429605152874087053884471459553034400822336*kd + 1362171197650720615916740156171465615738890689107575439360000*ki + 60538905764470062617174426996447380495905857364817563757838336*kp + 60538905764470062617174426996447380495905857364817563757838336*T*ki + 1709092917133541268429605152874087053884471459553034400822336*T*kp + 9527770640049606644179955720428217957821860374411197852680192)/(517217245042237489695942121663229349471324453736965734400000*T)"));
    pt_coefs.push_back(new Function("kp","ki","kd","T","(9527770640049606644179955720428217957821860374411197852680192*T + 60538905764470062617174426996447380495905857364817563757838336*kd + 15324955408658888583583470271503091836187391221836021760000*ki + 1362171197650720615916740156171465615738890689107575439360000*kp + 1362171197650720615916740156171465615738890689107575439360000*T*ki + 60538905764470062617174426996447380495905857364817563757838336*T*kp + 57735236822839624422708535872514782220796956928384476102262784)/(517217245042237489695942121663229349471324453736965734400000*T)"));
    pt_coefs.push_back(new Function("kp","ki","kd","T","(57735236822839624422708535872514782220796956928384476102262784*T + 1362171197650720615916740156171465615738890689107575439360000*kd + 15324955408658888583583470271503091836187391221836021760000*kp + 15324955408658888583583470271503091836187391221836021760000*T*ki + 1362171197650720615916740156171465615738890689107575439360000*T*kp + 2043188925176215615677889430516454170045046702214280110789361664)/(517217245042237489695942121663229349471324453736965734400000*T)"));
    pt_coefs.push_back(new Function("kp","ki","kd","T","(2043188925176215615677889430516454170045046702214280110789361664*T + 15324955408658888583583470271503091836187391221836021760000*kd + 15324955408658888583583470271503091836187391221836021760000*T*kp + 45973283667570099034270085511265165406205649806935123740852224)/(517217245042237489695942121663229349471324453736965734400000*T)"));
    pt_coefs.push_back(new Function("kp","ki","kd","T","(45973283667570099034270085511265165406205649806935123740852224*T + 517217245042237489695942121663229349471324453736965734400000)/(517217245042237489695942121663229349471324453736965734400000*T)"));
    pt_coefs.push_back(new Function("kp","ki","kd","T","1"));

    univar_polynomial upol(pt_coefs);
    pair<icomplex,Interval>* discs;1;
    int nbroot = upol.get_roots(stabox,discs);
//    return 0;
    pol_coef.clear();
    pol_coef.push_back(Function("kp","ki","kd","T","1"));
    pol_coef.push_back(Function("kp","ki","kd","T","(45973283667570099034270085511265165406205649806935123740852224*T + 517217245042237489695942121663229349471324453736965734400000)/(517217245042237489695942121663229349471324453736965734400000*T)"));
    pol_coef.push_back(Function("kp","ki","kd","T","(2043188925176215615677889430516454170045046702214280110789361664*T + 15324955408658888583583470271503091836187391221836021760000*kd + 15324955408658888583583470271503091836187391221836021760000*T*kp + 45973283667570099034270085511265165406205649806935123740852224)/(517217245042237489695942121663229349471324453736965734400000*T)"));
    pol_coef.push_back(Function("kp","ki","kd","T","(57735236822839624422708535872514782220796956928384476102262784*T + 1362171197650720615916740156171465615738890689107575439360000*kd + 15324955408658888583583470271503091836187391221836021760000*kp + 15324955408658888583583470271503091836187391221836021760000*T*ki + 1362171197650720615916740156171465615738890689107575439360000*T*kp + 2043188925176215615677889430516454170045046702214280110789361664)/(517217245042237489695942121663229349471324453736965734400000*T)"));
    pol_coef.push_back(Function("kp","ki","kd","T","(9527770640049606644179955720428217957821860374411197852680192*T + 60538905764470062617174426996447380495905857364817563757838336*kd + 15324955408658888583583470271503091836187391221836021760000*ki + 1362171197650720615916740156171465615738890689107575439360000*kp + 1362171197650720615916740156171465615738890689107575439360000*T*ki + 60538905764470062617174426996447380495905857364817563757838336*T*kp + 57735236822839624422708535872514782220796956928384476102262784)/(517217245042237489695942121663229349471324453736965734400000*T)"));
    pol_coef.push_back(Function("kp","ki","kd","T","(244937499167144058674511452919335327040734690096633144672256*T + 1709092917133541268429605152874087053884471459553034400822336*kd + 1362171197650720615916740156171465615738890689107575439360000*ki + 60538905764470062617174426996447380495905857364817563757838336*kp + 60538905764470062617174426996447380495905857364817563757838336*T*ki + 1709092917133541268429605152874087053884471459553034400822336*T*kp + 9527770640049606644179955720428217957821860374411197852680192)/(517217245042237489695942121663229349471324453736965734400000*T)"));
    pol_coef.push_back(Function("kp","ki","kd","T","(3439388105842363491494283302356616146467478076628569948160*T + 241926240942475369556027638570784249909947691164891995634641*kd + 60538905764470062617174426996447380495905857364817563757838336*ki + 1709092917133541268429605152874087053884471459553034400822336*kp + 1709092917133541268429605152874087053884471459553034400822336*T*ki + 241926240942475369556027638570784249909947691164891995634641*T*kp + 244937499167144058674511452919335327040734690096633144672256)/(517217245042237489695942121663229349471324453736965734400000*T)"));
    pol_coef.push_back(Function("kp","ki","kd","T","(6123437479178601467006569798917525860280205615571544779328*kd + 1709092917133541268429605152874087053884471459553034400822336*ki + 241926240942475369556027638570784249909947691164891995634641*kp + 241926240942475369556027638570784249909947691164891995634641*T*ki + 6123437479178601467006569798917525860280205615571544779328*T*kp + 3439388105842363491494283302356616146467478076628569948160)/(517217245042237489695942121663229349471324453736965734400000*T)"));
    pol_coef.push_back(Function("kp","ki","kd","T","(85984702646059087287357082558915403661686951915714248704*kd + 241926240942475369556027638570784249909947691164891995634641*ki + 6123437479178601467006569798917525860280205615571544779328*kp + 6123437479178601467006569798917525860280205615571544779328*T*ki + 85984702646059087287357082558915403661686951915714248704*T*kp)/(517217245042237489695942121663229349471324453736965734400000*T)"));
    pol_coef.push_back(Function("kp","ki","kd","T","(6123437479178601467006569798917525860280205615571544779328*ki + 85984702646059087287357082558915403661686951915714248704*kp + 85984702646059087287357082558915403661686951915714248704*T*ki)/(517217245042237489695942121663229349471324453736965734400000*T)"));
    pol_coef.push_back(Function("kp","ki","kd","T","(10537019111845495613241032739783*ki)/(63382530011411470074835160268800000*T)"));

    IntervalMatrix coef_mat(1,11),cvec(11,1);
//    icomplex c(Interval(- 0.84540181371938381285507169289191),Interval(0.34070063664340792059925292413312));
    icomplex c(Interval(0.34837408242172199618777614744032),Interval(0.46222767016820064334679430164482));
//    Interval pol_ceval(0);
    icomplex cpow(Interval(1),Interval(0));
    icomplex pol_eval(Interval(0),Interval(0));
    for(unsigned i=0;i<pol_coef.size();i++) {
        coef_mat[0][i] = pol_coef.at(i).eval_vector(stabox)[0];
//        cout<<"polev "<<i<<" : "<<coef_mat[0][i]<<endl;
//        cvec[i][0] = ibex::pow(Interval(c),Interval(pol_coef.size()-1-i));
//        cout<<"cev "<<i<<" : "<<cvec[i][0]<<endl;
//        cout<<"prod: "<<cvec[i][0]*coef_mat[0][i]<<endl;
//        pol_ceval+=cvec[i][0]*coef_mat[0][i];
        Interval pcoef = pol_coef.at(pol_coef.size()-1-i).eval_vector(stabox)[0];
//        cout<<"pcoef "<<pol_coef.size()-1-i<<": "<<pcoef<<endl;
        if(i!=0)
            cpow = cpow*c;
//        cout<<"cpow: "<<cpow<<endl;
//        cout<<"cpow*pcoef: "<<cpow*pcoef<<endl;
        pol_eval = pol_eval+cpow*pcoef;
//        cout<<"pol_eval: "<<pol_eval<<endl<<endl;

    }
//real [21.3329, 21.7977]

//    IntervalMatrix pol_ceval = coef_mat*cvec;

    cout<<"pol_eval: "<<pol_eval<<"*i" <<endl;
//    Interval R = ibex::pow(ibex::pow(Interval(2),19)*ibex::abs(pol_ceval),0.1);
    icomplex alpha_den( Interval(685.41),Interval(- 256.32));
    icomplex alphanu = pol_eval/alpha_den;
    icomplex rnu = alphanu/coef_mat[0][0];
    rnu = icomplex(rnu.real_part()*5,rnu.imag_part()*5);
    Interval R = 5*rnu.abs();
    cout<<"center: "<<c-rnu<<" ,radius: "<<R<<endl;


//    cout<<"function ok"<<endl;
    return 0;

//    staboxw[4] = w_dom[0];
//    cout<<"w_dom: "<<w_dom<<endl;
//    cout<<"evaluation: "<<real_part.eval_vector(staboxw)<<endl;
//    cout<<"init box: "<<stabox<<endl;
//    fex.contract(stabox);
//    cout<<"after contraction: "<<stabox<<endl;
//    return 0;




//    IntervalVector box1(4);box1[0]=Interval(-10,10);box1[1]=Interval(0,10);box1[2]=Interval(0,10);box1[3]=Interval(-3,3);
//    NumConstraint cst1(kp,ki,kd,Max12(kp,ki,kd,box1[0])<=100);
//    CtcFwdBwd ctc1(cst1);
//    ctc1.contract(box1);
//    cout<<"contraction: "<<box1<<endl;
//    return 0;


//    costf2 b;
//    Heap<heap_elem> heap(b);
//    double freqeps=0.00001;
//    double umin(-3),umax(3);
//    Interval freq(umin,umin+freqeps);
//    IntervalVector box(5),bbox(4),boxmat(5);
//    Interval bernres(-10000),cres;

//    box[0] = Interval(6.86104);box[1]=Interval(0.0122183);box[2]=Interval(8.37996);box[3]=Interval(3.93941);
//    boxmat[0] = Interval(7.953783690601573);boxmat[1]=Interval(0.078943640474045);boxmat[2]=Interval(2.863881179107744);boxmat[3]=Interval(2.263575234363492);
//    bbox[1] = Interval(-9.69235, -9.63094);bbox[2]=Interval(9.54305, 9.59969);bbox[3]=Interval(9.59853, 9.6546);
//    IntervalVector berndraw(2),cdraw(2);
//    vibes::beginDrawing();
//    vibes::newFigure("Twz2 go");
//    while(freq.ub()<umax) {
//        cdraw[0] = freq;
//        berndraw[0] = freq;
//        box[4]=freq;
//        boxmat[4]=freq;
//        bbox[0]=ibex::pow(Interval(10),freq.mid());
////        berndraw[1] = bernwfree.at(1).eval_bernstein(bbox)[0]*berncst.at(1).eval_bernstein(IntervalVector(1,bbox[0]))[0];

//        cdraw[1] = Max12.eval_vector(box)[0];
//        vibes::drawBox(cdraw,"blue[blue]");
//        cdraw[1] = Max12.eval_vector(boxmat)[0];
//        vibes::drawBox(cdraw,"red[red]");
////        vibes::drawBox(berndraw,"red[red]");
////        heap.push(new heap_elem(NULL,bbox,bernres));
//        freq+=Interval(freqeps);
//    }
////    cout<<"max of bernres: "<<heap.top()->eval<<endl;
////    heap.flush();
//    return 0;

//    -1.78007 ; -9.5839 ; -0.265944 ; 1.46923

//    string avar[4];avar[0] = "kp";avar[1] = "ki";avar[2] = "kd";avar[3] = "tf";
//    pair<polynomial*,polynomial*> rs1 = get_pol(&rS4,avar,4,avar,0);
//    [-10, -9.92527] ; [9.42821, 9.49672] ; [9.54245, 9.59852] ; [1.00253, 1.00912]
//    Interval Kp(-10, -9.92527),Ki(9.42821, 9.49672),Kd(9.54245, 9.59852),T(1.00253, 1.00912);
//    IntervalVector sp(4);sp[0] = Kp ;sp[1] = Ki;sp[2] = Kd;sp[3] =T;
//    cout<<"stability of point "<<sp<<" : "<<endl
//          << "  rS1: "<<rS1.eval_vector(sp)
//           << "  rS2: "<<rS2.eval_vector(sp)
//               << "  rS3: "<<rS3.eval_vector(sp)
//                   << "  rS4: "<<rS4.eval_vector(sp)
//                       << "  rS5: "<<rS5.eval_vector(sp)
//                            << "  rS6: "<<rS6.eval_vector(sp)
//                                << "  rS7: "<<rS7.eval_vector(sp)<<endl;
//    ctc_final->contract(sp);
//    cout<<"contraction compo: "<<sp<<endl;
////    cout<<"eval bern: "<<rs1.first->eval_bernstein(sp)[0]/rs1.second->eval_bernstein(sp)[0]<<endl;
//    return 0;


//    Function Max123(kp,ki,kd,u,ibex::max(F3u(kp,ki,kd,u),Max12(kp,ki,kd,u)));

    //other variables for algo

//    Function Prob("kp","ki","w","((kp+w-2)^6+0.2)*ln(1+(kp+w)^2)+((ki+w-2)^6+0.2)*ln(1+(ki+w)^2)"  );

    SetIntervalReg stab(IniboxK,0.1);

    vibes::beginDrawing();
    vibes::newFigure("stability");
    stab.contract(*ctc_final,ctc_routh);
    cout<<"contraction done"<<endl;
    visit_leaves(stab);
    vibes::endDrawing();
    return 0;


    LargestFirst lf;
    double lower_ub(POS_INFINITY);
//    double lower_ub(100);
    double upper_lb(0);

//    costf1 maxlb;
    costf2 maxub;
    Heap<heap_elem> heap(maxub);
    heap_elem * he = new heap_elem(Inilw,Interval::ALL_REALS);
    heap.push(he);


    cost_uplo cu;
    Heap<list_elem> list(cu);

    list_elem *elem = new list_elem(IniboxK, heap,IntervalVector(1,Interval(-1,1000000)));
//    list_elem *elem = new list_elem(Iniprob, tree,IntervalVector(1,Interval::ALL_REALS));
    list.push(elem);


    Vector respoint(4);
//    Vector respoint(2);

    optiw str(lower_ub,preclwmin,&Max12,bernwset,bernwfree,berncst,false);


    double vol_rejected(0),volout(0);
    list_elem * elemtmp;
//    Vector sol(4);
//    //matlab solution
//    sol[0] = 0.3998;sol[1] = 0.4806;sol[2] = -0.1952;
//    cout<<"evaluation of matlab result: "<<Sum.eval_vector(sol)<<endl;
//    return 0;

    bool lbfix(false);
    Timer::start();
    while(!list.empty()) {
        elemtmp = list.pop();
        if(!lbfix)
            upper_lb = elemtmp->fmax[0].lb();
        if(lower_ub-elemtmp->fmax[0].lb()<stop_criterion) { // stop creterion reached
            break;
        }
        cout<<"maxdiam: "<<elemtmp->box.max_diam()<<endl;
//        cout<<" fmax of current element: "<<elemtmp->fmax<<endl;
        if(elemtmp->fmax[0].lb()>lower_ub) { //better lower_ub may be computed since computation of the father, must check
            vol_rejected+= elemtmp->box.volume();
            cout<<"solution deletion,  "<<endl;
            cout<<"loup : "<<lower_ub<<" get for point: "<<respoint<<" uplo: "<<upper_lb<< " volume rejected: "<<vol_rejected/IniboxK.volume()*100<<endl;
            delete elemtmp;
            continue;
        }
//        double ratio = (elemtmp->box[0]).diam()/(IniboxK[0]).diam()+(elemtmp->box[1]).diam()/(IniboxK[1]).diam()+(elemtmp->box[2]).diam()/(IniboxK[2]).diam();
        double ratio = (elemtmp->box[0]).diam()/(IniboxK[0]).diam()+(elemtmp->box[1]).diam()/(IniboxK[1]).diam()+(elemtmp->box[2]).diam()/(IniboxK[2]).diam()
                +(elemtmp->box[2]).diam()/(IniboxK[2]).diam();
        str.preclw = ratio/Inilw.volume()>preclwmin?ratio/(10*Inilw.volume()):preclwmin;
//        cout<<"precision on w: "<<str.preclw<<endl;
        str.lower_ub = lower_ub;
        //******************** Routh contraction ***********************
        IntervalVector inib = elemtmp->box;
        ctc_final->contract(elemtmp->box);
//        ctc_routh.contract(elemtmp->box);
        if(elemtmp->box != inib) {
            vol_rejected += inib.volume()-elemtmp->box.volume();
            if(elemtmp->box.is_empty())
                cout<<"routn deletion "<<endl;
            else
                cout<<"routh contract  "<<endl;
            cout<<"loup : "<<lower_ub<<" get for point: "<<respoint<<" uplo: "<<upper_lb<< " volume rejected: "<<vol_rejected/IniboxK.volume()*100<<endl;
        }
        if(elemtmp->box.is_empty()) {
            delete elemtmp;
            continue;
        }


        // **************** compute (wmax,fmax) *******************
//        Function g(w,Max12(elemtmp->box[0],elemtmp->box[1],elemtmp->box[2],w));
//        Function g(w,Prob(elemtmp->box[0],elemtmp->box[1],w));
//        str.f = &g;

        Vector midp(4);
        bool stab(true);
        for(unsigned j=0;j<1000;j++)
        {
            midp = randpt(elemtmp->box,elemtmp->box.size());
//            cout<<"box: "<<elemtmp->box<<" midpt: "<<midp<<endl;
            stab = true;
            for(unsigned i=0;i<routh_func.size();i++) {
                stab &= routh_func.at(i)->eval_vector(midp)[0].lb()>0;//+(i==3)*(-0.1);
                //            cout<<"routh "<<i<<": "<<routh_func.at(i)->eval_vector(midp)[0]<<endl;
                if(!stab){
//                    cout<<"midp fail: routh("<<i<<"): "<<routh_func.at(i)->eval_vector(midp)[0]<<endl;
                    break;
                }
            }
            if(stab)
            {
                break;}

        }
//        elemtmp->computable = stab;
        if(elemtmp->computable)
            best_upper_bound_forall_w(&str,elemtmp,false,Vector(4));
        cout<<"res found for box "<<elemtmp->box<<" : "<<endl<<"        "<<elemtmp->fmax<<endl;
        if(elemtmp->fmax[0] == Interval::EMPTY_SET) { // maximum is guaranteed to be higher than the current lower_up for box elemtmp->box
            vol_rejected+= elemtmp->box.volume();
            cout<<"solution deletion,  "<<endl;
            cout<<"loup : "<<lower_ub<<" get for point: "<<respoint<<" uplo: "<<upper_lb<< " volume rejected: "<<vol_rejected/IniboxK.volume()*100<<endl;
            delete elemtmp;
            continue;
        }

        //*********** Contract K with constraints f(k,wmax)< lower_ub  ***********
//        if(lower_ub<POS_INFINITY){ // run contraction with every w that contains maximum
//            IntervalVector box(elemtmp->box);
////            cout<<"Init box befor contraction: "<<elemtmp->box<<endl;
////            cout<<"try to get boxes from main"<<endl;
//            vector<heap_elem*> stackheapelem;
//            heap_elem *htmp;
////            cout<<"get "<<node.size()<<" nodes and "<<nodebox.size()<<" boxes for contraction"<<endl;
//            while(!elemtmp->heap.empty()){ // means that exists w such as fk(w)<lower_up has no solution =>fk(w)>lower_ub
//                htmp = elemtmp->heap.pop();
//                stackheapelem.push_back(htmp);
//                NumConstraint cst(kp,ki,kd,Max12(kp,ki,kd,(htmp->box)[0])<=lower_ub);
////                NumConstraint cst(kp,ki,kd,Prob(kp,ki,(nodebox.back())[0])<=lower_ub);
//                CtcFwdBwd ctc(cst);
//                ctc.contract(elemtmp->box);
//                if(elemtmp->box.is_empty()){
////                    cout<<"empty for w_interest"<<endl;
//                    break;
//                }
//            }

////            cout<<"box after contraction: "<<elemtmp->box<<endl;
//            if(box != elemtmp->box) {
//                cout<<"contraction usefull, initbox: "<<box<<", after contraction: "<<elemtmp->box<<endl;
//                vol_rejected +=box.volume()-elemtmp->box.volume();
//                cout<<"loup : "<<lower_ub<<" get for point: "<<respoint<<" uplo: "<<upper_lb<< " volume rejected: "<<vol_rejected/IniboxK.volume()*100<<endl;
//            }
//            if(elemtmp->box.is_empty()){
//                elemtmp->heap.flush();
//                delete elemtmp;
//                continue;
//            }
//            else {
//                while(!stackheapelem.empty()){
//                    elemtmp->heap.push(stackheapelem.back());
//                    stackheapelem.pop_back();
//                }
//            }
//        }


//        *************** Mid point eval *********************

//        cout<<"stab: "<<midp<<endl;
//        stab = true;

//        6.86104);box[1]=Interval(0.0122183);box[2]=Interval(8.37996);box[3]=Interval(3.93941
//        midp[0] = 6.86104;midp[1]=0.0122183;midp[2]=8.37996;midp[3]=3.93941;

//        midp = sol;
//        str.preclw = 0.001;
//        if(elemtmp->box.max_diam()<0.1 && !stab) {
////            vol_rejected +=elemtmp->box.volume();
//            cout<<"No stab point found"<<endl;
////            cout<<"loup : "<<lower_ub<<" get for point: "<<respoint<<" uplo: "<<upper_lb<< " volume rejected: "<<vol_rejected/IniboxK.volume()*100<<endl;
////            delete elemtmp;
////            continue;
//        }

        if(stab) {// routh criterion ok for midpoint
//            Function max_mid(w,Max12(Interval(midp[0]),Interval(midp[1]),Interval(midp[2]),w));
//            Function max_mid(w,Prob(Interval(midp[0]),Interval(midp[1]),w));
//            str.f = &max_mid;
            str.preclw=0.00001;
            Heap<heap_elem> heap1(maxub);
            Heap<heap_elem> heap2(maxub);
            while(!elemtmp->heap.empty()) {
                heap_elem * he1=new heap_elem(*(elemtmp->heap.top()));
                heap1.push(he1);
                heap_elem * he2=new heap_elem(*(elemtmp->heap.top()));
                heap2.push(he2);
                elemtmp->heap.pop();
            }
            while(!heap2.empty()) {
                heap_elem * he2=new heap_elem(*(heap2.top()));
                elemtmp->heap.push(he2);
                heap2.pop();
            }
            double max = best_upper_bound_forall_w(&str,elemtmp,true,midp);
            while(!heap1.empty()) {
                heap_elem * he1=new heap_elem(*(heap1.top()));
                elemtmp->heap.push(he1);
                heap1.pop();
            }
            cout<<" max of mid point "<<midp<<": "<<max<<endl;
//            return 0;
            //            cout<<"fmax: "<<elemtmp->fmax<<"lower_ub: "<<lower_ub<<endl;

            if(max<lower_ub) {
                lower_ub = max;
                respoint = midp;
                cout<<"loup : "<<lower_ub<<" get for point: "<<respoint<<" uplo: "<<elemtmp->fmax.lb()<< " volume rejected: "<<vol_rejected/IniboxK.volume()*100<<endl;
            }
        }

        if(elemtmp->box.max_diam()>prec) {
            pair<IntervalVector,IntervalVector> boxes = lf.bisect(elemtmp->box);
            Heap<heap_elem> heap1(maxub);
            Heap<heap_elem> heap2(maxub);
            while(!elemtmp->heap.empty()) {
                heap_elem * he1=new heap_elem(*(elemtmp->heap.top()));
                heap1.push(he1);
                heap_elem * he2=new heap_elem(*(elemtmp->heap.top()));
                heap2.push(he2);
                elemtmp->heap.pop();
            }

            list.push(new list_elem(boxes.first,heap1,elemtmp->fmax));
            list.push(new list_elem(boxes.second,heap2,elemtmp->fmax));
            delete elemtmp;
        }
        else{
            cout<<"minprec reached! "<<endl;
            cout<<"box: "<<elemtmp->box<<endl;
            lbfix = true;
            volout += elemtmp->box.volume();
            cout<<"volume treated: "<<(volout+vol_rejected)/IniboxK.volume()*100<<endl;
            //cout<<"loup : "<<lower_ub<<" get for point: "<<respoint<<" uplo: "<<elemtmp->fmax.lb()<< " volume rejected: "<<vol_rejected/IniboxK.volume()*100<<endl;
            delete elemtmp;
//            return 0;
        }
    }
    Timer::stop();
    cout<<"maximum: "<<sqrt(lower_ub)<<"for point: "<<respoint<<endl;
    cout<<" uplo: "<< sqrt(upper_lb)<<endl;
    cout<<"result found in "<<Timer::VIRTUAL_TIMELAPSE()<<endl;

    delete elemtmp;
    list.flush();

    return 0;


}
/*
 * maximum: 2.16987for point: (0.470296 ; 0.468828 ; -0.426464)
 uplo: 1.75584
result found in 952.533
*/
/*
 maximum: 1.90531for point: (0.528457 ; 0.527102 ; -0.345843)
 uplo: 1.8792
result found in 806.181
*/


