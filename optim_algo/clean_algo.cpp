

#include "/home/kiwi/Tools/ibex-lib-ibex-2.1.18/OUT/include/ibex/ibex.h"
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
       res += modulus.at(i).first->eval_bernstein(box)[0]/modulus.at(i).second->eval_bernstein(box)[0];
    res = ibex::sqrt(res);
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

    /*Variable kp,ki,kd,u;
    NumConstraint infub(kp,ki,kd,u,(*(str->f))(kp,ki,kd,u)<=1000000);
    CtcFwdBwd ctc1(infub);*/


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
        if(midp)
            cout<<"iteration: "<<nbit<<endl;
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
            if (cres[0].lb()<str->lower_ub)
            {
//                ctc1.contract(box);
            }
            if(!ctcbox.is_empty())
            {
                //            cout<<"res for box "<<box<<" : "<<res<<endl;
//                IntervalVector round_box(elem->box.size()+1);
//                round_box[0] = ibex::pow(Interval(10),elemtmp->box[0]);
//                for(unsigned i=1;i<=elem->box.size();i++) {
//                    round_box[i] = box[i-1];
//                }
//                res[0] =str->bernwfree.at(0).eval_bernstein(round_box)[0];//-str->berncst.at(0).eval_bernstein(IntervalVector(1,round_box[0]))[0];
//                for(unsigned i=1;i<str->bernwfree.size();i++)
//                    res[0] = max(res[0],str->bernwfree.at(i).eval_bernstein(round_box)[0]);//-str->berncst.at(i).eval_bernstein(IntervalVector(1,round_box[0]))[0]);
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
//                ctc1.contract(ctcbox);
//                cout<<" after contraction: "<<ctcbox<<endl;
            }
            if(ctcbox!=box)
                cout<<"contraction hit! box: "<<box<<" ctcbox: "<<ctcbox<<endl;
            if(!ctcbox.is_empty())
            {
//                //****** mid w eval *******
                box[elem->box.size()] = elemtmp->box.mid()[0];
                IntervalVector cresmid = str->f->eval_vector(box);
//////                //            cout<<"resmid for box "<<box<<" : "<<cresmid[0]<<endl;
//                box[elem->box.size()] = ibex::pow(Interval(10),box[elem->box.size()]);
//                resmid[0] =str->bernwset.at(0).eval_bernstein(box)[0];//-str->berncst.at(0).eval_bernstein(IntervalVector(1,box[elem->box.size()]))[0];
//                for(unsigned i=1;i<str->bernwset.size();i++)
//                    resmid[0] = max(resmid[0],str->bernwset.at(i).eval_bernstein(box)[0]);//-str->berncst.at(i).eval_bernstein(IntervalVector(1,box[elem->box.size()]))[0]);
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

   // cout<<"lb: "<<lower_ub<<" ,ub: "<<upper_ub<<endl;
    return upper_ub;
}

Vector randpt(IntervalVector box, int size) {
    Vector pt(size);
    for(unsigned i=0;i<size;i++) {
        int rnd = std::rand();
        pt[i] = double(rnd)/double(RAND_MAX)*box[i].diam()+box[i].lb();
    }
    return pt;
}

bool pol_position_ok(pair<icomplex,Interval> * eigval, int nbeig){

    Interval zero(-0.01,0.01);
    Interval diameter;
    for(int i=0;i<nbeig;i++) {
        diameter = 2*eigval[i].second-eigval[i].second;
        //cout<<"center: ("<<eigval[i].first.real<<","<<eigval[i].first.real<<"), diameter: "<<eigval[i].second<<endl;
        if((eigval[i].first.real+diameter).is_subset(zero) && (eigval[i].first.imag+diameter).is_subset(zero)) { // zero pole
            //cout<<"zero rejection"<<endl;
            continue;}
        if(((((eigval[i].first.real+diameter) + (eigval[i].first.imag+diameter)).lb() > 0) || // out of pole area
                ((-(eigval[i].first.real+diameter) + (eigval[i].first.imag+diameter)).ub() < 0) ||
                (eigval[i].first.real -eigval[i].second).lb() >0) )
            return false;
    }
    //cout<<"pole not out"<<endl;
    return true;

}

int main() {

    std::srand(std::time(0));

    double prec(0.01);
    double preclwmin = 1.e-3;// dynamic initialization in loop
    double stop_criterion(0.1); // stop if distance between uplo and loup lower than stop_criterion

    // init boxes

    IntervalVector IniboxK(4,Interval(0,2));
   // IniboxK[0] = Interval(10,10);
   // IniboxK[1] = Interval(9.5,9.5);
   // IniboxK[2] = Interval(9.5,9.5);
    IniboxK[3] = Interval(0,5);
    IntervalVector Inilw(1,Interval(-3,3));

    //function definition
    Variable ki,kp,kd,tf,w,u;


    //pid tf free
    Function Twz1w("kp","ki","kd","tf","w","(w^2*(w^2 + 1.0)*(w^2 + 10000.0)*(25.0*w^4 - 1.0*w^2 + 25.0))/((10000.0*w^2 + 1.0)*(25.0*kd^2*w^4 + 25.0*kp^2*w^2 + 25.0*kp^2*w^4 - 1.0*ki*(50.0*kd*w^2 + 70.0*w^2 + 70.0*w^4) + ki^2*(25.0*w^2 + 25.0) + 120.0*kd*w^4 - 50.0*kd*w^6 + 50.0*kp*w^2 - 50.0*kp*w^6 + 25.0*w^2 + 24.0*w^4 + 24.0*w^6 + 25.0*w^8 + 50.0*kd*kp*w^4))");
    Function Twz2w("kp","ki","kd","tf","w", "((100.0*w^2 + 1.0)*(25.0*w^4 - 1.0*w^2 + 25.0)*(kd^2*w^4 - 2.0*kd*ki*w^2 + 2.0*kd*kp*w^4 + ki^2*w^2 + ki^2 + kp^2*w^4 + kp^2*w^2))/((w^2 + 100.0)*(25.0*kd^2*w^4 - 50.0*kd*ki*w^2 + 50.0*kd*kp*w^4 - 50.0*kd*w^6 + 120.0*kd*w^4 + 25.0*ki^2*w^2 + 25.0*ki^2 - 70.0*ki*w^4 - 70.0*ki*w^2 + 25.0*kp^2*w^4 + 25.0*kp^2*w^2 - 50.0*kp*w^6 + 50.0*kp*w^2 + 25.0*w^8 + 24.0*w^6 + 24.0*w^4 + 25.0*w^2))");
    Function Twz3w("kp","ki","kd","tf","w", "((kd*w^3-10*w*kp-11*w*ki+w^3*kp)^2+(10*kd*w^2-10*ki+w^2*ki+11*w^2*kp)^2)*1/((-(12.4)*w^2-kd*w^2+ki-10*w^2*ki-11*w^2*kp+(25.0)*w^4)^2+(10*kd*w^3+(26.4)*w^3-w*kp-w-11*w*ki+10*w^3*kp-10*w^5)^2)");
    cout<<"w func ok"<<endl;
   // Function Twz1u("kp","ki","kd","u","tf","(((0.2479)*exp(ln(10)*5*u)*tf-1.128068*exp(ln(10)*3*u)-0.828*exp(ln(10)*3*u)*tf)^2+(-1.128068*exp(ln(10)*4*u)*tf+0.828*exp(ln(10)*2*u)-0.2479*exp(ln(10)*4*u))^2)*1/(((0.0046)*exp(ln(10)*u)*kp-exp(ln(10)*3*u)*kd+0.4958*exp(ln(10)*5*u)*tf-1.80228068*exp(ln(10)*3*u)-exp(ln(10)*3*u)*tf*kp-0.0082800000000000000005*exp(ln(10)*3*u)*tf+0.0046*exp(ln(10)*u)*ki*tf+exp(ln(10)*u)*ki)^2+(-1.80228068*exp(ln(10)*4*u)*tf+exp(ln(10)*2*u)*ki*tf+0.0046*exp(ln(10)*2*u)*tf*kp-0.0046*ki+0.0082800000000000000005*exp(ln(10)*2*u)+0.0046*exp(ln(10)*2*u)*kd-0.4958*exp(ln(10)*4*u)+exp(ln(10)*2*u)*kp)^2)");

// linearization loop
/*    Function Twz1u("kp","ki","kd","tf","u","((-0.828*exp(ln(10)*3*u)*tf+0.2479*exp(ln(10)*5*u)*tf-1.128068*exp(ln(10)*3*u))^2+(-0.2479*exp(ln(10)*4*u)+0.828*exp(ln(10)*2*u)-1.128068*exp(ln(10)*4*u)*tf)^2)*1/((-0.0082800000000000000005*exp(ln(10)*3*u)*tf+exp(ln(10)*u)*ki+0.0046*exp(ln(10)*u)*ki*tf-exp(ln(10)*3*u)*tf*kp+0.0046*exp(ln(10)*u)*kp-exp(ln(10)*3*u)*kd+0.4958*exp(ln(10)*5*u)*tf-1.80228068*exp(ln(10)*3*u))^2+(exp(ln(10)*2*u)*kp-0.4958*exp(ln(10)*4*u)-0.0046*ki+0.0046*exp(ln(10)*2*u)*kd+exp(ln(10)*2*u)*ki*tf+0.0082800000000000000005*exp(ln(10)*2*u)-1.80228068*exp(ln(10)*4*u)*tf+0.0046*exp(ln(10)*2*u)*tf*kp)^2)");

Function Twz2u("kp","ki","kd","tf","u","1/((-1.8*exp(ln(10)*3*u)*tf+exp(ln(10)*u)*kp+exp(ln(10)*u)*ki*tf-0.4958*exp(ln(10)*3*u))^2+(exp(ln(10)*2*u)*tf*kp-ki-0.4958*exp(ln(10)*4*u)*tf+1.8*exp(ln(10)*2*u)+exp(ln(10)*2*u)*kd)^2)*(((0.017999999999999999999)*exp(ln(10)*2*u)*ki*tf-0.0049579999999999999997*exp(ln(10)*4*u)*kd+0.017999999999999999999*exp(ln(10)*2*u)*kp-0.0049579999999999999997*exp(ln(10)*4*u)*tf*kp+0.0049579999999999999997*exp(ln(10)*2*u)*ki)^2+(-0.0049579999999999999997*exp(ln(10)*3*u)*kp+0.017999999999999999999*exp(ln(10)*u)*ki-0.0049579999999999999997*exp(ln(10)*3*u)*ki*tf-0.017999999999999999999*exp(ln(10)*3*u)*tf*kp-0.017999999999999999999*exp(ln(10)*3*u)*kd)^2)");


    Function Twz3u("kp","ki","kd","tf","u","1/((ki-exp(ln(10)*2*u)*kd-1.8*exp(ln(10)*2*u)+0.4958*exp(ln(10)*4*u)*tf-exp(ln(10)*2*u)*tf*kp)^2+(-1.8*exp(ln(10)*3*u)*tf+exp(ln(10)*u)*ki*tf+exp(ln(10)*u)*kp-0.4958*exp(ln(10)*3*u))^2)*((100.0)*exp(ln(10)*4*u)*tf^2+100.0*exp(ln(10)*2*u))");
*/
    Function Twz1u("../symbolic/Twz1.txt");
    Function Twz2u("../symbolic/Twz2.txt");
    Function Twz3u("../symbolic/Twz3.txt");
    Function Twz4u("../symbolic/Twz4.txt");
    Function Twz1(kp,ki,kd,tf,u,ibex::sqrt(Twz1u(kp,ki,kd,tf,u)));
    Function Twz2(kp,ki,kd,tf,u,ibex::sqrt(Twz2u(kp,ki,kd,tf,u)));
    Function Twz3(kp,ki,kd,tf,u,ibex::sqrt(Twz3u(kp,ki,kd,tf,u)));
    Function Twz4(kp,ki,kd,tf,u,ibex::sqrt(Twz4u(kp,ki,kd,tf,u)));

// direct linearization 
    Function Twz13u("../symbolic/Twz13.txt");
    Function Twz23u("../symbolic/Twz23.txt");
    Function Twz33u("../symbolic/Twz33.txt");
    Function Twz43u("../symbolic/Twz43.txt");
    cout<<"ok 3"<<endl;
    Function Twz12u("../symbolic/Twz12.txt");
    Function Twz22u("../symbolic/Twz22.txt");
    Function Twz32u("../symbolic/Twz32.txt");
    Function Twz42u("../symbolic/Twz42.txt");
    cout<<"ok 2"<<endl;
    Function Twz11u("../symbolic/Twz11.txt");
    Function Twz21u("../symbolic/Twz21.txt");
    Function Twz31u("../symbolic/Twz31.txt");
    Function Twz41u("../symbolic/Twz41.txt");
    cout<<"ok 1"<<endl;
    Function Twz10u("../symbolic/Twz10.txt");
    Function Twz20u("../symbolic/Twz20.txt");
    Function Twz30u("../symbolic/Twz30.txt");
    Function Twz40u("../symbolic/Twz40.txt");
    cout<<"ok 0"<<endl;


    cout<<"u func ok"<<endl;
    Function Twz13(kp,ki,kd,tf,u,ibex::sqrt(Twz13u(kp,ki,kd,tf,u)));
    Function Twz23(kp,ki,kd,tf,u,ibex::sqrt(Twz23u(kp,ki,kd,tf,u)));
    Function Twz33(kp,ki,kd,tf,u,ibex::sqrt(Twz33u(kp,ki,kd,tf,u)));
    Function Twz43(kp,ki,kd,tf,u,ibex::sqrt(Twz43u(kp,ki,kd,tf,u)));
    Function Twz12(kp,ki,kd,tf,u,ibex::sqrt(Twz12u(kp,ki,kd,tf,u)));
    Function Twz22(kp,ki,kd,tf,u,ibex::sqrt(Twz22u(kp,ki,kd,tf,u)));
    Function Twz32(kp,ki,kd,tf,u,ibex::sqrt(Twz32u(kp,ki,kd,tf,u)));
    Function Twz42(kp,ki,kd,tf,u,ibex::sqrt(Twz42u(kp,ki,kd,tf,u)));
    Function Twz11(kp,ki,kd,tf,u,ibex::sqrt(Twz11u(kp,ki,kd,tf,u)));
    Function Twz21(kp,ki,kd,tf,u,ibex::sqrt(Twz21u(kp,ki,kd,tf,u)));
    Function Twz31(kp,ki,kd,tf,u,ibex::sqrt(Twz31u(kp,ki,kd,tf,u)));
    Function Twz41(kp,ki,kd,tf,u,ibex::sqrt(Twz41u(kp,ki,kd,tf,u)));
    Function Twz10(kp,ki,kd,tf,u,ibex::sqrt(Twz10u(kp,ki,kd,tf,u)));
    Function Twz20(kp,ki,kd,tf,u,ibex::sqrt(Twz20u(kp,ki,kd,tf,u)));
    Function Twz30(kp,ki,kd,tf,u,ibex::sqrt(Twz30u(kp,ki,kd,tf,u)));
    Function Twz40(kp,ki,kd,tf,u,ibex::sqrt(Twz40u(kp,ki,kd,tf,u)));


    Function Max1(kp,ki,kd,tf,u,ibex::max(Twz13(kp,ki,kd,tf,u),Twz23(kp,ki,kd,tf,u)));
    Function Max2(kp,ki,kd,tf,u,ibex::max(Max1(kp,ki,kd,tf,u),Twz33(kp,ki,kd,tf,u)));
    Function Max3(kp,ki,kd,tf,u,ibex::max(Max2(kp,ki,kd,tf,u),Twz12(kp,ki,kd,tf,u)));
    Function Max4(kp,ki,kd,tf,u,ibex::max(Max3(kp,ki,kd,tf,u),Twz22(kp,ki,kd,tf,u)));
    Function Max5(kp,ki,kd,tf,u,ibex::max(Max4(kp,ki,kd,tf,u),Twz32(kp,ki,kd,tf,u)));
    Function Max6(kp,ki,kd,tf,u,ibex::max(Max5(kp,ki,kd,tf,u),Twz11(kp,ki,kd,tf,u)));
    Function Max7(kp,ki,kd,tf,u,ibex::max(Max6(kp,ki,kd,tf,u),Twz21(kp,ki,kd,tf,u)));
    Function Max8(kp,ki,kd,tf,u,ibex::max(Max7(kp,ki,kd,tf,u),Twz31(kp,ki,kd,tf,u)));
    Function Max9(kp,ki,kd,tf,u,ibex::max(Max8(kp,ki,kd,tf,u),Twz10(kp,ki,kd,tf,u)));
    Function Max10(kp,ki,kd,tf,u,ibex::max(Max9(kp,ki,kd,tf,u),Twz20(kp,ki,kd,tf,u)));

    Function Max11(kp,ki,kd,tf,u,ibex::max(Max9(kp,ki,kd,tf,u),Twz43(kp,ki,kd,tf,u)));
    Function Max12(kp,ki,kd,tf,u,ibex::max(Max9(kp,ki,kd,tf,u),Twz42(kp,ki,kd,tf,u)));
    Function Max13(kp,ki,kd,tf,u,ibex::max(Max9(kp,ki,kd,tf,u),Twz41(kp,ki,kd,tf,u)));
    Function Max123(kp,ki,kd,tf,u,ibex::max(Max9(kp,ki,kd,tf,u),Twz40(kp,ki,kd,tf,u)));


    //Function Max123(kp,ki,kd,tf,u,ibex::max(Max10(kp,ki,kd,tf,u),Twz30(kp,ki,kd,tf,u)));




   // Function Max123(kp,ki,kd,tf,u,ibex::max(Max12(kp,ki,kd,tf,u),Twz3(kp,ki,kd,tf,u)));
  //  Function Max123(kp,ki,kd,tf,u,Twz1(kp,ki,kd,tf,u));

 Max123.eval_vector(IntervalVector(5,Interval(1)));
cout<<"max ok"<<endl;

    vector<bernfunc> bernwset;
    vector<bernfunc> bernwfree;
    vector<bernfunc> berncst;

   /* vector< pair<polynomial *,polynomial *> > frac_wfree;
    vector< pair<polynomial *,polynomial *> > frac_wset;

    string aff_var[4];aff_var[0] = "kp";aff_var[1] = "ki";aff_var[2] = "kd";aff_var[3] = "tf";
    string var[4];var[0]="w";
    frac_wset.push_back(get_pol(&Twz1w,aff_var,3,var,1));
    bernwset.push_back(bernfunc(frac_wset));
    frac_wset.clear();
    frac_wset.push_back(get_pol(&Twz2w,aff_var,3,var,1));
    bernwset.push_back(bernfunc(frac_wset));
    aff_var[0] = "w";
    var[0] = "kp";var[1] = "ki";var[2] = "kd";var[3] = "tf";
    frac_wfree.push_back(get_pol(&Twz1w,aff_var,1,var,3));
    bernwfree.push_back(bernfunc(frac_wfree));
    frac_wfree.clear();
    frac_wfree.push_back(get_pol(&Twz2w,aff_var,1,var,3));
    bernwfree.push_back(bernfunc(frac_wfree));*/
    cout<<"func ok"<<endl;

//    Function rfunc("../symbolic/routh.txt");


//    //Function rfunc(kp,ki,kd,tf,Return((0.4958)*tf,0.4958+(1.8)*tf,(0.4958+(1.8)*tf)*(0.89243999999999999994+(0.4958)*kd+(1.8)*kp*tf+(3.2399999999999999998)*tf-(0.4958)*ki*tf+(1.8)*kd*tf),((3.2399999999999999998)*kp*tf+(0.89243999999999999994)*kp-(0.4958)*kp*ki*tf+(1.8)*kp*kd*tf+(1.8)*kp*ki*tf-(0.4958)*ki*tf-(0.89243999999999999994)*ki*tf+(1.8)*kp*tf+(0.4958)*kd*ki*tf-(0.24581763999999999999)*ki+(0.4958)*kp*kd+(1.8)*kd*ki*tf)*(0.89243999999999999994+(0.4958)*kd+(1.8)*kp*tf+(3.2399999999999999998)*tf-(0.4958)*ki*tf+(1.8)*kd*tf),ki)  );
//    Vector mask(rfunc.image_dim());
//    mask[0] = 1;
//    Function rS1(kp,ki,kd,tf, rfunc(kp,ki,kd,tf)*mask);
//    mask[0] = 0;
//    mask[1] = 1;
//    Function rS2(kp,ki,kd,tf,rfunc(kp,ki,kd,tf)*mask);
//    mask[1] = 0;
//    mask[2] = 1;
//    Function rS3(kp,ki,kd,tf,rfunc(kp,ki,kd,tf)*mask);
//    mask[2] = 0;
//    mask[3] = 1;
//    Function rS4(kp,ki,kd,tf,rfunc(kp,ki,kd,tf)*mask);
//    mask[3] = 0;
//    mask[4] = 1;
//    Function rS5(kp,ki,kd,tf,rfunc(kp,ki,kd,tf)*mask);
//    cout<<rS5.eval(IntervalVector(4,Interval(1)))<<endl;
//    cout<<"vect mult"<<mask*mask<<endl;
//    //Function rS6(kp,ki,kd,tf,rfunc(kp,ki,kd,tf)[5]);
//    cout<<"rfunc ok"<<endl;

    Function rS1("kp","ki","kd","tf","tf");
    Function rS2("kp","ki","kd","tf","1+0.05*tf");
    Function rS3("kp","ki","kd","tf","(0.05-6.725*ki*tf^2+0.33625*kd*tf+6.725*kd+0.33625*kp*tf^2+0.0025000000000000000002*tf)*(1+0.05*tf)^(-1)");
    Function rS4("kp","ki","kd","tf","((0.33625)*kp+2.26128125*kd*ki*tf^2+45.225625*kp*kd+2.26128125*kp^2*tf^2-0.33625*ki*tf+45.225625*kd*ki*tf-6.725*ki-45.225625*ki^2*tf^3+2.26128125*kp*ki*tf^3+0.016812500000000000001*kp*tf-45.225625*kp*ki*tf^2+2.26128125*kp*kd*tf)*(0.05-6.725*ki*tf^2+0.33625*kd*tf+6.725*kd+0.33625*kp*tf^2+0.0025000000000000000002*tf)^(-1)");
    Function rS5("kp","ki","kd","tf","6.725*ki");
    Function rS6("kp","ki","kd","tf","1+tf");
    Function rS7("kp","ki","kd","tf","(1+tf)^(-1)*(1+6.725*kd+6.725*kp*tf^2+tf+6.725*kd*tf-6.725*ki*tf^2)");
    Function rS8("kp","ki","kd","tf","((45.225625)*kp*kd*tf+6.725*kp-45.225625*kp*ki*tf^2+6.725*kp*tf+45.225625*kp*ki*tf^3-45.225625*ki^2*tf^3+45.225625*kd*ki*tf-6.725*ki*tf+45.225625*kp^2*tf^2-6.725*ki+45.225625*kp*kd+45.225625*kd*ki*tf^2)*(1+6.725*kd+6.725*kp*tf^2+tf+6.725*kd*tf-6.725*ki*tf^2)^(-1)");
    Function rS9("kp","ki","kd","tf","1+3*tf");
    Function rS10("kp","ki","kd","tf","(1+3*tf)^(-1)*(3+20.175*kp*tf^2+6.725*kd-6.725*ki*tf^2+20.175*kd*tf+9*tf)");
    Function rS11("kp","ki","kd","tf","(3+20.175*kp*tf^2+6.725*kd-6.725*ki*tf^2+20.175*kd*tf+9*tf)^(-1)*((20.175)*kp+45.225625*kd*ki*tf+135.67687499999999999*kd*ki*tf^2+60.524999999999999998*kp*tf-45.225625*ki^2*tf^3+45.225625*kp*kd+135.67687499999999999*kp*ki*tf^3-6.725*ki-45.225625*kp*ki*tf^2+135.67687499999999999*kp*kd*tf+135.67687499999999999*kp^2*tf^2-20.175*ki*tf)");
    Function rS12("kp","ki","kd","tf","1+4*tf");
    Function rS13("kp","ki","kd","tf","(1+4*tf)^(-1)*(4+6.725*kd+26.9*kd*tf-6.725*ki*tf^2+16*tf+26.9*kp*tf^2)");
    Function rS14("kp","ki","kd","tf","((180.9025)*kp*ki*tf^3+26.9*kp+180.9025*kp^2*tf^2-26.9*ki*tf+180.9025*kp*kd*tf-45.225625*kp*ki*tf^2+45.225625*kp*kd+107.6*kp*tf-6.725*ki+45.225625*kd*ki*tf-45.225625*ki^2*tf^3+180.9025*kd*ki*tf^2)*(4+6.725*kd+26.9*kd*tf-6.725*ki*tf^2+16*tf+26.9*kp*tf^2)^(-1)");

    vector<Function *> routh_func;
    routh_func.push_back(&rS1);
    routh_func.push_back(&rS2);
    routh_func.push_back(&rS3);
    routh_func.push_back(&rS4);
    routh_func.push_back(&rS5);
    routh_func.push_back(&rS6);
    routh_func.push_back(&rS7);
    routh_func.push_back(&rS8);
    routh_func.push_back(&rS9);
    routh_func.push_back(&rS10);
    routh_func.push_back(&rS11);
    routh_func.push_back(&rS12);
    routh_func.push_back(&rS13);
    routh_func.push_back(&rS14);
    cout<<"contraintes ok"<<endl;

    NumConstraint *c1= new NumConstraint(kp,ki,kd,tf,rS1(kp,ki,kd,tf)>=0);
    NumConstraint *c2= new NumConstraint(kp,ki,kd,tf,rS2(kp,ki,kd,tf)>=0);
    NumConstraint *c3= new NumConstraint(kp,ki,kd,tf,rS3(kp,ki,kd,tf)>=0);
    NumConstraint *c4= new NumConstraint(kp,ki,kd,tf,rS4(kp,ki,kd,tf)>=0);
    NumConstraint *c5= new NumConstraint(kp,ki,kd,tf,rS5(kp,ki,kd,tf)>=0);
    NumConstraint *c6= new NumConstraint(kp,ki,kd,tf,rS6(kp,ki,kd,tf)>=0);
    NumConstraint *c7= new NumConstraint(kp,ki,kd,tf,rS7(kp,ki,kd,tf)>=0);
    NumConstraint *c8= new NumConstraint(kp,ki,kd,tf,rS8(kp,ki,kd,tf)>=0);
    NumConstraint *c9= new NumConstraint(kp,ki,kd,tf,rS9(kp,ki,kd,tf)>=0);
    NumConstraint *c10= new NumConstraint(kp,ki,kd,tf,rS10(kp,ki,kd,tf)>=0);
    NumConstraint *c11= new NumConstraint(kp,ki,kd,tf,rS11(kp,ki,kd,tf)>=0);
    NumConstraint *c12= new NumConstraint(kp,ki,kd,tf,rS12(kp,ki,kd,tf)>=0);
    NumConstraint *c13= new NumConstraint(kp,ki,kd,tf,rS13(kp,ki,kd,tf)>=0);
    NumConstraint *c14= new NumConstraint(kp,ki,kd,tf,rS14(kp,ki,kd,tf)>=0);

    cout<<"cst ok"<<endl;
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
    fac.add_ctr(*c8);
    fac.add_ctr(*c9);
    fac.add_ctr(*c10);
    fac.add_ctr(*c11);
    fac.add_ctr(*c12);
    fac.add_ctr(*c13);
    fac.add_ctr(*c14);

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
 //Ctc* ctc_final = ctc;
    cout<<"init ctc final ok"<<endl;
    string var[4]; var[0] = "kp"; var[1] = "ki"; var[2] = "kd"; var[3] = "tf";
    univar_polynomial cltf_pol("(6.725*((kp*tf+kd)*s*s+(kp+ki*tf)*s+ki)+s^4*tf+(0.5+tf)*s^3+0.5*s^2)*(s^4*tf+(0.5+tf)*s^3+0.5*s^2)","s",var,4);
    cout<<"univar polynomial ok"<<endl;

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

    optiw str(lower_ub,preclwmin,&Max123,bernwset,bernwfree,berncst,false);


    double vol_rejected(0),volout(0);
    list_elem * elemtmp;

   // double freqeps=0.001;
    //double umin(-5),umax(5);
   // Interval freq(umin,umin+freqeps);
   // IntervalVector box(5),boxmat(5);
////    box[0] = Interval(0.697786);box[1]=Interval(0.531368);box[2]=Interval(0.511734);
   // box[0] = 8.66596;box[1] = 3.15248;box[2] = 0.0138283;box[3] = 0.0492056;
   // boxmat[0] = Interval(0.505948825233638);boxmat[1]=Interval(1.002037363530295);boxmat[2]=Interval(-0.496914467315547);

   // IntervalVector cdraw(2);
   // vibes::beginDrawing();
   // vibes::newFigure("Twz1");
   // while(freq.ub()<umax) {
   //     cdraw[0] = freq;
   //     box[4]=freq;
//        boxmat[3]=freq;

   //     cdraw[1] = Twz1u.eval_vector(box)[0];
   //     vibes::drawBox(cdraw,"blue[blue]");
//        cdraw[1] = 20/ibex::log(10)*ibex::log(Twz2.eval_vector(box)[0]);
//        vibes::drawBox(cdraw,"red[red]");
//        cdraw[1] = 20/ibex::log(10)*ibex::log(Twz3.eval_vector(box)[0]);
//        vibes::drawBox(cdraw,"green[green]");
//        cdraw[1] = 20/ibex::log(10)*ibex::log(Max123.eval_vector(box)[0]);
//        vibes::drawBox(cdraw,"black[black]");
////        cdraw[1] = Max123.eval_vector(boxmat)[0];
////        vibes::drawBox(cdraw,"yellow[yellow]");
   //     freq+=Interval(freqeps);
   // }
   // return 0;
    int nbeig(8);
    pair<icomplex,Interval> roots[nbeig];


    bool lbfix(false);
    Timer::start();
    while(!list.empty()) {
        elemtmp = list.pop();
        if(!lbfix)
            upper_lb = elemtmp->fmax[0].lb();
        if(lower_ub-elemtmp->fmax[0].lb()<stop_criterion) { // stop creterion reached
            break;
        }
      //  cout<<"maxdiam: "<<elemtmp->box.max_diam()<<endl;
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
        // ********************* pole placement condition ***************
        if(cltf_pol.get_roots(elemtmp->box,roots,nbeig) == 1) {
            if(!pol_position_ok(roots,nbeig)) {
                vol_rejected += elemtmp->box.volume();
                cout<<"pole placement condition rejection  "<<endl;
                cout<<"loup : "<<lower_ub<<" get for point: "<<respoint<<" uplo: "<<upper_lb<< " volume rejected: "<<vol_rejected/IniboxK.volume()*100<<endl;
            }
        }
        else {
            cout<<"unable to compute eigenvalues"<<endl;
        }

        //******************* search for a stable point in the box **************
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
            if(cltf_pol.get_roots(IntervalVector(midp),roots,nbeig) == 1)
                stab&= pol_position_ok(roots,nbeig);
            else {
                cout<<"unable to compute eigenvalues"<<endl;
                stab = false;
            }
            //cout<<"point: "<<midp<<endl;
           /* for(int i=0;i<nbeig;i++) {
                cout<<"center: ("<<roots[i].first.real<<","<<roots[i].first.imag<<"), diameter: "<<roots[i].second<<endl;
            }*/
           // cout<<"*****************************"<<endl;
            if(stab)
            {

                break;}

        }
        // **************** compute (wmax,fmax) *******************
//        Function g(w,Max12(elemtmp->box[0],elemtmp->box[1],elemtmp->box[2],w));
//        Function g(w,Prob(elemtmp->box[0],elemtmp->box[1],w));
//        str.f = &g;


//        elemtmp->computable = stab;
        if(elemtmp->computable)
            best_upper_bound_forall_w(&str,elemtmp,false,Vector(4));
       // cout<<"res found for box "<<elemtmp->box<<" : "<<endl<<"        "<<elemtmp->fmax<<endl;
        if(elemtmp->fmax[0] == Interval::EMPTY_SET) { // maximum is guaranteed to be higher than the current lower_up for box elemtmp->box
            vol_rejected+= elemtmp->box.volume();
            cout<<"solution deletion,  "<<endl;
            cout<<"loup : "<<lower_ub<<" get for point: "<<respoint<<" uplo: "<<upper_lb<< " volume rejected: "<<vol_rejected/IniboxK.volume()*100<<endl;
            delete elemtmp;
            continue;
        }


        if(stab) {// routh criterion ok for midpoint
//            Function max_mid(w,Max12(Interval(midp[0]),Interval(midp[1]),Interval(midp[2]),w));
//            Function max_mid(w,Prob(Interval(midp[0]),Interval(midp[1]),w));
//            str.f = &max_mid;
//            midp[0] = 0.073613386821018;midp[1] = 0.096927953568474;midp[2] = 0.030496647684430; // matlab sol
   //         midp[0] = 8.66596;midp[1] = 3.15248;midp[2] = 0.0138283;midp[3] = 0.0492056; //our sol
            str.preclw=0.001;
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
         //   cout<<" max of mid point "<<midp<<": "<<max<<endl;
        //    return 0;
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
    cout<<"maximum: "<<lower_ub<<"for point: "<<respoint<<endl;
    cout<<" uplo: "<< upper_lb<<endl;
    cout<<"result found in "<<Timer::VIRTUAL_TIMELAPSE()<<endl;
    delete elemtmp;
    list.flush();

//    //****************** plot results **********
//    double freqeps=0.001;
//    double umin(-3),umax(3);
//    Interval freq(umin,umin+freqeps);
//    IntervalVector box(4),boxmat(4);

//    box[0] = Interval(0.697786);box[1]=Interval(0.531368);box[2]=Interval(0.511734);
//    box = respoint;
//    boxmat[0] = Interval(2.729926145077233);boxmat[1]=Interval(5.800988994294283e-07);boxmat[2]=Interval(2.889225437571279);

//    IntervalVector cdraw(2);
//    vibes::beginDrawing();
//    vibes::newFigure("Matlab solution in blue, Algo solution in red");
//    while(freq.ub()<umax) {
//        cdraw[0] = freq;
//        box[3]=freq;
//        boxmat[3]=freq;

//        cdraw[1] = Twz1u.eval_vector(box)[0];
//        vibes::drawBox(cdraw,"blue[blue]");
//        cdraw[1] = Twz2u.eval_vector(box)[0];
//        vibes::drawBox(cdraw,"red[red]");
//        cdraw[1] = Twz3u.eval_vector(box)[0];
//        vibes::drawBox(cdraw,"green[green]");
////        cdraw[1] = Max123.eval_vector(boxmat)[0];
////        vibes::drawBox(cdraw,"red[red]");
//        freq+=Interval(freqeps);
//    }

//    return 0;


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


