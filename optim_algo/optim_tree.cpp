

#include "/home/kiwi/Tools/ibex-lib-ibex-2.1.18/OUT/include/ibex/ibex.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "src/ibex_SetIntervalReg.h"
#include "src/polynomial.h"
#include "src/vibes.h"
#include <cstdlib>
#include <ctime>


/*
 * #ALGO4
 *
 * Differences with #ALGO2:
 * Instead of considering fixed sized slices of w, we run sivia on w for optim, with a precision that depends on the size of boxe K.
*/


using namespace ibex;
using namespace std;

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

// element of the heap used in slave B&B
class heap_elem {
public:
    SetNodeReg * node; // node of the tree that correspond to the box
    IntervalVector box; // box of slave B&B, correspond to a frequency interval
    Interval eval; // evaluation of objectif function with box
    bool pass;

    heap_elem(SetNodeReg * node, IntervalVector box, Interval eval,bool pass);
};
heap_elem::heap_elem(SetNodeReg * node, IntervalVector box, Interval eval,bool pass): node(node),box(box),eval(eval),pass(pass)
{}

class list_elem{ // element of the master B&B algo
public:
    IntervalVector box; // box of controller parameters
    IntervalVector fmax; // enclosure of fmax(box), see MACIS presentation for detail
    vector<SetNodeReg*> wint;
    list_elem(IntervalVector box,vector<SetNodeReg*> wint,IntervalVector fmax); // constructor

};
list_elem::list_elem(IntervalVector box,vector<SetNodeReg*> wint,IntervalVector fmax):box(box),
    wint(wint),fmax(fmax)
{} // constructor implementation


//*************** template function definition for heap use*********************

// penalize lower bound, use in master B&B to always deal with the lower fmax.lb that correspond to the uplo
class cost_uplo: public CostFunc<list_elem> {
public:
    virtual double cost(const list_elem& elem) const;
};

double cost_uplo::cost(const list_elem& elem) const {
    return elem.fmax[0].lb();
}

// penalize box size, not used
class cost_box: public CostFunc<list_elem> {
public:
    virtual double cost(const list_elem& elem) const;
};

double cost_box::cost(const list_elem& elem) const {
    return elem.box.volume();
}
// penalized lower ub, not used
class costf1 : public CostFunc<heap_elem> {
public:
    virtual double cost(const heap_elem& elem) const;

};
double costf1::cost(const heap_elem& elem) const {
    return elem.eval.ub();
}
// penalized upper ub, use in slave B&B
class costf2 : public CostFunc<heap_elem> {
public:
    virtual double cost(const heap_elem& elem) const;

};
double costf2::cost(const heap_elem& elem) const {
    return -elem.eval.ub();
}


// optimization structure gives by the master B&B to the slave, contains slave B&B initialization value and objective function.
class optiw {
public:

    double lower_ub;
    double preclw;
    Function *f;
    vector< bernfunc > bernwset;
    vector< bernfunc > bernwfree;
    vector< bernfunc> berncst;
    bool entire;
    SetIntervalReg tree;
    optiw(double lower_ub, double preclw,Function* f,vector< bernfunc> frac_wset,vector< bernfunc > frac_wfree,vector< bernfunc> frac_cst,bool entire,SetIntervalReg tree);

};
optiw::optiw(double lower_ub, double preclw,Function* f,vector< bernfunc> frac_wset,vector< bernfunc> frac_wfree,vector< bernfunc> frac_cst,bool entire,SetIntervalReg tree): lower_ub(lower_ub),preclw(preclw),f(f),bernwset(frac_wset),bernwfree(frac_wfree),berncst(frac_cst),entire(entire),tree(tree)
{}


Vector randpt(IntervalVector box, int size) {
    Vector pt(size);
    for(unsigned i=0;i<size;i++) {
        int rnd = std::rand();
        pt[i] = double(rnd)/double(RAND_MAX)*box[i].diam()+box[i].lb();
    }
    return pt;
}


// create pair of numer denom polynomial used for bernstein eval from ibex function
// arguments: f: ibex function, aff_var: variables in which f will be considered polynomial, nbvar_aff: number of aff_var,
//var: other variables in which f will not be considered polynomial, nbvar: number of var.
pair<polynomial*,polynomial*> get_pol(Function *f,string * aff_var,unsigned nbvar_aff,string * var,unsigned nbvar) {
    symtab tab;
    for(unsigned i=0;i<nbvar_aff;i++)
        tab[aff_var[i]] = symbol(aff_var[i]);
    for(unsigned i=0;i<nbvar;i++)
        tab[var[i]] = symbol(var[i]);
    stringstream expstring;
    expstring<<f->expr();
    cout<<"expstring: "<<expstring.str()<<endl;
    GiNaC::parser reader(tab);
    ex twz1ex = reader(expstring);
    ex twz1exnum = numer(twz1ex);
    stringstream ss;
    ss<<twz1exnum;
    cout<<"numerator expression: "<<ss.str()<<endl;
    polynomial *num = new polynomial(ss.str(),aff_var,nbvar_aff,var,nbvar);
    cout<<"numerator created"<<endl;
    ex twz1exden = denom(twz1ex);
    cout<<"get denominator"<<endl;
    stringstream ssden;
    ssden<<twz1exden;
    cout<<"denominator expression: "<<ssden.str()<<endl;
    polynomial *den = new polynomial(ssden.str(),aff_var,nbvar_aff,var,nbvar);
    cout<<"denominator created"<<endl;
    pair<polynomial*,polynomial*> frac;
    frac.first=num;
    frac.second=den;
    return frac;



}

//slave B&B, compute fmax and wmax, see MACIS
// argument: optiw: see class definition, list_elem: see class definition, midp: true if mind point on controller parameter require, false else.
double best_upper_bound_forall_w(optiw * str,list_elem * elem,bool midp,Vector midpt) {
    cout<<"in func"<<endl;
    costf2 b;
    Heap<heap_elem> heap(b); // focus on max ub
    stack<heap_elem*> save_stack; // save elements


    //list.push(*str.Inilw);
    IntervalVector res(1,Interval::EMPTY_SET),resmid(1,Interval::EMPTY_SET); // eval result
    heap_elem *elemtmp; // current element
    double upper_ub = NEG_INFINITY; // upper bound on max
    double lower_ub = elem->fmax[0].lb(); // lower bound on max, initialized with result inherited see MACIS
//    cout<<"preclw: "<<str->preclw<<endl;
    //cout<<"upper_ub initially: "<<upper_ub<<"fmax initially: "<<elem->fmax<<endl;
    //cout<<"number of w to check: "<<elem->w_interest.size()<<endl;
//    int nbit = elem->w_entire.size();
//    SetNodeReg* nodetmp;
//    IntervalVector boxtmp(1);

//    cout<<"get "<<node_vect.size()<<"nodes and: "<<node_box.size()<<"boxes"<<endl;

    // keep only wmax from the tree leaves and create heap elem from them
    cout<<"init stack"<<endl;

    IntervalVector box(elem->box.size()+1);
    for(unsigned i=0;i<elem->box.size();i++){
        if(midp)
            box[i] = midpt[i];
        else
            box[i] = elem->box[i];
    }
    for (unsigned i=0;i<elem->wint.size();i++) {
        IntervalVector nodebox = str->tree.findNodeBox(elem->wint.at(i));
        box[elem->box.size()] = nodebox[0];
        res = str->f->eval_vector(box);
        heap.push(new heap_elem(elem->wint.at(i),nodebox,res[0],true));
    }
    cout<<"stack ok"<<endl;

//    cout<<"**************************** box: "<<box<<endl;

    // B&B algo slave

    unsigned nbit(0),nbitmax;
    if(midp)
        nbitmax=1000000000;
    else
        nbitmax=10;
    cout<<"while start"<<endl;
    while(!heap.empty()) {

//    cout<<"nb box to treat: "<<count<<endl;
//    for(unsigned i=0;i<nbitmax;i++){
//        if(heap.empty()) {

//        }
//        cout<<"lower ub: "<<lower_ub<<endl;
        elemtmp = heap.pop();
//        cout<<"box: "<<elemtmp->box<<" nb it: "<<nbit<<endl;
//        cout<<"box: "<<elemtmp->box<<endl;
        box[elem->box.size()] = elemtmp->box[0];

//        box[elem->box.size()] = elemtmp->box.mid()[0];


//            cout<<"first res mid: "<<resmid<<" for box: "<<box<<endl;
        if(midp){
            //********box eval********
            IntervalVector cres = str->f->eval_vector(box);
//            cout<<"res for box "<<box<<" : "<<res<<endl;
//            IntervalVector round_box(elem->box.size()+1);
//            round_box[0] = ibex::pow(Interval(10),elemtmp->box[0]);
//            for(unsigned i=1;i<=elem->box.size();i++) {
//                round_box[i] = box[i-1];
//            }
//            res[0] =str->bernwfree.at(0).eval_bernstein(round_box)[0]*str->berncst.at(0).eval_bernstein(IntervalVector(1,round_box[0]))[0];
//            for(unsigned i=1;i<str->bernwfree.size();i++)
//                res[0] = max(res[0],str->bernwfree.at(i).eval_bernstein(round_box)[0]*str->berncst.at(i).eval_bernstein(IntervalVector(1,round_box[0]))[0]);
//////            cout<<"bernstein res for box "<<round_box<<" : "<<res<<endl;
////            res.operator &=(str->f->eval_vector(box));

//            if(res.is_disjoint(cres)) {
//                cout<<"ERROR: disjoint classic and bern eval"<<endl;
//                cout<<"cres: "<<cres<<" , res: "<<res<<endl;
//            }
//            else
//                res.operator &=(cres);
            res = cres;

            //****** mid w eval *******
            box[elem->box.size()] = elemtmp->box.mid()[0];
            resmid = str->f->eval_vector(box);
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
            //****** mid w eval *******
            box[elem->box.size()] = elemtmp->box.mid()[0];
            IntervalVector cresmid = str->f->eval_vector(box);
//            cout<<"resmid for box "<<box<<" : "<<cresmid[0]<<endl;
//            box[elem->box.size()] = ibex::pow(Interval(10),box[elem->box.size()]);
//            resmid[0] =str->bernwset.at(0).eval_bernstein(box)[0]*str->berncst.at(0).eval_bernstein(IntervalVector(1,box[elem->box.size()]))[0];
//            for(unsigned i=1;i<str->bernwset.size();i++)
//                resmid[0] = max(resmid[0],str->bernwset.at(i).eval_bernstein(box)[0]*str->berncst.at(i).eval_bernstein(IntervalVector(1,box[elem->box.size()]))[0]);
////           cout<<"bern resmid for box "<<box<<": "<<resmid<<endl;

//            if(resmid.is_disjoint(cresmid)) {
//                cout<<"ERROR: disjoint classic and bern eval"<<endl;
//                cout<<"cresmid: "<<cresmid<<" , resmid: "<<resmid<<endl;
//            }
//            else
//                resmid.operator &=(cresmid);
            resmid = cresmid;

        }

        if(resmid[0].lb()>str->lower_ub) {
            if(!midp){
                elem->fmax = IntervalVector(1,Interval::EMPTY_SET); // box K does not minimize the maximum
                delete elemtmp;
                heap.flush();
                return 0;
            }
            else{
                heap.flush();
                delete elemtmp;
//                cout<<"midpoint failure"<<endl;
                return POS_INFINITY;
            }
        }

        if(res[0].ub()<lower_ub) { // this box w does not contains the maximum for the box K
            delete elemtmp;
            continue;
        }

        if(!midp) {
//            cout<<"freq: "<<elemtmp->box[0]<<" res: "<<res[0]<<" resmid: "<<resmid[0]<<endl;
        }


        //elem->fmax = res[0].lb()>elem->fmax[0].lb() ? IntervalVector(1,Interval(res[0].lb(),elem->fmax[0].ub())): elem->fmax;  // increase upper_lb
        elemtmp->eval = res[0];

        lower_ub = resmid[0].lb()>lower_ub? resmid[0].lb(): lower_ub;
        if(elemtmp->box[elemtmp->node->var].diam() > str->preclw && nbit<nbitmax) { // bissect and push
            nbit ++;
            if(elemtmp->node->is_leaf())
                elemtmp->node->cut(elemtmp->box);
            heap.push(new heap_elem(elemtmp->node->left,elemtmp->node->left_box(elemtmp->box),elemtmp->eval,false));
            heap.push(new heap_elem(elemtmp->node->right,elemtmp->node->right_box(elemtmp->box),elemtmp->eval,false));
            delete elemtmp;
        }
        else { // min prec on box reached, keep heap elem
//            cout<<"max :"<< elemtmp->eval<<" for freq: "<<elemtmp->box<<endl;
            save_stack.push(elemtmp);
//            cout<<"update upper_ub to: "<<upper_ub<<endl;
        }
    }
    cout<<"while end"<<endl;
    while(!save_stack.empty()) {
        heap.push(save_stack.top());
        save_stack.pop();
    }
    upper_ub = heap.top()->eval.ub();
    while(!heap.empty()) {
        elemtmp = heap.pop();
        elem->wint.push_back(elemtmp->node);
    }

//            cout<<volw/4*100<<" % of w kept"<<endl;
    if(!midp)
        elem->fmax = IntervalVector(1,Interval(lower_ub,upper_ub));

    return upper_ub; // return ub on max, used as result for controller parameter midpoint eval
}


int main() {
    std::srand(std::time(0));

    IntervalVector Iniprob(2,Interval(-10,10));
    IntervalVector Iniprob2(1,Interval(-0.1,0.1));
    // precision:
    double prec(0.01);
    double preclwmin = 1.e-5;// dynamic initialization in loop
    double stop_criterion(0.1); // stop if distance between uplo and loup lower than stop_criterion

    // init boxes

    IntervalVector IniboxK(3,Interval(-10,10));
    IntervalVector Inilw(1,Interval(-3,3));


    //function definition
    Variable kp,ki,kd,w,u;
    Function Twz1("kp","ki","kd","u","(sqrt(((-0.1464*exp(ln(10)*4*u)+33.75*exp(ln(10)*6*u))^2+(-33.7524*exp(ln(10)*5*u)+0.14399999999999999999*exp(ln(10)*3*u))^2)*1/(((0.00366)*exp(ln(10)*2*u)*kp+0.0036*kd*exp(ln(10)*2*u)-0.1464*exp(ln(10)*4*u)+33.75*exp(ln(10)*6*u)-0.0036*ki-exp(ln(10)*4*u)*kp-kd*exp(ln(10)*4*u)+1.00006*exp(ln(10)*2*u)*ki)^2+(-0.00366*exp(ln(10)*u)*ki+(5.9999999999999999995*10^(-5))*kd*exp(ln(10)*3*u)-33.7524*exp(ln(10)*5*u)+1.00006*exp(ln(10)*3*u)*kp+0.14399999999999999999*exp(ln(10)*3*u)+exp(ln(10)*3*u)*ki-0.0036*exp(ln(10)*u)*kp)^2)) + sqrt((((0.036)*exp(ln(10)*u)-10.0006*exp(ln(10)*3*u))^2+(10*exp(ln(10)*4*u)-0.0366*exp(ln(10)*2*u))^2)*1/(((0.00366)*exp(ln(10)*2*u)*kp+0.0036*kd*exp(ln(10)*2*u)-0.1464*exp(ln(10)*4*u)+33.75*exp(ln(10)*6*u)-0.0036*ki-exp(ln(10)*4*u)*kp-kd*exp(ln(10)*4*u)+1.00006*exp(ln(10)*2*u)*ki)^2+(-0.00366*exp(ln(10)*u)*ki+(5.9999999999999999995*10^(-5))*kd*exp(ln(10)*3*u)-33.7524*exp(ln(10)*5*u)+1.00006*exp(ln(10)*3*u)*kp+0.14399999999999999999*exp(ln(10)*3*u)+exp(ln(10)*3*u)*ki-0.0036*exp(ln(10)*u)*kp)^2)))*sqrt(((-0.01248+exp(ln(10)*2*u))^2+0.024964*exp(ln(10)*2*u))*1/((0.001248-3.162*exp(ln(10)*2*u))^2+0.0078960996*exp(ln(10)*2*u)))");
    cout<<"func1 ok"<<endl;
    Function Twz2("kp","ki","kd","u","(sqrt(((-33.7524*exp(ln(10)*4*u)*ki-0.1464*exp(ln(10)*4*u)*kp+33.75*kd*exp(ln(10)*6*u)-0.14399999999999999999*kd*exp(ln(10)*4*u)+0.14399999999999999999*exp(ln(10)*2*u)*ki+33.75*exp(ln(10)*6*u)*kp)^2+((0.14399999999999999999)*exp(ln(10)*3*u)*kp-33.75*exp(ln(10)*5*u)*ki-33.7524*exp(ln(10)*5*u)*kp+0.1464*exp(ln(10)*3*u)*ki-0.0023999999999999999998*kd*exp(ln(10)*5*u))^2)*1/(((0.00366)*exp(ln(10)*2*u)*kp+0.0036*kd*exp(ln(10)*2*u)-0.1464*exp(ln(10)*4*u)+33.75*exp(ln(10)*6*u)-0.0036*ki-exp(ln(10)*4*u)*kp-kd*exp(ln(10)*4*u)+1.00006*exp(ln(10)*2*u)*ki)^2+(-0.00366*exp(ln(10)*u)*ki+(5.9999999999999999995*10^(-5))*kd*exp(ln(10)*3*u)-33.7524*exp(ln(10)*5*u)+1.00006*exp(ln(10)*3*u)*kp+0.14399999999999999999*exp(ln(10)*3*u)+exp(ln(10)*3*u)*ki-0.0036*exp(ln(10)*u)*kp)^2))+sqrt(((-0.0366*exp(ln(10)*2*u)*kp-0.036*kd*exp(ln(10)*2*u)+0.036*ki+10*exp(ln(10)*4*u)*kp+10*kd*exp(ln(10)*4*u)-10.0006*exp(ln(10)*2*u)*ki)^2+((0.0366)*exp(ln(10)*u)*ki-(5.9999999999999999995*10^(-4))*kd*exp(ln(10)*3*u)-10.0006*exp(ln(10)*3*u)*kp-10*exp(ln(10)*3*u)*ki+0.036*exp(ln(10)*u)*kp)^2)*1/(((0.00366)*exp(ln(10)*2*u)*kp+0.0036*kd*exp(ln(10)*2*u)-0.1464*exp(ln(10)*4*u)+33.75*exp(ln(10)*6*u)-0.0036*ki-exp(ln(10)*4*u)*kp-kd*exp(ln(10)*4*u)+1.00006*exp(ln(10)*2*u)*ki)^2+(-0.00366*exp(ln(10)*u)*ki+(5.9999999999999999995*10^(-5))*kd*exp(ln(10)*3*u)-33.7524*exp(ln(10)*5*u)+1.00006*exp(ln(10)*3*u)*kp+0.14399999999999999999*exp(ln(10)*3*u)+exp(ln(10)*3*u)*ki-0.0036*exp(ln(10)*u)*kp)^2)))*sqrt(1/((7896.0996)*exp(ln(10)*2*u)+(-3948+exp(ln(10)*2*u))^2)*((78.960996)*exp(ln(10)*2*u)+(-39.48+exp(ln(10)*2*u))^2))");
    //    Function Twz1("kp","ki","kd","u", "1/((-241.00000000000000001*exp(ln(10)*4*u)+kd*exp(ln(10)*2*u)+101*kp*exp(ln(10)*2*u)-ki+100*exp(ln(10)*2*u)*ki+102.4*exp(ln(10)*2*u))^2+(100*exp(ln(10)*5*u)+101*exp(ln(10)*u)*ki+kp*exp(ln(10)*u)-100*kd*exp(ln(10)*3*u)-100*kp*exp(ln(10)*3*u)+exp(ln(10)*u)-242.40000000000000001*exp(ln(10)*3*u))^2)*((exp(ln(10)*5*u)+100*exp(ln(10)*u)-242.4*exp(ln(10)*3*u))^2+(-102.4*exp(ln(10)*4*u)+241.0*exp(ln(10)*2*u))^2)");
    //    Function Twz2("kp","ki","kd","u","1/(((12.4)*exp(ln(10)*4*u)-10*kd*exp(ln(10)*2*u)-11*kp*exp(ln(10)*2*u)+10*ki-exp(ln(10)*2*u)*ki-25.0*exp(ln(10)*2*u))^2+(exp(ln(10)*5*u)+11*exp(ln(10)*u)*ki+10*kp*exp(ln(10)*u)-kd*exp(ln(10)*3*u)-kp*exp(ln(10)*3*u)+10*exp(ln(10)*u)-26.4*exp(ln(10)*3*u))^2)*(((12.4)*exp(ln(10)*u)*ki+kp*exp(ln(10)*u)-11.4*kd*exp(ln(10)*3*u)-26.4*kp*exp(ln(10)*3*u)+10*kd*exp(ln(10)*5*u)-25.0*exp(ln(10)*3*u)*ki+10*kp*exp(ln(10)*5*u))^2+(kd*exp(ln(10)*2*u)+12.4*kp*exp(ln(10)*2*u)-10*exp(ln(10)*4*u)*ki-25.0*kp*exp(ln(10)*4*u)-15.0*kd*exp(ln(10)*4*u)-ki+26.4*exp(ln(10)*2*u)*ki)^2)");
    cout<<"u func ok"<<endl;
    Function Tw1z1w("kp","ki","kd","w","1/((-(0.0036)*w*kp+(0.14399999999999999999)*w^3+w^3*ki-(33.7524)*w^5+(1.00006)*w^3*kp-(0.00366)*w*ki+(5.9999999999999999995e-5)*kd*w^3)^2+(kd*w^4-(1.00006)*w^2*ki+w^4*kp+(0.0036)*ki+(0.1464)*w^4-(33.75)*w^6-(0.00366)*w^2*kp-(0.0036)*kd*w^2)^2)*((-(0.1464)*w^4+(33.75)*w^6)^2+((0.14399999999999999999)*w^3-(33.7524)*w^5)^2)");
    //    Function Twz1w("kp","ki","kd","w","1/((101*w*ki+kp*w+100*w^5-100*kp*w^3-100*kd*w^3+w-(242.40000000000000001)*w^3)^2+((241.00000000000000001)*w^4-101*kp*w^2-kd*w^2+ki-100*w^2*ki-(102.4)*w^2)^2)*((w^5+100*w-(242.4)*w^3)^2+(-(102.4)*w^4+(241.0)*w^2)^2)");
    cout<<"11 ok"<<endl;
    Function Tw2z1w("kp","ki","kd","w","((-(0.0366)*w^2+10*w^4)^2+((0.036)*w-(10.0006)*w^3)^2)*1/((-(0.0036)*w*kp+(0.14399999999999999999)*w^3+w^3*ki-(33.7524)*w^5+(1.00006)*w^3*kp-(0.00366)*w*ki+(5.9999999999999999995e-5)*kd*w^3)^2+(kd*w^4-(1.00006)*w^2*ki+w^4*kp+(0.0036)*ki+(0.1464)*w^4-(33.75)*w^6-(0.00366)*w^2*kp-(0.0036)*kd*w^2)^2)");
    cout<<"out1 ok"<<endl;
    Function Tw1z2w("kp","ki","kd","w","(((33.75)*kd*w^6+(33.75)*w^6*kp+(0.14399999999999999999)*w^2*ki-(0.14399999999999999999)*kd*w^4-(0.1464)*w^4*kp-(33.7524)*w^4*ki)^2+(-(0.0023999999999999999998)*kd*w^5-(33.7524)*w^5*kp+(0.1464)*w^3*ki+(0.14399999999999999999)*w^3*kp-(33.75)*w^5*ki)^2)*1/((-(0.0036)*w*kp+w^3*ki-(33.7524)*w^5+(1.00006)*w^3*kp+(0.14399999999999999999)*w^3+(5.9999999999999999995e-5)*kd*w^3-(0.00366)*w*ki)^2+((1.00006)*w^2*ki-kd*w^4+(33.75)*w^6-(0.1464)*w^4-w^4*kp-(0.0036)*ki+(0.00366)*w^2*kp+(0.0036)*kd*w^2)^2)");
    //    Function Twz2w("kp","ki","kd","w","(((12.4)*w*ki+kp*w-(26.4)*kp*w^3-(11.4)*kd*w^3+10*kp*w^5-(25.0)*w^3*ki+10*kd*w^5)^2+((12.4)*kp*w^2+kd*w^2-10*w^4*ki-(15.0)*kd*w^4-(25.0)*kp*w^4-ki+(26.4)*w^2*ki)^2)*1/((-(12.4)*w^4+11*kp*w^2+10*kd*w^2-10*ki+w^2*ki+(25.0)*w^2)^2+(11*w*ki+10*kp*w+w^5-kp*w^3-kd*w^3+10*w-(26.4)*w^3)^2)");
    Function Tw2z2w("kp","ki","kd","w","1/(((0.14399999999999999999)*w^3+w^3*ki-(0.0036)*w*kp-(0.00366)*w*ki+(5.9999999999999999995e-5)*kd*w^3-(33.7524)*w^5+(1.00006)*w^3*kp)^2+(w^4*kp-(1.00006)*w^2*ki+kd*w^4+(0.0036)*ki-(0.00366)*w^2*kp-(0.0036)*kd*w^2+(0.1464)*w^4-(33.75)*w^6)^2)*((10*w^4*kp-(10.0006)*w^2*ki+10*kd*w^4+(0.036)*ki-(0.0366)*w^2*kp-(0.036)*kd*w^2)^2+(10*w^3*ki-(0.036)*w*kp-(0.0366)*w*ki+(5.9999999999999999995e-4)*kd*w^3+(10.0006)*w^3*kp)^2)");
    Function w1("kp","ki","kd","w","((-0.01248+w^2)^2+(0.024964)*w^2)*1/((0.001248-(3.162)*w^2)^2+(0.0078960996)*w^2)");
    //    Function w1("kp","ki","kd","w","1/(10000+w^2)*(1+10000*w^2)");
    Function w2("kp","ki","kd","w","1/((7896.0996)*w^2+(-3948+w^2)^2)*((78.960996)*w^2+(-39.48+w^2)^2)");
    //    Function w2("kp","ki","kd","w","(100+w^2)*1/(1+100*w^2)");


    vector<bernfunc> bernwset;
    vector<bernfunc> bernwfree;
    vector<bernfunc> berncst;

    vector< pair<polynomial *,polynomial *> > frac_wfree;
    vector< pair<polynomial *,polynomial *> > frac_wset;
    vector< pair<polynomial *,polynomial *> > frac_cst;

    string aff_var[4];aff_var[0] = "kp";aff_var[1] = "ki";aff_var[2] = "kd";
    string var[4];var[0]="w";
    frac_wset.push_back(get_pol(&Tw1z1w,aff_var,3,var,1));
    frac_wset.push_back(get_pol(&Tw2z1w,aff_var,3,var,1));
    bernwset.push_back(bernfunc(frac_wset));
    frac_wset.clear();
    frac_wset.push_back(get_pol(&Tw1z2w,aff_var,3,var,1));
    frac_wset.push_back(get_pol(&Tw2z2w,aff_var,3,var,1));
    bernwset.push_back(bernfunc(frac_wset));
    aff_var[0] = "w";
    var[0] = "kp";var[1] = "ki";var[2] = "kd";
    frac_wfree.push_back(get_pol(&Tw1z1w,aff_var,1,var,3));
    frac_wfree.push_back(get_pol(&Tw2z1w,aff_var,1,var,3));
    bernwfree.push_back(bernfunc(frac_wfree));
    frac_wfree.clear();
    frac_wfree.push_back(get_pol(&Tw1z2w,aff_var,1,var,3));
    frac_wfree.push_back(get_pol(&Tw2z2w,aff_var,1,var,3));
    bernwfree.push_back(bernfunc(frac_wfree));
    frac_cst.push_back(get_pol(&w1,aff_var,1,var,0));
    berncst.push_back(bernfunc(frac_cst));
    frac_cst.clear();
    frac_cst.push_back(get_pol(&w2,aff_var,1,var,0));
    berncst.push_back(bernfunc(frac_cst));


    //    Interval Kp(5.31125, 5.31125);
    //    Interval Ki(1.475, 1.475);
    //    Interval Kd(-7.975, -7.975);


    //    pair<polynomial *,polynomial *>  W2 = get_pol(&w2,aff_var,1,var,0);
    //    IntervalVector tbox(4);tbox[0]=Kp;tbox[1]=Ki;tbox[2]=Kd;tbox[3]=ibex::pow(Interval(10),Interval(1.96875, 1.98438));
    //    cout<<"classic evaluation: "<<Twz2w.eval_vector(tbox)<<endl;
    //    tbox[1]=Kp;tbox[2]=Ki;tbox[3]=Kd;tbox[0]=ibex::pow(Interval(10),Interval(1.96875, 1.98438));
    //    cout<<"polynum eval: "<<frac_wfree.at(1).first->eval_bernstein(tbox)[0]<<endl;
    //    cout<<"polyden eval: "<<frac_wfree.at(1).second->eval_bernstein(tbox)[0]<<endl;
    //    cout<<"bernstein evaluation: "<<(frac_wfree.at(1).first->eval_bernstein(tbox)[0]/frac_wfree.at(1).second->eval_bernstein(tbox)[0])<<endl;
    //    cout<<"W2 eval: "<<(W2.first->eval_bernstein(IntervalVector(1,tbox[0]))[0]/W2.second->eval_bernstein(IntervalVector(1,tbox[0]))[0])<<endl;
    ////    cout<<"bernstein evaluation frac: "<<(frac.first->eval_bernstein(tbox)[0]/frac.second->eval_bernstein(tbox)[0])<<endl;
    //    return 0;
     cout<<"objective function ok"<<endl;
    Function rS1("kp","ki","kd","0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki");
    Function rS2("kp","ki","kd","1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*((0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(33.7524)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))");
    Function rS3("kp","ki","kd","(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*1/((0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(33.7524)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*(1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*((0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(33.7524)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(((0.0036)*kp+(0.00366)*ki)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*1/(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))");
    Function rS4("kp","ki","kd","-1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*((0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(33.7524)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*((0.0036)*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*((0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(33.7524)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*ki-(((0.0036)*kp+(0.00366)*ki)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*1/((0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(33.7524)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*((0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(33.7524)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(((0.0036)*kp+(0.00366)*ki)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*1/(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)))*1/(1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*((0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(33.7524)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(((0.0036)*kp+(0.00366)*ki)*(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*1/(((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)-(0.1464+kp+kd)*((0.0036)*kp+(0.00366)*ki))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))");
    Function rS5("kp","ki","kd","(0.0036)*ki");
    Function rGS1("kp","ki","kd","0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki");
    Function rGS2("kp","ki","kd","((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)");
    Function rGS3("kp","ki","kd","(1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)");
    Function rGS4("kp","ki","kd","-1/(1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*(1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))*((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+(0.0036)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*ki*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd)))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)");
    Function rGS5("kp","ki","kd","(0.0036)*ki");
    Function rT1("kp","ki","kd","0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki");
    Function rT2("kp","ki","kd","((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)");
    Function rT3("kp","ki","kd","1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))");
    Function rT4("kp","ki","kd","-((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*(((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))+(0.0036)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)*ki)*1/(((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)-(33.7524)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*1/(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki)+((0.0036)*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*ki-((0.0036)*kp+(0.00366)*ki)*((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*1/(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki))*1/((0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*((0.00366)*kp+(0.0036)*kd+(1.00006)*ki)-((0.0036)*kp+(0.00366)*ki)*(0.1464+kp+kd))*(0.14399999999999999999+(1.00006)*kp+(5.9999999999999999995e-5)*kd+ki)*(0.002410239271874000076+(1.11103210438367897874e-5)*kp+(0.99994000426636328084)*kd-(0.999928893945319444)*ki))");
    Function rT5("kp","ki","kd", "(0.0036)*ki");


    //    Function rS1("kp","ki","kd","60*kd - 25*ki + 35*kp + 119");
    //    Function rS2("kp","ki","kd", "(5*(60*kd - 50*ki + 154*kp + 60*kd*ki + 60*kd*kp + 10*ki*kp - 25*ki^2 + 35*kp^2 + 119))/(60*kd - 25*ki + 35*kp + 119)");
    //    Function rS3("kp","ki","kd", "5*ki*(5*kd + 5*kp + 12) - (125*ki^2)/12");
    //    cout<<frac_wfree.at(1).first->eval_bernstein(IntervalVector(4,Interval(1)));
    //    cout<<frac_wfree.at(1).second->eval_bernstein(IntervalVector(4,Interval(1)));
    //    cout<<"eval ok"<<endl;


    vector<Function*> routh_func;
    vector<Ctc*> array_ctc;
    vector<Sep*> sep_array;
    NumConstraint *c1= new NumConstraint(kp,ki,kd,rS1(kp,ki,kd)>=0);
    array_ctc.push_back(new CtcFwdBwd(*c1));
    routh_func.push_back(&rS1);
    sep_array.push_back(new SepFwdBwd(rS1,GEQ));
    NumConstraint *c2= new NumConstraint(kp,ki,kd,rS2(kp,ki,kd)>=0);
    array_ctc.push_back(new CtcFwdBwd(*c2));
    routh_func.push_back(&rS2);
    sep_array.push_back(new SepFwdBwd(rS2,GEQ));
    NumConstraint *c3= new NumConstraint(kp,ki,kd,rS3(kp,ki,kd)>=0);
    array_ctc.push_back(new CtcFwdBwd(*c3));
    routh_func.push_back(&rS3);
    sep_array.push_back(new SepFwdBwd(rS3,GEQ));
    NumConstraint *c4= new NumConstraint(kp,ki,kd,rS4(kp,ki,kd)>=0);
    array_ctc.push_back(new CtcFwdBwd(*c4));
    routh_func.push_back(&rS4);
    sep_array.push_back(new SepFwdBwd(rS4,GEQ));
    NumConstraint *c5= new NumConstraint(kp,ki,kd,rS5(kp,ki,kd)>=0);
    array_ctc.push_back(new CtcFwdBwd(*c5));
    sep_array.push_back(new SepFwdBwd(rS5,GEQ));
    routh_func.push_back(&rS5);
    //    NumConstraint *c6= new NumConstraint(kp,ki,kd,rT1(kp,ki,kd)>=0);
    //    array_ctc.push_back(new CtcFwdBwd(*c6));
    //    NumConstraint *c7= new NumConstraint(kp,ki,kd,rT2(kp,ki,kd)>=0);
    //    array_ctc.push_back(new CtcFwdBwd(*c7));
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
    CtcCompo ctc_routh(array_ctc);
    //    SepUnion routh_sep(sep_array);
    //    SetIntervalReg stabpav(IniboxK,0.1);
    //    stabpav.contract(routh_sep);
    //    visit_leaves(stabpav);

    //    return 0;

    //    IntervalVector pointstab(3);pointstab[0] = Interval(7.916036169838760);pointstab[1] = Interval(0.074448101237470);pointstab[2] = Interval(1.509569215261232);
    //    cout<<"stability of point "<<pointstab<<" : "<<endl
    //          << "  rS1: "<<rS1.eval_vector(pointstab)
    //           << "  rS2: "<<rS2.eval_vector(pointstab)
    //               << "  rS3: "<<rS3.eval_vector(pointstab)
    //                   << "  rS4: "<<rS4.eval_vector(pointstab)
    //                       << "  rS5: "<<rS5.eval_vector(pointstab)<<endl;
    //    return 0;

    //    Function F3("kp","ki","kd","w","(w^4*(25*kd^2 + 50.0*kd*kp + 25*kp^2) + w^2*(25*ki^2 - 50.0*kd*ki + 25*kp^2) + 25*ki^2)/(25.0*kd^2*w^4 - 50.0*kd*ki*w^2 + 50.0*kd*kp*w^4 - 50.0*kd*w^6 + 120.0*kd*w^4 + 25.0*ki^2*w^2 + 25.0*ki^2 - 70.0*ki*w^4 - 70.0*ki*w^2 + 25.0*kp^2*w^4 + 25.0*kp^2*w^2 - 50.0*kp*w^6 + 50.0*kp*w^2 + 25.0*w^8 + 24.0*w^6 + 24.0*w^4 + 25.0*w^2)");
    //    Function F3u("kp","ki","kd","u","(25*(kd^2*exp(9.210340372*u) + ki^2*exp(4.605170186*u) + kp^2*exp(4.605170186*u) + kp^2*exp(9.210340372*u) + ki^2 - 2.0*kd*ki*exp(4.605170186*u) + 2.0*kd*kp*exp(9.210340372*u)))/(25.0*exp(4.605170186*u) + 24.0*exp(9.210340372*u) + 25.0*exp(18.42068074*u) + 24.0*exp(13.81551056*u) + 120.0*kd*exp(9.210340372*u) - 70.0*ki*exp(4.605170186*u) - 70.0*ki*exp(9.210340372*u) + 50.0*kp*exp(4.605170186*u) - 50.0*kd*exp(13.81551056*u) - 50.0*kp*exp(13.81551056*u) + 25.0*kd^2*exp(9.210340372*u) + 25.0*ki^2*exp(4.605170186*u) + 25.0*kp^2*exp(4.605170186*u) + 25.0*kp^2*exp(9.210340372*u) + 25.0*ki^2 - 50.0*kd*ki*exp(4.605170186*u) + 50.0*kd*kp*exp(9.210340372*u))");

    //    Function f(ki,kp,kd,w,ki+kp-kd+w);

    //    Function m1(ki,kp,ibex::max(ki,kp));
    //    IntervalVector test(2);
    //    test[0] = Interval(1,2);
    //    test[1] = Interval(3,4);
    //    cout<<"max of "<<test<<" = "<<m1.eval_vector(test);
    //    test[0] = Interval(1,3.5);
    //    test[1] = Interval(3,4);
    //    cout<<"max of "<<test<<" = "<<m1.eval_vector(test);

    //    return 0;

    Function Max12(kp,ki,kd,u,ibex::max(Twz1(kp,ki,kd,u),Twz2(kp,ki,kd,u)));
    Function Sum(kp,ki,kd,u,ibex::max(Twz1(kp,ki,kd,u),Twz2(kp,ki,kd,u)));



    //    costf2 b;
    //    Heap<heap_elem> heap(b);
    double freqeps=0.0001;
    //    Interval freq(-3,-3+freqeps);
    //    IntervalVector box(4),bbox(4);
    //    Interval bernres(-10000),cres;
    //    box[0] = Interval(-7.17513, -6.62382);box[1]=Interval( 6.2275, 6.56463);box[2]=Interval(6.22509, 6.97427);
    //    bbox[1] = Interval(-7.17513, -6.62382);bbox[2]=Interval( 6.2275, 6.56463);bbox[3]=Interval(6.22509, 6.97427);
    ////    box[0] = Interval(7.916);box[1]=Interval(0.0744);box[2]=Interval(1.5096);
    ////    bbox[1] = Interval(7.916);bbox[2]=Interval(0.0744);bbox[3]=Interval(1.5096);
    //    IntervalVector berndraw(2),cdraw(2);
    //    vibes::beginDrawing();
    //    vibes::newFigure("Twz2 go");
    //    while(freq.ub()<3) {
    //        cdraw[0] = freq;
    //        berndraw[0] = freq;
    //        box[3]=freq;
    //        bbox[0]=ibex::pow(Interval(10),freq);
    //        berndraw[1] = bernwfree.at(1).eval_bernstein(bbox)[0]*berncst.at(1).eval_bernstein(IntervalVector(1,bbox[0]))[0];

    //        cdraw[1] = Twz2.eval_vector(box)[0];
    //        vibes::drawBox(cdraw,"blue[blue]");
    //        vibes::drawBox(berndraw,"red[red]");
    ////        heap.push(new heap_elem(NULL,bbox,bernres));
    //        freq+=Interval(freqeps);
    //    }
    ////    cout<<"max of bernres: "<<heap.top()->eval<<endl;
    ////    heap.flush();
    //    return 0;


    //    Function Max123(kp,ki,kd,u,ibex::max(F3u(kp,ki,kd,u),Max12(kp,ki,kd,u)));

    //other variables for algo

    //    Function Prob("kp","ki","w","((kp+w-2)^6+0.2)*ln(1+(kp+w)^2)+((ki+w-2)^6+0.2)*ln(1+(ki+w)^2)"  );


    //    double lower_ub(POS_INFINITY);
    double lower_ub(13);
    double upper_lb(0);

    LargestFirst lf;
    //    SetIntervalReg tree(Inilw,preclwmin,false);
    //    SetIntervalReg tree(Iniprob2,preclwmin,false);
    //    tree.root->status = __IBEX_UNK__;

    //    costf1 maxlb;


    SetIntervalReg iset(Inilw,preclwmin);
    vector<SetNodeReg*> winterest;
    winterest.push_back(iset.root);

    cost_uplo cu;
    Heap<list_elem> list(cu);

    list_elem *elem = new list_elem(IniboxK,winterest,IntervalVector(1,Interval::ALL_REALS));
    //    list_elem *elem = new list_elem(Iniprob, tree,IntervalVector(1,Interval::ALL_REALS));
    list.push(elem);


    Vector respoint(3);
    //    Vector respoint(2);

    optiw str(lower_ub,preclwmin,&Max12,bernwset,bernwfree,berncst,false,iset);


    double vol_rejected(0);
    list_elem * elemtmp;
//    Vector sol(4);
//    //matlab solution
//    sol[0] = 0.3998;sol[1] = 0.4806;sol[2] = -0.1952;
//    cout<<"evaluation of matlab result: "<<Sum.eval_vector(sol)<<endl;
//    return 0;


//*********** Master B&B **************
    Timer::start();
    while(!list.empty()) {
        elemtmp = list.pop();
        upper_lb = elemtmp->fmax[0].lb();
        if(lower_ub-elemtmp->fmax[0].lb()<stop_criterion) { // stop creterion reached
            break;
        }
//        cout<<" fmax of current element: "<<elemtmp->fmax<<endl;
        if(elemtmp->fmax[0].lb()>lower_ub) { //better lower_ub may be computed since computation of the father, must check
            vol_rejected+= elemtmp->box.volume();
            cout<<"loup : "<<lower_ub<<" get for point: "<<respoint<<" uplo: "<<elemtmp->fmax.lb()<< " volume rejected: "<<vol_rejected/IniboxK.volume()*100<<endl;
            delete elemtmp;
            continue;
        }
        double ratio = (elemtmp->box[0]).diam()/(IniboxK[0]).diam()+(elemtmp->box[1]).diam()/(IniboxK[1]).diam()+(elemtmp->box[2]).diam()/(IniboxK[2]).diam();
        str.preclw = ratio/Inilw.volume()>preclwmin?ratio/(10*Inilw.volume()):preclwmin;
//        str.preclw = 0.1;
//        cout<<"precision on w: "<<str.preclw<<endl;
        str.lower_ub = lower_ub;
        //******************** Routh contraction ***********************
        IntervalVector inib = elemtmp->box;
        ctc_routh.contract(elemtmp->box);
        if(elemtmp->box != inib) {
            vol_rejected += inib.volume()-elemtmp->box.volume();
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

        best_upper_bound_forall_w(&str,elemtmp,false,Vector(3));
        cout<<"res found for box "<<elemtmp->box<<" : "<<elemtmp->fmax<<endl;
        if(elemtmp->fmax[0] == Interval::EMPTY_SET) { // maximum is guaranteed to be higher than the current lower_up for box elemtmp->box
            vol_rejected+= elemtmp->box.volume();
            cout<<"loup : "<<lower_ub<<" get for point: "<<respoint<<" uplo: "<<upper_lb<< " volume rejected: "<<vol_rejected/IniboxK.volume()*100<<endl;
            delete elemtmp;
            continue;
        }
        //********** cannot keep this part, k is not proved to respect routh criterion, useless because of midpoint *******
//        if(elemtmp->fmax[0].ub()< lower_ub) {
//            lower_ub = elemtmp->fmax[0].ub();
//            respoint = elemtmp->box.mid();
//            cout<<"New min found: "<<elemtmp->fmax<<"for box: "<<elemtmp->box<<endl;
//            cout<<"upper bound : "<<lower_ub<<" get for point: "<<respoint<<" volume rejected: "<<vol_rejected/IniboxK.volume()*100<<endl;
//        }

        //*********** Contract K with constraints f(k,wmax)< lower_ub  ***********
//        if(lower_ub<POS_INFINITY){ // run contraction with every w that contains maximum
//            IntervalVector box(elemtmp->box);
////            cout<<"Init box befor contraction: "<<elemtmp->box<<endl;
//            vector<SetNodeReg*>node;
//            vector<IntervalVector> nodebox;
////            cout<<"try to get boxes from main"<<endl;
//            elemtmp->tree.getLeaf(&node,&nodebox);
////            cout<<"get "<<node.size()<<" nodes and "<<nodebox.size()<<" boxes for contraction"<<endl;
//            while(!node.empty()){ // means that exists w such as fk(w)<lower_up has no solution =>fk(w)>lower_ub
//                if(node.back()->status == __IBEX_OUT__) {
//                    node.pop_back();
//                    nodebox.pop_back();
//                    continue;
//                }
//                NumConstraint cst(kp,ki,kd,Max12(kp,ki,kd,(nodebox.back())[0])<=lower_ub);
////                NumConstraint cst(kp,ki,kd,Prob(kp,ki,(nodebox.back())[0])<=lower_ub);
//                CtcFwdBwd ctc(cst);
//                ctc.contract(elemtmp->box);
//                if(elemtmp->box.is_empty()){
////                    cout<<"empty for w_interest"<<endl;
//                    break;
//                }
//                node.pop_back();
//                nodebox.pop_back();
//            }

////            cout<<"box after contraction: "<<elemtmp->box<<endl;
//            if(box != elemtmp->box) {
//                //cout<<"contraction usefull, initbox: "<<box<<"after contraction: "<<elemtmp->box<<endl<<elemtmp->w_interest.size()
//                  // <<" interest contraction run and "<<elemtmp->w_entire.size()<<" entire contraction run"<<endl;
//                vol_rejected +=box.volume()-elemtmp->box.volume();
//                cout<<"loup : "<<lower_ub<<" get for point: "<<respoint<<" uplo: "<<elemtmp->fmax.lb()<< " volume rejected: "<<vol_rejected/IniboxK.volume()*100<<endl;
//            }
//        }
//        if(elemtmp->box.is_empty()){
//            delete elemtmp;
//            continue;
//        }

//        *************** Mid point eval *********************
        Vector midp = elemtmp->box.mid();
        bool stab(true);
        for(unsigned j=0;j<1000;j++)
        {
            midp = randpt(elemtmp->box,elemtmp->box.size());
//            cout<<"box: "<<elemtmp->box<<" midpt: "<<midp<<endl;
            stab = true;
            for(unsigned i=0;i<routh_func.size();i++) {
                stab &= routh_func.at(i)->eval_vector(midp)[0].lb()>0;
                //            cout<<"routh "<<i<<": "<<routh_func.at(i)->eval_vector(midp)[0]<<endl;
                if(!stab){
//                    cout<<"midp fail: routh("<<i<<"): "<<routh_func.at(i)->eval_vector(midp)[0]<<endl;
                    break;
                }
            }
            if(stab)
            {
                break;
            }

        }
//        midp = sol;
//        str.preclw = 0.00001;
        // if couth criteria ok for midp
        if(stab) {// routh criterion ok for midpoint
//            Function max_mid(w,Max12(Interval(midp[0]),Interval(midp[1]),Interval(midp[2]),w));
//            Function max_mid(w,Prob(Interval(midp[0]),Interval(midp[1]),w));
//            str.f = &max_mid;
            if(str.preclw>0.001) str.preclw=0.001; // forces low precision to get tight enclosure
            double max = best_upper_bound_forall_w(&str,elemtmp,true,midp);
            cout<<" max of mid point "<<midp<<": "<<max<<endl;
//            return 0;
            //            cout<<"fmax: "<<elemtmp->fmax<<"lower_ub: "<<lower_ub<<endl;

            if(max<lower_ub) { // if better solution found
                lower_ub = max;
                respoint = midp;
                cout<<"loup : "<<lower_ub<<" get for point: "<<respoint<<" uplo: "<<upper_lb<< " volume rejected: "<<vol_rejected/IniboxK.volume()*100<<endl;
            }
        }

        if(elemtmp->box.max_diam()>prec) { // bissect and push
            pair<IntervalVector,IntervalVector> boxes = lf.bisect(elemtmp->box);

            list.push(new list_elem(boxes.first,elemtmp->wint,elemtmp->fmax));
            list.push(new list_elem(boxes.second,elemtmp->wint,elemtmp->fmax));
            delete elemtmp;
        }
        else if (elemtmp->fmax[0].ub() == lower_ub){
            //cout<<"loup : "<<lower_ub<<" get for point: "<<respoint<<" uplo: "<<elemtmp->fmax.lb()<< " volume rejected: "<<vol_rejected/IniboxK.volume()*100<<endl;
            delete elemtmp;
        }
    }
    Timer::stop();
    cout<<"maximum: "<<sqrt(lower_ub)<<"for point: "<<respoint<<endl;
    cout<<" uplo: "<< sqrt(elemtmp->fmax[0].lb())<<endl;
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
