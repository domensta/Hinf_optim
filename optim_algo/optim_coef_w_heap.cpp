
#include "/home/bumquist/Tools/ibex/OUT/include/ibex/ibex.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "src/ibex_SetIntervalReg.h"
#include "src/polynomial.h"


/*
 * #ALGO4
 *
 * Differences with #ALGO2:
 * Instead of considering fixed sized slices of w, we run sivia on w for optim, with a precision that depends on the size of boxe K.
*/


using namespace ibex;
using namespace std;



class list_elem{ // element of the master B&B algo
public:
    IntervalVector box; // box of controller parameters
    IntervalVector fmax; // enclosure of fmax(box), see MACIS presentation for detail
    SetIntervalReg tree; // UKN valued leaves represents interval of frequencies that may contains maximum, correspond to wmax, see MACIS for detail
    list_elem(IntervalVector box,SetIntervalReg tree,IntervalVector fmax); // constructor

};
list_elem::list_elem(IntervalVector box,SetIntervalReg tree,IntervalVector fmax):box(box),
    tree(tree),fmax(fmax)
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

    double lower_ub;  // uplo of master
    double preclw; // minimum size of box on w
    Function *f; // ibex objectif function, defined as max(|Twz1|,...|Twzn|)
    vector< pair<polynomial *,polynomial *> > frac_wset; // list of numer denom pair of |Twzi|, use for bernstein eval of w midpoint
    vector< pair<polynomial *,polynomial *> > frac_wfree; // list of numer denom pair of |Twzi|, use for bernstein eval of controller parameters midpoint
    vector< pair<polynomial *,polynomial *> > frac_cst; // list of numer denom pair of constraints |Wi|
    bool entire; // maximum cannot be evaluated ([ENTIRE]) if true, not used
    optiw(double lower_ub, double preclw,Function* f,vector< pair<polynomial *,polynomial *> > frac_wset,vector< pair<polynomial *,polynomial *> > frac_wfree,vector< pair<polynomial *,polynomial *> > frac_cst,bool entire);

};

optiw::optiw(double lower_ub, double preclw,Function* f,vector< pair<polynomial *,polynomial *> > frac_wset,vector< pair<polynomial *,polynomial *> > frac_wfree,vector< pair<polynomial *,polynomial *> > frac_cst,bool entire): lower_ub(lower_ub),preclw(preclw),f(f),frac_wset(frac_wset),frac_wfree(frac_wfree),frac_cst(frac_cst),entire(entire)
{} // constructor implementation


// element of the heap used in slave B&B
class heap_elem {
public:
    SetNodeReg * node; // node of the tree that correspond to the box
    IntervalVector box; // box of slave B&B, correspond to a frequency interval
    Interval eval; // evaluation of objectif function with box

    heap_elem(SetNodeReg * node, IntervalVector box, Interval eval);
};
heap_elem::heap_elem(SetNodeReg * node, IntervalVector box, Interval value): node(node),box(box),eval(eval)
{}


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
double best_upper_bound_forall_w(optiw * str,list_elem * elem,bool midp) {

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
    vector<SetNodeReg*> node_vect; // use to get wmax elements from the tree as node
    vector<IntervalVector> node_box; // use to get wmax elements from the tree as box
    elem->tree.getLeaf(&node_vect,&node_box); // get every leaves of the tree
//    SetNodeReg* nodetmp;
//    IntervalVector boxtmp(1);

//    cout<<"get "<<node_vect.size()<<"nodes and: "<<node_box.size()<<"boxes"<<endl;

    // keep only wmax from the tree leaves and create heap elem from them
    while (!node_vect.empty()) {
        if(node_vect.back()->status != __IBEX_OUT__) {
            heap.push(new heap_elem(node_vect.back(),node_box.back(),Interval::ALL_REALS));
        }

        node_vect.pop_back();
        node_box.pop_back();
    }
    IntervalVector box(elem->box.size()+1);
    for(unsigned i=0;i<elem->box.size();i++){
        if(midp)
            box[i] = elem->box.mid()[i];
        else
            box[i] = elem->box[i];
    }
//    cout<<"**************************** box: "<<box<<endl;

    // B&B algo slave
    while(!heap.empty()) {
        elemtmp = heap.pop();
//        cout<<"box: "<<elemtmp->box<<endl;
        box[elem->box.size()] = elemtmp->box[0]; // get box of current element
//        box[elem->box.size()] = elemtmp->box.mid()[0];


//            cout<<"first res mid: "<<resmid<<" for box: "<<box<<endl;
        if(midp){ // midpoint asked by master B&B, berstein evaluation of polynomial function in w.
            //********box eval********
         //   res = str->f->eval_vector(box);
         //   cout<<"res for box "<<box<<" : "<<res<<endl;
            IntervalVector round_box(elem->box.size()+1); // create box for bernstein evaluation
            round_box[0] = ibex::pow(Interval(10),elemtmp->box[0]); // need to convert u to w as function are polynomial in w for berstein eval
            for(unsigned i=1;i<=elem->box.size();i++) {
                round_box[i] = box[i-1];
            }
            // bernstein evaluation of |Twzi|-|Wi|, get the max wrt i
            res[0] = str->frac_wfree.at(0).first->eval_bernstein(round_box)[0]/str->frac_wfree.at(0).second->eval_bernstein(round_box)[0]
                    -str->frac_cst.at(0).first->eval_bernstein(IntervalVector(1,round_box[0]))[0]/str->frac_cst.at(0).second->eval_bernstein(IntervalVector(1,round_box[0]))[0];
            for(unsigned j=1;j<str->frac_wset.size();j++)
                res[0] = max(res[0],str->frac_wfree.at(j).first->eval_bernstein(round_box)[0]/str->frac_wfree.at(j).second->eval_bernstein(round_box)[0]
                        -str->frac_cst.at(j).first->eval_bernstein(IntervalVector(1,round_box[0]))[0]/str->frac_cst.at(j).second->eval_bernstein(IntervalVector(1,round_box[0]))[0]);
            cout<<"bernstein res for box "<<round_box<<" : "<<res<<endl;
            //****** mid w eval *******
            box[elem->box.size()] = elemtmp->box.mid()[0];
            resmid = str->f->eval_vector(box); // w midpoint eval to increase lb of max

        }
        else // if evaluation for a box of controller parameter required
        {
            //********box eval********
            res = str->f->eval_vector(box); // ibex evaluation, to heavy for bernstein eval
            //****** mid w eval *******
            box[elem->box.size()] = elemtmp->box.mid()[0]; // get mid of w box
//            resmid = str->f->eval_vector(box);
//            cout<<"resmid for box "<<box<<" : "<<str->f->eval_vector(box)<<endl;
            box[elem->box.size()] = ibex::pow(Interval(10),box[elem->box.size()]); // convert u to w for berstein eval
            // berstein evaluation of objective function, no more polynomial in w (because fixed), bernstein evaluation possible
            resmid[0] = str->frac_wset.at(0).first->eval_bernstein(box)[0]/str->frac_wset.at(0).second->eval_bernstein(box)[0]
                    -str->frac_cst.at(0).first->eval_bernstein(IntervalVector(1,box[elem->box.size()]))[0]/str->frac_cst.at(0).second->eval_bernstein(IntervalVector(1,box[elem->box.size()]))[0];
            for(unsigned i=1;i<str->frac_wset.size();i++)
                resmid[0] = max(resmid[0],str->frac_wset.at(i).first->eval_bernstein(box)[0]/str->frac_wset.at(i).second->eval_bernstein(box)[0]
                        -str->frac_cst.at(i).first->eval_bernstein(IntervalVector(1,box[elem->box.size()]))[0]/str->frac_cst.at(i).second->eval_bernstein(IntervalVector(1,box[elem->box.size()]))[0]);
//            cout<<"bern resmid for box "<<box<<": "<<resmid<<endl;
        }

        if(resmid[0].lb()>str->lower_ub && !midp) { // carefull, modification of master B&B element, only if !midpt
            elem->fmax = IntervalVector(1,Interval::EMPTY_SET); // box K does not minimize the maximum
            elem->tree.gather();
            delete elemtmp;
            heap.flush();
            return 0;}
        if(res[0].ub()<lower_ub) { // this box w does not contains the maximum for the box K
            if(!midp)
                elemtmp->node->status = __IBEX_OUT__; // modification of wmax, only if !midpt
            delete elemtmp;
            continue;
        }


        //elem->fmax = res[0].lb()>elem->fmax[0].lb() ? IntervalVector(1,Interval(res[0].lb(),elem->fmax[0].ub())): elem->fmax;  // increase upper_lb
        elemtmp->eval = res[0]; // store result for this element, useless as heap used? may be relic from previous version
//        if(!midp)
//            cout<<"midp: "<<elemtmp->box.mid()<<", eval mid: "<<resmid<<endl;

        lower_ub = resmid[0].lb()>lower_ub? resmid[0].lb(): lower_ub; // updax lb of max
        if(elemtmp->box[elemtmp->node->var].diam() > str->preclw) { // bissect and push 
            elemtmp->node->cut(elemtmp->box);
            heap.push(new heap_elem(elemtmp->node->left,elemtmp->node->left_box(elemtmp->box),elemtmp->eval));
            heap.push(new heap_elem(elemtmp->node->right,elemtmp->node->right_box(elemtmp->box),elemtmp->eval));
            delete elemtmp;
        }
        else { // min prec on box reached, keep heap elem
//            cout<<"max :"<< elemtmp->eval<<" for freq: "<<elemtmp->box<<endl;
            if(!midp)
                save_stack.push(elemtmp);
            upper_ub = res[0].ub()>upper_ub?  res[0].ub(): upper_ub;
        }
    }
    if(!midp) // modification of master B&B elem, only if !midpt, keep wmax from saved boxes
    {
        while(!save_stack.empty()) {
            if(save_stack.top()->eval.ub()< lower_ub) {
                save_stack.top()->node->status = __IBEX_OUT__;
            }
            delete save_stack.top();
            save_stack.pop();
        }

        elem->fmax = IntervalVector(1,Interval(lower_ub,upper_ub));
        elem->tree.gather();
    }
    return upper_ub; // return ub on max, used as result for controller parameter midpoint eval
}


int main() {

    // precision:
    double prec(0.01); 
    double preclwmin = 1.e-4;// dynamic initialization in loop
    double stop_criterion(0.1); // stop if distance between uplo and loup lower than stop_criterion

    // init boxes

    IntervalVector IniboxK(3,Interval(-10,10)); // controller parameters init
    IntervalVector Inilw(1,Interval(-2,2)); // u init, 
    vector< pair<polynomial *,polynomial *> > frac_wfree;
    vector< pair<polynomial *,polynomial *> > frac_wset;
    vector< pair<polynomial *,polynomial *> > frac_cst;

    //function definition
    Variable kp,ki,kd,w,u;
    Function Twz1("kp","ki","kd","u","(((33.75)*exp(ln(10)*5*u)+0.036*exp(ln(10)*u)-10.147*exp(ln(10)*3*u))^2+(-43.7524*exp(ln(10)*4*u)+0.1806*exp(ln(10)*2*u))^2)*1/((-33.7524*exp(ln(10)*4*u)+0.0036*kd*exp(ln(10)*2*u)+0.00366*kp*exp(ln(10)*2*u)-kp*exp(ln(10)*4*u)-kd*exp(ln(10)*4*u)-0.0036*ki+1.00006*exp(ln(10)*2*u)*ki+0.14399999999999999999*exp(ln(10)*2*u))^2+((33.75)*exp(ln(10)*5*u)+0.0036*kp*exp(ln(10)*u)+0.00366*exp(ln(10)*u)*ki-(5.9999999999999999995*10^(-5))*kd*exp(ln(10)*3*u)-1.00006*kp*exp(ln(10)*3*u)-exp(ln(10)*3*u)*ki-0.1464*exp(ln(10)*3*u))^2)-1/((-0.01248+exp(ln(10)*2*u))^2+0.024964*exp(ln(10)*2*u))*((0.001248-3.162*exp(ln(10)*2*u))^2+0.0078960996*exp(ln(10)*2*u))");
    Function Twz2("kp","ki","kd","u","(((33.75)*kd*exp(ln(10)*5*u)+33.75*kp*exp(ln(10)*5*u)-43.7524*exp(ln(10)*3*u)*ki-0.14459999999999999999*kd*exp(ln(10)*3*u)-10.1470000000000000005*kp*exp(ln(10)*3*u)+0.036*kp*exp(ln(10)*u)+0.1806*exp(ln(10)*u)*ki)^2+((10.1470000000000000005)*exp(ln(10)*2*u)*ki-43.7524*kp*exp(ln(10)*4*u)-10.0024*kd*exp(ln(10)*4*u)-33.75*exp(ln(10)*4*u)*ki+0.036*kd*exp(ln(10)*2*u)+0.1806*kp*exp(ln(10)*2*u)-0.036*ki)^2)*1/(((33.75)*exp(ln(10)*5*u)-exp(ln(10)*3*u)*ki-(5.9999999999999999995*10^(-5))*kd*exp(ln(10)*3*u)-0.1464*exp(ln(10)*3*u)-1.00006*kp*exp(ln(10)*3*u)+0.0036*kp*exp(ln(10)*u)+0.00366*exp(ln(10)*u)*ki)^2+((1.00006)*exp(ln(10)*2*u)*ki-kp*exp(ln(10)*4*u)-33.7524*exp(ln(10)*4*u)-kd*exp(ln(10)*4*u)+0.14399999999999999999*exp(ln(10)*2*u)+0.0036*kd*exp(ln(10)*2*u)+0.00366*kp*exp(ln(10)*2*u)-0.0036*ki)^2)-((7896.0996)*exp(ln(10)*2*u)+(-3948+exp(ln(10)*2*u))^2)*1/((78.960996)*exp(ln(10)*2*u)+(-39.48+exp(ln(10)*2*u))^2)");
    Function Twz1w("kp","ki","kd","w","1/((w^3*ki+(0.1464)*w^3-(33.75)*w^5-(0.00366)*w*ki-(0.0036)*kp*w+(1.00006)*kp*w^3+(5.9999999999999999995e-5)*kd*w^3)^2+(kd*w^4+kp*w^4-(0.14399999999999999999)*w^2-(1.00006)*w^2*ki+(33.7524)*w^4+(0.0036)*ki-(0.00366)*kp*w^2-(0.0036)*kd*w^2)^2)*(((0.1806)*w^2-(43.7524)*w^4)^2+(-(10.1470000000000000005)*w^3+(0.036)*w+(33.75)*w^5)^2)");
    Function Twz2w("kp","ki","kd","w","((-(10.0024)*kd*w^4-(43.7524)*kp*w^4+(10.1470000000000000005)*w^2*ki-(0.036)*ki+(0.1806)*kp*w^2+(0.036)*kd*w^2-(33.75)*w^4*ki)^2+((33.75)*kp*w^5-(43.7524)*w^3*ki+(33.75)*kd*w^5+(0.1806)*w*ki+(0.036)*kp*w-(10.147)*kp*w^3-(0.14459999999999999999)*kd*w^3)^2)*1/((w^3*ki+(0.1464)*w^3-(33.75)*w^5-(0.00366)*w*ki-(0.0036)*kp*w+(1.00006)*kp*w^3+(5.9999999999999999995e-5)*kd*w^3)^2+(kd*w^4+kp*w^4-(0.14399999999999999999)*w^2-(1.00006)*w^2*ki+(33.7524)*w^4+(0.0036)*ki-(0.00366)*kp*w^2-(0.0036)*kd*w^2)^2)");
    Function w2("kp","ki","kd","w","1/((-39.48+w^2)^2+(78.960996)*w^2)*((-3948+w^2)^2+(7896.0996)*w^2)");
    Function w1("kp","ki","kd","w","((0.0078960996)*w^2+(0.001248-(3.162)*w^2)^2)*1/((0.024964)*w^2+(-0.01248+w^2)^2)");
    
    
    // polynomial definition for berstein eval
    string aff_var[4];aff_var[0] = "kp";aff_var[1] = "ki";aff_var[2] = "kd";
    string var[4];var[0]="w";
    frac_wset.push_back(get_pol(&Twz1w,aff_var,3,var,1));
    frac_wset.push_back(get_pol(&Twz2w,aff_var,3,var,1));
    aff_var[0] = "w";
    var[0] = "kp";var[1] = "ki";var[2] = "kd";
    frac_wfree.push_back(get_pol(&Twz1w,aff_var,1,var,3));
    frac_wfree.push_back(get_pol(&Twz2w,aff_var,1,var,3));

    frac_cst.push_back(get_pol(&w1,aff_var,1,var,0));
    frac_cst.push_back(get_pol(&w2,aff_var,1,var,0));


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
//     cout<<"objective function ok"<<endl;

 // stability criteria 
    Function rS1("kp","ki","kd","33.7524+kp+kd");
    Function rS2("kp","ki","kd", "-0.14399999999999999999-1.00006*ki-0.0036600000000000000002*kp-0.0036*kd+0.02962962962962962963*(33.7524+kp+kd)*(0.1464+ki+1.00006*kp+(5.9999999999999999995*10^(-5))*kd)");
    Function rS3("kp","ki","kd","1/(33.7524+kp+kd)*((0.14399999999999999999+1.00006*ki+0.00366*kp+0.0036*kd)*(-0.14399999999999999999-1.00006*ki-0.0036600000000000000002*kp-0.0036*kd+0.02962962962962962963*(33.7524+kp+kd)*(0.1464+ki+1.00006*kp+(5.9999999999999999995*10^(-5))*kd))+(33.7524+kp+kd)*((0.0036)*ki*(0.1464+ki+1.00006*kp+(5.9999999999999999995*10^(-5))*kd)-((0.00366)*ki+0.0036*kp)*(0.14399999999999999999+1.00006*ki+0.00366*kp+0.0036*kd))*1/(0.1464+ki+1.00006*kp+(5.9999999999999999995*10^(-5))*kd))");

//    cout<<frac_wfree.at(1).first->eval_bernstein(IntervalVector(4,Interval(1)));
//    cout<<frac_wfree.at(1).second->eval_bernstein(IntervalVector(4,Interval(1)));
//    cout<<"eval ok"<<endl;
    
    // stability criteria contractor 

    vector<Ctc*> array_ctc;
    NumConstraint *c3= new NumConstraint(kp,ki,kd,rS1(kp,ki,kd)>=0);
    array_ctc.push_back(new CtcFwdBwd(*c3));
    NumConstraint *c4= new NumConstraint(kp,ki,kd,rS2(kp,ki,kd)>=0);
    array_ctc.push_back(new CtcFwdBwd(*c4));
    NumConstraint *c5= new NumConstraint(kp,ki,kd,rS3(kp,ki,kd)>=0);
    array_ctc.push_back(new CtcFwdBwd(*c5));
    CtcCompo ctc_routh(array_ctc);


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
    
    //ibex objective function definition
    Function Max12(kp,ki,kd,u,ibex::max(Twz1(kp,ki,kd,u),Twz2(kp,ki,kd,u)));
 //   Function Sum(kp,ki,kd,u,ibex::max(Twz1(kp,ki,kd,u),Twz2(kp,ki,kd,u)));




//    Timer::start();

//    for (unsigned i=0;i<1000;i++) {
//        frac_wfree.at(0).first->eval_bernstein(IntervalVector(1,Interval(-11,100)))[0]/frac_wfree.at(0).second->eval_bernstein(IntervalVector(1,Interval(-11,100)))[0];
//        frac_wfree.at(1).first->eval_bernstein(IntervalVector(1,Interval(-11,100)))[0]/frac_wfree.at(1).second->eval_bernstein(IntervalVector(1,Interval(-11,100)))[0];
//    }

//    Timer::stop();
//    cout<<"1000 eval in "<<Timer::VIRTUAL_TIMELAPSE()<<" second "<<endl;
//    Timer::start();
//    for (unsigned i=0;i<1000;i++) {
//        Max12.eval_vector(IntervalVector(1,Interval(-11,100)));
//    }
//    Timer::stop();
//    cout<<"1000 classic eval in "<<Timer::VIRTUAL_TIMELAPSE()<<" second "<<endl;
//    return 0;


//    Function Max123(kp,ki,kd,u,ibex::max(F3u(kp,ki,kd,u),Max12(kp,ki,kd,u)));

    //other variables for algo

//    Function Prob("kp","ki","w","((kp+w-2)^6+0.2)*ln(1+(kp+w)^2)+((ki+w-2)^6+0.2)*ln(1+(ki+w)^2)"  );


    double lower_ub(POS_INFINITY); // loup initialization

    LargestFirst lf;
    SetIntervalReg tree(Inilw,preclwmin,false); // first tree initialization that contain wmax
//    SetIntervalReg tree(Iniprob2,preclwmin,false);
    tree.root->status = __IBEX_UNK__; // wmax init
    cost_uplo cu; // heap policy 
    Heap<list_elem> list(cu);

    list_elem *elem = new list_elem(IniboxK, tree,IntervalVector(1,Interval::ALL_REALS)); // list of master B&B, first element init
//    list_elem *elem = new list_elem(Iniprob, tree,IntervalVector(1,Interval::ALL_REALS));
    list.push(elem);

    Vector respoint(3); // controller parameter midpoint
//    Vector respoint(2);

    optiw str(lower_ub,preclwmin,&Max12,frac_wset,frac_wfree,frac_cst,false); // slave B&B class init, see class definition


    double vol_rejected(0); // % of IniK proved to not contain solution
    list_elem * elemtmp; // current element of master B&B
//    Vector sol(4);
//    //matlab solution
//    sol[0] = 0.3998;sol[1] = 0.4806;sol[2] = -0.1952;
//    cout<<"evaluation of matlab result: "<<Sum.eval_vector(sol)<<endl;
//    return 0;

//*********** Master B&B **************
    Timer::start();
    while(!list.empty()) { 
        elemtmp = list.pop();
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
            cout<<"loup : "<<lower_ub<<" get for point: "<<respoint<<" uplo: "<<elemtmp->fmax.lb()<< " volume rejected: "<<vol_rejected/IniboxK.volume()*100<<endl;
        }
        if(elemtmp->box.is_empty()) {
            delete elemtmp;
            continue;
        }


        // **************** compute (wmax,fmax) *******************
//        Function g(w,Max12(elemtmp->box[0],elemtmp->box[1],elemtmp->box[2],w));
//        Function g(w,Prob(elemtmp->box[0],elemtmp->box[1],w));
//        str.f = &g;

        best_upper_bound_forall_w(&str,elemtmp,false);
//        cout<<"res found for box "<<elemtmp->box<<" : "<<elemtmp->fmax<<endl;
        if(elemtmp->fmax[0] == Interval::EMPTY_SET) { // maximum is guaranteed to be higher than the current lower_up for box elemtmp->box
            vol_rejected+= elemtmp->box.volume();
            cout<<"loup : "<<lower_ub<<" get for point: "<<respoint<<" uplo: "<<elemtmp->fmax[0].lb()<< " volume rejected: "<<vol_rejected/IniboxK.volume()*100<<endl;
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
//        midp = sol;
//        str.preclw = 0.00001;
        // if couth criteria ok for midp
        if((rS1.eval_vector(midp))[0].lb()>=0 && (rS2.eval_vector(midp))[0].lb()>=0 && (rS3.eval_vector(midp))[0].lb()>=0) {// routh criterion ok for midpoint
//            Function max_mid(w,Max12(Interval(midp[0]),Interval(midp[1]),Interval(midp[2]),w));
//            Function max_mid(w,Prob(Interval(midp[0]),Interval(midp[1]),w));
//            str.f = &max_mid;
            if(str.preclw>0.001) str.preclw=0.001; // forces low precision to get tight enclosure
            double max = best_upper_bound_forall_w(&str,elemtmp,true);
//            cout<<" max of mid point "<<midp<<": "<<max<<endl;
//            return 0;
            //            cout<<"fmax: "<<elemtmp->fmax<<"lower_ub: "<<lower_ub<<endl;

            if(max<lower_ub) { // if better solution found
                lower_ub = max;
                respoint = midp;
                cout<<"loup : "<<lower_ub<<" get for point: "<<respoint<<" uplo: "<<elemtmp->fmax.lb()<< " volume rejected: "<<vol_rejected/IniboxK.volume()*100<<endl;
            }
        }

        if(elemtmp->box.max_diam()>prec) { // bissect and push
            pair<IntervalVector,IntervalVector> boxes = lf.bisect(elemtmp->box);

            list.push(new list_elem(boxes.first,elemtmp->tree,elemtmp->fmax));
            list.push(new list_elem(boxes.second,elemtmp->tree,elemtmp->fmax));
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

