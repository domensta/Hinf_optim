#include "/home/kiwi/Tools/ibex-lib-ibex-2.1.18/OUT/include/ibex/ibex.h"
#include <stdio.h>
#include <stdlib.h>
#include "src/vibes.h"

using namespace ibex;
using namespace std;
int main() {
    Function Twz1("p","u","1/((-0.037950000000000000003-0.36*exp(ln(10)*4*u)*p-0.49662104479999999998*exp(ln(10)*4*u)+12.9076464*exp(ln(10)*2*u)+0.0046*exp(ln(10)*2*u)*p)^2+((8.309202)*exp(ln(10)*u)+0.17848800000000000001*exp(ln(10)*5*u)-8.186280680000000001*exp(ln(10)*3*u)-1.001656*exp(ln(10)*3*u)*p)^2)*(((0.089244000000000000004)*exp(ln(10)*5*u)-0.228068*exp(ln(10)*3*u)-0.6656*exp(ln(10)*3*u)*p)^2+(-0.18*exp(ln(10)*4*u)*p-0.33000448*exp(ln(10)*4*u)+0.46*exp(ln(10)*2*u)*p)^2)");
    Function Twz2("p","u","(((0.082499999999999999996)*exp(ln(10)*u)*p-0.063809459999999999995*exp(ln(10)*3*u)-0.08184*exp(ln(10)*3*u)*p)^2+(-0.040576271999999999997*exp(ln(10)*4*u)+0.0409035*exp(ln(10)*2*u)+0.1287*exp(ln(10)*2*u)*p)^2)*1/((-8.25-0.17848800000000000001*exp(ln(10)*4*u)+8.184*exp(ln(10)*2*u)+exp(ln(10)*2*u)*p)^2+((12.87)*exp(ln(10)*u)-0.4958*exp(ln(10)*3*u)-0.36*exp(ln(10)*3*u)*p)^2)");
    Function Twz3("p","u","((12.960000000000000001)*exp(ln(10)*4*u)+100.0*exp(ln(10)*2*u))*1/((-8.25-0.17848800000000000001*exp(ln(10)*4*u)+8.184*exp(ln(10)*2*u)+exp(ln(10)*2*u)*p)^2+((12.87)*exp(ln(10)*u)-0.4958*exp(ln(10)*3*u)-0.36*exp(ln(10)*3*u)*p)^2)");
    Function RS1("p","0.4958+(0.36)*p");
    Function RS2("p","1/(0.4958+(0.36)*p)*(1.7604866399999999998+(3.4420400000000000002)*p+(0.36)*p^2)");
    Function RS3("p","(20.629467526799999998+(41.3540028)*p+(3.5640000000000000002)*p^2)*1/(1.7604866399999999998+(3.4420400000000000002)*p+(0.36)*p^2)");
    Function RS4("p","1/(0.4958+(0.36)*p)*(4.0903499999999999998+(2.97)*p)");
    
    Interval p(0.8,1.2);
    Interval u(-3,3);
    IntervalVector ev(2,p);
    ev[1] = u;
    IntervalVector res1(1,Interval::EMPTY_SET);
    IntervalVector res2(1,Interval::EMPTY_SET);
    IntervalVector res3(1,Interval::EMPTY_SET);
    int step(7000),prm_step(60);
    for(unsigned i=0;i<step;i++) {
	for(unsigned j=0;j<prm_step;j++) {
	    ev[0] = Interval(p.diam()/prm_step*j+p.lb());
	    ev[1] = Interval(u.diam()/step*i+u.lb());
	    res1|=Twz1.eval_vector(ev);
	    res2|=Twz2.eval_vector(ev);
	    res3|=Twz3.eval_vector(ev);
	}
    }
    cout<<"Twz1: "<<res1<<endl;
    cout<<"Twz2: "<<res2<<endl;
    cout<<"Twz2: "<<res3<<endl;

    cout<<"RS1: "<<RS1.eval_vector(IntervalVector(1,p))<<endl;
    cout<<"RS2: "<<RS2.eval_vector(IntervalVector(1,p))<<endl;
    cout<<"RS3: "<<RS3.eval_vector(IntervalVector(1,p))<<endl;
    cout<<"RS4: "<<RS4.eval_vector(IntervalVector(1,p))<<endl;

    return 0;

}

