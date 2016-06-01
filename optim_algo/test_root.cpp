
#include "src/univar_polynomial.h"
#include "/home/kiwi/Tools/ibex-lib-ibex-2.1.18/OUT/include/ibex/ibex.h"
#include <stdio.h>
#include <stdlib.h>

using namespace ibex;
using namespace std;
int main() {

    string var[4]; var[0] = "kp"; var[1] = "ki"; var[2] = "kd"; var[3] = "tf";
    univar_polynomial cltf_pol("kp*s^3 + ki*s^2 +kd*s + tf","s",var,4);
    int nbeig(3);
    IntervalVector midp(nbeig+1);
    midp[0] = Interval(1,1.3);midp[1] = Interval(1.8,2);midp[2] = Interval(3);midp[2.9,3] = Interval(3.6,4);
    pair<icomplex,Interval> roots[nbeig];
    if(cltf_pol.get_roots(IntervalVector(midp),roots,nbeig) == 1) {
	for(int i=0;i<nbeig;i++) {
                cout<<"center: ("<<roots[i].first.real<<","<<roots[i].first.imag<<"), diameter: "<<roots[i].second<<endl;
	}
    }
   return 0;
}
