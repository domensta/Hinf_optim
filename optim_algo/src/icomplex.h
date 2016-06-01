#ifndef __ICOMPLEX_H__
#define __ICOMPLEX_H__

#include "ibex.h"
#include <stdio.h>
#include <stdlib.h>
using namespace ibex;
using namespace std;




class icomplex{
public:
//    template<typename T>
    Interval real;
    Interval imag;
    icomplex(Interval rp,Interval ip);
    icomplex();
    friend icomplex operator+(icomplex const &first, icomplex const &second);
    friend icomplex operator-(icomplex const &first, icomplex const &second);
    friend icomplex operator*(icomplex const &first, icomplex const &second);
    friend icomplex operator/(icomplex const &first, icomplex const &second);
    friend icomplex operator+(icomplex const &first, Interval const &second);
    friend icomplex operator-(icomplex const &first, Interval const &second);
    friend icomplex operator*(icomplex const &first, Interval const &second);
    friend icomplex operator/(icomplex const &first, Interval const &second);
//    friend icomplex operator+(icomplex const &first, const T &second);
//    friend icomplex operator-(icomplex const &first, const T &second);
//    friend icomplex operator*(icomplex const &first, const T &second);
//    friend icomplex operator/(icomplex const &first, const T &second);
    ibex::Interval abs();
    bool is_real();
    bool is_imag();
    Interval real_part();
    Interval imag_part();


};

std::ostream& operator<<(std::ostream& os, const icomplex& ic);







#endif
