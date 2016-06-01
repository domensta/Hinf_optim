#include "icomplex.h"
icomplex::icomplex():real(Interval::EMPTY_SET),imag(Interval::EMPTY_SET)
{}

icomplex::icomplex(Interval rp, Interval ip):real(rp),imag(ip)
{}

icomplex operator+(icomplex const &first, icomplex const &second) {
    return icomplex(first.real+second.real,first.imag+second.imag);
}

icomplex operator-(icomplex const &first, icomplex const &second) {
    return icomplex(first.real-second.real,first.imag-second.imag);
}

icomplex operator*(icomplex const &first, icomplex const &second) {
    return icomplex(first.real*second.real-first.imag*second.imag,first.real*second.imag+first.imag*second.real);
}

icomplex operator/(icomplex const &first, icomplex const &second) {
    Interval real_num = first.real*second.real+first.imag*second.imag;
    Interval imag_num = first.imag*second.real-first.real*second.imag;
    Interval den = ibex::pow(second.real,2)+ibex::pow(second.imag,2);
    return icomplex(real_num/den,imag_num/den);
}

icomplex operator+(icomplex const &first, Interval const &second) {
    return icomplex(first.real+second,first.imag);
}

icomplex operator-(icomplex const &first, Interval const &second) {
    return icomplex(first.real-second,first.imag);
}

icomplex operator*(icomplex const &first, Interval const &second) {
    return icomplex(first.real*second,first.imag*second);
}

icomplex operator/(icomplex const &first, Interval const &second) {
    return icomplex(first.real/second,first.imag/second);
}

//icomplex operator+(icomplex const &first, T const &second) {
//    return icomplex(first.real+T,first.imag);
//}

//icomplex operator-(icomplex const &first, T const &second) {
//    return icomplex(first.real-T,first.imag);
//}

//icomplex operator*(icomplex const &first, T const &second) {
//    return icomplex(first.real*T,first.imag*T);
//}

//icomplex operator/(icomplex const &first, T const &second) {
//    return icomplex(first.real/T,first.imag/T);
//}


std::ostream& operator<<(std::ostream& os, const icomplex& ic) {
    Interval r =ic.real;
    Interval i = ic.imag;
    os<<r<<" + "<<i<<"*i";
    return os;
}

bool icomplex::is_real() {
    return(imag==Interval(0)||imag.is_empty());
}

bool icomplex::is_imag(){
     return(real==Interval(0)||real.is_empty());
}

Interval icomplex::abs() {
    Interval res = ibex::sqrt(ibex::pow(real,2)+ibex::pow(imag,2));
    return res;
}

Interval icomplex::real_part(){
    return real;
}

Interval icomplex::imag_part(){
    return imag;
}
