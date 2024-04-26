#include "../../headers/math.hpp"
#include <cmath>
#include <cstdint>
#include <bit>
#include <iostream>

using namespace std;

double krcrand::signum(double x)
{
    uint64_t code = bit_cast<uint64_t>(x) & 0x8000000000000000;
    code |= bit_cast<uint64_t>(double(1));
    return bit_cast<double>(code);
}

int p3_solve_viete(long double &a, long double &b, long double &c, double &x1, double &x2, double &x3)
{
    long double q = (a*a - 3.0*b)/9.0;
    long double r = (2.0*a*a*a-9.0*a*b+27.0*c)/54.0;
    long double s = q*q*q-r*r;
    if(s > 0){
        //cout << "s > 0\n";//+
        long double phi = 1.0/3.0*acosl(r/(q*sqrt(q)));
        x1 = -2.0*sqrtl(q)*cosl(phi)-a/3.0;
        x2 = -2.0*sqrtl(q)*cosl(phi+2.0/3.0*M_PIl)-a/3.0;
        x3 = -2.0*sqrtl(q)*cosl(phi-2.0/3.0*M_PIl)-a/3.0;
        if(x1 > x3){
            swap(x1,x3);
        }
        if(x1 > x2){
            swap(x1,x2);
        }
        if(x2 > x3){
            swap(x2,x3);
        }
        return 3;
    } else if(s == 0){
        //cout << "s == 0\n";//+
        x1 = -2.0*cbrtl(r)-a/3.0;
        x2 = cbrtl(r)-a/3.0;
        if(x1 > x2){
            swap(x1,x2);
        }
        return 2;
    } else{//s<0
        //cout << "s < 0, ";
        if(q>0){
            //cout << "q > 0\n";//+
            long double phi = 1.0/3.0*acoshl(abs(r)/(q*sqrt(q)));
            x1 = -2*krcrand::signum(r)*sqrtl(q)*coshl(phi)-a/3.0;
            return 1;
        } else if(q==0){
            //cout << "q == 0\n";//+
            x1 = -cbrtl(c-a*a*a/27.0)-a/3.0;
            return 1;
        }else{//q<0
            //cout << "q < 0\n";//+
            long double phi = 1.0/3.0*asinhl(abs(r)/(abs(q)*sqrtl(abs(q))));
            x1 = -2*krcrand::signum(r)*sqrtl(abs(q))*sinhl(phi) - a/3.0;
            return 1;
        }
    }
}

int krcrand::p3_solve(double a, double b, double c, double d, double &x1, double &x2, double &x3)
{
    long double _a = b/a;
    long double _b = c/a;
    long double _c = d/a;
    return p3_solve_viete(_a, _b, _c, x1, x2, x3);
}


int krcrand::p2_solve(double a, double b, double c, double &x1, double &x2)
{
    long double a_= a, b_= b, c_= c;
    long double D = b_*b_-4*a_*c_;
    if(D < 0.0){
        return 0;
    }
    D = sqrtl(D);
    if(static_cast<double>(D) == 0.0){
        x1 = -0.5*b_/a_;
        return 1;
    } else{
        x1 = -0.5*(b_+D)/a_;
        x2 = (-b_+D)/a_*0.5;
        if(x1 > x2){
            swap(x1,x2);
        }
        return 2;
    }
}