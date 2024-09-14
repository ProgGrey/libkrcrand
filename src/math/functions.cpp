#include "../../headers/math.hpp"
#include <cmath>
#include <cstdint>
#include <bit>
#include <iostream>
#include <ccomplex>
#include <vector>
#include <algorithm>
#include <set>

using namespace std;

double krcrand::signum(double x)
{
    uint64_t code = bit_cast<uint64_t>(x) & 0x8000000000000000;
    code |= bit_cast<uint64_t>(double(1));
    return bit_cast<double>(code);
}

int p3_solve_viete(long double const &a, long double const &b, long double const &c, double &x1, double &x2, double &x3)
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

int krcrand::p3_solve_cardano(long double const &a, long double const &b, long double const &c, long double const &d, double &x1, double &x2, double &x3)
{
    long double p = (3.0*a*c-b*b)/(3.0*a*a);
    long double q = (2.0*b*b*b-9.0*a*b*c+27.0*a*a*d)/(27.0*a*a*a);
    long double Q = p*p*p/27.0 + q*q/4.0;
    long double alpha = cbrtl(-q*0.5+sqrt(Q));
    long double beta = cbrtl(-q*0.5-sqrt(Q));
    int solutions;
    long double y1 = 0, y2 = 0, y3 = 0;
    y1 = alpha + beta;
    if(Q > 0){
        solutions = 1;
    }else if ((Q == 0) && (q == 0)){
        solutions = 1;
    } else if (Q==0.0){
        solutions = 2;
        y2 = -0.5*(alpha+beta);
    }else{//Q < 0
        solutions = 3;
        long double sQ = sqrtl(-Q);
        //determine alpha and beta by using definition of cbrt for complex numbers:
        //cout << -2.0*sQ/q << ' ' << atanl(-2.0*sQ/q) << endl;
        long double phi_a = atanl(-2.0*sQ/q);
        long double phi_b = atanl(2.0*sQ/q);
        long double r = sqrt(q*q/4-sQ*sQ);//|alpha| = |beta|
        //long double apb = cbrtl(r)*(cosl(phi_a/3)+cosl(phi_b/3));// (alpha + beta)
        long double apb = cbrtl(r)*(cosl(phi_a/3)+cosl(phi_b/3));// (alpha + beta)
        y1 = apb;
        apb = apb*0.5;// 
        long double amb = cbrtl(r)*(sinl(phi_b/3)-sinl(phi_a/3))*0.5*sqrtl(3);// (alpha - beta)/2
        y2 = -apb+amb;
        y3 = -apb-amb;
        /*
        if(a > 0){
            y1 += -1.0/a+1.0/3.0/a;
            y2 += 1.0/a+1.0/3.0/a;
            y3 += -1.0/a+1.0/3.0/a;
        }//*/
        /*
        if(y1 > y3){
            swap(y1,y3);
        }
        if(y1 > y2){
            swap(y1,y2);
        }
        if(y2 > y3){
            swap(y2,y3);
        }//*/
    }
    x1 = y1 - b/(3*a);
    x2 = y2 - b/(3*a);
    x3 = y3 - b/(3*a);
    
    //cout << "res" << y1*y1*y1*a+y1*y1*b+y1*c+d << endl;
    return solutions;
}

int krcrand::p3_solve_modern(long double const &a, long double const &b, long double const &c, long double const &d, double &x1, double &x2, double &x3)
{
    long double tmp = a*d;
    tmp *= tmp;
    long double D = 27.0*tmp - ((18.0*a*c-4.0*b*b)*b)*d + 4.0*a*c*c*c;
    tmp = b*c;
    tmp *= tmp;
    D -= tmp;
    //cout << "D " << D << endl;
    //cout << endl << b*c << ' ' << -tmp << ' ' <<  D << endl;
    long double B = 4.0*(b*b-3.0*a*c);
    //cout << "B " << B << endl;
    if(D >= 0.0){
        D = sqrt(D);
        long double A = cbrtl(8*((1.5*sqrtl(3.0)*D-13.5*a*d+4.5*c*b)*a-b*b*b));
        //cout << "A " << A << endl;
        x1 = (B+(-2.0*b+A)*A)/(6.0*a*A);
        //long double A2B = A*A-B;
        //cout << A2B << endl;
        //if(abs(A2B) < 1e-15){
        tmp = A+2.0*b;
        tmp *= tmp;
        long double x_c = (12.0*a*c-tmp)/(12.0*a*A);
        if(abs(((a*x_c+b)*x_c+c)*x_c + d) < 1e-15){
            x2 = x_c;
            return 2;
        }
        return 1;
    } else{//D<0
        D = sqrt(-D);
        vector<double> results;
        long double A_r = (-13.5*a*a*d+(4.5*a*c-b*b)*b)*8.0;
        long double A_i = 12.0*sqrtl(3.0)*a*D;
        //cout << "D " << D << "i\n";
        //cout << "A " << A_r <<" + " << A_i << "i\n";
        //cout<< "k,m " << k << ' ' << mult << endl;
        long double phi = atanl(A_i/A_r);// + 2.0*M_PIl*k;
        long double r = sqrtl(A_i*A_i+A_r*A_r);
        long double A_r_ = cbrtl(r)*cosl(phi/3.0);
        long double A_i_ = cbrtl(r)*sinl(phi/3.0);
        // Fix cbrt for complex number: sometimes signs are incorrect:
        A_r = signum((A_r_*A_r_*A_r_ - 3.0*A_r_*A_i_*A_i_)/A_r)*A_r_;
        A_i = signum((3.0*A_r_*A_r_*A_i_-A_i_*A_i_*A_i_)/A_i)*A_i_;
        long double A2_r = A_r*A_r-A_i*A_i;
        long double A2_i = 2.0*A_r*A_i;
        //long double mult = -1;
            //A_r = -1.0*A_r;
            //A_i = -1.0*A_i;
        long double n_r = B-2.0*b*A_r+A2_r;
        long double n_i = -2.0*b*A_i+A2_i;
        long double d_r = 6.0*a*A_r;
        long double d_i = 6.0*a*A_i;
        long double r_r = (n_r*d_r+n_i*d_i)/(d_r*d_r + d_i*d_i);
        if(abs(((a*r_r+b)*r_r+c)*r_r + d) < 1e-15){
            results.push_back(r_r);
        }
        d_r *= 2.0;
        d_i *= 2.0;
        long double t_r = A_r + 2.0*b;
        //long double t_i = A_i;
        long double t2_r = t_r*t_r-A_i*A_i;
        //long double t2_i = 2.0*t_i*t_i;
        long double left_r = 12.0*a*c - t2_r;
        long double left_i = -2.0*A_i*(A_r+2*b);
        //cout << "L " << left_r <<" + " << left_i << "i\n";
        long double right_r = -(A2_i)*sqrtl(3.0);
        long double right_i = (A2_r-B)*sqrtl(3.0);
        //cout << "R " << right_r <<" + " << right_i << "i\n";
        n_r = left_r + right_r;
        n_i = left_i + right_i;
        //cout << "num " << n_r <<" + " << n_i << "i\n";
        r_r = (n_r*d_r+n_i*d_i)/(d_r*d_r + d_i*d_i);
        //cout << r_r << endl;
        if(abs(((a*r_r+b)*r_r+c)*r_r + d) < 1e-15){
            results.push_back(r_r);
        }
        n_r = left_r - right_r;
        n_i = left_i - right_i;
        r_r = (n_r*d_r+n_i*d_i)/(d_r*d_r + d_i*d_i);
        //cout << r_r << endl;
        if(abs(((a*r_r+b)*r_r+c)*r_r + d) < 1e-15){
            results.push_back(r_r);
        }
        if(results.size() >= 1){
            x1 = results[0];
        }
        if(results.size() >= 2){
            x2 = results[1];
        }
        if(results.size() >= 3){
            x3 = results[2];
        }
        return results.size();
    }
}

int krcrand::p2_solve(double a, double b, double c, double &x1, double &x2)
{
    long double a_= a, b_= b, c_= c;
    if(abs(a) <= 1e-15 && (abs(b) <= 1e-15)){
        return 0;
    }else if(abs(b)<= 1e-15 && abs(c)<= 1e-15){
        x1 = 0.0;
        return 1;
    }else if(abs(a) <= 1e-15){
        x1 =  -c/b;
        return 1;
    }else if(abs(c)<= 1e-15){
        x1 = min(0.0,-b/a);
        x2 = max(0.0,-b/a);
        return 2;
    } else if(abs(b)<= 1e-15){
        long double D = -c*a;
        if (D < 0.0){
            return 0;
        }
        D = sqrtl(D);
        x1 = a > 0.0 ? D/a : -D/a;
        x2 = a < 0.0 ? D/a : -D/a;
        return 2;
    }
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


int krcrand::p3_solve(double a, double b, double c, double d, double &x1, double &x2, double &x3)
{
    int res;
    if(abs(a) <= 1e-15){
        return p2_solve(b,c,d, x1,x2);
    } else if(abs(d) <= 1e-15){
        res = p2_solve(a,b,c, x1,x2) + 1;
        x3 = 0.0;
    }else{
        res = p3_solve_modern(static_cast<long double>(a), static_cast<long double>(b), static_cast<long double>(c), static_cast<long double>(d), x1, x2, x3);
    }
    if(res == 2){
        if(x1 > x2){
            swap(x1,x2);
        }
    }else if(res == 3){
        if(x1 > x3){
            swap(x1,x3);
        }
        if(x1 > x2){
            swap(x1,x2);
        }
        if(x2 > x3){
            swap(x2,x3);
        }
    }//*/
    return res;
}