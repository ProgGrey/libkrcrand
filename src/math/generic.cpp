#include "../../headers/math.hpp"
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>

double horner(double x, double *coefs, uint_fast8_t n)
{
    double val = coefs[0];
    for (uint_fast8_t k = 1; k < n; k++) {
        val = fma(x, val,  coefs[k]);
    }
    return val;
}

double horner1(double x, double* coefs, uint_fast8_t n)
{
    double val = x + coefs[0];
    for (uint_fast8_t k = 1; k < n; k++) {
        val = fma(x, val, coefs[k]);
    }
    return val;
}

namespace krcrand{

double uniform01_exclude0(uint64_t a)
{
    return fma(double(a), 5.4210108624275221703e-20, 5.4210108624275221700e-20);
}

// Based on algorithm from Cephes 2.8 math library by Stephen L. Moshier 1984, 1995, 2000
// http://www.netlib.org/cephes/
// Max relative error is 1.523918708144333e-16 for normal numbers and 1.039444891344354e-16
// for subnormals (10^10 test runs).
double unsafe_log(double x)
{
    //Subnormal numbers correction. Max rel. error is 1.039444891344354e-16 for 10^10 test runs.
    double addition = 0;
    if (x <= 2.2250738585072009e-308) {
        x *= 4503599627370496.0;
        addition = -36.0436533891171560897;
    }
    //x=2.2e-308
    // Log(2)
    const double loge2_p = 0.693359375;
    const double loge2_n = 2.121944400546905827679E-4;
    //const double loge2 = 0.693147180559945309417;
    //sqrt(2)/2
    const double sqrt2b2 = 7.07106781186547524400E-1;
    //Constants for polynom
    double a[6] = { 1.01875663804580931796E-4, 4.97494994976747001425E-1, 4.70579119878881725854E0,
                      1.44989225341610930846E1, 1.79368678507819816313E1,   7.70838733755885391666E0 };
    double b[6] = { 1.12873587189167450590E1, 4.52279145837532221105E1, 8.29875266912776603211E1, // first cooficient is 1 here
                    7.11544750618563894466E1, 2.31251620126765340583E1 };
    // unpack double
    uint64_t code;
    memcpy(&code, &x, sizeof(double));
    uint64_t man_code = (code & 0x000FFFFFFFFFFFFF) | 0x3FE0000000000000;
    uint64_t exp_code = code >> 52 | 0x4330000000000000;
    double man;
    double exp;
    memcpy(&man, &man_code, sizeof(double));
    memcpy(&exp, &exp_code, sizeof(double));
    exp = exp - 1022 - 4503599627370496.0;
    // Fix interval
    if (man < sqrt2b2) {
        man += man;
        exp--;
    }
    man -= 1;
    //Rational approximation
    double num = horner(man, a, 6);
    double den = horner1(man, b, 5);

    double man2 = man * man;

    double res = fma(man, (man2 * num / den), -loge2_n * exp);
    //Splitting increases precision:
    //res -= loge2_n * exp;
    res = fma(-0.5, man2, res);
    //res -= 0.5 * man2;
    res += man;
    //res += loge2 * exp; Using property of logarithms for powers.
    //res += loge2_p * exp;
    res = fma(loge2_p, exp, res);
    return res + addition;
}

}