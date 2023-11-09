#include "../../headers/math.hpp"
#include <cstdint>

#define uniform01_exclude0_function(name) \
d_type name(i_type a)\
{\
    return fma_d(conv_u_to_d(a), set1_d(5.4210108624275221703e-20), set1_d(5.4210108624275221700e-20));\
}

#define horner_function(name) \
d_type name(d_type x, double *coefs, uint_fast8_t n) \
{\
    d_type val = set1_d(coefs[0]); \
    for (uint_fast8_t k = 1; k < n; k++) { \
        val = fma_d(x, val,  set1_d(coefs[k])); \
    }\
    return val;\
}


#define horner1_function(name) \
d_type name(d_type x, double* coefs, uint_fast8_t n)\
{\
    d_type val = add_d(x, set1_d(coefs[0]));\
    for (uint_fast8_t k = 1; k < n; k++) {\
        val = fma_d(x, val, set1_d(coefs[k]));\
    }\
    return val;\
}

#ifndef __AVX512_BRANCH__
#define branch(flag_var, if_true, if_false) or_d(and_d(flag_var, if_true), notand_d(flag_var, if_false))
#define zero_branch(flag_var, if_true) and_d(flag_var, if_true)
#else
#define branch(flag_var, if_true, if_false) mov_mask_d(flag_var, if_true, if_false)
#define zero_branch(flag_var, if_true) mov_mask_d(flag_var, if_true, set1_d(0.0))
#endif
// Based on algorithm from Cephes 2.8 math library by Stephen L. Moshier 1984, 1995, 2000
// http://www.netlib.org/cephes/
// Max relative error is 1.523918708144333e-16 for normal numbers and 1.039444891344354e-16
// for subnormals (10^10 test runs).
#define log_function(name, horner_name, horner1_name) \
d_type name(d_type x)\
{\
    /*Subnormal numbers correction. Max rel. error is 1.039444891344354e-16 for 10^10 test runs.*/\
    auto sub_flag = le_d(x, set1_d(2.2250738585072009e-308));\
    d_type sub_norm = mull_d(x, set1_d(4503599627370496.0));\
    d_type addition = zero_branch(sub_flag, set1_d(-36.0436533891171560897));\
    /* Choose branch:*/\
    x = branch(sub_flag, sub_norm, x);\
    /* Log(2) */\
    const double loge2_p = 0.693359375;\
    const double loge2_n = 2.121944400546905827679E-4;\
    /*sqrt(2)/2*/\
    const double sqrt2b2 = 7.07106781186547524400E-1;\
    /*Constants for polynoms*/\
    double a[6] = { 1.01875663804580931796E-4, 4.97494994976747001425E-1, 4.70579119878881725854E0,\
                      1.44989225341610930846E1, 1.79368678507819816313E1,   7.70838733755885391666E0 };\
    double b[6] = { 1.12873587189167450590E1, 4.52279145837532221105E1, 8.29875266912776603211E1, /* first cooficient is 1 here */\
                    7.11544750618563894466E1, 2.31251620126765340583E1 };\
    /* unpack double*/\
    i_type code = to_i(x);\
    i_type man_code = or_i(and_i(code, set1_i(0x000FFFFFFFFFFFFF)), set1_i(0x3FE0000000000000));\
    i_type exp_code = or_i(shr_i64(code, 52), set1_i(0x4330000000000000));\
    d_type man = to_d(man_code);\
    d_type exp = to_d(exp_code);\
    exp = sub_d(exp, set1_d(4503599627370496.0+1022));\
    /* Fix interval*/\
    auto sqt2_flag = lt_d(man, set1_d(sqrt2b2));\
    d_type manc = add_d(man, man);\
    d_type expc = sub_d(exp, set1_d(1.0));\
    man = branch(sqt2_flag, manc, man);\
    exp = branch(sqt2_flag, expc, exp);\
    man = sub_d(man, set1_d(1.0));\
    /*Rational approximation*/\
    d_type num = horner_name(man, a, 6);\
    d_type den = horner1_name(man, b, 5);\
    d_type man2 = mull_d(man, man);\
    d_type tmp = div_d(mull_d(man2, num), den);\
    d_type res = fma_d(man, tmp, mull_d(set1_d(-loge2_n), exp));\
    /*Splitting increases precision:*/\
    /*res -= loge2_n * exp;*/\
    res = fma_d(set1_d(-0.5), man2, res);\
    /*res -= 0.5 * man2;*/\
    res = add_d(res, man);\
    /*res += loge2 * exp; Using property of logarithms for powers.\
     res += loge2_p * exp;*/\
    res = fma_d(set1_d(loge2_p), exp, res);\
    return add_d(res, addition);\
}

