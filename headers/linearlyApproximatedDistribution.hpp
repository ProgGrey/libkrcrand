#include "platform.hpp"
#include "math.hpp"
#include <typeinfo>
#include "generators.hpp"
#include <cmath>
#include <vector>

namespace krcrand{

template<typename GenType> class LAD{
    private:
    double *a = nullptr, *b = nullptr, *c = nullptr;//a b and c coefficients of 0.5*a*x^2+b*x+c
    double *x = nullptr;
    double left_s, left_b, left_a;
    double right_s, right_b, right_a;

    bool is_zero_derivative = false;

    GenType lad_generator;

    protected:
    double help_dist(double &rm)
    {
        uint64_t raw = lad_generator();
        unsigned int num = raw >> tbl_shift;
        if((num == 0) && (is_zero_derivative)){//left approximation
            double u = uniform01_exclude01(raw);
            double result = unsafe_lambertw0(left_a*left_s*sqrt(u)*exp(left_a*left_s*0.5)*0.5/sqrt(left_b))*2/left_a;
            rm = left_s*left_s/(left_b*exp(-left_a*(left_s-result))*result*(left_a*result+2.0)*max_val);
            return result;
        } else if(num == 255){// right approximation
            double u = uniform01_exclude01(raw);
            double result = right_s-unsafe_log((1.0-u)/right_b)/right_a;
            rm = 1.0/(right_b*right_a*exp((right_s-result)*right_a)*max_val);
            return result;
        } else{// main approximation
            double u = uniform01(raw);
            double x1,x2;
            double D = sqrt(b[num]*b[num]-2.0*a[num]*(c[num] - u));
            x1 = (-b[num]-D)/a[num];
            x2 = (-b[num]+D)/a[num];
            //p2_solve(0.5*a[num], b[num], c[num] - u, x1, x2);
            double result;
            double x1_dist = std::max(abs(x[num] - x1), abs(x[num+1] - x1));
            double x2_dist = std::max(abs(x[num] - x2), abs(x[num+1] - x2));
            if(x1_dist < x2_dist){
                result = x1;
            } else{
                result = x2;
            }
            rm = 1.0/((a[num]*result+b[num])*max_val);
            return result;
        }
    }

    double max_val;
    // Points where maximum of f(x)/g(x) is possible.
    virtual std::vector<double> left_max(double a, double b, double s) = 0;
    virtual std::vector<double> med_max(double a, double b) = 0;
    virtual std::vector<double> right_max(double a, double b, double s)  = 0;

    GenType::GeneratorStateType init_lad(GenType::GeneratorStateType gs, bool is_zero_derivative)
    {
        this->is_zero_derivative = is_zero_derivative;
        //build tables
        double p = 1.0/tbl_size;
        x = new double[tbl_size];
        double *yl = new double[tbl_size];
        for(unsigned int k = 0; k < tbl_size; k++){
            
            x[k] = qf(1.0/tbl_size*k);
            yl[k] = pdf(x[k]);
        }
        double *yr =  new double[tbl_size];
        //Fix squares
        for(unsigned int k = 0; k < tbl_size-1; k++){
            double h = (yl[k] + yl[k+1])*0.5+p/(x[k]-x[k+1]);
            double m = std::min(std::min(h, yl[k]),yl[k+1]);
            yr[k] = yl[k+1] - m;
            yl[k] -= m;
            if(h != m){
                //angle correction
                if(yr[k] > yl[k]){
                    h = yr[k]-yl[k]-2*p/(x[k+1]-x[k]);
                    yr[k] -= h;
                } else{
                    h = yl[k] - yr[k]+2*p/(x[k+1]-x[k]);
                    yl[k] -= h;
                }
            }
        }
        // Fix square for zero interval
        //compute constants
        a = new double[tbl_size-1];
        b = new double[tbl_size-1];
        c = new double[tbl_size-1];
        for(unsigned int k = 0; k < tbl_size-1; k++){
            a[k] = (yl[k]-yr[k])/(x[k]-x[k+1]);
            b[k] = (x[k]*yr[k]-yl[k]*x[k+1])/(x[k]-x[k+1]);
            c[k] = p*k-0.5*a[k]*x[k]*x[k]-b[k]*x[k];
        }
        // free memory
        delete [] yl;
        delete [] yr;
        // Compute constants for tails
        if(is_zero_derivative){
            left_s = x[1];
            left_b = p;
            left_a = (pdf(x[1])*x[1]-2.0*p)/(x[1]*p);
        }
        right_s = x[tbl_size-1];
        right_b = p;
        right_a = pdf(x[tbl_size-1])*tbl_size;
        //Compute maximum value
        max_val = 1.0;
        for(unsigned int k = (!is_zero_derivative ? 0 : 1); k < tbl_size-1; k++){
            std::vector<double> x_cands = med_max(a[k], b[k]);
            double pg;
            for(double val : x_cands){
                if(( val >= x[k]) && ( val <= x[k+1])){
                    pg = pdf(val);
                    //max_val = std::max(max_val, pg/(a[k]*val+b[k]));
                }
            }
            pg = pdf(x[k]);
            max_val = std::max(max_val, pg/(a[k]*x[k]+b[k]));
            pg = pdf(x[k+1]);
            max_val = std::max(max_val, pg/(a[k]*x[k+1]+b[k]));
        }
        //left tail
        if(is_zero_derivative){
            std::vector<double> x_cands = left_max(left_a, left_b, left_s);
            // pdf of aproximation distribution in 0ed interval
            auto appr = [=, this](double x){
                return left_b*exp((x-left_s)*left_a)*x*(left_a*x+2.0)/(left_s*left_s);
            };
            for(double val : x_cands){
                double pg = pdf(val);
                if((val > 0.0)&& (val <= x[1])){
                    max_val = std::max(max_val, pg/appr(val));
                }
            }
        }
        // right tail
        std::vector<double> x_cands = right_max(right_a, right_b, right_s);
        for(double val : x_cands){
            if(val > x[tbl_size-1]){
                double pg = pdf(val);
                max_val = std::max(max_val, pg/(right_b*right_a*exp(right_a*right_s-right_a*val)));
            }
        }
        // Init generator
        return lad_generator.set_state(gs);
    }

    public:
    const unsigned int tbl_size = 256;
    const unsigned int tbl_shift = 56;//64-log2(tbl_size)

    virtual double pdf(double x) = 0;
    virtual double qf(double x) = 0;

    ~LAD()
    {
        delete [] a;
        delete [] b;
        delete [] c;
        delete [] x;
        /*
        a = nullptr;
        b = nullptr;
        c = nullptr;
        x = nullptr;
        //*/
    }
};

}