#ifndef CHEB_HPP_
#define CHEB_HPP_

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>
#include <string>
#include <cstdio>
#include <variant>

#include <boost/math/quadrature/gauss_kronrod.hpp>
// add baobzi dependence
#include <baobzi_template.hpp>
#include "auxil.hpp"
#include "macros.hpp"

class Cheb {
    public:
        double half_length[2];
        double center[2]; 
        double param[3];
        baobzi_input_t input;
        const char* oname; 
    private:
        baobzi::Function<2,10,0,double> savefunc;
        
    public:
        Cheb() = default;
        ~Cheb() = default;

        /**
         * @brief constructor of Baobzi object
         */
        Cheb(double sbound1, double sbound2, double cen1, double cen2, int dim, 
        int output_dim, int order, double tol, double minimum_leaf_fraction, 
        int split_multi_eval, int min_depth, int max_depth, double M, double ell0,
        double D, const char* output_name) {
            // std::cout << "Construct Baobzi Object" << std::endl;
            half_length[0] = sbound1; half_length[1] = sbound2; 
            assert(half_length[0] <= 1); assert(half_length[1] <= 1); // half length, <= 1
            center[0] = cen1; center[1] = cen2; 
            // in case of out of range
            speak("half-length 1", half_length[0]); 
            speak("half-length 2", half_length[1]); 
            speak("center 1", center[0]);
            speak("center 2", center[1]); 

            input.dim = dim; 
            input.output_dim = output_dim; 
            input.order = order;
            input.tol = tol; 
            input.minimum_leaf_fraction = minimum_leaf_fraction;
            input.split_multi_eval = split_multi_eval;
            input.min_depth = min_depth;
            input.max_depth = max_depth; 
            // load external parameters
            // param[0] = M: exponential constant factor
            // param[1] = ell0: protein rest length
            param[0] = M; param[1] = ell0; param[2] = D; 
            input.data = &param; 

            oname = output_name; 
        }
        

        /**
         *  @brief approximate functions (subject to choose at different scanerios)
         *  1) carrier line cumulative distribution function 
         *  2) binding volume function
         *          
         *  @return integral value
         */
        inline std::function<void(const double*, double*, const void*)> conApproxFunc(const int choice) {
            // x[0] = r⊥/lm: perpendicular distance above rod
            // x[1] = s: upper limit of integral
            auto approxCDF = [](const double *x, double *y, const void *data) {
                // check upperbound limit
                if (x[1] <= 0) *y = 0; 
                else {
                    // unpack the external parameters
                    const double M = ((double *) data)[0];
                    const double ell0 = ((double *) data)[1];
                    const double D = ((double *) data)[2]; 
                    double error = 0; 

                    auto integrand = [&](double s) {
                        const double exponent = sqrt(s * s + x[0] * x[0]) - ell0;
                        return exp(-M * exponent * exponent);
                    };

                    *y = D * boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
                            integrand, 0, x[1] / D, 10, 1e-6, &error); // normalization with param
                }
            };

            // x[0] = r⊥/lm: perpendicular distance above rod
            // x[1] = integral value 
            auto reverApproxCDF = [](const double* x, double* y, const void* data) {
                if (x[1] <= 0)
                    *y = 0;
                else {
                    // unpack the external parameters
                    const double M = ((double*)data)[0];
                    const double ell0 = ((double*)data)[1];
                    const double D = ((double*)data)[2];
                    double error = 0;

                    double lowerbound = 0.0 + 1e-10;
                    double upperbound = 0.5;
                    boost::uintmax_t max_iter = 1000; // Maximum number of iterations
                    boost::math::tools::eps_tolerance<double> tolerance(20); // Desired tolerance

                    auto integrand = [&](double s) {
                        const double exponent = sqrt(s * s + x[0] * x[0]) - ell0;
                        return exp(-M * exponent * exponent);
                    };

                    auto solve_func = [&](double caluplimit) { 
                        // speak("caluplimit", caluplimit); 
                        double residue = boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
                            integrand, 0, caluplimit / D, 10, 1e-6, &error) - x[1] / D;
                        // speak("residue", residue); 
                        return residue;
                    }; 

                    std::pair<double,double> res = boost::math::tools::bisect(solve_func, lowerbound, upperbound, tolerance, max_iter);
                    // speak("value=======", res.first); 
                    *y = res.first;
                }
            };


            // x[0] = s: upper limit of integral
            // x[1] = M (manipulated, not public parameter of class)
            auto approxBindV = [](const double *x, double *y, const void *data) {
                if (x[1] <= 0) *y = 0; 
                else {
                    const double M = ((double *) data)[0];
                    const double ell0 = ((double *) data)[1];
                    const double D = ((double *) data)[2]; 
                    double error = 0; 

                    auto integrand = [&](double s) {
                        const double exponent = s - ell0; 
                        return s * s * exp(-x[1] * exponent * exponent);
                    };

                    double result = boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
                                    integrand, 0, x[0] / D, 10, 1e-6, &error); 

                    *y = 4. * M_PI * result;
                }
            };

            if (choice == 1) {
                return approxCDF;
            }
            else if (choice == 2) {
                return approxBindV; 
            }
            else if (choice == 3){
                return reverApproxCDF; 
            }
            else {
                throw std::invalid_argument("Invalid choice");
            }
        }
        
        /** 
         * @brief approximate function and save it to private parameter for further loading (time efficiency)
         * 
         * @return null
        */
        inline void approxFunc(const int choice) {
            std::cout << "Generate Function Approximator" << std::endl; 
            baobzi::Function<2,10,0,double> func_approx(&input, center, half_length, conApproxFunc(choice), {});
            func_approx.print_stats();
            savefunc = func_approx; 
            return;
        }

        /**
         * @brief evaluate function value at given input coordinates; notice that segmentation error would be generated if `inval[]`
         *        (after transformation, like *D or /D, not input here) is not in the predetermined domain of `half_length` and 
         *        `center`, so normalization of fitting is needed. 
         * 
         * @return apprxoximated value
        */
        inline double evalFunc(double inval[]) {
            // std::cout << "Evaluate Point at (" << inval[0] << "," << inval[1] << ")" << std::endl; 
            double res; 
            savefunc(inval, &res); 
            return res; 
        }
};


#endif
