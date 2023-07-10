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
        double param[4];
        baobzi_input_t input;
        const char* oname; 
        int indicator; 
        int printOrNot; 
    private:
        baobzi::Function<2,10,0,double> savefunc;
    public:
        Cheb() = default;
        ~Cheb() = default;

        /**
         * @brief constructor of Baobzi object
         */
        Cheb(double (&hl)[2], double (&cen)[2], double tol, double alpha, double freelength, 
        double D, const char* output_name, const int runind, const int porn, const double upperbound = 0, const double aconst = 0) {
            // std::cout << "Construct Baobzi Object" << std::endl;
            memcpy(&half_length, &hl, sizeof(hl)); 
            assert(half_length[0] <= 1); assert(half_length[1] <= 1); // half length, <= 1
            memcpy(&center, &cen, sizeof(cen)); 
            // whether to log output
            printOrNot = porn; 
            // in case of out of range
            if (printOrNot == 1) {
                speak("half-length 1", half_length[0]); 
                speak("half-length 2", half_length[1]); 
                speak("center 1", center[0]);
                speak("center 2", center[1]);
            } 
            // except tolerance, these parameters could be set as fixed and slight tuningWILL NOT 
            // give significant effect to the result
            input.dim = 2; 
            input.output_dim = 1; 
            input.order = 10;
            input.tol = tol; 
            input.minimum_leaf_fraction = 0.0;
            input.split_multi_eval = 1;
            input.min_depth = 0;
            input.max_depth = 40; 
            // load external parameters
            // param[0] = M: exponential constant factor -> exp_fact_
            // param[1] = ell0: protein rest length -> rest_length_
            // param[2] = D; diameter of rod crosslink is binding to  -> length_scale_
            // param[3] = upperbound; the cutoff radius for r⊥ and s
            // param[4] = constant; some arbitrarily defined constant
            double M = alpha * D * D; 
            double ell0 = freelength / D; 
            param[0] = M; param[1] = ell0; param[2] = D; param[3] = upperbound; param[4] = aconst; 
            input.data = &param; 
            // determines which function approximation is implemented (e.g. lookup, reverse lookup)
            indicator = runind; 
            oname = output_name; 
        }

        /**
         *  @brief approximate functions (subject to choose at different scanerios)
         *  1) carrier line cumulative distribution function 
         *  2) binding volume function
         *          
         *  @return integral value
         */
        inline std::function<void(const double*, double*, const void*)> conApproxFunc() {
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
                        /* if want to normalize r⊥, /D after each *x[0] */
                        const double exponent = sqrt(s * s + x[0] * x[0]) - ell0;
                        return exp(-M * exponent * exponent);
                    };

                    *y = D * boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
                            integrand, 0, x[1] / D, 10, 1e-6, &error); // normalization with param
                }
            };

            /* Regardless of inputs, approximate the constant value function based on the input parameter 
             * Used in saving computation cost through approximation in reverse lookup
             */
            auto approxConstantFunction = [](const double *x, double *y, const void *data) {
                if (x[1] <= 0) *y = 0; 
                else {
                    const double cst = ((double*)data)[4]; 
                    *y = cst; 
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
                    const double ub = ((double*)data)[3];

                    double error = 0;
                    double shift = 1e-20; 
                    double errortolerence = 1e-6; 

                    double lower_bound = 0.0 + shift; 
                    double upper_bound = 2; 
                    boost::uintmax_t max_iter = 1000; // Maximum number of iterations
                    boost::math::tools::eps_tolerance<double> tolerance(30); // Desired tolerance

                    auto integrand = [&](double s) {
                        const double exponent = sqrt(s * s + x[0] * x[0]) - ell0;
                        return exp(-M * exponent * exponent);
                    };

                    auto solve_func = [&](double caluplimit) { 
                        double residue = D * boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
                            integrand, 0, caluplimit / D, 10, 1e-6, &error) - x[1]; 
                        return residue;
                    }; 

                    try {
                        std::pair<double,double> res = boost::math::tools::bisect(solve_func, lower_bound, upper_bound, tolerance, max_iter);
                        *y = res.first;
                    } 
                    catch(...) {
                        // *y = 0; 
                        for (double i = lower_bound; i < upper_bound; i+=0.0001) {
                            double res = solve_func(i); 
                            if (ABS(res) < errortolerence) {
                                *y = i; 
                                break; 
                            }
                        }
                        // std::cout << "a" << std::endl; 
                    }
                    
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
            
            if (indicator == 1) return approxCDF;
            else if (indicator == 2) return approxBindV; 
            else if (indicator == 3) return reverApproxCDF; 
            else if (indicator == 4) return approxConstantFunction; 
            else throw std::invalid_argument("Invalid choice -> Parameter Setting Error");
        }
        
        /** 
         * @brief approximate function and save it to private parameter for further loading (time efficiency)
         * 
         * @return null
        */
        inline double approxFunc() {
            if (printOrNot == 1) std::cout << "Generate Function Approximator" << std::endl; 
            baobzi::Function<2,10,0,double> func_approx(&input, center, half_length, conApproxFunc(), {});
            double spaceTaken = func_approx.print_stats(printOrNot);
            savefunc = func_approx; 
            return spaceTaken;
        }

        /**
         * @brief evaluate function value at given input coordinates; notice that segmentation error would be generated if `inval[]`
         *        (after transformation, like *D or /D, not input here) is not in the predetermined domain of `half_length` and 
         *        `center`, so normalization of fitting is needed. 
         * 
         * @return apprxoximated value
        */
        inline double evalFunc(double inval[]) {
            if (printOrNot == 1) std::cout << "Evaluate Point at (" << inval[0] << "," << inval[1] << ")" << std::endl; 
            double res; 
            savefunc(inval, &res); 
            return res; 
        }

        /** 
         * @brief check whether the given 2D coordinate is included in the domain
         * 
         * @return binary value
         */
        inline int checkInclude(double (&pt)[2]) {
            try {
                // std::cout << center[0] << "," << half_length[0] << "," << pt[0] << std::endl; 
                // std::cout << center[1] << "," << half_length[1] << "," << pt[1] << std::endl; 
                if ((pt[0] >= center[0] - half_length[0]) && (pt[0] <= center[0] + half_length[0]) 
                && (pt[1] >= center[1] - half_length[1]) && (pt[1] <= center[1] + half_length[1])) {
                    return 1; 
                }
                else throw std::invalid_argument("NOT IN THE DOMAIN");
            }
            catch(...) {
                return 0; 
            }
        }
};


#endif