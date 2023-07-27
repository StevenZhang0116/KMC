/**
 * @file cheb.hpp
 * @author Zihan Zhang
 * @brief Create single Baobzi object
 *
 * @copyright Copyright (c) 2023
 *
 */

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
        double param[10];
        baobzi_input_t input;
        int indicator; 
        int printOrNot; 
        baobzi::Function<2, 10, 0, double> saveFunc;
    public:
        Cheb() = default;
        ~Cheb() = default;

        /**
         * @brief constructor of Baobzi object from scratch
         * @param[in]: hl: half length of search domain (2D)
         * @param[in]: cen: center of search domain (2D)
         * @param[in]: tol: maximum desired relative error between your function and the approximant
         * @param[in]: runind: which approximation function to be executed (in conApproxFunc())
         * @param[in]: porn: whether to have output log
         * -- other params' descriptions are included in inline remarks
         */
        Cheb(double (&hl)[2], double (&cen)[2], double tol, double alpha, double freelength, 
        double D, const int runind, const int porn, const double errtol = 1e-3, 
        const double upperbound = 40, const double smallbound = 1e-10, const double e_fact = 0, 
        const double fdep_length = 0, const double M1 = 0, const double M2 = 0, const double theconstant = 0) {
            // copy parameters
            memcpy(&half_length, &hl, sizeof(hl)); 
            // half length, <= 1, to satisfy Baobzi requirement
            assert(half_length[0] <= 1); 
            assert(half_length[1] <= 1); 
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

            // except tolerance, these parameters could be set as fixed and slight tuning
            // WILL NOT give significant effect to the result; 
            // since written as fixed, multiple places in [cheb.hpp] should be 
            // changed altogether, otherwise errors are triggered. 
            input.dim = 2; 
            input.output_dim = 1; 
            input.order = 10;
            input.tol = tol; 
            input.minimum_leaf_fraction = 0.0;
            input.split_multi_eval = 1;
            input.min_depth = 0;
            input.max_depth = 40; 

            // load external parameters to pass into Baobzi object
            // cross-used in multiple subfunctions in [conApproxFunc()]
            // 
            // param[0] = M: exponential constant factor -> exp_fact_
            // param[1] = ell0: protein rest length -> rest_length_
            // param[2] = D; dimensional diameter of rod crosslink is binding to  -> length_scale_
            // param[3] = upperbound; the cutoff radius for r⊥ and s
            // param[4] = theconstant; some arbitrarily defined constant for approximation
            // param[5] = errtol; error tolerance (threshold) in integral value apprximation
            // param[6] = e_fact; energy(load) sensitivity to unbinding
            // param[7] = fdep_length; characteristic length for force dependent unbinding
            // param[8] = M1; exponential constant factor when spring is compressed
            // param[9] = M2; exponential constant factor when spring is stretched
            double M = alpha * D * D; 
            double ell0 = freelength / D; 
            param[0] = M; 
            param[1] = ell0; 
            param[2] = D; 
            param[3] = upperbound; 
            param[4] = theconstant; 
            param[5] = errtol; 
            param[6] = e_fact; 
            param[7] = fdep_length; 
            param[8] = M1; 
            param[9] = M2; 
            param[10] = smallbound; 

            input.data = &param; 

            // determines which function approximation is implemented (e.g. lookup, reverse lookup)
            indicator = runind; 
        }


        /**
         * @brief reconstructed constructor of Baobzi object using precalculated approximated function
         * and other presaved parameters
         */
        Cheb(double (&hl)[2], double (&cen)[2], std::string fullFilePath, const int porn) {
            // define function object -- the dimensions should match with the fixed parameters 
            // declared above, otherwise runtime_error is generated
            baobzi::Function<2, 10, 0, double> readInFunc; 
            // load Baobzi function object
            try{
                baobzi::Function<2, 10, 0, double> tempFunc(fullFilePath.c_str());
                readInFunc = tempFunc; 
            }
            // unmatched dimension of Function object
            catch (const std::runtime_error& ex) {
                std::cout << "Caught std::runtime_error: " << ex.what() << std::endl;
            }
            // save to member variable
            saveFunc = readInFunc; 
            double spaceTaken = saveFunc.memory_usage();

            memcpy(&half_length, &hl, sizeof(hl)); 
            // half length, <= 1, to satisfy Baobzi requirement
            assert(half_length[0] <= 1); 
            assert(half_length[1] <= 1); 
            memcpy(&center, &cen, sizeof(cen)); 
            printOrNot = porn; 
        }

        /**
         *  @brief approximate functions (subject to choose at different scanerios)
         * -- input given in Cheb constructor
         * -- for whole list of explaining the integral functions, refer KMC_func/integrals.hpp
         * -- provides API for users to add their own energy/torque/angle-dependent functions.
         *  @return integral value
         */
        inline std::function<void(const double*, double*, const void*)> conApproxFunc() {
            // x[0] = r⊥/lm: perpendicular distance above rod
            // x[1] = s: upper limit of integral
            /* Normal (Energy Dependent) Integral CDF - 1 */
            auto approxDefaultCDF = [](const double *x, double *y, const void *data) {
                // check upperbound limit
                if (x[1] <= 0) *y = 0; 
                else {
                    // unpack the external parameters
                    const double M = ((double *) data)[0];
                    const double ell0 = ((double *) data)[1];
                    const double D = ((double *) data)[2]; 

                    auto integrand = [&](double s) {
                        /* if want to normalize r⊥, /D after each *x[0] */
                        const double exponent = sqrt(s * s + x[0] * x[0]) - ell0;
                        return exp(-M * exponent * exponent);
                    };

                    double error = 0; 
                    // normalization with param
                    *y = D * boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
                            integrand, 0, x[1] / D, 10, 1e-6, &error); 
                }
            };

            // x[0] = r⊥/lm: perpendicular distance above rod
            // x[1] = s: upper limit of integral
            /* Force Dependent Integral CDF - 6 */
            auto approxFdepCDF = [](const double *x, double *y, const void *data) {
                if (x[1] <= 0) *y = 0;
                else {
                    const double M = ((double *) data)[0];
                    const double ell0 = ((double *) data)[1];
                    const double D = ((double *) data)[2]; 
                    const double e_fact = ((double *) data)[6];
                    const double fdep_length = ((double *) data)[7]; 

                    auto integrand = [&](double s) {
                        const double rprime = sqrt(s * s + x[0] * x[0]) - ell0;
                        const double energy_term = 0.5 * (1.0 - e_fact) * rprime * rprime;
                        const double force_term = fdep_length * rprime;
                        return exp(-M * (energy_term - force_term));
                    };

                    double error = 0; 
                    *y = D * boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
                            integrand, 0, x[1] / D, 10, 1e-6, &error); 

                }
            };

            // x[0] = r⊥/lm: perpendicular distance above rod
            // x[1] = s: upper limit of integral
            /* Assymetric Integral CDF - 7 */
            auto approxAsymCDF = [](const double *x, double *y, const void *data) {
                if (x[1] <= 0) *y = 0;
                else {
                    const double ell0 = ((double *) data)[1];
                    const double D = ((double *) data)[2]; 
                    const double e_fact = ((double *) data)[6];
                    const double fdep_length = ((double *) data)[7];
                    const double M1 = ((double *) data)[8];
                    const double M2 = ((double *) data)[9];

                    auto integrand = [&](double s) {
                        const double rprime = sqrt(s * s + x[0] * x[0]) - ell0;
                        const double energy_term = 0.5 * (1.0 - e_fact) * rprime * rprime;
                        const double force_term = fdep_length * rprime;
                        const double M = rprime < 0. ? M1 : M2; 
                        return exp(-M * (energy_term - force_term));
                    };

                    double error = 0; 
                    *y = D * boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
                            integrand, 0, x[1] / D, 10, 1e-6, &error); 

                }
            };

            // x[0] = r⊥/lm: perpendicular distance above rod
            // x[1] = s: upper limit of integral
            /* Second Order Energy Dependent Integral CDF - 8 */
            auto approxEdep2ndCDF = [](const double *x, double *y, const void *data) {
                if (x[1] <= 0) *y = 0;
                else {
                    const double M = ((double *) data)[0];
                    const double ell0 = ((double *) data)[1];
                    const double D = ((double *) data)[2]; 

                    // what if x[0] < 1? 
                    auto integrand = [&](double s) {
                        const double exponent = ((x[0] - 1.0) * sqrt(1.0 + (s * s / (x[0] * x[0])))) - ell0; 
                        return exp(-M * exponent * exponent); 
                    };

                    double error = 0; 
                    *y = D * boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
                            integrand, 0, x[1] / D, 10, 1e-6, &error); 

                }
            };

            // x[0] = r⊥/lm: perpendicular distance above rod
            // x[1] = s': binding place on the rodll
            /* Normal PDf - 5 */
            auto approxPDF = [](const double *x, double *y, const void *data) {
                if (x[1] <= 0) *y = 0; 
                else {
                    const double M = ((double *) data)[0];
                    const double ell0 = ((double *) data)[1];
                    const double D = ((double *) data)[2];

                    const double threshold = 0.1; 

                    // normal integrand
                    const double exponent = sqrt(x[1] / D * x[1] / D + x[0] * x[0]) - ell0;
                    const double res = exp(-M * exponent * exponent); 

                    *y = D * res; 

                    /* this way would lead to singularity/non-derivative and is detrimental in approximation*/
                    // if (res > threshold) *y = res;
                    // else *y = 0; 
                }
            };

            /* Regardless of inputs, approximate the constant value function based on the input parameter 
             * Used in saving computation cost through approximation in reverse lookup */
            /* Constant Function Approximation - 4 */
            auto approxConstantFunction = [](const double *x, double *y, const void *data) {
                if (x[1] <= 0) *y = 0; 
                else {
                    const double theconstant = ((double*)data)[4]; 
                    *y = theconstant; 
                }
            }; 

            // x[0] = r⊥/lm: perpendicular distance above rod
            // x[1] = integral value 
            /* Reverse Lookup of Regular Integral CDF - 3 */
            auto reverApproxCDF = [](const double* x, double* y, const void* data) {
                if (x[1] <= 0)
                    *y = 0;
                else {
                    // unpack the external parameters
                    const double M = ((double*)data)[0];
                    const double ell0 = ((double*)data)[1];
                    const double D = ((double*)data)[2];
                    const double ub = ((double*)data)[3];
                    const double errtol = ((double*)data)[5]; 
                    const double sb = ((double*)data)[10];

                    double error = 0;
                    // define search range
                    double lower_bound = sb; 
                    double upper_bound = ub * D; 
                    // maximum number of iterations
                    boost::uintmax_t max_iter = 10000; 
                    // desired tolerance
                    boost::math::tools::eps_tolerance<double> tolerance(30); 

                    auto integrand = [&](double s) {
                        const double exponent = sqrt(s * s + x[0] * x[0]) - ell0;
                        return exp(-M * exponent * exponent);
                    };

                    // calculate residue towards the actual integral value (as input)
                    auto solve_func = [&](double caluplimit) { 
                        double residue; 
                        // if under certain tolerance, consider integral value (x[1]) as 0. 
                        if (ABS(x[1]) <= errtol) {
                            residue = D * boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
                                integrand, 0, caluplimit / D, 15, 1e-6, &error); 
                        }

                        else {
                            residue = D * boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
                                integrand, 0, caluplimit / D, 15, 1e-6, &error) - x[1]; 
                        }
                        
                        return residue;
                    }; 

                    try {
                        std::pair<double,double> res = boost::math::tools::bisect(solve_func, lower_bound, upper_bound, tolerance, max_iter);
                        *y = res.first;
                    } 
                    catch(...) {
                        int did = 1;   
                        if (did == 0) *y = 0; 
                        else if (did == 1) {
                            /* will significant decrease the performance */
                            for (double i = lower_bound; i < upper_bound; i += 1e-5) {
                                double res = solve_func(i); 
                                if (ABS(res) < errtol) {
                                    *y = i; 
                                    break; 
                                }
                            }
                        }
                        else throw std::invalid_argument("Root Seeking Error!");
                    }
                    
                }
            };

            // x[0] = s: upper limit of integral
            // x[1] = M (manipulated, not public parameter of class)
            /* Binding Volume Approximation in the normal CDF scanerio - 2 */
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
            
            if (indicator == 1) return approxDefaultCDF;
            else if (indicator == 2) return approxBindV; 
            else if (indicator == 3) return reverApproxCDF; 
            else if (indicator == 4) return approxConstantFunction; 
            else if (indicator == 5) return approxPDF; 
            else if (indicator == 6) return approxFdepCDF; 
            else if (indicator == 7) return approxAsymCDF; 
            else if (indicator == 8) return approxEdep2ndCDF; 
            else throw std::invalid_argument("Invalid choice -> Parameter Setting Error");
        }
        
        /** 
         * @brief approximate function and save it to private parameter for further loading (time efficiency)
         * @return space taken (in bytes) of this constructed function
        */
        inline size_t approxFunc() {
            baobzi::Function<2, 10, 0, double> func_approx(&input, center, half_length, conApproxFunc(), {});
            double spaceTaken = func_approx.memory_usage();
            saveFunc = func_approx; 
            return spaceTaken;
        }

        /**
         * @brief evaluate function value at given input coordinates
         * -- notice that segmentation error would be generated if [inval[]] (after transformation, 
         * like *D or /D, not input here) is not in the predetermined domain of [half_length] and 
         * [center], so normalization of fitting is needed. 
         * @param[in]: inval: 2D coordinates of evaluated point
         * @return apprxoximated value
        */
        inline double evalFunc(double inval[]) {
            if (printOrNot == 1) std::cout << "Evaluate Point at (" << inval[0] << "," << inval[1] << ")" << std::endl;
            // result 
            double res; 
            saveFunc(inval, &res); 
            // to avoid unreasonable approximation (both in normal and reverse lookup)
            return ABS(res); 
        }

        /** 
         * @brief check whether the given 2D coordinate is included in the domain
         * @param[in]: pt: input data (possibly dimensional or dimensionless based on construction) 
         * @param[in]: checkErrorIndex: whether to print log
         * @return binary indicator
         */
        inline int checkInclude(double (&pt)[2], const int checkErrorIndex) {
            try {
                if (checkErrorIndex == 1) {
                    std::cout << "==========" << std::endl; 
                    std::cout << center[0] << "," << half_length[0] << "," << pt[0] << std::endl; 
                    std::cout << center[1] << "," << half_length[1] << "," << pt[1] << std::endl; 
                    std::cout << "==========" << std::endl; 
                }
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

        /**
         * @brief save function(s) and auxiliary parameters to external file(s) that subject to be reloaded
         * using Cheb constructor  
         * @param[in]: folderName: output folder
         * @param[in]: objectIndex: iterative index to count current Baobzi object inside the whole family
         * @return void
         */
        inline void saveFunctionObject(const char* folderName, const int objectIndex) {
            // save Baobzi object
            std::string sL = std::string(folderName) + std::to_string(objectIndex);
            const char* saveLocation = sL.c_str();  
            int result = std::remove(saveLocation);
            if (result == 0) std::cout << "Baobzi Filename Detected and Deleted: " << saveLocation << std::endl; 
            double spaceTaken = saveFunc.memory_usage();
            saveFunc.save(saveLocation); 
            // save auxiliary information that related to tis Baobzi object
            std::ofstream myfile; 
            std::string fileSuffix = "-Res"; 
            std::string sLAux = sL + fileSuffix; 
            const char* saveLocation2 = sLAux.c_str(); 
            int result2 = std::remove(saveLocation2);
            if (result2 == 0) std::cout << "Baobzi Auxiliary File Detected and Deleted: " << saveLocation2 << std::endl; 
            myfile.open(saveLocation2);
            // ---- saving format ----
            // [center of 1D] -- [half length of 1D]
            // [center of 2D] -- [half length of 2Ds]
            myfile << center[0] << "," << half_length[0] << std::endl; 
            myfile << center[1] << "," << half_length[1] << std::endl; 
            // ---- saving format ----
            myfile.close(); 
        }
};


#endif