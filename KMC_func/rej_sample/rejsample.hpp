#ifndef REJSAMPLE_HPP_
#define REJSAMPLE_HPP_

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>
#include <string>
#include <cstdio>
#include <variant>
#include <filesystem>
#include <omp.h>
#include <random>
#include <fstream>
#include <functional>
#include <chrono>
#include <thread>

#include <sys/stat.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/tools/minima.hpp>

#include "auxiliary.hpp"

class RejSample {
    protected:
        double exp_fact_;
        double rest_length_;
        double length_scale_; 

    public: 
        double filterThres; 
        double tkindex; 
        double ttindex; 

    private:
        static constexpr double small_ = 1e-6;
        static constexpr double thres_epsilon_ = 1e-10; 

    public: 
        RejSample() = default;
        ~RejSample() = default;

        RejSample(double alpha, double freelength, double D, double thres, int kinetic_index, int thres_index = 0) {
            length_scale_ = D; 
            exp_fact_ = alpha * length_scale_ * length_scale_;
            rest_length_ = freelength / length_scale_;
            filterThres = thres; 
            tkindex = kinetic_index;
            ttindex = thres_index; 
        }

        inline double getUpperBound() const {
            return sqrt(-log(small_) / exp_fact_) + rest_length_;
        }

        // P_l distribution
        inline double approxPDF(const double r_perp, const double s, const double epsilon) {
            if (s <= 0) return 0; 
            else {
                const double exponent = sqrt(s * s + r_perp * r_perp) - rest_length_;
                const double res = exp(-exp_fact_ * exponent * exponent);
                return length_scale_ * res * epsilon ; 
            }
        }

        // binding volume calculation
        inline double approxBindVolume(const double sbound) {
            assert(sbound > 0);
            auto integrand = [&] (double s) {
                const double exponent = s - rest_length_;
                return s * s * exp(-exp_fact_ * exponent * exponent);
            };
            double error = 0;
            double result = boost::math::quadrature::gauss_kronrod<double, 21>::integrate(integrand, 0, sbound, 10, 1e-6, &error);
            return CUBE(length_scale_) * 4. * M_PI * result; 
        }

        inline double constFunction(const double threshold) {
            return threshold; 
        }

        inline std::tuple<double, double, std::vector<double>> evaluate_pdf(const double r_perp, double a, double b, const double epsilon) {
            // calcuate integral (consider [approxPDF] as 1D function with [r_perp] fixed)
            double error = 0; 
            auto inner_integral = [&](double s) {
                return approxPDF(r_perp, s, epsilon);
            };

            auto neg_inner_integral = [&](double s) {
                return -1 * approxPDF(r_perp, s, epsilon);
            }; 

            double integral = boost::math::quadrature::gauss_kronrod<double, 21>::integrate(inner_integral, a, b, 10, 1e-6, &error);

            // calculate maximum pointx
            const int double_bits = std::numeric_limits<double>::digits;

            auto result = boost::math::tools::brent_find_minima(neg_inner_integral, a, b, double_bits);
            double maxYval = -1 * result.second; // since we use minimum finding algorithm
            double maxXval = result.first;

            // calculate first & last index that beyond certain threshold
            // [just brute force]
            double threshold = 1e-4;
            double grid = 1e-5; 
            std::vector<double> lowupBounds; 
            for (double cnt = a; cnt <= b; cnt += grid) {
                double pdfVal = inner_integral(cnt);
                // low bound
                if ((lowupBounds.size() == 0) && (pdfVal >= threshold)) lowupBounds.push_back(cnt);
                // up bound
                if ((lowupBounds.size() == 1) && (pdfVal <= threshold)) lowupBounds.push_back(cnt);
            }

            return std::tie(integral, maxYval,lowupBounds);
        }

        inline std::tuple<std::vector<double>, double, double, double> doSampling(const int nsamples, const double r_perp, 
        double total_time = 1e2, double delta_t = 1) {
            std::cout << "====== START SAMPLING ======" << std::endl; 

            std::default_random_engine generator;
            std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);
        
            std::vector<double> samples;

            double bindingAccuracy = 0; 
            double bindUnbindRatio = 0; 

            double ub = getUpperBound(); 
            // domain of our target function
            double a = 0;
            double b = 1e2; // [TODO: if it should be a large number ]
            std::cout << "Upbound: " << ub << std::endl; 

            double bindVolume = approxBindVolume(ub);
            std::cout << "bindVolume: " << bindVolume << std::endl; 

            // physical parameters registeration [TODO: CHECK!]
            double ko = 1; // s-1, one-step bare off rate
            double epsilon = 0.25 * 1e1; // mu m-1, binding site density
            if (tkindex == 0) epsilon = 1; // not for monte carlo simulation setting
            double kos = 1; // s-1, off rate constant
            double Ke = 100; // dimensionless, effective association constant

            double integralVal;
            double calThres;
            std::vector<double> lowupBounds; 
            std::tie(integralVal,calThres,lowupBounds) = evaluate_pdf(r_perp, a, b, epsilon);

            std::cout << "Low Bound: " << lowupBounds[0] << "; " << "Up Bound: " << lowupBounds[1] << std::endl; 

            // a = lowupBounds[0]; 
            // b = lowupBounds[1]; 

            std::cout << "Integral of PDF: " << integralVal << std::endl; 
        
            // threhsold (of uniform distribution as the sampled distribution)
            double useThres; 
            if (ttindex == 0) {
                // normalized
                useThres = (calThres + thres_epsilon_);
                std::cout << "Use Calculated Threshold: " << useThres << std::endl; 
            }
            else if (ttindex == 1) {
                useThres = filterThres;
                std::cout << "Use User Input Threshold: " << useThres << std::endl;
            }

            double uniformVal = (b - a) * useThres; 
            std::cout << "Uniform Distribution Integral: " << uniformVal << std::endl; 


            int totalNum = 0;   
            // normalization using integral
            double normalFactor = 1; 

            // for large-scale resampling 
            if (tkindex == 0) {
                while (samples.size() < nsamples) {
                    double x = a + (b - a) * uniform_distribution(generator);
                    double acceptance_prob = (approxPDF(r_perp, x, epsilon)  / constFunction(useThres)) * normalFactor;
                    // std::cout << acceptance_prob << std::endl; 
                    totalNum += 1; 
                    if (uniform_distribution(generator) < acceptance_prob) {
                        samples.push_back(x);
                    }
                }
                // acceptance probability (can be calculated analytically as well)
                bindingAccuracy = static_cast<double>(nsamples) / totalNum;
            }
            // for monte carlo check
            else if (tkindex == 1) {
                double numSteps = total_time / delta_t; 
                std::cout << "## Total Number of Time Steps: " << numSteps << " ##" << std::endl; 
                // how many samples per timestep 
                double samplesPerStep = nsamples / numSteps;
                std::cout << "** Samples per Step: " << samplesPerStep << " **" << std::endl; 
                // calculate unbinding probability
                double unbindProb = 1 - exp(-1 * ko * delta_t / samplesPerStep); 
                std::cout << "Unbind Probability: " << unbindProb << std::endl; 
                // state index -- 0 for unbind, 1 for bind
                int stateIndex = 0; 

                // iterate each time step
                for (int j = 0; j < numSteps; j++) {
                    // std::cout << j << std::endl; 
                    // try to bind
                    int transitionIndex = 0; 
                    // iterate each sample in step
                    for (int i = 0; i < samplesPerStep; i++) {
                        if (transitionIndex == 0) {
                            double x = a + (b - a) * uniform_distribution(generator);
                            double acceptance_prob; 
                            double acceptance_prob2; 
                            if (stateIndex == 0) {
                                // sample P_l [TODO: normalized? ]
                                acceptance_prob = (approxPDF(r_perp, x, epsilon)  / constFunction(useThres)) * normalFactor;
                                // sample P_a [for future convenience, currently as dimensionless constant]
                                acceptance_prob2 = 1 - exp(-1 * delta_t / samplesPerStep * kos * Ke); 
                            }
                            else {
                                acceptance_prob = unbindProb; 
                                acceptance_prob2 = 1; 
                            }
                            // std::cout << acceptance_prob << "," << acceptance_prob2 << std::endl; 

                            // generate two random number to fulfill two distributions separately
                            // P_{on} = P_a * P_l (pointwisely)
                            if ((uniform_distribution(generator) < acceptance_prob) && (uniform_distribution(generator) < acceptance_prob2)) {
                                // successfully (un)bind
                                transitionIndex = 1;  
                                stateIndex = std::abs(stateIndex - 1); 
                                samples.push_back(stateIndex);
                            }
                            // unsuccessfully (un)bind
                            else samples.push_back(stateIndex); 
                        }
                        // if already transition, push current state
                        else samples.push_back(stateIndex); 
                    }
                }

                // calculate P_{on}^T/P_{off}^T
                bindUnbindRatio = Ke * integralVal;
                std::cout << "bindUnbindRatio: " << bindUnbindRatio << std::endl; 
            }
            std::cout << "====== END ======" << std::endl; 


            return std::tie(samples, integralVal, bindingAccuracy, bindUnbindRatio); 
        }

};

#endif