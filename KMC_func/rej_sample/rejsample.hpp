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

#include <sys/stat.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/tools/minima.hpp>


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

        inline double approxPDF(const double r_perp, const double s) {
            if (s <= 0) return 0; 
            else {
                const double exponent = sqrt(s * s + r_perp * r_perp) - rest_length_;
                const double res = exp(-exp_fact_ * exponent * exponent);
                return length_scale_ * res; 
            }
        }

        inline double constFunction(const double threshold) {
            return threshold; 
        }

        inline std::pair<double, double> evaluate_target(const double r_perp, double a, double b) {
            // calcuate integral (consider [approxPDF] as 1D function with [r_perp] fixed)
            double error = 0; 
            auto inner_integral = [&](double s) {
                return approxPDF(r_perp, s);
            };

            auto neg_inner_integral = [&](double s) {
                return -1 * approxPDF(r_perp, s);
            }; 

            double integral = boost::math::quadrature::gauss_kronrod<double, 21>::integrate(inner_integral, a, b, 10, 1e-6, &error);

            // calculate maximum pointx
            auto result = boost::math::tools::brent_find_minima(neg_inner_integral, a, b, 10);
            double maxYval = -1 * result.second; // since we use minimum finding algorithm
            double maxXval = result.first;

            return std::make_pair(integral, maxYval);
        }

        inline std::tuple<std::vector<double>, double, double> doSampling(const int nsamples, const double r_perp, int total_time = 1e2, int delta_t = 1) {

            std::default_random_engine generator;
            std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

            std::vector<double> samples;

            double bindingAccuracy; 

            double ub = getUpperBound(); 
            // domain of our target function
            double a = 0;
            double b = ub; 

            auto res = evaluate_target(r_perp, a, b);
            double integralVal = res.first; 
            double calThres = res.second; 

            // threhsold (of uniform distribution as the sampled distribution)
            double useThres; 
            if (ttindex == 0) {
                useThres = (calThres + thres_epsilon_) / integralVal;
                std::cout << "Use Calculated Threshold: " << useThres << std::endl; 
            }
            else if (ttindex == 1) {
                useThres = filterThres;
                std::cout << "Use User Input Threshold: " << useThres << std::endl;
            }

            double uniformVal = (b - a) * useThres; 
            int totalNum = 0;   

            std::cout << "Integral of PDF: " << integralVal << std::endl; 

            // normalization using integral
            double normalFactor = 1 / integralVal; 

            // for large-scale resampling 
            if (tkindex == 0) {
                while (samples.size() < nsamples) {
                    double x = a + (b - a) * uniform_distribution(generator);
                    double acceptance_prob = (approxPDF(r_perp, x)  / constFunction(useThres)) * normalFactor;
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
                int numSteps = total_time / delta_t; 
                // how many samples per timestep 
                int samplesPerStep = nsamples / numSteps;
                // calculate unbinding probability
                double k0 = 0.77; // sec-1, one-step bare off rate
                double unbindProb = 1 - exp(-k0 * delta_t); 
                // state index -- 0 for unbind, 1 for bind
                int stateIndex = 0; 
                // go over each step
                for (int j = 0; j < numSteps; j++) {
                    // std::cout << j << std::endl; 
                    // try to bind
                    int transitionIndex = 0; 
                    for (int i = 0; i < samplesPerStep; i++) {
                        if (transitionIndex == 0) {
                            double x = a + (b - a) * uniform_distribution(generator);
                            double acceptance_prob; 
                            if (stateIndex == 0) {
                                acceptance_prob = (approxPDF(r_perp, x)  / constFunction(useThres));
                            }
                            else {
                                acceptance_prob = unbindProb; 
                            }
                            // std::cout << acceptance_prob << std::endl; 
                            if (uniform_distribution(generator) < acceptance_prob) {
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

            }


            return std::tie(samples, integralVal, bindingAccuracy); 
        }

};

#endif