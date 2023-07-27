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

class RejSample {
    protected:
        double exp_fact_;
        double rest_length_;
        double length_scale_; 

    public: 

    private:
        static constexpr double small_ = 1e-6;

    public: 
        RejSample() = default;
        ~RejSample() = default;

        RejSample(double alpha, double freelength, double D) {
            length_scale_ = D; 
            exp_fact_ = alpha * length_scale_ * length_scale_;
            rest_length_ = freelength / length_scale_;
        }

        inline double getUpperBound() const {
            return sqrt(-log(small_) / exp_fact_) + rest_length_;
        }

        inline double approxPDF(const double r_perp, const double s) {
            if (s <= 0) return 0; 
            else {
                const double exponent = sqrt(s * s + r_perp * r_perp) - rest_length_;
                const double res = exp(-exp_fact_ * exponent * exponent);
                return length_scale_ * res; // normalization (?)
            }
        }

        inline double constFunction(const double threshold) {
            return threshold; 
        }

        inline std::pair<double, double> evaluate_target(const double r_perp, double a, double b, int num_points) {
            double step_size = (b - a) / num_points;
            double integral = 0.0;
            double maxpt = 0.0; 
            for (int i = 0; i < num_points; ++i) {
                double x0 = a + i * step_size;
                double x1 = a + (i + 1) * step_size;
                // calculate integral
                integral += (approxPDF(r_perp, x0) + approxPDF(r_perp, x1)) * 0.5 * step_size;
                // find maximum point
                if (approxPDF(r_perp, x0) >= maxpt) maxpt = approxPDF(r_perp, x0); 
            }
            return std::make_pair(integral, maxpt);
        }

        inline std::vector<double> doSampling(const int num_samples, const double r_perp) {
            std::default_random_engine generator;
            std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

            std::vector<double> samples;

            double ub = getUpperBound(); 
            double a = 0;
            double b = ub; 
            auto res = evaluate_target(r_perp, a, b, 1e7);

            double thres = 0.1;
            double constIntegral = thres * (b - a);


            while (samples.size() < num_samples) {
                double x = a + (b - a) * uniform_distribution(generator);
                double acceptance_prob = (approxPDF(r_perp, x) * constIntegral)  / (constFunction(thres) * res.first);
                if (uniform_distribution(generator) < acceptance_prob) {
                    samples.push_back(x);
                }
            }
            return samples; 
        }

};

#endif