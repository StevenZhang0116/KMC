/**
 * @file rejection_sampling_1d.cpp
 * @author Zihan Zhang
 * @brief Simple Demo code to test rejection sampling in 1D
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <cmath>
#include <filesystem>
#include <boost/math/tools/minima.hpp>

// Define the target normal distribution with a constant added
// TODO: what if I make constant nonzero? which other part of the code I should change? 
double target_distribution(double x) {
    double mean = 0.0;
    double variance = 1.0;
    double constant = 0.0; 
    return exp(-0.5 * ((x - mean) / variance) * ((x - mean) / variance)) / (variance * sqrt(2 * M_PI)) + constant;
}

// Define the proposal uniform distribution
// the value should be large or equal to the maximum of target_distribution
double proposal_distribution2(double x) {
    return 0.5;
}

// Numerical integration using the trapezoidal rule
double integrate_target(double a, double b, int num_points) {
    double step_size = (b - a) / num_points;
    double integral = 0.0;
    for (int i = 0; i < num_points; ++i) {
        double x0 = a + i * step_size;
        double x1 = a + (i + 1) * step_size;
        integral += (target_distribution(x0) + target_distribution(x1)) * 0.5 * step_size;
    }
    return integral;
}

int main() {
    int num_samples = 10000;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

    std::vector<double> samples;

    // Calculate the normalization constant for the target distribution
    double a = -10.0; // Lower bound of the uniform distribution
    double b = 15.0;  // Upper bound of the uniform distribution
    double normalization_constant = integrate_target(a, b, 10000);

    // Rejection sampling with uniform proposal distribution
    while (samples.size() < num_samples) {
        // Sample from the proposal distribution (uniformly within [a, b])
        double x = a + (b - a) * uniform_distribution(generator);

        // Calculate the acceptance probability with constant and normalization constant adjustment
        double acceptance_prob = target_distribution(x) / (proposal_distribution2(x) * normalization_constant);

        // Accept the sample with a probability of acceptance_prob
        if (uniform_distribution(generator) < acceptance_prob) {
            samples.push_back(x);
        }
    }

    // Save the samples to a file
    std::ofstream myfile;
    std::string searchfilename = "rej_sample_data_1d.txt";
    try {
        std::filesystem::remove(searchfilename);
    }
    catch (...) {}
    myfile.open(searchfilename);

    for (double sample : samples) {
        myfile << sample << "," << std::endl;
    }

    return 0;
}
