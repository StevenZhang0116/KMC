/**
 * @file mh_1d.cpp
 * @author Zihan Zhang
 * @brief Simple Demo code to test Metropolis-Hastings Algorithm in 1D
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

// Define the target normal distribution
double target_distribution(double x, double mean, double variance) {
    double exponent = -0.5 * pow((x - mean) / variance, 2);
    double normalization = 1.0 / (variance * sqrt(2 * M_PI));
    return normalization * exp(exponent);
}

int main() {
    int num_samples = 10000;
    double mean = 0.0;
    double variance = 1.0;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

    // Initialize the sample with an arbitrary value
    double current_sample = 0.0;

    // Initialize a vector to store the samples
    std::vector<double> samples;

    // Metropolis-Hastings algorithm
    for (int i = 0; i < num_samples; i++) {
        // Propose a new sample from a normal distribution with the current sample as mean
        std::normal_distribution<double> proposal_distribution(current_sample, 0.5);

        double proposed_sample = proposal_distribution(generator);

        // Calculate the acceptance probability
        double acceptance_prob = std::min(1.0, target_distribution(proposed_sample, mean, variance) /
                                              target_distribution(current_sample, mean, variance));

        // Accept the proposed sample with a probability of acceptance_prob
        if (uniform_distribution(generator) < acceptance_prob) {
            current_sample = proposed_sample;
        }

        // Store the accepted sample
        samples.push_back(current_sample);
    }

    std::ofstream myfile;
    std::string searchfilename = "mh_sample_data_1d.txt";
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
