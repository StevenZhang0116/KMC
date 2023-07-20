/**
 * @file rejection_sampling.cpp
 * @author Zihan Zhang
 * @brief Simple Demo code to test rejection sampling
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <iostream>
#include <random>
#include <cmath>

double normalPdf(double x, double mean, double variance) {
    return std::exp(-(x - mean) * (x - mean) / (2 * variance)) / std::sqrt(2 * M_PI * variance);
}

int main() {
    // Seed the random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Define the parameters of the target normal distribution
    double target_mean = 0.0;
    double target_variance = 1.0;

    // Define the parameters of the proposal exponential distribution
    double proposal_rate = 1.0;

    // Create distributions for the proposal and target distributions
    std::exponential_distribution<double> proposal_dist(proposal_rate);

    // Maximum value of the target distribution for scaling
    double max_target = normalPdf(target_mean, target_mean, target_variance);

    // Number of samples to generate
    int num_samples = 10000;

    // Perform rejection sampling
    int accepted_samples = 0;
    for (int i = 0; i < num_samples; ++i) {
        // Sample from the proposal distribution
        double x = proposal_dist(gen);

        // Compute the acceptance probability
        double acceptance_prob = normalPdf(x, target_mean, target_variance) / (proposal_rate * max_target);

        // Generate a uniform random number between 0 and 1
        double u = std::uniform_real_distribution<double>(0, 1)(gen);

        // Accept or reject the sample
        if (u <= acceptance_prob) {
            // Accept the sample
            accepted_samples++;
            std::cout << x << std::endl;
        }
    }

    // Calculate the acceptance rate
    double acceptance_rate = static_cast<double>(accepted_samples) / num_samples;
    std::cout << "Acceptance rate: " << acceptance_rate << std::endl;

    return 0;
}
