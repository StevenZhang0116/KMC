/**
 * @file rejection_sampling_2d.cpp
 * @author Zihan Zhang
 * @brief Simple Demo code to test rejection sampling in 2D
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <vector>
#include <filesystem>

// Define the target 2D normal distribution
double target_distribution(double x, double y, double mean_x, double mean_y, double var_x, double var_y) {
    double exponent = -0.5 * ((pow((x - mean_x) / var_x, 2)) + (pow((y - mean_y) / var_y, 2)));
    double normalization = 1.0 / (2.0 * M_PI * var_x * var_y);
    return normalization * exp(exponent);
}

// Define the proposal uniform distribution
double proposal_distribution(double x, double y, double x_min, double x_max, double y_min, double y_max) {
    return 1.0; // Uniform density within the specified rectangle
}

int main() {
    int num_samples = 1000000;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

    std::vector<double> samples_x;
    std::vector<double> samples_y;

    double mean_x = 0.0;
    double mean_y = 0.0;
    double var_x = 2.0;
    double var_y = 2.0;

    double x_min = -10.0; 
    double x_max = 10.0;  
    double y_min = -10.0; 
    double y_max = 10.0;  

    // Rejection sampling with uniform proposal distribution
    while (samples_x.size() < num_samples && samples_y.size() < num_samples) {
        // Sample from the proposal distribution (uniformly within the specified rectangle)
        double x = x_min + (x_max - x_min) * uniform_distribution(generator);
        double y = y_min + (y_max - y_min) * uniform_distribution(generator);

        // Calculate the acceptance probability
        double acceptance_prob = target_distribution(x, y, mean_x, mean_y, var_x, var_y) /
                                 proposal_distribution(x, y, x_min, x_max, y_min, y_max);

        // Accept the sample with a probability of acceptance_prob
        if (uniform_distribution(generator) < acceptance_prob) {
            samples_x.push_back(x);
            samples_y.push_back(y);
        }
    }

    // Save the samples to a file
    std::ofstream myfile;
    std::string searchfilename = "rej_sample_data_2d.txt";
    try {
        std::filesystem::remove(searchfilename);
    } catch (...) {
    }
    myfile.open(searchfilename);

    for (size_t i = 0; i < samples_x.size(); i++) {
        myfile << samples_x[i] << "," << samples_y[i] << std::endl;
    }

    return 0;
}
