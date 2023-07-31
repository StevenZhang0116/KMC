#include <math.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <fstream>
#include <iterator>
#include <array> 
#include <filesystem>

#include "rejsample.hpp"


int main() {
    const double D = 0.024;
    const double alpha = 0.1 / (2 * 0.00411);
    const double freelength = 0.05;
    // "redundent calculation" only for data recording purpose
    const double M = alpha * D * D;
    const double ell0 = freelength / D;

    double manualThreshold = 0.05; 

    // handle output
    std::string parentRoot = "./data/";
    std::ofstream myfile; std::ofstream myfile2; 
    std::string strD = std::to_string(D);
    std::string strAlpha = std::to_string(alpha);
    std::string strFreeLength = std::to_string(freelength);
    std::string searchfilename = parentRoot + "D=" + strD + "-" + "alpha=" + strAlpha + "-" + "fl=" + strFreeLength + "data.txt"; 
    std::string searchfilename2 = parentRoot + "D=" + strD + "-" + "alpha=" + strAlpha + "-" + "fl=" + strFreeLength + "samp.txt"; 
    try {
        std::filesystem::remove(searchfilename);
        std::filesystem::remove(searchfilename2);
    }

    catch (...) {}
    myfile.open(searchfilename); myfile2.open(searchfilename2);

    // construct Rejection Sampling constructor object
    RejSample rej1(alpha, freelength, D, manualThreshold, 1, 0);

    // record data
    double bdd = rej1.getUpperBound();
    double r_perp = 10 * D;
    for (double s = 0; s < bdd; s += 0.01) {
        double res = rej1.approxPDF(r_perp, s); 
        myfile << s << "," << res << std::endl;  
    }

    std::vector<double> sampledData;  
    double normalizationFactor;
    double requiredAcc;
    int total_time = 1e2;
    int step_time = 1; 
    int total_samples = 1e2; 
    std::tie(sampledData, normalizationFactor, requiredAcc) = rej1.doSampling(total_samples, r_perp, total_time, step_time);
    std::cout << "Accuracy: " << requiredAcc << std::endl; 
    
    for (double sample : sampledData) {
        myfile2 << sample << "," << normalizationFactor << "," << total_time << "," << step_time << std::endl;
    }
}