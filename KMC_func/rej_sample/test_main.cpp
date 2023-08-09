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
    // clear
    deleteFilesInDirectory("./data/");
    

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
    try {
        std::filesystem::remove(searchfilename);
    }
    catch (...) {}
    myfile.open(searchfilename); 

    // construct Rejection Sampling constructor object
    // 0 -- generate resampled histogram
    // 1 -- monte carlo simulation
    int tkindex = 1; 
    RejSample rej1(alpha, freelength, D, manualThreshold, tkindex, 0);

    double epsilon = 250; // mu m-1, binding site density
    if (tkindex == 0) epsilon = 1; 

    // record data
    double bdd = rej1.getUpperBound();
    double r_perp = 10 * D;
    for (double s = 0; s < bdd; s += 0.01) {
        double res = rej1.approxPDF(r_perp, s, epsilon); 
        myfile << s << "," << res << std::endl;  
    }

    std::vector<double> sampledData;  
    double normalizationFactor;
    double requiredAcc;
    double bindUnbindRatio; 

    double total_time = 1e4;
    double step_time = 1e-4; 
    // double total_samples = 1e5; 
    std::vector<double> total_samples_vec = pow_linspace(8,8);

    for(int i = 0; i < total_samples_vec.size(); i++) {
        double total_samples = total_samples_vec[i];

        assert((total_samples / (total_time/step_time)) >= 1); 

        std::string samplestep = std::to_string(total_samples/total_time); 

        std::string searchfilename2 = parentRoot + "D=" + strD + "-" + "alpha=" + strAlpha + "-" + "fl=" + strFreeLength + "-" + "st=" + samplestep + "samp.txt"; 
        try {
            std::filesystem::remove(searchfilename2);
        }
        catch(...) {}
        myfile2.open(searchfilename2);

        std::tie(sampledData, normalizationFactor, requiredAcc, bindUnbindRatio) = rej1.doSampling(total_samples, r_perp, total_time, step_time);
        std::cout << "Accuracy: " << requiredAcc << std::endl; 
        
        myfile2.rdbuf()->pubsetbuf(0, 0);

        int boundCnt = 0; 
        int unboundCnt = 0; 

        for (int j = 0; j < sampledData.size(); j++) {
            // for large value, normal way of writing into file is too costy
            // myfile2 << alpha << "," << freelength << "," << sampledData[j] << "," << normalizationFactor << "," << total_time << "," << step_time << "," << total_samples << "," << bindUnbindRatio << std::endl;
            if (tkindex == 1) {
                if (sampledData[j] == 1) boundCnt += 1;
                else unboundCnt += 1;
            }
        }

        double timeRatio = static_cast<double> (boundCnt)/unboundCnt; 
        std::cout << "Bound/Unbound Time Ratio: " << timeRatio << std::endl; 
        myfile2 << total_samples << "," << timeRatio << std::endl; 

        myfile2.close(); 

    }
}