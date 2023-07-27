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
    const double tol = 1e-2;
    const double D = 0.024;
    const double alpha = 0.1 / (2 * 0.00411);

    const double freelength = 0.05;
    const double M = alpha * D * D;
    const double ell0 = freelength / D;

    // handle output
    std::ofstream myfile; std::ofstream myfile2; 
    std::string strD = std::to_string(D);
    std::string strAlpha = std::to_string(alpha);
    std::string strFreeLength = std::to_string(freelength);
    std::string searchfilename = "D=" + strD + "-" + "alpha=" + strAlpha + "-" + "fl=" + strFreeLength + "data.txt"; 
    std::string searchfilename2 = "D=" + strD + "-" + "alpha=" + strAlpha + "-" + "fl=" + strFreeLength + "samp.txt"; 
    try {
        std::filesystem::remove(searchfilename);
        std::filesystem::remove(searchfilename2);
    }

    catch (...) {}
    myfile.open(searchfilename); myfile2.open(searchfilename2);

    RejSample rej1(alpha, freelength, D);
    double bdd = rej1.getUpperBound();
    double r_perp = 10 * D;
    for (double s = 0; s < bdd; s += 0.01) {
        double res = rej1.approxPDF(r_perp, s); 
        myfile << s << "," << res << std::endl;  
    }

    std::vector<double> sampledData = rej1.doSampling(1e6, r_perp);
    
    for (double sample : sampledData) {
        myfile2 << sample << std::endl;
    }
}