#ifndef CHEBCOLL_HPP_
#define CHEBCOLL_HPP_

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

#include <boost/math/quadrature/gauss_kronrod.hpp>
// add baobzi dependence
#include <baobzi_template.hpp>
#include "auxil.hpp"
#include "macros.hpp"
#include "cheb.hpp"

class Chebcoll {
    protected:
        double exp_fact_;
        double rest_length_;
        double length_scale_; 
        double the_upper_bound_; 
        double grid_size_magnitude_; 
        double integral_min_; 
        double integral_max_; 
    public: 
        double talpha; double tfreelength; 
        double tbbtol; int ri; const char* on; int tpon; 
        std::vector<Cheb> allCheb; 
        std::vector<double> breakPtsCollChange; 
        std::vector<double> breakPtsCollUnchange; 
        std::vector<double> chebSpaceTaken; 
    
    private:
        static constexpr double small_ = 1e-6;
        static constexpr double shift_small_ = 1e-10; 
    
    public:
        Chebcoll() = default;
        ~Chebcoll() = default;

        Chebcoll(double alpha, double freelength, double D, const int runind, double bbtol = 1e-4, 
        std::vector<double> integralMinMax = std::vector<double>(), const char* output_name = "func_approx.baobzi", 
        const int printOrNot = 0) {
            // initialize dimensionless 
            length_scale_ = D; 
            exp_fact_ = alpha * length_scale_ * length_scale_;
            rest_length_ = freelength / length_scale_;
            // (dummy) registeration of parameters
            talpha = alpha; 
            tfreelength = freelength;
            tbbtol = bbtol; 
            on = output_name; 
            ri = runind; 
            tpon = printOrNot; 
            speak("Tolerance of Bobazi Family", bbtol); 
            // only need to set up integral min/max value in the REVERSE LOOKUP 
            if (integralMinMax.size() > 0) {
                assert(runind == 3); 
                integral_min_ = integralMinMax[0];
                integral_max_ = integralMinMax[1]; 
                std::cout << "Setup integral max/min value in Constructor: " << integral_min_ << "," << integral_max_ << std::endl; 
            }
        }

        inline double calcBoltzmann(double dist_cent) const {
            double r = dist_cent / length_scale_;
            return exp(-exp_fact_ * SQR(r - rest_length_));
        }

        double getUpperBound() const {
            return sqrt(-log(small_) / exp_fact_) + rest_length_;
        }

        inline void createBaobziFamily(std::vector<std::vector<double>> tempkk = std::vector<std::vector<double>>()) {
            double upBound; double oneFixLength; 
            double oneFixCenter; double otherGrid; 
            std::vector<double> lengthVec; 
            std::vector<double> centerVec; 
            std::vector<double> gridVec;
            if (ri == 1){
                std::cout << "==== Create Family of Positive Checking ====" << std::endl; 
                upBound = getUpperBound();
                the_upper_bound_ = upBound; 
                speak("UpBound (for both dimensions)", upBound);  
                oneFixLength = upBound / 2 * length_scale_; 
                assert(oneFixLength <= 1); 
                oneFixCenter = oneFixLength; 
                // half-length
                otherGrid = find_order(oneFixLength); 
                while (the_upper_bound_ / otherGrid <= 1e2) {
                    otherGrid /= 10; 
                    std::cout << "=== MAKE GRID SMALLER === " << std::endl; 
                }
            }

            else if (ri == 3) {
                std::cout << "==== Create Family of Reverse Checking ====" << std::endl; 
                upBound = getUpperBound();
                the_upper_bound_ = upBound; 
                speak("UpBound (for one dimensions)", upBound);  
                // speak("UpBound (for another dimension)", rr); 
                otherGrid = 0.1;
                for (int i = 0; i < tempkk.size(); i++) {
                    double rr = tempkk[i][1];
                    double ll = tempkk[i][0];
                    oneFixLength = (rr - ll) / 2;
                    // speak("oneFixLength",oneFixLength); 
                    assert(oneFixLength <= 1); 
                    oneFixCenter = ll + oneFixLength; 
                    // normalization (?)
                    // oneFixLength *= length_scale_; 
                    // oneFixCenter *= length_scale_; 
                    // match order
                    double tGrid = find_order(oneFixLength);
                    tGrid = std::max(tGrid, 0.1); 
                    // should be integer, as division of 10's powers
                    double rrTimes = otherGrid / tGrid; 
                    // speak("rrTimes", rrTimes);
                    int kk = ceil(rrTimes);
                    if (kk > 1) kk = round10(kk); 
                    // speak("kk", kk);
                    for (int j = 0; j < kk; j++) {
                        gridVec.push_back(tGrid);
                        lengthVec.push_back(oneFixLength);
                        centerVec.push_back(oneFixCenter);
                    }
                }
                
                speak("Size of lengthVec", gridVec.size());
                speak("Sum of lengthVec", total_sum(gridVec)); 
                assert(total_sum(gridVec) <= upBound); 
            }

            std::vector<double> iterVec;
            if (ri == 1) {
                // declare bound for the fixed parameter, usually distPerp (vertical distance)
                double lbound = otherGrid; 
                double ubound = upBound - otherGrid; 
                double gg = 1; // range [1,2], inversely proportional to running time
                double gridSize = gg * otherGrid; 

                for (double iter = lbound; iter < ubound; iter += gridSize) iterVec.push_back(iter); 
                // iteratively create Baobzi object
                speak("Baobzi Objects need to be created", floor((ubound - lbound)/gridSize + 1)); 
                speak("Grid Size Magnitude", otherGrid); // magnitude of grid size along both dimension
                grid_size_magnitude_ = otherGrid; 
            }
            else if (ri == 3){
                iterVec = cumulativeSum(gridVec); 
            }

            int tri; 
            /* if want to normalize r⊥, change upper bound of iii to 1 (instead of iterVec.size()) */
            for (size_t iii = 0; iii < iterVec.size(); iii++) {
                double thisCenter = iterVec[iii]; 
                if (ri == 3) speak("thisCenter", thisCenter); 
                const auto st1 = get_wtime();
                if (ri == 3) {
                    oneFixLength = lengthVec[iii];
                    oneFixCenter = centerVec[iii]; 
                    otherGrid = gridVec[iii]; 
                }
                /* if want to normalize r⊥, change the first coordinates of hl[] and center[] 
                 * to the second coordinate */ 
                double hl[2] = {otherGrid, oneFixCenter};
                double center[2] = {thisCenter, oneFixCenter}; 
                if ((hl[1] == 0.0) && (center[1] == 0.0)) {
                    tri = 4; 
                    hl[1] = 1; 
                    center[1] = 1; 
                }
                else {
                    tri = ri; 
                }
                if (ri == 3) {
                    // std::cout << "hl: " << hl[0] << ";" << hl[1] << std::endl;
                    // std::cout << "center: " << center[0] << ";" << center[1] << std::endl;
                }
                Cheb theBaobzi(hl, center, tbbtol, talpha, tfreelength, length_scale_, on, tri, tpon, the_upper_bound_);
                double ssTaken = theBaobzi.approxFunc(); // taken space in Mb
                allCheb.push_back(theBaobzi); 
                breakPtsCollChange.push_back(thisCenter - otherGrid); 
                breakPtsCollUnchange.push_back(oneFixCenter - oneFixLength); 
                chebSpaceTaken.push_back(ssTaken);
                if ((iii % 10 == 0) && (ri == 3)) { // only print log in reverse lookup
                    const auto ft1 = get_wtime();
                    const double dt1 = get_wtime_diff(&st1, &ft1);
                    speak("Baobzi Created", iii); 
                    speak("Needed Time per object (s)", dt1); 
                }
            }
            speak("Total Baobzi Family Space (MiB)", total_sum(chebSpaceTaken)); 
        }

        inline double evalSinglePt(double (&ptCenter)[2], int bs = 1) {
            std::cout << ptCenter[0] << "," << ptCenter[1] << std::endl;
            // edge case detection
            if (ptCenter[0] > the_upper_bound_) {
                printf("Baobzi Family Warning: dist_perp %g very large, clamp to grid UB %g \n", ptCenter[0], the_upper_bound_);
                ptCenter[0] = the_upper_bound_ - grid_size_magnitude_; 
                bs = 0;  
            }

            if ((ptCenter[1] > integral_max_) && (ri == 3)) {
                printf("Baobzi Family Warning: integral %g very large, clamp to integral UB %g \n", ptCenter[1], integral_max_);
                ptCenter[1] = integral_max_ - integral_min_; 
                bs = 0;  
            }

            int gridNum = breakPtsCollChange.size();
            int pickPt = 0; 
            // binary search O(log n)
            if (bs == 1) {
                int kk = findIntervalIndex(breakPtsCollChange, ptCenter[0]); 
                if (breakPtsCollUnchange[kk] <= ptCenter[1]) {
                    pickPt = kk;
                } 
            }
            // brute force search O(n)
            else {
                for (int i = 0; i < gridNum - 1; i++) {
                    if ((breakPtsCollChange[i] <= ptCenter[0]) 
                        && (breakPtsCollChange[i+1] >= ptCenter[0]) 
                        // && (breakPtsCollUnchange[i] <= ptCenter[1])
                    ) {
                        pickPt = i; 
                        break;
                    }
                    pickPt = gridNum - 1; 
                }
            }            

            // speak("pickpt", pickPt); 
            Cheb pickBaobzi = allCheb[pickPt];
            int checkContain = pickBaobzi.checkInclude(ptCenter); 
            assert(checkContain == 1); 
            double approxVal = pickBaobzi.evalFunc(ptCenter);
            return approxVal; 
        }

        inline std::vector<std::vector<double>> findExtremeVal(int recorddata = 0, int calerror = 0) {
            std::vector<std::vector<double>> intSpecSaver; 
            std::vector<double> errorTrueSaver;
            // currently only allow normal lookup
            if (ri == 1) {
                std::ofstream myfile; 
                if (recorddata == 1) {
                    std::cout << "==START RECORD DATA==" << std::endl;
                    std::string rootpath = "3d-data/";
                    std::string strD = std::to_string(length_scale_);
                    std::string strAlpha = std::to_string(talpha);
                    std::string strFreeLength = std::to_string(tfreelength);
                    std::string searchfilename = rootpath + "D=" + strD + "-" + "alpha=" + strAlpha + "-" + "fl=" + strFreeLength + ".txt"; 
                    try{
                        std::filesystem::remove(searchfilename);
                    }
                    catch (...) {}
                    myfile.open(searchfilename);
                }
                int cnt = 0; 
                double prefac = 1; 
                for (double i = grid_size_magnitude_; i < the_upper_bound_ - grid_size_magnitude_; i += grid_size_magnitude_) {
                    std::vector<double> integralSaver; 
                    for (double j = grid_size_magnitude_; j < (the_upper_bound_ - grid_size_magnitude_) * length_scale_; j += grid_size_magnitude_ * length_scale_ / prefac){
                        double iter[2] = {i,j}; 
                        double intres = evalSinglePt(iter,0); 
                        integralSaver.push_back(intres); 
                        if (recorddata == 1){ 
                            // [vertical distance] << [scan length] << [lookup table result]
                            myfile << i << "," << j / length_scale_ << "," << intres; 
                        }
                        if (calerror == 1) {
                            double realres = length_scale_ * integral(i, 0, j / length_scale_, exp_fact_, rest_length_);
                            double relerr = ABS(intres - realres); 
                            errorTrueSaver.push_back(relerr); 
                            // << relative error with real integral value
                            myfile << "," << relerr; 
                        }
                        myfile << std::endl;
                    }
                    double maxval = *std::max_element(std::begin(integralSaver), std::end(integralSaver));
                    double minval = *std::min_element(std::begin(integralSaver), std::end(integralSaver));
                    // integral value must be positive
                    assert(maxval > 0); assert(minval > 0); 
                    if ((ABS(maxval) <= 1e-2) || (ABS(minval) <= 1e-2)) {
                        maxval = 0;
                        minval = 0; 
                    }
                    // std::cout << maxval << ";" << minval << std::endl; 
                    std::vector<double> perIntVal = {minval, maxval}; 
                    intSpecSaver.push_back(perIntVal); 
                    cnt++; 
                }
                if (recorddata == 1) {
                    std::cout << "==FINISH RECORD DATA==" << std::endl;
                    myfile.close(); 
                }
                if (calerror == 1){
                    std::cout << "==== Comparison Statistics [Inaccurate Time Measurement] ====" << std::endl;
                    speak("Mean Error of Cheb Family to Real Value over Whole Domain", mean_error(errorTrueSaver));
                    std::cout << "==== END ====" << std::endl; 
                }                

                speak("cnt",cnt); 
            }
            return intSpecSaver; 
        }

        inline std::vector<double> intMinMax(std::vector<std::vector<double>> intSpecSaver) {
            std::vector<double> res = findMinMaxVec(intSpecSaver);
            // load global max/min integral value [for reference]
            integral_min_ = res[0];
            integral_max_ = res[1]; 
            speak("integral_max_", integral_max_); 
            speak("integral_min_", integral_min_);
            return res;  

        }
};




#endif