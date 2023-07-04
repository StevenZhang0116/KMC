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
        double talpha; double tfreelength; double tD; 
        double tbbtol; int ri; const char* on; int tpon; 
        std::vector<Cheb> allCheb; 
        std::vector<double> breakPtsCollChange; 
        std::vector<double> breakPtsCollUnchange; 
        std::vector<double> chebSpaceTaken; 
    
    private:
        static constexpr double small_ = 1e-6;
    
    public:
        Chebcoll() = default;
        ~Chebcoll() = default;

        Chebcoll(double alpha, double freelength, double D, const int runind, double bbtol = 1e-4, const char* output_name = "func_approx.baobzi", 
             const int printOrNot = 0) {
            // initialize dimensionless 
            length_scale_ = D; 
            exp_fact_ = alpha * length_scale_ * length_scale_;
            rest_length_ = freelength / length_scale_;
            // (dummy) registeration
            talpha = alpha; tfreelength = freelength; tD = D; 
            tbbtol = bbtol; on = output_name; ri = runind; 
            tpon = printOrNot; 
            speak("Tolerance of Bobazi Family", bbtol); 
        }

        inline double calcBoltzmann(double dist_cent) const {
            double r = dist_cent / length_scale_;
            return exp(-exp_fact_ * SQR(r - rest_length_));
        }

        double getUpperBound() const {
            return sqrt(-log(small_) / exp_fact_) + rest_length_;
        }

        inline void createBaobziFamily(double ll = 0, double rr = 0) {
            double upBound; double oneFixLength; double oneFixCenter; double otherGrid; 
            if (ri == 1){
                std::cout << "==== Create Family of Positive Checking ====" << std::endl; 
                upBound = getUpperBound();
                the_upper_bound_ = upBound; 
                speak("UpBound (for both dimensions)", upBound);  
                oneFixLength = upBound / 2 * length_scale_; 
                assert(oneFixLength <= 1); 
                oneFixCenter = oneFixLength; 
                otherGrid = find_order(oneFixLength); // half-length
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
                speak("UpBound (for another dimension)", rr); 
                // rr = rr - 1e-2;
                double the_small_ = 1e-5; 
                oneFixLength = (rr - ll) / 2;
                assert(oneFixLength <= 1); 
                oneFixCenter = ll + oneFixLength; 
                otherGrid = find_order(oneFixLength);
                while (rr / otherGrid <= 1e2) {
                    otherGrid /= 10;
                    std::cout << "=== MAKE GRID SMALLER ===" << std::endl; 
                }
            }
            speak("OtherGrid", otherGrid); 
            speak("oneFixCenter", oneFixCenter); 
            speak("oneFixLength", oneFixLength); 

            // declare bound for the fixed parameter, usually distPerp (vertical distance)
            double lbound = otherGrid; 
            double ubound = upBound - otherGrid; 
            double gg = 1; // range [1,2], inversely proportional to running time
            double gridSize = gg * otherGrid; 

            // iteratively create Baobzi object
            speak("Baobzi Objects need to be created", (int)(ubound - lbound)/gridSize + 1); 
            speak("Grid Size Magnitude", otherGrid); // magnitude of grid size along both dimension
            grid_size_magnitude_ = otherGrid; 
            int cnt = 0; 
            for (double iter = lbound; iter < ubound; iter += gg * gridSize) {
                const auto st1 = get_wtime();
                double hl[2] = {otherGrid, oneFixLength};
                double center[2] = {iter, oneFixCenter}; 
                // speak("alpha",talpha);
                // speak("freelength",tfreelength);
                // speak("D", tD); 
                Cheb theBaobzi(hl, center, tbbtol, talpha, tfreelength, tD, on, ri, tpon, the_upper_bound_);
                double ssTaken = theBaobzi.approxFunc(); // taken space in Mb
                allCheb.push_back(theBaobzi); 
                breakPtsCollChange.push_back(iter - otherGrid); 
                breakPtsCollUnchange.push_back(oneFixCenter - oneFixLength); 
                chebSpaceTaken.push_back(ssTaken);
                cnt++; 
                if (cnt % 100 == 0) {
                    const auto ft1 = get_wtime();
                    const double dt1 = get_wtime_diff(&st1, &ft1);
                    speak("Baobzi Created", cnt); 
                    speak("Needed Time per object (s)", dt1); 
                }
            }
            speak("Total Baobzi Family Space (MiB)", total_sum(chebSpaceTaken)); 
        }

        inline double evalSinglePt(double (&ptCenter)[2], int bs = 1) {
            // std::cout << ptCenter[0] << "," << ptCenter[1] << std::endl;
            // edge case detectino
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
                    if ((breakPtsCollChange[i] <= ptCenter[0]) && 
                        (breakPtsCollChange[i+1] >= ptCenter[0]) &&
                        (breakPtsCollUnchange[i] <= ptCenter[1])
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
            /* To DO: For test purpose, add comparison with true value */
            return approxVal; 
        }

        inline std::vector<double> findExtremeVal(int recorddata = 0) {
            std::vector<double> integral_saver; 
            // currently only allow normal lookup
            if (ri == 1) {
                std::ofstream myfile; 
                if (recorddata == 1) {
                    std::cout << "==START RECORD DATA==" << std::endl;
                    std::string rootpath = "3d-data/";
                    std::string strD = std::to_string(tD);
                    std::string strAlpha = std::to_string(talpha);
                    std::string searchfilename = rootpath + "D=" + strD + "-" + "alpha=" + strAlpha + ".txt"; 
                    try{
                        std::filesystem::remove(searchfilename);
                    }
                    catch (...) {}
                    myfile.open(searchfilename);
                }
                for (double i = grid_size_magnitude_; i < the_upper_bound_ - grid_size_magnitude_; i += grid_size_magnitude_) {
                    for (double j = grid_size_magnitude_; j < (the_upper_bound_ - grid_size_magnitude_) * length_scale_; j += grid_size_magnitude_ * length_scale_){
                        double iter[2] = {i,j}; 
                        double intres = evalSinglePt(iter,0); 
                        integral_saver.push_back(intres); 
                        // [vertical distance] << [scan length] << [lookup table result]
                        myfile << i << "," << j / length_scale_ << "," << intres << std::endl;
                    }
                }
                if (recorddata == 1) {
                    std::cout << "==FINISH RECORD DATA==" << std::endl;
                    myfile.close(); 
                }
                /** TODO 
                  * what if change to min/max specific to every given interval (corresponding to normal lookup)
                  */ 
                double maxval = *std::max_element(std::begin(integral_saver), std::end(integral_saver));
                double minval = *std::min_element(std::begin(integral_saver), std::end(integral_saver));
                assert(maxval > 0); assert(minval > 0); // integral value must be positive
                integral_min_ = minval; integral_max_ = maxval; 
                std::vector<double> res = {minval, maxval}; 
                speakvec("Max/Min Integral", res); 
                return res; 
            }
            return integral_saver; 
        }

        inline void compareTrue() {
            std::vector<double> errorTrueSaver;
            for (double i = grid_size_magnitude_; i < the_upper_bound_ - grid_size_magnitude_; i += grid_size_magnitude_) {
                for (double j = grid_size_magnitude_; j < (the_upper_bound_ - grid_size_magnitude_) * length_scale_; j += grid_size_magnitude_ * length_scale_){
                    double iter[2] = {i,j};
                    double intres = evalSinglePt(iter,0);
                    double realres = length_scale_ * integral(i, 0, j / length_scale_, exp_fact_, rest_length_);
                    errorTrueSaver.push_back(ABS(intres - realres)); 
                }
            }
            speak("Mean Error of Cheb Family to Real Value", mean_error(errorTrueSaver));
            return; 
        }
};




#endif
