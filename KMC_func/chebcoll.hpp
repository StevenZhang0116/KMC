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

        Chebcoll(double alpha, double freelength, double D, double bbtol = 1e-5, const char* output_name = "func_approx.baobzi", 
            const int runind = 1, const int printOrNot = 0) {
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

        inline void createBaobziFamily() {
            double upBound = getUpperBound();
            the_upper_bound_ = upBound; 
            speak("UpBound (for both dimensions)", upBound);  
            double oneFixLength = upBound / 2 * length_scale_; 
            assert(oneFixLength <= 1); 
            double oneFixCenter = oneFixLength; 
            double otherGrid = find_order(oneFixLength); // half-length
            while (the_upper_bound_ / otherGrid <= 1e2) {
                otherGrid /= 10; 
                std::cout << "=== MAKE GRID SMALLER === " << std::endl; 
            }
            // declare bound 
            double lbound = otherGrid; 
            double ubound = upBound - otherGrid; 
            double gg = 1; // range [1,2], inversely proportional to running time
            double gridSize = gg * otherGrid; 
            // iteratively create Baobzi object
            speak("Baobzi Objects need to be created", (int)(ubound - lbound)/gridSize + 1); 
            speak("Grid Size Magnitude", otherGrid); // magnitude of grid size along both dimension
            grid_size_magnitude_ = otherGrid; 
            for (double iter = lbound; iter < ubound; iter += gridSize) {
                double hl[2] = {otherGrid, oneFixLength};
                double center[2] = {iter, oneFixCenter}; 
                // speak("alpha",talpha);
                // speak("freelength",tfreelength);
                // speak("D", tD); 
                Cheb theBaobzi(hl, center, tbbtol, talpha, tfreelength, tD, on, ri, tpon);
                double ssTaken = theBaobzi.approxFunc(); // taken space in Mb
                allCheb.push_back(theBaobzi); 
                breakPtsCollChange.push_back(iter - otherGrid); 
                breakPtsCollUnchange.push_back(oneFixCenter - oneFixLength); 
                chebSpaceTaken.push_back(ssTaken);
            }
            speak("Total Baobzi Family Space (MiB)", total_sum(chebSpaceTaken)); 
        }

        inline double evalSinglePt(double (&ptCenter)[2]) {
            // std::cout << ptCenter[0] << "," << ptCenter[1] << std::endl;
            if (ptCenter[0] > the_upper_bound_) {
                printf("Warning: dist_perp %g very large, clamp to grid UB %g \n", ptCenter[0], the_upper_bound_);
                ptCenter[0] = the_upper_bound_ - grid_size_magnitude_; 
            }
            int gridNum = breakPtsCollChange.size();
            int pickPt = 0; 
            for (int i = 0; i < gridNum - 1; i++) {
                if ((breakPtsCollChange[i] <= ptCenter[0]) && (breakPtsCollChange[i+1] >= ptCenter[0]) &&
                    (breakPtsCollUnchange[i] <= ptCenter[1])
                ) {
                    pickPt = i; 
                    break;
                }
                pickPt = gridNum - 1; 
            }
            // speak("pickpt", pickPt); 
            Cheb pickBaobzi = allCheb[pickPt];
            int checkContain = pickBaobzi.checkInclude(ptCenter); 
            assert(checkContain == 1); 
            double approxVal = pickBaobzi.evalFunc(ptCenter);
            return approxVal; 

        }
};




#endif
