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
    public: 
        double talpha; double tfreelength; double tD; 
        double tbbtol; int ri; const char* on; 
        std::vector<Cheb> allCheb; 
        std::vector<double> breakPtsCollChange; 
        std::vector<double> breakPtsCollUnchange; 
    
    private:
        static constexpr double small_ = 1e-4;
    
    public:
         Chebcoll() = default;
        ~Chebcoll() = default;

        Chebcoll(double alpha, double freelength, double D, double bbtol, const char* output_name, const int runind) {
            // initialize dimensionless 
            length_scale_ = D; 
            exp_fact_ = alpha * length_scale_ * length_scale_;
            rest_length_ = freelength / length_scale_;
            // (dummy) registeration
            talpha = alpha; tfreelength = freelength; tD = D; 
            tbbtol = bbtol; on = output_name; ri = runind; 
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
            double oneFixLength = upBound / 2 * length_scale_; 
            assert(oneFixLength <= 1); 
            double oneFixCenter = oneFixLength; 
            double otherGrid = find_order(oneFixLength); // half-length
            // iteratively create Baobzi object
            for (double iter = otherGrid; iter < upBound - otherGrid; iter += 1 * otherGrid) {
                double hl[2] = {otherGrid, oneFixLength};
                double center[2] = {iter, oneFixCenter}; 
                /* stupid parameter loading error */
                // speak("alpha",talpha);
                // speak("freelength",tfreelength);
                // speak("D", tD); 
                Cheb theBaobzi(hl,center,tbbtol,talpha,tfreelength,tD,on,ri);
                theBaobzi.approxFunc();
                allCheb.push_back(theBaobzi); 
                breakPtsCollChange.push_back(iter - otherGrid); 
                breakPtsCollUnchange.push_back(oneFixCenter - oneFixLength); 
            }
        }

        inline double evalSinglePt(double (&ptCenter)[2]) {
            // std::cout << ptCenter[0] << "," << ptCenter[1] << std::endl;
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
            assert(checkContain); 
            double approxVal = pickBaobzi.evalFunc(ptCenter);
            return approxVal; 

        }
};















#endif
