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
        double alpha; double freelength; double D; 
        double bbtol; int runind; const char* fn; 
        std::vector<Cheb> allCheb; 
    
    private:
        static constexpr double small_ = 1e-4;
    
    public:
         Chebcoll() = default;
        ~Chebcoll() = default;

        Chebcoll(double alpha, double freelength, double D, double bbtol, const char* output_name, int runind) {
            // initialize dimensionless 
            length_scale_ = D; 
            exp_fact_ = alpha * length_scale_ * length_scale_;
            rest_length_ = freelength / length_scale_;
            // registeration
            alpha = alpha; freelength = freelength; D = D; 
            bbtol = bbtol; runind = runind; fn = output_name; 
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
            for (double iter = 10; iter < upBound; iter += 2 * otherGrid) {
                speak("create one",iter); 
                double hl[2] = {iter, oneFixCenter};
                double center[2] = {otherGrid, oneFixLength}; 
                Cheb theBaobzi(hl,center,bbtol,alpha,freelength,D,fn,runind);
                theBaobzi.approxFunc();
                allCheb.push_back(theBaobzi); 


            }
        }
};















#endif
