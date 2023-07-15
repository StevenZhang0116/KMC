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
#include <sys/stat.h>

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
        double talpha; 
        double tfreelength; 
        double tbbtol; 
        int tpon; 
        double etl; 
        int bbcount; 
        const char* tfoldername; 
        int tifsave; 
        // key
        int workIndex; 

        std::vector<Cheb> allCheb; 
        std::vector<double> breakPtsCollChange; 
        std::vector<double> breakPtsCollUnchange; 
        std::vector<double> chebSpaceTaken; 
    
    private:
        static constexpr double small_ = 1e-6;
        const double buffer_distance_ = 1e-2; 
    
    public:
        Chebcoll() = default;
        ~Chebcoll() = default;

        Chebcoll(double alpha, double freelength, double D, const int runind, double bbtol = 1e-4, 
        std::vector<double> integralMinMax = std::vector<double>(), const double errortolerence = 1e-3, 
        const int ifsave = 0, const char* defaultfoldername = "./mytests/", const int printornot = 0) {
            // initialize dimensionless 
            length_scale_ = D; 
            exp_fact_ = alpha * length_scale_ * length_scale_;
            rest_length_ = freelength / length_scale_;
            // (dummy) registeration of parameters
            talpha = alpha; 
            tfreelength = freelength;
            tbbtol = bbtol; 
            workIndex = runind; 
            tpon = printornot; 
            tifsave = ifsave; 
            // generate folder name -- user input
            tfoldername = defaultfoldername; 
            // 
            speak("Tolerance of Bobazi Family", bbtol); 
            // only need to set up integral min/max value in the REVERSE LOOKUP 
            if (integralMinMax.size() > 0) {
                assert(runind == 3); 
                integral_min_ = integralMinMax[0];
                integral_max_ = integralMinMax[1]; 
                std::cout << "Setup integral max/min value in Constructor for Reverse Lookup: " << integral_min_ << "," << integral_max_ << std::endl; 
            }
            etl = errortolerence; 
            speak("Tolerance of Approximating Integral to 0", etl); 
        }

        /**
         * @brief calculate Boltzmann factor
         * @return the value
         */

        inline double calcBoltzmann(double dist_cent) const {
            double r = dist_cent / length_scale_;
            return exp(-exp_fact_ * SQR(r - rest_length_));
        }

        /** 
         * @brief calculate upperbound for s and r⊥; if value beyonds the limit, clamp to that
         * @return UB
         */

        double getUpperBound() const {
            return sqrt(-log(small_) / exp_fact_) + rest_length_;
        }

        /**
         * @brief create Baobzi family over the whole domain with one dimension normalized and other dimension linearly discretized
         * [the hint of how to normalize both dimensions are included in the comments, but not desirable] 
         * save all objects respectively in vectors (member variables) and will be used in evalSinglePt() function
         * @param[in]: tempkk: optional, only used in reverse lookup to input prior knowledge of integral range (workIndex == 3)
         * @return void 
        */

        inline void createBaobziFamily(std::vector<std::vector<double>> tempkk = std::vector<std::vector<double>>()) {
            double upBound; 
            double oneFixLength; 
            double oneFixCenter; 
            double otherGrid; 
            double gridLoader = 1e-10; 
            double smallBound; 
            
            // create folder to save Baobzi object files if needed
            if (tifsave == 1) {
                int result = mkdir(tfoldername, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                if (result == 0) std::cout << "Folder created successfully." << std::endl;
                else std::cout << "Failed to create the folder." << std::endl;
            }

            std::vector<double> lengthVec; 
            std::vector<double> centerVec; 
            std::vector<double> gridVec;
            // normal lookup scanerio
            if ((workIndex == 1) || (workIndex == 5)){
                std::cout << "==== Create Family of Positive Checking ====" << std::endl; 
                upBound = getUpperBound();
                the_upper_bound_ = upBound; 
                speak("UpBound (for both dimensions)", upBound);  
                oneFixLength = upBound / 2 * length_scale_; 
                assert(oneFixLength <= 1); 
                oneFixCenter = oneFixLength; 
                // half-length
                otherGrid = find_order(oneFixLength); 
                // magnitude of grid size along both dimension
                speak("Grid Size Magnitude - 1", otherGrid); 
                while (the_upper_bound_ / otherGrid <= 1e2) {
                    otherGrid /= 10; 
                    std::cout << "=== MAKE GRID SMALLER === " << std::endl; 
                }
            }
            // reverse lookup scanerio
            else if (workIndex == 3) {
                double roww = tempkk.size();
                std::cout << "==== Create Family of Reverse Checking ====" << std::endl; 
                upBound = getUpperBound();
                the_upper_bound_ = upBound; 
                speak("UpBound (for one dimensions)", upBound);  
                otherGrid = find_order(upBound / roww);
                speak("Grid Size Magnitude - 3", otherGrid); 
                for (int i = 0; i < tempkk.size(); i++) {
                    double rr = tempkk[i][1];
                    double ll = tempkk[i][0];
                    oneFixLength = (rr - ll) / 2;
                    assert(oneFixLength <= 1); 
                    oneFixCenter = ll + oneFixLength; 
                    // match order
                    double tGrid = find_order(oneFixLength);
                    // randomly chosen ~0.01 seems to be a reasonable choice
                    double smallBound = 0.001; 
                    gridLoader = smallBound; 

                    if (tbbtol < gridLoader) {
                        tbbtol = gridLoader; 
                        std::cout << "Change tolerance to " << gridLoader << std::endl; 
                    }
                    
                    // whether the grids are consistently unifrom
                    int evengridIndex = 1; 
                    if (evengridIndex == 1) tGrid = smallBound; 
                    else tGrid = std::max(tGrid, smallBound);
                    
                    // should be integer, as division of 10's powers
                    double rrTimes = otherGrid / tGrid; 
                    int kk = round(rrTimes);
                    if (kk > 1) kk = round10(kk); 
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
            if ((workIndex == 1) || (workIndex == 5)) {
                // declare bound for the fixed parameter, usually distPerp (vertical distance)
                double lbound = otherGrid; 
                double ubound = upBound - otherGrid; 
                double gg = 1; // range [1,2], inversely proportional to running time
                double gridSize = gg * otherGrid; 

                for (double iter = lbound; iter < ubound; iter += gridSize) iterVec.push_back(iter); 
                bbcount = floor((ubound - lbound)/gridSize + 1); 
                speak("Baobzi Objects need to be created", bbcount); 
                grid_size_magnitude_ = otherGrid; 
            }
            else if (workIndex == 3){
                iterVec = cumulativeSum(gridVec); 
                grid_size_magnitude_ = smallBound; 
            }

            // temporary variables
            std::vector<Cheb> theallCheb(iterVec.size()); 
            std::vector<double> thebreakPtsCollChange(iterVec.size()); 
            std::vector<double> thebreakPtsCollUnchange(iterVec.size()); 
            std::vector<double> thechebSpaceTaken(iterVec.size()); 

            int tri; 
            double totalTime; 
            /* if want to normalize r⊥, change upper bound of iii to 1 (instead of iterVec.size()) */
            // #pragma omp parallel for
            for (size_t iii = 0; iii < iterVec.size(); iii++) {
                double thisCenter = iterVec[iii]; 
                const auto st1 = get_wtime();
                if (workIndex == 3) {
                    oneFixLength = lengthVec[iii];
                    oneFixCenter = centerVec[iii]; 
                    otherGrid = gridVec[iii]; 
                }
                /* if want to normalize r⊥, change the first coordinates of hl[] and center[] 
                 * to the second coordinate */ 
                double hl[2] = {otherGrid, oneFixLength};
                double center[2] = {thisCenter, oneFixCenter}; 
                // only happens in reverse lookup where the small integral value is clamped to 0
                if ((hl[1] == 0.0) && (center[1] == 0.0)) {
                    tri = 4; 
                    hl[1] = 1; 
                    center[1] = 1; 
                }
                else {
                    tri = workIndex; 
                }

                if ((workIndex == 3)) {
                    std::cout << "hl: " << otherGrid << ";" << oneFixLength << std::endl;
                    std::cout << "center: " << thisCenter << ";" << oneFixCenter << std::endl;
                }

                /* Create Baobzi Object */
                Cheb theBaobzi(hl, center, tbbtol, talpha, tfreelength, length_scale_, tri, tpon, etl, the_upper_bound_);
                if (tifsave == 1) theBaobzi.saveFunctionObject(tfoldername, iii); 

                size_t ssTaken = theBaobzi.approxFunc(); // taken space in Mb
                theallCheb[iii] = theBaobzi; 
                thebreakPtsCollChange[iii] = thisCenter - otherGrid; 
                thebreakPtsCollUnchange[iii] = oneFixCenter - oneFixLength; 
                thechebSpaceTaken[iii] = ssTaken;
                const auto ft1 = get_wtime();
                const double dt1 = get_wtime_diff(&st1, &ft1);
                totalTime += dt1; 
                // only print log in reverse lookup
                if ((iii % 100 == 0) && (workIndex == 3)) { 
                    speak("Baobzi Created", iii); 
                    speak("Needed Time per object (s)", dt1); 
                }
            }
            // transfer to MB
            speak("Total Baobzi Family Space (MiB)", total_sum(thechebSpaceTaken) / (1024 * 1024)); 
            speak("Total Time to Build Up Baobzi Family", totalTime); 

            // save temporary variables to member variable of class
            allCheb = theallCheb; 
            breakPtsCollChange = thebreakPtsCollChange; 
            breakPtsCollUnchange = thebreakPtsCollUnchange;
            chebSpaceTaken = thechebSpaceTaken; 
        }

        /**
         * @brief evaluate (either normal or reverse) at given input using either brute-force search or binary search
         * @param[in]: ptCenter: 2D coordinate of evaluated points
         * @param[in]: bs: index of using binary search or brute force 
         * @return calculated result
         */

        inline double evalSinglePt(double (&ptCenter)[2], int bs = 1) {
            // std::cout << ptCenter[0] << "," << ptCenter[1] << std::endl;
            // edge case detection
            if (ptCenter[0] > the_upper_bound_) {
                // printf("Baobzi Family Warning: dist_perp %g very large, clamp to grid UB %g \n", ptCenter[0], the_upper_bound_);
                ptCenter[0] = the_upper_bound_ - grid_size_magnitude_; 
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
                        && (breakPtsCollUnchange[i] <= ptCenter[1])
                    ) {
                        pickPt = i; 
                        break;
                    }
                    pickPt = gridNum - 1; 
                }
            }    
            
            // desired Baobzi object after searching for calculation
            Cheb pickBaobzi = allCheb[pickPt]; 
            // check if integral is beyond the scope; only calculate after knowing which exact object is used
            double intOB = pickBaobzi.center[1] + pickBaobzi.half_length[1]; 
            if ((ptCenter[1] > intOB) && (workIndex == 3)) {
                printf("Baobzi Family Warning: integral %g very large, clamp to integral UB %g \n", ptCenter[1], intOB);
                ptCenter[1] = intOB - std::max(integral_min_, 1e-5);
            }       
            // final check of whether evaluated point is in the domain
            int checkContain = pickBaobzi.checkInclude(ptCenter, 0); 
            if (checkContain == 1) {
                double approxVal = pickBaobzi.evalFunc(ptCenter);
                return approxVal; 
            }
            else {
                throw std::invalid_argument("BAOBZI DOMAIN ERROR");
            }
        }

        /**
         * @brief (up to choice) to record the data and relative error to true integral with respect to the s and r⊥ over the domain
         * @param[in]: recorddata: whether to record data in all kinds of lookup
         * @param[in]: recorderror: whether to record error in all kinds of lookup
         * @param[in]: optional, only used in reverse lookup to set up grid bound
         * @return min/max integral value over each linearly discretized grid in the domain
         */

        inline std::vector<std::vector<double>> findExtremeVal(int recorddata = 0, int recorderror = 0, int ubound = 0, double prefactor = 1) {
            // TODO: calculate upper bound of reverse lookup internally in the class? 
            if (workIndex == 3) assert(ubound > 0); 

            std::vector<std::vector<double>> intSpecSaver; 
            std::vector<double> errorTrueSaver;
            std::ofstream myfile; 
            int cnt = 0; 
            // ostream setup
            if (recorddata == 1) {
                std::cout << "==START RECORD DATA==" << std::endl;
                std::string rootpath = "3d-data/";
                std::string strD = std::to_string(length_scale_);
                std::string strAlpha = std::to_string(talpha);
                std::string strFreeLength = std::to_string(tfreelength);
                std::string categoryName; 
                // set up category name
                if ((workIndex == 1) || (workIndex == 5)) categoryName = "CDF-"; 
                if (workIndex == 3) categoryName = "ReverseCDF-"; 
                std::string searchfilename = rootpath + categoryName + "D=" + strD + "-" + "alpha=" + strAlpha + "-" + "fl=" + strFreeLength + ".txt"; 
                try {
                    std::filesystem::remove(searchfilename);
                }
                catch (...) {}
                myfile.open(searchfilename);
            }

            if ((workIndex == 1) || (workIndex == 5)) {
                // for calculation on global domain for energy dependent first-order CDF/PDF
                for (double i = grid_size_magnitude_; i < the_upper_bound_ - grid_size_magnitude_; i += grid_size_magnitude_) {
                    std::vector<double> gridResultSaver; 
                    for (double j = 0; j < (the_upper_bound_ - grid_size_magnitude_) * length_scale_; j += grid_size_magnitude_ * length_scale_ / prefactor){
                        double iter[2] = {i,j}; 
                        double intres = evalSinglePt(iter, 0); 
                        gridResultSaver.push_back(intres); 
                        if (recorddata == 1){ 
                            // [vertical distance] << [scan length] << [lookup table result]
                            myfile << i << "," << j / length_scale_ << "," << intres; 
                        }
                        if (recorderror == 1) {
                            double realres = length_scale_ * integral(i, 0, j / length_scale_, exp_fact_, rest_length_);
                            double relerr = ABS(intres - realres); 
                            errorTrueSaver.push_back(relerr); 
                            // << relative error with real integral value
                            myfile << "," << relerr; 
                        }
                        myfile << std::endl;
                    }
                    double maxval = *std::max_element(std::begin(gridResultSaver), std::end(gridResultSaver));
                    double minval = *std::min_element(std::begin(gridResultSaver), std::end(gridResultSaver));
                    // integral value must be positive
                    if ((ABS(maxval) <= etl) && (ABS(minval) <= etl)) {
                        maxval = 0;
                        minval = 0; 
                    }
                    std::vector<double> perIntVal = {minval, maxval}; 
                    intSpecSaver.push_back(perIntVal); 
                    cnt++; 
                }
                assert(cnt == bbcount); 
            }
            // for calculation on global domain for reverse lookup
            else if (workIndex == 3) {
                double bdd = 0.01; 
                for (double i = prefactor * bdd; i < ubound - 0.2; i += prefactor * bdd) {
                    double distPerp = i * length_scale_; 
                    std::vector<double> gridResultSaver; 
                    for (double j = prefactor * bdd; j < ubound - 0.2; j += prefactor * bdd) {
                        double val = integral(distPerp / length_scale_, 0, j, exp_fact_, rest_length_); 
                        double inval[] = {distPerp / length_scale_, val * length_scale_}; 
                        double revres = evalSinglePt(inval, 0); 
                        double relerr = ABS(revres - j * length_scale_); 
                        gridResultSaver.push_back(revres); 
                        if (recorddata == 1) {
                            myfile << i << "," << j << "," << revres; 
                        }  
                        if (recorderror == 1) {
                            errorTrueSaver.push_back(relerr);
                            myfile << "," << relerr; 
                        }
                        myfile << std::endl; 
                    }
                    cnt++; 
                }
            }

            if (recorddata == 1) {
                std::cout << "== FINISH RECORD DATA ==" << std::endl;
                myfile.close(); 
            }

            if (recorderror == 1){
                speak("Mean Error of Cheb Family to Real Value over Whole Domain", mean_error(errorTrueSaver));
                std::cout << "==== END ====" << std::endl; 
            }                
            speak("cnt",cnt); 

            return intSpecSaver; 
        }

        /**
         * @brief return global min/max integral value over the calculated domain
         * @return vector record global min/max integral value
         */

        inline std::vector<double> intMinMax(std::vector<std::vector<double>> intSpecSaver) {
            std::vector<double> res = findMinMaxVec(intSpecSaver);
            // load global max/min integral value [for reference]
            integral_min_ = res[0];
            integral_max_ = res[1]; 
            return res;  
        }
};




#endif