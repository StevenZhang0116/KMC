/**
 * @file chebcoll.hpp
 * @author Zihan Zhang
 * @brief Create multiple Baobzi objects (family) on the global search domain
 *
 * @copyright Copyright (c) 2023
 *
 */

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
#include <iomanip> 

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
        double tsmallBound; 
        // key -- work mode
        // [1]: Normal Lookup (CDF)
        // [3]: Reverse Lookup
        // [5]: Normal Lookup (PDF)
        // rests are still under developed
        int tworkIndex; 

        std::vector<Cheb> allCheb; 
        std::vector<double> breakPtsCollChange; 
        std::vector<double> breakPtsCollUnchange; 
        std::vector<double> chebSpaceTaken; 
    
    private:
        static constexpr double SMALL_ = 1e-6;
        // declaration of function approximation
        const int NORMAL_CDF = 1;
        const int NORMAL_REVERSECDF = 3;
        const int CONST_ZERO = 4; 
        const int NORMAL_PDF = 5; 
    
    public:
        Chebcoll() = default;
        ~Chebcoll() = default;

        /**
         * Baobzi Family (colleection of objects) constructor from scratch
         */
        Chebcoll(const int ifsave, double alpha, double freelength, double D, const int runind, double bbtol = 1e-4, const double errortolerence = 1e-3, 
        std::vector<double> integralMinMax = std::vector<double>(), const double smallBound = 1e-10,
        const char* defaultfoldername = "./mytests/", const int printornot = 0) {
            std::cout << "************************************" << std::endl;
            std::cout << "=== Calculation (from scratch) Start ===" << std::endl; 
            std::cout << "************************************" << std::endl;
            // initialize dimensionless 
            length_scale_ = D; 
            exp_fact_ = alpha * length_scale_ * length_scale_;
            rest_length_ = freelength / length_scale_;
            // (dummy) registeration of parameters
            talpha = alpha; 
            tfreelength = freelength;
            tbbtol = bbtol; 
            tworkIndex = runind; 
            tpon = printornot; 
            tifsave = ifsave; 
            tsmallBound = smallBound; 
            // generate folder name -- user input
            tfoldername = defaultfoldername; 
            speak("Tolerance of Bobazi Family", bbtol); 
            // only need to set up integral min/max value in the REVERSE LOOKUP 
            if (integralMinMax.size() > 0) {
                assert(runind == NORMAL_REVERSECDF); 
                integral_min_ = integralMinMax[0];
                integral_max_ = integralMinMax[1]; 
                std::cout << "Setup integral max/min value in Constructor for Reverse Lookup: "; 
                std::cout << integral_min_ << "," << integral_max_ << std::endl; 
            }
            etl = errortolerence; 
            speak("Tolerance of Approximating Integral to 0", etl); 
        }

        /**
         * Read in already created Baobzi object into tbe constructor based on the given root
         * The current folder path in demo test is KMC/build/mytests
         * Only needed in reverse lookup, where the build up is time costly
         */
        Chebcoll(const int runind, std::string objpath = "./mytests/") {
            std::cout << "*******************************************" << std::endl;
            std::cout << "=== Calculation (from reconstruction) Start ===" << std::endl; 
            std::cout << "*******************************************" << std::endl;
            // load work mode
            tworkIndex = runind; 
            // all filenames are integer-formatted, so the following operations are eligible
            std::vector<int> saveNameIntForm; 
            
            for (const auto& entry : std::filesystem::directory_iterator(objpath)){
                if (std::filesystem::is_regular_file(entry)) {
                    std::string fileName = entry.path().filename().string();
                    try {
                        // will have error for [param] file, thus automatically ignored
                        saveNameIntForm.push_back(std::stoi(fileName)); 
                    }
                    catch (const std::exception& ex) {
                        std::cerr << "Error: " << ex.what() << "," << fileName << std::endl;
                    }
                }
            }

            speak("Num of Files in Folder", saveNameIntForm.size()); 
            // sort vector -- dependent on how the Baobzi objects are sequentially created
            std::sort(saveNameIntForm.begin(), saveNameIntForm.end()); 
            // delete duplicated filenames due to scanning over Baobzi + its auxiliary file
            std::vector<int> uniqNameIntForm = getUniqueElements(saveNameIntForm); 
            // num of objects that need to be created
            int totalNum = uniqNameIntForm.size(); 
            speak("totalNum", totalNum); 
            assert(totalNum > 0); 
            // create temporary variables 
            std::vector<Cheb> theallCheb(totalNum); 
            std::vector<double> thebreakPtsCollChange(totalNum); 
            std::vector<double> thebreakPtsCollUnchange(totalNum);
            // load parameters
            std::vector<double> readInParam = readDataFromFile(objpath + "param"); 
            
            // follow the format in [recordParameter()]
            exp_fact_ = readInParam[0]; 
            rest_length_ = readInParam[1]; 
            length_scale_ = readInParam[2]; 
            the_upper_bound_ = readInParam[3]; 
            grid_size_magnitude_ = readInParam[4]; 
            integral_min_ = readInParam[5]; 
            integral_max_ = readInParam[6]; 

            const auto st1 = get_wtime();
            for (size_t i = 0; i < totalNum; i++) {
                int temp = uniqNameIntForm[i]; 
                std::string resFileName = std::to_string(temp);
                std::string fullFilePath = objpath + resFileName; 
                std::string fileSuffix = "-Res"; 
                std::string fullAuxFilePath = fullFilePath + fileSuffix; 
                // readin auxiliary parameters
                std::vector<double> readInData = readDataFromFile(fullAuxFilePath);
                // following the save format defined in [cheb.hpp], function [saveFunctionObject]
                double hl[2] = {readInData[1], readInData[3]};
                double center[2] = {readInData[0], readInData[2]}; 
                const int porn = 0; 
                Cheb theCheb(hl, center, fullFilePath, porn); 
                // load variables, similar to createBaobziFamily()
                theallCheb[i] = theCheb; 
                thebreakPtsCollChange[i] = readInData[0] - readInData[1];
                thebreakPtsCollUnchange[i] = readInData[2] - readInData[3];
            }
            // save to member variable of class
            allCheb = theallCheb; 
            breakPtsCollChange = thebreakPtsCollChange; 
            breakPtsCollUnchange = thebreakPtsCollUnchange;
            // calculate reconstruction time
            const auto ft1 = get_wtime();
            const double dt1 = get_wtime_diff(&st1, &ft1);
            speak("Total Reconstruction Time (s)", dt1); 

            std::uintmax_t folderSize = calculateFolderSize(objpath);
            speak("Total Memory Size in Reconstructed Folder (MB)", folderSize / (1024 * 1024)); 
        }

        /**
         * @brief calculate Boltzmann factor
         * @return res
         */

        inline double calcBoltzmann(double dist_cent) const {
            double r = dist_cent / length_scale_;
            return exp(-exp_fact_ * SQR(r - rest_length_));
        }

        /** 
         * @brief calculate upperbound for s and r⊥; if value beyonds the limit, clamp to that
         * @return UB
         */

        inline double getUpperBound() const {
            return sqrt(-log(SMALL_) / exp_fact_) + rest_length_;
        }

        /**
         * @brief create Baobzi family over the whole domain with one dimension normalized and other dimension linearly discretized
         * [the hint of how to normalize both dimensions are included in the comments, but not desirable] 
         * save all objects respectively in vectors (member variables) and will be used in [evalSinglePt()] function
         * @param[in]: tempkk: optional, only used in reverse lookup to input prior knowledge of integral range (tworkIndex == NORMAL_REVERSECDF)
         * @param[in]: prefactor: optional, constant factor to manipulate linear grid [factor of 10]
         * -- under current tests, [0.05,1] is a reasonable range 
         * @return[out]: 1: required space (in MB) of this BF object
         * @return[out]: 2: required build time (in s) of this BF object
        */

        inline std::pair<double, double> createBaobziFamily(const double prefactor = 1, std::vector<std::vector<double>> tempkk = std::vector<std::vector<double>>()) {

            double upBound; 
            double oneFixLength; 
            double oneFixCenter; 
            double otherGrid; 
            double smallBound; 
            
            // create folder to save Baobzi object files if needed
            if (tifsave == 1) {
                int result = mkdir(tfoldername, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                if (result == 0) {
                    std::cout << "Folder created successfully." << std::endl;
                }
                else {
                    removeAllFilesInFolder(tfoldername);
                    speak("Remove All Files in", tfoldername); 
                }
            }

            std::vector<double> lengthVec; 
            std::vector<double> centerVec; 
            std::vector<double> gridVec;
            // normal lookup scanerio (CDF or PDF)
            // TODO: ENUM CLASS
            if ((tworkIndex == NORMAL_CDF) || (tworkIndex == NORMAL_PDF)){
                std::cout << "==== Create Family of Positive Checking ====" << std::endl; 
                upBound = getUpperBound();
                the_upper_bound_ = upBound; 
                speak("UpBound (for both dimensions)", upBound);  
                oneFixLength = upBound / 2 * length_scale_; 
                assert(oneFixLength <= 1); 
                oneFixCenter = oneFixLength; 
                // adaptive to user input prefactor and dependent to [oneFixLength] through calculation
                otherGrid = find_order(oneFixLength) * prefactor; 

                // magnitude of grid size along both dimension
                speak("Grid Size Magnitude - 1", otherGrid); 
                while (the_upper_bound_ / otherGrid <= 1e2) {
                    otherGrid /= 10; 
                    std::cout << "=== MAKE GRID SMALLER === " << std::endl; 
                }
            }
            // reverse lookup scanerio
            else if (tworkIndex == NORMAL_REVERSECDF) {
                double roww = tempkk.size();
                std::cout << "==== Create Family of Reverse Checking ====" << std::endl; 
                upBound = getUpperBound();
                the_upper_bound_ = upBound; 
                speak("UpBound (for one dimensions)", upBound);  
                otherGrid = find_order(upBound / roww);
                speak("Grid Size Magnitude - 3", otherGrid); 
                for (int i = 0; i < roww; i++) {
                    double rr = tempkk[i][1];
                    double ll = tempkk[i][0];
                    oneFixLength = (rr - ll) / 2;
                    assert(oneFixLength <= 1); 
                    oneFixCenter = ll + oneFixLength; 
                    // match order
                    double tGrid = find_order(oneFixLength);
                    // randomly chosen ~0.01 seems to be a reasonable choice
                    smallBound = 0.01 * prefactor; 
                    
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
            if ((tworkIndex == NORMAL_CDF) || (tworkIndex == NORMAL_PDF)) {
                // declare bound for the fixed parameter, usually distPerp (vertical distance)
                double lbound = otherGrid; 
                double ubound = upBound - otherGrid; 
                /* if want to normalize r⊥, set: 
                "" 
                ubound = lbound + 1e-10; 
                ""
                */
               ubound = lbound + 1e-10; 

                double gg = 1; // range [1,2], inversely proportional to running time
                double gridSize = gg * otherGrid; 
                for (double iter = lbound; iter < ubound; iter += gridSize) iterVec.push_back(iter); 
                grid_size_magnitude_ = gridSize; 
            }
            else if (tworkIndex == NORMAL_REVERSECDF){
                iterVec = cumulativeSum(gridVec); 
                grid_size_magnitude_ = smallBound; 
            }
            // How many Baobzi Object in current Chebcoll 
            bbcount = iterVec.size(); 
            speak("Baobzi Objects need to be created", bbcount); 


            // temporary variables
            std::vector<Cheb> theallCheb(bbcount); 
            std::vector<double> thebreakPtsCollChange(bbcount); 
            std::vector<double> thebreakPtsCollUnchange(bbcount); 
            std::vector<double> thechebSpaceTaken(bbcount); 

            int tri; 
            double totalTime; 
            // #pragma omp parallel for
            for (size_t iii = 0; iii < iterVec.size();  iii++) {
                double thisCenter = iterVec[iii]; 
                const auto st1 = get_wtime();
                if (tworkIndex == NORMAL_REVERSECDF) {
                    oneFixLength = lengthVec[iii];
                    oneFixCenter = centerVec[iii]; 
                    otherGrid = gridVec[iii]; 
                }
                /* if want to normalize r⊥, set: 
                ""
                thisCenter = oneFixCenter;
                otherGrid = oneFixCenter; 
                ""
                */

                thisCenter = oneFixCenter;
                otherGrid = oneFixCenter; 

                double hl[2] = {otherGrid, oneFixCenter};
                double center[2] = {thisCenter, oneFixCenter}; 
                // only happens in reverse lookup where the small integral value is clamped to 0
                if ((hl[1] == 0.0) && (center[1] == 0.0)) {
                    tri = CONST_ZERO;
                    // change domain specification 
                    hl[1] = 1; 
                    center[1] = 1; 
                }
                else {
                    tri = tworkIndex; 
                }

                // Create Baobzi Object
                Cheb theBaobzi(hl, center, tbbtol, talpha, tfreelength, length_scale_, tri, tpon, etl, the_upper_bound_, tsmallBound);
                // taken space in Byte -- Function Approximated (at this point)!
                size_t ssTaken = theBaobzi.approxFunc(); 
                // save Baobzi object to an external file
                if (tifsave == 1) {
                    theBaobzi.saveFunctionObject(tfoldername, iii); 
                }
                theallCheb[iii] = theBaobzi; 
                // important member variable that needed to be reconstructed when readin files
                // load variables
                thebreakPtsCollChange[iii] = thisCenter - otherGrid; 
                thebreakPtsCollUnchange[iii] = oneFixCenter - oneFixLength; 
                thechebSpaceTaken[iii] = ssTaken;
                const auto ft1 = get_wtime();
                const double dt1 = get_wtime_diff(&st1, &ft1);
                totalTime += dt1; 
                // only print log in reverse lookup
                if ((iii % 100 == 0) && (tworkIndex == NORMAL_REVERSECDF)) { 
                    speak("Baobzi Created", iii); 
                    speak("Needed Time per object (s)", dt1); 
                }
            }
            // record parameter list for reconstruction purpose
            if (tifsave == 1) recordParameter(); 
            // transfer to MB
            double requiredSpace = total_sum(thechebSpaceTaken) / (1024 * 1024); 
            speak("Total Baobzi Family Space (MiB)", requiredSpace); 
            speak("Total Time to Build Up Baobzi Family", totalTime); 

            // save temporary variables to member variable of class
            allCheb = theallCheb; 
            breakPtsCollChange = thebreakPtsCollChange; 
            breakPtsCollUnchange = thebreakPtsCollUnchange;
            chebSpaceTaken = thechebSpaceTaken; 

            return std::make_pair(requiredSpace, totalTime); 

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
                printf("Baobzi Family Warning: dist_perp %g very large, clamp to grid UB %g \n", ptCenter[0], the_upper_bound_);
                ptCenter[0] = the_upper_bound_ - grid_size_magnitude_; 
                bs = 0;  
            }

            int gridNum = breakPtsCollChange.size();
            int pickPt = 0; 

            /* if want to normalize r⊥, change ifsearch = 0 */
            // we know which Baobzi object to use (since only 1) so no need to search for
            int ifsearch = 0; 

            if (ifsearch == 1){
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
                        // interval search
                        // only use one dimension to filter [might be] enough
                        if ((breakPtsCollChange[i] <= ptCenter[0]) 
                            && (breakPtsCollChange[i + 1] >= ptCenter[0]) 
                            // && (breakPtsCollUnchange[i] <= ptCenter[1])
                        ) {
                            pickPt = i; 
                            break;
                        }
                        pickPt = gridNum - 1; 
                    }
                }    
            }
            
            // desired Baobzi object after searching for calculation
            Cheb pickBaobzi = allCheb[pickPt]; 
            // check if integral is beyond the scope; only calculate after knowing which exact object is used
            double intOB = pickBaobzi.center[1] + pickBaobzi.half_length[1]; 
            if ((ptCenter[1] > intOB) && (tworkIndex == NORMAL_REVERSECDF)) {
                printf("Baobzi Family Warning: integral %g very large, clamp to integral UB %g \n", ptCenter[1], intOB);
                ptCenter[1] = intOB - std::max(integral_min_, 1e-5);
            }       
            // final check of whether evaluated point is in the domain
            // might be redundant, 
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
         * @param[in]: recorddata [binary]: whether to record data in all kinds of lookup
         * @param[in]: recorderror [binary]: whether to record error in all kinds of lookup
         * @param[in]: ubound: optional, only used in reverse lookup to set up grid bound
         * @param[in]: prefactor [positive int]: grid of recording data, =1 -> =grid of linearly discretized Baobzi object
         *             <1 -> finer grid in evaluation
         * @param[in]: relOrAbs [binary]: whether to record relative error (0 <= [err] <= 1) or absolute error (with unit)
         * @return: 1: min/max integral value over each linearly discretized grid in the domain [only needed in regular lookup, empty otherwise (currently)]
         * @return: 2: average error over the calculated domain using unifrom grids (0 if [recorderror] = 0)
         * @return: 3: average evaluation time of uiformly distributed samples in the domain 
         */

        inline std::tuple<std::vector<std::vector<double>>, double, double> scanGlobalDomain(int recorddata = 0, int recorderror = 0, 
        int ubound = 0, double prefactor = 1, const int relOrAbs = 0) {
            // TODO: calculate upper bound of reverse lookup internally in the class? 
            if (tworkIndex == NORMAL_REVERSECDF) assert(ubound > 0); 

            // suffix of rel or abs error
            std::string errSuffix; 
            if (relOrAbs == 0) errSuffix = "abs-";
            else errSuffix = "rel-"; 

            std::vector<std::vector<double>> intSpecSaver; 
            std::vector<double> errorTrueSaver;
            std::ofstream myfile; 

            // row counter
            int cnt = 0; 
            // total case counter
            int cntTotal = 0;
            // total evaluation time
            double evalTime = 0; 
            // total error loader
            double totalMeanError = 0; 

            // ostream setup
            if (recorddata == 1) {
                std::cout << "==START RECORD DATA==" << std::endl;
                std::string rootPath = "3d-data/";
                std::string strD = std::to_string(round_up(length_scale_, 3));
                std::string strAlpha = std::to_string(round_up(talpha, 3));
                std::string strFreeLength = std::to_string(round_up(tfreelength, 3));
                std::string categoryName; 
                // set up category name
                if (tworkIndex == NORMAL_CDF) categoryName = "CDF-"; 
                else if (tworkIndex == NORMAL_REVERSECDF) categoryName = "ReverseCDF-"; 
                else if (tworkIndex == NORMAL_PDF) categoryName = "PDF-"; 
                std::string searchfilename = rootPath + categoryName + errSuffix + "D=" + strD + "-" + "alpha=" + strAlpha + "-" + "fl=" + strFreeLength + ".txt"; 
                try {
                    std::filesystem::remove(searchfilename);
                }
                catch (...) {}
                myfile.open(searchfilename);
            }

            if ((tworkIndex == NORMAL_CDF) || (tworkIndex == NORMAL_PDF)) {
                // in case too heavy workload
                // TODO: should it be [iterGrid = grid_size_magnitude_] instead？ 
                double iterGrid = std::max(grid_size_magnitude_, 0.1); 
                // for calculation on global domain for energy dependent first-order CDF/PDF
                /* if want to normalize r⊥, * length_scale_ for low bound, up bound, and grid size in for loop declaration */
                for (double i = iterGrid * length_scale_; i < (the_upper_bound_ - iterGrid) * length_scale_; i += iterGrid * length_scale_) {
                    // speak("curr", i); 
                    std::vector<double> gridResultSaver; 
                    for (double j = 0; j < (the_upper_bound_ - iterGrid) * length_scale_; j += iterGrid * length_scale_ * prefactor){
                        double iter[2] = {i, j}; 

                        // count time
                        const auto st = get_wtime();
                        // point evaluation
                        double intres = evalSinglePt(iter, 0); 
                        const auto ft = get_wtime();
                        // time difference of evaluation
                        const double dt = get_wtime_diff(&st, &ft);
                        // increment
                        evalTime += dt; 

                        gridResultSaver.push_back(intres); 
                        /* if want to normalize r⊥, / length_scale_ after i */
                        double realres = length_scale_ * integral(i / length_scale_, 0, j / length_scale_, exp_fact_, rest_length_);

                        // calculate relative or abselute error
                        double relerr; 
                        if (relOrAbs == 0) relerr = ABS(intres - realres); 
                        else relerr = ABS(intres - realres) / ABS(realres); 

                        if (recorddata == 1){ 
                            // [vertical distance] << [scan length] << [lookup table result]
                            /* if want to normalize r⊥, / length_scale_ after i */
                            myfile << i / length_scale_ << "," << j / length_scale_ << "," << intres; 
                        }
                        errorTrueSaver.push_back(relerr); 
                        if (recorderror == 1) {
                            // [error] 
                            myfile << "," << relerr; 
                        }
                        if ((recorddata == 1) || (recorderror == 1)) {
                            myfile << std::endl; 
                        }
                        cntTotal++; 
                    }
                    // formulate output (only useful in reverse lookup)
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
                // assert(cnt == bbcount); 
            }

            // for calculation on global domain for reverse lookup
            else if (tworkIndex == NORMAL_REVERSECDF) {
                // in case too heavy workload
                double iterGrid = std::max(grid_size_magnitude_, 0.01); 
                for (double i = prefactor * iterGrid; i < ubound - iterGrid; i += prefactor * iterGrid) {
                    double distPerp = i * length_scale_; 
                    std::vector<double> gridResultSaver; 
                    for (double j = prefactor * iterGrid; j < ubound - iterGrid; j += prefactor * iterGrid) {
                        double val = length_scale_ * integral(distPerp / length_scale_, 0, j, exp_fact_, rest_length_); 
                        double inval[] = {distPerp / length_scale_, val}; 

                        const auto st = get_wtime();
                        double revres = evalSinglePt(inval, 0); 
                        const auto ft = get_wtime();
                        const double dt = get_wtime_diff(&st, &ft);
                        evalTime += dt; 

                        gridResultSaver.push_back(revres); 
                        double trueval = j * length_scale_; 

                        // calculate relative or absolute error
                        double relerr; 
                        if (relOrAbs == 0) relerr = ABS(revres - trueval); 
                        else relerr = ABS(revres - trueval) / ABS(trueval); 
                        gridResultSaver.push_back(revres); 

                        if (recorddata == 1) {
                            // [vertical distance] << [scan length (dimensional) -- true value] << [reconstructed value (dimensionless) -- calculated value]
                            myfile << i << "," << j << "," << revres; 
                        }  
                        errorTrueSaver.push_back(relerr);
                        if (recorderror == 1) {
                            // [error]
                            myfile << "," << relerr; 
                        }
                        if ((recorddata == 1) || (recorderror == 1)){
                            // [integral ("correct" value)]
                            myfile << "," << val; 
                            myfile << std::endl; 
                        }
                        cntTotal++; 
                    }
                    cnt++; 
                }
            }

            if (recorddata == 1) {
                std::cout << "== FINISH RECORD DATA ==" << std::endl;
                myfile.close(); 
            }

            if (recorderror == 1){
                totalMeanError = mean_error(errorTrueSaver); 
                std::cout << errSuffix << std::endl; 
                speak("Mean Error of Cheb Family to Real Value over Whole Domain", totalMeanError);
                std::cout << "==== END ====" << std::endl; 
            }    

            double evalTimePerSample = evalTime / cntTotal; 
            return std::make_tuple(intSpecSaver, totalMeanError, evalTimePerSample); 
        }

        /**
         * @brief return global min/max integral value over the calculated domain
         * @param[in]: intSpecSaver: 2D matrix of all data in domain recorded by scanGlobalDomain()
         * @return vector record global min/max integral value
         */

        inline std::vector<double> intMinMax(std::vector<std::vector<double>> intSpecSaver) {
            std::vector<double> res = findMinMaxVec(intSpecSaver);
            // load global max/min integral value [for reference]
            integral_min_ = res[0];
            integral_max_ = res[1]; 
            return res;  
        }

        /**
         * @brief record important parameters of pre-calculated Chebcoll object to an external file
         */

        inline void recordParameter() {
            std::ofstream myfile;
            std::string suffix = "param"; 
            std::string paramFileName = std::string(tfoldername) + suffix; 
            const char* saveLocation = paramFileName.c_str();
            speak("saveLocation", saveLocation); 
            int result = std::remove(saveLocation);
            if (result == 0) std::cout << "Parameter File Detected and Deleted: " << saveLocation << std::endl; 
            myfile.open(saveLocation);
            // default value registeration
            if ((integral_min_ < 1e-30) || (integral_min_ > 1e30)) integral_min_ = 0; 
            if ((integral_max_ < 1e-30) || (integral_max_ > 1e30)) integral_max_ = 0; 

            myfile << std::setprecision(20);
            // ---- saving format ----
            myfile << exp_fact_ << "," << rest_length_ << "," << length_scale_ << "," << the_upper_bound_ << ","; 
            myfile << grid_size_magnitude_ << "," << integral_min_ << "," << integral_max_ << std::endl; 
            // ---- saving format ----
            myfile.close(); 
            std::cout << "Parameter List Finished Recording" << std::endl; 

        }
};



#endif