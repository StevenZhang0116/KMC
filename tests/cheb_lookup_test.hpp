#include "KMC_func/integrals.hpp"
#include "KMC_func/lookup_table.hpp"
#include "KMC_func/lut_filler_edep.hpp"
#include "KMC_func/macros.hpp"
#include "catch.hpp"
#include "test_helpers.hpp"
#include "KMC_func/cheb.hpp"
#include "KMC_func/auxil.hpp"

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

#include <boost/math/quadrature/gauss_kronrod.hpp>

/** 
 * Reverse look up search 
 * Problem with large stiffness -> CDF similar to step function
 * How to automatically set up search domain if no prior knowledge of it
 */

TEST_CASE("REVERSE Lookup table test (all kind) spring ", "[REVERSE lookup]") {
    // std::vector<double> hh {0.05,0.01,0.15,0.2};
    // for (size_t i = 0; i < hh.size(); i++) {
    std::cout << "==== REVERSE CHECK SPRING TEST ====" << std::endl; 
    const double tol = 1e-2;
    const double D = 0.024;
    const double alpha = 0.1 / (2 * 0.00411);
    const double freelength = 0.05;
    const double M = alpha * D * D;
    const double ell0 = freelength / D;

    static constexpr double small_ = 1e-30; 

    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init(alpha, freelength, D);
    LookupTable LUT(&lut_filler);

    double distPerp = 0;
    distPerp = 0.2;

    double testbound = LUT.getNonDsbound()/2; 
    double startbound = 0.1; 
    speak("testbound", testbound); 
    double boundgrid = 0.1; 
    size_t gridcnt = (testbound - startbound) / boundgrid; 
    std::vector<double> intval; 
    std::vector<double> rlerr; 
    std::vector<double> baobzierr; 

    // ("distPerp = 0.2 > D+ell0, single peaked")
    for (double sbound = startbound; sbound < testbound; sbound += boundgrid) {
        double val = integral(distPerp / D, 0, sbound, M, ell0);
        double a1 = LUT.ReverseLookup(distPerp, val * D); 
        // speak("val",val); 
        double err = ABS(a1 - sbound * D); 
        rlerr.push_back(err); 
        intval.push_back(val); 
        // CHECK(LUT.ReverseLookup(distPerp, val * D) == Approx(sbound * D).epsilon(tol));
    }
    speak("small", intval[0]); speak("large",intval[gridcnt]); 
    // speakvec(intval, gridcnt); 
    double intvaldiff = (intval[gridcnt] - intval[0]) / 2; 
    double midint = intval[0] + intvaldiff; 

    int dim = 2;
    int order = 10; 
    int odim = 1;
    double mlf = 0.0;
    int sme = 1;
    int mind = 0; 
    int maxd = 40;
    double hlfc = intvaldiff * D; 
    double hl[] = {powf(10.0f, floorf(log10f(hlfc))), (hlfc+ small_) * 1}; // half length
    double center[] = {distPerp / D, (midint * D + small_) * 1};  // center
    const char* fn = "func_approx.baobzi"; 

    double bbtol = 1e-6; 
    speak("tolerence",bbtol); 

    // std::string rootpath = "rvl-res/time/";
    // std::ofstream myfile; 
    // std::string strparam = std::to_string(distPerp);
    // std::string strTol = std::to_string(-1 * log10(bbtol));
    // std::string searchfilename = rootpath + strparam + "-" + strTol + ".txt"; 
    // try {
    //     std::filesystem::remove(searchfilename);
    // }
    // catch (...) {}
    // myfile.open(searchfilename);

    const auto st1 = get_wtime();
    Cheb theBaobzi(hl[0],hl[1],center[0],center[1],dim,odim,order,bbtol,mlf,sme,mind,maxd,M,ell0,D,fn);
    theBaobzi.approxFunc(3);

    const auto ft1 = get_wtime();
    const double dt1 = get_wtime_diff(&st1, &ft1);

    for (double sbound = startbound; sbound < testbound - startbound; sbound += boundgrid) {
        // speak("sbound", sbound); 
        double val = integral(distPerp / D, 0, sbound, M, ell0); 
        double inval[] = {distPerp / D, val * D}; 
        // speakvec(inval,2);
        double a1 = theBaobzi.evalFunc(inval); // calculate the Baobzi's upper limit of integral
        double bberr = ABS(a1 - sbound * D); 
        // speak("Baobzi Error", bberr); 
        baobzierr.push_back(bberr); 
        CHECK(a1 == Approx(sbound * D).epsilon(tol)); 

    }

    speak("Average Error for Reverse LookUP", mean_error(rlerr));
    speak("Randomly Chosen Error for Chebyshev", baobzierr[10]); 
    speak("Average Error for Chebyshev", mean_error(baobzierr)); 
    speak("Elapsed Time(s) for Chebyshev", dt1);

    // myfile << alpha << "," << bbtol << "," << dt1 << "," << mean_error(baobzierr) << "," << mean_error(rlerr) << std::endl;
    // myfile.close(); 
    // }
}

/** 
 * Spring with different stiffness and various parameter settings tests passed!  
 * Tests may not be inclusive and exhaustive enough -> subject to further examine
 * Problem: D > certain threshold -> need more adaptive formulation on prefactor 
 * Importance of wisely setting up domain -> Equation (47) in Lamson 2021 EPJ -> Currently easy when distPerp is fixed
*/

// TEST_CASE("Lookup table test (all kind) spring ", "[lookup_soft]") {
//     std::cout << "==== (ALL KIND) SPRING TEST ====" << std::endl; 

//     // Physical Parameters Setting
//     constexpr double errTol = 1e-3;
//     const double D = 0.024;
//     const double alpha = 10 / (2 * 0.00411);
//     const double freelength = 0.05;
//     const double M = alpha * D * D; speak("M",M); 
//     const double ell0 = freelength / D;

//     LUTFillerEdep lut_filler(256, 256);
//     lut_filler.Init(alpha, freelength, D);
//     LookupTable LUT(&lut_filler);

//     double distPerp = 0;
//     distPerp = 0.2; 
//     std::string rootpath = "int-res/";
//     std::string strPerp = std::to_string(distPerp);
//     std::ofstream myfile; 
//     std::string searchfilename = rootpath + strPerp + ".txt"; 
//     try {
//         std::filesystem::remove(searchfilename);
//     }
//     catch (...) {}

//     myfile.open(searchfilename);

//     double startbound = 0.1;
//     double testbound = 30; 
//     double boundgrid = 0.01; 
//     size_t gridcnt = floor((testbound - startbound) / boundgrid);
//     speak("total cases: ", gridcnt); 

//     std::vector<double> ludiff; std::vector<double> bbdiff; // error storer
//     std::vector<double> bbresl; // Baobzi calculation storer
//     std::vector<double> bbparm; // parameter storer

//     // Baobzi Parameter Setting
//     int dim = 2;
//     int order = 10; 
//     int odim = 1;
//     double bbtol = 1e-5;
//     double mlf = 0.0;
//     int sme = 1;
//     int mind = 0; 
//     int maxd = 40;
//     double hl[] = {1e-4, testbound / 2 * D}; // half length
//     double center[] = {distPerp / D, testbound / 2 * D};  // center
//     const char* fn = "func_approx.baobzi"; 

//     // Baobzi function approximator
//     Cheb theBaobzi(hl[0],hl[1],center[0],center[1],dim,odim,order,bbtol,mlf,sme,mind,maxd,M,ell0,D,fn);
//     theBaobzi.approxFunc(1);

//     // LOOKUP TABLE TEST
//     const auto st1 = get_wtime();
//     // ("distPerp = 0.2 > D+ell0, single peaked")
//     for (double sbound = startbound; sbound < testbound; sbound += boundgrid) {
//         double a1 = LUT.Lookup(distPerp, sbound * D);
//         double a2 = D * integral(distPerp / D, 0, sbound, M, ell0);
//         // CHECK(a1 == Approx(a2).epsilon(errTol));
//         ludiff.push_back(ABS(a1 - a2)); 
//     }
//     const auto ft1 = get_wtime();
//     const double dt1 = get_wtime_diff(&st1, &ft1);

//     // BAOBZI TEST
//     int cnt2 = 0; 
//     const auto st2 = get_wtime();
//     for (double sbound = startbound; sbound < testbound; sbound += boundgrid) {
//         // Baobzi test
//         double inval[] = {distPerp / D, sbound * D};
//         double a1 = theBaobzi.evalFunc(inval);  // baobzi result
//         double a2 = D * integral(distPerp / D, 0, sbound, M, ell0);  // integral for comparison
//         myfile << a1 << "," << sbound << std::endl;
//         CHECK(a1 == Approx(a2).epsilon(errTol));
//         // speak("Baobzi Error", ABS(a1 - a2)); 
//         bbdiff.push_back(ABS(a1 - a2)); 
//         bbresl.push_back(a1); 
//         bbparm.push_back(sbound); 
//     }
//     const auto ft2 = get_wtime();
//     const double dt2 = get_wtime_diff(&st2, &ft2);

//     speak("Average Error for Lookup", mean_error(ludiff));
//     speak("Average Error for Chebyshev", mean_error(bbdiff)); 
//     speak("Elapsed Time(s) for Lookup", dt1);
//     speak("Elapsed Time(s) for Chebyshev", dt2);

//     myfile.close(); 

//     // speakvec(bbresl, gridcnt); 
//     // speakvec(bbparm, gridcnt); 
// }

// /**
//  * binding volume test case passed! 
//  * The test 'Not using binding volume' has weird result (due to fixed standard), and not sure whether necessary to test on 
// */

// TEST_CASE("Test the calculation of binding volume.", "[bind volume]") {
//     std::cout << "==== BINDING VOLUME TEST ====" << std::endl; 

//     // Physical Parameters Setting
//     const double D = 0.024;
//     const double freelength = 0.05;
//     double alpha = 1. / (2 * 0.00411);
//     const double ell0 = freelength / D;
//     double M = alpha * D * D;
    
//     static constexpr double small_ = 1e-4;

//     double startbound = 0.1 * M;
//     double testbound = 1 * M; 
//     double boundgrid = 0.1 * M; 
//     int gridcnt = floor((testbound - startbound) / boundgrid);
//     speak("total cases: ", gridcnt); 

//     int dim = 2;
//     int order = 10; 
//     int odim = 1;
//     double bbtol = 1e-6;
//     double mlf = 0.0;
//     int sme = 1;
//     int mind = 0; 
//     int maxd = 40;
//     double hl[] = {0.5, testbound/2}; 
//     double center[] = {0.5, testbound/2};
//     const char* fn = "func_approx.baobzi"; 

//     Cheb theBaobzi(hl[0],hl[1],center[0],center[1],dim,odim,order,bbtol,mlf,sme,mind,maxd,M,ell0,D,fn);
//     theBaobzi.approxFunc(2);

//     LUTFillerEdep lut_filler(256, 256);
//     SECTION("Check LUTfiller for proper binding volume") {
//         for (double i = 0.1; i < 1.0; i += .1) {
//             double M1 = i * M; 
//             double uppbound = sqrt(-log(small_) / M1) + ell0; 
//             lut_filler.Init(alpha * i, freelength, D);
//             double bind_vol = bind_vol_integral(lut_filler.getUpperBound(), M1, ell0);
//             // speak("uppbound*D", uppbound  * D); 
//             // speak("i*M", M1); 
//             double inval[] = {uppbound * D, M1}; 
//             double bb_bind_vol = theBaobzi.evalFunc(inval); 
//             speak("Baobzi Error", ABS(bind_vol - bb_bind_vol)); 
//             CHECK(bb_bind_vol == Approx(bind_vol).epsilon(1e-8));
//         }
//     }


//     // SECTION("Check LUT for proper binding volume") {
//     //     for (double i = 0.1; i < 1.0; i += .1) {
//     //         lut_filler.Init(alpha * i, freelength, D);
//     //         LookupTable LUT(&lut_filler, true);
//     //         double bind_vol =
//     //             CUBE(D) *
//     //             bind_vol_integral(lut_filler.getUpperBound(), i * M, ell0);
//     //         CHECK(LUT.getBindVolume() == Approx(bind_vol).epsilon(1e-8));
//     //     }
//     // }


//     SECTION("Not using binding volume") {
//         lut_filler.Init(.5, freelength, D);
//         LookupTable LUT(&lut_filler);

//         double lookupval = LUT.getBindVolume();
//         double stdard = 1.; 

//         CHECK((lookupval) == Approx(stdard).epsilon(1e-12));

//         double convi = .5 / alpha;
//         double M1 = convi * M;
//         double comp = CUBE(D) * bind_vol_integral(lut_filler.getUpperBound(), M1, ell0);
//         speak("Integral Comp Result", comp); 

//         double uppbound = sqrt(-log(small_) / M1) + ell0;;
//         double inval1[] = {uppbound * D, M1};
//         speak("test uppbound",uppbound * D);
//         speak("test M1", M1); 
//         double hl1[] = {1e-3, M1/2}; 
//         double center1[] = {uppbound * D, M1};
//         Cheb theBaobzi1(hl1[0],hl1[1],center1[0],center1[1],dim,odim,order,bbtol,mlf,sme,mind,maxd,M,ell0,D,fn);
//         theBaobzi1.approxFunc(2);

//         double bb_bind_vol = CUBE(D) * theBaobzi1.evalFunc(inval1);
//         speak("Baobzi Result", bb_bind_vol);
//         speak("Baobzi Error", ABS(stdard - bb_bind_vol)); 
//     }
// }


// TEST_CASE("Lookup table test MEDIUM spring ", "[lookup_med]") {
//     constexpr double errTol = 1e-3;
//     const double D = 0.024;
//     const double alpha = 1.0 / (2 * 0.00411);
//     const double freelength = 0.05;
//     const double M = alpha * D * D;
//     const double ell0 = freelength / D;

//     LUTFillerEdep lut_filler(256, 256);
//     lut_filler.Init(alpha, freelength, D);
//     LookupTable LUT(&lut_filler);

//     double distPerp = 0;
//     // ("distPerp = 0.2 > D+ell0, single peaked")
//     distPerp = 0.2;
//     for (double sbound = 0; sbound < 20; sbound += 0.5) {
//         CHECK(LUT.Lookup(distPerp, sbound * D) ==
//               Approx(D * integral(distPerp / D, 0, sbound, M, ell0))
//                   .epsilon(errTol));
//     }
//     // ("distPerp = 0.1 > D+ell0, single peaked")
//     distPerp = 0.1;
//     for (double sbound = 0; sbound < 20; sbound += 0.5) {
//         CHECK(LUT.Lookup(distPerp, sbound * D) ==
//               Approx(D * integral(distPerp / D, 0, sbound, M, ell0))
//                   .epsilon(errTol));
//     }
//     // ("distPerp = 0.06 < D+ell0, double peaked")
//     distPerp = 0.06;
//     for (double sbound = 0; sbound < 20; sbound += 0.5) {
//         CHECK(LUT.Lookup(distPerp, sbound * D) ==
//               Approx(D * integral(distPerp / D, 0, sbound, M, ell0))
//                   .epsilon(errTol));
//     }
// }

// TEST_CASE("Lookup table test STIFF spring ", "[lookup_stiff]") {
//     constexpr double absTol = 1e-5;
//     const double D = 0.024;
//     const double alpha = 10.0 / (2 * 0.00411);
//     const double freelength = 0.05 + D;
//     const double M = alpha * D * D;
//     const double ell0 = freelength / D;

//     LUTFillerEdep lut_filler(256, 256);
//     lut_filler.Init(alpha, freelength, D);
//     LookupTable LUT(&lut_filler);

//     double distPerp = 0;
//     // ("distPerp = 0.2 > D+ell0, single peaked")
//     distPerp = 0.2;
//     for (double sbound = 0; sbound < 20; sbound += 0.5) {
//         // CHECK(errorPass(LUT.Lookup(distPerp, sbound * D),
//         // D * integral(distPerp / D, 0, sbound, M, ell0)));
//         CHECK(LUT.Lookup(distPerp, sbound * D) ==
//               Approx(D * integral(distPerp / D, 0, sbound, M, ell0))
//                   .margin(absTol));
//     }
//     // ("distPerp = 0.1 > D+ell0, single peaked")
//     distPerp = 0.1;
//     for (double sbound = 0; sbound < 20; sbound += 0.5) {
//         // CHECK(errorPass(LUT.Lookup(distPerp, sbound * D),
//         // D * integral(distPerp / D, 0, sbound, M, ell0)));
//         CHECK(LUT.Lookup(distPerp, sbound * D) ==
//               Approx(D * integral(distPerp / D, 0, sbound, M, ell0))
//                   .margin(absTol));
//     }
//     // ("distPerp = 0.06 < D+ell0, double peaked")
//     distPerp = 0.06;
//     for (double sbound = 0; sbound < 20; sbound += 0.5) {
//         CHECK(errorPass(LUT.Lookup(distPerp, sbound * D),
//                         D * integral(distPerp / D, 0, sbound, M, ell0)));
//         CHECK(LUT.Lookup(distPerp, sbound * D) ==
//               Approx(D * integral(distPerp / D, 0, sbound, M, ell0))
//                   .margin(absTol));
//     }
// }

// TEST_CASE("Lookup table test manual medium spring REL error", "[lookup]") {
//     // integrated by mathematica
//     const double D = 0.024;
//     constexpr double tol = 1e-4;

//     LUTFillerEdep lut_filler(256, 256);
//     lut_filler.Init(1.0 / (2 * 0.00411), 0.05 + D, D);
//     LookupTable LUT(&lut_filler);

//     double distPerp = 0;

//     distPerp = 0.2;
//     // ("distPerp = 0.2 > D+ell0, single peaked")
//     CHECK(LUT.Lookup(distPerp, 0.5 * D) / D == Approx(0.0722077).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 1.0 * D) / D == Approx(0.1428390).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 1.5 * D) / D == Approx(0.2104120).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 2.0 * D) / D == Approx(0.2736230).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 3.0 * D) / D == Approx(0.3830390).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 4.0 * D) / D == Approx(0.4663750).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 5.0 * D) / D == Approx(0.5238890).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 6.0 * D) / D == Approx(0.5596700).epsilon(tol));

//     distPerp = 0.08;
//     // "distPerp = 0.08 > D+ell0, single peaked"/D,
//     CHECK(LUT.Lookup(distPerp, 0.5 * D) / D == Approx(0.497588).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 1.0 * D) / D == Approx(0.993608).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 1.5 * D) / D == Approx(1.485540).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 2.0 * D) / D == Approx(1.969290).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 3.0 * D) / D == Approx(2.887840).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 4.0 * D) / D == Approx(3.692500).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 5.0 * D) / D == Approx(4.333200).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 6.0 * D) / D == Approx(4.789860).epsilon(tol));

//     distPerp = 0.06;
//     // "distPerp = 0.06 < D+ell0, double peaked"/D,
//     CHECK(LUT.Lookup(distPerp, 0.5 * D) / D == Approx(0.488864).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 1.0 * D) / D == Approx(0.981139).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 1.5 * D) / D == Approx(1.478150).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 2.0 * D) / D == Approx(1.977880).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 3.0 * D) / D == Approx(2.960520).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 4.0 * D) / D == Approx(3.858570).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 5.0 * D) / D == Approx(4.598640).epsilon(tol));
//     CHECK(LUT.Lookup(distPerp, 6.0 * D) / D == Approx(5.140580).epsilon(tol));
// }

// TEST_CASE("REVERSE Lookup table test manual medium spring REL error",
//           "[REVERSE lookup]") {
//     // integrated by mathematica
//     const double D = 0.024;

//     double distPerp = 0;
//     LUTFillerEdep lut_filler(256, 256);
//     lut_filler.Init(1.0 / (2 * 0.00411), 0.05 + D, D);
//     LookupTable LUT(&lut_filler);
//     // LUT.Init(&lut_filler);

//     // double tol = RELTOL * REVERSEFAC;
//     const double tol = 1e-4;

//     distPerp = 0.1;
//     // ("distPerp = 0.1 > D+ell0, single peaked")
//     CHECK(LUT.ReverseLookup(distPerp, D * 0) / D == Approx(0.0).epsilon(tol));
//     CHECK(LUT.ReverseLookup(distPerp, D * 0.0460519) / D ==
//           Approx(.05).epsilon(tol));
//     CHECK(LUT.ReverseLookup(distPerp, D * 0.3221280) / D ==
//           Approx(.35).epsilon(tol));
//     CHECK(LUT.ReverseLookup(distPerp, D * 0.4598240) / D ==
//           Approx(0.5).epsilon(tol));
//     CHECK(LUT.ReverseLookup(distPerp, D * 0.9153560) / D ==
//           Approx(1.0).epsilon(tol));
//     CHECK(LUT.ReverseLookup(distPerp, D * 1.3619600) / D ==
//           Approx(1.5).epsilon(tol));
//     CHECK(LUT.ReverseLookup(distPerp, D * 1.7944600) / D ==
//           Approx(2.0).epsilon(tol));
//     CHECK(LUT.ReverseLookup(distPerp, D * 2.2071800) / D ==
//           Approx(2.5).epsilon(tol));
//     CHECK(LUT.ReverseLookup(distPerp, D * 3.2701500) / D ==
//           Approx(4.0).epsilon(tol));
//     CHECK(LUT.ReverseLookup(distPerp, D * 3.7911500) / D ==
//           Approx(5.0).epsilon(tol));
//     CHECK(LUT.ReverseLookup(distPerp, D * 4.1524200) / D ==
//           Approx(6.0).epsilon(tol));
//     // CHECK(relError(LUT.ReverseLookup(distPerp / D, 4.37561), 7.0) < tol);
// }

// TEST_CASE("REVERSE Lookup table test stiff spring ", "[REVERSE lookup]") {
//     const double tol = 1e-2;

//     const double D = 0.024;
//     const double alpha = 10.0 / (2 * 0.00411);
//     const double freelength = 0.05;
//     const double M = alpha * D * D;
//     const double ell0 = freelength / D;

//     LUTFillerEdep lut_filler(256, 256);
//     lut_filler.Init(alpha, freelength, D);
//     LookupTable LUT(&lut_filler);

//     double distPerp = 0;
//     // ("distPerp = 0.2 > D+ell0, single peaked")
//     // WARNING: This reverse lookup fails because function is too flat
//     // distPerp = 0.2;
//     // for (double sbound = 0; sbound < LUT.getNonDsbound() / 8; sbound += 0.1)
//     // {
//     //    double val = integral(distPerp / D, 0, sbound, M, ell0);
//     //    double scalc = LUT.ReverseLookup(distPerp, val * D);
//     //    printf("scalc = %f\n", scalc);
//     //    printf("sbound = %f\n", sbound * D);
//     //    CHECK(errorPass(scalc, sbound * D, REVERSEFAC));
//     //}

//     // ("distPerp = 0.1 > D+ell0, single peaked")
//     distPerp = 0.1;
//     for (double sbound = 0; sbound < LUT.getNonDsbound() / 2; sbound += 0.1) {
//         double val = integral(distPerp / D, 0, sbound, M, ell0);
//         CHECK(LUT.ReverseLookup(distPerp, val * D) ==
//               Approx(sbound * D).margin(tol));
//         // CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
//         // REVERSEFAC));
//     }
//     // ("distPerp = 0.06 < D+ell0, double peaked")
//     distPerp = 0.06;
//     for (double sbound = 0; sbound < LUT.getNonDsbound() / 2; sbound += 0.1) {
//         double val = integral(distPerp / D, 0, sbound, M, ell0);
//         CHECK(LUT.ReverseLookup(distPerp, val * D) ==
//               Approx(sbound * D).epsilon(tol));
//         // CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
//         // REVERSEFAC));
//     }
// }

// TEST_CASE("REVERSE binary Lookup with different springs", "[REVERSE binary]") {

//     const double D = 0.024;
//     const double freelength = 0.05;
//     const double ell0 = freelength / D;
//     double alpha, M;

//     LUTFillerEdep lut_filler(256, 256);

//     double rowIndexMax = lut_filler.getDistPerpGridNum();
//     double colIndexMax = lut_filler.getDistParaGridNum();
//     double distPerpSpacing;

//     SECTION("Soft spring") {
//         alpha = 0.1 / (2 * 0.00411);
//         M = alpha * D * D;

//         lut_filler.Init(alpha, freelength, D);
//         LookupTable LUT(&lut_filler);

//         distPerpSpacing = LUT.perp_spacing_;
//         for (int i = 0; i < rowIndexMax - 2; ++i) {
//             double C0 = LUT.table_[LUT.getTableIndex(i, colIndexMax - 1)] * D;
//             double C1 =
//                 LUT.table_[LUT.getTableIndex(i + 1, colIndexMax - 1)] * D;
//             double Cavg = .5 * (C0 + C1);
//             double distPerpAvg = distPerpSpacing * (i + .5) * D;
//             double sbound = LUT.ReverseLookup(distPerpAvg, Cavg);
//             double Cintegral =
//                 integral(distPerpAvg / D, 0, sbound / D, M, ell0) * D;
//             CHECK(Cintegral == Approx(Cavg).epsilon(RELTOL));
//         }
//     }

//     SECTION("Medium spring") {
//         alpha = 1.0 / (2 * 0.00411);
//         M = alpha * D * D;

//         lut_filler.Init(alpha, freelength, D);
//         LookupTable LUT(&lut_filler);

//         distPerpSpacing = LUT.perp_spacing_;
//         for (int i = 0; i < rowIndexMax - 2; ++i) {
//             double C0 = LUT.table_[LUT.getTableIndex(i, colIndexMax - 1)] * D;
//             double C1 =
//                 LUT.table_[LUT.getTableIndex(i + 1, colIndexMax - 1)] * D;
//             double Cavg = .5 * (C0 + C1);
//             double distPerpAvg = distPerpSpacing * (i + .5) * D;
//             double sbound = LUT.ReverseLookup(distPerpAvg, Cavg);
//             double Cintegral =
//                 integral(distPerpAvg / D, 0, sbound / D, M, ell0) * D;
//             CHECK(Cintegral == Approx(Cavg).epsilon(RELTOL));
//         }
//     }

//     SECTION("Stiff spring") {
//         alpha = 10. / (2 * 0.00411);
//         M = alpha * D * D;

//         lut_filler.Init(alpha, freelength, D);
//         LookupTable LUT(&lut_filler);

//         distPerpSpacing = LUT.perp_spacing_;
//         for (int i = 0; i < rowIndexMax - 2; ++i) {
//             double C0 = LUT.table_[LUT.getTableIndex(i, colIndexMax - 1)] * D;
//             double C1 =
//                 LUT.table_[LUT.getTableIndex(i + 1, colIndexMax - 1)] * D;
//             double Cavg = .5 * (C0 + C1);
//             double distPerpAvg = distPerpSpacing * (i + .5) * D;
//             double sbound = LUT.ReverseLookup(distPerpAvg, Cavg);
//             double Cintegral =
//                 integral(distPerpAvg / D, 0, sbound / D, M, ell0) * D;
//             CHECK(Cintegral == Approx(Cavg).margin(1e-5));
//         }
//     }
// }
