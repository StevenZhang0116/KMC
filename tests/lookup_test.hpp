#include "KMC_func/integrals.hpp"
#include "KMC_func/lookup_table.hpp"
#include "KMC_func/lut_filler_edep.hpp"
#include "KMC_func/macros.hpp"
#include "catch.hpp"
#include "test_helpers.hpp"
#include <baobzi_template.hpp>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

#include <boost/math/quadrature/gauss_kronrod.hpp>

TEST_CASE("Lookup table test SOFT spring ", "[lookup_soft]") {
    constexpr double errTol = 1e-3;
    const double D = 0.024;
    const double alpha = 0.1 / (2 * 0.00411);
    const double freelength = 0.05;
    const double M = alpha * D * D;
    const double ell0 = freelength / D;

    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init(alpha, freelength, D);
    LookupTable LUT(&lut_filler);

    double distPerp = 0;
    distPerp = 0.01;
    // ("distPerp = 0.2 > D+ell0, single peaked")
    for (double sbound = 0; sbound < 20; sbound += 0.5) {
        CHECK(LUT.Lookup(distPerp, sbound * D) ==
              Approx(D * integral(distPerp / D, 0, sbound, M, ell0))
                  .epsilon(errTol));
        // test err
        double luval = LUT.Lookup(distPerp, sbound * D); 
        double caval = D * integral(distPerp / D, 0, sbound, M, ell0); 
        std::cout << ABS(luval-caval) << std::endl; 
        
    }

    // ("distPerp = 0.1 > D+ell0, single peaked")
    distPerp = 0.1;
    for (double sbound = 0; sbound < 20; sbound += 0.5) {
        CHECK(LUT.Lookup(distPerp, sbound * D) ==
              Approx(D * integral(distPerp / D, 0, sbound, M, ell0))
                  .epsilon(errTol));
    }

    // ("distPerp = 0.06 < D+ell0, double peaked")
    distPerp = 0.06;
    for (double sbound = 0; sbound < 20; sbound += 0.5) {
        CHECK(LUT.Lookup(distPerp, sbound * D) ==
              Approx(D * integral(distPerp / D, 0, sbound, M, ell0))
                  .epsilon(errTol));
    }
}

TEST_CASE("Lookup table test MEDIUM spring ", "[lookup_med]") {
    constexpr double errTol = 1e-3;
    const double D = 0.024;
    const double alpha = 1.0 / (2 * 0.00411);
    const double freelength = 0.05;
    const double M = alpha * D * D;
    const double ell0 = freelength / D;

    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init(alpha, freelength, D);
    LookupTable LUT(&lut_filler);

    double distPerp = 0;
    // ("distPerp = 0.2 > D+ell0, single peaked")
    distPerp = 0.2;
    for (double sbound = 0; sbound < 20; sbound += 0.5) {
        CHECK(LUT.Lookup(distPerp, sbound * D) ==
              Approx(D * integral(distPerp / D, 0, sbound, M, ell0))
                  .epsilon(errTol));
    }
    // ("distPerp = 0.1 > D+ell0, single peaked")
    distPerp = 0.1;
    for (double sbound = 0; sbound < 20; sbound += 0.5) {
        CHECK(LUT.Lookup(distPerp, sbound * D) ==
              Approx(D * integral(distPerp / D, 0, sbound, M, ell0))
                  .epsilon(errTol));
    }
    // ("distPerp = 0.06 < D+ell0, double peaked")
    distPerp = 0.06;
    for (double sbound = 0; sbound < 20; sbound += 0.5) {
        CHECK(LUT.Lookup(distPerp, sbound * D) ==
              Approx(D * integral(distPerp / D, 0, sbound, M, ell0))
                  .epsilon(errTol));
    }
}

TEST_CASE("Lookup table test STIFF spring ", "[lookup_stiff]") {
    constexpr double absTol = 1e-5;
    const double D = 0.024;
    const double alpha = 10.0 / (2 * 0.00411);
    const double freelength = 0.05 + D;
    const double M = alpha * D * D;
    const double ell0 = freelength / D;

    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init(alpha, freelength, D);
    LookupTable LUT(&lut_filler);

    double distPerp = 0;
    // ("distPerp = 0.2 > D+ell0, single peaked")
    distPerp = 0.2;
    for (double sbound = 0; sbound < 20; sbound += 0.5) {
        // CHECK(errorPass(LUT.Lookup(distPerp, sbound * D),
        // D * integral(distPerp / D, 0, sbound, M, ell0)));
        CHECK(LUT.Lookup(distPerp, sbound * D) ==
              Approx(D * integral(distPerp / D, 0, sbound, M, ell0))
                  .margin(absTol));
    }
    // ("distPerp = 0.1 > D+ell0, single peaked")
    distPerp = 0.1;
    for (double sbound = 0; sbound < 20; sbound += 0.5) {
        // CHECK(errorPass(LUT.Lookup(distPerp, sbound * D),
        // D * integral(distPerp / D, 0, sbound, M, ell0)));
        CHECK(LUT.Lookup(distPerp, sbound * D) ==
              Approx(D * integral(distPerp / D, 0, sbound, M, ell0))
                  .margin(absTol));
    }
    // ("distPerp = 0.06 < D+ell0, double peaked")
    distPerp = 0.06;
    for (double sbound = 0; sbound < 20; sbound += 0.5) {
        CHECK(errorPass(LUT.Lookup(distPerp, sbound * D),
                        D * integral(distPerp / D, 0, sbound, M, ell0)));
        CHECK(LUT.Lookup(distPerp, sbound * D) ==
              Approx(D * integral(distPerp / D, 0, sbound, M, ell0))
                  .margin(absTol));
    }
}

TEST_CASE("Lookup table test manual medium spring REL error", "[lookup]") {
    // integrated by mathematica
    const double D = 0.024;
    constexpr double tol = 1e-4;

    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init(1.0 / (2 * 0.00411), 0.05 + D, D);
    LookupTable LUT(&lut_filler);

    double distPerp = 0;

    distPerp = 0.2;
    // ("distPerp = 0.2 > D+ell0, single peaked")
    CHECK(LUT.Lookup(distPerp, 0.5 * D) / D == Approx(0.0722077).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 1.0 * D) / D == Approx(0.1428390).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 1.5 * D) / D == Approx(0.2104120).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 2.0 * D) / D == Approx(0.2736230).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 3.0 * D) / D == Approx(0.3830390).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 4.0 * D) / D == Approx(0.4663750).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 5.0 * D) / D == Approx(0.5238890).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 6.0 * D) / D == Approx(0.5596700).epsilon(tol));

    distPerp = 0.08;
    // "distPerp = 0.08 > D+ell0, single peaked"/D,
    CHECK(LUT.Lookup(distPerp, 0.5 * D) / D == Approx(0.497588).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 1.0 * D) / D == Approx(0.993608).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 1.5 * D) / D == Approx(1.485540).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 2.0 * D) / D == Approx(1.969290).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 3.0 * D) / D == Approx(2.887840).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 4.0 * D) / D == Approx(3.692500).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 5.0 * D) / D == Approx(4.333200).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 6.0 * D) / D == Approx(4.789860).epsilon(tol));

    distPerp = 0.06;
    // "distPerp = 0.06 < D+ell0, double peaked"/D,
    CHECK(LUT.Lookup(distPerp, 0.5 * D) / D == Approx(0.488864).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 1.0 * D) / D == Approx(0.981139).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 1.5 * D) / D == Approx(1.478150).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 2.0 * D) / D == Approx(1.977880).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 3.0 * D) / D == Approx(2.960520).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 4.0 * D) / D == Approx(3.858570).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 5.0 * D) / D == Approx(4.598640).epsilon(tol));
    CHECK(LUT.Lookup(distPerp, 6.0 * D) / D == Approx(5.140580).epsilon(tol));
}

TEST_CASE("REVERSE Lookup table test manual medium spring REL error",
          "[REVERSE lookup]") {
    // integrated by mathematica
    const double D = 0.024;

    double distPerp = 0;
    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init(1.0 / (2 * 0.00411), 0.05 + D, D);
    LookupTable LUT(&lut_filler);
    // LUT.Init(&lut_filler);

    // double tol = RELTOL * REVERSEFAC;
    const double tol = 1e-4;

    distPerp = 0.1;
    // ("distPerp = 0.1 > D+ell0, single peaked")
    CHECK(LUT.ReverseLookup(distPerp, D * 0) / D == Approx(0.0).epsilon(tol));
    CHECK(LUT.ReverseLookup(distPerp, D * 0.0460519) / D ==
          Approx(.05).epsilon(tol));
    CHECK(LUT.ReverseLookup(distPerp, D * 0.3221280) / D ==
          Approx(.35).epsilon(tol));
    CHECK(LUT.ReverseLookup(distPerp, D * 0.4598240) / D ==
          Approx(0.5).epsilon(tol));
    CHECK(LUT.ReverseLookup(distPerp, D * 0.9153560) / D ==
          Approx(1.0).epsilon(tol));
    CHECK(LUT.ReverseLookup(distPerp, D * 1.3619600) / D ==
          Approx(1.5).epsilon(tol));
    CHECK(LUT.ReverseLookup(distPerp, D * 1.7944600) / D ==
          Approx(2.0).epsilon(tol));
    CHECK(LUT.ReverseLookup(distPerp, D * 2.2071800) / D ==
          Approx(2.5).epsilon(tol));
    CHECK(LUT.ReverseLookup(distPerp, D * 3.2701500) / D ==
          Approx(4.0).epsilon(tol));
    CHECK(LUT.ReverseLookup(distPerp, D * 3.7911500) / D ==
          Approx(5.0).epsilon(tol));
    CHECK(LUT.ReverseLookup(distPerp, D * 4.1524200) / D ==
          Approx(6.0).epsilon(tol));
    // CHECK(relError(LUT.ReverseLookup(distPerp / D, 4.37561), 7.0) < tol);
}

TEST_CASE("REVERSE Lookup table test soft spring ", "[REVERSE lookup]") {
    const double tol = 1e-2;
    const double D = 0.024;
    const double alpha = 0.1 / (2 * 0.00411);
    const double freelength = 0.05;
    const double M = alpha * D * D;
    const double ell0 = freelength / D;

    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init(alpha, freelength, D);
    LookupTable LUT(&lut_filler);

    double distPerp = 0;
    distPerp = 0.2;
    // ("distPerp = 0.2 > D+ell0, single peaked")
    for (double sbound = 0; sbound < LUT.getNonDsbound() / 2; sbound += 0.2) {
        double val = integral(distPerp / D, 0, sbound, M, ell0);
        CHECK(LUT.ReverseLookup(distPerp, val * D) ==
              Approx(sbound * D).epsilon(tol));
    }
    // ("distPerp = 0.1 > D+ell0, single peaked")
    distPerp = 0.1;
    for (double sbound = 0; sbound < LUT.getNonDsbound() / 2; sbound += 0.2) {
        double val = integral(distPerp / D, 0, sbound, M, ell0);
        CHECK(LUT.ReverseLookup(distPerp, val * D) ==
              Approx(sbound * D).epsilon(tol));
    }
    // ("distPerp = 0.06 < D+ell0, double peaked")
    distPerp = 0.06;
    for (double sbound = 0; sbound < LUT.getNonDsbound() / 2; sbound += 0.2) {
        double val = integral(distPerp / D, 0, sbound, M, ell0);
        CHECK(LUT.ReverseLookup(distPerp, val * D) ==
              Approx(sbound * D).epsilon(tol));
    }
}

TEST_CASE("REVERSE Lookup table test medium spring ", "[REVERSE lookup]") {
    const double tol = 1e-2;
    const double D = 0.024;
    const double alpha = 1.0 / (2 * 0.00411);
    const double freelength = 0.05;
    const double M = alpha * D * D;
    const double ell0 = freelength / D;

    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init(alpha, freelength, D);
    LookupTable LUT(&lut_filler);

    double distPerp = 0;
    distPerp = 0.2;
    // ("distPerp = 0.2 > D+ell0, single peaked")
    for (double sbound = 0; sbound < LUT.getNonDsbound() / 3; sbound += 0.2) {
        double val = integral(distPerp / D, 0, sbound, M, ell0);
        // CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D));
        CHECK(LUT.ReverseLookup(distPerp, val * D) ==
              Approx(sbound * D).epsilon(tol));
    }
    // ("distPerp = 0.1 > D+ell0, single peaked")
    distPerp = 0.1;
    for (double sbound = 0; sbound < LUT.getNonDsbound() / 2; sbound += 0.2) {
        double val = integral(distPerp / D, 0, sbound, M, ell0);
        // CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D));
        CHECK(LUT.ReverseLookup(distPerp, val * D) ==
              Approx(sbound * D).epsilon(tol));
    }
    // ("distPerp = 0.06 < D+ell0, double peaked")
    distPerp = 0.06;
    for (double sbound = 0; sbound < LUT.getNonDsbound() / 2; sbound += 0.2) {
        double val = integral(distPerp / D, 0, sbound, M, ell0);
        // CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D));
        CHECK(LUT.ReverseLookup(distPerp, val * D) ==
              Approx(sbound * D).epsilon(tol));
    }
}

TEST_CASE("REVERSE Lookup table test stiff spring ", "[REVERSE lookup]") {
    const double tol = 1e-2;

    const double D = 0.024;
    const double alpha = 10.0 / (2 * 0.00411);
    const double freelength = 0.05;
    const double M = alpha * D * D;
    const double ell0 = freelength / D;

    LUTFillerEdep lut_filler(256, 256);
    lut_filler.Init(alpha, freelength, D);
    LookupTable LUT(&lut_filler);

    double distPerp = 0;
    // ("distPerp = 0.2 > D+ell0, single peaked")
    // WARNING: This reverse lookup fails because function is too flat
    // distPerp = 0.2;
    // for (double sbound = 0; sbound < LUT.getNonDsbound() / 8; sbound += 0.1)
    // {
    //    double val = integral(distPerp / D, 0, sbound, M, ell0);
    //    double scalc = LUT.ReverseLookup(distPerp, val * D);
    //    printf("scalc = %f\n", scalc);
    //    printf("sbound = %f\n", sbound * D);
    //    CHECK(errorPass(scalc, sbound * D, REVERSEFAC));
    //}

    // ("distPerp = 0.1 > D+ell0, single peaked")
    distPerp = 0.1;
    for (double sbound = 0; sbound < LUT.getNonDsbound() / 2; sbound += 0.1) {
        double val = integral(distPerp / D, 0, sbound, M, ell0);
        CHECK(LUT.ReverseLookup(distPerp, val * D) ==
              Approx(sbound * D).margin(tol));
        // CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
        // REVERSEFAC));
    }
    // ("distPerp = 0.06 < D+ell0, double peaked")
    distPerp = 0.06;
    for (double sbound = 0; sbound < LUT.getNonDsbound() / 2; sbound += 0.1) {
        double val = integral(distPerp / D, 0, sbound, M, ell0);
        CHECK(LUT.ReverseLookup(distPerp, val * D) ==
              Approx(sbound * D).epsilon(tol));
        // CHECK(errorPass(LUT.ReverseLookup(distPerp, val * D), sbound * D,
        // REVERSEFAC));
    }
}

TEST_CASE("REVERSE binary Lookup with different springs", "[REVERSE binary]") {

    const double D = 0.024;
    const double freelength = 0.05;
    const double ell0 = freelength / D;
    double alpha, M;

    LUTFillerEdep lut_filler(256, 256);

    double rowIndexMax = lut_filler.getDistPerpGridNum();
    double colIndexMax = lut_filler.getDistParaGridNum();
    double distPerpSpacing;

    SECTION("Soft spring") {
        alpha = 0.1 / (2 * 0.00411);
        M = alpha * D * D;

        lut_filler.Init(alpha, freelength, D);
        LookupTable LUT(&lut_filler);

        distPerpSpacing = LUT.perp_spacing_;
        for (int i = 0; i < rowIndexMax - 2; ++i) {
            double C0 = LUT.table_[LUT.getTableIndex(i, colIndexMax - 1)] * D;
            double C1 =
                LUT.table_[LUT.getTableIndex(i + 1, colIndexMax - 1)] * D;
            double Cavg = .5 * (C0 + C1);
            double distPerpAvg = distPerpSpacing * (i + .5) * D;
            double sbound = LUT.ReverseLookup(distPerpAvg, Cavg);
            double Cintegral =
                integral(distPerpAvg / D, 0, sbound / D, M, ell0) * D;
            CHECK(Cintegral == Approx(Cavg).epsilon(RELTOL));
        }
    }

    SECTION("Medium spring") {
        alpha = 1.0 / (2 * 0.00411);
        M = alpha * D * D;

        lut_filler.Init(alpha, freelength, D);
        LookupTable LUT(&lut_filler);

        distPerpSpacing = LUT.perp_spacing_;
        for (int i = 0; i < rowIndexMax - 2; ++i) {
            double C0 = LUT.table_[LUT.getTableIndex(i, colIndexMax - 1)] * D;
            double C1 =
                LUT.table_[LUT.getTableIndex(i + 1, colIndexMax - 1)] * D;
            double Cavg = .5 * (C0 + C1);
            double distPerpAvg = distPerpSpacing * (i + .5) * D;
            double sbound = LUT.ReverseLookup(distPerpAvg, Cavg);
            double Cintegral =
                integral(distPerpAvg / D, 0, sbound / D, M, ell0) * D;
            CHECK(Cintegral == Approx(Cavg).epsilon(RELTOL));
        }
    }

    SECTION("Stiff spring") {
        alpha = 10. / (2 * 0.00411);
        M = alpha * D * D;

        lut_filler.Init(alpha, freelength, D);
        LookupTable LUT(&lut_filler);

        distPerpSpacing = LUT.perp_spacing_;
        for (int i = 0; i < rowIndexMax - 2; ++i) {
            double C0 = LUT.table_[LUT.getTableIndex(i, colIndexMax - 1)] * D;
            double C1 =
                LUT.table_[LUT.getTableIndex(i + 1, colIndexMax - 1)] * D;
            double Cavg = .5 * (C0 + C1);
            double distPerpAvg = distPerpSpacing * (i + .5) * D;
            double sbound = LUT.ReverseLookup(distPerpAvg, Cavg);
            double Cintegral =
                integral(distPerpAvg / D, 0, sbound / D, M, ell0) * D;
            CHECK(Cintegral == Approx(Cavg).margin(1e-5));
        }
    }
}

TEST_CASE("Test the calculation of binding volume.", "[bind volume]") {
    const double D = 0.024;
    const double freelength = 0.05;
    double alpha = 1. / (2 * 0.00411);
    const double ell0 = freelength / D;
    double M = alpha * D * D;

    LUTFillerEdep lut_filler(256, 256);
    SECTION("Check LUTfiller for proper binding volume") {
        for (double i = 0.1; i < 1.0; i += .1) {
            lut_filler.Init(alpha * i, freelength, D);
            double bind_vol =
                bind_vol_integral(lut_filler.getUpperBound(), i * M, ell0);
            REQUIRE(lut_filler.getBindingVolume() ==
                    Approx(bind_vol).epsilon(1e-8));
        }
    }
    SECTION("Check LUT for proper binding volume") {
        for (double i = 0.1; i < 1.0; i += .1) {
            lut_filler.Init(alpha * i, freelength, D);
            LookupTable LUT(&lut_filler, true);
            double bind_vol =
                CUBE(D) *
                bind_vol_integral(lut_filler.getUpperBound(), i * M, ell0);
            REQUIRE(LUT.getBindVolume() == Approx(bind_vol).epsilon(1e-8));
        }
    }

    SECTION("Not using binding volume") {
        lut_filler.Init(.5, freelength, D);
        LookupTable LUT(&lut_filler);

        REQUIRE(LUT.getBindVolume() == Approx(1.).epsilon(1e-12));
    }
}

