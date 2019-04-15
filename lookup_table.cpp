#include "lookup_table.hpp"

#include <cassert>
#include <cstdio>
#include <iostream>
#include <string>

#include <boost/math/quadrature/gauss_kronrod.hpp>

/*! \brief Integrate exponential factor with the form of
 * e^{-M * (\sqrt{ s^2 + lm^2}-ell0 -1)^2}
 * from sbound0 to sbound1  with respect to the variable s.
 *
 * \param lm Physically, this is the perpendicular distance above rod
 * \param sbound lowerr limit of integral
 * \param sbound Upper limit of integral
 * \param M exponential constant factor. Physically, this is the product of
 (1-load_sensitivity)*spring_const/(k_B * Temperature)
 * \param ell0 Shift of the integrands mean. Physically, protein rest length
 * \return result The value of the integration

 */
double integral(double lm, double sbound0, double sbound1, double M,
                double ell0) {
    if (sbound0 >= sbound1) {
        return 0;
    }
    auto integrand = [&](double s) {
        // lambda capture variabls ell0 and M
        const double exponent = sqrt(s * s + lm * lm) - ell0 - 1;
        return exp(-M * exponent * exponent);
    };
    double error = 0;
    double result =
        boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
            integrand, sbound0, sbound1, 10, 1e-6, &error);
    return result;
}

void LookupTable::Init(double M_, double freeLength, double tubuleD) {

    // setup length scale and dimensionless lengths
    D = tubuleD;
    ell0 = freeLength / tubuleD;
    M = M_ * tubuleD * tubuleD;

    // truncate the integration when integrand < SMALL
    // Interation table in dimensionless lengths
    // lUB = sqrt(lm^2 + s^2), dimensionless scaled by D
    constexpr double SMALL = 1e-5;
    const double lUB = sqrt(-log(SMALL) / M) + 1 + ell0;

    // step 1 determine grid
    const double distPerpLB = 0; // allow min distPerp to be half D, some overlap
    const double distPerpUB = lUB;
    distPerpGridNumber = 48; // grid in dperp
    double distPerpGridSpacing =
        (distPerpUB - distPerpLB) / (distPerpGridNumber - 1);

    const double sboundLB = 0;
    const double sboundUB = sqrt(lUB * lUB - distPerpLB * distPerpLB);
    sboundGridNumber = 64; // grid in s bound
    double sboundGridSpacing = (sboundUB - sboundLB) / (sboundGridNumber - 1);

    // step 2 init grid
    distPerpGrid.resize(distPerpGridNumber);
    for (int i = 0; i < distPerpGridNumber; i++) {
        distPerpGrid[i] = i * distPerpGridSpacing + distPerpLB;
    }
    sboundGrid.resize(sboundGridNumber);
    for (int i = 0; i < sboundGridNumber; i++) {
        sboundGrid[i] = i * sboundGridSpacing + sboundLB;
    }
    distPerpGridSpacingInv = 1 / distPerpGridSpacing;
    sboundGridSpacingInv = 1 / sboundGridSpacing;

    // step 3 tabulate the matrix
    FillMatrix();
}

void LookupTable::FillMatrix() {
    table.resize(distPerpGridNumber * sboundGridNumber, 0);

    // boost integration parameters
    const int max_depth = 10; // maximum number of interval splittings
    const double tol = 1e-8;  // maximum relative error

    // row major
    // i is slow changing, should be distPerpGrid
    // j is fast changing, should be sboundGrid

    // IMPORTANT: Table stores dimensionless val
    for (int i = 0; i < distPerpGridNumber; i++) {
        for (int j = 0; j < sboundGridNumber; j++) {
            const double sbound = sboundGrid[j];
            assert(!(sbound < 0));
            double result = integral(distPerpGrid[i], 0, sbound, M, ell0);
            table[i * sboundGridNumber + j] = result;
        }
    }
}

void LookupTable::PrintTable() const {
    const int distPerpGridNumber = distPerpGrid.size();
    const int sboundGridNumber = sboundGrid.size();

    // IMPORTANT: Table stores dimensionless val
    for (int i = 0; i < distPerpGridNumber; i++) {
        for (int j = 0; j < sboundGridNumber; j++) {
            std::cout << "(distPerp, sbound) = (" << distPerpGrid[i] << ","
                      << sboundGrid[j]
                      << ") = " << table[i * sboundGridNumber + j] << std::endl;
        }
    }
}

double LookupTable::Lookup(double distPerp, double sbound) const {
    double rowFrac = 0, colFrac = 0;
    int rowIndex = getRowIndex(distPerp, rowFrac);
    int colIndex = getColIndex(sbound, colFrac);

    // clamp to grid bound
    if (rowIndex < 0) {
        printf("Error: distPerp too small: %g \n", distPerp);
        exit(1);
    }
    if (colIndex < 0) {
        printf("Error: sbound too small: %g \n", sbound);
        exit(1);
    }
    if (rowIndex > distPerpGridNumber - 2) {
#ifndef NDEBUG
        printf("Warning: distPerp %g very large, clamp to grid UB\n", distPerp);
#endif
        rowIndex = distPerpGridNumber - 2;
        rowFrac = 1;
    }
    if (colIndex > sboundGridNumber - 2) {
#ifndef NDEBUG
        printf("Warning: sbound %g very large, clamp to grid UB\n", sbound);
#endif
        colIndex = sboundGridNumber - 2;
        colFrac = 1;
    }

    double val = table[getTableIndex(rowIndex, colIndex)] * (1 - colFrac) *
                     (1 - rowFrac) //
                 + table[getTableIndex(rowIndex, colIndex + 1)] * colFrac *
                       (1 - rowFrac) //
                 + table[getTableIndex(rowIndex + 1, colIndex)] *
                       (1 - colFrac) * rowFrac //
                 + table[getTableIndex(rowIndex + 1, colIndex + 1)] * colFrac *
                       rowFrac; //

    return val; // IMPORTANT: table stores dimensionless values
}

double LookupTable::ReverseLookup(double distPerp, const double val) const {
    if (val == 0) {
        return 0;
    }

    double sbound;
    double rowFrac = 0;
    int rowIndex = getRowIndex(distPerp, rowFrac);

    // distPerp check
    if (rowIndex < 0) {
        printf("Error: distPerp too small: %g \n", distPerp);
        exit(1);
    }
    if (rowIndex > distPerpGridNumber - 2) {
#ifndef NDEBUG
        printf("Warning: distPerp %g very large, clamp to grid UB\n", distPerp);
#endif
        rowIndex = distPerpGridNumber - 2;
        rowFrac = 1;
    }

    // Every row of table starts from 0, for sbound=0
    // val should be positive
    if (val < 0) {
#ifndef NDEBUG
        printf("Error: val should be positive: %g \n", val);
#endif
        exit(1);
    }
    double colFrac = 0;

    // reverse lookup val on row rowIndex
    const int index1LB = getTableIndex(rowIndex, 0);
    const int index1UB = getTableIndex(rowIndex, sboundGridNumber - 1);
    const int index2LB = getTableIndex(rowIndex + 1, 0);
    const int index2UB = getTableIndex(rowIndex + 1, sboundGridNumber - 1);

    auto lower1 = std::lower_bound(table.begin() + index1LB,
                                   table.begin() + index1UB, val);
    // reverse lookup val on row rowIndex+1
    auto lower2 = std::lower_bound(table.begin() + index2LB,
                                   table.begin() + index2UB, val);

    if (lower1 == table.begin() + index1UB ||
        lower2 == table.begin() + index2UB) {
        sbound = sboundGrid.back();
#ifndef NDEBUG
        printf("Warning: val %g too large, setting sbound to max %g\n", val,
               sbound);
#endif
        return sbound;
    }

    // find cross point of distPerp and two points
    assert(lower1 - table.begin() >= index1LB);
    assert(lower2 - table.begin() >= index2LB);

    /**
     * value on table grids
     * rowIndex --------+--------+---
     *             colIndexm colIndexm+1
     *
     * rowIndex+1 -----------+--------+
     *             colIndexp colIndexp+1
     */
    const int colIndexm = lower1 - 1 - table.begin() - index1LB;
    const int colIndexp = lower2 - 1 - table.begin() - index2LB;
    double valmA = *(lower1 - 1);
    double valmB = *(lower1);
    double valpA = *(lower2 - 1);
    double valpB = *(lower2);
    double sboundGridSpacing = sboundGrid[1] - sboundGrid[0];
    double sboundm = sboundGrid[colIndexm] +
                     (val - valmA) / (valmB - valmA) * sboundGridSpacing;
    double sboundp = sboundGrid[colIndexp] +
                     (val - valpA) / (valpB - valpA) * sboundGridSpacing;
    sbound = sboundm * (1 - rowFrac) + sboundp * rowFrac;
    return sbound;
}
