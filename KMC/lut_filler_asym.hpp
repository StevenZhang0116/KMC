/**
 * @file lut_filler_fdep.hpp
 * @author Adam Lamson
 * @brief Class to fill lookup table with
 * @version 0.1
 * @date 2020-02-10
 *
 * @copyright Copyright (c) 2019
 *
 */
#ifndef LUT_FILLER_ASYM_HPP_
#define LUT_FILLER_ASYM_HPP_

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "integrals.hpp"
#include "lut_filler.hpp"
#include "macros.hpp"
/**
 * @brief 2D Lookup Table for the dimensionless functions
 *
 *
 * Both interface functions and internal caculations are dimensionless
 */
class LUTFillerAsym : public LUTFiller {
  private:
    double exp_fact1_;   ///< \beta * k_spring * D * D, dimensionless
    double exp_fact2_;   ///< \beta * k_spring * D * D, dimensionless
    double e_fact_;      ///< \lambda, dimensionless
    double fdep_length_; ///< x_c / D, dimensionless
    double rest_length_; ///< \ell_o / D, dimensionless

  public:
    static constexpr double fd_small_ = 1e-4;
    LUTFillerAsym(double dist_para_grid_num, double dist_perp_grid_num)
        : LUTFiller(dist_para_grid_num, dist_perp_grid_num) {}

    /*! \brief Initialize lookup table filler
     *
     * Stores precalculated values of the integral
     * \int_0^y exp{-exp_fact_ (sqrt[x^2 + l_min^2] - l_o - D)^2}dx
     * To do this, grid spacing are calculated based on the value of exp_fact_
     * and freeLength. The upper bound lUB is calculated so that the maximum
     * change at the end of the integral is less than 10^{-4}. Everything is
     * non-dimensionalized by rodD.
     *
     * \param exp_fact1 Exponential prefactor when spring is compressed with
     * the dimensions of L^{-2}.
     * \param exp_fact Exponential prefactor when spring is stretched with the
     * dimensions of L^{-2}.
     * \param freeLength Rest length of binding crosslink.
     * \param rodD Diameter of rod crosslink is binding to.
     *
     * \return void
     */
    void Init(double exp_fact1, double exp_fact2, double e_fact,
              double fdep_length, double rest_length, double rodD) {
        // constexpr double small = 1e-4;
        // setup length scale and dimensionless lengths
        length_scale_ = rodD;
        rest_length_ = rest_length / length_scale_;
        fdep_length_ = fdep_length / length_scale_;
        e_fact_ = e_fact;
        exp_fact1_ = exp_fact1 * length_scale_ * length_scale_;
        exp_fact2_ = exp_fact2 * length_scale_ * length_scale_;

        assert(exp_fact1_ >= 0);
        assert(exp_fact2_ > 0);
        assert(e_fact < 1);
        assert(rest_length_ > 0);

        LUTFiller::Init();
    }

    double getUpperBound() const {
        return rest_length_ +
               ((sqrt(SQR(fdep_length_) -
                      (2. * (1. - e_fact_) * log(fd_small_) / exp_fact2_)) +
                 fdep_length_) /
                (1. - e_fact_));
    }
    /*! \brief Get result from integrating Boltzmann factor along filament.
     *
     * \param dist_perp Minimum distance away from filament
     * \param dist_para_l Lower bound parallel distance along filament
     * \param dist_para_u Upper bound parallel distance along filament
     *
     * \return double Value of integral
     */
    double getIntegralResult(double dist_perp, double dist_para_l,
                             double dist_para_u) const {
        return asym_integral(dist_perp, dist_para_l, dist_para_u, exp_fact1_,
                             exp_fact2_, e_fact_, fdep_length_, rest_length_);
    }
    double getBindingVolume() const {
        return asym_bind_vol_integral(upper_bound_, exp_fact1_, exp_fact2_,
                                      e_fact_, fdep_length_, rest_length_);
    }
};

#endif
