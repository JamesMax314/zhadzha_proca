/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXCISIONPROCADIAGNOSTICS_HPP_
#define EXCISIONPROCADIAGNOSTICS_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Does excision for fixed BG BH solutions
//! Note that it is does not using simd so one must set disable_simd()
template <class matter_t, class background_t> class ExcisionProcaDiagnostics
{
    // Use matter_t class
    using Vars = typename matter_t::template Vars<double>;

  protected:
    const double m_dx;                              //!< The grid spacing
    const std::array<double, CH_SPACEDIM> m_center; //!< The BH center
    const FourthOrderDerivatives m_deriv;
    const background_t m_background;
    const double m_excision_width;

  public:
    ExcisionProcaDiagnostics(const double a_dx,
                             const std::array<double, CH_SPACEDIM> a_center,
                             background_t a_background,
                             double a_excision_width = 1.0)
        : m_dx(a_dx), m_deriv(m_dx), m_center(a_center),
          m_excision_width(a_excision_width), m_background(a_background)
    {
    }

    void compute(const Cell<double> current_cell) const
    {
        const Coordinates<double> coords(current_cell, m_dx, m_center);
        double horizon_distance = m_background.excise(current_cell);
        if (horizon_distance < m_excision_width)
        {
            current_cell.store_vars(0.0, c_rho);
            current_cell.store_vars(0.0, c_rhoJ);
        } // else do nothing
    }
};

#endif /* EXCISIONPROCADIAGNOSTICS_HPP_ */
