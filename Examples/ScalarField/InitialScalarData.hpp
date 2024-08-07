/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALSCALARDATA_HPP_
#define INITIALSCALARDATA_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4RHS.hpp"
#include "ScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "Potential.hpp"

#include "IntVect.H"
#include "MayDay.H"
#include <fstream>

//! Class which sets the initial scalar field matter config
class InitialScalarData
{
  public:
    //! A structure for the input params for scalar field properties and initial
    //! conditions
    struct params_t
    {
        double amplitude; //!< Amplitude of k=0 mode of initial SF
        double velocity;  //!< Amplitude of initial SF velocity
        double m;         //!< SF mass
        double m_pl = 1.; //!< Planck mass (units)
    };

    //! The constructor
    InitialScalarData(params_t a_params)
        : m_params(a_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // calculate and store the scalar field value
        const data_t phi = m_params.amplitude/m_params.m_pl;
        const data_t phidot = m_params.velocity/m_params.m_pl/m_params.m_pl;

        current_cell.store_vars(phi, c_phi);
        current_cell.store_vars(phidot, c_Pi);

        //calculate and store gauge variables
        data_t lapse = 1.0;
        data_t shift[3] = {0.0, 0.0, 0.0};

        current_cell.store_vars(lapse, c_lapse);
        current_cell.store_vars(shift[0], c_shift1);
        current_cell.store_vars(shift[1], c_shift2);
        current_cell.store_vars(shift[2], c_shift3);

        //calculate and store scalar metric variables
        data_t chi = 1.0; // a
        data_t K = -3.0*sqrt((8. * M_PI/3.)*(0.5*phidot*phidot + 0.5*pow(m_params.m/m_params.m_pl * phi, 2.0))); // K (Friedman's equations)

        current_cell.store_vars(chi, c_chi);
        current_cell.store_vars(K, c_K);
    }

  protected:
    const params_t m_params; //!< The matter initial condition params
};

#endif /* INITIALSCALARDATA_HPP_ */
