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
        double E;      //!< Energy scale [Mp]
    };

    //! The constructor
    InitialScalarData(params_t a_params)
        : m_params(a_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // calculate and store the scalar field and mass values
        const data_t Mp = 1./m_params.E;

        current_cell.store_vars(m_params.amplitude, c_phi);
        current_cell.store_vars(m_params.velocity, c_Pi);

        //calculate and store gauge variables
        data_t lapse = 1.0;
        data_t shift[3] = {0.0, 0.0, 0.0};

        current_cell.store_vars(lapse, c_lapse);
        current_cell.store_vars(shift[0], c_shift1);
        current_cell.store_vars(shift[1], c_shift2);
        current_cell.store_vars(shift[2], c_shift3);

        //calculate and store scalar metric variables
        data_t chi = 1.0; // a
        data_t K = -3.0*sqrt((8. * M_PI/3./pow(Mp, 2.))*(0.5*pow(m_params.velocity, 2.) + 0.5*pow(m_params.m * m_params.amplitude, 2.0))); // K (Friedman's equations)

        current_cell.store_vars(chi, c_chi);
        current_cell.store_vars(K, c_K);

        const data_t h11 = current_cell.load_vars(c_h11);
        const data_t A11 = current_cell.load_vars(c_A11);

        if(h11 == 0) { MayDay::Error("h11 has been overwritten."); }
        if(A11 == 0) { MayDay::Error("A11 has not been initialised."); }
    }

  protected:
    const params_t m_params; //!< The matter initial condition params
};

#endif /* INITIALSCALARDATA_HPP_ */
