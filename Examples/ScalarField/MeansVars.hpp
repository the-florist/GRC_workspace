/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

 //This calculates the mean and variance of cosmologically-relevant quantities
 //Specifically developed for inflation models.

#ifndef MEANSVARS_HPP_
#define MEANSVARS_HPP_

#include "Cell.hpp"
#include <array>
#include "Coordinates.hpp"
#include "VarsTools.hpp"
#include <cmath>

class MeansVars 
{
    public:
        template <class data_t> void compute(Cell<data_t> current_cell) const
        {
            CH_TIME("MeansVars::compute");

            data_t phi = current_cell.load_vars(c_phi);
            data_t chi = current_cell.load_vars(c_chi);
            data_t phisq = phi*phi;
            data_t chisq = chi*chi;

            //store class (Vars) variables as diagnostic variables on the grid
            current_cell.store_vars(phisq, c_sf2);
            current_cell.store_vars(chisq, c_ch2);
        }
};

 #endif /* MEANSVARS_HPP_ */