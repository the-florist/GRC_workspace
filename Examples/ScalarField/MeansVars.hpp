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
        MeansVars(int a_c_phi, int a_c_chi) : m_c_phi(a_c_phi), m_c_chi(a_c_chi) {}

    
        template <class data_t> void compute(Cell<data_t> current_cell) const
        {
            CH_TIME("MeansVars::compute");

            data_t phi = current_cell.load_vars(m_c_phi);
            data_t chi = current_cell.load_vars(m_c_chi);
            data_t phisq = phi*phi;
            data_t chisq = chi*chi;

            //store class (Vars) variables as diagnostic variables on the grid
            current_cell.store_vars(phisq, c_sf2);
            current_cell.store_vars(chisq, c_ch2);
        }

    protected:
        const int m_c_phi;
        const int m_c_chi;
};

 //#include "MeansVars.impl.hpp"
 #endif /* MEANSVARS_HPP_ */