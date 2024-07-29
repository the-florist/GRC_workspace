/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

 //This calculates the mean and variance of cosmologically-relevant quantities
 //Specifically developed for inflation models.

#ifndef CONSTRAINTSTATISTICS_HPP_
#define CONSTRAINTSTATISTICS_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "VarsTools.hpp"
#include <cmath>
#include <array>

class ConstraintStatistics
{
    public:
        ConstraintStatistics(double a_Ham_mean, double a_Ham_abs_mean,double a_Mom_mean) : 
            m_Ham_mean(a_Ham_mean), m_Ham_abs_mean(a_Ham_abs_mean), m_Mom_mean(a_Mom_mean) {}

    
        template <class data_t> void compute(Cell<data_t> current_cell) const
        {
            CH_TIME("ConstraintStatistics::compute");

            data_t Ham = current_cell.load_vars(c_Ham);
            data_t HamAbs = current_cell.load_vars(c_Ham_abs);
            data_t HamAbsAAD = abs(HamAbs - m_Ham_abs_mean);
            data_t HamVar = pow(Ham - m_Ham_mean, 2.);
	    data_t HamNorm = Ham/m_Ham_abs_mean;
	    data_t HamNormSq = HamNorm*HamNorm;

            data_t Mom = current_cell.load_vars(c_Mom);
            data_t MomAAD = abs(Mom - m_Mom_mean);

            //store class (Vars) variables as diagnostic variables on the grid
            current_cell.store_vars(HamNorm, c_Ham_norm);
	    current_cell.store_vars(HamNormSq, c_Ham_norm_sq);
	    current_cell.store_vars(HamVar, c_Ham_var);
            current_cell.store_vars(HamAbsAAD, c_Ham_abs_AAD);
            current_cell.store_vars(MomAAD, c_Mom_AAD);
        }

    protected:
        double m_Ham_mean;
        double m_Ham_abs_mean;
        double m_Mom_mean;
};

 #endif /* CONSTRAINTSTATISTICS_HPP_ */
