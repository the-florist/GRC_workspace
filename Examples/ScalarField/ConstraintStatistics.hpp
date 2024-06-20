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
        ConstraintStatistics(int a_c_Ham, int a_c_Ham_abs, int a_c_Ham_var, int a_c_Ham_abs_AAD,
            double a_Ham_mean, double a_Ham_abs_mean, int a_c_Mom, int a_c_Mom_AAD, double a_Mom_mean) : 
        m_c_Ham(a_c_Ham), m_c_Ham_abs(a_c_Ham_abs), m_c_Ham_var(a_c_Ham_var), m_c_Ham_abs_AAD(a_c_Ham_abs_AAD),
            m_Ham_mean(a_Ham_mean), m_Ham_abs_mean(a_Ham_abs_mean), m_c_Mom(a_c_Mom), m_c_Mom_AAD(a_c_Mom_AAD), m_Mom_mean(a_Mom_mean) {}

    
        template <class data_t> void compute(Cell<data_t> current_cell) const
        {
            CH_TIME("ConstraintStatistics::compute");

            data_t Ham = current_cell.load_vars(m_c_Ham);
            data_t HamAbs = current_cell.load_vars(m_c_Ham_abs);
            data_t HamAbsAAD = abs(HamAbs - m_Ham_abs_mean);
            data_t HamVar = pow(Ham - m_Ham_mean, 2.);

            data_t Mom = current_cell.load_vars(m_c_Mom);
            data_t MomAAD = abs(Mom - m_Mom_mean);

            //store class (Vars) variables as diagnostic variables on the grid
            current_cell.store_vars(HamVar, m_c_Ham_var);
            current_cell.store_vars(HamAbsAAD, m_c_Ham_abs_AAD);
            current_cell.store_vars(MomAAD, m_c_Mom_AAD);
        }

    protected:
        const int m_c_Ham;
        const int m_c_Ham_abs;
        const int m_c_Ham_var;
        const int m_c_Ham_abs_AAD;
        double m_Ham_mean;
        double m_Ham_abs_mean;
        const int m_c_Mom;
        const int m_c_Mom_AAD;
        double m_Mom_mean;
};

 #endif /* CONSTRAINTSTATISTICS_HPP_ */