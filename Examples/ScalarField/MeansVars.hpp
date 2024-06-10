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

class MeansVars 
{
    public:
        struct params_t
        {
            std::array<double, CH_SPACEDIM>
                center; 
        };

        /*template <class data_t>
        struct Vars 
        {
            data_t phi;
            data_t Pi;
            data_t chi;
            data_t K;
            
            template <typename mapping_function_t>
            void enum_mapping(mapping_function_t mapping_function);
        };*/

        MeansVars(double dx, params_t a_params, int a_c_Ham, int a_c_Ham_abs, int a_c_Ham_var, int a_c_Ham_abs_AAD,
         double a_Ham_mean, double a_Ham_abs_mean, int a_c_Mom, int a_c_Mom_AAD, double a_Mom_mean);

        template <class data_t>
        void compute(Cell<data_t> current_cell) const;

    protected:
        double m_dx;
        const params_t m_params;
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

 #include "MeansVars.impl.hpp"
 #endif /* MEANSVARS_HPP_ */