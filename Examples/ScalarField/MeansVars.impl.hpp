/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

 #if !defined(MEANSVARS_HPP_)
 #error "Only include this file through MeansVars.hpp"
 #endif

 #ifndef MEANSVARS_IMPL_HPP
 #define MEANSVARS_IMPL_HPP

#include "Coordinates.hpp"
#include "VarsTools.hpp"
#include <cmath>

inline
 MeansVars::MeansVars(double dx, params_t a_params, int a_c_Ham, int a_c_Ham_abs, int a_c_Ham_var, int a_c_Ham_abs_AAD,
         double a_Ham_mean, double a_Ham_abs_mean, int a_c_Mom, int a_c_Mom_AAD, double a_Mom_mean) : 
    m_dx (dx), m_params (a_params), m_c_Ham(a_c_Ham), m_c_Ham_abs(a_c_Ham_abs), m_c_Ham_var(a_c_Ham_var), m_c_Ham_abs_AAD(a_c_Ham_abs_AAD),
        m_Ham_mean(a_Ham_mean), m_Ham_abs_mean(a_Ham_abs_mean), m_c_Mom(a_c_Mom), m_c_Mom_AAD(a_c_Mom_AAD), m_Mom_mean(a_Mom_mean) {}

 template <class data_t>
 void MeansVars::compute(Cell<data_t> current_cell) const
 {
     CH_TIME("MeansVars::compute");

     Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

     /*const auto vars = current_cell.template load_vars<Vars>();

     data_t phisq = vars.phi*vars.phi;
     data_t chisq = vars.chi*vars.chi;
     data_t kin = vars.Pi*vars.Pi;*/

     data_t Ham = current_cell.load_vars(m_c_Ham);
     data_t HamAbs = current_cell.load_vars(m_c_Ham_abs);
     data_t HamAbsAAD = abs(HamAbs - m_Ham_abs_mean);
     data_t HamVar = pow(Ham - m_Ham_mean, 2.);

     data_t Mom = current_cell.load_vars(m_c_Mom);
     data_t MomAAD = abs(Mom - m_Mom_mean);

    //store class (Vars) variables as diagnostic variables on the grid
     /*current_cell.store_vars(vars.phi, c_sf);
     current_cell.store_vars(vars.Pi, c_sfd);
     current_cell.store_vars(vars.chi, c_a);
     current_cell.store_vars(vars.K, c_H);
     current_cell.store_vars(phisq, c_sf2);
     current_cell.store_vars(chisq, c_ch2);
     current_cell.store_vars(kin, c_kin);*/

     current_cell.store_vars(HamVar, m_c_Ham_var);
     current_cell.store_vars(HamAbsAAD, m_c_Ham_abs_AAD);
     current_cell.store_vars(MomAAD, m_c_Mom_AAD);
 }

 /*template <class data_t>
 template <typename mapping_function_t>
 void MeansVars::Vars<data_t>::enum_mapping(mapping_function_t mapping_function)
 {
     using namespace VarsTools;
     define_enum_mapping(mapping_function, c_phi, phi);
     define_enum_mapping(mapping_function, c_Pi, Pi);
     define_enum_mapping(mapping_function, c_chi, chi);
     define_enum_mapping(mapping_function, c_K, K);
 }*/

 #endif /* MEANSVARS_IMPL_HPP_ */