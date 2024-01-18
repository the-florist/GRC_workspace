/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

 #ifndef RANDOMFIELD_HPP_
 #error "This file should only be included via RandomField.hpp"
 #endif

 #ifndef RANDOMFIELD_IMPL_HPP_
 #define RANDOMFIELD_IMPL_HPP_

 inline RandomField::RandomField(params_t a_params, double a_dx)
    : m_dx(a_dx), m_params(a_params)
{
}

template <class data_t>
void RandomField::compute(Cell<data_t> current_cell, double m_dx) const
{
    InitialScalarData::params_t params;

    params.kstar = 16.*(2.*M_PI/params.L);
    params.epsilon = 2./params.L;
    params.H0 = -3.0*sqrt((8.0 * M_PI/3.0/params.m_pl/params.m_pl)*(0.5*params.velocity*params.velocity + 0.5*pow(params.m * params.amplitude, 2.0)));
    params.norm = pow(params.N, 3.);

    Coordinates<data_t> coords(current_cell, m_dx, params.center);
}

template <class data_t>
void RandomField::calc_spectrum() const
{

}

template <class data_t>
data_t