/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "InitialScalarData.hpp"
#include "RandomField.hpp"
#include "Potential.hpp"
#include "MeansVars.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp) 
    {
        pp.load("m_pl", initial_params.m_pl, 1.0);
        pp.load("G_Newton", G_Newton, pow(initial_params.m_pl, -2.0)); // natural units with m_pl mass units
        pp.load("scalar_mass", potential_params.scalar_mass, 0.01);

        // Initial scalar field data
        pp.load("scalar_amplitude", initial_params.amplitude, 10.0);
        pp.load("scalar_velocity", initial_params.velocity, -0.001);
        pp.load("scalar_mass", initial_params.m, 0.01);

        // Random fields initial data class
        random_field_params.center =
            center; // already read in SimulationParametersBase
        pp.load("N_full", random_field_params.N, 128);
        pp.load("N_fine", random_field_params.Nf, 128);
        pp.load("L_full", random_field_params.L, 4.);
        pp.load("tensor_amplitude", random_field_params.A, 1.e-6);
        pp.load("output_path", random_field_params.print_path);
    }

    void check_params()
    {
        warn_parameter("scalar_mass", potential_params.scalar_mass,
                       potential_params.scalar_mass <
                           0.2 / coarsest_dx / dt_multiplier,
                       "oscillations of scalar field do not appear to be "
                       "resolved on coarsest level");
        warn_parameter("N_fine", random_field_params.Nf,
                        random_field_params.Nf == random_field_params.N,
                        "initial conditions (if using RF class) will be coarse grained");
        warn_parameter("N_fine", random_field_params.Nf,
                        random_field_params.Nf > random_field_params.N,
                        "finest IC generation level appears to be less than the base resolution");
    }

    // Initial data for matter and potential and BH
    double G_Newton;
    InitialScalarData::params_t initial_params;
    RandomField::params_t random_field_params;
    Potential::params_t potential_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
