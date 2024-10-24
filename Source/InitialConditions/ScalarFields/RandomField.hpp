/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef RANDOMFIELD_HPP_
#define RANDOMFIELD_HPP_

#include "Cell.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp"
#include "InitialScalarData.hpp"
#include "Coordinates.hpp"
#include "VarsTools.hpp"
#include "fftw3.h"
#include <random> // needed for random number generator
#include <fstream>

class RandomField
{
    // Use the variable definition in CCZ4
    template <class data_t>
    using Vars = ADMConformalVars::VarsWithGauge<data_t>;
    
    public:
        struct params_t
        {
            std::array<double, CH_SPACEDIM>
                center;   //!< Centre of the grid
            int N;        //!< Box resolution (first level)
            int Nf;       //!< Finest resolution to generate ICs for
                                //! (used for convergence testing)
            double L;     //!< Box length in physical units
            double A;     //!< Amplitude applied to tensor field
            std::string print_path;
        };

        RandomField(params_t a_params, InitialScalarData::params_t a_bkgd_params, std::string a_spec_type);

        template <class data_t> void compute(Cell<data_t> current_cell) const; 

        void calc_spectrum();
        void clear_data();

    private:
        double** hx;

    protected:
        const params_t m_params;
        const InitialScalarData::params_t m_bkgd_params;
        const std::string m_spec_type;

        bool debug;
        double kstar;
        double epsilon;
        double H0;
        double norm;
        
        int flip_index(int I, int N);
        int invert_index(int I, int N);
        int invert_index_with_sign(int I, int N);
        double find_rayleigh_factor(double km, std::string spec_type, int comp, double rand_mod);
        void apply_symmetry_rules(int i, int j, int k, double field[][2], int N);
        void calc_transferse_vectors(int x, int y, int z, int N, double MHat[3], double NHat[3], double a = 0.);
        void Test_norm(double vec[]);
        void Test_orth(double vec1[], double vec2[]);
};

#include "RandomField.impl.hpp"

#endif /* RANDOMFIELD_HPP_ */
