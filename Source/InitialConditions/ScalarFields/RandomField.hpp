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
#include "fftw3.h"
#include <random> // needed for random number generator
#include <fstream>

class RandomField
{
    public:
        RandomField(InitialScalarData::params_t a_params, std::string a_spec_type);

        template <class data_t> void compute(Cell<data_t> current_cell) const; 

	    void clear_data();
        void calc_spectrum();

    private:
        double** hx;
        int N;

    protected:
        const InitialScalarData::params_t m_params;
        std::string m_spec_type;

        double kstar;
        double epsilon;
        double H0;
        double norm;
        
        int flip_index(int I, int N);
        int invert_index(int I, int N);
        int invert_index_with_sign(int I, int N);
        double find_rayleigh_factor(double km, std::string spec_type, double uniform_draw, int comp);
        void apply_symmetry_rules(int i, int j, int k, double field[][2], int N);
        void calc_transferse_vectors(int x, int y, int z, double MHat[3], double NHat[3], double a = 0.);
        void Test_norm(double vec[]);
        void Test_orth(double vec1[], double vec2[]);
};

#include "RandomField.impl.hpp"

#endif /* RANDOMFIELD_HPP_ */
