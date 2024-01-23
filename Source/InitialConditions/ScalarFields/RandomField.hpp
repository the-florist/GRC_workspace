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

        //template <class data_t>
        void calc_spectrum();

    private:
        fftw_complex** hk;
        fftw_complex** hx;

    protected:
        const InitialScalarData::params_t m_params;
        std::string m_spec_type;
        
        double find_rayleigh_factor(double km, double ks, double ep, std::string spec_type, double H0, double uniform_draw);
        int flip_index(int I, int N);
        int invert_index(int I, int N);
        void calc_transferse_vectors(int x, int y, int z, double MHat[3], double NHat[3], double a = 0.);
        void Test_norm(double vec[]);
        void Test_orth(double vec1[], double vec2[]);

};

#include "RandomField.impl.hpp"

#endif /* RANDOMFIELD_HPP_ */