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

class RandomField
{
    public:
        RandomField(InitialScalarData::params_t a_params);

        template <class data_t> void compute(Cell<data_t> current_cell) const; 

        //template <class data_t>
        void calc_spectrum(std::string spec_type);

    private:
        fftw_complex** hk;
        fftw_complex** hx;

    protected:
        const InitialScalarData::params_t m_params;
        double find_rayleigh_factor(double km, double ks, double ep, std::string spec_type, double H0, double uniform_draw);
        int flip_index(int I, int N);
        int invert_index(int I, int N);

};

#include "RandomField.impl.hpp"

#endif /* RANDOMFIELD_HPP_ */