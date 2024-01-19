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

        //ouble calc_h_real();
        void calc_spectrum();

    private:
        fftw_complex** hk;
        fftw_complex** hx;

    protected:
        const InitialScalarData::params_t m_params;

};

#include "RandomField.impl.hpp"

#endif /* RANDOMFIELD_HPP_ */