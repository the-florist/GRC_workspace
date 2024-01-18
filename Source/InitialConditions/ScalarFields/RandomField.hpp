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
        struct params_t
        {
            const unsigned int seed1 = 7586572;
            const unsigned int seed2 = 3539263;
        };

        RandomField(params_t a_params, double a_dx);

        template <class data_t> void compute(Cell<data_t> current_cell, double m_dx) const; 

    protected:
        const params_t m_params;
        double m_dx
        template <class data_t>
        void calc_spectrum(argtype arg) const;

}

#include RandomField.impl.hpp

#endif /* RANDOMFIELD_HPP_ */