/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GAMMACALCULATOR_HPP_
#define GAMMACALCULATOR_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

class GammaCalculator
{
    // Only variables needed are metric
    template <class data_t> struct Vars
    {
        Tensor<2, data_t> h;

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_symmetric_enum_mapping(
                mapping_function, GRInterval<c_h11, c_h33>(), h);
        }
    };

  protected:
    const SixthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables

  public:
    GammaCalculator(double a_dx) : m_deriv(a_dx) {}

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables, and calc 1st
        // derivs
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);

        using namespace TensorAlgebra;
        const auto h_UU = compute_inverse_sym(vars.h);
        const auto chris = compute_christoffel(d1.h, h_UU);
        MayDay::Error("Chris calc ended.");

        /*simd<double> tol(1e-8);
        if(simd_compare_gt(chris.contracted[0], tol))
        {
            std::cout << "Christoffel comp 0 is large here: " << chris.contracted[0] << "\n";
            MayDay::Error("Chris failed.");
        }*/

        /*simd<double> tol(1e-8);
        simd<double> chrs(0.);
        for(int s=0; s<3; s++) { chrs += chris.contracted[s]; }
        if(simd_compare_gt(chrs, tol))
        {
            std::cout << "Chris has this magnitude in GammaCalc: " << chrs << "\n";
            MayDay::Error("Check chris error.");
        }*/

        // assign values of Gamma^k = h_UU^ij * \tilde{Gamma}^k_ij in the output
        // FArrayBox
        current_cell.store_vars(chris.contracted,
                                GRInterval<c_Gamma1, c_Gamma3>());
    }
};

#endif /* GAMMACALCULATOR_HPP_ */
