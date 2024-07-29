/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_sf2,
    c_ch2,

    c_Ham,
    c_Mom,
    c_Ham_abs_terms,
    c_Mom_abs_terms,

    c_Ham_rescaled,
    c_Ham_var,
    c_Ham_abs_AAD,
    c_Mom_AAD,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "PhiSq",
    "ChiSq",
    
    "Ham",
    "Mom",
    "Ham_abs_terms",
    "Mom_abs_terms",
    
    "HamRescaled",
    "HamVar",
    "HamAbsAAD",
    "MomAAD"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
