/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_sf,
    c_a,
    c_H,
    c_sf2,

    c_Ham,
    c_Mom,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Phi",
    "ScaleFactor",
    "HubbleFactor",
    "PhiSq",
    "Ham",
    "Mom"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
