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
    c_sfd,
    c_sf2,
    c_ch2,
    c_Ham_pbp,
    c_Ham_pbp_norm,

    c_a,
    c_H,

    c_pot,
    c_kin,

    c_Ham,
    c_Mom,
    c_Ham_abs_terms,
    c_Mom_abs_terms,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Phi",
    "Pi",
    "PhiSq",
    "ChiSq",
    "HamPBP",
    "HamPBPNorm",

    "ScaleFactor",
    "HubbleFactor",

    "PotED",
    "KinED",
    
    "Ham",
    "Mom",
    "Ham_abs_terms",
    "Mom_abs_terms"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
