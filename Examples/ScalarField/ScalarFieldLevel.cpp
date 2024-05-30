/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "MatterCCZ4RHS.hpp"

// For constraints calculation
#include "NewMatterConstraints.hpp"

// For tag cells
#include "FixedGridsTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "GammaCalculator.hpp"
#include "InitialScalarData.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "SetValue.hpp"

//For printing out mean diagnostics
#include "AMRReductions.hpp"
#include "MeansVars.hpp"
#include <cmath>
#include <string>
#include <sstream>
#include <typeinfo>
#include <vector>
#include <fstream>

// Start time
#include <ctime>
#include <typeinfo>
#include "RandomField.hpp"

// Things to do at each advance step, after the RK4 is calculated 
void ScalarFieldLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void ScalarFieldLevel::initialData()
{
    time_t t;
    CH_TIME("ScalarFieldLevel::initialData");
    pout() << "ScalarFieldLevel::initialData " << m_level << endl;

    RandomField pfield(m_p.initial_params, "position");

    BoxLoops::loop(
        make_compute_pack(SetValue(0.),
                            pfield),
    m_state_new, m_state_new, INCLUDE_GHOST_CELLS, disable_simd());

    pfield.clear_data();
    pout() << "Calculating position ICs ended.\n";

    RandomField vfield(m_p.initial_params, "velocity");

    BoxLoops::loop(
        make_compute_pack(vfield),
    m_state_new, m_state_new, INCLUDE_GHOST_CELLS, disable_simd());

    vfield.clear_data();
    pout() << "Calculating velocity ICs ended.\n";

    BoxLoops::loop(
        make_compute_pack(InitialScalarData(m_p.initial_params)),
    m_state_new, m_state_new, INCLUDE_GHOST_CELLS,disable_simd());

    pout() << "IC set-up ended.\n";
    
    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);
}

#ifdef CH_USE_HDF5
// Things to do before outputting a checkpoint file
void ScalarFieldLevel::prePlotLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(
        MatterConstraints<ScalarFieldWithPotential>(
            scalar_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom, c_Mom), m_p.min_chi, 
	    c_Ham_abs_terms, Interval(c_Mom_abs_terms, c_Mom_abs_terms)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}
#endif

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    if (m_p.max_spatial_derivative_order == 4)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      FourthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      SixthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

// Things to do at ODE update, after soln + rhs
void ScalarFieldLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void ScalarFieldLevel::preTagCells()
{
    // we don't need any ghosts filled for the fixed grids tagging criterion
    // used here so don't fill any
}

void ScalarFieldLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state,
    const FArrayBox &current_state_diagnostics)
{
    BoxLoops::loop(
        FixedGridsTaggingCriterion(m_dx, m_level, 2.0 * m_p.L, m_p.center),
        current_state, tagging_criterion);
}

void ScalarFieldLevel::specificPostTimeStep()
{

    CH_TIME("ScalarFieldLevel::specificPostTimeStep");

    if (!FilesystemTools::directory_exists(m_p.data_path)) 
    {
        FilesystemTools::mkdir_recursive(m_p.data_path);
    }

    bool first_step = (m_time == 0.);
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(
        MatterConstraints<ScalarFieldWithPotential>(
            scalar_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom, c_Mom), m_p.min_chi, 
	    c_Ham_abs_terms, Interval(c_Mom_abs_terms, c_Mom_abs_terms)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

    double mass = m_p.potential_params.scalar_mass;//0.01;

    AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
    AMRReductions<VariableType::evolution> amr_reductions_evo(m_gr_amr);
    double vol = amr_reductions.get_domain_volume();

    double habsbar = amr_reductions.sum(c_Ham_abs_terms)/vol;

    BoxLoops::loop(MeansVars(m_dx, habsbar, m_p.grid_params, m_p.data_path), m_state_new, m_state_diagnostics, FILL_GHOST_CELLS);

    // Convergence testing only

    //double hampbpSum = amr_reductions.sum(c_Ham_pbp)/vol;
    double hamnormMax = amr_reductions.max(c_Ham_pbp_norm);
    //double mombar = amr_reductions.sum(c_Mom)/vol;

    // All other runs
    //Calculates means
    /*double phibar = amr_reductions.sum(c_sf)/vol;
    double pibar = amr_reductions.sum(c_sfd)/vol;

    double a = 1./sqrt(amr_reductions.sum(c_a)/vol);
    double H = -amr_reductions.sum(c_H)/vol/3.;

    double hambar = amr_reductions.sum(c_Ham)/vol;
    double hamabspbpSum = amr_reductions.sum(c_Ham_abs_pbp)/vol;
    double mombar = amr_reductions.sum(c_Mom)/vol;
    double habsbar = amr_reductions.sum(c_Ham_abs_terms)/vol;
    double mabsbar = amr_reductions.sum(c_Mom_abs_terms)/vol;

    //Calculates energy components and the slow-roll parameters
    double kinb = 0.5*pibar*pibar;
    double potb = 0.5*mass*mass*phibar*phibar;

    double epsilon = 3.*kinb/(kinb + potb);
    double delta = 3. + mass*mass*phibar/(pibar)/(H);

    //Calculates variances
    double phivar = amr_reductions.sum(c_sf2)/vol - phibar*phibar;
    double chivar = amr_reductions.sum(c_ch2)/vol - c_a*c_a;

    //Calculates gauge quantities
    double lapse = amr_reductions_evo.sum(c_lapse)/vol;*/

    //Prints all that out into the data/ directory
    SmallDataIO means_file(m_p.data_path+"means_file", m_dt, m_time, m_restart_time, SmallDataIO::APPEND, first_step, ".dat");
    means_file.remove_duplicate_time_data(); // removes any duplicate data from previous run (for checkpointing)

    if(first_step) 
    {
        means_file.write_header_line({"Ham Norm max"});
        /*means_file.write_header_line({"Scalar field mean","Scalar field variance","Pi mean","Scale factor","Conformal factor variance","Hubble factor",
            "Kinetic ED","Potential ED","First SRP","Second SRP","Avg Ham constr","Avg |Ham| constr (point by point)","Avg Mom constr",
            "Avg Ham abs term","Avg Mom abs term","Avg lapse"});*/
    }
    means_file.write_time_data_line({hamnormMax});
    //means_file.write_time_data_line({phibar, phivar, pibar, a, chivar, H, kinb, potb, epsilon, delta, hambar, hamabspbpSum, mombar, habsbar, mabsbar, lapse});
}
