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
#include "ConstraintStatistics.hpp"
#include <cmath>
#include <string>
#include <sstream>
#include <typeinfo>
#include <vector>
#include <fstream>

// Start time
#include <ctime>
#include <typeinfo>

// Initial conditions
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

    RandomField pfield(m_p.random_field_params, m_p.initial_params, "position");

    BoxLoops::loop(
        make_compute_pack(SetValue(0.), pfield),
    m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());

    pfield.clear_data();
    pout() << "Calculating position ICs ended.\n";

    RandomField vfield(m_p.random_field_params, m_p.initial_params, "velocity");

    BoxLoops::loop(
        make_compute_pack(vfield),
    m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());

    vfield.clear_data();
    pout() << "Calculating velocity ICs ended.\n";

    BoxLoops::loop(
        make_compute_pack(InitialScalarData(m_p.initial_params)),
    m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());

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

    double mass = m_p.potential_params.scalar_mass;

    AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
    AMRReductions<VariableType::evolution> amr_reductions_evo(m_gr_amr);
    double vol = amr_reductions.get_domain_volume();

    double hamBar = amr_reductions.sum(c_Ham)/vol;
    double hamAbsBar = amr_reductions.sum(c_Ham_abs_terms)/vol;
    double momBar = amr_reductions.sum(c_Mom)/vol;

    if (first_step) { pout() << "Domain volume: " << vol << "\n"; }

    BoxLoops::loop(
        ConstraintStatistics(hamBar, hamAbsBar, momBar),
        m_state_diagnostics, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

    // Convergence testing only
    double hamNormBar = amr_reductions.sum(c_Ham_norm)/vol;
    double hamNormVar = amr_reductions.sum(c_Ham_norm_var)/vol;
    double hamVar = amr_reductions.sum(c_Ham_var)/vol;
    double hamAbsAAD = amr_reductions.sum(c_Ham_abs_AAD)/vol;
    double momAAD = amr_reductions.sum(c_Mom_AAD)/vol;

    BoxLoops::loop(MeansVars(),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

    // All other runs
    //Calculates scalar field and FLRW metric means
    double phibar = amr_reductions_evo.sum(c_phi)/vol;
    double pibar = amr_reductions_evo.sum(c_Pi)/vol;

    double chibar = amr_reductions_evo.sum(c_chi)/vol;
    double a = 1./sqrt(chibar);
    double H = -amr_reductions_evo.sum(c_K)/vol/3.;

    //Calculates energy components and the slow-roll parameters
    double kinb = 0.5*pibar*pibar;
    double potb = 0.5*mass*mass*phibar*phibar;

    double epsilon = 3.*kinb/(kinb + potb);
    double delta = 3. + mass*mass*phibar/(pibar)/(H);

    //Calculates variances
    double phivar = amr_reductions.sum(c_sf2)/vol - phibar*phibar;
    double chivar = amr_reductions.sum(c_ch2)/vol - chibar*chibar;

    //Calculates gauge quantities
    double lapse = amr_reductions_evo.sum(c_lapse)/vol;

    //Prints all that out into the data/ directory
    SmallDataIO means_file(m_p.data_path+"means_file", m_dt, m_time, m_restart_time, SmallDataIO::APPEND, first_step, ".dat");
    means_file.remove_duplicate_time_data(); // removes any duplicate data from previous run (for checkpointing)

    SmallDataIO constrs_file(m_p.data_path+"constraints_file", m_dt, m_time, m_restart_time, SmallDataIO::APPEND, first_step, ".dat");
    constrs_file.remove_duplicate_time_data(); // removes any duplicate data from previous run (for checkpointing)

    if(first_step) 
    {
        constrs_file.write_header_line({"HamMean","HamSTD","HamAbsMean","HamNormMean","HamNormSTD","MomBar","MomAAD"});
        means_file.write_header_line({"PhiMean","PhiVar","PiMean","ScaleFact","ChiVar","HubbleFact",
            "KinED","PotED","SRP1","SRP2","LapseMean"});
    }
    constrs_file.write_time_data_line({hamBar, sqrt(hamVar), hamAbsBar, hamNormBar, sqrt(hamNormVar), momBar, momAAD});
    means_file.write_time_data_line({phibar, phivar, pibar, a, chivar, H, kinb, potb, epsilon, delta, lapse});
}
