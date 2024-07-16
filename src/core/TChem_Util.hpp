/* =====================================================================================
TChem-atm version 1.0
Copyright (2024) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.

This file is part of TChem-atm. TChem-atm is open source software: you can redistribute it
and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Contact Oscar Diaz-Ibarra at <odiazib@sandia.gov>, or
           Mike Schmidt at <mjschm@sandia.gov>, or
           Cosmin Safta at <csafta@sandia.gov>

Sandia National Laboratories, New Mexico/Livermore, NM/CA, USA
===================================================================================== */
#ifndef __TCHEM_UTIL_HPP__
#define __TCHEM_UTIL_HPP__

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <regex>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <cassert>
#include <cmath>
#include <ctime>
#include <limits>

#include <complex>

/// kokkos
#include "Kokkos_Complex.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_DualView.hpp"
#include "Kokkos_Timer.hpp"

/// blas
#include "Tines.hpp"

/// tchem configure header
#include "TChem_atm_Config.hpp"

#if defined(TCHEM_ATM_ENABLE_TPL_YAML_CPP)
#include "yaml-cpp/yaml.h"
#endif

namespace TChem {

/// exec spaces
using exec_space = Kokkos::DefaultExecutionSpace;
using host_exec_space = Kokkos::DefaultHostExecutionSpace;

/// this is for user interface
using interf_device_type = typename Tines::UseThisDevice<exec_space>::type;
using interf_host_device_type = typename Tines::UseThisDevice<host_exec_space>::type;

template <typename ExecSpace> struct UseThisTeamPolicy {
  using type = Kokkos::TeamPolicy<ExecSpace, Kokkos::Schedule<Kokkos::Static>>;
};

/// type defs
#if defined(TCHEM_ATM_ENABLE_REAL_TYPE_SINGLE_PRECISION)
using real_type = float;
#elif defined(TCHEM_ATM_ENABLE_REAL_TYPE_DOUBLE_PRECISION)
using real_type = double;
#endif
using ordinal_type = int;
using size_type = size_t;

/// control parameters
#if defined(TCHEM_ATM_ENABLE_VERBOSE)
static constexpr bool verboseEnabled = true;
#else
static constexpr bool verboseEnabled = false;
#endif

#if defined(TCHEM_ATM_ENABLE_DEBUG)
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
#define TCHEM_DEBUG_CHECK_ERROR(err, msg)                                                                              \
  if (err) {                                                                                                           \
    std::stringstream ss;                                                                                              \
    ss << msg << "\n" << __FILE__ << ", " << __LINE__ << "\n";                                                         \
    throw std::logic_error(ss.str().c_str());                                                                          \
  }
#else
#define TCHEM_DEBUG_CHECK_ERROR(err, msg)                                                                              \
  if (err)                                                                                                             \
    printf("%s\n%s, %d\n", msg, __FILE__, __LINE__);
#endif
#else
#define TCHEM_DEBUG_CHECK_ERROR(err, msg)
#endif

#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
#define TCHEM_CHECK_ERROR(err, msg)                                                                                    \
  if (err) {                                                                                                           \
    std::stringstream ss;                                                                                              \
    ss << msg << "\n" << __FILE__ << ", " << __LINE__ << "\n";                                                         \
    throw std::logic_error(ss.str().c_str());                                                                          \
  }
#else
#define TCHEM_CHECK_ERROR(err, msg)                                                                                    \
  if (err)                                                                                                             \
    printf("%s\n%s, %d\n", msg, __FILE__, __LINE__);
#endif

/// parameters

/**
 * \def LENGTHOFELEMNAME
 * Maximum number of characters for element names
 */
static constexpr ordinal_type LENGTHOFELEMNAME = 3;

/**
 * \def LENGTHOFSPECNAME
 *  Maximum number of characters for species names
 */
static constexpr ordinal_type LENGTHOFSPECNAME = 18;

/**
 * \def NUMBEROFELEMINSPEC
 * Maximum number of (different) elements that compose a species
 */
static constexpr ordinal_type NUMBEROFELEMINSPEC = 5;

/**
 * \def NTHRDBMAX
 * Maximum number of third body efficiencies
 */
constexpr ordinal_type NTHRDBMAX = 20;

/**
 * \def NPLOGMAX
 * Maximum number of interpolation ranges for PLOG
 */
static constexpr ordinal_type NPLOGMAX = 30;

/**
 * \def NTH9RNGMAX
 * Maximum number of temperature ranges for 9-coefficients NASA polynomials
 */
static constexpr ordinal_type NTH9RNGMAX = 5;

/**
 * \def NSPECREACMAX
 * Maximum number of reactant or product species in a reaction
 */
static constexpr ordinal_type NSPECREACMAX = 6;

/**
 * \def NAVOG
 * Avogadro's number
 */
static constexpr real_type NAVOG = 6.02214179E23;

/**
 * \def RUNIV
 * Universal gas constant \f$J/(mol\cdot K)\f$
 */
// constexpr real_type RUNIV = 8.3144621;
static constexpr real_type RUNIV = 8.314472;

// Avogadro's number (mole^{-1}).
static constexpr real_type AVAGADRO = 6.02214179E23;

// from camp converson #/cc -> ppm conversion prefactor
// CONV_PPM = AVAGADRO / RUNIV * 10^(-12)
static constexpr real_type CONV_PPM = 7.24296358205307E+10;

/**
 * \def CALJO
 * Conversion from calories to Joule
 */
static constexpr real_type CALJO = 4.184;

/**
 * \def KBOLT
 * Boltzmann's constant \f$(k_B)\f$ \f$[JK^{-1}]\f$
 */
static constexpr real_type KBOLT = 1.3806504E-23;

/**
 * \def EVOLT
 * electron volt (eV) unit \f$[J]\f$
 */
static constexpr real_type EVOLT = 1.60217653E-19;

/// callable from device
/// old cuda compilers do not support constexpr values well
/// need to make a constexpr function
/**
 * \def TCSMALL
 * Threshold for various numerical estimates
 */
static KOKKOS_FORCEINLINE_FUNCTION constexpr real_type TCSMALL() { return 1.e-12; }

static KOKKOS_FORCEINLINE_FUNCTION constexpr real_type PI() { return 3.141592653589793; }
// 2*pi
static KOKKOS_FORCEINLINE_FUNCTION constexpr real_type DPI() { return 6.283185307179586; }

/**
 * \def REACBALANCE
 * Treshold for checking reaction balance with
 *        real stoichiometric coefficients
 */
static KOKKOS_FORCEINLINE_FUNCTION constexpr real_type REACBALANCE() { return 1.e-4; }

/**
 * \def ATMPA
 * Standard atmospheric pressure \f$[Pa]\f$
 */
static KOKKOS_FORCEINLINE_FUNCTION constexpr real_type ATMPA() {
  return 101325.0;
}

struct LinozInputParameters {
  // Linoz climatological data
  real_type linoz_o3_clim;    //    fields(o3_clim_ndx)      %data(:,:,lchnk )
  real_type linoz_n2o_clim;   // fields(n2o_clim_ndx)     %data(:,:,lchnk )
  real_type linoz_noy_clim;   // fields(noy_clim_ndx)     %data(:,:,lchnk )
  real_type linoz_ch4_clim;   // fields(ch4_clim_ndx)     %data(:,:,lchnk )
  real_type linoz_h2o_clim;   // fields(h2o_clim_ndx)     %data(:,:,lchnk )
  real_type linoz_t_clim;     // fields(t_clim_ndx)       %data(:,:,lchnk )
  real_type linoz_o3col_clim; // fields(o3col_clim_ndx)   %data(:,:,lchnk )

  // OZONE P-L terms !unit vmr/sec
  real_type linoz_PmL_clim_no3; // fields(no3_PmL_clim_ndx)     %data(:,:,lchnk
                                // )    !unit vmr/sec
  real_type linoz_dPmL_dO3X; // fields(no3_dPmL_dO3_ndx)     %data(:,:,lchnk )
  real_type
      linoz_dPmL_dn2o_no3; // fields(no3_dPmL_dN2O_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dnoy_no3; // fields(no3_dPmL_dNOY_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dch4_no3; // fields(no3_dPmL_dCH4_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dh2o_no3;     // fields(no3_dPmL_dH2O_ndx)    %data(:,:,lchnk )
  real_type linoz_dPmL_dT_no3; // fields(no3_dPmL_dT_ndx)      %data(:,:,lchnk )
  real_type
      linoz_dPmL_dO3col_no3; // fields(no3_dPmL_dO3col_ndx)  %data(:,:,lchnk )

  // pointer to pn2o
  real_type linoz_PmL_clim_pn2o; // fields(pn2o_PmL_clim_ndx) %data(:,:,lchnk )
  real_type linoz_dPmL_dO3_pn2o; // fields(pn2o_dPmL_dO3_ndx) %data(:,:,lchnk )
  real_type
      linoz_dPmL_dn2o_pn2o; // fields(pn2o_dPmL_dN2O_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dnoy_pn2o; // fields(pn2o_dPmL_dNOY_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dch4_pn2o; // fields(pn2o_dPmL_dCH4_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dh2o_pn2o; // fields(pn2o_dPmL_dH2O_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dT_pn2o; // fields(pn2o_dPmL_dT_ndx)      %data(:,:,lchnk )
  real_type
      linoz_dPmL_dO3col_pn2o; // fields(pn2o_dPmL_dO3col_ndx)  %data(:,:,lchnk )

  // pointer to ln2o
  real_type linoz_PmL_clim_ln2o; // fields(ln2o_PmL_clim_ndx) %data(:,:,lchnk )
                                 // !1/sec not vrm/sec already scaled by n2o vmr
  real_type linoz_dPmL_dO3_ln2o; // fields(ln2o_dPmL_dO3_ndx) %data(:,:,lchnk )
  real_type
      linoz_dPmL_dn2o_ln2o; // fields(ln2o_dPmL_dN2O_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dnoy_ln2o; // fields(ln2o_dPmL_dNOY_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dch4_ln2o; // fields(ln2o_dPmL_dCH4_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dh2o_ln2o; // fields(ln2o_dPmL_dH2O_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dT_ln2o; // fields(ln2o_dPmL_dT_ndx)      %data(:,:,lchnk )
  real_type
      linoz_dPmL_dO3col_ln2o; // fields(ln2o_dPmL_dO3col_ndx)  %data(:,:,lchnk )

  // pointer to production of NOy
  real_type linoz_PmL_clim_pnoy; // fields(pnoy_PmL_clim_ndx) %data(:,:,lchnk )
  real_type linoz_dPmL_dO3_pnoy; // fields(pnoy_dPmL_dO3_ndx) %data(:,:,lchnk )
  real_type
      linoz_dPmL_dn2o_pnoy; // fields(pnoy_dPmL_dN2O_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dnoy_pnoy; // fields(pnoy_dPmL_dNOY_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dch4_pnoy; // fields(pnoy_dPmL_dCH4_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dh2o_pnoy; // fields(pnoy_dPmL_dH2O_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dT_pnoy; // fields(pnoy_dPmL_dT_ndx)      %data(:,:,lchnk )
  real_type
      linoz_dPmL_dO3col_pnoy; // fields(pnoy_dPmL_dO3col_ndx)  %data(:,:,lchnk )

  // pointer to loss of NOY
  real_type linoz_PmL_clim_lnoy; // fields(lnoy_PmL_clim_ndx) %data(:,:,lchnk )
  real_type linoz_dPmL_dO3_lnoy; // fields(lnoy_dPmL_dO3_ndx) %data(:,:,lchnk )
  real_type
      linoz_dPmL_dn2o_lnoy; // fields(lnoy_dPmL_dN2O_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dnoy_lnoy; // fields(lnoy_dPmL_dNOY_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dch4_lnoy; // fields(lnoy_dPmL_dCH4_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dh2o_lnoy; // fields(lnoy_dPmL_dH2O_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dT_lnoy; // fields(lnoy_dPmL_dT_ndx)      %data(:,:,lchnk )
  real_type
      linoz_dPmL_dO3col_lnoy; // fields(lnoy_dPmL_dO3col_ndx)  %data(:,:,lchnk )

  // pointer to loss of CH4
  real_type linoz_PmL_clim_nch4; // fields(nch4_PmL_clim_ndx) %data(:,:,lchnk )
  real_type linoz_dPmL_dO3_nch4; // fields(nch4_dPmL_dO3_ndx) %data(:,:,lchnk )
  real_type
      linoz_dPmL_dn2o_nch4; // fields(nch4_dPmL_dN2O_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dnoy_nch4; // fields(nch4_dPmL_dNOY_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dch4_nch4; // fields(nch4_dPmL_dCH4_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dh2o_nch4; // fields(nch4_dPmL_dH2O_ndx)    %data(:,:,lchnk )
  real_type
      linoz_dPmL_dT_nch4; // fields(nch4_dPmL_dT_ndx)      %data(:,:,lchnk )
  real_type
      linoz_dPmL_dO3col_nch4; // fields(nch4_dPmL_dO3col_ndx)  %data(:,:,lchnk )

  real_type linoz_cariolle_psc; // fields(cariolle_pscs_ndx)%data(:,:,lchnk )
};

using linoz_input_parameters_type = LinozInputParameters;

using linoz_input_parameters_1d_dual_view = Tines::value_type_1d_dual_view<linoz_input_parameters_type, exec_space>;
using linoz_input_parameters_0d_view = Tines::value_type_0d_view<linoz_input_parameters_type, interf_device_type>;
using linoz_input_parameters_1d_view = Tines::value_type_1d_view<linoz_input_parameters_type, interf_device_type>;
using linoz_input_parameters_1d_view_host = Tines::value_type_1d_view<linoz_input_parameters_type, interf_host_device_type>;
using linoz_input_parameters_0d_view_host = Tines::value_type_0d_view<linoz_input_parameters_type, interf_host_device_type>;

struct LinozVolumenMixinRatioIdx {
// input index
ordinal_type o3_ndx;
ordinal_type n2olnz_ndx;
ordinal_type noylnz_ndx;
ordinal_type ch4lnz_ndx;
ordinal_type h2olnz_ndx;
// output index
ordinal_type o3lnz_ndx;
ordinal_type ch4_ndx;
ordinal_type n2o_ndx;
ordinal_type no_ndx;
ordinal_type no2_ndx;
ordinal_type hno3_ndx;
};

using linoz_vmr_idx_type = LinozVolumenMixinRatioIdx;

struct ArrheniusReactionType {
  ordinal_type _reaction_index;
  real_type _E;
  real_type _A;
  real_type _B;
  real_type _C;
  real_type _D;
};

using arrhenius_reaction_type_1d_dual_view = Tines::value_type_1d_dual_view<ArrheniusReactionType, exec_space>;

struct CMAQ_H2O2ReactionType {
  ordinal_type _reaction_index;
  real_type _A1;
  real_type _B1;
  real_type _C1;
  real_type _A2;
  real_type _B2;
  real_type _C2;
};

using cmaq_h2o2_type_1d_dual_view = Tines::value_type_1d_dual_view<CMAQ_H2O2ReactionType, exec_space>;


struct JPL_ReactionType {
  ordinal_type _reaction_index;
  real_type _k0_A;
  real_type _k0_B;
  real_type _kinf_A;
  real_type _kinf_B;
  real_type _Fc;
  real_type _k0_C;
  real_type _kinf_C;
};

using jpl_reaction_type_1d_dual_view = Tines::value_type_1d_dual_view<JPL_ReactionType, exec_space>;

struct R_JPL_ArrheniusReactionType {
  ordinal_type _reaction_index;
  // parameters for JPL reaction
  real_type _k0_A;
  real_type _k0_B;
  real_type _kinf_A;
  real_type _kinf_B;
  real_type _Fc;
  real_type _k0_C;
  real_type _kinf_C;
  // parameters for JPL Arrhenius
  real_type _A;
  real_type _B;
  real_type _C;

};

using r_jpl_arrhenius_reaction_type_1d_dual_view = Tines::value_type_1d_dual_view<R_JPL_ArrheniusReactionType, exec_space>;

struct EMISSION_SourceType {
  ordinal_type _species_index;
  real_type _emissition_rate;
};

using emission_source_type_1d_dual_view = Tines::value_type_1d_dual_view<EMISSION_SourceType, exec_space>;

struct adjust_reactionType {
ordinal_type _reaction_index;
ordinal_type _species_index;
};

using adjust_reaction_type_1d_dual_view = Tines::value_type_1d_dual_view<adjust_reactionType, exec_space>;

// modifier for uci1, uci2, uci3 reaction types.
// factor = rate_photo_O1D/fc
// fc = k1[X1] + k2[X2] + k3[X3]
// X1=[N2], X2=[O2], X3[H2O]
// ki=Ai*exp(Ci/Temp)
struct prod_O1DType{
//NOTE: it only works for 3 reactions
ordinal_type _reaction_indices[3];
//
ordinal_type _photolysis_reaction_index;
ordinal_type _species_index_1;
real_type _A1;
real_type _C1;

ordinal_type _species_index_2;
real_type _A2;
real_type _C2;

ordinal_type _species_index_3;
real_type _A3;
real_type _C3;
};

using prod_O1D_type_1d_dual_view = Tines::value_type_1d_dual_view<prod_O1DType, exec_space>;


/// kinetic model data
struct KineticModelData;
using kmd_type = KineticModelData;
using kmd_type_1d_view_host = Tines::value_type_1d_view<kmd_type, interf_host_device_type>;

/// time marching data structure
struct TimeAdvance {
  real_type _tbeg, _tend;
  real_type _dt, _dtmin, _dtmax;
  ordinal_type _max_num_newton_iterations;
  ordinal_type _num_time_iterations_per_interval;
  ordinal_type _num_outer_time_iterations_per_interval;
  ordinal_type _jacobian_interval;
};

/// time tolerence; real_type_2d_view numberOfTimeODEs x 2 (atol,rtol)
/// newton tolerence; real_type_1d_view (atol, rtol)
using time_advance_type = TimeAdvance;
using time_advance_type_0d_dual_view = Tines::value_type_0d_dual_view<time_advance_type, interf_device_type>;
using time_advance_type_1d_dual_view = Tines::value_type_1d_dual_view<time_advance_type, interf_device_type>;

using time_advance_type_0d_view = Tines::value_type_0d_view<time_advance_type, interf_device_type>;
using time_advance_type_1d_view = Tines::value_type_1d_view<time_advance_type, interf_device_type>;

using time_advance_type_0d_view_host = Tines::value_type_0d_view<time_advance_type, interf_host_device_type>;
using time_advance_type_1d_view_host = Tines::value_type_1d_view<time_advance_type, interf_host_device_type>;

/// view
using real_type_0d_dual_view = Tines::value_type_0d_dual_view<real_type, interf_device_type>;
using real_type_1d_dual_view = Tines::value_type_1d_dual_view<real_type, interf_device_type>;
using real_type_2d_dual_view = Tines::value_type_2d_dual_view<real_type, interf_device_type>;
using real_type_3d_dual_view = Tines::value_type_3d_dual_view<real_type, interf_device_type>;

using ordinal_type_0d_dual_view = Tines::value_type_0d_dual_view<ordinal_type, interf_device_type>;
using ordinal_type_1d_dual_view = Tines::value_type_1d_dual_view<ordinal_type, interf_device_type>;
using ordinal_type_2d_dual_view = Tines::value_type_2d_dual_view<ordinal_type, interf_device_type>;
using ordinal_type_3d_dual_view = Tines::value_type_3d_dual_view<ordinal_type, interf_device_type>;

template <int S> using string_type_1d_dual_view = Tines::value_type_1d_dual_view<char[S], interf_device_type>;

using real_type_0d_view = Tines::value_type_0d_view<real_type, interf_device_type>;
using real_type_1d_view = Tines::value_type_1d_view<real_type, interf_device_type>;
using real_type_2d_view = Tines::value_type_2d_view<real_type, interf_device_type>;
using real_type_3d_view = Tines::value_type_3d_view<real_type, interf_device_type>;

using ordinal_type_0d_view = Tines::value_type_0d_view<ordinal_type, interf_device_type>;
using ordinal_type_1d_view = Tines::value_type_1d_view<ordinal_type, interf_device_type>;
using ordinal_type_2d_view = Tines::value_type_2d_view<ordinal_type, interf_device_type>;
using ordinal_type_3d_view = Tines::value_type_3d_view<ordinal_type, interf_device_type>;

template <int S> using string_type_1d_view = Tines::value_type_1d_dual_view<char[S], interf_device_type>;

using real_type_0d_view_host = Tines::value_type_0d_view<real_type, interf_host_device_type>;
using real_type_1d_view_host = Tines::value_type_1d_view<real_type, interf_host_device_type>;
using real_type_2d_view_host = Tines::value_type_2d_view<real_type, interf_host_device_type>;
using real_type_3d_view_host = Tines::value_type_3d_view<real_type, interf_host_device_type>;

using ordinal_type_0d_view_host = Tines::value_type_0d_view<ordinal_type, interf_host_device_type>;
using ordinal_type_1d_view_host = Tines::value_type_1d_view<ordinal_type, interf_host_device_type>;
using ordinal_type_2d_view_host = Tines::value_type_2d_view<ordinal_type, interf_host_device_type>;
using ordinal_type_3d_view_host = Tines::value_type_3d_view<ordinal_type, interf_host_device_type>;

template <int S> using string_type_1d_view_host = Tines::value_type_1d_dual_view<char[S], interf_host_device_type>;

using real_type_0d_const_view_host = typename real_type_0d_view_host::const_type;
using real_type_1d_const_view_host = typename real_type_1d_view_host::const_type;
using real_type_2d_const_view_host = typename real_type_2d_view_host::const_type;
using real_type_3d_const_view_host = typename real_type_3d_view_host::const_type;


/// utility function
using do_not_init_tag = std::string; // Kokkos::ViewAllocateWithoutInitializing;
                                     // // currently not working

template <typename T> using ats = Tines::ArithTraits<T>;

namespace Impl {
template <typename ViewType, typename MemoryTraitsType>
using ViewWithMemoryTraits = Kokkos::View<typename ViewType::data_type, typename ViewType::array_layout,
                                          typename ViewType::device_type, MemoryTraitsType>;

template <typename ViewType, typename MemoryTraitsType>
using ConstViewWithMemoryTraits = Kokkos::View<typename ViewType::const_data_type, typename ViewType::array_layout,
                                               typename ViewType::device_type, MemoryTraitsType>;

template <typename ViewType, typename MemoryTraitsType>
using ScratchViewWithMemoryTraits =
    Kokkos::View<typename ViewType::data_type, typename ViewType::array_layout,
                 typename ViewType::execution_space::scratch_memory_space, MemoryTraitsType>;

} // namespace Impl

template <typename ViewType>
using Unmanaged = Impl::ViewWithMemoryTraits<ViewType, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

template <typename ViewType> using Atomic = Impl::ViewWithMemoryTraits<ViewType, Kokkos::MemoryTraits<Kokkos::Atomic>>;

template <typename ViewType> using Const = Impl::ConstViewWithMemoryTraits<ViewType, typename ViewType::memory_traits>;

template <typename ViewType>
using ConstUnmanaged = Impl::ConstViewWithMemoryTraits<ViewType, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

template <typename ViewType>
using AtomicUnmanaged = Impl::ViewWithMemoryTraits<ViewType, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Atomic>>;

template <typename ViewType>
using Scratch = Impl::ScratchViewWithMemoryTraits<ViewType, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

template <typename T> KOKKOS_FORCEINLINE_FUNCTION T getValueInRange(const T &lo, const T &up, const T &val) {
  return val < lo ? lo : val > up ? up : val;
}


template <typename ViewType, typename T> void convertToKokkos(ViewType &out, const std::vector<T> &in) {
#if defined(TCHEM_ATM_ENABLE_DEBUG)
  std::cout << "WARNING: we highly discourage to use this utility function and "
               "recommend to use Kokkos::View directly\n";
#endif
  {
    static_assert(Kokkos::is_view<ViewType>::value, "Error: Output view is not Kokkos::View");
    static_assert(ViewType::rank == 1, "Error: Output view is not rank-1");
    static_assert(std::is_same<T, typename ViewType::non_const_value_type>::value,
                  "Error: std::vector value type does not match to Kokkos::View value "
                  "type");
    static_assert(std::is_same<typename ViewType::array_layout, Kokkos::LayoutRight>::value,
                  "Error: Output view is supposed to be layout right");
  }

  const ordinal_type n0 = in.size(), m0 = out.extent(0);
  if (n0 != m0)
    out = ViewType(Kokkos::ViewAllocateWithoutInitializing("convertStdVector1D"), n0);
  const auto out_host = Kokkos::create_mirror_view(out);
  Kokkos::parallel_for(Kokkos::RangePolicy<host_exec_space>(0, n0),
                       [&](const ordinal_type &i) { out_host(i) = in[i]; });
  Kokkos::deep_copy(out, out_host);
}


namespace Impl {

// default
#define TCHEM_ATM_LEAGUE_RANGE_DYNAMIC
//#define TCHEM_ATM_LEAGUE_RANGE_STATIC

/// decoupling team from nbatch
template <typename MemberType>
static KOKKOS_INLINE_FUNCTION void getLeagueRange(const MemberType &member, const ordinal_type n, ordinal_type &ibeg,
                                                  ordinal_type &iend, ordinal_type &iinc) {
#if defined(TCHEM_ATM_LEAGUE_RANGE_DYNAMIC)
  ibeg = member.league_rank();
  iend = n;
  iinc = member.league_size();
#elif defined(TCHEM_ATM_LEAGUE_RANGE_STATIC)
  const ordinal_type lsize = member.league_size();
  const ordinal_type istep = (n / lsize) + (n % lsize > 0);
  ibeg = member.league_rank() * istep;
  const ordinal_type itmp = ibeg + istep;
  iend = itmp < n ? itmp : n;
  iinc = 1;
#else
  TCHEM_ATM_CHECK_ERROR(member.league_size() != n, "Error: league size does not match to input");
  ibeg = member.league_rank();
  iend = ibeg + 1;
  iinc = 1;
#endif
}

/// state vector standard interface defining the internal structure
/// 0. Density (R)
/// 1. Pressure (P)
/// 2. Temperature (T)
/// 3. MassFractions (Yn)
template <typename RealType1DView> struct StateVector {
private:
  const ordinal_type _nSpec;
  const RealType1DView _v;

public:
  using real_value_type = typename RealType1DView::non_const_value_type;
  using range_type = Kokkos::pair<ordinal_type, ordinal_type>;

  KOKKOS_INLINE_FUNCTION StateVector() : _nSpec(), _v() {}
  KOKKOS_INLINE_FUNCTION StateVector(const StateVector &b) : _nSpec(b._nSpec), _v(b._v) {}
  KOKKOS_INLINE_FUNCTION StateVector(const ordinal_type nSpec, const RealType1DView &v) : _nSpec(nSpec), _v(v) {}

  /// validate input vector
  KOKKOS_INLINE_FUNCTION bool isValid() const {
    const bool is_valid_rank = (RealType1DView::rank == 1);
    const bool is_extent_valid = (_v.extent(0) <= (3 + _nSpec));
    return (is_valid_rank && is_extent_valid);
  }

  /// assign a pointer to update the vector in a batch fashion
  KOKKOS_INLINE_FUNCTION ordinal_type size() const { return _nSpec + 3; }
  KOKKOS_INLINE_FUNCTION void assign_data(real_value_type *ptr) { _v.assign_data(ptr); }

  /// copy access to private members
  KOKKOS_INLINE_FUNCTION RealType1DView KokkosView() const { return _v; }
  KOKKOS_INLINE_FUNCTION ordinal_type NumSpecies() const { return _nSpec; }

  /// interface to state vector
  KOKKOS_INLINE_FUNCTION real_value_type &Density() const { return _v(0); }
  KOKKOS_INLINE_FUNCTION real_value_type &Pressure() const { return _v(1); }
  KOKKOS_INLINE_FUNCTION real_value_type &Temperature() const { return _v(2); }
  KOKKOS_INLINE_FUNCTION auto MassFractions() const -> decltype(Kokkos::subview(_v, range_type(3, 3 + _nSpec))) {
    return Kokkos::subview(_v, range_type(3, 3 + _nSpec));
  }

  KOKKOS_INLINE_FUNCTION real_value_type *DensityPtr() const { return &_v(0); }
  KOKKOS_INLINE_FUNCTION real_value_type *PressurePtr() const { return &_v(1); }
  KOKKOS_INLINE_FUNCTION real_value_type *TemperaturePtr() const { return &_v(2); }
  KOKKOS_INLINE_FUNCTION real_value_type *MassFractionsPtr() const { return &_v(3); }
};
static KOKKOS_INLINE_FUNCTION ordinal_type getStateVectorSize(const ordinal_type nSpec) { return nSpec + 3; }
template <typename RealType1DView>
static KOKKOS_INLINE_FUNCTION StateVector<RealType1DView> wrapStateVectorView(const ordinal_type nSpec,
                                                                              const RealType1DView &v) {
  return StateVector<RealType1DView>(nSpec, v);
}

} // namespace Impl

///
/// These are only used by unit tests and examples
///
namespace Test {

template <typename ViewType> static inline void cloneView(const ViewType &v) {
  auto vp = v.data();
  if (ViewType::rank == 1) {
    const auto vs = v.stride(0);
    Kokkos::parallel_for(
        Kokkos::RangePolicy<typename ViewType::execution_space>(0, v.extent(0)),
        KOKKOS_LAMBDA(const ordinal_type i) { vp[i * vs] = vp[0]; });
  } else if (ViewType::rank == 2) {
    const auto vs0 = v.stride(0);
    const auto vs1 = v.stride(1);
    Kokkos::parallel_for(
        Kokkos::RangePolicy<typename ViewType::execution_space>(0, v.extent(0)), KOKKOS_LAMBDA(const ordinal_type i) {
          for (ordinal_type j = 0, jend = v.extent(1); j < jend; ++j)
            vp[i * vs0 + j * vs1] = vp[j * vs1];
        });
  } else {
    std::logic_error("Error: rank (>2) is not supported");
  }
}

template <typename RealType1DViewHostType>
static inline void readSiteFraction(const std::string &filename, const ordinal_type &nSpec,
                                    const RealType1DViewHostType &siteFraction) {
  std::ifstream file(filename);
  if (file.is_open()) {
    for (ordinal_type k = 0; k < nSpec; ++k)
      file >> siteFraction(k);
  } else {
    std::stringstream ss;
    ss << "Error: filename (" << filename << ") is not open\n";
    std::logic_error(ss.str());
  }
}

template <typename RealType1DViewHostType>
static inline void read1DVector(const std::string &filename, const ordinal_type &nSpec,
                                const RealType1DViewHostType &vector) {
  std::ifstream file(filename);
  if (file.is_open()) {
    printf("read1DVector: Reading data from %s\n", filename.c_str());
    for (ordinal_type k = 0; k < nSpec; ++k)
      file >> vector(k);
  } else {
    std::stringstream ss;
    ss << "Error: filename (" << filename << ") is not open\n";
    std::logic_error(ss.str());
    printf("read1DVector: Could not open %s -> Abort !\n", filename.c_str());
    exit(1);
  }
}

template <typename StringViewHostType, typename RealType2DViewHostType, typename KCMDRealType1DViewHostType>
static inline void readSample(const std::string &filename, const StringViewHostType &speciesNamesHost,
                              const KCMDRealType1DViewHostType &sMass, const ordinal_type &nSpec,
                              const ordinal_type &stateVecDim, RealType2DViewHostType &state_host, int &nBatch) {

  //
  std::ifstream file(filename);
  std::vector<std::string> varnames;
  std::vector<real_type> values;
  if (file.is_open()) {

    printf("readSample: Reading gas samples from %s\n", filename.c_str());

    // read header of file and save variable name in vector
    std::string line;
    std::getline(file, line);
    std::istringstream iss(line);
    std::string delimiter = " ";

    size_t pos = 0;
    std::string token;
    while ((pos = line.find(delimiter)) != std::string::npos) {
      varnames.push_back(line.substr(0, pos));
      line.erase(0, pos + delimiter.length());
    }
    varnames.push_back(line);
    // for (auto i = varnames.begin(); i != varnames.end(); ++i)
    //   std::cout << *i << "\n";
    // read values
    real_type value;
    while (file >> value) {
      values.push_back(value);
    }
  } else {
    printf("readSample : Could not open %s -> Abort !\n", filename.c_str());
    exit(1);
  }

  file.close();

  std::vector<ordinal_type> indx;
  for (ordinal_type sp = 2; sp < varnames.size(); sp++) {
    // convert to capital letters
    auto var_name = varnames[sp];
    std::transform(var_name.begin(), var_name.end(), var_name.begin(), ::toupper);
    // look for species name and save its index
    ordinal_type i = 0;
    for (i = 0; i < nSpec; i++) {
      if (strncmp(&speciesNamesHost(i, 0), var_name.c_str(), LENGTHOFSPECNAME) == 0) {
        indx.push_back(i);
        printf("species %s index %d \n", &speciesNamesHost(i, 0), i);
        break;
      } // end if
    }   // end for i
    // species name is not found
    if (i == nSpec) {
      printf("readSample : Species name is not part of kinetic model %s -> Abort !\n", (varnames[sp]).c_str());
      exit(1);
    } // end check last iter

  } // end sp

  nBatch = values.size() / varnames.size();
  const ordinal_type nVars = varnames.size();
  printf("Number of samples %d\n", nBatch);
  /// input: state vectors: temperature, pressure and concentration
  state_host = RealType2DViewHostType("StateVector", nBatch, stateVecDim);
  /// create a mirror view to store input from a file
  // state_host = Kokkos::create_mirror_view(state);

  Kokkos::parallel_for(Kokkos::RangePolicy<TChem::host_exec_space>(0, nBatch),
                       [state_host, nVars, &values, &indx](const ordinal_type &i) {
                         // density is set to zero

                         state_host(i, 1) = values[i * nVars + 1]; // pressure
                         state_host(i, 2) = values[i * nVars];     // temperatures

                         // mass fractions
                         for (ordinal_type sp = 0; sp < indx.size(); sp++) {
                           state_host(i, indx[sp] + 3) = values[i * nVars + 2 + sp];
                         }
                       });

  //
  // 3. compute density

  Kokkos::parallel_for(Kokkos::RangePolicy<TChem::host_exec_space>(0, nBatch), [&](const ordinal_type &i) {
    const real_type_1d_view_host state_at_i = Kokkos::subview(state_host, i, Kokkos::ALL());
    //
    const Impl::StateVector<real_type_1d_view_host> sv_at_i(nSpec, state_at_i);
    const auto Ys = sv_at_i.MassFractions();
    real_type Ysum(0);
    for (ordinal_type sp = 0; sp < indx.size(); sp++)
      Ysum += Ys(indx[sp]) / sMass(indx[sp]);

    const real_type Runiv = RUNIV * 1.0e3;
    sv_at_i.Density() = sv_at_i.Pressure() / (Runiv * Ysum * sv_at_i.Temperature());
  });
}

template <typename StringViewHostType, typename RealType2DViewHostType>
static inline void readSurfaceSample(const std::string &filename, const StringViewHostType &speciesNamesHost,
                                     const ordinal_type &nSpec, RealType2DViewHostType &sitefraction_host,
                                     int &nBatch) {

  //
  std::ifstream file(filename);
  std::vector<std::string> varnames;
  std::vector<real_type> values;
  if (file.is_open()) {

    printf("readSurfaceSample: Reading surface samples from %s\n", filename.c_str());

    // read header of file and save variable name in vector
    std::string line;
    std::getline(file, line);
    std::istringstream iss(line);
    std::string delimiter = " ";

    size_t pos = 0;
    std::string token;
    while ((pos = line.find(delimiter)) != std::string::npos) {
      varnames.push_back(line.substr(0, pos));
      line.erase(0, pos + delimiter.length());
    }
    varnames.push_back(line);
    // for (auto i = varnames.begin(); i != varnames.end(); ++i)
    //   std::cout << *i << "\n";
    // read values
    real_type value;
    while (file >> value) {
      values.push_back(value);
    }
  } else {
    printf("readSurfaceSample: Could not open %s -> Abort !\n", filename.c_str());
    exit(1);
  }
  file.close();

  std::vector<ordinal_type> indx;
  for (ordinal_type sp = 0; sp < varnames.size(); sp++) {

    auto var_name = varnames[sp];
    std::transform(var_name.begin(), var_name.end(), var_name.begin(), ::toupper);
    ordinal_type i = 0;
    for (i = 0; i < nSpec; i++) {
      if (strncmp(&speciesNamesHost(i, 0), var_name.c_str(), LENGTHOFSPECNAME) == 0) {
        indx.push_back(i);
        printf("species %s index %d \n", &speciesNamesHost(i, 0), i);
        break;
      } // end if
    }   // end for i
    // species name is not found
    if (i == nSpec) {
      printf("readSample : Species name is not part of kinetic model %s -> Abort !\n", (var_name).c_str());
      exit(1);
    } // end check last iter
  }   // end for sp

  nBatch = values.size() / varnames.size();
  const ordinal_type nVars = varnames.size();
  printf("Number of samples %d\n", nBatch);
  /// input: state vectors: temperature, pressure and concentration
  sitefraction_host = RealType2DViewHostType("Site fraction host", nBatch, nSpec);
  /// create a mirror view to store input from a file
  // sitefraction_host = Kokkos::create_mirror_view(state);

  Kokkos::parallel_for(Kokkos::RangePolicy<TChem::host_exec_space>(0, nBatch),
                       [sitefraction_host, nVars, &values, &indx](const ordinal_type &i) {
                         // mass fractions
                         for (ordinal_type sp = 0; sp < indx.size(); sp++) {
                           sitefraction_host(i, indx[sp]) = values[i * nVars + sp];
                         }
                       });
}

template <typename StringViewHostType, typename RealType1DViewHostType, typename KCMDRealType1DViewHostType>
static inline void readInput(const std::string &filename, const ordinal_type &nSpec,
                             const StringViewHostType &speciesNamesHost, const KCMDRealType1DViewHostType &sMass,
                             const RealType1DViewHostType &stateVector) {
  // read input file that contains: T, P and some species
  // T P IC8H18 O2 N2 AR
  // compute density and set state vector.

  // 1.  get data from file
  std::ifstream file(filename);
  std::vector<std::string> varnames;
  std::vector<real_type> values;
  if (file.is_open()) {

    // read header of file and save variable name in vector
    std::string line;
    std::getline(file, line);
    std::istringstream iss(line);
    std::string delimiter = " ";

    size_t pos = 0;
    std::string token;
    while ((pos = line.find(delimiter)) != std::string::npos) {
      varnames.push_back(line.substr(0, pos));
      line.erase(0, pos + delimiter.length());
    }
    varnames.push_back(line);
    // for (auto i = varnames.begin(); i != varnames.end(); ++i)
    //         std::cout << *i << "\n";
    // read values
    real_type value;
    while (file >> value) {
      values.push_back(value);
    }
  }
  file.close();

  std::vector<ordinal_type> indx;
  for (ordinal_type sp = 0; sp < varnames.size(); sp++) {
    for (ordinal_type i = 0; i < nSpec; i++) {
      if (strncmp(&speciesNamesHost(i, 0), (varnames[sp]).c_str(), LENGTHOFSPECNAME) == 0) {
        indx.push_back(i);
        printf("species %s index %d \n", &speciesNamesHost(i, 0), i);
      }
    }
  }

  const ordinal_type nBatch = values.size() / varnames.size();
  const ordinal_type nVars = varnames.size();
  printf("Number of samples %d\n", nBatch);

  // 2. fill state vector with previews data
  Impl::StateVector<RealType1DViewHostType> sv(nSpec, stateVector);

  if (sv.isValid()) {
    sv.Temperature() = values[0]; // temperature
    sv.Pressure() = values[1];    // pressure
    // sv.Density()
    const auto Ys = sv.MassFractions();
    for (ordinal_type sp = 0; sp < indx.size(); sp++)
      Ys(indx[sp]) = values[sp + 2];

  } else {
    std::logic_error("Error: stateVector is not valid");
  }
  // 3. compute density
  const auto Ys = sv.MassFractions();
  real_type Ysum(0);
  for (ordinal_type sp = 0; sp < indx.size(); sp++)
    Ysum += Ys(indx[sp]) / sMass(indx[sp]);

  const real_type Runiv = RUNIV * 1.0e3;
  sv.Density() = sv.Pressure() / (Runiv * Ysum * sv.Temperature());
}

template <typename RealType1DViewHostType>
static inline void readStateVector(const std::string &filename, const ordinal_type &nSpec,
                                   const RealType1DViewHostType &stateVector) {
  Impl::StateVector<RealType1DViewHostType> sv(nSpec, stateVector);
  if (sv.isValid()) {
    std::ifstream file(filename);
    if (file.is_open()) {
      file >> sv.Density();
      file >> sv.Pressure();
      file >> sv.Temperature();

      const auto Xc = sv.MassFractions();
      for (ordinal_type i = 0; i < nSpec; ++i)
        file >> Xc(i);
    } else {
      std::stringstream ss;
      ss << "Error: filename (" << filename << ") is not open\n";
      std::logic_error(ss.str());
      std::string error_ = "Error: filename (" + filename + ") is not open\n";
      printf("%s\n", error_.c_str());
      exit(-1);
    }
  } else {
    std::logic_error("Error: stateVector is not valid");
  }
}

template <typename RealType1DViewHostType>
static inline void readVector(const std::string &filename, const ordinal_type &nVariables,
                              const RealType1DViewHostType &vector) {
  std::ifstream file(filename);
  if (file.is_open()) {
    for (ordinal_type i = 0; i < nVariables; ++i) {
      file >> vector(i);
    }
  } else {
    std::stringstream ss;
    ss << "Error: filename (" << filename << ") is not open\n";
    std::logic_error(ss.str());
  }
}


template <typename RealType1DViewHostType>
static inline void writeReactionRates(const std::string &filename, const ordinal_type &nSpec,
                                      const RealType1DViewHostType &omega) {
  std::ofstream file(filename);
  if (file.is_open()) {
    file << std::scientific;
    file << "## nSpec(" << nSpec << ")    Omega(" << omega.extent(0) << ")\n";
    for (ordinal_type i = 0, iend = omega.extent(0); i < iend; ++i)
      file << i << "   " << omega(i) << "\n";
  } else {
    std::stringstream ss;
    ss << "Error: filename (" << filename << ") is not open\n";
    std::logic_error(ss.str());
  }
}

template <typename RealType3DViewHostType>
static inline void write3DMatrix(const std::string &filename, const RealType3DViewHostType &Matrix) {
  std::ofstream file(filename);
  if (file.is_open()) {
    file << std::scientific;
    file << "## dim 1: " << Matrix.extent(0) << "\n";
    file << "## dim 2: " << Matrix.extent(1) << "\n";
    file << "## dim 3: " << Matrix.extent(2) << "\n";

    const ordinal_type len0 = Matrix.extent(0);
    const ordinal_type len1 = Matrix.extent(1);
    const ordinal_type len2 = Matrix.extent(2);
    for (ordinal_type i = 0; i < len0; ++i) {
      for (ordinal_type j = 0; j < len1; ++j) {
        for (ordinal_type k = 0; k < len2; ++k)
          file << Matrix(i, j, k) << " ";
        // file << "\n";
      }
      file << "\n";
    }
  } else {
    std::stringstream ss;
    ss << "Error: filename (" << filename << ") is not open\n";
    std::logic_error(ss.str());
  }
}

template <typename RealType1DViewHostType>
static inline void write1DVector(const std::string &filename, const RealType1DViewHostType &Vector) {
  std::ofstream file(filename);
  if (file.is_open()) {
    file << std::scientific;
    file << "## nitems (" << Vector.extent(0) << ")\n";

    const ordinal_type len0 = Vector.extent(0);
    for (ordinal_type i = 0; i < len0; ++i) {
      file << Vector(i) << " ";
      file << "\n";
    }
  } else {
    std::stringstream ss;
    ss << "Error: filename (" << filename << ") is not open\n";
    std::logic_error(ss.str());
  }
}

template <typename RealType2DViewHostType>
static inline void write2DMatrix(const std::string &filename, const ordinal_type &nRow, const ordinal_type &nColumn,
                                 const RealType2DViewHostType &Matrix) {
  std::ofstream file(filename);
  if (file.is_open()) {
    file << std::scientific;
    file << "## nRows(" << nRow << ")    (" << Matrix.extent(0) << ")\n";
    file << "## nCols(" << nColumn << ")    (" << Matrix.extent(1) << ")\n";

    const ordinal_type len0 = Matrix.extent(0);
    const ordinal_type len1 = Matrix.extent(1);
    for (ordinal_type i = 0; i < len0; ++i) {
      for (ordinal_type j = 0; j < len1; ++j)
        file << Matrix(i, j) << " ";
      file << "\n";
    }
  } else {
    std::stringstream ss;
    ss << "Error: filename (" << filename << ") is not open\n";
    std::logic_error(ss.str());
  }
}



static inline bool compareFiles(const std::string &filename1, const std::string &filename2) {
  std::ifstream f1(filename1);
  std::string s1((std::istreambuf_iterator<char>(f1)), std::istreambuf_iterator<char>());
  std::ifstream f2(filename2);
  std::string s2((std::istreambuf_iterator<char>(f2)), std::istreambuf_iterator<char>());
  return (s1.compare(s2) == 0);
}

static inline bool compareFilesValues(const std::string &filename1, const std::string &filename2, const real_type& threshold_relative) {
  std::ifstream f1(filename1);
  std::string s1((std::istreambuf_iterator<char>(f1)), std::istreambuf_iterator<char>());
  std::ifstream f2(filename2);
  std::string s2((std::istreambuf_iterator<char>(f2)), std::istreambuf_iterator<char>());

  bool passTest(true);
  if (s1.compare(s2) != 0) {

    using ats = Tines::ats<real_type>;

    std::ifstream file1(filename1);
    std::ifstream file2(filename2);
    real_type max_relative_error(0);
    real_type max_absolute_error(0);
    real_type save_value1(0);
    real_type save_value2(0);

    if (file1.is_open() && file2.is_open()) {
      // header
      std::string line;
      std::getline(file1, line);
      std::istringstream iss1(line);

      std::string line2;
      std::getline(file2, line2);
      std::istringstream iss2(line2);

      real_type value1, value2, diff;
      while (file1 >> value1 && file2 >> value2) {

        diff = ats::abs(value1 - value2);
        if (diff > max_absolute_error) {
          max_absolute_error = diff;
          max_relative_error = diff / value2;
          save_value1 = value1;
          save_value2 = value2;
        }
      }

    } else {
      printf("test : Could not open %s -> Abort !\n", filename1.c_str());
      printf("test : Could not open %s -> Abort !\n", filename2.c_str());
      return (false);
    }
    printf("There are differences between the reference file and current file: reference: %s current: %s \n",
           filename2.c_str(), filename1.c_str());
    printf("Maximum relative error : %15.10e \n", max_relative_error);
    printf("Maximum absolute error : %15.10e \n", max_absolute_error);
    printf("Current value :  %15.10e Reference value : %15.10e \n", save_value1, save_value2);
    const real_type threshold_abs = ats::epsilon() * real_type(100);
    // check abs error
    if (max_absolute_error < threshold_abs) {
      printf("PASS with  absolute error threshold: %15.10e \n", threshold_abs);
      passTest = true;
    } else if (max_relative_error < threshold_relative)  {
      printf("PASS with relative error threshold: %15.10e \n", threshold_relative);
      passTest = true;
    } else {
      printf("FAIL with absolute error threshold: %15.10e  and relative error threshold: %15.10e \n", threshold_abs,  threshold_relative);
      passTest = false;
    }
  }

  return (passTest);
}

} // namespace Test

/** Calculate the mean speed of a gas-phase species
 * [ \f$\mbox{m}\,\mbox{s}^{-1}\f$ ] :
 * \f[
 *   v = \sqrt{\frac{8RT}{\pi MW}}
 * \f]
 * where R is the ideal gas constant [\f$\mbox{J}\, \mbox{K}^{-1}\,
 * \mbox{mol}^{-1}\f$], T is temperature [K], and MW is the molecular weight of
 * the gas-phase species
 * [\f$\mbox{kg}\, \mbox{mol}^{-1}\f$]
 * @param temperature__K Temperature [K]
 * @param mw__kg_mol Molecular weight of the gas-phase species [\f$\mbox{kg}\,
 * \mbox{mol}^{-1}\f$]
 */
 KOKKOS_INLINE_FUNCTION static
 real_type mean_speed__m_s(const real_type t,
                          const real_type mw) {
  return sqrt(8.0 * RUNIV * t / (PI() * mw));
}

  /** Calculate the mean free path of a gas-phase species [m]
 * \f[
 *   \lambda = 3.0 D_g / v
 * \f]
 * where \f$D_g\f$ is the gas-phase diffusion coefficient
 * [\f$\mbox{m}^2\,\mbox{s}^{-1}\f$] and
 * \f$v\f$ is the mean speed of the gas-phase molecules
 * [ \f$\mbox{m}\,\mbox{s}^{-1}\f$ ].
 *
 * @param diffusion_coeff__m2_s Diffusion coefficient of the gas species
 * [\f$\mbox{m}^2\, \mbox{s}^{-1}\f$]
 * @param temperature__K Temperature [K]
 * @param mw__kg_mol Molecular weight of the gas-phase species [\f$\mbox{kg}\,
 * \mbox{mol}^{-1}\f$]
 */
 KOKKOS_INLINE_FUNCTION static real_type mean_free_path_m(const real_type diffusion_coeff__m2_s,
                                       const real_type t,
                                       const real_type mw) {
  return 3.0 * diffusion_coeff__m2_s /
         mean_speed__m_s(t, mw);
}

} // namespace TChem

#endif
