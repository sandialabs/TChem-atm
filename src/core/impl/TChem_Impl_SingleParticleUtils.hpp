/* =====================================================================================
TChem-atm version 1.0
Copyright (2024) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC
(NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
Government retains certain rights in this software.

This file is part of TChem-atm. TChem-atm is open source software: you can redistribute
it and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Contact Oscar Diaz-Ibarra at <odiazib@sandia.gov>, or
           Mike Schmidt at <mjschm@sandia.gov>, or
           Cosmin Safta at <csafta@sandia.gov>

Sandia National Laboratories, New Mexico/Livermore, NM/CA, USA
=====================================================================================
*/
#ifndef __TCHEM_IMPL_SINGLE_PARTICLE_UTILS_HPP__
#define __TCHEM_IMPL_SINGLE_PARTICLE_UTILS_HPP__

/// Minimum aerosol-phase mass concentration [kg m-3]
#define MINIMUM_MASS_ 1.0e-25L
/// Minimum mass assumed molecular weight [kg mol-1]
#define MINIMUM_MW_ 0.1L
/// Minimum mass assumed density [kg m-3]
#define MINIMUM_DENSITY_ 1800.0L

#define I_SPEC_IN_PART amcd.nSpec_gas + i_spec + i_part*amcd.nSpec

#include "TChem_Util.hpp"

namespace TChem {
namespace Impl {
  template<typename ValueType, typename DeviceType>
struct Single_particle_utils
{
  using value_type = ValueType;
  using device_type = DeviceType;
  using scalar_type = typename ats<value_type>::scalar_type;

  using real_type = scalar_type;
  /// sacado is value type
  using value_type_1d_view_type = Tines::value_type_1d_view<value_type,device_type>;
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
  using ordinal_type_1d_view_type = Tines::value_type_1d_view<ordinal_type,device_type>;
  using aerosol_model_data_type= AerosolModel_ConstData<device_type>;

  using const_real_type_1d_view_type = Tines::value_type_1d_view<const real_type,device_type>;

 template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static void
  aero_phase_get_volume(const MemberType& member,
                        const ordinal_type i_part,
                        const value_type_1d_view_type& state,
                        value_type& volume,
                        const aerosol_model_data_type& amcd)
   {
     // Sum the mass and MW
   volume = MINIMUM_MASS_ / MINIMUM_DENSITY_;
  for (int i_spec = 0; i_spec < amcd.nSpec; i_spec++) {
      volume += state(I_SPEC_IN_PART) / amcd.aerosol_density(i_spec);
  }

  }
  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static void effective_radius(
    const MemberType& member,
    const ordinal_type i_part,
    const real_type& t,
    const real_type& p,
    const value_type_1d_view_type& state,
    value_type& radius,
    const aerosol_model_data_type& amcd
    )
    {
    radius=0.0;
    // for (int i_phase = 0; i_phase < num_of_phases; ++i_phase) {
    value_type volume=0.0;
    aero_phase_get_volume(member, i_part, state, volume, amcd );
    radius += volume;
    // }//
    radius = ats<value_type>::pow((radius * 3.0 / 4.0 / PI()), 1.0 / 3.0);
    }
  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static void
  aero_phase_get_mass(const MemberType& member,
                             const ordinal_type i_part,
                             const value_type_1d_view_type& state,
                             value_type&MW,
                             value_type&mass,
                             const aerosol_model_data_type& amcd
                             )
  {
     // Sum the mass and MW
  mass = MINIMUM_MASS_;
  value_type moles = MINIMUM_MASS_ / MINIMUM_MW_;
  for (int i_spec = 0; i_spec < amcd.nSpec; i_spec++) {
      mass += state(I_SPEC_IN_PART);
      moles += state(I_SPEC_IN_PART) / amcd.molecular_weights(i_spec);
    // }
  }
  MW = mass / moles;
  }

  /** Calculate the transition regime correction factor [unitless] \cite Fuchs1971
 * : \f[ f(K_n,\alpha) = \frac{0.75 \alpha ( 1 + K_n )}{K_n(1+K_n) + 0.283\alpha
 * K_n + 0.75 \alpha} \f] where the Knudsen Number \f$K_n = \lambda / r\f$
 * (where \f$\lambda\f$ is the mean free path [m] of the gas-phase species and
 * \f$r\f$ is the effective radius of the particles [m]), and \f$ \alpha \f$ is
 * the mass accomodation coefficient, which is typically assumed to equal 0.1
 * \cite Zaveri2008.
 *
 *  @param mean_free_path__m mean free path of the gas-phase species [m]
 *  @param radius__m Particle effective radius [m]
 *  @param alpha Mass accomodation coefficient [unitless]
 */

KOKKOS_INLINE_FUNCTION static
value_type
transition_regime_correction_factor(real_type mean_free_path__m,
                                    value_type radius__m,
                                    real_type alpha) {
  value_type K_n = mean_free_path__m / radius__m;
  return (0.75 * alpha * (1.0 + K_n)) /
         (K_n * K_n + (1.0 + 0.283 * alpha) * K_n + 0.75 * alpha);
}

/** Calculate the gas-aerosol reaction rate constant for the transition regime
 * [\f$\mbox{m}^3\, \mbox{particle}^{-1}\, \mbox{s}^{-1}\f$]
 *
 *  The rate constant \f$k_c\f$ is calculated according to \cite Zaveri2008 as:
 *  \f[
 *    k_c = 4 \pi r D_g f_{fs}( K_n, \alpha )
 *  \f]
 *  where \f$r\f$ is the radius of the particle(s) [m], \f$D_g\f$ is the
 * diffusion coefficient of the gas-phase species
 * [\f$\mbox{m}^2\,\mbox{s}^{-1}\f$], \f$f_{fs}( K_n, \alpha )\f$ is the Fuchs
 * Sutugin transition regime correction factor [unitless], \f$K_n\f$ is the
 * Knudsen Number [unitess], and \f$\alpha\f$ is the mass accomodation
 * coefficient.
 *
 *  Rates can be calculated as:
 *  \f[
 *    r_c = [G] N_a k_c
 *  \f]
 *  where \f$[G]\f$ is the gas-phase species concentration [ppm], \f$N_a\f$ is
 * the number concentration of particles [\f$\mbox{particle}\,\mbox{m}^{-3}\f$]
 * and the rate \f$r_c\f$ is in [\f$\mbox{ppm}\,\mbox{s}^{-1}\f$].
 *
 *  @param diffusion_coeff__m2_s Diffusion coefficent of the gas species
 *  [\f$\mbox{m}^2\, \mbox{s}^{-1}\f$]
 *  @param mean_free_path__m Mean free path of gas molecules [m]
 *  @param radius__m Particle radius [m]
 *  @param alpha Mass accomodation coefficient [unitless]
 */
KOKKOS_INLINE_FUNCTION
static value_type gas_aerosol_transition_rxn_rate_constant(
    real_type diffusion_coeff__m2_s, real_type mean_free_path__m,
    value_type radius__m, real_type alpha) {
  return 4.0 * PI() * radius__m * diffusion_coeff__m2_s *
         transition_regime_correction_factor(mean_free_path__m, radius__m,
                                             alpha);
}

};

} // namespace Impl
} // namespace TChem

#endif
