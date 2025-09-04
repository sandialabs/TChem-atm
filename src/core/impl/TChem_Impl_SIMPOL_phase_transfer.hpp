/* =====================================================================================
TChem-atm version 2.0.0
Copyright (2025) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2025 National Technology & Engineering Solutions of Sandia, LLC
(NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
Government retains certain rights in this software.

This file is part of TChem-atm. TChem-atm is open source software: you can redistribute
it and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Contact Oscar Diaz-Ibarra at <odiazib@sandia.gov>, or
           Cosmin Safta at <csafta@sandia.gov> or,
           Nicole Riemer at <nriemer@illinois.edu> or,
           Matthew West at <mwest@illinois.edu>

Sandia National Laboratories, New Mexico/Livermore, NM/CA, USA
=====================================================================================
*/
#ifndef __TCHEM_IMPL_SIMPOL_PHASE_TRANSFER_HPP__
#define __TCHEM_IMPL_SIMPOL_PHASE_TRANSFER_HPP__

#include "TChem_Util.hpp"
#include "TChem_ReactionTypes.hpp"
#include "TChem_Impl_SIMPOL_constant.hpp"
#include "TChem_Impl_SingleParticleUtils.hpp"

namespace TChem {
namespace Impl {
  template<typename ValueType, typename DeviceType>
struct SIMPOL_single_particle
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

  KOKKOS_INLINE_FUNCTION static ordinal_type getWorkSpaceSize()
  {
    ordinal_type workspace_size=0;
    return workspace_size;
  }
  // update RHS
  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static
  void team_invoke(const MemberType& member,
    const ordinal_type i_part,
    const ordinal_type i_simpol,
    const real_type& t,
    const real_type& p,
    const real_type_1d_view_type& number_conc,
    const value_type_1d_view_type& state,
    const value_type_1d_view_type& omega,
    const aerosol_model_data_type& amcd
    )
  {

  const SIMPOL_PhaseTransferType simpol_params = amcd.simpol_params(i_simpol);
  ordinal_type GAS_SPEC_=simpol_params.gas_sp_index;
  ordinal_type AERO_SPEC_i_phase=simpol_params.aerosol_sp_index;
  real_type DIFF_COEFF_=simpol_params.diffusion_coeff;
  real_type MW_=simpol_params.molecular_weight;

  using SIMPOL_constant_type
   = TChem::Impl::SIMPOL_constant<real_type, device_type >;

  using single_part_utils_type = Single_particle_utils<value_type, device_type >;

     // aerosol phase (m3/#/s)
  real_type mfp_m =0.0;
  real_type alpha= 0.0;
  real_type EQUIL_CONST_=0.0;
  real_type KGM3_TO_PPM_=0.0;
  SIMPOL_constant_type::team_invoke( member, t, p, alpha,
    mfp_m, KGM3_TO_PPM_, EQUIL_CONST_, simpol_params);

  // Compute radius
  value_type radius =1;
  single_part_utils_type::effective_radius(member, i_part, t, p,
                  state,
                  radius,
                  amcd);


   // Get the particle number concentration (#/m3)
   // Get the total mass of the aerosol phase (kg/mol)
   value_type aero_phase_avg_MW=0;
   // Get the total mass of the aerosol phase (kg/m3)
   value_type aero_phase_mass=0;
   single_part_utils_type::aero_phase_get_mass(member,
                       i_part,
                       state,
                       aero_phase_avg_MW,
                       aero_phase_mass,
                       amcd );


   // Calculate the rate constant for diffusion limited mass transfer to the

   value_type cond_rate =
   single_part_utils_type::gas_aerosol_transition_rxn_rate_constant(DIFF_COEFF_,
    mfp_m, radius, alpha);

  //FIXME: to be done
  // Get the activity coefficient (if one exists)
  value_type act_coeff = 1.0;
  // if (AERO_ACT_ID_(i_phase) > -1) {
  //   act_coeff = state[AERO_ACT_ID_(i_phase)];
  // }
  // Calculate the evaporation rate constant (ppm_x*m^3/kg_x/s)
  value_type evap_rate =
        cond_rate *  act_coeff * (EQUIL_CONST_ * aero_phase_avg_MW / aero_phase_mass);

  // Calculate the evaporation and condensation rates
  cond_rate *= state(GAS_SPEC_);
  evap_rate *= state(AERO_SPEC_i_phase+i_part*amcd.nSpec);
  // // per-particle mass concentration rates
  value_type diff = - evap_rate + cond_rate;

#if defined(TCHEM_ENABLE_SERIAL_TEST_OUTPUT)
  printf("radius %e \n", radius);
  printf("aero_phase_mass %e \n", aero_phase_mass);
  printf("aero_phase_avg_MW %e \n", aero_phase_avg_MW);
  printf("cond_rate %e \n", cond_rate);
  printf("act_coeff %e \n", act_coeff);
  printf("EQUIL_CONST_ %e \n", EQUIL_CONST_);
  printf("evap_rate %e \n", evap_rate);
  printf("state(GAS_SPEC_) %e \n", state(GAS_SPEC_));
  printf("state(%d) %e \n",AERO_SPEC_i_phase, state(AERO_SPEC_i_phase+i_part*amcd.nSpec));
  printf("diff %e \n",diff);
  printf("A evap_rate %e \n", evap_rate);
  printf("A cond_rate %e \n", cond_rate);
#endif
  // Change in the gas-phase is evaporation - condensation (ppm/s)
  Kokkos::atomic_add(&omega(GAS_SPEC_), -number_conc(i_part) * diff);
  Kokkos::atomic_add(&omega(AERO_SPEC_i_phase+i_part*amcd.nSpec), diff / KGM3_TO_PPM_);
  }

};

} // namespace Impl
} // namespace TChem

#endif
