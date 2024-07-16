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
#ifndef __TCHEM_IMPL_AEROSOL_RHS_HPP__
#define __TCHEM_IMPL_AEROSOL_RHS_HPP__

#include "TChem_Impl_SIMPOL_phase_transfer.hpp"
#include "TChem_Impl_ReactionRatesAerosol.hpp"
namespace TChem {
namespace Impl {
  template<typename ValueType, typename DeviceType>
struct Aerosol_RHS
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

  using kinetic_model_data_type= KineticModelNCAR_ConstData<device_type>;

  KOKKOS_INLINE_FUNCTION static ordinal_type getWorkSpaceSize(const kinetic_model_data_type& kmcd,
                                                              const aerosol_model_data_type& amcd)
  {
    using reaction_rates_gas = TChem::Impl::ReactionRatesAerosol<real_type, device_type >;
    ordinal_type workspace_size=reaction_rates_gas::getWorkSpaceSize(kmcd);
    return workspace_size;
  }
    // update RHS using real_type, i.e., no sacado
  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static
  void team_invoke(const MemberType& member,
    const real_type& t,
    const real_type& p,
    const real_type_1d_view_type& number_conc,
    const real_type_1d_view_type& state,
    const real_type_1d_view_type& const_X,
    const real_type_1d_view_type& omega,
    const real_type_1d_view_type& work,
    const kinetic_model_data_type& kmcd,
    const aerosol_model_data_type& amcd
    )
  {
   // set omega(rhs) to zero, because we are using Kokkos::atomic_add.
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(member, omega.extent(0)),
       [&](const ordinal_type& i) {
      omega(i) = 0.0;
    });
    // 1. update RHS of gas species
    using reaction_rates_gas = TChem::Impl::ReactionRatesAerosol<real_type, device_type >;
    reaction_rates_gas::team_invoke(member,
                         t,
                         p,
                         state,
                         const_X,
                         omega,
                         // work arrays
                         work,
                         kmcd);
    // 2. update RHS of gas and aerosol species
    member.team_barrier();
    using SIMPOL_single_particle_type = TChem::Impl::SIMPOL_single_particle<real_type, device_type >;
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(member, amcd.nParticles),
       [&](const ordinal_type& i_part) {
    for (size_t i_simpol = 0; i_simpol < amcd.nSimpol_tran; i_simpol++)
    {
    SIMPOL_single_particle_type
    ::team_invoke(member, i_part,i_simpol,
                  t, p, number_conc,
                  state, omega,
                  amcd);
    }// i_simpol
    });// i_part

#if defined(TCHEM_ENABLE_SERIAL_TEST_OUTPUT)
  printf("omega.extent(0) %d \n",omega.extent(0));
  printf("---RHSs--\n");
  printf("omega(%d) %e \n",0,omega(0));
  for (ordinal_type i_part = 0; i_part < amcd.nParticles; i_part++)
  {
    ordinal_type is = amcd.nSpec_gas + i_part*amcd.nSpec;
    for (ordinal_type i = 0; i < amcd.nSpec; i++)
    {
      printf("omega(%d) %e \n",is+i,omega(is+i));
    }
  }
#endif
  }



};
} // namespace Impl
} // namespace TChem

#endif
