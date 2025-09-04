/* =====================================================================================
TChem-atm version 2.0.0
Copyright (2025) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2025 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.

This file is part of TChem-atm. TChem-atm is open source software: you can redistribute it
and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Contact Oscar Diaz-Ibarra at <odiazib@sandia.gov>, or
           Cosmin Safta at <csafta@sandia.gov> or,
           Nicole Riemer at <nriemer@illinois.edu> or,
           Matthew West at <mwest@illinois.edu>

Sandia National Laboratories, New Mexico/Livermore, NM/CA, USA
===================================================================================== */
#ifndef __TCHEM_IMPL_RATE_OF_PROCESS_HPP__
#define __TCHEM_IMPL_RATE_OF_PROCESS_HPP__

#include "TChem_Util.hpp"
#include "TChem_KineticModelData.hpp"
#include "TChem_Impl_ReactionRates.hpp"
// #define TCHEM_ENABLE_SERIAL_TEST_OUTPUT
namespace TChem {
namespace Impl {

template<typename ValueType, typename DeviceType>
struct RateofProgress
{
  using value_type = ValueType;
  using device_type = DeviceType;
  using scalar_type = typename ats<value_type>::scalar_type;

  using real_type = scalar_type;
  /// sacado is value type
  using value_type_1d_view_type = Tines::value_type_1d_view<value_type,device_type>;
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
  using kinetic_model_type= KineticModelNCAR_ConstData<device_type>;

  KOKKOS_INLINE_FUNCTION static ordinal_type getWorkSpaceSize(
    const kinetic_model_type& kmcd)
  {
    return kmcd.nReac;
  }

  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static void team_invoke(
    const MemberType& member,
    /// input temperature
    const real_type& t,
    const real_type& p,
    //// inputs
    const value_type_1d_view_type& concX,
    /// photo rates
    const real_type_1d_view_type& photo_rates,
    /// output
    const value_type_1d_view_type& rate_of_progress, /// (kmcd.nSpec)
    // work
    const value_type_1d_view_type& kfor,
    /// const input from kinetic model
    const kinetic_model_type& kmcd)
  {


    // set photo rates:
    const ordinal_type n_photo_rates = photo_rates.extent(0);
     Kokkos::parallel_for(
      Tines::RangeFactory<value_type>::TeamVectorRange(member, n_photo_rates), [&](const ordinal_type& i) {
      kfor(i)=photo_rates(i);
     });
    member.team_barrier();

    // compute reaction constants
    using reation_rates_type = TChem::Impl::ReactionRates<value_type, device_type >;
    reation_rates_type::team_invoke( member, t, p, concX, kfor, kmcd);


   
    member.team_barrier();    
    // compute rate of progresses
    Kokkos::parallel_for(
      Tines::RangeFactory<value_type>::TeamVectorRange(member, kmcd.nReac), [&](const ordinal_type& i) {
        value_type ropFor_at_i = kfor(i);

        for (ordinal_type j = 0; j < kmcd.reacNreac(i); ++j) {
          const ordinal_type kspec = kmcd.reacSidx(i, j);
          const ordinal_type niup = ats<ordinal_type>::abs(kmcd.reacNuki(i, j));
          ropFor_at_i *= ats<value_type>::pow(concX(kspec), niup);
        }

        rate_of_progress(i) = ropFor_at_i;
        // printf("reaction i %d ropFor %e \n", i, ropFor(i));
    });

    member.team_barrier();

  }
};

} // namespace Impl
} // namespace TChem

#endif
