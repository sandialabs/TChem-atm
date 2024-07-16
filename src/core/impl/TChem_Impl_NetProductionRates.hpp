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
#ifndef __TCHEM_IMPL_NET_PRODUCTION_RATES_HPP__
#define __TCHEM_IMPL_NET_PRODUCTION_RATES_HPP__

#include "TChem_Util.hpp"
#include "TChem_KineticModelData.hpp"
#include "TChem_Impl_RateofProgress.hpp"
namespace TChem {
namespace Impl {

template<typename ValueType, typename DeviceType>
struct NetProductionRates
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
    return 2*kmcd.nReac + kmcd.nSpec;
  }

  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static void team_invoke_detail(
    const MemberType& member,
    /// input temperature
    const real_type& t,
    const real_type& p,
    //// inputs
    const value_type_1d_view_type& x,
    /// photo rates
    const real_type_1d_view_type& photo_rates,
    // external sources
    const real_type_1d_view_type& external_sources,
    /// output
    const value_type_1d_view_type& net_production_rate, /// (kmcd.nSpec)
    // work
    const value_type_1d_view_type& kfor,
    const value_type_1d_view_type& rate_of_progress,
    /// const input from kinetic model
    const kinetic_model_type& kmcd)
  {
  	const ordinal_type n_active_vars = kmcd.nSpec - kmcd.nConstSpec;

    // set net production rate to be equal to external sources.
    Kokkos::parallel_for(
      Tines::RangeFactory<value_type>::TeamVectorRange(member, n_active_vars), [&](const ordinal_type& i) {
      net_production_rate(i) = external_sources(i);
    });

    // compute rate of progress
  	using rateof_progress_type = TChem::Impl::RateofProgress<value_type, device_type >;
    rateof_progress_type::team_invoke(member, t, p, x, photo_rates, rate_of_progress, kfor, kmcd);
    member.team_barrier();

    // compute net production rate, i.e. RHS or omega
    auto rop = rate_of_progress;
    Kokkos::parallel_for(
      Tines::RangeFactory<value_type>::TeamVectorRange(member, kmcd.nReac), [=](const ordinal_type& i) {
        const value_type rop_at_i = rop(i);
        for (ordinal_type j = 0; j < kmcd.reacNreac(i); ++j) {
          const ordinal_type kspec = kmcd.reacSidx(i, j);
          // do not compute const-tracer species
          if (kspec < n_active_vars){
            const value_type val = kmcd.reacNuki(i, j) * rop_at_i;
            Kokkos::atomic_add(&net_production_rate(kspec), val);
          }

        }
        const ordinal_type joff = kmcd.reacSidx.extent(1) / 2;
        for (ordinal_type j = 0; j < kmcd.reacNprod(i); ++j) {
          const ordinal_type kspec = kmcd.reacSidx(i, j + joff);
          // do not compute cons-tracer species
          if (kspec < n_active_vars){
            const value_type val = kmcd.reacNuki(i, j + joff) * rop_at_i;
            Kokkos::atomic_add(&net_production_rate(kspec), val);
          }
        }
      });

    member.team_barrier();


  } // team_invoke

    template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static void team_invoke_sacado(
    const MemberType& member,
    /// input temperature
    const real_type& t,
    const real_type& p,
    //// inputs
    const value_type_1d_view_type& X,
    const real_type_1d_view_type& photo_rates,
    const real_type_1d_view_type& external_sources,
    const real_type_1d_view_type& const_X,
    /// output
    const value_type_1d_view_type& net_production_rate, /// (kmcd.nSpec)
    // work
    const real_type_1d_view_type& work,
    /// const input from kinetic model
    const kinetic_model_type& kmcd)
    {
      auto w = (real_type*)work.data();
      const ordinal_type len = ats<value_type>::sacadoStorageCapacity();
      // do not use either t or p, because they are real_type
      const ordinal_type sacadoStorageDimension = ats<value_type>::sacadoStorageDimension(X(0));

      auto ropFor = value_type_1d_view_type(w, kmcd.nReac, sacadoStorageDimension);
      w += kmcd.nReac*len;
      auto kfor = value_type_1d_view_type(w, kmcd.nReac, sacadoStorageDimension);
      w += kmcd.nReac*len;
      auto concX = value_type_1d_view_type(w, kmcd.nSpec, sacadoStorageDimension);
      w += kmcd.nSpec*len;

      const ordinal_type n_active_vars = kmcd.nSpec-kmcd.nConstSpec;

      Kokkos::parallel_for(
        Tines::RangeFactory<value_type>::TeamVectorRange(member, n_active_vars ),
         [&](const ordinal_type& i) {
        concX(i) = X(i);
      });
      // constant variables
      Kokkos::parallel_for(
        Tines::RangeFactory<value_type>::TeamVectorRange(member, kmcd.nConstSpec),
         [=](const ordinal_type& i) {
        concX(i + n_active_vars) = const_X(i);
      });
      member.team_barrier();

      team_invoke_detail(member,
                         t,
                         p,
                         concX,
                         photo_rates,
                         external_sources,
                         net_production_rate,
                         // work arrays
                         kfor,
                         ropFor,
                         kmcd);

    }

      template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static void team_invoke_sacado(
    const MemberType& member,
    /// input temperature
    const real_type& t,
    const real_type& p,
    //// inputs
    const value_type_1d_view_type& X,
    const real_type_1d_view_type& photo_rates,
    const real_type_1d_view_type& external_sources,
    /// output
    const value_type_1d_view_type& net_production_rate, /// (kmcd.nSpec)
    // work
    const real_type_1d_view_type& work,
    /// const input from kinetic model
    const kinetic_model_type& kmcd)
    {

      auto w = (real_type*)work.data();
      const ordinal_type len = ats<value_type>::sacadoStorageCapacity();
      // do not use either t or p, because they are real_type
      const ordinal_type sacadoStorageDimension = ats<value_type>::sacadoStorageDimension(X(0));

      auto ropFor = value_type_1d_view_type(w, kmcd.nReac, sacadoStorageDimension);
      w += kmcd.nReac*len;
      auto kfor = value_type_1d_view_type(w, kmcd.nReac, sacadoStorageDimension);
      w += kmcd.nReac*len;
      member.team_barrier();

      team_invoke_detail(member,
                         t,
                         p,
                         X,
                         photo_rates,
                         external_sources,
                         net_production_rate,
                         kfor,
                         ropFor,
                         kmcd);

    }
      template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static void team_invoke(
    const MemberType& member,
    /// input temperature
    const real_type& t,
    const real_type& p,
    //// inputs
    const real_type_1d_view_type& X, // molar concentration of active species
    const real_type_1d_view_type& photo_rates,
    const real_type_1d_view_type& external_sources,
    //// inputs
    const real_type_1d_view_type& const_X, // molar concentration of constant species
    /// output
    const real_type_1d_view_type& net_production_rate, /// (kmcd.nSpec)
    // work
    const real_type_1d_view_type& work,
    /// const input from kinetic model
    const kinetic_model_type& kmcd)
    {

      auto w = (real_type*)work.data();
      auto ropFor = real_type_1d_view_type(w, kmcd.nReac);
      w += kmcd.nReac;
      auto kfor = real_type_1d_view_type(w, kmcd.nReac);
      w += kmcd.nReac;

      auto concX = real_type_1d_view_type(w, kmcd.nSpec);
      w += kmcd.nSpec;

      const ordinal_type n_active_vars = kmcd.nSpec-kmcd.nConstSpec;

      Kokkos::parallel_for(
          Tines::RangeFactory<real_type>::TeamVectorRange(member, n_active_vars ),
          [&](const ordinal_type& i) {
            concX(i) = X(i);
      });
      // constant variables
      Kokkos::parallel_for(
        Tines::RangeFactory<real_type>::TeamVectorRange(member, kmcd.nConstSpec),
          [=](const ordinal_type& i) {
            const ordinal_type idx(i + n_active_vars);
            concX(idx) = const_X(i);
          });
      member.team_barrier();

      team_invoke_detail(member,
                         t,
                         p,
                         concX,
                         photo_rates,
                         external_sources,
                         net_production_rate,
                         kfor,
                         ropFor,
                         kmcd);

    }

  // source term assuming all species are not constant
  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static void team_invoke(
    const MemberType& member,
    /// input temperature
    const real_type& t,
    const real_type& p,
    //// inputs
    const real_type_1d_view_type& X, // molar concentration of active species
    const real_type_1d_view_type& photo_rates,
    const real_type_1d_view_type& external_sources,
    /// output
    const real_type_1d_view_type& net_production_rate, /// (kmcd.nSpec)
    // work
    const real_type_1d_view_type& work,
    /// const input from kinetic model
    const kinetic_model_type& kmcd)
    {

      auto w = (real_type*)work.data();
      auto ropFor = real_type_1d_view_type(w, kmcd.nReac);
      w += kmcd.nReac;
      auto kfor = real_type_1d_view_type(w, kmcd.nReac);
      w += kmcd.nReac;

      team_invoke_detail(member,
                         t,
                         p,
                         X,
                         photo_rates,
                         external_sources,
                         net_production_rate,
                         kfor,
                         ropFor,
                         kmcd);

    }

};//NetProductionRates

} // namespace Impl
} // namespace TChem
#endif
