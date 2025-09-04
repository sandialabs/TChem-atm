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
#ifndef __TCHEM_IMPL_REACTION_RATES_HPP__
#define __TCHEM_IMPL_REACTION_RATES_HPP__

#include "TChem_Util.hpp"
#include "TChem_KineticModelData.hpp"
#include "TChem_Impl_KForward.hpp"
#include "TChem_Impl_KForwardJPL.hpp" 
#include "TChem_Impl_AdjustReactions.hpp" 
// #define TCHEM_ENABLE_SERIAL_TEST_OUTPUT
namespace TChem {
namespace Impl {

template<typename ValueType, typename DeviceType>
struct ReactionRates
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
    return 0;
  }

  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static void team_invoke(
    const MemberType& member,
    /// input temperature
    const real_type& t,
    const real_type& p,
    //// inputs
    const value_type_1d_view_type& x,
    // input/output
    const value_type_1d_view_type& kfor,
    /// const input from kinetic model
    const kinetic_model_type& kmcd)
  {
    using kForward_type = TChem::Impl::KForward<value_type, device_type >;
    using AdjustReactions_type = TChem::Impl::AdjustReactions<value_type, device_type >;
    using kForwardJPL_type = TChem::Impl::KForwardJPL<value_type, device_type >;

    const value_type m = x(kmcd.M_index);

     // reactions constants
    kForward_type::team_invoke( member, t, p, kfor, kmcd);
    member.team_barrier();
    kForwardJPL_type::team_invoke(member, t, m, kfor, kmcd);
    member.team_barrier();
    AdjustReactions_type::team_invoke(member, t, x, kfor, kmcd);

  } // team_invoke_detail
};

} // namespace Impl
} // namespace TChem

#endif