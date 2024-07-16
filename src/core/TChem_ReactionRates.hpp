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
#ifndef __TCHEM_REACTION_RATES_HPP__
#define __TCHEM_REACTION_RATES_HPP__

#include "TChem_KineticModelData.hpp"
#include "TChem_Util.hpp"
#include "TChem_Impl_ReactionRates.hpp"

namespace TChem {

struct ReactionRates
{
  using host_device_type = typename Tines::UseThisDevice<host_exec_space>::type;
  using device_type      = typename Tines::UseThisDevice<exec_space>::type;

  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
  using real_type_2d_view_type = Tines::value_type_2d_view<real_type,device_type>;

  using real_type_1d_view_host_type = Tines::value_type_1d_view<real_type,host_device_type>;
  using real_type_2d_view_host_type = Tines::value_type_2d_view<real_type,host_device_type>;

  using kinetic_model_type = KineticModelNCAR_ConstData<device_type>;
  using kinetic_model_host_type = KineticModelNCAR_ConstData<host_device_type>;

  template<typename DeviceType>
  static inline ordinal_type getWorkSpaceSize(
    const KineticModelNCAR_ConstData<DeviceType>& kmcd)
  {
    return Impl::ReactionRates<real_type,DeviceType>::getWorkSpaceSize(kmcd);
  }

  static void runHostBatch( /// input
    const real_type_2d_view_host_type& state,
    /// output
    const real_type_2d_view_host_type& kfwd,
    /// const data from kinetic model
    const kinetic_model_host_type& kmcd);

  static void runDeviceBatch( /// input
    typename UseThisTeamPolicy<exec_space>::type& policy,
    const real_type_2d_view_type& state,
    /// output
    const real_type_2d_view_type& kfwd,
    /// const data from kinetic model
    const kinetic_model_type& kmcd);
  //
  static void runDeviceBatch( /// input
    const real_type_2d_view_type& state,
    /// output
    const real_type_2d_view_type& kfwd,
    /// const data from kinetic model
    const kinetic_model_type& kmcd);

  //

};

} // namespace TChem

#endif
