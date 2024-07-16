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
#ifndef __TCHEM_ATMOSPHERIC_CHEMISTRY_E3SM_HPP__
#define __TCHEM_ATMOSPHERIC_CHEMISTRY_E3SM_HPP__

#include "TChem_KineticModelData.hpp"
#include "TChem_Util.hpp"

#include "TChem_Impl_AtmosphericChemistryE3SM.hpp"

namespace TChem {

struct AtmosphericChemistryE3SM
{

  using host_device_type = typename Tines::UseThisDevice<host_exec_space>::type;
  using device_type      = typename Tines::UseThisDevice<exec_space>::type;

  using real_type_0d_view_type = Tines::value_type_0d_view<real_type,device_type>;
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
  using real_type_2d_view_type = Tines::value_type_2d_view<real_type,device_type>;

  using real_type_0d_view_host_type = Tines::value_type_0d_view<real_type,host_device_type>;
  using real_type_1d_view_host_type = Tines::value_type_1d_view<real_type,host_device_type>;
  using real_type_2d_view_host_type = Tines::value_type_2d_view<real_type,host_device_type>;

  template<typename DeviceType>
  static inline ordinal_type getWorkSpaceSize(
    const KineticModelNCAR_ConstData<DeviceType>& kmcd)
  {
    using device_type = DeviceType;
    using problem_type = Impl::AtmosphericChemistryE3SM_Problem<real_type, device_type>;
    const ordinal_type m = problem_type::getNumberOfEquations(kmcd) + 1;

    ordinal_type work_size(0);
#if defined(TCHEM_ATM_ENABLE_SACADO_JACOBIAN_ATMOSPHERIC_CHEMISTRY)
    if (m < 32) {
      using value_type = Sacado::Fad::SLFad<real_type,32>;
      work_size = Impl::AtmosphericChemistryE3SM<value_type, device_type>::getWorkSpaceSize(kmcd);
    } else if (m < 64) {
      using value_type = Sacado::Fad::SLFad<real_type,64>;
      work_size = Impl::AtmosphericChemistryE3SM<value_type, device_type>::getWorkSpaceSize(kmcd);
    } else if (m < 128) {
      using value_type = Sacado::Fad::SLFad<real_type,128>;
      work_size = Impl::AtmosphericChemistryE3SM<value_type, device_type>::getWorkSpaceSize(kmcd);
    } else if  (m < 256) {
      using value_type = Sacado::Fad::SLFad<real_type,256>;
      work_size = Impl::AtmosphericChemistryE3SM<value_type, device_type>::getWorkSpaceSize(kmcd);
    } else if  (m < 512) {
      using value_type = Sacado::Fad::SLFad<real_type,512>;
      work_size = Impl::AtmosphericChemistryE3SM<value_type, device_type>::getWorkSpaceSize(kmcd);
    } else if (m < 1024){
      using value_type = Sacado::Fad::SLFad<real_type,1024>;
      work_size = Impl::AtmosphericChemistryE3SM<value_type, device_type>::getWorkSpaceSize(kmcd);
    } else{
      TCHEM_CHECK_ERROR(0,
                        "Error: Number of equations is bigger than size of sacado fad type");
    }
#else
    {
      work_size = Impl::AtmosphericChemistryE3SM<real_type, device_type>::getWorkSpaceSize(kmcd);
    }
#endif

    return work_size;

  }

  static void runDeviceBatch( /// thread block size
    typename UseThisTeamPolicy<exec_space>::type& policy,
    /// global tolerence parameters that governs all samples
    const real_type_1d_view_type& tol_newton,
    const real_type_2d_view_type& tol_time,
    /// sample specific input
    const real_type_2d_view_type& fac,
    const time_advance_type_1d_view& tadv,
    const real_type_2d_view_type& state,
    const real_type_2d_view_type& photo_rates,
    const real_type_2d_view_type& external_sources,
    /// output
    const real_type_1d_view_type& t_out,
    const real_type_1d_view_type& dt_out,
    const real_type_2d_view_type& state_out,
    /// const data from kinetic model
    const KineticModelNCAR_ConstData<device_type >& kmcd);



  /// tadv - an input structure for time marching
  /// state (nSpec+3) - initial condition of the state vector
  /// work - work space sized by getWorkSpaceSize
  /// t_out - time when this code exits
  /// state_out - final condition of the state vector (the same input state can
  /// be overwritten) kmcd - const data for kinetic model
  static void runHostBatch( /// input
    typename UseThisTeamPolicy<host_exec_space>::type& policy,
    /// global tolerence parameters that governs all samples
    const real_type_1d_view_host_type& tol_newton,
    const real_type_2d_view_host_type& tol_time,
    /// sample specific input
    const real_type_2d_view_host_type& fac,
    const time_advance_type_1d_view_host& tadv,
    const real_type_2d_view_host_type& state,
    const real_type_2d_view_host_type& photo_rates,
    const real_type_2d_view_host_type& external_sources,
    /// output
    const real_type_1d_view_host_type& t_out,
    const real_type_1d_view_host_type& dt_out,
    const real_type_2d_view_host_type& state_out,
    /// const data from kinetic model
    const KineticModelNCAR_ConstData<host_device_type>& kmcd);



  static void runHostBatch( /// input
    typename UseThisTeamPolicy<host_exec_space>::type& policy,
    /// global tolerence parameters that governs all samples
    const real_type_1d_view_host_type& tol_newton,
    const real_type_2d_view_host_type& tol_time,
    /// sample specific input
    const real_type_2d_view_host_type& fac,
    const time_advance_type_1d_view_host& tadv,
    const real_type_2d_view_host_type& state,
    const real_type_2d_view_host_type& photo_rates,
    const real_type_2d_view_host_type& external_sources,
    /// output
    const real_type_1d_view_host_type& t_out,
    const real_type_1d_view_host_type& dt_out,
    const real_type_2d_view_host_type& state_out,
    /// const data from kinetic model
    const Kokkos::View<KineticModelNCAR_ConstData<host_device_type>*,host_device_type>& kmcds);
  //
  static void runDeviceBatch( /// thread block size
    typename UseThisTeamPolicy<exec_space>::type& policy,
    /// global tolerence parameters that governs all samples
    const real_type_1d_view_type& tol_newton,
    const real_type_2d_view_type& tol_time,
    /// sample specific input
    const real_type_2d_view_type& fac,
    const time_advance_type_1d_view& tadv,
    const real_type_2d_view_type& state,
    const real_type_2d_view_type& photo_rates,
    const real_type_2d_view_type& external_sources,
    /// output
    const real_type_1d_view_type& t_out,
    const real_type_1d_view_type& dt_out,
    const real_type_2d_view_type& state_out,
    /// const data from kinetic model
    const Kokkos::View<KineticModelNCAR_ConstData<device_type>*,device_type>& kmcds);

};

} // namespace TChem

#endif
