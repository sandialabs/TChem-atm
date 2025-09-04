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
#ifndef __TCHEM_ATMOSPHERIC_CHEMISTRY_E3SM_EXPLICIT_EULER_HPP__
#define __TCHEM_ATMOSPHERIC_CHEMISTRY_E3SM_EXPLICIT_EULER_HPP__

#include "TChem_KineticModelData.hpp"
#include "TChem_Util.hpp"
#include "TChem_Impl_AtmosphericChemistryE3SM_Problem.hpp"
#include "TChem_Impl_TimeIntegratorExplicitEuler.hpp"

namespace TChem {

struct AtmosphericChemistryE3SM_ExplicitEuler
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
    using time_integrator_type = Impl::TimeIntegratorExplicitEuler<real_type, DeviceType>;

    const ordinal_type m = problem_type::getNumberOfEquations(kmcd);

    ordinal_type work_problem(0);
    work_problem = problem_type::getWorkSpaceSize(kmcd);
    ordinal_type work_time_integration(0);
    time_integrator_type::workspace(m, work_time_integration);
    return work_problem + work_time_integration;

  }

	static void 
  runHostBatch( /// thread block size
           typename UseThisTeamPolicy<host_exec_space>::type& policy,
           /// input
           const time_advance_type_1d_view_host& tadv,
           const real_type_2d_view_host& state,
           const real_type_2d_view_host& photo_rates,
           const real_type_2d_view_host& external_sources,
           /// output
           const real_type_1d_view_host& t_out,
           const real_type_1d_view_host& dt_out,
           const real_type_2d_view_host& state_out,
           /// const data from kinetic model
           const KineticModelNCAR_ConstData<interf_host_device_type>& kmcd);

  static void 
  runDeviceBatch( /// thread block size
           typename UseThisTeamPolicy<exec_space>::type& policy,
           /// input
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


};

 
} // namespace TChem

#endif