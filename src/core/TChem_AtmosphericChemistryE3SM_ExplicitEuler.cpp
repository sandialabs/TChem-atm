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
#include "TChem_Util.hpp"
#include "TChem_AtmosphericChemistryE3SM_ExplicitEuler.hpp"

namespace TChem
{
	  template<typename PolicyType,
           typename DeviceType>
  void
  AtmosphericChemistryE3SM_ExplicitEuler_TemplateRun( /// required template arguments
          const std::string& profile_name,
          /// team size setting
          const PolicyType& policy,
          /// input
          const Tines::value_type_1d_view<time_advance_type, DeviceType>& tadv,
          const Tines::value_type_2d_view<real_type, DeviceType>& state,
          const Tines::value_type_2d_view<real_type, DeviceType>& photo_rates,
          const Tines::value_type_2d_view<real_type, DeviceType>& external_sources,
          /// output
          const Tines::value_type_1d_view<real_type, DeviceType>& t_out,
          const Tines::value_type_1d_view<real_type, DeviceType>& dt_out,
          const Tines::value_type_2d_view<real_type, DeviceType>& state_out,
          /// const data from kinetic model
          const KineticModelNCAR_ConstData<DeviceType>& kmcd) {
    Kokkos::Profiling::pushRegion(profile_name);

    using time_integrator_type = TChem::Impl::TimeIntegratorExplicitEuler<real_type, DeviceType>;
    using policy_type = PolicyType;
    using real_type_0d_view_type = Tines::value_type_0d_view<real_type, DeviceType>;
    using real_type_1d_view_type = Tines::value_type_1d_view<real_type, DeviceType>;
    using range_type = Kokkos::pair<ordinal_type, ordinal_type>;
    using problem_type = TChem::Impl::AtmosphericChemistryE3SM_Problem<real_type,DeviceType>;

    const ordinal_type level = 1;
    const ordinal_type per_team_extent = AtmosphericChemistryE3SM_ExplicitEuler::getWorkSpaceSize(kmcd);

    /// this assumes host space
    Kokkos::parallel_for
      (profile_name,
       policy,
       KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
  const ordinal_type i = member.league_rank();

  const auto tadv_at_i = tadv(i);
  const real_type t_end = tadv_at_i._tend;
  const real_type_0d_view_type t_out_at_i = Kokkos::subview(t_out, i);
  if (t_out_at_i() < t_end) {
    const real_type_1d_view_type state_at_i =
        Kokkos::subview(state, i, Kokkos::ALL());

      // Note: The number of photo reactions can be equal to zero.
      real_type_1d_view_type photo_rates_at_i;
      if(photo_rates.extent(0) > 0)
      {
        photo_rates_at_i =
        Kokkos::subview(photo_rates, i, Kokkos::ALL());
      }

      // Note: The number of external sources can be equal to zero.
      real_type_1d_view_type external_sources_at_i;
      if(external_sources.extent(0) > 0)
      {
        external_sources_at_i= Kokkos::subview(external_sources, i, Kokkos::ALL());
      }

    const real_type_1d_view_type state_out_at_i =
        Kokkos::subview(state_out, i, Kokkos::ALL());

    const real_type_0d_view_type dt_out_at_i = Kokkos::subview(dt_out, i);
      Scratch<real_type_1d_view_type> work(member.team_scratch(level),
                                       per_team_extent);
      Impl::StateVector<real_type_1d_view_type> sv_at_i(kmcd.nSpec, state_at_i);
      Impl::StateVector<real_type_1d_view_type> sv_out_at_i(kmcd.nSpec,
                                                        state_out_at_i);
      TCHEM_CHECK_ERROR(!sv_at_i.isValid(),
                        "Error: input state vector is not valid");
      TCHEM_CHECK_ERROR(!sv_out_at_i.isValid(),
                        "Error: input state vector is not valid");
    {
      const ordinal_type max_num_time_iterations = tadv_at_i._num_time_iterations_per_interval;
      const real_type
        dt_in = tadv_at_i._dt,
        dt_min = tadv_at_i._dtmin,
        dt_max = tadv_at_i._dtmax;
      const real_type t_beg = tadv_at_i._tbeg;

      const real_type temperature = sv_at_i.Temperature();
      const real_type pressure = sv_at_i.Pressure();
      const real_type density = sv_at_i.Density();
      const real_type_1d_view_type Ys = sv_at_i.MassFractions();
      const auto activeYs = Kokkos::subview(Ys,
          range_type(0, kmcd.nSpec - kmcd.nConstSpec));
      const auto constYs = Kokkos::subview(Ys,
          range_type(kmcd.nSpec - kmcd.nConstSpec, kmcd.nSpec));

      const real_type_0d_view_type temperature_out(sv_out_at_i.TemperaturePtr());
      const real_type_0d_view_type pressure_out(sv_out_at_i.PressurePtr());
      const real_type_1d_view_type Ys_out = sv_out_at_i.MassFractions();
      const auto activeYs_out = Kokkos::subview(Ys_out,
          range_type(0, kmcd.nSpec - kmcd.nConstSpec));

      const real_type_0d_view_type density_out(sv_out_at_i.DensityPtr());

      const ordinal_type m = problem_type::getNumberOfTimeODEs(kmcd);
      /// problem setup
      const ordinal_type problem_workspace_size = problem_type::getWorkSpaceSize(kmcd);

      problem_type problem;
      problem._kmcd = kmcd;

      /// problem workspace
      auto wptr = work.data();
      auto pw = real_type_1d_view_type(wptr, problem_workspace_size);
      wptr += problem_workspace_size;

      /// error check
      const ordinal_type workspace_used(wptr - work.data()), workspace_extent(work.extent(0));
      if (workspace_used > workspace_extent) {
        Kokkos::abort("Error Ignition ZeroD Sacado : workspace used is larger than it is provided\n");
      }

      /// time integrator workspace
      auto tw = real_type_1d_view_type(wptr, workspace_extent - workspace_used);

      /// initialize problem
      problem._work = pw;
      problem._temperature= temperature;
      problem._pressure =pressure;
      problem._const_concentration= constYs;
      problem._photo_rates=photo_rates_at_i;
      problem._external_sources=external_sources_at_i;

      const ordinal_type r_val =
      time_integrator_type::invoke(member,
                                  problem,
                                  max_num_time_iterations,
                                  dt_in,
                                  t_beg,
                                  t_end,
                                  activeYs,
                                  t_out_at_i,
                                  dt_out_at_i,
                                  activeYs_out,
                                  tw);
    }
  }
      });
    Kokkos::Profiling::popRegion();
  }



#define TCHEM_RUN_ATMOSPHERIC_CHEMISTRY_E3SM()     \
  AtmosphericChemistryE3SM_ExplicitEuler_TemplateRun(          \
          profile_name,       \
          policy,       \
          tadv,         \
          state,        \
          photo_rates,  \
          external_sources, \
          t_out,        \
          dt_out,       \
          state_out,        \
          kmcd)

  void
  AtmosphericChemistryE3SM_ExplicitEuler::runHostBatch( /// thread block size
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
           const KineticModelNCAR_ConstData<interf_host_device_type>& kmcd)
  {
    const std::string profile_name = "TChem::AtmosphericChemistryE3SM_ExpliciEuler::runHostBatch::kmcd array";
    TCHEM_RUN_ATMOSPHERIC_CHEMISTRY_E3SM();
  }// namespace TChem
  void
  AtmosphericChemistryE3SM_ExplicitEuler::runDeviceBatch( /// thread block size
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
           const KineticModelNCAR_ConstData<device_type >& kmcd)
  {
    const std::string profile_name = "TChem::AtmosphericChemistryE3SM_ExpliciEuler::runDeviceBatch::kmcd array";
    TCHEM_RUN_ATMOSPHERIC_CHEMISTRY_E3SM();
  }

}
