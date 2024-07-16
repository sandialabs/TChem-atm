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
#include "TChem_Util.hpp"

#include "TChem_AtmosphericChemistryE3SM_CVODE.hpp"

namespace TChem
{
	  template<typename PolicyType,
           typename ValueType,
           typename DeviceType>
  void
  AtmosphericChemistryE3SM_CVODE_TemplateRun( /// required template arguments
          const std::string& profile_name,
          const ValueType& dummyValueType,
          /// team size setting
          const PolicyType& policy,
          /// input
          const Tines::value_type_2d_view<real_type, DeviceType>& tol,
          const Tines::value_type_2d_view<real_type, DeviceType>& fac,
          const Tines::value_type_1d_view<time_advance_type, DeviceType>& tadv,
          const Tines::value_type_2d_view<real_type, DeviceType>& state,
          const Tines::value_type_2d_view<real_type, DeviceType>& photo_rates,
          const Tines::value_type_2d_view<real_type, DeviceType>& external_sources,
          /// output
          const Tines::value_type_1d_view<real_type, DeviceType>& t_out,
          const Tines::value_type_1d_view<real_type, DeviceType>& dt_out,
          const Tines::value_type_2d_view<real_type, DeviceType>& state_out,
          /// const data from kinetic model
          const Tines::value_type_1d_view<KineticModelNCAR_ConstData<DeviceType>,DeviceType>& kmcds,
          const Tines::value_type_1d_view<Tines::TimeIntegratorCVODE<real_type,DeviceType>,DeviceType>& cvodes) {
    Kokkos::Profiling::pushRegion(profile_name);
#if defined(TINES_ENABLE_TPL_SUNDIALS)
    using policy_type = PolicyType;

    using real_type_0d_view_type = Tines::value_type_0d_view<real_type, DeviceType>;
    using real_type_1d_view_type = Tines::value_type_1d_view<real_type, DeviceType>;
    // const auto kmcd_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Kokkos::subview(kmcds, 0));

    using problem_type = TChem::Impl::AtmosphericChemistryE3SM_Problem<real_type,DeviceType>;
    const ordinal_type level = 1;
    const ordinal_type per_team_extent = AtmosphericChemistryE3SM_CVODE::getWorkSpaceSize(kmcds(0));

    /// this assumes host space
    Kokkos::parallel_for
      (profile_name,
       policy,
       KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
  const ordinal_type i = member.league_rank();
  const auto kmcd = (kmcds.extent(0) == 1 ? kmcds(0) : kmcds(i));
        auto cvode = cvodes(i);

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
      // Note: The number of external source can be equal to zero.
      real_type_1d_view_type external_sources_at_i;
      if(external_sources.extent(0) > 0)
      {
        external_sources_at_i= Kokkos::subview(external_sources, i, Kokkos::ALL());
      }
      const real_type_1d_view_type state_out_at_i =
        Kokkos::subview(state_out, i, Kokkos::ALL());
      const real_type_1d_view_type fac_at_i =
        Kokkos::subview(fac, i, Kokkos::ALL());

      const real_type_0d_view_type dt_out_at_i = Kokkos::subview(dt_out, i);
      Scratch<real_type_1d_view_type> work(member.team_scratch(level),
                                       per_team_extent);
      auto wptr = work.data();
      const real_type_1d_view_type ww(wptr, work.extent(0));

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
      const auto activeYs = real_type_1d_view_type(Ys.data(),
                              kmcd.nSpec - kmcd.nConstSpec );
      const auto constYs  = real_type_1d_view_type(Ys.data()
                            + kmcd.nSpec - kmcd.nConstSpec,  kmcd.nSpec );

      const real_type_0d_view_type temperature_out(sv_out_at_i.TemperaturePtr());
      const real_type_0d_view_type pressure_out(sv_out_at_i.PressurePtr());
      const real_type_1d_view_type Ys_out = sv_out_at_i.MassFractions();
      const real_type_0d_view_type density_out(sv_out_at_i.DensityPtr());

      const ordinal_type m = problem_type::getNumberOfTimeODEs(kmcd);
      auto vals = cvode.getStateVector();

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
      problem._fac = fac_at_i;
      problem._work_cvode = tw; // time integrator workspace
      problem._work = pw;
      problem._temperature= temperature;
      problem._pressure =pressure;
      problem._const_concentration= constYs;
      problem._photo_rates=photo_rates_at_i;
      problem._external_sources=external_sources_at_i;

      for (ordinal_type i=0;i<m;++i)
        vals(i) = activeYs(i);

      real_type t = t_out_at_i(), dt = 0;
      cvode.initialize(t,
           dt_in, dt_min, dt_max,
           tol(0,0), tol(0,1),
           TChem::Impl::ProblemAtmosphericChemistryE3SM_ComputeFunctionCVODE,
           TChem::Impl::ProblemAtmosphericChemistryE3SM_ComputeJacobianCVODE);

      cvode.setProblem(problem);
      for (ordinal_type iter=0;iter<max_num_time_iterations && t<=t_end;++iter) {
        const real_type t_prev = t;
        cvode.advance(t_end, t, -1);
        dt = t - t_prev;
      }
      t_out_at_i() = t;
      dt_out_at_i() = dt;

      for (ordinal_type i=0;i<m;++i)
        Ys_out(i) = vals(i);
      temperature_out() =temperature;

    }
  }
      });
#endif
    Kokkos::Profiling::popRegion();
  }



#define TCHEM_RUN_ATMOSPHERIC_CHEMISTRY_E3SM()     \
  AtmosphericChemistryE3SM_CVODE_TemplateRun(          \
          profile_name,       \
          value_type(),       \
          policy,       \
          tol,          \
          fac,          \
          tadv,         \
          state,        \
          photo_rates,  \
          external_sources, \
          t_out,        \
          dt_out,       \
          state_out,        \
              kmcds,        \
          cvodes)

  void
  AtmosphericChemistryE3SM_CVODE::runHostBatch( /// thread block size
           typename UseThisTeamPolicy<host_exec_space>::type& policy,
           /// input
           const real_type_2d_view_host& tol,
           const real_type_2d_view_host& fac,
           const time_advance_type_1d_view_host& tadv,
           const real_type_2d_view_host& state,
           const real_type_2d_view_host& photo_rates,
           const real_type_2d_view_host& external_sources,
           /// output
           const real_type_1d_view_host& t_out,
           const real_type_1d_view_host& dt_out,
           const real_type_2d_view_host& state_out,
           /// const data from kinetic model
           const Tines::value_type_1d_view<KineticModelNCAR_ConstData<interf_host_device_type>, interf_host_device_type>& kmcds,
           const Tines::value_type_1d_view<Tines::TimeIntegratorCVODE<real_type, host_device_type>, host_device_type>& cvodes) {
#if defined(TINES_ENABLE_TPL_SUNDIALS)
    const std::string profile_name = "TChem::AtmosphericChemistryE3SM::runHostBatch::kmcd array";
    // Note: we do not support SACADO and CVODE. Thus, CVODE uses a numerical Jacobian.
    using value_type = real_type;
    TCHEM_RUN_ATMOSPHERIC_CHEMISTRY_E3SM();
#endif
  }// namespace TChem
}
