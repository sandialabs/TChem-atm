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
#include "TChem_Util.hpp"
#include "TChem_Impl_RateofProgress.hpp"
#include "TChem_RateofProgress.hpp"

namespace TChem {

  template<typename PolicyType,
           typename DeviceType>
void
RateofProgress_TemplateRun( /// input
  const std::string& profile_name,
  /// team size setting
  const PolicyType& policy,

  const Tines::value_type_2d_view<real_type, DeviceType>& state,
  const Tines::value_type_2d_view<real_type, DeviceType>& photo_rates,
  /// output
  const Tines::value_type_2d_view<real_type, DeviceType>& rate_of_progress,
  /// const data from kinetic model
  const KineticModelNCAR_ConstData<DeviceType >& kmcd)
{
  Kokkos::Profiling::pushRegion(profile_name);

  using policy_type = PolicyType;
  using device_type = DeviceType;

  using real_type_1d_view_type = Tines::value_type_1d_view<real_type, device_type>;

  const ordinal_type level = 1;
  const ordinal_type per_team_extent =
  TChem::RateofProgress::getWorkSpaceSize(kmcd);

  Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
      const ordinal_type i = member.league_rank();
      const real_type_1d_view_type state_at_i =
        Kokkos::subview(state, i, Kokkos::ALL());
      const real_type_1d_view_type rate_of_progress_at_i =
        Kokkos::subview(rate_of_progress, i, Kokkos::ALL());

      // Note: The number of photo reactions can be equal to zero.
      real_type_1d_view_type photo_rate_at_i;
      if(photo_rates.extent(0) > 0)
      {
        photo_rate_at_i =
        Kokkos::subview(photo_rates, i, Kokkos::ALL());
      }
      Scratch<real_type_1d_view_type> work(member.team_scratch(level),
                                      per_team_extent);
      const Impl::StateVector<real_type_1d_view_type> sv_at_i(kmcd.nSpec,
                                                         state_at_i);
      //
      auto wptr = work.data();
      const real_type_1d_view_type ww(wptr, work.extent(0));
      TCHEM_CHECK_ERROR(!sv_at_i.isValid(),
                        "Error: input state vector is not valid");
      {
        const real_type t = sv_at_i.Temperature();
        const real_type p = sv_at_i.Pressure();
        const real_type_1d_view_type Ys = sv_at_i.MassFractions();

        // const real_type density = sv_at_i.Density();
        Impl::RateofProgress<real_type, device_type>
        ::team_invoke(member, t, p, Ys, photo_rate_at_i, rate_of_progress_at_i, ww, kmcd);

        member.team_barrier();

        }
    });
  Kokkos::Profiling::popRegion();
}
void
RateofProgress::runHostBatch( /// input
  const real_type_2d_view_host_type& state,
  const real_type_2d_view_host_type& photo_rates,
  /// output
  const real_type_2d_view_host_type& rate_of_progress,
  /// const data from kinetic model
  const kinetic_model_host_type& kmcd)
{

  using policy_type = Kokkos::TeamPolicy<host_exec_space>;

  const ordinal_type level = 1;
  const ordinal_type per_team_extent = getWorkSpaceSize(kmcd);
  const ordinal_type per_team_scratch =
    Scratch<real_type_1d_view_host_type>::shmem_size(per_team_extent);

  policy_type policy(state.extent(0), Kokkos::AUTO()); // fine
  policy.set_scratch_size(level, Kokkos::PerTeam(per_team_scratch));

  RateofProgress_TemplateRun( /// input
    "TChem::RateofProgress::runHostBatch",
    policy,
    state,
    photo_rates,
    rate_of_progress,
    kmcd);

}

void
RateofProgress::runDeviceBatch( /// input
  typename UseThisTeamPolicy<exec_space>::type& policy,
  const real_type_2d_view_type& state,
  const real_type_2d_view_type& photo_rates,
  /// output
  const real_type_2d_view_type& rate_of_progress,
  /// const data from kinetic model
  const kinetic_model_type& kmcd)
{

  RateofProgress_TemplateRun( /// input
    "TChem::RateofProgress::runDeviceBatch",
    policy,
    state,
    photo_rates,
    rate_of_progress,
    kmcd);

}

void
RateofProgress::runDeviceBatch( /// input
  const real_type_2d_view_type& state,
  const real_type_2d_view_type& photo_rates,
  /// output
  const real_type_2d_view_type& rate_of_progress,
  /// const data from kinetic model
  const kinetic_model_type& kmcd)
{
  Kokkos::Profiling::pushRegion("TChem::RateofProgress::runDeviceBatch");
  using policy_type = Kokkos::TeamPolicy<exec_space>;

  const ordinal_type level = 1;
  const ordinal_type per_team_extent = getWorkSpaceSize(kmcd);
  const ordinal_type per_team_scratch =
    Scratch<real_type_1d_view_type>::shmem_size(per_team_extent);

  policy_type policy(state.extent(0), Kokkos::AUTO()); // fine
  policy.set_scratch_size(level, Kokkos::PerTeam(per_team_scratch));

  RateofProgress_TemplateRun( /// input
    "TChem::RateofProgress::runDeviceBatch",
    policy,
    state,
    photo_rates,
    rate_of_progress,
    kmcd);

}

} // namespace TChem
