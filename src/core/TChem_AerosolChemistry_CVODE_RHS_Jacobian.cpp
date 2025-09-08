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
#include "TChem_AerosolChemistry_CVODE_RHS_Jacobian.hpp"

namespace TChem
{

// Right hand side function dy/dt = f(t,y)
int AerosolChemistry_CVODE_K::f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  auto udata = static_cast<UserData*>(user_data);

  const auto nbatches  = udata->nbatches;
  const auto policy = udata->policy;
  const auto kmcd = udata->kmcd;
  const auto amcd = udata->amcd;
  const auto batchSize = udata->batchSize;
  const auto num_concentration = udata->num_concentration;
  const auto temperature= udata->temperature;
  const auto pressure= udata->pressure;
  const auto const_tracers= udata->const_tracers;

  real_type_2d_view_type vals(N_VGetDeviceArrayPointer(y), nbatches, batchSize);
  real_type_2d_view_type rhs(N_VGetDeviceArrayPointer(ydot), nbatches, batchSize);

  const ordinal_type level = 1;
  // const ordinal_type number_of_equations = batchSize;
  const ordinal_type per_team_extent
         = TChem::Impl::Aerosol_RHS<real_type, device_type>::getWorkSpaceSize(kmcd, amcd);
  const std::string profile_name = "TChem::AerosolChemistry::RHS_evaluation";

  TChemAerosolChemistryRHS rhs_tchem(rhs, vals, num_concentration,
     const_tracers, temperature, pressure, kmcd, amcd);

  Kokkos::Profiling::pushRegion(profile_name);
  Kokkos::parallel_for(profile_name, policy, rhs_tchem);
  Kokkos::Profiling::popRegion();

  udata->team_size_recommended_rhs = policy.team_size_recommended(
    rhs_tchem, Kokkos::ParallelForTag());
  udata->team_size_max_rhs = policy.team_size_max(
    rhs_tchem, Kokkos::ParallelForTag());

 return 0;
}

// Jacobian of f(t,y)
int AerosolChemistry_CVODE_K::Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data
, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  auto udata = static_cast<UserData*>(user_data);

  const auto nbatches  = udata->nbatches;
  const auto policy = udata->policy;
  const auto kmcd = udata->kmcd;
  const auto amcd = udata->amcd;
  const auto num_concentration = udata->num_concentration;
  const auto temperature= udata->temperature;
  const auto pressure= udata->pressure;
  const auto const_tracers= udata->const_tracers;
  const auto batchSize = udata->batchSize;
  const auto fac= udata->fac;
#if defined(TCHEM_ATM_ENABLE_GPU)
  const auto JacRL= udata->JacRL;
#endif
  auto J_data = sundials::kokkos::GetDenseMat<MatType>(J)->View();

  real_type_2d_view_type y2d(N_VGetDeviceArrayPointer(y), nbatches, batchSize);
  real_type_2d_view_type vals(N_VGetDeviceArrayPointer(y), nbatches, batchSize);

  const ordinal_type level = 1;
  const ordinal_type number_of_equations = problem_type::getNumberOfTimeODEs(kmcd, amcd);
  const ordinal_type per_team_extent = problem_type::getWorkSpaceSize(kmcd,amcd)
    + number_of_equations;
  const std::string profile_name = "TChem::AerosolChemistry::Jacobian_evaluation";
  auto jac_lambda =        KOKKOS_LAMBDA(const typename policy_type::member_type& member) {

    const ordinal_type i = member.league_rank();

    const ordinal_type m = problem_type::getNumberOfTimeODEs(kmcd,amcd);
    const real_type_1d_view_type vals_at_i =
    Kokkos::subview(vals, i, Kokkos::ALL());

    const real_type_1d_view_type fac_at_i =
    Kokkos::subview(fac, i, Kokkos::ALL());

    const real_type_2d_view_type jacobian_at_i =
#if defined(TCHEM_ATM_ENABLE_GPU)
    Kokkos::subview(JacRL, i, Kokkos::ALL(), Kokkos::ALL());
#else
    Kokkos::subview(J_data, i, Kokkos::ALL(), Kokkos::ALL());
#endif

    const real_type_1d_view_type number_conc_at_i =
    Kokkos::subview(num_concentration, i, Kokkos::ALL());

    TChem::Scratch<real_type_1d_view_type> work(member.team_scratch(level),
                                   per_team_extent);

    auto wptr = work.data();

    const real_type_1d_view_type constYs  = Kokkos::subview(const_tracers, i, Kokkos::ALL());

    const ordinal_type problem_workspace_size = problem_type::getWorkSpaceSize(kmcd,amcd);
    auto pw = real_type_1d_view_type(wptr, problem_workspace_size);
    wptr +=problem_workspace_size;
    problem_type problem;
    problem._kmcd = kmcd;
    problem._amcd = amcd;
    /// initialize problem
    problem._fac = fac_at_i;
    problem._work = pw;
    problem._temperature= temperature(i);
    problem._pressure =pressure(i);
    problem._const_concentration= constYs;
    problem._number_conc =number_conc_at_i;
    problem.computeNumericalJacobian(member,vals_at_i,jacobian_at_i);
#if defined(TCHEM_ATM_ENABLE_GPU)
    //NOTE: Sundials uses a left layout on the device, while we are using a right layout.
    // This incompatible layouts is producing a runtime error.
    for (int k=0;k<batchSize;++k)
      for (int j=0;j<batchSize;++j)
       J_data(i, k, j) = jacobian_at_i(k,j);
#endif
    };
  Kokkos::Profiling::pushRegion(profile_name);
  Kokkos::parallel_for(profile_name, policy, jac_lambda);
  Kokkos::Profiling::popRegion();

  udata->team_size_recommended_jac = policy.team_size_recommended(
    jac_lambda, Kokkos::ParallelForTag());
  udata->team_size_max_jac = policy.team_size_max(
    jac_lambda, Kokkos::ParallelForTag());

  return 0;

}

}// namespace TChem
