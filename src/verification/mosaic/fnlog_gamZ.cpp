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
#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using oridinal_type_1d_view = TChem::ordinal_type_1d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void fnlog_gamZ(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto jA_arr = input.get_array("jA");
    const auto jE_arr = input.get_array("jE");
    const auto aH2O_arr = input.get_array("aH2O");

    real_type_1d_view jA("jA", 1);
    verification::convert_1d_vector_to_1d_view_device(jA_arr, jA);

    real_type_1d_view jE("jE", 1);
    verification::convert_1d_vector_to_1d_view_device(jE_arr, jE);

    real_type_1d_view aH2O("aH20", 1);
    verification::convert_1d_vector_to_1d_view_device(aH2O_arr, aH2O);

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    // Prepare variables for output

    // Reals or int that are defined outside of the parallel_for region are passed as const.
    real_type_1d_view outputs_fnlog_gamZ("outputs_fnlog_gamZ", 1);

    std::string profile_name ="Verification_test_fnlog_gamZ";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    const auto host_exec_space = TChem::host_exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    // Check this routines on GPUs.
     Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {

      Real& log_gamZ = outputs_fnlog_gamZ(0);

    // Perform the adjustment calculation
    TChem::Impl::MOSAIC<real_type, device_type>::fnlog_gamZ(
      mmd,
      jA(0)-1,
      jE(0)-1,
      aH2O(0),
      log_gamZ);
    });

    const auto outputs_fnlog_gamZ_h = Kokkos::create_mirror_view_and_copy(host_exec_space, outputs_fnlog_gamZ);

    Real log_gamZ = outputs_fnlog_gamZ_h(0);

    // Assuming the outputs are scalar and can be directly set in the ensemble
    output.set("log_gamZ", log_gamZ);

  });
}
