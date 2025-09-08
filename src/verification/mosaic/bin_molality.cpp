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

void bin_molality(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto je_arr = input.get_array("je");
    const auto aH2O_a_arr = input.get_array("aH2O_a");

    real_type_1d_view je("je", 1);
    verification::convert_1d_vector_to_1d_view_device(je_arr, je);

    real_type_1d_view aH2O_a("aH2O_a", 1);
    verification::convert_1d_vector_to_1d_view_device(aH2O_a_arr, aH2O_a);

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    // Prepare variables for output

    // Reals or int that are defined outside of the parallel_for region are passed as const.
    real_type_1d_view outputs_bin_molality("outputs_bin_molality", 1);

    std::string profile_name ="Verification_test_bin_molality";
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

      Real& molality = outputs_bin_molality(0);

    // Perform the adjustment calculation
    TChem::Impl::MOSAIC<real_type, device_type>::bin_molality(
      mmd,
      je(0)-1,
      aH2O_a(0),
      molality);
    });

    const auto outputs_bin_molality_h = Kokkos::create_mirror_view_and_copy(host_exec_space, outputs_bin_molality);

    Real molality = outputs_bin_molality_h(0);

    // Assuming the outputs are scalar and can be directly set in the ensemble
    output.set("molality", molality);

  });
}
