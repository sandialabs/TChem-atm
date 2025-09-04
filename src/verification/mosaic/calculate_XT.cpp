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
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void calculate_XT(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {
    // Assuming the input arrays are 1D and of the same length

    auto aer_db = input.get_array("aer");

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    real_type_1d_view aer("aer", aer_db.size());

    // Convert input data from std::vector or similar structure to Kokkos views
    verification::convert_1d_vector_to_1d_view_device(aer_db, aer);

    // Prepare variables for output

    // Reals or int that are defined outside of the parallel_for region are passed as const.
    real_type_1d_view outputs_calculate_XT("outputs_calculate_XT", 1);

    std::string profile_name ="Verification_test_calculate_XT";
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

      Real& XT = outputs_calculate_XT(0);

    // Perform the adjustment calculation
    TChem::Impl::MOSAIC<real_type, device_type>::calculate_XT(
      mmd,
      aer,
      XT);
    });

    const auto outputs_calculate_XT_h = Kokkos::create_mirror_view_and_copy(host_exec_space, outputs_calculate_XT);

    Real XT = outputs_calculate_XT_h(0);

    // Assuming the outputs are scalar and can be directly set in the ensemble
    output.set("XT", XT);

    // Convert input data from Kokkos views to std::vector
    verification::convert_1d_view_device_to_1d_vector(aer, aer_db);

    output.set("aer", aer_db);
  });
}
