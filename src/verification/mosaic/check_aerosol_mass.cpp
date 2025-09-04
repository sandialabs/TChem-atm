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

void check_aerosol_mass(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto aer_db = input.get_array("aer");
    auto mass_dry_a_arr = input.get_array("mass_dry_a");
    auto jphase_arr = input.get_array("jphase");
    auto jaerosolstate_arr = input.get_array("jaerosolstate");
    auto num_a_arr = input.get_array("num_a");

    real_type_1d_view mass_dry_a("mass_dry_a", 1);
    verification::convert_1d_vector_to_1d_view_device(mass_dry_a_arr, mass_dry_a);

    real_type_1d_view jphase("jphase", 1);
    verification::convert_1d_vector_to_1d_view_device(jphase_arr, jphase);

    real_type_1d_view jaerosolstate("jaerosolstate", 1);
    verification::convert_1d_vector_to_1d_view_device(jaerosolstate_arr, jaerosolstate);

    real_type_1d_view num_a("num_a", 1);
    verification::convert_1d_vector_to_1d_view_device(num_a_arr, num_a);

    const ordinal_type n = 3;
    const auto nsize_aero = static_cast<ordinal_type>(aer_db.size())/n;

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    real_type_1d_view aer_solid("aer_solid", nsize_aero);
    real_type_1d_view aer_liquid("aer_liquid", nsize_aero);
    real_type_1d_view aer_total("aer_total", nsize_aero);

    std::vector<std::vector<real_type>> aer_db_2d;
    for (size_t i = 0; i < n*nsize_aero; i += nsize_aero) {
        // Create a small vector from a segment of long_vector
        std::vector<real_type> aer_temp(aer_db.begin() + i, aer_db.begin() + i + nsize_aero);
        aer_db_2d.push_back(aer_temp);
    }

    verification::convert_1d_vector_to_1d_view_device(aer_db_2d[0], aer_solid);
    verification::convert_1d_vector_to_1d_view_device(aer_db_2d[1], aer_liquid);
    verification::convert_1d_vector_to_1d_view_device(aer_db_2d[2], aer_total);

    std::string profile_name ="Verification_test_check_aerosol_mass";
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
    // Perform the adjustment calculation
    TChem::Impl::MOSAIC<real_type, device_type>::check_aerosol_mass(
      mmd,
      aer_total,
      mass_dry_a(0),
      jphase(0),
      jaerosolstate(0),
      num_a(0));
    });

    // Convert modified variables data from Kokkos views to std::vector for output
    verification::convert_1d_view_device_to_1d_vector(aer_solid, aer_db_2d[0]);
    verification::convert_1d_view_device_to_1d_vector(aer_liquid, aer_db_2d[1]);
    verification::convert_1d_view_device_to_1d_vector(aer_total, aer_db_2d[2]);

    std::vector<real_type> aer_flattened;
    for (const auto& row : aer_db_2d) {
        aer_flattened.insert(aer_flattened.end(), row.begin(), row.end());
    }
    output.set("aer", aer_flattened);

    verification::convert_1d_view_device_to_1d_vector(mass_dry_a, mass_dry_a_arr);
    verification::convert_1d_view_device_to_1d_vector(jphase, jphase_arr);
    verification::convert_1d_view_device_to_1d_vector(jaerosolstate, jaerosolstate_arr);
    verification::convert_1d_view_device_to_1d_vector(num_a, num_a_arr);

    output.set("mass_dry_a", mass_dry_a_arr);
    output.set("jphase", jphase_arr);
    output.set("jaerosolstate", jaerosolstate_arr);
    output.set("num_a", num_a_arr);
  });
}
