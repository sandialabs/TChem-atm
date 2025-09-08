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

void adjust_solid_aerosol(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {
    // Assuming the input arrays are 1D and of the same length

    const auto aer_db = input.get_array("aer");
    const auto electrolyte_db = input.get_array("electrolyte");
    const auto epercent_db = input.get_array("epercent");

    // Assuming the size of arrays is known and uniform across all inputs
    const ordinal_type n = 3;
    const auto nsize_aero = static_cast<ordinal_type>(aer_db.size())/n;
    const auto nsize_electrolyte = static_cast<ordinal_type>(electrolyte_db.size())/n;
    const auto nsize_epercent = static_cast<ordinal_type>(epercent_db.size())/n;

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    real_type_1d_view aer_solid("aer_solid", nsize_aero);
    real_type_1d_view aer_liquid("aer_liquid", nsize_aero);
    real_type_1d_view aer_total("aer_total", nsize_aero);
    real_type_1d_view electrolyte_solid("electrolyte_solid", nsize_electrolyte);
    real_type_1d_view electrolyte_liquid("electrolyte_liquid", nsize_electrolyte);
    real_type_1d_view electrolyte_total("electrolyte_total", nsize_electrolyte);
    real_type_1d_view epercent_solid("epercent_solid", nsize_epercent);
    real_type_1d_view epercent_liquid("epercent_liquid", nsize_epercent);
    real_type_1d_view epercent_total("epercent_total", nsize_epercent);

    std::vector<std::vector<real_type>> aer_db_2d;
    for (size_t i = 0; i < n*nsize_aero; i += nsize_aero) {
        // Create a small vector from a segment of long_vector
        std::vector<real_type> aer_temp(aer_db.begin() + i, aer_db.begin() + i + nsize_aero);
        aer_db_2d.push_back(aer_temp);
    }

    // Convert input data from std::vector or similar structure to Kokkos views
    verification::convert_1d_vector_to_1d_view_device(aer_db_2d[0], aer_solid);
    verification::convert_1d_vector_to_1d_view_device(aer_db_2d[1], aer_liquid);
    verification::convert_1d_vector_to_1d_view_device(aer_db_2d[2], aer_total);

    std::vector<std::vector<real_type>> epercent_db_2d;
    for (size_t i = 0; i < n*nsize_epercent; i += nsize_epercent) {
        // Create a small vector from a segment of long_vector
        std::vector<real_type> epercent_temp(epercent_db.begin() + i, epercent_db.begin() + i + nsize_epercent);
        epercent_db_2d.push_back(epercent_temp);
    }
    // Convert input data from std::vector or similar structure to Kokkos views
    verification::convert_1d_vector_to_1d_view_device(epercent_db_2d[0], epercent_solid);
    verification::convert_1d_vector_to_1d_view_device(epercent_db_2d[1], epercent_liquid);
    verification::convert_1d_vector_to_1d_view_device(epercent_db_2d[2], epercent_total);


    std::vector<std::vector<real_type>> electrolyte_db_2d;
    for (size_t i = 0; i < n*nsize_electrolyte; i += nsize_electrolyte) {
        // Create a small vector from a segment of long_vector
        std::vector<real_type> electrolyte_temp(electrolyte_db.begin() + i, electrolyte_db.begin() + i + nsize_electrolyte);
        electrolyte_db_2d.push_back(electrolyte_temp);
    }

    // Convert input data from std::vector or similar structure to Kokkos views
    verification::convert_1d_vector_to_1d_view_device(electrolyte_db_2d[0], electrolyte_solid);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_db_2d[1], electrolyte_liquid);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_db_2d[2], electrolyte_total);

    // Prepare variables for output

    // Reals or int that are defined outside of the parallel_for region are passed as const.
    real_type_1d_view outputs_adjust_solid_aerosol("outputs_adjust_solid_aerosol", 3);

    std::string profile_name ="Verification_test_adjust_solid_aerosol";
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

      Real& water_a = outputs_adjust_solid_aerosol(0);
      Real& jphase = outputs_adjust_solid_aerosol(1);
      Real& jhyst_leg = outputs_adjust_solid_aerosol(2);

    // Perform the adjustment calculation
    TChem::Impl::MOSAIC<real_type, device_type>::adjust_solid_aerosol(
      mmd,
      aer_solid, aer_liquid, aer_total,
      electrolyte_solid, electrolyte_liquid, electrolyte_total,
      epercent_solid, epercent_liquid, epercent_total,
      water_a, jphase, jhyst_leg);
    });

    const auto outputs_adjust_solid_aerosol_h = Kokkos::create_mirror_view_and_copy(host_exec_space, outputs_adjust_solid_aerosol);

    Real water_a = outputs_adjust_solid_aerosol_h(0);
    Real jphase = outputs_adjust_solid_aerosol_h(1);
    Real jhyst_leg = outputs_adjust_solid_aerosol_h(2);


    // Assuming the outputs are scalar and can be directly set in the ensemble
    output.set("water_a", water_a);
    output.set("jphase", jphase);
    output.set("jhyst_leg", jhyst_leg);

    // Convert input data from Kokkos views to std::vector
    verification::convert_1d_view_device_to_1d_vector(aer_solid, aer_db_2d[0]);
    verification::convert_1d_view_device_to_1d_vector(aer_liquid, aer_db_2d[1]);
    verification::convert_1d_view_device_to_1d_vector(aer_total, aer_db_2d[2]);

    std::vector<real_type> aer_flattened;
    for (const auto& row : aer_db_2d) {
        aer_flattened.insert(aer_flattened.end(), row.begin(), row.end());
    }
    output.set("aer", aer_flattened);

    verification::convert_1d_view_device_to_1d_vector(epercent_solid, epercent_db_2d[0]);
    verification::convert_1d_view_device_to_1d_vector(epercent_liquid, epercent_db_2d[1]);
    verification::convert_1d_view_device_to_1d_vector(epercent_total, epercent_db_2d[2]);

    std::vector<real_type> epercent_flattened;
    for (const auto& row : epercent_db_2d) {
        epercent_flattened.insert(epercent_flattened.end(), row.begin(), row.end());
    }
    output.set("epercent", epercent_flattened);

    verification::convert_1d_view_device_to_1d_vector(electrolyte_solid, electrolyte_db_2d[0]);
    verification::convert_1d_view_device_to_1d_vector(electrolyte_liquid, electrolyte_db_2d[1]);
    verification::convert_1d_view_device_to_1d_vector(electrolyte_total, electrolyte_db_2d[2]);

    std::vector<real_type> electrolyte_flattened;
    for (const auto& row : electrolyte_db_2d) {
        electrolyte_flattened.insert(electrolyte_flattened.end(), row.begin(), row.end());
    }
    output.set("electrolyte", electrolyte_flattened);

  });
}
