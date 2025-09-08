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

void aerosol_water_up(Ensemble *ensemble) {
    ensemble->process([=](const Input &input, Output &output) {

        const auto electrolyte_db = input.get_array("electrolyte");

        const ordinal_type n = 3;
        const auto nsize_electrolyte = static_cast<ordinal_type>(electrolyte_db.size())/n;

        const auto mmd = TChem::Impl::MosaicModelData<device_type>();

        real_type_1d_view electrolyte_solid("electrolyte_solid", nsize_electrolyte);
        real_type_1d_view electrolyte_liquid("electrolyte_liquid", nsize_electrolyte);
        real_type_1d_view electrolyte_total("electrolyte_total", nsize_electrolyte);

        std::vector<std::vector<real_type>> electrolyte_db_2d;
        for (size_t i = 0; i < n*nsize_electrolyte; i += nsize_electrolyte) {
            std::vector<real_type> electrolyte_temp(electrolyte_db.begin() + i, electrolyte_db.begin() + i + nsize_electrolyte);
            electrolyte_db_2d.push_back(electrolyte_temp);
        }

        verification::convert_1d_vector_to_1d_view_device(electrolyte_db_2d[0], electrolyte_solid);
        verification::convert_1d_vector_to_1d_view_device(electrolyte_db_2d[1], electrolyte_liquid);
        verification::convert_1d_vector_to_1d_view_device(electrolyte_db_2d[2], electrolyte_total);

        real_type_1d_view outputs_aerosol_water_up("outputs_aerosol_water_up", 1);

        std::string profile_name ="Verification_test_aerosol_water_up";
        using policy_type =
            typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
        const auto exec_space_instance = TChem::exec_space();
        const auto host_exec_space = TChem::host_exec_space();
        policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

        Kokkos::parallel_for(
        profile_name,
        policy,
        KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
          Real& aerosol_water = outputs_aerosol_water_up(0);

        TChem::Impl::MOSAIC<real_type, device_type>::aerosol_water_up(
          mmd,
          electrolyte_total,
          aerosol_water);
        });

        const auto outputs_aerosol_water_up_h = Kokkos::create_mirror_view_and_copy(host_exec_space, outputs_aerosol_water_up);

        Real aerosol_water = outputs_aerosol_water_up_h(0);

        output.set("aerosol_water", aerosol_water);
    });
}
