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
#include "TChem_CommandLineParser.hpp"
#include "TChem_Linv3StratosphereSolver.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using ordinal_type = TChem::ordinal_type;
using real_type = TChem::real_type;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_2d_view = TChem::real_type_2d_view;
using real_type_2d_view_host = TChem::real_type_2d_view_host;
using real_type_1d_view_host = TChem::real_type_1d_view_host;
using linoz_input_parameters_1d_view = TChem::linoz_input_parameters_1d_view;
using linoz_input_parameters_0d_view = TChem::linoz_input_parameters_0d_view;
using linoz_input_parameters_1d_view_host = TChem::linoz_input_parameters_1d_view_host;
using linoz_input_parameters_0d_view_host = TChem::linoz_input_parameters_0d_view_host;
using linoz_vmr_idx_type = TChem::linoz_vmr_idx_type;
using linoz_input_parameters_type = TChem::linoz_input_parameters_type;
using ordinal_type_1d_view_host = TChem::ordinal_type_1d_view_host;
using ordinal_type_1d_view = TChem::ordinal_type_1d_view;

using namespace skywalker;
using namespace TChem::verification;


void linv3_stratosphere_solver_test( const int team_size,
                                    const int vector_size, Ensemble *ensemble) {
  ensemble->process([=]( const Input &input, Output &output) {
    const Real zero = 0;

    const real_type psc_T = input.get_array("psc_T")[0];
    const real_type dt = input.get_array("delta_t")[0];
    const real_type lats = input.get_array("rlats")[0];
    const real_type chlorine_loading = input.get_array("chlorine_loading")[0];
    const real_type sza = input.get_array("sza")[0];

    const auto temp_db = input.get_array("temp");
    const int pver = temp_db.size();
    auto temp_host = real_type_1d_view_host((Real *)temp_db.data(), temp_db.size());
    const auto temperature = real_type_1d_view("temperature", pver);
    Kokkos::deep_copy(temperature, temp_host);

    const auto pmid_db = input.get_array("pmid");
    auto pmid_host = real_type_1d_view_host((Real *)pmid_db.data(), pmid_db.size());
    const auto pressure = real_type_1d_view("pressure", pver);
    Kokkos::deep_copy(pressure, pmid_host);

    const auto o3col_db = input.get_array("o3col");
    auto o3col_host = real_type_1d_view_host((Real *)o3col_db.data(), o3col_db.size());
    const auto o3col = real_type_1d_view("o3col", pver);
    Kokkos::deep_copy(o3col, o3col_host);

    // const real_type  = input.get_array("")[0];

    // output.set("scavratenum", std::vector(1, scavratenum));
    // output.set("scavratevol", std::vector(1, scavratevol));

    const auto xvmr_db = input.get_array("xvmr");
    const int gas_pcnst = xvmr_db.size()/pver;
    real_type_2d_view volume_mixing_ratio("volume_mixing_ratio", pver, gas_pcnst);

    convert_1d_vector_to_2d_view_device(xvmr_db,volume_mixing_ratio );

    std::vector<std::string> linoz_list = {
   "linoz_o3_clim",
   "linoz_n2o_clim",
   "linoz_noy_clim",
   "linoz_ch4_clim",
   "linoz_h2o_clim",
   "linoz_t_clim",
   "linoz_o3col_clim",
   "linoz_PmL_clim_no3",
   "linoz_dPmL_dO3X",
   "linoz_dPmL_dn2o_no3",
   "linoz_dPmL_dnoy_no3",
   "linoz_dPmL_dch4_no3",
   "linoz_dPmL_dh2o_no3",
   "linoz_dPmL_dT_no3",
   "linoz_dPmL_dO3col_no3",
   "linoz_PmL_clim_pn2o",
   "linoz_dPmL_dO3_pn2o",
   "linoz_dPmL_dn2o_pn2o",
   "linoz_dPmL_dnoy_pn2o",
   "linoz_dPmL_dch4_pn2o",
   "linoz_dPmL_dh2o_pn2o",
   "linoz_dPmL_dT_pn2o",
   "linoz_dPmL_dO3col_pn2o",
   "linoz_PmL_clim_ln2o",
   "linoz_dPmL_dO3_ln2o",
   "linoz_dPmL_dn2o_ln2o",
   "linoz_dPmL_dnoy_ln2o",
   "linoz_dPmL_dch4_ln2o",
   "linoz_dPmL_dh2o_ln2o",
   "linoz_dPmL_dT_ln2o",
   "linoz_dPmL_dO3col_ln2o",
   "linoz_PmL_clim_pnoy",
   "linoz_dPmL_dO3_pnoy",
   "linoz_dPmL_dn2o_pnoy",
   "linoz_dPmL_dnoy_pnoy",
   "linoz_dPmL_dch4_pnoy",
   "linoz_dPmL_dh2o_pnoy",
   "linoz_dPmL_dT_pnoy",
   "linoz_dPmL_dO3col_pnoy",
   "linoz_PmL_clim_lnoy",
   "linoz_dPmL_dO3_lnoy",
   "linoz_dPmL_dn2o_lnoy",
   "linoz_dPmL_dnoy_lnoy",
   "linoz_dPmL_dch4_lnoy",
   "linoz_dPmL_dh2o_lnoy",
   "linoz_dPmL_dT_lnoy",
   "linoz_dPmL_dO3col_lnoy",
   "linoz_PmL_clim_nch4",
   "linoz_dPmL_dO3_nch4",
   "linoz_dPmL_dn2o_nch4",
   "linoz_dPmL_dnoy_nch4",
   "linoz_dPmL_dch4_nch4",
   "linoz_dPmL_dh2o_nch4",
   "linoz_dPmL_dT_nch4",
   "linoz_dPmL_dO3col_nch4",
   "linoz_cariolle_psc",
   };
  const int n_linoz_vars = linoz_list.size();

      TChem::linoz_input_parameters_1d_dual_view linoz_inputs_dual("linoz inputs", pver);
      auto linoz_inputs_host = linoz_inputs_dual.view_host();

    {
      std::vector<std::vector<real_type>> db;
      for (int ivar = 0; ivar < n_linoz_vars; ++ivar)
      {
      const auto var_name = linoz_list[ivar];
      db.push_back(input.get_array(var_name));
      } // var

      for (int kk = 0; kk < pver; ++kk)
      {
      auto &linix_at_kk =  linoz_inputs_host(kk);
      for (int ivar = 0; ivar < n_linoz_vars; ++ivar)
      {
        if (linoz_list[ivar] =="linoz_o3_clim" )
        {
          linix_at_kk.linoz_o3_clim=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_n2o_clim" ) {
          linix_at_kk.linoz_n2o_clim=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_noy_clim" ) {
          linix_at_kk.linoz_noy_clim=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_ch4_clim" ) {
          linix_at_kk.linoz_ch4_clim=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_h2o_clim" ) {
          linix_at_kk.linoz_h2o_clim=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_t_clim" ) {
          linix_at_kk.linoz_t_clim=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_o3col_clim" ) {
          linix_at_kk.linoz_o3col_clim=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_PmL_clim_no3" ) {
          linix_at_kk.linoz_PmL_clim_no3=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dO3X" ) {
          linix_at_kk.linoz_dPmL_dO3X=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dn2o_no3" ) {
          linix_at_kk.linoz_dPmL_dn2o_no3=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dnoy_no3" ) {
          linix_at_kk.linoz_dPmL_dnoy_no3=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dch4_no3" ) {
          linix_at_kk.linoz_dPmL_dch4_no3=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dh2o_no3" ) {
          linix_at_kk.linoz_dPmL_dh2o_no3=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dT_no3" ) {
          linix_at_kk.linoz_dPmL_dT_no3=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dO3col_no3" ) {
          linix_at_kk.linoz_dPmL_dO3col_no3=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_PmL_clim_pn2o" ) {
          linix_at_kk.linoz_PmL_clim_pn2o=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dO3_pn2o" ) {
          linix_at_kk.linoz_dPmL_dO3_pn2o=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dn2o_pn2o" ) {
          linix_at_kk.linoz_dPmL_dn2o_pn2o=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dnoy_pn2o" ) {
          linix_at_kk.linoz_dPmL_dnoy_pn2o=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dch4_pn2o" ) {
          linix_at_kk.linoz_dPmL_dch4_pn2o=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dh2o_pn2o" ) {
          linix_at_kk.linoz_dPmL_dh2o_pn2o=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dT_pn2o" ) {
          linix_at_kk.linoz_dPmL_dT_pn2o=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dO3col_pn2o" ) {
          linix_at_kk.linoz_dPmL_dO3col_pn2o=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_PmL_clim_ln2o" ) {
          linix_at_kk.linoz_PmL_clim_ln2o=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dO3_ln2o" ) {
          linix_at_kk.linoz_dPmL_dO3_ln2o=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dn2o_ln2o" ) {
          linix_at_kk.linoz_dPmL_dn2o_ln2o=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dnoy_ln2o" ) {
          linix_at_kk.linoz_dPmL_dnoy_ln2o=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dch4_ln2o" ) {
          linix_at_kk.linoz_dPmL_dch4_ln2o=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dh2o_ln2o" ) {
          linix_at_kk.linoz_dPmL_dh2o_ln2o=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dT_ln2o" ) {
          linix_at_kk.linoz_dPmL_dT_ln2o=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dO3col_ln2o" ) {
          linix_at_kk.linoz_dPmL_dO3col_ln2o=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_PmL_clim_pnoy" ) {
          linix_at_kk.linoz_PmL_clim_pnoy=db[ivar][kk];
        }  else if (linoz_list[ivar] =="linoz_dPmL_dO3_pnoy" ) {
          linix_at_kk.linoz_dPmL_dO3_pnoy=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dn2o_pnoy" ) {
          linix_at_kk.linoz_dPmL_dn2o_pnoy=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dnoy_pnoy" ) {
          linix_at_kk.linoz_dPmL_dnoy_pnoy=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dch4_pnoy" ) {
          linix_at_kk.linoz_dPmL_dch4_pnoy=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dh2o_pnoy" ) {
          linix_at_kk.linoz_dPmL_dh2o_pnoy=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dT_pnoy" ) {
          linix_at_kk.linoz_dPmL_dT_pnoy=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dO3col_pnoy" ) {
          linix_at_kk.linoz_dPmL_dO3col_pnoy=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_PmL_clim_lnoy" ) {
          linix_at_kk.linoz_PmL_clim_lnoy=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dO3_lnoy" ) {
          linix_at_kk.linoz_dPmL_dO3_lnoy=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dn2o_lnoy" ) {
          linix_at_kk.linoz_dPmL_dn2o_lnoy=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dnoy_lnoy" ) {
          linix_at_kk.linoz_dPmL_dnoy_lnoy=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dch4_lnoy" ) {
          linix_at_kk.linoz_dPmL_dch4_lnoy=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dh2o_lnoy" ) {
          linix_at_kk.linoz_dPmL_dh2o_lnoy=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dT_lnoy" ) {
          linix_at_kk.linoz_dPmL_dT_lnoy=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dO3col_lnoy" ) {
          linix_at_kk.linoz_dPmL_dO3col_lnoy=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_PmL_clim_nch4" ) {
          linix_at_kk.linoz_PmL_clim_nch4=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dO3_nch4" ) {
          linix_at_kk.linoz_dPmL_dO3_nch4=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dn2o_nch4" ) {
          linix_at_kk.linoz_dPmL_dn2o_nch4=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dnoy_nch4" ) {
          linix_at_kk.linoz_dPmL_dnoy_nch4=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dch4_nch4" ) {
          linix_at_kk.linoz_dPmL_dch4_nch4=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dh2o_nch4" ) {
          linix_at_kk.linoz_dPmL_dh2o_nch4=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dT_nch4" ) {
          linix_at_kk.linoz_dPmL_dT_nch4=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_dPmL_dO3col_nch4" ) {
          linix_at_kk.linoz_dPmL_dO3col_nch4=db[ivar][kk];
        } else if (linoz_list[ivar] =="linoz_cariolle_psc" ) {
          linix_at_kk.linoz_cariolle_psc=db[ivar][kk];
        } else {
          printf("var does not exits %s \n",linoz_list[ivar].c_str() );
          exit(1);
        }

      }
      }// kk

    }
    linoz_inputs_dual.modify_host();
    linoz_inputs_dual.sync_device();
    const auto& linoz_inputs=linoz_inputs_dual.view_device();

    linoz_vmr_idx_type linoz_vmr_idx{};
    // Fortran to C++ indexing
    linoz_vmr_idx.o3_ndx = ordinal_type(input.get_array("o3_ndx")[0]) -1;
    linoz_vmr_idx.n2olnz_ndx = ordinal_type(input.get_array("n2olnz_ndx")[0]) -1;
    linoz_vmr_idx.noylnz_ndx = ordinal_type(input.get_array("noylnz_ndx")[0]) -1;
    linoz_vmr_idx.ch4lnz_ndx = ordinal_type(input.get_array("ch4lnz_ndx")[0]) -1;
    linoz_vmr_idx.h2olnz_ndx = ordinal_type(input.get_array("h2olnz_ndx")[0]) -1;
    linoz_vmr_idx.o3lnz_ndx= ordinal_type(input.get_array("o3lnz_ndx")[0]) -1;
    linoz_vmr_idx.ch4_ndx= ordinal_type(input.get_array("ch4_ndx")[0]) -1;
    linoz_vmr_idx.n2o_ndx= ordinal_type(input.get_array("n2o_ndx")[0]) -1;
    linoz_vmr_idx.no_ndx= ordinal_type(input.get_array("no_ndx")[0]) -1;
    linoz_vmr_idx.no2_ndx= ordinal_type(input.get_array("no2_ndx")[0]) -1;
    linoz_vmr_idx.hno3_ndx= ordinal_type(input.get_array("hno3_ndx")[0]) -1;


    const auto tropFlag_db = input.get_array("tropFlag");


    ordinal_type_1d_view tropFlag("tropFlag",pver);
    auto tropFlag_host = Kokkos::create_mirror_view(tropFlag);

    for (int i = 0; i < pver; ++i)
    {
      tropFlag_host(i) = ordinal_type(tropFlag_db[i]);
    }// end for i
    Kokkos::deep_copy(tropFlag, tropFlag_host);

    const auto h2ovmr_db = input.get_array("h2ovmr");
    auto h2ovm_host = real_type_1d_view_host((Real *)h2ovmr_db.data(), h2ovmr_db.size());
    const auto water_vapor_volume_mixing_ratio =real_type_1d_view("water_vapor_volume_mixing_ratio",pver);
    Kokkos::deep_copy(water_vapor_volume_mixing_ratio, h2ovm_host);

    const auto exec_space_instance = TChem::exec_space();

    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;

    /// team policy
    const int nBatch = pver;
    policy_type policy(exec_space_instance, nBatch, Kokkos::AUTO());

    if (team_size > 0 && vector_size > 0) {
        policy = policy_type(exec_space_instance,  nBatch, team_size, vector_size);
    } else if (team_size > 0 && vector_size < 0) {
      // only set team size
       policy = policy_type(exec_space_instance, nBatch,  team_size);
    }

    TChem::Linv3StratosphereSolver::runDeviceBatch( /// input
        policy,
        temperature, pressure, volume_mixing_ratio, dt, lats, psc_T, sza,
        chlorine_loading, o3col, tropFlag, water_vapor_volume_mixing_ratio,
        linoz_inputs, linoz_vmr_idx, volume_mixing_ratio);

    std::vector<Real> volume_mixing_ratio_out(pver*gas_pcnst, zero);
    convert_2d_view_device_to_1d_vector(volume_mixing_ratio,
                                                          volume_mixing_ratio_out);

    output.set("xvmr", volume_mixing_ratio_out);


  });
}


int main(int argc, char *argv[]) {

  int team_size(-1), vector_size(-1);
  std::string input_file("None");
  std::string output_file("None");
  /// parse command line arguments
  TChem::CommandLineParser opts(
      "validation test for linv3 stratosphere solver using Skywalker");
  opts.set_option<int>("team_thread_size", "time thread size ", &team_size);
  // opts.set_option<int>("batch_size", " number of batches or samples, e.g. 10  ", &nBatch);
  opts.set_option<int>("vector_thread_size", "vector thread size ",
                       &vector_size);
  opts.set_option<std::string>("input_file", "Chem file name e.g., chem.inp",
                               &input_file);
  opts.set_option<std::string>(
      "output_file", "Output file name e.g., IgnSolution.dat", &output_file);


  const bool r_parse = opts.parse(argc, argv);
  if (r_parse)
    return 0; // print help return
  Kokkos::initialize(argc, argv);
  {

    if (output_file =="None"){
     output_file = output_name(input_file);
    }
    std::cout << argv[0] << ": reading " << input_file << std::endl;

    Ensemble *ensemble = skywalker::load_ensemble(input_file, "mam4xx");

    // the settings.
    Settings settings = ensemble->settings();
    if (!settings.has("function")) {
      std::cerr << "No function specified in mam4xx.settings!" << std::endl;
      exit(1);
    }

    auto func_name = settings.get("function");

    linv3_stratosphere_solver_test(team_size, vector_size, ensemble);
    std::cout << argv[0] << ": writing " << output_file << std::endl;
    ensemble->write(output_file);

    // using Linv3StratosphereSolver
  }
  Kokkos::finalize();

  return 0;
}
