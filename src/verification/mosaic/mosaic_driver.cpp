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
#include <iostream>
#include <verification.hpp>

#if defined(TCHEM_ATM_ENABLE_SKYWALKER)
 #include "skywalker.hpp"
 using namespace skywalker;
 using namespace TChem;
#endif
void usage() {
  std::cerr << "mosaic_driver: a Skywalker driver for validating the "
               "mosaic routines."
            << std::endl;
  std::cerr << "mosaic_driver: usage:" << std::endl;
  std::cerr << "mosaic_driver <input.yaml>" << std::endl;
  exit(0);
}

void check_aerosol_mass(Ensemble *ensemble);

void adjust_liquid_aerosol(Ensemble *ensemble);

void adjust_solid_aerosol(Ensemble *ensemble);

void do_full_deliquescence(Ensemble *ensemble);

void calculate_XT(Ensemble *ensemble);

void fnlog_gamZ(Ensemble *ensemble);

void fuchs_sutugin(Ensemble *ensemble);

void gas_diffusivity(Ensemble *ensemble);

void mean_molecular_speed(Ensemble *ensemble);

void fn_Keq(Ensemble *ensemble);

void drh_mutual(Ensemble *ensemble);

void fn_Po(Ensemble *ensemble);

void molality_0(Ensemble *ensemble);

void bin_molality(Ensemble *ensemble);

void bin_molality_60(Ensemble *ensemble);

void MTEM_compute_log_gamZ(Ensemble *ensemble);

void aerosol_water_up(Ensemble *ensemble);

int main(int argc, char **argv) {
  if (argc == 1) {
    usage();
  }

  verification::initialize(argc, argv);
  std::string input_file = argv[1];
  std::string output_file = verification::output_name(input_file);
  std::cout << argv[0] << ": reading " << input_file << std::endl;

  // Load the ensemble. Any error encountered is fatal.
  Ensemble *ensemble = skywalker::load_ensemble(input_file, "TChem-atm");

  // the settings.
  Settings settings = ensemble->settings();
  if (!settings.has("function")) {
    std::cerr << "No function specified in TChem-atm.settings!" << std::endl;
    exit(1);
  }

  // Dispatch to the requested function.
  auto func_name = settings.get("function");
  try {
    if (func_name == "check_aerosol_mass") {
      check_aerosol_mass(ensemble);
    } else if (func_name == "adjust_liquid_aerosol") {
      adjust_liquid_aerosol(ensemble);
    } else if (func_name == "adjust_solid_aerosol") {
      adjust_solid_aerosol(ensemble);
    } else if (func_name == "do_full_deliquescence") {
      do_full_deliquescence(ensemble);
    } else if (func_name == "calculate_XT") {
      calculate_XT(ensemble);
    } else if (func_name == "fnlog_gamZ") {
      fnlog_gamZ(ensemble);
    } else if (func_name == "fuchs_sutugin") {
      fuchs_sutugin(ensemble);
    } else if (func_name == "gas_diffusivity") {
      gas_diffusivity(ensemble);
    } else if (func_name == "mean_molecular_speed") {
      mean_molecular_speed(ensemble);
    } else if (func_name == "fn_Keq") {
      fn_Keq(ensemble);
    } else if (func_name == "drh_mutual") {
      drh_mutual(ensemble);
    } else if (func_name == "fn_Po") {
      fn_Po(ensemble);
    } else if (func_name == "molality_0") {
      molality_0(ensemble);
    } else if (func_name == "bin_molality") {
      bin_molality(ensemble);
    } else if (func_name == "bin_molality_60") {
      bin_molality_60(ensemble);
    } else if (func_name == "MTEM_compute_log_gamZ") {
      MTEM_compute_log_gamZ(ensemble);
    } else if (func_name == "aerosol_water_up") {
      aerosol_water_up(ensemble);
    } else {
      std::cerr << "Error: Function name '" << func_name
                << "' does not have an implemented test!" << std::endl;
      exit(1);
    }
  } catch (std::exception &e) {
    std::cerr << argv[0] << ": Error: " << e.what() << std::endl;
  }

  // Write out a Python module.
  std::cout << argv[0] << ": writing " << output_file << std::endl;
  ensemble->write(output_file);

  // Clean up.
  delete ensemble;
  verification::finalize();
}
