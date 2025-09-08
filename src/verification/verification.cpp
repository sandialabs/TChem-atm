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
#include "verification.hpp"

namespace TChem {
namespace verification {
void initialize(int argc, char **argv) { Kokkos::initialize(argc, argv); }

void finalize() {
  Kokkos::finalize();
  }

std::string output_name(const std::string &input_file) {
  std::string output_file;
  size_t slash = input_file.find_last_of('/');
  size_t dot = input_file.find_last_of('.');
  if ((dot == std::string::npos) and (slash == std::string::npos)) {
    dot = input_file.length();
  }
  if (slash == std::string::npos) {
    slash = 0;
  } else {
    slash += 1;
    dot -= slash;
  }
  return std::string("TChem_") + input_file.substr(slash, dot) +
         std::string(".py");
}

void convert_1d_vector_to_2d_view_device(const std::vector<real_type> &var_std,
                                         const real_type_2d_view &var_device) {
  auto host = Kokkos::create_mirror_view(var_device);
  int count = 0;
  for (int d2 = 0; d2 < var_device.extent(1); ++d2) {
    for (int d1 = 0; d1 < var_device.extent(0); ++d1) {
      host(d1, d2) = var_std[count];
      count++;
    }
  }
  Kokkos::deep_copy(var_device, host);
}

void convert_2d_view_device_to_1d_vector(const real_type_2d_view &var_device,
                                         std::vector<real_type> &var_std)
{
  auto host = Kokkos::create_mirror_view(var_device);
  Kokkos::deep_copy(host, var_device);
  int count = 0;
  for (int d2 = 0; d2 < var_device.extent(1); ++d2) {
    for (int d1 = 0; d1 < var_device.extent(0); ++d1) {
      var_std[count] = host(d1, d2);
      count++;
    }
  }
}

void convert_1d_vector_to_1d_view_device(const std::vector<real_type> &var_std,
                                         const real_type_1d_view &var_device)
{
  auto var_host = real_type_1d_view_host((real_type *)var_std.data(), var_std.size());
  Kokkos::deep_copy(var_device, var_host);
}

void convert_1d_view_device_to_1d_vector(const real_type_1d_view &var_device,
                                         std::vector<real_type> &var_std)
{
  auto var_host = real_type_1d_view_host((real_type *)var_std.data(), var_std.size());
  Kokkos::deep_copy(var_host, var_device);
}

}   // namespace verification
} // namespace TChem
