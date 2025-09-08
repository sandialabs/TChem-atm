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
#ifndef TCHEM_ATM_VERIFICATION_HPP
#define TCHEM_ATM_VERIFICATION_HPP

#include "TChem.hpp"

namespace TChem {
using real_type = TChem::real_type;
using real_type_2d_view = TChem::real_type_2d_view;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_1d_view_host = TChem::real_type_1d_view_host;
namespace verification {

/// Call this function to initialize a validation driver.
void initialize(int argc, char **argv);

/// Call this function to finalize a validation driver.
void finalize();

std::string output_name(const std::string &input_file);

void convert_1d_vector_to_2d_view_device(const std::vector<real_type> &var_std,
                                         const real_type_2d_view &var_device);

void convert_2d_view_device_to_1d_vector(const real_type_2d_view &var_device,
                                         std::vector<real_type> &var_std);

void convert_1d_vector_to_1d_view_device(const std::vector<real_type> &var_std,
                                         const real_type_1d_view &var_device);

void convert_1d_view_device_to_1d_vector(const real_type_1d_view &var_device,
                                         std::vector<real_type> &var_std);
}// namespace verification
} // namespace TChem

#endif
