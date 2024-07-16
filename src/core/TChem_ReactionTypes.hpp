/* =====================================================================================
TChem-atm version 1.0
Copyright (2024) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC
(NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
Government retains certain rights in this software.

This file is part of TChem-atm. TChem-atm is open source software: you can redistribute
it and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Contact Oscar Diaz-Ibarra at <odiazib@sandia.gov>, or
           Mike Schmidt at <mjschm@sandia.gov>, or
           Cosmin Safta at <csafta@sandia.gov>

Sandia National Laboratories, New Mexico/Livermore, NM/CA, USA
=====================================================================================
*/
#ifndef __TCHEM_REACTION_TYPES_HPP__
#define __TCHEM_REACTION_TYPES_HPP__

#include "TChem_Util.hpp"
namespace TChem {


struct SIMPOL_PhaseTransferType{
  // index in state vector
  ordinal_type aerosol_sp_index;
  // mass transfer parameters
  real_type B1;
  real_type B2;
  real_type B3;
  real_type B4;
  // gas species info
  ordinal_type gas_sp_index;
  real_type diffusion_coeff;
  real_type N_star;
  bool compute_alpha;
  real_type molecular_weight;
};

using simpol_phase_transfer_type_1d_dual_view =
 Tines::value_type_1d_dual_view<SIMPOL_PhaseTransferType, exec_space>;

} // namespace TChem

#endif
