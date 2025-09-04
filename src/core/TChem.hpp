/* =====================================================================================
TChem-atm version 2.0.0
Copyright (2025) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2025 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.

This file is part of TChem-atm. TChem-atm is open source software: you can redistribute it
and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Contact Oscar Diaz-Ibarra at <odiazib@sandia.gov>, or
           Cosmin Safta at <csafta@sandia.gov> or,
           Nicole Riemer at <nriemer@illinois.edu> or,
           Matthew West at <mwest@illinois.edu>

Sandia National Laboratories, New Mexico/Livermore, NM/CA, USA
===================================================================================== */
#ifndef __TCHEM_HPP__
#define __TCHEM_HPP__

#include "TChem_Util.hpp"
#include "TChem_KineticModelData.hpp"
#include "TChem_AerosolModelData.hpp"
#include "TChem_AerosolChemistry_CVODE.hpp"
#include "TChem_AerosolChemistry_KokkosKernels.hpp"
#include "TChem_AerosolChemistry.hpp"
#include "TChem_AtmosphericChemistry.hpp"
#include "TChem_AtmosphericChemistryE3SM.hpp"
#include "TChem_AtmosphericChemistryE3SM_CVODE.hpp"
#include "TChem_AtmosphericChemistryE3SM_ExplicitEuler.hpp"
#include "TChem_AtmosphericChemistryE3SM_ImplicitEuler.hpp"
#include "TChem_ReactionRates.hpp"
#include "TChem_ReactionTypes.hpp"
#include "TChem_RateofProgress.hpp"
#include "TChem_NetProductionRates.hpp"
#include "TChem_Linv3StratosphereSolver.hpp"
#include "TChem_AerosolChemistry_CVODE_RHS_Jacobian.hpp"




#endif
