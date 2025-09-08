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
#ifndef __TCHEM_TEST_SIMPOL_RHS_HOST_HPP__
#define __TCHEM_TEST_SIMPOL_RHS_HOST_HPP__

#include "TChem.hpp"
#include "TChem_Impl_SIMPOL_constant.hpp"
#include "TChem_Impl_SIMPOL_phase_transfer.hpp"
using ordinal_type = TChem::ordinal_type;
using real_type = TChem::real_type;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_2d_view = TChem::real_type_2d_view;
using real_type_2d_view_host = TChem::real_type_2d_view_host;

TEST(SimpolRHS, single_host)
{
    using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;

    std::string chemFile ="../examples/runs/atmospheric_chemistry/Simpol/config_gas.yaml";
    std::string aerochemFile="../examples/runs/atmospheric_chemistry/Simpol/test_SIMPOL_phase_transfer.yaml";

    // construct kmd and use the view for testing
    TChem::KineticModelData kmd(chemFile);
    const auto kmcd = TChem::createNCAR_KineticModelConstData<host_device_type>(kmd);

    TChem::AerosolModelData amd = TChem::AerosolModelData(aerochemFile, kmd);
    const auto amcd = TChem::create_AerosolModelConstData<host_device_type>(amd);

    using SIMPOL_single_particle_type = TChem::Impl::SIMPOL_single_particle<real_type, host_device_type >;
    using SIMPOL_constant_type = TChem::Impl::SIMPOL_constant<real_type, host_device_type >;

    const auto member = Tines::HostSerialTeamMember();
    // inputs
    real_type t = 272.5; // K
    real_type p = 101253.3; // pa
    // outputs
    real_type alpha=0;
    real_type mfp_m=0;
    real_type KGM3_TO_PPM_=0;
    real_type equil_constant=0;

    SIMPOL_constant_type::team_invoke(member, t,
     p, alpha, mfp_m, KGM3_TO_PPM_, equil_constant, amcd.simpol_params(0));
#if defined(TCHEM_ENABLE_SERIAL_TEST_OUTPUT)
    printf("amcd.simpol_params(0).B1 %e \n", amcd.simpol_params(0).B1);
    printf("alpha %e \n", alpha);
    printf("mfp_m %e \n", mfp_m);
    printf("KGM3_TO_PPM_ %e \n", KGM3_TO_PPM_);
    printf("equil_constant %e \n", equil_constant);
#endif
    using value_type_1d_view_type = typename SIMPOL_single_particle_type::value_type_1d_view_type;
    using real_type_1d_view_type = typename SIMPOL_single_particle_type::real_type_1d_view_type;
    real_type_1d_view_type number_conc("number_conc", amcd.nParticles);
    Kokkos::deep_copy(number_conc, 1.3e6); // particle number concentration (#/cc)
    // initial concentration
    real_type ethanol=0.1;
    real_type ethanol_aq = 1.0e-8 ;
    real_type H2O_aq = 1.4e-2;
    ordinal_type ntotal_species = amcd.nSpec_gas + amcd.nSpec*amcd.nParticles;

    value_type_1d_view_type state("state", ntotal_species);
    // gas species
    state(0) = ethanol;
    // aerosol species
    for (int i_part = 0; i_part < amcd.nParticles; i_part++)
    {
      state(1+amcd.nSpec*i_part) = ethanol_aq/number_conc(i_part);
      state(2+amcd.nSpec*i_part) = H2O_aq/number_conc(i_part);
    }

    value_type_1d_view_type omega("omega", ntotal_species);

    ordinal_type i_simpol=0;
    for (int i_part = 0; i_part < amcd.nParticles; i_part++)
    {
#if defined(TCHEM_ENABLE_SERIAL_TEST_OUTPUT)
    printf("----Working on particle No %d ---\n", i_part);
#endif
    SIMPOL_single_particle_type
    ::team_invoke(member, i_part,i_simpol, t, p, number_conc, state, omega, amcd);
    }

  //save results to a file.
  std::string outputFile ="Simpol_RHS_HOST.txt";
  FILE *fout = fopen(outputFile.c_str(), "w");
  TChem::Test::writeReactionRates(outputFile, omega.extent(0), omega);
  fclose(fout);
}
TEST(SimpolRHS, verification_host)
{
 /// pass test with a relative error of
  const real_type threshold =1e-12;
  EXPECT_TRUE(TChem::Test::compareFilesValues("Simpol_RHS_HOST.txt",
          "references/partmc_integration/simpol/Simpol_RHS_HOST.txt",threshold)
        );
}
#endif
