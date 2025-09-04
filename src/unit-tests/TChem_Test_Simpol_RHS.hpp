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
#ifndef __TCHEM_TEST_SIMPOL_RHS_HPP__
#define __TCHEM_TEST_SIMPOL_RHS_HPP__

#include "TChem.hpp"
#include "TChem_Impl_SIMPOL_constant.hpp"
#include "TChem_Impl_SIMPOL_phase_transfer.hpp"
using ordinal_type = TChem::ordinal_type;
using real_type = TChem::real_type;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_2d_view = TChem::real_type_2d_view;
using real_type_2d_view_host = TChem::real_type_2d_view_host;

#define TCHEM_TEST_IMPL_SACADO
namespace TChem {
namespace Test {
void static SimpolRHS_test()
{
    using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
    // input files need to located in same directory as the executable.
    std::string chemFile ="../examples/runs/atmospheric_chemistry/Simpol/config_gas.yaml";
    std::string aerochemFile="../examples/runs/atmospheric_chemistry/Simpol/test_SIMPOL_phase_transfer.yaml";

    // construct kmd; gas phase
    TChem::KineticModelData kmd(chemFile);
    const auto kmcd = TChem::createNCAR_KineticModelConstData<device_type>(kmd);
    // construct amd; aerosols
    TChem::AerosolModelData amd = TChem::AerosolModelData(aerochemFile, kmd);
    const auto amcd = TChem::create_AerosolModelConstData<device_type>(amd);

#if defined(TCHEM_TEST_IMPL_SACADO)
    // length of value type is needed at compilation time.
    using value_type = Sacado::Fad::SLFad<real_type,128>;
#else
   using value_type = real_type;
#endif
    using SIMPOL_single_particle_type = TChem::Impl::SIMPOL_single_particle<value_type, device_type >;


    using value_type_1d_view_type = typename SIMPOL_single_particle_type::value_type_1d_view_type;
    using real_type_1d_view_type = typename SIMPOL_single_particle_type::real_type_1d_view_type;

    real_type_1d_view_type number_conc("number_conc", amcd.nParticles);
    // assuming constant number concentration
    Kokkos::deep_copy(number_conc, 1.3e6); // particle number concentration (#/cc)

    ordinal_type ntotal_species = amcd.nSpec_gas + amcd.nSpec*amcd.nParticles;

    value_type_1d_view_type state("state", ntotal_species);
    // by default view are initialized with zeros.
    value_type_1d_view_type omega("omega", ntotal_species);

#if defined(TCHEM_TEST_IMPL_SACADO)
    // this view is used to copy value of omega; omegas has value and its derivatives.
    real_type_1d_view_type omega_out("omega", ntotal_species);
#else
    auto omega_out = omega;
#endif
    const auto exec_space_instance = TChem::exec_space();

    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    // FIXME:
    // is this a bug in nvcc?
    //https://github.com/google/googletest/issues/4104
 #if  1
    // team policy
    ordinal_type nBatch =1;
    policy_type policy(exec_space_instance, nBatch, Kokkos::AUTO());
    Kokkos::parallel_for
      ("SIMPOL RHS",
       policy,
       KOKKOS_LAMBDA(const typename policy_type::member_type& member) {

       real_type t = 272.5; // K
       real_type p = 101253.3; // pa

       // initial concentration
       real_type ethanol=0.1;
       real_type ethanol_aq = 1.0e-8 ;
       real_type H2O_aq = 1.4e-2;

       // gas species
       state(0) = ethanol;
       // aerosol species
       for (int i_part = 0; i_part < amcd.nParticles; i_part++)
       {
        state(1+amcd.nSpec*i_part) = ethanol_aq/number_conc(i_part);
        state(2+amcd.nSpec*i_part) = H2O_aq/number_conc(i_part);
       }// i_part
       member.team_barrier();

       ordinal_type i_simpol=0;
       Kokkos::parallel_for(
      Kokkos::TeamThreadRange(member, amcd.nParticles),
       [&](const ordinal_type& i_part) {
           SIMPOL_single_particle_type
          ::team_invoke(member, i_part, i_simpol, t, p, number_conc, state, omega, amcd);
        });
        member.team_barrier();
#if defined(TCHEM_TEST_IMPL_SACADO)
        //copy values of omega
        for (ordinal_type i = 0; i < ntotal_species; i++)
          omega_out(i) = omega(i).val();

#endif
  });

#endif
  // we need to copy data from device to host.
  auto omega_host = Kokkos::create_mirror_view_and_copy(TChem::host_exec_space(), omega_out);
  // save results to a file.
  std::string outputFile ="Simpol_RHS_DEVICE.txt";
  FILE *fout = fopen(outputFile.c_str(), "w");
  TChem::Test::writeReactionRates(outputFile, omega_host.extent(0), omega_host);
  fclose(fout);

}
}// namespace Test

}// namespace TChem
TEST(SimpolRHS, Device)
{
  TChem::Test::SimpolRHS_test();
}

TEST(SimpolRHS, verification_device)
{
 /// pass test with a relative error of
  const real_type threshold =1e-12;
  EXPECT_TRUE(TChem::Test::compareFilesValues("Simpol_RHS_DEVICE.txt",
          "references/partmc_integration/simpol/Simpol_RHS_HOST.txt",threshold)
        );
}

#endif
