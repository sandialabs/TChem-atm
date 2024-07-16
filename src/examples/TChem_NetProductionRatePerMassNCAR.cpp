/* =====================================================================================
TChem-atm version 1.0
Copyright (2024) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.

This file is part of TChem-atm. TChem-atm is open source software: you can redistribute it
and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Contact Oscar Diaz-Ibarra at <odiazib@sandia.gov>, or
           Mike Schmidt at <mjschm@sandia.gov>, or
           Cosmin Safta at <csafta@sandia.gov>

Sandia National Laboratories, New Mexico/Livermore, NM/CA, USA
===================================================================================== */

#include "TChem_CommandLineParser.hpp"
#include "TChem.hpp"
#include "TChem_Impl_ReactionRates.hpp"
#include "TChem_Impl_KForward.hpp"
#include "TChem_Impl_KForwardJPL.hpp"
#include "TChem_Impl_AdjustReactions.hpp"
#include "TChem_Impl_RateofProgress.hpp"
#include "TChem_Impl_NetProductionRates.hpp"
#include "TChem_AtmosphericChemistry.hpp"

using ordinal_type = TChem::ordinal_type;
using real_type = TChem::real_type;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_2d_view = TChem::real_type_2d_view;
using real_type_2d_view_host = TChem::real_type_2d_view_host;

int
main(int argc, char* argv[])
{
  #if defined(TCHEM_ATM_ENABLE_TPL_YAML_CPP)
  /// default inputs
  std::string prefixPath("");
  std::string chemFile(prefixPath + "chem.yaml");
  std::string inputFile(prefixPath + "input.dat");
  std::string outputFile(prefixPath + "omega.dat");
  std::string thermFile(prefixPath + "therm.dat");
  int nBatch(1);
  bool verbose(true);
  bool useYaml(true);
  bool use_sample_format(false);
  std::string test("arrhenius");

  /// parse command line arguments
  TChem::CommandLineParser opts(
    "This example computes reaction rates with a given state vector");
  opts.set_option<std::string>(
    "chemfile", "Chem file name e.g., chem.yaml", &chemFile);
  opts.set_option<std::string>(
    "inputfile", "Input state file name e.g., input.dat", &inputFile);
  opts.set_option<std::string>(
    "outputfile", "Output omega file name e.g., omega.dat", &outputFile);
  opts.set_option<std::string>(
    "thermfile", "Therm file name e.g., therm.dat", &thermFile);
  //
  opts.set_option<std::string>(
    "unit-test", "unit test name e.g., arrhenius", &test);
  opts.set_option<int>(
    "batchsize",
    "Batchsize the same state vector described in statefile is cloned",
    &nBatch);
  opts.set_option<bool>(
    "verbose", "If true, printout the first omega values", &verbose);
  opts.set_option<bool>(
    "useYaml", "If true, use yaml to parse input file", &useYaml);
  opts.set_option<bool>(
    "use_sample_format", "If true, input file does not header or format", &use_sample_format);

  const bool r_parse = opts.parse(argc, argv);
  if (r_parse)
    return 0; // print help return

  Kokkos::initialize(argc, argv);
  {
    const bool detail = false;

    TChem::exec_space().print_configuration(std::cout, detail);
    TChem::host_exec_space().print_configuration(std::cout, detail);

    using host_device_type      = typename Tines::UseThisDevice<TChem::host_exec_space>::type;

    /// construct kmd and use the view for testing
    TChem::KineticModelData kmd = TChem::KineticModelData(chemFile);
    const auto kmcd = TChem::createNCAR_KineticModelConstData<host_device_type>(kmd);
    const auto member = Tines::HostSerialTeamMember();

    const auto speciesNamesHost = Kokkos::create_mirror_view(kmcd.speciesNames);
    Kokkos::deep_copy(speciesNamesHost, kmcd.speciesNames);


    // forward constant
    using kForward_type = TChem::Impl::KForward<real_type, host_device_type >;
    using value_type_1d_view_type = typename kForward_type::value_type_1d_view_type;

    value_type_1d_view_type kfor("forward reate constant", kmcd.nReac);

    using reaction_rates_type = TChem::Impl::ReactionRatesAerosol<real_type, host_device_type >;

    value_type_1d_view_type omega("omega", kmcd.nSpec);

    value_type_1d_view_type work("work", 3*kmcd.nReac);
    const ordinal_type stateVecDim = kmcd.nSpec +3;

    // read scenario condition from yaml file
    real_type_2d_view_host state_host;
    TChem::AtmChemistry
         ::setScenarioConditions(chemFile, speciesNamesHost,
                                 kmcd.nSpec, stateVecDim, state_host, nBatch );


    // TODO:
    // create a map or dic that returns the index of a variable in the state vector.
    // for example "temperature" -> 2 or "M" -> 76+3
    // FIXME
    const auto t = state_host(0,2);
    const auto p = state_host(0,1);

    using range_type = Kokkos::pair<ordinal_type, ordinal_type>;
    const auto x =
        Kokkos::subview(state_host, 0, range_type(3, 3 + kmcd.nSpec));

    const auto m = x(kmcd.M_index);

    printf("Using value of m = %e \n", m);

    kForward_type::team_invoke(member, t, p, kfor, kmcd);
    using kForwardJPL_type = TChem::Impl::KForwardJPL<real_type, host_device_type >;
    kForwardJPL_type::team_invoke(member, t, m, kfor, kmcd);
    std::string output_file_kfor(prefixPath + "kfwd.dat");

    TChem::Test::writeReactionRates(output_file_kfor, kmcd.nReac, kfor);

    std::string output_file_kfor_wadjust_reaction(prefixPath + "kfwd_wadjust_reactions.dat");
    using AdjustReactions_type = TChem::Impl::AdjustReactions<real_type, host_device_type >;

    printf("x(%d)  %e \n", kmcd.M_index,m);

    AdjustReactions_type::team_invoke(member, t, x, kfor, kmcd);

    TChem::Test::writeReactionRates(output_file_kfor_wadjust_reaction, kmcd.nReac, kfor);

    using ReactionRates_type = TChem::Impl::ReactionRates<real_type, host_device_type >;

    for (int i = 0; i < kmcd.nReac; ++i)
    {
      kfor(i)=0;
    }

    std::string output(prefixPath + "reaction_rates.dat");
    ReactionRates_type::team_invoke(member, t, p, x, kfor, kmcd);
    TChem::Test::writeReactionRates(output, kmcd.nReac, kfor);

    for (int i = 0; i < kmcd.nReac; ++i)
    {
      kfor(i)=0;
    }

    value_type_1d_view_type rate_of_progress("rate_of_progress", kmcd.nReac);

    // for (int i = 73; i < kmcd.nSpec; ++i)
    // {
    //   x(i)=1;
    // }

    ordinal_type n_photo_rates = 0;
    value_type_1d_view_type photo_rates("photo_rates", n_photo_rates);

    using rateof_progress_type = TChem::Impl::RateofProgress<real_type, host_device_type >;
    rateof_progress_type::team_invoke(member, t, p, x, photo_rates, rate_of_progress, kfor, kmcd);



    std::string output_rop(prefixPath + "rate_of_progress.dat");
    TChem::Test::writeReactionRates(output_rop, kmcd.nReac, rate_of_progress);

    // for (int i = 0; i < 73; ++i)
    // {
    //   printf("%.15e, & ! (%d) %s \n",x(i),i,&speciesNamesHost(i, 0));
    // }
    // //  invariants
    // printf("Invariants \n");
    // for (int i = 73; i < kmcd.nSpec; ++i)
    // {
    //   printf("%.15e, & ! (%d) %s \n",x(i),i,&speciesNamesHost(i, 0));
    // }

    for (int i = 0; i < kmcd.nReac; ++i)
    {
      kfor(i)=0;
      rate_of_progress(i)=0;
    }

    value_type_1d_view_type net_production_rates("net_production_rate", kmcd.nSpec);
    const ordinal_type n_active_vars = kmcd.nSpec - kmcd.nConstSpec;
    value_type_1d_view_type external_sources("external_sources", n_active_vars);

    using net_production_rate_type = TChem::Impl::NetProductionRates<real_type, host_device_type >;
    net_production_rate_type::team_invoke_detail(member, t, p, x, photo_rates, external_sources, net_production_rates, rate_of_progress, kfor, kmcd);

    std::string output_net_production_rates(prefixPath + "net_production_rate.dat");
    TChem::Test::writeReactionRates(output_net_production_rates, kmcd.nSpec, net_production_rates);


  }
  Kokkos::finalize();

  #else
   printf("This example requires Yaml ...\n" );
  #endif

  return 0;
}
