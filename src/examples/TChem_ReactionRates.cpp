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
#include "TChem.hpp"
#include "TChem_CommandLineParser.hpp"

using ordinal_type = TChem::ordinal_type;
using real_type = TChem::real_type;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_2d_view = TChem::real_type_2d_view;
using real_type_2d_view_host= TChem::real_type_2d_view_host;

int main(int argc, char *argv[]) {

  std::string prefix_path("");
  std::string chemFile(prefix_path + "chem.yaml");
  std::string inputFile(prefix_path + "input.dat");
  std::string outputFile(prefix_path + "reaction_rates.dat");
  std::string outputFileTimes(prefix_path + "wall_times.json");
  int nBatch(1), team_size(-1), vector_size(-1);
  bool verbose(true);
  bool use_cloned_samples(false);

  /// parse command line arguments
  TChem::CommandLineParser opts("This example computes reaction rates with a given state vector");
  opts.set_option<std::string>("chemfile", "Chem file name e.g., chem.yaml", &chemFile);
  opts.set_option<std::string>("inputfile", "Input state file name e.g., input.yaml", &inputFile);
  opts.set_option<std::string>("outputfile", "Output omega file name e.g., omega.dat", &outputFile);
  opts.set_option<std::string>("outputfile_times", "Wal times file name e.g., times.json", &outputFileTimes);
  opts.set_option<int>("team_thread_size", "time thread size ", &team_size);
  opts.set_option<int>("vector_thread_size", "vector thread size ", &vector_size);
  opts.set_option<int>("batchsize", "Batchsize the same state vector described in statefile is cloned", &nBatch);
  opts.set_option<bool>("verbose", "If true, printout the first omega values", &verbose);
  opts.set_option<int>("batch_size", " number of batches or samples, e.g. 10  ", &nBatch);
  opts.set_option<bool>(
      "use_cloned_samples", "If true, one state vector will be cloned.", &use_cloned_samples);

  const bool r_parse = opts.parse(argc, argv);
  if (r_parse)
    return 0; // print help return

  Kokkos::initialize(argc, argv);
  {

  	const bool detail = false;

    using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;

    TChem::exec_space().print_configuration(std::cout, detail);
    TChem::host_exec_space().print_configuration(std::cout, detail);
    const auto exec_space_instance = TChem::exec_space();

    /// construct kmd and use the view for testing
    TChem::KineticModelData kmd = TChem::KineticModelData(chemFile);
    const auto kmcd = TChem::createNCAR_KineticModelConstData<device_type>(kmd);

    const ordinal_type stateVecDim = TChem::Impl::getStateVectorSize(kmcd.nSpec);
    const auto speciesNamesHost = kmd.sNames_.view_host();


    // read scenario condition from yaml file
    real_type_2d_view_host state_scenario_host;
    ordinal_type nbacth_files=0;
    TChem::AtmChemistry::setScenarioConditions(chemFile,
     speciesNamesHost, kmcd.nSpec, stateVecDim, state_scenario_host, nbacth_files);

    real_type_2d_view state;

    if (nbacth_files == 1 && use_cloned_samples && nBatch > 1) {
      // only clone samples if nbacth_files is 1
      printf("-------------------------------------------------------\n");
      printf("--------------------Warning----------------------------\n");
      printf("Using cloned samples ... only for numerical experiments\n");
      state = real_type_2d_view("StateVector Devices", nBatch, stateVecDim);
      auto state_host = Kokkos::create_mirror_view(state);
      auto state_scenario_host_at_0 = Kokkos::subview(state_scenario_host, 0, Kokkos::ALL);
      auto state_host_at_0 = Kokkos::subview(state_host, 0, Kokkos::ALL);
      Kokkos::deep_copy(state_host_at_0, state_scenario_host_at_0);
      TChem::Test::cloneView(state_host);
      Kokkos::deep_copy(state, state_host);
    } else {
      nBatch = nbacth_files;
      state = real_type_2d_view("StateVector Devices", nBatch, stateVecDim);
      Kokkos::deep_copy(state, state_scenario_host);
    }

    // output
    real_type_2d_view kfwd("kfwd", nBatch, kmcd.nReac);

    using policy_type = typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    Kokkos::Timer timer;

    FILE *fout_times = fopen(outputFileTimes.c_str(), "w");
    fprintf(fout_times, "{\n");
    fprintf(fout_times, " \"Reaction Rates\": \n {\n");
    // fprintf(fout_times, " Computation type, total time [sec], time per sample [sec/sample]\n ");

    policy_type policy(exec_space_instance, nBatch, Kokkos::AUTO());

    if (team_size > 0 && vector_size > 0) {
        policy = policy_type(exec_space_instance,  nBatch, team_size, vector_size);
    } else if (team_size > 0 && vector_size < 0) {
      // only set team size
       policy = policy_type(exec_space_instance, nBatch,  team_size);
    }

    Kokkos::fence();
    /// skip the first run
    TChem::ReactionRates::runDeviceBatch(policy, state, kfwd, kmcd);
    exec_space_instance.fence();
    timer.reset();
    TChem::ReactionRates::runDeviceBatch(policy, state, kfwd, kmcd);
    exec_space_instance.fence();
    const real_type t_device_batch = timer.seconds();
    fprintf(fout_times, "%s: %20.14e, \n","\"wall_time\"", t_device_batch);
    fprintf(fout_times, "%s: %20.14e, \n","\"wall_time_per_sample\"", t_device_batch / real_type(nBatch));
    fprintf(fout_times, "%s: %d \n","\"number_of_samples\"", nBatch);
    fprintf(fout_times, "} \n ");// reaction rates


    if (verbose) {
    	auto kfwd_host = Kokkos::create_mirror_view(kfwd);
        Kokkos::deep_copy(kfwd_host, kfwd);

     {
       FILE *fout = fopen(outputFile.c_str(), "w");
       auto kfwd_host_at_0 = Kokkos::subview(kfwd_host, 0, Kokkos::ALL());
       TChem::Test::writeReactionRates(outputFile, kmcd.nReac, kfwd_host_at_0);
       fclose(fout);
     }
    }

    fprintf(fout_times, "}\n ");// end index time
    fclose(fout_times);




  }
  Kokkos::finalize();

  return 0;
}
