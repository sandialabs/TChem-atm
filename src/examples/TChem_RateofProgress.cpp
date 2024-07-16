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
  std::string inputFile("None");
  std::string outputFile(prefix_path + "reaction_rates.dat");
  std::string outputFileTimes(prefix_path + "wall_times.json");
  int nBatch(1), team_size(-1), vector_size(-1);
  bool verbose(true);
  bool use_cloned_samples(false);

  /// parse command line arguments
  TChem::CommandLineParser opts("This example computes the rate of progress with a given state vector");
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

    if (inputFile == "None") {
       inputFile=chemFile;
    }

    printf("Using this input file: %s \n", inputFile.c_str());

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
    TChem::AtmChemistry::setScenarioConditions(inputFile,
     speciesNamesHost, kmcd.nSpec, stateVecDim, state_scenario_host, nbacth_files);

    // read photolysis reaction values
    // we assume photolysis reaction  are computed by another tool.
    real_type_2d_view_host photo_rates_scenario_host;
    ordinal_type n_photo_rates = 0;
    TChem::AtmChemistry::setScenarioConditionsPhotolysisReactions(inputFile,
             nbacth_files,
             // output
             photo_rates_scenario_host,
             n_photo_rates
             );

    real_type_2d_view state;
    real_type_2d_view photo_rates;

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

      if (n_photo_rates > 0 )
      {
        photo_rates = real_type_2d_view("photo rates", nBatch, n_photo_rates);
        auto photo_rates_host = Kokkos::create_mirror_view(photo_rates);
        auto photo_rates_scenario_host_at_0 = Kokkos::subview(photo_rates_scenario_host, 0, Kokkos::ALL);
        auto photo_rates_host_at_0 = Kokkos::subview(photo_rates_host, 0, Kokkos::ALL);
        Kokkos::deep_copy(photo_rates_host_at_0, photo_rates_scenario_host_at_0);
        TChem::Test::cloneView(photo_rates_host);
        Kokkos::deep_copy(photo_rates, photo_rates_host);
      } // n_photo_rates

    } else {
      nBatch = nbacth_files;
      state = real_type_2d_view("StateVector Devices", nBatch, stateVecDim);
      Kokkos::deep_copy(state, state_scenario_host);

      if (n_photo_rates > 0 )
      {
        nBatch = nbacth_files;
        photo_rates = real_type_2d_view("StateVector Devices", nBatch, n_photo_rates);
        Kokkos::deep_copy(photo_rates, photo_rates_scenario_host);

      }  // n_photo_rates
    }

    // output
    real_type_2d_view rate_of_progress("rate_of_progress", nBatch, kmcd.nReac);

    using policy_type = typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    Kokkos::Timer timer;

    FILE *fout_times = fopen(outputFileTimes.c_str(), "w");
    fprintf(fout_times, "{\n");
    fprintf(fout_times, " \"Rate of Progress\": \n {\n");
    // fprintf(fout_times, " Computation type, total time [sec], time per sample [sec/sample]\n ");

    policy_type policy(exec_space_instance, nBatch, Kokkos::AUTO());

    if (team_size > 0 && vector_size > 0) {
        policy = policy_type(exec_space_instance,  nBatch, team_size, vector_size);
    } else if (team_size > 0 && vector_size < 0) {
      // only set team size
       policy = policy_type(exec_space_instance, nBatch,  team_size);
    }

    const ordinal_type level = 1;
    const ordinal_type per_team_extent = TChem::RateofProgress::getWorkSpaceSize(kmcd);
    const ordinal_type per_team_scratch =
    TChem::Scratch<real_type_1d_view>::shmem_size(per_team_extent);
    policy.set_scratch_size(level, Kokkos::PerTeam(per_team_scratch));

    Kokkos::fence();
    /// skip the first run
    TChem::RateofProgress::runDeviceBatch(policy, state, photo_rates, rate_of_progress, kmcd);
    exec_space_instance.fence();
    timer.reset();
    TChem::RateofProgress::runDeviceBatch(policy, state, photo_rates, rate_of_progress, kmcd);
    exec_space_instance.fence();
    const real_type t_device_batch = timer.seconds();
    fprintf(fout_times, "%s: %20.14e, \n","\"wall_time\"", t_device_batch);
    fprintf(fout_times, "%s: %20.14e, \n","\"wall_time_per_sample\"", t_device_batch / real_type(nBatch));
    fprintf(fout_times, "%s: %d \n","\"number_of_samples\"", nBatch);
    fprintf(fout_times, "} \n ");// reaction rates


    if (verbose) {
    	auto rate_of_progress_host = Kokkos::create_mirror_view(rate_of_progress);
        Kokkos::deep_copy(rate_of_progress_host, rate_of_progress);

     if (use_cloned_samples) {
       FILE *fout = fopen(outputFile.c_str(), "w");
       auto rate_of_progress_host_at_0 = Kokkos::subview(rate_of_progress_host, 0, Kokkos::ALL());
       TChem::Test::writeReactionRates(outputFile, kmcd.nReac, rate_of_progress_host_at_0);
       fclose(fout);
     } else {
       FILE *fout = fopen(outputFile.c_str(), "w");
       TChem::Test::write2DMatrix(outputFile, nBatch, kmcd.nReac,rate_of_progress_host);
       fclose(fout);
     }

    }

    fprintf(fout_times, "}\n ");// end index time
    fclose(fout_times);




  }
  Kokkos::finalize();

  return 0;
}
