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
#include "TChem.hpp"
#include "TChem_CommandLineParser.hpp"
#include "TChem_Impl_AerosolChemistry.hpp"

using ordinal_type = TChem::ordinal_type;
using real_type = TChem::real_type;

int main(int argc, char *argv[]) {
  std::string chemFile( "chem.yaml");
  std::string aeroFile("aero.yaml");
  std::string outputFile("AerosolChemistry.dat");
  std::string inputFile("None");
  std::string input_file_particles("None");
  std::string outputFileTimes("wall_times.json");
  int nBatch(1), team_size(-1), vector_size(-1);
  bool verbose(true);
  bool use_cloned_samples(false);
  int number_of_particles(-1);
  bool do_jac(true);
  bool do_rhs(true);

  /// parse command line arguments
  TChem::CommandLineParser opts("This example computes the net production rates with a given state vector");
  opts.set_option<std::string>("chemfile", "Chem file name e.g., chem.inp",
                               &chemFile);
  opts.set_option<std::string>("aerofile", "Aerosol Chem file name e.g., aero.yaml",
                               &aeroFile);
  opts.set_option<std::string>("inputfile", "Chem file name e.g., chem.inp",
                               &inputFile);
  opts.set_option<std::string>("inputfile_particles", "input file name e.g., chem.inp",
                               &input_file_particles);
  opts.set_option<std::string>("outputfile", "Output omega file name e.g., omega.dat", &outputFile);
  opts.set_option<std::string>("outputfile_times", "Wal times file name e.g., times.json", &outputFileTimes);
  opts.set_option<int>("team_thread_size", "time thread size ", &team_size);
  opts.set_option<int>("vector_thread_size", "vector thread size ", &vector_size);
  opts.set_option<int>("number_of_particles", "Set the number of particles; this will overwrite the value from the input file. ", &number_of_particles);
  opts.set_option<bool>("verbose", "If true, printout the first omega values", &verbose);
  opts.set_option<int>("batch_size", " number of batches or samples, e.g. 10  ", &nBatch);
  opts.set_option<bool>(
      "use_cloned_samples", "If true, one state vector will be cloned.", &use_cloned_samples);
  opts.set_option<bool>("do_jac", "Evaluate Jacobian matrix", &do_jac);
  opts.set_option<bool>("do_rhs", "Evaluate RHS", &do_rhs);


  const bool r_parse = opts.parse(argc, argv);
  if (r_parse)
    return 0; // print help return

  Kokkos::initialize(argc, argv);
  {

    // by default use input condition from chemFile
    if (inputFile == "None") {
       inputFile=chemFile;
    }
    // check first in chem file
    if (input_file_particles == "None") {
       input_file_particles=chemFile;
    }

  	const bool detail = false;
    constexpr real_type zero =0;

    using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;

    using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
    using problem_type = TChem::Impl::AerosolChemistry_Problem<real_type,device_type>;

    using real_type_0d_view_type = Tines::value_type_0d_view<real_type, device_type>;
    using real_type_1d_view_type = Tines::value_type_1d_view<real_type, device_type>;
    using real_type_2d_view_type = Tines::value_type_2d_view<real_type, device_type>;
    using real_type_3d_view_type = Tines::value_type_3d_view<real_type, device_type>;
    using real_type_2d_view_host_type = Tines::value_type_2d_view<real_type, host_device_type>;


    TChem::exec_space().print_configuration(std::cout, detail);
    TChem::host_exec_space().print_configuration(std::cout, detail);
    const auto exec_space_instance = TChem::exec_space();

    /// construct kmd and use the view for testing
    printf("kmd parsing ...\n");
    TChem::KineticModelData kmd(chemFile);
    const auto kmcd = TChem::createNCAR_KineticModelConstData<device_type>(kmd);

    printf("amd parsing ...\n");
    TChem::AerosolModelData amd(aeroFile, kmd);
    if(number_of_particles > 0) {
      amd.setNumberofParticles(number_of_particles);
    }
    const auto amcd = TChem::create_AerosolModelConstData<device_type>(amd);

    const ordinal_type total_n_species =kmcd.nSpec + amcd.nSpec * amcd.nParticles;
    const ordinal_type number_of_equations = problem_type::getNumberOfTimeODEs(kmcd, amcd);
    const ordinal_type stateVecDim =
        TChem::Impl::getStateVectorSize(total_n_species);
    const auto speciesNamesHost = Kokkos::create_mirror_view(kmcd.speciesNames);
    Kokkos::deep_copy(speciesNamesHost, kmcd.speciesNames);

    printf("Number of Gas Species %d \n", kmcd.nSpec);
    printf("Number of Gas Reactions %d \n", kmcd.nReac);

    printf("Number of Aerosol Species %d \n", amcd.nSpec);
    printf("Number of Aerosol Particles %d \n", amcd.nParticles);

    printf("stateVecDim %d \n", stateVecDim);

    // read scenario condition from yaml file
    real_type_2d_view_host_type state_scenario_host;
    ordinal_type nbatch_files=1;
    printf("conditions parsing ...\n");
    TChem::AtmChemistry::setScenarioConditions(inputFile,
     speciesNamesHost, kmcd.nSpec, stateVecDim, state_scenario_host, nbatch_files);
    printf("Number of const species %d \n", kmcd.nConstSpec);

    real_type_2d_view_host_type num_concentration_host;
    real_type_2d_view_host_type state_host;

    // scenario particles
    amd.scenarioConditionParticles(input_file_particles,
                                   nbatch_files,
                                   num_concentration_host,
                                   state_scenario_host);

    real_type_2d_view_type num_concentration;
    real_type_2d_view_type state;

    if (nbatch_files == 1 && use_cloned_samples && nBatch > 1) {
      // only clone samples if nbatch_files is 1
      printf("-------------------------------------------------------\n");
      printf("--------------------Warning----------------------------\n");
      printf("Using cloned samples ... only for numerical experiments\n");
      state = real_type_2d_view_type("StateVector Devices", nBatch, stateVecDim);
      state_host = Kokkos::create_mirror_view(state);
      auto state_scenario_host_at_0 = Kokkos::subview(state_scenario_host, 0, Kokkos::ALL);
      auto state_at_0 = Kokkos::subview(state, 0, Kokkos::ALL);
      Kokkos::deep_copy(state_at_0, state_scenario_host_at_0);
      TChem::Test::cloneView(state);
      //make sure to copy state to host view.
      Kokkos::deep_copy(state_host, state);
      //scenario particles
      auto num_concentration_host_at_0 = Kokkos::subview(num_concentration_host, 0, Kokkos::ALL);
      num_concentration = real_type_2d_view_type("num_concentration", nBatch, amd.nParticles_);
      auto num_concentration_at_0 = Kokkos::subview(num_concentration, 0, Kokkos::ALL);
      Kokkos::deep_copy(num_concentration_at_0, num_concentration_host_at_0);
      TChem::Test::cloneView(num_concentration);

    } else {
      nBatch = nbatch_files;
      state = real_type_2d_view_type("StateVector Devices", nBatch, stateVecDim);
      Kokkos::deep_copy(state, state_scenario_host);
      state_host = state_scenario_host;

      // scenario particles
      num_concentration = real_type_2d_view_type("num_concentration", nBatch, amd.nParticles_);
      Kokkos::deep_copy(num_concentration, num_concentration_host);
    }
    // output
    real_type_2d_view_type fac("fac", nBatch, number_of_equations);
    real_type_2d_view_type rhs("rhs", nBatch, number_of_equations);

    using policy_type = typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;

    policy_type policy(exec_space_instance, nBatch, Kokkos::AUTO());

    if (team_size > 0 && vector_size > 0) {
        policy = policy_type(exec_space_instance,  nBatch, team_size, vector_size);
    } else if (team_size > 0 && vector_size < 0) {
      // only set team size
       policy = policy_type(exec_space_instance, nBatch,  team_size);
    }

    using range_type = Kokkos::pair<ordinal_type, ordinal_type>;
    const ordinal_type level = 1;
    ordinal_type per_team_extent=0.0;
    if (do_jac){
      per_team_extent = problem_type::getWorkSpaceSize(kmcd,amcd);
    }

    if (do_rhs){
      per_team_extent
         = TChem::Impl::Aerosol_RHS<real_type, device_type>::getWorkSpaceSize(kmcd, amcd);
    }
    const ordinal_type per_team_scratch =
    TChem::Scratch<real_type_1d_view_type>::shmem_size(per_team_extent);
    policy.set_scratch_size(level, Kokkos::PerTeam(per_team_scratch));

        real_type_2d_view_type y2d("y2d", nBatch, number_of_equations);

    Kokkos::parallel_for(
      "fill_y", Kokkos::RangePolicy<TChem::exec_space>(0, nBatch),
    KOKKOS_LAMBDA(const int i) {
        const ordinal_type n_active_gas_species = kmcd.nSpec - kmcd.nConstSpec;
        const real_type_1d_view_type state_at_i =
        Kokkos::subview(state, i, Kokkos::ALL());
        const ordinal_type total_n_species = kmcd.nSpec + amcd.nParticles*amcd.nSpec;
        TChem::Impl::StateVector<real_type_1d_view_type> sv_at_i(total_n_species, state_at_i);
        const real_type_1d_view_type Ys = sv_at_i.MassFractions();

        const auto activeYs = Kokkos::subview(Ys, range_type(0, n_active_gas_species));
        const real_type_1d_view_type partYs = Kokkos::subview(Ys, range_type(kmcd.nSpec, total_n_species));

        for (ordinal_type j=0;j<n_active_gas_species;++j){
          y2d(i, j) = activeYs(j);
        }

        for (ordinal_type j=n_active_gas_species;j<total_n_species- kmcd.nConstSpec;++j)
        {
          y2d(i, j) = partYs(j-n_active_gas_species);
        }

    });




    Kokkos::Timer timer;
    FILE *fout_times = fopen(outputFileTimes.c_str(), "w");
    fprintf(fout_times, "{\n");

    if (do_rhs){
      printf("..evaluating RHS\n");
      fprintf(fout_times, " \"Aerosol RHSs\": \n {\n");
      const ordinal_type level = 1;
      const std::string profile_name = "TChem::AerosolChemistry::RHS_evaluation";
      Kokkos::fence();
      timer.reset();
      Kokkos::Profiling::pushRegion(profile_name);
      Kokkos::parallel_for
      (profile_name,
       policy,
       KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
        const ordinal_type i = member.league_rank();

        const ordinal_type m = problem_type::getNumberOfTimeODEs(kmcd,amcd);
        const real_type_1d_view_type rhs_at_i =
        Kokkos::subview(rhs, i, Kokkos::ALL());

        const real_type_1d_view_type state_at_i =
        Kokkos::subview(state, i, Kokkos::ALL());

        const real_type_1d_view_type vals =
        Kokkos::subview(y2d, i, Kokkos::ALL());

        const real_type_1d_view_type number_conc_at_i =
        Kokkos::subview(num_concentration, i, Kokkos::ALL());
        TChem::Scratch<real_type_1d_view_type> work(member.team_scratch(level),
                                       per_team_extent);
        auto wptr = work.data();
        const ordinal_type total_n_species = kmcd.nSpec + amcd.nParticles*amcd.nSpec;
        TChem::Impl::StateVector<real_type_1d_view_type> sv_at_i(total_n_species, state_at_i);

        const real_type temperature = sv_at_i.Temperature();
        const real_type pressure = sv_at_i.Pressure();
        const real_type_1d_view_type Ys = sv_at_i.MassFractions();
        const ordinal_type n_active_gas_species = kmcd.nSpec - kmcd.nConstSpec;
        const real_type_1d_view_type constYs = Kokkos::subview(Ys,
              range_type(n_active_gas_species, kmcd.nSpec));

        /// problem workspace
        /// problem setup
        const ordinal_type problem_workspace_size
         = TChem::Impl::Aerosol_RHS<real_type, device_type>::getWorkSpaceSize(kmcd, amcd);
        auto pw = real_type_1d_view_type(wptr, problem_workspace_size);
        wptr +=problem_workspace_size;
        TChem::Impl::Aerosol_RHS<real_type, device_type>
        ::team_invoke(member, temperature, pressure,
                      number_conc_at_i, vals,
                      constYs,  rhs_at_i, pw, kmcd, amcd);
      });

      Kokkos::Profiling::popRegion();
      exec_space_instance.fence();
      const real_type t_device_batch = timer.seconds();
      fprintf(fout_times, "%s: %20.14e, \n","\"wall_time\"", t_device_batch);
      fprintf(fout_times, "%s: %20.14e, \n","\"wall_time_per_sample\"", t_device_batch / real_type(nBatch));
      fprintf(fout_times, "%s: %d \n","\"number_of_samples\"", nBatch);

      if (do_jac){
          fprintf(fout_times, "}, \n ");// reaction rates
      }
      else {
          fprintf(fout_times, "} \n ");// reaction rates
      }

    }

    if (do_jac){
      printf("..evaluating Jacobian\n");
      real_type_3d_view_type jacobian("jacobian", nBatch, number_of_equations, number_of_equations);
      const ordinal_type level = 1;
      fprintf(fout_times, " \"Aerosol Numerical Jacobian\": \n {\n");
      const std::string profile_name = "TChem::AerosolChemistry::NumericalJacobian_evaluation";
      Kokkos::fence();
      timer.reset();
      Kokkos::Profiling::pushRegion(profile_name);
      Kokkos::parallel_for
      (profile_name,
       policy,
       KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
        const ordinal_type i = member.league_rank();

        const ordinal_type m = problem_type::getNumberOfTimeODEs(kmcd,amcd);
        const real_type_2d_view_type jacobian_at_i =
        Kokkos::subview(jacobian, i, Kokkos::ALL(), Kokkos::ALL());

        const real_type_1d_view_type state_at_i =
        Kokkos::subview(state, i, Kokkos::ALL());

        const real_type_1d_view_type fac_at_i =
        Kokkos::subview(fac, i, Kokkos::ALL());

        const real_type_1d_view_type vals =
        Kokkos::subview(y2d, i, Kokkos::ALL());

        const real_type_1d_view_type number_conc_at_i =
        Kokkos::subview(num_concentration, i, Kokkos::ALL());
        TChem::Scratch<real_type_1d_view_type> work(member.team_scratch(level),
                                       per_team_extent);
        auto wptr = work.data();
        const ordinal_type total_n_species = kmcd.nSpec + amcd.nParticles*amcd.nSpec;
        TChem::Impl::StateVector<real_type_1d_view_type> sv_at_i(total_n_species, state_at_i);

        const real_type temperature = sv_at_i.Temperature();
        const real_type pressure = sv_at_i.Pressure();
        const real_type_1d_view_type Ys = sv_at_i.MassFractions();
        const ordinal_type n_active_gas_species = kmcd.nSpec - kmcd.nConstSpec;
        const real_type_1d_view_type constYs = Kokkos::subview(Ys,
              range_type(n_active_gas_species, kmcd.nSpec));
        /// problem workspace
        /// problem setup
        const ordinal_type problem_workspace_size = problem_type::getWorkSpaceSize(kmcd,amcd);
        auto pw = real_type_1d_view_type(wptr, problem_workspace_size);
        wptr +=problem_workspace_size;
        problem_type problem;
        problem._kmcd = kmcd;
        problem._amcd = amcd;

        /// initialize problem
        problem._fac = fac_at_i;
        problem._work = pw;
        problem._temperature= temperature;
        problem._pressure =pressure;
        problem._const_concentration= constYs;
        problem._number_conc =number_conc_at_i;
        // active gas species
        problem.computeNumericalJacobian(member,vals,jacobian_at_i);
      });
      Kokkos::Profiling::popRegion();
      exec_space_instance.fence();
      const real_type t_device_batch = timer.seconds();
      fprintf(fout_times, "%s: %20.14e, \n","\"wall_time\"", t_device_batch);
      fprintf(fout_times, "%s: %20.14e, \n","\"wall_time_per_sample\"", t_device_batch / real_type(nBatch));
      fprintf(fout_times, "%s: %d \n","\"number_of_samples\"", nBatch);
      fprintf(fout_times, "} \n ");// reaction rates
    }

    if (verbose & do_rhs) {


     if (use_cloned_samples) {
       auto rhs_at_0 = Kokkos::subview(rhs, 0, Kokkos::ALL());
       auto rhs_host_at_0 = Kokkos::create_mirror_view(rhs_at_0);
       Kokkos::deep_copy(rhs_host_at_0, rhs_at_0);
       FILE *fout = fopen(outputFile.c_str(), "w");
       TChem::Test::writeReactionRates(outputFile, number_of_equations, rhs_host_at_0);
       fclose(fout);
     } else {
      auto rhs_host= Kokkos::create_mirror_view(rhs);
      Kokkos::deep_copy(rhs_host, rhs);
      FILE *fout = fopen(outputFile.c_str(), "w");
      TChem::Test::write2DMatrix(outputFile, nBatch, number_of_equations,rhs_host);
      fclose(fout);
     }

    }

    fprintf(fout_times, "}\n ");// end index time
    fclose(fout_times);




  }
  Kokkos::finalize();

  return 0;
}
