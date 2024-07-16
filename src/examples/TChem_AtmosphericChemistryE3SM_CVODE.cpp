/* =====================================================================================
TChem_atm version 2.0
Copyright (2024) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC
(NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
Government retains certain rights in this software.

This file is part of TChem-atm_atm. TChem_atm is open source software: you can
redistribute it and/or modify it under the terms of BSD 2-Clause License
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
using time_advance_type = TChem::time_advance_type;

using real_type_0d_view = TChem::real_type_0d_view;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_2d_view = TChem::real_type_2d_view;

using time_advance_type_0d_view = TChem::time_advance_type_0d_view;
using time_advance_type_1d_view = TChem::time_advance_type_1d_view;

using real_type_0d_view_host = TChem::real_type_0d_view_host;
using real_type_1d_view_host = TChem::real_type_1d_view_host;
using real_type_2d_view_host = TChem::real_type_2d_view_host;

using time_advance_type_0d_view_host = TChem::time_advance_type_0d_view_host;
using time_advance_type_1d_view_host = TChem::time_advance_type_1d_view_host;

int main(int argc, char *argv[]) {
// FIXME: remove TPL YAML
// FIXME: use consistent name format
// FIXME: I will keep tchem name space for now; but I want to change tchem to
// tchem_atm  (or other name that we decided to keep)
#if defined(TCHEM_ATM_ENABLE_TPL_YAML_CPP) && defined(TINES_ENABLE_TPL_SUNDIALS)

  /// default inputs
  std::string prefixPath("");

  const real_type zero(0);
  real_type tbeg(0), tend(1);
  real_type dtmin(1e-8), dtmax(1e-1);
  real_type rtol_time(1e-4), atol_newton(1e-10), rtol_newton(1e-6);
  int num_time_iterations_per_interval(1e1), max_num_time_iterations(1e3),
      max_num_newton_iterations(100), jacobian_interval(1);
  real_type atol_time(1e-12);
  int nBatch(1), team_size(-1), vector_size(-1);
  std::string outputFileTimes("wall_times.json");
  bool verbose(true);
  std::string chemFile("chem.yaml");
  std::string outputFile("AtmosphericChemistryE3SM.dat");
  std::string inputFile("None");
  bool use_cloned_samples(false);
  bool use_cvode(false);
  /// parse command line arguments
  TChem::CommandLineParser opts(
      "This example computes the solution of gas chemistry problem");
  opts.set_option<std::string>(
      "inputsPath", "path to input files e.g., data/inputs", &prefixPath);
  opts.set_option<std::string>("chemfile", "Chem file name e.g., chem.inp",
                               &chemFile);
  opts.set_option<std::string>("inputfile", "Chem file name e.g., chem.inp",
                               &inputFile);
  opts.set_option<real_type>("tbeg", "Time begin", &tbeg);
  opts.set_option<real_type>("tend", "Time end", &tend);
  opts.set_option<real_type>("dtmin", "Minimum time step size", &dtmin);
  opts.set_option<real_type>("dtmax", "Maximum time step size", &dtmax);
  opts.set_option<real_type>(
      "atol-newton", "Absolute tolerance used in newton solver", &atol_newton);
  opts.set_option<real_type>(
      "rtol-newton", "Relative tolerance used in newton solver", &rtol_newton);
  opts.set_option<std::string>(
      "outputfile", "Output file name e.g., IgnSolution.dat", &outputFile);
  opts.set_option<real_type>(
      "tol-time", "Tolerence used for adaptive time stepping", &rtol_time);
  opts.set_option<int>("time-iterations-per-interval",
                       "Number of time iterations per interval to store qoi",
                       &num_time_iterations_per_interval);
  opts.set_option<std::string>("outputfile_times",
                               "Wal times file name e.g., times.json",
                               &outputFileTimes);
  opts.set_option<int>("max-time-iterations",
                       "Maximum number of time iterations",
                       &max_num_time_iterations);
  opts.set_option<int>("jacobian-interval",
                       "Jacobians are evaluated once in this interval",
                       &jacobian_interval);
  opts.set_option<int>("max-newton-iterations",
                       "Maximum number of newton iterations",
                       &max_num_newton_iterations);
  opts.set_option<bool>(
      "verbose", "If true, printout the first Jacobian values", &verbose);
  opts.set_option<int>("team_thread_size", "time thread size ", &team_size);
  opts.set_option<int>("batch_size", " number of batches or samples, e.g. 10  ", &nBatch);
  opts.set_option<int>("vector_thread_size", "vector thread size ",
                       &vector_size);
  opts.set_option<real_type>("atol-time", "Absolute tolerance used for adaptive time stepping", &atol_time);
  opts.set_option<bool>(
      "use_cloned_samples", "If true, one state vector will be cloned.", &use_cloned_samples);
  opts.set_option<bool>("use-cvode",
      "if true code runs ignition zero d reactor with cvode; otherwise, it uses TrBDF2",
      &use_cvode);

  const bool r_parse = opts.parse(argc, argv);
  if (r_parse)
    return 0; // print help return

  Kokkos::initialize(argc, argv);
  {
    const bool detail = false;

    // by default use input condition from chemFile
    if (inputFile == "None") {
       inputFile=chemFile;
    }


    printf("   TChem is running ATM model box with CVODE\n");

    TChem::exec_space().print_configuration(std::cout, detail);
    TChem::host_exec_space().print_configuration(std::cout, detail);

    using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;

    /// construct kmd and use the view for testing
    TChem::KineticModelData kmd(chemFile);
    const auto kmcd = TChem::createNCAR_KineticModelConstData<host_device_type>(kmd);

    const ordinal_type stateVecDim =
        TChem::Impl::getStateVectorSize(kmcd.nSpec);

    printf("Number of Species %d \n", kmcd.nSpec);
    printf("Number of Reactions %d \n", kmcd.nReac);

    const auto speciesNamesHost = Kokkos::create_mirror_view(kmcd.speciesNames);
    Kokkos::deep_copy(speciesNamesHost, kmcd.speciesNames);

    FILE *fout = fopen(outputFile.c_str(), "w");
    // read scenario condition from yaml file
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
    // read external forcing
    ordinal_type count_ext_forcing=0;

    const ordinal_type n_active_vars = kmcd.nSpec - kmcd.nConstSpec;
    printf("Number of species %d \n", kmcd.nSpec);
    printf("Number of const species %d \n", kmcd.nConstSpec);
    // Note: We allocate external_sources_scenario_host
    real_type_2d_view_host external_sources_scenario_host("external_sources_host",nbacth_files,n_active_vars);
    TChem::AtmChemistry::setScenarioConditionsExternalForcing(inputFile,
             speciesNamesHost,
             // output
             external_sources_scenario_host,
             count_ext_forcing);

    real_type_2d_view_host state;
    real_type_2d_view_host photo_rates;
    real_type_2d_view_host external_sources;

    if (nbacth_files == 1 && use_cloned_samples && nBatch > 1) {
      // only clone samples if nbacth_files is 1
      printf("-------------------------------------------------------\n");
      printf("--------------------Warning----------------------------\n");
      printf("Using cloned samples ... only for numerical experiments\n");
      state = real_type_2d_view_host("StateVector Devices", nBatch, stateVecDim);
      auto state_scenario_host_at_0 = Kokkos::subview(state_scenario_host, 0, Kokkos::ALL);
      auto state_host_at_0 = Kokkos::subview(state, 0, Kokkos::ALL);
      Kokkos::deep_copy(state_host_at_0, state_scenario_host_at_0);
      TChem::Test::cloneView(state);

      if (n_photo_rates > 0 )
      {
        photo_rates = real_type_2d_view_host("photo rates", nBatch, n_photo_rates);
        auto photo_rates_scenario_host_at_0 = Kokkos::subview(photo_rates_scenario_host, 0, Kokkos::ALL);
        auto photo_rates_host_at_0 = Kokkos::subview(photo_rates, 0, Kokkos::ALL);
        Kokkos::deep_copy(photo_rates_host_at_0, photo_rates_scenario_host_at_0);
        TChem::Test::cloneView(photo_rates);
      } // n_photo_rates

      external_sources = real_type_2d_view("external_sources", nBatch, n_active_vars);
      if (count_ext_forcing >  0) {
        auto external_sources_scenario_host_at_0 = Kokkos::subview(external_sources_scenario_host, 0, Kokkos::ALL);
        auto external_sources_host_at_0 = Kokkos::subview(external_sources, 0, Kokkos::ALL);
        Kokkos::deep_copy(external_sources_host_at_0, external_sources_scenario_host_at_0);
        TChem::Test::cloneView(external_sources);
      }else {
        // make sure that external source are zero.
       Kokkos::deep_copy(external_sources,zero);
      }

    } else {
      nBatch = nbacth_files;
      state = state_scenario_host;

      if (n_photo_rates > 0 )
      {
        photo_rates = photo_rates_scenario_host;
      }  // n_photo_rates


      if (count_ext_forcing >  0) {
        external_sources= external_sources_scenario_host;
      } else {
        // make sure that external source are zero.
        external_sources = real_type_2d_view_host("external_sources", nBatch, n_active_vars);
       Kokkos::deep_copy(external_sources,zero);
      } //   count_ext_forcing

    }
    auto writeState = [](const ordinal_type iter,
                         const real_type_1d_view_host _t,
                         const real_type_1d_view_host _dt,
                         const real_type_2d_view_host _state_at_i,
                         FILE *fout) { // sample, time, density, pressure,
                                       // temperature, mass fraction
      for (size_t sp = 0; sp < _state_at_i.extent(0); sp++) {
        fprintf(fout, "%d \t %15.10e \t  %15.10e \t ", iter, _t(sp), _dt(sp));
        //
        for (ordinal_type k = 0, kend = _state_at_i.extent(1); k < kend; ++k)
          fprintf(fout, "%15.10e \t", _state_at_i(sp, k));

        fprintf(fout, "\n");
      }

    };

    auto printState = [](const time_advance_type _tadv, const real_type _t,
                         const real_type_1d_view_host _state_at_i) {
      /// iter, t, dt, rho, pres, temp, Ys ....
      printf("%e %e %e %e %e", _t, _t - _tadv._tbeg, _state_at_i(0),
             _state_at_i(1), _state_at_i(2));
      for (ordinal_type k = 3, kend = _state_at_i.extent(0); k < kend; ++k)
        printf(" %e", _state_at_i(k));
      printf("\n");
    };

    FILE *fout_times = fopen(outputFileTimes.c_str(), "w");
    fprintf(fout_times, "{\n");
    fprintf(fout_times, " \"Atmospheric Chemistry E3SM\": \n {\n");

    {
      const auto exec_space_instance = TChem::exec_space();

      using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;

      /// team policy
      policy_type policy(exec_space_instance, nBatch, Kokkos::AUTO());

      // if (team_size > 0 && vector_size > 0) {
      //   policy = policy_type(exec_space_instance,  nBatch, team_size, vector_size);
      // } else if (team_size > 0 && vector_size < 0) {
      // // only set team size
      //  policy = policy_type(exec_space_instance, nBatch,  team_size);
      // }
          using time_integrator_cvode_type = Tines::TimeIntegratorCVODE<real_type,host_device_type>;
    Tines::value_type_1d_view<time_integrator_cvode_type,host_device_type> cvodes;

    if (use_cvode) {
      cvodes = Tines::value_type_1d_view<time_integrator_cvode_type,host_device_type>("cvodes", nBatch);
      for (ordinal_type i=0;i<nBatch;++i)
        cvodes(i).create(n_active_vars);
    }

        ordinal_type number_of_equations(0);

        using problem_type =
            TChem::Impl::AtmosphericChemistryE3SM_Problem<real_type,
                                                          host_device_type>;
        number_of_equations = problem_type::getNumberOfTimeODEs(kmcd);


      const ordinal_type level = 1;
      ordinal_type per_team_extent(0);

      per_team_extent = TChem::AtmosphericChemistryE3SM_CVODE::getWorkSpaceSize(kmcd);

      const ordinal_type per_team_scratch =
          TChem::Scratch<real_type_1d_view_host>::shmem_size(per_team_extent);
      policy.set_scratch_size(level, Kokkos::PerTeam(per_team_scratch));

      { /// time integration
        real_type_1d_view_host t("time", nBatch);
        Kokkos::deep_copy(t, tbeg);
        real_type_1d_view_host dt("delta time", nBatch);
        Kokkos::deep_copy(dt, dtmax);

        real_type_2d_view_host tol_time("tol time", number_of_equations, 2);
        real_type_1d_view_host tol_newton("tol newton", 2);

        real_type_2d_view_host fac("fac", nBatch, number_of_equations);



        /// tune tolerence
        {
          for (ordinal_type i = 0, iend = tol_time.extent(0); i < iend; ++i) {
            tol_time(i, 0) = atol_time;
            tol_time(i, 1) = rtol_time;
          }
          tol_newton(0) = atol_newton;
          tol_newton(1) = rtol_newton;
        }

        time_advance_type tadv_default;
        tadv_default._tbeg = tbeg;
        tadv_default._tend = tend;
        tadv_default._dt = dtmin;
        tadv_default._dtmin = dtmin;
        tadv_default._dtmax = dtmax;
        tadv_default._max_num_newton_iterations = max_num_newton_iterations;
        tadv_default._num_time_iterations_per_interval =
            num_time_iterations_per_interval;
        tadv_default._jacobian_interval = jacobian_interval;

        time_advance_type_1d_view_host tadv("tadv", nBatch);
        Kokkos::deep_copy(tadv, tadv_default);

        ordinal_type iter = 0;

        fprintf(fout, "%s \t %s \t %s \t ", "iter", "t", "dt");
        fprintf(fout, "%s \t %s \t %s \t", "Density[kg/m3]", "Pressure[Pascal]",
                "Temperature[K]");

        for (ordinal_type k = 0; k < kmcd.nSpec; k++)
          fprintf(fout, "%s \t", &speciesNamesHost(k, 0));
        fprintf(fout, "\n");
        writeState(-1, t, dt, state, fout);

        auto kmds = kmd.clone(nBatch);
        auto kmcds = TChem::createNCAR_KineticModelConstData<host_device_type>(kmds);


        real_type tsum(0);
        Kokkos::Timer timer;

        for (; iter < max_num_time_iterations && tsum <= tend * 0.9999;
             ++iter) {


      timer.reset();
       TChem::AtmosphericChemistryE3SM_CVODE::runHostBatch(
              policy, tol_time, fac, tadv, state, photo_rates, external_sources, t, dt, state,
              kmcds, cvodes);
      exec_space_instance.fence();

          const real_type t_device_batch = timer.seconds();
          fprintf(fout_times, "\"%s%d\": %20.14e, \n", "wall_time_iter_", iter,
                  t_device_batch);

          writeState(iter, t, dt, state, fout);

          /// carry over time and dt computed in this step
          tsum = zero;
          Kokkos::parallel_reduce(
              Kokkos::RangePolicy<TChem::host_exec_space>(0, nBatch),
              KOKKOS_LAMBDA(const ordinal_type &i, real_type &update) {
                tadv(i)._tbeg = t(i);
                tadv(i)._dt = dt(i);
                // printf("t %e, dt %e\n", t(i), dt(i));
                // printf("tadv t %e, tadv dt %e\n", tadv(i)._tbeg, tadv(i)._dt
                // );
                update += t(i);
              },
              tsum);
          Kokkos::fence();
          tsum /= nBatch;
        }// end iter
        fprintf(fout_times, "%s: %d, \n", "\"number_of_time_iters\"", iter);

      }
    }

    fprintf(fout_times, "%s: %d \n", "\"number_of_samples\"", nBatch);
    fprintf(fout_times, "} \n "); // reaction rates


    fclose(fout);
    fprintf(fout_times, "}\n "); // end index time
    fclose(fout_times);
  }
  Kokkos::finalize();

#else
  printf("This example requires Yaml ...\n");
#endif

  return 0;
}
