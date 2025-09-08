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

using ordinal_type = TChem::ordinal_type;
using real_type = TChem::real_type;
// User data structure available in user-supplied callback functions
using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using problem_type = TChem::Impl::AerosolChemistry_Problem<real_type,device_type>;
using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
using real_type_0d_view_type = Tines::value_type_0d_view<real_type, device_type>;
using real_type_1d_view_type = Tines::value_type_1d_view<real_type, device_type>;
using real_type_2d_view_type = Tines::value_type_2d_view<real_type, device_type>;
using real_type_3d_view_type = Tines::value_type_3d_view<real_type, device_type>;
using real_type_1d_view_host_type = Tines::value_type_1d_view<real_type, host_device_type>;
using real_type_2d_view_host_type = Tines::value_type_2d_view<real_type, host_device_type>;
using ordinal_type_type_1d_view_type = Tines::value_type_1d_view<ordinal_type, device_type>;

using VecType   = sundials::kokkos::Vector<TChem::exec_space>;
using MatType   = sundials::kokkos::DenseMatrix<TChem::exec_space>;
using LSType    = sundials::kokkos::DenseLinearSolver<TChem::exec_space>;
using SizeType  = VecType::size_type;

// Check for an unrecoverable (negative) return value from a SUNDIALS function
static int check_flag(const int flag, const std::string funcname)
{
  if (flag < 0)
  {
    std::cerr << "ERROR: " << funcname << " returned " << flag << std::endl;
    return 1;
  }
  return 0;
}

// Check if a function returned a NULL pointer
static int check_ptr(const void* ptr, const std::string funcname)
{
  if (ptr) { return 0; }
  std::cerr << "ERROR: " << funcname << " returned NULL" << std::endl;
  return 1;
}



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
  int number_of_particles(-1), max_num_time_iterations(-1);
  real_type rtol_time(1e-4), atol_time(1e-12);
  real_type tbeg(0), tend(1);
  real_type dtmin(1e-8);
  bool write_time_profiles(true);
  // , dtmax(1e-1);
  // Linear solver type
  int solver_type = 1;



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
  opts.set_option<real_type>(
      "rtol-time", "Tolerence used for adaptive time stepping", &rtol_time);
  opts.set_option<real_type>(
      "atol-time", "Tolerence used for adaptive time stepping", &atol_time);
  opts.set_option<real_type>("tbeg", "Time begin", &tbeg);
  opts.set_option<real_type>("tend", "Time end", &tend);
  opts.set_option<real_type>("dtmin", "Minimum time step size", &dtmin);
  opts.set_option<int>("solver_type", "solver type. ", &solver_type);
  opts.set_option<bool>(
      "write-time-profiles", "If true, this example will write the time profile output to a file.", &write_time_profiles);
  opts.set_option<int>("max-time-iterations",
                       "Maximum number of time iterations",
                       &max_num_time_iterations);
  const bool r_parse = opts.parse(argc, argv);
  if (r_parse)
    return 0; // print help return
  // Create the SUNDIALS context
  sundials::Context sunctx;
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
#if 1
    FILE *fout;
    if (write_time_profiles) {
     fout = fopen(outputFile.c_str(), "w");
    }
    FILE *fout_times = fopen(outputFileTimes.c_str(), "w");
    fprintf(fout_times, "{\n");
    fprintf(fout_times, " \"Aerosol Chemistry\": \n {\n");

    auto writeState = [](const ordinal_type iter,
                         const real_type _t,
                         const real_type _dt,
                         const real_type_1d_view_host_type density,
                         const real_type_1d_view_host_type pressure,
                         const real_type_1d_view_host_type temperature,
                         const real_type_2d_view_host_type _const_tracers_at_i,
                         const real_type_2d_view_host_type _sol_at_i,
                         const ordinal_type n_active_species,
                         FILE *fout) { // sample, time, density, pressure,
                                       // temperature, mass fraction

      for (size_t sp = 0; sp < _sol_at_i.extent(0); sp++) {
        fprintf(fout, "%d \t %15.10e \t  %15.10e \t ", iter, _t, _dt);
        //
        fprintf(fout, "%15.10e \t", density(sp));
        fprintf(fout, "%15.10e \t", pressure(sp));
        fprintf(fout, "%15.10e \t", temperature(sp));
        for (ordinal_type k = 0, kend = n_active_species; k < kend; ++k)
          fprintf(fout, "%15.10e \t", _sol_at_i(sp, k));
        for (ordinal_type k = 0, kend = _const_tracers_at_i.extent(1); k < kend; ++k)
          fprintf(fout, "%15.10e \t", _const_tracers_at_i(sp, k));
        for (ordinal_type k = n_active_species, kend = _sol_at_i.extent(1); k < kend; ++k)
          fprintf(fout, "%15.10e \t", _sol_at_i(sp, k));
        fprintf(fout, "\n");
      }

    };
#endif
  	const bool detail = false;
    // constexpr real_type zero =0;

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
    real_type_2d_view_type fac;

    using policy_type = typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;

    // fprintf(fout_times, " Computation type, total time [sec], time per sample [sec/sample]\n ");

    policy_type policy(exec_space_instance, nBatch, Kokkos::AUTO());

    if (team_size > 0 && vector_size > 0) {
        policy = policy_type(exec_space_instance,  nBatch, team_size, vector_size);
    } else if (team_size > 0 && vector_size < 0) {
      // only set team size
       policy = policy_type(exec_space_instance, nBatch,  team_size, Kokkos::AUTO());
    }

    using range_type = Kokkos::pair<ordinal_type, ordinal_type>;
    const ordinal_type level = 1;
    // Create UserData
    TChem::UserData udata;

    udata.nbatches = nBatch;
    udata.num_concentration = num_concentration;
    udata.batchSize=number_of_equations;
    udata.kmcd=kmcd;
    udata.amcd=amcd;
    // Create vector with the initial condition
    const sunrealtype T0 = SUN_RCONST(tbeg);

    SizeType length{static_cast<SizeType>(nBatch * number_of_equations)};
    VecType y{length, sunctx};
    real_type_2d_view_type y2d((y.View()).data(), nBatch, number_of_equations);

    const ordinal_type n_active_gas_species = kmcd.nSpec - kmcd.nConstSpec;
    real_type_2d_view_type const_tracers("const_tracers",nBatch, kmcd.nConstSpec);

    real_type_1d_view_type temperature("temperature",nBatch);
    real_type_1d_view_type pressure("pressure",nBatch);

    Kokkos::parallel_for(
      "fill_y", Kokkos::RangePolicy<TChem::exec_space>(0, nBatch),
      KOKKOS_LAMBDA(const SizeType i) {
        const real_type_1d_view_type state_at_i =
        Kokkos::subview(state, i, Kokkos::ALL());
        const ordinal_type total_n_species = kmcd.nSpec + amcd.nParticles*amcd.nSpec;
        TChem::Impl::StateVector<real_type_1d_view_type> sv_at_i(total_n_species, state_at_i);
        temperature(i) = sv_at_i.Temperature();
        pressure(i) = sv_at_i.Pressure();
        const real_type_1d_view_type Ys = sv_at_i.MassFractions();

        const auto activeYs = Kokkos::subview(Ys, range_type(0, n_active_gas_species));
        const auto constYs  = Kokkos::subview(Ys, range_type(n_active_gas_species, kmcd.nSpec));
        const real_type_1d_view_type partYs = Kokkos::subview(Ys, range_type(kmcd.nSpec, total_n_species));

        for (ordinal_type j=0;j<n_active_gas_species;++j){
          y2d(i, j) = activeYs(j);
        }

        for (ordinal_type j=n_active_gas_species;j<total_n_species- kmcd.nConstSpec;++j)
        {
          y2d(i, j) = partYs(j-n_active_gas_species);
        }

        for (ordinal_type j=0;j<kmcd.nConstSpec;++j){
          const_tracers(i, j) = constYs(j);
        }

     });

    udata.temperature = temperature;
    udata.pressure = pressure;
    udata.const_tracers =const_tracers;
    if (write_time_profiles) {
        fprintf(fout, "%s \t %s \t %s \t ", "iter", "t", "dt");
        fprintf(fout, "%s \t %s \t %s \t", "Density[kg/m3]", "Pressure[Pascal]",
                "Temperature[K]");

        for (ordinal_type k = 0; k < kmcd.nSpec; k++)
          fprintf(fout, "%s \t", &speciesNamesHost(k, 0));
        //given aero idx return species name
        std::map<int, std::string> aero_idx_sp_name;
        for (std::map<std::string, int>::iterator
           i = amd.aerosol_sp_name_idx_.begin();
          i != amd.aerosol_sp_name_idx_.end(); ++i)
          aero_idx_sp_name[i->second] = i->first;

        for (ordinal_type ipart = 0; ipart < amd.nParticles_; ipart++)
        {
          for (ordinal_type isp = 0; isp < amd.nSpec_; isp++)
          {
            // std::cout << "species Name : "<< aero_idx_sp_name[i] << "\n";
            auto aero_sp_name = aero_idx_sp_name[isp]+"_p"+std::to_string(ipart);
            fprintf(fout, "%s \t", aero_sp_name.c_str());
          }// isp
        }// ipar
      fprintf(fout, "\n");
    }


    // Create vector of absolute tolerances
    VecType abstol{length, sunctx};
    N_VConst(SUN_RCONST(atol_time), abstol);

    // Create CVODE using Backward Differentiation Formula methods
    void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
    if (check_ptr(cvode_mem, "CVodeCreate")) { return 1; }


    // Initialize the integrator and set the ODE right-hand side function
    int retval = CVodeInit(cvode_mem, TChem::AerosolChemistry_CVODE_K::f, T0, y);
    if (check_flag(retval, "CVodeInit")) { return 1; }


    // Attach the user data structure
    retval = CVodeSetUserData(cvode_mem, &udata);
    if (check_flag(retval, "CVodeSetUserData")) { return 1; }

    // Specify the scalar relative tolerance and vector absolute tolerances
    retval = CVodeSVtolerances(cvode_mem, SUN_RCONST(rtol_time), abstol);
    if (check_flag(retval, "CVodeSVtolerances")) { return 1; }

     // Create the matrix and linear solver objects
    std::unique_ptr<sundials::ConvertibleTo<SUNMatrix>> A;
    std::unique_ptr<sundials::ConvertibleTo<SUNLinearSolver>> LS;

    ordinal_type per_team_extent=0;

    if (solver_type == 0)
    {
      // Create Kokkos dense block diagonal matrix
      A = std::make_unique<MatType>(udata.nbatches, udata.batchSize, udata.batchSize, sunctx);

      // Create Kokkos batched dense linear solver
      LS = std::make_unique<LSType>(sunctx);

      // Attach the matrix and linear solver to CVODE
      retval = CVodeSetLinearSolver(cvode_mem, LS->Convert(), A->Convert());
      if (check_flag(retval, "CVodeSetLinearSolver")) { return 1; }

      // Set the user-supplied Jacobian function
      retval = CVodeSetJacFn(cvode_mem, TChem::AerosolChemistry_CVODE_K::Jac);
      if (check_flag(retval, "CVodeSetJacFn")) { return 1; }
      fac = real_type_2d_view_type("fac", nBatch, number_of_equations);
      udata.fac=fac;
      per_team_extent = problem_type::getWorkSpaceSize(kmcd,amcd)
        + number_of_equations;
#if defined(TCHEM_ATM_ENABLE_GPU)
    real_type_3d_view_type JacRL("JacRL", nBatch, number_of_equations, number_of_equations);
    udata.JacRL=JacRL;
#endif
    }
    else
    {
      // Create matrix-free GMRES linear solver
      LS = std::make_unique<sundials::experimental::SUNLinearSolverView>(
        SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, sunctx));

      // Attach the linear solver to CVODE
      retval = CVodeSetLinearSolver(cvode_mem, LS->Convert(), nullptr);
      if (check_flag(retval, "CVodeSetLinearSolver")) { return 1; }
      per_team_extent
         = TChem::Impl::Aerosol_RHS<real_type, device_type>::getWorkSpaceSize(kmcd, amcd);

    }
    const ordinal_type per_team_scratch =
      TChem::Scratch<real_type_1d_view_type>::shmem_size(per_team_extent);
      policy.set_scratch_size(level, Kokkos::PerTeam(per_team_scratch));
    udata.policy = policy;

    // Final time and time between outputs

    const sunrealtype Tf    = SUN_RCONST(tend);
    const sunrealtype dTout = SUN_RCONST(dtmin);

    // Number of output times
    const int Nt_p = static_cast<int>(ceil(Tf / dTout));

    const int Nt =  max_num_time_iterations > 0 ? max_num_time_iterations: Nt_p;

    // Current time and first output time
    sunrealtype t    = T0;
    sunrealtype tout = T0 + dTout;

    // Initial output
    real_type_2d_view_host_type y2d_h((y.HostView()).data(), udata.nbatches, udata.batchSize);
    sundials::kokkos::CopyFromDevice(y);
    Kokkos::fence();

    const auto density_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace() , temperature);
    const auto temperature_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace() , temperature);
    const auto pressure_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pressure);
    const auto const_tracers_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),const_tracers);
    if (write_time_profiles) {
    writeState(-1, tbeg, dTout,
     density_host, pressure_host, temperature_host,
     const_tracers_host, y2d_h, n_active_gas_species, fout);
    }

    // TChem::TChemAerosolChemistryRHS rhs_tchem(rhs, vals, num_concentration,
    //   const_tracers, temperature, pressure, kmcd, amcd);

    // Loop over output times
    Kokkos::Timer timer;
    int iout = 0;
    for (iout = 0; iout < Nt && tout <= tend * 0.9999; iout++)
    {
      // Advance in time
      timer.reset();
      retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
      exec_space_instance.fence();
      const real_type t_device_batch = timer.seconds();
      fprintf(fout_times, "\"%s%d\": %20.14e, \n", "wall_time_iter_", iout,
                  t_device_batch);
      if (check_flag(retval, "CVode")) { break; }

      // // Output solution from some batches
      sundials::kokkos::CopyFromDevice(y);
      Kokkos::fence();

      if (write_time_profiles) {
      writeState(iout, tout, dTout,
       density_host, pressure_host, temperature_host,
       const_tracers_host, y2d_h,
        n_active_gas_species, fout);
      }
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }

   fprintf(fout_times, "%s: %d, \n", "\"number_of_time_iters\"", iout);
   fprintf(fout_times, "%s: %d \n", "\"number_of_samples\"", nBatch);
   fprintf(fout_times, "} \n "); // reaction rates

   fprintf(fout_times, ", \n ");
   fprintf(fout_times, " \"Kokkos Settings\": \n  {\n");
   fprintf(fout_times, "%s: %d, \n", "\"team_size\"", team_size);
   fprintf(fout_times, "%s: %d, \n", "\"vector_size\"", vector_size);
   fprintf(fout_times, "%s: %d, \n", "\"vector_length_max\"", policy.vector_length_max());
   fprintf(fout_times, "%s: %d, \n", "\"impl_vector_length\"", policy.impl_vector_length());
   fprintf(fout_times, "%s: %d, \n", "\"team_size_recommended_rhs\"", udata.team_size_recommended_rhs);
   fprintf(fout_times, "%s: %d, \n", "\"team_size_recommended_jac\"", udata.team_size_recommended_jac);
   fprintf(fout_times, "%s: %d, \n", "\"team_size_max_rhs\"", udata.team_size_max_rhs);
   fprintf(fout_times, "%s: %d \n", "\"team_size_max_jac\"", udata.team_size_max_jac);
   fprintf(fout_times, "} \n "); // reaction rates
        // Print some final statistics
    long int nst, nfe, nsetups, nje, nni, ncfn, netf;

    retval = CVodeGetNumSteps(cvode_mem, &nst);
    check_flag(retval, "CVodeGetNumSteps");
    retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    check_flag(retval, "CVodeGetNumRhsEvals");
    retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
    check_flag(retval, "CVodeGetNumLinSolvSetups");
    retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
    check_flag(retval, "CVodeGetNumErrTestFails");
    retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    check_flag(retval, "CVodeGetNumNonlinSolvIters");
    retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    check_flag(retval, "CVodeGetNumNonlinSolvConvFails");
    retval = CVodeGetNumJacEvals(cvode_mem, &nje);
    check_flag(retval, "CVodeGetNumJacEvals");

   fprintf(fout_times, ", \n ");
   fprintf(fout_times, " \"Sundials Final Statistics\": \n  {\n");
   fprintf(fout_times, "%s: %d, \n", "\"steps\"", nst);
   fprintf(fout_times, "%s: %d, \n", "\"RHS_evals\"", nfe);
   fprintf(fout_times, "%s: %d, \n", "\"LS_setups\"", nsetups);
   fprintf(fout_times, "%s: %d, \n", "\"Jac_evals\"", nje);
   fprintf(fout_times, "%s: %d, \n", "\"NLS_iters\"", nni);
   fprintf(fout_times, "%s: %d, \n", "\"NLS_fails\"", ncfn);
   fprintf(fout_times, "%s: %d \n", "\"Error_test_fails\"", netf);
   fprintf(fout_times, "} \n "); // reaction rates

    std::cout << "\nFinal Statistics:\n"
              << "  Steps            = " << nst << "\n"
              << "  RHS evals        = " << nfe << "\n"
              << "  LS setups        = " << nsetups << "\n"
              << "  Jac evals        = " << nje << "\n"
              << "  NLS iters        = " << nni << "\n"
              << "  NLS fails        = " << ncfn << "\n"
              << "  Error test fails = " << netf << "\n";

    // Free objects
    CVodeFree(&cvode_mem);


    // fprintf(fout_times, "}\n ");// end index time
    // fclose(fout_times);
    if (write_time_profiles) {
    fclose(fout);
    }
    fprintf(fout_times, "}\n "); // end index time
    fclose(fout_times);




  }
  Kokkos::finalize();

  return 0;
}
