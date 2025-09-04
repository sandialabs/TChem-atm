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
#include "TChem_Driver.hpp"
#include "TChem.hpp"
#include "TChem_KineticModelNCAR_ConstData.hpp"
#include "TChem_CommandLineParser.hpp"

using real_type = TChem::real_type;
using ordinal_type = TChem::ordinal_type;

using exec_space = Kokkos::DefaultExecutionSpace;
using host_exec_space = Kokkos::DefaultHostExecutionSpace;

using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
using device_type = typename Tines::UseThisDevice<exec_space>::type;

static TChem::Driver *g_tchem = nullptr;

static int check_flag(const int flag, const std::string funcname)
{
  if (flag < 0) {
    std::cerr << "ERROR: " << funcname << " returned " << flag << std::endl;
    return 1;
  }
  return 0;
}
/**
 * Output statistics from CVODE solver.
 */
void print_cvode_statistics(void* cvode_mem) {
  long int nst, nfe, nsetups, nje, nni, ncfn, netf;
  int retval;

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

  std::cout << "\nFinal Statistics:\n"
            << "  Steps            = " << nst << "\n"
            << "  RHS evals        = " << nfe << "\n"
            << "  LS setups        = " << nsetups << "\n"
            << "  Jac evals        = " << nje << "\n"
            << "  NLS iters        = " << nni << "\n"
            << "  NLS fails        = " << ncfn << "\n"
            << "  Error test fails = " << netf << "\n";
}

/**
 * Initialize the TChem model given input YAML files.
 */
void initialize(const char* chemFile, const char* aeroFile, const char* numericsFile,
  const ordinal_type nBatch){

  g_tchem = new TChem::Driver();

  Kokkos::InitializationSettings settings;
  settings.set_device_id(0);
  Kokkos::initialize(settings);

  using exec_space = Kokkos::DefaultExecutionSpace;
  using host_exec_space = Kokkos::DefaultHostExecutionSpace;

  const bool detail = false;
  TChem::exec_space().print_configuration(std::cout, detail);
  TChem::host_exec_space().print_configuration(std::cout, detail);

  using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
  using device_type = typename Tines::UseThisDevice<exec_space>::type;

  using time_integrator_cvode_type = Tines::TimeIntegratorCVODE<real_type,host_device_type>;
  Tines::value_type_1d_view<time_integrator_cvode_type,host_device_type> cvodes;
  using real_type_2d_view_host_type = Tines::value_type_2d_view<real_type,host_device_type>;
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
  using real_type_2d_view_type = Tines::value_type_2d_view<real_type,device_type>;

  using real_type_1d_view_host_type = Tines::value_type_1d_view<real_type,host_device_type>;
  using real_type_2d_view_host_type = Tines::value_type_2d_view<real_type,host_device_type>;

  using kinetic_model_type = TChem::KineticModelNCAR_ConstData<device_type>;
  using kinetic_model_host_type = TChem::KineticModelNCAR_ConstData<host_device_type>;

  using ordinal_type = TChem::ordinal_type;

  g_tchem->setBatchSize(nBatch);
  g_tchem->createGasKineticModel(chemFile);
  g_tchem->createGasKineticModelConstData();
  g_tchem->createAerosolModel(aeroFile);
  g_tchem->createAerosolModelConstData();
  g_tchem->createStateVector(nBatch);
  g_tchem->getLengthOfStateVector();
  g_tchem->createNumerics(numericsFile);
  g_tchem->createNumberConcentrationVector(nBatch);

}

/**
 * Finalize the model by freeing memory and finalizing Kokkos.
 */
void finalize(){
  g_tchem->freeAll();
  delete g_tchem;
  Kokkos::finalize();
}

/**
 * Set the batch size of the problem.
 */
void TChem::Driver::setBatchSize(const ordinal_type nBatch) {
   _nBatch = nBatch;
}

/**
 * Read in the solver settings from a YAML file.
 */
void TChem::Driver::createNumerics(const std::string &numerics_file) {

  YAML::Node root = YAML::LoadFile(numerics_file);
  YAML::Node solver_info = root["solver_info"];

  auto atol_newton = solver_info["atol_newton"];
  auto rtol_newton = solver_info["rtol_newton"];
  auto dtmin = solver_info["dtmin"];
  auto atol_time = solver_info["atol_time"];
  auto rtol_time = solver_info["rtol_time"];
  auto max_num_newton_iterations = solver_info["max_newton_iterations"];
  auto max_num_time_iterations = solver_info["max_num_time_iterations"];

  auto team_size = solver_info["team_size"];
  auto vector_size = solver_info["vector_size"];

  auto verbose = solver_info["verbose"];

  _atol_newton = atol_newton.as<real_type>(1e-10);
  _rtol_newton = rtol_newton.as<real_type>(1e-6);
  _dtmin = dtmin.as<real_type>(1e-8);
  _atol_time = atol_time.as<real_type>(1e-12);
  _rtol_time = rtol_time.as<real_type>(1e-4);
  _max_num_newton_iterations = max_num_newton_iterations.as<ordinal_type>(100);
  _max_num_time_iterations = max_num_time_iterations.as<ordinal_type>(1e3);

  // If team_size and vector_size are not specified, default to -1
  _team_size = team_size.as<ordinal_type>(-1);
  _vector_size = vector_size.as<ordinal_type>(-1);

  _verbose = verbose.as<bool>(false);

}

/**
 * Create the gas kinetic model from a YAML file.
 */
void TChem::Driver::createGasKineticModel(const std::string &chem_file) {
  _chem_file = chem_file;
  _kmd = KineticModelData(_chem_file);
}

/**
 * Create the gas kinetic constant model data.
 */
void TChem::Driver::createGasKineticModelConstData() {
  printf("Creating kmcd \n");
  using interf_host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
  _kmcd_host = TChem::createNCAR_KineticModelConstData<interf_host_device_type>(_kmd);
  _kmcd_device = TChem::createNCAR_KineticModelConstData<device_type>(_kmd);
  printf("Number of Species %d \n", _kmcd_host.nSpec);
  printf("Number of Reactions %d \n", _kmcd_host.nReac);
  printf("End creating kmcd \n");

}

/**
 * Free model.
 */
void TChem::Driver::freeAll() {
  g_tchem->freeGasKineticModel();
  g_tchem->freeAerosolModel();
}

void TChem::Driver::freeGasKineticModel() {
  _chem_file = std::string();
}

void TChem::Driver::freeAerosolModel() {
  _aero_file = std::string();
}

/**
 * Create the aerosol model from a YAML file.
 */
void TChem::Driver::createAerosolModel(const std::string &aero_file) {
  _aero_file = aero_file;
  _amd = AerosolModelData(_aero_file, _kmd);
}

/**
 * Create the aerosol constant model data
 */
void TChem::Driver::createAerosolModelConstData() {
  printf("Creating amcd \n");
  using interf_host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
  _amcd_host = TChem::create_AerosolModelConstData<interf_host_device_type>(_amd);
  _amcd_device = TChem::create_AerosolModelConstData<device_type>(_amd);
  printf("Number of aerosol species %d \n", _amcd_host.nSpec);
  printf("Maximum number of particles %d \n", _amcd_host.nParticles);
  printf("End creating amcd \n");
}

/**
 * Create the state vector.
 */
void TChem::Driver::createStateVector(ordinal_type nBatch) {
  const ordinal_type len = TChem::Impl::getStateVectorSize(_kmcd_host.nSpec + _amcd_host.nSpec * _amcd_host.nParticles);
  _state = real_type_2d_view_host("state dev", nBatch, len);
}

/**
 * Create number concentration vector.
 */
void TChem::Driver::createNumberConcentrationVector(ordinal_type nBatch) {
  const ordinal_type len = _amcd_host.nParticles;
  _number_concentration = real_type_2d_view_host("number_concentration", nBatch, len);
}

/**
 * Set the values of the state vector for a given batch.
 */
void TChem_setNumberConcentrationVector(double *array, const ordinal_type iBatch){
  g_tchem->setNumberConcentrationVector(array, iBatch);
}

/**
 * Set the number concentration vector for a given batch.
 */
void TChem::Driver::setNumberConcentrationVector(double *array, const ordinal_type iBatch) {
  auto len = _amcd_host.nParticles;
  for (ordinal_type k = 0; k < len; k++){
     _number_concentration(iBatch,k) = array[k];
  }
}

/**
 * Get the state vector for a given batch.
 */
auto TChem::Driver::getStateVector(const ordinal_type iBatch) {
  auto state_at_i_batch = Kokkos::subview(_state, iBatch, Kokkos::ALL);
  return state_at_i_batch;
}

/**
 * Get the state vector for a given batch.
 */
void TChem_getStateVector(TChem::real_type *state, const ordinal_type iBatch){
  auto q = g_tchem->getStateVector(iBatch);
  auto len = TChem_getLengthOfStateVector();
  for (ordinal_type k = 0; k < len; k++) {
     state[k] = q[k];
  }
}

/**
 * Set the values of the state vector for a given batch.
 */
void TChem_setStateVector(double *array, const ordinal_type iBatch){
  g_tchem->setStateVector(array, iBatch);
}

/**
 *
 */
void TChem::Driver::setStateVector(double *array, const ordinal_type iBatch) {
  auto len = TChem_getLengthOfStateVector();
  for (ordinal_type k = 0; k < len; k++){
     _state(iBatch,k) = array[k];
  }
}

/**
 * Return species name at a given index.
 */
int TChem_getSpeciesName(int * index, char* result, const std::string::size_type buffer_size){
  std::string specName = g_tchem->getSpeciesName(index);
  specName.copy(result, buffer_size);
  result[specName.length()] = '\0';
  return specName.length();
}

/**
 * Return species name at a given index.
 */
std::string TChem::Driver::getSpeciesName(int *index){
  const auto speciesNamesHost = Kokkos::create_mirror_view(_kmcd_host.speciesNames);
  int k = *index;
  std::string species_name = &_kmcd_host.speciesNames(k,0);
  return species_name;
}

/**
 * Return number of gas species in the mechanism.
 */
ordinal_type TChem_getNumberOfSpecies(){
  return g_tchem->getNumberOfSpecies();
}

ordinal_type TChem::Driver::getNumberOfSpecies() { return _kmcd_host.nSpec;}

/**
 * Return the number of aerosol species in the mechanism.
 */
ordinal_type TChem_getNumberOfAeroSpecies(){
  ordinal_type nSpec = g_tchem->getNumberOfAeroSpecies();
  return nSpec;
}

ordinal_type TChem::Driver::getNumberOfAeroSpecies() { return _amcd_host.nSpec;}

/**
 * Get the density of the aerosol species at a given index.
 */
real_type TChem_getAerosolSpeciesDensity(int *index) {
  auto density = g_tchem->getAerosolSpeciesDensity(index);
  return density;
}

real_type TChem::Driver::getAerosolSpeciesDensity(int *index) {
  return (_amcd_host.aerosol_density(*index));
}

/**
 * Get the molecular weight of the aerosol species at a given index.
 */
real_type TChem_getAerosolSpeciesMW(int *index) {
  auto mw = g_tchem->getAerosolSpeciesMW(index);
  return mw;
}

real_type TChem::Driver::getAerosolSpeciesMW(int *index) {
  return (_amcd_host.molecular_weights(*index));
}

/**
 * Get the hygroscopicity parameter kappa of the aerosol species at a given index
 */
real_type TChem_getAerosolSpeciesKappa(int *index) {
  auto kappa = g_tchem->getAerosolSpeciesKappa(index);
  return kappa;
}

real_type TChem::Driver::getAerosolSpeciesKappa(int *index) {
  return (_amcd_host.aerosol_kappa(*index));
}

/**
 * Return aerosol species name at a given index.
 */
int TChem_getAerosolSpeciesName(int * index, char* result, const std::string::size_type buffer_size){
  std::string specName = g_tchem->getAerosolSpeciesName(index);
  specName.copy(result, buffer_size);
  result[specName.length()] = '\0';
  return specName.length();
}

std::string TChem::Driver::getAerosolSpeciesName(int *index){

  std::map<int, std::string> aero_idx_sp_name;
  for (std::map<std::string, int>::iterator
        i = _amd.aerosol_sp_name_idx_.begin();
        i != _amd.aerosol_sp_name_idx_.end(); ++i)
        aero_idx_sp_name[i->second] = i->first;

  std::string species_name = aero_idx_sp_name[*index];

  return species_name;
}

/**
 * Return length of the state vector.
 */
int TChem_getLengthOfStateVector() {
  return g_tchem == nullptr ? -1 : g_tchem->getLengthOfStateVector();
}

ordinal_type TChem::Driver::getLengthOfStateVector() const {
  return Impl::getStateVectorSize(_kmcd_host.nSpec + _amcd_host.nSpec * _amcd_host.nParticles);
}

/**
 * Return the length of the number concentration vector.
 */
int TChem_getNumberConcentrationVectorSize() {
  return g_tchem == nullptr ? -1 : g_tchem->getNumberConcentrationVectorSize();
}

ordinal_type TChem::Driver::getNumberConcentrationVectorSize() const{
  return _amcd_host.nParticles;
}

/**
 * Integrate a time step of chemistry.
 */
void TChem_doTimestep(const double &del_t){
  g_tchem->doTimestep(del_t);
}

/* Return the state vector from host */
void TChem::Driver::getStateVectorHost(real_type_2d_const_view_host &view) {
  TCHEM_CHECK_ERROR(_state.span() == 0, "State vector should be constructed");
  view = real_type_2d_const_view_host(&_state(0,0), _state.extent(0), _state.extent(1));
}

/* Integrate a time step */
void TChem::Driver::doTimestep(const double del_t){

  using VecType   = sundials::kokkos::Vector<TChem::exec_space>;
  using MatType   = sundials::kokkos::DenseMatrix<TChem::exec_space>;
  using LSType    = sundials::kokkos::DenseLinearSolver<TChem::exec_space>;
  using SizeType  = VecType::size_type;

  const auto exec_space_instance = TChem::exec_space();
  using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
  using interf_host_device_type =
      typename Tines::UseThisDevice<TChem::host_exec_space>::type;
  using problem_type =
      TChem::Impl::AerosolChemistry_Problem<real_type, interf_host_device_type>;
  using policy_type =
      typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;

  policy_type policy(exec_space_instance, _nBatch, Kokkos::AUTO());

  if (_team_size > 0 && _vector_size > 0) {
    policy = policy_type(exec_space_instance,  _nBatch, _team_size, _vector_size);
  } else if (_team_size > 0 && _vector_size < 0) {
    policy = policy_type(exec_space_instance, _nBatch, _team_size);
  }

  auto number_of_equations = problem_type::getNumberOfTimeODEs(_kmcd_host, _amcd_host);

  using real_type_1d_view_type =
      Tines::value_type_1d_view<real_type, device_type>;
  using real_type_2d_view_type =
      Tines::value_type_2d_view<real_type, device_type>;
  using real_type_2d_view_host_type =
      Tines::value_type_2d_view<real_type, host_device_type>;

  // Create UserData
  TChem::UserData udata;

  udata.nbatches = _nBatch;
  real_type_2d_view number_conc("NumberConcentration", _nBatch,
      _amcd_host.nParticles);
  Kokkos::deep_copy(number_conc, _number_concentration);
  udata.num_concentration = number_conc;
  udata.batchSize = number_of_equations;
  udata.kmcd = _kmcd_device;
  udata.amcd = _amcd_device;

  // Create the SUNDIALS context
  sundials::Context sunctx;

  // Set the initial time to be zero.
  const sunrealtype T0 = SUN_RCONST(0.0);

  SizeType length{static_cast<SizeType>(_nBatch * number_of_equations)};
  VecType y{length, sunctx};
  real_type_2d_view_type y2d((y.View()).data(), _nBatch, number_of_equations);

  const ordinal_type n_active_gas_species = _kmcd_host.nSpec - _kmcd_host.nConstSpec;
  real_type_2d_view_type const_tracers("const_tracers", _nBatch, _kmcd_host.nConstSpec);

  real_type_1d_view_type temperature("temperature", _nBatch);
  real_type_1d_view_type pressure("pressure", _nBatch);
  using range_type = Kokkos::pair<ordinal_type, ordinal_type>;
  const ordinal_type level = 1;

  auto stateVecDim = TChem_getLengthOfStateVector();
  real_type_2d_view state_device("StateVector Devices", _nBatch, stateVecDim);
  Kokkos::deep_copy(state_device, _state);

  const ordinal_type n_gas_spec = _kmcd_device.nSpec;
  const ordinal_type n_gas_spec_constant = _kmcd_device.nConstSpec;
  const ordinal_type total_n_species = _kmcd_device.nSpec  +
      _amcd_device.nParticles * _amcd_device.nSpec;

  Kokkos::parallel_for(
    "fill_y", Kokkos::RangePolicy<TChem::exec_space>(0, _nBatch),
    KOKKOS_LAMBDA(const SizeType i) {
      const real_type_1d_view_type state_at_i =
            Kokkos::subview(state_device, i, Kokkos::ALL());
      TChem::Impl::StateVector<real_type_1d_view_type> sv_at_i(total_n_species, state_at_i);
      temperature(i) = sv_at_i.Temperature();
      pressure(i) = sv_at_i.Pressure();
      const real_type_1d_view_type Ys = sv_at_i.MassFractions();

      const auto activeYs = Kokkos::subview(Ys, range_type(0, n_active_gas_species));
      const auto constYs  = Kokkos::subview(Ys, range_type(n_active_gas_species, n_gas_spec));
      const real_type_1d_view_type partYs = Kokkos::subview(Ys,
            range_type(n_gas_spec, total_n_species));

      for (ordinal_type j = 0; j < n_active_gas_species; ++j){
        y2d(i, j) = activeYs(j);
      }

      for (ordinal_type j = n_active_gas_species; j < total_n_species - n_gas_spec_constant; ++j)
      {
        y2d(i, j) = partYs(j - n_active_gas_species);
      }

      for (ordinal_type j = 0; j < n_gas_spec_constant; ++j){
        const_tracers(i, j) = constYs(j);
      }
  });

  udata.temperature = temperature;
  udata.pressure = pressure;
  udata.const_tracers = const_tracers;

  // Create vector of absolute tolerances
  VecType abstol{length, sunctx};
  N_VConst(SUN_RCONST(_atol_time), abstol);

  // Create CVODE using Backward Differentiation Formula methods
  void* cvode_mem = CVodeCreate(CV_BDF, sunctx);

  // Initialize the integrator and set the ODE right-hand side function
  int retval = CVodeInit(cvode_mem, TChem::AerosolChemistry_CVODE_K::f, T0, y);

  retval = CVodeSetMaxNumSteps(cvode_mem, _max_num_time_iterations);
  retval = CVodeSetMinStep(cvode_mem, _dtmin);
  // Attach the user data structure
  retval = CVodeSetUserData(cvode_mem, &udata);

  // Specify the scalar relative tolerance and vector absolute tolerances
  retval = CVodeSVtolerances(cvode_mem, SUN_RCONST(_rtol_time), abstol);

   // Create the matrix and linear solver objects
  std::unique_ptr<sundials::ConvertibleTo<SUNMatrix>> A;
  std::unique_ptr<sundials::ConvertibleTo<SUNLinearSolver>> LS;

  ordinal_type per_team_extent = 0;

  // Create matrix-free GMRES linear solver
  LS = std::make_unique<sundials::experimental::SUNLinearSolverView>(
      SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, sunctx));

  // Attach the linear solver to CVODE
  retval = CVodeSetLinearSolver(cvode_mem, LS->Convert(), nullptr);

  per_team_extent
       = TChem::Impl::Aerosol_RHS<real_type, device_type>::getWorkSpaceSize(_kmcd_device, _amcd_device);

  const ordinal_type per_team_scratch =
    TChem::Scratch<real_type_1d_view_type>::shmem_size(per_team_extent);
  policy.set_scratch_size(level, Kokkos::PerTeam(per_team_scratch));
  udata.policy = policy;

  // Final time and time between outputs
  const sunrealtype Tf    = SUN_RCONST(del_t);
  const sunrealtype dTout = SUN_RCONST(del_t);

  // Number of output times
  const int Nt_p = static_cast<int>(ceil(Tf / dTout));

  const int Nt = _max_num_time_iterations > 0 ? _max_num_time_iterations : Nt_p;

  // Current time and first output time
  sunrealtype t    = T0;
  sunrealtype tout = T0 + dTout;

  // Initial output
  real_type_2d_view_host_type y2d_h((y.HostView()).data(), udata.nbatches, udata.batchSize);

  // Time stepping
  int iout = 0;
  for (iout = 0; iout < Nt; iout++)
  {
    // Advance in time
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    exec_space_instance.fence();
    if (check_flag(retval, "CVode")) { break; }

    sundials::kokkos::CopyFromDevice(y);
    Kokkos::fence();

    tout += dTout;
    tout = (tout > Tf) ? Tf : tout;
  }

  // Copy results back
  int i = 0;
  for (ordinal_type j = 0; j < n_active_gas_species; ++j){
    _state(i, j + 3) = y2d_h(i, j);
  }

  for (ordinal_type j = n_active_gas_species; j < total_n_species - _kmcd_host.nConstSpec; ++j){
    _state(i, j + 3 + _kmcd_host.nConstSpec) = y2d_h(i, j);
  }

  if (_verbose) {
    print_cvode_statistics(cvode_mem);
  };

  CVodeFree(&cvode_mem);

}
