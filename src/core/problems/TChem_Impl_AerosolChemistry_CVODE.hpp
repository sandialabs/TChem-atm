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
#ifndef __TCHEM_IMPL_AEROSOL_CHEMISTRY_CVODE_HPP__
#define __TCHEM_IMPL_AEROSOL_CHEMISTRY_CVODE_HPP__

#include "Tines_Internal.hpp"

#include "TChem_KineticModelData.hpp"
#include "TChem_AerosolModelData.hpp"
#include "TChem_Impl_AerosolChemistry_Problem.hpp"


namespace TChem {
namespace Impl {

#if defined(TINES_ENABLE_TPL_SUNDIALS)
#include "Tines_Interface.hpp"

    static int ProblemAerosolChemistry_ComputeFunctionCVODE(sunrealtype t,
               N_Vector u,
               N_Vector f,
               void *user_data) {
      using host_device_type = Tines::UseThisDevice<Kokkos::Serial>::type;
      using problem_type = AerosolChemistry_Problem<sunrealtype,host_device_type>;
      using realtype_1d_view_type = Tines::value_type_1d_view<sunrealtype, host_device_type>;

      problem_type * problem = (problem_type*)(user_data);
      TINES_CHECK_ERROR(problem == nullptr, "user data is failed to cast to problem type");

      int m = problem->getNumberOfEquations();
      const auto member = Tines::HostSerialTeamMember();

      sunrealtype * u_data = N_VGetArrayPointer_Serial(u);
      sunrealtype * f_data = N_VGetArrayPointer_Serial(f);

      realtype_1d_view_type uu(u_data, m);
      realtype_1d_view_type ff(f_data, m);

      problem->computeFunction(member, uu, ff);
      return 0;
    }

    static int ProblemAerosolChemistry_ComputeJacobianCVODE(sunrealtype t,
               N_Vector u,
               N_Vector f,
               SUNMatrix J,
               void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
      using host_device_type = Tines::UseThisDevice<Kokkos::Serial>::type;
      using problem_type = AerosolChemistry_Problem<sunrealtype,host_device_type>;
      using realtype_1d_view_type = Tines::value_type_1d_view<sunrealtype, host_device_type>;
      using realtype_2d_view_type = Tines::value_type_2d_view<sunrealtype, host_device_type>;

      problem_type * problem = (problem_type*)(user_data);;
      TINES_CHECK_ERROR(problem == nullptr, "user data is failed to cast to problem type");

      int m = problem->getNumberOfEquations();
      const auto member = Tines::HostSerialTeamMember();

      sunrealtype * u_data = N_VGetArrayPointer_Serial(u);

      realtype_1d_view_type uu(u_data, m);
      realtype_2d_view_type JJ(problem->_work_cvode.data(), m, m);
      problem->computeNumericalJacobian(member, uu, JJ);

      /// prevent potential mismatch of data layout between kokkos and sundials
      for (int i=0;i<m;++i)
  for (int j=0;j<m;++j)
    SM_ELEMENT_D(J, i, j) = JJ(i,j);
      return 0;
    }
#endif

} // Impl
}// TChem

#endif
