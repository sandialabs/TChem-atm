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
#ifndef __TCHEM_IMPL_AEROSOL_CHEMISTRY_KOKKOSKERNELS_HPP__
#define __TCHEM_IMPL_AEROSOL_CHEMISTRY_KOKKOSKERNELS_HPP__

#include "Tines_Internal.hpp"

#include "TChem_KineticModelData.hpp"
#include "TChem_AerosolModelData.hpp"
#include "TChem_Impl_AerosolChemistry_Problem.hpp"
#include "KokkosODE_BDF.hpp"

namespace TChem {
namespace Impl {

template<typename MemberType, typename ValueType, typename DeviceType>
struct KokkosKernelsODE {  

  using value_type = ValueType;
  using device_type = DeviceType;
  using member_type = MemberType;
  using scalar_type = typename ats<value_type>::scalar_type;
  using problem_type = TChem::Impl::AerosolChemistry_Problem<value_type, device_type>;
  using value_type_1d_view_type = Tines::value_type_1d_view<value_type,device_type>;

  problem_type problem;
  ordinal_type neqs;
  member_type member;

  KOKKOS_FUNCTION
  KokkosKernelsODE(const ordinal_type& neqs_,
                   const problem_type& problem_,
                   const member_type& member_)
      : neqs(neqs_), problem(problem_), member(member_) {}


  KOKKOS_FUNCTION void evaluate_function(const double /*t*/,
                                         const double /*dt*/,
                                         const value_type_1d_view_type& y,
                                         const value_type_1d_view_type& f) const {
  problem.computeFunction(member,
                             y,
                             f);
  }

    template <class vec_type, class mat_type>
  KOKKOS_FUNCTION void evaluate_jacobian(const double /*t*/,
                                         const double /*dt*/, const vec_type& y,
                                         const mat_type& jac) const {
   
  problem.computeJacobian(member,
             y,
             jac);
              
  }


};

} // Impl
}// TChem

#endif    