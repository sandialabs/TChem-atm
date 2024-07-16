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
#ifndef __TCHEM_IMPL_TIME_INTEGRATOR_EXPLICIT_EULER_HPP__
#define __TCHEM_IMPL_TIME_INTEGRATOR_EXPLICIT_EULER_HPP__

#include "Tines.hpp"
using namespace Tines;

namespace TChem {
namespace Impl {

template <typename ValueType, typename DeviceType>
  struct TimeIntegratorExplicitEuler {
    using value_type = ValueType;
    using device_type = DeviceType;
    using scalar_type = typename ats<value_type>::scalar_type;

    using real_type = scalar_type;
    using real_type_0d_view_type = value_type_0d_view<real_type, device_type>;
    using real_type_1d_view_type = value_type_1d_view<real_type, device_type>;

    KOKKOS_INLINE_FUNCTION
    static void workspace(const int m, int &wlen) {
    // rhs	
      wlen = 2*m;
    }

template <typename MemberType,
              template <typename, typename> class ProblemType>
    KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType &member,
      /// problem
      const ProblemType<value_type, device_type> &problem,
      /// constant time step 
      const int &max_num_time_iterations,
      const real_type &dt_in,
      const real_type &t_beg,const real_type &t_end,
      /// input (initial condition)
      const real_type_1d_view_type &vals,
      /// output (final output conditions)
      const real_type_0d_view_type &t_out,
      const real_type_0d_view_type& dt_out,
      const real_type_1d_view_type &vals_out,
      /// workspace
      const real_type_1d_view_type &work) 
    {

     
     constexpr real_type zero(0);	
     /// early return
      if (dt_in < zero)
        return 3;
     
      /// workspace
      auto wptr = work.data();
      const int m = problem.getNumberOfEquations();
      auto fn = real_type_1d_view_type(wptr, m);
      wptr += m;
      auto un = real_type_1d_view_type(wptr, m);
      wptr += m;

      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m),
                           [&](const int &k) {
                             un(k) = vals(k);
                           });
      member.team_barrier();

      /// time integration
      real_type t(t_beg), dt(dt_in);
      for (int iter = 0; iter < max_num_time_iterations && dt != zero; ++iter)
      {
       // compute rhs 	
       // printf("dt %e t %e iter %d \n", dt, t, iter);	
       problem.computeFunction(member, un, fn);
       // time advance with explicit euler 
       Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m),
                                 [&](const int &k) { 
                                 	un(k) += fn(k)*dt;	
                                 	});
       t += dt;
       dt = ((t + dt) > t_end) ? t_end - t : dt;
      }
      member.team_barrier();

      /// finalize with output for next iterations of time solutions
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m),
                               [&](const int &k) {
                                 vals_out(k) = un(k);
                                 if (k == 0) {
                                   t_out()  = t;
                                   dt_out() = dt;
                                 }
      });

     return 0;
}
    };	
} // namespace Impl
} // namespace TChem

#endif	