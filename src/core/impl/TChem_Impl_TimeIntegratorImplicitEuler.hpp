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
#ifndef __TCHEM_IMPL_TIME_INTEGRATOR_IMPLICIT_EULER_HPP__
#define __TCHEM_IMPL_TIME_INTEGRATOR_IMPLICIT_EULER_HPP__

#include "Tines.hpp"
using namespace Tines;

namespace TChem {
namespace Impl {


  template <typename ValueType, typename DeviceType,
            template <typename, typename> class ProblemType>
  struct ImplicitEuler {
    using value_type = ValueType;
    using device_type = DeviceType;
    using scalar_type = typename ats<value_type>::scalar_type;

    using real_type = scalar_type;
    using real_type_1d_view_type = value_type_1d_view<real_type, device_type>;
    using real_type_2d_view_type = value_type_2d_view<real_type, device_type>;

    using problem_type = ProblemType<value_type, device_type>;

    problem_type _problem;
    real_type _dt;
    real_type_1d_view_type _un, _fn;

    KOKKOS_INLINE_FUNCTION
    ImplicitEuler()
      : _problem(),_dt(),
        _un(), _fn() {}


    KOKKOS_INLINE_FUNCTION
    void setWorkspace(real_type_1d_view_type &work) {
      const int m = _problem.getNumberOfEquations();
      assert(2 * m <= int(work.extent(0)) &&
             "Error: workspace is smaller than required");
      _un = real_type_1d_view_type(work.data() + 0 * m, m);
      _fn = real_type_1d_view_type(work.data() + 1 * m, m);
    }

    KOKKOS_INLINE_FUNCTION
    int getNumberOfTimeODEs() const { return _problem.getNumberOfTimeODEs(); }

    KOKKOS_INLINE_FUNCTION
    int getNumberOfConstraints() const {
      return _problem.getNumberOfConstraints();
    }

    KOKKOS_INLINE_FUNCTION
    int getNumberOfEquations() const { return _problem.getNumberOfEquations(); }

    template <typename MemberType>
    KOKKOS_INLINE_FUNCTION void
    computeInitValues(const MemberType &member,
                      const real_type_1d_view_type &u) const {
      const int m = _problem.getNumberOfEquations();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m),
                           [&](const int &i) { u(i) = _un(i); });
      member.team_barrier();
    }

    template <typename MemberType>
    KOKKOS_INLINE_FUNCTION void
    computeJacobian(const MemberType &member, const real_type_1d_view_type &u,
                    const real_type_2d_view_type &J) const {
      const real_type one(1), zero(0), half(0.5);
      const int m = _problem.getNumberOfTimeODEs(),
                n = _problem.getNumberOfEquations();

      /// evaluate problem Jacobian (n x n)
      _problem.computeJacobian(member, u, J);

      /// modify time ODE parts for the trapezoidal rule
      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(member, m), [&](const int &i) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, n),
                               [&](const int &j) {
                                 const real_type val = _dt*J(i, j);
                                 J(i, j) = (i == j ? one : zero) - val;
                               });
        });
      member.team_barrier();
    }

    template <typename MemberType>
    KOKKOS_INLINE_FUNCTION void
    computeFunction(const MemberType &member, const real_type_1d_view_type &u,
                    const real_type_1d_view_type &f) const {
      const int m = _problem.getNumberOfTimeODEs();

      /// evaluate problem function (n x 1)
      _problem.computeFunction(member, u, f);

      /// modify time ODE parts for the trapezoidal rule
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m),
                           [&](const int &i) {
                             const real_type val = f(i);
                             f(i) = (u(i) - _un(i)) - _dt * val;
                           });
      member.team_barrier();
    }
  };

 template <typename ValueType, typename DeviceType>
  struct TimeIntegratorImplicitEuler {
    using value_type = ValueType;
    using device_type = DeviceType;
    using scalar_type = typename ats<value_type>::scalar_type;

    using real_type = scalar_type;
    using real_type_0d_view_type = value_type_0d_view<real_type, device_type>;
    using real_type_1d_view_type = value_type_1d_view<real_type, device_type>;
    using real_type_2d_view_type = value_type_2d_view<real_type, device_type>;

    KOKKOS_INLINE_FUNCTION
    static void workspace(const int m, int &wlen) {
      using newton_solver_type = NewtonSolver<value_type, device_type>;
      // using trbdf2_type = TrBDF2<value_type, device_type>;
      /// problem.setWorkspace should be invoked before
      int wlen_newton(0);
      newton_solver_type::workspace(m, wlen_newton); /// utv workspace
      int wlen_implicit_euler=2*m;
      // trbdf2_type::workspace(m, wlen_trbdf); /// un, unr, fn
      const int wlen_this = (2 * m /* u, fnr */ + 2 * m + m * m /* dx, f, J */);

      wlen = (wlen_newton + wlen_implicit_euler + wlen_this);
    }

    template <typename MemberType,
              template <typename, typename> class ProblemType>
    KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType &member,
      /// problem
      const ProblemType<value_type, device_type> &problem,
      /// input iteration and qoi index to store
      const int &jacobian_interval,
      const int &max_num_newton_iterations,
      const int &max_num_time_iterations,
      const real_type_1d_view_type &tol_newton,
      const real_type_2d_view_type &tol_time,
      /// input time step and time range
      const real_type &dt_in,const real_type &dt_min,
      const real_type &t_beg, const real_type &t_end,
      /// input (initial condition)
      const real_type_1d_view_type &vals,
      /// output (final output conditions)
      const real_type_0d_view_type &t_out, const real_type_0d_view_type &dt_out,
      const real_type_1d_view_type &vals_out,
      /// workspace
      const real_type_1d_view_type &work) {
      using newton_solver_type = NewtonSolver<value_type, device_type>;

      using implicit_euler_type =
        ImplicitEuler<value_type, device_type, ProblemType>;
      /// return value; when it fails it return non-zero value
      int r_val(0);

      /// const values
      const real_type zero(0),half(0.5), minus_one(-1);

      /// early return
      if (dt_in < zero)
        return 3;

      /// data structure here is temperature, mass fractions of species...
      const int m = problem.getNumberOfEquations(),
                m_ode = problem.getNumberOfTimeODEs();

      /// time stepping object
      implicit_euler_type implicit_euler;

      /// to compute workspace correctly the problem information should be given first
      /// assign the problem to trbdf
      implicit_euler._problem = problem;
      /// workspace
      auto wptr = work.data();

      int wlen_newton(0);
      newton_solver_type::workspace(m, wlen_newton); /// utv workspace
      auto work_newton = real_type_1d_view_type(wptr, wlen_newton);
      wptr += wlen_newton;

      int wlen_implicit_euler=2*m;
      // trbdf2_type::workspace(m, wlen_trbdf); /// un, unr, fn
      auto work_implicti_euler = real_type_1d_view_type(wptr, wlen_implicit_euler);
      wptr += wlen_implicit_euler;
      implicit_euler.setWorkspace(work_implicti_euler);

      auto un = implicit_euler._un;
      auto fn = implicit_euler._fn;

      /// newton workspace
      auto dx = real_type_1d_view_type(wptr, m);
      wptr += m;
      auto f = real_type_1d_view_type(wptr, m);
      wptr += m;
      auto J = real_type_2d_view_type(wptr, m, m);
      wptr += m * m;

      auto u = real_type_1d_view_type(wptr, m);
      wptr += m;

      /// error check
      const int workspace_used(wptr - work.data()),
        workspace_extent(work.extent(0));

      assert(workspace_used <= workspace_extent &&
             "Error: workspace is used more than allocated");

      /// initial conditions
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m),
                           [&](const int &k) {
                             un(k) = vals(k);
                           });
      member.team_barrier();

      /// time integration
      real_type t(t_beg), dt(dt_in);
      for (int iter = 0; iter < max_num_time_iterations && dt != zero; ++iter) {
        {
          int converge(0);
          for (int i = 0; i < 4 && converge == 0; ++i) {
              dt = (dt > dt_min ? dt : dt_min);
              implicit_euler._dt = dt;
              int newton_iteration_count(0);
              newton_solver_type::invoke(
                member, implicit_euler, jacobian_interval,
                tol_newton(0), tol_newton(1),
                max_num_newton_iterations,
                u, dx, f, J, work_newton,
                newton_iteration_count, converge);

              if (!converge) {
                /// try again with half time step
                dt *= half;
                continue;
              }
          }// i

          if (converge) {
            t += dt;
            dt = ((t + dt) > t_end) ? t_end - t : dt;
            Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m),
                                 [&](const int &k) { un(k) = u(k); });
          } else {
            Kokkos::single(Kokkos::PerTeam(member), [&]() {
              printf("Warning: TimeIntegrator, sample (%d) trbdf fails to "
                     "converge with current time step %e\n",
                     int(member.league_rank()), dt);
            });
            r_val = 1;
            break;
          }
        }
        member.team_barrier();
      }

      {
        /// finalize with output for next iterations of time solutions
        if (r_val == 0) {
          Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m),
                               [&](const int &k) {
                                 vals_out(k) = u(k);
                                 if (k == 0) {
                                   t_out() = t;
                                   dt_out() = dt;
                                 }
                               });
        } else {
          /// if newton fails,
          /// - set values with zero
          /// - t_out becomes t_end
          /// - dt_out is minus one
          Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m),
                               [&](const int &k) {
                                 vals_out(k) = zero;
                                 if (k == 0) {
                                   t_out() = t_end;
                                   dt_out() = minus_one;
                                 }
                               });
        }
      }

      return r_val;
    }

    template <typename MemberType,
              template <typename, typename> class ProblemType>
    KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType &member,
      /// problem
      const ProblemType<value_type, device_type> &problem,
      /// input iteration and qoi index to store
      const int &max_num_newton_iterations,
      const int &max_num_time_iterations,
      const real_type_1d_view_type &tol_newton,
      const real_type_2d_view_type &tol_time,
      /// input time step and time range
      const real_type &dt_in,const real_type &dt_min,
      const real_type &t_beg, const real_type &t_end,
      /// input (initial condition)
      const real_type_1d_view_type &vals,
      /// output (final output conditions)
      const real_type_0d_view_type &t_out, const real_type_0d_view_type &dt_out,
      const real_type_1d_view_type &vals_out,
      /// workspace
      const real_type_1d_view_type &work) {

      const int jacobian_interval(1);
      return invoke(member, problem,
                    jacobian_interval,
                    max_num_newton_iterations,
                    max_num_time_iterations,
                    tol_newton,
                    tol_time,
                    dt_in, dt_min,
                    t_beg, t_end,
                    vals,
                    t_out, dt_out,
                    vals_out,
                    work);
    }
  };


} // namespace Impl
} // namespace TChem

#endif	
