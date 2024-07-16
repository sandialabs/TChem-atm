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
using real_type = TChem::real_type;
using ordinal_type = TChem::ordinal_type;
#include "KokkosBatched_Gemm_Decl.hpp"
#include "KokkosBatched_Gemm_TeamVector_Impl.hpp"

int main(int argc, char *argv[]) {

  Kokkos::initialize(argc, argv);
  {

  const bool detail = false;
  TChem::exec_space().print_configuration(std::cout, detail);
  const auto exec_space_instance = TChem::exec_space();
  using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
  using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;

  using policy_type = typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;

  using ats = Tines::ats<real_type>;
  using Trans = Tines::Trans;

  const ordinal_type np = 10000, m = 10;
  Tines::value_type_3d_view<real_type, device_type> A("A", np, m, m);
  Tines::value_type_3d_view<real_type, device_type> B("B", np, m, m);
  Tines::value_type_3d_view<real_type, device_type> C("C", np, m, m);
  Tines::value_type_3d_view<real_type, host_device_type> CC("CC", np, m, m);

  const real_type one(1), zero(.5);

  Kokkos::Random_XorShift64_Pool<device_type> random(13718);
  Kokkos::fill_random(A, random, real_type(1.0));
  Kokkos::fill_random(B, random, real_type(1.0));
  Kokkos::fill_random(C, random, real_type(1.0));
  Kokkos::deep_copy(CC, C);

  const real_type flops = double(np)*double(m)*double(m)*2/1e9;

  typedef typename policy_type::member_type member_type;
  policy_type policy(exec_space_instance, np, Kokkos::AUTO());

  double t_gemm(0);
  {
    Kokkos::Timer timer;
    Kokkos::parallel_for(
    "test GEMM tines",
    policy,
    KOKKOS_LAMBDA(const member_type& member) {
          const int k = member.league_rank();
    auto aa = Kokkos::subview(A, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb = Kokkos::subview(B, k, Kokkos::ALL(), Kokkos::ALL());
    auto cc = Kokkos::subview(C, k, Kokkos::ALL(), Kokkos::ALL());
        Tines::Gemm<Trans::NoTranspose, Trans::NoTranspose>::invoke(member, one, aa, bb,
                                                              zero, cc);
        });
    Kokkos::fence();
    t_gemm = timer.seconds();
  }

  printf("Time per device problem (Tines) %e (s), %e (gflop)\n", t_gemm / double(np), flops/t_gemm);

  double t_gemm_kernels(0);
  {
    Kokkos::Timer timer;
    Kokkos::parallel_for(
    "test GEMM",
    policy,
    KOKKOS_LAMBDA(const member_type& member) {
          const int k = member.league_rank();
    auto aa = Kokkos::subview(A, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb = Kokkos::subview(B, k, Kokkos::ALL(), Kokkos::ALL());
    auto cc = Kokkos::subview(C, k, Kokkos::ALL(), Kokkos::ALL());
    KokkosBatched::TeamVectorGemm<member_type, KokkosBatched::Trans::NoTranspose,
                                  KokkosBatched::Trans::NoTranspose,
                                  KokkosBatched::Algo::Gemm::Unblocked>::invoke(member, one, aa, bb,
                                                       zero, cc);
    });
    Kokkos::fence();
    t_gemm_kernels = timer.seconds();
  }

  printf("Time per device problem (Kokkos-kernels)%e (s), %e (gflop)\n", t_gemm_kernels / double(np), flops/t_gemm_kernels);


  }
  Kokkos::finalize();

  return 0;
}
