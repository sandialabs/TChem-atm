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
#ifndef __TCHEM_IMPL_ADJUST_REACTIONS_JPL_HPP__
#define __TCHEM_IMPL_ADJUST_REACTIONS_JPL_HPP__

#include "TChem_KineticModelData.hpp"
#include "TChem_Util.hpp"
// #define TCHEM_ENABLE_SERIAL_TEST_OUTPUT
namespace TChem {
namespace Impl {

template <typename ValueType, typename DeviceType> struct AdjustReactions {
  using value_type = ValueType;
  using device_type = DeviceType;
  using scalar_type = typename ats<value_type>::scalar_type;

  using real_type = scalar_type;
  /// sacado is value type
  using value_type_1d_view_type =
      Tines::value_type_1d_view<value_type, device_type>;
  using kinetic_model_type = KineticModelNCAR_ConstData<device_type>;

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION static void
  team_invoke(const MemberType &member,
              /// input temperature
              const real_type &t,
              const value_type_1d_view_type &x,
              /// output
              const value_type_1d_view_type &kfor,
              /// const input from kinetic model
              const kinetic_model_type &kmcd) {

    const ordinal_type n_adjust_reaction = kmcd.adjust_reaction.extent(0);

    const real_type t_1 = 1.0 / t;

    // Kokkos::parallel_for(
    //     Tines::RangeFactory<value_type>::TeamVectorRange(member,
    //                                                      n_adjust_reaction),
    //     [&](const ordinal_type &i) {
    // FIXME: will Kokkos::single slow down this computation?
      Kokkos::single(Kokkos::PerTeam(member), [&]() {
    for (ordinal_type i = 0; i < n_adjust_reaction; i++)
    {
      const auto param = kmcd.adjust_reaction(i);
        kfor(param._reaction_index) *= x(param._species_index);
    }
    });

    member.team_barrier();

    // adjust reaction in uci1, uci2, and uc3
    Kokkos::parallel_for(
        Tines::RangeFactory<value_type>::TeamVectorRange(member,
                                                         kmcd.number_of_prod_O1D),
        [&](const ordinal_type &i) {

         const auto param = kmcd.prod_O1D(i);
         const real_type k1 = param._A1 * ats<real_type>::exp(param._C1 * t_1);
         const real_type k2 = param._A2 * ats<real_type>::exp(param._C2 * t_1);
         const real_type k3 = param._A3 * ats<real_type>::exp(param._C3 * t_1);
         const value_type fc = k1*x(param._species_index_1) + k2*x(param._species_index_2) + k3*x(param._species_index_3);
         // printf("fc %e \n", fc);
         // printf("k1 %e x %e \n", k1,x(param._species_index_1));
         // printf("k2 %e x %e \n", k2,x(param._species_index_2));
         // printf("k3 %e x %e \n", k3,x(param._species_index_3));
         // printf("t_1 %e \n", t_1);
         // printf("A1 %e C1 %e \n", param._A1,param._C1);

         for (int ireac = 0; ireac < 3; ++ireac)
         {
          // printf("param._reaction_indices[%d] %d \n", ireac, param._reaction_indices[ireac]);
           kfor(param._reaction_indices[ireac]) *= kfor(param._photolysis_reaction_index)/fc;
         }


    });



  }
  }; // AdjustReactions
} // namespace Impl
} // namespace TChem

#endif
