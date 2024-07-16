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
#ifndef __TCHEM_IMPL_KFORWARD_JPL_HPP__
#define __TCHEM_IMPL_KFORWARD_JPL_HPP__

#include "TChem_KineticModelData.hpp"
#include "TChem_Util.hpp"
// #define TCHEM_ENABLE_SERIAL_TEST_OUTPUT
namespace TChem {
namespace Impl {

template <typename ValueType, typename DeviceType> struct KForwardJPL {
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
              const real_type &t, const value_type &m,
              /// output
              const value_type_1d_view_type &kfor,
              /// const input from kinetic model
              const kinetic_model_type &kmcd) {

    constexpr real_type zero = 0.0;
    constexpr real_type one = 1.0;
    constexpr real_type two = 2.0;
    constexpr real_type three_hundred = 300.0;
    const real_type t_1 = real_type(1) / t;

    // aux factor
    const ordinal_type n_troe_jlp_reac = kmcd.JPL_Coef.extent(0);
    // jpls
    Kokkos::parallel_for(
        Tines::RangeFactory<value_type>::TeamVectorRange(member,
                                                         n_troe_jlp_reac),
        [&](const ordinal_type &i) {
          const auto param = kmcd.JPL_Coef(i);

          // printf("k0_A %e k0_B %e   \n", param._k0_A, param._k0_B);
          // printf("kinf_A %e kinf_B %e   \n", param._kinf_A, param._kinf_B);
          // printf("Fc %e m %e   \n", param._Fc, m);

          // const value_type k0 =
          //     param._k0_A *
          //     ats<value_type>::pow(t / three_hundred, param._k0_B);

          const value_type k0 = param._k0_A * ( param._k0_C == zero ? one : ats<real_type>::exp(param._k0_C * t_1) ) *
              ( param._k0_B == zero ? one : ats<value_type>::pow (t/three_hundred, param._k0_B) );
        

          const value_type kinf =  param._kinf_A * ( param._kinf_C == zero ? one :ats<real_type>::exp(param._kinf_C * t_1) )*
                                  ( param._kinf_B == zero ? one : ats<real_type>::pow (t/three_hundred, param._kinf_B) );

          // const value_type kinf =
          //     param._kinf_A *
          //     ats<value_type>::pow(t / three_hundred, param._kinf_B);



          const value_type xpo = k0 * m / kinf;
                      
          const value_type factor = param._Fc == one ? value_type(1) : ats<value_type>::pow(param._Fc,
                                   one / (one + ats<value_type>::log10(xpo) * ats<value_type>::log10(xpo)));                        

          kfor(param._reaction_index) = k0 / (one + xpo) * factor;
          // printf("k0 %e  \n", k0);
          // printf("kinf %e  \n", kinf);
          // printf("kfor_t %e  \n", kfor(param._reaction_index));
        });


        const ordinal_type n_ratio_troe_arrhenius_reac = kmcd.R_JPL_ArrheniusCoef.extent(0);


        Kokkos::parallel_for(
        Tines::RangeFactory<value_type>::TeamVectorRange(member,
                                                         n_ratio_troe_arrhenius_reac),
        [&](const ordinal_type &i) {
          const auto param = kmcd.R_JPL_ArrheniusCoef(i);

          // printf("k0_A %e k0_B %e   \n", param._k0_A, param._k0_B);
          // printf("kinf_A %e kinf_B %e   \n", param._kinf_A, param._kinf_B);
          // printf("Fc %e m %e   \n", param._Fc, m);

          // const value_type k0 =
          //     param._k0_A *
          //     ats<value_type>::pow(t / three_hundred, param._k0_B);

          const real_type k0 = param._k0_A * ( param._k0_C == zero ? one : ats<real_type>::exp(param._k0_C * t_1) ) *
              ( param._k0_B == zero ? one : ats<real_type>::pow (t/three_hundred, param._k0_B) );
        

          const real_type kinf =  param._kinf_A * ( param._kinf_C == zero ? one :ats<real_type>::exp(param._kinf_C * t_1) )*
                                  ( param._kinf_B == zero ? one : ats<real_type>::pow (t/three_hundred, param._kinf_B) );

          // const value_type kinf =
          //     param._kinf_A *
          //     ats<value_type>::pow(t / three_hundred, param._kinf_B);

          const value_type xpo = k0 * m / kinf;
                      
          const value_type factor = param._Fc == one ? value_type(1) : ats<value_type>::pow(param._Fc,
                                   one / (one + ats<value_type>::log10(xpo) * ats<value_type>::log10(xpo)));        


          const value_type r_num =  k0 / (one + xpo) * factor;   
          

          const value_type k_deno = param._A * ( param._C == zero ? one : ats<real_type>::exp(param._C * t_1) ) *
              ( param._B == zero ? one : ats<value_type>::pow (t/three_hundred, param._B) );

          kfor(param._reaction_index) = r_num/k_deno;
          // printf("k0 %e  \n", k0);
          // printf("kinf %e  \n", kinf);
          // printf("kfor_t %e  \n", kfor(param._reaction_index));
        });

  }
  }; // KForwardJPL
} // namespace Impl
} // namespace TChem

#endif
