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
#include "TChem_Util.hpp"
#include "TChem_Impl_linv3_strat_chem_solve.hpp"
#include "TChem_Linv3StratosphereSolver.hpp"

namespace TChem {
#if 1
  template<typename PolicyType,
           typename DeviceType>
void
Linv3StratosphereSolver_TemplateRun( /// input
  const std::string& profile_name,
  /// team size setting
  const PolicyType& policy,
  const Tines::value_type_1d_view<real_type, DeviceType>& temperature,
  const Tines::value_type_1d_view<real_type, DeviceType>& pressure,

  const Tines::value_type_2d_view<real_type, DeviceType>& volume_mixing_ratio,
  const real_type dt, const real_type rlats, const real_type psc_T, const real_type  sza,
  const real_type chlorine_loading, 
  const Tines::value_type_1d_view<real_type, DeviceType>& o3col,
  const Tines::value_type_1d_view<ordinal_type, DeviceType>& tropFlag,
  const Tines::value_type_1d_view<real_type, DeviceType>& water_vapor_volume_mixing_ratio,
  const Tines::value_type_1d_view<linoz_input_parameters_type, DeviceType>& linoz_inputs, 
  const linoz_vmr_idx_type& linoz_vmr_idx,
  // output 
  const Tines::value_type_2d_view<real_type, DeviceType>& volume_mixing_ratio_out
  /// output
  )
{
  Kokkos::Profiling::pushRegion(profile_name);

  using policy_type = PolicyType;
  using device_type = DeviceType;
  constexpr real_type pressure_threshold = 237.1375e+2; //last vaild linoz pressure layer
  // convert lats from radians to degrees
  const real_type lats = rlats * 180./PI();

  Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
      const ordinal_type i = member.league_rank();

      if (!tropFlag(i) &&  pressure (i) <= pressure_threshold )
      {
      const real_type o3_vmr = volume_mixing_ratio(i, linoz_vmr_idx.o3_ndx);
      const real_type n2o_vmr = volume_mixing_ratio(i, linoz_vmr_idx.n2olnz_ndx);
      const real_type noy_vmr = volume_mixing_ratio(i, linoz_vmr_idx.noylnz_ndx);
      const real_type ch4_vmr = volume_mixing_ratio(i, linoz_vmr_idx.ch4lnz_ndx); 
      const real_type h2o_vmr = volume_mixing_ratio(i, linoz_vmr_idx.h2olnz_ndx);
  
      real_type o3_new;
      real_type n2o_new;
      real_type noy_new;
      real_type ch4_new;
      real_type h2o_new;
      real_type o3_value_ndx;
      real_type n2o_value_ndx;
      real_type ch4_value_ndx;
      real_type no_value_ndx;
      real_type no2_value_ndx;
      real_type hno3_value_ndx;
          // outputs diagnostics 
      real_type do3_linoz_psc;
        Impl::Linv3StratosphereSolver<real_type, device_type>
        ::team_invoke(member,
              temperature(i),
              pressure(i),
              dt, 
              rlats, 
              psc_T,
              sza, 
              chlorine_loading,
              o3col(i),
              linoz_inputs(i), 
              // inputs
              o3_vmr, 
              n2o_vmr,
              noy_vmr,
              ch4_vmr, 
              h2o_vmr,
              // outputs
              o3_new,
              n2o_new,
              noy_new,
              ch4_new,
              h2o_new,
              o3_value_ndx,
              n2o_value_ndx,
              ch4_value_ndx,
              no_value_ndx,
              no2_value_ndx,
              hno3_value_ndx,
          // outputs diagnostics 
          do3_linoz_psc);

        member.team_barrier();

        if (linoz_vmr_idx.o3lnz_ndx > -1) volume_mixing_ratio_out(i,  linoz_vmr_idx.o3lnz_ndx) =   o3_new;
        if (linoz_vmr_idx.n2olnz_ndx > -1) volume_mixing_ratio_out(i, linoz_vmr_idx.n2olnz_ndx)   = n2o_new;
        if (linoz_vmr_idx.noylnz_ndx > -1) volume_mixing_ratio_out(i, linoz_vmr_idx.noylnz_ndx)   = noy_new;
        if (linoz_vmr_idx.ch4lnz_ndx > -1) volume_mixing_ratio_out(i, linoz_vmr_idx.ch4lnz_ndx)   = ch4_new;
        if (linoz_vmr_idx.h2olnz_ndx > -1) volume_mixing_ratio_out(i, linoz_vmr_idx.h2olnz_ndx)   = h2o_new;

        //update real o3, ch4, n2o      
        if(linoz_vmr_idx.o3_ndx  > -1) volume_mixing_ratio_out(i, linoz_vmr_idx.o3_ndx ) =  o3_value_ndx +  volume_mixing_ratio_out(i, linoz_vmr_idx.o3_ndx );
        if(linoz_vmr_idx.ch4_ndx > -1) volume_mixing_ratio_out(i, linoz_vmr_idx.ch4_ndx) =  ch4_value_ndx  +  volume_mixing_ratio_out(i, linoz_vmr_idx.ch4_ndx);
        if(linoz_vmr_idx.n2o_ndx > -1) volume_mixing_ratio_out(i, linoz_vmr_idx.n2o_ndx) =  n2o_value_ndx +  volume_mixing_ratio_out(i, linoz_vmr_idx.n2o_ndx);
        if(linoz_vmr_idx.no_ndx >-1)  volume_mixing_ratio_out(i, linoz_vmr_idx.no_ndx)   =  no_value_ndx + volume_mixing_ratio_out(i, linoz_vmr_idx.no_ndx);
        if(linoz_vmr_idx.no2_ndx>-1)  volume_mixing_ratio_out(i, linoz_vmr_idx.no2_ndx)  =  no2_value_ndx+ volume_mixing_ratio_out(i, linoz_vmr_idx.no2_ndx);
        if(linoz_vmr_idx.hno3_ndx>-1) volume_mixing_ratio_out(i,linoz_vmr_idx.hno3_ndx)  =  hno3_value_ndx + volume_mixing_ratio_out(i, linoz_vmr_idx.hno3_ndx);
       }// end if
       



      member.team_barrier();
      if (tropFlag(i)) {
              volume_mixing_ratio_out(i, linoz_vmr_idx.h2olnz_ndx)   = water_vapor_volume_mixing_ratio(i);
      }

    });

  Kokkos::Profiling::popRegion();
}
void
Linv3StratosphereSolver::runHostBatch( /// input
  const real_type_1d_view_host_type& temperature,
  const real_type_1d_view_host_type& pressure,

  const real_type_2d_view_host_type& volume_mixing_ratio,
  const real_type dt, const real_type rlats, const real_type psc_T, const real_type  sza,
  const real_type chlorine_loading, 
  const real_type_1d_view_host_type& o3col,
  const ordinal_type_1d_view_host_type& tropFlag,
  const real_type_1d_view_host_type& water_vapor_volume_mixing_ratio,
  const linoz_input_parameters_1d_view_host& linoz_inputs, 
  const linoz_vmr_idx_type& linoz_vmr_idx,
  const real_type_2d_view_host_type& volume_mixing_ratio_out
  )
{
  using policy_type = Kokkos::TeamPolicy<host_exec_space>;
  policy_type policy(volume_mixing_ratio.extent(0), Kokkos::AUTO()); // fine
 
  Linv3StratosphereSolver_TemplateRun( /// input
    "TChem::Linv3StratosphereSolver::runHostBatch",
    policy,
    temperature,
    pressure,
    volume_mixing_ratio,
    dt,
    rlats, 
    psc_T, 
    sza,
    chlorine_loading, 
    o3col,
    tropFlag,
    water_vapor_volume_mixing_ratio,
    linoz_inputs, 
    linoz_vmr_idx,
    volume_mixing_ratio_out);

}

void
Linv3StratosphereSolver::runDeviceBatch( /// input
  typename UseThisTeamPolicy<exec_space>::type& policy,
   const real_type_1d_view_type& temperature,
  const real_type_1d_view_type& pressure,

 const real_type_2d_view_type& volume_mixing_ratio,
  const real_type dt, const real_type rlats, const real_type psc_T, const real_type  sza,
  const real_type chlorine_loading, 
  const real_type_1d_view_type& o3col,
  const ordinal_type_1d_view_type& tropFlag,
  const real_type_1d_view_type& water_vapor_volume_mixing_ratio,
  const linoz_input_parameters_1d_view& linoz_inputs, 
  const linoz_vmr_idx_type& linoz_vmr_idx,
  const real_type_2d_view_type& volume_mixing_ratio_out
  )
{

  Linv3StratosphereSolver_TemplateRun( /// input
    "TChem::Linv3StratosphereSolver::runDeviceBatch",
    policy,
    temperature,
    pressure,
    volume_mixing_ratio,
    dt,
    rlats, 
    psc_T, 
    sza,
    chlorine_loading, 
    o3col,
    tropFlag,
    water_vapor_volume_mixing_ratio,
    linoz_inputs, 
    linoz_vmr_idx,
    volume_mixing_ratio_out);

}

void
Linv3StratosphereSolver::runDeviceBatch( /// input
  const real_type_1d_view_type& temperature,
  const real_type_1d_view_type& pressure,

  const real_type_2d_view_type& volume_mixing_ratio,
  const real_type dt, const real_type rlats, const real_type psc_T, const real_type  sza,
  const real_type chlorine_loading, 
  const real_type_1d_view_type& o3col,
  const ordinal_type_1d_view_type& tropFlag,
  const real_type_1d_view_type& water_vapor_volume_mixing_ratio,
  const linoz_input_parameters_1d_view& linoz_inputs, 
  const linoz_vmr_idx_type& linoz_vmr_idx,
  const real_type_2d_view_type& volume_mixing_ratio_out
  )
{
  using policy_type = Kokkos::TeamPolicy<exec_space>;

  policy_type policy(volume_mixing_ratio.extent(0), Kokkos::AUTO()); // fine
  Linv3StratosphereSolver_TemplateRun( /// input
    "TChem::Linv3StratosphereSolver::runDeviceBatch",
    policy,
    temperature,
    pressure,
    volume_mixing_ratio,
    dt,
    rlats, 
    psc_T, 
    sza,
    chlorine_loading, 
    o3col,
    tropFlag,
    water_vapor_volume_mixing_ratio,
    linoz_inputs, 
    linoz_vmr_idx,
    volume_mixing_ratio_out);

}

#endif
} // namespace TChem
