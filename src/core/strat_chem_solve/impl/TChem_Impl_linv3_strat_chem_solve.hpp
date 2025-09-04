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
#ifndef __TCHEM_IMPL_LINV3_START_CHEM_SOLVE_HPP__
#define __TCHEM_IMPL_LINV3_START_CHEM_SOLVE_HPP__

#include "TChem_Util.hpp"

/*--------------------------------------------------------------------
 linearized ozone chemistry linoz-v3 (o3-n2o-noy-ch4 prognostic equations, h2o diagnosed from ch4) 
 from Hsu and Prather, grl, 2009   (https://doi.org/10.1029/2009GL042243)
 
written by Juno Hsu (junoh@uci.edu), 09/2020 */

namespace TChem {
namespace Impl {

template <typename ValueType, typename DeviceType> struct Linv3StratosphereSolver{

  using value_type = ValueType;
  using device_type = DeviceType;
  using scalar_type = typename ats<value_type>::scalar_type;

  using real_type = scalar_type;
  /// sacado is value type
  using value_type_1d_view_type =
      Tines::value_type_1d_view<value_type, device_type>;
  
  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION static void
  team_invoke_detail(const MemberType &member,
  	          const real_type &temp,
  	          const real_type &pmid,
  	          const real_type &delta_t, 
  	          const real_type &lats, 
  	          const real_type &psc_T,
  	          const real_type &sza, 
  	          const real_type &chlorine_loading, 

  	          // inputs
  	          const real_type &o3_vmr, 
  	          const real_type &n2o_vmr,
  	          const real_type &noy_vmr,
  	          const real_type &ch4_vmr, 
  	          const real_type &h2o_vmr,
  	          const real_type &linoz_o3_clim,
  	          const real_type &linoz_dPmL_dO3X,  
  	          const real_type &linoz_cariolle_psc, 
  	          const real_type &PL_O3, 
  	          const real_type &P_n2o, 
  	          const real_type &Lfreq_n2o,
  	          const real_type &Pfreq_noy,
  	          const real_type &Lfreq_noy,
  	          const real_type &Lfreq_ch4,

  	          // outputs
  	          real_type &o3_new,
              real_type &n2o_new,
              real_type &noy_new,
              real_type &ch4_new,
              real_type &h2o_new,
              real_type &o3_value_ndx,
              real_type &n2o_value_ndx,
              real_type &ch4_value_ndx,
              real_type &no_value_ndx,
   			  real_type &no2_value_ndx,
   			  real_type &hno3_value_ndx,
   			  // outputs diagnostics 
   			  real_type &do3_linoz_psc
  	          ) 
  {

  	constexpr real_type one(1.0);
  	constexpr real_type two(2.0);
  	constexpr real_type zero(0.0);
  	//The  local maxval calculation causes NBFB when ncol varies due to threading or pe-layout change
    // Tentatively uses a fixed value
     constexpr real_type ch4max =  1.8e-6;
     constexpr real_type pw= two * ch4max + 3.65e-6;

    constexpr real_type radians_to_degrees = 180./PI(); 
    constexpr real_type chlorine_loading_1987    = 2.5977; 
    constexpr real_type chlorine_loading_bgnd    = 0.0000; 

  // inputs : o3_vmr, n2o_vmr, noy_vmr, ch4_vmr, h2o_vmr
  // linoz_o3_clim, PL_O3, linoz_dPml_dO3X, delta_t, P_n2o,Lfreq_n2o
  // Pfreq_noy, Lfreq_noy, Lfreq_ch4, lats, chlorine_loading, chlorine_loading_bgnd, ,psc_T	
  	// sza, linoz_cariolle_psc, chlorine_loading_1987

  	// outputs: do3_linoz_psc
  // current mixing ratio
   real_type o3_old  =   o3_vmr;
   real_type n2o_old =   n2o_vmr;
   real_type noy_old =   noy_vmr;
   real_type ch4_old =   ch4_vmr;
   real_type h2o_old =   h2o_vmr;
   // climatological ozone
   real_type o3_clim = linoz_o3_clim;

   // o3          
   real_type ss_x= o3_clim - PL_O3 / linoz_dPmL_dO3X;            
   real_type delo3 = (ss_x - o3_old ) * (one - ats<real_type>::exp(linoz_dPmL_dO3X*delta_t));
   o3_new = o3_old + delo3;

   //n2o  
   real_type dn2op   = P_n2o * delta_t < zero ? P_n2o * delta_t : zero ;//  (dp/dt *dt in vmr)
   real_type Lfreq   = Lfreq_n2o < zero ? Lfreq_n2o : zero;//  n2o loss frequency          
   real_type dn2ol   = n2o_old*(ats<real_type>::exp(-Lfreq * delta_t) - one);
   n2o_new = n2o_old + dn2op + dn2ol;

   // noy
   const real_type value = Pfreq_noy * n2o_old* delta_t;
   real_type dnoyp =   value <  zero ? value : zero ;//  scaled by fn2o in Pfreq so multiply it back to get mr/sec
   Lfreq =   Lfreq_noy < zero ? Lfreq_noy : zero; //n2o loss frequency          
   real_type dnoyl =   noy_old*(exp(-Lfreq*delta_t) - one);
   noy_new = noy_old + dnoyp + dnoyl; 

   // 
   // ch4   
   Lfreq   =   Lfreq_ch4 < zero ? Lfreq_ch4 : zero;                    
   real_type delch4  =    ch4_old*(ats<real_type>::exp(-Lfreq*delta_t) - one);
   ch4_new =    ch4_old + delch4;
   h2o_new =    pw - two * ch4_new;

/*alternative h2o, -2*delch4 gain as the loss of delch4             
             h2o_new  =    h2o_old - 2._r8* del1ch4          
           PSC activation (follows Cariolle et al 1990.)
   use only if abs(latitude) > 40. */

	
   o3_value_ndx = zero; 
   if ( ats<real_type>::abs(lats) > 40.  ) {
   if ( (chlorine_loading-chlorine_loading_bgnd) > zero )  {
   if ( temp <= psc_T ) {
   	// define maximum SZA for PSC loss (= tangent height at sunset)
    const real_type value = 16.*ats<real_type>::log10(100000./pmid);
   	real_type max_sza = 90. + ats<real_type>::sqrt( value < zero ? value : zero);
   	if ( (sza*radians_to_degrees) <= max_sza )  {
        const real_type ratio_chlorine = chlorine_loading/chlorine_loading_1987;
   		const real_type psc_loss = ats<real_type>::exp(-linoz_cariolle_psc * ratio_chlorine * ratio_chlorine * delta_t );
   		// update from last o3_new due to gas chemisry change
        // use o3_new to prevent negative value                     
        const real_type delo3_psc = o3_new*(psc_loss -one);
        o3_new += delo3_psc;  //o3_new:update from gas chem loss
        do3_linoz_psc = delo3_psc/delta_t;
        o3_value_ndx = delo3   + delo3_psc + o3_old;

   	}// end sza*radians_to_degrees) <= max_sza

   } // temp <= psc_T
   
   } // end if chlorine_loading-chlorine_loading_bgnd) > zero

   } // end if ats<real_type>::abs(lats) > 40. 

   // update vmr
   // o3_new, n2o_new, noy_new, ch4_new, h2o_new

   // update real o3, ch4, n2o    
   // o3_value_ndx
   // if(o3_ndx  > 0) xvmr(i,k, o3_ndx ) =  delo3   + delo3_psc +  xvmr(i,k, o3_ndx )
   ch4_value_ndx = delch4 + ch4_old; 
   // if(ch4_ndx > 0) xvmr(i,k, ch4_ndx) =  delch4  +  xvmr(i,k, ch4_ndx)
   n2o_value_ndx = dn2op + dn2ol + n2o_old; 
   // if(n2o_ndx > 0) xvmr(i,k, n2o_ndx) =  (dn2op + dn2ol)  +  xvmr(i,k, n2o_ndx)
   no_value_ndx  = 0.05 *(dnoyp + dnoyl); 
   // if(no_ndx >0)  xvmr(i,k, no_ndx)   =  0.05 *(dnoyp + dnoyl) + xvmr(i,k, no_ndx)
   no2_value_ndx = 0.05 *(dnoyp + dnoyl);
   // if(no2_ndx>0)  xvmr(i,k, no2_ndx)  =  0.05 *(dnoyp + dnoyl) + xvmr(i,k, no2_ndx)
   hno3_value_ndx = 0.90 *(dnoyp + dnoyl); 
   // if(hno3_ndx>0) xvmr(i,k,hno3_ndx)  =  0.90 *(dnoyp + dnoyl) + xvmr(i,k, hno3_ndx)  

  }


  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION static void
  team_invoke(const MemberType &member,
  	          const real_type &temp,
  	          const real_type &pmid,
  	          const real_type &delta_t, 
  	          const real_type &lats, 
  	          const real_type &psc_T,
  	          const real_type &sza, 
  	          const real_type &chlorine_loading,
  	          const real_type &o3col,
  	          const LinozInputParameters& linoz_inputs, 
  	          // inputs
  	          const real_type &o3_vmr, 
  	          const real_type &n2o_vmr,
  	          const real_type &noy_vmr,
  	          const real_type &ch4_vmr, 
  	          const real_type &h2o_vmr,
              // outputs
  	          real_type &o3_new,
              real_type &n2o_new,
              real_type &noy_new,
              real_type &ch4_new,
              real_type &h2o_new,
              real_type &o3_value_ndx,
              real_type &n2o_value_ndx,
              real_type &ch4_value_ndx,
              real_type &no_value_ndx,
   			  real_type &no2_value_ndx,
   			  real_type &hno3_value_ndx,
   			  // outputs diagnostics 
   			  real_type &do3_linoz_psc

  	          )
  {

    constexpr real_type one(1.0);
    constexpr real_type convert_to_du = one/2.687e16;//      ! convert ozone column from mol/cm^2 to DU

    const real_type dO3     =   o3_vmr  - linoz_inputs.linoz_o3_clim;
    const real_type dN2O    =  n2o_vmr  - linoz_inputs.linoz_n2o_clim;
    const real_type dNOY    =  noy_vmr  - linoz_inputs.linoz_noy_clim; 
    const real_type dCH4    =  ch4_vmr  - linoz_inputs.linoz_ch4_clim;
    const real_type dH2O    =  h2o_vmr  - linoz_inputs.linoz_h2o_clim; 
    const real_type dTemp   =  temp - linoz_inputs.linoz_t_clim;
    const real_type dCOL     =  o3col*convert_to_du - linoz_inputs.linoz_o3col_clim;

   // for steady-state sol. plug-in, no need for dPmL_dO3 term 
    const real_type PL_O3 =  linoz_inputs.linoz_PmL_clim_no3        
            +  linoz_inputs.linoz_dPmL_dn2o_no3   * dN2O  
            +  linoz_inputs.linoz_dPmL_dnoy_no3   * dNOY  
            +  linoz_inputs.linoz_dPmL_dch4_no3   * dCH4  
            +  linoz_inputs.linoz_dPmL_dh2o_no3    * dH2O  
            +  linoz_inputs.linoz_dPmL_dT_no3       * dTemp 
            +  linoz_inputs.linoz_dPmL_dO3col_no3  * dCOL;

// Taylor-expansion of P (mr/sec)
    const real_type P_n2o =  linoz_inputs.linoz_PmL_clim_pn2o
            +  linoz_inputs.linoz_dPmL_dO3_pn2o   * dO3
            +  linoz_inputs.linoz_dPmL_dn2o_pn2o   * dN2O
            +  linoz_inputs.linoz_dPmL_dnoy_pn2o   * dNOY
            +  linoz_inputs.linoz_dPmL_dch4_pn2o   * dCH4
            +  linoz_inputs.linoz_dPmL_dh2o_pn2o   * dH2O
            +  linoz_inputs.linoz_dPmL_dT_pn2o     * dTemp
            +  linoz_inputs.linoz_dPmL_dO3col_pn2o * dCOL;

   //  Taylor expanding loss freq (1/sec) of N2O
    const real_type Lfreq_n2o =  linoz_inputs.linoz_PmL_clim_ln2o           
            +  linoz_inputs.linoz_dPmL_dO3_ln2o    * dO3    
            +  linoz_inputs.linoz_dPmL_dn2o_ln2o   * dN2O   
            +  linoz_inputs.linoz_dPmL_dnoy_ln2o   * dNOY   
            +  linoz_inputs.linoz_dPmL_dch4_ln2o   * dCH4   
            +  linoz_inputs.linoz_dPmL_dh2o_ln2o   * dH2O   
            +  linoz_inputs.linoz_dPmL_dT_ln2o     * dTemp  
            +  linoz_inputs.linoz_dPmL_dO3col_ln2o * dCOL;  

//Taylar expanding production of noy divided by fn2o (so unit 1/sec)      
    const real_type Pfreq_noy =  linoz_inputs.linoz_PmL_clim_pnoy            
            +  linoz_inputs.linoz_dPmL_dO3_pnoy  * dO3    
            +  linoz_inputs.linoz_dPmL_dn2o_pnoy  * dN2O   
            +  linoz_inputs.linoz_dPmL_dnoy_pnoy  * dNOY   
            +  linoz_inputs.linoz_dPmL_dch4_pnoy  * dCH4   
            +  linoz_inputs.linoz_dPmL_dh2o_pnoy  * dH2O   
            +  linoz_inputs.linoz_dPmL_dT_pnoy   * dTemp  
            +  linoz_inputs.linoz_dPmL_dO3col_pnoy* dCOL;   

// Taylar expanding loss of noy divided by fn2o (so unit 1/sec)      
    const real_type Lfreq_noy =  linoz_inputs.linoz_PmL_clim_lnoy    
            +  linoz_inputs.linoz_dPmL_dO3_lnoy  * dO3  
            +  linoz_inputs.linoz_dPmL_dn2o_lnoy * dN2O 
            +  linoz_inputs.linoz_dPmL_dnoy_lnoy * dNOY 
            +  linoz_inputs.linoz_dPmL_dch4_lnoy * dCH4 
            +  linoz_inputs.linoz_dPmL_dh2o_lnoy * dH2O 
            +  linoz_inputs.linoz_dPmL_dT_lnoy   * dTemp
            +  linoz_inputs.linoz_dPmL_dO3col_lnoy* dCOL;

// Taylor expanding loss frequency of CH4
    const real_type Lfreq_ch4 =  linoz_inputs.linoz_PmL_clim_nch4             
            +  linoz_inputs.linoz_dPmL_dO3_nch4   * dO3    
            +  linoz_inputs.linoz_dPmL_dn2o_nch4  * dN2O   
            +  linoz_inputs.linoz_dPmL_dnoy_nch4  * dNOY   
            +  linoz_inputs.linoz_dPmL_dch4_nch4  * dCH4   
            +  linoz_inputs.linoz_dPmL_dh2o_nch4  * dH2O   
            +  linoz_inputs.linoz_dPmL_dT_nch4    * dTemp  
            +  linoz_inputs.linoz_dPmL_dO3col_nch4* dCOL;                            

  	team_invoke_detail(member,
  	          temp,
  	          pmid,
  	          delta_t, 
  	          lats, 
  	          psc_T,
  	          sza, 
  	          chlorine_loading, 
  	          // inputs
  	          o3_vmr, 
  	          n2o_vmr,
  	          noy_vmr,
  	          ch4_vmr, 
  	          h2o_vmr,
  	          linoz_inputs.linoz_o3_clim,
  	          linoz_inputs.linoz_dPmL_dO3X,  
  	          linoz_inputs.linoz_cariolle_psc, 
  	          PL_O3, 
  	          P_n2o, 
  	          Lfreq_n2o,
  	          Pfreq_noy,
  	          Lfreq_noy,
  	          Lfreq_ch4, 
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
   			  do3_linoz_psc
  	          );
  }
}; // 

} // namespace Impl
} // namespace TChem

#endif

