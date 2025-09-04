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
#ifndef __TCHEM_IMPL_MOSAIC_HPP__
#define __TCHEM_IMPL_MOSAIC_HPP__

/*
Code based on the Model for Simulating Aerosol
Interactions and Chemistry (MOSAIC)

Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008).
Model for simulating aerosol interactions and chemistry (MOSAIC).
Journal of Geophysical Research: Atmospheres.
*/

#include "TChem.hpp"
#include "TChem_Util.hpp"

namespace TChem {
namespace Impl {

template<typename DeviceType>
struct MosaicModelData {

  using device_type = DeviceType;

  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
  using real_type_2d_view_type = Tines::value_type_2d_view<real_type,device_type>;
  using real_type_3d_view_type = Tines::value_type_3d_view<real_type,device_type>;

  public:
    const ordinal_type iso4_a  = 0;   // <-> ih2so4_g
    const ordinal_type ino3_a  = 1;   // <-> ihno3_g
    const ordinal_type icl_a   = 2;   // <-> ihcl_g
    const ordinal_type inh4_a  = 3;   // <-> inh3_g
    const ordinal_type imsa_a  = 4;   // <-> imsa_g
    const ordinal_type iaro1_a = 5;   // <-> iaro1_g
    const ordinal_type iaro2_a = 6;   // <-> iaro2_g
    const ordinal_type ialk1_a = 7;   // <-> ialk1_g
    const ordinal_type iole1_a = 8;   // <-> iole1_g
    const ordinal_type iapi1_a = 9;   // <-> iapi1_g
    const ordinal_type iapi2_a = 10;  // <-> iapi2_g
    const ordinal_type ilim1_a = 11;  // <-> ilim1_g
    const ordinal_type ilim2_a = 12;  // <-> ilim2_g
    const ordinal_type ico3_a  = 13;  // <-> ico2_g
    const ordinal_type ina_a   = 14;
    const ordinal_type ica_a   = 15;
    const ordinal_type ioin_a  = 16;
    const ordinal_type ioc_a   = 17;
    const ordinal_type ibc_a   = 18;

    //electrolyte indices (used for water content calculations)
    // these indices are order sensitive

    const ordinal_type jnh4so4  = 0;  // soluble
    const ordinal_type jlvcite  = 1;  // soluble
    const ordinal_type jnh4hso4 = 2;  // soluble
    const ordinal_type jnh4msa  = 3;  // soluble: new
    const ordinal_type jnh4no3  = 4;  // soluble
    const ordinal_type jnh4cl   = 5;  // soluble
    const ordinal_type jna2so4  = 6;  // soluble
    const ordinal_type jna3hso4 = 7;  // soluble
    const ordinal_type jnahso4  = 8;  // soluble
    const ordinal_type jnamsa   = 9;  // soluble: new
    const ordinal_type jnano3   = 10; // soluble
    const ordinal_type jnacl    = 11; // soluble
    const ordinal_type jcano3   = 12; // soluble
    const ordinal_type jcacl2   = 13; // soluble
    const ordinal_type jcamsa2  = 14; // soluble
    const ordinal_type jh2so4   = 15; // soluble
    const ordinal_type jmsa     = 16; // soluble
    const ordinal_type jhno3    = 17; // soluble
    const ordinal_type jhcl     = 18; // soluble
    const ordinal_type jhhso4   = 19; // soluble
    const ordinal_type jcaso4   = 20; // insoluble
    const ordinal_type jcaco3   = 21; // insoluble
    const ordinal_type joc      = 22; // insoluble - part of naercomp
    const ordinal_type jbc      = 23; // insoluble - part of naercomp
    const ordinal_type join     = 24; // insoluble - part of naercomp
    const ordinal_type jaro1    = 25; // insoluble - part of naercomp
    const ordinal_type jaro2    = 26; // insoluble - part of naercomp
    const ordinal_type jalk1    = 27; // insoluble - part of naercomp
    const ordinal_type jole1    = 28; // insoluble - part of naercomp
    const ordinal_type japi1    = 29; // insoluble - part of naercomp
    const ordinal_type japi2    = 30; // insoluble - part of naercomp
    const ordinal_type jlim1    = 31; // insoluble - part of naercomp
    const ordinal_type jlim2    = 32; // insoluble - part of naercomp
    const ordinal_type jh2o     = 33; // water - part of naercomp

    const ordinal_type naer         = 19;
    const ordinal_type nelectrolyte = 22;
    const ordinal_type nsalt        = 15;
    const ordinal_type nsoluble     = 20;
    const ordinal_type ncation      =  4;
    const ordinal_type nanion       =  5;

    // FIXME: make enum
    const ordinal_type jsolid  = 1;
    const ordinal_type jliquid = 2;
    const ordinal_type jtotal  = 3;

    const ordinal_type no_aerosol = 0;
    const ordinal_type all_solid  = 1;
    const ordinal_type all_liquid = 2;
    const ordinal_type mixed      = 3;

    const ordinal_type soluble   = 1;
    const ordinal_type insoluble = 2;

    const ordinal_type jhyst_up = 1;
    const ordinal_type jhyst_lo = 0;

    const real_type mass_cutoff = 0.0;

    // polynomial coefficients for binary molality (used in ZSR)
       // for aw < 0.97
    real_type_2d_dual_view a_zsr;  // ("a_zsr", 6, nelectrolyte);
    real_type_1d_dual_view aw_min; // ("aw_min", nelectrolyte);
    real_type_1d_dual_view b_zsr;  // ("b_zsr", nelectrolyte)

    // local aerosol ions
    // cations
    const ordinal_type jc_h   = 0;
    const ordinal_type jc_nh4 = 1;
    const ordinal_type jc_na  = 2;
    const ordinal_type jc_ca  = 3;

    // anions
    const ordinal_type ja_hso4 = 0;
    const ordinal_type ja_so4  = 1;
    const ordinal_type ja_no3  = 2;
    const ordinal_type ja_cl   = 3;
    const ordinal_type ja_msa  = 4;
    // const ordinal_type ja_co3  = 6; (note: this is commented out in MOSAIC as well)

    // MESA parameters
    const ordinal_type Nmax_MESA = 80;
    const real_type rtol_mesa    =  0.01;

    // Temperature-dependent thermodynamic parameters
    // liquid-liquid
    const ordinal_type nrxn_aer_ll = 3;
    real_type_1d_dual_view Keq_298_ll; // ("Keq_298_ll_", nrxn_aer_ll);
    real_type_1d_dual_view Keq_a_ll; // ("Keq_a_ll_", nrxn_aer_ll);
    real_type_1d_dual_view Keq_b_ll; // ("Keq_b_ll_", nrxn_aer_ll);

    // parameters for MTEM mixing rule (Zaveri, Easter, and Wexler, 2005)
    real_type_3d_dual_view b_mtem; // ("b_mtem_", 6, nelectrolyte, nelectrolyte);

    // Coefficients for %MDRH(T) = d1 + d2*T + d3*T^2 + d4*T^3    (T in Kelvin)
    // valid Temperature Range: 240 - 320 K
    real_type_2d_dual_view d_mdrh; // ("d_mdrh_", 4, 63);

    // molecular weights for aerosol species
    real_type_1d_dual_view mw_aer_mac; // ("mw_aer_mac_", naer)

    // densities for aerosol species
    real_type_1d_dual_view dens_aer_mac; // ("dens_aer_mac_", naer)

    // molecular weights for electrolytes
    real_type_1d_dual_view mw_electrolyte; // ("mw_electrolyte_", nelectrolyte)

    // magnitude of charges on ions
    ordinal_type_1d_dual_view zc; // ("zc_", ncation)
    ordinal_type_1d_dual_view za; // ("za_", nanion)

    // molecular weights of ions [g/mol]
    real_type_1d_dual_view mw_c; // ("MW_c_", ncation)
    real_type_1d_dual_view mw_a; // ("MW_a_", anion)

    public:
      MosaicModelData() {

        a_zsr = real_type_2d_dual_view(do_not_init_tag("MosaicModelData::a_zsr"), 6, nelectrolyte);
        aw_min = real_type_1d_dual_view(do_not_init_tag("MosaicModelData::aw_min"), nelectrolyte);
        b_zsr = real_type_1d_dual_view(do_not_init_tag("MosaicModelData::b_zsr"), nelectrolyte);
        Keq_298_ll = real_type_1d_dual_view(do_not_init_tag("MosaicModelData::Keq_298_ll"), nrxn_aer_ll);
        Keq_a_ll = real_type_1d_dual_view(do_not_init_tag("MosaicModelData::Keq_a_ll"), nrxn_aer_ll);
        Keq_b_ll = real_type_1d_dual_view(do_not_init_tag("MosaicModelData::Keq_b_ll"), nrxn_aer_ll);
        b_mtem = real_type_3d_dual_view(do_not_init_tag("MosaicModelData::b_mtem"), 6, nelectrolyte, nelectrolyte);
        d_mdrh = real_type_2d_dual_view(do_not_init_tag("MosaicModelData::d_mdrh"), 63, 4);
        mw_aer_mac = real_type_1d_dual_view(do_not_init_tag("MosaicModelData::mw_aer_mac"), naer);
        dens_aer_mac = real_type_1d_dual_view(do_not_init_tag("MosaicModelData::dens_aer_mac"), naer);
        mw_electrolyte = real_type_1d_dual_view(do_not_init_tag("MosaicModelData::mw_electrolyte"), nelectrolyte);
        zc = ordinal_type_1d_dual_view(do_not_init_tag("MosaicModelData::zc"), ncation);
        za = ordinal_type_1d_dual_view(do_not_init_tag("MosaicModelData::za"), ncation);
        mw_c = real_type_1d_dual_view(do_not_init_tag("MosaicModelData::mw_c"), ncation);
        mw_a = real_type_1d_dual_view(do_not_init_tag("MosaicModelData::mw_a"), nanion);

        auto a_zsr_host = a_zsr.view_host();
        auto aw_min_host = aw_min.view_host();
        auto b_zsr_host = b_zsr.view_host();
        auto Keq_298_ll_host = Keq_298_ll.view_host();
        auto Keq_a_ll_host = Keq_a_ll.view_host();
        auto Keq_b_ll_host = Keq_b_ll.view_host();
        auto b_mtem_host = b_mtem.view_host();
        auto d_mdrh_host = d_mdrh.view_host();
        auto mw_aer_mac_host = mw_aer_mac.view_host();
        auto dens_aer_mac_host = dens_aer_mac.view_host();
        auto mw_electrolyte_host = mw_electrolyte.view_host();
        auto zc_host = zc.view_host();
        auto za_host = za.view_host();
        auto mw_c_host = mw_c.view_host();
        auto mw_a_host = mw_a.view_host();

        a_zsr_host(0, jnh4so4) =   1.30894;
        a_zsr_host(1, jnh4so4) =  -7.09922;
        a_zsr_host(2, jnh4so4) =  20.62831;
        a_zsr_host(3, jnh4so4) = -32.19965;
        a_zsr_host(4, jnh4so4) =  25.17026;
        a_zsr_host(5, jnh4so4) =  -7.81632;
        aw_min_host(jnh4so4)   =       0.1;

        a_zsr_host(0, jlvcite)  =   1.10725;
        a_zsr_host(1, jlvcite)  =  -5.17978;
        a_zsr_host(2, jlvcite)  =  12.29534;
        a_zsr_host(3, jlvcite)  = -16.32545;
        a_zsr_host(4, jlvcite)  =  11.29274;
        a_zsr_host(5, jlvcite)  =  -3.19164;
        aw_min_host(jlvcite)    =       0.1;

        a_zsr_host(0, jnh4hso4) =  1.15510;
        a_zsr_host(1, jnh4hso4) = -3.20815;
        a_zsr_host(2, jnh4hso4) =  2.71141;
        a_zsr_host(3, jnh4hso4) =  2.01155;
        a_zsr_host(4, jnh4hso4) = -4.71014;
        a_zsr_host(5, jnh4hso4) =  2.04616;
        aw_min_host(jnh4hso4)   =      0.1;

        a_zsr_host(0, jnh4msa) =  1.15510;
        a_zsr_host(1, jnh4msa) = -3.20815;
        a_zsr_host(2, jnh4msa) =  2.71141;
        a_zsr_host(3, jnh4msa) =  2.01155;
        a_zsr_host(4, jnh4msa) = -4.71014;
        a_zsr_host(5, jnh4msa) =  2.04616;
        aw_min_host(jnh4msa)   =      0.1;

        a_zsr_host(0, jnh4no3) =   0.43507;
        a_zsr_host(1, jnh4no3) =   6.38220;
        a_zsr_host(2, jnh4no3) = -30.19797;
        a_zsr_host(3, jnh4no3) =  53.36470;
        a_zsr_host(4, jnh4no3) = -43.44203;
        a_zsr_host(5, jnh4no3) =  13.46158;
        aw_min_host(jnh4no3)   =       0.1;

        a_zsr_host(0, jnh4cl) =  0.45309;
        a_zsr_host(1, jnh4cl) =  2.65606;
        a_zsr_host(2, jnh4cl) = -14.7730;
        a_zsr_host(3, jnh4cl) =  26.2936;
        a_zsr_host(4, jnh4cl) = -20.5735;
        a_zsr_host(5, jnh4cl) =  5.94255;
        aw_min_host(jnh4cl)   =      0.1;

        a_zsr_host(0, jnacl) =  0.42922;
        a_zsr_host(1, jnacl) = -1.17718;
        a_zsr_host(2, jnacl) =  2.80208;
        a_zsr_host(3, jnacl) = -4.51097;
        a_zsr_host(4, jnacl) =  3.76963;
        a_zsr_host(5, jnacl) = -1.31359;
        aw_min_host(jnacl)   =      0.1;

        a_zsr_host(0, jnano3) =   1.34966;
        a_zsr_host(1, jnano3) =  -5.20116;
        a_zsr_host(2, jnano3) =  11.49011;
        a_zsr_host(3, jnano3) = -14.41380;
        a_zsr_host(4, jnano3) =   9.07037;
        a_zsr_host(5, jnano3) =  -2.29769;
        aw_min_host(jnano3)   =       0.1;

        a_zsr_host(0, jna2so4) =  0.39888;
        a_zsr_host(1, jna2so4) = -1.27150;
        a_zsr_host(2, jna2so4) =  3.42792;
        a_zsr_host(3, jna2so4) = -5.92632;
        a_zsr_host(4, jna2so4) =  5.33351;
        a_zsr_host(5, jna2so4) = -1.96541;
        aw_min_host(jna2so4)   =      0.1;

        a_zsr_host(0, jna3hso4) =  0.31480;
        a_zsr_host(1, jna3hso4) = -1.01087;
        a_zsr_host(2, jna3hso4) =  2.44029;
        a_zsr_host(3, jna3hso4) = -3.66095;
        a_zsr_host(4, jna3hso4) =  2.77632;
        a_zsr_host(5, jna3hso4) = -0.86058;
        aw_min_host(jna3hso4)   =      0.1;

        a_zsr_host(0, jnahso4)  =   0.62764;
        a_zsr_host(1, jnahso4)  =  -1.63520;
        a_zsr_host(2, jnahso4)  =   4.62531;
        a_zsr_host(3, jnahso4)  = -10.06925;
        a_zsr_host(4, jnahso4)  =  10.33547;
        a_zsr_host(5, jnahso4)  =  -3.88729;
        aw_min_host(jnahso4)    =       0.1;

        a_zsr_host(0, jnamsa)  =   0.62764;
        a_zsr_host(1, jnamsa)  =  -1.63520;
        a_zsr_host(2, jnamsa)  =   4.62531;
        a_zsr_host(3, jnamsa)  = -10.06925;
        a_zsr_host(4, jnamsa)  =  10.33547;
        a_zsr_host(5, jnamsa)  =  -3.88729;
        aw_min_host(jnamsa)    =       0.1;

        a_zsr_host(0, jcano3) =  0.38895;
        a_zsr_host(1, jcano3) = -1.16013;
        a_zsr_host(2, jcano3) =  2.16819;
        a_zsr_host(3, jcano3) = -2.23079;
        a_zsr_host(4, jcano3) =  1.00268;
        a_zsr_host(5, jcano3) = -0.16923;
        aw_min_host(jcano3)   =      0.1;

        a_zsr_host(0, jcacl2) =  0.29891;
        a_zsr_host(1, jcacl2) = -1.31104;
        a_zsr_host(2, jcacl2) =  3.68759;
        a_zsr_host(3, jcacl2) = -5.81708;
        a_zsr_host(4, jcacl2) =  4.67520;
        a_zsr_host(5, jcacl2) = -1.53223;
        aw_min_host(jcacl2)   =      0.1;

        a_zsr_host(0, jh2so4) =  0.32751;
        a_zsr_host(1, jh2so4) = -1.00692;
        a_zsr_host(2, jh2so4) =  2.59750;
        a_zsr_host(3, jh2so4) = -4.40014;
        a_zsr_host(4, jh2so4) =  3.88212;
        a_zsr_host(5, jh2so4) = -1.39916;
        aw_min_host(jh2so4)   =      0.1;

        a_zsr_host(0, jmsa) =  0.32751;
        a_zsr_host(1, jmsa) = -1.00692;
        a_zsr_host(2, jmsa) =  2.59750;
        a_zsr_host(3, jmsa) = -4.40014;
        a_zsr_host(4, jmsa) =  3.88212;
        a_zsr_host(5, jmsa) = -1.39916;
        aw_min_host(jmsa)   =      0.1;

        a_zsr_host(0, jhhso4) =  0.32751;
        a_zsr_host(1, jhhso4) = -1.00692;
        a_zsr_host(2, jhhso4) =  2.59750;
        a_zsr_host(3, jhhso4) = -4.40014;
        a_zsr_host(4, jhhso4) =  3.88212;
        a_zsr_host(5, jhhso4) = -1.39916;
        aw_min_host(jhhso4)   =     1.0;

        a_zsr_host(0, jhno3) =   0.75876;
        a_zsr_host(1, jhno3) =  -3.31529;
        a_zsr_host(2, jhno3) =   9.26392;
        a_zsr_host(3, jhno3) = -14.89799;
        a_zsr_host(4, jhno3) =  12.08781;
        a_zsr_host(5, jhno3) =  -3.89958;
        aw_min_host(jhno3)   =       0.1;

        a_zsr_host(0, jhcl) =  0.31133;
        a_zsr_host(1, jhcl) = -0.79688;
        a_zsr_host(2, jhcl) =  1.93995;
        a_zsr_host(3, jhcl) = -3.31582;
        a_zsr_host(4, jhcl) =  2.93513;
        a_zsr_host(5, jhcl) = -1.07268;
        aw_min_host(jhcl)   =      0.1;

        a_zsr_host(0, jcaso4)  =  0.0;
        a_zsr_host(1, jcaso4)  =  0.0;
        a_zsr_host(2, jcaso4)  =  0.0;
        a_zsr_host(3, jcaso4)  =  0.0;
        a_zsr_host(4, jcaso4)  =  0.0;
        a_zsr_host(5, jcaso4)  =  0.0;
        aw_min_host(jcaso4)    =  1.0;

        a_zsr_host(0, jcamsa2) =  0.38895;
        a_zsr_host(1, jcamsa2) = -1.16013;
        a_zsr_host(2, jcamsa2) =  2.16819;
        a_zsr_host(3, jcamsa2) = -2.23079;
        a_zsr_host(4, jcamsa2) =  1.00268;
        a_zsr_host(5, jcamsa2) = -0.16923;
        aw_min_host(jcamsa2)   =      0.1;

        a_zsr_host(0, jcaco3)  =  0.0;
        a_zsr_host(1, jcaco3)  =  0.0;
        a_zsr_host(2, jcaco3)  =  0.0;
        a_zsr_host(3, jcaco3)  =  0.0;
        a_zsr_host(4, jcaco3)  =  0.0;
        a_zsr_host(5, jcaco3)  =  0.0;
        aw_min_host(jcaco3)    =  1.0;

        b_zsr_host(jnh4so4)  = 28.0811;
        b_zsr_host(jlvcite)  = 14.7178;
        b_zsr_host(jnh4hso4) = 29.4779;
        b_zsr_host(jnh4msa)  = 29.4779;
        b_zsr_host(jnh4no3)  = 33.4049;
        b_zsr_host(jnh4cl)   = 30.8888;
        b_zsr_host(jnacl)    = 29.8375;
        b_zsr_host(jnano3)   = 32.2756;
        b_zsr_host(jna2so4)  = 27.6889;
        b_zsr_host(jna3hso4) = 14.2184;
        b_zsr_host(jnahso4)  = 28.3367;
        b_zsr_host(jnamsa)   = 28.3367;
        b_zsr_host(jcano3)   = 18.3661;
        b_zsr_host(jcacl2)   = 20.8792;
        b_zsr_host(jh2so4)   = 26.7347;
        b_zsr_host(jhhso4)   = 26.7347;
        b_zsr_host(jhno3)    = 28.8257;
        b_zsr_host(jhcl)     = 27.7108;
        b_zsr_host(jmsa)     = 26.7347;
        b_zsr_host(jcaso4)   = 0.0;
        b_zsr_host(jcamsa2)  = 18.3661;
        b_zsr_host(jcaco3)   = 0.0;

        Keq_298_ll_host(0) = 1.0502e-2; // HSO4- <=> SO4= + H+
        Keq_298_ll_host(1) = 1.805e-5;  // NH3(l) + H2O = NH4+ + OH-
        Keq_298_ll_host(2) = 1.01e-14;  // H2O(l) <=> H+ + OH-

        ///////////////////////////////////////
        // (NH4)2SO4 in (NH4)2SO4
        b_mtem_host(0,jnh4so4,jnh4so4) =   -2.94685;
        b_mtem_host(1,jnh4so4,jnh4so4) =   17.3328;
        b_mtem_host(2,jnh4so4,jnh4so4) =  -64.8441;
        b_mtem_host(3,jnh4so4,jnh4so4) =  122.7070;
        b_mtem_host(4,jnh4so4,jnh4so4) = -114.4373;
        b_mtem_host(5,jnh4so4,jnh4so4) =   41.6811;

        // (NH4)2SO4 in NH4NO3
        b_mtem_host(0,jnh4so4,jnh4no3) = -2.7503;
        b_mtem_host(1,jnh4so4,jnh4no3) =  4.3806;
        b_mtem_host(2,jnh4so4,jnh4no3) = -1.1110;
        b_mtem_host(3,jnh4so4,jnh4no3) = -1.7005;
        b_mtem_host(4,jnh4so4,jnh4no3) = -4.4207;
        b_mtem_host(5,jnh4so4,jnh4no3) =  5.1990;

        //(NH4)2SO4 in NH4Cl(revised on 11/15/2003)
        b_mtem_host(0,jnh4so4,jnh4cl) =  -2.06952;
        b_mtem_host(1,jnh4so4,jnh4cl) =   7.1240;
        b_mtem_host(2,jnh4so4,jnh4cl) = -24.4274;
        b_mtem_host(3,jnh4so4,jnh4cl) =  51.1458;
        b_mtem_host(4,jnh4so4,jnh4cl) = -54.2056;
        b_mtem_host(5,jnh4so4,jnh4cl) =  22.0606;

        //(NH4)2SO4 in Na2SO4
        b_mtem_host(0,jnh4so4,jna2so4) =   -2.17361;
        b_mtem_host(1,jnh4so4,jna2so4) =   15.9919;
        b_mtem_host(2,jnh4so4,jna2so4) =  -69.0952;
        b_mtem_host(3,jnh4so4,jna2so4) =  139.8860;
        b_mtem_host(4,jnh4so4,jna2so4) = -134.9890;
        b_mtem_host(5,jnh4so4,jna2so4) =   49.8877;

        //(NH4)2SO4 in NaNO3
        b_mtem_host(0,jnh4so4,jnano3) =   -4.4370;
        b_mtem_host(1,jnh4so4,jnano3) =   24.0243;
        b_mtem_host(2,jnh4so4,jnano3) =  -76.2437;
        b_mtem_host(3,jnh4so4,jnano3) =  128.6660;
        b_mtem_host(4,jnh4so4,jnano3) = -110.0900;
        b_mtem_host(5,jnh4so4,jnano3) =   37.7414;

        //(NH4)2SO4 in NaCl
        b_mtem_host(0,jnh4so4,jnacl) =  -1.5394;
        b_mtem_host(1,jnh4so4,jnacl) =   5.8671;
        b_mtem_host(2,jnh4so4,jnacl) = -22.7726;
        b_mtem_host(3,jnh4so4,jnacl) =  47.0547;
        b_mtem_host(4,jnh4so4,jnacl) = -47.8266;
        b_mtem_host(5,jnh4so4,jnacl) =  18.8489;

        //(NH4)2SO4 in HNO3
        b_mtem_host(0,jnh4so4,jhno3) =  -0.35750;
        b_mtem_host(1,jnh4so4,jhno3) =  -3.82466;
        b_mtem_host(2,jnh4so4,jhno3) =   4.55462;
        b_mtem_host(3,jnh4so4,jhno3) =   5.05402;
        b_mtem_host(4,jnh4so4,jhno3) = -14.7476;
        b_mtem_host(5,jnh4so4,jhno3) =   8.8009;

        //(NH4)2SO4 in HCl
        b_mtem_host(0,jnh4so4,jhcl) =  -2.15146;
        b_mtem_host(1,jnh4so4,jhcl) =   5.50205;
        b_mtem_host(2,jnh4so4,jhcl) = -19.1476;
        b_mtem_host(3,jnh4so4,jhcl) =  39.1880;
        b_mtem_host(4,jnh4so4,jhcl) = -39.9460;
        b_mtem_host(5,jnh4so4,jhcl) =  16.0700;

        //(NH4)2SO4 in H2SO4
        b_mtem_host(0,jnh4so4,jh2so4) =  -2.52604;
        b_mtem_host(1,jnh4so4,jh2so4) =   9.76022;
        b_mtem_host(2,jnh4so4,jh2so4) = -35.2540;
        b_mtem_host(3,jnh4so4,jh2so4) =  71.2981;
        b_mtem_host(4,jnh4so4,jh2so4) = -71.8207;
        b_mtem_host(5,jnh4so4,jh2so4) =  28.0758;

        //(NH4)2SO4 in NH4HSO4
        b_mtem_host(0,jnh4so4,jnh4hso4) =  -4.13219;
        b_mtem_host(1,jnh4so4,jnh4hso4) =  13.8863;
        b_mtem_host(2,jnh4so4,jnh4hso4) = -34.5387;
        b_mtem_host(3,jnh4so4,jnh4hso4) =  56.5012;
        b_mtem_host(4,jnh4so4,jnh4hso4) = -51.8702;
        b_mtem_host(5,jnh4so4,jnh4hso4) =  19.6232;

        //(NH4)2SO4 in(NH4)3H_host(SO4)2
        b_mtem_host(0,jnh4so4,jlvcite) =  -2.53482;
        b_mtem_host(1,jnh4so4,jlvcite) =  12.3333;
        b_mtem_host(2,jnh4so4,jlvcite) = -46.1020;
        b_mtem_host(3,jnh4so4,jlvcite) =  90.4775;
        b_mtem_host(4,jnh4so4,jlvcite) = -88.1254;
        b_mtem_host(5,jnh4so4,jlvcite) =  33.4715;

        //(NH4)2SO4 in NaHSO4
        b_mtem_host(0,jnh4so4,jnahso4) =   -3.23425;
        b_mtem_host(1,jnh4so4,jnahso4) =   18.7842;
        b_mtem_host(2,jnh4so4,jnahso4) =  -78.7807;
        b_mtem_host(3,jnh4so4,jnahso4) =  161.517;
        b_mtem_host(4,jnh4so4,jnahso4) = -154.940;
        b_mtem_host(5,jnh4so4,jnahso4) =   56.2252;

        //(NH4)2SO4  in Na3H_host(SO4)2
        b_mtem_host(0,jnh4so4,jna3hso4) =  -1.25316;
        b_mtem_host(1,jnh4so4,jna3hso4) =   7.40960;
        b_mtem_host(2,jnh4so4,jna3hso4) = -34.8929;
        b_mtem_host(3,jnh4so4,jna3hso4) =  72.8853;
        b_mtem_host(4,jnh4so4,jna3hso4) = -72.4503;
        b_mtem_host(5,jnh4so4,jna3hso4) =  27.7706;

        ///////////////////////////////////////

        ///////////////////////////////////////
        // NH4NO3 in(NH4)2SO4
        b_mtem_host(0,jnh4no3,jnh4so4) =   -3.5201;
        b_mtem_host(1,jnh4no3,jnh4so4) =   21.6584;
        b_mtem_host(2,jnh4no3,jnh4so4) =  -72.1499;
        b_mtem_host(3,jnh4no3,jnh4so4) =  126.7000;
        b_mtem_host(4,jnh4no3,jnh4so4) = -111.4550;
        b_mtem_host(5,jnh4no3,jnh4so4) =   38.5677;

        // NH4NO3 in NH4NO3
        b_mtem_host(0,jnh4no3,jnh4no3) =  -2.2630;
        b_mtem_host(1,jnh4no3,jnh4no3) =  -0.1518;
        b_mtem_host(2,jnh4no3,jnh4no3) =  17.0898;
        b_mtem_host(3,jnh4no3,jnh4no3) = -36.7832;
        b_mtem_host(4,jnh4no3,jnh4no3) =  29.8407;
        b_mtem_host(5,jnh4no3,jnh4no3) =  -7.9314;

        // NH4NO3 in NH4Cl(revised on 11/15/2003)
        b_mtem_host(0,jnh4no3,jnh4cl) =  -1.3851;
        b_mtem_host(1,jnh4no3,jnh4cl) =  -0.4462;
        b_mtem_host(2,jnh4no3,jnh4cl) =   8.4567;
        b_mtem_host(3,jnh4no3,jnh4cl) = -11.5988;
        b_mtem_host(4,jnh4no3,jnh4cl) =   2.9802;
        b_mtem_host(5,jnh4no3,jnh4cl) =   1.8132;

        // NH4NO3 in Na2SO4
        b_mtem_host(0,jnh4no3,jna2so4) =  -1.7602;
        b_mtem_host(1,jnh4no3,jna2so4) =  10.4044;
        b_mtem_host(2,jnh4no3,jna2so4) = -35.5894;
        b_mtem_host(3,jnh4no3,jna2so4) =  64.3584;
        b_mtem_host(4,jnh4no3,jna2so4) = -57.8931;
        b_mtem_host(5,jnh4no3,jna2so4) =  20.2141;

        // NH4NO3 in NaNO3
        b_mtem_host(0,jnh4no3,jnano3) =  -3.24346;
        b_mtem_host(1,jnh4no3,jnano3) =  16.2794;
        b_mtem_host(2,jnh4no3,jnano3) = -48.7601;
        b_mtem_host(3,jnh4no3,jnano3) =  79.2246;
        b_mtem_host(4,jnh4no3,jnano3) = -65.8169;
        b_mtem_host(5,jnh4no3,jnano3) =  22.1500;

        // NH4NO3 in NaCl
        b_mtem_host(0,jnh4no3,jnacl) =  -1.75658;
        b_mtem_host(1,jnh4no3,jnacl) =   7.71384;
        b_mtem_host(2,jnh4no3,jnacl) = -22.7984;
        b_mtem_host(3,jnh4no3,jnacl) =  39.1532;
        b_mtem_host(4,jnh4no3,jnacl) = -34.6165;
        b_mtem_host(5,jnh4no3,jnacl) =  12.1283;

        // NH4NO3 in Ca_host(NO3)2
        b_mtem_host(0,jnh4no3,jcano3) =  -0.97178;
        b_mtem_host(1,jnh4no3,jcano3) =   6.61964;
        b_mtem_host(2,jnh4no3,jcano3) = -26.2353;
        b_mtem_host(3,jnh4no3,jcano3) =  50.5259;
        b_mtem_host(4,jnh4no3,jcano3) = -47.6586;
        b_mtem_host(5,jnh4no3,jcano3) =  17.5074;

        // NH4NO3 in CaCl2 added on 12/22/2003
        b_mtem_host(0,jnh4no3,jcacl2) =  -0.41515;
        b_mtem_host(1,jnh4no3,jcacl2) =   6.44101;
        b_mtem_host(2,jnh4no3,jcacl2) = -26.4473;
        b_mtem_host(3,jnh4no3,jcacl2) =  49.0718;
        b_mtem_host(4,jnh4no3,jcacl2) = -44.2631;
        b_mtem_host(5,jnh4no3,jcacl2) =  15.3771;

        // NH4NO3 in HNO3
        b_mtem_host(0,jnh4no3,jhno3) =  -1.20644;
        b_mtem_host(1,jnh4no3,jhno3) =   5.70117;
        b_mtem_host(2,jnh4no3,jhno3) = -18.2783;
        b_mtem_host(3,jnh4no3,jhno3) =  31.7199;
        b_mtem_host(4,jnh4no3,jhno3) = -27.8703;
        b_mtem_host(5,jnh4no3,jhno3) =   9.7299;

        // NH4NO3 in HCl
        b_mtem_host(0,jnh4no3,jhcl) =  -0.680862;
        b_mtem_host(1,jnh4no3,jhcl) =   3.59456;
        b_mtem_host(2,jnh4no3,jhcl) = -10.7969;
        b_mtem_host(3,jnh4no3,jhcl) =  17.8434;
        b_mtem_host(4,jnh4no3,jhcl) = -15.3165;
        b_mtem_host(5,jnh4no3,jhcl) =   5.17123;

        ///////////////////////////////////////

        ///////////////////////////////////////
        // NH4Cl in(NH4)2SO4
        b_mtem_host(0,jnh4cl,jnh4so4) =   -2.8850;
        b_mtem_host(1,jnh4cl,jnh4so4) =   20.6970;
        b_mtem_host(2,jnh4cl,jnh4so4) =  -70.6810;
        b_mtem_host(3,jnh4cl,jnh4so4) =  124.3690;
        b_mtem_host(4,jnh4cl,jnh4so4) = -109.2880;
        b_mtem_host(5,jnh4cl,jnh4so4) =   37.5831;

        // NH4Cl in NH4NO3
        b_mtem_host(0,jnh4cl,jnh4no3) =  -1.9386;
        b_mtem_host(1,jnh4cl,jnh4no3) =   1.3238;
        b_mtem_host(2,jnh4cl,jnh4no3) =  11.8500;
        b_mtem_host(3,jnh4cl,jnh4no3) = -28.1168;
        b_mtem_host(4,jnh4cl,jnh4no3) =  21.8543;
        b_mtem_host(5,jnh4cl,jnh4no3) =  -5.1671;

        // NH4Cl in NH4Cl(revised on 11/15/2003)
        b_mtem_host(0,jnh4cl,jnh4cl) = -0.9559;
        b_mtem_host(1,jnh4cl,jnh4cl) =  0.8121;
        b_mtem_host(2,jnh4cl,jnh4cl) =  4.3644;
        b_mtem_host(3,jnh4cl,jnh4cl) = -8.9258;
        b_mtem_host(4,jnh4cl,jnh4cl) =  4.2362;
        b_mtem_host(5,jnh4cl,jnh4cl) =  0.2891;

        // NH4Cl in Na2SO4
        b_mtem_host(0,jnh4cl,jna2so4) =  0.0377;
        b_mtem_host(1,jnh4cl,jna2so4) =  6.0752;
        b_mtem_host(2,jnh4cl,jna2so4) = -30.8641;
        b_mtem_host(3,jnh4cl,jna2so4) =  63.3095;
        b_mtem_host(4,jnh4cl,jna2so4) = -61.0070;
        b_mtem_host(5,jnh4cl,jna2so4) =  22.1734;

        // NH4Cl in NaNO3
        b_mtem_host(0,jnh4cl,jnano3) =  -1.8336;
        b_mtem_host(1,jnh4cl,jnano3) =  12.8160;
        b_mtem_host(2,jnh4cl,jnano3) = -42.3388;
        b_mtem_host(3,jnh4cl,jnano3) =  71.1816;
        b_mtem_host(4,jnh4cl,jnano3) = -60.5708;
        b_mtem_host(5,jnh4cl,jnano3) =  20.5853;

        // NH4Cl in NaCl
        b_mtem_host(0,jnh4cl,jnacl) =  -0.1429;
        b_mtem_host(1,jnh4cl,jnacl) =   2.3561;
        b_mtem_host(2,jnh4cl,jnacl) = -10.4425;
        b_mtem_host(3,jnh4cl,jnacl) =  20.8951;
        b_mtem_host(4,jnh4cl,jnacl) = -20.7739;
        b_mtem_host(5,jnh4cl,jnacl) =   7.9355;

        // NH4Cl in Ca_host(NO3)2
        b_mtem_host(0,jnh4cl,jcano3) =   0.76235;
        b_mtem_host(1,jnh4cl,jcano3) =   3.08323;
        b_mtem_host(2,jnh4cl,jcano3) = -23.6772;
        b_mtem_host(3,jnh4cl,jcano3) =  53.7415;
        b_mtem_host(4,jnh4cl,jcano3) = -55.4043;
        b_mtem_host(5,jnh4cl,jcano3) =  21.2944;

        // NH4Cl in CaCl2(revised on 11/27/2003)
        b_mtem_host(0,jnh4cl,jcacl2) =   1.13864;
        b_mtem_host(1,jnh4cl,jcacl2) =  -0.340539;
        b_mtem_host(2,jnh4cl,jcacl2) =  -8.67025;
        b_mtem_host(3,jnh4cl,jcacl2) =  22.8008;
        b_mtem_host(4,jnh4cl,jcacl2) = -24.5181;
        b_mtem_host(5,jnh4cl,jcacl2) =   9.3663;

        // NH4Cl in HNO3
        b_mtem_host(0,jnh4cl,jhno3) =   2.42532;
        b_mtem_host(1,jnh4cl,jhno3) = -14.1755;
        b_mtem_host(2,jnh4cl,jhno3) =  38.804;
        b_mtem_host(3,jnh4cl,jhno3) = -58.2437;
        b_mtem_host(4,jnh4cl,jhno3) =  43.5431;
        b_mtem_host(5,jnh4cl,jhno3) = -12.5824;

        // NH4Cl in HCl
        b_mtem_host(0,jnh4cl,jhcl) =  0.330337;
        b_mtem_host(1,jnh4cl,jhcl) =  0.0778934;
        b_mtem_host(2,jnh4cl,jhcl) = -2.30492;
        b_mtem_host(3,jnh4cl,jhcl) =  4.73003;
        b_mtem_host(4,jnh4cl,jhcl) = -4.80849;
        b_mtem_host(5,jnh4cl,jhcl) =  1.78866;

        ///////////////////////////////////////

        ///////////////////////////////////////
        // Na2SO4 in(NH4)2SO4
        b_mtem_host(0,jna2so4,jnh4so4) =   -2.6982;
        b_mtem_host(1,jna2so4,jnh4so4) =   22.9875;
        b_mtem_host(2,jna2so4,jnh4so4) =  -98.9840;
        b_mtem_host(3,jna2so4,jnh4so4) =  198.0180;
        b_mtem_host(4,jna2so4,jnh4so4) = -188.7270;
        b_mtem_host(5,jna2so4,jnh4so4) =   69.0548;

        // Na2SO4 in NH4NO3
        b_mtem_host(0,jna2so4,jnh4no3) =  -2.4844;
        b_mtem_host(1,jna2so4,jnh4no3) =   6.5420;
        b_mtem_host(2,jna2so4,jnh4no3) =  -9.8998;
        b_mtem_host(3,jna2so4,jnh4no3) =  11.3884;
        b_mtem_host(4,jna2so4,jnh4no3) = -13.6842;
        b_mtem_host(5,jna2so4,jnh4no3) =   7.7411;

        // Na2SO4 in NH4Cl(revised on 11/15/2003)
        b_mtem_host(0,jna2so4,jnh4cl) =  -1.3325;
        b_mtem_host(1,jna2so4,jnh4cl) =  13.0406;
        b_mtem_host(2,jna2so4,jnh4cl) = -56.1935;
        b_mtem_host(3,jna2so4,jnh4cl) = 107.1170;
        b_mtem_host(4,jna2so4,jnh4cl) = -97.3721;
        b_mtem_host(5,jna2so4,jnh4cl) =  34.3763;

        // Na2SO4 in Na2SO4
        b_mtem_host(0,jna2so4,jna2so4) =  -1.2832;
        b_mtem_host(1,jna2so4,jna2so4) =  12.8526;
        b_mtem_host(2,jna2so4,jna2so4) =  -62.2087;
        b_mtem_host(3,jna2so4,jna2so4) =  130.3876;
        b_mtem_host(4,jna2so4,jna2so4) = -128.2627;
        b_mtem_host(5,jna2so4,jna2so4) =  48.0340;

        // Na2SO4 in NaNO3
        b_mtem_host(0,jna2so4,jnano3) =   -3.5384;
        b_mtem_host(1,jna2so4,jnano3) =   21.3758;
        b_mtem_host(2,jna2so4,jnano3) =  -70.7638;
        b_mtem_host(3,jna2so4,jnano3) =  121.1580;
        b_mtem_host(4,jna2so4,jnano3) = -104.6230;
        b_mtem_host(5,jna2so4,jnano3) =   36.0557;

        // Na2SO4 in NaCl
        b_mtem_host(0,jna2so4,jnacl) =   0.2175;
        b_mtem_host(1,jna2so4,jnacl) =  -0.5648;
        b_mtem_host(2,jna2so4,jnacl) =  -8.0288;
        b_mtem_host(3,jna2so4,jnacl) =  25.9734;
        b_mtem_host(4,jna2so4,jnacl) = -32.3577;
        b_mtem_host(5,jna2so4,jnacl) =  14.3924;

        // Na2SO4 in HNO3
        b_mtem_host(0,jna2so4,jhno3) =  -0.309617;
        b_mtem_host(1,jna2so4,jhno3) =  -1.82899;
        b_mtem_host(2,jna2so4,jhno3) =  -1.5505;
        b_mtem_host(3,jna2so4,jhno3) =  13.3847;
        b_mtem_host(4,jna2so4,jhno3) = -20.1284;
        b_mtem_host(5,jna2so4,jhno3) =   9.93163;

        // Na2SO4 in HCl
        b_mtem_host(0,jna2so4,jhcl) =  -0.259455;
        b_mtem_host(1,jna2so4,jhcl) =  -0.819366;
        b_mtem_host(2,jna2so4,jhcl) =  -4.28964;
        b_mtem_host(3,jna2so4,jhcl) =  16.4305;
        b_mtem_host(4,jna2so4,jhcl) = -21.8546;
        b_mtem_host(5,jna2so4,jhcl) =  10.3044;

        // Na2SO4 in H2SO4
        b_mtem_host(0,jna2so4,jh2so4) =  -1.84257;
        b_mtem_host(1,jna2so4,jh2so4) =   7.85788;
        b_mtem_host(2,jna2so4,jh2so4) = -29.9275;
        b_mtem_host(3,jna2so4,jh2so4) =  61.7515;
        b_mtem_host(4,jna2so4,jh2so4) = -63.2308;
        b_mtem_host(5,jna2so4,jh2so4) =  24.9542;

        // Na2SO4 in NH4HSO4
        b_mtem_host(0,jna2so4,jnh4hso4) =  -1.05891;
        b_mtem_host(1,jna2so4,jnh4hso4) =   2.84831;
        b_mtem_host(2,jna2so4,jnh4hso4) = -21.1827;
        b_mtem_host(3,jna2so4,jnh4hso4) =  57.5175;
        b_mtem_host(4,jna2so4,jnh4hso4) = -64.8120;
        b_mtem_host(5,jna2so4,jnh4hso4) =  26.1986;

        // Na2SO4 in(NH4)3H_host(SO4)2
        b_mtem_host(0,jna2so4,jlvcite) =  -1.16584;
        b_mtem_host(1,jna2so4,jlvcite) =   8.50075;
        b_mtem_host(2,jna2so4,jlvcite) = -44.3420;
        b_mtem_host(3,jna2so4,jlvcite) =  97.3974;
        b_mtem_host(4,jna2so4,jlvcite) = -98.4549;
        b_mtem_host(5,jna2so4,jlvcite) =  37.6104;

        // Na2SO4 in NaHSO4
        b_mtem_host(0,jna2so4,jnahso4) =  -1.95805;
        b_mtem_host(1,jna2so4,jnahso4) =   6.62417;
        b_mtem_host(2,jna2so4,jnahso4) = -31.8072;
        b_mtem_host(3,jna2so4,jnahso4) =  77.8603;
        b_mtem_host(4,jna2so4,jnahso4) = -84.6458;
        b_mtem_host(5,jna2so4,jnahso4) =  33.4963;

        // Na2SO4 in Na3H_host(SO4)2
        b_mtem_host(0,jna2so4,jna3hso4) =  -0.36045;
        b_mtem_host(1,jna2so4,jna3hso4) =   3.55223;
        b_mtem_host(2,jna2so4,jna3hso4) = -24.0327;
        b_mtem_host(3,jna2so4,jna3hso4) =  54.4879;
        b_mtem_host(4,jna2so4,jna3hso4) = -56.6531;
        b_mtem_host(5,jna2so4,jna3hso4) =  22.4956;

        ///////////////////////////////////////

        ///////////////////////////////////////
        // NaNO3 in(NH4)2SO4
        b_mtem_host(0,jnano3,jnh4so4) =   -2.5888;
        b_mtem_host(1,jnano3,jnh4so4) =   17.6192;
        b_mtem_host(2,jnano3,jnh4so4) =  -63.2183;
        b_mtem_host(3,jnano3,jnh4so4) =  115.3520;
        b_mtem_host(4,jnano3,jnh4so4) = -104.0860;
        b_mtem_host(5,jnano3,jnh4so4) =   36.7390;

        // NaNO3 in NH4NO3
        b_mtem_host(0,jnano3,jnh4no3) =  -2.0669;
        b_mtem_host(1,jnano3,jnh4no3) =   1.4792;
        b_mtem_host(2,jnano3,jnh4no3) =  10.5261;
        b_mtem_host(3,jnano3,jnh4no3) = -27.0987;
        b_mtem_host(4,jnano3,jnh4no3) =  23.0591;
        b_mtem_host(5,jnano3,jnh4no3) =  -6.0938;

        // NaNO3 in NH4Cl(revised on 11/15/2003)
        b_mtem_host(0,jnano3,jnh4cl) =  -0.8325;
        b_mtem_host(1,jnano3,jnh4cl) =   3.9933;
        b_mtem_host(2,jnano3,jnh4cl) = -15.3789;
        b_mtem_host(3,jnano3,jnh4cl) =  30.4050;
        b_mtem_host(4,jnano3,jnh4cl) = -29.4204;
        b_mtem_host(5,jnano3,jnh4cl) =  11.0597;

        // NaNO3 in Na2SO4
        b_mtem_host(0,jnano3,jna2so4) =  -1.1233;
        b_mtem_host(1,jnano3,jna2so4) =   8.3998;
        b_mtem_host(2,jnano3,jna2so4) = -31.9002;
        b_mtem_host(3,jnano3,jna2so4) =  60.1450;
        b_mtem_host(4,jnano3,jna2so4) = -55.5503;
        b_mtem_host(5,jnano3,jna2so4) =  19.7757;

        // NaNO3 in NaNO3
        b_mtem_host(0,jnano3,jnano3) =  -2.5386;
        b_mtem_host(1,jnano3,jnano3) =  13.9039;
        b_mtem_host(2,jnano3,jnano3) = -42.8467;
        b_mtem_host(3,jnano3,jnano3) =  69.7442;
        b_mtem_host(4,jnano3,jnano3) = -57.8988;
        b_mtem_host(5,jnano3,jnano3) =  19.4635;

        // NaNO3 in NaCl
        b_mtem_host(0,jnano3,jnacl) =  -0.4351;
        b_mtem_host(1,jnano3,jnacl) =   2.8311;
        b_mtem_host(2,jnano3,jnacl) = -11.4485;
        b_mtem_host(3,jnano3,jnacl) =  22.7201;
        b_mtem_host(4,jnano3,jnacl) = -22.4228;
        b_mtem_host(5,jnano3,jnacl) =   8.5792;

        // NaNO3 in Ca_host(NO3)2
        b_mtem_host(0,jnano3,jcano3) = -0.72060;
        b_mtem_host(1,jnano3,jcano3) =  5.64915;
        b_mtem_host(2,jnano3,jcano3) = -23.5020;
        b_mtem_host(3,jnano3,jcano3) =  46.0078;
        b_mtem_host(4,jnano3,jcano3) = -43.8075;
        b_mtem_host(5,jnano3,jcano3) =  16.1652;

        // NaNO3 in CaCl2
        b_mtem_host(0,jnano3,jcacl2) =   0.003928;
        b_mtem_host(1,jnano3,jcacl2) =   3.54724;
        b_mtem_host(2,jnano3,jcacl2) = -18.6057;
        b_mtem_host(3,jnano3,jcacl2) =  38.1445;
        b_mtem_host(4,jnano3,jcacl2) = -36.7745;
        b_mtem_host(5,jnano3,jcacl2) =  13.4529;

        // NaNO3 in HNO3
        b_mtem_host(0,jnano3,jhno3) =  -1.1712;
        b_mtem_host(1,jnano3,jhno3) =   7.20907;
        b_mtem_host(2,jnano3,jhno3) = -22.9215;
        b_mtem_host(3,jnano3,jhno3) =  38.1257;
        b_mtem_host(4,jnano3,jhno3) = -32.0759;
        b_mtem_host(5,jnano3,jhno3) =  10.6443;

        // NaNO3 in HCl
        b_mtem_host(0,jnano3,jhcl) =  0.738022;
        b_mtem_host(1,jnano3,jhcl) = -1.14313;
        b_mtem_host(2,jnano3,jhcl) =  0.32251;
        b_mtem_host(3,jnano3,jhcl) =  0.838679;
        b_mtem_host(4,jnano3,jhcl) = -1.81747;
        b_mtem_host(5,jnano3,jhcl) =  0.873986;
        ///////////////////////////////////////

        ///////////////////////////////////////
        // NaCl in(NH4)2SO4
        b_mtem_host(0,jnacl,jnh4so4) =   -1.9525;
        b_mtem_host(1,jnacl,jnh4so4) =   16.6433;
        b_mtem_host(2,jnacl,jnh4so4) =  -61.7090;
        b_mtem_host(3,jnacl,jnh4so4) =  112.9910;
        b_mtem_host(4,jnacl,jnh4so4) = -101.9370;
        b_mtem_host(5,jnacl,jnh4so4) =   35.7760;

        // NaCl in NH4NO3
        b_mtem_host(0,jnacl,jnh4no3) =  -1.7525;
        b_mtem_host(1,jnacl,jnh4no3) =   3.0713;
        b_mtem_host(2,jnacl,jnh4no3) =   4.8063;
        b_mtem_host(3,jnacl,jnh4no3) = -17.5334;
        b_mtem_host(4,jnacl,jnh4no3) =  14.2872;
        b_mtem_host(5,jnacl,jnh4no3) =  -3.0690;

        // NaCl in NH4Cl(revised on 11/15/2003)
        b_mtem_host(0,jnacl,jnh4cl) =  -0.4021;
        b_mtem_host(1,jnacl,jnh4cl) =   5.2399;
        b_mtem_host(2,jnacl,jnh4cl) = -19.4278;
        b_mtem_host(3,jnacl,jnh4cl) =  33.0027;
        b_mtem_host(4,jnacl,jnh4cl) = -28.1020;
        b_mtem_host(5,jnacl,jnh4cl) =   9.5159;

        // NaCl in Na2SO4
        b_mtem_host(0,jnacl,jna2so4) =   0.6692;
        b_mtem_host(1,jnacl,jna2so4) =   4.1207;
        b_mtem_host(2,jnacl,jna2so4) = -27.3314;
        b_mtem_host(3,jnacl,jna2so4) =  59.3112;
        b_mtem_host(4,jnacl,jna2so4) = -58.7998;
        b_mtem_host(5,jnacl,jna2so4) =  21.7674;

        // NaCl in NaNO3
        b_mtem_host(0,jnacl,jnano3) =  -1.17444;
        b_mtem_host(1,jnacl,jnano3) =  10.9927;
        b_mtem_host(2,jnacl,jnano3) = -38.9013;
        b_mtem_host(3,jnacl,jnano3) =  66.8521;
        b_mtem_host(4,jnacl,jnano3) = -57.6564;
        b_mtem_host(5,jnacl,jnano3) =  19.7296;

        // NaCl in NaCl
        b_mtem_host(0,jnacl,jnacl) =  1.17679;
        b_mtem_host(1,jnacl,jnacl) = -2.5061;
        b_mtem_host(2,jnacl,jnacl) =  0.8508;
        b_mtem_host(3,jnacl,jnacl) =  4.4802;
        b_mtem_host(4,jnacl,jnacl) = -8.4945;
        b_mtem_host(5,jnacl,jnacl) =  4.3182;

        // NaCl in Ca_host(NO3)2
        b_mtem_host(0,jnacl,jcano3) =   1.01450;
        b_mtem_host(1,jnacl,jcano3) =   2.10260;
        b_mtem_host(2,jnacl,jcano3) = -20.9036;
        b_mtem_host(3,jnacl,jcano3) =  49.1481;
        b_mtem_host(4,jnacl,jcano3) = -51.4867;
        b_mtem_host(5,jnacl,jcano3) =  19.9301;

        // NaCl in CaCl2(PSC92: revised on 11/27/2003)
        b_mtem_host(0,jnacl,jcacl2) =   1.55463;
        b_mtem_host(1,jnacl,jcacl2) =  -3.20122;
        b_mtem_host(2,jnacl,jcacl2) =  -0.957075;
        b_mtem_host(3,jnacl,jcacl2) =  12.103;
        b_mtem_host(4,jnacl,jcacl2) = -17.221;
        b_mtem_host(5,jnacl,jcacl2) =   7.50264;

        // NaCl in HNO3
        b_mtem_host(0,jnacl,jhno3) =   2.46187;
        b_mtem_host(1,jnacl,jhno3) = -12.6845;
        b_mtem_host(2,jnacl,jhno3) =  34.2383;
        b_mtem_host(3,jnacl,jhno3) = -51.9992;
        b_mtem_host(4,jnacl,jhno3) =  39.4934;
        b_mtem_host(5,jnacl,jhno3) = -11.7247;

        // NaCl in HCl
        b_mtem_host(0,jnacl,jhcl) =   1.74915;
        b_mtem_host(1,jnacl,jhcl) =  -4.65768;
        b_mtem_host(2,jnacl,jhcl) =   8.80287;
        b_mtem_host(3,jnacl,jhcl) = -12.2503;
        b_mtem_host(4,jnacl,jhcl) =   8.668751;
        b_mtem_host(5,jnacl,jhcl) =  -2.50158;

        ///////////////////////////////////////

        ///////////////////////////////////////
        // Ca_host(NO3)2 in NH4NO3
        b_mtem_host(0,jcano3,jnh4no3) =  -1.86260;
        b_mtem_host(1,jcano3,jnh4no3) =  11.6178;
        b_mtem_host(2,jcano3,jnh4no3) = -30.9069;
        b_mtem_host(3,jcano3,jnh4no3) =  41.7578;
        b_mtem_host(4,jcano3,jnh4no3) = -33.7338;
        b_mtem_host(5,jcano3,jnh4no3) =  12.7541;

        // Ca_host(NO3)2 in NH4Cl(revised on 11/15/2003)
        b_mtem_host(0,jcano3,jnh4cl) =   -1.1798;
        b_mtem_host(1,jcano3,jnh4cl) =   25.9608;
        b_mtem_host(2,jcano3,jnh4cl) =  -98.9373;
        b_mtem_host(3,jcano3,jnh4cl) =  160.2300;
        b_mtem_host(4,jcano3,jnh4cl) = -125.9540;
        b_mtem_host(5,jcano3,jnh4cl) =   39.5130;

        // Ca_host(NO3)2 in NaNO3
        b_mtem_host(0,jcano3,jnano3) =  -1.44384;
        b_mtem_host(1,jcano3,jnano3) =  13.6044;
        b_mtem_host(2,jcano3,jnano3) = -54.4300;
        b_mtem_host(3,jcano3,jnano3) = 100.582;
        b_mtem_host(4,jcano3,jnano3) = -91.2364;
        b_mtem_host(5,jcano3,jnano3) =  32.5970;

        // Ca_host(NO3)2 in NaCl
        b_mtem_host(0,jcano3,jnacl) =  -0.099114;
        b_mtem_host(1,jcano3,jnacl) =   2.84091;
        b_mtem_host(2,jcano3,jnacl) = -16.9229;
        b_mtem_host(3,jcano3,jnacl) =  37.4839;
        b_mtem_host(4,jcano3,jnacl) = -39.5132;
        b_mtem_host(5,jcano3,jnacl) =  15.8564;

        // Ca_host(NO3)2 in Ca_host(NO3)2
        b_mtem_host(0,jcano3,jcano3) =   0.055116;
        b_mtem_host(1,jcano3,jcano3) =   4.58610;
        b_mtem_host(2,jcano3,jcano3) = -27.6629;
        b_mtem_host(3,jcano3,jcano3) =  60.8288;
        b_mtem_host(4,jcano3,jcano3) = -61.4988;
        b_mtem_host(5,jcano3,jcano3) =  23.3136;

        // Ca_host(NO3)2 in CaCl2(PSC92: revised on 11/27/2003)
        b_mtem_host(0,jcano3,jcacl2) =   1.57155;
        b_mtem_host(1,jcano3,jcacl2) =  -3.18486;
        b_mtem_host(2,jcano3,jcacl2) =  -3.35758;
        b_mtem_host(3,jcano3,jcacl2) =  18.7501;
        b_mtem_host(4,jcano3,jcacl2) = -24.5604;
        b_mtem_host(5,jcano3,jcacl2) =  10.3798;

        // Ca_host(NO3)2 in HNO3
        b_mtem_host(0,jcano3,jhno3) =  1.04446;
        b_mtem_host(1,jcano3,jhno3) = -3.19066;
        b_mtem_host(2,jcano3,jhno3) =  2.44714;
        b_mtem_host(3,jcano3,jhno3) =  2.07218;
        b_mtem_host(4,jcano3,jhno3) = -6.43949;
        b_mtem_host(5,jcano3,jhno3) =  3.66471;

        // Ca_host(NO3)2 in HCl
        b_mtem_host(0,jcano3,jhcl) =  1.05723;
        b_mtem_host(1,jcano3,jhcl) = -1.46826;
        b_mtem_host(2,jcano3,jhcl) = -1.0713;
        b_mtem_host(3,jcano3,jhcl) =  4.64439;
        b_mtem_host(4,jcano3,jhcl) = -6.32402;
        b_mtem_host(5,jcano3,jhcl) =  2.78202;

        ///////////////////////////////////////

        ///////////////////////////////////////
        // CaCl2 in NH4NO3(PSC92: revised on 12/22/2003)
        b_mtem_host(0,jcacl2,jnh4no3) = -1.43626;
        b_mtem_host(1,jcacl2,jnh4no3) =  13.6598;
        b_mtem_host(2,jcacl2,jnh4no3) = -38.2068;
        b_mtem_host(3,jcacl2,jnh4no3) =  53.9057;
        b_mtem_host(4,jcacl2,jnh4no3) = -44.9018;
        b_mtem_host(5,jcacl2,jnh4no3) =  16.6120;

        // CaCl2 in NH4Cl(PSC92: revised on 11/27/2003)
        b_mtem_host(0,jcacl2,jnh4cl) =   -0.603965;
        b_mtem_host(1,jcacl2,jnh4cl) =   27.6027;
        b_mtem_host(2,jcacl2,jnh4cl) = -104.258;
        b_mtem_host(3,jcacl2,jnh4cl) =  163.553;
        b_mtem_host(4,jcacl2,jnh4cl) = -124.076;
        b_mtem_host(5,jcacl2,jnh4cl) =   37.4153;

        // CaCl2 in NaNO3(PSC92: revised on 12/22/2003)
        b_mtem_host(0,jcacl2,jnano3) =   0.44648;
        b_mtem_host(1,jcacl2,jnano3) =   8.8850;
        b_mtem_host(2,jcacl2,jnano3) = -45.5232;
        b_mtem_host(3,jcacl2,jnano3) =  89.3263;
        b_mtem_host(4,jcacl2,jnano3) = -83.8604;
        b_mtem_host(5,jcacl2,jnano3) =  30.4069;

        // CaCl2 in NaCl(PSC92: revised on 11/27/2003)
        b_mtem_host(0,jcacl2,jnacl) =   1.61927;
        b_mtem_host(1,jcacl2,jnacl) =   0.247547;
        b_mtem_host(2,jcacl2,jnacl) = -18.1252;
        b_mtem_host(3,jcacl2,jnacl) =  45.2479;
        b_mtem_host(4,jcacl2,jnacl) = -48.6072;
        b_mtem_host(5,jcacl2,jnacl) =  19.2784;

        // CaCl2 in Ca_host(NO3)2(PSC92: revised on 11/27/2003)
        b_mtem_host(0,jcacl2,jcano3) =  2.36667;
        b_mtem_host(1,jcacl2,jcano3) = -0.123309;
        b_mtem_host(2,jcacl2,jcano3) = -24.2723;
        b_mtem_host(3,jcacl2,jcano3) =  65.1486;
        b_mtem_host(4,jcacl2,jcano3) = -71.8504;
        b_mtem_host(5,jcacl2,jcano3) =  28.3696;

        // CaCl2 in CaCl2(PSC92: revised on 11/27/2003)
        b_mtem_host(0,jcacl2,jcacl2) =   3.64023;
        b_mtem_host(1,jcacl2,jcacl2) = -12.1926;
        b_mtem_host(2,jcacl2,jcacl2) =  20.2028;
        b_mtem_host(3,jcacl2,jcacl2) = -16.0056;
        b_mtem_host(4,jcacl2,jcacl2) =   1.52355;
        b_mtem_host(5,jcacl2,jcacl2) =   2.44709;

        // CaCl2 in HNO3
        b_mtem_host(0,jcacl2,jhno3) =    5.88794;
        b_mtem_host(1,jcacl2,jhno3) =  -29.7083;
        b_mtem_host(2,jcacl2,jhno3) =   78.6309;
        b_mtem_host(3,jcacl2,jhno3) = -118.037;
        b_mtem_host(4,jcacl2,jhno3) =   88.932;
        b_mtem_host(5,jcacl2,jhno3) =  -26.1407;

        // CaCl2 in HCl
        b_mtem_host(0,jcacl2,jhcl) =   2.40628;
        b_mtem_host(1,jcacl2,jhcl) =  -6.16566;
        b_mtem_host(2,jcacl2,jhcl) =  10.2851;
        b_mtem_host(3,jcacl2,jhcl) = -12.9035;
        b_mtem_host(4,jcacl2,jhcl) =   7.7441;
        b_mtem_host(5,jcacl2,jhcl) =  -1.74821;

        ///////////////////////////////////////

        ///////////////////////////////////////
        // HNO3 in(NH4)2SO4
        b_mtem_host(0,jhno3,jnh4so4) =   -3.57598;
        b_mtem_host(1,jhno3,jnh4so4) =   21.5469;
        b_mtem_host(2,jhno3,jnh4so4) =  -77.4111;
        b_mtem_host(3,jhno3,jnh4so4) =  144.136;
        b_mtem_host(4,jhno3,jnh4so4) = -132.849;
        b_mtem_host(5,jhno3,jnh4so4) =   47.9412;

        // HNO3 in NH4NO3
        b_mtem_host(0,jhno3,jnh4no3) =  -2.00209;
        b_mtem_host(1,jhno3,jnh4no3) =  -3.48399;
        b_mtem_host(2,jhno3,jnh4no3) =  34.9906;
        b_mtem_host(3,jhno3,jnh4no3) = -68.6653;
        b_mtem_host(4,jhno3,jnh4no3) =  54.0992;
        b_mtem_host(5,jhno3,jnh4no3) = -15.1343;

        // HNO3 in NH4Cl revised on 12/22/2003
        b_mtem_host(0,jhno3,jnh4cl) =  -0.63790;
        b_mtem_host(1,jhno3,jnh4cl) =  -1.67730;
        b_mtem_host(2,jhno3,jnh4cl) =  10.1727;
        b_mtem_host(3,jhno3,jnh4cl) = -14.9097;
        b_mtem_host(4,jhno3,jnh4cl) =   7.67410;
        b_mtem_host(5,jhno3,jnh4cl) =  -0.79586;

        // HNO3in NaCl
        b_mtem_host(0,jhno3,jnacl) =  1.3446;
        b_mtem_host(1,jhno3,jnacl) = -2.5578;
        b_mtem_host(2,jhno3,jnacl) =  1.3464;
        b_mtem_host(3,jhno3,jnacl) =  2.90537;
        b_mtem_host(4,jhno3,jnacl) = -6.53014;
        b_mtem_host(5,jhno3,jnacl) =  3.31339;

        // HNO3 in NaNO3
        b_mtem_host(0,jhno3,jnano3) =  -0.546636;
        b_mtem_host(1,jhno3,jnano3) =  10.3127;
        b_mtem_host(2,jhno3,jnano3) = -39.9603;
        b_mtem_host(3,jhno3,jnano3) =  71.4609;
        b_mtem_host(4,jhno3,jnano3) = -63.4958;
        b_mtem_host(5,jhno3,jnano3) =  22.0679;

        // HNO3 in Na2SO4
        b_mtem_host(0,jhno3,jna2so4) =   1.35059;
        b_mtem_host(1,jhno3,jna2so4) =   4.34557;
        b_mtem_host(2,jhno3,jna2so4) = -35.8425;
        b_mtem_host(3,jhno3,jna2so4) =  80.9868;
        b_mtem_host(4,jhno3,jna2so4) = -81.6544;
        b_mtem_host(5,jhno3,jna2so4) =  30.4841;

        // HNO3 in Ca_host(NO3)2
        b_mtem_host(0,jhno3,jcano3) =   0.869414;
        b_mtem_host(1,jhno3,jcano3) =   2.98486;
        b_mtem_host(2,jhno3,jcano3) = -22.255;
        b_mtem_host(3,jhno3,jcano3) =  50.1863;
        b_mtem_host(4,jhno3,jcano3) = -51.214;
        b_mtem_host(5,jhno3,jcano3) =  19.2235;

        // HNO3 in CaCl2(KM) revised on 12/22/2003
        b_mtem_host(0,jhno3,jcacl2) =   1.42800;
        b_mtem_host(1,jhno3,jcacl2) =  -1.78959;
        b_mtem_host(2,jhno3,jcacl2) =  -2.49075;
        b_mtem_host(3,jhno3,jcacl2) =  10.1877;
        b_mtem_host(4,jhno3,jcacl2) = -12.1948;
        b_mtem_host(5,jhno3,jcacl2) =   4.64475;

        // HNO3 in HNO3(added on 12/06/2004)
        b_mtem_host(0,jhno3,jhno3) =  0.22035;
        b_mtem_host(1,jhno3,jhno3) =  2.94973;
        b_mtem_host(2,jhno3,jhno3) = -12.1469;
        b_mtem_host(3,jhno3,jhno3) =  20.4905;
        b_mtem_host(4,jhno3,jhno3) = -17.3966;
        b_mtem_host(5,jhno3,jhno3) =   5.70779;

        // HNO3 in HCl(added on 12/06/2004)
        b_mtem_host(0,jhno3,jhcl) =  1.55503;
        b_mtem_host(1,jhno3,jhcl) = -3.61226;
        b_mtem_host(2,jhno3,jhcl) =  6.28265;
        b_mtem_host(3,jhno3,jhcl) = -8.69575;
        b_mtem_host(4,jhno3,jhcl) =  6.09372;
        b_mtem_host(5,jhno3,jhcl) = -1.80898;

        // HNO3 in H2SO4
        b_mtem_host(0,jhno3,jh2so4) =  1.10783;
        b_mtem_host(1,jhno3,jh2so4) = -1.3363;
        b_mtem_host(2,jhno3,jh2so4) = -1.83525;
        b_mtem_host(3,jhno3,jh2so4) =  7.47373;
        b_mtem_host(4,jhno3,jh2so4) = -9.72954;
        b_mtem_host(5,jhno3,jh2so4) =  4.12248;

        // HNO3 in NH4HSO4
        b_mtem_host(0,jhno3,jnh4hso4) =  -0.851026;
        b_mtem_host(1,jhno3,jnh4hso4) =  12.2515;
        b_mtem_host(2,jhno3,jnh4hso4) = -49.788;
        b_mtem_host(3,jhno3,jnh4hso4) =  91.6215;
        b_mtem_host(4,jhno3,jnh4hso4) = -81.4877;
        b_mtem_host(5,jhno3,jnh4hso4) =  28.0002;

        // HNO3 in(NH4)3H_host(SO4)2
        b_mtem_host(0,jhno3,jlvcite) =  -3.09464;
        b_mtem_host(1,jhno3,jlvcite) =  14.9303;
        b_mtem_host(2,jhno3,jlvcite) = -43.0454;
        b_mtem_host(3,jhno3,jlvcite) =  72.6695;
        b_mtem_host(4,jhno3,jlvcite) = -65.2140;
        b_mtem_host(5,jhno3,jlvcite) =  23.4814;

        // HNO3 in NaHSO4
        b_mtem_host(0,jhno3,jnahso4) =   1.22973;
        b_mtem_host(1,jhno3,jnahso4) =   2.82702;
        b_mtem_host(2,jhno3,jnahso4) = -17.5869;
        b_mtem_host(3,jhno3,jnahso4) =  28.9564;
        b_mtem_host(4,jhno3,jnahso4) = -23.5814;
        b_mtem_host(5,jhno3,jnahso4) =   7.91153;

        // HNO3 in Na3H_host(SO4)2
        b_mtem_host(0,jhno3,jna3hso4) =   1.64773;
        b_mtem_host(1,jhno3,jna3hso4) =   0.94188;
        b_mtem_host(2,jhno3,jna3hso4) = -19.1242;
        b_mtem_host(3,jhno3,jna3hso4) =  46.9887;
        b_mtem_host(4,jhno3,jna3hso4) = -50.9494;
        b_mtem_host(5,jhno3,jna3hso4) =  20.2169;

        ///////////////////////////////////////

        ///////////////////////////////////////
        // HCl in(NH4)2SO4
        b_mtem_host(0,jhcl,jnh4so4) =   -2.93783;
        b_mtem_host(1,jhcl,jnh4so4) =   20.5546;
        b_mtem_host(2,jhcl,jnh4so4) =  -75.8548;
        b_mtem_host(3,jhcl,jnh4so4) =  141.729;
        b_mtem_host(4,jhcl,jnh4so4) = -130.697;
        b_mtem_host(5,jhcl,jnh4so4) =   46.9905;

        // HCl in NH4NO3
        b_mtem_host(0,jhcl,jnh4no3) =  -1.69063;
        b_mtem_host(1,jhcl,jnh4no3) =  -1.85303;
        b_mtem_host(2,jhcl,jnh4no3) =  29.0927;
        b_mtem_host(3,jhcl,jnh4no3) = -58.7401;
        b_mtem_host(4,jhcl,jnh4no3) =  44.999;
        b_mtem_host(5,jhcl,jnh4no3) = -11.9988;

        // HCl in NH4Cl(revised on 11/15/2003)
        b_mtem_host(0,jhcl,jnh4cl) =  -0.2073;
        b_mtem_host(1,jhcl,jnh4cl) =  -0.4322;
        b_mtem_host(2,jhcl,jnh4cl) =   6.1271;
        b_mtem_host(3,jhcl,jnh4cl) = -12.3146;
        b_mtem_host(4,jhcl,jnh4cl) =   8.9919;
        b_mtem_host(5,jhcl,jnh4cl) =  -2.3388;

        // HCl in NaCl
        b_mtem_host(0,jhcl,jnacl) =  2.95913;
        b_mtem_host(1,jhcl,jnacl) = -7.92254;
        b_mtem_host(2,jhcl,jnacl) =  13.736;
        b_mtem_host(3,jhcl,jnacl) = -15.433;
        b_mtem_host(4,jhcl,jnacl) =  7.40386;
        b_mtem_host(5,jhcl,jnacl) = -0.918641;

        // HCl in NaNO3
        b_mtem_host(0,jhcl,jnano3) =   0.893272;
        b_mtem_host(1,jhcl,jnano3) =   6.53768;
        b_mtem_host(2,jhcl,jnano3) = -32.3458;
        b_mtem_host(3,jhcl,jnano3) =  61.2834;
        b_mtem_host(4,jhcl,jnano3) = -56.4446;
        b_mtem_host(5,jhcl,jnano3) =  19.9202;

        // HCl in Na2SO4
        b_mtem_host(0,jhcl,jna2so4) =   3.14484;
        b_mtem_host(1,jhcl,jna2so4) =   0.077019;
        b_mtem_host(2,jhcl,jna2so4) = -31.4199;
        b_mtem_host(3,jhcl,jna2so4) =  80.5865;
        b_mtem_host(4,jhcl,jna2so4) = -85.392;
        b_mtem_host(5,jhcl,jna2so4) =  32.6644;

        // HCl in Ca_host(NO3)2
        b_mtem_host(0,jhcl,jcano3) =   2.60432;
        b_mtem_host(1,jhcl,jcano3) =  -0.55909;
        b_mtem_host(2,jhcl,jcano3) = -19.6671;
        b_mtem_host(3,jhcl,jcano3) =  53.3446;
        b_mtem_host(4,jhcl,jcano3) = -58.9076;
        b_mtem_host(5,jhcl,jcano3) =  22.9927;

        // HCl in CaCl2(KM) revised on 3/13/2003 and again on 11/27/2003
        b_mtem_host(0,jhcl,jcacl2) =   2.98036;
        b_mtem_host(1,jhcl,jcacl2) =  -8.55365;
        b_mtem_host(2,jhcl,jcacl2) =  15.2108;
        b_mtem_host(3,jhcl,jcacl2) = -15.9359;
        b_mtem_host(4,jhcl,jcacl2) =   7.41772;
        b_mtem_host(5,jhcl,jcacl2) =  -1.32143;

        // HCl in HNO3(added on 12/06/2004)
        b_mtem_host(0,jhcl,jhno3) =   3.8533;
        b_mtem_host(1,jhcl,jhno3) = -16.9427;
        b_mtem_host(2,jhcl,jhno3) =  45.0056;
        b_mtem_host(3,jhcl,jhno3) = -69.6145;
        b_mtem_host(4,jhcl,jhno3) =  54.1491;
        b_mtem_host(5,jhcl,jhno3) = -16.6513;

        // HCl in HCl(added on 12/06/2004)
        b_mtem_host(0,jhcl,jhcl) =   2.56665;
        b_mtem_host(1,jhcl,jhcl) =  -7.13585;
        b_mtem_host(2,jhcl,jhcl) =  14.8103;
        b_mtem_host(3,jhcl,jhcl) = -21.8881;
        b_mtem_host(4,jhcl,jhcl) =  16.6808;
        b_mtem_host(5,jhcl,jhcl) =  -5.22091;

        // HCl in H2SO4
        b_mtem_host(0,jhcl,jh2so4) =   2.50179;
        b_mtem_host(1,jhcl,jh2so4) =  -6.69364;
        b_mtem_host(2,jhcl,jh2so4) =  11.6551;
        b_mtem_host(3,jhcl,jh2so4) = -13.6897;
        b_mtem_host(4,jhcl,jh2so4) =   7.36796;
        b_mtem_host(5,jhcl,jh2so4) =  -1.33245;

        // HCl in NH4HSO4
        b_mtem_host(0,jhcl,jnh4hso4) =   0.149955;
        b_mtem_host(1,jhcl,jnh4hso4) =  11.8213;
        b_mtem_host(2,jhcl,jnh4hso4) = -53.9164;
        b_mtem_host(3,jhcl,jnh4hso4) = 101.574;
        b_mtem_host(4,jhcl,jnh4hso4) = -91.4123;
        b_mtem_host(5,jhcl,jnh4hso4) =  31.5487;

        // HCl in(NH4)3H(SO4)2
        b_mtem_host(0,jhcl,jlvcite) =  -2.36927;
        b_mtem_host(1,jhcl,jlvcite) =  14.8359;
        b_mtem_host(2,jhcl,jlvcite) = -44.3443;
        b_mtem_host(3,jhcl,jlvcite) =  73.6229;
        b_mtem_host(4,jhcl,jlvcite) = -65.3366;
        b_mtem_host(5,jhcl,jlvcite) =  23.3250;

        // HCl in NaHSO4
        b_mtem_host(0,jhcl,jnahso4) =   2.72993;
        b_mtem_host(1,jhcl,jnahso4) =  -0.23406;
        b_mtem_host(2,jhcl,jnahso4) = -10.4103;
        b_mtem_host(3,jhcl,jnahso4) =  13.1586;
        b_mtem_host(4,jhcl,jnahso4) =  -7.79925;
        b_mtem_host(5,jhcl,jnahso4) =   2.30843;

        // HCl in Na3H(SO4)2
        b_mtem_host(0,jhcl,jna3hso4) =   3.51258;
        b_mtem_host(1,jhcl,jna3hso4) =  -3.95107;
        b_mtem_host(2,jhcl,jna3hso4) = -11.0175;
        b_mtem_host(3,jhcl,jna3hso4) =  38.8617;
        b_mtem_host(4,jhcl,jna3hso4) = -48.1575;
        b_mtem_host(5,jhcl,jna3hso4) =  20.4717;

        ///////////////////////////////////////

        ///////////////////////////////////////
        // H2SO4 in H2SO4
        b_mtem_host(0,jh2so4,jh2so4) =   0.76734;
        b_mtem_host(1,jh2so4,jh2so4) =  -1.12263;
        b_mtem_host(2,jh2so4,jh2so4) =  -9.08728;
        b_mtem_host(3,jh2so4,jh2so4) =  30.3836;
        b_mtem_host(4,jh2so4,jh2so4) = -38.4133;
        b_mtem_host(5,jh2so4,jh2so4) =  17.0106;

        // H2SO4 in NH4HSO4
        b_mtem_host(0,jh2so4,jnh4hso4) =   -2.03879;
        b_mtem_host(1,jh2so4,jnh4hso4) =   15.7033;
        b_mtem_host(2,jh2so4,jnh4hso4) =  -58.7363;
        b_mtem_host(3,jh2so4,jnh4hso4) =  109.242;
        b_mtem_host(4,jh2so4,jnh4hso4) = -102.237;
        b_mtem_host(5,jh2so4,jnh4hso4) =   37.5350;

        // H2SO4 in(NH4)3H(SO4)2
        b_mtem_host(0,jh2so4,jlvcite) =   -3.10228;
        b_mtem_host(1,jh2so4,jlvcite) =   16.6920;
        b_mtem_host(2,jh2so4,jlvcite) =  -59.1522;
        b_mtem_host(3,jh2so4,jlvcite) =  113.487;
        b_mtem_host(4,jh2so4,jlvcite) = -110.890;
        b_mtem_host(5,jh2so4,jlvcite) =   42.4578;

        // H2SO4 in(NH4)2SO4
        b_mtem_host(0,jh2so4,jnh4so4) =   -3.43885;
        b_mtem_host(1,jh2so4,jnh4so4) =   21.0372;
        b_mtem_host(2,jh2so4,jnh4so4) =  -84.7026;
        b_mtem_host(3,jh2so4,jnh4so4) =  165.324;
        b_mtem_host(4,jh2so4,jnh4so4) = -156.101;
        b_mtem_host(5,jh2so4,jnh4so4) =   57.3101;

        // H2SO4 in NaHSO4
        b_mtem_host(0,jh2so4,jnahso4) =   0.33164;
        b_mtem_host(1,jh2so4,jnahso4) =   6.55864;
        b_mtem_host(2,jh2so4,jnahso4) = -33.5876;
        b_mtem_host(3,jh2so4,jnahso4) =  65.1798;
        b_mtem_host(4,jh2so4,jnahso4) = -63.2046;
        b_mtem_host(5,jh2so4,jnahso4) =  24.1783;

        // H2SO4 in Na3H(SO4)2
        b_mtem_host(0,jh2so4,jna3hso4) =   3.06830;
        b_mtem_host(1,jh2so4,jna3hso4) =  -3.18408;
        b_mtem_host(2,jh2so4,jna3hso4) = -19.6332;
        b_mtem_host(3,jh2so4,jna3hso4) =  61.3657;
        b_mtem_host(4,jh2so4,jna3hso4) = -73.4438;
        b_mtem_host(5,jh2so4,jna3hso4) =  31.2334;

        // H2SO4 in Na2SO4
        b_mtem_host(0,jh2so4,jna2so4) =    2.58649;
        b_mtem_host(1,jh2so4,jna2so4) =    0.87921;
        b_mtem_host(2,jh2so4,jna2so4) =  -39.3023;
        b_mtem_host(3,jh2so4,jna2so4) =  101.603;
        b_mtem_host(4,jh2so4,jna2so4) = -109.469;
        b_mtem_host(5,jh2so4,jna2so4) =   43.0188;

        // H2SO4 in HNO3
        b_mtem_host(0,jh2so4,jhno3) =   1.54587;
        b_mtem_host(1,jh2so4,jhno3) =  -7.50976;
        b_mtem_host(2,jh2so4,jhno3) =  12.8237;
        b_mtem_host(3,jh2so4,jhno3) = -10.1452;
        b_mtem_host(4,jh2so4,jhno3) =  -0.541956;
        b_mtem_host(5,jh2so4,jhno3) =   3.34536;

        // H2SO4 in HCl
        b_mtem_host(0,jh2so4,jhcl) =   0.829757;
        b_mtem_host(1,jh2so4,jhcl) =  -4.11316;
        b_mtem_host(2,jh2so4,jhcl) =   3.67111;
        b_mtem_host(3,jh2so4,jhcl) =   3.6833;
        b_mtem_host(4,jh2so4,jhcl) = -11.2711;
        b_mtem_host(5,jh2so4,jhcl) =   6.71421;

        ///////////////////////////////////////

        ///////////////////////////////////////
        // H.HSO4 in H2SO4
        b_mtem_host(0,jhhso4,jh2so4) =   2.63953;
        b_mtem_host(1,jhhso4,jh2so4) =  -6.01532;
        b_mtem_host(2,jhhso4,jh2so4) =  10.0204;
        b_mtem_host(3,jhhso4,jh2so4) = -12.4840;
        b_mtem_host(4,jhhso4,jh2so4) =   7.78853;
        b_mtem_host(5,jhhso4,jh2so4) =  -2.12638;

        // H.HSO4 in NH4HSO4
        b_mtem_host(0,jhhso4,jnh4hso4) =  -0.77412;
        b_mtem_host(1,jhhso4,jnh4hso4) =  14.1656;
        b_mtem_host(2,jhhso4,jnh4hso4) = -53.4087;
        b_mtem_host(3,jhhso4,jnh4hso4) =  93.2013;
        b_mtem_host(4,jhhso4,jnh4hso4) = -80.5723;
        b_mtem_host(5,jhhso4,jnh4hso4) =  27.1577;

        // H.HSO4 in(NH4)3H(SO4)2
        b_mtem_host(0,jhhso4,jlvcite) =  -2.98882;
        b_mtem_host(1,jhhso4,jlvcite) =  14.4436;
        b_mtem_host(2,jhhso4,jlvcite) = -40.1774;
        b_mtem_host(3,jhhso4,jlvcite) =  67.5937;
        b_mtem_host(4,jhhso4,jlvcite) = -61.5040;
        b_mtem_host(5,jhhso4,jlvcite) =  22.3695;

        // H.HSO4 in(NH4)2SO4
        b_mtem_host(0,jhhso4,jnh4so4) =  -1.15502;
        b_mtem_host(1,jhhso4,jnh4so4) =   8.12309;
        b_mtem_host(2,jhhso4,jnh4so4) = -38.4726;
        b_mtem_host(3,jhhso4,jnh4so4) =  80.8861;
        b_mtem_host(4,jhhso4,jnh4so4) = -80.1644;
        b_mtem_host(5,jhhso4,jnh4so4) =  30.4717;

        // H.HSO4 in NaHSO4
        b_mtem_host(0,jhhso4,jnahso4) =   1.99641;
        b_mtem_host(1,jhhso4,jnahso4) =  -2.96061;
        b_mtem_host(2,jhhso4,jnahso4) =   5.54778;
        b_mtem_host(3,jhhso4,jnahso4) = -14.5488;
        b_mtem_host(4,jhhso4,jnahso4) =  14.8492;
        b_mtem_host(5,jhhso4,jnahso4) =  -5.1389;

        // H.HSO4 in Na3H(SO4)2
        b_mtem_host(0,jhhso4,jna3hso4) =   2.23816;
        b_mtem_host(1,jhhso4,jna3hso4) =  -3.20847;
        b_mtem_host(2,jhhso4,jna3hso4) =  -4.82853;
        b_mtem_host(3,jhhso4,jna3hso4) =  20.9192;
        b_mtem_host(4,jhhso4,jna3hso4) = -27.2819;
        b_mtem_host(5,jhhso4,jna3hso4) =  11.8655;

        // H.HSO4 in Na2SO4
        b_mtem_host(0,jhhso4,jna2so4) =   2.56907;
        b_mtem_host(1,jhhso4,jna2so4) =   1.13444;
        b_mtem_host(2,jhhso4,jna2so4) = -34.6853;
        b_mtem_host(3,jhhso4,jna2so4) =  87.9775;
        b_mtem_host(4,jhhso4,jna2so4) = -93.2330;
        b_mtem_host(5,jhhso4,jna2so4) =  35.9260;

        // H.HSO4 in HNO3
        b_mtem_host(0,jhhso4,jhno3) =   2.00024;
        b_mtem_host(1,jhhso4,jhno3) =  -4.80868;
        b_mtem_host(2,jhhso4,jhno3) =   8.29222;
        b_mtem_host(3,jhhso4,jhno3) = -11.0849;
        b_mtem_host(4,jhhso4,jhno3) =   7.51262;
        b_mtem_host(5,jhhso4,jhno3) =  -2.07654;

        // H.HSO4 in HCl
        b_mtem_host(0,jhhso4,jhcl) =   2.8009;
        b_mtem_host(1,jhhso4,jhcl) =  -6.98416;
        b_mtem_host(2,jhhso4,jhcl) =  14.3146;
        b_mtem_host(3,jhhso4,jhcl) = -22.0068;
        b_mtem_host(4,jhhso4,jhcl) =  17.5557;
        b_mtem_host(5,jhhso4,jhcl) =  -5.84917;

        ///////////////////////////////////////

        ///////////////////////////////////////
        // NH4HSO4 in H2SO4
        b_mtem_host(0,jnh4hso4,jh2so4) =   0.169160;
        b_mtem_host(1,jnh4hso4,jh2so4) =   2.15094;
        b_mtem_host(2,jnh4hso4,jh2so4) =  -9.62904;
        b_mtem_host(3,jnh4hso4,jh2so4) =  18.2631;
        b_mtem_host(4,jnh4hso4,jh2so4) = -17.3333;
        b_mtem_host(5,jnh4hso4,jh2so4) =   6.19835;

        // NH4HSO4 in NH4HSO4
        b_mtem_host(0,jnh4hso4,jnh4hso4) =  -2.34457;
        b_mtem_host(1,jnh4hso4,jnh4hso4) =  12.8035;
        b_mtem_host(2,jnh4hso4,jnh4hso4) = -35.2513;
        b_mtem_host(3,jnh4hso4,jnh4hso4) =  53.6153;
        b_mtem_host(4,jnh4hso4,jnh4hso4) = -42.7655;
        b_mtem_host(5,jnh4hso4,jnh4hso4) =  13.7129;

        // NH4HSO4 in(NH4)3H(SO4)2
        b_mtem_host(0,jnh4hso4,jlvcite) =  -2.56109;
        b_mtem_host(1,jnh4hso4,jlvcite) =  11.1414;
        b_mtem_host(2,jnh4hso4,jlvcite) = -30.2361;
        b_mtem_host(3,jnh4hso4,jlvcite) =  50.0320;
        b_mtem_host(4,jnh4hso4,jlvcite) = -44.1586;
        b_mtem_host(5,jnh4hso4,jlvcite) =  15.5393;

        // NH4HSO4 in(NH4)2SO4
        b_mtem_host(0,jnh4hso4,jnh4so4) =  -0.97315;
        b_mtem_host(1,jnh4hso4,jnh4so4) =   7.06295;
        b_mtem_host(2,jnh4hso4,jnh4so4) = -29.3032;
        b_mtem_host(3,jnh4hso4,jnh4so4) =  57.6101;
        b_mtem_host(4,jnh4hso4,jnh4so4) = -54.9020;
        b_mtem_host(5,jnh4hso4,jnh4so4) =  20.2222;

        // NH4HSO4 in NaHSO4
        b_mtem_host(0,jnh4hso4,jnahso4) =  -0.44450;
        b_mtem_host(1,jnh4hso4,jnahso4) =   3.33451;
        b_mtem_host(2,jnh4hso4,jnahso4) = -15.2791;
        b_mtem_host(3,jnh4hso4,jnahso4) =  30.1413;
        b_mtem_host(4,jnh4hso4,jnahso4) = -26.7710;
        b_mtem_host(5,jnh4hso4,jnahso4) =   8.78462;

        // NH4HSO4 in Na3H(SO4)2
        b_mtem_host(0,jnh4hso4,jna3hso4) =  -0.99780;
        b_mtem_host(1,jnh4hso4,jna3hso4) =   4.69200;
        b_mtem_host(2,jnh4hso4,jna3hso4) = -16.1219;
        b_mtem_host(3,jnh4hso4,jna3hso4) =  29.3100;
        b_mtem_host(4,jnh4hso4,jna3hso4) = -26.3383;
        b_mtem_host(5,jnh4hso4,jna3hso4) =   9.20695;

        // NH4HSO4 in Na2SO4
        b_mtem_host(0,jnh4hso4,jna2so4) =  -0.52694;
        b_mtem_host(1,jnh4hso4,jna2so4) =   7.02684;
        b_mtem_host(2,jnh4hso4,jna2so4) = -33.7508;
        b_mtem_host(3,jnh4hso4,jna2so4) =  70.0565;
        b_mtem_host(4,jnh4hso4,jna2so4) = -68.3226;
        b_mtem_host(5,jnh4hso4,jna2so4) =  25.2692;

        // NH4HSO4 in HNO3
        b_mtem_host(0,jnh4hso4,jhno3) =  0.572926;
        b_mtem_host(1,jnh4hso4,jhno3) = -2.04791;
        b_mtem_host(2,jnh4hso4,jhno3) =  2.1134;
        b_mtem_host(3,jnh4hso4,jhno3) =  0.246654;
        b_mtem_host(4,jnh4hso4,jhno3) = -3.06019;
        b_mtem_host(5,jnh4hso4,jhno3) =  1.98126;

        // NH4HSO4 in HCl
        b_mtem_host(0,jnh4hso4,jhcl) =  0.56514;
        b_mtem_host(1,jnh4hso4,jhcl) =  0.22287;
        b_mtem_host(2,jnh4hso4,jhcl) = -2.76973;
        b_mtem_host(3,jnh4hso4,jhcl) =  4.54444;
        b_mtem_host(4,jnh4hso4,jhcl) = -3.86549;
        b_mtem_host(5,jnh4hso4,jhcl) =  1.13441;

        ///////////////////////////////////////

        ///////////////////////////////////////
        // (NH4)3H(SO4)2 in H2SO4
        b_mtem_host(0,jlvcite,jh2so4) =  -1.44811;
        b_mtem_host(1,jlvcite,jh2so4) =   6.71815;
        b_mtem_host(2,jlvcite,jh2so4) = -25.0141;
        b_mtem_host(3,jlvcite,jh2so4) =  50.1109;
        b_mtem_host(4,jlvcite,jh2so4) = -50.0561;
        b_mtem_host(5,jlvcite,jh2so4) =  19.3370;

        // (NH4)3H(SO4)2 in NH4HSO4
        b_mtem_host(0,jlvcite,jnh4hso4) =  -3.41707;
        b_mtem_host(1,jlvcite,jnh4hso4) =  13.4496;
        b_mtem_host(2,jlvcite,jnh4hso4) = -34.8018;
        b_mtem_host(3,jlvcite,jnh4hso4) =  55.2987;
        b_mtem_host(4,jlvcite,jnh4hso4) = -48.1839;
        b_mtem_host(5,jlvcite,jnh4hso4) =  17.2444;

        // (NH4)3H(SO4)2 in(NH4)3H(SO4)2
        b_mtem_host(0,jlvcite,jlvcite) =  -2.54479;
        b_mtem_host(1,jlvcite,jlvcite) =  11.8501;
        b_mtem_host(2,jlvcite,jlvcite) = -39.7286;
        b_mtem_host(3,jlvcite,jlvcite) =  74.2479;
        b_mtem_host(4,jlvcite,jlvcite) = -70.4934;
        b_mtem_host(5,jlvcite,jlvcite) =  26.2836;

        // (NH4)3H(SO4)2 in(NH4)2SO4
        b_mtem_host(0,jlvcite,jnh4so4) =  -2.30561;
        b_mtem_host(1,jlvcite,jnh4so4) =  14.5806;
        b_mtem_host(2,jlvcite,jnh4so4) = -55.1238;
        b_mtem_host(3,jlvcite,jnh4so4) = 103.451;
        b_mtem_host(4,jlvcite,jnh4so4) = -95.2571;
        b_mtem_host(5,jlvcite,jnh4so4) =  34.2218;

        // (NH4)3H(SO4)2 in NaHSO4
        b_mtem_host(0,jlvcite,jnahso4) =   -2.20809;
        b_mtem_host(1,jlvcite,jnahso4) =   13.6391;
        b_mtem_host(2,jlvcite,jnahso4) =  -57.8246;
        b_mtem_host(3,jlvcite,jnahso4) =  117.907;
        b_mtem_host(4,jlvcite,jnahso4) = -112.154;
        b_mtem_host(5,jlvcite,jnahso4) =   40.3058;

        // (NH4)3H(SO4)2 in Na3H(SO4)2
        b_mtem_host(0,jlvcite,jna3hso4) =  -1.15099;
        b_mtem_host(1,jlvcite,jna3hso4) =   6.32269;
        b_mtem_host(2,jlvcite,jna3hso4) = -27.3860;
        b_mtem_host(3,jlvcite,jna3hso4) =  55.4592;
        b_mtem_host(4,jlvcite,jna3hso4) = -54.0100;
        b_mtem_host(5,jlvcite,jna3hso4) =  20.3469;

        // (NH4)3H(SO4)2 in Na2SO4
        b_mtem_host(0,jlvcite,jna2so4) =  -1.15678;
        b_mtem_host(1,jlvcite,jna2so4) =   8.28718;
        b_mtem_host(2,jlvcite,jna2so4) = -37.3231;
        b_mtem_host(3,jlvcite,jna2so4) =  76.6124;
        b_mtem_host(4,jlvcite,jna2so4) = -74.9307;
        b_mtem_host(5,jlvcite,jna2so4) =  28.0559;

        // (NH4)3H(SO4)2 in HNO3
        b_mtem_host(0,jlvcite,jhno3) =  0.01502;
        b_mtem_host(1,jlvcite,jhno3) = -3.1197;
        b_mtem_host(2,jlvcite,jhno3) =  3.61104;
        b_mtem_host(3,jlvcite,jhno3) =  3.05196;
        b_mtem_host(4,jlvcite,jhno3) = -9.98957;
        b_mtem_host(5,jlvcite,jhno3) =  6.04155;

        // (NH4)3H(SO4)2 in HCl
        b_mtem_host(0,jlvcite,jhcl) =  -1.06477;
        b_mtem_host(1,jlvcite,jhcl) =   3.38801;
        b_mtem_host(2,jlvcite,jhcl) = -12.5784;
        b_mtem_host(3,jlvcite,jhcl) =  25.2823;
        b_mtem_host(4,jlvcite,jhcl) = -25.4611;
        b_mtem_host(5,jlvcite,jhcl) =  10.0754;

        ///////////////////////////////////////

        ///////////////////////////////////////
        // NaHSO4 in H2SO4
        b_mtem_host(0,jnahso4,jh2so4) =   0.68259;
        b_mtem_host(1,jnahso4,jh2so4) =   0.71468;
        b_mtem_host(2,jnahso4,jh2so4) =  -5.59003;
        b_mtem_host(3,jnahso4,jh2so4) =  11.0089;
        b_mtem_host(4,jnahso4,jh2so4) = -10.7983;
        b_mtem_host(5,jnahso4,jh2so4) =   3.82335;

        // NaHSO4 in NH4HSO4
        b_mtem_host(0,jnahso4,jnh4hso4) =  -0.03956;
        b_mtem_host(1,jnahso4,jnh4hso4) =   4.52828;
        b_mtem_host(2,jnahso4,jnh4hso4) = -25.2557;
        b_mtem_host(3,jnahso4,jnh4hso4) =  54.4225;
        b_mtem_host(4,jnahso4,jnh4hso4) = -52.5105;
        b_mtem_host(5,jnahso4,jnh4hso4) =  18.6562;

        // NaHSO4 in(NH4)3H(SO4)2
        b_mtem_host(0,jnahso4,jlvcite) =  -1.53503;
        b_mtem_host(1,jnahso4,jlvcite) =   8.27608;
        b_mtem_host(2,jnahso4,jlvcite) = -28.9539;
        b_mtem_host(3,jnahso4,jlvcite) =  55.2876;
        b_mtem_host(4,jnahso4,jlvcite) = -51.9563;
        b_mtem_host(5,jnahso4,jlvcite) =  18.6576;

        // NaHSO4 in(NH4)2SO4
        b_mtem_host(0,jnahso4,jnh4so4) =  -0.38793;
        b_mtem_host(1,jnahso4,jnh4so4) =   7.14680;
        b_mtem_host(2,jnahso4,jnh4so4) = -38.7201;
        b_mtem_host(3,jnahso4,jnh4so4) =  84.3965;
        b_mtem_host(4,jnahso4,jnh4so4) = -84.7453;
        b_mtem_host(5,jnahso4,jnh4so4) =  32.1283;

        // NaHSO4 in NaHSO4
        b_mtem_host(0,jnahso4,jnahso4) =  -0.41982;
        b_mtem_host(1,jnahso4,jnahso4) =   4.26491;
        b_mtem_host(2,jnahso4,jnahso4) = -20.2351;
        b_mtem_host(3,jnahso4,jnahso4) =  42.6764;
        b_mtem_host(4,jnahso4,jnahso4) = -40.7503;
        b_mtem_host(5,jnahso4,jnahso4) =  14.2868;

        // NaHSO4 in Na3H(SO4)2
        b_mtem_host(0,jnahso4,jna3hso4) =  -0.32912;
        b_mtem_host(1,jnahso4,jna3hso4) =   1.80808;
        b_mtem_host(2,jnahso4,jna3hso4) =  -8.01286;
        b_mtem_host(3,jnahso4,jna3hso4) =  15.5791;
        b_mtem_host(4,jnahso4,jna3hso4) = -14.5494;
        b_mtem_host(5,jnahso4,jna3hso4) =   5.27052;

        // NaHSO4 in Na2SO4
        b_mtem_host(0,jnahso4,jna2so4) =   0.10271;
        b_mtem_host(1,jnahso4,jna2so4) =   5.09559;
        b_mtem_host(2,jnahso4,jna2so4) = -30.3295;
        b_mtem_host(3,jnahso4,jna2so4) =  66.2975;
        b_mtem_host(4,jnahso4,jna2so4) = -66.3458;
        b_mtem_host(5,jnahso4,jna2so4) =  24.9443;

        // NaHSO4 in HNO3
        b_mtem_host(0,jnahso4,jhno3) =  0.608309;
        b_mtem_host(1,jnahso4,jhno3) = -0.541905;
        b_mtem_host(2,jnahso4,jhno3) = -2.52084;
        b_mtem_host(3,jnahso4,jhno3) =  6.63297;
        b_mtem_host(4,jnahso4,jhno3) = -7.24599;
        b_mtem_host(5,jnahso4,jhno3) =  2.88811;

        // NaHSO4 in HCl
        b_mtem_host(0,jnahso4,jhcl) =   1.98399;
        b_mtem_host(1,jnahso4,jhcl) =  -4.51562;
        b_mtem_host(2,jnahso4,jhcl) =   8.36059;
        b_mtem_host(3,jnahso4,jhcl) = -12.4948;
        b_mtem_host(4,jnahso4,jhcl) =   9.67514;
        b_mtem_host(5,jnahso4,jhcl) =  -3.18004;

        ///////////////////////////////////////

        ///////////////////////////////////////
        // Na3H(SO4)2 in H2SO4
        b_mtem_host(0,jna3hso4,jh2so4) =  -0.83214;
        b_mtem_host(1,jna3hso4,jh2so4) =   4.99572;
        b_mtem_host(2,jna3hso4,jh2so4) = -20.1697;
        b_mtem_host(3,jna3hso4,jh2so4) =  41.4066;
        b_mtem_host(4,jna3hso4,jh2so4) = -42.2119;
        b_mtem_host(5,jna3hso4,jh2so4) =  16.4855;

        // Na3H(SO4)2 in NH4HSO4
        b_mtem_host(0,jna3hso4,jnh4hso4) =  -0.65139;
        b_mtem_host(1,jna3hso4,jnh4hso4) =   3.52300;
        b_mtem_host(2,jna3hso4,jnh4hso4) = -22.8220;
        b_mtem_host(3,jna3hso4,jnh4hso4) =  56.2956;
        b_mtem_host(4,jna3hso4,jnh4hso4) = -59.9028;
        b_mtem_host(5,jna3hso4,jnh4hso4) =  23.1844;

        // Na3H(SO4)2 in(NH4)3H(SO4)2
        b_mtem_host(0,jna3hso4,jlvcite) =  -1.31331;
        b_mtem_host(1,jna3hso4,jlvcite) =   8.40835;
        b_mtem_host(2,jna3hso4,jlvcite) = -38.1757;
        b_mtem_host(3,jna3hso4,jlvcite) =  80.5312;
        b_mtem_host(4,jna3hso4,jlvcite) = -79.8346;
        b_mtem_host(5,jna3hso4,jlvcite) =  30.0219;

        // Na3H(SO4)2 in(NH4)2SO4
        b_mtem_host(0,jna3hso4,jnh4so4) =  -1.03054;
        b_mtem_host(1,jna3hso4,jnh4so4) =   8.08155;
        b_mtem_host(2,jna3hso4,jnh4so4) = -38.1046;
        b_mtem_host(3,jna3hso4,jnh4so4) =  78.7168;
        b_mtem_host(4,jna3hso4,jnh4so4) = -77.2263;
        b_mtem_host(5,jna3hso4,jnh4so4) =  29.1521;

        // Na3H(SO4)2 in NaHSO4
        b_mtem_host(0,jna3hso4,jnahso4) =   -1.90695;
        b_mtem_host(1,jna3hso4,jnahso4) =   11.6241;
        b_mtem_host(2,jna3hso4,jnahso4) =  -50.3175;
        b_mtem_host(3,jna3hso4,jnahso4) =  105.884;
        b_mtem_host(4,jna3hso4,jnahso4) = -103.258;
        b_mtem_host(5,jna3hso4,jnahso4) =   37.6588;

        // Na3H(SO4)2 in Na3H(SO4)2
        b_mtem_host(0,jna3hso4,jna3hso4) =  -0.34780;
        b_mtem_host(1,jna3hso4,jna3hso4) =   2.85363;
        b_mtem_host(2,jna3hso4,jna3hso4) = -17.6224;
        b_mtem_host(3,jna3hso4,jna3hso4) =  38.9220;
        b_mtem_host(4,jna3hso4,jna3hso4) = -39.8106;
        b_mtem_host(5,jna3hso4,jna3hso4) =  15.6055;

        // Na3H(SO4)2 in Na2SO4
        b_mtem_host(0,jna3hso4,jna2so4) =   -0.75230;
        b_mtem_host(1,jna3hso4,jna2so4) =   10.0140;
        b_mtem_host(2,jna3hso4,jna2so4) =  -50.5677;
        b_mtem_host(3,jna3hso4,jna2so4) =  106.941;
        b_mtem_host(4,jna3hso4,jna2so4) = -105.534;
        b_mtem_host(5,jna3hso4,jna2so4) =   39.5196;

        // Na3H(SO4)2 in HNO3
        b_mtem_host(0,jna3hso4,jhno3) =   0.057456;
        b_mtem_host(1,jna3hso4,jhno3) =  -1.31264;
        b_mtem_host(2,jna3hso4,jhno3) =  -1.94662;
        b_mtem_host(3,jna3hso4,jhno3) =  10.7024;
        b_mtem_host(4,jna3hso4,jhno3) = -14.9946;
        b_mtem_host(5,jna3hso4,jhno3) =   7.12161;

        // Na3H(SO4)2 in HCl
        b_mtem_host(0,jna3hso4,jhcl) =  0.637894;
        b_mtem_host(1,jna3hso4,jhcl) = -2.29719;
        b_mtem_host(2,jna3hso4,jhcl) =  0.765361;
        b_mtem_host(3,jna3hso4,jhcl) =  4.8748;
        b_mtem_host(4,jna3hso4,jhcl) = -9.25978;
        b_mtem_host(5,jna3hso4,jhcl) =  4.91773;

        /////////////////////////////////////

        // MDRH coefficients
        // SULFATE-POOR SYSTEMS
        // AC
        d_mdrh_host(0,0) = -58.00268351;
        d_mdrh_host(0,1) = 2.031077573;
        d_mdrh_host(0,2) = -0.008281218;
        d_mdrh_host(0,3) = 1.00447E-05;

        // AN
        d_mdrh_host(1,0) = 1039.137773;
        d_mdrh_host(1,1) = -11.47847095;
        d_mdrh_host(1,2) = 0.047702786;
        d_mdrh_host(1,3) = -6.77675E-05;

        // AS
        d_mdrh_host(2,0) = 115.8366357;
        d_mdrh_host(2,1) = 0.491881663;
        d_mdrh_host(2,2) = -0.00422807;
        d_mdrh_host(2,3) = 7.29274E-06;

        // SC
        d_mdrh_host(3,0) = 253.2424151;
        d_mdrh_host(3,1) = -1.429957864;
        d_mdrh_host(3,2) = 0.003727554;
        d_mdrh_host(3,3) = -3.13037E-06;

        // SN
        d_mdrh_host(4,0) = -372.4306506;
        d_mdrh_host(4,1) = 5.3955633;
        d_mdrh_host(4,2) = -0.019804438;
        d_mdrh_host(4,3) = 2.25662E-05;

        // SS
        d_mdrh_host(5,0) = 286.1271416;
        d_mdrh_host(5,1) = -1.670787758;
        d_mdrh_host(5,2) = 0.004431373;
        d_mdrh_host(5,3) = -3.57757E-06;

        // CC
        d_mdrh_host(6,0) = -1124.07059;
        d_mdrh_host(6,1) = 14.26364209;
        d_mdrh_host(6,2) = -0.054816822;
        d_mdrh_host(6,3) = 6.70107E-05;

        // CN
        d_mdrh_host(7,0) = 1855.413934;
        d_mdrh_host(7,1) = -20.29219473;
        d_mdrh_host(7,2) = 0.07807482;
        d_mdrh_host(7,3) = -1.017887858e-4;

        // AN + AC
        d_mdrh_host(8,0) = 1761.176886;
        d_mdrh_host(8,1) = -19.29811062;
        d_mdrh_host(8,2) = 0.075676987;
        d_mdrh_host(8,3) = -1.0116959e-4;

        // AS + AC
        d_mdrh_host(9,0) = 122.1074303;
        d_mdrh_host(9,1) = 0.429692122;
        d_mdrh_host(9,2) = -0.003928277;
        d_mdrh_host(9,3) = 6.43275E-06;

        // AS + AN
        d_mdrh_host(10,0) = 2424.634678;
        d_mdrh_host(10,1) = -26.54031307;
        d_mdrh_host(10,2) = 0.101625387;
        d_mdrh_host(10,3) = -1.31544547798e-4;

        // AS + AN + AC
        d_mdrh_host(11,0) = 2912.082599;
        d_mdrh_host(11,1) = -31.8894185;
        d_mdrh_host(11,2) = 0.121185849;
        d_mdrh_host(11,3) = -1.556534623e-4;

        // SC + AC
        d_mdrh_host(12,0) = 172.2596493;
        d_mdrh_host(12,1) = -0.511006195;
        d_mdrh_host(12,2) = 4.27244597e-4;
        d_mdrh_host(12,3) = 4.12797E-07;

        // SN + AC
        d_mdrh_host(13,0) = 1596.184935;
        d_mdrh_host(13,1) = -16.37945565;
        d_mdrh_host(13,2) = 0.060281218;
        d_mdrh_host(13,3) = -7.6161E-05;

        // SN + AN
        d_mdrh_host(14,0) = 1916.072988;
        d_mdrh_host(14,1) = -20.85594868;
        d_mdrh_host(14,2) = 0.081140141;
        d_mdrh_host(14,3) = -1.07954274796e-4;

        // SN + AN + AC
        d_mdrh_host(15,0) = 1467.165935;
        d_mdrh_host(15,1) = -16.01166196;
        d_mdrh_host(15,2) = 0.063505582;
        d_mdrh_host(15,3) = -8.66722E-05;

        // SN + SC
        d_mdrh_host(16,0) = 158.447059;
        d_mdrh_host(16,1) = -0.628167358;
        d_mdrh_host(16,2) = 0.002014448;
        d_mdrh_host(16,3) = -3.13037E-06;

        // SN + SC + AC
        d_mdrh_host(17,0) = 1115.892468;
        d_mdrh_host(17,1) = -11.76936534;
        d_mdrh_host(17,2) = 0.045577399;
        d_mdrh_host(17,3) = -6.05779E-05;

        // SS + AC
        d_mdrh_host(18,0) = 269.5432407;
        d_mdrh_host(18,1) = -1.319963885;
        d_mdrh_host(18,2) = 0.002592363;
        d_mdrh_host(18,3) = -1.44479E-06;

        // SS + AN
        d_mdrh_host(19,0) = 2841.334784;
        d_mdrh_host(19,1) = -31.1889487;
        d_mdrh_host(19,2) = 0.118809274;
        d_mdrh_host(19,3) = -1.53007e-4;

        // SS + AN + AC
        d_mdrh_host(20,0) = 2199.36914;
        d_mdrh_host(20,1) = -24.11926569;
        d_mdrh_host(20,2) = 0.092932361;
        d_mdrh_host(20,3) = -1.21774e-4;

        // SS + AS
        d_mdrh_host(21,0) = 395.0051604;
        d_mdrh_host(21,1) = -2.521101657;
        d_mdrh_host(21,2) = 0.006139319;
        d_mdrh_host(21,3) = -4.43756E-06;

        //SS + AS + AC
        d_mdrh_host(22,0) = 386.5150675;
        d_mdrh_host(22,1) = -2.4632138;
        d_mdrh_host(22,2) = 0.006139319;
        d_mdrh_host(22,3) = -4.98796E-06;

        // SS + AS + AN
        d_mdrh_host(23,0) = 3101.538491;
        d_mdrh_host(23,1) = -34.19978105;
        d_mdrh_host(23,2) = 0.130118605;
        d_mdrh_host(23,3) = -1.66873e-4;

        // SS + AS + AN + AC
        d_mdrh_host(24,0) = 2307.579403;
        d_mdrh_host(24,1) = -25.43136774;
        d_mdrh_host(24,2) = 0.098064728;
        d_mdrh_host(24,3) = -1.28301e-4;

        // SS + SC
        d_mdrh_host(25,0) = 291.8309602;
        d_mdrh_host(25,1) = -1.828912974;
        d_mdrh_host(25,2) = 0.005053148;
        d_mdrh_host(25,3) = -4.57516E-06;

        // SS + SC + AC
        d_mdrh_host(26,0) = 188.3914345;
        d_mdrh_host(26,1) = -0.631345031;
        d_mdrh_host(26,2) = 0.000622807;
        d_mdrh_host(26,3) = 4.47196E-07;

        // SS + SN
        d_mdrh_host(27,0) = -167.1252839;
        d_mdrh_host(27,1) = 2.969828002;
        d_mdrh_host(27,2) = -0.010637255;
        d_mdrh_host(27,3) = 1.13175E-05;

        // SS + SN + AC
        d_mdrh_host(28,0) = 1516.782768;
        d_mdrh_host(28,1) = -15.7922661;
        d_mdrh_host(28,2) = 0.058942209;
        d_mdrh_host(28,3) = -7.5301E-05;

        // SS + SN + AN
        d_mdrh_host(29,0) = 1739.963163;
        d_mdrh_host(29,1) = -19.06576022;
        d_mdrh_host(29,2) = 0.07454963;
        d_mdrh_host(29,3) = -9.94302E-05;

        // SS + SN + AN + AC
        d_mdrh_host(30,0) = 2152.104877;
        d_mdrh_host(30,1) = -23.74998008;
        d_mdrh_host(30,2) = 0.092256654;
        d_mdrh_host(30,3) = -1.21953e-4;

        // SS + SN + SC
        d_mdrh_host(31,0) = 221.9976265;
        d_mdrh_host(31,1) = -1.311331272;
        d_mdrh_host(31,2) = 0.004406089;
        d_mdrh_host(31,3) = -5.88235E-06;

        // SS + SN + SC + AC
        d_mdrh_host(32,0) = 1205.645615;
        d_mdrh_host(32,1) = -12.71353459;
        d_mdrh_host(32,2) = 0.048803922;
        d_mdrh_host(32,3) = -6.41899E-05;

        // CC + AC
        d_mdrh_host(33,0) = 506.6737879;
        d_mdrh_host(33,1) = -3.723520818;
        d_mdrh_host(33,2) = 0.010814242;
        d_mdrh_host(33,3) = -1.21087E-05;

        // CC + SC
        d_mdrh_host(34,0) = -1123.523841;
        d_mdrh_host(34,1) = 14.08345977;
        d_mdrh_host(34,2) = -0.053687823;
        d_mdrh_host(34,3) = 6.52219E-05;

        // CC + SC + AC
        d_mdrh_host(35,0) = -1159.98607;
        d_mdrh_host(35,1) = 14.44309169;
        d_mdrh_host(35,2) = -0.054841073;
        d_mdrh_host(35,3) = 6.64259E-05;

        // CN + AC
        d_mdrh_host(36,0) = 756.0747916;
        d_mdrh_host(36,1) = -8.546826257;
        d_mdrh_host(36,2) = 0.035798677;
        d_mdrh_host(36,3) = -5.06629E-05;

        // CN + AN
        d_mdrh_host(37,0) = 338.668191;
        d_mdrh_host(37,1) = -2.971223403;
        d_mdrh_host(37,2) = 0.012294866;
        d_mdrh_host(37,3) = -1.87558E-05;

        // CN + AN + AC
        d_mdrh_host(38,0) = -53.18033508;
        d_mdrh_host(38,1) = 0.663911748;
        d_mdrh_host(38,2) = 9.16326e-4;
        d_mdrh_host(38,3) = -6.70354E-06;

        // CN + SC
        d_mdrh_host(39,0) = 3623.831129;
        d_mdrh_host(39,1) = -39.27226457;
        d_mdrh_host(39,2) = 0.144559515;
        d_mdrh_host(39,3) = -1.78159e-4;

        // CN + SC + AC
        d_mdrh_host(40,0) = 3436.656743;
        d_mdrh_host(40,1) = -37.16192684;
        d_mdrh_host(40,2) = 0.136641377;
        d_mdrh_host(40,3) = -1.68262e-4;

        // CN + SN
        d_mdrh_host(41,0) = 768.608476;
        d_mdrh_host(41,1) = -8.051517149;
        d_mdrh_host(41,2) = 0.032342332;
        d_mdrh_host(41,3) = -4.52224E-05;

        // CN + SN + AC
        d_mdrh_host(42,0) = 33.58027951;
        d_mdrh_host(42,1) = -0.308772182;
        d_mdrh_host(42,2) = 0.004713639;
        d_mdrh_host(42,3) = -1.19658E-05;

        // CN + SN + AN
        d_mdrh_host(43,0) = 57.80183041;
        d_mdrh_host(43,1) = 0.215264604;
        d_mdrh_host(43,2) = 4.11406e-4;
        d_mdrh_host(43,3) = -4.30702E-06;

        // CN + SN + AN + AC
        d_mdrh_host(44,0) = -234.368984;
        d_mdrh_host(44,1) = 2.721045204;
        d_mdrh_host(44,2) = -0.006688341;
        d_mdrh_host(44,3) = 2.31729E-06;

        // CN + SN + SC
        d_mdrh_host(45,0) = 3879.080557;
        d_mdrh_host(45,1) = -42.13562874;
        d_mdrh_host(45,2) = 0.155235005;
        d_mdrh_host(45,3) = -1.91387e-4;

        // CN + SN + SC + AC
        d_mdrh_host(46,0) = 3600.576985;
        d_mdrh_host(46,1) = -39.0283489;
        d_mdrh_host(46,2) = 0.143710316;
        d_mdrh_host(46,3) = -1.77167e-4;

        // CN + CC
        d_mdrh_host(47,0) = -1009.729826;
        d_mdrh_host(47,1) = 12.9145339;
        d_mdrh_host(47,2) = -0.049811146;
        d_mdrh_host(47,3) = 6.09563E-05;

        // CN + CC + AC
        d_mdrh_host(48,0) = -577.0919514;
        d_mdrh_host(48,1) = 8.020324227;
        d_mdrh_host(48,2) = -0.031469556;
        d_mdrh_host(48,3) = 3.82181E-05;

        // CN + CC + SC
        d_mdrh_host(49,0) = -728.9983499;
        d_mdrh_host(49,1) = 9.849458215;
        d_mdrh_host(49,2) = -0.03879257;
        d_mdrh_host(49,3) = 4.78844E-05;

        // CN + CC + SC + AC
        d_mdrh_host(50,0) = -803.7026845;
        d_mdrh_host(50,1) = 10.61881494;
        d_mdrh_host(50,2) = -0.041402993;
        d_mdrh_host(50,3) = 5.08084E-05;

        // SULFATE-RICH SYSTEMS
        // AB
        d_mdrh_host(51,0) = -493.6190458;
        d_mdrh_host(51,1) = 6.747053851;
        d_mdrh_host(51,2) = -0.026955267;
        d_mdrh_host(51,3) = 3.45118E-05;

        // LV
        d_mdrh_host(52,0) = 53.37874093;
        d_mdrh_host(52,1) = 1.01368249;
        d_mdrh_host(52,2) = -0.005887513;
        d_mdrh_host(52,3) = 8.94393E-06;

        // SB
        d_mdrh_host(53,0) = 206.619047;
        d_mdrh_host(53,1) = -1.342735684;
        d_mdrh_host(53,2) = 0.003197691;
        d_mdrh_host(53,3) = -1.93603E-06;

        // AB + LV
        d_mdrh_host(54,0) = -493.6190458;
        d_mdrh_host(54,1) = 6.747053851;
        d_mdrh_host(54,2) = -0.026955267;
        d_mdrh_host(54,3) = 3.45118E-05;

        // AS + LV
        d_mdrh_host(55,0) = 53.37874093;
        d_mdrh_host(55,1) = 1.01368249;
        d_mdrh_host(55,2) = -0.005887513;
        d_mdrh_host(55,3) = 8.94393E-06;

        // SS + SB
        d_mdrh_host(56,0) = 206.619047;
        d_mdrh_host(56,1) = -1.342735684;
        d_mdrh_host(56,2) = 0.003197691;
        d_mdrh_host(56,3) = -1.93603E-06;

        // SS + LV
        d_mdrh_host(57,0) = 41.7619047;
        d_mdrh_host(57,1) = 1.303872053;
        d_mdrh_host(57,2) = -0.007647908;
        d_mdrh_host(57,3) = 1.17845E-05;

        // SS + AS + LV
        d_mdrh_host(58,0) = 41.7619047;
        d_mdrh_host(58,1) = 1.303872053;
        d_mdrh_host(58,2) = -0.007647908;
        d_mdrh_host(58,3) = 1.17845E-05;

        // SS + AB
        d_mdrh_host(59,0) = -369.7142842;
        d_mdrh_host(59,1) = 5.512878771;
        d_mdrh_host(59,2) = -0.02301948;
        d_mdrh_host(59,3) = 3.0303E-05;

        // SS + LV + AB
        d_mdrh_host(60,0) = -369.7142842;
        d_mdrh_host(60,1) = 5.512878771;
        d_mdrh_host(60,2) = -0.02301948;
        d_mdrh_host(60,3) = 3.0303E-05;

        // SB + AB
        d_mdrh_host(61,0) = -162.8095232;
        d_mdrh_host(61,1) = 2.399326592;
        d_mdrh_host(61,2) = -0.009336219;
        d_mdrh_host(61,3) = 1.17845E-05;

        // SS + SB + AB
        d_mdrh_host(62,0) = -735.4285689;
        d_mdrh_host(62,1) = 8.885521857;
        d_mdrh_host(62,2) = -0.033488456;
        d_mdrh_host(62,3) = 4.12458E-05;

        // molecular weights for aerosol species
        mw_aer_mac_host(iso4_a) =  96.0;
        mw_aer_mac_host(ino3_a) =  62.0;
        mw_aer_mac_host(icl_a)  =  35.5;
        mw_aer_mac_host(imsa_a) =  95.0;    // CH3SO3
        mw_aer_mac_host(ico3_a) =  60.0;
        mw_aer_mac_host(inh4_a) =  18.0;
        mw_aer_mac_host(ina_a)  =  23.0;
        mw_aer_mac_host(ica_a)  =  40.0;
        mw_aer_mac_host(ioin_a) =   1.0;		// not used
        mw_aer_mac_host(ibc_a)  =   1.0;		// not used
        mw_aer_mac_host(ioc_a)  =   1.0;	  // 200 assumed for primary organics
        mw_aer_mac_host(iaro1_a)= 150.0;
        mw_aer_mac_host(iaro2_a)= 150.0;
        mw_aer_mac_host(ialk1_a)= 140.0;
        mw_aer_mac_host(iole1_a)= 140.0;
        mw_aer_mac_host(iapi1_a)= 184.0;
        mw_aer_mac_host(iapi2_a)= 184.0;
        mw_aer_mac_host(ilim1_a)= 200.0;
        mw_aer_mac_host(ilim2_a)= 200.0;

        // densities for aerosol species
        dens_aer_mac_host(iso4_a) = 1.8;	// used
        dens_aer_mac_host(ino3_a) = 1.8;	// used
        dens_aer_mac_host(icl_a)  = 2.2;	// used
        dens_aer_mac_host(imsa_a) = 1.8;	// used
        dens_aer_mac_host(ico3_a) = 2.6;	// used
        dens_aer_mac_host(inh4_a) = 1.8;	// used
        dens_aer_mac_host(ina_a)  = 2.2;	// used
        dens_aer_mac_host(ica_a)  = 2.6;	// used
        dens_aer_mac_host(ioin_a) = 2.6;	// used
        dens_aer_mac_host(ioc_a)  = 1.0;	// used
        dens_aer_mac_host(ibc_a)  = 1.8;	// used
        dens_aer_mac_host(iaro1_a)= 1.4;
        dens_aer_mac_host(iaro2_a)= 1.4;
        dens_aer_mac_host(ialk1_a)= 1.4;
        dens_aer_mac_host(iole1_a)= 1.4;
        dens_aer_mac_host(iapi1_a)= 1.4;
        dens_aer_mac_host(iapi2_a)= 1.4;
        dens_aer_mac_host(ilim1_a)= 1.4;
        dens_aer_mac_host(ilim2_a)= 1.4;

        // molecular weights for electrolytes
        mw_electrolyte_host(jnh4so4) = 132.0;
        mw_electrolyte_host(jlvcite) = 247.0;
        mw_electrolyte_host(jnh4hso4)= 115.0;
        mw_electrolyte_host(jnh4msa) = 113.0;
        mw_electrolyte_host(jnh4no3) =  80.0;
        mw_electrolyte_host(jnh4cl)  =  53.5;
        mw_electrolyte_host(jnacl)   =  58.5;
        mw_electrolyte_host(jnano3)  =  85.0;
        mw_electrolyte_host(jna2so4) = 142.0;
        mw_electrolyte_host(jna3hso4)= 262.0;
        mw_electrolyte_host(jnahso4) = 120.0;
        mw_electrolyte_host(jnamsa)  = 118.0;
        mw_electrolyte_host(jcaso4)  = 136.0;
        mw_electrolyte_host(jcamsa2) = 230.0;
        mw_electrolyte_host(jcano3)  = 164.0;
        mw_electrolyte_host(jcacl2)  = 111.0;
        mw_electrolyte_host(jcaco3)  = 100.0;
        mw_electrolyte_host(jh2so4)  =  98.0;
        mw_electrolyte_host(jhno3)   =  63.0;
        mw_electrolyte_host(jhcl)    =  36.5;
        mw_electrolyte_host(jmsa)    =  96.0;

        // magnitude of charges on ions
        zc_host(jc_h)   = 1;
        zc_host(jc_nh4) = 1;
        zc_host(jc_na)  = 1;
        zc_host(jc_ca)  = 2;

        za_host(ja_hso4) = 1;
        za_host(ja_so4)  = 2;
        za_host(ja_no3)  = 1;
        za_host(ja_cl)   = 1;
        za_host(ja_msa)  = 1;

        // molecular weights of ions [g/mol]
        // cations
        mw_c_host(jc_h)   = 1.0;
        mw_c_host(jc_nh4) = 18.0;
        mw_c_host(jc_na)  = 23.0;
        mw_c_host(jc_ca)  = 40.0;

        // anions
        mw_a_host(ja_so4)  = 96.0;
        mw_a_host(ja_hso4) = 97.0;
        mw_a_host(ja_no3)  = 62.0;
        mw_a_host(ja_cl)   = 35.5;
        mw_a_host(ja_msa)  = 95.0;

        a_zsr.modify_host();
        aw_min.modify_host();
        b_zsr.modify_host();
        Keq_298_ll.modify_host();
        Keq_a_ll.modify_host();
        Keq_b_ll.modify_host();
        b_mtem.modify_host();
        d_mdrh.modify_host();
        mw_aer_mac.modify_host();
        dens_aer_mac.modify_host();
        mw_electrolyte.modify_host();
        zc.modify_host();
        za.modify_host();
        mw_c.modify_host();
        mw_a.modify_host();

        a_zsr.sync_device();
        aw_min.sync_device();
        b_zsr.sync_device();
        Keq_298_ll.sync_device();
        Keq_a_ll.sync_device();
        Keq_b_ll.sync_device();
        b_mtem.sync_device();
        d_mdrh.sync_device();
        mw_aer_mac.sync_device();
        dens_aer_mac.sync_device();
        mw_electrolyte.sync_device();
        zc.sync_device();
        za.sync_device();
        mw_c.sync_device();
        mw_a.sync_device();
      }

      ~MosaicModelData() = default;
};

template<typename ValueType, typename DeviceType>
struct MOSAIC{

  using value_type = ValueType;
  using device_type = DeviceType;
  using scalar_type = typename ats<value_type>::scalar_type;

  using real_type = scalar_type;

  using value_type_1d_view_type = Tines::value_type_1d_view<value_type,device_type>;
  using value_type_2d_view_type = Tines::value_type_2d_view<value_type,device_type>;
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
  using real_type_2d_view_type = Tines::value_type_2d_view<real_type,device_type>;

  KOKKOS_INLINE_FUNCTION static ordinal_type getWorkSpaceSize() {
    ordinal_type workspace_size=0;
    return workspace_size;
  }

  KOKKOS_INLINE_FUNCTION static
  void check_aerosol_mass(const MosaicModelData<DeviceType>& mosaic,
                          const real_type_1d_view_type& aer_total,
                          real_type& mass_dry_a,
                          real_type& jphase,
                          real_type& jaerosolstate,
                          real_type& num_a) {

    mass_dry_a = 0.0;

    real_type aer_H = (2.0*aer_total(mosaic.iso4_a)  +
                           aer_total(mosaic.ino3_a)  +
                           aer_total(mosaic.icl_a)   +
                           aer_total(mosaic.imsa_a)  +
                       2.0*aer_total(mosaic.ico3_a)) -
                      (2.0*aer_total(mosaic.ica_a)   +
                           aer_total(mosaic.ina_a)   +
                           aer_total(mosaic.inh4_a));
    aer_H = max(aer_H, 0.0);

    auto mw_aer_mac = mosaic.mw_aer_mac.template view<DeviceType>();
    for (ordinal_type iaer = 0; iaer < mosaic.naer; iaer++) {
      mass_dry_a = mass_dry_a + aer_total(iaer)*mw_aer_mac(iaer); // ng/m^3(air)
    }
    mass_dry_a = mass_dry_a + aer_H;

    const real_type drymass = mass_dry_a; // ng/m^3(air)
    mass_dry_a = mass_dry_a*1.0e-15; // g/cc(air)

    if (drymass < mosaic.mass_cutoff) {
      jaerosolstate = mosaic.no_aerosol;
      jphase = 0;
      if (drymass == 0.0) {
        num_a = 0.0;
      }
    }
   } // check_aerosol_mass

  KOKKOS_INLINE_FUNCTION static
  void adjust_liquid_aerosol(const MosaicModelData<DeviceType>& mosaic,
                             const real_type_1d_view_type& aer_solid,
                             const real_type_1d_view_type& aer_liquid,
                             const real_type_1d_view_type& aer_total,
                             const real_type_1d_view_type& electrolyte_solid,
                             const real_type_1d_view_type& electrolyte_liquid,
                             const real_type_1d_view_type& electrolyte_total,
                             const real_type_1d_view_type& epercent_solid,
                             const real_type_1d_view_type& epercent_liquid,
                             const real_type_1d_view_type& epercent_total,
                             real_type& jphase, real_type& jhyst_leg) {

    jphase    = mosaic.jliquid;
    jhyst_leg = mosaic.jhyst_up; // upper curve

    // Partition all electrolytes into liquid phase
    for (ordinal_type je = 0; je < mosaic.nelectrolyte; ++je) {
      electrolyte_solid(je) = 0.0;
      epercent_solid(je) = 0.0;
      electrolyte_liquid(je) = electrolyte_total(je);
      epercent_liquid(je) = epercent_total(je);
    }

    // Except these electrolytes, which always remain in the solid phase
    electrolyte_solid(mosaic.jcaco3) = electrolyte_total(mosaic.jcaco3);
    electrolyte_solid(mosaic.jcaso4) = electrolyte_total(mosaic.jcaso4);
    epercent_solid(mosaic.jcaco3) = epercent_total(mosaic.jcaco3);
    epercent_solid(mosaic.jcaso4) = epercent_total(mosaic.jcaso4);
    electrolyte_liquid(mosaic.jcaco3) = 0.0;
    electrolyte_liquid(mosaic.jcaso4) = 0.0;
    epercent_liquid(mosaic.jcaco3) = 0.0;
    epercent_liquid(mosaic.jcaso4) = 0.0;

    // Partition all the aer species into solid and liquid phases
    // Solid phase
    aer_solid(mosaic.iso4_a) = electrolyte_solid(mosaic.jcaso4);
    aer_solid(mosaic.ino3_a) = 0.0;
    aer_solid(mosaic.icl_a) = 0.0;
    aer_solid(mosaic.inh4_a) = 0.0;
    aer_solid(mosaic.ioc_a) = aer_total(mosaic.ioc_a);
    aer_solid(mosaic.imsa_a) = 0.0;
    aer_solid(mosaic.ico3_a) = aer_total(mosaic.ico3_a);
    aer_solid(mosaic.ina_a) = 0.0;
    aer_solid(mosaic.ica_a) = electrolyte_solid(mosaic.jcaco3) + electrolyte_solid(mosaic.jcaso4);
    aer_solid(mosaic.ibc_a) = aer_total(mosaic.ibc_a);
    aer_solid(mosaic.ioin_a) = aer_total(mosaic.ioin_a);
    aer_solid(mosaic.iaro1_a) = aer_total(mosaic.iaro1_a);
    aer_solid(mosaic.iaro2_a) = aer_total(mosaic.iaro2_a);
    aer_solid(mosaic.ialk1_a) = aer_total(mosaic.ialk1_a);
    aer_solid(mosaic.iole1_a) = aer_total(mosaic.iole1_a);
    aer_solid(mosaic.iapi1_a) = aer_total(mosaic.iapi1_a);
    aer_solid(mosaic.iapi2_a) = aer_total(mosaic.iapi2_a);
    aer_solid(mosaic.ilim1_a) = aer_total(mosaic.ilim1_a);
    aer_solid(mosaic.ilim2_a) = aer_total(mosaic.ilim2_a);

    // Liquid phase
    aer_liquid(mosaic.iso4_a) = max(0.0, aer_total(mosaic.iso4_a) - aer_solid(mosaic.iso4_a));
    aer_liquid(mosaic.ino3_a) = aer_total(mosaic.ino3_a);
    aer_liquid(mosaic.icl_a) = aer_total(mosaic.icl_a);
    aer_liquid(mosaic.inh4_a) = aer_total(mosaic.inh4_a);
    aer_liquid(mosaic.ioc_a) = 0.0;
    aer_liquid(mosaic.imsa_a) = aer_total(mosaic.imsa_a);
    aer_liquid(mosaic.ico3_a) = 0.0;
    aer_liquid(mosaic.ina_a) = aer_total(mosaic.ina_a);
    aer_liquid(mosaic.ica_a) = max(0.0, aer_total(mosaic.ica_a) - aer_solid(mosaic.ica_a));
    aer_liquid(mosaic.ibc_a) = 0.0;
    aer_liquid(mosaic.ioin_a) = 0.0;
    aer_liquid(mosaic.iaro1_a) = 0.0;
    aer_liquid(mosaic.iaro2_a) = 0.0;
    aer_liquid(mosaic.ialk1_a) = 0.0;
    aer_liquid(mosaic.iole1_a) = 0.0;
    aer_liquid(mosaic.iapi1_a) = 0.0;
    aer_liquid(mosaic.iapi2_a) = 0.0;
    aer_liquid(mosaic.ilim1_a) = 0.0;
    aer_liquid(mosaic.ilim2_a) = 0.0;
  } // adjust_liquid_aerosol

  KOKKOS_INLINE_FUNCTION static
  void adjust_solid_aerosol(const MosaicModelData<DeviceType>& mosaic,
                            const real_type_1d_view_type& aer_solid,
                            const real_type_1d_view_type& aer_liquid,
                            const real_type_1d_view_type& aer_total,
                            const real_type_1d_view_type& electrolyte_solid,
                            const real_type_1d_view_type& electrolyte_liquid,
                            const real_type_1d_view_type& electrolyte_total,
                            const real_type_1d_view_type& epercent_solid,
                            const real_type_1d_view_type& epercent_liquid,
                            const real_type_1d_view_type& epercent_total,
                            real_type& water_a, real_type& jphase, real_type& jhyst_leg) {

    jphase = mosaic.jsolid;
    jhyst_leg = mosaic.jhyst_lo;
    water_a = 0.0;

    // Transfer aer(total) to aer(jsolid) and set aer(liquid) to 0
    for (ordinal_type iaer = 0; iaer < mosaic.naer; iaer++) {
        aer_solid(iaer) = aer_total(iaer);
        aer_liquid(iaer) = 0.0;
    }

    // Transfer electrolyte(total) to electrolyte(solid) and set electrolyte(liquid) to 0
    for (ordinal_type je = 0; je < mosaic.nelectrolyte; je++) {
      electrolyte_liquid(je) = 0.0;
      epercent_liquid(je) = 0.0;
      electrolyte_solid(je) = electrolyte_total(je);
      epercent_solid(je) = epercent_total(je);
    }

    // Update aer(total) that may have been affected above
    aer_total(mosaic.inh4_a) = aer_solid(mosaic.inh4_a);
    aer_total(mosaic.ino3_a) = aer_solid(mosaic.ino3_a);
    aer_total(mosaic.icl_a) = aer_solid(mosaic.icl_a);

    // Update electrolyte(total)
    for (ordinal_type je = 0; je < mosaic.nelectrolyte; je++) {
      electrolyte_total(je) = electrolyte_solid(je);
      epercent_total(je) = epercent_solid(je);
    }
  } // adjust_solid_aerosol

  KOKKOS_INLINE_FUNCTION static
  void do_full_deliquescence(const MosaicModelData<DeviceType>& mosaic,
                             const real_type_1d_view_type& electrolyte_solid,
                             const real_type_1d_view_type& electrolyte_liquid,
                             const real_type_1d_view_type& electrolyte_total,
                             const real_type_1d_view_type& aer_solid,
                             const real_type_1d_view_type& aer_liquid,
                             const real_type_1d_view_type& aer_total) {

    // Partition all electrolytes to the liquid phase
    for (ordinal_type js = 0; js < mosaic.nelectrolyte; js++) {
        electrolyte_solid(js) = 0.0;
        electrolyte_liquid(js) = electrolyte_total(js);
    }

    // Except these electrolytes, which always remain in the solid phase
    electrolyte_solid(mosaic.jcaco3) = electrolyte_total(mosaic.jcaco3);
    electrolyte_solid(mosaic.jcaso4) = electrolyte_total(mosaic.jcaso4);
    electrolyte_liquid(mosaic.jcaco3) = 0.0;
    electrolyte_liquid(mosaic.jcaso4) = 0.0;

    // Partition all the generic aer species into solid and liquid phases
    // Solid phase
    aer_solid(mosaic.iso4_a) = electrolyte_solid(mosaic.jcaso4);
    aer_solid(mosaic.ino3_a) = 0.0;
    aer_solid(mosaic.icl_a) = 0.0;
    aer_solid(mosaic.inh4_a) = 0.0;
    aer_solid(mosaic.ioc_a) = aer_total(mosaic.ioc_a);
    aer_solid(mosaic.imsa_a) = 0.0;
    aer_solid(mosaic.ico3_a) = aer_total(mosaic.ico3_a);
    aer_solid(mosaic.ina_a) = 0.0;
    aer_solid(mosaic.ica_a) = electrolyte_solid(mosaic.jcaco3) + electrolyte_solid(mosaic.jcaso4);
    aer_solid(mosaic.ibc_a) = aer_total(mosaic.ibc_a);
    aer_solid(mosaic.ioin_a) = aer_total(mosaic.ioin_a);
    aer_solid(mosaic.iaro1_a) = aer_total(mosaic.iaro1_a);
    aer_solid(mosaic.iaro2_a) = aer_total(mosaic.iaro2_a);
    aer_solid(mosaic.ialk1_a) = aer_total(mosaic.ialk1_a);
    aer_solid(mosaic.iole1_a) = aer_total(mosaic.iole1_a);
    aer_solid(mosaic.iapi1_a) = aer_total(mosaic.iapi1_a);
    aer_solid(mosaic.iapi2_a) = aer_total(mosaic.iapi2_a);
    aer_solid(mosaic.ilim1_a) = aer_total(mosaic.ilim1_a);
    aer_solid(mosaic.ilim2_a) = aer_total(mosaic.ilim2_a);

    // Liquid phase
    aer_liquid(mosaic.iso4_a) = aer_total(mosaic.iso4_a) - electrolyte_solid(mosaic.jcaso4);
    aer_liquid(mosaic.ino3_a) = aer_total(mosaic.ino3_a);
    aer_liquid(mosaic.icl_a) = aer_total(mosaic.icl_a);
    aer_liquid(mosaic.inh4_a) = aer_total(mosaic.inh4_a);
    aer_liquid(mosaic.ioc_a) = 0.0;
    aer_liquid(mosaic.imsa_a) = aer_total(mosaic.imsa_a);
    aer_liquid(mosaic.ico3_a) = 0.0;
    aer_liquid(mosaic.ina_a) = aer_total(mosaic.ina_a);
    aer_liquid(mosaic.ica_a) = electrolyte_total(mosaic.jcano3) + electrolyte_total(mosaic.jcacl2);
    aer_liquid(mosaic.ibc_a) = 0.0;
    aer_liquid(mosaic.ioin_a) = 0.0;
    aer_liquid(mosaic.iaro1_a) = 0.0;
    aer_liquid(mosaic.iaro2_a) = 0.0;
    aer_liquid(mosaic.ialk1_a) = 0.0;
    aer_liquid(mosaic.iole1_a) = 0.0;
    aer_liquid(mosaic.iapi1_a) = 0.0;
    aer_liquid(mosaic.iapi2_a) = 0.0;
    aer_liquid(mosaic.ilim1_a) = 0.0;
    aer_liquid(mosaic.ilim2_a) = 0.0;
  } // do_full_deliquescence

  KOKKOS_INLINE_FUNCTION static
  void calculate_XT(const MosaicModelData<DeviceType>& mosaic,
                    const real_type_1d_view_type& aer,
                    real_type &XT) {

    //aer nmol/m^3
    if ((aer(mosaic.iso4_a) + aer(mosaic.imsa_a)) > 0.0) {
        XT = (aer(mosaic.inh4_a) + aer(mosaic.ina_a) +
             2.0 * aer(mosaic.ica_a)) /
            (aer(mosaic.iso4_a) + 0.5 * aer(mosaic.imsa_a));
    } else {
        XT = -1.0;
    }
  }// calculate_XT

  KOKKOS_INLINE_FUNCTION static
  void fnlog_gamZ(const MosaicModelData<DeviceType>& mosaic,
                  const ordinal_type& jA,
                  const ordinal_type& jE,
                  const real_type& aH2O,
                  real_type& log_gamZ_) {

    // FIXME: aH2O should not be local; make sure updated with RH

    auto b_mtem = mosaic.b_mtem.template view<DeviceType>();
    auto aw_min = mosaic.aw_min.template view<DeviceType>();

    const real_type aw = max(aH2O, aw_min(jE));

    log_gamZ_ = b_mtem(0,jA,jE) + aw *
               (b_mtem(1,jA,jE) + aw *
               (b_mtem(2,jA,jE) + aw *
               (b_mtem(3,jA,jE) + aw *
               (b_mtem(4,jA,jE) + aw *
                b_mtem(5,jA,jE) ))));
  } // fnlog_gamZ

  KOKKOS_INLINE_FUNCTION static
  void fuchs_sutugin(const real_type& rkn,
                     const real_type& a,
                     real_type& fuchs_sut) {

    const real_type rnum  = 0.75*a*(1. + rkn);
    const real_type denom = ats<real_type>::pow(rkn, 2.0) + rkn + 0.283*rkn*a + 0.75*a;
    fuchs_sut = rnum/denom;
  } // fuchs_sutugin

  KOKKOS_INLINE_FUNCTION static
  void gas_diffusivity(const real_type& T_K,
                       const real_type& P,
                       const real_type& MW,
                       const real_type& Vm,
                       real_type& gas_diff) {
    gas_diff = (1.0e-3 * ats<real_type>::pow(T_K, 1.75) * ats<real_type>::sqrt(1.0/MW + 0.035))/
                    (P * ats<real_type>::pow(ats<real_type>::pow(Vm, 0.333333) + 2.7189, 2.0));
  } // gas_diffusivity

  KOKKOS_INLINE_FUNCTION static
  void mean_molecular_speed(const real_type& T_K,
                            const real_type& MW,
                            real_type& mean_molec_speed) {

    mean_molec_speed = 1.455e4 * ats<real_type>::sqrt(T_K/MW);
  } // mean_molecular_speed

  KOKKOS_INLINE_FUNCTION static
  void fn_Keq(const real_type& Keq_298,
              const real_type& a,
              const real_type& b,
              const real_type& T,
              real_type& Keq) {

    real_type tt = 298.15/T;

    Keq = Keq_298*ats<real_type>::exp(a*(tt - 1.0) + b*(1.0 + ats<real_type>::log(tt) - tt));
  } // fn_Keq

  KOKKOS_INLINE_FUNCTION static
  void fn_Po(const real_type& Po_298,
             const real_type& DH,
             const real_type& T_K,
             real_type& Po) {

    // Van't Hoff Equation
    Po = Po_298*ats<real_type>::exp(-(DH/(RUNIV/1000))*(1.0/T_K - (1/298.15)));
  } // fn_Po

  KOKKOS_INLINE_FUNCTION static
  void drh_mutual(const MosaicModelData<DeviceType>& mosaic,
                  const ordinal_type& j_index,
                  const real_type& T_K,
                  real_type& drh_mutual) {

    if (j_index == 7 || j_index == 8 || (j_index >= 34 && j_index <= 51)) {
      drh_mutual = 10.0; // CaNO3 or CaCl2 containing mixtures
    } else {
      auto d_mdrh = mosaic.d_mdrh.template view<DeviceType>();
      drh_mutual = d_mdrh(j_index,0) + T_K *
                  (d_mdrh(j_index,1) + T_K *
                  (d_mdrh(j_index,2) + T_K *
                   d_mdrh(j_index,3) )) + 1.0;
    }
  } // drh_mutual

  KOKKOS_INLINE_FUNCTION static
  void molality_0(const MosaicModelData<DeviceType>& mosaic,
                  const ordinal_type& je,
                  real_type& aw,
                  real_type& molality) {

    auto aw_min = mosaic.aw_min.template view<DeviceType>();
    auto a_zsr = mosaic.a_zsr.template view<DeviceType>();

    aw = max(aw, aw_min(je));
    aw = min(aw, 0.999999);

    if (aw < 0.97) {

      real_type xm = a_zsr(0,je) +
                 aw*(a_zsr(1,je) +
                 aw*(a_zsr(2,je) +
                 aw*(a_zsr(3,je) +
                 aw*(a_zsr(4,je) +
                 aw* a_zsr(5,je) ))));

        molality = 55.509*xm/(1.0 - xm);
    } else {
      auto b_zsr = mosaic.b_zsr.template view<DeviceType>();
      molality = -b_zsr(je)*ats<real_type>::log(aw);
    }
  } // molality_0

  KOKKOS_INLINE_FUNCTION static
  void bin_molality(const MosaicModelData<DeviceType>& mosaic,
                    const ordinal_type& je,
                    const real_type& aH2O_a,
                    real_type& molality) {

    auto aw_min = mosaic.aw_min.template view<DeviceType>();
    auto a_zsr = mosaic.a_zsr.template view<DeviceType>();

    real_type aw = max(aH2O_a, aw_min(je));
    aw = min(aw, 0.999999);

    if (aw < 0.97) {

      real_type xm = a_zsr(0,je) +
                 aw*(a_zsr(1,je) +
                 aw*(a_zsr(2,je) +
                 aw*(a_zsr(3,je) +
                 aw*(a_zsr(4,je) +
                 aw* a_zsr(5,je) ))));

        molality = 55.509*xm/(1.0 - xm);
    } else {
      auto b_zsr = mosaic.b_zsr.template view<DeviceType>();
      molality = -b_zsr(je)*ats<real_type>::log(aw);
    }
  } // bin_molality

  KOKKOS_INLINE_FUNCTION static
  void bin_molality_60(const MosaicModelData<DeviceType>& mosaic,
                    const ordinal_type& je,
                    real_type& molality) {

    auto a_zsr = mosaic.a_zsr.template view<DeviceType>();

    const real_type aw = 0.6;

    real_type xm = a_zsr(0,je) +
               aw*(a_zsr(1,je) +
               aw*(a_zsr(2,je) +
               aw*(a_zsr(3,je) +
               aw*(a_zsr(4,je) +
               aw* a_zsr(5,je) ))));

    molality = 55.509*xm/(1.0 - xm);
  } // bin_molality_60

  KOKKOS_INLINE_FUNCTION static
  void MTEM_compute_log_gamZ(const MosaicModelData<DeviceType>& mosaic,
                             const real_type& aH2O,
                             real_type_2d_view_type log_gamZ) {

    real_type log_gamZ_ = 0.0;

  // sulfate-poor species
    ordinal_type jA = mosaic.jhno3;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;

    jA = mosaic.jhcl;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;

    jA = mosaic.jnh4so4;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;

    jA = mosaic.jnh4no3;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jnh4cl;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jna2so4;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;

    jA = mosaic.jnano3;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jnacl;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jcano3;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jcacl2;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

  // sulfate-rich species
    jA = mosaic.jh2so4;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jhhso4;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jnh4hso4;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jlvcite;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jnahso4;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jna3hso4;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;
  } // MTEM_compute_log_gamZ

  KOKKOS_INLINE_FUNCTION static
  void aerosol_water_up(const MosaicModelData<DeviceType>& mosaic,
                        const real_type_1d_view_type& electrolyte_total,
                        real_type& aerosol_water) {

    real_type dum = 0.0;

    // TODO: change this for loop to a parallel_reduce
    for (ordinal_type je = 0; je < mosaic.nsalt + 4; je++) {
      real_type molality = 0.0;
      bin_molality_60(mosaic, je, molality);
      dum += 1.e-9*electrolyte_total(je)/molality;
    }

    aerosol_water = dum;
  } // aerosol_water_up

};

} // namespace Impl
} // namespace TChem

#endif
