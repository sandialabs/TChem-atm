/* =====================================================================================
TChem-atm version 1.0
Copyright (2024) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.

This file is part of TChem-atm. TChem-atm is open source software: you can redistribute it
and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Contact Oscar Diaz-Ibarra at <odiazib@sandia.gov>, or
           Mike Schmidt at <mjschm@sandia.gov>, or
           Cosmin Safta at <csafta@sandia.gov>

Sandia National Laboratories, New Mexico/Livermore, NM/CA, USA
===================================================================================== */
#ifndef __TCHEM_KINETIC_MODELMDATA_HPP__
#define __TCHEM_KINETIC_MODELMDATA_HPP__

#include "TChem_Util.hpp"

#if defined(TCHEM_ENABLE_TPL_YAML_CPP)
#include "yaml-cpp/yaml.h"
#include <iostream>
#endif

namespace TChem {
  struct KineticModelData {
  public:

    // aerosol chemistry
    // number of species
    ordinal_type nSpec_;

    /* Reaction data */
    /* is reaction reversible ? */
    ordinal_type nReac_{0};

    /** \var ordinal_type maxSpecInReac_
     *  \ingroup maxpar
     *  \brief Maximum number of species in a reaction */
    ordinal_type maxSpecInReac_;

     /* Pressure dependent reactions */
    ordinal_type nFallReac_;

    /* Arrhenius parameters, forward and reverse */
    real_type_2d_dual_view reacArhenFor_, reacArhenRev_;

    /* stoichiometric coeffs and reactants and product indices */
    ordinal_type_2d_dual_view reacSidx_;
    real_type_2d_dual_view reacNuki_;

    // FIXME: add a description of these types
    real_type_2d_dual_view reacPpar_;
    ordinal_type_1d_dual_view reacPfal_;


    ordinal_type_1d_dual_view isRev_;

    /* no. of reac+prod, no. of reac only, no. of prod only, stoichiom.coef.
     * indicator */
    ordinal_type_1d_dual_view reacNrp_, reacNreac_, reacNprod_, reacScoef_;

    string_type_1d_dual_view<LENGTHOFSPECNAME + 1> sNames_;

    arrhenius_reaction_type_1d_dual_view ArrheniusCoef_;
    cmaq_h2o2_type_1d_dual_view CMAQ_H2O2Coef_;
    jpl_reaction_type_1d_dual_view JPL_Coef_;

    r_jpl_arrhenius_reaction_type_1d_dual_view R_JPL_ArrheniusCoef_;
    adjust_reaction_type_1d_dual_view adjust_reaction_;

    // modifier for uci1, uci2, and uci3
    prod_O1D_type_1d_dual_view prod_O1D_;
    ordinal_type number_of_prod_O1D_{0};

    emission_source_type_1d_dual_view EmissionCoef_;
    ordinal_type  nConstSpec_;
    real_type CONV_PPM_;
    ordinal_type_1d_dual_view convExponent_;
    // index of M species
    ordinal_type M_index_;
    // return index of gas species
    std::map<std::string, int> species_indx_;

  public:
      /**
       * constructor that allows use of a yaml file and custom echo and error streams
       * @param mechfile
       * @param echofile
       */
      KineticModelData(const std::string& mechfile, std::ostream& echofile);
      KineticModelData(const std::string& mechfile);

    /// constructor and destructor
    KineticModelData() = default;
    KineticModelData(const KineticModelData& b) = default;
    ~KineticModelData() = default;

#if defined(TCHEM_ATM_ENABLE_TPL_YAML_CPP)
    void initYamlFile(const std::string &mechfile, std::ostream& echofile);
    ordinal_type initChemNCAR(YAML::Node& doc, std::ostream& echofile);
#endif
    /// create multiple models sharing the data in this model
    kmd_type_1d_view_host clone(const int n_models);

  };



} // namespace TChem
#endif

#include "TChem_KineticModelNCAR_ConstData.hpp"
