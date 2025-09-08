/* =====================================================================================
TChem-atm version 2.0.0
Copyright (2025) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2025 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.

This file is part of TChem-atm. TChem-atm is open source software: you can redistribute it
and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Contact Cosmin Safta at <csafta@sandia.gov> or,, or
           Nicole Riemer at <nriemer@illinois.edu> or,
           Matthew West at <mwest@illinois.edu>
         Oscar Diaz-Ibarra at <odiazib@sandia.gov>

Sandia National Laboratories, New Mexico/Livermore, NM/CA, USA
===================================================================================== */
#include "TChem_KineticModelData.hpp"

#if defined(TCHEM_ATM_ENABLE_TPL_YAML_CPP)
#include <cstdarg>
#include <algorithm>
#endif

namespace TChem {
#if defined(TCHEM_ATM_ENABLE_TPL_YAML_CPP)

/**
 * Helper function to print to the stream.
 * @param stream
 * @param buffer
 * @param format
 * @param ...
 */
void StreamPrint(std::ostream& stream, const char* format, ...)
{
    va_list args;
    va_start(args, format);
    auto len = vsnprintf(nullptr,0, format, args);
    va_end(args);
    std::vector<char> vec(len + 1);
    va_start(args, format);
    std::vsnprintf(&vec[0], len + 1, format, args);
    va_end(args);
    va_end(args);
    stream << vec.data();
}

KineticModelData::KineticModelData(const std::string &mechfile, std::ostream& echofile) {
    initYamlFile(mechfile, echofile);
}

KineticModelData::KineticModelData(const std::string &mechfile) {
    std::ofstream echofile;
    echofile.open("kmod.echo");
    initYamlFile(mechfile,echofile);
    echofile.close();
}

void KineticModelData::initYamlFile(const std::string &mechfile,
                                    std::ostream& echofile){

  YAML::Node doc = YAML::LoadFile(mechfile);
  // FIXME: add error checking in yaml parser
  if (doc["NCAR-version"]) {
    if (verboseEnabled) {
      printf("Using parser for NCAR atmospheric chemistry\n");
    }
    initChemNCAR(doc, echofile);
  } else {
    printf("we do not have a parser for this file\n");
    // exit
  }
}

int KineticModelData::initChemNCAR(YAML::Node &root, std::ostream& echofile) {

#define DASHLINE(stream)                                                                                                 \
  stream << "------------------------------------------------------------"                                         \
                "-------------\n";

  CONV_PPM_ = CONV_PPM;

  auto species_names = root["species"];
  auto reactions = root["reactions"];

  YAML::Node const_species_names;
  nConstSpec_ = 0;
  if (root["constant_species"]) {
    const_species_names = root["constant_species"];
    nConstSpec_ = const_species_names.size();
  }
  // total number of species includes constant species
  // constant species are tracers
  nSpec_ = species_names.size() + nConstSpec_;
  nReac_ = reactions.size();

  /* Species' name */
  sNames_ = string_type_1d_dual_view<LENGTHOFSPECNAME + 1>(do_not_init_tag("KMD::sNames"), nSpec_);

  auto sNamesHost = sNames_.view_host();

  // std::map<std::string, int> species_indx_;
  int sp_i(0);
  // non constant species
  for (auto const &sp : species_names) {
    std::string sp_name = sp["name"].as<std::string>();
    // std::transform(sp_name.begin(),
    // sp_name.end(),sp_name.begin(), ::toupper);
    species_indx_.insert(std::pair<std::string, int>(sp_name, sp_i));
    char *specNm = &*sp_name.begin();
    strncat(&sNamesHost(sp_i, 0), specNm, LENGTHOFSPECNAME);
    sp_i++;
  }
  // do not reset sp_i, add const species at end of species list
  //
  // constant species
  for (auto const &sp : const_species_names) {

    std::string sp_name = sp["name"].as<std::string>();
    // std::transform(sp_name.begin(),
    // sp_name.end(),sp_name.begin(), ::toupper);
    species_indx_.insert(std::pair<std::string, int>(sp_name, sp_i));
    char *specNm = &*sp_name.begin();
    strncat(&sNamesHost(sp_i, 0), specNm, LENGTHOFSPECNAME);
    sp_i++;
  }

  // M species index
  {
    //Note: Name of M species is hard-coded to M
  auto it = species_indx_.find("M");

  if (it != species_indx_.end()) {
    M_index_ = it->second;
    // reacSidxHost(i, count) = it->second;
  } else {
    printf("Yaml : Error when interpreting kinetic model  !!!");
    printf("M species does not exit; Note: Name of M species is hard-coded to M \n");
      exit(1);
  }
  }

  StreamPrint(echofile, "kmod.list : # of species                                      : %d\n", nSpec_);
  StreamPrint(echofile, "kmod.list : # of species with constant concentration          : %d\n", nConstSpec_);
  StreamPrint(echofile, "kmod.list : # of reactions                                    : %d\n", nReac_);
  //
  StreamPrint(echofile, "No. \t Species\n");
  for (int i = 0; i < nSpec_; i++)
    StreamPrint(echofile, "%-3d\t%-32s\n", i + 1, &sNamesHost(i, 0));
  DASHLINE(echofile);
  // fflush(echofile);

  /* reaction info */
  if (nReac_ > 0) {
    isRev_ = ordinal_type_1d_dual_view(do_not_init_tag("KMD::isRev"), nReac_);
    reacNrp_ = ordinal_type_1d_dual_view(do_not_init_tag("KMD::reacNrp"), nReac_);
    reacNreac_ = ordinal_type_1d_dual_view(do_not_init_tag("KMD::reacNrp"), nReac_);
    reacNprod_ = ordinal_type_1d_dual_view(do_not_init_tag("KMD::reacNrp"), nReac_);
    convExponent_ = ordinal_type_1d_dual_view(do_not_init_tag("KMD::convExponent"), nReac_);
    // reacArhenFor_ =
    //   real_type_2d_dual_view(do_not_init_tag("KMD::reacArhenFor"), nReac_, 3);
    // isDup_ = ordinal_type_1d_dual_view(do_not_init_tag("KMD::isDup"), nReac_)
  }

  auto isRevHost = isRev_.view_host();
  auto reacNrpHost = reacNrp_.view_host();
  auto reacNreacHost = reacNreac_.view_host();
  auto reacNprodHost = reacNprod_.view_host();
  /* Arrhenius parameters */
  // auto reacArhenForHost = reacArhenFor_.view_host();

  ordinal_type countReac(0);
  ordinal_type countArrheniusReac(0); // count index for arrhenius pressure parameter E.
  ordinal_type countCMAQ_H2O2Reac(0);
  ordinal_type countJPL_Reac(0);
  ordinal_type countR_JPL_ArrheniusReac(0);
  ordinal_type count_adjust_reaction(0);
  ordinal_type countR_Phot_Reac(0);
  maxSpecInReac_ = 0;
  nFallReac_ = 0; // troe type reactions
  std::map<std::string, int>::iterator it;
  // species index and stoicoeff
  std::vector<std::map<ordinal_type, real_type>> productsInfo, reactantsInfo;
  std::map<ordinal_type, real_type> reactants_sp, products_sp; //, aux_parameter_arhenius
  for (auto const &reaction : reactions) {
    // get stoichimetric coefficients
    if (reaction["reactants"]) {
      auto reactants = reaction["reactants"];
      /* no of reactants only */
      reacNreacHost(countReac) = reactants.size();
      for (auto const &reac : reactants) {
        // get species index
        it = species_indx_.find(reac.first.as<std::string>());

        if (it != species_indx_.end()) {
          reactants_sp.insert(std::pair<ordinal_type, real_type>(it->second, reac.second.as<real_type>()));
          // reacSidxHost(i, count) = it->second;
        } else {
          printf("Yaml : Error when interpreting kinetic model  !!!");
          printf("species does not exit %s in reactants of reaction No %d\n", reac.first.as<std::string>().c_str(),
                 countReac);
          exit(1);
        }
      }
    } else {
      printf("Yaml : Error when interpreting kinetic model  !!!");
      printf("Reactions No %d does not have reactants\n", countReac);
      exit(1);
    }

    if (reaction["products"]) {
      auto products = reaction["products"];
      /* no of products */
      reacNprodHost(countReac) = products.size();
      for (auto const &prod : products) {
        // get species index
        it = species_indx_.find(prod.first.as<std::string>());

        if (it != species_indx_.end()) {
          products_sp.insert(std::pair<ordinal_type, real_type>(it->second, prod.second.as<real_type>()));
        } else {
          // NOTE: it is okay if a species is not register as long the reaction where this specie is involve is only a fwd reactions.
          printf("Yaml : Warning when interpreting kinetic model  !!!");
          printf("species does not exit %s in products of reaction No %d\n", prod.first.as<std::string>().c_str(),
                 countReac);
          // exit(1);
        }
      }

    } else {
      printf("Yaml : Error when interpreting kinetic model  !!!");
      printf("Reactions No %d does not have products\n", countReac);
      exit(1);
    }

    isRevHost(countReac) = 1; // all reaction are ireversible

    reacNrpHost(countReac) = reacNreacHost(countReac) + reacNprodHost(countReac);

    // check max in reactansts
    maxSpecInReac_ = maxSpecInReac_ > reacNreacHost(countReac) ? maxSpecInReac_ : reacNreacHost(countReac);

    // check max in products
    maxSpecInReac_ = maxSpecInReac_ > reacNprodHost(countReac) ? maxSpecInReac_ : reacNprodHost(countReac);

    productsInfo.push_back(products_sp);
    reactantsInfo.push_back(reactants_sp);
    products_sp.clear();
    reactants_sp.clear();

    // get Arrhenius constants
    auto reaction_type = reaction["type"].as<std::string>();

    if (reaction_type == "ARRHENIUS") {
      countArrheniusReac++;
    } else if (reaction_type == "TROE") {
      nFallReac_++;
    } else if (reaction_type == "CMAQ_H2O2") {
      countCMAQ_H2O2Reac++;
    } else if (reaction_type == "JPL") {
      countJPL_Reac++;
    } else if (reaction_type == "R_JPL_ARRHENIUS") {
      countR_JPL_ArrheniusReac++;
    } else if (reaction_type == "PHOTO_RATE") {
      countR_Phot_Reac++;
    } else {
      printf("Error: reaction type does not exit \n ");
      printf("Reaction No: %d type %s \n", countReac, reaction_type.c_str());
      exit(1);
    }
    // adjust reaction section
    if (reaction["adjust_reaction"]){
      count_adjust_reaction += reaction["adjust_reaction"].size();
    }

    countReac++;
  }

  // twice because we only consider max (products, reactants)
  maxSpecInReac_ *= 2;

  StreamPrint(echofile, "kmod.list : Max # of species in a reaction                    : %d\n", maxSpecInReac_);
  //
  StreamPrint(echofile, "kmod.list :  # of arrhenius type reactions                    : %d\n", countArrheniusReac);
  //
  StreamPrint(echofile, "kmod.list :  # of CMAQ_H2O2 type reactions                    : %d\n", countCMAQ_H2O2Reac);
  //
  StreamPrint(echofile, "kmod.list :  # of Troe type reactions                         : %d\n", nFallReac_);

  StreamPrint(echofile, "kmod.list :  # of JPL-Troe type reactions                     : %d\n", countJPL_Reac);

  StreamPrint(echofile, "kmod.list :  # of Ratio JPL arrhenius type reactions          : %d\n", countR_JPL_ArrheniusReac);
  StreamPrint(echofile, "kmod.list :  # of adjustment pairs                            : %d\n", count_adjust_reaction);
  StreamPrint(echofile, "kmod.list :  # of Photolysis rates                            : %d\n", countR_Phot_Reac);

  //
  if (nReac_ > 0) {
    reacSidx_ = ordinal_type_2d_dual_view(do_not_init_tag("KMD::reacNrp"), nReac_, maxSpecInReac_);
    reacNuki_ = real_type_2d_dual_view(do_not_init_tag("KMD::reacNrp"), nReac_, maxSpecInReac_);
    // reacScoef_ =
    // ordinal_type_1d_dual_view(do_not_init_tag("KMD::reacScoef"), nReac_);
  }

  auto reacNukiHost = reacNuki_.view_host();
  auto reacSidxHost = reacSidx_.view_host();
  // auto reacScoefHost = reacScoef_.view_host();

  /* Stoichiometric coefficients */
  for (int i = 0; i < nReac_; i++) {
    /* by default reaction has integer stoichiometric coefficients */
    // reacScoefHost(i) = -1;

    int count(0);
    auto reactants_sp = reactantsInfo[i];

    for (auto &reac : reactants_sp) {
      // set species index
      reacSidxHost(i, count) = reac.first;
      // set stoichiometric coefficient
      reacNukiHost(i, count) = -reac.second;
      count++;
    }

    count = maxSpecInReac_ / 2;
    auto products_sp = productsInfo[i];

    for (auto &prod : products_sp) {
      // set species Index
      reacSidxHost(i, count) = prod.first;
      reacNukiHost(i, count) = prod.second;
      count++;
    }
  }

  if (countArrheniusReac > 0) {
    ArrheniusCoef_ =
        arrhenius_reaction_type_1d_dual_view(do_not_init_tag("KMD::auxParamReacArhenFor"), countArrheniusReac);
  }

  if (countCMAQ_H2O2Reac > 0) {
    CMAQ_H2O2Coef_ = cmaq_h2o2_type_1d_dual_view(do_not_init_tag("KMD::coefCMAQ_H2O2ReacFor"), countCMAQ_H2O2Reac);
  }

  auto ArrheniusCoefHost = ArrheniusCoef_.view_host();
  auto CMAQ_H2O2CoefHost = CMAQ_H2O2Coef_.view_host();

  if (countJPL_Reac > 0) {
    JPL_Coef_ = jpl_reaction_type_1d_dual_view(do_not_init_tag("KMD::coefJPLReacFor"), countJPL_Reac);
  }

  auto JPL_CoefHost = JPL_Coef_.view_host();

  if (countR_JPL_ArrheniusReac > 0) {
    R_JPL_ArrheniusCoef_ = r_jpl_arrhenius_reaction_type_1d_dual_view(do_not_init_tag("KMD::coefRJPLArrheniusReacFor"), countR_JPL_ArrheniusReac);
  }

  if (count_adjust_reaction > 0) {
    adjust_reaction_ = adjust_reaction_type_1d_dual_view(do_not_init_tag("KMD::adjustReaction"), count_adjust_reaction);
  }
  auto adjust_reaction_Host = adjust_reaction_.view_host();

  auto R_JPL_ArrheniusCoefHost = R_JPL_ArrheniusCoef_.view_host();

  //// k0_A, k0_B, k0_C kinf_A,  kinf_B,  kinf_C, Fc, N
  ordinal_type number_of_param_troe(8);

  // troe type reactions
  if (nFallReac_ > 0) {
    // reactions index
    reacPfal_ = ordinal_type_1d_dual_view(do_not_init_tag("KMD::reacPfal"), nFallReac_);
    reacPpar_ = real_type_2d_dual_view(do_not_init_tag("KMD::reacPpar"), nFallReac_, number_of_param_troe);
  }

  auto reacPfalHost = reacPfal_.view_host();
  auto reacPparHost = reacPpar_.view_host();

  ordinal_type count_troe(0);
  ordinal_type count_arrhen_aux_param(0);
  ordinal_type count_cmaq_h2o2(0);
  ordinal_type ireac(0);
  ordinal_type count_jpl(0);
  ordinal_type count_r_jpl_arrhen(0);
  ordinal_type icount_adjust_reaction(0);


  for (auto const &reaction : reactions) {
    auto reaction_type = reaction["type"].as<std::string>();

    if (reaction_type == "ARRHENIUS") {
      auto rate_coefficients = reaction["coefficients"];

      // pre_exponential
      real_type A_coef(1.0);
      if (rate_coefficients["A"]) {
        A_coef = rate_coefficients["A"].as<real_type>();
        if (reaction["time_unit"]) {
          if (reaction["time_unit"].as<std::string>() == "min") {
            A_coef /= real_type(60.0);
          }
        }
      }
      // temperature coefficient
      real_type D_coef(300.0);
      if (rate_coefficients["D"]) {
        D_coef = rate_coefficients["D"].as<real_type>();
      }
      // temperature coefficient
      real_type B_coef(0);
      if (rate_coefficients["B"]) {
        B_coef = rate_coefficients["B"].as<real_type>();
      }
      // activation energy
      real_type Ea_coef(0);
      if (rate_coefficients["Ea"]) {
        // Boltzmann's constant (k_B) [JK^{-1}]$
        Ea_coef = -rate_coefficients["Ea"].as<real_type>() / KBOLT;
        if (rate_coefficients["C"]) {
          printf("Ea and C are presented in Yaml input file\n");
          exit(1);
        }
      } else if (rate_coefficients["C"]) {
        Ea_coef = rate_coefficients["C"].as<real_type>();
      }

      real_type E_coef(0);
      if (rate_coefficients["E"]) {
        E_coef = rate_coefficients["E"].as<real_type>();
      }

      ArrheniusReactionType arrhenius_reaction_type;

      // pre_exponential A_tchem = A/D^B;
      arrhenius_reaction_type._A = A_coef;
      arrhenius_reaction_type._B = B_coef;
      arrhenius_reaction_type._C = Ea_coef;
      arrhenius_reaction_type._D = D_coef;
      arrhenius_reaction_type._E = E_coef;
      arrhenius_reaction_type._reaction_index = ireac;
      ArrheniusCoefHost(count_arrhen_aux_param) = arrhenius_reaction_type;
      count_arrhen_aux_param++;
      // printf("Kforward A_ %e B_ %e C_ %e D_ %e E_ %e \n",A_coef, B_coef,  Ea_coef,  D_coef, E_coef );
    }

    if (reaction_type == "CMAQ_H2O2") {
      auto rate_coefficients = reaction["coefficients"];

      real_type k1_A(1.0);
      if (rate_coefficients["k1_A"]) {
        k1_A = rate_coefficients["k1_A"].as<real_type>();
        if (reaction["time_unit"]) {
          if (reaction["time_unit"].as<std::string>() == "min") {
            k1_A /= real_type(60.0);
          }
        }
      }

      real_type k1_B(0.0);
      if (rate_coefficients["k1_B"]) {
        k1_B = rate_coefficients["k1_B"].as<real_type>();
      }

      real_type k1_C(0.0);
      if (rate_coefficients["k1_C"]) {
        k1_C = rate_coefficients["k1_C"].as<real_type>();
      }

      real_type k2_A(1.0);
      if (rate_coefficients["k2_A"]) {
        k2_A = rate_coefficients["k2_A"].as<real_type>();
        if (reaction["time_unit"]) {
          if (reaction["time_unit"].as<std::string>() == "min") {
            k2_A /= real_type(60.0);
          }
        }
      }

      real_type k2_B(0.0);
      if (rate_coefficients["k2_B"]) {
        k2_B = rate_coefficients["k2_B"].as<real_type>();
      }

      real_type k2_C(0.0);
      if (rate_coefficients["k2_C"]) {
        k2_C = rate_coefficients["k2_C"].as<real_type>();
      }

      CMAQ_H2O2ReactionType cqmaq_h2o2_reaction;

      cqmaq_h2o2_reaction._A1 = k1_A;
      cqmaq_h2o2_reaction._B1 = k1_B;
      cqmaq_h2o2_reaction._C1 = k1_C;
      // [M]*k2
      cqmaq_h2o2_reaction._A2 = k2_A * real_type(1e6);
      cqmaq_h2o2_reaction._B2 = k2_B;
      cqmaq_h2o2_reaction._C2 = k2_C;
      cqmaq_h2o2_reaction._reaction_index = ireac;
      CMAQ_H2O2CoefHost(count_cmaq_h2o2) = cqmaq_h2o2_reaction;
      count_cmaq_h2o2++;

    } // end "CMAQ_H2O2"

    // parser of JPL, which is the troe implementation in e3sm code.

    if (reaction_type == "JPL") {
      auto rate_coefficients = reaction["coefficients"];

      real_type k0_A(1.0);
      if (rate_coefficients["k0_A"]) {
        k0_A = rate_coefficients["k0_A"].as<real_type>();
        if (reaction["time_unit"]) {
          if (reaction["time_unit"].as<std::string>() == "min") {
            k0_A /= real_type(60.0);
          }
        }
      }

      // temperature coefficient
      real_type k0_B(0);
      if (rate_coefficients["k0_B"]) {
        k0_B = rate_coefficients["k0_B"].as<real_type>();
      }

      // activation energy
      real_type k0_C(0);
      if (rate_coefficients["k0_C"]) {
        // Boltzmann's constant (k_B) [JK^{-1}]$
        k0_C = rate_coefficients["k0_C"].as<real_type>();
      }

      // kinf_A  kinf_B  kinf_C Fc N
      real_type kinf_A(1);
      if (rate_coefficients["kinf_A"]) {
        kinf_A = rate_coefficients["kinf_A"].as<real_type>();
        if (reaction["time_unit"]) {
          if (reaction["time_unit"].as<std::string>() == "min") {
            kinf_A /= real_type(60.0);
          }
        }
      }

      real_type kinf_B(0);
      if (rate_coefficients["kinf_B"]) {
        kinf_B = rate_coefficients["kinf_B"].as<real_type>();
      }

      // activation energy
      real_type kinf_C(0);
      if (rate_coefficients["kinf_C"]) {
        // Boltzmann's constant (k_B) [JK^{-1}]$
        kinf_C = rate_coefficients["kinf_C"].as<real_type>();
      }



      real_type Fc(0.6);
      if (rate_coefficients["Fc"]) {
        Fc = rate_coefficients["Fc"].as<real_type>();
      }

      JPL_ReactionType jpl_reaction;

      jpl_reaction._k0_A = k0_A;
      jpl_reaction._k0_B = k0_B;
      jpl_reaction._k0_C = k0_C;
      jpl_reaction._kinf_A =kinf_A ;
      jpl_reaction._kinf_B = kinf_B;
      jpl_reaction._kinf_C = kinf_C;
      jpl_reaction._Fc = Fc;
      jpl_reaction._reaction_index = ireac;
      JPL_CoefHost(count_jpl) = jpl_reaction;
      count_jpl++;

    } // end "JPL"

    // ratio JPL ARRHENIUS

    if (reaction_type == "R_JPL_ARRHENIUS") {

      // JPL reaction
      auto rate_coefficients = reaction["coefficients"];

      real_type k0_A(1.0);
      if (rate_coefficients["k0_A"]) {
        k0_A = rate_coefficients["k0_A"].as<real_type>();
        if (reaction["time_unit"]) {
          if (reaction["time_unit"].as<std::string>() == "min") {
            k0_A /= real_type(60.0);
          }
        }
      }

      // temperature coefficient
      real_type k0_B(0);
      if (rate_coefficients["k0_B"]) {
        k0_B = rate_coefficients["k0_B"].as<real_type>();
      }

      // activation energy
      real_type k0_C(0);
      if (rate_coefficients["k0_C"]) {
        // Boltzmann's constant (k_B) [JK^{-1}]$
        k0_C = rate_coefficients["k0_C"].as<real_type>();
      }

      // kinf_A  kinf_B  kinf_C Fc N
      real_type kinf_A(1);
      if (rate_coefficients["kinf_A"]) {
        kinf_A = rate_coefficients["kinf_A"].as<real_type>();
        if (reaction["time_unit"]) {
          if (reaction["time_unit"].as<std::string>() == "min") {
            kinf_A /= real_type(60.0);
          }
        }
      }

      real_type kinf_B(0);
      if (rate_coefficients["kinf_B"]) {
        kinf_B = rate_coefficients["kinf_B"].as<real_type>();
      }

      // activation energy
      real_type kinf_C(0);
      if (rate_coefficients["kinf_C"]) {
        // Boltzmann's constant (k_B) [JK^{-1}]$
        kinf_C = rate_coefficients["kinf_C"].as<real_type>();
      }

      real_type Fc(0.6);
      if (rate_coefficients["Fc"]) {
        Fc = rate_coefficients["Fc"].as<real_type>();
      }

      //ARRHENIUS reaction

      real_type A(1);
      if (rate_coefficients["A"]) {
        A = rate_coefficients["A"].as<real_type>();
        if (reaction["time_unit"]) {
          if (reaction["time_unit"].as<std::string>() == "min") {
            A /= real_type(60.0);
          }
        }
      }

      real_type B(0);
      if (rate_coefficients["B"]) {
        B = rate_coefficients["B"].as<real_type>();
      }

      // activation energy
      real_type C(0);
      if (rate_coefficients["C"]) {
        // Boltzmann's constant (k_B) [JK^{-1}]$
        C = rate_coefficients["C"].as<real_type>();
      }


      R_JPL_ArrheniusReactionType r_jpl_arrhenius_reaction;

      r_jpl_arrhenius_reaction._k0_A = k0_A;
      r_jpl_arrhenius_reaction._k0_B = k0_B;
      r_jpl_arrhenius_reaction._k0_C = k0_C;
      r_jpl_arrhenius_reaction._kinf_A =kinf_A ;
      r_jpl_arrhenius_reaction._kinf_B = kinf_B;
      r_jpl_arrhenius_reaction._kinf_C = kinf_C;
      r_jpl_arrhenius_reaction._Fc = Fc;
      r_jpl_arrhenius_reaction._reaction_index = ireac;

      r_jpl_arrhenius_reaction._A =A ;
      r_jpl_arrhenius_reaction._B = B;
      r_jpl_arrhenius_reaction._C = C;


      R_JPL_ArrheniusCoefHost(count_r_jpl_arrhen) = r_jpl_arrhenius_reaction;
      count_r_jpl_arrhen++;

    } // end  R JPL ARRHENIUS


    if (reaction_type == "TROE") {
      auto rate_coefficients = reaction["coefficients"];

      // k0_A, k0_B, k0_C
      // pre_exponential
      real_type A_coef(1.0);
      if (rate_coefficients["k0_A"]) {
        A_coef = rate_coefficients["k0_A"].as<real_type>();
        if (reaction["time_unit"]) {
          if (reaction["time_unit"].as<std::string>() == "min") {
            A_coef /= real_type(60.0);
          }
        }
      }

      // temperature coefficient
      real_type B_coef(0);
      if (rate_coefficients["k0_B"]) {
        B_coef = rate_coefficients["k0_B"].as<real_type>();
      }
      // activation energy
      real_type Ea_coef(0);
      if (rate_coefficients["k0_C"]) {
        // Boltzmann's constant (k_B) [JK^{-1}]$
        Ea_coef = rate_coefficients["k0_C"].as<real_type>();
      }
      // kinf_A  kinf_B  kinf_C Fc N
      real_type kinf_A(1);
      if (rate_coefficients["kinf_A"]) {
        kinf_A = rate_coefficients["kinf_A"].as<real_type>();
        if (reaction["time_unit"]) {
          if (reaction["time_unit"].as<std::string>() == "min") {
            kinf_A /= real_type(60.0);
          }
        }
      }

      real_type kinf_B(0);
      if (rate_coefficients["kinf_B"]) {
        kinf_B = rate_coefficients["kinf_B"].as<real_type>();
      }

      real_type kinf_C(0);
      if (rate_coefficients["kinf_C"]) {
        kinf_C = rate_coefficients["kinf_C"].as<real_type>();
      }

      real_type Fc(0.6);
      if (rate_coefficients["Fc"]) {
        Fc = rate_coefficients["Fc"].as<real_type>();
      }

      real_type N(1.0);
      if (rate_coefficients["N"]) {
        N = rate_coefficients["N"].as<real_type>();
      }
      // reaction index
      reacPfalHost(count_troe) = ireac;
      // auxiliary parameters
      // From camp ! Include [M] in K0_A_
      // K0_A_ = K0_A_ * real(1.0d6, kind=dp)
      reacPparHost(count_troe, 0) = A_coef * real_type(1e6);
      reacPparHost(count_troe, 1) = B_coef;
      reacPparHost(count_troe, 2) = Ea_coef;
      //
      reacPparHost(count_troe, 3) = kinf_A;
      reacPparHost(count_troe, 4) = kinf_B;
      reacPparHost(count_troe, 5) = kinf_C;
      reacPparHost(count_troe, 6) = Fc;
      reacPparHost(count_troe, 7) = N;
      count_troe++;

    } // end troe type

    // adjust reaction
    if (reaction["adjust_reaction"]){

      for (auto const &iadjust_reaction : reaction["adjust_reaction"]) {
        auto sp_name = iadjust_reaction.as<std::string>();
        it = species_indx_.find(sp_name);
        if (it != species_indx_.end()) {
          adjust_reactionType adjust_reaction;
          adjust_reaction._species_index = it->second;
          adjust_reaction._reaction_index= ireac;
          // printf("ireac %d, name %s , index %d \n ",  ireac, sp_name.c_str(), it->second );
          adjust_reaction_Host(icount_adjust_reaction) = adjust_reaction;
          icount_adjust_reaction++;
        } else {
          printf("Yaml : Error when interpreting kinetic model  !!!");
          printf("species does not exit %s in adjust_reaction for reaction No %d\n", sp_name.c_str(), ireac);
          exit(1);
        }
      } //

    }  // adjust_reaction


    ireac++;
  } // end reactions

  if (root["modifier_prod_O1D"]) {

    // NOTE: only one item in this type of modifier.
    number_of_prod_O1D_=1;
    prod_O1D_ = prod_O1D_type_1d_dual_view(do_not_init_tag("KMD::prod_O1D_"), 1);


    auto modifier_prod_O1D = root["modifier_prod_O1D"];
    auto prod_O1D_host = prod_O1D_.view_host();


    auto rate_coefficients = modifier_prod_O1D["coefficients"];
    real_type A1(1);
      if (rate_coefficients["A1"]) {
        A1 = rate_coefficients["A1"].as<real_type>();
      }

      // activation energy
      real_type C1(0);
      if (rate_coefficients["C1"]) {
        C1 = rate_coefficients["C1"].as<real_type>();
      }

      real_type A2(1);
      if (rate_coefficients["A2"]) {
        A2 = rate_coefficients["A2"].as<real_type>();
      }

      // activation energy
      real_type C2(0);
      if (rate_coefficients["C2"]) {
        C2 = rate_coefficients["C2"].as<real_type>();
      }

      real_type A3(1);
      if (rate_coefficients["A3"]) {
        A3 = rate_coefficients["A3"].as<real_type>();
      }

      // activation energy
      real_type C3(0);
      if (rate_coefficients["C3"]) {
        C3 = rate_coefficients["C3"].as<real_type>();
      }

      prod_O1DType prod_O1D;
      prod_O1D._A1 = A1;
      prod_O1D._C1 = C1;
      prod_O1D._A2 = A2;
      prod_O1D._C2 = C2;
      prod_O1D._A3 = A3;
      prod_O1D._C3 = C3;

      auto sp_name1_y = modifier_prod_O1D["species_name_1"];

      auto sp_name1 = sp_name1_y.as<std::string>();
      it = species_indx_.find(sp_name1);
      if (it != species_indx_.end()) {
          prod_O1D._species_index_1 = it->second;
          // printf("In modifier_prod_O1D; name %s , index %d \n ", sp_name1.c_str(), it->second );
      } else {
          printf("Yaml : Error when interpreting kinetic model  !!!");
          printf("species does not exit %s in modifier_prod_O1D\n", sp_name1.c_str());
          exit(1);
      }

      auto sp_name2_y = modifier_prod_O1D["species_name_2"];

      auto sp_name2 = sp_name2_y.as<std::string>();
      it = species_indx_.find(sp_name2);
      if (it != species_indx_.end()) {
          prod_O1D._species_index_2 = it->second;
          // printf("In modifier_prod_O1D; name %s , index %d \n ", sp_name2.c_str(), it->second );
      } else {
          printf("Yaml : Error when interpreting kinetic model  !!!");
          printf("species does not exit %s in modifier_prod_O1D\n", sp_name2.c_str());
          exit(1);
      }

      auto sp_name3_y = modifier_prod_O1D["species_name_3"];

      auto sp_name3 = sp_name3_y.as<std::string>();
      it = species_indx_.find(sp_name3);
      if (it != species_indx_.end()) {
          prod_O1D._species_index_3 = it->second;
          // printf("In modifier_prod_O1D; name %s , index %d \n ", sp_name3.c_str(), it->second );
      } else {
          printf("Yaml : Error when interpreting kinetic model  !!!");
          printf("species does not exit %s in modifier_prod_O1D\n", sp_name3.c_str());
          exit(1);
      }
      // FIXME: pass a reaction identifier, and then get the reaction index.
      prod_O1D._photolysis_reaction_index = modifier_prod_O1D["photolysis_reaction_index"].as<ordinal_type>() ;
      int count(0);
      for (auto const &ireac : modifier_prod_O1D["reaction_list"]) {
        prod_O1D._reaction_indices[count]= ireac.as<ordinal_type>();
        count++;
      }

      prod_O1D_host(0) = prod_O1D;

      prod_O1D_.modify_host();
      prod_O1D_.sync_device();

  }  // prod_O1D

  ordinal_type countEmissionSources(0);
  for (auto const &source : root["sources"]) {
    if (source["type"].as<std::string>() == "EMISSION") {
      countEmissionSources++;
    } else {
      printf("Error: source type does not exit \n ");
      printf("Source No: %d type %s \n", countEmissionSources, source["type"].as<std::string>().c_str());
      exit(1);
    }
  } // end source terms

  if (countEmissionSources > 0) {
    EmissionCoef_ = emission_source_type_1d_dual_view(do_not_init_tag("KMD::EmissionCoefHost_"), countEmissionSources);
  }

  auto EmissionCoefHost = EmissionCoef_.view_host();
  ordinal_type count_emission_coef(0);
  for (auto const &source : root["sources"]) {
    if (source["type"].as<std::string>() == "EMISSION") {
      EMISSION_SourceType emission_source;
      if (source["coefficients"]["emission_rate"]) {
        emission_source._emissition_rate = source["coefficients"]["emission_rate"].as<real_type>();

      } else {
        printf("Error EMISSION, emission_rate is not provided \n");
        exit(1);
      }

      if (source["species"]) {
        // get species index
        it = species_indx_.find(source["species"].as<std::string>());
        if (it != species_indx_.end()) {
          emission_source._species_index = it->second;
        } else {
          printf("Yaml : Error when interpreting kinetic model  !!!");
          printf("species does not exit %s\n", source["species"].as<std::string>().c_str());
          exit(1);
        }

      } else {
        printf("Error EMISSION, species name  is not provided \n");
        exit(1);
      }
      EmissionCoefHost(count_emission_coef) = emission_source;
      count_emission_coef++;
    } // end source type emission
  }   // end sources

  auto convExponentHost = convExponent_.view_host();

  for (ordinal_type i = 0; i < nReac_; i++) {
    ordinal_type nu_sum = 0;
    for (ordinal_type j = 0; j < reacNreacHost(i); ++j) {
      nu_sum += ats<ordinal_type>::abs(reacNukiHost(i, j));
    }
    convExponentHost(i) = nu_sum - ordinal_type(1);
  }

  StreamPrint(echofile, "kmod.list :  # sources e.g., EMISSION                         : %d\n", countEmissionSources);

  StreamPrint(echofile, "Reaction data : species and Arrhenius pars\n");
  for (int i = 0; i < nReac_; i++) {
    StreamPrint(echofile, "%-5d\t%1d\t%2d\t%2d | ", i + 1, isRevHost(i), reacNreacHost(i), reacNprodHost(i));
    //
    for (int j = 0; j < reacNreacHost(i); j++)
      StreamPrint(echofile, "%f*%s | ", reacNukiHost(i, j), &sNamesHost(reacSidxHost(i, j), 0));

    /// KJ why do we do this way ?
    const int joff = maxSpecInReac_ / 2;
    for (int j = 0; j < reacNprodHost(i); j++)
      StreamPrint(echofile, "%f*%s | ", reacNukiHost(i, j + joff), &sNamesHost(reacSidxHost(i, j + joff), 0));

    // StreamPrint(echofile,
    //         "%16.8e\t%16.8e\t%16.8e",
    //         reacArhenForHost(i, 0),
    //         reacArhenForHost(i, 1),
    //         reacArhenForHost(i, 2));

    StreamPrint(echofile, "\n");

    if (verboseEnabled)
      printf("KineticModelData::initChem() : Done reading reaction data\n");
  }

  sNames_.modify_host();
  isRev_.modify_host();
  reacNrp_.modify_host();
  reacNreac_.modify_host();
  reacNprod_.modify_host();
  reacNuki_.modify_host();
  reacSidx_.modify_host();
  // reacScoef_.modify_host();
  // reacArhenFor_.modify_host();
  ArrheniusCoef_.modify_host();
  CMAQ_H2O2Coef_.modify_host();
  EmissionCoef_.modify_host();
  JPL_Coef_.modify_host();
  R_JPL_ArrheniusCoef_.modify_host();

  reacPfal_.modify_host();
  reacPpar_.modify_host();

  convExponent_.modify_host();
  adjust_reaction_.modify_host();

  /* Species' name and weights */
  sNames_.sync_device();
  isRev_.sync_device();
  reacNrp_.sync_device();
  reacNreac_.sync_device();
  reacNprod_.sync_device();
  reacNuki_.sync_device();
  reacSidx_.sync_device();
  // reacArhenFor_.sync_device();
  ArrheniusCoef_.sync_device();
  CMAQ_H2O2Coef_.sync_device();
  JPL_Coef_.sync_device();
  R_JPL_ArrheniusCoef_.sync_device();
  EmissionCoef_.sync_device();
  reacPfal_.sync_device();
  reacPpar_.sync_device();
  convExponent_.sync_device();
  adjust_reaction_.sync_device();


  return (0);
}
kmd_type_1d_view_host KineticModelData::clone(const int n_models) {
  kmd_type_1d_view_host r_val(do_not_init_tag("KMD::cloned models"), n_models);
  Kokkos::deep_copy(r_val, *this);
  return r_val;
}

#endif

} // namespace TChem
