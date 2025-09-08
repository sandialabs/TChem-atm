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
#include "TChem_AerosolModelData.hpp"

#include <cstdarg>
#include <algorithm>

namespace TChem {


AerosolModelData::AerosolModelData(const std::string &mechfile,
                                   const KineticModelData& kmd,
 std::ostream& echofile) {
    setGasParameters(kmd);
    initFile(mechfile, echofile);
}

AerosolModelData::AerosolModelData(const std::string &mechfile,
                                   const KineticModelData& kmd) {
    std::ofstream echofile;
    echofile.open("kmod.echo");
    setGasParameters(kmd);
    initFile(mechfile ,echofile);
    echofile.close();
}

void AerosolModelData::setGasParameters(const KineticModelData& kmd){
  // FIXME: I maybe need to add a map of gas species as an input in initChem
  const ordinal_type n_active_species = kmd.nSpec_- kmd.nConstSpec_;
  nSpec_gas_=n_active_species;
  nConstSpec_gas_=kmd.nConstSpec_;//
  // maps from kmd has order of gas species.
  gas_sp_name_idx_ = kmd.species_indx_;
  is_gas_parameters_set_=true;
  printf("Setting gas parameters... nSpec_gas_ %d  nConstSpec_gas_%d \n", nSpec_gas_, nConstSpec_gas_);
 }
void AerosolModelData::initFile(const std::string &mechfile,
                                std::ostream& echofile){

  YAML::Node doc = YAML::LoadFile(mechfile);
  // FIXME: add error checking in yaml parser
  if (doc["NCAR-version"]) {
    if (verboseEnabled) {
      printf("Using TChem parser for atmospheric chemistry\n");
    }
    initChem(doc, echofile);
  } else {
    TCHEM_CHECK_ERROR(true,"we do not have a parser for this file." );
    // exit
  }
}

int AerosolModelData::initChem(YAML::Node &root,
                               std::ostream& echofile) {

    TCHEM_CHECK_ERROR(!is_gas_parameters_set_,"Error: gas parameters are not set use: setGasParameters(n_active_gas) ." );

    //gas_idx_sp_name : given index return spacies name
    std::map<int, std::string> gas_idx_sp_name;
    for (std::map<std::string, int>::iterator
         i = gas_sp_name_idx_.begin();
         i != gas_sp_name_idx_.end(); ++i)
    gas_idx_sp_name[i->second] = i->first;

    // 1. let's make a map of aerosol species and gas species.
    // given the name of species return the YAML::NODE
    std::map<std::string, YAML::Node> gas_sp_info;
    // 2. get molecular weitghs and density of aerosol_species
    std::vector<real_type> mw_aerosol_sp, density_aero_sp, kappa_aero_sp;
    int i_aero_sp=0;
    // loops over species, only make map from aerosol species.
    // we assume that map of gas species was previously created in kmd.
    for (auto const &item : root["NCAR-version"] ) {
      auto type =item["type"].as<std::string>();
      if (type=="CHEM_SPEC"){
        // std::cout <<"sp_name " << item["name"]<<"\n";
        const auto sp_name = item["name"].as<std::string>();
        if (item["phase"]){
           auto phase_name =item["phase"].as<std::string>();
          //  std::cout <<"sp_name " << sp_name<< item["phase"]<<"\n";
           if (phase_name == "AEROSOL" ){
            // std::cout <<"sp_name " << sp_name<< item["phase"]<<"\n";
             aerosol_sp_name_idx_[sp_name] = i_aero_sp;
             i_aero_sp++;
            //  std::cout <<"molecular weight" << item["molecular weight [kg mol-1]"]<<"\n";
             mw_aerosol_sp.push_back(item["molecular weight [kg mol-1]"].as<real_type>());
             density_aero_sp.push_back(item["density [kg m-3]"].as<real_type>());
             kappa_aero_sp.push_back(item["kappa"].as<real_type>());
           }
        } else
        {
          // save YAML::NODE, we need it in simpol.
          // order of species can be different in gas conf file and aero mech file.
          gas_sp_info.insert(std::pair<std::string, YAML::Node>(sp_name,item));
        }// phase
      } // if CHEM_SPEC

    // get number of particles
      if (type=="AERO_REP_SINGLE_PARTICLE"){
        nParticles_= item["maximum computational particles"].as<int>();
      }// AERO_REP_SINGLE_PARTICLE

    } // item loop
    printf("Done with species...\n");
    std::vector<SIMPOL_PhaseTransferType> simpol_info;
    for (auto const &item : root["NCAR-version"]) {
      auto type =item["type"].as<std::string>();
      if (type=="MECHANISM"){
        // loop over reactions
       auto reactions = item["reactions"];
       for (auto const &ireac : reactions) {
        // std::cout <<"reactions " << ireac<<"\n";
         auto reaction_type = ireac["type"].as<std::string>();
        //  std::cout <<"reactions " << reaction_type<<"\n";
         if (reaction_type=="SIMPOL_PHASE_TRANSFER")
         {
          SIMPOL_PhaseTransferType simpol_info_at;

          simpol_info_at.B1 = ireac["B"][0].as<real_type>();
          simpol_info_at.B2 = ireac["B"][1].as<real_type>();
          simpol_info_at.B3 = ireac["B"][2].as<real_type>();
          simpol_info_at.B4 = ireac["B"][3].as<real_type>();
          // getting aerosol and gas species index.
          const auto aero_sp = ireac["aerosol-phase species"].as<std::string>();
          // std::cout <<"aero_sp " << aero_sp<<"\n";
          auto it_aero = aerosol_sp_name_idx_.find(aero_sp);
          if (it_aero != aerosol_sp_name_idx_.end()) {
           // aerosol species are place after gas species.
           //Note: that we do not use number of particles here.
           simpol_info_at.aerosol_sp_index = it_aero->second + nSpec_gas_;
           } else {
           printf("species does not exit %s in reactants of SIMPOL_PHASE_TRANSFER reaction No \n", aero_sp.c_str());
           TCHEM_CHECK_ERROR(true,"Yaml : Error when interpreting kinetic model" );
          }

           const auto gas_sp = ireac["gas-phase species"].as<std::string>();
           auto it_gas = gas_sp_name_idx_.find(gas_sp);
           if (it_gas != gas_sp_name_idx_.end()) {
            //
            simpol_info_at.gas_sp_index = it_gas->second;
           } else {
            printf("Error : species does not exit %s in reactants of SIMPOL_PHASE_TRANSFER\n", gas_sp.c_str());
            TCHEM_CHECK_ERROR(true,"Yaml : Error when interpreting aerosol model" );
          }
          simpol_info.push_back(simpol_info_at);
         } else {
        printf("Warning : TChem (amd) did not parse this reaction type %s \n", reaction_type.c_str());
       }// ireac
      }
    }// item loop

    } // item
    printf("Done with reactions/phase trans...\n");

    // number of aerosol species
    nSpec_= density_aero_sp.size();

    nSimpol_tran_=simpol_info.size();
    printf("Number of Simpol phase transfer %d \n",nSimpol_tran_ );
    // update simpol
    simpol_params_ = simpol_phase_transfer_type_1d_dual_view(do_not_init_tag("AMD::simpol_params_"), nSimpol_tran_);
    const auto simpol_params_host = simpol_params_.view_host();

    for (ordinal_type isimpol = 0; isimpol < nSimpol_tran_; isimpol++)
    {
      // SIMPOL_PhaseTransferType simpol_params_host_at_i;
      auto simpol_params_host_at_i = simpol_info[isimpol];
      // get name of gas species
      auto sp_name = gas_idx_sp_name[simpol_params_host_at_i.gas_sp_index];
      // get YAML::NODE using sp_name gas species
      auto gas_sp_at_i = gas_sp_info[sp_name];
      // std::cout <<"gas_sp_index " << simpol_params_host_at_i.gas_sp_index<<"\n";
      // std::cout <<"sp_name " << sp_name<<"\n";
      // std::cout <<"gas_sp_at_i " << gas_sp_at_i<<"\n";
      simpol_params_host_at_i.diffusion_coeff=gas_sp_at_i["diffusion coeff [m2 s-1]"].as<real_type>(); // m2 s-;
      if (gas_sp_at_i["N star"]) {
        simpol_params_host_at_i.N_star=gas_sp_at_i["N star"].as<real_type>();
        simpol_params_host_at_i.compute_alpha=true;
      } else {
        simpol_params_host_at_i.compute_alpha=false;
        simpol_params_host_at_i.N_star=999;
      }
      simpol_params_host_at_i.molecular_weight=gas_sp_at_i["molecular weight [kg mol-1]"].as<real_type>();
      // save simpol_params_host_at_i in kokkos array
      simpol_params_host(isimpol)=simpol_params_host_at_i;
    } // nsimpol


    molecular_weights_ = real_type_1d_dual_view(do_not_init_tag("AMD::molecular_weights"), nSpec_);
    auto molecular_weights_host = molecular_weights_.view_host();
    aerosol_density_= real_type_1d_dual_view(do_not_init_tag("AMD::aerosol_density_"), nSpec_);
    auto aerosol_density_host = aerosol_density_.view_host();
    aerosol_kappa_ = real_type_1d_dual_view(do_not_init_tag("AMD::kappa"), nSpec_);
    auto aerosol_kappa_host = aerosol_kappa_.view_host();
    for (int i = 0; i < nSpec_; i++)
    {
      molecular_weights_host(i) = mw_aerosol_sp[i];
      aerosol_density_host(i) = density_aero_sp[i];
      aerosol_kappa_host(i) = kappa_aero_sp[i];
    }

    molecular_weights_.modify_host();
    aerosol_density_.modify_host();
    simpol_params_.modify_host();

    molecular_weights_.sync_device();
    aerosol_density_.sync_device();
    simpol_params_.sync_device();

    return 0;
}

void AerosolModelData::setNumberofParticles(const ordinal_type number_of_particles)
{
  printf("-------------------------------------------------------\n");
  printf("--------------------Warning----------------------------\n");
  printf("Setting number of particles\n");
  printf("Old value : %d \n", nParticles_);
  nParticles_=number_of_particles;
  printf("Current value : %d \n", nParticles_);
  printf("-------------------------------------------------------\n");
}

void AerosolModelData::scenarioConditionParticles(const std::string &mechfile,
                                                  const ordinal_type nBatch,
                                                  real_type_2d_view_host& num_concentration,
                                                  real_type_2d_view_host& state)
{
  YAML::Node root = YAML::LoadFile(mechfile);

  if (root["particles"]){

    auto num_conc_info = root["particles"]["num_concentration"]["initial_value"];
    num_concentration = real_type_2d_view_host("num_concentration",nBatch, nParticles_ );
    // loop over batch
    for (ordinal_type ibatch = 0; ibatch < nBatch; ibatch++)
    {
      // Note: we assume all particles have equal num_concentration
      const auto value = num_conc_info[ibatch].as<real_type>();
      auto num_con_ibatch = Kokkos::subview(num_concentration, ibatch, Kokkos::ALL());
      Kokkos::deep_copy(num_con_ibatch,value);
    }

    // read initial conditions for part concentraion
    auto initial_state = root["particles"]["initial_state"];
    // density, pressure, temperature, active gas species, and invariant gas species
    const ordinal_type initial_index = 3 +nSpec_gas_+nConstSpec_gas_;
    for (auto const& sp_cond : initial_state) {
      auto species_name = sp_cond.first.as<std::string>();
      auto it = aerosol_sp_name_idx_.find(species_name);
      if (it == aerosol_sp_name_idx_.end()) {
          printf("species does not exit %s in reactants of SIMPOL_PHASE_TRANSFER reaction No \n", species_name.c_str());
          TCHEM_CHECK_ERROR(true,"Yaml : Error when interpreting kinetic model" );
      }
      const ordinal_type species_idx = it->second;
      std::cout << species_idx << "\n";
      auto species_info = initial_state[species_name]["initial_value"];

      for (ordinal_type ibatch = 0; ibatch < nBatch; ibatch++)
      {
       const auto value = species_info[ibatch].as<real_type>();
       auto state_ibatch = Kokkos::subview(state, ibatch, Kokkos::ALL());
       auto num_con_ibatch = Kokkos::subview(num_concentration, ibatch, Kokkos::ALL());
       for (int ipart = 0; ipart < nParticles_; ipart++)
      {
      state_ibatch(initial_index+nSpec_*ipart+species_idx) = value/num_con_ibatch(ipart);
    }
    } // batches

  }





  }

}




} // namespace TChem
