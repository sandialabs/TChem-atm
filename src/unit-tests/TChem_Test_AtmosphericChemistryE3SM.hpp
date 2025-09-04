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

Questions? Contact Oscar Diaz-Ibarra at <odiazib@sandia.gov>, or
           Cosmin Safta at <csafta@sandia.gov> or,
           Nicole Riemer at <nriemer@illinois.edu> or,
           Matthew West at <mwest@illinois.edu>

Sandia National Laboratories, New Mexico/Livermore, NM/CA, USA
===================================================================================== */
#ifndef __TCHEM_TEST_ATMOSPHERIC_CHEMISTRY_E3SM_HPP__
#define __TCHEM_TEST_ATMOSPHERIC_CHEMISTRY_E3SM_HPP__

#include "TChem_KineticModelData.hpp"
#include "TChem_AtmosphericChemistryE3SM.hpp"
// FIXME:
// Test is not passing, but we are using a old version of chem file. 
#if 0
TEST(AtmosphericChemistry, single)
{
  std::string exec="../examples/TChem_AtmosphericChemistryE3SM.x";

  std::string prefixPath="../examples/runs/atmospheric_chemistry/uci_e3sm_v3/";
  std::string chemFile(prefixPath + "uci_v2.yaml");
  std::string outputFile("atmospheric-chemistry-e3sm.dat");
  std::string invoke=(exec+" "+
		      "--chemfile="+chemFile+" "+
		      "--outputfile="+outputFile+" "+		
          "--max-time-iterations=3 " +      
		      "--batch_size=10 --use_cloned_samples=true verbose=true");
  const auto invoke_c_str = invoke.c_str();
  printf("testing : %s\n", invoke_c_str);
  std::system(invoke_c_str);

  /// pass test with a relative error of
  const real_type threshold =1e-8; 
  EXPECT_TRUE(TChem::Test::compareFilesValues("atmospheric-chemistry-e3sm.dat",
					"references/atmospheric_chemistry_e3sm/atmospheric-chemistry-e3sm.dat",threshold)
	      );

  EXPECT_TRUE(TChem::Test::compareFiles("kmod.echo",
          "references/atmospheric_chemistry_e3sm/kmod.echo")
        );
}
#endif
TEST(AtmosphericChemistry_uci, single)
{
  std::string exec="../examples/TChem_AtmosphericChemistryE3SM.x";

  std::string prefixPath="../examples/runs/atmospheric_chemistry/uci_col/";
  std::string chemFile(prefixPath + "uci_v2_test3.yaml");
  std::string inputFile(prefixPath + "input_conditions_col.yaml");
  std::string outputFile("full_gas.dat");

  std::string invoke=(exec+" "+
          "--chemfile="+chemFile+" "+
          "--outputfile="+outputFile+" "+  
          "--inputfile="+inputFile+" "+ 
          "--time-iterations-per-interval=100 " +  
          "--tol-time=1e-6 " +
          "--dtmin=1e-20 " +
          "--dtmax=3600 " +
          "--tend=1800 " +
          "--atol-newton=1e-18 "+
          "--rtol-newton=1e-8 " +
          "--max-newton-iterations=20 "+
          "--max-time-iterations=20000 " +
          "--verbose=true");
  const auto invoke_c_str = invoke.c_str();
  printf("testing : %s\n", invoke_c_str);
  std::system(invoke_c_str);

  /// pass test with a relative error of
  const real_type threshold =1e-8; 
  EXPECT_TRUE(TChem::Test::compareFilesValues("full_gas.dat",
          "references/atmospheric_chemistry_e3sm/uci_col/full_gas.dat",threshold)
        );

  EXPECT_TRUE(TChem::Test::compareFiles("kmod.echo",
          "references/atmospheric_chemistry_e3sm/uci_col/kmod.echo")
        );
}


#endif
