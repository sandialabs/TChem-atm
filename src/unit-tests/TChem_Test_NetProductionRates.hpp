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
#ifndef __TCHEM_TEST_NET_PRODUCTION_RATES_HPP__
#define __TCHEM_TEST_NET_PRODUCTION_RATES_HPP__

#include "TChem_KineticModelData.hpp"
#include "TChem_NetProductionRates.hpp"

TEST(NetProductionRates, single)
{
  std::string exec="../examples/TChem_NetProductionRates.x";

  std::string prefixPath="../examples/runs/atmospheric_chemistry/uci_e3sm_v3/";
  std::string chemFile(prefixPath + "uci_v2.yaml");
  std::string outputFile("net-production-rates.dat");
  std::string invoke=(exec+" "+
		      "--chemfile="+chemFile+" "+
		      "--outputfile="+outputFile+" "+		      
		      "--batch_size=10 --use_cloned_samples=true verbose=true");
  const auto invoke_c_str = invoke.c_str();
  printf("testing : %s\n", invoke_c_str);
  std::system(invoke_c_str);

  /// compare with ref
  EXPECT_TRUE(TChem::Test::compareFiles("net-production-rates.dat",
					"references/net_production_rates/net-production-rates.dat")
	      );

  EXPECT_TRUE(TChem::Test::compareFiles("kmod.echo",
          "references/net_production_rates/kmod.echo")
        );
}


TEST(NetProductionRates_uci, single)
{
  std::string exec="../examples/TChem_NetProductionRates.x";

  std::string prefixPath="../examples/runs/atmospheric_chemistry/uci_col/";
  std::string chemFile(prefixPath + "uci_v2_test3.yaml");
  std::string inputFile(prefixPath + "input_conditions_col.yaml");
  std::string outputFile("net-production-rates.dat");
  std::string invoke=(exec+" "+
          "--chemfile="+chemFile+" "+
          "--inputfile="+inputFile+" "+
          "--outputfile="+outputFile+" "+         
          "--verbose=true");
  const auto invoke_c_str = invoke.c_str();
  printf("testing : %s\n", invoke_c_str);
  std::system(invoke_c_str);

  /// compare with ref
  EXPECT_TRUE(TChem::Test::compareFiles("net-production-rates.dat",
          "references/net_production_rates/uci_col/net_production_rates.dat")
        );

  EXPECT_TRUE(TChem::Test::compareFiles("kmod.echo",
          "references/net_production_rates/uci_col/kmod.echo")
        );
}

#endif
