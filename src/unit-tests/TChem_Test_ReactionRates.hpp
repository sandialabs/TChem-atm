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
#ifndef __TCHEM_TEST_REACTIONRATES_HPP__
#define __TCHEM_TEST_REACTIONRATES_HPP__

#include "TChem_KineticModelData.hpp"
#include "TChem_ReactionRates.hpp"

TEST(ReactionRates, single)
{
  std::string exec="../examples/TChem_ReactionRates.x";

  std::string prefixPath="../examples/runs/atmospheric_chemistry/uci_e3sm_v3/";
  std::string chemFile(prefixPath + "uci.yaml");
  std::string outputFile("reaction-rates.dat");
  std::string invoke=(exec+" "+
		      "--chemfile="+chemFile+" "+
		      "--outputfile="+outputFile+" "+		      
		      "--batch_size=10 --use_cloned_samples=true verbose=true");
  const auto invoke_c_str = invoke.c_str();
  printf("testing : %s\n", invoke_c_str);
  std::system(invoke_c_str);

  /// compare with ref
  EXPECT_TRUE(TChem::Test::compareFiles("reaction-rates.dat",
					"references/reaction_rates/reaction-rates.dat")
	      );

  EXPECT_TRUE(TChem::Test::compareFiles("kmod.echo",
          "references/reaction_rates/kmod.echo")
        );
}

#endif
