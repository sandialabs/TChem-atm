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
#ifndef __TCHEM_TEST_RATE_OF_PROGRESS_HPP__
#define __TCHEM_TEST_RATE_OF_PROGRESS_HPP__

#include "TChem_KineticModelData.hpp"
#include "TChem_RateofProgress.hpp"

TEST(RateofProgress, single)
{
  std::string exec="../examples/TChem_RateofProgress.x";

  std::string prefixPath="../examples/runs/atmospheric_chemistry/uci_e3sm_v3/";
  std::string chemFile(prefixPath + "uci_v2.yaml");
  std::string outputFile("rate-of-progress.dat");
  std::string invoke=(exec+" "+
		      "--chemfile="+chemFile+" "+
		      "--outputfile="+outputFile+" "+		      
		      "--batch_size=10 --use_cloned_samples=true verbose=true");
  const auto invoke_c_str = invoke.c_str();
  printf("testing : %s\n", invoke_c_str);
  std::system(invoke_c_str);

  /// compare with ref
  EXPECT_TRUE(TChem::Test::compareFiles("rate-of-progress.dat",
					"references/rate_of_progress/rate-of-progress.dat")
	      );

  EXPECT_TRUE(TChem::Test::compareFiles("kmod.echo",
          "references/rate_of_progress/kmod.echo")
        );
}

#endif
