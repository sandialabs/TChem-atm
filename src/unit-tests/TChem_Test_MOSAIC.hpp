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
#ifndef __TCHEM_TEST_MOSAIC_HPP__
#define __TCHEM_TEST_MOSAIC_HPP__

#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"

namespace TChem {
namespace Test {
void static MOSAIC_Test() {

    using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
    const auto host_exec_space = TChem::host_exec_space();
    const auto mmd = TChem::Impl::MosaicModelData<device_type>();
    auto b_mtem = mmd.b_mtem.template view<device_type>();
    const auto b_mtem_h = Kokkos::create_mirror_view_and_copy(host_exec_space, b_mtem);
    EXPECT_TRUE(b_mtem_h(0,mmd.jnh4so4,mmd.jnh4hso4) == -4.13219);
}

}// namespace Test
}// namespace TChem

TEST(MOSAIC, device) {
  TChem::Test::MOSAIC_Test();
}

#endif
