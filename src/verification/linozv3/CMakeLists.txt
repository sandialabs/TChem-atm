SET(DRIVER_NAME "TChem_Linozv3_StratosphereChemSolver")
add_executable(${DRIVER_NAME}.x
               ${DRIVER_NAME}.cpp)
TARGET_LINK_LIBRARIES(${DRIVER_NAME}.x verification;${TCHEM_ATM_LINK_LIBRARIES})
