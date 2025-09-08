SET(TCHEM_VERIFICATION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/data_sets)

ADD_LIBRARY(verification verification.cpp)
target_link_libraries(verification ${TCHEM_ATM_LINK_LIBRARIES})

TARGET_INCLUDE_DIRECTORIES(verification
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
# Skywalker drivers
#--------------------------
add_subdirectory(mosaic)
add_subdirectory(linozv3)
