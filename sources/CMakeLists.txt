add_library(CryoSimLib SHARED cryostat.cxx )
add_executable(CryoSim CryoSim.cxx)
target_link_libraries(CryoSimLib ${GSL_LIBRARIES})
target_link_libraries(CryoSim CryoSimLib)
