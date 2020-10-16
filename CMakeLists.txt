cmake_minimum_required (VERSION 2.6)
project (CryoSim)
# The version number.
set (CRYOSIM_VERSION_MAJOR 0)
set (CRYOSIM_VERSION_MINOR 1)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_COMPILER clang++)
cmake_policy(SET CMP0042 NEW)

if (GSL_PATH)
  set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${GSL_PATH})
  set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${GSL_PATH})
endif (GSL_PATH) 

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package (GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  add_compile_options (-fdiagnostics-color=always)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  add_compile_options (-fcolor-diagnostics)
endif ()

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set(CMAKE_INCLUDE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/include )
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )

# add lib
# set(CRYOSIM_SUB_DIRECTORIES sources test)
set(CRYOSIM_SUB_DIRECTORIES sources)

# Macro that add the sub dir and the include file install path
foreach(subdir ${CRYOSIM_SUB_DIRECTORIES})
  add_subdirectory (${subdir})  
  include_directories(${subdir}) 

  file(GLOB incfiles "${subdir}/*.h")
  foreach(incfile ${incfiles})
    install (FILES ${incfile} DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})
  endforeach(incfile)
endforeach(subdir)

