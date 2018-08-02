cmake_minimum_required(VERSION 2.8)

project(Fastmech_iRICv3 Fortran)

enable_language(C)
include(GNUInstallDirs)

SET(SOURCES
"src/betaCDFandINV.f90"
"src/betavol.f90"
"src/elevdt.f90"
"src/functionrough.f90"
"src/initdt.f90"
"src/initHdt.f90"
"src/initimeta.f90"
"src/rough.f90"
"src/roughnessinv.f90"
"src/volfrac.f90"
"src/CalcCondMod.F90"
"src/GridCondMod.F90"
"src/GridCoordMod.F90"
"src/iriclibf.f90"
"src/ReadCGNS.F90"
"src/WriteCGNS.F90"
"src/ParmSet.f90"
"src/StbBkwOutput.f90"
"src/StpBkwInput.f90"
"src/Support.f90"
"src/xAnlzFlow.f90"
"src/xFArea2.f90"
"src/xFChPreis.f90"
"src/xFInput.f90"
"src/xFlowSupport.F90"
"src/xFMatInv.f90"
"src/xFNetInitConds.f90"
"src/xFNetPreisMain.f90"
"src/xFTSBoundCond.f90"
"src/CsedMod.f90"
"src/CsedMod_DT.f90"
"src/Einit.f90"
"src/GeometryMod.F90"
"src/Median.f90"
"src/RivCalcInitCond.F90"
"src/RivConnectivityMod.F90"
"src/RivRoughnessMod.F90"
"src/RivStagrMod_jmn.F90"
"src/RivVarMod.f90"
"src/RivVarTimeMod.f90"
"src/RivVarVertMod.f90"
"src/RivVarWMod.f90"
"src/RivVertMod2.f90"
"src/Stagr_Console.f90"
"src/stressDivMod.f90"
"src/tridag.f90"
"src/Uinit.f90"
"src/Welcome.f90"
)

# executable
add_executable(Fastmech_iRICv3 ${SOURCES})

# includes
find_path(CGNS_INC cgnslib_f.h)
include_directories("${CGNS_INC}")
find_path(IRIC_INC iriclib_f.h)
include_directories("${IRIC_INC}")

# library dependencies
if(WIN32)
  find_library(CGNS_LIBRARY cgnsdll)
  find_library(HDF5_LIBRARY hdf5)
  find_library(IRICLIB_LIBRARY iriclib)
else()
  find_library(CGNS_LIBRARY libcgns)
  find_library(HDF5_LIBRARY libhdf5)
  find_library(IRICLIB_LIBRARY libiriclib)
endif()
#find_library(HDF5_LIBRARY hdf5)
#find_library(IRICLIB_LIBRARY iriclib)
SET(EXTRA_LIBS ${EXTRA_LIBS} ${IRICLIB_LIBRARY} ${CGNS_LIBRARY} ${HDF5_LIBRARY})

# link
target_link_libraries(Fastmech_iRICv3 ${EXTRA_LIBS})

# default real kind 8
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  if(WIN32)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /real-size:64")
  else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -real-size 64")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
endif()

# set stack on windows
if(WIN32)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /HEAP:100000000,100000000 /STACK:1000000000,1000000000")
else()
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fmax-stack-var-size=100000000")
endif()

# install
install(TARGETS Fastmech_iRICv3 DESTINATION "${CMAKE_INSTALL_BINDIR}")
