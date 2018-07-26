cmake_minimum_required(VERSION 2.8)

project(Fastmech_iRICv3 Fortran)

enable_language(C)
include(GNUInstallDirs)

SET(DT_Wilcox_Kenworthy_SOURCES
"Daniele Tonina/betaCDFandINV.f90"
"Daniele Tonina/betavol.f90"
"Daniele Tonina/elevdt.f90"
"Daniele Tonina/functionrough.f90"
"Daniele Tonina/initdt.f90"
"Daniele Tonina/initHdt.f90"
"Daniele Tonina/initimeta.f90"
"Daniele Tonina/rough.f90"
"Daniele Tonina/roughnessinv.f90"
"Daniele Tonina/volfrac.f90"
)


SET(iRIC_SOURCES
CalcCondMod.F90
GridCondMod.F90
GridCoordMod.F90
iriclibf.f90
ReadCGNS.F90
WriteCGNS.F90
)

SET(Jim1D_Code_SOURCES
"Jim1D Code/ParmSet.f90"
"Jim1D Code/StbBkwOutput.f90"
"Jim1D Code/StpBkwInput.f90"
"Jim1D Code/Support.f90"
"Jim1D Code/xAnlzFlow.f90"
"Jim1D Code/xFArea2.f90"
"Jim1D Code/xFChPreis.f90"
"Jim1D Code/xFInput.f90"
"Jim1D Code/xFlowSupport.F90"
"Jim1D Code/xFMatInv.f90"
"Jim1D Code/xFNetInitConds.f90"
"Jim1D Code/xFNetPreisMain.f90"
"Jim1D Code/xFTSBoundCond.f90"
)

SET(Fastmech_iRICv3_SOURCES
CsedMod.f90
CsedMod_DT.f90
Einit.f90
falveld.F90
GeometryMod.F90
Median.f90
RivCalcInitCond.F90
RivConnectivityMod.F90
RivRoughnessMod.F90
RivStagrMod_jmn.F90
RivVarMod.for
RivVarTimeMod.f90
RivVarVertMod.f90
RivVarWMod.f90
RivVertMod2.f90
Stagr_Console.f90
stressDivMod.f90
tridag.f90
Uinit.f90
Welcome.f90
)

# executable
add_executable(Fastmech_iRICv3 ${DT_Wilcox_Kenworthy_SOURCES} ${iRIC_SOURCES} ${Jim1D_Code_SOURCES} ${Fastmech_iRICv3_SOURCES})

# includes
find_path(CGNS_INC cgnslib_f.h)
include_directories("${CGNS_INC}")
find_path(IRIC_INC iriclib_f.h)
include_directories("${IRIC_INC}")

# library dependencies
if(WIN32)
  find_library(CGNS_LIBRARY cgnsdll)
else()
  find_library(CGNS_LIBRARY cgns)
endif()
find_library(HDF5_LIBRARY hdf5)
find_library(IRICLIB_LIBRARY iriclib)
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
endif()

# install
install(TARGETS Fastmech_iRICv3 DESTINATION "${CMAKE_INSTALL_BINDIR}")
