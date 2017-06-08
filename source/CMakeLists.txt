cmake_minimum_required(VERSION 2.6)

# Define project
PROJECT (predictive_STS)

SET(BUILD_BINARY_DIR ${CMAKE_BINARY_DIR}
    CACHE PATH 
    "The Predictive_STS build (not the install) puts all the libraries and executables together here (with /Release, etc. appended on some platforms).")

# Make everything go in the same binary directory. (These are CMake-defined
# variables.)
SET(EXECUTABLE_OUTPUT_PATH ${BUILD_BINARY_DIR})
SET(LIBRARY_OUTPUT_PATH ${BUILD_BINARY_DIR})

INCLUDE_DIRECTORIES(${OpenSim_SOURCE_DIR} ${OpenSim_SOURCE_DIR}/Vendors Common)

SET(OPENSIM_INSTALL_DIR $ENV{OPENSIM_HOME} 
		CACHE PATH "Top-level directory of OpenSim install")

# To add Debug feature add ";Debug" after Release on the line below
SET(CMAKE_CONFIGURATION_TYPES "RelWithDebInfo;Release;Debug" 
      CACHE STRING "Semicolon separated list of supported configuration types, only supports Debug, Release, MinSizeRel, and RelWithDebInfo, anything else will be ignored." FORCE )
	
# Location of headers
SET(SIMTK_HEADERS_DIR ${OPENSIM_INSTALL_DIR}/sdk/include/SimTK/include)
SET(OPENSIM_HEADERS_DIR ${OPENSIM_INSTALL_DIR}/sdk/include)
INCLUDE_DIRECTORIES(${SIMTK_HEADERS_DIR} ${OPENSIM_HEADERS_DIR})
# Libraries and dlls
SET(OPENSIM_LIBS_DIR ${OPENSIM_INSTALL_DIR}/sdk/lib ${OPENSIM_INSTALL_DIR}/lib)
SET(OPENSIM_DLLS_DIR ${OPENSIM_INSTALL_DIR}/bin)

LINK_DIRECTORIES(${OPENSIM_LIBS_DIR} ${OPENSIM_DLLS_DIR})

# Namespace
SET(NameSpace "OpenSim_" CACHE STRING "Prefix for simtk lib names, includes trailing '_'. Leave empty to use stock SimTK libraries.")

SET(OSIM_LIBRARIES
	debug osimSimulation_d	optimized osimSimulation
	debug osimActuators_d	optimized osimActuators
	debug osimCommon_d		optimized osimCommon
	debug osimAnalyses_d	optimized osimAnalyses
	debug osimTools_d		optimized osimTools
	debug ${NameSpace}SimTKcommon_d optimized ${NameSpace}SimTKcommon
	debug ${NameSpace}SimTKmath_d optimized  ${NameSpace}SimTKmath
	debug ${NameSpace}SimTKsimbody_d optimized ${NameSpace}SimTKsimbody
	${PLATFORM_LIBS}
)

IF(WIN32)
	SET(PLATFORM_LIBS  pthreadVC2)
ELSE (WIN32)
  SET(NameSpace "")
  IF(APPLE)
	SET(PLATFORM_LIBS  SimTKAtlas)
  ELSE(APPLE)
	SET(PLTAFORM_LIBS SimTKAtlas_Lin_generic)
  ENDIF(APPLE)
ENDIF (WIN32)

# This block copies the additional files into the running directory
# For example vtp, obj files. Add to the end for more extensions
FILE(GLOB DATA_FILES *.vtp *.obj *.osim)
FOREACH (dataFile ${DATA_FILES})
 ADD_CUSTOM_COMMAND(
    TARGET ${TARGET}
    COMMAND ${CMAKE_COMMAND}
    ARGS -E copy
    ${dataFile}
    ${OptimizationExample_BINARY_DIR})
 ENDFOREACH (dataFile) 

MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX)
MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH)
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)

FILE(GLOB COMMON_FILES Common/*.h Common/*.cpp Common/*.c)
#SET_SOURCE_FILES_PROPERTIES( ${CMAKE_SOURCE_DIR}/Common/cmaes.c PROPERTIES LANGUAGE CXX )
#SET(SOURCE ${COMMON_FILES})


SET(ENABLE_MPI OFF CACHE BOOL "Enable the use of MPI libraries for parallel computing if installed.)")

SET(MPI_INSTALL_DIR "C:/Program Files/Microsoft HPC Pack 2012" 
		CACHE PATH "Top-level directory of MPI install")
		

IF(ENABLE_MPI)
	IF(UNIX)
		INCLUDE_DIRECTORIES(/usr/include/mpi)
		LINK_DIRECTORIES(/usr/bin /usr/lib)
		add_definitions(-DUSE_MPI)
		SET(MPI_LIBRARIES mpich)
	ENDIF(UNIX)
	IF(WIN32)
		INCLUDE_DIRECTORIES(${MPI_INSTALL_DIR}/Inc)
		LINK_DIRECTORIES(${MPI_INSTALL_DIR}/Bin ${MPI_INSTALL_DIR}/Lib/amd64)
		add_definitions(-DUSE_MPI)
		SET(MPI_LIBRARIES msmpi)
	ENDIF(WIN32)
ENDIF(ENABLE_MPI)


# Other than Windows we can debug without debuggable SimTK libraries
IF(WIN32)
SET(CMAKE_DEBUG_POSTFIX "_d" CACHE INTERNAL "" FORCE)
ELSE(WIN32)
SET(CMAKE_DEBUG_POSTFIX "_d" CACHE STRING "Suffix for debug libraries")
ENDIF(WIN32)

ADD_SUBDIRECTORY(Common)
ADD_SUBDIRECTORY(nominal)