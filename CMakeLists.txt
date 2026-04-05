##################
# General settings
##################

cmake_minimum_required(VERSION 3.15)
project(nemaktis)

# generate json file for LSP
set(CMAKE_EXPORT_COMPILE_COMMANDS 1) 

# compile flags
if(WIN32)
	set(CMAKE_CXX_FLAGS_RELEASE "/O2")
else()
	set(CMAKE_CXX_FLAGS_RELEASE "-O3")
	set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -g")
endif()

# fix pixi dev env (omp path for clangd and numpy path)
if(DEFINED ENV{OMP_INCLUDE_DIR})
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I$ENV{OMP_INCLUDE_DIR}")
endif()
if(DEFINED ENV{NUMPY_INCLUDE_DIR})
	set(Python_NumPy_INCLUDE_DIR $ENV{NUMPY_INCLUDE_DIR})
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Necessary libraries
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
find_package(MuParser REQUIRED)
find_package(Eigen3 REQUIRED)

if(WIN32)
	set(FFTW_USE_STATIC_LIBS ON)
endif()
find_package(FFTW REQUIRED)

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.85.0 CONFIG EXACT REQUIRED COMPONENTS filesystem regex system)
if (WIN32)
	add_definitions( -DBOOST_ALL_NO_LIB )
	add_definitions( -DBOOST_ALL_DYN_LINK )
endif()

find_package(OpenMP REQUIRED COMPONENTS CXX)

find_package(VTK REQUIRED COMPONENTS CommonCore CommonDataModel IOXML)

find_package(Python 3.12 EXACT REQUIRED COMPONENTS Interpreter Development.Module NumPy)

# cmake_policy(SET CMP0078 NEW)
# cmake_policy(SET CMP0086 NEW)
find_package(SWIG REQUIRED)
include(UseSWIG)
set(SWIG_USE_SWIG_DEPENDENCIES true)
if (WIN32)
	set(SWIG_MODULE_SUFFIX .${Python_SOABI}.pyd)
else()
	set(SWIG_MODULE_SUFFIX .${Python_SOABI}${CMAKE_SHARED_MODULE_SUFFIX})
endif()


##########################################
# Include CMakeLists.txt from subdirectory
##########################################

add_subdirectory(BeamPropagationSolver)
add_subdirectory(RayTracingSolver)
