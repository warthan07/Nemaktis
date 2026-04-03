##################
# General settings
##################

cmake_minimum_required(VERSION 3.15)
project(nemaktis)

# generate json file for completion with YCM
set(CMAKE_EXPORT_COMPILE_COMMANDS 1) 

# compile flags
if(WIN32)
	set(CMAKE_CXX_FLAGS_RELEASE "/O2")
else()
	set(CMAKE_CXX_FLAGS_RELEASE "-O3")
	set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -g")
endif()

if(DEFINED ENV{Python_NumPy_INCLUDE_DIR})
	set(Python_NumPy_INCLUDE_DIR $ENV{Python_NumPy_INCLUDE_DIR})
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

find_package(OpenMP REQUIRED)

find_package(VTK REQUIRED COMPONENTS CommonCore CommonDataModel IOXML)

cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)
find_package(SWIG REQUIRED)
include(UseSWIG)

find_package(Python 3.12 EXACT REQUIRED COMPONENTS Interpreter Development.Module NumPy)


##########################################
# Include CMakeLists.txt from subdirectory
##########################################

add_subdirectory(BeamPropagationSolver)
add_subdirectory(RayTracingSolver)
