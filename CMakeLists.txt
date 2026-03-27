##################
# General settings
##################

cmake_minimum_required(VERSION 3.15)
project(nemaktis)

# generate json file for completion with YCM
set(CMAKE_EXPORT_COMPILE_COMMANDS 1) 

if(DEFINED ENV{MYROOT})
	link_directories($ENV{MYROOT}/lib)
endif()
# if(DEFINED ENV{PREFIX})
# 	set(CMAKE_INSTALL_PREFIX $ENV{PREFIX})
# endif()
# if(DEFINED ENV{CONDA_PREFIX})
# 	set(CMAKE_PREFIX_PATH $ENV{CONDA_PREFIX})
# 	include_directories($ENV{CONDA_PREFIX}/include)
# 	link_directories($ENV{CONDA_PREFIX}/lib)
# endif()

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


###########################################################################
# Custom targets to ease the switching between debug or release compilation
###########################################################################

if(WIN32)
	set(CMAKE_CXX_FLAGS_RELEASE "/O2")
else()
	set(CMAKE_CXX_FLAGS_RELEASE "-O3")
	set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -g")
endif()
add_custom_target(debug
	COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
	COMMAND +${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
	COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
)
add_custom_target(release
	COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
	COMMAND +${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
	COMMENT "Switch CMAKE_BUILD_TYPE to Release"
)


##########################################
# Include CMakeLists.txt from subdirectory
##########################################

add_subdirectory(BeamPropagationSolver)
add_subdirectory(RayTracingSolver)
add_subdirectory(HighLevelPythonInterface)
