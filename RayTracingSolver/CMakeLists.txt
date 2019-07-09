##################
# Project settings
##################

project(rt-solver)
cmake_minimum_required(VERSION 2.8.8)

# generate json file for completion with YCM
set(CMAKE_EXPORT_COMPILE_COMMANDS 1) 

if(DEFINED ENV{MYROOT})
	link_directories($ENV{MYROOT}/lib)
endif()
if(DEFINED ENV{PREFIX})
	set(CMAKE_INSTALL_PREFIX $ENV{PREFIX})
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Necessary libraries
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/scripts)
find_package(Eigen3 REQUIRED)

if(WIN32)
	set(FFTW_USE_STATIC_LIBS ON)
endif()
find_package(FFTW REQUIRED)

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS filesystem regex)
if (WIN32)
  add_definitions( -DBOOST_ALL_NO_LIB )
  add_definitions( -DBOOST_ALL_DYN_LINK )
endif()

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

find_package(VTK REQUIRED COMPONENTS vtkCommonCore vtkCommonDataModel vtkIOXML)
include(${VTK_USE_FILE})


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


#################################################
# Custom target for building the doc with doxygen
#################################################

add_custom_target(doc
	COMMAND doxygen Doxyfile
	COMMAND ln -sf html/index.html index.html
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
	COMMENT "Building documentation with Doxygen"
)


#############################################
# Custom target for showing the project stats
#############################################

file(GLOB_RECURSE PROJECTFILES src/*.cpp src/*.h)
file(GLOB_RECURSE EXCLUDE json.h)
list(REMOVE_ITEM PROJECTFILES ${EXCLUDE})

add_custom_target(stats
	COMMAND wc ${PROJECTFILES}
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src
	COMMENT "Statistics of the project:"
)


##############################
# Definition of dirs and files
##############################

set(UTILITIESLIB_DIR src/UtilitiesLib)
file(GLOB UTILITIESLIB_CPP ${UTILITIESLIB_DIR}/*.cpp)

set(PROPAGATIONLIB_DIR src/PropagationLib)
file(GLOB PROPAGATIONLIB_CPP ${PROPAGATIONLIB_DIR}/*.cpp)

set(MAPPINGLIB_DIR src/MappingLib)
file(GLOB MAPPINGLIB_CPP ${MAPPINGLIB_DIR}/*.cpp)

set(SIMULATIONLIB_DIR src/SimulationLib)
file(GLOB SIMULATIONLIB_CPP ${SIMULATIONLIB_DIR}/*.cpp)

INCLUDE_DIRECTORIES(
	${UTILITIESLIB_DIR}
	${PROPAGATIONLIB_DIR}
	${MAPPINGLIB_DIR}
	${SIMULATIONLIB_DIR}
	${EIGEN3_INCLUDE_DIR}
	${Boost_INCLUDE_DIRS}
)


#################
# The main target
#################

add_executable(${PROJECT_NAME} src/main.cpp)

set_property(
	TARGET ${PROJECT_NAME}
	PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
)
add_custom_target(run
	COMMAND ./${PROJECT_NAME}
	DEPENDS ${PROJECT_NAME}
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	COMMENT "Run with ${CMAKE_BUILD_TYPE} configuration"
)


#####################
# The library targets
#####################

add_library(Utilities ${UTILITIESLIB_CPP})
add_library(Propagation ${PROPAGATIONLIB_CPP})
add_library(Mapping ${MAPPINGLIB_CPP})
add_library(Simulation ${SIMULATIONLIB_CPP})


###########################################
# We link the main target and the libraries
###########################################

target_link_libraries(${PROJECT_NAME} Simulation)
target_link_libraries(Simulation Propagation Mapping  ${Boost_LIBRARIES} ${FFTW_LIBRARIES})
target_link_libraries(Mapping Utilities)
target_link_libraries(Propagation Utilities)
target_link_libraries(Utilities ${VTK_LIBRARIES})


#################
# Install targets
#################

install(TARGETS ${PROJECT_NAME}
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
)
