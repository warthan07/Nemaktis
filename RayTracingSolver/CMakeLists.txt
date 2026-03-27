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
	${UTILITIESLIB_DIR} ${PROPAGATIONLIB_DIR}
    ${MAPPINGLIB_DIR} ${SIMULATIONLIB_DIR})


#################
# The main target
#################

add_executable(rt-solver src/main.cpp)
set_property(
	TARGET rt-solver
	PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
add_custom_target(run-rt
	COMMAND ./rt-solver
	DEPENDS rt-solver
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	COMMENT "Run with ${CMAKE_BUILD_TYPE} configuration")


#####################
# The library targets
#####################

add_library(rt-lib STATIC
    ${UTILITIESLIB_CPP} ${PROPAGATIONLIB_CPP}
    ${MAPPINGLIB_CPP} ${SIMULATIONLIB_CPP})


###########################################
# We link the main target and the libraries
###########################################

target_link_libraries(rt-solver PRIVATE rt-lib)
target_link_libraries(rt-lib PUBLIC
	VTK::CommonCore VTK::CommonDataModel VTK::IOXML
	Boost::filesystem Boost::regex Eigen3::Eigen ${FFTW_LIBRARIES})


#################
# Install targets
#################

install(TARGETS rt-solver DESTINATION ${SKBUILD_SCRIPTS_DIR})