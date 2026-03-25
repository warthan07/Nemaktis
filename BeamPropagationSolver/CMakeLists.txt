##################
# General settings
##################

cmake_minimum_required(VERSION 3.15)
project(bpm-solver)

# generate json file for completion with YCM
set(CMAKE_EXPORT_COMPILE_COMMANDS 1) 

if(DEFINED ENV{MYROOT})
	link_directories($ENV{MYROOT}/lib)
endif()
# if(DEFINED ENV{PREFIX})
# 	set(CMAKE_INSTALL_PREFIX $ENV{PREFIX})
# endif()
if(DEFINED ENV{CONDA_PREFIX})
	set(CMAKE_PREFIX_PATH $ENV{CONDA_PREFIX})
	include_directories($ENV{CONDA_PREFIX}/include)
	link_directories($ENV{CONDA_PREFIX}/lib)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Necessary libraries
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/scripts)
find_package(MuParser REQUIRED)
find_package(Eigen3 3.3 REQUIRED)

if(WIN32)
	set(FFTW_USE_STATIC_LIBS ON)
endif()
find_package(FFTW REQUIRED)

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.82.0 CONFIG EXACT REQUIRED COMPONENTS filesystem regex system)
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

find_package(Python 3.9 EXACT REQUIRED COMPONENTS Interpreter Development.Module NumPy)


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

set(BPMLIB_DIR src/BeamPropagationLib)
file(GLOB BPMLIB_CPP ${BPMLIB_DIR}/*.cpp)

set(PHYSICSLIB_DIR src/PhysicsLib)
file(GLOB PHYSICSLIB_CPP ${PHYSICSLIB_DIR}/*.cpp)

set(POSTPROCESSORLIB_DIR src/PostprocessorLib)
file(GLOB POSTPROCESSORLIB_CPP ${POSTPROCESSORLIB_DIR}/*.cpp)

set(SETTINGSLIB_DIR src/SettingsLib)
file(GLOB SETTINGSLIB_CPP ${SETTINGSLIB_DIR}/*.cpp)

set(UTILITIESLIB_DIR src/UtilitiesLib)
file(GLOB UTILITIESLIB_CPP ${UTILITIESLIB_DIR}/*.cpp)

set(PYTHONBACKEND_BASENAME bpm_backend)
set(PYTHONBACKEND_DIR src/PythonBackend)
set(PYTHONBACKEND_CPP ${PYTHONBACKEND_DIR}/${PYTHONBACKEND_BASENAME}.cpp)
set(PYTHONBACKEND_SWIG ${PYTHONBACKEND_DIR}/${PYTHONBACKEND_BASENAME}.i)
set_property(SOURCE ${PYTHONBACKEND_SWIG} PROPERTY CPLUSPLUS ON)

include_directories(
	${BPMLIB_DIR}
	${PHYSICSLIB_DIR}
	${POSTPROCESSORLIB_DIR}
	${SETTINGSLIB_DIR}
	${UTILITIESLIB_DIR}
	${PYTHONBACKEND_DIR}
	${MUPARSER_INCLUDE_DIRS}
	${FFTW_INCLUDES}
	${Boost_INCLUDE_DIRS}
	${Python_INCLUDE_DIRS}
	${Python_NumPy_INCLUDE_DIRS}
	${EIGEN3_INCLUDE_DIR}
)


##################
# The main targets
##################

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

add_library(bpm-lib STATIC ${BPMLIB_CPP} ${PHYSICSLIB_CPP} ${POSTPROCESSORLIB_CPP}
	${SETTINGSLIB_CPP} ${UTILITIESLIB_CPP})

set(SWIG_USE_SWIG_DEPENDENCIES true)
swig_add_library(${PYTHONBACKEND_BASENAME} LANGUAGE python
	SOURCES ${PYTHONBACKEND_SWIG} ${PYTHONBACKEND_CPP})
get_property(PYTHONBACKEND_PY TARGET ${PYTHONBACKEND_BASENAME} PROPERTY SWIG_SUPPORT_FILES)
message(${PYTHONBACKEND_PY})
if(WIN32)
	set_property(TARGET ${PYTHONBACKEND_BASENAME}
		PROPERTY SUFFIX ".${Python_SOABI}.pyd")
else()
	set_property(TARGET ${PYTHONBACKEND_BASENAME}
		PROPERTY SUFFIX ".${Python_SOABI}${CMAKE_SHARED_MODULE_SUFFIX}")
endif()


############################################
# We link the main targets and the libraries
############################################

target_link_libraries(${PROJECT_NAME}
	bpm-lib OpenMP::OpenMP_CXX Eigen3::Eigen)
target_link_libraries(${PYTHONBACKEND_BASENAME}
	bpm-lib Python::NumPy OpenMP::OpenMP_CXX Eigen3::Eigen)
target_link_libraries(bpm-lib 
	${Boost_LIBRARIES} ${VTK_LIBRARIES} ${FFTW_LIBRARIES} ${MUPARSER_LIBRARIES})

install(FILES ${PYTHONBACKEND_PY} DESTINATION .)
install(TARGETS ${PYTHONBACKEND_BASENAME} DESTINATION .)
install(TARGETS ${PROJECT_NAME} bpm-lib 
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
)