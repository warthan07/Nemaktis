##################
# General settings
##################

project(bpm-solver)
cmake_minimum_required(VERSION 3.11.2)

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
find_package(MuParser REQUIRED)

if(WIN32)
	set(FFTW_USE_STATIC_LIBS ON)
endif()
find_package(FFTW REQUIRED)

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.67 EXACT REQUIRED COMPONENTS filesystem regex system)
if (WIN32)
  add_definitions( -DBOOST_ALL_NO_LIB )
  add_definitions( -DBOOST_ALL_DYN_LINK )
endif()
message(STATUS "Boost libs: ${Boost_LIBRARIES}")

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

find_package(VTK REQUIRED COMPONENTS vtkCommonCore vtkCommonDataModel vtkIOXML)
include(${VTK_USE_FILE})

cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

set(Python3_FIND_REGISTRY LAST)
find_package(Python3 REQUIRED COMPONENTS Development NumPy)


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

set(PYTHONBACKEND_DIR src/PythonBackend)
set(PYTHONBACKEND_BASENAME bpm_backend)
set(PYTHONBACKEND_CPP ${PYTHONBACKEND_DIR}/${PYTHONBACKEND_BASENAME}.cpp)

if(WIN32)
	set(PYTHON_INSTALL_FILES
        	${CMAKE_CURRENT_BINARY_DIR}/${PYTHONBACKEND_BASENAME}.py
	        ${CMAKE_CURRENT_BINARY_DIR}/_${PYTHONBACKEND_BASENAME}.pyd)
else()
	set(PYTHON_INSTALL_FILES
        	${CMAKE_CURRENT_BINARY_DIR}/${PYTHONBACKEND_BASENAME}.py
	        ${CMAKE_CURRENT_BINARY_DIR}/_${PYTHONBACKEND_BASENAME}.so)
endif()
set(SETUP_PY_IN ${PYTHONBACKEND_DIR}/setup.py.in)
set(SETUP_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
configure_file(${SETUP_PY_IN} ${SETUP_PY_OUT})

INCLUDE_DIRECTORIES(
	${BPMLIB_DIR}
	${PHYSICSLIB_DIR}
	${POSTPROCESSORLIB_DIR}
	${SETTINGSLIB_DIR}
	${UTILITIESLIB_DIR}
	${PYTHONBACKEND_DIR}
	${MUPARSER_INCLUDE_DIRS}
	${FFTW_INCLUDES}
	${Boost_INCLUDE_DIRS}
	${Python3_INCLUDE_DIRS}
	${Python3_NumPy_INCLUDE_DIRS}
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

add_library(BPM ${BPMLIB_CPP})
add_library(Physics ${PHYSICSLIB_CPP})
add_library(Postprocessor ${POSTPROCESSORLIB_CPP})
add_library(Settings ${SETTINGSLIB_CPP})
add_library(Utilities ${UTILITIESLIB_CPP})

add_library(bpm-lib ${BPMLIB_CPP} ${PHYSICSLIB_CPP} ${POSTPROCESSORLIB_CPP}
	${SETTINGSLIB_CPP} ${UTILITIESLIB_CPP} ${PYTHONBACKEND_CPP})
set_target_properties(bpm-lib PROPERTIES
	OUTPUT_NAME ${PROJECT_NAME}
	PUBLIC_HEADER ${PYTHONBACKEND_DIR}/${PYTHONBACKEND_BASENAME}.h)

set_property(SOURCE ${PYTHONBACKEND_DIR}/${PYTHONBACKEND_BASENAME}.i PROPERTY CPLUSPLUS ON)
swig_add_library(${PYTHONBACKEND_BASENAME} LANGUAGE python
	SOURCES ${PYTHONBACKEND_DIR}/${PYTHONBACKEND_BASENAME}.i)


############################################
# We link the main targets and the libraries
############################################

target_link_libraries(${PROJECT_NAME} Settings)
target_link_libraries(Settings BPM Postprocessor)
target_link_libraries(Postprocessor Physics
	${Boost_LIBRARIES} ${VTK_LIBRARIES} ${FFTW_LIBRARIES})
target_link_libraries(BPM Utilities Physics)
target_link_libraries(Physics ${MUPARSER_LIBRARIES})

target_link_libraries(bpm-lib
	${Boost_LIBRARIES} ${VTK_LIBRARIES} ${FFTW_LIBRARIES} ${MUPARSER_LIBRARIES})
swig_link_libraries(${PYTHONBACKEND_BASENAME} bpm-lib ${Python3_LIBRARIES})


#################
# Install targets
#################

add_custom_target(python-install
	DEPENDS bpm_backend
	COMMAND python ${SETUP_PY_OUT} install
)
install(TARGETS ${PROJECT_NAME} bpm-lib
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
)
