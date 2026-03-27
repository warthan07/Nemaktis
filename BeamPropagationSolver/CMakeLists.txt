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
set(PYTHONBACKEND_CPP ${PYTHONBACKEND_DIR}/bpm_backend.cpp)
set(PYTHONBACKEND_SWIG ${PYTHONBACKEND_DIR}/bpm_backend.i)
set_property(SOURCE ${PYTHONBACKEND_SWIG} PROPERTY CPLUSPLUS ON)

include_directories(
	${BPMLIB_DIR} ${PHYSICSLIB_DIR} ${POSTPROCESSORLIB_DIR} 
	${SETTINGSLIB_DIR} ${UTILITIESLIB_DIR} ${PYTHONBACKEND_DIR}
	${MUPARSER_INCLUDE_DIRS} ${FFTW_INCLUDES}
)


##################
# The main targets
##################

add_executable(bpm-solver src/main.cpp)
set_property(
	TARGET bpm-solver
	PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
)
add_custom_target(run-bpm
	COMMAND ./bpm-solver
	DEPENDS bpm-solver
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	COMMENT "Run with ${CMAKE_BUILD_TYPE} configuration"
)

#####################
# The library targets
#####################

add_library(bpm-lib STATIC ${BPMLIB_CPP} ${PHYSICSLIB_CPP} ${POSTPROCESSORLIB_CPP}
	${SETTINGSLIB_CPP} ${UTILITIESLIB_CPP})
set_property(TARGET bpm-lib PROPERTY POSITION_INDEPENDENT_CODE ON)

set(SWIG_USE_SWIG_DEPENDENCIES true)
swig_add_library(bpm_backend LANGUAGE python
	SOURCES ${PYTHONBACKEND_SWIG} ${PYTHONBACKEND_CPP})
get_property(PYTHONBACKEND_PY TARGET bpm_backend PROPERTY SWIG_SUPPORT_FILES)
if(WIN32)
	set_property(TARGET bpm_backend
		PROPERTY SUFFIX ".${Python_SOABI}.pyd")
else()
	set_property(TARGET bpm_backend
		PROPERTY SUFFIX ".${Python_SOABI}${CMAKE_SHARED_MODULE_SUFFIX}")
endif()


############################################
# We link the main targets and the libraries
############################################

target_link_libraries(bpm-solver PRIVATE bpm-lib)
target_link_libraries(bpm_backend PRIVATE bpm-lib Python::NumPy)
target_link_libraries(bpm-lib PUBLIC
	VTK::CommonCore VTK::CommonDataModel VTK::IOXML
	Boost::filesystem Boost::regex Eigen3::Eigen OpenMP::OpenMP_CXX 
	${FFTW_LIBRARIES} ${MUPARSER_LIBRARIES})

install(FILES ${PYTHONBACKEND_PY} DESTINATION ${SKBUILD_PLATLIB_DIR})
install(TARGETS bpm_backend DESTINATION ${SKBUILD_PLATLIB_DIR})
install(TARGETS bpm-solver DESTINATION ${SKBUILD_SCRIPTS_DIR})
