# - Find MuParser
# muParser is an extensible high performance math expression parser library written in C++
# http://muparser.sourceforge.net
#
# The module defines the following variables:
#  MUPARSER_FOUND        - True if MuParser found.
#  MUPARSER_INCLUDE_DIRS - where to find muParser.h, etc.
#  MUPARSER_LIBRARIES    - List of libraries when using MuParser.
#
#=============================================================================
# Copyright (C) 2005-2013 EDF-EADS-Phimeca
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)

if ( MUPARSER_INCLUDE_DIR AND MUPARSER_LIBRARIES )
  # Already in cache, be silent
  set ( MuParser_FIND_QUIETLY TRUE )
endif ()

find_path ( MUPARSER_INCLUDE_DIR muParser.h
            HINTS
            ${MUPARSER_ROOT_DIR}/include
            PATHS
            C:/muparser/include
            PATH_SUFFIXES muParser )            
set ( MUPARSER_INCLUDE_DIRS ${MUPARSER_INCLUDE_DIR} )

find_library ( MUPARSER_LIBRARY
               NAMES muparser
               HINTS
               ${MUPARSER_ROOT_DIR}/lib
               PATH_SUFFIXES muparser )    
set ( MUPARSER_LIBRARIES ${MUPARSER_LIBRARY} )

# root dir
# try to guess root dir from include dir
if ( MUPARSER_INCLUDE_DIR )
  string ( REGEX REPLACE "(.*)/include.*" "\\1" MUPARSER_ROOT_DIR ${MUPARSER_INCLUDE_DIR} )

# try to guess root dir from library dir
elseif ( MUPARSER_LIBRARY )
  string ( REGEX REPLACE "(.*)/lib[/|32|64].*" "\\1" MUPARSER_ROOT_DIR ${MUPARSER_LIBRARY} )
endif ()

# handle REQUIRED and QUIET options
include ( FindPackageHandleStandardArgs )
find_package_handle_standard_args ( MuParser DEFAULT_MSG MUPARSER_LIBRARY
  MUPARSER_LIBRARIES
  MUPARSER_INCLUDE_DIR
  MUPARSER_INCLUDE_DIRS
  MUPARSER_ROOT_DIR
)

mark_as_advanced (
  MUPARSER_LIBRARY
  MUPARSER_LIBRARIES
  MUPARSER_INCLUDE_DIR
  MUPARSER_INCLUDE_DIRS
  MUPARSER_ROOT_DIR
)
