# - Try to find Petsc lib
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(Petsc)
#
# Once done this will define
#
#  PETSC_FOUND - system has Petsc lib with correct version
#  PETSC_INCLUDE_DIRS - the Petsc include directory
#  PETSC_LIBRARIES   - the Petsc library.

# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Copyright (c) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
# Redistribution and use is allowed according to the terms of the 2-clause BSD license.

# find out the size of a pointer. this is required to only search for
# libraries in the directories relevant for the architecture
if (CMAKE_SIZEOF_VOID_P)
  math (EXPR _BITS "8 * ${CMAKE_SIZEOF_VOID_P}")
endif (CMAKE_SIZEOF_VOID_P)

# look for a system-wide BLAS library
find_package(BLAS QUIET)
set(PETSC_BLAS_LIBRARY "")
if (BLAS_FOUND)
    list(APPEND PETSC_BLAS_LIBRARY "${BLAS_LIBRARIES}")
elseif(PETSC_ROOT)
    find_library(PETST_BLAS_LIBRARY
        NAME "blas"
        PATH ${PETSC_ROOT}
        PATH_SUFFIXES "lib" "lib${_BITS}" "lib/${CMAKE_LIBRARY_ARCHITECTURE}"
        NO_DEFAULT_PATH)
endif()
# print message if there was still no blas found!
if(NOT BLAS_FOUND AND NOT PETSC_BLAS_LIBRARY)
  message(STATUS "BLAS not found but required for PETSC")
  return()
endif()
list(APPEND CMAKE_REQUIRED_LIBRARIES "${PETSC_BLAS_LIBRARY}")
find_package(LAPACK QUIET)

set(PETSC_LAPACK_LIBRARY "")
if (LAPACK_FOUND)
    list(APPEND PETSC_LAPACK_LIBRARY "${LAPACK_LIBRARIES}")
elseif(PETSC_ROOT)
    find_library(PETST_LAPACK_LIBRARY
        NAME "lapack"
        PATH ${PETSC_ROOT}
        PATH_SUFFIXES "lib" "lib${_BITS}" "lib/${CMAKE_LIBRARY_ARCHITECTURE}"
        NO_DEFAULT_PATH)
endif()
# print message if there was still no blas found!
if(NOT LAPACK_FOUND AND NOT PETSC_LAPACK_LIBRARY)
  message(STATUS "LAPACK not found but required for PETSC")
  return()
endif()
list(APPEND CMAKE_REQUIRED_LIBRARIES "${PETSC_LAPACK_LIBRARY}")

find_package(X11 QUIET)
if (X11_FOUND)
    list(APPEND PETSC_X11_LIBRARY "${X11_LIBRARIES}")
endif()
list(APPEND CMAKE_REQUIRED_LIBRARIES "${PETSC_X11_LIBRARY}")
# only probe if we haven't a path in our cache
if (Petsc_ROOT)
 set (PETSC_ROOT "${Petsc_ROOT}")
endif (Petsc_ROOT)
if (NOT PETSC_NORMAL_INCLUDE_DIR)
	find_path (PETSC_NORMAL_INCLUDE_DIR
	  NAMES "petsc.h"
	  PATHS ${PETSC_ROOT}
	  PATH_SUFFIXES "petsc-3.4.4" "include" "petsc"
	  )
endif (NOT PETSC_NORMAL_INCLUDE_DIR)

# if parallel computing is explicitly enabled - reuse the paths and links from
# OpmMainLib + OpmFind
# this needs to be called explicitly as FindPetsc runs before OpmMainLib starts
# looking for MPI (for some reason). Ideally this isn't necessary
if(USE_MPI)
    find_package(MPI)
endif()

set(PETSC_MPI_FOUND ${MPI_FOUND})
set(PETSC_MPI_INCLUDE_DIRS ${MPI_INCLUDE_PATH})
set(PETSC_MPI_LIBRARIES ${MPI_LIBRARIES})

# fallback - use the petsc provided implementation of serial MPI.
if (NOT PETSC_MPI_INCLUDE_DIRS)
    message(STATUS "Building without MPI support - looking for PETSc provided serial implementation")
    find_path (PETSC_MPI_INCLUDE_DIRS
        NAMES "mpi.h"
        PATHS ${PETSC_ROOT}/include
        PATH_SUFFIXES "mpiuni"
        )

    if(PETSC_MPI_INCLUDE_DIRS)
        # not setting special linkage
        set(PETSC_MPI_FOUND 1)
    endif(PETSC_MPI_INCLUDE_DIRS)

endif(NOT PETSC_MPI_INCLUDE_DIRS)

# couldn't find any usable mpi implementation - abort
if(NOT PETSC_MPI_FOUND)
    message(STATUS "Could not find any suitable MPI implementation. Is PETSC_ROOT set?")
    return()
endif()

# look for actual Petsc library
if (NOT PETSC_LIBRARY)
  find_library(PETSC_LIBRARY
    NAMES "petsc-3.4.3" "petsc-3.4.4" "petsc" 
    PATHS ${PETSC_ROOT}
    PATH_SUFFIXES "lib" "lib${_BITS}" "lib/${CMAKE_LIBRARY_ARCHITECTURE}"
    )
endif()
if(NOT PETSC_LIBRARY)
  message(STATUS "Could not find the PETSc library")
  return()
endif()

list(APPEND CMAKE_REQUIRED_LIBRARIES "${PETSC_LIBRARY}" "${PETSC_MPI_LIBRARIES}")

if (PETSC_MPI_INCLUDE_DIRS AND PETSC_NORMAL_INCLUDE_DIR)
    list (APPEND PETSC_INCLUDE_DIR ${PETSC_MPI_INCLUDE_DIRS} ${PETSC_NORMAL_INCLUDE_DIR})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Petsc DEFAULT_MSG PETSC_INCLUDE_DIR PETSC_LIBRARY)
mark_as_advanced(PETSC_INCLUDE_DIR PETSC_LIBRARY)

# if both headers and library are found, store results
if(PETSC_FOUND)
  set(PETSC_INCLUDE_DIRS ${PETSC_INCLUDE_DIR})

  set(PETSC_LIBRARIES ${PETSC_LIBRARY})

  if (PETSC_BLAS_LIBRARY)
    list(APPEND PETSC_LIBRARIES ${PETSC_BLAS_LIBRARY})
  endif()
  if (PETSC_LAPACK_LIBRARY)
    list(APPEND PETSC_LIBRARIES ${PETSC_LAPACK_LIBRARY})
  endif()
  if (PETSC_X11_LIBRARY)
    list(APPEND PETSC_LIBRARIES ${PETSC_X11_LIBRARY})
  endif()
endif()

