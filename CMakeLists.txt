# CMakeLists.txt for compiling MYSTRAN with gfortran in a GNU environment
# based on an older CMakeLists.txt by ceanwang@gmail.com
# made to work again by Bruno Borges Paschoalinoto (2020)

# set up basic project info
cmake_minimum_required(VERSION 3.18)
include(CheckFunctionExists)

project(Mystran)
enable_language(Fortran)

if(WIN32)
  enable_language(RC)
endif()

# basic compiler and output options
set(CMAKE_SOURCE_DIR "${PROJECT_SOURCE_DIR}/Source")
set(BLAS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/BLAS")
set(PROJECT_BINARY_DIR "${PROJECT_SOURCE_DIR}/Binaries")
set(CMAKE_FLAGS "-c -fbacktrace")

# set(CMAKE_BUILD_TYPE Debug)

# set some dirs
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/mod")
include_directories("${CMAKE_SOURCE_DIR}/INCLUDE")

# uncomment this to debug
# set(CMAKE_VERBOSE_MAKEFILE true)

# register the Profile build type, just in case
set(CMAKE_CONFIGURATION_TYPES
  Debug
  Release
  RelWithDebInfo
  MinSizeRel
  Profiling
  Deterministic
  CACHE STRING "" FORCE
)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

# suppress cmake warnings for superlu
if(NOT DEFINED CMAKE_SUPPRESS_DEVELOPER_WARNINGS)
  set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE INTERNAL "No dev warnings")
endif()

# build static libraries throughout this tree
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static libraries" FORCE)

# recommend the appropriate make command
if(WIN32)
  set(RECOMMENDED_MAKE "mingw32-make")
else()
  set(RECOMMENDED_MAKE "make")
endif()

# submodules (i.e. SuperLU)
find_package(Git QUIET)

if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  option(GIT_SUBMODULE "Check submodules during build" ON)

  if(GIT_SUBMODULE)
    message(STATUS "Updating submodules")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive --force
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT
    )

    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(
        FATAL_ERROR
        "git submodule update --init --recursive --force failed"
      )
    endif()
  endif()
endif()

set(SUPERLU_DIR "${PROJECT_SOURCE_DIR}/superlu")
set(SUPERLU_PATCHES_DIR "${PROJECT_SOURCE_DIR}/superlu_patches")
set(SUPERLU_MT_DIR "${PROJECT_SOURCE_DIR}/superlu_mt")
set(SUPERLU_MT_PATCHES_DIR "${PROJECT_SOURCE_DIR}/superlu_mt_patches")
set(METIS_DIR "${PROJECT_SOURCE_DIR}/metis/METIS")
set(GKLIB_DIR "${PROJECT_SOURCE_DIR}/metis/GKlib")
set(METIS_PATCHES_DIR "${PROJECT_SOURCE_DIR}/metis/metis_patches")
set(GKLIB_PATCHES_DIR "${PROJECT_SOURCE_DIR}/metis/gklib_patches")

if(NOT DEFINED TPL_ENABLE_METISLIB)
  set(TPL_ENABLE_METISLIB ON CACHE BOOL "Enable METIS for SuperLU")
endif()

if(TPL_ENABLE_METISLIB)
  if(NOT EXISTS "${GKLIB_DIR}/CMakeLists.txt")
    message(FATAL_ERROR "The GKlib submodule was not downloaded.")
  endif()

  if(NOT EXISTS "${METIS_DIR}/CMakeLists.txt")
    message(FATAL_ERROR "The METIS submodule was not downloaded.")
  endif()

  message(STATUS "Applying patches to GKlib...")
  file(GLOB_RECURSE GKLIB_PATCH_FILES RELATIVE "${GKLIB_PATCHES_DIR}" "${GKLIB_PATCHES_DIR}/*")

  foreach(file IN LISTS GKLIB_PATCH_FILES)
    set(src "${GKLIB_PATCHES_DIR}/${file}")
    set(dst "${GKLIB_DIR}/${file}")

    get_filename_component(dst_dir "${dst}" DIRECTORY)
    file(MAKE_DIRECTORY "${dst_dir}")

    message(STATUS "  Patching: ${file}")
    file(COPY "${src}" DESTINATION "${dst_dir}")
  endforeach()

  message(STATUS "GKlib patched.")

  message(STATUS "Applying patches to METIS...")
  file(GLOB_RECURSE METIS_PATCH_FILES RELATIVE "${METIS_PATCHES_DIR}" "${METIS_PATCHES_DIR}/*")

  foreach(file IN LISTS METIS_PATCH_FILES)
    set(src "${METIS_PATCHES_DIR}/${file}")
    set(dst "${METIS_DIR}/${file}")

    get_filename_component(dst_dir "${dst}" DIRECTORY)
    file(MAKE_DIRECTORY "${dst_dir}")

    message(STATUS "  Patching: ${file}")
    file(COPY "${src}" DESTINATION "${dst_dir}")
  endforeach()

  message(STATUS "METIS patched.")

  set(SHARED OFF CACHE BOOL "Build vendored GKlib/METIS as static libraries" FORCE)
  set(GKLIB_BUILD_APPS OFF CACHE BOOL "Build GKlib applications" FORCE)
  set(METIS_BUILD_PROGRAMS OFF CACHE BOOL "Build METIS command line programs" FORCE)
  set(GKLIB_PATH "${GKLIB_DIR}" CACHE PATH "GKlib source directory" FORCE)

  # Use a binary directory that is disjoint from the source tree. The MYSTRAN
  # build is performed in-source (CMAKE_BINARY_DIR == CMAKE_SOURCE_DIR), so
  # using "${CMAKE_BINARY_DIR}/metis/METIS" would make METIS' binary dir alias
  # its source dir and cause the patched CMakeLists.txt to overwrite the source
  # include/metis.h on every reconfigure (forcing a full METIS rebuild).
  add_subdirectory("${GKLIB_DIR}" "${PROJECT_BINARY_DIR}/GKlib")
  add_subdirectory("${METIS_DIR}" "${PROJECT_BINARY_DIR}/METIS")

  set(TPL_METIS_LIBRARIES "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/libmetis.a;${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/libGKlib.a" CACHE STRING "METIS static libraries" FORCE)
  set(TPL_METIS_INCLUDE_DIRS "${PROJECT_BINARY_DIR}/METIS/include" CACHE PATH "METIS include directory" FORCE)

  message(STATUS "Applying patches to SUPERLU...")
  file(GLOB_RECURSE SUPERLU_PATCH_FILES RELATIVE "${SUPERLU_PATCHES_DIR}" "${SUPERLU_PATCHES_DIR}/*")

  foreach(file IN LISTS SUPERLU_PATCH_FILES)
    set(src "${SUPERLU_PATCHES_DIR}/${file}")
    set(dst "${SUPERLU_DIR}/${file}")

    get_filename_component(dst_dir "${dst}" DIRECTORY)
    file(MAKE_DIRECTORY "${dst_dir}")

    message(STATUS "  Patching: ${file}")
    file(COPY "${src}" DESTINATION "${dst_dir}")
  endforeach()

  message(STATUS "SuperLU patched.")
endif()

set(enable_examples OFF FORCE)
set(enable_tests OFF FORCE)
set(enable_doc OFF FORCE)

if(USE_SUPERLU_MT)
  message(STATUS "Will link against faster SuperLU_MT.")
  set(PLAT "_OPENMP" CACHE STRING "threading flavor _PTHREAD/_OPENMP" FORCE)

  # message(STATUS "Applying patches to SuperLU_MT...")
  # file(GLOB_RECURSE PATCH_FILES RELATIVE "${SUPERLU_MT_PATCHES_DIR}" "${SUPERLU_MT_PATCHES_DIR}/*")

  # foreach(file IN LISTS PATCH_FILES)
  # set(src "${SUPERLU_MT_PATCHES_DIR}/${file}")
  # set(dst "${SUPERLU_MT_DIR}/${file}")

  # # ensure the destination directory exists
  # get_filename_component(dst_dir "${dst}" DIRECTORY)
  # file(MAKE_DIRECTORY "${dst_dir}")
  # # print the file being copied
  # message(STATUS "  Patching: ${file}")

  # # copy the file
  # file(COPY "${src}" DESTINATION "${dst_dir}")
  # endforeach()

  # message(STATUS "SuperLU_MT patched.")
  if(NOT EXISTS "${SUPERLU_MT_DIR}/CMakeLists.txt")
    message(
      FATAL_ERROR
      "The SuperLU_MT submodule was not downloaded. You have to do it manually!"
    )
  endif()

  if(WIN32)
    # force static linking to OpenMP library
    # temporarily set CMAKE_FIND_LIBRARY_SUFFIXES to prefer .a over .dll.a
    set(CMAKE_FIND_LIBRARY_SUFFIXES_BACKUP ${CMAKE_FIND_LIBRARY_SUFFIXES})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")

    set(OpenMP_C_FLAGS "-fopenmp" CACHE STRING "" FORCE)
    set(OpenMP_C_LIB_NAMES "gomp" CACHE STRING "" FORCE)

    # find the static gomp library explicitly
    find_library(GOMP_STATIC_LIB NAMES libgomp.a gomp PATHS /mingw64/lib NO_DEFAULT_PATH)

    if(GOMP_STATIC_LIB)
      set(OpenMP_gomp_LIBRARY "${GOMP_STATIC_LIB}" CACHE FILEPATH "" FORCE)
    endif()

    find_package(OpenMP REQUIRED)

    # restore original library suffixes
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_BACKUP})
  endif()

  include_directories("${SUPERLU_MT_DIR}/SRC")
  add_subdirectory(${SUPERLU_MT_DIR})

  if(TPL_ENABLE_METISLIB)
    if(PLAT STREQUAL "_PTHREAD")
      add_dependencies(superlu_mt_PTHREAD metis GKlib)
    elseif(PLAT STREQUAL "_OPENMP")
      add_dependencies(superlu_mt_OPENMP metis GKlib)
    endif()
  endif()

  set(SLU_DRIVER "${SUPERLU_MT_DIR}/SRC/c_fortran_pdgssv.c")
else()
  message(STATUS "Will link against regular SuperLU.")

  if(NOT EXISTS "${SUPERLU_DIR}/CMakeLists.txt")
    message(
      FATAL_ERROR
      "The SuperLU submodule was not downloaded. You have to do it manually!"
    )
  endif()

  include_directories("${SUPERLU_DIR}/SRC")
  add_subdirectory(${SUPERLU_DIR})

  if(TPL_ENABLE_METISLIB)
    add_dependencies(superlu metis GKlib)
  endif()

  set(SLU_DRIVER "${SUPERLU_DIR}/FORTRAN/c_fortran_dgssv.c")
endif()

# f2c stuff
set(F2C_DIR "${PROJECT_SOURCE_DIR}/f2c")
set(F2C_INCLUDE_DIR "${F2C_DIR}/include")
set(F2C_FN "${F2C_DIR}/libf2c.zip")
set(F2C_URL "https://www.netlib.org/f2c/libf2c.zip")

# download f2c
if(NOT EXISTS ${F2C_DIR})
  message(STATUS "Downloading libf2c source from ${F2C_URL}...")
  make_directory("${F2C_DIR}")
  file(DOWNLOAD ${F2C_URL} ${F2C_FN} TIMEOUT 60 STATUS DOWNLOAD_STATUS)

  # Check if download was successful.
  list(GET DOWNLOAD_STATUS 0 STATUS_CODE)
  list(GET DOWNLOAD_STATUS 1 ERROR_MESSAGE)

  if(${STATUS_CODE} EQUAL 0)
    message(STATUS "Done downloading libf2c.")
  else()
    # Exit CMake if the download failed, printing the error message.
    file(REMOVE_RECURSE ${F2C_DIR})
    message(FATAL_ERROR "Error downloading libf2c: ${ERROR_MESSAGE}")
  endif()
endif()

# extract libf2c source
file(ARCHIVE_EXTRACT INPUT ${F2C_FN} DESTINATION ${F2C_DIR})

# prepare libf2c header files
file(GLOB_RECURSE F2C_PREHEADERS "${F2C_DIR}/*.h0")

foreach(H0 ${F2C_PREHEADERS})
  string(REGEX REPLACE "[.]h0$" ".h" H0_R ${H0})
  file(RENAME "${H0}" "${H0_R}")
  file(COPY "${H0_R}" DESTINATION "${F2C_INCLUDE_DIR}")
endforeach()

# get a load of this: f2c generates its own "arith.h" on the fly
# so we gotta compile arithchk and run it
set(F2C_ARITHCHK_SRC "${F2C_DIR}/arithchk.c")
set(F2C_ARITHCHK_BIN "${F2C_DIR}/arithchk")

if(WIN32)
  set(F2C_ARITHCHK_BIN "${F2C_ARITHCHK_BIN}.exe")
endif()

set(F2C_ARITH_H "${F2C_INCLUDE_DIR}/arith.h")
set_source_files_properties(
  ${F2C_ARITHCHK_SRC} PROPERTIES COMPILE_FLAGS "-DNO_LONG_LONG -DNO_FPINIT"
)
add_executable(arithchk ${F2C_ARITHCHK_SRC})
target_link_libraries(arithchk m)
set_target_properties(
  arithchk PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${F2C_DIR}
)
add_custom_command(
  OUTPUT ${F2C_ARITH_H}
  COMMAND ${F2C_ARITHCHK_BIN} > ${F2C_ARITH_H}
  DEPENDS ${F2C_ARITHCHK_BIN}
)

# add libf2c to the compilation procedures
include_directories(${F2C_INCLUDE_DIR})
file(GLOB_RECURSE F2C_CFILES "${F2C_DIR}/*.c")
add_definitions(-DINTEGER_STAR_8)
add_library(f2c ${F2C_CFILES} ${F2C_ARITH_H})

# add some extra win32 flags for libf2c
if(WIN32)
  add_definitions(-DUSE_CLOCK -DMSDOS)
endif()

# set some extra vars for MSYS builds to make the binary portable
if(WIN32)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")
endif()

# collect all fortran source files
file(GLOB_RECURSE ALL_FORTRAN_FILES
  "${CMAKE_SOURCE_DIR}/*.f"
  "${CMAKE_SOURCE_DIR}/*.F"
  "${CMAKE_SOURCE_DIR}/*.f90"
  "${CMAKE_SOURCE_DIR}/*.F90"
  "${CMAKE_SOURCE_DIR}/*.f95"
  "${CMAKE_SOURCE_DIR}/*.F95"
  "${CMAKE_SOURCE_DIR}/*.f03"
  "${CMAKE_SOURCE_DIR}/*.F03"
)

# same BLAS-finding subroutine as SuperLU
if(NOT enable_internal_blaslib)
  if(TPL_BLAS_LIBRARIES)
    set(BLAS_FOUND TRUE)
  else()
    find_package(BLAS)

    if(BLAS_FOUND)
      set(TPL_BLAS_LIBRARIES "${BLAS_LIBRARIES}" CACHE FILEPATH
        "Set from FindBLAS.cmake BLAS_LIBRARIES." FORCE)
    endif()
  endif()
endif()

if(NOT BLAS_FOUND)
  message(STATUS "BLAS not found, building local BLAS")
  list(
    APPEND blas_fns dgemm dgemv dlamch dlanst dscal dsteqr dsterf dswap dtrsm
    dtrtri ilaenv lsame xerbla
  )

  foreach(fname IN LISTS blas_fns)
    # message(STATUS "Checking for BLAS subr ${fname}")
    set(CMAKE_REQUIRED_QUIET TRUE)
    check_function_exists("${fname}" BLAS_FN_EXISTS)
    unset(CMAKE_REQUIRED_QUIET)
    set(BLAS_FN_EXISTS CACHE "1" STRING)

    if(NOT BLAS_FN_EXISTS)
      # message(STATUS "BLAS subr ${fname} not found")
      string(TOUPPER ${fname} fname_upper)
      list(APPEND missing_blas_src "${BLAS_SOURCE_DIR}/${fname_upper}.f")
      list(APPEND missing_blas_fns ${fname})
    endif()

    unset(BLAS_FN_EXISTS CACHE)
  endforeach()

  # if any subroutines have bene found, create an inner blas library
  list(LENGTH missing_blas_fns MISSING_FNS_TOTAL)

  if(MISSING_FNS_TOTAL GREATER 0)
    if(MISSING_FNS_TOTAL GREATER 1)
      message(
        STATUS
        "BLAS subrs (${missing_blas_fns}) are absent and will be built locally."
      )
    else()
      message(
        STATUS
        "BLAS subr ${missing_blas_fns} is absent and will be built locally."
      )
    endif()
  endif()
else()
  message(STATUS "Using system BLAS (${TPL_BLAS_LIBRARIES})")
endif()

# prepare the main executable, linked against the specifics and the m
# it appears utils used to be a module, but that is no longer the case?
# file(GLOB UTIL_FILES "${CMAKE_SOURCE_DIR}/UTIL/*.f*")
# file(GLOB MAIN_FILES "${CMAKE_SOURCE_DIR}/MAIN/*.[fF]*")
add_executable(
  mystran
  ${ALL_FORTRAN_FILES}
  ${SLU_DRIVER}
  ${missing_blas_src}
)

# determine which superlu variant to link against
if(USE_SUPERLU_MT)
  if(PLAT STREQUAL "_PTHREAD")
    message(WARNING "We recommend using OpenMP, not pthread!")
    target_link_libraries(mystran superlu_mt_PTHREAD f2c)
  elseif(PLAT STREQUAL "_OPENMP")
    target_link_libraries(mystran superlu_mt_OPENMP f2c)
  endif()

  target_compile_definitions(mystran PRIVATE USE_SUPERLU_MT)

  if(TPL_ENABLE_METISLIB)
    target_compile_definitions(mystran PRIVATE HAVE_METIS)
  endif()
else()
  target_link_libraries(mystran superlu f2c)
endif()

set_target_properties(
  mystran PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
)

if(WIN32)
  target_sources(mystran PRIVATE mystran.rc)
endif()

# -----------------------------------------------------------------------
# Static library tracking
#
# For every library that is statically linked into mystran we:
# 1. Record a compile-time preprocessor symbol (_STATIC_<NAME>)
# 2. Record a human-readable filename for the status message
#
# After populating the lists we apply the definitions, print a summary,
# and warn for any symbol not referenced in PRINT_BUILD_INFO.F90.
# -----------------------------------------------------------------------
set(_MYSTRAN_STATIC_DEFS "")
set(_MYSTRAN_STATIC_NAMES "")

# f2c is always built and linked as a static archive
list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_F2C)
list(APPEND _MYSTRAN_STATIC_NAMES libf2c.a)

if(USE_SUPERLU_MT)
  list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_SUPERLU_MT)

  if(WIN32)
    list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_WINPTHREAD)
    list(APPEND _MYSTRAN_STATIC_NAMES libwinpthread.a)
  else()
    list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_PTHREAD)

    # no name because it's part of glibc
  endif()

  if(PLAT STREQUAL "_OPENMP")
    list(APPEND _MYSTRAN_STATIC_NAMES libsuperlu_mt_OPENMP.a)

    if(WIN32)
      # On Windows the gomp runtime is force-linked as a static archive
      list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_GOMP)
      list(APPEND _MYSTRAN_STATIC_NAMES libgomp.a)
    endif()
  elseif(PLAT STREQUAL "_PTHREAD")
    list(APPEND _MYSTRAN_STATIC_NAMES libsuperlu_mt_PTHREAD.a)
  endif()
else()
  list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_SUPERLU)
  list(APPEND _MYSTRAN_STATIC_NAMES libsuperlu.a)
endif()

if(TPL_ENABLE_METISLIB)
  list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_METIS)
  list(APPEND _MYSTRAN_STATIC_NAMES libmetis.a)
  list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_GKLIB)
  list(APPEND _MYSTRAN_STATIC_NAMES libGKlib.a)

  if(NOT WIN32)
    list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_REGEX)

    # only on linux because on windows we build embedded gkegex.
    # no name because it's part of glibc.
  endif()
endif()

# Local BLAS functions are compiled directly into the binary when the
# system does not provide them; treat them as a bundled static object
if(MISSING_FNS_TOTAL GREATER 0)
  list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_BLAS)
  list(APPEND _MYSTRAN_STATIC_NAMES libblas_local.a)
endif()

# On Windows the -static linker flag pulls in Fortran and C runtimes too
if(WIN32)
  list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_GFORTRAN)
  list(APPEND _MYSTRAN_STATIC_NAMES libgfortran.a)
  list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_GCC)
  list(APPEND _MYSTRAN_STATIC_NAMES libgcc.a)
endif()

# Apply preprocessor definitions to the mystran compilation
foreach(_def IN LISTS _MYSTRAN_STATIC_DEFS)
  target_compile_definitions(mystran PRIVATE ${_def})
endforeach()

# Print a human-readable summary at configure time
list(JOIN _MYSTRAN_STATIC_NAMES ", " _STATIC_NAMES_STR)
message(STATUS "Static libraries linked in: ${_STATIC_NAMES_STR}")

# force some libraries
list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_LAPACK)
list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_ARPACK)

# -----------------------------------------------------------------------
# License file associations
#
# Map each _STATIC_<NAME> suffix to a license file path (relative to
# PROJECT_SOURCE_DIR).  CMake will:
# - define  _STATIC_<NAME>_LICENSE  for the Fortran preprocessor
# - generate a SUBROUTINE PRINT_LICENSE_<NAME>(IUNIT) that writes the
# license text line-by-line to whatever unit number is passed in
#
# Add new entries here whenever a new statically-linked library is added.
# -----------------------------------------------------------------------
set(_MYSTRAN_LICENSE_MAP
  "F2C:f2c/Notice:libf2c:Fortran-to-C runtime support library"
  "SUPERLU:superlu/License.txt:SuperLU:sparse solver"
  "SUPERLU_MT:superlu_mt/License.txt:SuperLU_MT:multi-threaded sparse solver"
  "METIS:metis/METIS/LICENSE:METIS:graph partitioner for SuperLU"
  "GKLIB:metis/GKlib/LICENSE.txt:GKlib:utility library for METIS"
  "BLAS:extra_licenses/blas-lapack.txt:BLAS:embedded subset of reference BLAS"
  "LAPACK:extra_licenses/blas-lapack.txt:LAPACK:linear algebra library"
  "ARPACK:extra_licenses/arpack.txt:ARPACK:Lanczos eigensolver implementation"
  "GOMP:extra_licenses/gpl-3.0-with-gcc-exception.txt:libgomp:GNU OpenMP runtime"
  "GFORTRAN:extra_licenses/gpl-3.0-with-gcc-exception.txt:libgfortran:GNU Fortran runtime"
  "GCC:extra_licenses/gpl-3.0-with-gcc-exception.txt:libgcc:GNU C runtime"
  "WINPTHREAD:extra_licenses/libwinpthread.txt:libwinpthread:Windows POSIX threads, used by SuperLU_MT"
  "PTHREAD:extra_licenses/lgpl-2.1.txt:libpthread:POSIX threads, used by SuperLU_MT"
  "REGEX:extra_licenses/lgpl-2.1.txt:libregex:POSIX regexes, used by GKlib"
)

# Generate STATIC_LICENSES_GENERATED.F90 in the build directory
set(_GEN_LIC_FILE "${CMAKE_SOURCE_DIR}/INCLUDE/STATIC_LICENSES_GENERATED.F90.inc")
set(_GEN_LIC_CONTENT "! AUTO-GENERATED by CMakeLists.txt -- DO NOT EDIT\n")
set(_MYSTRAN_LICENSED_NAMES "")

# ---- MYSTRAN's own license (always emitted, no #ifdef guard) ----
set(_MYSTRAN_OWN_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE.txt")

if(EXISTS "${_MYSTRAN_OWN_LICENSE}")
  file(STRINGS "${_MYSTRAN_OWN_LICENSE}" _mlic_lines NEWLINE_CONSUME)
  string(REPLACE "\n" ";" _mlic_lines "${_mlic_lines}")
  string(APPEND _GEN_LIC_CONTENT
    "\nSUBROUTINE PRINT_LICENSE_MYSTRAN(IUNIT)\n"
    "  IMPLICIT NONE\n"
    "  INTEGER, INTENT(IN) :: IUNIT\n"
  )

  foreach(_line IN LISTS _mlic_lines)
    string(REPLACE "'" "''" _line_esc "${_line}")
    string(REPLACE "\\" "\\\\" _line_esc "${_line_esc}")
    string(APPEND _GEN_LIC_CONTENT
      "  WRITE(IUNIT,'(A)') '${_line_esc}'\n"
    )
  endforeach()

  string(APPEND _GEN_LIC_CONTENT
    "END SUBROUTINE PRINT_LICENSE_MYSTRAN\n"
  )
else()
  message(WARNING "MYSTRAN license file not found: ${_MYSTRAN_OWN_LICENSE}")
endif()

foreach(_entry IN LISTS _MYSTRAN_LICENSE_MAP)
  string(REPLACE ":" ";" _pair "${_entry}")
  list(GET _pair 0 _lname)
  list(GET _pair 1 _lpath)
  list(LENGTH _pair _pair_len)

  if(_pair_len GREATER 2)
    list(GET _pair 2 _llibname)
  else()
    set(_llibname "${_lname}")
  endif()

  if(_pair_len GREATER 3)
    list(GET _pair 3 _ldesc)
  else()
    set(_ldesc "")
  endif()

  set(_MYSTRAN_LIBNAME_${_lname} "${_llibname}")
  set(_MYSTRAN_DESC_${_lname} "${_ldesc}")
  set(_lfile "${PROJECT_SOURCE_DIR}/${_lpath}")

  if(EXISTS "${_lfile}")
    # Read the license file line by line
    file(STRINGS "${_lfile}" _lic_lines NEWLINE_CONSUME)

    # file(STRINGS) may return the whole file as one item when using
    # NEWLINE_CONSUME; split on actual newlines ourselves
    string(REPLACE "\n" ";" _lic_lines "${_lic_lines}")

    string(APPEND _GEN_LIC_CONTENT
      "\nSUBROUTINE PRINT_LICENSE_${_lname}(IUNIT)\n"
      "  IMPLICIT NONE\n"
      "  INTEGER, INTENT(IN) :: IUNIT\n"
    )

    foreach(_line IN LISTS _lic_lines)
      # Escape Fortran string delimiters: ' -> ''
      string(REPLACE "'" "''" _line_esc "${_line}")

      # Escape backslashes that CMake would otherwise swallow
      string(REPLACE "\\" "\\\\" _line_esc "${_line_esc}")
      string(APPEND _GEN_LIC_CONTENT
        "  WRITE(IUNIT,'(A)') '${_line_esc}'\n"
      )
    endforeach()

    string(APPEND _GEN_LIC_CONTENT
      "END SUBROUTINE PRINT_LICENSE_${_lname}\n"
    )

    list(APPEND _MYSTRAN_LICENSED_NAMES "${_lname}")
  else()
    message(WARNING
      "License file for _STATIC_${_lname} not found: ${_lfile}\n"
      "  '_STATIC_${_lname}_LICENSE' will NOT be defined."
    )
  endif()
endforeach()

# Append PRINT_ALL_LICENSES: MYSTRAN first, then delegates to each PRINT_LICENSE_* under its #ifdef guard
string(APPEND _GEN_LIC_CONTENT
  "\nSUBROUTINE PRINT_ALL_LICENSES(IUNIT)\n"
  "  IMPLICIT NONE\n"
  "  INTEGER, INTENT(IN) :: IUNIT\n"

  # "  WRITE(IUNIT,'(A)') '========================================================'\n"
  # "  WRITE(IUNIT,'(A)') ' MYSTRAN '\n"
  "  WRITE(IUNIT,'(A)') '========================================================'\n"
  "  CALL PRINT_LICENSE_MYSTRAN(IUNIT)\n"
  "  WRITE(IUNIT,'(A)') ' '\n"
)

foreach(_lname IN LISTS _MYSTRAN_LICENSED_NAMES)
  set(_llibname "${_MYSTRAN_LIBNAME_${_lname}}")
  set(_ldesc "${_MYSTRAN_DESC_${_lname}}")
  string(APPEND _GEN_LIC_CONTENT
    "#ifdef _STATIC_${_lname}\n"
    "  WRITE(IUNIT,'(A)') '========================================================'\n"
    "  WRITE(IUNIT,'(A)') ' ${_llibname} -- ${_ldesc}'\n"
    "  WRITE(IUNIT,'(A)') '========================================================'\n"
    "  CALL PRINT_LICENSE_${_lname}(IUNIT)\n"
    "  WRITE(IUNIT,'(A)') ' '\n"
    "#endif\n"
  )
endforeach()

string(APPEND _GEN_LIC_CONTENT
  "END SUBROUTINE PRINT_ALL_LICENSES\n"
)

file(WRITE "${_GEN_LIC_FILE}" "${_GEN_LIC_CONTENT}")
message(STATUS "Generated license subroutines: ${_GEN_LIC_FILE}")

# Warn if any statically-linked library lacks a license mapping.
# Runtime and utility libs with no separate OSS license file are exempt.
set(_LICENSE_EXEMPT_LIBS)

foreach(_def IN LISTS _MYSTRAN_STATIC_DEFS)
  string(REPLACE "_STATIC_" "" _lname "${_def}")

  if(NOT _lname IN_LIST _MYSTRAN_LICENSED_NAMES)
    if(NOT _lname IN_LIST _LICENSE_EXEMPT_LIBS)
      message(WARNING
        "No license file is mapped for '${_def}'. "
        "Add an entry to _MYSTRAN_LICENSE_MAP in CMakeLists.txt."
      )
    endif()
  endif()
endforeach()

if(CMAKE_COMPILER_IS_GNUCC)
  set(CMAKE_Fortran_FLAGS "-Wall -Wextra \
    -Wno-unused-variable -Wno-unused-label -Wno-unused-parameter -Wno-tabs \
    -Wno-compare-reals -Wno-character-truncation -Wno-unused-dummy-argument \
    -Wmaybe-uninitialized -Wrealloc-lhs -fbounds-check -fcheck=all \
    -ffree-line-length-none")
  set(CMAKE_Fortran_FLAGS_DEBUG
    "-g -O0 -fcheck=all -fbacktrace -fbounds-check \
    -fno-inline -fno-ipa-sra -fno-ipa-cp -fno-optimize-sibling-calls")
  set(CMAKE_Fortran_FLAGS_DETERMINISTIC
    "-O0 -g \
     -fno-fast-math \
     -ffp-contract=off \
     -fno-unsafe-math-optimizations \
     -fno-associative-math \
     -fno-reciprocal-math \
     -frounding-math"
  )
  set(CMAKE_C_FLAGS_DETERMINISTIC
    "-O0 -g \
     -fno-fast-math \
     -ffp-contract=off \
     -frounding-math"
  )
  set(CMAKE_CXX_FLAGS_DETERMINISTIC
    "${CMAKE_C_FLAGS_DETERMINISTIC}"
  )
  set(CMAKE_PROFILING_FLAGS "-O2 -g -fno-inline -fno-omit-frame-pointer -fno-inline-functions")
  set(CMAKE_C_FLAGS_PROFILING "${CMAKE_PROFILING_FLAGS}")
  set(CMAKE_CXX_FLAGS_PROFILING "${CMAKE_PROFILING_FLAGS}")
  set(CMAKE_Fortran_FLAGS_PROFILING "${CMAKE_PROFILING_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS_PROFILING "")
  set(CMAKE_SHARED_LINKER_FLAGS_PROFILING "")
endif()

# -----------------------------------------------------------------------
# Build metadata generated file
#
# Generates BUILD_INFO_GENERATED.F90.inc with:
# PRINT_BUILD_CONSTANTS(IUNIT) - timestamp, compiler, flags, solver info
# PRINT_STATIC_LIB_LIST(IUNIT) - names of all statically linked libraries
# -----------------------------------------------------------------------

# Build timestamp captured at configure time (local clock)
string(TIMESTAMP _BUILD_TIMESTAMP "%Y-%m-%d %H:%M:%S UTC" UTC)

set(_BUILD_TYPE_DISPLAY "${CMAKE_BUILD_TYPE}")

if(NOT _BUILD_TYPE_DISPLAY)
  set(_BUILD_TYPE_DISPLAY "(not specified)")
endif()

# Effective Fortran flags = base flags + build-type-specific flags
string(TOUPPER "${CMAKE_BUILD_TYPE}" _BT_UPPER)
set(_FC_FLAGS "${CMAKE_Fortran_FLAGS}")

if(_BT_UPPER AND DEFINED CMAKE_Fortran_FLAGS_${_BT_UPPER})
  string(APPEND _FC_FLAGS " ${CMAKE_Fortran_FLAGS_${_BT_UPPER}}")
endif()

string(REPLACE "\n" " " _FC_FLAGS "${_FC_FLAGS}")
string(REGEX REPLACE " +" " " _FC_FLAGS "${_FC_FLAGS}")
string(STRIP "${_FC_FLAGS}" _FC_FLAGS)
string(REPLACE "'" "''" _FC_FLAGS "${_FC_FLAGS}")

# Effective C flags
set(_CC_FLAGS "${CMAKE_C_FLAGS}")

if(_BT_UPPER AND DEFINED CMAKE_C_FLAGS_${_BT_UPPER})
  string(APPEND _CC_FLAGS " ${CMAKE_C_FLAGS_${_BT_UPPER}}")
endif()

string(REPLACE "\n" " " _CC_FLAGS "${_CC_FLAGS}")
string(REGEX REPLACE " +" " " _CC_FLAGS "${_CC_FLAGS}")
string(STRIP "${_CC_FLAGS}" _CC_FLAGS)
string(REPLACE "'" "''" _CC_FLAGS "${_CC_FLAGS}")

# Compiler identity strings
set(_FC_DISPLAY "${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
set(_CC_DISPLAY "${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")
string(REPLACE "'" "''" _FC_PATH "${CMAKE_Fortran_COMPILER}")
string(REPLACE "'" "''" _CC_PATH "${CMAKE_C_COMPILER}")

# Sparse solver variant
if(USE_SUPERLU_MT)
  if(PLAT STREQUAL "_OPENMP")
    set(_SLU_VARIANT "SuperLU_MT (OpenMP)")
  else()
    set(_SLU_VARIANT "SuperLU_MT (pthreads)")
  endif()
else()
  set(_SLU_VARIANT "SuperLU (single-threaded)")
endif()

# BLAS source
if(BLAS_FOUND)
  set(_BLAS_INFO "System BLAS: ${TPL_BLAS_LIBRARIES}")
  string(REPLACE "'" "''" _BLAS_INFO "${_BLAS_INFO}")
elseif(MISSING_FNS_TOTAL GREATER 0)
  set(_BLAS_INFO "CBLAS + built-in subset (${MISSING_FNS_TOTAL} routine(s) compiled in)")
else()
  set(_BLAS_INFO "System BLAS")
endif()

# METIS
if(TPL_ENABLE_METISLIB)
  set(_METIS_INFO "Enabled")
else()
  set(_METIS_INFO "Disabled")
endif()

set(_GEN_BLD_FILE "${CMAKE_SOURCE_DIR}/INCLUDE/BUILD_INFO_GENERATED.F90.inc")
set(_GEN_BLD_CONTENT "! AUTO-GENERATED by CMakeLists.txt -- DO NOT EDIT\n")

# PRINT_BUILD_CONSTANTS subroutine
string(APPEND _GEN_BLD_CONTENT
  "\nSUBROUTINE PRINT_BUILD_CONSTANTS(IUNIT)\n"
  "  IMPLICIT NONE\n"
  "  INTEGER, INTENT(IN) :: IUNIT\n"
  "  WRITE(IUNIT,'(A)') ' '\n"
  "  WRITE(IUNIT,'(A)') ' Build information:'\n"
  "  WRITE(IUNIT,'(A)') ' ------------------'\n"
  "  WRITE(IUNIT,'(A)') '  Build date (UTC)    : ${_BUILD_TIMESTAMP}'\n"
  "  WRITE(IUNIT,'(A)') '  CMake build type    : ${_BUILD_TYPE_DISPLAY}'\n"
  "  WRITE(IUNIT,'(A)') ' '\n"
  "  WRITE(IUNIT,'(A)') '  Fortran compiler    : ${_FC_DISPLAY}'\n"
  "  WRITE(IUNIT,'(A)') '  Fortran path        : ${_FC_PATH}'\n"
  "  WRITE(IUNIT,'(A)') '  Fortran flags       : ${_FC_FLAGS}'\n"
  "  WRITE(IUNIT,'(A)') ' '\n"
  "  WRITE(IUNIT,'(A)') '  C compiler          : ${_CC_DISPLAY}'\n"
  "  WRITE(IUNIT,'(A)') '  C compiler path     : ${_CC_PATH}'\n"
  "  WRITE(IUNIT,'(A)') '  C flags             : ${_CC_FLAGS}'\n"
  "  WRITE(IUNIT,'(A)') ' '\n"
  "  WRITE(IUNIT,'(A)') '  Sparse solver       : ${_SLU_VARIANT}'\n"
  "  WRITE(IUNIT,'(A)') '  BLAS                : ${_BLAS_INFO}'\n"
  "  WRITE(IUNIT,'(A)') '  METIS partitioner   : ${_METIS_INFO}'\n"
  "END SUBROUTINE PRINT_BUILD_CONSTANTS\n"
)

# PRINT_STATIC_LIB_LIST subroutine
string(APPEND _GEN_BLD_CONTENT
  "\nSUBROUTINE PRINT_STATIC_LIB_LIST(IUNIT)\n"
  "  IMPLICIT NONE\n"
  "  INTEGER, INTENT(IN) :: IUNIT\n"
  "  WRITE(IUNIT,'(A)') ' Static libraries compiled into this MYSTRAN binary:'\n"
  "  WRITE(IUNIT,'(A)') ' ----------------------------------------------------'\n"
)

foreach(_def IN LISTS _MYSTRAN_STATIC_DEFS)
  string(REPLACE "_STATIC_" "" _suffix "${_def}")
  set(_libname "${_MYSTRAN_LIBNAME_${_suffix}}")
  set(_desc "${_MYSTRAN_DESC_${_suffix}}")

  if(NOT _libname)
    set(_libname "${_def}")
    set(_desc "(no description)")
  endif()

  string(REPLACE "'" "''" _libname_esc "${_libname}")
  string(REPLACE "'" "''" _desc_esc "${_desc}")

  string(APPEND _GEN_BLD_CONTENT
    "  WRITE(IUNIT,'(3X,A,T20,A)') '${_libname_esc}', '(${_desc_esc})'\n"
  )
endforeach()

string(APPEND _GEN_BLD_CONTENT
  "  WRITE(IUNIT,'(A)') ' '\n"
  "END SUBROUTINE PRINT_STATIC_LIB_LIST\n"
)

file(WRITE "${_GEN_BLD_FILE}" "${_GEN_BLD_CONTENT}")
message(STATUS "Generated build-info subroutines: ${_GEN_BLD_FILE}")

# issue a couple messages about compilation
include(ProcessorCount)
ProcessorCount(NCPU)
message(STATUS "You can now compile MYSTRAN with `cmake --build .`")

if(NOT NCPU EQUAL 0)
  message(STATUS "Compile faster by passing the -j${NCPU} flag.")
endif()
