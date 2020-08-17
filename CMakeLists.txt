# CMakeTutorial

# This minimum version is mostly set in order to get a newer version
# of the FindMPI check.  Note that you can easily install a newer cmake version
# using conda or pip.
cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

foreach(policy
    CMP0048
    CMP0074
    CMP0077
  )
  if(POLICY ${policy})
    cmake_policy(SET ${policy} NEW)
  endif()
endforeach()

project(CMakeTutorial VERSION 1.0.0 LANGUAGES C)

# set cflags

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -ffast-math -Wall -O3 -ftree-vectorize -std=gnu99")

# External packages

find_package(GSL REQUIRED)


# Set git version header file
message(STATUS "Resolving GIT Version")

find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    WORKING_DIRECTORY "${local_dir}"
    OUTPUT_VARIABLE _build_version
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  message( STATUS "GIT hash: ${_build_version}")
else()
  message(STATUS "GIT not found")
endif()


# The version number.
set (GITVERSION "\"${_build_version}\"")

configure_file (
  "${PROJECT_SOURCE_DIR}/gitversion.h.in"
  "${PROJECT_SOURCE_DIR}/gbmcmc/src/gitversion.h"
  )

install (FILES "${PROJECT_SOURCE_DIR}/gbmcmc/src/gitversion.h"
         DESTINATION include)

message(STATUS "Resolving GIT Version")

# The version number.
set (GITVERSION "\"${_build_version}\"")

# Internal products

add_subdirectory(tools)
add_subdirectory(gbmcmc)
add_subdirectory(gbfisher)
