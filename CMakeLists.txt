cmake_minimum_required(VERSION 3.10)

project(BootCMatch VERSION 0.9
									 DESCRIPTION "Bootstrap AMG based on Compatible weighted Matching"
									 HOMEPAGE_URL "https://psctoolkit.github.io/"
									 LANGUAGES C Fortran)

# Directory with CMAKE modules, used to search for libraries
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

#-----------------------------------------------------------------------------#
# We look for required libraries
#-----------------------------------------------------------------------------#
find_package(BLAS REQUIRED)
if(BLAS_FOUND)
	message(STATUS "Found BLAS Library.")
endif(BLAS_FOUND)

find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
	message(STATUS "Found LAPACK Library.")
endif(LAPACK_FOUND)

find_package(HSL REQUIRED)
if(HSL_FOUND)
	message(STATUS "Found HSL Library.")
endif(HSL_FOUND)

find_package(SPRAL REQUIRED)
if(SPRAL_FOUND)
	message(STATUS "Found SPRAL Library.")
endif(SPRAL_FOUND)

#-----------------------------------------------------------------------------#
# We look for the auxiliary libraries
#-----------------------------------------------------------------------------#
option(USE_SUPERLU "Link with SuperLU (Version 4 or 5)" OFF)
find_package(SUPERLU)
if(SUPERLU_FOUND)
	message(STATUS "Found SUPERLU Library.")
endif(SUPERLU_FOUND)

#-----------------------------------------------------------------------------#
# Start The compilation of the Library
#-----------------------------------------------------------------------------#
add_subdirectory(SRC)

#-----------------------------------------------------------------------------#
# Do I need to compile the tests?
#-----------------------------------------------------------------------------#
option(COMPILE_TEST "Compile the tests?" ON)
if(COMPILE_TEST)
	add_subdirectory(Test)
endif(COMPILE_TEST)

#-----------------------------------------------------------------------------#
# Settings for the installation
#-----------------------------------------------------------------------------#
install(TARGETS BCM
        LIBRARY
          DESTINATION lib
          COMPONENT Libraries
          NAMELINK_COMPONENT Development
        PUBLIC_HEADER
          DESTINATION include
          COMPONENT Development
       )
