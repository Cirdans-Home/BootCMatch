cmake_minimum_required(VERSION 3.10)

project(BootCMatch VERSION 0.9
									 DESCRIPTION "Bootstrap AMG based on Compatible weighted Matching"
									 HOMEPAGE_URL "https://psctoolkit.github.io/"
									 LANGUAGES C CXX Fortran)

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

option(INSTALL_HSL "Use HSL Library" OFF)
if(INSTALL_HSL)
	find_package(HSL REQUIRED)
	if(HSL_FOUND)
		message(STATUS "Found HSL Library.")
	endif(HSL_FOUND)
endif(INSTALL_HSL)

option(INSTALL_SPRAL "Use SPRAL Library" OFF)
if(INSTALL_SPRAL)
	find_package(SPRAL REQUIRED)
	if(SPRAL_FOUND)
		message(STATUS "Found SPRAL Library.")
	endif(SPRAL_FOUND)
endif(INSTALL_SPRAL)

find_package(GFortran REQUIRED)

#-----------------------------------------------------------------------------#
# We look for the auxiliary libraries
#-----------------------------------------------------------------------------#
option(USE_SUPERLU "Link with SuperLU (Version 4 or 5)" OFF)
find_package(SUPERLU)
if(SUPERLU_FOUND)
	message(STATUS "Found SUPERLU Library.")
endif(SUPERLU_FOUND)

# TwoThird Approximate matching with λ parameter choice
option(USE_AMGMatch "Use matchLambdaTwoThirdeps and matchLambdaOpt Algorithms" OFF)
if(USE_AMGMatch)
	add_subdirectory(AUXIL/Matchbox)	# Matchbox Library
	# We look for an external installation of the Lemon Libary, it is required
	find_package(LEMON REQUIRED)
	if(LEMON_FOUND)
		message(STATUS "Found LEMON Library.")
	endif(LEMON_FOUND)
	add_subdirectory(AUXIL/AMGMatch)
	# Set Variables for compilation
	set(AUXIL_FLAGS HAVE_AMGMATCH)
	set(AUXIL_LIBS matchbox amgmatch)
endif(USE_AMGMatch)

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
