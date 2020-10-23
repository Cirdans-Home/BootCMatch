# Compiling the tests

set(testprogram "testbootstrap" "testkbootsolvers" "testsolving")

foreach(testname IN LISTS testprogram)
	add_executable(${testname} 		${testname}
														 		ioutil.c
														 		ioutil.h)
	target_include_directories( ${testname} PUBLIC  ${CMAKE_SOURCE_DIR}/SRC/utilities
																				 	${CMAKE_SOURCE_DIR}/SRC/MATVEC
																					${CMAKE_SOURCE_DIR}/SRC/AMG
																					${CMAKE_SOURCE_DIR}/SRC/BOOTAMG
																					${CMAKE_SOURCE_DIR}/SRC/LINSOLVERS
																					${SUPERLU_INCLUDE_DIR}
																					${HSL_INCDIR}
																					${SPRAL_INCDIR})
	target_link_libraries( ${testname} PUBLIC BCM)

	# Install the executables
	install(TARGETS ${testname}
		RUNTIME
			DESTINATION samples/RUNS
			COMPONENT Runtime)
endforeach()

# Directory with the scripts
install(DIRECTORY RUNS DESTINATION samples)
