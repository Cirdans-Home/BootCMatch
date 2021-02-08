# Compiling the tests

set(testprogram "testbootstrap" "testkbootsolvers" "testsolving"
								"testkbootsolvers3dlap")

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

# Install Makefile from template and source files
string(REPLACE ";" " " testprogram "${testprogram}")
configure_file(
	"${CMAKE_CURRENT_SOURCE_DIR}/Makefile.in"
	"${CMAKE_INSTALL_PREFIX}/samples/Makefile"
	@ONLY
)
install(FILES ioutil.c testbootstrap.c testkbootsolvers.c
							testkbootsolvers3dlap.c testsolving.c ioutil.h
	DESTINATION samples
	COMPONENT Runtime)
