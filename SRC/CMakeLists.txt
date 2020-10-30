# We build here the libBCM library

add_subdirectory(utilities)
add_subdirectory(MATVEC)
add_subdirectory(AMG)
add_subdirectory(BOOTAMG)
add_subdirectory(LINSOLVERS)

add_library(BCM STATIC ${utilities_sources}
											 ${matvec_sources}
											 ${amg_sources}
											 ${bootamg_sources}
											 ${linsolvers_sources})

target_include_directories( BCM PUBLIC  ${CMAKE_SOURCE_DIR}/SRC/utilities
																			 	${CMAKE_SOURCE_DIR}/SRC/MATVEC
																				${CMAKE_SOURCE_DIR}/SRC/AMG
																				${CMAKE_SOURCE_DIR}/SRC/BOOTAMG
																				${CMAKE_SOURCE_DIR}/SRC/LINSOLVERS
																				${SUPERLU_INCLUDE_DIR}
																				${HSL_INCDIR}
																				${SPRAL_INCDIR})

target_link_libraries( BCM PUBLIC -lm ${BLAS_LIBRARIES}
																			${LAPACK_LIBRARIES}
																			${SPRAL_LIBS}
																			${HSL_LIBS}
																			${SUPERLU_LIBRARIES}
																			${AUXIL_LIBS}
																			${GFORTRAN_LIBRARY})

target_link_directories( BCM PUBLIC ${SPRAL_LIBDIR} ${HSL_LIBDIR} ${LEMON_LIBDIR})

target_compile_definitions( BCM PUBLIC ${HSL_FLAGS}
																			 ${SPRAL_FLAGS}
																			 ${SUPERLU_FLAGS}
																			 ${AUXIL_FLAGS})

set_property(TARGET BCM PROPERTY
		PUBLIC_HEADER ${amg_header}
									${bootamg_header}
									${linsolvers_header}
									${matvec_header}
									${utilities_header})

if(BUILD_SHARED_LIBS)
	add_library(BCMshared SHARED ${utilities_sources}
												 			 ${matvec_sources}
												 		 	 ${amg_sources}
												 		 	 ${bootamg_sources}
												 		 	 ${linsolvers_sources})

	target_include_directories( BCMshared PUBLIC  ${CMAKE_SOURCE_DIR}/SRC/utilities
																				 	${CMAKE_SOURCE_DIR}/SRC/MATVEC
																					${CMAKE_SOURCE_DIR}/SRC/AMG
																					${CMAKE_SOURCE_DIR}/SRC/BOOTAMG
																					${CMAKE_SOURCE_DIR}/SRC/LINSOLVERS
																					${SUPERLU_INCLUDE_DIR}
																					${HSL_INCDIR}
																					${SPRAL_INCDIR})

	target_link_libraries( BCMshared PUBLIC -lm ${BLAS_LIBRARIES}
																				${LAPACK_LIBRARIES}
																				${SPRAL_LIBS}
																				${HSL_LIBS}
																				${SUPERLU_LIBRARIES}
																				${AUXIL_LIBS}
																				${GFORTRAN_LIBRARY})

	target_link_directories( BCMshared PUBLIC ${SPRAL_LIBDIR} ${HSL_LIBDIR} ${LEMON_LIBDIR})

	target_compile_definitions( BCMshared PUBLIC ${HSL_FLAGS}
																				 ${SPRAL_FLAGS}
																				 ${SUPERLU_FLAGS}
																				 ${AUXIL_FLAGS})

	set_target_properties(BCMshared PROPERTIES OUTPUT_NAME BCM)

	install(TARGETS BCMshared
	        LIBRARY
	          DESTINATION lib
	          COMPONENT Libraries
	          NAMELINK_COMPONENT Development
	       )

endif(BUILD_SHARED_LIBS)
