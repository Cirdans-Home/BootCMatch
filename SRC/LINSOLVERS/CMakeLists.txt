# LINSOLVERS

set(linsolvers_sources 	LINSOLVERS/bcm_boot_solver.c
 											 	LINSOLVERS/bcm_bootpcg.c
 												LINSOLVERS/bcm_linsolvers.h
												LINSOLVERS/bcm.h PARENT_SCOPE)

set(linsolvers_header ${CMAKE_SOURCE_DIR}/SRC/LINSOLVERS/bcm.h
											${CMAKE_SOURCE_DIR}/SRC/LINSOLVERS/bcm_linsolvers.h
											PARENT_SCOPE)
