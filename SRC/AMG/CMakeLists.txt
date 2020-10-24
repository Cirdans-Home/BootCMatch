# AMG

set(amg_sources 	AMG/bcm_amg.c
								 	AMG/bcm_linmatch.c
									AMG/bcm_weightmat.c
									AMG/bcm_adaptivecoarsening.c
									AMG/bcm_callcxxfunction.cpp
									AMG/bcm_amg.h PARENT_SCOPE)

set(amg_header ${CMAKE_SOURCE_DIR}/SRC/AMG/bcm_amg.h PARENT_SCOPE)
