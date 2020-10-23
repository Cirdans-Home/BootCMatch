# MATVEC

set(matvec_sources 	MATVEC/bcm_csr_matop.c
                  	MATVEC/bcm_csr_matrix.c
                    MATVEC/bcm_csr_matvec.c
                    MATVEC/bcm_csr_relax.c
                    MATVEC/bcm_vector.c
										MATVEC/bcm_matvec.h PARENT_SCOPE)

set(matvec_header ${CMAKE_SOURCE_DIR}/SRC/MATVEC/bcm_matvec.h PARENT_SCOPE)
