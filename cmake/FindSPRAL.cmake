# Looking for SPRAL Library
set(SPRAL_DIR /opt/hsl
						CACHE PATH "/opt/hsl")

find_library(SPRAL_LIBRARIES libspral.a
						 PATHS ${SPRAL_DIR}/lib REQUIRED)

set(SPRAL_INCDIR ${SPRAL_DIR}/include)
set(SPRAL_LIBDIR ${SPRAL_DIR}/lib)
set(SPRAL_LIBS -lspral)
set(SPRAL_FLAGS HAVE_SPRAL)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
	"SPRAL"
	FOUND_VAR SPRAL_FOUND
	REQUIRED_VARS SPRAL_INCDIR SPRAL_LIBDIR SPRAL_LIBS SPRAL_FLAGS
	FAIL_MESSAGE "Couldn't find SPRAL"
	)
