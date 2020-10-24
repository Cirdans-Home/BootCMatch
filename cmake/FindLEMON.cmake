# Looking for LEMON Library
set(LEMON_DIR /opt/hsl
						CACHE PATH "/opt/hsl")

find_library(LEMON_LIBRARIES libemon.a
						 PATHS ${LEMON_DIR}/lib REQUIRED)

set(LEMON_INCDIR ${LEMON_DIR}/include)
set(LEMON_LIBDIR ${LEMON_DIR}/lib)
set(LEMON_LIBS -lemon)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
	"LEMON"
	FOUND_VAR LEMON_FOUND
	REQUIRED_VARS LEMON_INCDIR LEMON_LIBDIR LEMON_LIBS
	FAIL_MESSAGE "Couldn't find HSL"
	)
