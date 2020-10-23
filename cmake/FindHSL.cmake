# Looking for HSL MC64 Library
set(HSL_DIR /opt/hsl
						CACHE PATH "/opt/hsl")

find_library(HSL_LIBRARIES libhsl_mc64.a
						 PATHS ${HSL_DIR}/lib REQUIRED)

set(HSL_INCDIR ${HSL_DIR}/include)
set(HSL_LIBDIR ${HSL_DIR}/lib)
set(HSL_LIBS -lhsl_mc64)
set(HSL_FLAGS HAVE_HSL)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
	"HSL"
	FOUND_VAR HSL_FOUND
	REQUIRED_VARS HSL_INCDIR HSL_LIBDIR HSL_LIBS HSL_FLAGS
	FAIL_MESSAGE "Couldn't find HSL"
	)
