# Make an automatic library - will be static or dynamic based on user setting
add_library(amgmatch
    matchTwoThirdeps.cpp
    matchOpt.cpp
    amgmatchUtil.cpp
    )

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(amgmatch PUBLIC OpenMP::OpenMP_CXX)
endif()

#link light_graph
target_link_libraries(amgmatch PRIVATE matchbox ${LEMON_LIBS})
target_link_directories(amgmatch PRIVATE ${LEMON_LIBDIR})

target_include_directories(amgmatch PUBLIC
    "${LEMON_INCDIR}"
    "${PROJECT_SOURCE_DIR}/include"
    "${PROJECT_SOURCE_DIR}/include/Matchbox"
    )

# All users of this library will need at least C++11
target_compile_options(amgmatch PUBLIC -DSTATS -DMKL_ILP64 -m64 -O3 -ffast-math)

# Install AMGMatch library
install(TARGETS amgmatch
        LIBRARY
          DESTINATION lib
          COMPONENT Libraries
          NAMELINK_COMPONENT Development
        PUBLIC_HEADER
          DESTINATION include
          COMPONENT Development
       )

set_property(TARGET matchbox PROPERTY
                PUBLIC_HEADER ${PROJECT_SOURCE_DIR}/AUXIL/include/AMGMatch/amgmatch.h)
