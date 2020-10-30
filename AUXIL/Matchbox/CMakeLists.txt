# Make an automatic library - will be static or dynamic based on user setting
add_library(matchbox STATIC
    BptGraph.cc
    BptMatchingEngine.cc
    AdjGraphEngine.cc
    BptCanonicalEngine.cc
    BptEdgCovEngine.cc
    BptIndSetEngine.cc
    BptVtxCovEngine.cc
    ConvertEngine.cc
    Graph.cc
    MatchingEngine.cc
    Stats.cc
    Timer.cc
    )

# We need this directory, and users of our library will need it too
target_include_directories(matchbox
        PUBLIC "${CMAKE_SOURCE_DIR}/AUXIL/include"
               "${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox" )

# All users of this library will need at least C++11
target_compile_features(matchbox PUBLIC cxx_std_11 )
target_compile_options(matchbox PUBLIC -DSTATS)
target_compile_options(matchbox PRIVATE -ansi)

# Set the PUBLIC_HEADER to install
set_property(TARGET matchbox PROPERTY
            PUBLIC_HEADER ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/AdjGraphEngine.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/BptEdgCovEngine.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/BptVtxCovEngine.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/Graph.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/LstItmAndValIdxdQue.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/Timer.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/VecItmQue.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/BhpItmAndValIdxdMaxPriQue.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/BptGraph.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/ConvertEngine.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/gtest.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/LstItmIdxdQue.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/MatchingEngine.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/Types.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/VecItmStk.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/BhpItmAndValIdxdMinPriQue.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/BptIndSetEngine.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/EdgeSet.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/LstItmAndValIdxdMaxPriQue.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/LstItmQue.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/MatchingEngine.h.gch
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/Utility.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/VertexSet.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/BptCanonicalEngine.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/BptMatchingEngine.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/Functional.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/LstItmAndValIdxdMinPriQue.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/LstItmStk.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/Stats.h
                    ${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox/Utility.h.gch)

install(TARGETS matchbox
        LIBRARY
          DESTINATION lib
          COMPONENT Libraries
          NAMELINK_COMPONENT Development
        PUBLIC_HEADER
          DESTINATION include
          COMPONENT Development
       )

if(BUILD_SHARED_LIBS)
    add_library(matchboxshared SHARED
        BptGraph.cc
        BptMatchingEngine.cc
        AdjGraphEngine.cc
        BptCanonicalEngine.cc
        BptEdgCovEngine.cc
        BptIndSetEngine.cc
        BptVtxCovEngine.cc
        ConvertEngine.cc
        Graph.cc
        MatchingEngine.cc
        Stats.cc
        Timer.cc
    )

    # We need this directory, and users of our library will need it too
    target_include_directories(matchboxshared
            PUBLIC "${CMAKE_SOURCE_DIR}/AUXIL/include"
                   "${CMAKE_SOURCE_DIR}/AUXIL/include/Matchbox" )

    # All users of this library will need at least C++11
    target_compile_features(matchboxshared PUBLIC cxx_std_11 )
    target_compile_options(matchboxshared PUBLIC -fPIC -DSTATS)
    target_compile_options(matchboxshared PRIVATE -ansi)

    set_target_properties(matchboxshared PROPERTIES OUTPUT_NAME matchbox)

    install(TARGETS matchboxshared
            LIBRARY
              DESTINATION lib
              COMPONENT Libraries
              NAMELINK_COMPONENT Development
           )
endif(BUILD_SHARED_LIBS)
