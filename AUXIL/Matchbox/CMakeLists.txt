# Make an automatic library - will be static or dynamic based on user setting
add_library(matchbox
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
        PUBLIC "${PROJECT_SOURCE_DIR}/AUXIL/include"
               "${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox" )

# All users of this library will need at least C++11
target_compile_features(matchbox PUBLIC cxx_std_11 )
target_compile_options(matchbox PUBLIC -ansi -DSTATS)

# Set the PUBLIC_HEADER to install
set_property(TARGET matchbox PROPERTY
            PUBLIC_HEADER ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/AdjGraphEngine.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/BptEdgCovEngine.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/BptVtxCovEngine.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/Graph.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/LstItmAndValIdxdQue.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/Timer.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/VecItmQue.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/BhpItmAndValIdxdMaxPriQue.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/BptGraph.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/ConvertEngine.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/gtest.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/LstItmIdxdQue.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/MatchingEngine.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/Types.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/VecItmStk.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/BhpItmAndValIdxdMinPriQue.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/BptIndSetEngine.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/EdgeSet.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/LstItmAndValIdxdMaxPriQue.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/LstItmQue.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/MatchingEngine.h.gch
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/Utility.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/VertexSet.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/BptCanonicalEngine.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/BptMatchingEngine.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/Functional.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/LstItmAndValIdxdMinPriQue.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/LstItmStk.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/Stats.h
                    ${PROJECT_SOURCE_DIR}/AUXIL/include/Matchbox/Utility.h.gch)

install(TARGETS matchbox
        LIBRARY
          DESTINATION lib
          COMPONENT Libraries
          NAMELINK_COMPONENT Development
        PUBLIC_HEADER
          DESTINATION include
          COMPONENT Development
       )
