add_library(debugGenPrecoInterface src/G4DebugGeneratorPrecompoundInterface.cc)
target_link_libraries(debugGenPrecoInterface
    ${Geant4_LIBRARIES}
)
target_include_directories(debugGenPrecoInterface
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_executable(testFTF src/testFTF.cc)

target_include_directories(testFTF
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(testFTF
    ${Geant4_LIBRARIES}
    debugGenPrecoInterface
)
