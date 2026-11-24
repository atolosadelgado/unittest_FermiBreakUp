add_library(mydebuglibs src/G4DebugGeneratorPrecompoundInterface.cc src/G4DebugPreCompoundModel.cc)
target_link_libraries(mydebuglibs
    ${Geant4_LIBRARIES}
)
target_include_directories(mydebuglibs
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
    mydebuglibs
)
