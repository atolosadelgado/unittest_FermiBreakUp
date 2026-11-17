add_executable(testFTF src/testFTF.cc)

target_include_directories(testFTF
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(testFTF
    ${Geant4_LIBRARIES}
)
