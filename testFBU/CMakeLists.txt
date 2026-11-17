add_executable(testFermiBreakUp src/testFermiBreakUp.cc)

target_include_directories(testFermiBreakUp
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(testFermiBreakUp
    ${Geant4_LIBRARIES}
)
