add_executable(testG4analysis testG4analysis.cc)

target_link_libraries(testG4analysis
    ${Geant4_LIBRARIES}
)
