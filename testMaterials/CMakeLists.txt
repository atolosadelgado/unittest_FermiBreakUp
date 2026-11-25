add_executable(NistMaterialDump NistMaterialDump.cc)

target_link_libraries(NistMaterialDump
    ${Geant4_LIBRARIES}
)
