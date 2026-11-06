cmake_minimum_required(VERSION 3.16)

project(testFermiBreakUp)

find_package(ROOT REQUIRED COMPONENTS Core Hist RIO)

find_package(Geant4 REQUIRED)

include(${Geant4_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR})

add_executable(testFermiBreakUp testFermiBreakUp.cc)

target_link_libraries(testFermiBreakUp
    ${Geant4_LIBRARIES}
)

add_executable(testFTFP testFTFP.cc)

target_link_libraries(testFTFP
    ${Geant4_LIBRARIES}
)
target_link_libraries(testFTFP
    ${Geant4_LIBRARIES}
    ROOT::Core
    ROOT::Hist
    ROOT::RIO
)
