cmake_minimum_required(VERSION 3.16)

project(testFermiBreakUp)

find_package(Geant4 REQUIRED)

include(${Geant4_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR})

# Enable UndefinedBehaviorSanitizer
option(ENABLE_UBSAN "Enable Undefined Behavior Sanitizer" OFF)

if(ENABLE_UBSAN)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(STATUS "Undefined Behavior Sanitizer enabled")

        add_compile_options(
            -fsanitize=undefined
            -fno-sanitize-recover=undefined
            -fno-omit-frame-pointer
            -g
        )

        add_link_options(
            -fsanitize=undefined
        )
    else()
        message(WARNING "UBSan is not supported with this compiler")
    endif()
endif()



add_subdirectory(testFBU)

add_subdirectory(testFTF)

add_subdirectory(testAnalysis)


