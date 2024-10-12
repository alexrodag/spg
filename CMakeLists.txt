cmake_minimum_required(VERSION 3.18)

# Set required C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

find_package(OpenMP REQUIRED)

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type specified, defaulting to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

if (MSVC)
    add_definitions(/MP) # Parallel compilation
    # add_definitions(/W4) # Show level 4 warnings
    # add_definitions(/WX) # Treat warnings as errors
endif ()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # add_definitions(-Wall) # Enable all warnings
    # add_definitions(-Werror) # Treat warnings as errors (optional)
endif()

add_subdirectory(spg)

add_subdirectory(apps)