cmake_minimum_required(VERSION 3.10)
project(RgAnalysis VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Include directories
include_directories(${PROJECT_SOURCE_DIR})

# Source files
set(RGANALYSIS_SOURCES
    RgAnalysis.cpp
    StaticStep.cpp
    DynamicImplicitStep.cpp
    ModalStep.cpp
)

# Header files
set(RGANALYSIS_HEADERS
    RgAnalysis.h
    AnalysisStep.h
    StaticStep.h
    DynamicImplicitStep.h
    ModalStep.h
    StepControl.h
)

# Create library
add_library(RgAnalysis STATIC ${RGANALYSIS_SOURCES})

# Optional: Build examples
option(BUILD_EXAMPLES "Build example programs" ON)

if(BUILD_EXAMPLES)
    # Example executable
    add_executable(ExampleUsage ExampleUsage.cpp)
    target_link_libraries(ExampleUsage RgAnalysis)
endif()

# Installation
install(TARGETS RgAnalysis
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

install(FILES ${RGANALYSIS_HEADERS}
    DESTINATION include/RgAnalysis
)

# Print configuration
message(STATUS "")
message(STATUS "RgAnalysis Configuration:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build examples: ${BUILD_EXAMPLES}")
message(STATUS "")
