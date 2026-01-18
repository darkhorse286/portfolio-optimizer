cmake_minimum_required(VERSION 3.15)
project(PortfolioOptimizer VERSION 1.0.0 LANGUAGES CXX)

# C++17 is required
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Intellisense fix for VS Code
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Find required packages
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(nlohmann_json 3.11.3 REQUIRED)

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
    ${nlohmann_json_INCLUDE_DIRS}
)

# Collect source files
file(GLOB_RECURSE SOURCES 
    "src/data/*.cpp"
    "src/optimizer/*.cpp"
    "src/risk/*.cpp"
    "src/backtest/*.cpp"
)

# Main executable
add_executable(portfolio_optimizer
    src/main.cpp
    ${SOURCES}
)

target_link_libraries(portfolio_optimizer
    Eigen3::Eigen
    nlohmann_json::nlohmann_json
)

# Enable testing
enable_testing()

# Download Catch2 for unit tests
include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.0
)
FetchContent_MakeAvailable(Catch2)

# Test executable
file(GLOB_RECURSE TEST_SOURCES
    "tests/*.cpp"
)

if(TEST_SOURCES)
    add_executable(run_tests ${TEST_SOURCES} ${SOURCES})
    target_link_libraries(run_tests
        Eigen3::Eigen
        nlohmann_json::nlohmann_json
        Catch2::Catch2WithMain
    )
    # Register tests with CTest
    include(CTest)
    include(Catch)
    catch_discover_tests(run_tests)
endif()

# Installation rules
install(TARGETS portfolio_optimizer RUNTIME DESTINATION bin)

# Print configuration summary
message(STATUS "")
message(STATUS "=== Portfolio Optimizer Configuration ===")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Eigen3 version: ${EIGEN3_VERSION}")
message(STATUS "nlohmann_json version: ${nlohmann_json_VERSION}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX} ")
message(STATUS "=========================================")
message(STATUS "")