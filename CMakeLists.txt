# Portfolio Optimizer - CMakeLists.txt
cmake_minimum_required(VERSION 3.15)
project(PortfolioOptimizer VERSION 1.0.0 LANGUAGES CXX)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# Export compile commands for IntelliSense
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Find required packages
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# nlohmann/json - header-only library
find_package(nlohmann_json 3.2.0 QUIET)
if(NOT nlohmann_json_FOUND)
    message(STATUS "nlohmann/json not found, using system headers")
    # Assume it's installed in system include path
endif()

# =====================
# Data Layer Library
# =====================
file(GLOB_RECURSE DATA_SOURCES "src/data/*.cpp")

add_library(data ${DATA_SOURCES})
target_include_directories(data PUBLIC 
    ${CMAKE_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
)
target_link_libraries(data PUBLIC 
    Eigen3::Eigen
    risk
)

if(nlohmann_json_FOUND)
    target_link_libraries(data PUBLIC nlohmann_json::nlohmann_json)
endif()

# ====================
# Risk Model Library
# ====================
file(GLOB_RECURSE RISK_SOURCES "src/risk/*.cpp")

add_library(risk ${RISK_SOURCES})
target_include_directories(risk PUBLIC 
    ${CMAKE_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
)
target_link_libraries(risk PUBLIC 
    Eigen3::Eigen
)

if(nlohmann_json_FOUND)
    target_link_libraries(risk PUBLIC nlohmann_json::nlohmann_json)
endif()

# ===================
# Main Executable
# ===================
add_executable(portfolio_optimizer src/main.cpp)
target_include_directories(portfolio_optimizer PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
)
target_link_libraries(portfolio_optimizer PRIVATE 
    data
    risk
    Eigen3::Eigen
)

if(nlohmann_json_FOUND)
    target_link_libraries(portfolio_optimizer PRIVATE nlohmann_json::nlohmann_json)
endif()

# =======================
# Testing with Catch2
# =======================
include(FetchContent)

# Download Catch2
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.0
)
FetchContent_MakeAvailable(Catch2)

# Enable testing
enable_testing()

# Data layer tests
add_executable(test_data_loader tests/test_data_loader.cpp)
target_link_libraries(test_data_loader PRIVATE 
    data
    risk
    Catch2::Catch2WithMain
)
add_test(NAME DataLoader COMMAND test_data_loader)

# Risk model tests
add_executable(test_risk_models tests/test_risk_models.cpp)
target_link_libraries(test_risk_models PRIVATE 
    risk
    Catch2::Catch2WithMain
)
add_test(NAME RiskModels COMMAND test_risk_models)

# Convenience target to run all tests
add_custom_target(run_tests
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    DEPENDS test_data_loader test_risk_models
    COMMENT "Running all tests..."
)

# ========================
# Installation (optional)
# ========================
install(TARGETS portfolio_optimizer DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)

# ============================
# Print configuration summary
# ============================
message(STATUS "")
message(STATUS "=== Configuration Summary ===")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Eigen3 version: ${Eigen3_VERSION}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "=============================")
message(STATUS "")