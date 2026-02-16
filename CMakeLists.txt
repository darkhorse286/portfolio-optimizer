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

find_library(OSQP_LIBRARY 
    NAMES osqp 
    PATHS /usr/local/lib /usr/lib
    REQUIRED
)
find_path(OSQP_INCLUDE_DIR 
    NAMES osqp/osqp.h 
    PATHS /usr/local/include /usr/include
    REQUIRED
)

message(STATUS "Found OSQP library: ${OSQP_LIBRARY}")
message(STATUS "Found OSQP include: ${OSQP_INCLUDE_DIR}")

# nlohmann/json - header-only library
find_package(nlohmann_json 3.2.0 QUIET)
if(NOT nlohmann_json_FOUND)
    message(STATUS "nlohmann/json not found, using system headers")
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

# =======================
# Optimizer Library
# =======================
file(GLOB_RECURSE OPTIMIZER_SOURCES "src/optimizer/*.cpp")

add_library(optimizer ${OPTIMIZER_SOURCES})
target_include_directories(optimizer PUBLIC 
    ${CMAKE_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
    ${OSQP_INCLUDE_DIR}  
)
target_link_libraries(optimizer PUBLIC 
    Eigen3::Eigen
    risk
    ${OSQP_LIBRARY}  
)

if(nlohmann_json_FOUND)
    target_link_libraries(optimizer PUBLIC nlohmann_json::nlohmann_json)
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
    optimizer
    Eigen3::Eigen
)

if(nlohmann_json_FOUND)
    target_link_libraries(portfolio_optimizer PRIVATE nlohmann_json::nlohmann_json)
endif()

# =======================
# Backtest Library
# =======================
file(GLOB_RECURSE BACKTEST_SOURCES "src/backtest/*.cpp")
if(BACKTEST_SOURCES)
    add_library(backtest ${BACKTEST_SOURCES})
    target_include_directories(backtest PUBLIC
        ${CMAKE_SOURCE_DIR}/include
        ${EIGEN3_INCLUDE_DIR}
    )
    target_link_libraries(backtest PUBLIC
        Eigen3::Eigen
        data
        risk
        optimizer
    )
    if(nlohmann_json_FOUND)
        target_link_libraries(backtest PUBLIC nlohmann_json::nlohmann_json)
    endif()
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

# Discover and register all tests in tests/ directory
file(GLOB TEST_SOURCES "${CMAKE_SOURCE_DIR}/tests/*.cpp")
set(ALL_TEST_TARGETS "")

# Common test libraries (conditionally append existing targets)
set(TEST_LIBS Catch2::Catch2WithMain Eigen3::Eigen)
if(TARGET data)
    list(APPEND TEST_LIBS data)
endif()
if(TARGET risk)
    list(APPEND TEST_LIBS risk)
endif()
if(TARGET optimizer)
    list(APPEND TEST_LIBS optimizer)
endif()
if(TARGET backtest)
    list(APPEND TEST_LIBS backtest)
endif()

foreach(test_src ${TEST_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)
    set(target_name test_${test_name})
    add_executable(${target_name} ${test_src})
    target_include_directories(${target_name} PRIVATE ${CMAKE_SOURCE_DIR}/include ${EIGEN3_INCLUDE_DIR})
    target_link_libraries(${target_name} PRIVATE ${TEST_LIBS})
    add_test(NAME ${test_name} COMMAND ${target_name})
    list(APPEND ALL_TEST_TARGETS ${target_name})
endforeach()

# Aggregate run_tests target
add_custom_target(run_tests
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    DEPENDS ${ALL_TEST_TARGETS}
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
message(STATUS "Libraries:")
message(STATUS "  - Data layer: ENABLED")
message(STATUS "  - Risk models: ENABLED")
message(STATUS "  - Optimizer: ENABLED")
message(STATUS "=============================")
message(STATUS "")