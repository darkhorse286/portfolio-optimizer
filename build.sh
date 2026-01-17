#!/bin/bash

# Portfolio Optimizer Build Script
# Usage: ./build.sh [clean|test|run|help]

set -e # Exit on error

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$PROJECT_DIR/build"
BUILD_TYPE="${BUILD_TYPE:-Release}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_header(){
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE} Portfolio Optimizer Build Script ${NC} "
    echo -e "${BLUE}========================================${NC}"
    echo ""
}

print_success(){
    echo -e "${GREEN}[SUCCESS]${NC} $1${NC}"
}

print_error(){
    echo -e "${RED}[ERROR]${NC} $1${NC}"
}

print_info(){
    echo -e "${YELLOW}[INFO]${NC} $1${NC}"
}

clean_build(){
    print_info "Cleaning build directory..."
    if [ -d "$BUILD_DIR" ]; then
        rm -rf "$BUILD_DIR"
        print_success "Build directory cleaned."
    else
        print_info "Build directory does not exist. Nothing to clean."
    fi
    
}

configure_cmake(){
    print_info "Configuring CMake..."
    mkdir -p "${BUILD_DIR}"
    cd "${BUILD_DIR}"
    cmake -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
            -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
            "${PROJECT_DIR}"

    print_success "CMake configuration completed."
}

build_project(){
    print_info "Building project..."
    cd "${BUILD_DIR}"
    
    # Get number of CPU cores
    if ["$(uname)" == "Darwin"]; then
        NPROC=$(sysctl -n hw.ncpu)
    else
        NPROC=$(nproc)
    fi
    
    make -j${NPROC}

    print_success "Build complete."
}

run_tests(){
    print_info "Running tests..."
    cd "${BUILD_DIR}"
    
    if [ ! -f "bin/run_tests"]; then
        print_error "Test executable not found. Please build the project first."
        exit 1
    fi

    .bin/run_tests
    print_success "All tests passed."
}

run_optimizer(){
    print_info "Running Portfolio Optimizer..."
    cd "${BUILD_DIR}"
    
    if [ ! -f "bin/portfolio_optimizer"]; then
        print_error "Optimizer executable not found. Please build the project first."
        exit 1
    fi

    "${BUILD_DIR}/bin/portfolio_optimizer" \
    --config data/config/portfolio_config.json \
    --verbose
}

show_usage(){
    cat << EOF
Usage: ./build.sh [COMMAND]

Commands:
    clean       Remove build directory
    config      Configure CMake
    build       Build the project (default)
    test        Build and run tests
    run         Build and run optimizer
    all         Clean, configure, build, test, and run
    help        Show this help message

Environment Variables:
    BUILD_TYPE  Build type: Release (default), Debug, RelWithDebInfo

Examples:
    ./build.sh              # Just build
    ./build.sh clean build  # Clean and build
    ./build.sh test         # Build and run tests
    ./build.sh run          # Build and run optimizer
    BUILD_TYPE=Debug ./build.sh build  # Debug build

EOF
}

# Main script logic
print_header
if [ $# -eq 0 ]; then
    # Default action is to build
    if [ ! -d "${BUILD_DIR}" ]; then
        configure_cmake
    fi
    build_project
    exit 0
fi

for cmd in "$@"; do
    case "$cmd" in
        clean)
            clean_build
            ;;
        config)
            configure_cmake
            ;;
        build)
            if [ ! -d "${BUILD_DIR}" ]; then
                configure_cmake
            fi
            build_project
            ;;
        test)
            if [ ! -d "${BUILD_DIR}" ]; then
                configure_cmake
            fi
            build_project
            run_tests
            ;;
        run)
            if [ ! -d "${BUILD_DIR}" ]; then
                configure_cmake
            fi
            build_project
            run_optimizer
            ;;
        all)
            clean_build
            configure_cmake
            build_project
            run_tests
            run_optimizer
            ;;
        help|--help|-h)
            show_usage
            exit 0
            ;;
        *)
            print_error "Unknown command: $cmd"
            show_usage
            exit 1
            ;;
    esac
done

print_success "All requested operations completed."