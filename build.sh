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