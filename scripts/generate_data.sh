#!/bin/bash
# Generate synthetic market data for portfolio optimizer

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

echo "=== Synthetic Data Generator ==="
echo ""

# Compile the data generator
echo "Compiling data generator..."
g++ -std=c++17 \
    -I"${PROJECT_DIR}/include" \
    -I/usr/include/eigen3 \
    "${SCRIPT_DIR}/generate_synthetic_data.cpp" \
    "${PROJECT_DIR}/src/data/market_data.cpp" \
    "${PROJECT_DIR}/src/data/data_loader.cpp" \
    -o "${PROJECT_DIR}/build/generate_data" \
    -O2

if [ $? -eq 0 ]; then
    echo "✓ Compilation successful"
else
    echo "✗ Compilation failed"
    exit 1
fi

echo ""

# Run the generator
echo "Running data generator..."
"${PROJECT_DIR}/build/generate_data" "$@"

echo ""
echo "Done! Data saved to data/market/historical_prices.csv"