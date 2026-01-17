# Portfolio Optimizer - Docker Container
FROM ubuntu:24.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    gdb \
    git \
    wget \
    libeigen3-dev \
    nlohmann-json3-dev \
    python3 \
    python3-pip \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages for visualization
RUN pip3 install --no-cache-dir \
    numpy==1.24.3 \
    pandas==2.0.3 \
    matplotlib==3.7.2 \
    seaborn==0.12.2 \
    plotly==5.15.0 \
    kaleido==0.2.1 \
    tabulate==0.9.0 \
    jinja2==3.1.2

# Set working directory
WORKDIR /app

# Copy project files
COPY . /app

# Create build directory and compile the project
RUN mkdir -p build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    make -j$(nproc)

# Create results directory
RUN mkdir -p results

# Set environment variables
ENV PYTHONPATH=/app/scripts:$PYTHONPATH
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

# Default command to run when container starts
CMD ["./build/bin/portfolio_optimizer", "--config", "data/config/portfolio_config.json"]
