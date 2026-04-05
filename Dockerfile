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
    libblas-dev \
    liblapack-dev \
    python3 \
    python3-pip \
    python3-dev \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages for visualization
# Create an isolated virtual environment to install Python packages
RUN python3 -m venv /opt/venv && \
    /opt/venv/bin/python -m pip install --upgrade pip setuptools wheel && \
    /opt/venv/bin/pip install --no-cache-dir \
        numpy \
        pandas \
        matplotlib \
        seaborn \
        plotly \
        kaleido \
        tabulate \
        jinja2

# Ensure the virtualenv is used by default in the container
ENV PATH="/opt/venv/bin:$PATH"

# Set working directory
WORKDIR /app

# Copy project files
COPY . /app

# Create build directory and compile the project
# Remove any existing build directory (copied from host) to avoid stale CMake cache
# Build and install OSQP from source (libosqp-dev may not be available on all base images)
RUN git clone --depth 1 https://github.com/osqp/osqp.git /tmp/osqp && \
    cd /tmp/osqp && mkdir -p build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    cmake --build . -j$(nproc) && \
    cmake --install . && \
    rm -rf /tmp/osqp

# Build the project
RUN rm -rf build && mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    make -j$(nproc)

# Create results directory
RUN mkdir -p results

# Set environment variables
ENV PYTHONPATH=/app/scripts:$PYTHONPATH
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

# Default command to run when container starts
CMD ["./build/bin/portfolio_optimizer", "--config", "data/config/portfolio_config.json"]
