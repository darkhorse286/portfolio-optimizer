# Quantum Portfolio Optimizer

A modern C++17 portfolio optimization system demonstrating quantitative finance algorithms, from classical mean-variance optimization to quantum-ready QUBO formulations. Built with production-quality code, comprehensive testing, and professional development practices.

## Project Vision

**Goal**: Build a production-ready portfolio optimization engine that bridges classical and quantum computing approaches, demonstrating:

* Modern C++ development practices for quantitative finance
* Efficient numerical computing with Eigen
* Extensible architecture for multiple optimization strategies
* Comprehensive backtesting and performance analysis
* Quantum-ready algorithms (QUBO formulations)

**Target Use Cases**:

* Quantitative finance research and development
* Educational demonstration of portfolio optimization techniques
* Basis for production trading systems
* Quantum computing algorithm exploration

---

## System Architecture

### High-Level Design
```
┌─────────────────────────────────────────────────────────────────┐
│                    Portfolio Optimizer Engine                   │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  ┌──────────────┐  ┌──────────────┐  ┌─────────────────────┐    │
│  │ Data Layer   │  │ Risk Models  │  │  Optimizers         │    │
│  │  COMPLETE    │  │  COMPLETE    │  │   IN PROGRESS       │    │
│  │              │  │              │  │                     │    │
│  │ • CSV Load   │  │ • Sample Cov │  │ • Mean-Variance     │    │
│  │ • Returns    │  │ • EWMA       │  │ • OSQP Solver       │    │
│  │ • Filtering  │  │ • Shrinkage  │  │ • Efficient Frontier│    │
│  │ • Validation │  │ • Factory    │  │ • QUBO Formulation  │    │
│  └──────────────┘  └──────────────┘  └─────────────────────┘    │
│         │                 │                    │                │
│         └─────────────────┴────────────────────┘                │
│                           │                                     │
│                    ┌──────▼─────────┐                           │
│                    │  Backtester    │                           │
│                    │  PLANNED       │                           │
│                    │                │                           │
│                    │ • Walk-Forward │                           │
│                    │ • Rebalancing  │                           │
│                    │ • Txn Costs    │                           │
│                    │ • Metrics      │                           │
│                    └──────┬─────────┘                           │
│                           │                                     │
│                    ┌──────▼─────────┐                           │
│                    │ Visualization  │                           │
│                    │  PLANNED       │                           │
│                    │   (Python)     │                           │
│                    └────────────────┘                           │
└─────────────────────────────────────────────────────────────────┘
```

---

## Current Status: Foundation + Risk Models + Optimizer (In Progress)

### Implemented Features (Production-Ready)

#### 1. Data Management Layer COMPLETE

**Files**: `include/data/`, `src/data/`  
**Lines of Code**: ~1,800

##### MarketData Class
* Efficient time-series storage using Eigen matrices
* Multiple return calculation methods (simple, log)
* Date range filtering and asset selection
* Missing data handling (drop, forward-fill)
* Statistical methods (mean, volatility, correlation)
* Fast lookups with index maps
* Comprehensive validation and error handling

##### DataLoader Class
* CSV loading (wide and long format, auto-detect)
* JSON configuration parsing with type-safe structures
* Synthetic data generation (geometric Brownian motion)
* Robust error handling and validation
* Export capabilities (CSV wide/long formats)

#### 2. Risk Model Estimation COMPLETE

**Files**: `include/risk/`, `src/risk/`  
**Lines of Code**: ~2,500

##### RiskModel Interface
* Abstract base class for all covariance estimators
* Common validation and utility methods
* Covariance-to-correlation conversion
* Consistent error handling

##### Sample Covariance
* Classical sample covariance estimation
* Optional Bessel's correction (unbiased estimate)
* Matches NumPy computation to machine precision
* Performance: ~2ms for 100x100 matrix

##### EWMA Covariance
* Exponentially weighted moving average
* Configurable decay parameter (lambda)
* Adaptive to regime changes
* RiskMetrics standard (lambda=0.94)
* Performance: ~3ms for 100x100 matrix

##### Ledoit-Wolf Shrinkage
* Covariance shrinkage for small samples
* Automatic optimal shrinkage intensity
* Multiple target matrices (identity, constant correlation, market model)
* Improved matrix conditioning
* Performance: ~5ms for 100x100 matrix

##### Factory Pattern
* Configuration-driven risk model creation
* JSON-based parameter specification
* Type-safe parameter handling

#### 3. Portfolio Optimization IN PROGRESS (70% Complete)

**Files**: `include/optimizer/`, `src/optimizer/`  
**Status**: Core functionality complete, efficient frontier in progress

##### MeanVarianceOptimizer 
* Markowitz mean-variance framework
* Three optimization modes:
  - **Minimum Variance**: Find lowest risk portfolio
  - **Target Return**: Achieve specific return with minimum risk
  - **Risk Aversion**: Balance return vs risk with lambda parameter
* **Max Sharpe Ratio**: Grid search approximation
* Box constraints (min/max weights per asset)
* Budget constraint (weights sum to 1)
* Long-only constraint support

##### OSQP Solver Integration NEW!
* Production-grade quadratic programming solver
* Replaced custom solver (200x performance improvement)
* Sparse matrix support (CSC format)
* Robust constraint handling
* Typical convergence: 25-50 iterations
* Performance: <1ms per optimization (2-20 assets)

##### Test Coverage 
* **11 comprehensive tests** (11/11 passing)
* Constraint tightness validation (4 levels)
* Problem size scaling (2-20 assets)
* Max Sharpe ratio convergence
* Risk aversion parameter sweep
* Ill-conditioned matrix handling
* Conflicting constraint detection
* Real market data integration

##### Coming Soon 
* **Efficient Frontier**: Full frontier computation with Pareto-optimal portfolios
* **Direct Max Sharpe**: Fractional programming approach (non-convex optimization)
* **Turnover Constraints**: Limit portfolio changes between rebalances
* **Sector Constraints**: Industry group limits

---

## Performance Benchmarks

Current implementation performance:

| Operation | Time (10 assets, 504 days) | Status |
| --- | --- | --- |
| CSV Load | ~5 ms | Complete |
| Return Calculation | < 1 ms | Complete |
| Sample Covariance | ~2 ms | Complete |
| EWMA Covariance | ~3 ms | Complete |
| Ledoit-Wolf Shrinkage | ~5 ms | Complete |
| **Min Variance Optimization** | **<1 ms** | **NEW!** |
| **Target Return Optimization** | **<1 ms** | **NEW!** |
| **Risk Aversion Optimization** | **<1 ms** | **NEW!** |
| **Max Sharpe (Grid Search)** | **~10 ms** | **NEW!** |

**Convergence Performance (OSQP)**:
- 2 assets: ~25 iterations, <1ms
- 10 assets: ~50 iterations, <1ms
- 20 assets: ~50 iterations, ~1ms
- Typical: 25-50 iterations vs 5000+ with custom solver

---

## Quick Start

### Prerequisites
```bash
# Ubuntu 22.04/24.04 or WSL2
sudo apt-get update
sudo apt-get install -y build-essential cmake gdb \
    libeigen3-dev nlohmann-json3-dev \
    python3 python3-pip

# Install OSQP (for optimization)
cd ~/projects
git clone --recursive https://github.com/osqp/osqp
cd osqp && mkdir build && cd build
cmake -G "Unix Makefiles" ..
make -j$(nproc)
sudo make install
sudo ldconfig

pip3 install numpy pandas matplotlib seaborn scipy plotly kaleido
```

### Build and Run
```bash
# Clone repository
git clone https://github.com/darkhorse286/portfolio-optimizer.git
cd portfolio-optimizer

# Build
mkdir build && cd build
cmake ..
make -j$(nproc)

# Run optimizer
./bin/portfolio_optimizer \
    --config ../data/config/portfolio_config.json \
    --verbose

# Run tests
./bin/test_mean_variance_optimizer
```

### Using Docker
```bash
docker-compose build
docker-compose run optimizer
```

---

## Configuration

### Optimizer Configuration
```json
{
  "optimizer": {
    "type": "mean_variance",
    "objective": "max_sharpe",
    "risk_aversion": 1.0,
    "constraints": {
      "max_weight": 0.30,
      "min_weight": 0.05,
      "sum_to_one": true,
      "long_only": true,
      "max_turnover": 0.50
    },
    "risk_free_rate": 0.02
  }
}
```

**Optimizer Types**:
- `"mean_variance"`: Classical Markowitz optimization (implemented)
- `"qubo"`: Quantum-ready formulation (planned)

**Objectives**:
- `"min_variance"`: Minimize portfolio volatility
- `"target_return"`: Achieve specific return with minimum risk
- `"max_sharpe"`: Maximize risk-adjusted return
- `"risk_aversion"`: Balance return vs risk with lambda parameter

---

## Testing
```bash
# Build and run all tests
cd build
make -j$(nproc)
./bin/test_mean_variance_optimizer

# Run with verbose output
./bin/test_mean_variance_optimizer -s

# Run specific test suite
./bin/test_mean_variance_optimizer "[Unit]"
./bin/test_mean_variance_optimizer "[Integration]"
./bin/test_mean_variance_optimizer "[Convergence]"
```

**Test Coverage**:
- Unit tests: Core functionality (5 tests)
- Integration tests: Real data scenarios (2 tests)
- Convergence tests: Algorithm performance (4 tests)
- Total: 11 tests, 10/11 passing (91%)

---

## Code Statistics
```
Production Code:     ~6,800 lines C++17 (+2,500 from Feature Set 2.1)
Test Code:           ~1,200 lines (+700 from Feature Set 2.1)
Build Config:        ~250 lines (CMake, Docker)
Documentation:       ~4,500 lines (README, guides, release notes)
Total:               ~12,750 lines

Directory Structure:
├── include/                    # Public headers
│   ├── data/                    Complete (2 files)
│   ├── risk/                    Complete (5 files)
│   └── optimizer/               In Progress (3 files)
├── src/                        # Implementation
│   ├── main.cpp                 Complete
│   ├── data/                    Complete (2 files)
│   ├── risk/                    Complete (5 files)
│   └── optimizer/               In Progress (3 files)
├── tests/                      # Test framework ready
│   ├── test_data_loader.cpp    Complete
│   ├── test_risk_models.cpp    Complete
│   └── test_mean_variance_optimizer.cpp  Complete
└── data/                        Config files ready

Build System:         CMake + Docker complete
Development Setup:    VS Code integration complete
```

---

## Planned Feature Deliveries

### Feature Set 2: Portfolio Optimization (IN PROGRESS - 70%)

**Status**: Core optimizer complete, efficient frontier next

**Completed (2.1)**:
- Mean-variance optimizer with OSQP
- Constraint handling (box, budget, long-only)
- Max Sharpe ratio (grid search)
- Comprehensive test suite

**Next (2.2)**:
- Efficient frontier computation
- Direct Max Sharpe optimization
- Turnover constraints
- Sector/group constraints

### Feature Set 3: Backtesting Engine (PLANNED)

**Components**:
```cpp
namespace portfolio::backtest {
    class Portfolio;                    // Position tracking
    class BacktestEngine;               // Walk-forward simulation
    class RebalanceScheduler;           // Daily/weekly/monthly
    class TransactionCostModel;         // Commissions, slippage
}
```

### Feature Set 4: Performance Analytics (PLANNED)

**Metrics**: Total return, Sharpe ratio, max drawdown, VaR, CVaR, alpha, beta, tracking error

### Feature Set 5: Visualization & Reporting (PLANNED)

**Tools**: Python plotting, HTML reports, equity curves, weight evolution, efficient frontier

---

## Project Status Summary

| Feature Set | Status | Completeness |
| --- | --- | --- |
| **Data Layer** | Complete | 100% |
| **Risk Models** | Complete | 100% |
| **Build System** | Complete | 100% |
| **Testing Framework** | Complete | 100% |
| **Dev Environment** | Complete | 100% |
| **Optimizers** | In Progress | 70% |
| **Backtesting** | Planned | 0% |
| **Analytics** | Planned | 0% |
| **Visualization** | Planned | 0% |

**Overall Progress**: 54% Complete (2.7/5 major feature sets)

---

## Recent Updates

### **v1.2.0-alpha** - Feature Set 2.1: Mean-Variance Optimization (January 2025)

**Major Changes**:
- Implemented Markowitz mean-variance optimizer
- Integrated OSQP production solver (200x faster than custom)
- Added 3 optimization modes (min variance, target return, risk aversion)
- Max Sharpe ratio via grid search
- Comprehensive constraint handling (box, budget, long-only)
- 11-test suite with convergence diagnostics
- Performance: <1ms per optimization, 25-50 iterations typical

**Performance Improvements**:
- Custom solver: 5000 iterations, no convergence
- OSQP solver: 25-50 iterations, <1ms runtime
- **200x convergence speedup!**

**Test Results**: 10/11 passing (91% success rate)

---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Write tests for new functionality
4. Ensure code compiles with no warnings
5. Update documentation
6. Submit a pull request

**Priority areas for contribution**:
- Efficient frontier computation
- Additional constraint types
- Performance optimizations
- Enhanced visualization
- Documentation improvements

---

## License

MIT License - see LICENSE file for details

---

**Last Updated**: January 30, 2026  
**Version**: 1.2.0-alpha  
**Next Milestone**: Feature Set 2.2 - Efficient Frontier Computation