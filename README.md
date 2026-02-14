# Portfolio Optimizer

A high-performance C++ portfolio optimization library implementing modern portfolio theory with production-grade solvers.

## Features

### Completed

- **Data Layer**: CSV/JSON loading, return calculation, statistical analysis
- **Risk Models**: Sample covariance, EWMA, Ledoit-Wolf shrinkage
- **Mean-Variance Optimization**: Markowitz portfolio optimization with OSQP solver
- **Efficient Frontier**: Complete Pareto-optimal frontier computation
- **Turnover Constraints**: Per-asset position change limits
 - **Tracking Error (soft penalty)**: Benchmark-relative penalty-based constraint
 - **Tracking Error (soft penalty)**: Benchmark-relative penalty-based constraint
 - **Direct Max Sharpe (Schaible's refinement)**: Fractional programming refinement for direct Max-Sharpe
 - **Sector / Group Constraints**: `SectorMapping` helper and `GroupConstraint` support
 - **Inequality Constraints**: Partial support for inequality constraints via `A_ineq` in QP
- **Build System**: CMake with Docker support
- **Testing**: Comprehensive test suite with Catch2 (65 tests, 100% passing)
 - **Testing**: Comprehensive test suite with Catch2 (coverage expanded for FS 2.3)

### Feature Set 2.3 (Completed)

- Tracking Error soft-constraint (penalty method) for benchmark-relative portfolios
- Integration of tracking-error penalty into all mean-variance objectives
- Auto-tuning heuristic for penalty coefficient
- Unit tests covering tracking-error behavior and validation

### Planned

- Feature Set 3: Backtesting engine with transaction costs
- Feature Set 4: Performance analytics and attribution
- Feature Set 5: Visualization and reporting

## Quick Start

### Prerequisites

- C++17 compiler (GCC 13+ or Clang 14+)
- CMake 3.15+
- Eigen 3.3+
- OSQP 0.6.3+
- nlohmann/json 3.2.0+

### Build

```bash
# Clone repository
git clone https://github.com/darkhorse286/portfolio-optimizer.git
cd portfolio-optimizer

# Build
mkdir build && cd build
cmake ..
make -j$(nproc)

# Run tests
make test

# Run optimizer
./bin/portfolio_optimizer --config ../data/config/portfolio_config.json
```

### Docker Build

```bash
docker build -t portfolio-optimizer .
docker run -it portfolio-optimizer
```

## Usage Examples

### Basic Mean-Variance Optimization

```cpp
#include "optimizer/mean_variance_optimizer.hpp"
#include "data/data_loader.hpp"
#include "risk/sample_covariance.hpp"

// Load market data
auto data = DataLoader::load_csv("prices.csv");

// Calculate returns
auto returns = data.calculate_returns(ReturnType::LOG);

// Estimate covariance
SampleCovariance cov_estimator;
auto covariance = cov_estimator.estimate_covariance(returns);

// Set up optimizer
MeanVarianceOptimizer optimizer(ObjectiveType::MAX_SHARPE, 0.02);

// Configure constraints
OptimizationConstraints constraints;
constraints.min_weight = 0.0;
constraints.max_weight = 0.3;
constraints.long_only = true;
constraints.sum_to_one = true;

// Optimize
auto result = optimizer.optimize(
    data.get_mean_returns(),
    covariance,
    constraints
);

// Display results
result.print_summary();
```

### Efficient Frontier Computation

```cpp
#include "optimizer/efficient_frontier.hpp"

// Configure frontier
EfficientFrontier frontier;
frontier.set_num_points(20);
frontier.set_auto_range(true);

// Compute frontier
auto result = frontier.compute(
    expected_returns,
    covariance,
    constraints,
    0.02  // risk-free rate
);

// Export and display
result.export_to_csv("frontier.csv");
result.print_summary();

std::cout << "Max Sharpe Portfolio:\n";
std::cout << "  Return: " << result.max_sharpe_portfolio.expected_return * 100 << "%\n";
std::cout << "  Risk: " << result.max_sharpe_portfolio.volatility * 100 << "%\n";
std::cout << "  Sharpe: " << result.max_sharpe_portfolio.sharpe_ratio << "\n";
```

### Turnover Constraints

```cpp
// Current portfolio
Eigen::VectorXd current_weights(n_assets);
current_weights << 0.2, 0.3, 0.5;

// Limit position changes to 10% per asset
OptimizationConstraints constraints;
constraints.max_turnover = 0.10;

// Rebalance with turnover limit
auto result = optimizer.optimize(
    expected_returns,
    covariance,
    constraints,
    current_weights
);

// Check actual turnover
for (int i = 0; i < n_assets; ++i) {
    double change = std::abs(result.weights(i) - current_weights(i));
    std::cout << "Asset " << i << " turnover: " << change * 100 << "%\n";
}
```

## Configuration

Portfolio configuration is specified in JSON format:

```json
{
  "data": {
    "source": "csv",
    "filepath": "data/prices/market_data.csv",
    "tickers": ["AAPL", "GOOGL", "MSFT", "AMZN", "TSLA"],
    "start_date": "2020-01-01",
    "end_date": "2024-12-31"
  },
  "risk_model": {
    "type": "ledoit_wolf",
    "shrinkage_target": "constant_correlation",
    "bias_correction": true
  },
  "optimizer": {
    "type": "mean_variance",
    "objective": "max_sharpe",
    "constraints": {
      "min_weight": 0.0,
      "max_weight": 0.3,
      "long_only": true,
      "sum_to_one": true,
      "max_turnover": 0.15
    },
    "risk_free_rate": 0.02
  }
}
```

**Optimizer Types**:
- `mean_variance`: Markowitz optimization (implemented)
- `qubo`: Quantum-ready formulation (planned for 2.3)

**Objectives**:
- `min_variance`: Minimize portfolio volatility
- `target_return`: Achieve specific return with minimum risk
- `max_sharpe`: Maximize risk-adjusted return
- `risk_aversion`: Balance return vs risk with lambda parameter

**Constraints**:
- `min_weight`, `max_weight`: Box constraints per asset
- `long_only`: No short positions
- `sum_to_one`: Full investment (weights sum to 1)
- `max_turnover`: Maximum position change per asset

## Testing

```bash
# Build and run all tests
cd build
make -j$(nproc)
make test

# Run specific test suites
./bin/test_data_loader
./bin/test_risk_models
./bin/test_mean_variance_optimizer
./bin/test_efficient_frontier
./bin/test_turnover_constraints

# Run with verbose output
./bin/test_efficient_frontier -s

# Run specific test tags
./bin/test_efficient_frontier "[Basic]"
./bin/test_turnover_constraints "[Validation]"
```

**Test Coverage**:
- Data Layer: 8 tests
- Risk Models: 10 tests
- Mean-Variance Optimizer: 11 tests
- Efficient Frontier: 17 tests
- Turnover Constraints: 10 tests
- Extended Integration: 9 tests
- **Total**: 65 tests, 100% passing

## Project Structure

```
portfolio-optimizer/
├── include/                        # Public headers
│   ├── data/                       Complete (2 files)
│   │   ├── market_data.hpp
│   │   └── data_loader.hpp
│   ├── risk/                       Complete (5 files)
│   │   ├── risk_model.hpp
│   │   ├── sample_covariance.hpp
│   │   ├── ewma_covariance.hpp
│   │   ├── ledoit_wolf_shrinkage.hpp
│   │   └── risk_model_factory.hpp
│   └── optimizer/                  Complete (5 files)
│       ├── optimizer_interface.hpp
│       ├── mean_variance_optimizer.hpp
│       ├── efficient_frontier.hpp
│       ├── osqp_solver.hpp
│       └── quadratic_solver.hpp
├── src/                            # Implementation
│   ├── main.cpp
│   ├── data/                       Complete (2 files)
│   ├── risk/                       Complete (5 files)
│   └── optimizer/                  Complete (5 files)
├── tests/                          # Test suite
│   ├── test_data_loader.cpp        8 tests
│   ├── test_risk_models.cpp        10 tests
│   ├── test_mean_variance_optimizer.cpp  11 tests
│   ├── test_efficient_frontier.cpp       17 tests
│   └── test_turnover_constraints.cpp     10 tests
├── data/                           # Configuration and sample data
│   ├── config/
│   └── prices/
├── docs/                           # Documentation
└── CMakeLists.txt                  # Build configuration
```

## Project Status

| Component | Status | Completeness |
|-----------|--------|--------------|
| Data Layer | Complete | 100% |
| Risk Models | Complete | 100% |
| Build System | Complete | 100% |
| Testing Framework | Complete | 100% |
| Mean-Variance Optimizer | Complete | 100% |
| Efficient Frontier | Complete | 100% |
| Turnover Constraints | Complete | 100% |
| Advanced Constraints | Planned | 0% |
| Backtesting | Planned | 0% |
| Analytics | Planned | 0% |
| Visualization | Planned | 0% |

**Overall Progress**: 65% Complete (3.25/5 major feature sets)

## Recent Updates

### v1.3.0-beta - Feature Set 2.2 (February 2026)

**Completed**:
- Efficient frontier computation with dual methods (risk aversion, target return)
- Automatic feasible return range detection
- Special portfolio identification (min variance, max Sharpe)
- CSV export and summary statistics
- Per-asset turnover constraints across all optimization modes
- 27 new tests (17 frontier, 10 turnover), 100% passing
- Command-line integration with --frontier flag
- Performance: 20-point frontier in 20ms (10 assets)

### v1.2.0-alpha - Feature Set 2.1 (January 2026)

**Completed**:
- Mean-variance optimizer with OSQP integration
- Three optimization modes (min variance, target return, risk aversion)
- Max Sharpe ratio via grid search
- Comprehensive constraint handling (box, budget, long-only)
- 200x performance improvement over custom solver
- 11-test suite with convergence diagnostics

### v1.1.0 - Feature Set 1 (January 2026)

**Completed**:
- Risk model estimation framework
- Sample covariance, EWMA, Ledoit-Wolf shrinkage estimators
- Risk model factory pattern
- 30+ risk model tests

## Performance

| Operation | Time | Details |
|-----------|------|---------|
| CSV loading | 5ms | 1000 rows, 10 assets |
| Return calculation | 1ms | 1000 periods |
| Covariance estimation | 2-5ms | Sample/EWMA/Ledoit-Wolf |
| Portfolio optimization | <1ms | Single optimization, 2-50 assets |
| Efficient frontier | 20ms | 20 points, 10 assets |
| Complete pipeline | <150ms | Data load + risk + frontier |

All performance targets exceeded.

## Roadmap

### Feature Set 2.3: Advanced Constraints and Optimization (Planned)

**Scope**:
- Direct Max Sharpe optimization (fractional programming)
- Sector and group constraints (industry exposure limits)
- Tracking error constraints (vs benchmark)
- Cardinality constraints (max number of positions)

**Status**: Design phase
**Estimated Effort**: 2-3 weeks

### Feature Set 3: Backtesting Engine (Planned)

**Components**:
```cpp
namespace portfolio::backtest {
    class Portfolio;                // Position tracking
    class BacktestEngine;           // Walk-forward simulation
    class RebalanceScheduler;       // Daily/weekly/monthly
    class TransactionCostModel;     // Commissions, slippage
    class TradeLogger;              // Trade history
}
```

**Key Features**:
- Walk-forward backtesting with rolling windows
- Transaction cost modeling (commissions, slippage, market impact)
- Flexible rebalancing schedules (calendar and threshold-based)
- Portfolio state tracking (positions, cash, P&L)
- Trade attribution and analysis

**Status**: Requirements gathering
**Estimated Effort**: 4-6 weeks

### Feature Set 4: Performance Analytics (Planned)

**Metrics**:
- Return metrics: Total return, CAGR, annualized returns
- Risk metrics: Volatility, max drawdown, VaR, CVaR, downside deviation
- Risk-adjusted: Sharpe, Sortino, Calmar, Information ratios
- Relative: Alpha, beta, tracking error, R-squared
- Attribution: Brinson-Fachler decomposition

**Status**: Planned
**Estimated Effort**: 3-4 weeks

### Feature Set 5: Visualization and Reporting (Planned)

**Tools**:
- Python integration for plotting (matplotlib/plotly)
- HTML report generation
- Equity curve visualization
- Weight evolution charts
- Risk decomposition plots
- Efficient frontier plots

**Status**: Planned
**Estimated Effort**: 3-4 weeks

### Feature Set 6: Quantum: Quantum Optimization (Planned)

**Key Features**
- Portfolio optimization as QUBO problem
- Portfolio optimization as QUOPLE problem
- Portfolio optimization as QAOA problem
- Portfolio optimization as QAMOO problem

## Code Statistics

```
Production Code:     ~8,500 lines C++17
Test Code:           ~2,800 lines
Build Config:        ~300 lines (CMake, Docker)
Documentation:       ~6,000 lines (README, guides, release notes)
Total:               ~17,600 lines

Test Coverage:
  Unit tests:        45 tests
  Integration tests: 15 tests
  Convergence tests: 5 tests
  Total:            65 tests (100% passing)
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Write tests for new functionality
4. Ensure code compiles with no warnings (-Wall -Wextra -Wpedantic)
5. Update documentation
6. Submit a pull request

**Current Priority Areas**:
- Direct Max Sharpe optimization implementation
- Sector constraint formulation
- Backtesting engine components
- Performance optimizations
- Documentation improvements

## License

MIT License - see LICENSE file for details

## References

1. Markowitz, H. (1952). "Portfolio Selection", The Journal of Finance
2. Stellato, B. et al. (2020). "OSQP: An Operator Splitting Solver for Quadratic Programs", Mathematical Programming Computation
3. Boyd, S. & Vandenberghe, L. (2004). Convex Optimization, Cambridge University Press
4. Ledoit, O. & Wolf, M. (2004). "Honey, I Shrunk the Sample Covariance Matrix", Journal of Portfolio Management

## Contact

For questions, issues, or contributions, please visit the GitHub repository:
https://github.com/darkhorse286/portfolio-optimizer

## Version

Current Version: 1.3.0
Last Updated: February 4, 2026
Next Milestone: Feature Set 2.3 - Advanced Constraints and Direct Optimization