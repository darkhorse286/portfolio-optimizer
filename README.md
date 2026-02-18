# Portfolio Optimizer

A high-performance C++ portfolio optimization library implementing modern portfolio theory with production-grade solvers. Designed for extensibility toward quantum optimization algorithms.

## Features

### Completed

- **Data Layer**: CSV/JSON loading, return calculation, statistical analysis
- **Risk Models**: Sample covariance, EWMA, Ledoit-Wolf shrinkage
- **Mean-Variance Optimization**: Markowitz portfolio optimization with OSQP solver
- **Efficient Frontier**: Complete Pareto-optimal frontier computation
- **Turnover Constraints**: Per-asset position change limits
- **Tracking Error (soft penalty)**: Benchmark-relative penalty-based constraint
- **Direct Max Sharpe (Schaible's refinement)**: Fractional programming refinement for direct Max-Sharpe
- **Sector / Group Constraints**: `SectorMapping` helper and `GroupConstraint` support
- **Inequality Constraints**: Partial support for inequality constraints via `A_ineq` in QP
- **Backtesting Engine**: Walk-forward simulation with transaction costs, rebalancing, and trade logging
- **Performance Analytics**: Comprehensive return, risk, risk-adjusted, benchmark, rolling, drawdown, and attribution analysis
- **Build System**: CMake with Docker support
- **Testing**: Comprehensive test suite with Catch2 (157+ tests, 100% passing)

### Feature Set 4: Performance Analytics (Completed)

- Core metrics: total return, CAGR, rolling returns, monthly return table, annualized volatility, downside deviation, skewness, kurtosis
- Risk metrics: maximum drawdown with full event detail, Value at Risk (historical and parametric), Conditional VaR / Expected Shortfall
- Risk-adjusted metrics: Sharpe ratio, Sortino ratio, Calmar ratio, Omega ratio
- Drawdown analysis: discrete event identification (peak-trough-recovery cycles), top-N extraction, aggregate statistics, underwater curve
- Benchmark analysis: CAPM regression (alpha, beta, R-squared), tracking error, information ratio, active return, up/down capture ratios
- Rolling statistics: rolling window computations for volatility, Sharpe, Sortino, max drawdown, beta, tracking error, skewness, kurtosis, plus generic `apply()` for custom functions
- Performance attribution: Brinson-Fachler single-period decomposition (allocation, selection, interaction effects), multi-period Carino logarithmic linking
- Integration with `BacktestResult` via `compute_analytics()` for seamless end-to-end analysis
- Export to JSON and CSV from both standalone analytics and backtest results

### Planned

- Feature Set 5: Visualization and reporting
- Feature Set 6: Quantum optimization (classical vs quantum algorithm comparison)

## Architecture

```
portfolio-optimizer/
├── include/
│   ├── data/                    # Data layer (complete)
│   │   ├── market_data.hpp
│   │   └── data_loader.hpp
│   ├── risk/                    # Risk models (complete)
│   │   ├── risk_model.hpp
│   │   ├── sample_covariance.hpp
│   │   ├── ewma_covariance.hpp
│   │   ├── ledoit_wolf_shrinkage.hpp
│   │   └── risk_model_factory.hpp
│   ├── optimizer/               # Portfolio optimization (complete)
│   │   ├── optimizer_interface.hpp
│   │   ├── mean_variance_optimizer.hpp
│   │   ├── efficient_frontier.hpp
│   │   ├── osqp_solver.hpp
│   │   └── quadratic_solver.hpp
│   ├── backtest/                # Backtesting engine (complete)
│   │   ├── backtest_engine.hpp
│   │   ├── backtest_result.hpp
│   │   ├── portfolio.hpp
│   │   ├── rebalance_scheduler.hpp
│   │   ├── trade_logger.hpp
│   │   └── transaction_cost_model.hpp
│   └── analytics/               # Performance analytics (complete)
│       ├── performance_metrics.hpp
│       ├── drawdown_analysis.hpp
│       ├── benchmark_analysis.hpp
│       ├── rolling_statistics.hpp
│       └── attribution.hpp
├── src/
│   ├── data/
│   ├── risk/
│   ├── optimizer/
│   ├── backtest/
│   └── analytics/
├── tests/
├── data/config/
├── docs/
└── CMakeLists.txt
```

### Namespace Structure

```cpp
namespace portfolio {
    namespace data { }        // Data loading and market data management
    namespace risk { }        // Risk model estimation
    namespace optimizer { }   // Portfolio optimization
    namespace backtest { }    // Backtesting engine
    namespace analytics { }   // Performance analytics
}
```

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
make run_tests

# Run optimizer
./bin/portfolio_optimizer --config ../data/config/portfolio_config.json
```

### Docker Build

```bash
docker build -t portfolio-optimizer .
docker run -it portfolio-optimizer
```

## Usage

### End-to-End Backtest with Analytics

```cpp
#include "backtest/backtest_engine.hpp"
#include "data/data_loader.hpp"

using namespace portfolio;
using namespace portfolio::backtest;

// Load market data
auto data = DataLoader::load_csv("data/prices/market_data.csv");

// Configure backtest parameters
BacktestParams params;
params.initial_capital = 1000000.0;
params.lookback_window = 252;
params.min_history = 60;
params.risk_model_type = "ewma";
params.risk_free_rate = 0.02;

// Configure transaction costs
params.transaction_costs.commission_rate = 0.001;
params.transaction_costs.slippage_bps = 5.0;
params.transaction_costs.market_impact_model = "linear";
params.transaction_costs.market_impact_coeff = 0.002;

// Configure rebalancing
params.rebalance.frequency = RebalanceFrequency::MONTHLY;
params.rebalance.drift_threshold = 0.05;
params.rebalance.min_days_between = 5;

// Run backtest
BacktestEngine engine(params);
auto result = engine.run(data);

// Full analytics summary (uses analytics layer automatically)
result.print_summary();

// Access the complete analytics suite
auto metrics = result.compute_analytics();
std::cout << "Sharpe Ratio:  " << metrics.sharpe_ratio() << "\n";
std::cout << "Sortino Ratio: " << metrics.sortino_ratio() << "\n";
std::cout << "Calmar Ratio:  " << metrics.calmar_ratio() << "\n";
std::cout << "Max Drawdown:  " << metrics.max_drawdown() << "\n";
std::cout << "VaR (95%):     " << metrics.value_at_risk(0.95) << "\n";
std::cout << "CVaR (95%):    " << metrics.conditional_var(0.95) << "\n";
std::cout << "Skewness:      " << metrics.skewness() << "\n";
std::cout << "Kurtosis:      " << metrics.kurtosis() << "\n";

// Export analytics
result.export_nav_to_csv("results/nav_series.csv");
result.export_analytics_csv("results/analytics_timeseries.csv");
std::string json = result.export_analytics_json();
```

### Standalone Performance Analytics

```cpp
#include "analytics/performance_metrics.hpp"
#include "analytics/drawdown_analysis.hpp"
#include "analytics/benchmark_analysis.hpp"
#include "analytics/rolling_statistics.hpp"
#include "analytics/attribution.hpp"

using namespace portfolio::analytics;

// From raw vectors (no backtest needed)
std::vector<double> nav = { /* daily NAV values */ };
std::vector<double> returns = { /* daily returns */ };
std::vector<std::string> dates = { /* YYYY-MM-DD dates */ };

PerformanceMetrics metrics(nav, returns, dates, 0.02, 252);
std::cout << metrics.summary();

// Drawdown event analysis
DrawdownAnalysis dd(nav, dates);
auto worst = dd.worst_drawdown();
auto top5 = dd.top_drawdowns(5);
std::cout << dd.report();

// Benchmark-relative analysis
std::vector<double> benchmark_returns = { /* benchmark daily returns */ };
BenchmarkAnalysis bench(returns, benchmark_returns, 0.02, 252);
std::cout << "Alpha: " << bench.alpha() << "\n";
std::cout << "Beta:  " << bench.beta() << "\n";
std::cout << "R^2:   " << bench.r_squared() << "\n";
std::cout << "IR:    " << bench.information_ratio() << "\n";

// Rolling statistics
RollingStatistics rolling(returns, 63);  // 63-day (quarterly) window
auto vol = rolling.volatility();
auto sharpe = rolling.sharpe_ratio();
auto beta = rolling.beta(benchmark_returns);

// Performance attribution (Brinson-Fachler)
std::vector<SectorAllocation> port_sectors = {
    {"Technology", 0.40, 0.06},
    {"Healthcare", 0.35, 0.03},
    {"Finance",    0.25, -0.01}
};
std::vector<SectorAllocation> bench_sectors = {
    {"Technology", 0.30, 0.05},
    {"Healthcare", 0.40, 0.02},
    {"Finance",    0.30, 0.01}
};
auto attrib = Attribution::single_period(port_sectors, bench_sectors);
// attrib.total_allocation, attrib.total_selection, attrib.total_interaction
```

### Portfolio Management Example

```cpp
#include "backtest/portfolio.hpp"

using namespace portfolio::backtest;

std::vector<std::string> tickers = {"AAPL", "GOOGL", "MSFT"};
Portfolio portfolio(1000000.0, tickers);

// Update prices
Eigen::VectorXd prices(3);
prices << 150.0, 140.0, 380.0;
portfolio.update_prices("2025-01-15", prices);

// Set target allocation
Eigen::VectorXd target(3);
target << 0.4, 0.3, 0.3;
auto shares_traded = portfolio.set_target_weights(target, prices);
```

## Testing

```bash
# Build and run all tests
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)
make run_tests

# Run individual test suites
./bin/test_test_data_loader
./bin/test_test_risk_models
./bin/test_test_mean_variance_optimizer
./bin/test_test_efficient_frontier
./bin/test_test_turnover_constraints
./bin/test_test_tracking_error
./bin/test_test_portfolio
./bin/test_test_transaction_costs
./bin/test_test_rebalance_scheduler
./bin/test_test_trade_logger
./bin/test_test_backtest_engine
./bin/test_test_performance_metrics
./bin/test_test_drawdown_analysis
./bin/test_test_benchmark_analysis
./bin/test_test_rolling_statistics
./bin/test_test_attribution

# Run with verbose output
./bin/test_test_performance_metrics -s
```

## Performance Benchmarks

| Operation | Problem Size | Time |
|-----------|--------------|------|
| CSV Load | 500 dates, 10 assets | ~5ms |
| Return Calculation | 500 dates, 10 assets | <1ms |
| Sample Covariance | 100x100 matrix | ~2ms |
| EWMA Covariance | 100x100 matrix | ~3ms |
| Ledoit-Wolf Shrinkage | 100x100 matrix | ~5ms |
| Min Variance Optimization | 10 assets | <1ms |
| Max Sharpe Optimization | 10 assets | ~25ms |
| Efficient Frontier (20 pts) | 10 assets | ~20ms |
| Risk Model Estimation | 100x100 matrix | <5ms |
| Portfolio Optimization | 50 assets | <50ms |
| Backtest (2 years, monthly) | 10 assets | <100ms |
| All analytics metrics | 252-day series | <2ms |
| Rolling Sharpe (252-day window) | 1000 days | <5ms |
| Drawdown analysis | 2500 days | <2ms |
| Full analytics report | 252-day series | <10ms |

## Roadmap

### Feature Set 5: Visualization and Reporting (Planned)

**Tools**:
- Python integration for plotting (matplotlib/plotly)
- HTML report generation
- Equity curve visualization
- Weight evolution charts
- Risk decomposition plots
- Efficient frontier plots
- Drawdown and underwater charts
- Rolling metrics visualization

**Status**: Planned
**Estimated Effort**: 3-4 weeks

### Feature Set 6: Quantum Optimization (Planned)

**Scope**: Compare current classical optimizations against emerging quantum and quantum-inspired algorithms for portfolio optimization. The goal is to benchmark solution quality, convergence speed, and scalability across classical and quantum approaches on identical problem instances.

**Algorithms under evaluation**:
- QUBO (Quadratic Unconstrained Binary Optimization): Reformulate the portfolio problem as a binary optimization suitable for quantum annealers and classical QUBO solvers (simulated annealing, tabu search). Compare against the Markowitz QP baseline.
- QAOA (Quantum Approximate Optimization Algorithm): Gate-model quantum algorithm for combinatorial optimization. Evaluate circuit depth, parameter landscape, and approximation quality versus classical solutions.
- QAMO (Quantum Alternating Mean-field Optimization): Variational mean-field approach to portfolio selection. Assess convergence behavior and solution quality relative to EWMA-based classical pipelines.
- QAMOO (Quantum Alternating Multi-Objective Optimization): Multi-objective extension for simultaneously optimizing return, risk, and other portfolio objectives. Compare Pareto frontier quality against classical efficient frontier computation.

**Planned deliverables**:
- Classical QUBO solver (simulated annealing) as a reference baseline
- Quantum annealer interface (D-Wave) for QUBO problems
- QAOA implementation targeting gate-model simulators
- QAMO and QAMOO implementations with variational circuits
- Side-by-side benchmarking framework: solution quality, wall-clock time, and scaling characteristics
- Penalty-based constraint encoding compatible with the existing tracking-error penalty infrastructure

**Status**: Research and design
**Estimated Effort**: 8-12 weeks

## Code Statistics

```
Production Code:     ~13,900 lines C++17
Test Code:           ~5,750 lines
Build Config:        ~350 lines (CMake, Docker)
Documentation:       ~9,500 lines (README, guides, release notes, changelogs)
Total:               ~29,500 lines

Test Coverage:
  Unit tests:        85 tests
  Integration tests: 52 tests
  Convergence tests: 5 tests
  Edge case tests:   15 tests
  Total:             157 tests (100% passing)
```

## Dependencies

### Runtime Dependencies
- **Eigen** 3.3+ (linear algebra)
- **nlohmann/json** 3.2.0+ (configuration)
- **OSQP** 0.6.3+ (quadratic programming)

### Build Dependencies
- **CMake** 3.15+
- **GCC** 13+ or **Clang** 14+
- **Catch2** v3 (testing, auto-downloaded)

### Development Dependencies
- **VS Code** with C/C++ extension (recommended)
- **gdb** (debugging)
- **Python** 3.8+ (visualization, planned)

## Contributing

Contributions are welcome. Please:

1. Fork the repository
2. Create a feature branch
3. Write tests for new functionality
4. Ensure code compiles with no warnings (`-Wall -Wextra -Wpedantic`)
5. Update documentation
6. Submit a pull request

## Version History

- **1.6.0** (2026-02-17) - Performance Analytics (Feature Set 4)
- **1.5.0** (2026-02-16) - Backtesting Engine (Feature Set 3)
- **1.4.0** (2026-02-12) - Advanced Constraints and Optimization (Feature Set 2.3)
- **1.3.0** (2026-02-04) - Advanced Portfolio Optimization (Feature Set 2.2)
- **1.2.0** (2026-01-30) - Mean-Variance Optimization (Feature Set 2.1)
- **1.1.0** (2026-01-19) - Risk Model Estimation (Feature Set 1)
- **1.0.0** (2026-01-17) - Foundation and Data Layer

## License

MIT License. See LICENSE for details.

## References

1. Markowitz, H. (1952). "Portfolio Selection", The Journal of Finance
2. J.P. Morgan (1996). "RiskMetrics Technical Document"
3. Ledoit, O. & Wolf, M. (2004). "Honey, I Shrunk the Sample Covariance Matrix"
4. Brinson, G., Hood, L. & Beebower, G. (1986). "Determinants of Portfolio Performance", Financial Analysts Journal
5. Carino, D. (1999). "Combining Attribution Effects Over Time", Journal of Performance Measurement
6. Beasley, J., Springer, S. & Moro, B. (1996). "Algorithm for the Percentage Points of the Normal Distribution"
7. Stellato, B. et al. (2020). "OSQP: An Operator Splitting Solver for Quadratic Programs"