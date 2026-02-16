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
- **Build System**: CMake with Docker support
- **Testing**: Comprehensive test suite with Catch2 (87+ tests, 100% passing)

### Feature Set 3: Backtesting Engine (Completed)

- Walk-forward backtesting engine with configurable lookback windows
- Portfolio state management with position tracking, cash accounting, and NAV computation
- Transaction cost modeling: commissions, slippage (bps), and market impact (linear, sqrt, none)
- Flexible rebalancing: calendar-based (daily, weekly, monthly, quarterly, annually) and drift-threshold triggers
- Minimum days between rebalances to prevent excessive trading
- Trade logging with per-trade and per-rebalance recording, CSV export, and summary statistics
- Backtest result analytics: total return, annualized return, annualized volatility, Sharpe ratio, max drawdown
- NAV series export to CSV for external analysis
- Integration with risk model and optimizer layers for end-to-end backtesting

### Planned

- Feature Set 4: Performance analytics and attribution
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
│   └── backtest/                # Backtesting engine (complete)
│       ├── backtest_engine.hpp
│       ├── portfolio.hpp
│       ├── rebalance_scheduler.hpp
│       ├── trade_logger.hpp
│       └── transaction_cost_model.hpp
├── src/
│   ├── data/
│   ├── risk/
│   ├── optimizer/
│   └── backtest/
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
make test

# Run optimizer
./bin/portfolio_optimizer --config ../data/config/portfolio_config.json
```

### Docker Build

```bash
docker build -t portfolio-optimizer .
docker run -it portfolio-optimizer
```

## Usage

### Backtesting Example

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

// Analyze results
result.print_summary();
std::cout << "Total return: " << result.total_return() << "\n";
std::cout << "Sharpe ratio: " << result.sharpe_ratio() << "\n";
std::cout << "Max drawdown: " << result.max_drawdown() << "\n";

// Export results
result.export_nav_to_csv("results/nav_series.csv");
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

// Inspect state
std::cout << "NAV: " << portfolio.nav() << "\n";
std::cout << "Cash: " << portfolio.cash() << "\n";
auto weights = portfolio.current_weights();
```

### Transaction Cost Modeling Example

```cpp
#include "backtest/transaction_cost_model.hpp"

using namespace portfolio::backtest;

TransactionCostConfig cfg;
cfg.commission_rate = 0.001;       // 10 bps
cfg.slippage_bps = 5.0;           // 5 bps
cfg.market_impact_model = "sqrt";  // Sub-linear impact
cfg.market_impact_coeff = 1e-4;

TransactionCostModel model(cfg);

TradeOrder order{"AAPL", 1000.0, 150.0};
TradeCost cost = model.calculate_cost(order);
// cost.commission, cost.slippage, cost.market_impact, cost.total()
```

### Rebalance Scheduling Example

```cpp
#include "backtest/rebalance_scheduler.hpp"

using namespace portfolio::backtest;

RebalanceConfig cfg;
cfg.frequency = RebalanceFrequency::MONTHLY;
cfg.drift_threshold = 0.05;    // 5% drift triggers rebalance
cfg.min_days_between = 5;      // At least 5 days between rebalances

RebalanceScheduler scheduler(cfg);

Eigen::VectorXd current(3), target(3);
current << 0.55, 0.25, 0.20;
target << 0.40, 0.30, 0.30;

if (scheduler.should_rebalance("2025-02-03", current, target)) {
    // Execute rebalance
    scheduler.record_rebalance("2025-02-03");
}
```

### Optimization Example

```cpp
#include "optimizer/mean_variance_optimizer.hpp"
#include "risk/risk_model_factory.hpp"

using namespace portfolio;

// Estimate risk
risk::RiskModelConfig risk_cfg;
risk_cfg.type = "ledoit_wolf";
auto model = risk::RiskModelFactory::create(risk_cfg);
auto cov = model->estimate_covariance(returns);

// Optimize
optimizer::MeanVarianceOptimizer opt(optimizer::ObjectiveType::MAX_SHARPE, 0.02);
optimizer::OptimizationConstraints constraints;
constraints.min_weight = 0.0;
constraints.max_weight = 0.3;

auto result = opt.optimize(expected_returns, cov, constraints);
```

### Turnover-Constrained Optimization Example

```cpp
optimizer::OptimizationConstraints constraints;
constraints.min_weight = 0.0;
constraints.max_weight = 0.30;
constraints.max_turnover = 0.10; // 10% per-asset limit

auto result = optimizer.optimize(
    expected_returns,
    covariance,
    constraints,
    current_weights
);
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
  },
  "backtest": {
    "initial_capital": 1000000.0,
    "rebalance_frequency": "monthly",
    "lookback_window": 252,
    "min_history": 60,
    "transaction_costs": {
      "commission_rate": 0.001,
      "slippage_bps": 5.0,
      "market_impact": "linear"
    },
    "cash_management": {
      "allow_cash": true,
      "cash_return": 0.0
    }
  }
}
```

## Testing

```bash
# Build and run all tests
cd build
cmake ..
make -j$(nproc)
make test

# Run specific test suites
./bin/test_data_loader
./bin/test_risk_models
./bin/test_mean_variance_optimizer
./bin/test_efficient_frontier
./bin/test_turnover_constraints
./bin/test_tracking_error
./bin/test_test_portfolio
./bin/test_test_transaction_costs
./bin/test_test_rebalance_scheduler
./bin/test_test_trade_logger
./bin/test_test_backtest_engine

# Run with verbose output
./bin/test_test_backtest_engine -s
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

## Roadmap

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
Production Code:     ~10,500 lines C++17
Test Code:           ~3,800 lines
Build Config:        ~350 lines (CMake, Docker)
Documentation:       ~8,000 lines (README, guides, release notes, changelogs)
Total:               ~22,650 lines

Test Coverage:
  Unit tests:        55 tests
  Integration tests: 22 tests
  Convergence tests: 5 tests
  Edge case tests:   5 tests
  Total:             87 tests (100% passing)
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
4. Stellato, B. et al. (2020). "OSQP: An Operator Splitting Solver for Quadratic Programs"
5. Boyd, S. & Vandenberghe, L. (2004). Convex Optimization, Cambridge University Press
6. Farhi, E. et al. (2014). "A Quantum Approximate Optimization Algorithm"
7. Lucas, A. (2014). "Ising formulations of many NP problems", Frontiers in Physics