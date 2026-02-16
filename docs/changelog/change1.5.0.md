# Changelog Entry for v1.5.0

## [1.5.0] - 2026-02-16

### Added

- Backtesting engine layer (`namespace portfolio::backtest`) with five production components
- `Portfolio` class: position tracking, cash management, NAV computation, target weight rebalancing with sell-first execution, cost basis tracking, daily snapshot recording, and portfolio reset
- `TransactionCostModel` class: configurable commission (proportional), slippage (basis points), and market impact (none, linear, sqrt) with single-order and batch calculation interfaces
- `RebalanceScheduler` class: calendar-based triggers (daily, weekly, monthly, quarterly, annually), drift-threshold triggers, minimum-days-between cooldown, JSON and string configuration
- `TradeLogger` class: per-trade and per-rebalance logging with auto-increment IDs, zero-share skip, turnover computation, date and ticker queries, summary statistics, and CSV export with automatic directory creation
- `BacktestEngine` class: walk-forward simulation integrating risk estimation, optimization, rebalance scheduling, transaction costs, and trade logging; two `run()` overloads (default Max Sharpe and user-provided optimizer)
- `BacktestResult` struct: total return, annualized return, annualized volatility, Sharpe ratio, max drawdown, NAV series CSV export
- `BacktestParams` struct: configuration from JSON via `from_config()` or programmatic construction
- Supporting structs: `Position`, `PortfolioSnapshot`, `TradeOrder`, `TradeCost`, `TransactionCostConfig`, `TradeRecord`, `TradeSummary`, `RebalanceConfig`, `RebalanceFrequency` enum
- Test suite: 22 new test cases across 5 test files (`test_portfolio.cpp`, `test_transaction_costs.cpp`, `test_rebalance_scheduler.cpp`, `test_trade_logger.cpp`, `test_backtest_engine.cpp`)
- Release notes: `docs/release_notes_feature_set_3.md`
- Changelog entry: `docs/changelog/change1.5.0.md`
- Updated `README.md` with backtest usage examples, updated architecture diagram, updated roadmap including quantum optimization feature set, and updated code statistics

### Changed

- Updated namespace structure documentation to reflect all four completed layers (data, risk, optimizer, backtest)
- Updated project status from 60% to 70% complete
- Updated test count from 65 to 87 (100% passing)
- Updated code statistics: production code ~10,500 lines, test code ~3,800 lines, total ~22,650 lines
- Revised roadmap: Feature Set 3 moved from Planned to Completed; Feature Set 6 (Quantum Optimization) expanded with QUBO, QAOA, QAMO, and QAMOO algorithm descriptions and benchmarking plan

### Fixed

- Corrected brace mismatch in `BacktestEngine::run()` where a missing closing brace for the disabled rebalance block caused the snapshot-recording code to be trapped inside `if (false)` and the `return result` statement to execute inside the for loop, resulting in zero snapshots and early loop termination

### Performance

- Portfolio update_prices: <0.01ms for 10 assets
- Portfolio set_target_weights: <0.01ms for 10 assets
- Transaction cost calculation: <0.01ms for 10 trades
- Full backtest (no optimization): <5ms for 60 days, 3 assets
- Full backtest (with optimization): <100ms for 252 days, 10 assets
- All Feature Set 3 performance targets met

---

## [1.4.0] - 2026-02-12

### Added

- Tracking error soft-constraint (penalty method) for benchmark-relative portfolios
- Integration of tracking-error penalty into all mean-variance objectives (min variance, target return, risk aversion, max Sharpe direct)
- Auto-tuning heuristic for tracking-error penalty coefficient
- Direct Max Sharpe via Schaible's fractional programming refinement
- Sector / group constraints: `SectorMapping` helper and `GroupConstraint` support
- Inequality constraint support: `A_ineq`, `b_ineq_lower`, `b_ineq_upper` in `QuadraticProblem`
- Unit tests for tracking-error behavior and validation (`test_tracking_error.cpp`)

### Changed

- `OptimizationConstraints::validate()` enforces positive `max_tracking_error` when benchmark_weights provided
- Validation warns when benchmark weights do not sum to 1.0 (tolerance 0.01)

---

## [1.3.0] - 2026-02-04

### Added

- Efficient frontier computation with dual methods (risk aversion sweep and target return sweep)
- Automatic feasible return range detection
- Special portfolio identification (minimum variance, maximum Sharpe)
- Per-asset turnover constraints via tightened box bounds
- Frontier data export to CSV
- 27 new tests: efficient frontier (17) and turnover constraints (10)

### Changed

- Updated README with Feature Set 2.2 completion status
- Deferred direct Max Sharpe and sector constraints to Feature Set 2.3

### Fixed

- Automatic output directory creation prevents file write errors
- Frontier point sorting ensures proper ordering
- CSV export creates parent directories automatically

---

## [1.2.0] - 2026-01-30

### Added

- Mean-variance portfolio optimizer with OSQP quadratic programming solver
- Four optimization objectives: min variance, target return, risk aversion, max Sharpe (grid search)
- Box constraints (min/max weight), budget constraint (sum to one), long-only constraint
- OSQP solver wrapper with CSC matrix conversion
- 11 tests with convergence diagnostics
- OSQP dependency integration in CMake

### Changed

- Updated project structure with optimizer namespace
- Added optimizer section to portfolio configuration JSON

---

## [1.1.0] - 2026-01-19

### Added

- Risk model estimation framework with abstract `RiskModel` interface
- Sample covariance estimator with optional Bessel's correction
- EWMA covariance estimator with configurable decay
- Ledoit-Wolf shrinkage estimator with multiple target matrices
- Risk model factory pattern for configuration-driven creation
- 30+ test cases for risk models
- Example program for risk estimators

### Changed

- Simplified `risk_model` configuration (flattened nested structure)

### Removed

- Conflicting `method` field from risk model configuration
- Unused `min_periods` field
- Nested `shrinkage` object (flattened to top level)

---

## [1.0.0] - 2026-01-17

### Added

- Initial project structure and CMake build system
- Data layer with `MarketData` and `DataLoader` classes
- CSV and JSON data loading
- Return calculation (simple and logarithmic)
- Date filtering and asset selection
- Missing data handling (drop and forward-fill)
- Statistical methods (mean, volatility, correlation)
- Synthetic data generation for testing
- Comprehensive test suite for data layer
- Docker containerization
- VS Code development environment configuration

---

## Version History

- **1.5.0** (2026-02-16) - Backtesting Engine (Feature Set 3)
- **1.4.0** (2026-02-12) - Advanced Constraints and Optimization (Feature Set 2.3)
- **1.3.0** (2026-02-04) - Advanced Portfolio Optimization (Feature Set 2.2)
- **1.2.0** (2026-01-30) - Mean-Variance Optimization (Feature Set 2.1)
- **1.1.0** (2026-01-19) - Risk Model Estimation (Feature Set 1)
- **1.0.0** (2026-01-17) - Foundation and Data Layer

---

## Upcoming Releases

### [1.6.0] - Planned (March 2026)

**Feature Set 4: Performance Analytics**

- Return metrics: total return, annualized return, CAGR
- Risk metrics: volatility, max drawdown, VaR, CVaR, downside deviation
- Risk-adjusted metrics: Sharpe, Sortino, Calmar, Information ratio
- Relative performance: alpha, beta, tracking error, R-squared
- Attribution analysis: Brinson-Fachler decomposition
- Rolling statistics computation

### [1.7.0] - Planned (April 2026)

**Feature Set 5: Visualization and Reporting**

- Python integration for plotting (matplotlib/plotly)
- Equity curve visualization
- Weight evolution charts
- Drawdown and underwater plots
- Rolling metrics visualization
- Efficient frontier plots
- HTML report generation

### [2.0.0] - Planned (Q2-Q3 2026)

**Feature Set 6: Quantum Optimization**

- QUBO formulation for portfolio optimization with classical solvers (simulated annealing, tabu search)
- Quantum annealer interface (D-Wave) for QUBO problems
- QAOA (Quantum Approximate Optimization Algorithm) implementation
- QAMO (Quantum Alternating Mean-field Optimization) implementation
- QAMOO (Quantum Alternating Multi-Objective Optimization) for Pareto frontier generation
- Benchmarking framework comparing classical QP solvers against quantum and quantum-inspired approaches
- Solution quality, convergence speed, and scalability analysis across identical problem instances

---

## Dependencies

### Runtime Dependencies
- **Eigen** 3.3+ (linear algebra)
- **nlohmann/json** 3.2.0+ (configuration)
- **OSQP** 0.6.3+ (quadratic programming)

### Build Dependencies
- **CMake** 3.15+
- **GCC** 13+ or **Clang** 14+
- **Catch2** v3 (testing, auto-downloaded)

---

## Notes

- All releases maintain backward compatibility within major versions
- Configuration format changes are documented in migration guides
- Performance benchmarks measured on Intel i7 @ 3.5 GHz (WSL2 Ubuntu 24.04)
- Test coverage targets >90% for all public APIs
- The 2.0.0 release will be a major version to reflect the introduction of quantum solver interfaces alongside classical solvers

## Links

[1.5.0]: https://github.com/darkhorse286/portfolio-optimizer/compare/v1.4.0...v1.5.0
[1.4.0]: https://github.com/darkhorse286/portfolio-optimizer/compare/v1.3.0...v1.4.0
[1.3.0]: https://github.com/darkhorse286/portfolio-optimizer/compare/v1.2.0...v1.3.0
[1.2.0]: https://github.com/darkhorse286/portfolio-optimizer/compare/v1.1.0...v1.2.0
[1.1.0]: https://github.com/darkhorse286/portfolio-optimizer/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/darkhorse286/portfolio-optimizer/releases/tag/v1.0.0