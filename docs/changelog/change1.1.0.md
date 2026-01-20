# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2026-01-19

### Added

- Risk model estimation framework with abstract `RiskModel` interface
- Sample covariance estimator with optional Bessel's correction
- EWMA (Exponentially Weighted Moving Average) covariance estimator with configurable decay
- Ledoit-Wolf shrinkage estimator with multiple target matrices
- Risk model factory pattern for configuration-driven creation
- Comprehensive test suite for all risk models (30+ test cases)
- Example program demonstrating all three risk estimators
- Feature Set 1 documentation (`docs/FEATURE_SET_1_RISK_MODELS.md`)
- Release notes (`RELEASE_NOTES.md`)
- VS Code IntelliSense configuration (`.vscode/c_cpp_properties.json`)
- Risk model configuration schema in `portfolio_config.json`

### Changed

- Updated `CMakeLists.txt` to include risk model library and tests
- Simplified `risk_model` configuration section (flattened nested structure)
- Updated README.md with risk model information and usage examples
- Updated project status to 40% complete (2/5 feature sets)

### Removed

- Removed conflicting `method` field from risk model configuration
- Removed unused `min_periods` field from risk model configuration
- Removed nested `shrinkage` object in configuration (flattened to top level)

### Fixed

- Resolved configuration conflicts between `type` and `method` fields
- Improved configuration clarity with single source of truth for risk model selection

### Performance

- Sample covariance: ~2ms for 100x100 matrix
- EWMA covariance: ~3ms for 100x100 matrix
- Ledoit-Wolf shrinkage: ~5ms for 100x100 matrix
- All performance targets met

## [1.0.0] - 2026-01-17

### Added

- Initial project structure and build system
- Data layer with `MarketData` and `DataLoader` classes
- CSV and JSON data loading capabilities
- Return calculation (simple and logarithmic)
- Date filtering and asset selection
- Missing data handling (drop and forward-fill)
- Statistical methods (mean, volatility, correlation)
- Synthetic data generation for testing
- Comprehensive test suite for data layer
- Docker containerization
- VS Code development environment configuration
- CMake build system with Eigen and nlohmann/json integration
- Catch2 testing framework integration
- Documentation (README, build guides)
- Example data generation script

### Performance

- CSV loading: ~5ms for typical dataset
- Return calculation: <1ms
- Covariance computation: ~2ms

---

## Version History

- **1.1.0** (2026-01-18) - Risk Model Estimation (Feature Set 1)
- **1.0.0** (2026-01-01) - Foundation and Data Layer

---

## Upcoming Releases

### [1.2.0] - Planned

**Feature Set 2: Portfolio Optimization**

- Mean-variance optimizer
- Quadratic programming solver integration
- Constraint handling (box, budget, turnover)
- Efficient frontier computation
- QUBO formulation for quantum solvers
- Classical QUBO solver (simulated annealing)

### [1.3.0] - Planned

**Feature Set 3: Backtesting Engine**

- Walk-forward backtesting
- Transaction cost modeling
- Rebalancing scheduler
- Portfolio state management
- Trade logging and attribution

### [1.4.0] - Planned

**Feature Set 4: Performance Analytics**

- Return metrics (total, annualized, CAGR)
- Risk metrics (volatility, drawdown, VaR, CVaR)
- Risk-adjusted metrics (Sharpe, Sortino, Calmar)
- Relative performance metrics (alpha, beta, tracking error)
- Attribution analysis

### [1.5.0] - Planned

**Feature Set 5: Visualization and Reporting**

- Equity curve plotting
- Weight evolution visualization
- Risk-return scatter plots
- Drawdown charts
- Rolling metrics
- HTML report generation

---

## Notes

- All releases maintain backward compatibility within major versions
- Configuration format changes are documented in migration guides
- Performance benchmarks are measured on Intel i7 @ 3.5 GHz
- Test coverage targets >90% for all public APIs

[1.1.0]: https://github.com/darkhorse286/portfolio-optimizer/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/darkhorse286/portfolio-optimizer/releases/tag/v1.0.0