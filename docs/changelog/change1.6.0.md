# Changelog Entry for v1.6.0

## [1.6.0] - 2026-02-17

### Added

- Performance analytics layer (`namespace portfolio::analytics`) with five production components
- `PerformanceMetrics` class: total return, CAGR, rolling returns, monthly return table, annualized volatility, downside deviation, max drawdown with full detail, VaR/CVaR (historical and parametric), skewness, kurtosis, Sharpe ratio, Sortino ratio, Calmar ratio, Omega ratio, with JSON/CSV export and formatted summary output
- `DrawdownAnalysis` class: discrete drawdown event identification (peak-trough-recovery cycles), top-N extraction by depth, aggregate summary statistics (average depth, average duration, time in drawdown percentage), underwater curve, formatted report
- `BenchmarkAnalysis` class: CAPM regression via OLS (alpha, beta, R-squared), tracking error, information ratio, active return, up/down capture ratios, with JSON export
- `RollingStatistics` class: rolling window computations for volatility, mean return, Sharpe ratio, Sortino ratio, max drawdown, beta, tracking error, skewness, kurtosis, plus generic `apply()` for custom rolling functions; configurable via `RollingConfig` struct
- `Attribution` class: Brinson-Fachler single-period performance attribution with allocation, selection, and interaction effects; multi-period linking via Carino logarithmic smoothing method; formatted report output
- `BacktestResult::compute_analytics()` method bridging backtest output to the full analytics suite
- `BacktestResult::export_analytics_json()` and `BacktestResult::export_analytics_csv()` for analytics-powered export
- `backtest_result.hpp` forwarding header for clean analytics-to-backtest dependency
- Beasley-Springer-Moro inverse normal CDF implementation for parametric VaR/CVaR
- Fisher bias-corrected skewness and kurtosis (matching scipy.stats definitions)
- Comprehensive test suite: 5 test files with 70+ test cases covering all public APIs
- Release notes: `docs/release_notes_feature_set_4.md`
- Changelog entry: `docs/changelog/change1.6.0.md`

### Changed

- **Activated optimizer loop** in `BacktestEngine::run()`: removed the `if (false && ...)` debug guard that disabled portfolio rebalancing; the walk-forward backtest now performs live optimization, risk estimation, rebalance scheduling, trade execution, and cost deduction on each eligible date
- Enhanced `BacktestResult::print_summary()` to use the full analytics layer when sufficient data is available, with graceful fallback to inline metrics
- Updated `CMakeLists.txt`: added `analytics` library target linked to `backtest`, added `analytics` to test library list, linked `analytics` to main executable, updated project version to 1.6.0, updated configuration summary
- Updated project status to 100% complete (5/5 feature sets)

### Architecture

- Analytics layer uses `std::vector<double>` throughout (no Eigen dependency) for simpler consumer interface; data volumes (single time series) do not benefit from SIMD
- Forward declaration of `BacktestResult` in `performance_metrics.hpp` avoids circular include dependency
- Lazy drawdown cache via `mutable` members for efficient repeated access
- Pointer+length internal rolling apply avoids per-window vector copies
- Carino linking handles edge cases (zero active return, zero single-period active) via L'Hopital limits

### Performance

All performance targets met:

| Component | Target | Achieved |
|-----------|--------|----------|
| All metrics on 252-day series | < 5ms | < 2ms |
| Rolling Sharpe (252-day window, 1000 days) | < 10ms | < 5ms |
| Full analytics report | < 20ms | < 10ms |
| Drawdown analysis on 2500 days | < 5ms | < 2ms |

### Code Quality

- Production code (analytics): ~3,400 lines (1,238 headers + 2,165 implementations)
- Test code: ~1,944 lines across 5 test files
- Integration code: ~400 lines (updated CMakeLists.txt, backtest_engine.hpp/cpp, backtest_result.hpp, changelog)
- Total new code: ~5,750 lines
- Compiler warnings: 0 (with -Wall -Wextra -Wpedantic)
- Doxygen comments on all public APIs
- Numerical accuracy validated against NumPy/scipy reference implementations

### Design Decisions

- VaR/CVaR returns positive loss magnitude (risk management convention)
- Drawdown series uses non-positive values (0.0 at peak, -0.10 = 10% below)
- Alpha annualized via simple scaling (daily_alpha * 252), not compounding
- Rolling output aligned to window end (pandas convention)
- Downside deviation uses full count in denominator (standard semi-deviation)
- Bias-corrected moments match scipy.stats.skew and scipy.stats.kurtosis

### Dependencies

No new external dependencies. The analytics layer uses only the existing C++17 standard library and nlohmann/json (already in the project).

### Breaking Changes

- `BacktestResult::print_summary()` now produces richer output when analytics are available. The basic metrics are still present; the output is a superset of the previous format.
- Activating the optimizer loop changes backtest behavior: results now reflect actual portfolio optimization and rebalancing. Backtests that previously showed flat equal-weight returns will now show optimized allocation performance.

---

## Previous Releases

- **1.5.0** (2026-02-16) - Backtesting Engine (Feature Set 3)
- **1.4.0** (2026-02-12) - Advanced Constraints and Optimization (Feature Set 2.3)
- **1.3.0** (2026-02-04) - Advanced Portfolio Optimization (Feature Set 2.2)
- **1.2.0** (2026-01-30) - Mean-Variance Optimization (Feature Set 2.1)
- **1.1.0** (2026-01-19) - Risk Model Estimation (Feature Set 1)
- **1.0.0** (2026-01-17) - Foundation and Data Layer