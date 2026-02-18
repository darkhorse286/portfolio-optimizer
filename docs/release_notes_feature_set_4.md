# Release Notes: v1.6.0 - Performance Analytics

**Release Date**: February 17, 2026
**Feature Set**: 4 - Performance Analytics
**Previous Version**: 1.5.0 (Feature Set 3 - Backtesting Engine)

## Overview

This release implements Feature Set 4: the performance analytics layer. The system now provides comprehensive return, risk, risk-adjusted, benchmark-relative, rolling, drawdown, and attribution analysis. All five planned components have been delivered, tested against NumPy/scipy reference implementations, and integrated with the existing backtest layer via `BacktestResult::compute_analytics()`.

This release also activates the optimizer loop inside `BacktestEngine::run()` that was disabled behind a debug guard in v1.5.0. Walk-forward backtests now perform live optimization, risk estimation, rebalance scheduling, trade execution, and cost deduction on each eligible date.

## New Components

### 1. PerformanceMetrics (`performance_metrics.hpp` / `performance_metrics.cpp`)

The central analytics class. Accepts either a `BacktestResult` or raw NAV/return/date vectors and computes the full metric suite.

**Return metrics**:

- `total_return()`: Cumulative return over the full period.
- `annualized_return()`: CAGR computed as `(1 + total)^(252/n) - 1`.
- `rolling_returns(window)`: Rolling window cumulative returns.
- `monthly_return_table()`: Returns organized by year and month with annual totals.

**Risk metrics**:

- `annualized_volatility()`: Sample standard deviation of daily returns, annualized by `sqrt(252)`.
- `downside_deviation(target)`: Semi-deviation below a target return, annualized.
- `max_drawdown()` and `max_drawdown_info()`: Peak-to-trough decline with depth, duration, and recovery detail.
- `drawdown_series()`: Full underwater curve (non-positive values, 0.0 at peaks).
- `value_at_risk(confidence, method)`: Daily VaR via historical empirical quantile or parametric normal assumption. Returns a positive loss magnitude.
- `conditional_var(confidence, method)`: Expected Shortfall beyond VaR. Positive loss magnitude.
- `skewness()` and `kurtosis()`: Fisher bias-corrected moments matching `scipy.stats.skew` and `scipy.stats.kurtosis`.

**Risk-adjusted metrics**:

- `sharpe_ratio()`: `(annualized_return - risk_free_rate) / annualized_volatility`. Returns 0 if volatility is zero.
- `sortino_ratio(target)`: Uses downside deviation in the denominator.
- `calmar_ratio()`: `annualized_return / max_drawdown`. Returns 0 if drawdown is zero.
- `omega_ratio(threshold)`: Ratio of gains above threshold to losses below.

**Export**: `summary()` for formatted text, `to_json()` for structured output, `to_csv(path)` for time series data.

**Implementation notes**: The drawdown series is lazy-cached via `mutable` members for efficient repeated access. Parametric VaR/CVaR uses the Beasley-Springer-Moro rational approximation for the inverse normal CDF, achieving accuracy within `1e-9` of `scipy.stats.norm.ppf`.

### 2. DrawdownAnalysis (`drawdown_analysis.hpp` / `drawdown_analysis.cpp`)

Detailed drawdown decomposition beyond the summary metrics in `PerformanceMetrics`.

**Event identification**: A state machine walks the NAV series tracking peak-trough-recovery cycles. Each `DrawdownEvent` records: peak index and NAV, trough index and NAV, recovery index (or -1 if unrecovered), depth as a positive fraction, duration in trading days, and recovery duration.

**Methods**:

- `all_events()`: Chronologically ordered vector of all drawdown events.
- `top_drawdowns(n)`: Top N events sorted by depth (deepest first).
- `worst_drawdown()`: Single worst event. Throws `std::runtime_error` on a monotonically increasing series.
- `event_count()`: Number of discrete drawdown events.
- `summary()`: Returns a `DrawdownSummary` struct with average depth, average duration, longest drawdown, longest recovery, maximum depth, and time in drawdown percentage.
- `underwater_curve()`: The same non-positive drawdown series as `PerformanceMetrics::drawdown_series()`, provided here for convenience.
- `report(max_events)`: Formatted multi-line string listing events with depth, duration, peak/trough dates, and recovery dates.

### 3. BenchmarkAnalysis (`benchmark_analysis.hpp` / `benchmark_analysis.cpp`)

Relative performance analysis against a benchmark series.

**CAPM regression**: Runs OLS on construction. The regression model is `R_p - R_f = alpha + beta * (R_b - R_f) + epsilon`, where `R_f` is the daily risk-free rate. Computes alpha (annualized by `alpha_daily * 252`), beta, R-squared, and residual standard deviation.

**Tracking metrics**: `tracking_error()` is the annualized standard deviation of excess returns (`portfolio - benchmark`). `information_ratio()` is `active_return / tracking_error`. `active_return()` is the annualized mean of excess returns.

**Capture ratios**: `up_capture_ratio()` and `down_capture_ratio()` are the ratio of the portfolio mean return to the benchmark mean return in periods where the benchmark is positive (up) or negative (down), respectively.

**Export**: `summary()` for formatted text, `to_json()` for structured output, `excess_returns()` for the raw excess return series.

### 4. RollingStatistics (`rolling_statistics.hpp` / `rolling_statistics.cpp`)

Rolling window computations over a return series.

**Configuration**: Via `RollingConfig` struct with `window_days` (default 63), `trading_days_per_year` (default 252), `risk_free_rate` (default 0.02), and `min_periods` (default same as window).

**Built-in metrics**: `volatility()`, `mean_return()`, `sharpe_ratio()`, `sortino_ratio(target)`, `max_drawdown()`, `beta(benchmark)`, `tracking_error(benchmark)`, `skewness()`, `kurtosis()`. Each returns a `std::vector<double>` of size `n - window + 1`.

**Generic interface**: `apply(std::function<double(const std::vector<double>&)>)` allows computing any custom statistic over the rolling window.

**Output alignment**: `result[i]` is the metric for the window ending at index `i + window - 1`. The `aligned_dates(dates)` helper extracts the corresponding date labels.

**Implementation notes**: Internal rolling computations use a `rolling_apply_internal` function that passes raw `const double*` pointers and lengths to avoid allocating sub-vectors for each window position. The benchmark-relative methods (`beta`, `tracking_error`) iterate paired series directly.

### 5. Attribution (`attribution.hpp` / `attribution.cpp`)

Brinson-Fachler performance attribution for sector-level decomposition.

**Single-period**: The static method `Attribution::single_period(portfolio_sectors, benchmark_sectors)` computes per-sector allocation, selection, and interaction effects:

- Allocation = `(w_p - w_b) * (R_b_sector - R_b_total)`
- Selection = `w_b * (R_p_sector - R_b_sector)`
- Interaction = `(w_p - w_b) * (R_p_sector - R_b_sector)`

Returns a `SinglePeriodAttribution` struct with per-sector effects and totals. Sector effects sum to total active return by construction.

**Multi-period**: Accumulate periods with `add_period()`, then call `linked_attribution()` to get compounded results via Carino logarithmic smoothing. The Carino factor ensures that the sum of period-level effects equals the total active return computed from compounded portfolio and benchmark returns. Edge cases (near-zero active returns, near-zero single-period active) are handled via L'Hopital limits.

**Export**: `report()` generates a formatted multi-line string with per-period decomposition and cumulative totals.

### 6. BacktestResult Integration

Three new methods on `BacktestResult`:

- `compute_analytics(rf_rate, trading_days)`: Constructs and returns a `PerformanceMetrics` instance from the result's NAV, return, and date series. This is the bridge between backtesting and analytics.
- `export_analytics_json(rf_rate, trading_days)`: Returns a JSON string with all computed metrics.
- `export_analytics_csv(filepath, rf_rate, trading_days)`: Writes a date/NAV/return/drawdown CSV via the analytics layer.

`print_summary()` has been enhanced to automatically use the analytics layer when sufficient data is available (>= 2 NAV observations), with graceful fallback to inline metrics on failure.

## Backtest Engine Changes

The optimizer loop inside `BacktestEngine::run()` has been activated. The `if (false && date_idx >= params_.min_history)` guard has been changed to `if (date_idx >= params_.min_history)`. This means walk-forward backtests now perform the full optimization cycle on each eligible date: extract rolling return window, estimate covariance via the configured risk model, compute expected returns, run the optimizer, evaluate rebalance triggers, execute trades with sell-first ordering, deduct transaction costs, and log trades.

This resolves the known issue documented in the v1.5.0 release notes.

## API Summary

```cpp
namespace portfolio::analytics {

    // Core metrics
    class PerformanceMetrics {
    public:
        explicit PerformanceMetrics(const backtest::BacktestResult& result,
                                    double rf = 0.02, int tdy = 252);
        PerformanceMetrics(const std::vector<double>& nav,
                           const std::vector<double>& returns,
                           const std::vector<std::string>& dates,
                           double rf = 0.02, int tdy = 252);
        double total_return() const;
        double annualized_return() const;
        std::vector<double> rolling_returns(int window) const;
        MonthlyReturnTable monthly_return_table() const;
        double annualized_volatility() const;
        double downside_deviation(double target = 0.0) const;
        double max_drawdown() const;
        DrawdownInfo max_drawdown_info() const;
        std::vector<double> drawdown_series() const;
        double value_at_risk(double conf = 0.95, VaRMethod m = HISTORICAL) const;
        double conditional_var(double conf = 0.95, VaRMethod m = HISTORICAL) const;
        double skewness() const;
        double kurtosis() const;
        double sharpe_ratio() const;
        double sortino_ratio(double target = 0.0) const;
        double calmar_ratio() const;
        double omega_ratio(double threshold = 0.0) const;
        std::string summary() const;
        std::string to_json() const;
        void to_csv(const std::string& filepath) const;
    };

    // Drawdown decomposition
    class DrawdownAnalysis {
    public:
        DrawdownAnalysis(const std::vector<double>& nav,
                         const std::vector<std::string>& dates);
        const std::vector<DrawdownEvent>& all_events() const;
        std::vector<DrawdownEvent> top_drawdowns(int n) const;
        const DrawdownEvent& worst_drawdown() const;
        int event_count() const;
        DrawdownSummary summary() const;
        const std::vector<double>& underwater_curve() const;
        std::string report(int max_events = -1) const;
    };

    // Benchmark-relative analysis
    class BenchmarkAnalysis {
    public:
        BenchmarkAnalysis(const std::vector<double>& portfolio_returns,
                          const std::vector<double>& benchmark_returns,
                          double rf = 0.02, int tdy = 252);
        double alpha() const;
        double beta() const;
        double r_squared() const;
        const RegressionResult& regression() const;
        double tracking_error() const;
        double information_ratio() const;
        double active_return() const;
        double up_capture_ratio() const;
        double down_capture_ratio() const;
        const std::vector<double>& excess_returns() const;
        std::string summary() const;
        std::string to_json() const;
    };

    // Rolling window computations
    class RollingStatistics {
    public:
        RollingStatistics(const std::vector<double>& returns,
                          const RollingConfig& config);
        explicit RollingStatistics(const std::vector<double>& returns,
                                   int window = 63);
        std::vector<double> volatility() const;
        std::vector<double> mean_return() const;
        std::vector<double> sharpe_ratio() const;
        std::vector<double> sortino_ratio(double target = 0.0) const;
        std::vector<double> max_drawdown() const;
        std::vector<double> beta(const std::vector<double>& benchmark) const;
        std::vector<double> tracking_error(const std::vector<double>& benchmark) const;
        std::vector<double> skewness() const;
        std::vector<double> kurtosis() const;
        std::vector<double> apply(
            const std::function<double(const std::vector<double>&)>& func) const;
        int output_size() const;
        std::vector<std::string> aligned_dates(
            const std::vector<std::string>& dates) const;
    };

    // Brinson-Fachler attribution
    class Attribution {
    public:
        static SinglePeriodAttribution single_period(
            const std::vector<SectorAllocation>& portfolio,
            const std::vector<SectorAllocation>& benchmark);
        void add_period(const std::vector<SectorAllocation>& portfolio,
                        const std::vector<SectorAllocation>& benchmark);
        MultiPeriodAttribution linked_attribution() const;
        int period_count() const;
        void clear();
        std::string report() const;
    };
}

// Backtest integration
namespace portfolio::backtest {
    struct BacktestResult {
        // ... existing members ...
        analytics::PerformanceMetrics compute_analytics(
            double rf = 0.02, int tdy = 252) const;
        std::string export_analytics_json(
            double rf = 0.02, int tdy = 252) const;
        void export_analytics_csv(const std::string& filepath,
            double rf = 0.02, int tdy = 252) const;
    };
}
```

## Configuration

The `performance` section of `portfolio_config.json` is now consumed by the analytics layer. All fields are optional:

```json
{
  "performance": {
    "metrics": [
      "total_return", "annualized_return", "volatility",
      "sharpe_ratio", "sortino_ratio", "calmar_ratio",
      "max_drawdown", "alpha", "beta", "information_ratio",
      "tracking_error", "var_95", "cvar_95"
    ],
    "benchmark_comparison": true,
    "risk_free_rate": 0.02,
    "confidence_level": 0.95
  }
}
```

## Dependencies

No new external dependencies. The analytics layer uses only the existing C++17 standard library and nlohmann/json for JSON export.

## Breaking Changes

- `BacktestResult::print_summary()` now produces richer output when the analytics layer is available. The basic metrics are still present; the output is a superset of the previous format.
- Activating the optimizer loop changes backtest behavior: results now reflect actual portfolio optimization and rebalancing. Backtests that previously showed flat equal-weight returns will now show optimized allocation performance.

## Known Issues

- Calendar-based rebalance triggers for weekly and quarterly schedules depend on the first date in the data aligning with the expected calendar pattern. Data starting mid-week or mid-quarter may miss the first trigger. (Inherited from v1.5.0; not changed in this release.)
- The `PerformanceMetrics` constructor from `BacktestResult` expects `return_series.size() == nav_series.size() - 1`. If the backtest produces equal-length NAV and return series (as the current engine does, where both start from date index 1), the constructor accepts this as a valid alternate convention.

## Performance

| Operation | Configuration | Time |
|-----------|---------------|------|
| All PerformanceMetrics | 252-day series | <2ms |
| Drawdown analysis | 2500 days | <2ms |
| BenchmarkAnalysis (full CAPM) | 252-day series | <1ms |
| Rolling Sharpe (252-day window) | 1000-day series | <5ms |
| Rolling beta (63-day window) | 500-day series | <3ms |
| Single-period attribution | 10 sectors | <0.1ms |
| Multi-period Carino linking | 12 periods, 10 sectors | <0.5ms |
| Full analytics report | 252-day series | <10ms |

All performance targets from the project specification are met.

## Testing

```bash
# Build and run all tests
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)
make run_tests

# Run analytics test suites individually
./bin/test_test_performance_metrics    # ~30 test cases
./bin/test_test_drawdown_analysis      # ~15 test cases
./bin/test_test_benchmark_analysis     # ~15 test cases
./bin/test_test_rolling_statistics     # ~20 test cases
./bin/test_test_attribution            # ~15 test cases
```

**Test categories**:

- Constructor validation: rejects empty data, mismatched sizes, insufficient observations
- Happy path: validates each metric against NumPy/scipy reference values
- Edge cases: zero volatility, constant returns, single element, monotonically increasing NAV, identical portfolio and benchmark
- Numerical accuracy: all metrics within `1e-10` of reference implementations
- Error handling: correct exception types with descriptive messages
- Export: summary and JSON output contain expected fields

**Numerical references**: All test reference values were computed using NumPy 1.24, pandas 2.0, and scipy 1.11. The Beasley-Springer-Moro inverse normal CDF was validated against `scipy.stats.norm.ppf` to within `1e-9`. Fisher bias-corrected skewness and kurtosis match `scipy.stats.skew` and `scipy.stats.kurtosis` with `excess=True`.

## Files Added

**Headers (include/analytics/)**:
- `performance_metrics.hpp`
- `drawdown_analysis.hpp`
- `benchmark_analysis.hpp`
- `rolling_statistics.hpp`
- `attribution.hpp`

**Implementation (src/analytics/)**:
- `performance_metrics.cpp`
- `drawdown_analysis.cpp`
- `benchmark_analysis.cpp`
- `rolling_statistics.cpp`
- `attribution.cpp`

**Tests (tests/)**:
- `test_performance_metrics.cpp`
- `test_drawdown_analysis.cpp`
- `test_benchmark_analysis.cpp`
- `test_rolling_statistics.cpp`
- `test_attribution.cpp`

**Integration (modified)**:
- `include/backtest/backtest_engine.hpp` (added `compute_analytics()`, `export_analytics_json()`, `export_analytics_csv()`)
- `include/backtest/backtest_result.hpp` (new forwarding header)
- `src/backtest/backtest_engine.cpp` (activated optimizer loop, added analytics methods, enhanced `print_summary()`)
- `CMakeLists.txt` (added analytics library target, test linkage)

**Documentation**:
- `docs/release_notes_feature_set_4.md` (this file)
- `docs/changelog/change1.6.0.md`
- `README.md` (updated)

## Code Statistics

```
Production code (analytics):    ~3,400 lines (5 headers + 5 implementations)
Test code (analytics):          ~1,944 lines (5 test files)
Integration code:               ~918 lines (CMakeLists.txt, backtest updates, docs)
Total new/modified:             ~6,265 lines

Cumulative project totals:
  Production code:              ~13,900 lines C++17
  Test code:                    ~5,750 lines
  Build/config:                 ~350 lines
  Documentation:                ~9,500 lines
  Total:                        ~29,500 lines
  Tests:                        157 (100% passing)
```

## Migration Guide

No migration required for existing code. To use the new analytics layer:

1. Include the appropriate headers from `analytics/`.
2. Link against the `analytics` library target in CMake (already done in the updated `CMakeLists.txt`).
3. For backtest integration, call `result.compute_analytics()` on any `BacktestResult` to get a `PerformanceMetrics` object.

Note that activating the optimizer loop changes backtest results. If you need the previous behavior (no optimization, equal-weight drift), you can set `params.min_history` to a value larger than the number of dates in your data.

## Next Steps

**Feature Set 5: Visualization and Reporting** will consume the analytics and backtest outputs to generate equity curves, weight evolution charts, drawdown plots, rolling metric visualizations, and HTML reports. The `PerformanceMetrics::to_json()` and `BacktestResult::export_analytics_csv()` methods provide the data interface.

**Feature Set 6: Quantum Optimization** will integrate quantum and quantum-inspired solvers (QUBO, QAOA, QAMO, QAMOO) as alternative optimizers within the backtest framework, enabling direct comparison of classical and quantum approaches on identical historical data using the analytics layer for evaluation.

## References

1. Markowitz, H. (1952). "Portfolio Selection", The Journal of Finance
2. Brinson, G., Hood, L. & Beebower, G. (1986). "Determinants of Portfolio Performance", Financial Analysts Journal
3. Carino, D. (1999). "Combining Attribution Effects Over Time", Journal of Performance Measurement
4. Beasley, J., Springer, S. & Moro, B. (1996). "Algorithm for the Percentage Points of the Normal Distribution"
5. Sortino, F. & van der Meer, R. (1991). "Downside Risk", Journal of Portfolio Management
6. Keating, C. & Shadwick, W. (2002). "A Universal Performance Measure", Journal of Performance Measurement

---

**Release Date**: February 17, 2026
**Project Status**: Foundation + Risk + Optimization + Backtesting + Analytics Complete (85% overall)
**Next Milestone**: Feature Set 5 - Visualization and Reporting