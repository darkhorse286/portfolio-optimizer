# Changelog Entry for v1.7.0

## [1.7.0] - 2026-04-04

### Added

- `scripts/generate_report.py`: self-contained Python 3.8+ report generator that
  - loads `results/analytics_timeseries.csv` and `results/analytics_metrics.json`
  - produces `results/equity_curve.png`, `results/drawdown.png`, and a single-file `results/report.html`
  - CLI flags: `--input-dir`, `--output-dir`, `--title`
- Plotting functions implemented (initial set):
  - `plot_equity_curve()` — NAV normalization, underwater shading, worst-drawdown annotation
  - `plot_drawdown()` — drawdown fill, worst-trough marker, percent formatting
- HTML report builder: `build_html_report()` renders an inline Jinja2 template and embeds PNGs as base64 data URIs
- C++ integration: `--report` flag in `src/main.cpp` to optionally invoke the report generator after backtest completion; safe fallback if `python3` is not available
- Python test suite (pytest): unit tests for plotting/report functions and an integration test that runs the report script end-to-end
- C++ smoke test: verifies the `--report` behavior when `python3` is not present
- Test runner: `run_all_tests.sh` creates `./.venv`, runs pytest suites, runs C++ `ctest`, and performs the smoke test
- Docker improvements: create a Python virtualenv inside the image for visualization deps, add `.dockerignore`, build OSQP from source when system packages are unavailable
- Documentation and release materials:
  - `docs/release_notes_feature_set_5.md`
  - `docs/changelog/CHANGELOG.md`
  - Updated `README.md` to reflect Feature Set 5 deliverables

### Changed

- Feature Set 5 moved from "Planned" to "Initial Delivery" in `README.md` with an updated list of deliverables and notes about future work
- Dockerfile: use a contained Python venv for pip installs, ensure clean build directory during image creation, add BLAS/LAPACK build deps and OSQP source build step
- `docker-compose.yml` guidance: corrected `visualizer` flags and recommended `docker compose run` usage when passing additional command-line args
- Test harness: improved fake writer utilities, ensured integration tests use project venv Python when available

### Architecture

- Reporting pipeline is implemented as a thin Python layer that consumes the existing C++ analytics outputs. This keeps the C++ analytics logic unchanged while enabling richer visualization and reporting in an interpreted environment.
- The HTML report is produced server-side (no external JS/CSS), making it fully portable and shareable as a single file.
- Docker builds now isolate Python dependencies inside `/opt/venv` to avoid host/system packaging constraints and to maintain reproducible image builds.

### Performance

- Typical report generation (252-day daily series) — chart generation + HTML render: < 0.5s on developer workstation (PNG write at 150 DPI, 10x4 & 10x3 figures)
- Equity/drawdown plotting functions are lightweight; for 252-day series they complete in < 50ms each on a modern laptop

### Code Quality

- Python unit tests: coverage for plotting stubs and HTML builder logic (5 unit tests passing)
- Python integration tests: end-to-end script invocation with temporary results files (1 integration test passing)
- C++ tests: existing Catch2 test suite (13 targets) run via `ctest` as part of `run_all_tests.sh`
- Test runner: `run_all_tests.sh` standardizes test execution and creates an isolated Python venv for test dependencies

### Design Decisions

- Reports are intentionally single-file HTML artifacts (base64-embedded PNGs) to simplify sharing and archival without external dependencies.
- Drawdown series follows the non-positive convention (0.0 at peak, negative during drawdown) used across the analytics layer.
- The `--report` flag in `main.cpp` uses `std::system()` to call the Python script; the program checks for `python3` first and emits a clear warning if unavailable rather than failing.

### Dependencies

- New runtime dependencies for reporting (Python): `pandas`, `numpy`, `matplotlib`, `jinja2` (installed into image venv or expected in developer venv)
- Docker build dependencies: `git`, `cmake`, `build-essential`, `libblas-dev`, `liblapack-dev`, and other base build tools used to compile OSQP

### Breaking Changes

- None in the public API. The `--report` flag is optional and default behavior is unchanged.

---

## Previous Releases

- **1.6.0** (2026-02-17) - Performance analytics and DrawdownAnalysis (Feature Set 4)
- **1.5.0** (2026-02-16) - Backtesting Engine (Feature Set 3)
- **1.4.0** (2026-02-12) - Advanced Constraints and Optimization (Feature Set 2.3)
- **1.3.0** (2026-02-04) - Advanced Portfolio Optimization (Feature Set 2.2)
- **1.2.0** (2026-01-30) - Mean-Variance Optimization (Feature Set 2.1)
- **1.1.0** (2026-01-19) - Risk Model Estimation (Feature Set 1)
- **1.0.0** (2026-01-17) - Foundation and Data Layer
