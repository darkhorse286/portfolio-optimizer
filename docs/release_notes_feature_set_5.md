# Release Notes: Feature Set 5 — Visualization & Reporting

Release: Feature Set 5 (Initial Delivery)
Date: 2026-04-04

Summary
-------
Feature Set 5 delivers the first iteration of visualization and reporting for the Portfolio Optimizer project. The release focuses on a reliable, reproducible reporting pipeline that integrates with the existing C++ backtest and analytics outputs.

Highlights
----------
- Added `scripts/generate_report.py` — a Python 3.8+ script that:
  - Loads `results/analytics_timeseries.csv` and `results/analytics_metrics.json`
  - Outputs `results/equity_curve.png`, `results/drawdown.png`, and a self-contained `results/report.html`
  - Uses `pandas`, `matplotlib`, `numpy`, and `jinja2`
  - Provides CLI flags: `--input-dir`, `--output-dir`, `--title`
- Implemented plotting functions:
  - `plot_equity_curve()` — normalizes NAV, shades underwater regions, and annotates worst drawdown
  - `plot_drawdown()` — fills drawdown area and marks worst trough
  - Plot implementations follow a consistent styling and save PNGs at 150 DPI
- HTML report builder:
  - `build_html_report()` embeds PNGs as base64 data URIs and formats a metrics table
  - Inline CSS template for a single-file, shareable report
- C++ integration:
  - `--report` flag added to `src/main.cpp` to optionally invoke the Python reporter after a backtest
  - Safe behavior when `python3` is not available (prints a warning, continues execution)
- Docker improvements for reproducible builds:
  - Use a dedicated Python virtualenv inside the image (`/opt/venv`) to install visualization deps
  - Added `.dockerignore` to avoid copying host build artifacts and virtualenvs
  - Build OSQP from source in the image when system package is unavailable
- Testing & CI readiness:
  - Python unit tests and integration tests (pytest)
  - C++ tests (Catch2) and a C++ smoke test for `--report` behavior
  - `run_all_tests.sh` creates a local venv, runs Python tests, C++ `ctest`, and the smoke test

Backward Compatibility
----------------------
- The `--report` flag is optional; existing workflows remain unchanged unless the flag is specified.
- The Python reporter is a separate component and does not affect core C++ analytics outputs.

Known Limitations & Next Steps
------------------------------
- Report plotting currently implements a core subset of charts; additional charts (weight evolution, efficient frontier) will be added in subsequent iterations.
- Interactive HTML/Plotly exports are planned for Feature Set 5.x.
- Add CI integration (GitHub Actions) to run `./run_all_tests.sh` and build the Docker image.

Credits
-------
- Implementation and tests: repository contributors
- Plotting template and Jinja2 integration: internal tooling

