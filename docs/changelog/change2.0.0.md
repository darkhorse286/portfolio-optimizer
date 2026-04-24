# Changelog Entry for v2.0.0

## [2.0.0] - 2026-04-24

### Added

- Quantum optimization layer (`include/quantum/`, `src/quantum/`):
  - `QuantumOptimizer` abstract base, parallel hierarchy to classical `Optimizer`
  - `QUBOFormulation` class: mean-variance → binary quadratic unconstrained form with three-pass matrix construction (risk, return, budget penalty)
  - `SimulatedAnnealingSolver` classical QUBO baseline with O(n) delta_E and multi-start support
  - `BenchmarkRunner` for uniform execution of multiple solvers on identical data with `SolverResult` standardization
  - `QiskitSolver` C++ subprocess wrapper
  - `QuantumSolverFactory` registration pattern

- Qiskit integration (`scripts/quantum/`):
  - `qiskit_submit.py`: circuit construction (QAOA, QAMO, QAMOO modes), backend submission (Aer, IBM Quantum), job queue management
  - `qiskit_collect.py`: async job polling, bitstring decoding, portfolio metric computation via static backtest
  - `qiskit_solver.py`: QAMO algorithm (mean-field mixer with self-consistent equations), QAMOO sweep with Pareto deduplication, COBYLA random restarts
  - `ibm_credentials_check.py`: credential validation, backend enumeration, system property queries
  - `benchmark_viz.py` (extended from Phase 5): frontier comparison chart, 8-color solver palette (Markowitz, SA, QAOA×2, QAMO×2, QAMOO×2)
  - JSON file protocol: `quantum_jobs.json`, `quantum_result_<id>.json` with circuit metadata, execution times, readout correction status

- Frontier comparison visualization:
  - `plot_frontier_comparison()` in benchmark_viz.py overlays QAMOO points on classical `EfficientFrontier` CSV
  - Annualized axes (return × 252, volatility × √252)
  - NISQ overhead annotation and signal quality callout

- IBM Quantum integration:
  - Readout error correction via mthree (M3) with graceful fallback
  - Circuit depth recording (pre- and post-transpilation)
  - Backend calibration date extraction
  - Signal quality flag (low if top bitstring probability < 5%)
  - Execution time recording (microseconds)
  - `--list-backends` flag in qiskit_submit.py

- Configuration:
  - `portfolio_config.json` section `quantum` with `qamo` (mf_iterations, tol) and `qamoo` (n_frontier_points, lambda_min, lambda_max) subsections

- Command-line flags:
  - `--quantum-submit` and `--quantum-collect` in src/main.cpp for async job management

- C++ tests:
  - `test_qubo_formulation.cpp` — matrix construction, symmetry, bit encoding validation
  - `test_simulated_annealing.cpp` — baseline solver correctness, restarts
  - `test_benchmark_runner.cpp` — multi-solver orchestration, result aggregation
  - `test_qiskit_solver_integration.cpp` — subprocess wrapper, job queue protocol

- Python tests:
  - `tests/python/quantum/test_qiskit_solver.py` (19 tests): mean-field solver, QAMO circuit generation, QUBO matrix construction, QAMOO frontier deduplication, frontier chart rendering, HTML report section generation

- Documentation:
  - `docs/release_notes_feature_set_6.md`
  - `docs/changelog/change2.0.0.md`
  - Updated `README.md` with Feature Set 6 completion and v2.0.0 status

### Changed

- `src/backtest/backtest_engine.hpp` / `src/backtest/backtest_engine.cpp`: added `run()` overload accepting `QuantumOptimizer` for quantum solver backtests
- `src/main.cpp`: added `--quantum-submit` and `--quantum-collect` flag handling, quantum submit phase before benchmarking
- `CMakeLists.txt`: added `quantum` library target, test linkage, Python virtualenv configuration for CI
- `scripts/generate_report.py` → `scripts/quantum/benchmark_viz.py`: extended with solver palette, frontier comparison, quantum_report.html generation
- `run_all_tests.sh`: added quantum test suite execution, C++ and Python test integration, TEST_IBM conditional block for hardware runs

### Architecture

- Quantum solvers do not inherit from classical `Optimizer`; instead, both register with `BenchmarkRunner` that unifies execution via `SolverResult` outputs. This separation makes architectural intent explicit and avoids type-conversion overhead.

- Two-script async model (qiskit_submit.py + qiskit_collect.py) allows asynchronous IBM Quantum job management while supporting synchronous Aer reruns, decoupling job submission from result collection via JSON protocol.

- Fair frontier comparison: both classical and quantum solvers use identical annualized expected_returns (simple daily mean × 252) and annualized EWMA covariance (λ=0.94, daily cov × 252), ensuring mathematical comparability.

- QAMOO frontier noisy relative to classical: gap reflects depth-1 NISQ hardware overhead, not algorithm failure. Circuit depth ∝ problem size; deeper circuits require error mitigation (ZNE, M3) for improvement.

### Performance

| Operation | Configuration | Time |
|-----------|---------------|------|
| QAOA circuit (p=1, 20 qubits) | Aer simulator, 1024 shots | ~21s per run |
| QAMO optimization (COBYLA) | 50 iterations, Aer backtest | ~2.5s per lambda |
| QAMOO sweep | 20 lambda points, COBYLA | ~50s total (Aer) |
| Frontier comparison chart | 6 quantum + 1 classical point | <100ms |
| Full quantum_report.html | 6 solvers, 3 embedded PNGs | <0.5s |

All timings include Python startup and Qiskit overhead.

### Code Quality

- C++ tests: 4 test files covering QUBO construction, SA baseline, multi-solver orchestration, Qiskit subprocess interface
- Python tests: 19 pytest test cases covering mean-field solver, circuit generation, matrix construction, frontier deduplication, HTML rendering
- Test runner: `run_all_tests.sh` integrates C++ (Catch2), Python (pytest), and optional IBM hardware tests via TEST_IBM guard
- Numerical validation: QUBO matrix verified symmetric; mean-field equations converge within tolerance; bitstring weights sum to 1.0

### Design Decisions

1. Parallel abstract base `QuantumOptimizer` (not a subclass of `Optimizer`): Quantum and classical solvers use different problem formulations and methods; separate interfaces prevent type confusion.

2. Two-script async model (qiskit_submit.py + qiskit_collect.py): IBM Quantum jobs are asynchronous with unpredictable queue waits; separation allows overnight submission and delayed collection.

3. Python for all circuit code: Qiskit SDK is Python-only; JSON file handoff avoids FFI complexity and reuses existing `--report` design pattern.

4. Fair frontier comparison: annualized expected_returns (daily mean × 252) and EWMA covariance (λ=0.94, daily cov × 252) applied to both classical and quantum, ensuring mathematical comparability.

5. QAMOO frontier gap as a scientific finding: depth-1 circuits on NISQ hardware introduce overhead; not a failure but a measurement of current hardware capability. Future work: deeper circuits, zero-noise extrapolation, Q matrix sparsification.

### Dependencies

- New Python: `qiskit`, `qiskit-aer`, `qiskit-ibm-runtime`, `mthree`, `scipy`
- No new C++ dependencies
- IBM Quantum: free tier (10 min/month QPU) or paid plan; credentials via environment variables

### Breaking Changes

- None. Quantum layer is opt-in; classical code paths unchanged. All new APIs are additive.

### Known Issues

- QAMOO frontier noisy at depth-1: NISQ hardware overhead. Future: Q sparsification, ZNE, deeper circuits.
- IBM run time includes queue wait (minutes to hours). Future: Kafka job queue for parallel submissions.
- Phase 5 deferred items (weight evolution charts, Plotly exports) available for future minor release.

### Test Results

- C++ tests: 4 new test files; all passing
- Python tests: 19 pytest cases; all passing
- Integration test: quantum_report.html generation; all passing
- Full test suite: 176 tests (100% passing)

---

## Previous Releases

- **1.7.0** (2026-04-04) - Visualization and Reporting (Feature Set 5)
- **1.6.0** (2026-02-17) - Performance Analytics (Feature Set 4)
- **1.5.0** (2026-02-16) - Backtesting Engine (Feature Set 3)
- **1.4.0** (2026-02-12) - Advanced Constraints and Optimization (Feature Set 2.3)
- **1.3.0** (2026-02-04) - Advanced Portfolio Optimization (Feature Set 2.2)
- **1.2.0** (2026-01-30) - Mean-Variance Optimization (Feature Set 2.1)
- **1.1.0** (2026-01-19) - Risk Model Estimation (Feature Set 1)
- **1.0.0** (2026-01-17) - Foundation and Data Layer
