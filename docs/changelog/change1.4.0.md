# Changelog Entry for v1.4.0

## [1.4.0] - 2026-02-12

### Added

- **Tracking Error (soft penalty)**
  - Add tracking-error penalty term to objective: penalizes deviation from a benchmark via \(\lambda_{te} (w - w_b)^T \Sigma (w - w_b)\).
  - New `OptimizationConstraints` fields:
    - `benchmark_weights` (Eigen::VectorXd)
    - `max_tracking_error` (double)
    - `tracking_error_penalty` (double)
  - Auto-tuning method: `OptimizationConstraints::compute_tracking_error_penalty(const Eigen::MatrixXd&)`.
  - Integrated into `optimize_min_variance`, `optimize_target_return`, `optimize_risk_aversion`, and `optimize_max_sharpe_direct` (Schaible iteration).

  - **Direct Max Sharpe (Schaible refinement)**: Fractional programming iteration implemented in `MeanVarianceOptimizer::schaibles_iteration()` and exposed via `optimize_max_sharpe_direct`.

  - **Sector / Group Constraints**: `SectorMapping` utility to produce `GroupConstraint` objects from asset-to-sector mappings; `GroupConstraint` validation and integration points documented.

  - **Inequality Constraints**: `QuadraticProblem` now documents `A_ineq`, `b_ineq_lower`, and `b_ineq_upper` fields and unit tests verify basic inequality handling.

- **Unit tests**
  - `tests/test_tracking_error.cpp` covering heuristic, validation, warning, and integration behaviors.

### Changed

- `OptimizationConstraints::validate()` now enforces positive `max_tracking_error` when `benchmark_weights` provided and warns if benchmark weights do not sum to 1.0 (tolerance 0.01).
- README updated to document Feature Set 2.3 and new configuration fields.

### Fixed

- N/A (minor validation fixes included)

### Testing

- Added tests for tracking-error behavior; all tests pass locally.

### Notes

- Penalty-based TE enables fair comparison with quantum penalty encodings (QUBO) and keeps classical QP pipeline compatible.
- Auto-tuning heuristic: `penalty = 1.0 / (max_tracking_error^2 * trace(Σ) / n)`.

---

Previous: v1.3.0 - Feature Set 2.2 (Efficient Frontier, Turnover Constraints)
