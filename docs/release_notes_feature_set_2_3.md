# Release Notes: v1.4.0 - Tracking Error Penalty & Enhancements

**Release Date**: February 12, 2026  
**Feature Set**: 2.3 - Tracking Error (Penalty Method)

## Overview

Feature Set 2.3 introduces a production-grade, soft tracking-error constraint implemented via a penalty term in the objective. The penalty-based formulation enables fair comparison with quantum-ready QUBO encodings (which support soft constraints) while remaining fully compatible with classical QP solvers.

This release also includes validation, auto-tuning for the penalty coefficient, and unit tests covering the new behavior.

## New Features

### 1. Tracking Error as a Soft Penalty

- Formulation: add to objective

  \(\lambda_{te} (w - w_b)^T \Sigma (w - w_b)\)

  Expands to: \(\lambda_{te} (w^T \Sigma w - 2 w^T \Sigma w_b + w_b^T \Sigma w_b)\)

- Implementation details:
  - New fields on `OptimizationConstraints`:
    - `Eigen::VectorXd benchmark_weights` (empty = no TE penalty)
    - `double max_tracking_error` (target TE threshold; must be > 0 when benchmark provided)
    - `double tracking_error_penalty` (penalty coefficient; if <= 0 it is auto-tuned)
  - Auto-tuning heuristic implemented as:

    \[ penalty = 1.0 / (max_tracking_error^2 * trace(\Sigma) / n) \]

    where \(n\) is the number of assets and `trace(Σ)/n` is the average asset variance.

  - Integration points: the penalty term modifies the QP objective (both quadratic and linear terms) in all optimization modes:
    - `optimize_min_variance`
    - `optimize_target_return`
    - `optimize_risk_aversion`
    - `optimize_max_sharpe_direct` (Schaible iteration included)

### 2. Direct Max Sharpe (Schaible refinement)

- Implemented Schaible's fractional programming refinement for direct maximization of the Sharpe ratio as an alternative to grid search. The method iterates by solving a sequence of QP subproblems (see `schaibles_iteration()` in `MeanVarianceOptimizer`). It is exposed via `optimize_max_sharpe_direct` and invoked as a HYBRID/DIRECT method. Unit tests validate convergence on simple problems.

### 3. Sector / Group Constraints

- Added `SectorMapping` helper (include/data/sector_mapper.hpp) to convert asset-to-sector mappings into `GroupConstraint` entries consumable by the optimizer. `GroupConstraint` supports min/max exposure per group and is validated by `GroupConstraint::validate()`.

### 4. Inequality Constraint Support

- The QP specification (`QuadraticProblem`) includes `A_ineq`, `b_ineq_lower`, and `b_ineq_upper` to represent linear inequality constraints. OSQP-based tests (`test_osqp_inequality.cpp`) verify support for simple inequality constraints (e.g., sector exposure caps expressed as A_ineq rows). Full production integration for complex group-based inequalities is available as needed.

  - Objective adjustments follow the QP form used by our solver: for a QP written as

    \( (1/2) x^T P x + q^T x \), adding \(\lambda_{te} w^T \Sigma w\) corresponds to adding \(2 \lambda_{te} \Sigma\) to `P` (implementation uses consistent factors per subproblem) and adding the linear term \(-2 \lambda_{te} \Sigma w_b\) to `q`.

### 2. Validation and Warnings

- `OptimizationConstraints::validate()` now:
  - Throws `std::invalid_argument` if `benchmark_weights` is provided but `max_tracking_error <= 0`.
  - Emits a warning (stderr) when `benchmark_weights` do not sum to 1.0 (tolerance 0.01); this is non-fatal.

### 3. Unit Tests

- Added `tests/test_tracking_error.cpp` covering:
  - Penalty auto-tuning correctness
  - Validation error when `max_tracking_error` missing
  - Warning behavior when benchmark weights do not sum to 1
  - Integration test showing strong penalty pulls min-variance solution toward benchmark

## API Changes

- `OptimizationConstraints` (header: include/optimizer/optimizer_interface.hpp)
  - New public fields: `benchmark_weights`, `max_tracking_error`, `tracking_error_penalty`
  - New method: `double compute_tracking_error_penalty(const Eigen::MatrixXd&) const`

- No changes to public `OptimizerInterface` signatures. All changes are additive and backward compatible.

## Configuration

Add the following optional fields in optimizer constraints JSON (example):

```json
"constraints": {
  "min_weight": 0.0,
  "max_weight": 0.3,
  "long_only": true,
  "sum_to_one": true,
  "benchmark_weights": [0.5, 0.5],
  "max_tracking_error": 0.05,
  "tracking_error_penalty": 0.0
}
```

- If `tracking_error_penalty` is `0.0` or negative, the system will auto-tune using `max_tracking_error` and the covariance matrix.
- If `benchmark_weights` is omitted or empty, no tracking-error penalty is applied.

## Testing

- New tests added in `tests/test_tracking_error.cpp` (integration + unit checks)
- Existing test suites remain passing
- Full test run confirms all tests pass locally

## Performance

- Penalty integration is O(n^2) for dense `Σ` updates (same complexity as the base QP construction)
- Auto-tuning requires only `trace(Σ)` (O(n)) and is negligible relative to solver cost
- No measurable runtime regression for typical problem sizes (N &lt; 200)

## Backwards Compatibility

- All existing code remains functional. New fields are optional.
- Behavior is unchanged unless `benchmark_weights` is set.

## Notes and Rationale

- Using a penalty (soft) formulation ensures fair comparison with quantum approaches that typically express constraints as energy penalties in QUBO form.
- The auto-tuning heuristic provides a reasonable default penalty scale so users do not need to hand-tune `lambda_te` in most cases.

## Contributors

- Feature implementation and tests: Feature Set 2.3 team

## Next Steps

- Expand documentation examples to demonstrate tracking-error tuning
- Consider adding a CLI flag to enable/disable auto-tuning
- Investigate alternate tuning heuristics (e.g., using portfolio-level variance)
