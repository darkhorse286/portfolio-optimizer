#!/usr/bin/env python3
"""QAMO and QAMOO quantum optimization algorithms for portfolio optimization.

QAMO: Quantum Alternating Mean-field Optimization — replaces the QAOA mixer
with per-qubit RX angles derived from mean-field self-consistent equations.

QAMOO: Quantum Alternating Multi-Objective Optimization — sweeps over a range
of risk-aversion (lambda) values using QAMO as the inner solver, producing a
set of (risk, return) Pareto-frontier points for comparison against the
classical Markowitz efficient frontier.
"""

from __future__ import annotations

import numpy as np
from typing import Any, Dict, List, Tuple

from qiskit import QuantumCircuit, transpile


# ---------------------------------------------------------------------------
# Mean-field solver
# ---------------------------------------------------------------------------

def solve_mean_field(
    Q: np.ndarray,
    beta: float,
    max_iterations: int,
    tol: float,
) -> Tuple[np.ndarray, bool, int]:
    """Solve the mean-field self-consistent equations for a QUBO problem.

    Args:
        Q: QUBO matrix (n x n).
        beta: Inverse temperature parameter.
        max_iterations: Maximum self-consistent iterations.
        tol: Convergence tolerance on max |delta_m|.

    Returns:
        m: Magnetisations in [-0.5, 0.5] where m_i = <x_i> - 0.5.
        converged: Whether the equations converged.
        iterations: Number of iterations performed.
    """
    n = Q.shape[0]
    m = np.zeros(n)
    for iteration in range(max_iterations):
        m_new = np.tanh(-beta * Q @ m) / 2.0
        delta = np.max(np.abs(m_new - m))
        m = m_new
        if delta < tol:
            return m, True, iteration + 1
    return m, False, max_iterations


# ---------------------------------------------------------------------------
# QAMO circuit builder
# ---------------------------------------------------------------------------

def build_qamo_circuit(
    Q: np.ndarray,
    gamma: float,
    beta: float,
    max_mf_iters: int,
    mf_tol: float,
) -> Tuple[QuantumCircuit, np.ndarray]:
    """Build a QAMO circuit for the given QUBO matrix.

    Identical to build_qaoa_circuit() except the mixer step uses per-qubit
    RX angles derived from the mean-field solution instead of a uniform angle.

    Args:
        Q: QUBO matrix (n x n).
        gamma: Cost unitary angle.
        beta: Global mean-field inverse temperature.
        max_mf_iters: Maximum mean-field iterations.
        mf_tol: Mean-field convergence tolerance.

    Returns:
        circuit: QuantumCircuit with n qubits and n classical bits.
        theta: Per-qubit RX angles (for diagnostics), shape (n,).
    """
    n = Q.shape[0]
    m, _converged, _iters = solve_mean_field(Q, beta, max_mf_iters, mf_tol)
    # Map magnetisation in [-0.5, 0.5] to a rotation angle
    theta = 2.0 * np.arctan(m + 0.5)

    circuit = QuantumCircuit(n, n)
    circuit.h(range(n))  # uniform superposition

    # Cost unitary — identical structure to QAOA (RZ for diagonal, RZZ for off-diagonal)
    for i in range(n):
        for j in range(i, n):
            val = Q[i, j]
            if abs(val) < 1e-10:
                continue
            angle = 2.0 * gamma * val
            if i == j:
                circuit.rz(angle, i)
            else:
                circuit.rzz(angle, i, j)

    # Mean-field mixer: per-qubit rotation angles
    for i in range(n):
        circuit.rx(theta[i], i)

    circuit.measure(range(n), range(n))
    return circuit, theta


# ---------------------------------------------------------------------------
# Python QUBO matrix builder (mirrors QUBOFormulation::build() in C++)
# ---------------------------------------------------------------------------

def build_qubo_matrix(
    covariance: np.ndarray,
    expected_returns: np.ndarray,
    num_bits: int,
    lambda_: float,
    penalty: float,
) -> np.ndarray:
    """Build the QUBO matrix, matching QUBOFormulation::build() in C++ exactly.

    Three-pass construction:
      Pass 1 — Risk term: 0.5 * cov[i,j] * (2^k * 2^l) / denom^2
      Pass 2 — Return term (negated): -lambda * mu[i] * 2^k / denom  (diagonal)
      Pass 3 — Budget penalty: P*(sum w - 1)^2, expanded into quadratic + linear

    Args:
        covariance: Asset covariance matrix (n x n).
        expected_returns: Expected returns per asset (n,).
        num_bits: Number of bits per asset (B).
        lambda_: Risk-aversion / return scaling parameter (plays role of risk_aversion_).
        penalty: Budget penalty coefficient.

    Returns:
        Q: QUBO matrix (n*B x n*B), symmetrized.
    """
    n = len(expected_returns)
    B = num_bits
    Nvar = n * B
    Q = np.zeros((Nvar, Nvar))
    denom = float((1 << B) - 1)

    # Pre-compute bit weights
    pow2 = np.array([float(1 << k) for k in range(B)])

    # Pass 1: Risk term
    bit_outer = np.outer(pow2, pow2) / (denom * denom)  # (B, B)
    for i in range(n):
        for j in range(n):
            pi = slice(i * B, (i + 1) * B)
            pj = slice(j * B, (j + 1) * B)
            Q[pi, pj] += 0.5 * covariance[i, j] * bit_outer

    # Pass 2: Return term (diagonal only, negated)
    for i in range(n):
        for k in range(B):
            p = i * B + k
            Q[p, p] -= lambda_ * expected_returns[i] * (pow2[k] / denom)

    # Pass 3: Budget penalty P*(sum_p c_p x_p - 1)^2
    c = np.array([pow2[k] / denom for i in range(n) for k in range(B)])
    Q += penalty * np.outer(c, c)                         # quadratic part
    Q[np.arange(Nvar), np.arange(Nvar)] -= 2.0 * penalty * c  # linear part

    # Symmetrize (matches C++ final step)
    Q = 0.5 * (Q + Q.T)
    return Q


# ---------------------------------------------------------------------------
# Helpers shared with QAMOO
# ---------------------------------------------------------------------------

def _decode_bitstring_np(bitstring: str, num_assets: int, num_bits: int) -> np.ndarray:
    """Decode a Qiskit bitstring into portfolio weights as a numpy array.

    Args:
        bitstring: Measurement bitstring in Qiskit convention (reversed LSB).
        num_assets: Number of assets.
        num_bits: Bits per asset.

    Returns:
        weights: Portfolio weights in [0, 1], shape (num_assets,).
    """
    bits = [int(c) for c in reversed(bitstring.strip())]
    denom = float((1 << num_bits) - 1)
    weights = np.zeros(num_assets)
    for i in range(num_assets):
        start = i * num_bits
        value = sum(bits[start + k] * (1 << k) for k in range(num_bits))
        weights[i] = value / denom
    return weights


def _objective_np(Q: np.ndarray, bitstring: str) -> float:
    """Compute QUBO objective x^T Q x for a bitstring (Qiskit convention)."""
    bits = np.array([int(c) for c in reversed(bitstring.strip())], dtype=float)
    return float(bits @ Q @ bits)


def _execute_aer(circuit: QuantumCircuit, backend_instance: Any, shots: int) -> Dict[str, int]:
    """Transpile and run a circuit on an Aer backend, returning measurement counts."""
    transpiled = transpile(circuit, backend_instance)
    result = backend_instance.run(transpiled, shots=shots).result()
    return result.get_counts()


# ---------------------------------------------------------------------------
# COBYLA objective for QAMO
# ---------------------------------------------------------------------------

def evaluate_qamo(
    Q: np.ndarray,
    gamma: float,
    beta: float,
    config: Dict[str, Any],
    backend_instance: Any,
    shots: int,
) -> float:
    """Evaluate the mean QUBO energy for a QAMO circuit parameterised by (gamma, beta).

    Used as the objective function for COBYLA parameter optimisation.

    Args:
        Q: QUBO matrix.
        gamma: Cost unitary angle.
        beta: Mean-field inverse temperature.
        config: Must contain 'max_mf_iterations' and 'mf_convergence_tol'.
        backend_instance: Aer simulator backend.
        shots: Number of measurement shots per evaluation.

    Returns:
        Mean QUBO energy E[x^T Q x] over the measurement distribution.
    """
    circuit, _ = build_qamo_circuit(
        Q, gamma, beta,
        config["max_mf_iterations"],
        config["mf_convergence_tol"],
    )
    counts = _execute_aer(circuit, backend_instance, shots)
    total = sum(counts.values())
    if total == 0:
        return 0.0
    return sum(
        (cnt / total) * _objective_np(Q, bs)
        for bs, cnt in counts.items()
    )


# ---------------------------------------------------------------------------
# COBYLA optimiser with random restarts
# ---------------------------------------------------------------------------

def _optimize_qamo_params(
    Q: np.ndarray,
    config: Dict[str, Any],
    backend_instance: Any,
    shots: int,
    n_restarts: int = 4,
    rng: "np.random.Generator | None" = None,
) -> Tuple[float, float, Dict[str, int]]:
    """Run COBYLA for QAMO with random restarts when the result is all-zeros.

    Starts from [0.1, 0.1]. If the most-probable bitstring is the trivial
    all-zeros solution, retries up to n_restarts times from random starting
    points drawn uniformly from (0.05, pi) x (0.05, pi), stopping as soon as
    a non-trivial bitstring is found.

    Args:
        Q: QUBO matrix (n x n).
        config: Must contain 'max_mf_iterations' and 'mf_convergence_tol'.
        backend_instance: Aer backend.
        shots: Shots per COBYLA evaluation and per final circuit run.
        n_restarts: Maximum number of random-restart attempts.
        rng: Optional seeded numpy Generator for reproducibility.

    Returns:
        gamma_opt: Best gamma found.
        beta_opt: Best beta found.
        counts: Measurement counts from the final circuit run.
    """
    from scipy.optimize import minimize as sp_minimize

    if rng is None:
        rng = np.random.default_rng()

    zero_bs = "0" * Q.shape[0]

    def _attempt(x0: List[float]) -> Tuple[float, float, Dict[str, int]]:
        opt = sp_minimize(
            lambda p: evaluate_qamo(Q, p[0], p[1], config, backend_instance, shots),
            x0=x0,
            method="COBYLA",
            options={"maxiter": 50, "rhobeg": 0.1},
        )
        g, b = float(opt.x[0]), float(opt.x[1])
        circuit, _ = build_qamo_circuit(
            Q, g, b,
            config["max_mf_iterations"],
            config["mf_convergence_tol"],
        )
        return g, b, _execute_aer(circuit, backend_instance, shots)

    gamma_opt, beta_opt, counts = _attempt([0.1, 0.1])
    best_bs = max(counts, key=counts.get) if counts else zero_bs

    for _ in range(n_restarts):
        if best_bs != zero_bs:
            break
        x0 = rng.uniform(0.05, np.pi, size=2).tolist()
        g, b, c = _attempt(x0)
        if c:
            bs = max(c, key=c.get)
            if bs != zero_bs:
                gamma_opt, beta_opt, counts, best_bs = g, b, c, bs

    return gamma_opt, beta_opt, counts


# ---------------------------------------------------------------------------
# Pareto deduplication
# ---------------------------------------------------------------------------

def deduplicate_frontier(
    points: List[Dict[str, Any]],
    tol: float,
) -> List[Dict[str, Any]]:
    """Remove dominated and near-duplicate Pareto frontier points.

    A point p is dominated if another point q satisfies:
      q.return >= p.return AND q.volatility <= p.volatility
      (with at least one strict inequality).

    Near-duplicates within ``tol`` in both dimensions keep the higher-return point.

    Args:
        points: List of frontier point dicts with 'portfolio_return' and
            'portfolio_volatility' keys.
        tol: Absolute tolerance for near-duplicate detection (e.g. 0.005 = 0.5%).

    Returns:
        Deduplicated and filtered list (may be shorter than input).
    """
    if not points:
        return points

    # Step 1: Remove dominated points
    non_dominated: List[Dict[str, Any]] = []
    for p in points:
        dominated = any(
            q["portfolio_return"] >= p["portfolio_return"]
            and q["portfolio_volatility"] <= p["portfolio_volatility"]
            and (
                q["portfolio_return"] > p["portfolio_return"]
                or q["portfolio_volatility"] < p["portfolio_volatility"]
            )
            for q in points
            if q is not p
        )
        if not dominated:
            non_dominated.append(p)

    # Step 2: Merge near-duplicates (keep higher return within each cluster)
    result: List[Dict[str, Any]] = []
    for p in non_dominated:
        merged = False
        for idx, q in enumerate(result):
            if (
                abs(p["portfolio_return"] - q["portfolio_return"]) <= tol
                and abs(p["portfolio_volatility"] - q["portfolio_volatility"]) <= tol
            ):
                if p["portfolio_return"] > q["portfolio_return"]:
                    result[idx] = p
                merged = True
                break
        if not merged:
            result.append(p)

    return result


# ---------------------------------------------------------------------------
# QAMOO sweep
# ---------------------------------------------------------------------------

def run_qamoo_sweep(
    problem: Dict[str, Any],
    config: Dict[str, Any],
    backend_instance: Any,
    shots: int,
) -> Tuple[List[Dict[str, Any]], np.ndarray]:
    """Run a QAMOO sweep across lambda values to produce a quantum Pareto frontier.

    For each lambda value, rebuilds the QUBO matrix and runs QAMO with COBYLA
    parameter optimisation. Collects (risk, return) portfolio solutions, then
    deduplicates and sorts by volatility.

    Args:
        problem: Loaded from quantum_problem.json; must include 'covariance',
            'expected_returns', 'penalty_budget', 'num_bits_per_asset',
            'num_assets'.
        config: Dict with keys: n_frontier_points, lambda_min, lambda_max,
            max_mf_iterations, mf_convergence_tol.
        backend_instance: Aer backend (COBYLA optimisation runs synchronously).
        shots: Measurement shots per COBYLA evaluation.

    Returns:
        frontier_points: Sorted (ascending volatility) list of Pareto points, each
            containing lambda, portfolio_volatility, portfolio_return, weights,
            bitstring, objective_value.
        best_weights: Weights of the highest Sharpe-ratio solution from the sweep.
    """
    rng = np.random.default_rng(42)

    lambdas = np.logspace(
        np.log10(float(config["lambda_min"])),
        np.log10(float(config["lambda_max"])),
        int(config["n_frontier_points"]),
    )

    covariance = np.array(problem["covariance"])
    expected_returns = np.array(problem["expected_returns"])
    penalty_budget = float(problem["penalty_budget"])
    num_bits = int(problem["num_bits_per_asset"])
    num_assets = int(problem["num_assets"])

    frontier_points: List[Dict[str, Any]] = []

    for lam in lambdas:
        Q = build_qubo_matrix(
            covariance, expected_returns, num_bits, float(lam), penalty_budget
        )

        gamma_opt, beta_opt, counts = _optimize_qamo_params(
            Q, config, backend_instance, shots, rng=rng
        )

        if not counts:
            continue

        # Select most-probable bitstring (spec: max(counts, key=counts.get))
        best_bs = max(counts, key=counts.get)
        weights = _decode_bitstring_np(best_bs, num_assets, num_bits)

        port_return = float(expected_returns @ weights)
        cov_w = covariance @ weights
        port_vol = float(np.sqrt(max(weights @ cov_w, 0.0)))
        obj_val = _objective_np(Q, best_bs)

        frontier_points.append({
            "lambda": float(lam),
            "portfolio_volatility": port_vol,
            "portfolio_return": port_return,
            "weights": weights.tolist(),
            "bitstring": best_bs,
            "objective_value": obj_val,
        })

    frontier_points = deduplicate_frontier(frontier_points, tol=0.005)
    frontier_points.sort(key=lambda p: p["portfolio_volatility"])

    if not frontier_points:
        return frontier_points, np.ones(num_assets) / num_assets

    best = max(
        frontier_points,
        key=lambda p: p["portfolio_return"] / max(p["portfolio_volatility"], 1e-10),
    )
    return frontier_points, np.array(best["weights"])
