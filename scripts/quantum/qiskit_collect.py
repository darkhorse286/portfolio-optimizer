#!/usr/bin/env python3
"""Collect QAOA job results for portfolio optimization using Qiskit."""

from __future__ import annotations

import argparse
import json
import sys
import time
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
from qiskit import QuantumCircuit, transpile


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Collect QAOA job results from Qiskit backends."
    )
    parser.add_argument(
        "--jobs-file",
        default="results/quantum_jobs.json",
        help="Path to quantum jobs JSON file.",
    )
    parser.add_argument(
        "--results-dir",
        default="results/",
        help="Directory to write result JSON files.",
    )
    parser.add_argument(
        "--timeout-min",
        type=int,
        default=30,
        help="Timeout in minutes for IBM job collection.",
    )
    parser.add_argument(
        "--poll-interval",
        type=int,
        default=30,
        help="Seconds between status checks.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print job status for each pending job without waiting or writing results.",
    )
    return parser.parse_args()


def load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        print(f"Job file not found: {path}")
        sys.exit(0)

    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def save_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def decode_bitstring(bitstring: str, num_bits_per_asset: int, num_assets: int) -> List[float]:
    bitstring = bitstring.strip()
    if not bitstring:
        raise ValueError("Empty bitstring")

    bits = [int(c) for c in reversed(bitstring)]
    total_bits = num_bits_per_asset * num_assets
    if len(bits) != total_bits:
        raise ValueError(
            f"Bitstring length {len(bits)} does not match expected {total_bits}"
        )

    denom = float((1 << num_bits_per_asset) - 1)
    weights: List[float] = []

    for asset_idx in range(num_assets):
        start = asset_idx * num_bits_per_asset
        value = 0
        for bit_idx in range(num_bits_per_asset):
            value += bits[start + bit_idx] << bit_idx
        weights.append(value / denom)

    return weights


def objective_value(q_matrix: List[List[float]], bitstring: str) -> float:
    bits = [int(c) for c in reversed(bitstring.strip())]
    vector = [float(b) for b in bits]
    total = 0.0
    for i, xi in enumerate(vector):
        for j, xj in enumerate(vector):
            total += xi * q_matrix[i][j] * xj
    return total


def build_qaoa_circuit(
    num_qubits: int,
    q_matrix: List[List[float]],
    qaoa_depth: int,
    gamma: float,
    beta: float,
) -> QuantumCircuit:
    circuit = QuantumCircuit(num_qubits, num_qubits)

    for qubit in range(num_qubits):
        circuit.h(qubit)

    for _ in range(qaoa_depth):
        for i in range(num_qubits):
            for j in range(i, num_qubits):
                value = q_matrix[i][j]
                if abs(value) < 1e-10:
                    continue
                angle = 2.0 * gamma * value
                if i == j:
                    circuit.rz(angle, i)
                else:
                    circuit.rzz(angle, i, j)

        for qubit in range(num_qubits):
            circuit.rx(2.0 * beta, qubit)

    circuit.measure(range(num_qubits), range(num_qubits))
    return circuit


def select_best_bitstring(counts: Dict[str, int], q_matrix: List[List[float]]) -> str:
    """Return the bitstring with the lowest QUBO objective value among measured samples."""
    return min(counts.keys(), key=lambda bs: objective_value(q_matrix, bs))


def _qaoa_expected_energy(
    params: List[float],
    num_qubits: int,
    q_matrix: List[List[float]],
    qaoa_depth: int,
    backend_instance: Any,
    shots: int,
) -> float:
    """Expected QUBO energy E[H] over the measurement distribution for given (gamma, beta)."""
    gamma, beta = float(params[0]), float(params[1])
    circuit = build_qaoa_circuit(num_qubits, q_matrix, qaoa_depth, gamma, beta)
    transpiled = transpile(circuit, backend_instance)
    result = backend_instance.run(transpiled, shots=shots).result()
    counts = result.get_counts()
    total = sum(counts.values())
    return sum((cnt / total) * objective_value(q_matrix, bs) for bs, cnt in counts.items())


def optimize_qaoa_parameters(
    num_qubits: int,
    q_matrix: List[List[float]],
    qaoa_depth: int,
    backend_instance: Any,
    shots: int = 512,
    max_iter: int = 100,
) -> Tuple[float, float]:
    """Optimize QAOA variational parameters (gamma, beta) via COBYLA.

    Minimizes the expected QUBO energy over the measurement distribution.
    Returns (gamma_opt, beta_opt).
    """
    try:
        from scipy.optimize import minimize
    except ImportError:
        print("Warning: scipy not available; using default QAOA parameters.", file=sys.stderr)
        return 0.5, 0.5

    result = minimize(
        _qaoa_expected_energy,
        x0=[0.5, 0.5],
        args=(num_qubits, q_matrix, qaoa_depth, backend_instance, shots),
        method="COBYLA",
        options={"maxiter": max_iter, "rhobeg": 0.5},
    )
    return float(result.x[0]), float(result.x[1])


def extract_execution_us(result: Any, job: Any) -> Tuple[float, str]:
    """Extract QPU execution time in microseconds.

    Tries multiple known metadata keys in priority order across different
    IBM Runtime versions, then falls back to job-level timestamps.

    Args:
        result: PrimitiveResult or similar result object from the job.
        job: IBM runtime job object (for fallback metrics).

    Returns:
        Tuple of (time_us, key_name_used).
    """
    # Try execution_spans first (SamplerV2, current API)
    try:
        spans = result.metadata['execution']['execution_spans']
        if spans:
            delta = spans[0].stop - spans[0].start
            us = delta.total_seconds() * 1e6
            return us, 'execution_spans'
    except (KeyError, AttributeError, IndexError, TypeError):
        pass

    # Try time_taken (older backends, returns seconds)
    try:
        return result.metadata['time_taken'] * 1e6, 'time_taken'
    except (KeyError, TypeError):
        pass

    # Try job-level timing as fallback
    try:
        metrics = job.metrics()
        if 'timestamps' in metrics:
            t = metrics['timestamps']
            if 'running' in t and 'finished' in t:
                fmt = "%Y-%m-%dT%H:%M:%S.%fZ"
                start = datetime.strptime(t['running'], fmt)
                end = datetime.strptime(t['finished'], fmt)
                return (end - start).total_seconds() * 1e6, 'job_metrics'
    except Exception:
        pass

    return -1.0, 'unavailable'


def extract_counts(pub_result: Any) -> Dict[str, int]:
    """Extract bitstring counts from a SamplerV2 PubResult.

    Handles varying classical register names across IBM backend versions.

    Args:
        pub_result: SamplerV2 PubResult object.

    Returns:
        Dict mapping bitstring to count.

    Raises:
        ValueError: If no BitArray field is found in PubResult data.
    """
    # Try named register 'meas' first
    try:
        return pub_result.data.meas.get_counts()
    except AttributeError:
        pass

    # Try register 'c'
    try:
        return pub_result.data.c.get_counts()
    except AttributeError:
        pass

    # Iterate all fields in data to find a BitArray
    try:
        from qiskit.primitives.containers.bit_array import BitArray
        for field_name in pub_result.data.__dataclass_fields__:
            field = getattr(pub_result.data, field_name)
            if isinstance(field, BitArray):
                return field.get_counts()
    except Exception:
        pass

    raise ValueError(
        "Could not find BitArray in PubResult data fields: "
        + str(list(getattr(pub_result.data, "__dataclass_fields__", {}).keys()))
    )


def _compute_signal_quality(counts: Dict[str, int]) -> Tuple[float, str]:
    """Compute top-bitstring fraction and signal quality label.

    Args:
        counts: Bitstring measurement counts.

    Returns:
        Tuple of (top_fraction, quality) where quality is 'ok' or 'low'.
    """
    total = sum(counts.values())
    top = max(counts.values()) if counts else 0
    fraction = top / total if total > 0 else 0.0
    quality = "low" if fraction < 0.05 else "ok"
    return fraction, quality


def collect_aer_job(job: Dict[str, Any], results_dir: Path) -> Dict[str, Any]:
    mode = job.get("mode", "qaoa")
    if mode == "qamo":
        return collect_aer_qamo_job(job, results_dir)
    if mode == "qamoo":
        return collect_aer_qamoo_job(job, results_dir)

    problem_file = Path(job["problem_file"])
    if not problem_file.exists():
        raise FileNotFoundError(f"Problem file not found: {problem_file}")

    with problem_file.open("r", encoding="utf-8") as f:
        problem = json.load(f)

    qaoa_depth = int(job.get("circuit_params", {}).get("qaoa_depth", job.get("qaoa_depth", 1)))
    num_qubits = int(problem["num_variables"])
    num_assets = int(problem["num_assets"])
    num_bits_per_asset = int(problem["num_bits_per_asset"])
    q_matrix = problem["Q"]
    shots = int(job["shots"])
    solver_name = f"qaoa_p{qaoa_depth}_{job['backend']}"

    try:
        from qiskit_aer import AerSimulator
        backend_instance = AerSimulator()

        # Optimize variational parameters via COBYLA
        print(f"  Optimizing QAOA parameters ({num_qubits} qubits, p={qaoa_depth})...")
        gamma_opt, beta_opt = optimize_qaoa_parameters(
            num_qubits, q_matrix, qaoa_depth, backend_instance,
            shots=min(shots, 512), max_iter=100,
        )
        print(f"  Optimal params: gamma={gamma_opt:.4f} beta={beta_opt:.4f}")

        # Final measurement with optimized parameters
        circuit = build_qaoa_circuit(num_qubits, q_matrix, qaoa_depth, gamma_opt, beta_opt)
        pre_depth = int(circuit.depth())
        transpiled = transpile(circuit, backend_instance)
        post_depth = int(transpiled.depth())

        start_time = time.time()
        result = backend_instance.run(transpiled, shots=shots).result()
        end_time = time.time()
        counts = result.get_counts()

        top_fraction, signal_quality = _compute_signal_quality(counts)

        # Select the sample with the lowest QUBO objective value
        bitstring = select_best_bitstring(counts, q_matrix)
        raw_weights = decode_bitstring(bitstring, num_bits_per_asset, num_assets)

        # Normalize so weights sum to 1; fall back to equal-weight if all-zero bitstring
        total_weight = sum(raw_weights)
        if total_weight > 1e-9:
            weights = [w / total_weight for w in raw_weights]
        else:
            weights = [1.0 / num_assets] * num_assets

        obj = objective_value(q_matrix, bitstring)
        convergence_info = (
            f"gamma_opt={gamma_opt:.4f} beta_opt={beta_opt:.4f} "
            f"obj={obj:.6f} bitstring_count={counts[bitstring]}/{shots}"
        )

        return {
            "schema_version": "1.0",
            "job_id": job["job_id"],
            "problem_id": problem["problem_id"],
            "solver_name": solver_name,
            "execution_backend": job["backend"],
            "weights": weights,
            "bitstring": bitstring,
            "objective_value": obj,
            "pre_transpilation_depth": pre_depth,
            "circuit_depth": post_depth,
            "circuit_execution_us": -1.0,
            "metadata_key_used": "n/a",
            "shots": shots,
            "solve_time_ms": float((end_time - start_time) * 1000.0),
            "convergence_info": convergence_info,
            "backend_calibration_date": "n/a",
            "error_mitigation_method": "none",
            "physical_qubits": job.get("physical_qubits", list(range(num_qubits))),
            "top_bitstring_fraction": top_fraction,
            "signal_quality": signal_quality,
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "status": "COMPLETED",
        }
    except Exception as exc:
        raise RuntimeError(f"Aer collection failed: {exc}")


def collect_aer_qamo_job(job: Dict[str, Any], results_dir: Path) -> Dict[str, Any]:
    """Collect an Aer QAMO job: optimise with COBYLA, return a single portfolio."""
    problem_file = Path(job["problem_file"])
    if not problem_file.exists():
        raise FileNotFoundError(f"Problem file not found: {problem_file}")

    with problem_file.open("r", encoding="utf-8") as f:
        problem = json.load(f)

    num_qubits = int(problem["num_variables"])
    num_assets = int(problem["num_assets"])
    num_bits_per_asset = int(problem["num_bits_per_asset"])
    q_matrix = problem["Q"]
    shots = int(job["shots"])

    try:
        from qiskit_aer import AerSimulator
        from qiskit_solver import (
            build_qamo_circuit,
            _decode_bitstring_np,
            _objective_np,
            _execute_aer,
            _optimize_qamo_params,
            solve_mean_field,
        )

        backend_instance = AerSimulator()
        Q_np = np.array(q_matrix)
        mf_config: Dict[str, Any] = {"max_mf_iterations": 50, "mf_convergence_tol": 1e-6}

        print(f"  Optimising QAMO parameters ({num_qubits} qubits, up to 4 restarts)...", file=sys.stderr)
        gamma_opt, beta_opt, _opt_counts = _optimize_qamo_params(
            Q_np, mf_config, backend_instance, min(shots, 512)
        )

        print(f"  Optimal params: gamma={gamma_opt:.4f} beta={beta_opt:.4f}", file=sys.stderr)

        circuit, theta = build_qamo_circuit(
            Q_np, gamma_opt, beta_opt,
            mf_config["max_mf_iterations"],
            mf_config["mf_convergence_tol"],
        )
        pre_depth = int(circuit.depth())
        transpiled = transpile(circuit, backend_instance)
        post_depth = int(transpiled.depth())

        start_time = time.time()
        counts = _execute_aer(circuit, backend_instance, shots)
        end_time = time.time()

        top_fraction, signal_quality = _compute_signal_quality(counts)
        bitstring = select_best_bitstring(counts, q_matrix)
        raw_weights = _decode_bitstring_np(bitstring, num_assets, num_bits_per_asset)
        total_weight = float(raw_weights.sum())
        if total_weight > 1e-9:
            weights = (raw_weights / total_weight).tolist()
        else:
            weights = [1.0 / num_assets] * num_assets

        # Mean-field diagnostics
        m, mf_converged, mf_iters = solve_mean_field(
            Q_np, beta_opt,
            mf_config["max_mf_iterations"],
            mf_config["mf_convergence_tol"],
        )
        import numpy as _np
        per_qubit_angles = (2.0 * _np.arctan(m + 0.5)).tolist()

        obj = _objective_np(Q_np, bitstring)
        solver_name = f"qamo_p1_{job['backend']}"

        return {
            "schema_version": "1.0",
            "job_id": job["job_id"],
            "problem_id": problem["problem_id"],
            "solver_name": solver_name,
            "execution_backend": job["backend"],
            "mode": "qamo",
            "weights": weights,
            "bitstring": bitstring,
            "objective_value": obj,
            "pre_transpilation_depth": pre_depth,
            "circuit_depth": post_depth,
            "circuit_execution_us": -1.0,
            "metadata_key_used": "n/a",
            "shots": shots,
            "solve_time_ms": float((end_time - start_time) * 1000.0),
            "convergence_info": (
                f"gamma_opt={gamma_opt:.4f} beta_opt={beta_opt:.4f} "
                f"obj={obj:.6f}"
            ),
            "backend_calibration_date": "n/a",
            "error_mitigation_method": "none",
            "physical_qubits": job.get("physical_qubits", list(range(num_qubits))),
            "top_bitstring_fraction": top_fraction,
            "signal_quality": signal_quality,
            "mf_converged": bool(mf_converged),
            "mf_iterations": int(mf_iters),
            "per_qubit_angles": per_qubit_angles,
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "status": "COMPLETED",
        }
    except Exception as exc:
        raise RuntimeError(f"Aer QAMO collection failed: {exc}")


def collect_aer_qamoo_job(job: Dict[str, Any], results_dir: Path) -> Dict[str, Any]:
    """Collect an Aer QAMOO job: sweep lambdas, return Pareto frontier."""
    problem_file = Path(job["problem_file"])
    if not problem_file.exists():
        raise FileNotFoundError(f"Problem file not found: {problem_file}")

    with problem_file.open("r", encoding="utf-8") as f:
        problem = json.load(f)

    if "covariance" not in problem:
        raise RuntimeError(
            "Problem JSON is missing 'covariance' field required for QAMOO. "
            "Re-submit with --mode qamoo to augment the problem file."
        )

    num_qubits = int(problem["num_variables"])
    num_assets = int(problem["num_assets"])
    shots = int(job["shots"])
    n_pts = int(job.get("n_frontier_points", 20))
    lambda_min = float(job.get("lambda_min", 0.1))
    lambda_max = float(job.get("lambda_max", 20.0))

    try:
        from qiskit_aer import AerSimulator
        from qiskit_solver import run_qamoo_sweep

        backend_instance = AerSimulator()
        qamoo_config: Dict[str, Any] = {
            "n_frontier_points": n_pts,
            "lambda_min": lambda_min,
            "lambda_max": lambda_max,
            "max_mf_iterations": 50,
            "mf_convergence_tol": 1e-6,
        }

        lambdas = list(
            np.logspace(np.log10(lambda_min), np.log10(lambda_max), n_pts).tolist()
        )
        print(
            f"  Running QAMOO sweep: {n_pts} lambda points on {num_qubits} qubits...",
            file=sys.stderr,
        )
        start_time = time.time()
        frontier_points, best_weights = run_qamoo_sweep(
            problem, qamoo_config, backend_instance, shots
        )
        end_time = time.time()

        # Backward-compat: best single solution (highest Sharpe from sweep)
        total_w = float(best_weights.sum())
        if total_w > 1e-9:
            best_weights_list = (best_weights / total_w).tolist()
        else:
            best_weights_list = [1.0 / num_assets] * num_assets

        best_bitstring = (
            frontier_points[0]["bitstring"] if frontier_points else ""
        )
        best_obj = frontier_points[0]["objective_value"] if frontier_points else 0.0

        solver_name = f"qamoo_{n_pts}pts_{job['backend']}"

        return {
            "schema_version": "1.0",
            "job_id": job["job_id"],
            "problem_id": problem["problem_id"],
            "solver_name": solver_name,
            "execution_backend": job["backend"],
            "mode": "qamoo",
            "weights": best_weights_list,
            "bitstring": best_bitstring,
            "objective_value": best_obj,
            "pre_transpilation_depth": job.get("pre_transpilation_depth", 0),
            "circuit_depth": job.get("circuit_depth", 0),
            "circuit_execution_us": -1.0,
            "metadata_key_used": "n/a",
            "shots": shots,
            "solve_time_ms": float((end_time - start_time) * 1000.0),
            "convergence_info": (
                f"QAMOO sweep: {n_pts} points requested, "
                f"{len(frontier_points)} returned after deduplication"
            ),
            "backend_calibration_date": "n/a",
            "error_mitigation_method": "none",
            "physical_qubits": job.get("physical_qubits", list(range(num_qubits))),
            "top_bitstring_fraction": 0.0,
            "signal_quality": "ok" if frontier_points else "low",
            "frontier_points": frontier_points,
            "n_frontier_points_requested": n_pts,
            "n_frontier_points_returned": len(frontier_points),
            "sweep_lambdas": lambdas,
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "status": "COMPLETED",
        }
    except Exception as exc:
        raise RuntimeError(f"Aer QAMOO collection failed: {exc}")


def collect_ibm_job(job: Dict[str, Any], timeout_min: int, poll_interval: int) -> Dict[str, Any]:
    token = os.environ.get("IBM_QUANTUM_TOKEN")
    instance = os.environ.get("IBM_QUANTUM_INSTANCE", "")
    if not token:
        raise RuntimeError("IBM_QUANTUM_TOKEN is not set. Cannot collect IBM jobs.")

    with Path(job["problem_file"]).open("r", encoding="utf-8") as f:
        problem = json.load(f)

    mode = job.get("mode", "qaoa")
    qaoa_depth = int(job.get("qaoa_depth", 1))
    if mode == "qamo":
        solver_name = f"qamo_p1_{job['backend']}"
    elif mode == "qamoo":
        n_pts = int(job.get("n_frontier_points", 20))
        solver_name = f"qamoo_{n_pts}pts_{job['backend']}"
    else:
        solver_name = f"qaoa_p{qaoa_depth}_{job['backend']}"
    circuit_params = job.get("circuit_params", {})
    pre_depth = job.get("pre_transpilation_depth", 0)
    post_depth = job.get("circuit_depth", 0)
    physical_qubits = job.get("physical_qubits", list(range(int(problem["num_variables"]))))
    error_mitigation_requested = bool(circuit_params.get("error_mitigation", False))

    _phase4_defaults: Dict[str, Any] = {
        "pre_transpilation_depth": pre_depth,
        "circuit_depth": post_depth,
        "metadata_key_used": "unavailable",
        "backend_calibration_date": "unavailable",
        "error_mitigation_method": "none",
        "physical_qubits": physical_qubits,
        "top_bitstring_fraction": 0.0,
        "signal_quality": "low",
    }

    try:
        from qiskit_ibm_runtime import QiskitRuntimeService

        service = (
            QiskitRuntimeService(channel="ibm_quantum_platform", token=token, instance=instance)
            if instance
            else QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
        )
        backend = service.backend(job["backend"])
        runner_job = service.job(job["job_id"])
        print(f"  IBM job {job['job_id']} retrieved from {job['backend']}", file=sys.stderr)
    except Exception as exc:
        raise RuntimeError(f"IBM job lookup failed: {exc}")

    # Retrieve calibration date
    cal_date = "unavailable"
    try:
        cal_date = backend.properties().last_update_date.isoformat()
        print(f"  Backend calibration date: {cal_date}", file=sys.stderr)
    except Exception as exc:
        print(f"  Warning: could not retrieve calibration date: {exc}", file=sys.stderr)

    # Poll until terminal status
    start_time = time.time()
    deadline = start_time + timeout_min * 60
    status = None
    wait = poll_interval
    while time.time() < deadline:
        try:
            status_obj = runner_job.status()
            status = getattr(status_obj, "name", str(status_obj))
        except Exception:
            status = str(runner_job.status())

        print(f"  Job status: {status}", file=sys.stderr)
        if status in ("DONE", "COMPLETED"):
            break
        if status in ("ERROR", "FAILED", "CANCELLED"):
            break
        time.sleep(wait)
        wait = min(wait * 1.5, 300)

    if status is None:
        status = "TIMEOUT"

    if status in ("ERROR", "FAILED"):
        return {
            "schema_version": "1.0",
            "job_id": job["job_id"],
            "problem_id": problem["problem_id"],
            "solver_name": solver_name,
            "execution_backend": job["backend"],
            "weights": [],
            "bitstring": "",
            "objective_value": 0.0,
            "circuit_execution_us": -1.0,
            "shots": job["shots"],
            "solve_time_ms": 0.0,
            "convergence_info": f"Job failed with status {status}",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "status": "FAILED",
            **_phase4_defaults,
            "backend_calibration_date": cal_date,
        }

    if status == "TIMEOUT":
        return {
            "schema_version": "1.0",
            "job_id": job["job_id"],
            "problem_id": problem["problem_id"],
            "solver_name": solver_name,
            "execution_backend": job["backend"],
            "weights": [],
            "bitstring": "",
            "objective_value": 0.0,
            "circuit_execution_us": -1.0,
            "shots": job["shots"],
            "solve_time_ms": 0.0,
            "convergence_info": "Timeout reached while polling IBM job.",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "status": "TIMEOUT",
            **_phase4_defaults,
            "backend_calibration_date": cal_date,
        }

    # DONE — extract counts
    try:
        result = runner_job.result()
    except Exception as exc:
        return {
            "schema_version": "1.0",
            "job_id": job["job_id"],
            "problem_id": problem["problem_id"],
            "solver_name": solver_name,
            "execution_backend": job["backend"],
            "weights": [],
            "bitstring": "",
            "objective_value": 0.0,
            "circuit_execution_us": -1.0,
            "shots": job["shots"],
            "solve_time_ms": 0.0,
            "convergence_info": f"IBM result retrieval failed: {exc}",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "status": "FAILED",
            **_phase4_defaults,
            "backend_calibration_date": cal_date,
        }

    # Robust counts extraction (SamplerV2 PrimitiveResult or older API)
    counts: Dict[str, int] = {}
    try:
        pub_result = result[0]
        counts = extract_counts(pub_result)
        print(f"  Extracted {sum(counts.values())} shots from PubResult", file=sys.stderr)
    except Exception:
        try:
            counts = result.get_counts() if hasattr(result, "get_counts") else {}
        except Exception:
            counts = {}
        if counts:
            print(f"  Extracted {sum(counts.values())} shots via legacy API", file=sys.stderr)
        else:
            print("  Warning: could not extract counts from result", file=sys.stderr)

    # mthree readout error correction
    error_mitigation_method = "none"
    if error_mitigation_requested and counts:
        try:
            import mthree
            mit = mthree.M3Mitigation(backend)
            print(f"  Calibrating mthree for qubits {physical_qubits}...", file=sys.stderr)
            mit.cals_from_system(physical_qubits)
            quasi_probs = mit.apply_correction(counts, physical_qubits)
            shots_total = job["shots"]
            corrected = {k: max(0, int(v * shots_total)) for k, v in quasi_probs.items()}
            # Keep only non-zero entries
            corrected = {k: v for k, v in corrected.items() if v > 0}
            if corrected:
                counts = corrected
                error_mitigation_method = "mthree"
                print("  mthree correction applied.", file=sys.stderr)
            else:
                print("  Warning: mthree produced empty corrected counts; using raw.", file=sys.stderr)
        except Exception as exc:
            print(f"  Warning: mthree correction failed: {exc}", file=sys.stderr)
            error_mitigation_method = "failed"

    # Signal quality
    top_fraction, signal_quality = _compute_signal_quality(counts)
    if signal_quality == "low":
        print(
            f"  Warning: low signal quality (top bitstring = {top_fraction:.1%} of shots). "
            f"Results reflect NISQ hardware noise.",
            file=sys.stderr,
        )

    # Decode best solution
    if counts:
        bitstring = select_best_bitstring(counts, problem["Q"])
        raw_weights = decode_bitstring(
            bitstring, int(problem["num_bits_per_asset"]), int(problem["num_assets"])
        )
        total_weight = sum(raw_weights)
        weights = (
            [w / total_weight for w in raw_weights]
            if total_weight > 1e-9
            else [1.0 / int(problem["num_assets"])] * int(problem["num_assets"])
        )
        convergence_info = (
            f"min-energy bitstring count: {counts.get(bitstring, 0)}/{job['shots']} "
            f"signal={signal_quality} em={error_mitigation_method}"
        )
        obj = objective_value(problem["Q"], bitstring)
    else:
        bitstring = ""
        weights = []
        convergence_info = "IBM job completed but no counts returned."
        obj = 0.0

    # Execution timing
    execution_us, metadata_key = extract_execution_us(result, runner_job)
    print(f"  Execution time: {execution_us:.1f} µs (key={metadata_key})", file=sys.stderr)

    return {
        "schema_version": "1.0",
        "job_id": job["job_id"],
        "problem_id": problem["problem_id"],
        "solver_name": solver_name,
        "execution_backend": job["backend"],
        "weights": weights,
        "bitstring": bitstring,
        "objective_value": obj,
        "pre_transpilation_depth": pre_depth,
        "circuit_depth": post_depth,
        "circuit_execution_us": execution_us,
        "metadata_key_used": metadata_key,
        "shots": job["shots"],
        "solve_time_ms": float((time.time() - start_time) * 1000.0),
        "convergence_info": convergence_info,
        "backend_calibration_date": cal_date,
        "error_mitigation_method": error_mitigation_method,
        "physical_qubits": physical_qubits,
        "top_bitstring_fraction": top_fraction,
        "signal_quality": signal_quality,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETED",
    }


def main() -> int:
    args = parse_args()
    jobs_path = Path(args.jobs_file)
    results_dir = Path(args.results_dir)
    data = load_json(jobs_path)

    jobs = data.get("jobs", [])
    pending = [job for job in jobs if job.get("status") not in {"COMPLETED", "FAILED", "TIMEOUT"}]

    if not pending:
        print("All jobs already resolved.")
        return 0

    if args.dry_run:
        for job in pending:
            if job["backend"] == "aer_simulator":
                print(f"Job {job['job_id']}: Aer simulator - COMPLETED")
            elif job["backend"].startswith("ibm_"):
                token = os.environ.get("IBM_QUANTUM_TOKEN")
                instance = os.environ.get("IBM_QUANTUM_INSTANCE", "")
                if not token:
                    print(f"Job {job['job_id']}: IBM - TOKEN MISSING")
                    continue
                try:
                    from qiskit_ibm_runtime import QiskitRuntimeService
                    service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token, instance=instance)
                    runner_job = service.job(job["job_id"])
                    status = runner_job.status().name
                    print(f"Job {job['job_id']}: IBM - {status}")
                except Exception as exc:
                    print(f"Job {job['job_id']}: IBM - ERROR: {exc}")
        return 0

    collected = 0
    failed = 0

    for job in pending:
        try:
            if job["backend"] == "aer_simulator":
                result_payload = collect_aer_job(job, results_dir)
            elif job["backend"].startswith("ibm_"):
                result_payload = collect_ibm_job(job, args.timeout_min, args.poll_interval)
            else:
                raise ValueError(f"Unsupported backend: {job['backend']}")

            result_file = results_dir / f"quantum_result_{job['job_id']}.json"
            save_json(result_file, result_payload)
            job["status"] = result_payload["status"]
            if result_payload["status"] != "COMPLETED":
                failed += 1
            collected += 1
        except Exception as exc:
            job["status"] = "FAILED"
            job["convergence_info"] = str(exc)
            failed += 1
            collected += 1

    data["jobs"] = jobs
    save_json(jobs_path, data)

    print(f"Collected {collected}/{len(jobs)} jobs. {failed} failed or timed out.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
