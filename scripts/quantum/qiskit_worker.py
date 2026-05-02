#!/usr/bin/env python3
"""Persistent Qiskit worker process for portfolio optimization requests."""

from __future__ import annotations

import json
import os
import sys
import time
import traceback
from pathlib import Path
from typing import Any, Dict, List

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import numpy as np
from qiskit import transpile

import qiskit_submit
import qiskit_collect


def normalize_weights(weights: List[float]) -> List[float]:
    total = sum(weights)
    if total > 1e-9:
        return [w / total for w in weights]
    if weights:
        return [1.0 / len(weights)] * len(weights)
    return []


def select_best_bitstring(counts: Dict[str, int], q_matrix: List[List[float]]) -> str:
    if not counts:
        raise ValueError("No measurement counts available.")
    return min(counts.keys(), key=lambda bs: qiskit_collect.objective_value(q_matrix, bs))


def _resolve_ibm_backend(request: Dict[str, Any]) -> str:
    """Resolve 'ibm_auto' to the best available IBM backend by queue length."""
    token = os.environ.get("IBM_QUANTUM_TOKEN")
    instance = os.environ.get("IBM_QUANTUM_INSTANCE", "")
    if not token:
        raise RuntimeError("IBM_QUANTUM_TOKEN is not set.")
    from qiskit_ibm_runtime import QiskitRuntimeService
    service = (
        QiskitRuntimeService(channel="ibm_quantum_platform", token=token, instance=instance)
        if instance
        else QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
    )
    min_qubits = int(request.get("ibm_min_qubits", 20))
    candidates = [
        b for b in service.backends()
        if b.num_qubits >= min_qubits and b.status().operational
    ]
    if not candidates:
        raise RuntimeError(
            f"No operational IBM backends with >= {min_qubits} qubits"
        )
    selected = min(candidates, key=lambda b: b.status().pending_jobs)
    print(f"  ibm_auto selected: {selected.name}", file=sys.stderr)
    return selected.name


def _run_qamoo_sweep(request: Dict[str, Any], problem: Dict[str, Any], shots: int) -> Dict[str, Any]:
    """Run a QAMOO lambda sweep locally on Aer, regardless of target backend.

    IBM hardware QAMOO would require N separate submissions (one per lambda),
    making it impractical. The sweep always runs synchronously on Aer.
    """
    from qiskit_aer import AerSimulator
    from qiskit_solver import run_qamoo_sweep

    if "covariance" not in problem:
        raise RuntimeError(
            "QAMOO requires 'covariance' in the problem file. "
            "Set augment_problem_data: true in the solver config."
        )

    n_frontier_points = int(request.get("n_frontier_points", 20))
    qamoo_config: Dict[str, Any] = {
        "n_frontier_points": n_frontier_points,
        "lambda_min": 0.1,
        "lambda_max": 20.0,
        "max_mf_iterations": 50,
        "mf_convergence_tol": 1e-6,
    }

    num_assets = int(problem["num_assets"])
    backend_instance = AerSimulator()

    start_time = time.time()
    frontier_points, best_weights = run_qamoo_sweep(problem, qamoo_config, backend_instance, shots)
    end_time = time.time()

    total_w = float(best_weights.sum())
    weights = (best_weights / total_w).tolist() if total_w > 1e-9 else [1.0 / num_assets] * num_assets

    return {
        "status": "ok",
        "job_id": request.get("request_id", ""),
        "execution_backend": request.get("backend", "aer_simulator"),
        "weights": weights,
        "solve_time_ms": float((end_time - start_time) * 1000.0),
        "circuit_execution_us": -1.0,
        "convergence_info": (
            f"QAMOO sweep: {n_frontier_points} points requested, "
            f"{len(frontier_points)} returned after deduplication"
        ),
        "signal_quality": "ok" if frontier_points else "low",
    }


def run_aer_request(request: Dict[str, Any], problem: Dict[str, Any]) -> Dict[str, Any]:
    from qiskit_aer import AerSimulator

    q_matrix = problem["Q"]
    num_variables = int(problem["num_variables"])
    num_assets = int(problem["num_assets"])
    num_bits_per_asset = int(problem["num_bits_per_asset"])
    shots = int(request["shots"])
    mode = request["mode"]
    qaoa_depth = int(request["qaoa_depth"])
    optimization_level = int(request.get("optimization_level", 1))

    if mode == "qamoo":
        return _run_qamoo_sweep(request, problem, shots)

    if mode == "qamo":
        from qiskit_solver import build_qamo_circuit

        Q_np = np.array(q_matrix)
        circuit, _ = build_qamo_circuit(Q_np, 0.1, 0.1, 50, 1e-6)
    else:
        circuit = qiskit_submit.build_qaoa_circuit(num_variables, q_matrix, qaoa_depth, 0.1, 0.1)

    backend = AerSimulator()
    transpiled = transpile(circuit, optimization_level=optimization_level)

    start_time = time.time()
    job = backend.run(transpiled, shots=shots)
    result = job.result()
    end_time = time.time()

    counts = result.get_counts()
    bitstring = select_best_bitstring(counts, q_matrix)
    raw_weights = qiskit_collect.decode_bitstring(bitstring, num_bits_per_asset, num_assets)
    weights = normalize_weights(raw_weights)
    top_fraction, signal_quality = qiskit_collect._compute_signal_quality(counts)
    convergence_info = (
        f"bitstring_count={counts[bitstring]}/{shots} "
        f"signal={signal_quality}"
    )

    return {
        "status": "ok",
        "job_id": request.get("request_id", ""),
        "execution_backend": request["backend"],
        "weights": weights,
        "solve_time_ms": float((end_time - start_time) * 1000.0),
        "circuit_execution_us": -1.0,
        "convergence_info": convergence_info,
        "signal_quality": signal_quality,
    }


def run_ibm_request(request: Dict[str, Any], problem: Dict[str, Any]) -> Dict[str, Any]:
    # QAMOO lambda sweep always runs on Aer — N IBM submissions per sweep is impractical.
    if request.get("mode") == "qamoo":
        return _run_qamoo_sweep(request, problem, int(request["shots"]))

    token = os.environ.get("IBM_QUANTUM_TOKEN")
    instance = os.environ.get("IBM_QUANTUM_INSTANCE", "")
    if not token:
        raise RuntimeError("IBM_QUANTUM_TOKEN is not set.")

    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2

    q_matrix = problem["Q"]
    num_variables = int(problem["num_variables"])
    num_assets = int(problem["num_assets"])
    num_bits_per_asset = int(problem["num_bits_per_asset"])
    shots = int(request["shots"])
    mode = request["mode"]
    qaoa_depth = int(request["qaoa_depth"])
    backend_name = request["backend"]
    optimization_level = int(request.get("optimization_level", 1))

    if mode == "qamo":
        from qiskit_solver import build_qamo_circuit

        Q_np = np.array(q_matrix)
        circuit, _ = build_qamo_circuit(Q_np, 0.1, 0.1, 50, 1e-6)
    else:
        circuit = qiskit_submit.build_qaoa_circuit(num_variables, q_matrix, qaoa_depth, 0.1, 0.1)

    service = (
        QiskitRuntimeService(channel="ibm_quantum_platform", token=token, instance=instance)
        if instance
        else QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
    )
    backend = service.backend(backend_name)
    transpiled = transpile(circuit, backend=backend, optimization_level=optimization_level)
    sampler = SamplerV2(mode=backend)

    start_time = time.time()
    job = sampler.run([transpiled], shots=shots)
    result = job.result()
    end_time = time.time()

    pub_result = result[0]
    counts = qiskit_collect.extract_counts(pub_result)
    bitstring = select_best_bitstring(counts, q_matrix)
    raw_weights = qiskit_collect.decode_bitstring(bitstring, num_bits_per_asset, num_assets)
    weights = normalize_weights(raw_weights)
    top_fraction, signal_quality = qiskit_collect._compute_signal_quality(counts)
    convergence_info = (
        f"bitstring_count={counts[bitstring]}/{shots} "
        f"signal={signal_quality}"
    )
    execution_us, _ = qiskit_collect.extract_execution_us(pub_result, job)

    return {
        "status": "ok",
        "job_id": job.job_id(),
        "execution_backend": backend_name,
        "weights": weights,
        "solve_time_ms": float((end_time - start_time) * 1000.0),
        "circuit_execution_us": execution_us,
        "convergence_info": convergence_info,
        "signal_quality": signal_quality,
    }


def handle_request(request: Dict[str, Any]) -> Dict[str, Any]:
    try:
        mode = request.get("mode", "qaoa")
        backend = request.get("backend", "aer_simulator")
        problem_file = Path(request["problem_file"])
        if not problem_file.exists():
            raise FileNotFoundError(f"Problem file not found: {problem_file}")

        problem = qiskit_submit.load_problem(problem_file)
        if request.get("augment_problem_data", False):
            qiskit_submit._augment_problem_data(problem_file, problem)

        if backend == "ibm_auto":
            # QAMOO always runs on Aer regardless of backend — skip IBM resolution.
            if mode == "qamoo":
                return _run_qamoo_sweep(request, problem, int(request.get("shots", 1024)))
            backend = _resolve_ibm_backend(request)
            request = {**request, "backend": backend}

        if backend == "aer_simulator":
            return run_aer_request(request, problem)
        if backend.startswith("ibm_"):
            return run_ibm_request(request, problem)

        raise RuntimeError(f"Unsupported backend: {backend}")
    except SystemExit as exc:
        raise RuntimeError(f"Worker request aborted: {exc}")


def main() -> int:
    sys.stdout.reconfigure(line_buffering=True)
    for raw_line in sys.stdin:
        line = raw_line.strip()
        if not line:
            continue
        try:
            request = json.loads(line)
            response = handle_request(request)
        except Exception as exc:
            response = {
                "status": "error",
                "error": str(exc),
                "traceback": traceback.format_exc(),
                "request_id": request.get("request_id") if isinstance(request, dict) else None,
            }
        sys.stdout.write(json.dumps(response) + "\n")
        sys.stdout.flush()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
