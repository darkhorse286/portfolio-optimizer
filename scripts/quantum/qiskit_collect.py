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
from typing import Any, Dict, List

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


def select_top_bitstring(counts: Dict[str, int]) -> str:
    return max(counts.items(), key=lambda item: item[1])[0]


def result_metadata_time_us(result: Any) -> float:
    metadata = getattr(result, "metadata", {}) or {}
    time_taken = metadata.get("time_taken")
    if isinstance(time_taken, (int, float)):
        return float(time_taken) * 1e6
    return -1.0


def collect_aer_job(job: Dict[str, Any], results_dir: Path) -> Dict[str, Any]:
    problem_file = Path(job["problem_file"])
    if not problem_file.exists():
        raise FileNotFoundError(f"Problem file not found: {problem_file}")

    with problem_file.open("r", encoding="utf-8") as f:
        problem = json.load(f)

    params = job.get("circuit_params", {})
    gamma = float(params.get("gamma", 0.1))
    beta = float(params.get("beta", 0.1))
    qaoa_depth = int(params.get("qaoa_depth", job.get("qaoa_depth", 1)))

    num_qubits = int(problem["num_variables"])
    q_matrix = problem["Q"]
    shots = int(job["shots"])
    solver_name = f"qaoa_p{qaoa_depth}_{job['backend']}"

    try:
        from qiskit_aer import AerSimulator

        backend_instance = AerSimulator()
        circuit = build_qaoa_circuit(num_qubits, q_matrix, qaoa_depth, gamma, beta)
        transpiled = transpile(circuit, backend_instance)
        start_time = time.time()
        aer_job = backend_instance.run(transpiled, shots=shots)
        result = aer_job.result()
        end_time = time.time()
        counts = result.get_counts()
        bitstring = select_top_bitstring(counts)
        weights = decode_bitstring(
            bitstring, int(problem["num_bits_per_asset"]), int(problem["num_assets"])
        )
        objective = objective_value(q_matrix, bitstring)
        result_payload = {
            "schema_version": "1.0",
            "job_id": job["job_id"],
            "problem_id": problem["problem_id"],
            "solver_name": solver_name,
            "execution_backend": job["backend"],
            "weights": weights,
            "bitstring": bitstring,
            "objective_value": objective,
            "circuit_depth": int(transpiled.depth()),
            "circuit_execution_us": -1.0,
            "shots": shots,
            "solve_time_ms": float((end_time - start_time) * 1000.0),
            "convergence_info": f"top bitstring count: {counts[bitstring]}/{shots}",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "status": "COMPLETED",
        }
        return result_payload
    except Exception as exc:
        raise RuntimeError(f"Aer collection failed: {exc}")


def collect_ibm_job(job: Dict[str, Any], timeout_min: int, poll_interval: int) -> Dict[str, Any]:
    token = os.environ.get("IBM_QUANTUM_TOKEN")
    instance = os.environ.get("IBM_QUANTUM_INSTANCE", "")
    if not token:
        raise RuntimeError("IBM_QUANTUM_TOKEN is not set. Cannot collect IBM jobs.")

    with Path(job["problem_file"]).open("r", encoding="utf-8") as f:
        problem = json.load(f)

    try:
        from qiskit_ibm_runtime import QiskitRuntimeService

        service = (
            QiskitRuntimeService(channel="ibm_quantum", token=token, instance=instance)
            if instance
            else QiskitRuntimeService(channel="ibm_quantum", token=token)
        )
        backend = service.backend(job["backend"])
        runner_job = service.job(job["job_id"])
    except Exception as exc:
        raise RuntimeError(f"IBM job lookup failed: {exc}")

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

        if status in ("DONE", "COMPLETED"):
            break
        if status in ("ERROR", "FAILED", "CANCELLED"):
            break
        time.sleep(wait)
        wait = min(wait * 1.5, 300)  # cap at 5 minutes

    if status is None:
        status = "TIMEOUT"

    if status in ("ERROR", "FAILED"):
        return {
            "schema_version": "1.0",
            "job_id": job["job_id"],
            "problem_id": problem["problem_id"],
            "solver_name": f"qaoa_p{job['qaoa_depth']}_{job['backend']}",
            "execution_backend": job["backend"],
            "weights": [],
            "bitstring": "",
            "objective_value": 0.0,
            "circuit_depth": 0,
            "circuit_execution_us": -1.0,
            "shots": job["shots"],
            "solve_time_ms": 0.0,
            "convergence_info": f"Job failed with status {status}",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "status": "FAILED",
        }
    if status == "TIMEOUT":
        return {
            "schema_version": "1.0",
            "job_id": job["job_id"],
            "problem_id": problem["problem_id"],
            "solver_name": f"qaoa_p{job['qaoa_depth']}_{job['backend']}",
            "execution_backend": job["backend"],
            "weights": [],
            "bitstring": "",
            "objective_value": 0.0,
            "circuit_depth": 0,
            "circuit_execution_us": -1.0,
            "shots": job["shots"],
            "solve_time_ms": 0.0,
            "convergence_info": "Timeout reached while polling IBM job.",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "status": "TIMEOUT",
        }

    try:
        result = runner_job.result()
        counts = result.get_counts() if hasattr(result, "get_counts") else {}
        if not counts:
            convergence_info = "IBM job completed but no counts returned."
            bitstring = ""
            weights = []
        else:
            bitstring = select_top_bitstring(counts)
            weights = decode_bitstring(
                bitstring, int(problem["num_bits_per_asset"]), int(problem["num_assets"])
            )
            convergence_info = f"top bitstring count: {counts[bitstring]}/{job['shots']}"

        execution_us = result_metadata_time_us(result)
        circuit = build_qaoa_circuit(
            int(problem["num_variables"]),
            problem["Q"],
            int(job["qaoa_depth"]),
            float(job["circuit_params"]["gamma"]),
            float(job["circuit_params"]["beta"]),
        )
        transpiled = transpile(circuit, backend)

        return {
            "schema_version": "1.0",
            "job_id": job["job_id"],
            "problem_id": problem["problem_id"],
            "solver_name": f"qaoa_p{job['qaoa_depth']}_{job['backend']}",
            "execution_backend": job["backend"],
            "weights": weights,
            "bitstring": bitstring,
            "objective_value": objective_value(problem["Q"], bitstring),
            "circuit_depth": int(transpiled.depth()),
            "circuit_execution_us": execution_us,
            "shots": job["shots"],
            "solve_time_ms": float((time.time() - start_time) * 1000.0),
            "convergence_info": convergence_info,
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "status": "COMPLETED",
        }
    except Exception as exc:
        return {
            "schema_version": "1.0",
            "job_id": job["job_id"],
            "problem_id": problem["problem_id"],
            "solver_name": f"qaoa_p{job['qaoa_depth']}_{job['backend']}",
            "execution_backend": job["backend"],
            "weights": [],
            "bitstring": "",
            "objective_value": 0.0,
            "circuit_depth": 0,
            "circuit_execution_us": -1.0,
            "shots": job["shots"],
            "solve_time_ms": 0.0,
            "convergence_info": f"IBM collection exception: {exc}",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "status": "FAILED",
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
                    service = QiskitRuntimeService(channel="ibm_quantum", token=token, instance=instance)
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
