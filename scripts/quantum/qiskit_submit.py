#!/usr/bin/env python3
"""Submit a QAOA job for portfolio optimization using Qiskit."""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

from qiskit import QuantumCircuit, transpile


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Submit a QAOA job for portfolio optimization."
    )
    parser.add_argument(
        "--problem-file",
        default="results/quantum_problem.json",
        help="Path to the quantum problem JSON file.",
    )
    parser.add_argument(
        "--jobs-file",
        default="results/quantum_jobs.json",
        help="Path to the quantum jobs JSON file.",
    )
    parser.add_argument(
        "--backend",
        default="aer_simulator",
        help="Qiskit backend name (default: aer_simulator).",
    )
    parser.add_argument(
        "--shots",
        type=int,
        default=1024,
        help="Number of measurement shots.",
    )
    parser.add_argument(
        "--qaoa-depth",
        type=int,
        default=1,
        help="QAOA circuit depth p.",
    )
    return parser.parse_args()


def load_problem(problem_file: Path) -> Dict[str, Any]:
    if not problem_file.exists():
        print(f"Error: problem file not found: {problem_file}", file=sys.stderr)
        sys.exit(1)

    with problem_file.open("r", encoding="utf-8") as f:
        data = json.load(f)

    if data.get("schema_version") != "1.0":
        print(
            f"Error: unsupported schema_version in {problem_file}; expected 1.0.",
            file=sys.stderr,
        )
        sys.exit(1)

    required = ["problem_id", "num_assets", "num_bits_per_asset", "num_variables", "Q"]
    for key in required:
        if key not in data:
            print(f"Error: missing required field '{key}' in {problem_file}", file=sys.stderr)
            sys.exit(1)

    return data


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


def safe_job_id(metadata: Dict[str, Any], backend_name: str) -> str:
    if isinstance(metadata, dict):
        job_id = metadata.get("job_id") or metadata.get("id")
        if isinstance(job_id, str) and job_id:
            return job_id
    timestamp = int(datetime.now(timezone.utc).timestamp() * 1_000_000)
    return f"{backend_name}_{timestamp}"


def append_job_entry(jobs_file: Path, job_entry: Dict[str, Any]) -> None:
    if jobs_file.exists():
        with jobs_file.open("r", encoding="utf-8") as f:
            data = json.load(f)
        if not isinstance(data, dict) or "jobs" not in data:
            print(f"Error: invalid jobs file format: {jobs_file}", file=sys.stderr)
            sys.exit(1)
    else:
        data = {
            "schema_version": "1.0",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "jobs": [],
        }

    data["jobs"].append(job_entry)
    jobs_file.parent.mkdir(parents=True, exist_ok=True)
    with jobs_file.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def main() -> int:
    args = parse_args()
    problem_file = Path(args.problem_file)
    jobs_file = Path(args.jobs_file)
    problem = load_problem(problem_file)

    num_variables = problem["num_variables"]
    q_matrix = problem["Q"]
    backend_name = args.backend
    shots = args.shots
    qaoa_depth = args.qaoa_depth

    circuit = build_qaoa_circuit(num_variables, q_matrix, qaoa_depth, 0.1, 0.1)

    if backend_name == "aer_simulator":
        try:
            from qiskit_aer import AerSimulator

            backend_instance = AerSimulator()
            transpiled = transpile(circuit, backend_instance)
            job = backend_instance.run(transpiled, shots=shots)
            result = job.result()
            metadata = getattr(result, "metadata", {}) or {}
            job_id = safe_job_id(metadata, backend_name)
        except Exception as exc:
            print(f"Error: Aer submission failed: {exc}", file=sys.stderr)
            return 1
    elif backend_name.startswith("ibm_"):
        token = os.environ.get("IBM_QUANTUM_TOKEN")
        instance = os.environ.get("IBM_QUANTUM_INSTANCE", "")
        if not token:
            print(
                "Error: IBM_QUANTUM_TOKEN is not set. Cannot submit to IBM backends.",
                file=sys.stderr,
            )
            return 1

        try:
            from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2

            service = (
                QiskitRuntimeService(token=token, instance=instance)
                if instance
                else QiskitRuntimeService(token=token)
            )
            backend_instance = service.backend(backend_name)
            transpiled = transpile(circuit, backend_instance)
            sampler = SamplerV2(backend_instance)
            job = sampler.run([transpiled], shots=shots)
            job_id = job.job_id()
        except Exception as exc:
            print(f"Error: IBM submission failed: {exc}", file=sys.stderr)
            return 1
    else:
        print(f"Error: unsupported backend '{backend_name}'", file=sys.stderr)
        return 1

    job_entry = {
        "job_id": job_id,
        "backend": backend_name,
        "problem_file": str(problem_file),
        "problem_size": problem["num_assets"],
        "qaoa_depth": qaoa_depth,
        "shots": shots,
        "submitted_at": datetime.now(timezone.utc).isoformat(),
        "status": "QUEUED",
        "circuit_params": {
            "gamma": 0.1,
            "beta": 0.1,
            "qaoa_depth": qaoa_depth,
        },
    }

    append_job_entry(jobs_file, job_entry)
    print(f"Job {job_id} submitted to {backend_name}. Run --quantum-collect to retrieve results.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
