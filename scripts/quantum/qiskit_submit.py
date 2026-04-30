#!/usr/bin/env python3
"""Submit a QAOA job for portfolio optimization using Qiskit."""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
from qiskit import QuantumCircuit, transpile


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Submit a QAOA/QAMO/QAMOO job for portfolio optimization."
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
    parser.add_argument(
        "--mode",
        default="qaoa",
        choices=["qaoa", "qamo", "qamoo"],
        help="Algorithm mode: qaoa (default), qamo, or qamoo sweep.",
    )
    parser.add_argument(
        "--n-frontier-points",
        type=int,
        default=20,
        help="Number of lambda sweep points for QAMOO mode.",
    )
    parser.add_argument(
        "--lambda-min",
        type=float,
        default=0.1,
        help="Minimum lambda value for QAMOO sweep.",
    )
    parser.add_argument(
        "--lambda-max",
        type=float,
        default=20.0,
        help="Maximum lambda value for QAMOO sweep.",
    )
    parser.add_argument(
        "--augment-problem-data",
        action="store_true",
        help="Overwrite expected_returns and covariance in the problem file with "
             "values computed from historical prices before circuit submission.",
    )
    parser.add_argument(
        "--list-backends",
        action="store_true",
        help="Authenticate with IBM, print available backends with qubit counts and queue depths, then exit 0.",
    )
    parser.add_argument(
        "--select-backend",
        action="store_true",
        help="Print the name of the operational IBM backend with the shortest queue "
             "(fewest pending jobs) and at least --min-qubits qubits, then exit 0. "
             "Designed for use with command substitution: BACKEND=$(... --select-backend).",
    )
    parser.add_argument(
        "--min-qubits",
        type=int,
        default=20,
        help="Minimum qubit count when selecting a backend with --select-backend (default: 20).",
    )
    parser.add_argument(
        "--error-mitigation",
        action="store_true",
        help="Enable mthree readout error correction when collecting IBM jobs.",
    )
    return parser.parse_args()


def validate_ibm_credentials() -> bool:
    """Validate IBM Quantum credentials and return True on success, False with warning on failure."""
    token = os.environ.get("IBM_QUANTUM_TOKEN")
    instance = os.environ.get("IBM_QUANTUM_INSTANCE", "")
    if not token:
        print("FAIL: IBM_QUANTUM_TOKEN environment variable is not set.", file=sys.stderr)
        return False

    try:
        from qiskit_ibm_runtime import QiskitRuntimeService
        service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token, instance=instance)
        # Try to list backends to verify
        backends = service.backends()
        if backends:
            print("PASS: IBM Quantum credentials are valid.", file=sys.stderr)
            return True
        else:
            print("FAIL: No backends available.", file=sys.stderr)
            return False
    except Exception as exc:
        print(f"FAIL: IBM Quantum authentication failed: {exc}", file=sys.stderr)
        return False


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


def _ewma_covariance(returns: "np.ndarray", lam: float) -> "np.ndarray":
    """EWMA covariance matching C++ EWMACovariance::estimate_covariance.

    Demeaned recursive update: Cov_t = λ·Cov_{t-1} + (1-λ)·r_t·r_t^T,
    initialised with the outer product of the first centred return row.
    Returns a daily covariance matrix (not annualised).
    """
    means = returns.mean(axis=0)
    centered = returns - means
    cov = np.outer(centered[0], centered[0])
    for t in range(1, len(centered)):
        r = centered[t]
        cov = lam * cov + (1.0 - lam) * np.outer(r, r)
    return 0.5 * (cov + cov.T)


def _augment_problem_data(problem_file: Path, problem: Dict[str, Any]) -> None:
    """Overwrite expected_returns and covariance in problem dict and file.

    Computes both fields from historical prices using the same estimators as
    the C++ classical optimizer so that QAMOO and Markowitz frontiers are
    directly comparable:
      - expected_returns: simple arithmetic daily mean × 252 (annualised)
      - covariance:       EWMA (λ=0.94) daily covariance × 252 (annualised)

    Falls back to zero expected returns and a 20%-vol diagonal covariance if
    historical prices cannot be loaded.
    """
    import pandas as pd

    repo_root = Path(__file__).resolve().parent.parent.parent
    config_path = repo_root / "data" / "config" / "portfolio_config.json"
    num_assets = int(problem.get("num_assets", 0))

    try:
        with config_path.open("r", encoding="utf-8") as f:
            cfg = json.load(f)

        ewma_lambda: float = cfg.get("risk_model", {}).get("ewma_lambda", 0.94)
        prices_path = repo_root / "data" / "market" / "historical_prices.csv"
        if not prices_path.exists():
            raise FileNotFoundError(f"Prices not found: {prices_path}")

        # Resolve asset names: C++ writes placeholder "asset_0" names, so fall
        # back to the universe from config.
        asset_names: List[str] = cfg.get("data", {}).get("universe", [])
        prices = pd.read_csv(prices_path, index_col="date", parse_dates=True)
        cols = [c for c in asset_names if c in prices.columns]
        if len(cols) != num_assets:
            raise ValueError(
                f"Expected {num_assets} assets in prices, found {len(cols)}"
            )

        daily_ret = prices[cols].pct_change().dropna().values

        expected_returns = (daily_ret.mean(axis=0) * 252).tolist()
        cov = (_ewma_covariance(daily_ret, ewma_lambda) * 252).tolist()

        print(
            f"  Augmented problem: annualised expected_returns + "
            f"EWMA covariance (λ={ewma_lambda}) for {num_assets} assets.",
            file=sys.stderr,
        )
    except Exception as exc:
        print(
            f"  Warning: could not compute market statistics ({exc}); "
            "using zero expected_returns and diagonal covariance fallback.",
            file=sys.stderr,
        )
        expected_returns = [0.0] * num_assets
        cov = (0.04 * np.eye(num_assets)).tolist()  # 20% vol diagonal

    problem["expected_returns"] = expected_returns
    problem["covariance"] = cov
    _atomic_json_write(problem_file, problem)


def _atomic_json_write(target_file: Path, obj: Any) -> None:
    target_file.parent.mkdir(parents=True, exist_ok=True)
    tmp_fd, tmp_path_str = tempfile.mkstemp(
        dir=target_file.parent,
        prefix=target_file.stem + "_",
        suffix=".tmp",
    )
    try:
        with os.fdopen(tmp_fd, "w", encoding="utf-8") as f:
            json.dump(obj, f, indent=2)
        Path(tmp_path_str).replace(target_file)
    except Exception:
        try:
            os.unlink(tmp_path_str)
        except OSError:
            pass
        raise


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
    _atomic_json_write(jobs_file, data)


def main() -> int:
    args = parse_args()

    if args.list_backends or args.select_backend:
        if not validate_ibm_credentials():
            sys.exit(1)
        from qiskit_ibm_runtime import QiskitRuntimeService
        token = os.environ.get("IBM_QUANTUM_TOKEN")
        instance = os.environ.get("IBM_QUANTUM_INSTANCE", "")
        service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token, instance=instance)
        backends = service.backends()

        if args.list_backends:
            print("Available IBM Quantum backends:")
            for backend in backends:
                status = backend.status()
                pending_jobs = getattr(status, 'pending_jobs', 'N/A')
                print(f"  {backend.name}: {backend.num_qubits} qubits, "
                      f"operational={status.operational}, pending_jobs={pending_jobs}")
            sys.exit(0)

        # --select-backend: pick the operational backend with the fewest pending jobs.
        candidates = []
        for backend in backends:
            status = backend.status()
            if not status.operational:
                continue
            if backend.num_qubits < args.min_qubits:
                continue
            pending = getattr(status, 'pending_jobs', float('inf'))
            candidates.append((pending, backend.name))

        if not candidates:
            print(
                f"Error: no operational IBM backend with >= {args.min_qubits} qubits found.",
                file=sys.stderr,
            )
            sys.exit(1)

        candidates.sort()
        # Print only the backend name — clean for shell command substitution.
        print(candidates[0][1])
        sys.exit(0)

    problem_file = Path(args.problem_file)
    jobs_file = Path(args.jobs_file)
    problem = load_problem(problem_file)

    num_variables = problem["num_variables"]
    q_matrix = problem["Q"]
    backend_name = args.backend
    shots = args.shots
    qaoa_depth = args.qaoa_depth
    mode = args.mode

    if args.augment_problem_data:
        _augment_problem_data(problem_file, problem)

    # Build initial circuit for submit-time validation / job-id generation
    if mode in ("qamo", "qamoo"):
        from qiskit_solver import build_qamo_circuit
        _mf_cfg: Dict[str, Any] = {"max_mf_iterations": 50, "mf_convergence_tol": 1e-6}
        Q_np = np.array(q_matrix)
        circuit, _ = build_qamo_circuit(Q_np, 0.1, 0.1,
                                        _mf_cfg["max_mf_iterations"],
                                        _mf_cfg["mf_convergence_tol"])
    else:
        circuit = build_qaoa_circuit(num_variables, q_matrix, qaoa_depth, 0.1, 0.1)

    pre_depth: int = 0
    circuit_depth: int = 0
    physical_qubits: List[int] = []

    if backend_name == "aer_simulator":
        try:
            from qiskit_aer import AerSimulator

            # Preflight: statevector memory grows as 2^n × 16 bytes.
            # 26 qubits = 1 GB; beyond that simulation will likely OOM.
            if num_variables > 25:
                print(
                    f"Error: Circuit requires {num_variables} qubits. Aer statevector "
                    f"simulation needs 2^{num_variables} × 16 B ≈ "
                    f"{(2 ** num_variables * 16) // (1024 ** 2):,} MB of RAM, which "
                    f"will exceed available memory on most machines. "
                    f"Reduce num_bits_per_asset (currently implied by problem size).",
                    file=sys.stderr,
                )
                return 1

            backend_instance = AerSimulator()
            pre_depth = int(circuit.depth())
            transpiled_aer = transpile(circuit, backend_instance, optimization_level=1)
            circuit_depth = int(transpiled_aer.depth())
            physical_qubits = list(range(num_variables))
            job = backend_instance.run(transpiled_aer, shots=shots)
            result = job.result()
            if not result.success:
                status = getattr(result, "status", "simulation failed")
                print(f"Error: Aer submission failed: {status}", file=sys.stderr)
                return 1
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
                QiskitRuntimeService(channel="ibm_quantum_platform", token=token, instance=instance)
                if instance
                else QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
            )
            backend = service.backend(backend_name)
            pre_depth = int(circuit.depth())
            transpiled = transpile(circuit, backend=backend, optimization_level=1)
            circuit_depth = int(transpiled.depth())
            try:
                layout = transpiled.layout
                if layout is not None and getattr(layout, "final_layout", None) is not None:
                    physical_qubits = sorted(layout.final_layout.get_physical_bits().keys())
                else:
                    physical_qubits = list(range(num_variables))
            except Exception:
                physical_qubits = list(range(num_variables))
            print(f"  Circuit depth: {pre_depth} (logical) → {circuit_depth} (post-transpile)", file=sys.stderr)
            print(f"  Physical qubits: {physical_qubits}", file=sys.stderr)
            sampler = SamplerV2(mode=backend)
            job = sampler.run([transpiled], shots=shots)
            job_id = job.job_id()
        except Exception as exc:
            print(f"Error: IBM submission failed: {exc}", file=sys.stderr)
            return 1
    else:
        print(f"Error: unsupported backend '{backend_name}'", file=sys.stderr)
        return 1

    job_entry: Dict[str, Any] = {
        "job_id": job_id,
        "backend": backend_name,
        "problem_file": str(problem_file),
        "problem_size": problem["num_assets"],
        "qaoa_depth": qaoa_depth,
        "shots": shots,
        "submitted_at": datetime.now(timezone.utc).isoformat(),
        "status": "QUEUED",
        "mode": mode,
        "circuit_params": {
            "gamma": 0.1,
            "beta": 0.1,
            "qaoa_depth": qaoa_depth,
            "error_mitigation": args.error_mitigation,
        },
        "pre_transpilation_depth": pre_depth,
        "circuit_depth": circuit_depth,
        "physical_qubits": physical_qubits,
    }
    if mode == "qamoo":
        job_entry["n_frontier_points"] = args.n_frontier_points
        job_entry["lambda_min"] = args.lambda_min
        job_entry["lambda_max"] = args.lambda_max

    append_job_entry(jobs_file, job_entry)
    print(f"Job {job_id} submitted to {backend_name}. Run --quantum-collect to retrieve results.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
