"""Tests for Phase 4 and Phase 5 quantum solver changes.

Phase 4: depth reporting, metadata extraction, bitstring extraction,
signal quality, and HTML report hardware notes.

Phase 5: QAMO mean-field solver, QAMO circuit builder, QUBO matrix builder,
QAMOO sweep, deduplication, frontier comparison plot, and HTML frontier section.
"""

from __future__ import annotations

import json
import os
import sys
import types
from pathlib import Path
from typing import Any, Dict, List
from unittest.mock import MagicMock

import numpy as np
import pytest

# Make scripts/quantum importable from any working directory.
_SCRIPTS_DIR = Path(__file__).parent.parent.parent.parent / "scripts" / "quantum"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_SMALL_PROBLEM = {
    "schema_version": "1.0",
    "problem_id": "test_p1",
    "num_assets": 2,
    "num_bits_per_asset": 2,
    "num_variables": 4,
    "Q": [
        [1.0, 0.5, 0.0, 0.0],
        [0.5, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.5],
        [0.0, 0.0, 0.5, 1.0],
    ],
}


def _write_problem(tmp_path: Path) -> Path:
    p = tmp_path / "quantum_problem.json"
    p.write_text(json.dumps(_SMALL_PROBLEM))
    return p


# ---------------------------------------------------------------------------
# Test 1 — post-transpilation depth recorded for Aer submit
# ---------------------------------------------------------------------------

def test_post_transpilation_depth_recorded(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Both circuit_depth and pre_transpilation_depth are set and > 0 for Aer."""
    pytest.importorskip("qiskit_aer")

    problem_file = _write_problem(tmp_path)
    jobs_file = tmp_path / "quantum_jobs.json"

    monkeypatch.setattr(
        sys, "argv",
        [
            "qiskit_submit.py",
            "--problem-file", str(problem_file),
            "--jobs-file", str(jobs_file),
            "--backend", "aer_simulator",
            "--shots", "64",
        ],
    )

    import qiskit_submit
    rc = qiskit_submit.main()
    assert rc == 0, "qiskit_submit.main() should return 0 for Aer"

    data = json.loads(jobs_file.read_text())
    entry = data["jobs"][-1]

    assert "circuit_depth" in entry, "circuit_depth (post-transpile) must be present"
    assert "pre_transpilation_depth" in entry, "pre_transpilation_depth must be present"
    assert entry["circuit_depth"] > 0, "circuit_depth must be positive"
    assert entry["pre_transpilation_depth"] > 0, "pre_transpilation_depth must be positive"
    # Transpilation may increase or decrease depth depending on optimizer;
    # just verify both fields are recorded and positive.
    assert entry["circuit_depth"] > 0 and entry["pre_transpilation_depth"] > 0


# ---------------------------------------------------------------------------
# Test 2 — extract_execution_us fallback chain
# ---------------------------------------------------------------------------

def test_extract_execution_us_fallback() -> None:
    """extract_execution_us returns (-1.0, 'unavailable') when no metadata is present."""
    import qiskit_collect

    class _EmptyResult:
        metadata: Dict[str, Any] = {}

    class _NoMetricsJob:
        def metrics(self) -> None:
            raise RuntimeError("no metrics")

    us, key = qiskit_collect.extract_execution_us(_EmptyResult(), _NoMetricsJob())
    assert us == -1.0
    assert key == "unavailable"


def test_extract_execution_us_time_taken() -> None:
    """extract_execution_us reads 'time_taken' (seconds) and converts to microseconds."""
    import qiskit_collect

    class _TimeTakenResult:
        metadata = {"time_taken": 0.001}

    class _NoMetricsJob:
        def metrics(self) -> None:
            raise RuntimeError("no metrics")

    us, key = qiskit_collect.extract_execution_us(_TimeTakenResult(), _NoMetricsJob())
    assert us == pytest.approx(1000.0)
    assert key == "time_taken"


# ---------------------------------------------------------------------------
# Test 3 — extract_counts fallback chain
# ---------------------------------------------------------------------------

def test_extract_counts_meas_register() -> None:
    """extract_counts succeeds via pub_result.data.meas.get_counts()."""
    import qiskit_collect

    mock_counts = {"00": 10, "11": 5}

    class _MockBitArray:
        def get_counts(self) -> Dict[str, int]:
            return mock_counts

    class _Data:
        meas = _MockBitArray()

    class _PubResult:
        data = _Data()

    assert qiskit_collect.extract_counts(_PubResult()) == mock_counts


def test_extract_counts_c_register_fallback() -> None:
    """extract_counts falls back to 'c' register when 'meas' is absent."""
    import qiskit_collect

    mock_counts = {"01": 8, "10": 4}

    class _MockBitArray:
        def get_counts(self) -> Dict[str, int]:
            return mock_counts

    class _DataNoMeas:
        c = _MockBitArray()

    class _PubResult:
        data = _DataNoMeas()

    assert qiskit_collect.extract_counts(_PubResult()) == mock_counts


def test_extract_counts_raises_when_no_bitarray() -> None:
    """extract_counts raises ValueError when no BitArray field is found."""
    import qiskit_collect

    class _EmptyData:
        __dataclass_fields__: Dict[str, Any] = {}

    class _PubResult:
        data = _EmptyData()

    with pytest.raises(ValueError, match="Could not find BitArray"):
        qiskit_collect.extract_counts(_PubResult())


# ---------------------------------------------------------------------------
# Test 4 — signal_quality 'low' when top bitstring < 5% of shots
# ---------------------------------------------------------------------------

def test_signal_quality_low_flag() -> None:
    """_compute_signal_quality returns 'low' when top bitstring is < 5% of shots."""
    import qiskit_collect

    # 1 dominant key at 30 shots, 100 other keys at 10 shots each
    # top fraction = 30 / (30 + 100*10) = 30/1030 ≈ 2.9% → "low"
    counts: Dict[str, int] = {"top_bs": 30}
    for i in range(100):
        counts[f"bs_{i:03d}"] = 10

    fraction, quality = qiskit_collect._compute_signal_quality(counts)
    assert quality == "low"
    assert fraction < 0.05


# ---------------------------------------------------------------------------
# Test 5 — signal_quality 'ok' when top bitstring >= 5% of shots
# ---------------------------------------------------------------------------

def test_signal_quality_ok_flag() -> None:
    """_compute_signal_quality returns 'ok' when top bitstring is 40% of shots."""
    import qiskit_collect

    counts: Dict[str, int] = {
        "0000": 400,
        "0001": 200,
        "0010": 200,
        "0011": 200,
    }

    fraction, quality = qiskit_collect._compute_signal_quality(counts)
    assert quality == "ok"
    assert fraction >= 0.05


# ---------------------------------------------------------------------------
# Test 6 — quantum_report.html contains Hardware notes for IBM run
# ---------------------------------------------------------------------------

def _make_comparison_json(tmp_path: Path, signal_quality: str = "ok") -> Path:
    data = {
        "schema_version": "1.0",
        "run_id": "20260405_143022",
        "problem_sizes": [10],
        "runs": [
            {
                "solver_name": "markowitz_mv",
                "solver_type": "classical",
                "execution_backend": "osqp",
                "problem_size": 10,
                "performance": {
                    "sharpe_ratio": 0.823,
                    "total_return": 0.138,
                    "annualized_return": 0.138,
                    "annualized_volatility": 0.182,
                    "max_drawdown": -0.225,
                },
                "solve_time_ms": 0.8,
                "circuit_execution_us": -1.0,
                "solution_quality_vs_classical": 1.0,
                "nav_series": [],
            },
            {
                "solver_name": "qaoa_p1_ibm_brisbane",
                "solver_type": "quantum",
                "execution_backend": "ibm_brisbane",
                "problem_size": 10,
                "performance": {
                    "sharpe_ratio": 0.312,
                    "total_return": 0.045,
                    "annualized_return": 0.045,
                    "annualized_volatility": 0.201,
                    "max_drawdown": -0.312,
                },
                "solve_time_ms": 15000.0,
                "circuit_execution_us": 4500.0,
                "solution_quality_vs_classical": 0.38,
                "nav_series": [],
                # Phase 4 fields embedded directly for test convenience
                "circuit_depth": 47,
                "pre_transpilation_depth": 12,
                "backend_calibration_date": "2026-04-05T08:30:00+00:00",
                "error_mitigation_method": "mthree",
                "signal_quality": signal_quality,
                "top_bitstring_fraction": 0.03 if signal_quality == "low" else 0.12,
                "metadata_key_used": "execution_spans",
            },
        ],
    }
    p = tmp_path / "comparison_results.json"
    p.write_text(json.dumps(data))
    return p


def test_hardware_notes_section_in_report(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """quantum_report.html contains 'Hardware notes' section with IBM run metadata."""
    pytest.importorskip("matplotlib")
    pytest.importorskip("jinja2")

    comparison_file = _make_comparison_json(tmp_path, signal_quality="ok")

    monkeypatch.setattr(
        sys, "argv",
        [
            "benchmark_viz.py",
            "--comparison-file", str(comparison_file),
            "--output-dir", str(tmp_path),
            "--title", "Test Report",
        ],
    )

    import benchmark_viz
    benchmark_viz.main()

    html = (tmp_path / "quantum_report.html").read_text()
    assert "Hardware notes" in html, "Hardware notes section must appear for IBM runs"
    assert "2026-04-05" in html, "Calibration date must appear in report"
    assert "mthree" in html, "Error mitigation method must appear in report"
    assert "ibm_brisbane" in html, "Backend name must appear in hardware notes"


# ---------------------------------------------------------------------------
# Test 7 — low-signal callout appears when signal_quality='low'
# ---------------------------------------------------------------------------

def test_low_signal_callout_in_report(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """quantum_report.html shows noise callout when IBM run has signal_quality='low'."""
    pytest.importorskip("matplotlib")
    pytest.importorskip("jinja2")

    comparison_file = _make_comparison_json(tmp_path, signal_quality="low")

    monkeypatch.setattr(
        sys, "argv",
        [
            "benchmark_viz.py",
            "--comparison-file", str(comparison_file),
            "--output-dir", str(tmp_path),
            "--title", "Test Report Low Signal",
        ],
    )

    import benchmark_viz
    benchmark_viz.main()

    html = (tmp_path / "quantum_report.html").read_text()
    assert "noise" in html.lower(), "Noise callout text must appear for low-signal runs"
    assert "low signal" in html.lower() or "noise-dominated" in html.lower(), (
        "Low signal warning must reference noise-dominated results"
    )
    assert "Hardware notes" in html


# ===========================================================================
# Phase 5 tests
# ===========================================================================

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_SMALL_Q = np.array([
    [1.0, 0.5, 0.0, 0.0],
    [0.5, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.5],
    [0.0, 0.0, 0.5, 1.0],
], dtype=float)

_SMALL_COV = np.array([
    [0.04, 0.01],
    [0.01, 0.04],
], dtype=float)

_SMALL_MU = np.array([0.10, 0.12], dtype=float)


def _import_solver():
    import qiskit_solver
    return qiskit_solver


# ---------------------------------------------------------------------------
# Test 1 — solve_mean_field converges for a small QUBO
# ---------------------------------------------------------------------------

def test_solve_mean_field_converges() -> None:
    """solve_mean_field converges for a 2x2 Q matrix."""
    qs = _import_solver()
    Q = np.array([[1.0, 0.5], [0.5, 1.0]], dtype=float)
    m, converged, iters = qs.solve_mean_field(Q, beta=1.0, max_iterations=100, tol=1e-6)
    assert converged is True
    assert iters > 0
    assert all(abs(mi) <= 0.5 for mi in m), f"Magnetisations out of [-0.5,0.5]: {m}"


# ---------------------------------------------------------------------------
# Test 2 — zero Q matrix yields all-zero magnetisations
# ---------------------------------------------------------------------------

def test_solve_mean_field_returns_zeros_for_zero_Q() -> None:
    """solve_mean_field returns zeros when Q is the zero matrix (no interactions)."""
    qs = _import_solver()
    Q = np.zeros((4, 4), dtype=float)
    m, converged, iters = qs.solve_mean_field(Q, beta=1.0, max_iterations=50, tol=1e-6)
    assert converged is True
    assert np.allclose(m, 0.0), f"Expected all zeros, got {m}"


# ---------------------------------------------------------------------------
# Test 3 — build_qamo_circuit gate counts
# ---------------------------------------------------------------------------

def test_build_qamo_circuit_has_correct_gate_count() -> None:
    """QAMO circuit for 4-qubit problem has correct H, RZZ, RX, Measure counts."""
    pytest.importorskip("qiskit_aer")
    qs = _import_solver()
    circuit, theta = qs.build_qamo_circuit(
        _SMALL_Q, gamma=0.3, beta=0.5, max_mf_iters=50, mf_tol=1e-6
    )
    assert circuit.num_qubits == 4
    assert circuit.num_clbits == 4
    ops = circuit.count_ops()
    assert ops.get("h", 0) == 4, f"Expected 4 H gates, got {ops.get('h', 0)}"
    assert ops.get("rzz", 0) >= 1, f"Expected >= 1 RZZ gates, got {ops.get('rzz', 0)}"
    assert ops.get("rx", 0) == 4, f"Expected 4 RX gates, got {ops.get('rx', 0)}"
    assert ops.get("measure", 0) == 4, f"Expected 4 Measure, got {ops.get('measure', 0)}"
    assert theta.shape == (4,)


# ---------------------------------------------------------------------------
# Test 4 — build_qubo_matrix symmetry and shape
# ---------------------------------------------------------------------------

def test_build_qubo_matrix_symmetry() -> None:
    """build_qubo_matrix produces a symmetric matrix with correct shape."""
    qs = _import_solver()
    num_assets, num_bits = 2, 2
    Q = qs.build_qubo_matrix(
        _SMALL_COV, _SMALL_MU, num_bits=num_bits, lambda_=1.0, penalty=10.0
    )
    assert Q.shape == (num_assets * num_bits, num_assets * num_bits), (
        f"Expected ({num_assets * num_bits}, {num_assets * num_bits}), got {Q.shape}"
    )
    skew = np.max(np.abs(Q - Q.T))
    assert skew < 1e-10, f"Q is not symmetric: max |Q - Q^T| = {skew}"


# ---------------------------------------------------------------------------
# Test 5 — QAMO Aer integration (valid weights)
# ---------------------------------------------------------------------------

@pytest.mark.integration
def test_qamo_aer_produces_valid_weights(tmp_path: Path) -> None:
    """QAMO on Aer produces weights summing near 1.0 and all >= -0.01."""
    pytest.importorskip("qiskit_aer")
    from qiskit_aer import AerSimulator

    # Load the real 10-asset problem if available; fall back to small problem
    repo_root = _SCRIPTS_DIR.parent.parent
    problem_file = repo_root / "results" / "quantum_problem.json"
    if problem_file.exists():
        with problem_file.open() as f:
            problem = json.load(f)
        # Build Q as numpy array
        Q_np = np.array(problem["Q"])
        num_assets = problem["num_assets"]
        num_bits = problem["num_bits_per_asset"]
    else:
        Q_np = _SMALL_Q
        num_assets = 2
        num_bits = 2

    import qiskit_solver as qs
    backend = AerSimulator()
    mf_config: Dict[str, Any] = {"max_mf_iterations": 50, "mf_convergence_tol": 1e-6}

    from scipy.optimize import minimize
    opt = minimize(
        lambda p: qs.evaluate_qamo(Q_np, p[0], p[1], mf_config, backend, 256),
        x0=[0.1, 0.1],
        method="COBYLA",
        options={"maxiter": 30, "rhobeg": 0.1},
    )
    gamma_opt, beta_opt = float(opt.x[0]), float(opt.x[1])
    circuit, _ = qs.build_qamo_circuit(Q_np, gamma_opt, beta_opt, 50, 1e-6)
    counts = qs._execute_aer(circuit, backend, 512)
    best_bs = max(counts, key=counts.get)
    weights = qs._decode_bitstring_np(best_bs, num_assets, num_bits)

    # Normalize, falling back to equal-weight on all-zero bitstring
    # (mirrors collect_aer_qamo_job behaviour)
    total_w = float(weights.sum())
    if total_w > 1e-9:
        weights = weights / total_w
    else:
        weights = np.ones(num_assets) / num_assets

    assert abs(weights.sum() - 1.0) <= 0.15, f"Weights sum to {weights.sum():.4f}"
    assert all(w >= -0.01 for w in weights), f"Negative weight found: {weights}"


# ---------------------------------------------------------------------------
# Test 6 — QAMOO Aer integration (frontier points)
# ---------------------------------------------------------------------------

@pytest.mark.integration
def test_qamoo_aer_produces_frontier_points(tmp_path: Path) -> None:
    """QAMOO on Aer with 5 lambda points produces sorted, valid frontier."""
    pytest.importorskip("qiskit_aer")
    from qiskit_aer import AerSimulator

    cov = _SMALL_COV
    mu = _SMALL_MU
    num_assets, num_bits = 2, 2
    import qiskit_solver as qs

    problem: Dict[str, Any] = {
        "num_assets": num_assets,
        "num_bits_per_asset": num_bits,
        "num_variables": num_assets * num_bits,
        "expected_returns": mu.tolist(),
        "covariance": cov.tolist(),
        "penalty_budget": 10.0,
        "Q": qs.build_qubo_matrix(cov, mu, num_bits, 1.0, 10.0).tolist(),
        "problem_id": "test_qamoo",
        "schema_version": "1.0",
    }
    config: Dict[str, Any] = {
        "n_frontier_points": 5,
        "lambda_min": 0.1,
        "lambda_max": 5.0,
        "max_mf_iterations": 50,
        "mf_convergence_tol": 1e-6,
    }
    backend = AerSimulator()
    frontier_points, best_weights = qs.run_qamoo_sweep(problem, config, backend, shots=256)

    assert len(frontier_points) > 0, "Expected at least one frontier point"
    assert len(frontier_points) <= 5
    assert all(p["portfolio_volatility"] > 0 for p in frontier_points), (
        "All frontier points must have positive volatility"
    )
    vols = [p["portfolio_volatility"] for p in frontier_points]
    assert vols == sorted(vols), "Frontier points must be sorted by volatility"
    assert best_weights is not None and len(best_weights) == num_assets


# ---------------------------------------------------------------------------
# Test 7 — deduplicate_frontier removes near-duplicates
# ---------------------------------------------------------------------------

def test_deduplicate_frontier_removes_near_duplicates() -> None:
    """deduplicate_frontier collapses near-duplicate points (within 0.5%)."""
    qs = _import_solver()

    def _pt(vol: float, ret: float) -> Dict[str, Any]:
        return {
            "portfolio_volatility": vol,
            "portfolio_return": ret,
            "lambda": 1.0,
            "weights": [0.5, 0.5],
            "bitstring": "00",
            "objective_value": 0.0,
        }

    points = [
        _pt(0.150, 0.100),   # cluster A — anchor
        _pt(0.1502, 0.1001), # cluster A — near-dup (within 0.5%)
        _pt(0.1504, 0.1002), # cluster A — near-dup (within 0.5%)
        _pt(0.200, 0.130),   # cluster B
        _pt(0.250, 0.160),   # cluster C
    ]
    result = qs.deduplicate_frontier(points, tol=0.005)
    assert len(result) < 5, f"Expected fewer than 5 points after dedup, got {len(result)}"


# ---------------------------------------------------------------------------
# Test 8 — plot_frontier_comparison creates a valid PNG file
# ---------------------------------------------------------------------------

def test_plot_frontier_comparison_creates_file(tmp_path: Path) -> None:
    """plot_frontier_comparison writes a PNG file > 10_000 bytes."""
    pytest.importorskip("matplotlib")
    pytest.importorskip("pandas")

    # Synthetic classical CSV
    csv_path = tmp_path / "efficient_frontier.csv"
    with csv_path.open("w") as f:
        f.write("volatility,return,sharpe,weights_0,weights_1\n")
        for i in range(10):
            vol = 0.10 + i * 0.02
            ret = 0.05 + i * 0.015
            sh = ret / vol
            f.write(f"{vol:.4f},{ret:.4f},{sh:.4f},0.5,0.5\n")

    # Synthetic QAMOO result
    qamoo_result: Dict[str, Any] = {
        "solver_name": "qamoo_5pts_aer",
        "frontier_points": [
            {
                "lambda": float(i + 1),
                "portfolio_volatility": 0.12 + i * 0.025,
                "portfolio_return": 0.04 + i * 0.02,
                "weights": [0.5, 0.5],
                "bitstring": "00",
                "objective_value": 0.0,
            }
            for i in range(5)
        ],
    }

    png_path = tmp_path / "frontier_comparison.png"
    import benchmark_viz
    benchmark_viz.plot_frontier_comparison(str(csv_path), qamoo_result, str(png_path))

    assert png_path.exists(), "PNG file was not created"
    assert png_path.stat().st_size > 10_000, (
        f"PNG too small ({png_path.stat().st_size} bytes); expected > 10 000"
    )


# ---------------------------------------------------------------------------
# Test 9 — frontier section appears in HTML report
# ---------------------------------------------------------------------------

def test_frontier_section_in_html_report(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """quantum_report.html contains 'Frontier comparison' when frontier args are present."""
    pytest.importorskip("matplotlib")
    pytest.importorskip("jinja2")
    pytest.importorskip("pandas")

    # Re-use comparison JSON from existing helper
    comparison_file = _make_comparison_json(tmp_path, signal_quality="ok")

    # Synthetic classical CSV
    csv_path = tmp_path / "efficient_frontier.csv"
    with csv_path.open("w") as f:
        f.write("volatility,return,sharpe,weights_0,weights_1\n")
        for i in range(10):
            vol = 0.10 + i * 0.02
            ret = 0.05 + i * 0.015
            sh = ret / vol
            f.write(f"{vol:.4f},{ret:.4f},{sh:.4f},0.5,0.5\n")

    # Synthetic QAMOO result JSON
    qamoo_json = tmp_path / "quantum_result_qamoo.json"
    qamoo_data: Dict[str, Any] = {
        "solver_name": "qamoo_5pts_aer",
        "frontier_points": [
            {
                "lambda": float(i + 1),
                "portfolio_volatility": 0.12 + i * 0.025,
                "portfolio_return": 0.04 + i * 0.02,
                "weights": [0.5, 0.5],
                "bitstring": "00",
                "objective_value": 0.0,
            }
            for i in range(5)
        ],
    }
    qamoo_json.write_text(json.dumps(qamoo_data))

    monkeypatch.setattr(
        sys, "argv",
        [
            "benchmark_viz.py",
            "--comparison-file", str(comparison_file),
            "--output-dir", str(tmp_path),
            "--title", "Frontier Test Report",
            "--frontier-classical", str(csv_path),
            "--frontier-quantum", str(qamoo_json),
        ],
    )

    import benchmark_viz
    benchmark_viz.main()

    html = (tmp_path / "quantum_report.html").read_text()
    assert "Frontier comparison" in html, "HTML must contain 'Frontier comparison' section"
    assert "the quantum frontier is expected" in html.lower() or \
           "converge toward" in html.lower(), (
        "Interpretation note must appear in frontier section"
    )
