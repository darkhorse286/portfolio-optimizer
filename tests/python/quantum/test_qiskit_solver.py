"""Tests for Phase 4 quantum solver changes: depth reporting, metadata extraction,
bitstring extraction, signal quality, and HTML report hardware notes."""

from __future__ import annotations

import json
import sys
import types
from pathlib import Path
from typing import Any, Dict
from unittest.mock import MagicMock

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
