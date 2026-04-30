#!/usr/bin/env python3
"""Aggregate statistics across multiple benchmark runs.

Discovers comparison_results_<timestamp>.json files in the results directory,
normalises solver names, groups metrics across runs, computes mean/std, flags
high-variance solvers, and writes aggregated_benchmark.json and
aggregated_benchmark.md.

Single source of truth: comparison_results_[0-9]*.json only.
quantum_result_*.json, quantum_jobs_*.json, and quantum_problem_*.json are
never read here.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import statistics
from pathlib import Path
from typing import Any, Dict, List, Optional


# ---------------------------------------------------------------------------
# Solver name normalisation — maps legacy internal names to display names so
# older and newer comparison_results files consolidate into the same rows.
# ---------------------------------------------------------------------------
SOLVER_NAME_MAP: Dict[str, str] = {
    "markowitz_mv":            "Markowitz",
    "sa_classical":            "QUBO (Sim. Annealing)",
    "qaoa_p1_aer_simulator":   "QAOA Uninformed (Aer)",
    "qamo_p1_aer_simulator":   "QAMO Uninformed (Aer)",
    "qaoa_p1_ibm_fez":         "QAOA (IBM fez)",
    "qamo_p1_ibm_fez":         "QAMO (IBM fez)",
    "qamoo_20pts_ibm_fez":     "QAMOO (IBM fez)",
    "qamoo_5pts_ibm_fez":      "QAMOO 5pt (IBM fez)",
    "qamoo_20pts_aer_simulator": "QAMOO (Aer Simulator)",
    "qaoa_p1_ibm_kingston":        "QAOA (IBM kingston)",
    "qamo_p1_ibm_kingston":        "QAMO (IBM kingston)",
    "qamoo_20pts_ibm_kingston":    "QAMOO (IBM kingston)",
    "qamoo_5pts_ibm_kingston":     "QAMOO 5pt (IBM kingston)",
    "qaoa_p1_ibm_marrakesh":       "QAOA (IBM marrakesh)",
    "qamo_p1_ibm_marrakesh":       "QAMO (IBM marrakesh)",
    "qamoo_20pts_ibm_marrakesh":   "QAMOO (IBM marrakesh)",
    "qamoo_5pts_ibm_marrakesh":    "QAMOO 5pt (IBM marrakesh)",
}


def normalize_name(name: str) -> str:
    return SOLVER_NAME_MAP.get(name, name)


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------

def _mean(values: List[float]) -> float:
    return statistics.mean(values) if values else 0.0


def _stdev(values: List[float]) -> float:
    return statistics.stdev(values) if len(values) >= 2 else 0.0


def _is_high_variance(mean: float, std: float) -> bool:
    """Return True when std / |mean| > 0.5."""
    return abs(mean) > 1e-9 and (std / abs(mean)) > 0.5


def _fmt(value: float, precision: int = 3) -> str:
    return f"{value:.{precision}f}"


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

def discover_files(results_dir: str | Path) -> List[Path]:
    """Return timestamped comparison_results_[0-9]*.json files, sorted."""
    pattern = os.path.join(str(results_dir), "comparison_results_[0-9]*.json")
    files = sorted(
        Path(f) for f in glob.glob(pattern)
        if not f.endswith("comparison_results_latest.json")
    )
    return files


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _get_float(d: dict, key: str) -> Optional[float]:
    v = d.get(key)
    try:
        return float(v) if v is not None else None
    except (TypeError, ValueError):
        return None


def load_runs(files: List[Path]) -> Dict[str, List[Dict[str, Any]]]:
    """Load and group solver runs by normalised solver_name.

    Fields are read from the actual schema produced by main.cpp:
      - Sharpe, return, volatility, drawdown: run["performance"][field]
      - Turnover: run["trade_summary"]["turnover"]
    """
    by_solver: Dict[str, List[Dict[str, Any]]] = {}

    for path in files:
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:
            print(f"Warning: could not read {path.name}: {exc}")
            continue

        for run in data.get("runs", []):
            raw_name = run.get("solver_name", "")
            if not raw_name:
                continue
            name = normalize_name(raw_name)

            perf = run.get("performance", {}) or {}
            ts   = run.get("trade_summary", {}) or {}

            entry: Dict[str, Any] = {
                "sharpe_ratio":          _get_float(perf, "sharpe_ratio"),
                "total_return":          _get_float(perf, "total_return"),
                "annualized_volatility": _get_float(perf, "annualized_volatility"),
                "max_drawdown":          _get_float(perf, "max_drawdown"),
                "turnover":              _get_float(ts,   "turnover"),
                "solve_time_ms":         _get_float(run,  "solve_time_ms"),
                "solver_type":           run.get("solver_type", ""),
                "execution_backend":     run.get("execution_backend", ""),
            }
            by_solver.setdefault(name, []).append(entry)

    return by_solver


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def aggregate_solver(entries: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Compute aggregate stats for one solver across multiple run entries."""
    def _collect(key: str) -> List[float]:
        return [e[key] for e in entries if e.get(key) is not None]

    def _stat(values: List[float]) -> Dict[str, Any]:
        mean = _mean(values)
        std  = _stdev(values)
        return {"mean": mean, "std": std, "n": len(values),
                "high_variance": _is_high_variance(mean, std)}

    return {
        "num_runs":          len(entries),
        "solver_type":       entries[0].get("solver_type", "") if entries else "",
        "execution_backend": entries[0].get("execution_backend", "") if entries else "",
        "sharpe":            _stat(_collect("sharpe_ratio")),
        "total_return":      _stat(_collect("total_return")),
        "volatility":        _stat(_collect("annualized_volatility")),
        "max_drawdown":      _stat(_collect("max_drawdown")),
        "turnover":          _stat(_collect("turnover")),
        "solve_time_ms":     _stat(_collect("solve_time_ms")),
    }


# ---------------------------------------------------------------------------
# Markdown output
# ---------------------------------------------------------------------------

def build_md_table(stats: Dict[str, Dict[str, Any]]) -> str:
    header = (
        "| Solver | Runs | Mean Sharpe | ±σ | Mean Return | ±σ | "
        "Mean Turnover | ±σ | Variance |"
    )
    sep = (
        "|--------|------|-------------|-----|-------------|-----|"
        "--------------|-----|----------|"
    )
    rows = [header, sep]

    for solver_name, agg in stats.items():
        sharpe   = agg["sharpe"]
        ret      = agg["total_return"]
        turnover = agg["turnover"]

        high_var = (
            sharpe["high_variance"]
            or ret["high_variance"]
            or turnover["high_variance"]
        )
        variance_flag = "⚠ high" if high_var else "ok"
        turn_mean = _fmt(turnover["mean"]) if turnover["n"] > 0 else "n/a"
        turn_std  = _fmt(turnover["std"])  if turnover["n"] > 0 else "n/a"

        rows.append(
            f"| {solver_name} "
            f"| {agg['num_runs']} "
            f"| {_fmt(sharpe['mean'])} "
            f"| {_fmt(sharpe['std'])} "
            f"| {_fmt(ret['mean'], 4)} "
            f"| {_fmt(ret['std'], 4)} "
            f"| {turn_mean} "
            f"| {turn_std} "
            f"| {variance_flag} |"
        )

    return "\n".join(rows)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Aggregate statistics across timestamped benchmark runs."
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results"),
        help="Directory containing comparison_results_[0-9]*.json files (default: results/).",
    )
    args = parser.parse_args()

    results_dir: Path = args.results_dir
    if not results_dir.exists():
        raise FileNotFoundError(f"Results directory not found: {results_dir}")

    files = discover_files(results_dir)
    if not files:
        print(
            f"No timestamped comparison_results_[0-9]*.json files found in {results_dir}. "
            "Run --benchmark first to generate archived results."
        )
        return

    print(f"Found {len(files)} run file(s):")
    for f in files:
        print(f"  {f.name}")

    by_solver = load_runs(files)
    if not by_solver:
        print("No solver runs found across discovered files.")
        return

    stats: Dict[str, Dict[str, Any]] = {
        name: aggregate_solver(entries)
        for name, entries in sorted(by_solver.items())
    }

    # Write JSON
    json_path = results_dir / "aggregated_benchmark.json"
    json_path.write_text(
        json.dumps({"num_files": len(files), "solvers": stats}, indent=2),
        encoding="utf-8",
    )
    print(f"\nAggregated stats written to: {json_path}")

    # Write and print Markdown
    md_table = build_md_table(stats)
    md_path = results_dir / "aggregated_benchmark.md"
    md_path.write_text(md_table + "\n", encoding="utf-8")
    print(f"Markdown table written to: {md_path}\n")
    print(md_table)


if __name__ == "__main__":
    main()
