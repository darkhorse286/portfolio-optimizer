#!/usr/bin/env python3
"""Benchmark visualization script for quantum portfolio optimizers."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import base64
import csv

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from jinja2 import Template

SOLVER_TYPE_MAP: dict[str, str] = {
    "classical": "Classical",
    "quantum_inspired": "QUBO",
    "quantum": "Quantum",
}

# Complete color palette — one entry per (algorithm-family × backend) combination.
# Algorithm families: markowitz, sa_classical, qaoa, qamo, qamoo
# Backends: aer (simulator) or ibm (real hardware)
# Classical solvers have no backend distinction so use a single key.
#
# Hue families:
#   Markowitz  — blue
#   SA/QUBO    — amber
#   QAOA       — violet / rose
#   QAMO       — teal
#   QAMOO      — pink / fuchsia
SOLVER_COLOR_MAP: dict[str, str] = {
    # Classical
    "markowitz":    "#2563eb",  # blue
    # QUBO / simulated-annealing
    "sa_classical": "#d97706",  # amber
    # QAOA
    "qaoa_aer":     "#7c3aed",  # violet
    "qaoa_ibm":     "#be123c",  # rose
    # QAMO
    "qamo_aer":     "#0d9488",  # teal
    "qamo_ibm":     "#0f766e",  # dark teal
    # QAMOO
    "qamoo_aer":    "#db2777",  # pink
    "qamoo_ibm":    "#9d174d",  # dark pink
}


def _display_name(solver_name: str) -> str:
    """Return a human-readable label for a solver name."""
    if solver_name.startswith("markowitz"):
        return "Markowitz"
    if solver_name.startswith("sa_classical"):
        return "QUBO (Sim. Annealing)"
    if solver_name.startswith("qaoa"):
        if "ibm_" in solver_name:
            backend = solver_name.split("ibm_", 1)[1]
            return f"QAOA (IBM {backend})"
        if "aer" in solver_name:
            return "QAOA (Aer Simulator)"
        return "QAOA"
    if solver_name.startswith("qamoo"):
        if "ibm_" in solver_name:
            backend = solver_name.split("ibm_", 1)[1]
            return f"QAMOO (IBM {backend})"
        return "QAMOO (Aer Simulator)"
    if solver_name.startswith("qamo"):
        if "ibm_" in solver_name:
            backend = solver_name.split("ibm_", 1)[1]
            return f"QAMO (IBM {backend})"
        return "QAMO (Aer Simulator)"
    return solver_name


def _solver_color(run: Dict[str, Any]) -> str:
    """Return the plot colour for a run.

    Dispatches on solver_name prefix first so every known algorithm family
    gets a unique, stable colour regardless of what solver_type is set to.
    Falls back to solver_type only for algorithms not yet in SOLVER_COLOR_MAP.
    """
    solver_name = run.get("solver_name", "")
    backend = run.get("execution_backend", "")
    suffix = "ibm" if backend.startswith("ibm_") else "aer"

    # Check longest prefixes first (qamoo before qamo)
    if solver_name.startswith("qamoo"):
        return SOLVER_COLOR_MAP[f"qamoo_{suffix}"]
    if solver_name.startswith("qamo"):
        return SOLVER_COLOR_MAP[f"qamo_{suffix}"]
    if solver_name.startswith("qaoa"):
        return SOLVER_COLOR_MAP[f"qaoa_{suffix}"]
    if solver_name.startswith("markowitz"):
        return SOLVER_COLOR_MAP["markowitz"]
    if solver_name.startswith("sa_classical"):
        return SOLVER_COLOR_MAP["sa_classical"]

    # Unknown algorithm — fall back to solver_type for forward-compatibility
    solver_type = run.get("solver_type", "")
    return SOLVER_COLOR_MAP.get(solver_type, "#6b7280")


def _display_type(solver_type: str) -> str:
    """Return a human-readable label for a solver type."""
    return SOLVER_TYPE_MAP.get(solver_type, solver_type)


def _compute_static_portfolio_metrics(
    weights: List[float],
    universe: List[str],
    prices_csv: Path,
    risk_free_rate: float = 0.02,
) -> Optional[Tuple[Dict[str, Any], List[Dict[str, Any]]]]:
    """Compute backtest metrics for a fixed-weight portfolio using historical prices.

    Returns (performance_dict, nav_series_list) or None on failure.
    """
    try:
        prices = pd.read_csv(prices_csv, index_col="date", parse_dates=True)
    except Exception:
        return None

    cols = [c for c in universe if c in prices.columns]
    if not cols:
        return None

    w = {c: weights[universe.index(c)] for c in cols}
    total_w = sum(w.values())
    if total_w <= 0:
        return None
    w = {c: v / total_w for c, v in w.items()}

    p = prices[cols].dropna()
    if len(p) < 2:
        return None

    daily_ret = p.pct_change().dropna()
    port_ret = sum(daily_ret[c] * w[c] for c in cols)

    nav = (1 + port_ret).cumprod()
    total_return = float(nav.iloc[-1] - 1.0)
    n_days = len(port_ret)
    ann_return = float((1 + total_return) ** (252.0 / n_days) - 1)
    ann_vol = float(port_ret.std() * math.sqrt(252))
    sharpe = (ann_return - risk_free_rate) / ann_vol if ann_vol > 0 else 0.0

    peak = nav.cummax()
    drawdown = (nav - peak) / peak
    max_drawdown = float(drawdown.min())

    nav_series = [
        {"date": str(d.date()), "nav": float(v)}
        for d, v in zip(nav.index, nav.values)
    ]

    performance = {
        "sharpe_ratio": round(sharpe, 4),
        "total_return": round(total_return, 4),
        "annualized_return": round(ann_return, 4),
        "annualized_volatility": round(ann_vol, 4),
        "max_drawdown": round(max_drawdown, 4),
    }
    return performance, nav_series


def _try_augment_runs_from_result_files(
    runs: List[Dict[str, Any]], results_dir: Path
) -> None:
    """Merge Phase 4 fields from quantum_result_*.json files into runs in-place.

    For each run, looks for a matching quantum_result file by solver_name and
    merges any Phase 4 fields not already present in the run dict.

    Args:
        runs: List of run dicts from comparison_results.json (mutated in-place).
        results_dir: Directory to scan for quantum_result_*.json files.
    """
    result_files = list(results_dir.glob("quantum_result_*.json"))
    existing_names = {run.get("solver_name") for run in runs}

    for result_file in result_files:
        try:
            with result_file.open("r", encoding="utf-8") as f:
                qr = json.load(f)
        except Exception:
            continue

        solver_name = qr.get("solver_name", "")
        backend = qr.get("execution_backend", "")

        # Case 1: merge Phase 4 metadata into an existing matching run
        for run in runs:
            if run.get("solver_name") == solver_name:
                for key in (
                    "circuit_depth",
                    "pre_transpilation_depth",
                    "backend_calibration_date",
                    "error_mitigation_method",
                    "signal_quality",
                    "top_bitstring_fraction",
                    "metadata_key_used",
                    "physical_qubits",
                ):
                    if key in qr and key not in run:
                        run[key] = qr[key]

        # Case 2: result not yet in comparison_results.json — inject it.
        # Covers IBM hardware runs and any Aer QAMO/QAMOO runs added after
        # the C++ comparison was last generated.
        # Attempt a static-weight backtest so real metrics are shown.
        if solver_name not in existing_names and qr.get("status") == "COMPLETED":
            weights: List[float] = qr.get("weights", [])
            universe: List[str] = qr.get("universe", [])

            # Resolve market data relative to this script's repo root
            repo_root = Path(__file__).resolve().parent.parent.parent
            config_path = repo_root / "data" / "config" / "portfolio_config.json"
            prices_csv = repo_root / "data" / "market" / "historical_prices.csv"
            risk_free_rate = 0.02
            if config_path.exists():
                try:
                    cfg = json.loads(config_path.read_text())
                    risk_free_rate = cfg.get("performance", {}).get("risk_free_rate", 0.02)
                    if not universe:
                        universe = cfg.get("data", {}).get("universe", [])
                except Exception:
                    pass

            static_result = None
            if weights and universe and prices_csv.exists():
                static_result = _compute_static_portfolio_metrics(
                    weights, universe, prices_csv, risk_free_rate
                )

            if static_result is not None:
                perf, nav_series = static_result
                has_full_backtest = True
            else:
                perf = {
                    "sharpe_ratio": 0.0,
                    "total_return": 0.0,
                    "annualized_return": 0.0,
                    "annualized_volatility": 0.0,
                    "max_drawdown": 0.0,
                }
                nav_series = []
                has_full_backtest = False

            injected: Dict[str, Any] = {
                "solver_name": solver_name,
                "solver_type": "quantum",
                "execution_backend": backend,
                "problem_size": qr.get("problem_size", 0),
                "performance": perf,
                "solve_time_ms": qr.get("solve_time_ms", 0.0),
                "circuit_execution_us": qr.get("circuit_execution_us", -1.0),
                "solution_quality_vs_classical": 0.0,
                "nav_series": nav_series,
                # Phase 4 hardware metadata
                "circuit_depth": qr.get("circuit_depth", 0),
                "pre_transpilation_depth": qr.get("pre_transpilation_depth", 0),
                "backend_calibration_date": qr.get("backend_calibration_date", "unavailable"),
                "error_mitigation_method": qr.get("error_mitigation_method", "none"),
                "signal_quality": qr.get("signal_quality", "low"),
                "top_bitstring_fraction": qr.get("top_bitstring_fraction", 0.0),
                "metadata_key_used": qr.get("metadata_key_used", "unavailable"),
                "_hardware_only": not has_full_backtest,
            }
            runs.append(injected)
            existing_names.add(solver_name)


def _fill_solution_quality(runs: List[Dict[str, Any]]) -> None:
    """Compute solution_quality_vs_classical for any run where it is 0.0.

    Uses the best classical Sharpe ratio as the baseline (same denominator the
    C++ comparison layer uses).  Mutates runs in-place; safe to call multiple
    times.  Leaves the value as 0.0 only when there is genuinely no usable
    classical baseline or no backtest data for the run.
    """
    classical_sharpe: Optional[float] = None
    for run in runs:
        if run.get("solver_type") == "classical" and not run.get("_hardware_only"):
            sharpe = run.get("performance", {}).get("sharpe_ratio", 0.0)
            if classical_sharpe is None or sharpe > classical_sharpe:
                classical_sharpe = sharpe

    if classical_sharpe is None or abs(classical_sharpe) < 1e-9:
        return  # no usable classical baseline — leave values unchanged

    for run in runs:
        if run.get("solution_quality_vs_classical", 0.0) != 0.0:
            continue  # already set (including the 1.0 on the Markowitz run itself)
        if run.get("_hardware_only"):
            continue  # no backtest data — can't compute a meaningful ratio
        sharpe = run.get("performance", {}).get("sharpe_ratio")
        if sharpe is not None:
            run["solution_quality_vs_classical"] = sharpe / classical_sharpe


def plot_comparison_equity_curves(runs: List[Dict[str, Any]], output_path: str) -> None:
    """Plot normalized equity curves for all solvers.

    Args:
        runs: List of run dictionaries from comparison_results.json
        output_path: Path to save the PNG file
    """
    plt.figure(figsize=(12, 5))

    for run in runs:
        solver_name = run["solver_name"]
        if run.get("_hardware_only"):
            continue  # no backtest NAV — hardware notes only
        nav_series = run.get("nav_series", [])
        if not nav_series:
            continue

        # nav_series may be a list of floats or a list of {"nav", "date"} dicts
        if nav_series and isinstance(nav_series[0], dict):
            nav_values = [point["nav"] for point in nav_series]
            x_axis = [point.get("date", i) for i, point in enumerate(nav_series)]
        else:
            nav_values = [float(v) for v in nav_series]
            x_axis = list(range(len(nav_values)))

        if nav_values:
            first_nav = nav_values[0]
            normalized = [v / first_nav for v in nav_values]
            plt.plot(x_axis, normalized, label=_display_name(solver_name), color=_solver_color(run))

    plt.xlabel("Date")
    plt.ylabel("Portfolio Value (Normalized to 1.0)")
    plt.title("Solver Comparison — Equity Curves")
    plt.legend()
    plt.grid(True, alpha=0.3)
    ax = plt.gca()
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=10, integer=True))
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_scaling_table(results: List[Dict[str, Any]], output_path: str) -> None:
    """Plot scaling table as matplotlib table.

    Args:
        results: List of result dictionaries
        output_path: Path to save the PNG file
    """
    n_rows = len(results)
    fig_height = max(2.0, 0.55 * (n_rows + 1) + 0.5)
    fig, ax = plt.subplots(figsize=(16, fig_height))
    ax.axis('off')

    data = []
    for res in results:
        backend = res.get("execution_backend", "")
        is_ibm = backend.startswith("ibm_")
        is_hw_only = res.get("_hardware_only", False)
        circuit_depth_str = str(res["circuit_depth"]) if is_ibm and res.get("circuit_depth") else "-"
        quality = res.get("solution_quality_vs_classical", 0.0)
        if is_hw_only:
            vs_markowitz_str = "hw only"
        elif quality == 0.0:
            vs_markowitz_str = "N/A"
        else:
            vs_markowitz_str = f"{quality:.2f}x"
        row = [
            _display_name(res["solver_name"]),
            _display_type(res["solver_type"]),
            backend,
            "hw only" if is_hw_only else f"{res['performance']['sharpe_ratio']:.3f}",
            "hw only" if is_hw_only else f"{res['performance']['total_return'] * 100:.1f}%",
            vs_markowitz_str,
            f"{res['solve_time_ms']:.1f}ms",
            circuit_depth_str,
        ]
        data.append(row)

    columns = ["Solver", "Type", "Backend", "Sharpe", "Total Return", "vs Markowitz", "Solve Time (ms)", "Circuit Depth"]

    table = ax.table(cellText=data, colLabels=columns, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.auto_set_column_width(list(range(len(columns))))
    for (row, _col), cell in table.get_celld().items():
        cell.set_height(0.12 if row == 0 else 0.10)
        cell.PAD = 0.05

    # Color vs Markowitz column (index 5)
    for i, res in enumerate(results):
        quality = res.get("solution_quality_vs_classical", 0.0)
        cell = table[(i + 1, 5)]
        if quality == 0.0:
            pass  # N/A or hw only — leave default white
        elif quality > 0.95:
            cell.set_facecolor("#dcfce7")
        elif 0.70 <= quality <= 0.95:
            cell.set_facecolor("#fef3c7")
        else:
            cell.set_facecolor("#fee2e2")

    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_frontier_comparison(
    classical_csv: str,
    qamoo_result: Dict[str, Any],
    output_path: str,
) -> None:
    """Plot quantum vs classical Pareto frontier comparison.

    Args:
        classical_csv: Path to results/efficient_frontier.csv
            (columns: volatility, return, sharpe, weights_*).
        qamoo_result: Loaded from a QAMOO quantum_result JSON; must contain
            'frontier_points' list.
        output_path: Where to write the PNG.
    """
    fig, ax = plt.subplots(figsize=(10, 7))

    def _pct(ax_obj, axis: str) -> None:
        fmt = mticker.FuncFormatter(lambda x, _: f"{x * 100:.0f}%")
        if axis == "x":
            ax_obj.xaxis.set_major_formatter(fmt)
        else:
            ax_obj.yaxis.set_major_formatter(fmt)

    # --- Classical frontier ---
    classical_loaded = False
    classical_vols: List[float] = []
    classical_rets: List[float] = []
    classical_sharpes: List[float] = []
    try:
        df = pd.read_csv(classical_csv)
        if "volatility" in df.columns and "return" in df.columns:
            df = df.sort_values("volatility")
            # CSV stores raw daily values; annualize for display.
            classical_vols = (df["volatility"] * np.sqrt(252)).tolist()
            classical_rets = (df["return"] * 252).tolist()
            classical_sharpes = df.get("sharpe", pd.Series(dtype=float)).tolist()
            ax.plot(
                classical_vols, classical_rets,
                color="#2563eb", linewidth=2, label="Markowitz efficient frontier",
                zorder=3,
            )
            if classical_vols:
                min_idx = int(df["volatility"].idxmin())
                ax.scatter(
                    classical_vols[min_idx], classical_rets[min_idx],
                    marker="v", color="#2563eb", s=100, zorder=5,
                    label="Classical min-variance",
                )
            if classical_sharpes:
                max_sh_idx = int(pd.Series(classical_sharpes).idxmax())
                ax.scatter(
                    classical_vols[max_sh_idx], classical_rets[max_sh_idx],
                    marker="*", color="#2563eb", s=200, zorder=5,
                    label="Classical max-Sharpe",
                )
            classical_loaded = True
    except Exception as exc:
        ax.annotate(
            f"Classical frontier unavailable: {exc}",
            xy=(0.5, 0.95), xycoords="axes fraction",
            ha="center", fontsize=9, color="#6b7280",
        )

    # --- QAMOO frontier ---
    # Filter out degenerate all-zero solutions (zero volatility = zero-weight bitstring).
    frontier_points: List[Dict[str, Any]] = [
        p for p in qamoo_result.get("frontier_points", [])
        if p.get("portfolio_volatility", 0.0) > 1e-9
    ]
    qamoo_vols = [p["portfolio_volatility"] for p in frontier_points]
    qamoo_rets = [p["portfolio_return"] for p in frontier_points]

    if not frontier_points:
        ax.annotate(
            "QAMOO: no valid frontier points (circuit returned trivial all-zero solutions)",
            xy=(0.5, 0.5), xycoords="axes fraction",
            ha="center", fontsize=9, color="#6b7280",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#f3f4f6",
                      edgecolor="#d1d5db", alpha=0.9),
        )

    if qamoo_vols:
        ax.scatter(
            qamoo_vols, qamoo_rets,
            color="#7c3aed", s=60, alpha=0.8, zorder=4,
            label="QAMOO quantum frontier",
        )
        # Mark max-Sharpe QAMOO point
        rf = 0.0
        sharpes_q = [
            r / max(v, 1e-10) for v, r in zip(qamoo_vols, qamoo_rets)
        ]
        best_idx = int(max(range(len(sharpes_q)), key=lambda i: sharpes_q[i]))
        ax.scatter(
            qamoo_vols[best_idx], qamoo_rets[best_idx],
            marker="*", color="#7c3aed", s=200, zorder=5,
            label="QAMOO max-Sharpe",
        )

    # NISQ annotation when mean QAMOO volatility exceeds classical by > 10%
    if classical_loaded and qamoo_vols and classical_vols:
        mean_qamoo_vol = sum(qamoo_vols) / len(qamoo_vols)
        mean_classical_vol = sum(classical_vols) / len(classical_vols)
        if mean_classical_vol > 0 and mean_qamoo_vol > mean_classical_vol * 1.10:
            ax.annotate(
                "QAMOO results reflect current NISQ hardware noise at this circuit depth",
                xy=(0.5, 0.03), xycoords="axes fraction",
                ha="center", fontsize=9, color="#92400e",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="#fef3c7",
                          edgecolor="#f59e0b", alpha=0.9),
            )

    ax.set_xlabel("Annualized volatility", fontsize=12)
    ax.set_ylabel("Annualized return", fontsize=12)
    ax.set_title("Quantum vs Classical Efficient Frontier", fontsize=13)
    _pct(ax, "x")
    _pct(ax, "y")
    ax.grid(True, color="#d1d5db", alpha=0.3)
    ax.legend(loc="upper left", fontsize=10)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def build_quantum_report(
    metrics_by_solver: Dict[str, Dict[str, Any]],
    comparison_png: str,
    scaling_png: str,
    output_path: str,
    title: str,
    frontier_png: Optional[str] = None,
    frontier_classical_csv: Optional[str] = None,
    frontier_qamoo_result: Optional[Dict[str, Any]] = None,
) -> None:
    """Build self-contained HTML report.

    Args:
        metrics_by_solver: Metrics grouped by solver (may include Phase 4 fields
            such as circuit_depth, backend_calibration_date, error_mitigation_method,
            and signal_quality when augmented from quantum_result files).
        comparison_png: Path to comparison equity curves PNG.
        scaling_png: Path to scaling table PNG.
        output_path: Path to write HTML.
        title: Report title.
        frontier_png: Optional path to frontier comparison PNG.
        frontier_classical_csv: Optional path to efficient_frontier.csv.
        frontier_qamoo_result: Optional loaded QAMOO result dict.
    """
    def encode_png(path: str) -> str:
        try:
            with open(path, "rb") as fh:
                b64 = base64.b64encode(fh.read()).decode("ascii")
                return f"data:image/png;base64,{b64}"
        except Exception:
            return ""

    equity_data_uri = encode_png(comparison_png)
    scaling_data_uri = encode_png(scaling_png)
    frontier_data_uri = encode_png(frontier_png) if frontier_png else ""

    # Frontier comparison summary table
    has_frontier = bool(frontier_png and frontier_data_uri)
    frontier_summary: List[List[str]] = []
    frontier_note = (
        "Points above the classical frontier are not achievable with the given assets. "
        "Points below it indicate the current cost of using quantum methods on NISQ hardware. "
        "As hardware improves and circuit depths increase, the quantum frontier is expected "
        "to converge toward the classical result."
    )
    if has_frontier and frontier_classical_csv and frontier_qamoo_result:
        try:
            df_cls = pd.read_csv(frontier_classical_csv)
            if "volatility" in df_cls.columns and "return" in df_cls.columns:
                min_vol_row = df_cls.loc[df_cls["volatility"].idxmin()]
                max_sh_row = df_cls.loc[df_cls["sharpe"].idxmax()] if "sharpe" in df_cls.columns else min_vol_row
                # CSV stores raw daily values; annualize for display.
                frontier_summary.append([
                    "Classical (Markowitz)",
                    f"{float(min_vol_row['return']) * 252 * 100:.1f}%",
                    f"{float(min_vol_row['volatility']) * np.sqrt(252) * 100:.1f}%",
                    f"{float(min_vol_row.get('sharpe', 0)):.3f}",
                    f"{float(max_sh_row['return']) * 252 * 100:.1f}%",
                    f"{float(max_sh_row['volatility']) * np.sqrt(252) * 100:.1f}%",
                    f"{float(max_sh_row.get('sharpe', 0)):.3f}",
                ])
        except Exception:
            pass
        fp_list: List[Dict[str, Any]] = frontier_qamoo_result.get("frontier_points", [])
        if fp_list:
            min_v = min(fp_list, key=lambda p: p["portfolio_volatility"])
            best_sh = max(
                fp_list,
                key=lambda p: p["portfolio_return"] / max(p["portfolio_volatility"], 1e-10),
            )
            frontier_summary.append([
                "QAMOO (Quantum)",
                f"{float(min_v['portfolio_return']) * 100:.1f}%",
                f"{float(min_v['portfolio_volatility']) * 100:.1f}%",
                f"{float(min_v['portfolio_return']) / max(float(min_v['portfolio_volatility']), 1e-10):.3f}",
                f"{float(best_sh['portfolio_return']) * 100:.1f}%",
                f"{float(best_sh['portfolio_volatility']) * 100:.1f}%",
                f"{float(best_sh['portfolio_return']) / max(float(best_sh['portfolio_volatility']), 1e-10):.3f}",
            ])

    # Per-solver metrics table rows
    rows = []
    for solver, metrics in metrics_by_solver.items():
        perf = metrics.get("performance", {})
        rows.append([
            _display_name(solver),
            _display_type(metrics.get("solver_type", "")),
            metrics.get("execution_backend", ""),
            f"{perf.get('sharpe_ratio', 0):.3f}",
            f"{perf.get('total_return', 0) * 100:.1f}%",
            f"{perf.get('annualized_volatility', 0) * 100:.1f}%",
            f"{perf.get('max_drawdown', 0) * 100:.1f}%"
        ])

    # Hardware notes: collect IBM runs that have Phase 4 metadata
    ibm_runs = [
        (solver, metrics)
        for solver, metrics in metrics_by_solver.items()
        if metrics.get("execution_backend", "").startswith("ibm_")
    ]
    has_ibm_runs = bool(ibm_runs)

    hw_rows = []
    has_low_signal = False
    for solver, metrics in ibm_runs:
        sq = metrics.get("signal_quality", "")
        if sq == "low":
            has_low_signal = True
        hw_rows.append([
            _display_name(solver),
            metrics.get("execution_backend", ""),
            str(metrics.get("circuit_depth", "N/A")),
            metrics.get("backend_calibration_date", "N/A"),
            metrics.get("error_mitigation_method", "N/A"),
            "ok" if sq == "ok" else ("low — results may be noise-dominated" if sq == "low" else sq or "N/A"),
        ])

    template_str = """<!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>{{ title }}</title>
    <style>
        *, *::before, *::after { box-sizing: border-box; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            font-size: 14px;
            line-height: 1.5;
            background: #f8fafc;
            color: #111827;
            margin: 0;
            padding: 24px;
        }
        .page { max-width: 1100px; margin: 0 auto; }
        h1 {
            font-size: 1.6rem;
            font-weight: 700;
            margin: 0 0 4px;
            color: #0f172a;
        }
        .subtitle {
            font-size: 0.85rem;
            color: #6b7280;
            margin-bottom: 32px;
        }
        h2 {
            font-size: 1rem;
            font-weight: 600;
            color: #374151;
            margin: 0 0 12px;
            padding-bottom: 6px;
            border-bottom: 2px solid #e5e7eb;
        }
        .section {
            background: #ffffff;
            border: 1px solid #e5e7eb;
            border-radius: 8px;
            padding: 20px 24px;
            margin-bottom: 24px;
        }
        .callout {
            background: #fef3c7;
            border: 1px solid #f59e0b;
            border-radius: 6px;
            padding: 12px 16px;
            margin-bottom: 20px;
            font-size: 13px;
            color: #92400e;
        }
        .chart-wrap { overflow-x: auto; }
        .chart-wrap img { display: block; max-width: 100%; height: auto; }
        .table-wrap { overflow-x: auto; }
        table.metrics {
            border-collapse: collapse;
            width: 100%;
            min-width: 560px;
            font-size: 13px;
        }
        table.metrics th {
            background: #f1f5f9;
            color: #374151;
            font-weight: 600;
            text-align: left;
            padding: 9px 14px;
            border-bottom: 2px solid #e2e8f0;
            white-space: nowrap;
        }
        table.metrics td {
            padding: 8px 14px;
            border-bottom: 1px solid #f1f5f9;
            white-space: nowrap;
        }
        table.metrics tbody tr:hover { background: #f8fafc; }
    </style>
</head>
<body>
<div class="page">
    <h1>{{ title }}</h1>
    <p class="subtitle">Quantum vs classical portfolio optimization benchmark</p>

    {% if has_low_signal %}
    <div class="callout">
        <strong>Notice:</strong> One or more IBM hardware results show low signal quality
        (top bitstring &lt; 5% of shots). These results reflect current NISQ hardware noise
        rather than algorithmic performance. Aer simulation results provide the noise-free baseline.
    </div>
    {% endif %}

    <div class="section">
        <h2>Solver Comparison — Equity Curves</h2>
        <div class="chart-wrap">
        {% if equity_data_uri %}
            <img alt="Equity Curves" src="{{ equity_data_uri }}">
        {% else %}
            <p>Equity chart not available.</p>
        {% endif %}
        </div>
    </div>

    <div class="section">
        <h2>Scaling Table</h2>
        <div class="chart-wrap">
        {% if scaling_data_uri %}
            <img alt="Scaling Table" src="{{ scaling_data_uri }}">
        {% else %}
            <p>Scaling table not available.</p>
        {% endif %}
        </div>
    </div>

    <div class="section">
        <h2>Per-Solver Metrics</h2>
        <div class="table-wrap">
        <table class="metrics">
            <thead>
                <tr>
                    <th>Solver</th>
                    <th>Type</th>
                    <th>Backend</th>
                    <th>Sharpe</th>
                    <th>Total Return</th>
                    <th>Volatility</th>
                    <th>Max Drawdown</th>
                </tr>
            </thead>
            <tbody>
                {% for row in rows %}
                <tr>
                    <td>{{ row[0] }}</td>
                    <td>{{ row[1] }}</td>
                    <td>{{ row[2] }}</td>
                    <td>{{ row[3] }}</td>
                    <td>{{ row[4] }}</td>
                    <td>{{ row[5] }}</td>
                    <td>{{ row[6] }}</td>
                </tr>
                {% endfor %}
            </tbody>
        </table>
        </div>
    </div>

    {% if has_ibm_runs %}
    <div class="section">
        <h2>Hardware notes</h2>
        <div class="table-wrap">
        <table class="metrics">
            <thead>
                <tr>
                    <th>Solver</th>
                    <th>Backend</th>
                    <th>Circuit Depth</th>
                    <th>Calibration Date</th>
                    <th>Error Mitigation</th>
                    <th>Signal Quality</th>
                </tr>
            </thead>
            <tbody>
                {% for row in hw_rows %}
                <tr>
                    <td>{{ row[0] }}</td>
                    <td>{{ row[1] }}</td>
                    <td>{{ row[2] }}</td>
                    <td>{{ row[3] }}</td>
                    <td>{{ row[4] }}</td>
                    <td>{{ row[5] }}</td>
                </tr>
                {% endfor %}
            </tbody>
        </table>
        </div>
    </div>
    {% endif %}

    {% if has_frontier %}
    <div class="section">
        <h2>Frontier comparison</h2>
        <div class="chart-wrap">
        {% if frontier_data_uri %}
            <img alt="Frontier Comparison" src="{{ frontier_data_uri }}">
        {% else %}
            <p>Frontier chart not available.</p>
        {% endif %}
        </div>
        {% if frontier_summary %}
        <div class="table-wrap" style="margin-top:16px">
        <table class="metrics">
            <thead>
                <tr>
                    <th>Method</th>
                    <th>Min-Vol Return</th>
                    <th>Min-Vol Volatility</th>
                    <th>Min-Vol Sharpe</th>
                    <th>Max-Sharpe Return</th>
                    <th>Max-Sharpe Volatility</th>
                    <th>Max-Sharpe Ratio</th>
                </tr>
            </thead>
            <tbody>
                {% for row in frontier_summary %}
                <tr>
                    {% for cell in row %}<td>{{ cell }}</td>{% endfor %}
                </tr>
                {% endfor %}
            </tbody>
        </table>
        </div>
        {% endif %}
        <p style="margin-top:12px; font-size:13px; color:#374151;">{{ frontier_note }}</p>
    </div>
    {% endif %}

</div>
</body>
</html>"""

    tmpl = Template(template_str)
    html = tmpl.render(
        title=title,
        rows=rows,
        equity_data_uri=equity_data_uri,
        scaling_data_uri=scaling_data_uri,
        has_ibm_runs=has_ibm_runs,
        hw_rows=hw_rows,
        has_low_signal=has_low_signal,
        has_frontier=has_frontier,
        frontier_data_uri=frontier_data_uri,
        frontier_summary=frontier_summary,
        frontier_note=frontier_note,
    )

    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(html)


def load_comparison_results(file_path: Path) -> Dict[str, Any]:
    """Load comparison results JSON."""
    with file_path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate quantum benchmark report")
    parser.add_argument("--comparison-file", type=Path, default=Path("results/comparison_results.json"),
                        help="Path to comparison_results.json")
    parser.add_argument("--output-dir", type=Path, default=Path("results"),
                        help="Output directory")
    parser.add_argument("--title", type=str, default="Quantum Benchmark Report",
                        help="Report title")
    parser.add_argument(
        "--frontier-classical",
        type=Path,
        default=None,
        help="Path to efficient_frontier.csv for frontier comparison.",
    )
    parser.add_argument(
        "--frontier-quantum",
        type=Path,
        default=None,
        help="Path to QAMOO result JSON for frontier comparison.",
    )

    args = parser.parse_args()

    if not args.comparison_file.exists():
        raise FileNotFoundError(f"Comparison file not found: {args.comparison_file}")

    data = load_comparison_results(args.comparison_file)
    runs = data.get("runs", [])

    # Augment runs with Phase 4 fields from quantum_result_*.json files
    _try_augment_runs_from_result_files(runs, args.output_dir)

    # Fill any missing solution_quality_vs_classical values from Sharpe ratios
    _fill_solution_quality(runs)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Generate comparison plots
    equity_png = args.output_dir / "comparison_equity_curves.png"
    plot_comparison_equity_curves(runs, str(equity_png))

    scaling_png = args.output_dir / "scaling_table.png"
    plot_scaling_table(runs, str(scaling_png))

    # Resolve default paths for frontier comparison
    classical_csv = args.frontier_classical
    quantum_json = args.frontier_quantum
    if classical_csv is None:
        default_cls = args.output_dir.parent / "results" / "efficient_frontier.csv"
        if default_cls.exists():
            classical_csv = default_cls
    if quantum_json is None:
        default_q = args.output_dir / "quantum_result_latest_qamoo.json"
        if default_q.exists():
            quantum_json = default_q

    frontier_png: Optional[str] = None
    qamoo_result: Optional[Dict[str, Any]] = None
    if classical_csv is not None and classical_csv.exists() \
            and quantum_json is not None and quantum_json.exists():
        try:
            with quantum_json.open("r", encoding="utf-8") as fh:
                qamoo_result = json.load(fh)
            frontier_png_path = args.output_dir / "frontier_comparison.png"
            plot_frontier_comparison(str(classical_csv), qamoo_result, str(frontier_png_path))
            frontier_png = str(frontier_png_path)
            print(f"Wrote frontier comparison chart to: {frontier_png_path}")
        except Exception as exc:
            print(f"Warning: frontier comparison failed: {exc}")

    # Group metrics by solver
    metrics_by_solver = {run["solver_name"]: run for run in runs}

    # Build HTML
    html_path = args.output_dir / "quantum_report.html"
    build_quantum_report(
        metrics_by_solver,
        str(equity_png),
        str(scaling_png),
        str(html_path),
        args.title,
        frontier_png=frontier_png,
        frontier_classical_csv=str(classical_csv) if classical_csv else None,
        frontier_qamoo_result=qamoo_result,
    )

    print(f"Wrote report to: {html_path}")


if __name__ == "__main__":
    main()
