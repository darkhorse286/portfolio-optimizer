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
import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from jinja2 import Template

# Import aggregation helpers from the sibling script without requiring a package.
sys.path.insert(0, str(Path(__file__).parent))
try:
    from benchmark_aggregate import (
        discover_files, load_runs, aggregate_solver,
        normalize_name, build_md_table,
    )
    _HAS_AGGREGATE = True
except ImportError:
    _HAS_AGGREGATE = False

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


# Walk-forward solver names registered in main.cpp.
# Used to suppress stale single-submission Aer result files.
WALK_FORWARD_SOLVER_NAMES: set[str] = {
    "QAOA Uninformed (Aer)",
    "QAOA Informed (Aer)",
    "QAMO Uninformed (Aer)",
    "QAMO Informed (Aer)",
}

# Legacy single-submission solver names written by older main.cpp builds.
# These are superseded by the walk-forward variants above and must not be
# injected into the report even when found in quantum_result_*.json files.
_LEGACY_AER_SOLVER_NAMES: set[str] = {
    "qaoa_p1_aer_simulator",
    "qamo_p1_aer_simulator",
    # Old QAMOO naming scheme from pre-Phase-5 result files.
    # These are used only as frontier_points sources; never inject as solver rows.
    "qamoo_20pts_aer_simulator",
    "qamoo_5pts_ibm_fez",
    "qamoo_20pts_ibm_fez",
    "qamoo_20pts_ibm_kingston",
    "qamoo_20pts_ibm_marrakesh",
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

    # Case-insensitive prefix check so both legacy snake_case names (e.g.
    # qaoa_p1_aer_simulator) and the new human-readable names (e.g.
    # "QAOA Uninformed (Aer)") resolve to the same colour family.
    # Check longest prefixes first (qamoo before qamo).
    solver_lower = solver_name.lower()
    if solver_lower.startswith("qamoo"):
        return SOLVER_COLOR_MAP[f"qamoo_{suffix}"]
    if solver_lower.startswith("qamo"):
        return SOLVER_COLOR_MAP[f"qamo_{suffix}"]
    if solver_lower.startswith("qaoa"):
        return SOLVER_COLOR_MAP[f"qaoa_{suffix}"]
    if solver_lower.startswith("markowitz"):
        return SOLVER_COLOR_MAP["markowitz"]
    if solver_lower.startswith("sa_classical"):
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


def _get_sector_mapping() -> Dict[str, str]:
    return {
        "AAPL": "Technology",
        "MSFT": "Technology",
        "GOOGL": "Technology",
        "JPM": "Financials",
        "BAC": "Financials",
        "JNJ": "Healthcare",
        "PFE": "Healthcare",
        "XOM": "Energy",
        "CVX": "Energy",
        "WMT": "Consumer",
    }


def _compute_sector_allocations(
    universe: List[str],
    weights: List[float],
    prices_csv: Path,
) -> Optional[List[Dict[str, Any]]]:
    if not universe or not weights or len(universe) != len(weights) or not prices_csv.exists():
        return None

    try:
        prices = pd.read_csv(prices_csv, index_col="date", parse_dates=True)
    except Exception:
        return None

    cols = [c for c in universe if c in prices.columns]
    if not cols:
        return None

    p = prices[cols].dropna()
    if len(p) < 2:
        return None

    mapping = _get_sector_mapping()
    sector_weights: dict[str, float] = {}
    sector_returns: dict[str, float] = {}
    sector_asset_weights: dict[str, list[tuple[str, float]]] = {}

    for ticker, weight in zip(universe, weights):
        sector = mapping.get(ticker)
        if sector is None:
            continue
        sector_weights[sector] = sector_weights.get(sector, 0.0) + weight
        sector_asset_weights.setdefault(sector, []).append((ticker, weight))

    if not sector_weights:
        return None

    # Use total return per asset over the available history window.
    asset_returns: dict[str, float] = {}
    for ticker in cols:
        first = float(p[ticker].iloc[0])
        last = float(p[ticker].iloc[-1])
        asset_returns[ticker] = (last / first - 1.0) if first > 0 else 0.0

    allocations: List[Dict[str, Any]] = []
    total_weight = sum(weights)
    if total_weight <= 0.0:
        total_weight = 1.0

    for sector, sector_w in sector_weights.items():
        asset_list = sector_asset_weights.get(sector, [])
        sector_ret = 0.0
        if sector_w > 0 and asset_list:
            for ticker, w in asset_list:
                asset_return = asset_returns.get(ticker, 0.0)
                sector_ret += asset_return * (w / sector_w)
        allocations.append({
            "sector_name": sector,
            "weight": sector_w / total_weight,
            "return_value": sector_ret,
        })

    return allocations


def _compute_benchmark_sector_allocations(
    universe: List[str],
    prices_csv: Path,
) -> Optional[List[Dict[str, Any]]]:
    if not universe or not prices_csv.exists():
        return None

    allocations = _compute_sector_allocations(universe, [1.0 / len(universe)] * len(universe), prices_csv)
    if allocations is None:
        return None

    # Re-weight sectors equally for the benchmark comparison.
    sector_list = [alloc["sector_name"] for alloc in allocations]
    sector_set = sorted(set(sector_list))
    if not sector_set:
        return None

    sector_return_map = {alloc["sector_name"]: alloc["return_value"] for alloc in allocations}
    return [
        {
            "sector_name": sector,
            "weight": 1.0 / len(sector_set),
            "return_value": sector_return_map.get(sector, 0.0),
        }
        for sector in sector_set
    ]


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

    # Case 2 deduplication: track best (most-recent) candidate per solver_name.
    # submitted_at is often absent; filename (lexicographic order) serves as proxy.
    case2_candidates: Dict[str, Dict[str, Any]] = {}
    case2_sort_keys: Dict[str, str] = {}

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
                    if key in qr and not run.get(key):
                        run[key] = qr[key]
                # Resolve "ibm_auto" to the actual backend from the result file
                if qr.get("execution_backend") and run.get("execution_backend") == "ibm_auto":
                    run["execution_backend"] = qr["execution_backend"]

        # Accumulate Case 2 candidates — skip legacy/walk-forward names and
        # solvers already present in comparison_results.json.
        if (solver_name
                and solver_name not in _LEGACY_AER_SOLVER_NAMES
                and solver_name not in WALK_FORWARD_SOLVER_NAMES
                and solver_name not in existing_names
                and qr.get("status") == "COMPLETED"):
            sort_key = qr.get("submitted_at") or result_file.name
            if sort_key > case2_sort_keys.get(solver_name, ""):
                case2_candidates[solver_name] = qr
                case2_sort_keys[solver_name] = sort_key

    # Case 2: inject one entry per new solver_name using the most-recent file.
    for solver_name, qr in case2_candidates.items():
        backend = qr.get("execution_backend", "")
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

        sector_allocations = None
        benchmark_sector_allocations = None
        trade_summary = None
        if weights and universe and prices_csv.exists():
            sector_allocations = _compute_sector_allocations(universe, weights, prices_csv)
            benchmark_sector_allocations = _compute_benchmark_sector_allocations(universe, prices_csv)
            trade_summary = {
                "rebalance_count": None,
                "turnover": None,
                "avg_cost_per_trade": None,
                "total_costs": None,
            }

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
            "completed_at": qr.get("completed_at", ""),
            "_hardware_only": not has_full_backtest,
            "_synthetic_nav": True,  # Mark quantum_result injections as synthetic so they render dashed
        }
        if trade_summary is not None:
            injected["trade_summary"] = trade_summary
        if sector_allocations is not None:
            injected["sector_weights"] = sector_allocations
        if benchmark_sector_allocations is not None:
            injected["benchmark_sector_weights"] = benchmark_sector_allocations

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


def _fill_synthetic_nav_for_runs(
    runs: List[Dict[str, Any]],
    prices_csv: Path,
    risk_free_rate: float = 0.02,
) -> None:
    """Compute buy-and-hold NAV series for runs with avg_weights but empty nav_series."""
    for run in runs:
        if run.get("nav_series"):
            continue
        weights = run.get("avg_weights", [])
        tickers = run.get("tickers", [])
        if not weights or not tickers or len(weights) != len(tickers):
            continue
        result = _compute_static_portfolio_metrics(weights, tickers, prices_csv, risk_free_rate)
        if result is None:
            continue
        perf, nav_series = result
        run["nav_series"] = nav_series
        run["performance"] = perf
        run["_synthetic_nav"] = True


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
            is_synthetic = run.get("_synthetic_nav", False)
            plt.plot(
                x_axis, normalized,
                label=_display_name(solver_name),
                color=_solver_color(run),
                linestyle="--" if is_synthetic else "-",
            )

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

    # Determine whether the NISQ disclaimer should be shown
    show_nisq = False
    if classical_loaded and qamoo_vols and classical_vols:
        mean_qamoo_vol = sum(qamoo_vols) / len(qamoo_vols)
        mean_classical_vol = sum(classical_vols) / len(classical_vols)
        if mean_classical_vol > 0 and mean_qamoo_vol > mean_classical_vol * 1.10:
            show_nisq = True

    ax.set_xlabel("Annualized volatility", fontsize=12)
    ax.set_ylabel("Annualized return", fontsize=12)
    ax.set_title("Quantum vs Classical Efficient Frontier", fontsize=13)
    _pct(ax, "x")
    _pct(ax, "y")
    ax.grid(True, color="#d1d5db", alpha=0.3)
    legend = ax.legend(loc="upper left", fontsize=10)
    plt.tight_layout()

    # Place NISQ disclaimer just below the legend box
    if show_nisq:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        leg_bb = legend.get_window_extent(renderer=renderer)
        ax_inv = ax.transAxes.inverted()
        x_left = ax_inv.transform((leg_bb.x0, leg_bb.y0))[0]
        y_below = ax_inv.transform((leg_bb.x0, leg_bb.y0 - 6))[1]
        ax.annotate(
            "NISQ: results reflect current hardware noise at this circuit depth",
            xy=(x_left, y_below), xycoords="axes fraction",
            ha="left", va="top", fontsize=9, color="#92400e",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#fef3c7",
                      edgecolor="#f59e0b", alpha=0.9),
        )

    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def _parse_sector_allocations(sectors: Any) -> Optional[List[Dict[str, Any]]]:
    """Parse a sector allocation structure from run metadata."""
    if not sectors or not isinstance(sectors, list):
        return None
    parsed: List[Dict[str, Any]] = []
    for item in sectors:
        if not isinstance(item, dict):
            continue
        name = item.get("sector_name") or item.get("sector") or item.get("name")
        weight = item.get("weight")
        ret = item.get("return_value")
        if name is None or weight is None or ret is None:
            continue
        try:
            parsed.append({
                "sector_name": str(name),
                "weight": float(weight),
                "return_value": float(ret),
            })
        except (TypeError, ValueError):
            continue
    return parsed if parsed else None


def _compute_brinson_fachler(
    portfolio_sectors: List[Dict[str, Any]],
    benchmark_sectors: List[Dict[str, Any]],
) -> Optional[Dict[str, float]]:
    """Compute Brinson-Fachler single-period allocation, selection, interaction."""
    if not portfolio_sectors or not benchmark_sectors:
        return None

    bench_map = {s["sector_name"]: s for s in benchmark_sectors}
    for sector in portfolio_sectors:
        if sector["sector_name"] not in bench_map:
            return None

    total_bench_return = sum(s["weight"] * s["return_value"] for s in benchmark_sectors)
    total_port_return = sum(s["weight"] * s["return_value"] for s in portfolio_sectors)

    allocation = 0.0
    selection = 0.0
    interaction = 0.0
    for port_sector in portfolio_sectors:
        bench_sector = bench_map[port_sector["sector_name"]]
        w_p = port_sector["weight"]
        w_b = bench_sector["weight"]
        r_p = port_sector["return_value"]
        r_b = bench_sector["return_value"]

        allocation += (w_p - w_b) * (r_b - total_bench_return)
        selection += w_b * (r_p - r_b)
        interaction += (w_p - w_b) * (r_p - r_b)

    total_active = allocation + selection + interaction
    return {
        "allocation": allocation,
        "selection": selection,
        "interaction": interaction,
        "total_active": total_active,
        "portfolio_return": total_port_return,
        "benchmark_return": total_bench_return,
    }


def _augment_cpp_sector_data(runs: List[Dict[str, Any]]) -> None:
    """Convert dict-format sector_weights from C++ JSON runs to list format with return values.

    The C++ benchmark serializes sector_weights as {"Technology": 0.4, ...}.
    _parse_sector_allocations expects [{sector_name, weight, return_value}, ...].
    Reads historical prices to fill in per-sector returns. Mutates runs in-place.
    """
    repo_root = Path(__file__).resolve().parent.parent.parent
    prices_csv = repo_root / "data" / "market" / "historical_prices.csv"

    sector_map = _get_sector_mapping()
    sector_tickers: dict[str, list[str]] = {}
    for ticker, sector in sector_map.items():
        sector_tickers.setdefault(sector, []).append(ticker)

    sector_returns: dict[str, float] = {s: 0.0 for s in sector_tickers}
    if prices_csv.exists():
        try:
            prices = pd.read_csv(prices_csv, index_col="date", parse_dates=True)
            for sector, tickers in sector_tickers.items():
                cols = [t for t in tickers if t in prices.columns]
                if not cols:
                    continue
                p = prices[cols].dropna()
                if len(p) < 2:
                    continue
                rets = []
                for ticker in cols:
                    first = float(p[ticker].iloc[0])
                    last = float(p[ticker].iloc[-1])
                    if first > 0:
                        rets.append(last / first - 1.0)
                if rets:
                    sector_returns[sector] = sum(rets) / len(rets)
        except Exception:
            pass

    for run in runs:
        sw = run.get("sector_weights")
        bsw = run.get("benchmark_sector_weights")
        avg_weights = run.get("avg_weights", [])
        tickers = run.get("tickers", [])

        # Compute per-solver portfolio sector returns using actual asset weights
        portfolio_sector_returns: dict[str, float] = {}
        portfolio_sector_weights_within: dict[str, list[tuple[str, float]]] = {}

        if avg_weights and tickers and len(avg_weights) == len(tickers) and prices_csv.exists():
            try:
                prices = pd.read_csv(prices_csv, index_col="date", parse_dates=True)
                asset_returns: dict[str, float] = {}
                for ticker in tickers:
                    if ticker in prices.columns:
                        col = prices[ticker].dropna()
                        if len(col) >= 2:
                            asset_returns[ticker] = (col.iloc[-1] / col.iloc[0]) - 1.0

                for ticker, weight in zip(tickers, avg_weights):
                    sector = sector_map.get(ticker)
                    if sector is None:
                        continue
                    portfolio_sector_weights_within.setdefault(sector, []).append((ticker, weight))

                for sector, asset_list in portfolio_sector_weights_within.items():
                    total_sector_w = sum(w for _, w in asset_list)
                    if total_sector_w > 0:
                        portfolio_sector_returns[sector] = sum(
                            asset_returns.get(t, 0.0) * (w / total_sector_w)
                            for t, w in asset_list
                        )
            except Exception:
                pass

        if isinstance(sw, dict) and sw:
            total_w = sum(sw.values()) or 1.0
            run["sector_weights"] = [
                {
                    "sector_name": sector,
                    "weight": weight / total_w,
                    "return_value": portfolio_sector_returns.get(sector, sector_returns.get(sector, 0.0)),
                }
                for sector, weight in sw.items()
            ]

        if isinstance(bsw, dict) and bsw:
            total_bw = sum(bsw.values()) or 1.0
            run["benchmark_sector_weights"] = [
                {
                    "sector_name": sector,
                    "weight": weight / total_bw,
                    "return_value": sector_returns.get(sector, 0.0),  # benchmark stays equal-weighted
                }
                for sector, weight in bsw.items()
            ]


def _format_percentage(value: Any, precision: int = 2) -> str:
    try:
        return f"{float(value) * 100:.{precision}f}%"
    except (TypeError, ValueError):
        return "n/a"


def _format_currency(value: Any) -> str:
    try:
        return f"{float(value):.2f}"
    except (TypeError, ValueError):
        return "n/a"


def _build_agg_rows(results_dir: Path) -> List[Dict[str, Any]]:
    """Return per-solver aggregate stats for the HTML report section.

    Source 1: timestamped comparison_results_[0-9]*.json archives (walk-forward Aer).
    Source 2: quantum_result_*.json files (IBM hardware, QAMOO Aer) — metrics are
              computed via static-weight backtest because these files carry weights
              only, with no embedded performance stats.

    Also writes aggregated_benchmark.json and aggregated_benchmark.md to results_dir
    as a side effect so running benchmark_viz.py is the only command needed.

    Returns an empty list when no data is available or the aggregate module
    could not be imported, so the section is simply omitted.
    """
    if not _HAS_AGGREGATE:
        return []

    # --- Part 1: walk-forward results from timestamped archives ---
    files = discover_files(results_dir)
    by_solver: Dict[str, List[Dict[str, Any]]] = load_runs(files) if files else {}
    # Snapshot names from archives so Part 2 can skip them without also
    # blocking accumulation of multiple files for the same IBM solver.
    archive_solvers: set[str] = set(by_solver.keys())

    # --- Part 2: IBM hardware and QAMOO Aer from quantum_result_*.json ---
    # Each file is one submission. Accumulate all files per solver so σ reflects
    # genuine run-to-run variance. Metrics are computed via static-weight backtest.
    repo_root = Path(__file__).resolve().parent.parent.parent
    config_path = repo_root / "data" / "config" / "portfolio_config.json"
    prices_csv  = repo_root / "data" / "market" / "historical_prices.csv"
    risk_free_rate = 0.02
    config_universe: List[str] = []
    if config_path.exists():
        try:
            cfg = json.loads(config_path.read_text(encoding="utf-8"))
            risk_free_rate  = cfg.get("performance", {}).get("risk_free_rate", 0.02)
            config_universe = cfg.get("data", {}).get("universe", [])
        except Exception:
            pass

    for qr_file in sorted(results_dir.glob("quantum_result_*.json")):
        try:
            qr = json.loads(qr_file.read_text(encoding="utf-8"))
        except Exception:
            continue

        raw_name = qr.get("solver_name", "")
        if not raw_name:
            continue
        if raw_name in _LEGACY_AER_SOLVER_NAMES:
            continue
        if raw_name in WALK_FORWARD_SOLVER_NAMES:
            continue
        if qr.get("status") != "COMPLETED":
            continue

        solver_name = normalize_name(raw_name)
        if solver_name in archive_solvers:
            continue  # already covered by walk-forward archives

        weights  = qr.get("weights", [])
        universe = qr.get("universe") or config_universe
        if not weights or not universe or not prices_csv.exists():
            continue

        result = _compute_static_portfolio_metrics(weights, universe, prices_csv, risk_free_rate)
        if result is None:
            continue
        perf, _ = result

        by_solver.setdefault(solver_name, []).append({
            "sharpe_ratio":          perf.get("sharpe_ratio"),
            "total_return":          perf.get("total_return"),
            "annualized_volatility": perf.get("annualized_volatility"),
            "max_drawdown":          perf.get("max_drawdown"),
            "turnover":              None,
            "solve_time_ms":         qr.get("solve_time_ms"),
            "solver_type":           "quantum",
            "execution_backend":     qr.get("execution_backend", ""),
        })

    if not by_solver:
        return []

    stats: Dict[str, Dict[str, Any]] = {
        name: aggregate_solver(entries)
        for name, entries in sorted(by_solver.items())
    }

    # Write side-effect output files so this is the only command needed.
    try:
        (results_dir / "aggregated_benchmark.json").write_text(
            json.dumps({"num_files": len(files), "solvers": stats}, indent=2),
            encoding="utf-8",
        )
        (results_dir / "aggregated_benchmark.md").write_text(
            build_md_table(stats) + "\n", encoding="utf-8",
        )
    except Exception as exc:
        print(f"Warning: could not write aggregate output files: {exc}")

    def _f(v: float, p: int = 3) -> str:
        return f"{v:.{p}f}"

    rows = []
    for solver_name, agg in stats.items():
        sharpe   = agg["sharpe"]
        ret      = agg["total_return"]
        turnover = agg["turnover"]

        high_var = sharpe["high_variance"] or ret["high_variance"] or turnover["high_variance"]
        variance_cell = "⚠ high" if high_var else "ok"

        rows.append({
            "solver":      solver_name,
            "num_runs":    agg["num_runs"],
            "sharpe_mean": _f(sharpe["mean"]),
            "sharpe_std":  _f(sharpe["std"]),
            "return_mean": _f(ret["mean"], 4),
            "return_std":  _f(ret["std"],  4),
            "turn_mean":   _f(turnover["mean"]) if turnover["n"] > 0 else "n/a",
            "turn_std":    _f(turnover["std"])  if turnover["n"] > 0 else "n/a",
            "variance":    variance_cell,
            "high_var":    high_var,
        })
    return rows


def build_quantum_report(
    metrics_by_solver: Dict[str, Dict[str, Any]],
    comparison_png: str,
    scaling_png: str,
    output_path: str,
    title: str,
    frontier_png: Optional[str] = None,
    frontier_classical_csv: Optional[str] = None,
    frontier_qamoo_result: Optional[Dict[str, Any]] = None,
    agg_rows: Optional[List[Dict[str, Any]]] = None,
    has_synthetic_nav: bool = False,
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
    activity_rows = []
    attribution_rows = []
    turnover_warnings = []

    _INFORMED_TIP = (
        "Informed: problem file augmented with EWMA expected_returns and covariance "
        "computed from historical prices before circuit submission."
    )
    _UNINFORMED_TIP = (
        "Uninformed: raw QUBO Q matrix only, no market data augmentation."
    )

    for solver, metrics in metrics_by_solver.items():
        perf = metrics.get("performance", {})
        solver_name = _display_name(solver)
        backend = metrics.get("execution_backend", "")

        # Wrap informed/uninformed names with a tooltip so the distinction is
        # visible on hover without cluttering the table column.
        if "Informed" in solver_name and "Uninformed" not in solver_name:
            solver_name = f'<span title="{_INFORMED_TIP}">{solver_name}</span>'
        elif "Uninformed" in solver_name:
            solver_name = f'<span title="{_UNINFORMED_TIP}">{solver_name}</span>'

        rows.append([
            solver_name,
            _display_type(metrics.get("solver_type", "")),
            backend,
            f"{perf.get('sharpe_ratio', 0):.3f}",
            f"{perf.get('total_return', 0) * 100:.1f}%",
            f"{perf.get('annualized_volatility', 0) * 100:.1f}%",
            f"{perf.get('max_drawdown', 0) * 100:.1f}%"
        ])

        trade_summary = metrics.get("trade_summary", {}) or {}
        rebalance_count = trade_summary.get("rebalance_count")
        skipped_rebalances = trade_summary.get("skipped_rebalances")
        turnover = trade_summary.get("turnover")
        avg_cost = trade_summary.get("avg_cost_per_trade")
        total_costs = trade_summary.get("total_costs")

        is_ibm = backend.startswith("ibm_")
        _hw_na = (
            '<span title="Single hardware submission — re-running the circuit at each '
            'rebalancing date is not feasible at current queue constraints.">'
            'n/a — single hardware submission, not walk-forward</span>'
        )

        warning_row = (
            (isinstance(turnover, (int, float)) and turnover > 1.0) or
            (isinstance(skipped_rebalances, int) and skipped_rebalances > 0)
        )
        if isinstance(turnover, (int, float)) and turnover > 1.0:
            turnover_warnings.append(solver_name)

        if is_ibm:
            activity_rows.append([
                solver_name,
                backend,
                _hw_na, _hw_na, _hw_na, _hw_na, _hw_na,
                False,
            ])
        else:
            activity_rows.append([
                solver_name,
                backend,
                str(rebalance_count) if rebalance_count is not None else "n/a",
                str(skipped_rebalances) if skipped_rebalances is not None else "n/a",
                _format_percentage(turnover, precision=1) if turnover is not None else "n/a",
                _format_currency(avg_cost),
                _format_currency(total_costs),
                warning_row,
            ])

        portfolio_sectors = _parse_sector_allocations(metrics.get("sector_weights"))
        benchmark_sectors = _parse_sector_allocations(metrics.get("benchmark_sector_weights"))
        attribution = None
        if not is_ibm and portfolio_sectors is not None and benchmark_sectors is not None:
            attribution = _compute_brinson_fachler(portfolio_sectors, benchmark_sectors)

        if is_ibm:
            attribution_rows.append([
                solver_name,
                backend,
                _hw_na, _hw_na, _hw_na, _hw_na,
            ])
        elif attribution is None:
            attribution_rows.append([
                solver_name,
                backend,
                '<span title="Attribution unavailable — weights may be outside feasible region.">n/a</span>',
                '<span title="Attribution unavailable — weights may be outside feasible region.">n/a</span>',
                '<span title="Attribution unavailable — weights may be outside feasible region.">n/a</span>',
                '<span title="Attribution unavailable — weights may be outside feasible region.">n/a</span>',
            ])
        else:
            attribution_rows.append([
                solver_name,
                backend,
                _format_percentage(attribution["allocation"], precision=2),
                _format_percentage(attribution["selection"], precision=2),
                _format_percentage(attribution["interaction"], precision=2),
                _format_percentage(attribution["total_active"], precision=2),
            ])

    turnover_warning_note = None
    if turnover_warnings:
        if len(turnover_warnings) == 1:
            turnover_warning_note = (
                f"Warning: turnover exceeded 100% for {turnover_warnings[0]}. "
                "This may indicate the constraint violation path was triggered."
            )
        else:
            solver_list = ", ".join(turnover_warnings)
            turnover_warning_note = (
                f"Warning: turnover exceeded 100% for {solver_list}. "
                "This may indicate the constraint violation path was triggered."
            )

    # Hardware notes: collect all IBM runs, ordered by completed_at ascending
    # so readers see the chronological progression of hardware submissions.
    ibm_runs = sorted(
        [
            (solver, metrics)
            for solver, metrics in metrics_by_solver.items()
            if metrics.get("execution_backend", "").startswith("ibm_")
        ],
        key=lambda t: t[1].get("completed_at", ""),
    )
    has_ibm_runs = bool(ibm_runs)

    hw_rows = []
    has_low_signal = False
    for solver, metrics in ibm_runs:
        sq = metrics.get("signal_quality", "")
        if sq == "low":
            has_low_signal = True
        completed = metrics.get("completed_at", "")
        submitted_display = completed[:10] if completed else "N/A"
        hw_rows.append([
            _display_name(solver),
            metrics.get("execution_backend", ""),
            str(metrics.get("circuit_depth", "N/A")),
            metrics.get("backend_calibration_date", "N/A"),
            metrics.get("error_mitigation_method", "N/A"),
            "ok" if sq == "ok" else ("low — results may be noise-dominated" if sq == "low" else sq or "N/A"),
            submitted_display,
            str(metrics.get("rounds_averaged", "N/A")),
        ])

    has_agg = bool(agg_rows)

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
        .warning-row td { background: #fef3c7; }
        .footnote {
            margin-top: 12px;
            font-size: 13px;
            color: #4b5563;
        }
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

    {% if has_synthetic_nav %}
    <div class="callout">
        <strong>Notice:</strong> IBM hardware solvers produced a single set of optimized weights
        rather than a full walk-forward backtest. Their equity curves (shown as dashed lines) are
        computed by holding those weights constant over the entire backtest period.
        Solid lines represent monthly rebalanced walk-forward backtests.
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
                    <td>{{ row[0] | safe }}</td>
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

    <div class="section">
        <h2>Transaction Activity</h2>
        {% if turnover_warning_note %}
        <div class="callout">
            <strong>Warning:</strong> {{ turnover_warning_note }}
        </div>
        {% endif %}
        <div class="table-wrap">
        <table class="metrics">
            <thead>
                <tr>
                    <th>Solver</th>
                    <th>Backend</th>
                    <th>Rebalance Events</th>
                    <th>Skipped Rebalances</th>
                    <th>Total Turnover</th>
                    <th>Avg Cost Per Trade</th>
                    <th>Total Transaction Costs</th>
                </tr>
            </thead>
            <tbody>
                {% for row in activity_rows %}
                <tr{% if row[7] %} class="warning-row"{% endif %}>
                    <td>{{ row[0] | safe }}</td>
                    <td>{{ row[1] }}</td>
                    <td>{{ row[2] | safe }}</td>
                    <td>{{ row[3] | safe }}</td>
                    <td>{{ row[4] | safe }}</td>
                    <td>{{ row[5] | safe }}</td>
                    <td>{{ row[6] | safe }}</td>
                </tr>
                {% endfor %}
            </tbody>
        </table>
        </div>
    </div>

    <div class="section">
        <h2>Performance Attribution (Brinson-Fachler)</h2>
        <div class="table-wrap">
        <table class="metrics">
            <thead>
                <tr>
                    <th>Solver</th>
                    <th>Backend</th>
                    <th>Allocation Effect</th>
                    <th>Selection Effect</th>
                    <th>Interaction Effect</th>
                    <th>Total Active Return</th>
                </tr>
            </thead>
            <tbody>
                {% for row in attribution_rows %}
                <tr>
                    <td>{{ row[0] | safe }}</td>
                    <td>{{ row[1] }}</td>
                    <td>{{ row[2] | safe }}</td>
                    <td>{{ row[3] | safe }}</td>
                    <td>{{ row[4] | safe }}</td>
                    <td>{{ row[5] | safe }}</td>
                </tr>
                {% endfor %}
            </tbody>
        </table>
        </div>
        <p class="footnote">Single-period decomposition. Multi-period Carino linking available via export_analytics_json().</p>
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
                    <th>Submitted</th>
                    <th>Rounds Averaged</th>
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
                    <td>{{ row[6] }}</td>
                    <td>{{ row[7] }}</td>
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

    {% if has_agg %}
    <div class="section">
        <h2>Aggregate Statistics (across {{ agg_rows|length }} solver × run combinations)</h2>
        <div class="table-wrap">
        <table class="metrics">
            <thead>
                <tr>
                    <th>Solver</th>
                    <th>Runs</th>
                    <th>Mean Sharpe</th>
                    <th title="Standard deviation of Sharpe across runs">±σ</th>
                    <th>Mean Return</th>
                    <th title="Standard deviation of total return across runs">±σ</th>
                    <th>Mean Turnover</th>
                    <th title="Standard deviation of turnover across runs">±σ</th>
                    <th title="⚠ when std / |mean| &gt; 0.5 for any metric">Variance</th>
                </tr>
            </thead>
            <tbody>
                {% for row in agg_rows %}
                <tr{% if row.high_var %} style="background:#fef3c7"{% endif %}>
                    <td>{{ row.solver }}</td>
                    <td>{{ row.num_runs }}</td>
                    <td>{{ row.sharpe_mean }}</td>
                    <td>{{ row.sharpe_std }}</td>
                    <td>{{ row.return_mean }}</td>
                    <td>{{ row.return_std }}</td>
                    <td>{{ row.turn_mean }}</td>
                    <td>{{ row.turn_std }}</td>
                    <td>{{ row.variance }}</td>
                </tr>
                {% endfor %}
            </tbody>
        </table>
        </div>
        <p class="footnote">Classical solvers show σ=0.000 because they are deterministic — this is the control baseline. ⚠ rows exceed 50% coefficient of variation on at least one metric.</p>
    </div>
    {% endif %}

</div>
</body>
</html>"""

    tmpl = Template(template_str)
    html = tmpl.render(
        title=title,
        rows=rows,
        activity_rows=activity_rows,
        attribution_rows=attribution_rows,
        turnover_warning_note=turnover_warning_note,
        equity_data_uri=equity_data_uri,
        scaling_data_uri=scaling_data_uri,
        has_ibm_runs=has_ibm_runs,
        hw_rows=hw_rows,
        has_low_signal=has_low_signal,
        has_synthetic_nav=has_synthetic_nav,
        has_frontier=has_frontier,
        frontier_data_uri=frontier_data_uri,
        frontier_summary=frontier_summary,
        frontier_note=frontier_note,
        has_agg=has_agg,
        agg_rows=agg_rows or [],
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

    # Convert C++ dict-format sector_weights to list format with return values
    _augment_cpp_sector_data(runs)

    # Fill any missing solution_quality_vs_classical values from Sharpe ratios
    _fill_solution_quality(runs)

    # Compute buy-and-hold equity curves for IBM runs that have weights but no nav_series
    _repo_root = Path(__file__).resolve().parent.parent.parent
    _prices_csv = _repo_root / "data" / "market" / "historical_prices.csv"
    _risk_free_rate = 0.02
    _cfg_path = _repo_root / "data" / "config" / "portfolio_config.json"
    if _cfg_path.exists():
        try:
            with _cfg_path.open() as _f:
                _cfg = json.load(_f)
            _risk_free_rate = float(_cfg.get("optimizer", {}).get("risk_free_rate", 0.02))
        except Exception:
            pass
    _fill_synthetic_nav_for_runs(runs, _prices_csv, _risk_free_rate)

    # Update solution_quality now that synthetic NAV has populated performance metrics
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
        for _candidate in [
            args.output_dir / "efficient_frontier.csv",
            args.output_dir.parent / "results" / "efficient_frontier.csv",
            args.output_dir.parent / "efficient_frontier.csv",
        ]:
            if _candidate.exists():
                classical_csv = _candidate
                break
    if quantum_json is None:
        # Scan for any quantum_result file that contains frontier_points;
        # pick the most recent by submitted_at, falling back to filename.
        best_fp_file: Optional[Path] = None
        best_fp_key: str = ""
        for _qr_path in args.output_dir.glob("quantum_result_*.json"):
            try:
                _qr = json.loads(_qr_path.read_text(encoding="utf-8"))
            except Exception:
                continue
            if not _qr.get("frontier_points"):
                continue
            _key = _qr.get("submitted_at") or _qr_path.name
            if _key > best_fp_key:
                best_fp_file = _qr_path
                best_fp_key = _key
        if best_fp_file is not None:
            quantum_json = best_fp_file

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

    has_synthetic_nav = any(run.get("_synthetic_nav") for run in runs)

    # Aggregate stats across timestamped run archives
    agg_rows = _build_agg_rows(args.output_dir)
    if agg_rows:
        print(f"Aggregated {sum(r['num_runs'] for r in agg_rows)} solver-run(s) "
              f"across {len(discover_files(args.output_dir))} archive file(s).")

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
        agg_rows=agg_rows,
        has_synthetic_nav=has_synthetic_nav,
    )

    print(f"Wrote report to: {html_path}")


if __name__ == "__main__":
    main()
