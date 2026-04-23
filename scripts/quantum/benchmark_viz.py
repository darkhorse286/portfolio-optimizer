#!/usr/bin/env python3
"""Benchmark visualization script for quantum portfolio optimizers."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
from jinja2 import Template

SOLVER_NAME_MAP: dict[str, str] = {
    "markowitz": "Markowitz",
    "sa_classical": "QUBO (Sim. Annealing)",
    "qaoa": "QAOA",
}

SOLVER_TYPE_MAP: dict[str, str] = {
    "classical": "Classical",
    "quantum_inspired": "QUBO",
    "quantum": "Quantum (QAOA)",
}


def _display_name(solver_name: str) -> str:
    """Return a human-readable label for a solver name."""
    for key, label in SOLVER_NAME_MAP.items():
        if solver_name.startswith(key):
            return label
    return solver_name


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

        # Case 2: IBM hardware result not yet in comparison_results.json — inject it.
        # Attempt a static-weight backtest so real metrics are shown.
        if solver_name not in existing_names and backend.startswith("ibm_") \
                and qr.get("status") == "COMPLETED":
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


def plot_comparison_equity_curves(runs: List[Dict[str, Any]], output_path: str) -> None:
    """Plot normalized equity curves for all solvers.

    Args:
        runs: List of run dictionaries from comparison_results.json
        output_path: Path to save the PNG file
    """
    plt.figure(figsize=(12, 5))

    color_map = {
        "classical": "#2563eb",
        "quantum_inspired": "#d97706",
        "quantum": "#7c3aed"
    }

    for run in runs:
        solver_type = run["solver_type"]
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
            plt.plot(x_axis, normalized, label=_display_name(solver_name), color=color_map.get(solver_type, "#6b7280"))

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
        row = [
            _display_name(res["solver_name"]),
            _display_type(res["solver_type"]),
            backend,
            "hw only" if is_hw_only else f"{res['performance']['sharpe_ratio']:.3f}",
            "hw only" if is_hw_only else f"{res['performance']['total_return'] * 100:.1f}%",
            "hw only" if is_hw_only else f"{res['solution_quality_vs_classical']:.2f}x",
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
        quality = res["solution_quality_vs_classical"]
        cell = table[(i + 1, 5)]
        if quality > 0.95:
            cell.set_facecolor("#dcfce7")
        elif 0.70 <= quality <= 0.95:
            cell.set_facecolor("#fef3c7")
        else:
            cell.set_facecolor("#fee2e2")

    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def build_quantum_report(metrics_by_solver: Dict[str, Dict[str, Any]],
                         comparison_png: str, scaling_png: str,
                         output_path: str, title: str) -> None:
    """Build self-contained HTML report.

    Args:
        metrics_by_solver: Metrics grouped by solver (may include Phase 4 fields
            such as circuit_depth, backend_calibration_date, error_mitigation_method,
            and signal_quality when augmented from quantum_result files).
        comparison_png: Path to comparison equity curves PNG
        scaling_png: Path to scaling table PNG
        output_path: Path to write HTML
        title: Report title
    """
    def encode_png(path: str) -> str:
        try:
            with open(path, "rb") as fh:
                import base64
                b64 = base64.b64encode(fh.read()).decode("ascii")
                return f"data:image/png;base64,{b64}"
        except Exception:
            return ""

    equity_data_uri = encode_png(comparison_png)
    scaling_data_uri = encode_png(scaling_png)

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

    args = parser.parse_args()

    if not args.comparison_file.exists():
        raise FileNotFoundError(f"Comparison file not found: {args.comparison_file}")

    data = load_comparison_results(args.comparison_file)
    runs = data.get("runs", [])

    # Augment runs with Phase 4 fields from quantum_result_*.json files
    _try_augment_runs_from_result_files(runs, args.output_dir)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Generate plots
    equity_png = args.output_dir / "comparison_equity_curves.png"
    plot_comparison_equity_curves(runs, str(equity_png))

    scaling_png = args.output_dir / "scaling_table.png"
    plot_scaling_table(runs, str(scaling_png))

    # Group metrics by solver
    metrics_by_solver = {run["solver_name"]: run for run in runs}

    # Build HTML
    html_path = args.output_dir / "quantum_report.html"
    build_quantum_report(metrics_by_solver, str(equity_png), str(scaling_png), str(html_path), args.title)

    print(f"Wrote report to: {html_path}")


if __name__ == "__main__":
    main()
