#!/usr/bin/env python3
"""Benchmark visualization script for quantum portfolio optimizers."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Any

import matplotlib.pyplot as plt
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
    plt.ylabel("Normalized NAV")
    plt.title("Solver Comparison - Equity Curves")
    plt.legend()
    plt.grid(True, alpha=0.3)
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
    fig, ax = plt.subplots(figsize=(14, fig_height))
    ax.axis('off')

    data = []
    for res in results:
        row = [
            _display_name(res["solver_name"]),
            _display_type(res["solver_type"]),
            res["execution_backend"],
            f"{res['performance']['sharpe_ratio']:.3f}",
            f"{res['performance']['total_return'] * 100:.1f}%",
            f"{res['solution_quality_vs_classical']:.2f}x",
            f"{res['solve_time_ms']:.1f}ms",
        ]
        data.append(row)

    columns = ["Solver", "Type", "Backend", "Sharpe", "Return", "vs Markowitz", "Solve ms"]

    table = ax.table(cellText=data, colLabels=columns, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    # Let matplotlib size columns to content, then add a little vertical padding
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
        metrics_by_solver: Metrics grouped by solver
        comparison_png: Path to comparison equity curves PNG
        scaling_png: Path to scaling table PNG
        output_path: Path to write HTML
        title: Report title
    """
    # Encode images
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

    # Build metrics table
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
            f"{perf.get('max_drawdown', 0) * 100:.1f}%"  # already negative
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
        .chart-wrap {
            overflow-x: auto;
        }
        .chart-wrap img {
            display: block;
            max-width: 100%;
            height: auto;
        }
        .table-wrap {
            overflow-x: auto;
        }
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
</div>
</body>
</html>"""

    tmpl = Template(template_str)
    html = tmpl.render(title=title, rows=rows, equity_data_uri=equity_data_uri, scaling_data_uri=scaling_data_uri)

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