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

        # Normalize to 1.0
        nav_values = [point["nav"] for point in nav_series]
        if nav_values:
            first_nav = nav_values[0]
            normalized = [v / first_nav for v in nav_values]
            dates = [point["date"] for point in nav_series]
            plt.plot(dates, normalized, label=solver_name, color=color_map.get(solver_type, "#6b7280"))

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
    fig, ax = plt.subplots(figsize=(10, 3))
    ax.axis('off')

    # Prepare data
    data = []
    for res in results:
        row = [
            res["solver_name"],
            res["solver_type"],
            res["execution_backend"],
            f"{res['performance']['sharpe_ratio']:.3f}",
            f"{res['performance']['total_return'] * 100:.1f}%",
            f"{res['solution_quality_vs_classical']:.2f}x",
            f"{res['solve_time_ms']:.1f}ms"
        ]
        data.append(row)

    columns = ["Solver", "Type", "Backend", "Sharpe", "Return", "vs Markowitz", "Solve ms"]

    table = ax.table(cellText=data, colLabels=columns, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.2)

    # Color vs Markowitz cells
    for i, res in enumerate(results):
        quality = res["solution_quality_vs_classical"]
        cell = table[(i+1, 5)]  # +1 for header
        if quality > 0.95:
            cell.set_facecolor("#dcfce7")  # green
        elif 0.70 <= quality <= 0.95:
            cell.set_facecolor("#fef3c7")  # amber
        else:
            cell.set_facecolor("#fee2e2")  # red

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
            solver,
            metrics.get("solver_type", ""),
            metrics.get("execution_backend", ""),
            f"{perf.get('sharpe_ratio', 0):.3f}",
            f"{perf.get('total_return', 0) * 100:.1f}%",
            f"{perf.get('annualized_volatility', 0) * 100:.1f}%",
            f"{perf.get('max_drawdown', 0) * 100:.1f}%"
        ])

    template_str = """<!doctype html>
<html>
<head>
    <meta charset="utf-8">
    <title>{{ title }}</title>
    <style>
        body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: #ffffff; color: #111827; margin: 20px; }
        h1 { font-size: 1.5rem; margin-bottom: 0.5rem; }
        .section { margin-bottom: 2rem; }
        .metrics { border-collapse: collapse; width: 100%; margin-bottom: 16px; }
        .metrics td, .metrics th { padding: 8px 10px; border-bottom: 1px solid #e5e7eb; }
        .metrics tr:nth-child(even) { background: #f9fafb; }
        .chart { margin-bottom: 18px; }
    </style>
</head>
<body>
    <h1>{{ title }}</h1>

    <div class="section">
        <h2>Solver Comparison</h2>
        {% if equity_data_uri %}
        <img alt="Equity Curves" src="{{ equity_data_uri }}" style="max-width: 100%;">
        {% else %}
        <div>Equity chart not available.</div>
        {% endif %}
    </div>

    <div class="section">
        <h2>Scaling Table</h2>
        {% if scaling_data_uri %}
        <img alt="Scaling Table" src="{{ scaling_data_uri }}" style="max-width: 100%;">
        {% else %}
        <div>Scaling table not available.</div>
        {% endif %}
    </div>

    <div class="section">
        <h2>Per-Solver Metrics</h2>
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