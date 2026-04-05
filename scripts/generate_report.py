#!/usr/bin/env python3
"""Generate an HTML report from analytics output.

This script loads time series and metrics produced by the C++ backtest
and wires them into plotting and report-building functions (stubs).

Dependencies: pandas, matplotlib, numpy, jinja2
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import base64
from jinja2 import Template

import pandas as pd


def plot_equity_curve(df: pd.DataFrame, output_path: str, title: str = "Equity Curve") -> None:
    """Plot the equity (NAV) curve and save to `output_path`.

    The function normalizes NAV to start at 1.0, shades underwater periods,
    annotates the worst drawdown trough, and saves a PNG image.

    Args:
        df: DataFrame with columns [`date`, `nav`, `return`, `drawdown`].
        output_path: File path where the PNG will be saved.
        title: Title for the plot.

    Returns:
        None
    """
    # Make a defensive copy
    data = df.copy()

    # Ensure date column is datetime
    if "date" in data.columns:
        data["date"] = pd.to_datetime(data["date"])  # type: ignore[index]
    else:
        raise ValueError("Input DataFrame must contain a 'date' column")

    if "nav" not in data.columns:
        raise ValueError("Input DataFrame must contain a 'nav' column")

    # Sort by date
    data = data.sort_values(by="date")

    # Normalize NAV to start at 1.0
    first_nav = float(data["nav"].iat[0])
    if first_nav == 0:
        raise ValueError("First NAV value is zero, cannot normalize")
    nav_norm = data["nav"] / first_nav

    dates = data["date"].to_numpy()

    # Running peak
    running_peak = np.maximum.accumulate(nav_norm.to_numpy())

    # Prepare figure
    fig, ax = plt.subplots(figsize=(10, 4))

    # Plot NAV
    ax.plot(dates, nav_norm, color="#2563eb", linewidth=1.5)

    # Horizontal reference at 1.0
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8)

    # Shade underwater areas where nav_norm < running_peak
    underwater = nav_norm.to_numpy() < running_peak
    if np.any(underwater):
        ax.fill_between(dates, nav_norm, running_peak, where=underwater,
                        facecolor="#ef4444", alpha=0.15, interpolate=True)

    # Annotate worst drawdown trough
    if "drawdown" in data.columns:
        dd = data["drawdown"].to_numpy()
        # choose worst as min if negative values present, otherwise max
        if np.nanmin(dd) < 0:
            worst_idx = int(np.nanargmin(dd))
        else:
            worst_idx = int(np.nanargmax(dd))

        worst_date = dates[worst_idx]
        worst_nav = nav_norm.iat[worst_idx]
        worst_dd_val = float(dd[worst_idx])

        # Interpret drawdown as fraction if reasonable
        if abs(worst_dd_val) <= 1.5:
            pct = worst_dd_val * 100.0
        else:
            pct = worst_dd_val

        label = f"{pct:.1f}%"

        # Marker and annotation
        ax.plot(worst_date, worst_nav, marker="v", color="#ef4444", markersize=8)
        ax.annotate(label,
                    xy=(worst_date, worst_nav),
                    xytext=(0, -12),
                    textcoords="offset points",
                    ha="center",
                    va="top",
                    color="#ef4444",
                    fontsize=9)

    # X-axis formatting
    locator = mdates.AutoDateLocator(minticks=3, maxticks=8)
    formatter = mdates.DateFormatter("%b %Y")
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    fig.autofmt_xdate(rotation=30)

    # Y-axis as multipliers
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:.1f}x"))

    # Grid on y-axis only
    ax.grid(axis="y", color="lightgray", linestyle="-", linewidth=0.7)

    # Title
    ax.set_title(title)

    # No legend by default

    # Save figure
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_drawdown(df: pd.DataFrame, output_path: str) -> None:
    """Plot drawdown series and save to `output_path`.

    Args:
        df: DataFrame with columns [`date`, `nav`, `return`, `drawdown`].
        output_path: File path where the PNG will be saved.

    Returns:
        None
    """
    data = df.copy()

    # Ensure date column is datetime
    if "date" in data.columns:
        data["date"] = pd.to_datetime(data["date"])  # type: ignore[index]
    else:
        raise ValueError("Input DataFrame must contain a 'date' column")

    if "drawdown" not in data.columns:
        raise ValueError("Input DataFrame must contain a 'drawdown' column")

    data = data.sort_values(by="date")
    dates = data["date"].to_numpy()
    dd = data["drawdown"].to_numpy()

    # Identify maximum drawdown (most negative value)
    if np.all(np.isnan(dd)):
        # Nothing to plot
        return

    worst_idx = int(np.nanargmin(dd))
    worst_date = dates[worst_idx]
    worst_val = float(dd[worst_idx])

    # Format label as percentage (e.g. -18.3%)
    label = f"{worst_val * 100:.1f}%"

    fig, ax = plt.subplots(figsize=(10, 3))

    # Fill area between drawdown and 0
    ax.fill_between(dates, dd, 0.0, where=~np.isnan(dd), facecolor="#dc2626", alpha=0.5)

    # Draw drawdown line
    ax.plot(dates, dd, color="#dc2626", linewidth=1)

    # Horizontal line at y=0
    ax.axhline(0.0, color="gray", linewidth=0.5)

    # Mark maximum drawdown point
    ax.plot(worst_date, worst_val, marker="o", color="#dc2626", markersize=6)
    ax.annotate(label,
                xy=(worst_date, worst_val),
                xytext=(0, 8),
                textcoords="offset points",
                ha="center",
                va="bottom",
                color="#dc2626",
                fontsize=9)

    # X-axis formatting
    locator = mdates.AutoDateLocator(minticks=3, maxticks=8)
    formatter = mdates.DateFormatter("%b %Y")
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    fig.autofmt_xdate(rotation=30)

    # Y-axis as percentages
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y * 100:.0f}%"))

    # Grid on y-axis only
    ax.grid(axis="y", color="lightgray", linestyle="-", linewidth=0.7)

    ax.set_title("Drawdown")

    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def build_html_report(metrics: Dict[str, float], equity_png: str, drawdown_png: str, output_path: str, title: str) -> None:
        """Build a self-contained HTML report embedding metrics and PNG charts.

        Args:
                metrics: Dictionary of analytics metrics.
                equity_png: Path to the equity curve PNG file.
                drawdown_png: Path to the drawdown PNG file.
                output_path: Path where the HTML file will be written.
                title: Title to use in the HTML report.

        Returns:
                None. Writes a self-contained HTML file to `output_path`.
        """

        def fmt_pct(key: str, decimals: int = 1) -> str:
                val = metrics.get(key)
                if val is None:
                        return "N/A"
                try:
                        return f"{float(val) * 100:.{decimals}f}%"
                except Exception:
                        return "N/A"

        def fmt_float(key: str, decimals: int) -> str:
                val = metrics.get(key)
                if val is None:
                        return "N/A"
                try:
                        return f"{float(val):.{decimals}f}"
                except Exception:
                        return "N/A"

        # Prepare excess kurtosis: prefer explicit key, otherwise derive from 'kurtosis'
        if "excess_kurtosis" in metrics:
                kurt_key = "excess_kurtosis"
        else:
                kurt_key = "kurtosis"

        # Build rows in the required order
        rows = [
                ("Total return", fmt_pct("total_return", 1)),
                ("Annualized return", fmt_pct("annualized_return", 1)),
                ("Annualized volatility", fmt_pct("annualized_volatility", 1)),
                ("Sharpe ratio", fmt_float("sharpe_ratio", 2)),
                ("Sortino ratio", fmt_float("sortino_ratio", 2)),
                ("Calmar ratio", fmt_float("calmar_ratio", 2)),
                ("Max drawdown", fmt_pct("max_drawdown", 1)),
                ("VaR 95%", fmt_pct("var_95", 1)),
                ("CVaR 95%", fmt_pct("cvar_95", 1)),
                ("Skewness", fmt_float("skewness", 3)),
                ("Excess kurtosis", fmt_float(kurt_key, 3)),
        ]

        # Load and base64-encode images (if available)
        def encode_png(path: str) -> str:
                try:
                        with open(path, "rb") as fh:
                                b64 = base64.b64encode(fh.read()).decode("ascii")
                        return f"data:image/png;base64,{b64}"
                except Exception:
                        return ""

        equity_data_uri = encode_png(equity_png)
        drawdown_data_uri = encode_png(drawdown_png)

        template_str = """<!doctype html>
<html>
    <head>
        <meta charset="utf-8">
        <title>{{ title }}</title>
        <style>
            body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; background: #ffffff; color: #111827; margin: 20px; }
            h1 { font-size: 1.5rem; margin-bottom: 0.5rem; }
            .metrics { border-collapse: collapse; width: 360px; margin-bottom: 16px; }
            .metrics td, .metrics th { padding: 8px 10px; border-bottom: 1px solid #e5e7eb; }
            .metrics tr:nth-child(even) { background: #f9fafb; }
            .chart { margin-bottom: 18px; }
            .charts { display: flex; flex-direction: column; gap: 12px; }
            .center { display: flex; gap: 20px; align-items: center; }
            .metric-label { color: #374151; }
            .metric-value { font-weight: 600; }
            @media (min-width: 900px) { .charts { flex-direction: row; } .chart img { max-width: 100%; height: auto; } }
        </style>
    </head>
    <body>
        <h1>{{ title }}</h1>
        <div class="center">
            <table class="metrics">
                <tbody>
                {% for label, value in rows %}
                    <tr>
                        <td class="metric-label">{{ label }}</td>
                        <td class="metric-value">{{ value }}</td>
                    </tr>
                {% endfor %}
                </tbody>
            </table>
        </div>

        <div class="charts">
            <div class="chart">
                {% if equity_data_uri %}
                <img alt="Equity Curve" src="{{ equity_data_uri }}">
                {% else %}
                <div>Equity chart not available.</div>
                {% endif %}
            </div>

            <div class="chart">
                {% if drawdown_data_uri %}
                <img alt="Drawdown" src="{{ drawdown_data_uri }}">
                {% else %}
                <div>Drawdown chart not available.</div>
                {% endif %}
            </div>
        </div>
    </body>
</html>
"""

        tmpl = Template(template_str)
        html = tmpl.render(title=title, rows=rows, equity_data_uri=equity_data_uri, drawdown_data_uri=drawdown_data_uri)

        # Write to output_path
        with open(output_path, "w", encoding="utf-8") as fh:
                fh.write(html)



def _validate_input_files(input_dir: Path) -> None:
    """Validate that required input files exist in `input_dir`.

    Args:
        input_dir: Directory expected to contain analytics files.

    Raises:
        FileNotFoundError: If any required file is missing.
    """
    csv_path = input_dir / "analytics_timeseries.csv"
    json_path = input_dir / "analytics_metrics.json"

    missing = []
    if not csv_path.exists():
        missing.append(str(csv_path))
    if not json_path.exists():
        missing.append(str(json_path))

    if missing:
        raise FileNotFoundError(
            "Required analytics file(s) not found: " + ", ".join(missing)
        )


def load_inputs(input_dir: Path) -> (pd.DataFrame, Dict[str, float]):
    """Load the analytics time series CSV and metrics JSON.

    Args:
        input_dir: Directory where `analytics_timeseries.csv` and
            `analytics_metrics.json` are located.

    Returns:
        A tuple of (timeseries DataFrame, metrics dict).

    Raises:
        FileNotFoundError: If input files are missing.
    """
    _validate_input_files(input_dir)

    csv_path = input_dir / "analytics_timeseries.csv"
    json_path = input_dir / "analytics_metrics.json"

    # Load CSV (parse `date` if present)
    try:
        ts = pd.read_csv(csv_path, parse_dates=["date"])  # type: pd.DataFrame
    except Exception:
        # Fallback: load without parsing dates
        ts = pd.read_csv(csv_path)

    # Load JSON metrics
    with open(json_path, "r", encoding="utf-8") as fh:
        metrics = json.load(fh)

    return ts, metrics


def main() -> None:
    """CLI entry point for report generation.

    Parses command-line arguments, loads inputs, calls plotting/report
    functions, and writes an HTML placeholder report to the output
    directory as `report.html`.
    """
    parser = argparse.ArgumentParser(description="Generate portfolio HTML report")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("results"),
        help="Directory containing analytics_timeseries.csv and analytics_metrics.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
        help="Directory to write report and plots",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="Portfolio Report",
        help="Title to show in the generated report",
    )

    args = parser.parse_args()
    input_dir: Path = args.input_dir
    output_dir: Path = args.output_dir
    title: str = args.title

    # Ensure output dir exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load inputs (will raise FileNotFoundError with a clear message if missing)
    timeseries, metrics = load_inputs(input_dir)

    # Call plotting stubs (equity curve and drawdown write PNGs to output_dir)
    plot_equity_curve(timeseries, str(output_dir / "equity_curve.png"), title="Equity Curve")
    plot_drawdown(timeseries, str(output_dir / "drawdown.png"))

    # Build and write HTML report, embedding generated PNGs
    out_path = output_dir / "report.html"
    build_html_report(metrics, str(output_dir / "equity_curve.png"), str(output_dir / "drawdown.png"), str(out_path), title)
    print(f"Wrote report to: {out_path}")


if __name__ == "__main__":
    main()
