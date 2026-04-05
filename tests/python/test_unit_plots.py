import io
import builtins
import matplotlib
import pytest

from scripts.generate_report import plot_equity_curve, plot_drawdown


def _capture_fig_save(monkeypatch):
    """Monkeypatch Figure.savefig to capture bytes into a buffer."""
    buf = io.BytesIO()

    def fake_save(self, *args, **kwargs):
        # write a small PNG header to indicate content was produced
        buf.write(b"\x89PNG\r\n\x1a\n")

    monkeypatch.setattr(matplotlib.figure.Figure, "savefig", fake_save)
    return buf


def test_plot_equity_curve_basic(timeseries_df, monkeypatch):
    buf = _capture_fig_save(monkeypatch)
    # Should not raise
    plot_equity_curve(timeseries_df, "/nonexistent/path/equity.png", title="Equity Curve")
    # Ensure buffer got written to by fake_save
    assert buf.getvalue().startswith(b"\x89PNG")


def test_plot_equity_curve_flat_nav(flat_nav_df, monkeypatch):
    buf = _capture_fig_save(monkeypatch)
    # Flat NAV should not error (no drawdown trough)
    plot_equity_curve(flat_nav_df, "/nonexistent/path/equity.png")
    assert buf.getvalue().startswith(b"\x89PNG")


def test_plot_drawdown_basic(timeseries_df, monkeypatch):
    buf = _capture_fig_save(monkeypatch)
    plot_drawdown(timeseries_df, "/nonexistent/path/drawdown.png")
    assert buf.getvalue().startswith(b"\x89PNG")


def test_plot_drawdown_single_row(single_row_df, monkeypatch):
    buf = _capture_fig_save(monkeypatch)
    # Single row should not crash
    plot_drawdown(single_row_df, "/nonexistent/path/drawdown.png")
    assert buf.getvalue().startswith(b"\x89PNG")
