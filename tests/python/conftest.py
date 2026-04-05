import io
import json
import pandas as pd
import numpy as np
import pytest
from pathlib import Path


@pytest.fixture
def timeseries_df() -> pd.DataFrame:
    """Return a simple synthetic timeseries DataFrame.

    Columns: date (YYYY-MM-DD), nav, return, drawdown
    """
    dates = pd.date_range("2023-01-01", periods=10, freq="D")
    rets = np.array([0.0, 0.01, -0.005, 0.02, -0.01, 0.0, 0.005, -0.02, 0.01, 0.0])
    nav = (1 + rets).cumprod()
    running_max = np.maximum.accumulate(nav)
    drawdown = (nav - running_max) / running_max
    df = pd.DataFrame({
        "date": dates.strftime("%Y-%m-%d"),
        "nav": nav,
        "return": rets,
        "drawdown": drawdown,
    })
    return df


@pytest.fixture
def flat_nav_df() -> pd.DataFrame:
    """Return a timeseries DataFrame with flat NAV (no drawdown)."""
    dates = pd.date_range("2023-01-01", periods=5, freq="D")
    nav = np.ones(5)
    rets = np.zeros(5)
    drawdown = np.zeros(5)
    return pd.DataFrame({"date": dates.strftime("%Y-%m-%d"), "nav": nav, "return": rets, "drawdown": drawdown})


@pytest.fixture
def single_row_df() -> pd.DataFrame:
    dates = pd.date_range("2023-01-01", periods=1, freq="D")
    return pd.DataFrame({"date": dates.strftime("%Y-%m-%d"), "nav": [1.0], "return": [0.0], "drawdown": [0.0]})


@pytest.fixture
def metrics_dict() -> dict:
    return {
        "total_return": 0.083,
        "annualized_return": 0.083,
        "annualized_volatility": 0.154,
        "sharpe_ratio": 0.54321,
        "sortino_ratio": 0.54321,
        "calmar_ratio": 0.6,
        "max_drawdown": -0.139,
        "var_95": -0.015,
        "cvar_95": -0.018,
        "skewness": 0.302,
        "kurtosis": 0.615,
    }


@pytest.fixture
def tmp_results_dir(tmp_path: Path) -> Path:
    d = tmp_path / "results"
    d.mkdir()
    return d
