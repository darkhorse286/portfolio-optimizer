import subprocess
import json
from pathlib import Path
import pytest


def test_end_to_end_report(tmp_path: Path):
    results = tmp_path / 'results'
    results.mkdir()

    # Create analytics files (small)
    ts = results / 'analytics_timeseries.csv'
    ts.write_text('date,nav,return,drawdown\n2023-01-01,1.0,0.0,0.0\n')

    metrics = results / 'analytics_metrics.json'
    metrics.write_text(json.dumps({'total_return': 0.05, 'annualized_return': 0.05}))

    # Run the report generator
    # Prefer the project's virtualenv python if available (created by run_all_tests.sh)
    venv_py = Path.cwd() / '.venv' / 'bin' / 'python'
    python_exec = str(venv_py) if venv_py.exists() else 'python3'
    cmd = [python_exec, 'scripts/generate_report.py', '--input-dir', str(results), '--output-dir', str(results), '--title', 'IT Report']
    proc = subprocess.run(cmd, cwd=Path.cwd(), capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr

    report = results / 'report.html'
    assert report.exists()
    content = report.read_text()
    assert 'IT Report' in content
    assert 'Total return' in content
