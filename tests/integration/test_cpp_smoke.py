import subprocess
import os
from pathlib import Path


def test_cpp_smoke_no_python():
    # Run the compiled binary with an empty PATH so python3 is unavailable
    bin_path = Path('build/bin/portfolio_optimizer')
    assert bin_path.exists(), 'Compiled binary not found; build project first.'

    env = os.environ.copy()
    env['PATH'] = ''

    cmd = [str(bin_path), '--config', 'data/config/portfolio_config.json', '--report']
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True)

    # The program should exit cleanly (0) even if python3 isn't available
    assert proc.returncode == 0, f"Binary exited with {proc.returncode}: {proc.stderr}"
