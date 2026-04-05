#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")" && pwd)
VENV="$ROOT/.venv"
PY="$VENV/bin/python"
PIP="$VENV/bin/pip"
PYTEST="$VENV/bin/pytest"

echo "Ensuring virtualenv exists at $VENV"
if [ ! -x "$PY" ]; then
  python3 -m venv "$VENV"
fi

echo "Installing test dependencies into virtualenv"
"$PIP" install --upgrade pip setuptools wheel >/dev/null
"$PIP" install --upgrade pytest >/dev/null

echo "Running Python unit tests..."
"$PY" -m pytest -q tests/python || { echo 'Python unit tests failed'; exit 1; }

echo "Running Python integration tests..."
"$PY" -m pytest -q tests/integration -k 'not cpp_smoke' || { echo 'Python integration tests failed'; exit 1; }

echo "Running C++ tests (ctest)..."
if [ -d build ]; then
  (cd build && ctest --output-on-failure) || { echo 'C++ tests failed'; exit 1; }
else
  echo 'Build directory not found; please run cmake and build first.'
  exit 1
fi

echo "Running C++ smoke test (no python)..."
"$PY" -m pytest -q tests/integration/test_cpp_smoke.py || { echo 'C++ smoke test failed'; exit 1; }

echo "All tests passed."
exit 0
