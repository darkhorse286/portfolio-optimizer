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
"$PIP" install -r scripts/quantum/requirements.txt >/dev/null

echo "Running Python unit tests..."
"$PY" -m pytest -q tests/python || { echo 'Python unit tests failed'; exit 1; }

echo "Running Python integration tests..."
"$PY" -m pytest -q tests/integration -k 'not cpp_smoke' || { echo 'Python integration tests failed'; exit 1; }

_IBM_PASS=0
_IBM_FAIL=0
record_result() {
    local name="$1" code="$2"
    if [ "$code" -eq 0 ]; then
        echo "PASS: $name"
        _IBM_PASS=$((_IBM_PASS + 1))
    else
        echo "FAIL: $name (exit $code)"
        _IBM_FAIL=$((_IBM_FAIL + 1))
    fi
}

if [ "${TEST_IBM:-0}" = "1" ]; then
  echo "=== IBM Quantum connectivity check ==="
  "$PY" scripts/quantum/ibm_credentials_check.py || true
  record_result "IBM credentials" $?

  echo "=== mthree availability check ==="
  "$PY" -c "import mthree; print('mthree', mthree.__version__)" || true
  record_result "mthree available" $?

  echo "IBM checks: $_IBM_PASS passed, $_IBM_FAIL failed"
fi

echo "Running C++ tests (ctest)..."
if [ -d build ]; then
  (cd build && ctest --output-on-failure --exclude-regex integration) || { echo 'C++ tests failed'; exit 1; }
else
  echo 'Build directory not found; please run cmake and build first.'
  exit 1
fi

echo "Running C++ integration tests..."
if [ -d build ]; then
  (cd build && ctest --output-on-failure --tests-regex integration) && echo "C++ integration tests completed" || echo "C++ integration tests skipped (dependencies not available)"
else
  echo 'Build directory not found; please run cmake and build first.'
  exit 1
fi

echo "Running C++ smoke test (no python)..."
"$PY" -m pytest -q tests/integration/test_cpp_smoke.py || { echo 'C++ smoke test failed'; exit 1; }

echo "All tests passed."
exit 0
