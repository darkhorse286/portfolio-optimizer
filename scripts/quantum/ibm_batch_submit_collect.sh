#!/usr/bin/env bash
# Submit 3 rounds of QAOA, QAMO, and QAMOO to IBM hardware, collect all
# results, then regenerate the benchmark report.
#
# Usage:
#   bash scripts/quantum/ibm_batch_submit_collect.sh
#
# Environment variables required:
#   IBM_QUANTUM_TOKEN    — IBM Quantum API token
#   IBM_QUANTUM_INSTANCE — IBM Quantum instance (optional, uses default if unset)

set -euo pipefail

# ── Configuration ────────────────────────────────────────────────────────────
ROUNDS=3
BACKEND="ibm_fez"
SHOTS=1024
N_FRONTIER=20          # lambda sweep points for QAMOO
TIMEOUT_MIN=120        # per-algorithm collect timeout
POLL_INTERVAL=30       # seconds between IBM status checks
RESULTS_DIR="results"

QAOA_PROBLEM="$RESULTS_DIR/quantum_problem_qaoa.json"
QAMO_PROBLEM="$RESULTS_DIR/quantum_problem_qamo.json"

# Dedicated jobs files for this batch (separate from Aer jobs files)
QAOA_JOBS="$RESULTS_DIR/ibm_batch_qaoa.json"
QAMO_JOBS="$RESULTS_DIR/ibm_batch_qamo.json"
QAMOO_JOBS="$RESULTS_DIR/ibm_batch_qamoo.json"

# ── Pre-flight checks ────────────────────────────────────────────────────────
if [[ -z "${IBM_QUANTUM_TOKEN:-}" ]]; then
    echo "ERROR: IBM_QUANTUM_TOKEN is not set." >&2
    exit 1
fi

for f in "$QAOA_PROBLEM" "$QAMO_PROBLEM"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Problem file not found: $f" >&2
        echo "  Run './build/bin/portfolio_optimizer --config data/config/portfolio_config.json --output results --benchmark' first." >&2
        exit 1
    fi
done

# Start fresh — remove any leftover batch jobs files from a prior run
rm -f "$QAOA_JOBS" "$QAMO_JOBS" "$QAMOO_JOBS"

# ── Submit ───────────────────────────────────────────────────────────────────
echo "=== Submitting $ROUNDS rounds × 3 algorithms to $BACKEND ==="

for round in $(seq 1 $ROUNDS); do
    echo ""
    echo "--- Round $round / $ROUNDS ---"

    echo "[Round $round] QAOA"
    python3 scripts/quantum/qiskit_submit.py \
        --problem-file "$QAOA_PROBLEM" \
        --jobs-file    "$QAOA_JOBS"    \
        --backend      "$BACKEND"      \
        --mode         qaoa            \
        --shots        "$SHOTS"

    echo "[Round $round] QAMO"
    python3 scripts/quantum/qiskit_submit.py \
        --problem-file "$QAMO_PROBLEM" \
        --jobs-file    "$QAMO_JOBS"    \
        --backend      "$BACKEND"      \
        --mode         qamo            \
        --shots        "$SHOTS"

    echo "[Round $round] QAMOO ($N_FRONTIER frontier points)"
    python3 scripts/quantum/qiskit_submit.py \
        --problem-file      "$QAMO_PROBLEM"  \
        --jobs-file         "$QAMOO_JOBS"    \
        --backend           "$BACKEND"       \
        --mode              qamoo            \
        --shots             "$SHOTS"         \
        --n-frontier-points "$N_FRONTIER"
done

TOTAL=$((ROUNDS * 3))
echo ""
echo "=== $TOTAL jobs submitted. Polling IBM until complete... ==="
echo "    Timeout: ${TIMEOUT_MIN} min per algorithm  |  Poll interval: ${POLL_INTERVAL} s"

# ── Collect ──────────────────────────────────────────────────────────────────
collect() {
    local label="$1"
    local jobs_file="$2"
    echo ""
    echo "--- Collecting $label ---"
    python3 scripts/quantum/qiskit_collect.py \
        --jobs-file     "$jobs_file"    \
        --results-dir   "$RESULTS_DIR"  \
        --timeout-min   "$TIMEOUT_MIN"  \
        --poll-interval "$POLL_INTERVAL"
}

collect "QAOA"  "$QAOA_JOBS"
collect "QAMO"  "$QAMO_JOBS"
collect "QAMOO" "$QAMOO_JOBS"

# ── Report ───────────────────────────────────────────────────────────────────
echo ""
echo "=== Regenerating benchmark report ==="
python3 scripts/quantum/benchmark_viz.py --output-dir "$RESULTS_DIR"

echo ""
echo "=== Done. Report written to $RESULTS_DIR/quantum_report.html ==="
