#!/usr/bin/env python3
"""IBM Quantum credentials validation script."""

from __future__ import annotations

import os
import sys
from typing import Any, Dict


def check_token() -> bool:
    """Check if IBM_QUANTUM_TOKEN is set."""
    token = os.environ.get("IBM_QUANTUM_TOKEN")
    if token:
        print("PASS: IBM_QUANTUM_TOKEN is set.")
        return True
    else:
        print("FAIL: IBM_QUANTUM_TOKEN environment variable is not set.", file=sys.stderr)
        print("      Please set it with: export IBM_QUANTUM_TOKEN=your_token_here", file=sys.stderr)
        return False


def check_instance() -> bool:
    """Check if IBM_QUANTUM_INSTANCE is set or defaults."""
    instance = os.environ.get("IBM_QUANTUM_INSTANCE", "ibm-q/open/main")
    print(f"PASS: IBM_QUANTUM_INSTANCE is set to '{instance}'.")
    return True


def check_authentication() -> bool:
    """Attempt authentication with IBM Quantum."""
    token = os.environ.get("IBM_QUANTUM_TOKEN")
    instance = os.environ.get("IBM_QUANTUM_INSTANCE", "ibm-q/open/main")
    if not token:
        return False

    try:
        from qiskit_ibm_runtime import QiskitRuntimeService
        service = QiskitRuntimeService(channel="ibm_quantum", token=token, instance=instance)
        backends = service.backends()
        if backends:
            print("PASS: Successfully authenticated with IBM Quantum.")
            return True
        else:
            print("FAIL: Authentication succeeded but no backends available.", file=sys.stderr)
            return False
    except Exception as exc:
        print(f"FAIL: IBM Quantum authentication failed: {exc}", file=sys.stderr)
        return False


def list_backends() -> None:
    """List first 5 available backends."""
    token = os.environ.get("IBM_QUANTUM_TOKEN")
    instance = os.environ.get("IBM_QUANTUM_INSTANCE", "ibm-q/open/main")
    if not token:
        return

    try:
        from qiskit_ibm_runtime import QiskitRuntimeService
        service = QiskitRuntimeService(channel="ibm_quantum", token=token, instance=instance)
        backends = service.backends()[:5]
        print("Available backends (first 5):")
        for backend in backends:
            name = backend.name
            num_qubits = backend.num_qubits
            status = backend.status()
            operational = status.operational
            print(f"  {name}: {num_qubits} qubits, operational={operational}")
    except Exception as exc:
        print(f"Error listing backends: {exc}", file=sys.stderr)


def main() -> int:
    """Main entry point."""
    print("=== IBM Quantum Credentials Check ===")

    checks = [check_token, check_instance, check_authentication]
    results = [check() for check in checks]

    if all(results):
        print("\nAll checks passed! IBM Quantum is ready.")
        list_backends()
        return 0
    else:
        print("\nSome checks failed. Please resolve the issues above.", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())