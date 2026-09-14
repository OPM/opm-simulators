#!/bin/bash
# Shell wrapper for the unified State Rollback & Timestep Replay test suite
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec python3 "${SCRIPT_DIR}/run-all-state-rollback-tests.py" "$@"
