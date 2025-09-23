# Python scripts for testing the generated wheels

## Installation

Compatibility note (Windows checkout):
The directory names under `docker/test_wheels/dockerfiles` no longer use `:`.
Mappings:

- `debian:11` -> `debian-11`
- `ubuntu:22.04` -> `ubuntu-22.04`
- `ubuntu:24.04` -> `ubuntu-24.04`

The CLI accepts both forms for `--docker-os` (e.g., `ubuntu:22.04` or `ubuntu-22.04`) but
uses hyphenated names for filesystem paths and image tags (e.g., `test-opm-wheels-ubuntu-22.04`).

**Important:** This package requires access to the `python_versions.json` configuration file in the opm-simulators repository. Choose one of these installation methods:

### Option 1: Install within repository (Recommended)
```bash
# Create venv inside the opm-simulators repository
$ cd /path/to/opm-simulators/python/docker/test_wheels/scripts
$ python -m venv .venv
$ source .venv/bin/activate
$ pip install .    # Installs package and dependencies
```

### Option 2: Global installation with environment variable
```bash
# Install anywhere, but set environment variable
$ python -m venv .venv
$ source .venv/bin/activate
$ pip install .
$ export OPM_SIMULATORS_ROOT=/path/to/opm-simulators
```

### For development (editable install with locked dependency versions):
```bash
$ cd /path/to/opm-simulators/python/docker/test_wheels/scripts
$ uv venv
$ uv sync          # Installs package in editable mode with exact versions from uv.lock
$ source .venv/bin/activate
```

## Build the PyPI wheels

Build the wheels from the `python` folder of the `opm-simulators` repository using
the provided helper script for simplified building with automatic logging:

```bash
# Static build with wheel extraction
./docker/run-docker-build.sh static --output wheelhouse
# Creates: docker/wheelhouse/ with wheels

# Shared libraries build
./docker/run-docker-build.sh shared --libtype shared --output wheelhouse-shared
# Creates: docker/wheelhouse-shared/ with wheels
```
see `.docker/run-docker-build.sh` for more examples.

## Build a test container

Build a Docker container for testing the generated wheels:

```bash
# Basic usage (Ubuntu 24.04)
$ opm-wheels build-docker-image --docker-os="ubuntu:24.04"

# Build with custom repository and branch
$ opm-wheels build-docker-image \
    --docker-os="ubuntu:22.04" \
    --opm-simulators-branch="pr-branch" \
    --opm-common-branch="development" \
    --docker-progress

# Build for Debian
$ opm-wheels build-docker-image --docker-os="debian:11"
```

### Available options for `build-docker-image`:
- `--docker-os`, `-d` - OS for Docker image (default: ubuntu-22.04; legacy colon form accepted)
  - Supported: `ubuntu-22.04`, `ubuntu-24.04`, `debian-11`
- `--python-versions`, `-p` - Python versions to install in container
  - Supported: `3.8,3.9,3.10,3.11,3.12,3.13`
  - Format: comma-separated list without spaces (e.g., `3.11,3.12`)
  - Default: installs all supported versions
- `--opm-simulators-repo`, `-osr` - OPM simulators repository URL
- `--opm-common-repo`, `-ocr` - OPM common repository URL
- `--opm-simulators-branch`, `-osb` - Simulators branch (default: master)
- `--opm-common-branch`, `-ocb` - Common branch (default: master)
- `--docker-progress`, `-dp` - Show Docker build progress output

## Run tests

Run unit tests for `opm-common` and `opm-simulators` wheels in the test container:

```bash
# Test all Python versions with default wheel directory
$ opm-wheels run-docker-image --docker-os="ubuntu:24.04"

# Test specific Python versions
$ opm-wheels run-docker-image \
    --docker-os="ubuntu:22.04" \
    --python-versions="3.11,3.12,3.13"

# Test with custom wheel directory
$ opm-wheels run-docker-image \
    --docker-os="debian:11" \
    --wheel-dir="/path/to/wheelhouse-shared" \
    --python-versions="3.12"

# Test using host test directories (CI/CD friendly)
$ opm-wheels run-docker-image \
    --docker-os="ubuntu:24.04" \
    --wheel-dir="/path/to/wheelhouse" \
    --host-tests-dir="/path/to/opm-repos"
```

### Available options for `run-docker-image`:
- `--docker-os`, `-d` - OS for Docker image (default: ubuntu-22.04; legacy colon form accepted)
- `--wheel-dir` - Directory containing wheel files (default: python/wheelhouse)
- `--python-versions`, `-p` - Python versions for testing
  - Supported: `3.8,3.9,3.10,3.11,3.12,3.13`
  - Format: comma-separated list without spaces (e.g., `3.11,3.12`)
- `--host-tests-dir` - Use test directories from host instead of cloned repos
  - Expects: `host-tests-dir/opm-common/python` and `host-tests-dir/opm-simulators/python`
  - Eliminates need for repository cloning and nightly container rebuilds
- `--test-cases-common` - Run specific opm-common test cases
  - Format: comma-separated test patterns (e.g., `eclfile,esmry`)
- `--test-cases-simulators` - Run specific opm-simulators test cases
  - Format: comma-separated test patterns (e.g., `basic,mpi`)
- `--stop-on-error` - Stop execution on first test failure
  - Default: continue running all tests even if some fail

## CI/CD Usage

### Jenkins Pipeline Example

For Jenkins environments with network restrictions, use `--host-tests-dir` to avoid repository cloning:

```bash
#!/bin/bash
# Jenkins script - assumes workspace contains checked-out opm repositories

# Build test container once (can be cached/reused)
opm-wheels build-docker-image --docker-os="ubuntu:24.04"

# Run tests using host directories (no network access needed)
opm-wheels run-docker-image \
    --docker-os="ubuntu:24.04" \
    --wheel-dir="$WORKSPACE/wheels" \
    --host-tests-dir="$WORKSPACE" \
    --python-versions="3.11,3.12,3.13"
```

### Directory Structure for Host Tests

When using `--host-tests-dir`, ensure your directory structure matches:

```
your-repos-dir/
├── opm-common/
│   └── python/           # Contains Python test files
│       ├── tests/
│       └── test_*.py
└── opm-simulators/
    └── python/           # Contains Python test files
        ├── tests/
        └── test_*.py
```

### Benefits for CI/CD

- **No nightly rebuilds**: Container image can be built once and reused
- **Network isolation**: No git clones during test execution
- **Version consistency**: Tests run against exact workspace codebase
- **Faster execution**: Skip repository cloning step
- **Debugging friendly**: Test failures map directly to workspace files

## Debugging

The `opm-wheels` tool provides several options specifically designed to streamline debugging workflows:

### Debugging Features Overview

1. **Python Version Selection** - Build and test only specific Python versions
2. **Focused Test Selection** - Run specific test cases instead of entire test suites
3. **Continue-on-Error Mode** - Run all tests to see all failures at once
4. **Host Test Directories** - Test against local code without container rebuilds

### Python Version Selection

Speed up debugging iterations by building/testing only the Python versions you need:

```bash
# Build container with only Python 3.12 (much faster than all versions)
$ opm-wheels build-docker-image --docker-os="ubuntu-24.04" --python-versions="3.12"

# Test only Python 3.11 and 3.12
$ opm-wheels run-docker-image --docker-os="ubuntu-24.04" --python-versions="3.11,3.12"

# Debug single version in detail
$ opm-wheels run-docker-image --docker-os="ubuntu-24.04" --python-versions="3.13"
```

**Benefits:**
- Faster Docker builds when changing entrypoint.py or test configuration
- Reduced test execution time during development
- Focus debugging on specific Python versions showing issues

### Focused Test Selection

Target specific test cases to isolate issues:

```bash
# Run only specific opm-common tests
$ opm-wheels run-docker-image \
    --docker-os="ubuntu-24.04" \
    --python-versions="3.12" \
    --test-cases-common="eclfile,esmry"

# Run only specific opm-simulators tests
$ opm-wheels run-docker-image \
    --docker-os="ubuntu-24.04" \
    --python-versions="3.12" \
    --test-cases-simulators="basic,mpi"

# Combine both for targeted debugging
$ opm-wheels run-docker-image \
    --docker-os="ubuntu-24.04" \
    --python-versions="3.12" \
    --test-cases-common="eclfile" \
    --test-cases-simulators="basic"
```

**Test Case Patterns:**
- opm-common: `eclfile`, `esmry`, `restart_runs_ext` (runs `test_<pattern>.py`)
- opm-simulators: `basic`, `mpi`, `fluidstate_variables` (runs `test/test_<pattern>.py`)

### Continue-on-Error Mode

By default, tests continue running even if some fail, showing all issues at once:

```bash
# Default: Continue running all tests (recommended for debugging)
$ opm-wheels run-docker-image --docker-os="ubuntu-24.04" --python-versions="3.12"
# Output shows: opm-simulators failures + opm-common results

# Stop on first error (legacy behavior)
$ opm-wheels run-docker-image \
    --docker-os="ubuntu-24.04" \
    --python-versions="3.12" \
    --stop-on-error
# Output: stops after first test suite failure
```

**Default Behavior:**
- ✅ Run opm-simulators tests (may fail)
- ✅ Run opm-common tests regardless of previous failures
- ✅ Show comprehensive failure summary
- ✅ Continue with additional Python versions if specified

**With `--stop-on-error`:**
- ❌ Stop at first test failure
- ❌ No comprehensive failure view
- ❌ Must fix issues one by one

### Complete Debugging Workflow

Here's a recommended debugging workflow combining all features:

```bash
# 1. Build lightweight container for single Python version
$ opm-wheels build-docker-image --docker-os="ubuntu-24.04" --python-versions="3.12"

# 2. Get comprehensive test results without stopping on errors
$ opm-wheels run-docker-image \
    --docker-os="ubuntu-24.04" \
    --wheel-dir="wheelhouse-static2" \
    --python-versions="3.12" \
    --host-tests-dir="/path/to/opm-repos"

# 3. Focus on specific failing tests
$ opm-wheels run-docker-image \
    --docker-os="ubuntu-24.04" \
    --python-versions="3.12" \
    --test-cases-simulators="fluidstate_variables" \
    --host-tests-dir="/path/to/opm-repos"

# 4. Verify fixes with full test suite
$ opm-wheels run-docker-image \
    --docker-os="ubuntu-24.04" \
    --python-versions="3.12" \
    --host-tests-dir="/path/to/opm-repos"
```

## Additional commands

### Show wheel contents
```bash
# List files inside a wheel
$ opm-wheels show-wheel-files opm-2025.10-cp312-cp312-manylinux_2_27_x86_64.whl

# Use custom wheel directory
$ opm-wheels show-wheel-files \
    --wheel-dir="wheelhouse-shared" \
    opm_simulators-2025.10-cp312-cp312-manylinux_2_27_x86_64.whl
```

## Troubleshooting

### FileNotFoundError: Could not locate python_versions.json

This error occurs when the CLI cannot find the configuration file. Solutions:

1. **Run from within the repository:**
   ```bash
   cd /path/to/opm-simulators
   cd python/docker/test_wheels/scripts
   source .venv/bin/activate
   opm-wheels run-docker-image --help
   ```

2. **Set environment variable for global installations:**
   ```bash
   export OPM_SIMULATORS_ROOT=/path/to/opm-simulators
   opm-wheels run-docker-image --help
   ```

3. **Verify the configuration file exists:**
   ```bash
   ls /path/to/opm-simulators/python/docker/python_versions.json
   ```

### RuntimeError: Not a valid Git repository

This error means you're not inside a Git repository. Either:
- Navigate to the opm-simulators repository directory, or
- Use the `OPM_SIMULATORS_ROOT` environment variable method

### FileNotFoundError: Expected opm-*/python directory not found

When using `--host-tests-dir`, this error means the directory structure is incorrect. Ensure your directory contains:
```
host-tests-dir/
├── opm-common/python/     # Must exist with test files
└── omp-simulators/python/ # Must exist with test files
```
