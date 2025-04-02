#! /usr/bin/env python3

# This script is called as ENTRYPOINT from the Dockerfiles, e.g. docker/test_wheels/dockerfiles/.../Dockerfile
# For each python version in PYTHON_VERSIONS (with the values in docker/python_versions.json as default), it:
#  - installs the wheels for the python version
#  - runs the opm-common and opm-simulators python unit tests
# Returns 0 if all tests passed, 1 otherwise

import glob
import json
import os
import subprocess
import sys

def run_opm_common_tests(python_bin: str) -> None:
    """Run all unittests in opm-common/python/tests"""
    subprocess.run(
        [
            python_bin, "-m", "unittest", "discover",
            "-s", "/test/opm-common/python",  # path to the tests
            "-v"  # verbose output
        ],
        check=True
    )

def run_opm_simulators_tests(python_bin: str) -> None:
    """Run all unittests in opm-simulators/python/test"""
    subprocess.run(
        [
            python_bin, "-m", "unittest", "discover",
            "-s", "/test/opm-simulators/python",  # path to the tests
            "-v"  # verbose output
        ],
        check=True
    )


def install_and_test(python_version, pyenv_full_version):
    """Install wheels for a single Python version and run a basic import test.
       Return True if everything succeeded, or False if it failed."""

    # Determine the correct wheel tags
    if tuple(map(int, python_version.split('.'))) < (3, 8):
        py_tag = f"cp{python_version.replace('.', '')}-cp{python_version.replace('.', '')}m"
    else:
        py_tag = f"cp{python_version.replace('.', '')}-cp{python_version.replace('.', '')}"

    # Find wheel files dynamically instead of hardcoding version
    wheel_path = f"/test/wheels"
    simulators_wheels = glob.glob(f"{wheel_path}/opm_simulators-*-{py_tag}-*manylinux*x86_64.whl")
    common_wheels = glob.glob(f"{wheel_path}/opm-*-{py_tag}-*manylinux*x86_64.whl")

    if not simulators_wheels:
        raise FileNotFoundError(f"No opm_simulators wheel found for Python {python_version} in {wheel_path}")
    if not common_wheels:
        raise FileNotFoundError(f"No opm-common wheel found for Python {python_version} in {wheel_path}")

    wheel_file_simulators = simulators_wheels[0]  # Use first match
    wheel_file_common = common_wheels[0]  # Use first match

    # Absolute path to the specific Python version
    python_bin = f"/root/.pyenv/versions/{pyenv_full_version}/bin/python"

    print(f"\n--- Installing wheels and running tests with Python {python_version} ---")
    print(f"Found opm-common wheel: {os.path.basename(wheel_file_common)}")
    print(f"Found opm_simulators wheel: {os.path.basename(wheel_file_simulators)}")

    try:
        # 1) Upgrade pip
        subprocess.run(
            [python_bin, "-m", "pip", "install", "--upgrade", "pip"],
            check=True
        )

        # 2) Install the wheels
        subprocess.run(
            [python_bin, "-m", "pip", "install",
             wheel_file_common,
             wheel_file_simulators],
            check=True
        )

        # 3) Run all unittests
        run_opm_simulators_tests(python_bin)
        run_opm_common_tests(python_bin)

        return True  # All good!

    except subprocess.CalledProcessError:
        print(f"Installation or test failed for Python {python_version}", file=sys.stderr)
        return False

def main():
    # Try to load version map, with fallback handling
    version_map_file = "/test/common/python_version_map.json"

    try:
        with open(version_map_file) as f:
            version_map = json.load(f)
        print(f"Loaded version map from {version_map_file}")
    except FileNotFoundError:
        print(f"Error: Version map file not found at {version_map_file}")
        print("This file should be created during Docker image build by install_python_versions.py")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in {version_map_file}: {e}")
        sys.exit(1)

    default_versions = ",".join(version_map.keys())
    python_versions_env = os.environ.get("PYTHON_VERSIONS", default_versions)

    print(f"Running tests.. Requested Python versions: {python_versions_env}")

    python_versions = [v.strip() for v in python_versions_env.split(",") if v.strip()]

    # Validate requested versions against available versions
    invalid_versions = []
    for ver in python_versions:
        if ver not in version_map:
            invalid_versions.append(ver)

    if invalid_versions:
        print(f"Error: Invalid Python version(s) requested: {', '.join(invalid_versions)}")
        print(f"   Available versions: {', '.join(sorted(version_map.keys()))}")
        sys.exit(1)

    # Track how many pass/fail
    results = []
    for ver in python_versions:
        full_version = version_map[ver]
        passed = install_and_test(ver, full_version)
        results.append((ver, passed))

    # Summarize
    total_tests = len(results)
    passed_tests = sum(1 for _, passed in results if passed)
    failed_tests = total_tests - passed_tests

    print("\n--- Summary of tests ---")
    if failed_tests > 0:
        print(f"{failed_tests}/{total_tests} tests failed.")
        print("result: FAILED")
        sys.exit(1)
    else:
        print(f"{passed_tests}/{total_tests} tests passed.")
        print("result: OK")
        sys.exit(0)

if __name__ == "__main__":
    main()
