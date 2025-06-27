#! /usr/bin/env python3

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

    wheel_name_simulators = f"opm_simulators-2025.4-{py_tag}-manylinux_2_28_x86_64.whl"
    wheel_name_common = f"opm-2025.4-{py_tag}-manylinux_2_28_x86_64.whl"
    wheel_path = f"/test/wheels"

    # Absolute path to the specific Python version
    python_bin = f"/root/.pyenv/versions/{pyenv_full_version}/bin/python"

    print(f"\n--- Installing wheels and running tests with Python {python_version} ---")

    try:
        # 1) Upgrade pip
        subprocess.run(
            [python_bin, "-m", "pip", "install", "--upgrade", "pip"],
            check=True
        )

        # 2) Install the wheels
        subprocess.run(
            [python_bin, "-m", "pip", "install",
             f"{wheel_path}/{wheel_name_common}",
             f"{wheel_path}/{wheel_name_simulators}"],
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
    with open("/test/common/python_version_map.json") as f:
        version_map = json.load(f)
    default_versions = ",".join(version_map.keys())
    python_versions_env = os.environ.get("PYTHON_VERSIONS", default_versions)

    print(f"Running tests.. Requested Python versions: {python_versions_env}")

    python_versions = [v.strip() for v in python_versions_env.split(",") if v.strip()]

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
