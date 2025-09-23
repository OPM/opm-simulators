#! /usr/bin/env python3

# This script is called as ENTRYPOINT from the Dockerfiles, e.g. docker/test_wheels/dockerfiles/.../Dockerfile
# For each python version in PYTHON_VERSIONS (with the values in docker/python_versions.json as default), it:
#  - installs the wheels for the python version
#  - runs the opm-common and opm-simulators python unit tests
# Returns 0 if all tests passed, 1 otherwise

import glob
import json
import os
import shutil
import subprocess
import sys

def should_stop_on_error() -> bool:
    """Check if we should stop on test errors based on environment variable."""
    return os.environ.get("STOP_ON_ERROR") == "1"

def copy_essential_directories(source_dir: str, dest_dir: str, essential_items: list) -> None:
    """Copy essential directories from source to destination, avoiding problematic items.

    Ensures clean slate by removing any existing destination directory first.
    """
    # Clean slate approach - delete and recreate for predictable, reproducible behavior
    if os.path.exists(dest_dir):
        shutil.rmtree(dest_dir)
    os.makedirs(dest_dir)

    for item in essential_items:
        item_path = os.path.join(source_dir, item)
        dest_path = os.path.join(dest_dir, item)

        if os.path.exists(item_path):
            if os.path.isdir(item_path):
                shutil.copytree(item_path, dest_path)  # No dirs_exist_ok needed - directory is guaranteed empty
            else:
                shutil.copy2(item_path, dest_path)

def run_opm_common_tests(python_bin: str) -> bool:
    """Run all unittests in opm-common/python/tests"""
    # Check for specific test selection
    test_cases = os.environ.get("TEST_CASES_COMMON")
    stop_on_error = should_stop_on_error()

    if os.environ.get("USE_HOST_TESTS"):
        # Copy strategy for host-mounted directories:
        # When using --host-tests-dir, the host directories contain 'opm/' source code
        # that conflicts with installed wheel packages. The local 'opm/' directories
        # shadow the installed packages, causing import failures because the local
        # source code doesn't have compiled extensions that the wheels provide.
        # Solution: Copy only test files (not conflicting opm/ directories) to /tmp
        # and run tests from there, allowing imports to find the installed wheel packages.
        test_dir = "/tmp/opm-common-tests"

        # Copy essential directories to clean location, avoiding conflicting opm/ directory
        copy_essential_directories("/test/opm-common/python", test_dir, ["tests", "test", "examples"])

        # Run tests from clean location (preserves package structure for relative imports)
        if test_cases:
            # Run specific test files
            all_passed = True
            for test_case in test_cases.split(","):
                test_case = test_case.strip()
                print(f"Running opm-common test: {test_case}")
                result = subprocess.run([
                    python_bin, "-m", "unittest", f"tests.test_{test_case}",
                    "-v"
                ], cwd=test_dir, check=False)
                if result.returncode != 0:
                    print(f"Test {test_case} failed with return code {result.returncode}", file=sys.stderr)
                    all_passed = False
                    if stop_on_error:
                        raise subprocess.CalledProcessError(result.returncode, result.args)
            if not all_passed:
                print("Some opm-common specific tests failed", file=sys.stderr)
            return all_passed
        else:
            # Run all tests
            result = subprocess.run([
                python_bin, "-m", "unittest", "discover",
                "-s", test_dir,
                "-p", "test_*.py",
                "-v"
            ], cwd=test_dir, check=False)
            if result.returncode != 0:
                print(f"opm-common tests failed with return code {result.returncode}", file=sys.stderr)
                if stop_on_error:
                    raise subprocess.CalledProcessError(result.returncode, result.args)
                return False
            return True
    else:
        if test_cases:
            # Run specific test files in cloned repo environment
            all_passed = True
            for test_case in test_cases.split(","):
                test_case = test_case.strip()
                print(f"Running opm-common test: {test_case}")
                result = subprocess.run([
                    python_bin, "-m", "unittest", f"test_{test_case}",
                    "-v"
                ], cwd="/test/opm-common/python", check=False)
                if result.returncode != 0:
                    print(f"Test {test_case} failed with return code {result.returncode}", file=sys.stderr)
                    all_passed = False
                    if stop_on_error:
                        raise subprocess.CalledProcessError(result.returncode, result.args)
            if not all_passed:
                print("Some opm-common specific tests failed", file=sys.stderr)
            return all_passed
        else:
            # Run from cloned repos (opm directory was moved to opm_original during build)
            result = subprocess.run([
                python_bin, "-m", "unittest", "discover",
                "-s", "/test/opm-common/python",  # path to the tests
                "-v"  # verbose output
            ], check=False)
            if result.returncode != 0:
                print(f"opm-common tests failed with return code {result.returncode}", file=sys.stderr)
                if stop_on_error:
                    raise subprocess.CalledProcessError(result.returncode, result.args)
                return False
            return True

def run_opm_simulators_tests(python_bin: str) -> bool:
    """Run all unittests in opm-simulators/python/test"""
    # Check for specific test selection
    test_cases = os.environ.get("TEST_CASES_SIMULATORS")
    stop_on_error = should_stop_on_error()

    if os.environ.get("USE_HOST_TESTS"):
        # Copy strategy for host-mounted directories:
        # Same issue as opm-common - host 'opm/' directories conflict with installed wheels.
        # Copy entire python directory structure except conflicting opm/ directory
        test_dir = "/tmp/opm-simulators-tests"

        # Copy essential directories to clean location, avoiding conflicting opm/ directory
        copy_essential_directories("/test/opm-simulators/python", test_dir, ["test", "test_data"])

        # Run tests from clean location (preserves package structure for relative imports)
        # Use parent directory as start and specify test pattern to make relative imports work
        if test_cases:
            # Run specific test files
            all_passed = True
            for test_case in test_cases.split(","):
                test_case = test_case.strip()
                print(f"Running opm-simulators test: {test_case}")
                result = subprocess.run([
                    python_bin, "-m", "unittest", f"test.test_{test_case}",
                    "-v"
                ], cwd=test_dir, check=False)
                if result.returncode != 0:
                    print(f"Test {test_case} failed with return code {result.returncode}", file=sys.stderr)
                    all_passed = False
                    if stop_on_error:
                        raise subprocess.CalledProcessError(result.returncode, result.args)
            if not all_passed:
                print("Some opm-simulators specific tests failed", file=sys.stderr)
            return all_passed
        else:
            # Run all tests
            result = subprocess.run([
                python_bin, "-m", "unittest", "discover",
                "-s", test_dir,
                "-p", "test_*.py",
                "-v"
            ], cwd=test_dir, check=False)
            if result.returncode != 0:
                print(f"opm-simulators tests failed with return code {result.returncode}", file=sys.stderr)
                if stop_on_error:
                    raise subprocess.CalledProcessError(result.returncode, result.args)
                return False
            return True
    else:
        if test_cases:
            # Run specific test files in cloned repo environment
            all_passed = True
            for test_case in test_cases.split(","):
                test_case = test_case.strip()
                print(f"Running opm-simulators test: {test_case}")
                result = subprocess.run([
                    python_bin, "-m", "unittest", f"test.test_{test_case}",
                    "-v"
                ], cwd="/test/opm-simulators/python", check=False)
                if result.returncode != 0:
                    print(f"Test {test_case} failed with return code {result.returncode}", file=sys.stderr)
                    all_passed = False
                    if stop_on_error:
                        raise subprocess.CalledProcessError(result.returncode, result.args)
            if not all_passed:
                print("Some opm-simulators specific tests failed", file=sys.stderr)
            return all_passed
        else:
            # Run from cloned repos (opm directory was moved to opm_original during build)
            result = subprocess.run([
                python_bin, "-m", "unittest", "discover",
                "-s", "/test/opm-simulators/python",  # path to the tests
                "-v"  # verbose output
            ], check=False)
            if result.returncode != 0:
                print(f"opm-simulators tests failed with return code {result.returncode}", file=sys.stderr)
                if stop_on_error:
                    raise subprocess.CalledProcessError(result.returncode, result.args)
                return False
            return True


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

        # 3) Run unittests based on what's requested
        test_cases_common = os.environ.get("TEST_CASES_COMMON")
        test_cases_simulators = os.environ.get("TEST_CASES_SIMULATORS")

        # Track overall test success
        overall_success = True

        # If no specific test cases are specified, run all tests
        if not test_cases_common and not test_cases_simulators:
            simulators_success = run_opm_simulators_tests(python_bin)
            common_success = run_opm_common_tests(python_bin)
            overall_success = simulators_success and common_success
        else:
            # Run only the requested test suites
            if test_cases_simulators or not test_cases_common:
                # Run simulators tests if explicitly requested OR if common tests not specified
                simulators_success = run_opm_simulators_tests(python_bin)
                overall_success = overall_success and simulators_success
            if test_cases_common or not test_cases_simulators:
                # Run common tests if explicitly requested OR if simulators tests not specified
                common_success = run_opm_common_tests(python_bin)
                overall_success = overall_success and common_success

        # Report overall result
        if overall_success:
            return True
        elif should_stop_on_error():
            # When stop_on_error is True, failures should stop execution
            return False
        else:
            # When stop_on_error is False, continue even if tests failed
            print(f"Tests failed for Python {python_version}, but continuing due to continue-on-error mode", file=sys.stderr)
            return True

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

    # Auto-detect actually installed Python versions instead of using all from JSON
    # This allows runtime to use only the versions that were installed during build
    try:
        result = subprocess.run(["pyenv", "versions", "--bare"],
                              capture_output=True, text=True, check=True)
        installed_versions = []
        for line in result.stdout.strip().split('\n'):
            if line.strip():
                # Extract short version (e.g., "3.12.9" -> "3.12")
                full_version = line.strip()
                short_version = '.'.join(full_version.split('.')[:2])
                if short_version in version_map:
                    installed_versions.append(short_version)

        if installed_versions:
            default_versions = ",".join(sorted(set(installed_versions)))
            print(f"Auto-detected installed Python versions: {default_versions}")
        else:
            # Fallback to JSON versions if detection fails
            default_versions = ",".join(version_map.keys())
            print(f"Could not detect installed versions, using all from JSON: {default_versions}")
    except Exception as e:
        # Fallback to JSON versions if detection fails
        default_versions = ",".join(version_map.keys())
        print(f"Could not detect installed versions ({e}), using all from JSON: {default_versions}")

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
