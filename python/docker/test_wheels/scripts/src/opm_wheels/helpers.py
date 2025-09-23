import logging
import subprocess
import sys

from pathlib import Path

from opm_wheels.constants import Directories, PythonVersion


def canonicalize_docker_os(name: str) -> str:
    """Return a canonical filesystem/tag-friendly OS identifier.

    Accepts legacy values like "ubuntu:22.04" and returns "ubuntu-22.04".
    If already hyphenated (e.g., "ubuntu-22.04"), it is returned unchanged.
    """
    return name.replace(":", "-")

def get_docker_tag(docker_os: str) -> str:
    """Get the Docker tag for the given OS.

    The tag is always produced from the canonicalized OS name, e.g.:
    - input: "ubuntu:22.04" or "ubuntu-22.04" => tag: "test-opm-wheels-ubuntu-22.04"

    Supported values (canonical): ubuntu-22.04, ubuntu-24.04, debian-11.
    """
    canonical = canonicalize_docker_os(docker_os)
    return f"test-opm-wheels-{canonical}"

def get_git_root() -> Path:
    """Return the absolute path of the repository's root."""
    try:
        output = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            stderr=subprocess.STDOUT
        )
    except subprocess.CalledProcessError:
        # Handle the error if we're not inside a Git repo, etc.
        raise RuntimeError("Not a valid Git repository or other error occurred.")
    root = output.decode("utf-8").strip()
    # Extra check to ensure we're in the correct repository, check that there is a
    # directory called python/docker/test_wheels in the root directory

    if not Path(root).joinpath("python/docker/test_wheels").exists():
        raise RuntimeError("Current cwd is not inside repository.")
    return Path(root)

def get_testing_root_dir() -> Path:
    """Return the absolute path of the testing directory."""
    root = get_git_root()
    testing_root = root / "python" / "docker" / "test_wheels"
    if not testing_root.exists():
        raise RuntimeError("Testing directory does not exist.")
    return testing_root

def get_wheel_default_dir() -> Path:
    """Return the default wheel directory as a relative path."""
    wheel_default_dir = Path(Directories.python) / Directories.wheelhouse
    return wheel_default_dir

def get_wheel_abs_dir(wheel_dir: Path|str) -> Path:
    """Return the absolute path of the wheel directory."""
    # Check if the wheel_dir is an absolute path
    wheel_dir = Path(wheel_dir)
    if wheel_dir.is_absolute():
        return wheel_dir
    # If it's not absolute, resolve it relative to the python/docker directory
    git_root = get_git_root()
    docker_root = git_root / "python" / "docker"
    wheel_dir = docker_root / wheel_dir
    if not wheel_dir.exists():
        raise RuntimeError(f"Wheel directory {wheel_dir} does not exist.")
    return wheel_dir

def get_wheel_file(wheel_dir: str, name: str) -> Path:
    """Get the path to the wheel file. NAME is the name of a wheel file, relative to the
    python/wheelhouse directory of the Git repository."""
    wheel_dir = get_wheel_abs_dir(wheel_dir)
    wheel_file = wheel_dir / name
    if not wheel_file.exists():
        raise FileNotFoundError(f"Wheel file {wheel_file} does not exist")
    return wheel_file

def run_docker_build(
    docker_tag: str,
    opm_simulators_repo: str,
    opm_common_repo: str,
    opm_simulators_branch: str,
    opm_common_branch: str,
    build_ctxt_dir: Path,
    docker_file: Path,
    docker_progress: bool = False,
    python_versions: list[PythonVersion] = None,
) -> None:
    if docker_progress:
        progress="plain"
    else:
        progress="auto"
    # Build Docker command with base arguments
    docker_cmd = [
        "docker", "build",
        "--progress", progress,
        "--build-arg", "opm_simulators_repo=" + opm_simulators_repo,
        "--build-arg", "opm_common_repo=" + opm_common_repo,
        "--build-arg", "opm_simulators_branch=" + opm_simulators_branch,
        "--build-arg", "opm_common_branch=" + opm_common_branch,
    ]

    # Add Python versions build argument if specified
    if python_versions:
        python_versions_str = ",".join([str(v) for v in python_versions])
        docker_cmd.extend(["--build-arg", f"BUILD_PYTHON_VERSIONS={python_versions_str}"])
        logging.info(f"Building Docker image with Python versions: {python_versions_str}")

    # Add tag and context
    docker_cmd.extend([
        "-t", docker_tag,
        "-f", str(docker_file),
        str(build_ctxt_dir)
    ])

    try:
        result = subprocess.run(docker_cmd, check=True)
        logging.info(f"Docker image with tag {docker_tag} built successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error building Docker image for tag {docker_tag}: {e}")
        sys.exit(1)
    return

def run_docker_run(
    docker_tag: str,
    wheel_dir: Path,
    python_versions: list[PythonVersion] | None,
    host_tests_dir: Path = None,
    test_cases_common: str = None,
    test_cases_simulators: str = None,
    stop_on_error: bool = False
) -> None:
    # Build Docker run command
    docker_cmd = [
        "docker", "run",
        "--rm",
        "-v", str(wheel_dir) + ":/test/wheels/",
    ]

    # Only set PYTHON_VERSIONS environment variable if explicitly specified
    # If not specified (None), let the container auto-detect installed versions
    if python_versions is not None:
        docker_cmd.extend(["-e", "PYTHON_VERSIONS=" + ",".join([str(v) for v in python_versions])])
        logging.info(f"Using specified Python versions: {','.join([str(v) for v in python_versions])}")
    else:
        logging.info("Using auto-detected Python versions from container")

    # Add host test directory mounts if specified
    if host_tests_dir:
        # Validate that the expected directories exist
        opm_common_python = host_tests_dir / "opm-common" / "python"
        opm_simulators_python = host_tests_dir / "opm-simulators" / "python"

        if not opm_common_python.exists():
            raise FileNotFoundError(f"Expected opm-common/python directory not found at: {opm_common_python}")
        if not opm_simulators_python.exists():
            raise FileNotFoundError(f"Expected opm-simulators/python directory not found at: {opm_simulators_python}")

        # Add volume mounts for test directories
        docker_cmd.extend([
            "-v", str(opm_common_python) + ":/test/opm-common/python",
            "-v", str(opm_simulators_python) + ":/test/opm-simulators/python"
        ])

        logging.info(f"Using host test directories from: {host_tests_dir}")

        # Signal to entrypoint that host directories are mounted
        docker_cmd.extend(["-e", "USE_HOST_TESTS=1"])

    # Add test selection environment variables if specified
    if test_cases_common:
        docker_cmd.extend(["-e", f"TEST_CASES_COMMON={test_cases_common}"])
        logging.info(f"Running only opm-common tests: {test_cases_common}")

    if test_cases_simulators:
        docker_cmd.extend(["-e", f"TEST_CASES_SIMULATORS={test_cases_simulators}"])
        logging.info(f"Running only opm-simulators tests: {test_cases_simulators}")

    # Add stop on error environment variable
    if stop_on_error:
        docker_cmd.extend(["-e", "STOP_ON_ERROR=1"])
        logging.info("Running with stop-on-error enabled")
    else:
        logging.info("Running with continue-on-error (default behavior)")

    # Add the Docker image tag
    docker_cmd.append(docker_tag)

    try:
        # Debug: log the full Docker command
        logging.info(f"DEBUG: Executing Docker command: {' '.join(docker_cmd)}")
        result = subprocess.run(
            docker_cmd,
            check=True,  # Raise an exception if the process returns a non-zero exit code
        )
        logging.info(f"Docker image with tag {docker_tag} run successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running Docker image for tag {docker_tag}: {e}")
        sys.exit(1)
    return
