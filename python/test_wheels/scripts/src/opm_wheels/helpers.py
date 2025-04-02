import logging
import subprocess
import sys

from pathlib import Path

from opm_wheels.constants import Directories, PythonVersion

def get_docker_tag(docker_os: str) -> str:
    """Get the Docker tag for the given OS. Currently only ubuntu:22.04 is supported."""
    return f"test-opm-wheels-{docker_os}"

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
    # directory called wheels in the root directory
    if not Path(root).joinpath("wheels").exists():
        raise RuntimeError("Current cwd is not inside repository.")
    return Path(root)

def get_testing_root_dir() -> Path:
    """Return the absolute path of the testing directory."""
    root = get_git_root()
    testing_root = root / "python" / "test_wheels"
    if not testing_root.exists():
        raise RuntimeError("Testing directory does not exist.")
    return testing_root

def get_wheel_dir() -> Path:
    """Return the absolute path of the wheel directory."""
    testing_root = get_testing_root_dir()
    wheel_dir = testing_root / Directories.wheelhouse
    if not wheel_dir.exists():
        raise RuntimeError("Wheel directory does not exist.")
    return wheel_dir

def get_wheel_file(name: str) -> Path:
    """Get the path to the wheel file. NAME is the name of a wheel file, relative to the
    python/wheelhouse directory of the Git repository."""
    wheel_dir = get_wheel_dir()
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
) -> None:
    if docker_progress:
        progress="plain"
    else:
        progress="auto"
    try:
        result = subprocess.run(
            [
                "docker", "build",
                "--progress", progress,
                "--build-arg", "opm_simulators_repo=" + opm_simulators_repo,
                "--build-arg", "opm_common_repo=" + opm_common_repo,
                "--build-arg", "opm_simulators_branch=" + opm_simulators_branch,
                "--build-arg", "opm_common_branch=" + opm_common_branch,
                "-t", docker_tag,
                "-f", str(docker_file),
                str(build_ctxt_dir)
            ],
            check=True,  # Raise an exception if the process returns a non-zero exit code
        )
        logging.info(f"Docker image with tag {docker_tag} built successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error building Docker image for tag {docker_tag}: {e}")
        sys.exit(1)
    return

def run_docker_run(
    docker_tag: str,
    wheel_dir: Path,
    python_versions: list[PythonVersion]
) -> None:
    try:
#-v $WHEELDIR_HOST:${WHEELDIR_CONTAINER}
        result = subprocess.run(
            [
                "docker", "run",
                "--rm",
                "-v", str(wheel_dir) + ":/test/wheels/",
                "-e", "PYTHON_VERSIONS=" + ",".join([str(v) for v in python_versions]),
                docker_tag
            ],
            check=True,  # Raise an exception if the process returns a non-zero exit code
        )
        logging.info(f"Docker image with tag {docker_tag} run successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running Docker image for tag {docker_tag}: {e}")
        sys.exit(1)
    return
