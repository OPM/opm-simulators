import logging
import zipfile

from pathlib import Path

import click

from opm_wheels import click_helpers, helpers
from opm_wheels.click_helpers import ClickOptions
from opm_wheels.constants import Directories, PythonVersion
from opm_wheels.config_reader import get_default_versions

@click.group()
@click.option("--verbose", "-v", is_flag=True, help="Show verbose output")
@click.pass_context
def main(ctx: click.Context, verbose: bool) -> None:
    ctx.ensure_object(dict)
    ctx.obj["VERBOSE"] = verbose
    if verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
        # logging.basicConfig(level=logging.WARNING)


@main.command()
@click.argument("name", type=str)
@ClickOptions.wheel_dir
def show_wheel_files(name: str, wheel_dir: str) -> None:
    """Show the files in wheel file. NAME is the name of a wheel file, relative to the
    python/wheelhouse directory."""
    if not name.endswith(".whl"):
        raise click.BadParameter("NAME must end with .whl")
    wheel_file = helpers.get_wheel_file(wheel_dir, name)
    # Check that
    with zipfile.ZipFile(wheel_file, "r") as zf:
        for name in zf.namelist():
            print(name)

@main.command()
@click.option(
    "--docker-os", "-d",
    type=str,
    default="ubuntu-22.04",
    help=(
        "The OS to use in the Docker image. Supported (canonical): "
        "ubuntu-22.04, ubuntu-24.04, debian-11. Legacy colon form accepted."
    )
)
@click.option(
    "--opm-simulators-repo", "-osr",
    type=str,
    default="https://github.com/OPM/opm-simulators.git",
    help="The URL of the opm-simulators repository."
)
@click.option(
    "--opm-common-repo", "-ocr",
    type=str,
    default="https://github.com/OPM/opm-common.git",
    help="The URL of the opm-common repository."
)
@click.option(
    "--opm-simulators-branch", "-osb",
    type=str,
    default="master",
    help="The branch of the opm-simulators repository."
)
@click.option(
    "--opm-common-branch", "-ocb",
    type=str,
    default="master",
    help="The branch of the opm-common repository."
)
@click.option(
    "--docker-progress", "-dp",
    is_flag=True,
    default=False,
    help="Show Docker progress output."
)
@click.option(
    "--python-versions", "-p",
    type=str,
    callback=click_helpers.validate_python_versions,
    help="Python versions to install in the Docker image. Valid values are: "
            f"{', '.join(PythonVersion.valid_versions())}. If not specified, all supported versions are installed."
)
def build_docker_image(
    docker_os: str,
    opm_simulators_repo: str,
    opm_common_repo: str,
    opm_simulators_branch: str,
    opm_common_branch: str,
    docker_progress: bool,
    python_versions: list[PythonVersion] = None
) -> None:
    """Build the Docker image. The option "--docker-os specifies the name of the OS to use
    in the Docker image. Supported: ubuntu:22.04, ubuntu:24.04, debian:11. If option --docker-progress
    is given, the Docker build will use option "--progress=plain". The default value is False."""

    testing_root = helpers.get_testing_root_dir()
    # Use docker directory as build context so python_versions.json is accessible
    docker_root = testing_root.parent  # python/docker/test_wheels -> python/docker
    canonical_os = helpers.canonicalize_docker_os(docker_os)
    docker_tag = helpers.get_docker_tag(canonical_os)
    # Prefer canonical directory, fall back to legacy colon directory (pre-rename)
    dockerfile = testing_root / Directories.docker_files / canonical_os / "Dockerfile"
    if not dockerfile.exists():
        legacy_os = canonical_os.replace("-", ":", 1)
        candidate = testing_root / Directories.docker_files / legacy_os / "Dockerfile"
        if candidate.exists():
            dockerfile = candidate
        else:
            raise FileNotFoundError(f"Dockerfile {dockerfile} does not exist")
    # Build the Docker image
    helpers.run_docker_build(
        docker_tag,
        opm_simulators_repo,
        opm_common_repo,
        opm_simulators_branch,
        opm_common_branch,
        build_ctxt_dir=docker_root,
        docker_file=dockerfile,
        docker_progress=docker_progress,
        python_versions=python_versions
    )

@main.command()
@click.option(
    "--docker-os", "-d",
    type=str,
    default="ubuntu-22.04",
    help=(
        "The OS for the Docker image. Supported (canonical): "
        "ubuntu-22.04, ubuntu-24.04, debian-11. Legacy colon form accepted."
    )
)
@ClickOptions.wheel_dir
@click.option(
    "--python-versions", "-p",
    type=str,
    callback=click_helpers.validate_python_versions,
    help="The Python versions to use for testing. Valid values are: "
            f"{', '.join(PythonVersion.valid_versions())}"
)
@click.option(
    "--host-tests-dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
    help="Directory containing opm-common and opm-simulators repositories from host. "
         "If specified, tests will be run from host directories instead of cloned repos. "
         "Expected structure: host-tests-dir/opm-common/python and host-tests-dir/opm-simulators/python"
)
@click.option(
    "--test-cases-common",
    type=str,
    help="Comma-separated list of specific test case patterns to run for opm-common (e.g., 'eclfile,esmry'). "
         "Will run test_<pattern>.py files. If not specified, all tests are run."
)
@click.option(
    "--test-cases-simulators",
    type=str,
    help="Comma-separated list of specific test case patterns to run for opm-simulators (e.g., 'basic,mpi'). "
         "Will run test_<pattern>.py files. If not specified, all tests are run."
)
@click.option(
    "--stop-on-error",
    is_flag=True,
    default=False,
    help="Stop execution on first test failure. If not specified, all tests will run even if some fail."
)
def run_docker_image(
    docker_os: str,
    wheel_dir: str,
    python_versions: list[PythonVersion] | None,
    host_tests_dir: Path = None,
    test_cases_common: str = None,
    test_cases_simulators: str = None,
    stop_on_error: bool = False,
) -> None:
    f"""Run the Docker image. The option "--docker-tag" specifies the tag of the Docker
    image to run. The option "--python-versions" specifies the Python versions to
    use for testing. Valid values are: {', '.join(PythonVersion.valid_versions())}. The value should
    be a comma-separated list of versions without separating space. For example: 3.8,3.9,3.10.
    If not given, the default value is {','.join(get_default_versions())}."""
    wheel_dir = helpers.get_wheel_abs_dir(wheel_dir)
    canonical_os = helpers.canonicalize_docker_os(docker_os)
    docker_tag = helpers.get_docker_tag(canonical_os)
    helpers.run_docker_run(docker_tag, wheel_dir, python_versions, host_tests_dir, test_cases_common, test_cases_simulators, stop_on_error)
