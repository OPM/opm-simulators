import logging
import zipfile

from pathlib import Path

import click

from opm_wheels import click_helpers, helpers
from opm_wheels.click_helpers import ClickOptions
from opm_wheels.constants import Directories, PythonVersion

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
    default="ubuntu:22.04",
    help="The OS to use in the Docker image. Currently only ubuntu:22.04 is supported."
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
def build_docker_image(
    docker_os: str,
    opm_simulators_repo: str,
    opm_common_repo: str,
    opm_simulators_branch: str,
    opm_common_branch: str,
    docker_progress: bool
) -> None:
    """Build the Docker image. The option "--docker-os specifies the name of the OS to use
    in the Docker image. Currently only ubuntu:22.04 is supported. If option --docker-progress
    is given, the Docker build will use option "--progress=plain". The default value is False."""
    testing_root = helpers.get_testing_root_dir()
    docker_tag = helpers.get_docker_tag(docker_os)
    dockerfile = testing_root / Directories.docker_files / docker_os / "Dockerfile"
    if not dockerfile.exists():
        raise FileNotFoundError(f"Dockerfile {dockerfile} does not exist")
    # Build the Docker image
    helpers.run_docker_build(
        docker_tag,
        opm_simulators_repo,
        opm_common_repo,
        opm_simulators_branch,
        opm_common_branch,
        build_ctxt_dir=testing_root,
        docker_file=dockerfile,
        docker_progress=docker_progress
    )

@main.command()
@click.option(
    "--docker-os", "-d",
    type=str,
    default="ubuntu:22.04",
    help="The OS for the Docker image. Currently only ubuntu:22.04 is supported."
)
@ClickOptions.wheel_dir
@click.option(
    "--python-versions", "-p",
    type=str,
    callback=click_helpers.validate_python_versions,
    help="The Python versions to use for testing. Valid values are: "
            f"{', '.join(PythonVersion.valid_versions())}"
)
def run_docker_image(
    docker_os: str,
    wheel_dir: str,
    python_versions: list[PythonVersion],
) -> None:
    """Run the Docker image. The option "--docker-tag" specifies the tag of the Docker
    image to run. The option "--python-versions" specifies the Python versions to
    use for testing. Valid values are: 3.6, 3.7, 3.8, 3.9, 3.10, 3.11, 3.12. The value should
    be a comma-separated list of versions without separating space. For example: 3.6,3.7,3.8.
    If not given, the default value is 3.6,3.7,3.8,3.9,3.10,3.11,3.12."""
    wheel_dir = helpers.get_wheel_abs_dir(wheel_dir)
    docker_tag = helpers.get_docker_tag(docker_os)
    helpers.run_docker_run(docker_tag, wheel_dir, python_versions)
