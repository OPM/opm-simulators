from opm_wheels.constants import PythonVersion
from typing import Union
import click

from opm_wheels import helpers

def validate_python_versions(
    ctx: click.Context, param: Union[click.Option, click.Parameter], value: str
) -> list[PythonVersion]:
    if value is None:
        value = ",".join(PythonVersion.valid_versions())
    versions = value.split(",")
    versions = [PythonVersion.from_str(str(v)) for v in versions]
    valid_versions = [PythonVersion.from_str(v) for v in PythonVersion.values()]
    for version in versions:
        if version not in valid_versions:
            raise click.BadParameter(
                f"Invalid python version '{version}'. Must be one of {', '.join([str(v) for v in valid_versions])}."
            )
    return versions

class ClickOptions:
    wheel_dir = lambda func: click.option(
        "--wheel-dir", "-wd",
        type=str,
        default=str(helpers.get_wheel_default_dir()),
        help=f"The directory containing the wheels. Can be relative to the root of the opm-simulators repository, or an absolute path. The default value is {str(helpers.get_wheel_default_dir())} directory."
    )(func)
