from opm_wheels.constants import LibType, PythonVersion
from typing import Union
import click

def validate_lib_type(
    ctx: click.Context, param: Union[click.Option, click.Parameter], value: str
) -> str:
    if value is None:
        value = LibType.SHARED.value
    if value not in LibType.values():
        raise click.BadParameter(
            f"Invalid library type '{value}'. Must be one of {', '.join(LibType.values())}."
        )
    return value

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
