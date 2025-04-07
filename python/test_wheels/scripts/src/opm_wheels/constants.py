import enum

class BaseEnum(enum.Enum):
    """Base class for enums with keys and values."""

    @classmethod
    def keys(cls) -> set[str]:
        return set(cls.__members__.keys())

    @classmethod
    def values(cls) -> set[str]:
        return {member.value for member in cls}


class Directories:
    docker_files = "dockerfiles"
    python = "python"
    wheelhouse = "wheelhouse"


class PythonVersion(BaseEnum):
    """Python version enum."""
    v3_6 = "3.6"
    v3_7 = "3.7"
    v3_8 = "3.8"
    v3_9 = "3.9"
    v3_10 = "3.10"
    v3_11 = "3.11"
    v3_12 = "3.12"

    @classmethod
    def from_str(cls, version_str):
        try:
            return cls[f"v{version_str.replace('.', '_')}"]
        except KeyError:
            valid_versions = ', '.join([str(version) for version in cls])
            raise ValueError(f"Invalid python version: {version_str}, valid versions are {valid_versions}")

    def __str__(self):
        # Return the version in the desired format (e.g., "3.8" instead of "v3_8")
        return self.name[1:].replace('_', '.')

    @classmethod
    def valid_versions(cls) -> list[str]:
        """Return a list of valid python versions."""
        return cls.values()
