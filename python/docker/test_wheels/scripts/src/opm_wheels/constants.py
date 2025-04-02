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


# Load Python versions from JSON to create enum dynamically
def _get_python_versions_for_enum():
    """Get Python versions from JSON for enum creation."""
    try:
        from pathlib import Path
        import json

        # Get path to python_versions.json
        current_dir = Path(__file__).parent
        config_path = current_dir / ".." / ".." / ".." / ".." / "python_versions.json"
        config_path = config_path.resolve()

        with open(config_path) as f:
            config = json.load(f)

        return list(config["supported_versions"].keys())
    except Exception:
        # Fallback to hardcoded values
        return ["3.8", "3.9", "3.10", "3.11", "3.12", "3.13"]

# Dynamically create PythonVersion enum using functional API
_python_versions = _get_python_versions_for_enum()

# Create enum members dict
enum_members = {}
for version in _python_versions:
    enum_name = f'v{version.replace(".", "_")}'
    enum_members[enum_name] = version

# Create PythonVersion enum dynamically
PythonVersion = enum.Enum('PythonVersion', enum_members, type=BaseEnum)

# Add class methods
@classmethod
def _from_str(cls, version_str):
    try:
        return cls[f"v{version_str.replace('.', '_')}"]
    except KeyError:
        valid_versions = ', '.join([member.value for member in cls])
        raise ValueError(f"Invalid python version: {version_str}, valid versions are {valid_versions}")

def _str(self):
    # Return the version in the desired format (e.g., "3.8" instead of "v3_8")
    return self.value

@classmethod
def _valid_versions(cls) -> list[str]:
    """Return a list of valid python versions."""
    return [member.value for member in cls]

# Bind methods to class
PythonVersion.from_str = _from_str
PythonVersion.__str__ = _str
PythonVersion.valid_versions = _valid_versions
