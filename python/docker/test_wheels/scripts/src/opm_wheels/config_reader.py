import json
import os
from pathlib import Path
from typing import Dict, List

from opm_wheels.constants import PythonVersion
from opm_wheels.helpers import get_git_root


def get_config_path() -> Path:
    """Get the path to the python_versions.json config file."""
    try:
        # Try to find git root (works when inside repo)
        repo_root = get_git_root()
        return repo_root / "python/docker/python_versions.json"
    except RuntimeError:
        # Fallback to environment variable for global installations
        env_root = os.environ.get('OPM_SIMULATORS_ROOT')
        if env_root:
            repo_root = Path(env_root)
            config_path = repo_root / "python/docker/python_versions.json"
            if config_path.exists():
                return config_path
        
        raise FileNotFoundError(
            "Could not locate python_versions.json. Either run from within "
            "opm-simulators repository or set OPM_SIMULATORS_ROOT environment variable."
        )


def load_python_config() -> Dict:
    """Load the python versions configuration from JSON."""
    config_path = get_config_path()
    if not config_path.exists():
        raise FileNotFoundError(f"Python config file not found at {config_path}")

    with open(config_path, 'r') as f:
        return json.load(f)


def get_supported_versions() -> Dict[str, str]:
    """Get the mapping of short versions to full versions."""
    config = load_python_config()
    return config["supported_versions"]


def get_default_versions() -> List[str]:
    """Get the list of default Python versions (derived from supported_versions)."""
    config = load_python_config()
    return list(config["supported_versions"].keys())


def get_python_version_objects() -> List[PythonVersion]:
    """Get PythonVersion objects for all default versions."""
    default_versions = get_default_versions()
    return [PythonVersion.from_str(version) for version in default_versions]