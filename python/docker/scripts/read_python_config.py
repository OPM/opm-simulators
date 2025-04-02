#!/usr/bin/env python3
"""Helper script to read Python configuration from JSON config file.

This script is called from:
- docker/scripts/generate-pypi-package.sh
- docker/scripts/sync_versions.sh

The script is used to read the Python versions from the docker/python_versions.json config file.
"""

import json
import sys
from pathlib import Path

def get_config_path():
    """Get the path to the python_versions.json config file."""
    script_dir = Path(__file__).parent
    # Config is now at docker root
    config_path = script_dir.parent / "python_versions.json"
    return config_path

def main():
    if len(sys.argv) != 2:
        print("Usage: python read_python_config.py <key>", file=sys.stderr)
        print("Available keys: default_versions (derived), supported_versions", file=sys.stderr)
        sys.exit(1)

    key = sys.argv[1]
    config_path = get_config_path()

    if not config_path.exists():
        print(f"Error: Config file not found at {config_path}", file=sys.stderr)
        sys.exit(1)

    try:
        with open(config_path, 'r') as f:
            config = json.load(f)

        # Handle default_versions by deriving from supported_versions
        if key == "default_versions":
            if "supported_versions" not in config:
                print("Error: supported_versions not found in config", file=sys.stderr)
                sys.exit(1)
            result = list(config["supported_versions"].keys())
            print(",".join(result))
        elif key in config:
            result = config[key]
            if key == "supported_versions":
                # Output as JSON for shell parsing
                print(json.dumps(result))
            else:
                print(result)
        else:
            print(f"Error: Key '{key}' not found in config", file=sys.stderr)
            sys.exit(1)

    except Exception as e:
        print(f"Error reading config: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
