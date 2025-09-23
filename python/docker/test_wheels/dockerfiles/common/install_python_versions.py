import os
import json
import subprocess

# This script is called from the Dockerfiles, e.g. docker/test_wheels/dockerfiles/.../Dockerfile
# It is used to install the Python versions specified in docker/python_versions.json
try:
    # Use the file copied by Docker build
    config_file = "/test/python_versions.json"
    with open(config_file) as f:
        config = json.load(f)
    full_version_map = config["supported_versions"]

    # Check if specific Python versions were requested via build argument
    build_versions = os.environ.get("BUILD_PYTHON_VERSIONS")
    if build_versions:
        requested_versions = [v.strip() for v in build_versions.split(",")]
        version_map = {v: full_version_map[v] for v in requested_versions if v in full_version_map}
        print(f"📋 Using requested Python versions: {list(version_map.keys())}")

        # Validate that all requested versions are supported
        missing_versions = [v for v in requested_versions if v not in full_version_map]
        if missing_versions:
            print(f"Error: Unsupported Python version(s): {missing_versions}")
            print(f"Supported versions: {list(full_version_map.keys())}")
            exit(1)
    else:
        version_map = full_version_map
        print(f"📋 Using all Python versions from python_versions.json: {list(version_map.keys())}")

except FileNotFoundError:
    print("Error: python_versions.json not found!")
    exit(1)
except KeyError as e:
    print(f"Error: Missing key {e} in python_versions.json")
    exit(1)

# Initialize pyenv
subprocess.run('eval "$(pyenv init -)"', shell=True, executable="/bin/bash")

# Install all full versions
for short, full in version_map.items():
    print(f"Installing Python {full} (for {short})...")
    subprocess.check_call(["pyenv", "install", full])

# Set global version to latest
latest_full_version = version_map[max(version_map)]
subprocess.check_call(["pyenv", "global", latest_full_version])
subprocess.run("pyenv rehash", shell=True)

# Create the version map file that entrypoint.py expects
version_map_output = "/test/common/python_version_map.json"
try:
    # Ensure the directory exists
    os.makedirs(os.path.dirname(version_map_output), exist_ok=True)

    # Write the version map to the expected location
    with open(version_map_output, "w") as f:
        json.dump(version_map, f, indent=2)
    print(f"Created version map file at {version_map_output}")
except Exception as e:
    print(f"Warning: Could not create version map file: {e}")

print("Python versions installed and pyenv configured.")
