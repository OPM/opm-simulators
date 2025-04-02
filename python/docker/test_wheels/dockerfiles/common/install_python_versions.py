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
    version_map = config["supported_versions"]
    print(f"📋 Using Python versions from python_versions.json: {list(version_map.keys())}")
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
