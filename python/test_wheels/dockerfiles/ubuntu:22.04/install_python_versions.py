import os
import json
import subprocess

version_map = json.loads(os.environ["PYTHON_VERSION_MAP"])

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

print("✅ Python versions installed and pyenv configured.")
