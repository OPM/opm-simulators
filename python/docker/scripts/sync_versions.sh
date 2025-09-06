#!/bin/bash

# - Sync Python versions from python_versions.json to static locations
# - This script updates locations that cannot dynamically read from JSON
# - Should be run after editing python_versions.json
#
# NOTE: This script in itself introduces a maintainability overhead; the alternative is to
#       manually update all the locations that are updated by this script. However, it is easy to
#       forget some of the locations, and the script also has documentation purpose giving the
#       necessary manual steps to update the locations.
set -e

DOCKER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPTS_DIR="$DOCKER_DIR/scripts"
PYTHON_VERSIONS_FILE="$DOCKER_DIR/python_versions.json"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

# Check if python_versions.json exists
if [ ! -f "$PYTHON_VERSIONS_FILE" ]; then
    log_error "python_versions.json not found at $PYTHON_VERSIONS_FILE"
    exit 1
fi

log_info "🔄 Syncing Python versions from python_versions.json to static locations..."

# Read versions from JSON (default_versions is now derived from supported_versions)
if ! command -v python3 >/dev/null 2>&1; then
    log_error "python3 not found - required to read JSON config"
    exit 1
fi

DEFAULT_VERSIONS=$(python3 "$SCRIPTS_DIR/read_python_config.py" default_versions 2>/dev/null)
if [ $? -ne 0 ] || [ -z "$DEFAULT_VERSIONS" ]; then
    log_error "Could not read default versions from JSON config"
    exit 1
fi

SUPPORTED_VERSIONS=$(python3 "$SCRIPTS_DIR/read_python_config.py" supported_versions 2>/dev/null)
if [ $? -ne 0 ] || [ -z "$SUPPORTED_VERSIONS" ]; then
    log_error "Could not read supported versions from JSON config"
    exit 1
fi

MANYLINUX_PLATFORM=$(python3 "$SCRIPTS_DIR/read_python_config.py" manylinux_platform 2>/dev/null)
if [ $? -ne 0 ] || [ -z "$MANYLINUX_PLATFORM" ]; then
    log_error "Could not read manylinux platform from JSON config"
    exit 1
fi

# Extract manylinux tag for Docker image (e.g., manylinux_2_28 from manylinux_2_28_x86_64)
MANYLINUX_TAG=$(echo "$MANYLINUX_PLATFORM" | sed 's/_x86_64$//')

log_info "📋 Versions from JSON (derived from supported_versions): $DEFAULT_VERSIONS"
log_info "📋 Manylinux platform: $MANYLINUX_PLATFORM (Docker tag: $MANYLINUX_TAG)"

# Extract minimum Python version (first when sorted)
MIN_VERSION=$(echo "$DEFAULT_VERSIONS" | tr ',' '\n' | sort -V | head -n1)
log_info "📋 Minimum Python version: $MIN_VERSION"

# Update Dockerfile default ARG and FROM line
DOCKERFILE="$DOCKER_DIR/Dockerfile"
if [ -f "$DOCKERFILE" ]; then
    log_info "📝 Updating $DOCKERFILE..."
    # Update the python_versions ARG line
    sed -i "s/^ARG python_versions=.*$/ARG python_versions=\"$DEFAULT_VERSIONS\"/" "$DOCKERFILE"
    log_success "Updated Dockerfile ARG python_versions"

    # Update the FROM line with the correct manylinux tag
    sed -i "s|^FROM quay.io/pypa/manylinux_[0-9_]*x86_64 AS stage1$|FROM quay.io/pypa/${MANYLINUX_TAG}_x86_64 AS stage1|" "$DOCKERFILE"
    log_success "Updated Dockerfile FROM line to use ${MANYLINUX_TAG}_x86_64"
else
    log_warning "Dockerfile not found at $DOCKERFILE"
fi

# Update test_wheels README.md
README_FILE="$DOCKER_DIR/test_wheels/scripts/README.md"
if [ -f "$README_FILE" ]; then
    log_info "📝 Updating $README_FILE..."
    # Update the supported versions line
    sed -i "s/- Supported: \`[^*]*\`/- Supported: \`$DEFAULT_VERSIONS\`/" "$README_FILE"
    log_success "Updated test_wheels README.md supported versions"
else
    log_warning "README.md not found at $README_FILE"
fi

# Update hardcoded version mappings in generate-pypi-package.sh
GEN_SCRIPT="$SCRIPTS_DIR/generate-pypi-package.sh"
if [ -f "$GEN_SCRIPT" ]; then
    log_info "📝 Updating version mappings in $GEN_SCRIPT..."

    # This is more complex as it involves the associative array
    # For now, just add a comment that these should match the JSON
    if ! grep -q "# NOTE: These version mappings should match python_versions.json" "$GEN_SCRIPT"; then
        sed -i '/^declare -A all_python_versions/i # NOTE: These version mappings should match python_versions.json - run sync_versions.sh after JSON changes' "$GEN_SCRIPT"
        log_success "Added maintenance comment to generate-pypi-package.sh"
    fi
else
    log_warning "generate-pypi-package.sh not found at $GEN_SCRIPT"
fi

# Update setup.py.in python_requires
SETUP_PY="$DOCKER_DIR/../setup.py.in"
if [ -f "$SETUP_PY" ]; then
    log_info "📝 Updating $SETUP_PY..."
    # Update the python_requires line
    sed -i "s/python_requires='>=.*'/python_requires='>=$MIN_VERSION'/" "$SETUP_PY"
    log_success "Updated setup.py.in python_requires to >=$MIN_VERSION"
else
    log_warning "setup.py.in not found at $SETUP_PY"
fi

# Summary
echo
log_success "🎉 Version sync completed!"
log_info "📍 Updated locations:"
log_info "  ✅ docker/Dockerfile (ARG python_versions, FROM line)"
log_info "  ✅ docker/test_wheels/scripts/README.md (supported versions)"
log_info "  ✅ python/setup.py.in (python_requires minimum version)"
log_info "  ⚠️  docker/scripts/generate-pypi-package.sh (manual update needed for version mappings)"
log_info "  ⚠️  docker/Dockerfile line 16 (manual update needed for manylinux compatibility comment)"
log_info "  ⚠️  docker/README.md line 166 (manual update needed for platform compatibility description)"
echo
log_info "📝 To add/remove Python versions:"
log_info "  1. Edit docker/python_versions.json"
log_info "  2. Run this script: ./docker/scripts/sync_versions.sh"
log_info "  3. Manually update generate-pypi-package.sh version mappings (lines 33-38)"
