#!/bin/bash

# Docker Build Helper Script
# Simplifies building PyPI wheels with standardized commands and logging
#
# Usage: ./run-docker-build.sh <tag> [options]
# Examples:
#   ./run-docker-build.sh static --output wheelhouse
#   ./run-docker-build.sh shared --libtype shared --output wheelhouse-shared
#   ./run-docker-build.sh py312 --python-versions "3.12" --output test-wheels

set -e

# Script directory and paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_DIR="$(dirname "$SCRIPT_DIR")"
DOCKER_DIR="$SCRIPT_DIR"
CONFIG_FILE="$DOCKER_DIR/python_versions.json"

# Default values
DEFAULT_LIBTYPE="static"
DEFAULT_JOBS="4"
DEFAULT_PYTHON_VERSIONS=""  # Will be read from config

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
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

# Show usage information
show_help() {
    cat << EOF
OPM Docker Build Helper Script

USAGE:
    $(basename "$0") <tag> [options]

ARGUMENTS:
    <tag>                   Docker tag name (e.g., static, shared, py312)

OPTIONS:
    --output <dir>          Extract wheels to directory (creates if needed)
    --libtype <type>        Library type: static or shared (default: static)
    --python-versions <versions>  Comma-separated Python versions (default: from config)
    --jobs <n>              Number of build jobs (default: 4)
    --version-common <ref>  Git ref for opm-common (branch/tag/PR) (default: master)
    --version-grid <ref>    Git ref for opm-grid (branch/tag/PR) (default: master)
    --version-simulators <ref>  Git ref for opm-simulators (branch/tag/PR) (default: master)
    --target-common <targets>   CMake targets for opm-common (default: opmcommon_python)
    --target-simulators <targets>   CMake targets for opm-simulators (default: simulators)
    --version-tag <tag>     Version tag for wheels (optional)
    --no-buildx             Use DOCKER_BUILDKIT=0 for detailed output
    --help                  Show this help message

EXAMPLES:
    $(basename "$0") static --output wheelhouse
        # Creates docker/wheelhouse/ with wheels
    $(basename "$0") shared --libtype shared --output wheelhouse-shared
        # Creates docker/wheelhouse-shared/ with wheels
    $(basename "$0") py312 --python-versions "3.12" --output test-wheels
        # Creates docker/test-wheels/ with wheels
    $(basename "$0") pr-test --version-simulators "pull/1234/head" --output pr-wheels
        # Test PR #1234 from opm-simulators
    $(basename "$0") pr6075 --version-simulators "pull/6075/head" --target-simulators "GasWater BlackOil"
        # Test PR #6075 with custom targets
    $(basename "$0") multi-pr --version-common "pull/567/head" --version-simulators "pull/1234/head"
        # Test multiple PRs together
    $(basename "$0") debug --libtype shared --no-buildx --jobs 8 --output /tmp/wheels
        # Creates /tmp/wheels/ with wheels (absolute path)
    $(basename "$0") release --version-tag "rc1" --output release-wheels
        # Creates docker/release-wheels/ with wheels tagged as version rc1

NOTES:
    - Script automatically runs from python/ directory
    - Relative output paths are created inside docker/ directory
    - Absolute output paths (starting with /) are used as-is
    - Logs are saved to docker/build-<tag>-<timestamp>.log
    - Default Python versions are read from docker/python_versions.json
    - Requires Docker to be installed and running

EOF
}

# Validate Docker is available
check_docker() {
    if ! command -v docker >/dev/null 2>&1; then
        log_error "Docker is not installed or not in PATH"
        exit 1
    fi

    if ! docker info >/dev/null 2>&1; then
        log_error "Docker daemon is not running or not accessible"
        exit 1
    fi
}

# Read default Python versions from config (derived from supported_versions)
get_default_python_versions() {
    if [ -f "$CONFIG_FILE" ] && command -v python3 >/dev/null 2>&1; then
        python3 -c "
import json, sys
try:
    with open('$CONFIG_FILE') as f:
        config = json.load(f)
    print(','.join(config['supported_versions'].keys()))
except:
    print('0.0')  # This should not happen
" 2>/dev/null || echo "0.0"  # This should never happen
    else
        echo "0.0"   # This should never happen
    fi
}

# Parse command line arguments
parse_arguments() {
    if [ $# -eq 0 ] || [ "$1" = "--help" ] || [ "$1" = "-h" ]; then
        show_help
        exit 0
    fi

    TAG="$1"
    shift

    # Initialize variables
    OUTPUT_DIR=""
    LIBTYPE="$DEFAULT_LIBTYPE"
    PYTHON_VERSIONS=""
    BUILD_JOBS="$DEFAULT_JOBS"
    VERSION_COMMON="master"
    VERSION_GRID="master"
    VERSION_SIMULATORS="master"
    TARGET_COMMON="opmcommon_python"
    TARGET_SIMULATORS="simulators"
    VERSION_TAG=""
    NO_BUILDX=false

    # Parse options
    while [[ $# -gt 0 ]]; do
        case $1 in
            --output)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            --libtype)
                LIBTYPE="$2"
                if [[ "$LIBTYPE" != "static" && "$LIBTYPE" != "shared" ]]; then
                    log_error "Invalid libtype: $LIBTYPE. Must be 'static' or 'shared'"
                    exit 1
                fi
                shift 2
                ;;
            --python-versions)
                PYTHON_VERSIONS="$2"
                shift 2
                ;;
            --jobs)
                BUILD_JOBS="$2"
                if ! [[ "$BUILD_JOBS" =~ ^[0-9]+$ ]]; then
                    log_error "Invalid jobs value: $BUILD_JOBS. Must be a positive integer"
                    exit 1
                fi
                shift 2
                ;;
            --version-common)
                VERSION_COMMON="$2"
                shift 2
                ;;
            --version-grid)
                VERSION_GRID="$2"
                shift 2
                ;;
            --version-simulators)
                VERSION_SIMULATORS="$2"
                shift 2
                ;;
            --target-common)
                TARGET_COMMON="$2"
                shift 2
                ;;
            --target-simulators)
                TARGET_SIMULATORS="$2"
                shift 2
                ;;
            --version-tag)
                VERSION_TAG="$2"
                shift 2
                ;;
            --no-buildx)
                NO_BUILDX=true
                shift
                ;;
            --help|-h)
                show_help
                exit 0
                ;;
            *)
                log_error "Unknown option: $1"
                echo
                show_help
                exit 1
                ;;
        esac
    done

    # Set default Python versions if not specified
    if [ -z "$PYTHON_VERSIONS" ]; then
        PYTHON_VERSIONS=$(get_default_python_versions)
        if [ "$PYTHON_VERSIONS" = "0.0" ]; then
            log_error "No default Python versions found in $CONFIG_FILE"
            exit 1
        fi
        log_info "Using default Python versions: $PYTHON_VERSIONS"
    fi

    # Validate tag name
    if [ -z "$TAG" ]; then
        log_error "Tag name is required"
        exit 1
    fi
}

# Main execution starts here
main() {
    # Parse arguments
    parse_arguments "$@"

    # Validate environment
    check_docker

    # Create timestamp for logging
    TIMESTAMP=$(date +%Y%m%d-%H%M%S)
    LOG_FILE="$DOCKER_DIR/build-${TAG}-${TIMESTAMP}.log"

    # Change to python directory
    cd "$PYTHON_DIR"

    log_info "Starting Docker build for tag: $TAG"
    log_info "Working directory: $(pwd)"
    log_info "Build configuration:"
    echo "  - Tag: $TAG"
    echo "  - Library type: $LIBTYPE"
    echo "  - Python versions: $PYTHON_VERSIONS"
    echo "  - Build jobs: $BUILD_JOBS"
    echo "  - OPM versions:"
    echo "    - opm-common: $VERSION_COMMON"
    echo "    - opm-grid: $VERSION_GRID"
    echo "    - opm-simulators: $VERSION_SIMULATORS"
    echo "  - Build targets:"
    echo "    - opm-common: $TARGET_COMMON"
    echo "    - opm-simulators: $TARGET_SIMULATORS"
    if [ -n "$OUTPUT_DIR" ]; then
        if [[ "$OUTPUT_DIR" = /* ]]; then
            echo "  - Output directory: $OUTPUT_DIR (absolute)"
        else
            echo "  - Output directory: docker/$OUTPUT_DIR (relative to docker/)"
        fi
    else
        echo "  - Output directory: none"
    fi
    echo "  - Log file: $LOG_FILE"
    echo

    # Construct Docker build command
    DOCKER_CMD=(docker build -t "manylinux2014_opm:$TAG" . -f docker/Dockerfile)

    # Add build arguments
    DOCKER_CMD+=(--build-arg "libtype=$LIBTYPE")
    DOCKER_CMD+=(--build-arg "python_versions=$PYTHON_VERSIONS")
    DOCKER_CMD+=(--build-arg "build_jobs=$BUILD_JOBS")
    DOCKER_CMD+=(--build-arg "version_common=$VERSION_COMMON")
    DOCKER_CMD+=(--build-arg "version_grid=$VERSION_GRID")
    DOCKER_CMD+=(--build-arg "version_simulators=$VERSION_SIMULATORS")
    DOCKER_CMD+=(--build-arg "target_common=$TARGET_COMMON")
    DOCKER_CMD+=(--build-arg "target_simulators=$TARGET_SIMULATORS")

    # Add version_tag if specified
    if [ -n "$VERSION_TAG" ]; then
        DOCKER_CMD+=(--build-arg "version_tag=$VERSION_TAG")
    fi

    # Add output option if specified
    if [ -n "$OUTPUT_DIR" ]; then
        # Resolve output directory relative to docker/ directory for relative paths
        if [[ "$OUTPUT_DIR" = /* ]]; then
            # Absolute path - use as-is
            RESOLVED_OUTPUT_DIR="$OUTPUT_DIR"
        else
            # Relative path - make relative to docker/ directory
            RESOLVED_OUTPUT_DIR="docker/$OUTPUT_DIR"
        fi

        DOCKER_CMD+=(--output "$RESOLVED_OUTPUT_DIR")
        log_info "Wheels will be extracted to: $RESOLVED_OUTPUT_DIR"
    fi

    # Set build environment
    if [ "$NO_BUILDX" = true ]; then
        export DOCKER_BUILDKIT=0
        log_info "Using legacy Docker build (detailed output)"
    fi

    # Show command being executed
    log_info "Executing: ${DOCKER_CMD[*]}"
    echo

    # Execute build with logging
    "${DOCKER_CMD[@]}" 2>&1 | tee "$LOG_FILE"

    # Check the exit status of the Docker command (not tee)
    DOCKER_EXIT_CODE=${PIPESTATUS[0]}

    if [ $DOCKER_EXIT_CODE -eq 0 ]; then
        echo
        log_success "Docker build completed successfully!"
        log_success "Tag: manylinux2014_opm:$TAG"

        if [ -n "$RESOLVED_OUTPUT_DIR" ] && [ -d "$RESOLVED_OUTPUT_DIR" ]; then
            WHEEL_COUNT=$(find "$RESOLVED_OUTPUT_DIR" -name "*.whl" | wc -l)
            log_success "Generated $WHEEL_COUNT wheel files in $RESOLVED_OUTPUT_DIR"
        fi

        log_info "Build log saved to: $LOG_FILE"
    else
        echo
        log_error "Docker build failed! (exit code: $DOCKER_EXIT_CODE)"
        log_error "Check the log file for details: $LOG_FILE"
        exit $DOCKER_EXIT_CODE
    fi
}

# Execute main function
main "$@"
