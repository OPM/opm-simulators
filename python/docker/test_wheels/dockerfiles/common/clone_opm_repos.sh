#!/bin/bash

# This script is called from the Dockerfiles, e.g. docker/test_wheels/dockerfiles/.../Dockerfile

set -eux

# These should be passed as environment variables from the Dockerfile.
: "${OPM_COMMON_REPO:?}"
: "${OPM_SIMULATORS_REPO:?}"
: "${OPM_COMMON_BRANCH:?}"
: "${OPM_SIMULATORS_BRANCH:?}"

# Function to clone repository with PR support (similar to generate-pypi-package.sh)
clone_repo_for_testing() {
    local repo_url=$1
    local branch=$2
    local target_dir=$3

    if [[ "$branch" == pull/* ]]; then
        # Handle PR references
        echo "Cloning $(basename "$repo_url") and fetching PR: $branch"
        git clone --depth=1 "$repo_url" "$target_dir"
        cd "$target_dir"
        git fetch origin "$branch"
        git checkout FETCH_HEAD
        cd ..
    else
        # Handle regular branches/tags
        echo "Cloning $(basename "$repo_url") at branch: $branch"
        git clone --depth=1 "$repo_url" -b "$branch" "$target_dir"
    fi
}

# Clone each repository with PR support
clone_repo_for_testing "$OPM_COMMON_REPO" "$OPM_COMMON_BRANCH" "opm-common"
clone_repo_for_testing "$OPM_SIMULATORS_REPO" "$OPM_SIMULATORS_BRANCH" "opm-simulators"

mv opm-simulators/python/opm opm-simulators/python/opm_original
mv opm-common/python/opm opm-common/python/opm_original
