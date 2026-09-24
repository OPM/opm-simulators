#!/bin/bash

# This runs before opm-common's shared Jenkins helpers are available.
if [[ ${ghprbCommentBody:-} =~ (^|[[:space:]])opm-common=([^[:space:]]*) ]]
then
    requested_revision=${BASH_REMATCH[2]}
    remainder=${ghprbCommentBody#*"${BASH_REMATCH[0]}"}
    if grep -qiF 'opm-common=' <<< "$remainder"
    then
        echo "Multiple opm-common revisions specified" >&2
        exit 1
    fi
    if [[ -n ${absolute_revisions:-} ]]
    then
        if [[ -z $requested_revision ]]
        then
            echo "Invalid opm-common revision: a ref is required after opm-common=" >&2
            exit 1
        fi
        if [[ $requested_revision == -* ]] || ! git check-ref-format --branch "$requested_revision" >/dev/null 2>&1
        then
            echo "Invalid opm-common revision '$requested_revision': expected a branch, tag, or commit" >&2
            exit 1
        fi
        upstreamRev[opm-common]=$requested_revision
    else
        if [[ ! $requested_revision =~ ^[1-9][0-9]*$ ]]
        then
            echo "Invalid opm-common PR number '$requested_revision': expected a positive integer" >&2
            exit 1
        fi
        upstreamRev[opm-common]=pull/$requested_revision/merge
    fi
elif grep -qiF 'opm-common=' <<< "${ghprbCommentBody:-}"
then
    echo "Invalid opm-common trigger: use opm-common= as a standalone lowercase token" >&2
    exit 1
fi

repo_dir="$WORKSPACE/deps/opm-common"
repo_root=${OPM_REPO_ROOT:-git@github.com:OPM}
mkdir -p "$repo_dir" || exit 1
pushd "$repo_dir" || exit 1
if ! test -e .git
then
    git init . || exit 1
fi
if git remote get-url origin >/dev/null 2>&1
then
    git remote set-url origin "$repo_root/opm-common" || exit 1
else
    git remote add origin "$repo_root/opm-common" || exit 1
fi
if ! git fetch --depth 1 -- origin "${upstreamRev[opm-common]}"
then
    echo "Failed to fetch opm-common revision '${upstreamRev[opm-common]}'; check that the PR or ref exists" >&2
    exit 1
fi
if ! git checkout -B branch_to_build FETCH_HEAD
then
    echo "Failed to check out opm-common revision '${upstreamRev[opm-common]}'" >&2
    exit 1
fi
popd
