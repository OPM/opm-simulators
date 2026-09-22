#!/usr/bin/env bash
# Invoke only from build/opm-simulators/bin inside compile_and_test.sh.
set -euo pipefail
if [[ ! -f /.dockerenv ]]; then
    echo 'Run this through bash dockerrunscript.sh, inside the container workflow.' >&2
    exit 2
fi
prefix=${1:-gpu_newton_20260922_profile}
profiler=$(command -v rocprofv3 || true)
if [[ -z $profiler && -x /opt/rocm/bin/rocprofv3 ]]; then
    profiler=/opt/rocm/bin/rocprofv3
fi
if [[ -z $profiler ]]; then
    echo 'rocprofv3 is unavailable in the active Docker image' >&2
    exit 2
fi
for mode in disabled enabled; do
    output="${prefix}_${mode}"
    trace="${output}_trace"
    if [[ -e $output || -e $trace ]]; then
        echo "Refusing to overwrite existing output: $output or $trace" >&2
        exit 2
    fi
    enabled=true
    if [[ $mode == disabled ]]; then enabled=false; fi
    "$profiler" --hip-trace --memory-copy-trace --kernel-trace \
        --collection-period 5:2:1 --output-format csv --output-directory "$trace" \
        -- ./flow_gpu \
        --threads-per-process=16 \
        --newton-min-iterations=1 \
        --matrix-add-well-contributions=true \
        --linear-solver-accelerator=gpu \
        --linear-solver=dilu \
        --enable-storage-cache=false \
        --experimental-compute-properties-on-gpu=true \
        --experimental-gpu-newton-update="$enabled" \
        --output-dir="$output" \
        /home/tobias/opm/RSC_cases/SPE11C_cases/RSC_16/deck/RSC_16_nodisp_nodiff
done
