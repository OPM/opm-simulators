#!/usr/bin/env bash
# Invoke from build/opm-simulators/bin, inside compile_and_test.sh's container.
set -euo pipefail
if [[ ! -f /.dockerenv ]]; then
    echo 'Run this through bash dockerrunscript.sh, inside the container workflow.' >&2
    exit 2
fi
mode=${1:?Expected disabled, shadow, enabled, lifecycle, rollback, timing, timing-profile, or acceptance}
prefix=${2:-gpu_newton_20260922}
args=(
    --threads-per-process=16
    --newton-min-iterations=1
    --matrix-add-well-contributions=true
    --linear-solver-accelerator=gpu
    --linear-solver=dilu
    --enable-storage-cache=false
    --experimental-compute-properties-on-gpu=true
)
deck=/home/tobias/opm/RSC_cases/SPE11C_cases/RSC_16/deck/RSC_16_nodisp_nodiff
run_case()
{
    local name=$1
    shift
    if [[ -e ${prefix}_${name} ]]; then
        echo "Refusing to overwrite existing output: ${prefix}_${name}" >&2
        exit 2
    fi
    ./flow_gpu "${args[@]}" --output-dir="${prefix}_${name}" "$@" "$deck"
}
case "$mode" in
    disabled)
        run_case disabled --experimental-gpu-newton-update=false
        ;;
    shadow)
        run_case shadow --experimental-gpu-newton-validation=true --enable-opm-rst-file=true
        ;;
    enabled)
        run_case enabled
        ;;
    lifecycle)
        run_case lifecycle_disabled --experimental-gpu-newton-update=false \
            --enable-opm-rst-file=true
        run_case lifecycle_enabled --enable-opm-rst-file=true
        ;;
    rollback)
        run_case rollback_disabled --experimental-gpu-newton-update=false \
            --experimental-gpu-newton-reject-once=true --enable-opm-rst-file=true
        run_case rollback_enabled --experimental-gpu-newton-reject-once=true \
            --experimental-gpu-newton-validation=true --enable-opm-rst-file=true
        ;;
    timing)
        for repetition in 1 2 3; do
            run_case "timing_disabled_${repetition}" --experimental-gpu-newton-update=false
            run_case "timing_enabled_${repetition}"
        done
        ;;
    timing-profile)
        bash "${BASH_SOURCE[0]}" timing "$prefix"
        script_directory="$(dirname "${BASH_SOURCE[0]}")"
        for repetition in 1 2 3; do
            python3 "$script_directory/compare_spe11c.py" \
                "${prefix}_timing_disabled_${repetition}" "${prefix}_timing_enabled_${repetition}" \
                > "${prefix}_timing_comparison_${repetition}.json"
        done
        python3 "$script_directory/summarize_spe11c_timings.py" . "$prefix" \
            > "${prefix}_timings.json"
        bash "$script_directory/profile_spe11c_bridge.sh" "${prefix}_profile"
        python3 "$script_directory/compare_spe11c.py" \
            "${prefix}_profile_disabled" "${prefix}_profile_enabled" \
            > "${prefix}_profile_comparison.json"
        ;;
    acceptance)
        profiler_path=$(command -v rocprofv3 || true)
        if [[ -z $profiler_path && -x /opt/rocm/bin/rocprofv3 ]]; then
            profiler_path=/opt/rocm/bin/rocprofv3
        fi
        if [[ -n $profiler_path ]]; then
            "$profiler_path" --help > "${prefix}_rocprofv3_help.txt" 2>&1
        fi
        for phase in shadow lifecycle; do
            bash "${BASH_SOURCE[0]}" "$phase" "$prefix"
        done
        comparison_script="$(dirname "${BASH_SOURCE[0]}")/compare_spe11c.py"
        python3 "$comparison_script" "${prefix}_lifecycle_disabled" "${prefix}_shadow" \
            > "${prefix}_shadow_comparison.json"
        python3 "$comparison_script" "${prefix}_lifecycle_disabled" "${prefix}_lifecycle_enabled" \
            > "${prefix}_lifecycle_comparison.json"
        bash "${BASH_SOURCE[0]}" rollback "$prefix"
        python3 "$comparison_script" "${prefix}_rollback_disabled" "${prefix}_rollback_enabled" \
            > "${prefix}_rollback_comparison.json"
        bash "${BASH_SOURCE[0]}" timing "$prefix"
        for repetition in 1 2 3; do
            python3 "$comparison_script" "${prefix}_timing_disabled_${repetition}" \
                "${prefix}_timing_enabled_${repetition}" \
                > "${prefix}_timing_comparison_${repetition}.json"
        done
        for profiler in rocprofv3 rocprofv2 rocprof; do
            if command -v "$profiler"; then
                "$profiler" --help > "${prefix}_${profiler}_help.txt" 2>&1 || true
            fi
        done
        ;;
    *)
        echo "Unknown mode: $mode" >&2
        exit 2
        ;;
esac
