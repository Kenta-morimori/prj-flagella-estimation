#!/usr/bin/env bash
set -euo pipefail

mode="${1:-}"
if [[ "${mode}" != "probe" && "${mode}" != "full" ]]; then
  printf 'Usage: %s {probe|full}\n' "$0" >&2
  exit 2
fi

campaign_config="conf/phase2_multi_run/flagella_count_duration_3s_r1.yaml"
study_config="conf/phase4/duration_seed_study_v1_r1.yaml"
timestamp="$(TZ=Asia/Tokyo date +'%Y-%m-%d/%H%M%S')"
run_dir="${RUN_DIR:-outputs/${timestamp}/phase4_duration_seed_study}"
dataset_dir="${run_dir}/dataset/v1_r1_duration_3s"
study_dir="${run_dir}/analysis/phase4_duration_seed_study"

if [[ "${mode}" == "probe" ]]; then
  uv run python scripts/01_simulate_swimming/run_multi_run.py \
    config="${campaign_config}" \
    dry_run=true \
    sample_limit=3
  exit 0
fi

overwrite_requested="${OVERWRITE:-false}"
if [[ -e "${run_dir}" ]]; then
  if [[ "${overwrite_requested}" != "true" ]]; then
    printf 'Existing output found. Re-run with OVERWRITE=true after inspection.\n' >&2
    exit 2
  fi
fi

campaign_overrides=(
  "output.base_dir=${run_dir}"
  "output.timestamp_subdir=false"
  "dataset.output_dir=${dataset_dir}"
)

if [[ "${overwrite_requested}" == "true" ]]; then
  uv run python scripts/01_simulate_swimming/run_multi_run.py \
    config="${campaign_config}" \
    "${campaign_overrides[@]}" \
    overwrite=true
  uv run python scripts/02_phase2_analysis/build_dataset.py \
    config="${campaign_config}" \
    run_dir="${run_dir}" \
    "${campaign_overrides[@]}" \
    overwrite=true
else
  uv run python scripts/01_simulate_swimming/run_multi_run.py \
    config="${campaign_config}" \
    "${campaign_overrides[@]}"
  uv run python scripts/02_phase2_analysis/build_dataset.py \
    config="${campaign_config}" \
    run_dir="${run_dir}" \
    "${campaign_overrides[@]}"
fi

uv run python scripts/04_phase4/evaluate_duration_seed_study.py \
  config="${study_config}" \
  dataset_dir="${dataset_dir}" \
  output_dir="${study_dir}"

printf 'Campaign root: %s\n' "${run_dir}"
printf 'Study output: %s\n' "${study_dir}"
