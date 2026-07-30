#!/usr/bin/env bash
set -euo pipefail

mode="${1:-}"
if [[ "${mode}" != "probe" && "${mode}" != "full" ]]; then
  printf 'Usage: %s {probe|full}\n' "$0" >&2
  exit 2
fi

campaign_config="conf/phase2_multi_run/flagella_count_duration_3s_r1.yaml"
study_config="conf/phase4/duration_seed_study_v1_r1.yaml"
run_dir="outputs/phase2_multi_run/flagella_count_duration_3s_r1"
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
if [[ -e "${run_dir}" || -e "${dataset_dir}" ]]; then
  if [[ "${overwrite_requested}" != "true" ]]; then
    printf 'Existing output found. Re-run with OVERWRITE=true after inspection.\n' >&2
    exit 2
  fi
fi

if [[ "${overwrite_requested}" == "true" ]]; then
  uv run python scripts/01_simulate_swimming/run_multi_run.py \
    config="${campaign_config}" \
    overwrite=true
  uv run python scripts/02_phase2_analysis/build_dataset.py \
    config="${campaign_config}" \
    run_dir="${run_dir}" \
    overwrite=true
else
  uv run python scripts/01_simulate_swimming/run_multi_run.py \
    config="${campaign_config}"
  uv run python scripts/02_phase2_analysis/build_dataset.py \
    config="${campaign_config}" \
    run_dir="${run_dir}"
fi

uv run python scripts/04_phase4/evaluate_duration_seed_study.py \
  config="${study_config}"

printf 'Study output: %s\n' "${study_dir}"
