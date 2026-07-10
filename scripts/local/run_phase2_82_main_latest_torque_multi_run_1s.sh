#!/usr/bin/env bash
set -euo pipefail

# Requires a clean git worktree because run_multi_run.py checks git status.

CONFIG="conf/phase2_multi_run/latest_model_torque_shape_stability.yaml"
TORQUE_VALUES="[1.5e-20,2.0e-20,2.5e-20]"

echo "[1/3] phase2_82 latest-model torque multi-run (1.0 s)"
summary_csv="$(
  uv run python scripts/01_simulate_swimming/run_multi_run.py \
    config="${CONFIG}" \
    time.duration_s=1.0 \
    sweep.axes.torque.values="${TORQUE_VALUES}"
)"
summary_csv="$(printf '%s\n' "${summary_csv}" | tail -n 1)"
run_root="$(dirname "${summary_csv}")"

echo "summary_csv=${summary_csv}"
echo "run_root=${run_root}"

echo "[2/3] phase2_82 latest-model torque summary plot"
uv run python scripts/01_simulate_swimming/plot_heatmap.py \
  config="${CONFIG}" \
  summary_csv="${summary_csv}" \
  output_dir="${run_root}/plots"

echo "[3/3] phase2_82 latest-model torque replay"
uv run python scripts/01_simulate_swimming/render_shape_stability_grid_replay.py \
  --input-dir "${run_root}" \
  --mode both \
  --fps-out-3d 10 \
  --output-dir "${run_root}/replay" \
  --overwrite

echo "Outputs:"
echo "  ${summary_csv}"
echo "  ${run_root}/run_manifest.json"
echo "  ${run_root}/plots/plot_data.csv"
echo "  ${run_root}/plots/first_fail_t_s_vs_torque.png"
echo "  ${run_root}/plots/hook_len_rel_err_max_vs_torque.png"
echo "  ${run_root}/plots/max_flag_bond_rel_err_vs_torque.png"
echo "  ${run_root}/plots/axis_center_to_body_roll_ratio_mean_vs_torque.png"
echo "  ${run_root}/replay/shape_stability_metrics.csv"
echo "  ${run_root}/replay/shape_stability_metrics.png"
echo "  ${run_root}/replay/grid_swim3d.mp4"
