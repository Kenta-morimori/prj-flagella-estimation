#!/usr/bin/env bash
# Create the isolated cs10 runtime environment. Run from the repository root.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

python_path="${CS10_PYTHON:-$HOME/.local/bin/python3.11}"
venv_path="${CS10_VENV:-$repo_root/.venv-cs10}"
ffmpeg_version="7.0.2"
ffmpeg_asset="ffmpeg-${ffmpeg_version}-amd64-static.tar.xz"
ffmpeg_url="${CS10_FFMPEG_URL:-https://johnvansickle.com/ffmpeg/releases/${ffmpeg_asset}}"
ffmpeg_sha256="${CS10_FFMPEG_SHA256:-abda8d77ce8309141f83ab8edf0596834087c52467f6badf376a6a2a4c87cf67}"
ffmpeg_root="${CS10_FFMPEG_ROOT:-$HOME/.local/opt/ffmpeg-${ffmpeg_version}-amd64-static}"

install_ffmpeg() {
  if [[ -x "$ffmpeg_root/ffmpeg" && -x "$ffmpeg_root/ffprobe" ]]; then
    return
  fi
  local temporary_dir archive archive_sha
  temporary_dir="$(mktemp -d)"
  trap 'rm -rf "$temporary_dir"' RETURN
  archive="$temporary_dir/$ffmpeg_asset"
  curl --fail --location --retry 3 --output "$archive" "$ffmpeg_url"
  archive_sha="$(sha256sum "$archive" | awk '{print $1}')"
  if [[ "$archive_sha" != "$ffmpeg_sha256" ]]; then
    echo "FFmpeg archive SHA-256 mismatch: $archive_sha" >&2
    exit 2
  fi
  tar -xJf "$archive" -C "$temporary_dir"
  local extracted_dir
  extracted_dir="$(find "$temporary_dir" -mindepth 1 -maxdepth 1 -type d -name 'ffmpeg-*-amd64-static' -print -quit)"
  if [[ -z "$extracted_dir" || ! -x "$extracted_dir/ffmpeg" || ! -x "$extracted_dir/ffprobe" ]]; then
    echo "FFmpeg archive did not contain static ffmpeg and ffprobe binaries" >&2
    exit 2
  fi
  mkdir -p "$(dirname "$ffmpeg_root")" "$HOME/.local/bin"
  rm -rf "$ffmpeg_root"
  mv "$extracted_dir" "$ffmpeg_root"
  ln -sfn "$ffmpeg_root/ffmpeg" "$HOME/.local/bin/ffmpeg"
  ln -sfn "$ffmpeg_root/ffprobe" "$HOME/.local/bin/ffprobe"
}

if [[ ! -x "$python_path" ]]; then
  echo "cs10 CPython 3.11 not found: $python_path" >&2
  echo "Set CS10_PYTHON to the uv-managed CPython 3.11 executable." >&2
  exit 2
fi

glibc_line="$(ldd --version | head -n 1 || true)"
if [[ "$glibc_line" != *"2.17"* ]]; then
  echo "Expected glibc 2.17-compatible cs10 host; found: $glibc_line" >&2
  exit 2
fi

install_ffmpeg
"$ffmpeg_root/ffmpeg" -hide_banner -encoders | grep -q 'libx264'
"$ffmpeg_root/ffmpeg" -hide_banner -version | head -n 1

uv venv --clear --python "$python_path" "$venv_path"
uv pip install --python "$venv_path/bin/python" --only-binary :all: -r requirements/cs10.txt
uv pip check --python "$venv_path/bin/python"
probe_dir="outputs/$(TZ=Asia/Tokyo date +%F)/$(TZ=Asia/Tokyo date +%H%M%S)/cs10_setup_probe"
"$venv_path/bin/python" scripts/cs10/probe_runtime.py --output-dir "$probe_dir"

echo "cs10 environment ready: $venv_path"
