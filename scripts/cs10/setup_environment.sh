#!/usr/bin/env bash
# Create the isolated cs10 runtime environment. Run from the repository root.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

python_path="${CS10_PYTHON:-$HOME/.local/bin/python3.11}"
venv_path="${CS10_VENV:-$repo_root/.venv-cs10}"

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

uv venv --clear --python "$python_path" "$venv_path"
uv pip install --python "$venv_path/bin/python" --only-binary :all: -r requirements/cs10.txt
uv pip check --python "$venv_path/bin/python"
probe_dir="outputs/$(TZ=Asia/Tokyo date +%F)/$(TZ=Asia/Tokyo date +%H%M%S)/cs10_setup_probe"
"$venv_path/bin/python" scripts/cs10/probe_runtime.py --output-dir "$probe_dir"

echo "cs10 environment ready: $venv_path"
