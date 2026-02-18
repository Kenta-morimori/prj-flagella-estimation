#!/usr/bin/env sh
# Configure repository-local git hooks.

set -eu

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

git config core.hooksPath .githooks
echo "Configured: core.hooksPath=$(git config --get core.hooksPath)"
echo "To disable: git config --unset core.hooksPath"
