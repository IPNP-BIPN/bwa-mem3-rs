#!/usr/bin/env bash
# Local CI gate: format, lint, test. Mirrors what CI runs.
set -euo pipefail
cd "$(dirname "$0")/.."

echo "cargo: $(cargo --version)"
echo "== rustfmt =="
cargo fmt --check
echo "== clippy (-D warnings) =="
cargo clippy --workspace --all-targets -- -D warnings
echo "== test =="
cargo test --workspace
echo "OK"
