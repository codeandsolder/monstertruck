#!/usr/bin/env bash
set -euo pipefail

# Runtime limits apply only to the built validator, never its compiler.
# Explicit executable path supports managed ephemeral Cargo artifact export.
if (( $# < 4 )); then
    echo "Usage: $0 VALIDATOR TIMEOUT_SECONDS VM_KIB COMMAND [ARGS...]" >&2
    exit 2
fi
validator=$(realpath -- "$1")
seconds=$2
vm_kib=$3
shift 3
[[ -x "$validator" && "$seconds" =~ ^[1-9][0-9]*$ && "$vm_kib" =~ ^[1-9][0-9]*$ ]]
ulimit -v "$vm_kib"
exec timeout --signal=KILL "$seconds" "$validator" "$@"
