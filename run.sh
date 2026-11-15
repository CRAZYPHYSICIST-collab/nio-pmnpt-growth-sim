#!/usr/bin/env bash
# Usage:  ./run.sh [test|exact|virtual|nonlinear|all] [output_dir]
# Requires: Mono (mcs + mono) and libgdiplus (for the figures)
set -e
cd "$(dirname "$0")"
MODE="${1:-all}"; OUT="${2:-output}"
mcs -optimize+ -r:System.Drawing.dll -r:System.Core.dll -out:NiOSim.exe NiOSim.cs
mono NiOSim.exe "$MODE" "$OUT"
