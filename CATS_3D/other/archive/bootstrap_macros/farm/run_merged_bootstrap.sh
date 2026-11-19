#!/usr/bin/env bash
set -euo pipefail
# SUBDIR must match the requested plane's actual subsamples; a label cannot rotate a plane.
plane=${1:?FT0A or FT0C}
mt=${2:?merged mT 0 or 1}
qn=${3:?qn 0..2}
phi=${4:?phi 0..5}
subdir=${5:?absolute merged subsample directory}
jobs=${6:-2}
[[ "$subdir" = /* && -d "$subdir" ]] || { echo 'SUBDIR must be an existing absolute directory'; exit 2; }
[[ "$jobs" =~ ^[1-9][0-9]*$ ]] || exit 2
cd "$(dirname "$0")"
mkdir -p logs
export OMP_NUM_THREADS=1
export plane mt qn phi subdir
seq 0 99 | xargs -P "$jobs" -I{} bash -c '
  cd build
  exec ./fit_CATS3D_farm merged_bootstrap "$plane" "$mt" "$qn" "$phi" "$1" "$subdir" > "../logs/merged_boot_${plane}_mt${mt}_qn${qn}_phi${phi}_${1}.log" 2>&1
' _ {}
