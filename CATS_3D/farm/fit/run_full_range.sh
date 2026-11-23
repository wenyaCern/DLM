#!/usr/bin/env bash
# One 3D CATS fit cell, or a batch that reuses one loaded kernel.
set -euo pipefail

batch=0
if [[ ${1:-} == --batch ]]; then
  if [[ $# -ne 6 && ! ( $# -eq 7 && $7 == --check ) ]]; then
    echo "Usage: $0 --batch MT QN_FIRST QN_LAST PHI_FIRST PHI_LAST [--check]" >&2
    exit 2
  fi
  batch=1
  mt=$2; qn_first=$3; qn_last=$4; phi_first=$5; phi_last=$6
  qn=$qn_first; phi=$phi_first
  check_flag=${7:-}
  values=("$mt" "$qn_first" "$qn_last" "$phi_first" "$phi_last")
else
  if [[ $# -ne 3 && ! ( $# -eq 4 && $4 == --check ) ]]; then
    echo "Usage: $0 MT QN PHI [--check]" >&2
    exit 2
  fi
  mt=$1; qn=$2; phi=$3
  check_flag=${4:-}
  values=("$mt" "$qn" "$phi")
fi
for value in "${values[@]}"; do
  [[ $value =~ ^[0-9]+$ ]] || { echo "Indices must be nonnegative integers" >&2; exit 2; }
done
(( mt <= 2 && qn <= 4 && phi <= 5 )) || {
  echo "Allowed indices: mT 0..2, qn 0..4, phi 0..5" >&2
  exit 2
}
if (( batch && (qn_first > qn_last || qn_last > 4 || phi_first > phi_last || phi_last > 5) )); then
  echo "Invalid batch range" >&2
  exit 2
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export CATS3D_ROOT=$(cd "$script_dir/../.." && pwd)
if [[ -f $script_dir/paths.sh ]]; then
  # Optional farm-only path overrides; this file is ignored by Git.
  source "$script_dir/paths.sh"
fi
: "${CATS3D_DERIVED_CF:=$CATS3D_ROOT/inputs/data/CF3D_derivedFull_Bminus_cent0to10_anyCPR.root}"
: "${CATS3D_LAMBDA_FILE:=$CATS3D_ROOT/inputs/data/PbPb_cent0to10_lambda_gen_vs_kstar_mt.root}"
: "${CATS3D_SMEARING_FILE:=$CATS3D_ROOT/inputs/data/AnalysisResults_lambdaPara_smearing_mTdiff.root}"
: "${CATS3D_KERNEL_FILE:=/scratch8/Wenya/pipi_coulomb_q24_r72_nonuniform.root}"
: "${CATS3D_FIT_OUTPUT_DIR:=$script_dir/farm_output/cent0to10_peak100_afterDCA_weak155}"
: "${CATS3D_LAMBDA_SELECTION:=aftercut}"
: "${CATS3D_FIT_QMAX_MEV:=120}"
: "${CATS3D_FOLD_THREADS:=4}"
: "${CATS3D_USE_SIMPLEX:=0}"
: "${CATS3D_OVERWRITE:=1}"
export CATS3D_FOLD_THREADS
case $CATS3D_USE_SIMPLEX in
  0) simplex_args=(no-simplex) ;;
  1) simplex_args=() ;;
  *) echo "CATS3D_USE_SIMPLEX must be 0 or 1" >&2; exit 2 ;;
esac
if (( mt == 0 )); then
  : "${CATS3D_WEAK_FRACTION:=0.155}"
else
  : "${CATS3D_WEAK_FRACTION:?Set the mT-dependent CATS3D_WEAK_FRACTION for mT index $mt}"
fi
rinit=(8.2 5.5 4.5)
: "${CATS3D_RINIT:=${rinit[$mt]}}"
export CATS3D_DERIVED_CF CATS3D_LAMBDA_FILE CATS3D_FIT_OUTPUT_DIR

case $CATS3D_LAMBDA_SELECTION in
  aftercut) lambda_args=(lambda-after-dca-cut) ;;
  fitwindow) lambda_args=() ;;
  *) echo "CATS3D_LAMBDA_SELECTION must be aftercut or fitwindow" >&2; exit 2 ;;
esac
for name in CATS3D_DERIVED_CF CATS3D_LAMBDA_FILE CATS3D_SMEARING_FILE CATS3D_KERNEL_FILE; do
  path=${!name}
  [[ -s $path ]] || { echo "Missing or empty $name: $path" >&2; exit 1; }
done
exe="$CATS3D_ROOT/other/build/fit_CATS3D_derived_full"
[[ -x $exe ]] || { echo "Build first: make -C CATS_3D/farm/fit gaussian" >&2; exit 1; }

out_dir="$CATS3D_FIT_OUTPUT_DIR/mt$mt/qn$qn"
out_file="$out_dir/fit_phi$phi.root"
if (( batch )); then
  cell_args=(--batch "$mt" "$qn_first" "$qn_last" "$phi_first" "$phi_last")
  output_args=()
else
  cell_args=("$mt" "$qn" "$phi")
  output_args=("output=$out_file")
fi
echo "Data: $CATS3D_DERIVED_CF"
echo "Lambda: $CATS3D_LAMBDA_FILE :: $CATS3D_LAMBDA_SELECTION"
echo "Smearing: $CATS3D_SMEARING_FILE"
echo "Kernel: $CATS3D_KERNEL_FILE"
echo "Kernel fold threads: $CATS3D_FOLD_THREADS"
echo "Minuit SIMPLEX pre-search: $CATS3D_USE_SIMPLEX"
echo "Weak-decay fraction per pion: $CATS3D_WEAK_FRACTION"
echo "Fit box: |q_out|, |q_side|, |q_long| < $CATS3D_FIT_QMAX_MEV MeV/c"
echo "Outer-box baseline and output: up to 200 MeV/c"
if (( batch )); then
  echo "Cells: mT=$mt, qn=$qn_first..$qn_last, phi=$phi_first..$phi_last"
  echo "Output directory: $CATS3D_FIT_OUTPUT_DIR"
else
  echo "Output: $out_file"
fi
[[ $check_flag == --check ]] && exit 0

if (( ! batch )); then
  mkdir -p "$out_dir"
elif [[ $CATS3D_OVERWRITE != 1 ]]; then
  for (( iqn=qn_first; iqn<=qn_last; ++iqn )); do
    for (( iphi=phi_first; iphi<=phi_last; ++iphi )); do
      existing="$CATS3D_FIT_OUTPUT_DIR/mt$mt/qn$iqn/fit_phi$iphi.root"
      [[ ! -e $existing ]] || { echo "Output exists: $existing" >&2; exit 1; }
    done
  done
fi
[[ $batch == 1 || ! -e $out_file || ${CATS3D_OVERWRITE:-0} == 1 ]] || {
  echo "Output exists; set CATS3D_OVERWRITE=1 to replace it: $out_file" >&2
  exit 1
}
cd "$CATS3D_ROOT/fit/gaussian"
"$exe" "${cell_args[@]}" \
  "rinit=$CATS3D_RINIT" "fitqmax=$CATS3D_FIT_QMAX_MEV" "weakfraction=$CATS3D_WEAK_FRACTION" diagonal-band band-fixed kernel-fitbox-only \
  "${lambda_args[@]}" \
  "adaptive-kernel=$CATS3D_KERNEL_FILE" \
  "smearfile=$CATS3D_SMEARING_FILE" \
  "${simplex_args[@]}" \
  "${output_args[@]}"
