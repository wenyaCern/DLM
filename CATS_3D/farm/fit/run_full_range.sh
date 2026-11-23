#!/usr/bin/env bash
# One 3D CATS fit cell using the completed farm kernel.
set -euo pipefail

if [[ $# -ne 3 && ! ( $# -eq 4 && $4 == --check ) ]]; then
  echo "Usage: $0 MT_INDEX QN_INDEX PHI_INDEX [--check]" >&2
  exit 2
fi
mt=$1
qn=$2
phi=$3
for value in "$mt" "$qn" "$phi"; do
  [[ $value =~ ^[0-9]+$ ]] || { echo "Indices must be nonnegative integers" >&2; exit 2; }
done
(( mt <= 2 && qn <= 4 && phi <= 5 )) || {
  echo "Allowed indices: mT 0..2, qn 0..4, phi 0..5" >&2
  exit 2
}

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
: "${CATS3D_OVERWRITE:=1}"
export CATS3D_FOLD_THREADS
if (( mt == 0 )); then
  : "${CATS3D_WEAK_FRACTION:=0.155}"
else
  : "${CATS3D_WEAK_FRACTION:?Set the mT-dependent CATS3D_WEAK_FRACTION for mT index $mt}"
fi
rinit=(8.2 5.5 4.5)
: "${CATS3D_RINIT:=${rinit[$mt]}}"
export CATS3D_DERIVED_CF CATS3D_LAMBDA_FILE

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
echo "Data: $CATS3D_DERIVED_CF"
echo "Lambda: $CATS3D_LAMBDA_FILE :: $CATS3D_LAMBDA_SELECTION"
echo "Smearing: $CATS3D_SMEARING_FILE"
echo "Kernel: $CATS3D_KERNEL_FILE"
echo "Kernel fold threads: $CATS3D_FOLD_THREADS"
echo "Weak-decay fraction per pion: $CATS3D_WEAK_FRACTION"
echo "Fit box: |q_out|, |q_side|, |q_long| < $CATS3D_FIT_QMAX_MEV MeV/c"
echo "Outer-box baseline and output: up to 200 MeV/c"
echo "Output: $out_file"
[[ ${4:-} == --check ]] && exit 0

mkdir -p "$out_dir"
[[ ! -e $out_file || ${CATS3D_OVERWRITE:-0} == 1 ]] || {
  echo "Output exists; set CATS3D_OVERWRITE=1 to replace it: $out_file" >&2
  exit 1
}
cd "$CATS3D_ROOT/fit/gaussian"
"$exe" "$mt" "$qn" "$phi" \
  "rinit=$CATS3D_RINIT" "fitqmax=$CATS3D_FIT_QMAX_MEV" "weakfraction=$CATS3D_WEAK_FRACTION" diagonal-band band-fixed kernel-fitbox-only \
  "${lambda_args[@]}" \
  "adaptive-kernel=$CATS3D_KERNEL_FILE" \
  "smearfile=$CATS3D_SMEARING_FILE" \
  "output=$out_file"
