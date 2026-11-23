#!/usr/bin/env bash
# Optional overrides for run_full_range.sh. Copy to paths.sh on the farm.
# CATS3D_ROOT is set by the runner before this file is sourced.
# Upload the three ROOT inputs yourself; they are ignored by Git.
export CATS3D_DERIVED_CF="$CATS3D_ROOT/inputs/data/CF3D_derivedFull_Bminus_cent0to10_anyCPR.root"
export CATS3D_LAMBDA_FILE="$CATS3D_ROOT/inputs/data/PbPb_cent0to10_lambda_gen_vs_kstar_mt.root"
export CATS3D_SMEARING_FILE="$CATS3D_ROOT/inputs/data/AnalysisResults_lambdaPara_smearing_mTdiff.root"
export CATS3D_KERNEL_FILE="/scratch8/Wenya/pipi_coulomb_q24_r72_nonuniform.root"
export CATS3D_FIT_OUTPUT_DIR="$CATS3D_ROOT/farm/fit/farm_output/cent0to10_peak100_afterDCA_weak155"
export CATS3D_LAMBDA_SELECTION="aftercut"  # or fitwindow
# export CATS3D_FIT_QMAX_MEV="120"
# mT index 0 uses 0.155 by default. Set an mT-specific value for indices 1, 2.
# export CATS3D_WEAK_FRACTION="0.155"
