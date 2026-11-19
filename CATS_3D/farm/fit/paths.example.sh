#!/usr/bin/env bash
# Copy to paths.sh on the farm and replace filenames with the uploaded inputs.
# These paths are resolved while running from CATS_3D/other/build/.
# Inputs are centralized under CATS_3D/inputs/. A large farm kernel may live
# on scratch storage, with a symlink in inputs/OutputFiles/.

export CATS3D_DERIVED_CF="../../inputs/data/CF3D_derivedFull_Bminus_cent0to10_anyCPR.root"
export CATS3D_LAMBDA_FILE="../../inputs/data/PbPb_cent0to10_lambda_gen_derived_fullPt.root"
export CATS3D_SMEARING_FILE="../../inputs/data/AnalysisResults_lambdaPara_smearing3DLCMS_cent0to10.root"
export CATS3D_KERNEL_FILE="../../inputs/OutputFiles/dataKernel_qData_r51_rmax20_coulomb_byQ.root"
export CATS3D_FIT_OUTPUT_DIR="../farm_output/derivedFull_Bminus_cent0to10_anyCPR_outerbox_weak14_fit"
