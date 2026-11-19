# Data-fit inputs

All active fit inputs are read through `CATS_3D/inputs/`. On this laptop,
`inputs/data/*.root` are local symlinks to the large ROOT files listed below;
on the farm, place or link the corresponding files under `inputs/data/`.
The `*.root` files and links are excluded from Git.

- SAME/MIXED CF: `/Users/wuwenya/Documents/code/femto_flow/CATS/3D_pion_CATS/output/derivedFull_Bminus/cent0to10_anyCPR/CF3D_derivedFull_Bminus_cent0to10_anyCPR.root`. Histograms are `mt{m}_qn{q}_phi{p}/same3d` and `mixed3d`.
- Lambda: `/Users/wuwenya/Documents/code/femto_flow/CATS/3D_pion_CATS/ROOT/lambdaPara/results_lambda_gen_projection/TPCOnlyMC_fraction_mergedFirst3/derived_test/PbPb_cent0to10_lambda_gen_derived_fullPt.root`, histogram `lambdaFitWindow`.
- Current 0–10% q response: `/Users/wuwenya/Documents/code/femto_flow/CATS/3D_pion_CATS/ROOT/lambdaPara/AnalysisResults_lambdaPara_smearing3DLCMS_test_cent0to10_cent10to30_mT02to04.root`, folder `femto-pair-track-track/TrackTrack/SE/MC/`.
- Earlier inclusive-centrality PID q response: `/Users/wuwenya/Documents/code/femto_flow/CATS/3D_pion_CATS/ROOT/lambdaPara/AnalysisResults_lambdaPara_smearing3DLCMS_test_PID.root`.
- Older input ROOT files and lambda variants are read from `inputs/File/`.

The current q-response histograms in the first new smearing file are divided
by centrality but not actually gated by the displayed `hMt` axis range.
