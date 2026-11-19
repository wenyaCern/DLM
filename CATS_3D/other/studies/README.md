# Numerical studies and diagnostics

This folder holds standalone comparisons and scans of q/r grid choices,
source range, fixed-radius slices, and resonance checks. They are not required
to build the adaptive kernel generator, direct CF entry, or current data fitter.
The `other/build_config/CMakeLists.txt` still exposes their individual targets. Recheck the
input/output paths before running older studies; some write into `other/OutputFiles/`.

ROOT drawing macros and plotting scripts live in `other/analysis/`.
