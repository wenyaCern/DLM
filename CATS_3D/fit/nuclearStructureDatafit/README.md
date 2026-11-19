# Nuclear structure data fits

O-O analysis lives here: `fit_CATS3D_oo_main.cpp`, plotting scripts,
`prepare_oo_pbpb_lambda.C`, `inputs/`, and `results/`.

The OO fit target is built from the workspace root with
`cmake --build other/build --target fit_CATS3D_oo`. Run the executable
from `other/build/` so its relative paths resolve.
The plotting scripts may be run from any directory;
they resolve the parent workspace from their file path.

Ne-Ne data and fit code can be added here alongside the O-O analysis. No
Ne-Ne-specific files were present during this reorganization.
