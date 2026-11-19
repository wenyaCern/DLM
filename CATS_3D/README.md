# CATS3D analysis workspace

The top level has six directories. The CATS3D physics class itself lives in
`../CATS/CATS3D.{h,cpp}` and is not modified by this workspace layout.

| Directory | Contents |
| --- | --- |
| `kernel_generation/` | Wavefunction kernel generation and checks |
| `direct_cf/` | Direct CATS3D correlation-function calculation |
| `fit/` | Gaussian and resonance Data fits, shared code, O-O/Ne-Ne fits |
| `inputs/` | Local ROOT inputs and their path manifest |
| `farm/` | Farm scripts, snapshot, and local farm output |
| `other/` | CMake configuration, analysis plots, studies, archives, and generated builds |

The large ROOT datasets, generated kernels, fit results, figures, archives,
and build products are kept locally and excluded from Git. Obtain the required
ROOT files separately, place them in `inputs/EPOS_data/`, `inputs/File/`, and
`inputs/OutputFiles/`, then run `python3 inputs/setup_local_paths.py`. The
current derived-data input paths are recorded in `inputs/data/README.md`.

Build the main targets from this directory:

```sh
make -C kernel_generation all
make -C direct_cf run_CATS3D
make -C fit gaussian
make -C fit resonance
```

The CMake graph is in `other/build_config/CMakeLists.txt`; generated binaries
are in `other/build/`. Legacy fit configurations use relative paths such as
`../File/`, so run their binaries with `other/build/` as the working directory.
ROOT plotting macros that use `File/`, `farm_output/`, or `figures/` should run
with `other/` as the working directory. New Python plotting scripts resolve
the workspace root from their own file location.

The old ROOT runner macros are retained in `other/archive/legacy_macros/`.
The source for bootstrap analysis is in `other/analysis/bootstrap/`.
