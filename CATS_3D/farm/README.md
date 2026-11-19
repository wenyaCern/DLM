# Farm production workspace

The farm has the same three calculation paths as local work:

| Folder | Role |
| --- | --- |
| `kernel/` | Build the shared kernel generator; farm uses a finer, larger grid |
| `fit/` | Build the shared fit code and run production nominal, bootstrap, and systematic jobs |
| `direct_cf/` | Build the shared direct-CF code |

Local runs are for workflow tests. Farm production must use the same source and
fit model, with only the kernel grid made finer. The exact farm grid, bootstrap
inputs, systematic variations, and submission scripts are being defined one at
a time. The previous independent farm fit snapshot is preserved locally in
`other/archive/farm_fit_snapshot_2026-11-19/`; it is not the active source.

## Installation on the farm

`quick_install_farm.py` is the original farm installer, kept unchanged. It
uses the directory containing the script as the DLM repository root. To run
it, copy it to the repository root first (without replacing `quick_install.py`):

```sh
cd /path/to/DLM
cp CATS_3D/farm/quick_install_farm.py ./quick_install_farm.py
python3 ./quick_install_farm.py basic
```

After CATS is installed under `install/`, build the three shared modules with
`make -C CATS_3D/farm/kernel`, `make -C CATS_3D/farm/fit`, and
`make -C CATS_3D/farm/direct_cf`. The input ROOT files and kernels are not in
Git; put them under `CATS_3D/inputs/` or link them there from scratch storage.
