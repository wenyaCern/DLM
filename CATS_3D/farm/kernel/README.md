# Farm pion kernel

This uses the same Coulomb + quantum-statistics wavefunction generator as the
local workflow. The `--farm` preset changes only its nonuniform q and r grids.

| Positive q interval (MeV/c) | Width (MeV/c) |
| --- | ---: |
| 0–5 | 1 |
| 5–15 | 2.5 |
| 15–40 | 5 |
| 40–80 | 10 |
| 80–200 | 20 |

| Positive r interval (fm) | Width (fm) |
| --- | ---: |
| 0–10 | 0.5 |
| 10–20 | 5/3 (approximately 1.667) |
| 20–40 | 2 |

The negative r axis mirrors the positive axis. There are 24 q cells per
positive axis and 72 r cells per full axis: `24^3 * 72^3` float values,
20.639 GB before ROOT compression. The 15–40 MeV/c region has five 5 MeV/c
cells. The finer low-q cells address the peak; the 2 fm outer-r cells avoid
an abrupt change to much coarser spatial spacing. A small three-point check
of the folded wavefunction supports this choice, but it does not establish a
global grid-convergence bound.

From the DLM repository root on the farm:

```sh
make -C CATS_3D/farm/kernel estimate
make -C CATS_3D/farm/kernel smoke
df -h /scratch8/Wenya
sbatch CATS_3D/farm/kernel/run_kernel.sbatch
```

After the job completes, run `make -C CATS_3D/farm/kernel check`. This scans
all q-cell indices and samples payloads without loading the full kernel into
RAM.

The default output is `/scratch8/Wenya/pipi_coulomb_q24_r72_nonuniform.root`.
The SLURM script refuses to overwrite an existing file. Its log is written
to `/scratch8/Wenya/cats3d_kernel_<jobid>.log`. Generation uses little RAM,
but loading the complete kernel for a fit requires at least 20.6 GB plus ROOT
and fit overhead; request a node with at least 48 GB for subsequent fits.
The output is not committed to Git. After generation, link it from
`CATS_3D/inputs/OutputFiles/` for the fit workflow.
