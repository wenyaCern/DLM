# Farm pion kernel

This uses the same Coulomb + quantum-statistics wavefunction generator as the
local workflow. The `--farm` preset changes only its nonuniform q and r grids.

| Positive q interval (MeV/c) | Width (MeV/c) |
| --- | ---: |
| 0–30 | 2.5 |
| 30–80 | 5 |
| 80–120 | 10 |
| 120–200 | 20 |

| Positive r interval (fm) | Width (fm) |
| --- | ---: |
| 0–12 | 0.75 |
| 12–24 | 2 |
| 24–40 | 3.2 |

The negative r axis mirrors the positive axis. There are 30 q cells per
positive axis and 54 r cells per full axis: `30^3 * 54^3` float values,
17.006 GB before ROOT compression. This leaves roughly 3 GB below the 20 GB
file target. The finer q cells address the rapidly changing low-q peak; the
finer r cells address the short-range Coulomb part. The high-q baseline and
large-r region use coarser cells. Check grid convergence before production
fits; this allocation is a numerical choice, not a demonstrated error bound.

From the DLM repository root on the farm:

```sh
make -C CATS_3D/farm/kernel estimate
make -C CATS_3D/farm/kernel smoke
df -h /scratch8/Wenya
sbatch CATS_3D/farm/kernel/run_kernel.sbatch
```

After the job completes, run `make -C CATS_3D/farm/kernel check`. This scans
all q-cell indices and samples payloads without loading the 17 GB kernel into
RAM.

The default output is `/scratch8/Wenya/pipi_coulomb_q30_r54_nonuniform.root`.
The SLURM script refuses to overwrite an existing file. Its log is written
to `/scratch8/Wenya/cats3d_kernel_<jobid>.log`. Generation uses little RAM,
but loading the complete kernel for a fit requires at least 17 GB plus ROOT
and fit overhead; request a node with at least 32 GB for subsequent fits.
The output is not committed to Git. After generation, link it from
`CATS_3D/inputs/OutputFiles/` for the fit workflow.
