# Farm fit

This folder builds the same fit source used in `../../fit/`. Local single-bin
runs test the workflow; farm runs will provide production nominal, bootstrap,
and systematic inputs to that shared model. The old independent fit snapshot
and submission scripts are archived under
`../../other/archive/farm_fit_snapshot_2026-11-19/` for reference.

Local farm results remain in `farm_output/` and are excluded from Git.

The Gaussian runner is `run_full_range.sh`. It fits one
`(mT, qn, phi)` cell at a time in the Cartesian
`|q_out|, |q_side|, |q_long| < 0.12 GeV/c` peak box by default. The baseline
is prefitted in the outer box and the saved prediction extends to 0.2 GeV/c.
Set `CATS3D_FIT_QMAX_MEV` to change the peak fit box. It uses the same baseline
prefit, diagonal band, k*-dependent lambda, factorized 3D smearing, and
Gaussian source as the local single-cell test. `lambdaAfterDCACut` and a
single-pion weak-decay fraction of 0.155 are the defaults for mT index 0.
The alternative `lambdaFitWindow` can be selected with
`CATS3D_LAMBDA_SELECTION=fitwindow`. For mT indices 1 and 2, set their
own `CATS3D_WEAK_FRACTION` before running.

Upload the CF, lambda, and smearing ROOT inputs into `../../inputs/data/`
with the names in `paths.example.sh`, or set absolute paths in a local
`paths.sh`. The kernel defaults to the completed farm file at
`/scratch8/Wenya/pipi_coulomb_q24_r72_nonuniform.root`. No ROOT data,
kernel, output, or `paths.sh` is committed to Git.

From the DLM repository root on the farm:

```sh
make -C CATS_3D/farm/fit gaussian
bash CATS_3D/farm/fit/run_full_range.sh 0 0 0 --check
mkdir -p CATS_3D/farm/fit/logs
screen -S cats3d_fit
bash CATS_3D/farm/fit/run_full_range.sh 0 0 0 \
  > CATS_3D/farm/fit/logs/fit_mt0_qn0_phi0.log 2>&1
```

Detach screen with Ctrl-A, D. Output is written under
`CATS_3D/farm/fit/farm_output/cent0to10_peak100_afterDCA_weak155/`.
The directory name is retained for overwriting the previous fit; the actual
fit limit is 120 MeV/c. The runner overwrites an existing cell by default.
The full kernel is about 20.6 GB on disk and is kept at full spatial
resolution in memory. The fold uses four independent accumulators per thread
and defaults to four threads (`CATS3D_FOLD_THREADS=4`) on the farm; no
spatial cells are merged. Set the variable to 1 if only one CPU is available.
Benchmark this one cell before
launching further cells. The preflight `--check` verifies input paths and
executable presence without loading the kernel. The old farm path convention
is preserved only in the archived snapshot.

The runner skips the SIMPLEX pre-search by default and keeps MIGRAD and HESSE.
The local low-mT reference cell converged to the same radii without SIMPLEX.
Set `CATS3D_USE_SIMPLEX=1` if another cell does not converge.

Run adjacent cells in one process to load the 20.6 GB kernel once:

```sh
bash CATS_3D/farm/fit/run_full_range.sh --batch 0 0 4 0 5 \
  > CATS_3D/farm/fit/logs/fit_mt0_batch.log 2>&1
```

This example fits all qn and phi cells for mT index 0. Each cell still has its
own data, baseline prefit, Minuit fit, and ROOT output. The in-memory kernel
is reused until the batch process exits. Keep the batch on one mT index because
its weak-decay fraction and radius initial value are set per mT.

The Gaussian derived-data runner uses the adaptive kernel directly. The old
external fast-kernel implementation remains available only to legacy
resonance and O-O fit targets that still reference its format.
