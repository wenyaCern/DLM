# Farm fit

This folder builds the same fit source used in `../../fit/`. Local single-bin
runs test the workflow; farm runs will provide production nominal, bootstrap,
and systematic inputs to that shared model. The old independent fit snapshot
and submission scripts are archived under
`../../other/archive/farm_fit_snapshot_2026-11-19/` for reference.

Local farm results remain in `farm_output/` and are excluded from Git.

The Gaussian full-range runner is `run_full_range.sh`. It fits one
`(mT, qn, phi)` cell at a time over the entire Cartesian
`|q_out|, |q_side|, |q_long| < 0.2 GeV/c` box. It uses the same baseline
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
`CATS_3D/farm/fit/farm_output/cent0to10_fullQ_afterDCA_weak155/`.
The full kernel is over 20 GB in memory; benchmark this one cell before
launching further cells. The preflight `--check` verifies input paths and
executable presence without loading the kernel. The old farm path convention
is preserved only in the archived snapshot.
