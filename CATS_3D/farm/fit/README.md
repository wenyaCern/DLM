# Farm fit

This folder builds the same fit source used in `../../fit/`. Local single-bin
runs test the workflow; farm runs will provide production nominal, bootstrap,
and systematic inputs to that shared model. The old independent fit snapshot
and submission scripts are archived under
`../../other/archive/farm_fit_snapshot_2026-11-19/` for reference.

Local farm results remain in `farm_output/` and are excluded from Git.

All active input reads are under `../../inputs/`: `data/` for CF, lambda, and
smearing files, `OutputFiles/` for kernels, and `EPOS_data/` for resonance
inputs. A large farm kernel can physically live on scratch storage and be
linked from `inputs/OutputFiles/`. Copy `paths.example.sh` to a local
`paths.sh`, update its filenames, and source it before running from
`../../other/build/`. The old farm path convention is preserved only in the
archived snapshot. The final fine kernel path will be set when its grid is
chosen.
