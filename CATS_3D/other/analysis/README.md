# Analysis scripts

- `bootstrap/`: bootstrap summaries and plots.
- `plots/`: fit, baseline, lambda, and resolution figures.
- `diagnostics/`: small model, smearing, and performance studies.
- `preparation/`: data preparation and resolution macros.

Run ROOT macros from `other/` because their default input and output paths
resolve through its local compatibility links. Python plots resolve the
workspace root from their file paths.
