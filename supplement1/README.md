# Supplement S1 — the read-out on known networks

Small simulation study behind Supplement S1 of the bridge-node paper: *given
the true network, how does the bridge score do?* Networks are drawn from the
paper's planted MMSBM, the thesis curvature algorithm is run on each known
network, and the corrected (and raw) bridge scores are compared with the
continuous target b_true (manuscript Eq. B1–B4). No bgms, no data generation —
estimation is the main study's job, not S1's.

## Design

The grid mirrors the thesis's simulation conditions, 200 networks per cell:

- **Part A — separability × bridge strength** (balanced 10+10): the thesis's
  five separability levels (p_within, p_between) ∈ {(.8, 0), (.8, .05),
  (.75, .10), (.7, .15), (.6, .30)} × π ∈ {(.5,.5), (.7,.3), (.9,.1)} = 15 cells.
- **Part B — community imbalance** at baseline separability (.8, .05):
  splits {13+7, 15+5} × the three π = 6 cells.

21 conditions, 4,200 networks. Detection uses T_iter = 8, δ_prune = 0.1,
automatic δ_cut (the thesis settings). The whole-count discretization floor
(Appendix B's claim) is computed at the base condition (10+10, π = (.5,.5),
baseline separability).

**Convention:** community 1 is the larger community and π[1] is the bridge's
membership in it, so an asymmetric bridge leans toward the larger community.
Swap the entries of `pi_bridge` in `s1_run.R` to flip this.

## Files

| file | contents |
|---|---|
| `curvature_CD.R` | thesis algorithm, **verbatim** (line-graph pipeline from `algorithm_implementations/curvature_CD.R` in the thesis repo); do not edit |
| `ricci_flow_tests.R` | thesis test suite, **verbatim** (`Tests/testthat/ricci_flow_tests.R`; only the `source()` path adapted) — 21 hand-verified ORC / Ricci-flow / detection tests |
| `s1_functions.R` | generator, b_true, corrected + raw scores, whole-count floor |
| `s1_run.R` | the 21 × 200 grid; per-condition CSVs, resumable |
| `s1_figures.R` | 5 figures + summary tables: tracking, separability curve, invariance from the true structure (no detection), invariance under detection, whole-count floor |
| `s1_tests.R` | correctness anchors for the S1-specific code (Table B1, worked example, floor) |
| `s1_sensitivity.R` | flow-parameter sensitivity at the base condition: T ∈ {4, 16, 30}, δ_prune ∈ {.05, .2}, same 200 networks; T = 8 reference read from `results/cond_2.csv` |

## Run

```sh
Rscript ricci_flow_tests.R         # thesis algorithm tests (21x "Test passed")
Rscript s1_tests.R                 # S1 anchors — must print "all tests passed"
Rscript s1_run.R                   # full grid, ~5–7 h single core
Rscript s1_figures.R               # results/fig_s1_*.pdf, s1_summary.csv
Rscript s1_sensitivity.R           # optional: flow-parameter sensitivity
                                   # (needs results/cond_2.csv from s1_run.R)
```

`S1_CORES=8 Rscript s1_run.R` parallelizes over replications (~1 h on an
M-series MacBook; results are bit-identical — every replication has its own
seed). `S1_REPS=20 Rscript s1_run.R` for a quick pass. Each condition writes
`results/cond_<i>.csv` and is skipped when the file exists; delete a file to
re-run that condition. Seeds are `10000*condition + rep`, fully deterministic.

## Dependencies

R with `igraph`, `lpSolve`, `transport` (what the thesis code loads),
`testthat` (thesis tests only) and `ggplot2` (figures only).
