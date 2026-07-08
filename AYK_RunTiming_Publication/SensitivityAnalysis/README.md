# Sensitivity Analysis

Two sensitivity analyses for the Kuskokwim Chinook natal-origin assignment, both
built directly on the paper's assignment engine (`R/00_setup.R` `perform_assignment`)
and its six regional groups (`R/06_cluster_analysis.R`).

The assignment engine is reimplemented in Python (`scripts/sa_core.py`) using the
exact same math and the same stream-edge attributes as the R code. It was validated
to machine precision (max abs. difference ≈ 1e-15) against the existing paper outputs:

- `outputs/Management_River_Analysis/management_river_analysis_tidy.csv` (per-unit production), and
- `outputs/Cluster Analysis Results/CSV_Exports/analysis4_cluster_closure_protection_by_year.csv` (closure protection).

So every number here is reproducible from, and consistent with, the published pipeline.

---

## 1. Individual assignment concentration  →  `1_Individual_Assignment_Concentration/`

**Question:** After the 0.7 cutoff, how much of each *individual fish's* assignment
falls in a single regional group? (i.e. is each fish resolved to one group, or smeared
across several?)

**How:** For every fish we take its 0.7-thresholded assignment across all stream edges,
map edges to the six regional groups, sum the assignment mass in each group, and take the
largest group's share of the fish's total assignment.

**Headline result (n = 1,465 fish):** on average **62%** of an individual's assignment
falls in its single dominant regional group (median 62%; interquartile range 51–69%).
**79.7%** of fish have ≥50% of their assignment in one group; 22% have ≥70%. Almost none
of the assignment mass (mean 0.24%) falls outside the six groups, so this is genuinely a
statement about concentration *among* the groups, not about lower-river leakage. On
average each fish's assignment still touches ~4.3 of the 6 groups at some level, which is
why the dominant-group share sits near 60% rather than near 100% — individual otolith
assignments are informative but not point-precise, which is exactly why the paper
aggregates to regional groups and to run-timing quartiles.

## 2. Threshold sensitivity  →  `2_Threshold_Sensitivity/`

**Question:** Would the final results change if a different cutoff than 0.7 were used?

**How:** The full cumulative-distribution and front-end-closure-protection analysis is
re-run at thresholds **0.5, 0.6, 0.7 (baseline), 0.8, 0.9**, and the per-group results
are compared.

**Headline result:** the conclusions are robust. Mean closure protection changes by at
most **~2 percentage points** for any group across the whole 0.5–0.9 range, and the core
finding is unchanged at every threshold: the early upper-river groups (Group 1 ≈ 22–23%,
Group 2 ≈ 15–17%) receive far more front-end-closure protection than the later
lower/middle groups (Groups 4–6 ≈ 5–6%). The most-to-least-protected ranking is stable,
and the group cumulative run-timing curves are nearly identical across thresholds.

---

## Reproduce

```
cd scripts
python3 part1.py && python3 part1_fig.py     # analysis 1
python3 part2.py && python3 part2_fig.py      # analysis 2
```
Requires Python with `numpy`, `pandas`, `matplotlib`, and `pyshp`. Paths point at the
repository's `data/` and `outputs/` folders.
