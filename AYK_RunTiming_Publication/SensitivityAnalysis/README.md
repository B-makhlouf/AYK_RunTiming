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

**Question:** After the 0.7 cutoff, how much of the total stream length (km) of each
*individual fish's* assignment lies within a single regional group?

**How:** A stream edge belongs to a fish's assignment if it survives the 0.7 threshold.
Each retained edge contributes its full length (km). For every fish we sum the assigned
length in each of the six groups and take the largest group's share of the fish's total
assigned length. This is length-weighted, not posterior-weighted.

**Headline result (n = 1,465 fish):** on average **64%** of an individual's assigned
stream length falls in its single dominant regional group (median 66%; IQR 55–72%).
**84% of fish** have ≥50% of their assigned kilometres in one group, and 30% have ≥70%.
Only ~0.2% of assigned length lies outside the six groups, and the average assignment
spans ~4.3 of the 6 groups — individual assignments concentrate in one region but still
cover neighbouring stream network, which is why the paper aggregates to regional groups.

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
