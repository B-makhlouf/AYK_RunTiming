# 2. Threshold sensitivity of the final results

Re-runs the cumulative-distribution and front-end closure-protection analysis at
assignment thresholds 0.5, 0.6, 0.7 (paper baseline), 0.8, 0.9.

## Figures
- `figures/closure_protection_vs_threshold.png` — mean closure protection per group vs
  threshold (shaded = year-to-year min–max). Lines are nearly flat → robust.
- `figures/closure_protection_bars_by_threshold.png` — same numbers as grouped bars
  (baseline 0.7 outlined in black), easiest side-by-side read.
- `figures/group_cumulative_distribution_by_threshold.png` — group cumulative run-timing
  curves; the five threshold curves sit almost on top of each other.

## Data
- `data/closure_protection_mean_by_group_x_threshold.csv` — compact group × threshold
  table of mean protection %.
- `data/closure_protection_summary_by_group_and_threshold.csv` — mean/min/max per group/threshold.
- `data/closure_protection_by_year_and_threshold.csv` — full per-year detail.
- `data/group_cumulative_distribution_by_threshold.csv` — cumulative % by group/DOY/threshold.

## Reading it
Mean protection moves ≤ ~2 percentage points across the full 0.5–0.9 range for every
group, and the ranking (Group 1 > Group 2 > … > Groups 4–6) is stable. The paper's
conclusion does not depend on the choice of 0.7.
