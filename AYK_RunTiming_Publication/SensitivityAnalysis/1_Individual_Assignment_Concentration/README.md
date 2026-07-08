# 1. Individual assignment concentration — stream length (threshold = 0.7)

For each individual fish: of the **total stream length (km)** its assignment covers after
the 0.7 threshold, what fraction lies within its single dominant regional group.

A stream edge is part of a fish's assignment if it survives the 0.7 cutoff (rescaled
posterior ≥ 0.7 of that fish's peak). Each retained edge contributes its full length
(`Shape_Leng`, Alaska Albers metres → km). The metric is length-weighted, not
posterior-weighted: it simply asks how much of the assigned river network sits in one group.

## Headline (n = 1,465 fish)
- Mean **64.3%**, median **65.6%** of assigned stream length is in a single group.
- IQR 54.6–72.4%; **84% of fish** have ≥50% of their assigned km in one group; 30% have ≥70%.
- Mean assigned footprint ≈ 9,836 km across ≈2,179 edges, spanning ~4.3 of the 6 groups.
- Only ~0.2% of assigned length falls outside the six groups.

## Figures
- `figures/individual_concentration_distribution.png` — histogram + cumulative distribution.
- `figures/concentration_by_group.png` — split by dominant group.

## Data (`data/`)
- `per_individual_group_concentration.csv` — one row per fish. Key columns:
  - `prop_km_in_dominant_group_of_all` — **primary metric**: dominant-group km ÷ total assigned km.
  - `prop_km_in_dominant_group_of_grouped` — same, over grouped km only (nearly identical).
  - `total_assigned_km`, `dominant_group_km`, `prop_km_outside_regional_groups`.
  - `n_groups_covered`, `n_stream_edges_retained`.
- `summary_statistics.csv`, `by_dominant_group.csv`, `by_year.csv`.
