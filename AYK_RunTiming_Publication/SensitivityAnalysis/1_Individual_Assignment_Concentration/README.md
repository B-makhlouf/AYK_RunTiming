# 1. Individual assignment concentration (threshold = 0.7)

For each individual fish: the proportion of its 0.7-thresholded assignment that falls
in its single dominant regional group.

## Figures
- `figures/individual_concentration_distribution.png` — histogram + cumulative distribution
  across all 1,465 fish. **Mean 62%, median 62%.**
- `figures/concentration_by_group.png` — the same metric split by which group is dominant.

## Data
- `data/per_individual_group_concentration.csv` — one row per fish. Key columns:
  - `prop_in_dominant_group_of_all` — **primary metric**: dominant-group share of the fish's
    total assignment mass.
  - `prop_in_dominant_group_of_grouped` — same, but as a share of mass within the six groups
    only (nearly identical, because <0.3% of mass lands outside the groups).
  - `prop_outside_regional_groups` — share of assignment outside the six groups.
  - `n_groups_with_mass` — how many of the 6 groups received any assignment mass.
  - `n_stream_edges_retained` — number of stream edges surviving the 0.7 cutoff.
- `data/summary_statistics.csv` — overall summary.
- `data/by_dominant_group.csv`, `data/by_year.csv` — breakdowns.

## Reading it
Higher = each fish is more cleanly resolved to one region. A mean of 62% (with 80% of
fish above 50%) means individual assignments carry clear regional signal but are not
point-precise — the justification for aggregating to regional groups.
