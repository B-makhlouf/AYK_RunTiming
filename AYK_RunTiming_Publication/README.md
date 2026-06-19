# Kuskokwim Chinook Run-Timing — Reproducible Analysis

Code to reproduce the figures and underlying results in:

> *Strategic fisheries closure balances harvest opportunity and stock diversity in a large Alaska river basin.* Makhlouf, Schindler, Cline, Fernandez, Esquible, Hoffman. **Canadian Journal of Fisheries and Aquatic Sciences (CanJFS).**

This is a cleaned, self-contained version of the original analysis. The statistical
code (the Bayesian natal-origin assignment, the k-means regional grouping, and the
six cluster analyses) is carried over **verbatim** from the working repository — only
file paths were centralized and dead code removed (see *Relationship to the original
code* below). It is intended to produce the same results as the paper when run against
the same input data.

---

## What it does

The pipeline reconstructs where and when Chinook salmon return to the Kuskokwim
basin from otolith strontium isotopes, groups tributaries into six regional stock
groups, and evaluates how the June 1–11 front-end fishery closure protects each
group.

## Requirements

R (≥ 4.1) and these packages:

```r
install.packages(c(
  "sf", "dplyr", "tidyr", "readr", "ggplot2", "scales", "RColorBrewer",
  "viridis", "glue", "gridExtra", "cowplot", "patchwork",
  "cluster", "here"
))
```

## Repository layout

```
AYK_RunTiming_Publication/
├── AYK_RunTiming_Publication.Rproj   open this in RStudio first
├── config.R          ← ALL paths + shared parameters (the one file you edit)
├── run_all.R         ← runs the whole pipeline end to end
├── R/
│   ├── 00_setup.R                      parameters + Bayesian assignment engine
│   ├── 01_plot_functions.R             shared plotting functions
│   ├── 02_tributary_management.R       natal origin → tributary/management production
│   ├── 03_cumulative_distribution.R    seasonal cumulative distribution + run duration
│   ├── 04_cpue_remaining.R             daily CPUE & remaining proportion
│   ├── 05_clustering_regional_groups.R k-means → 6 regional groups (+ diagnostics)
│   ├── 06_cluster_analysis.R           cluster analyses 1–6          (Figs 3–6)
│   └── 07_export_manuscript_figures.R  copy final figures → Figures/ in paper order
├── data/
│   ├── natal_origins/   six per-fish CSVs (2017–2022) — INCLUDED
│   └── spatial/         place the two shapefiles here — YOU SUPPLY (see below)
├── Figures/            final deliverables (filled when you run the pipeline)
│   ├── Figure3_*.png … Figure6_*.png      manuscript figures, paper-ordered
│   └── Clustering/                        cluster-selection + isotope boxplot figures
└── outputs/            all working outputs (figures, CSV tables, intermediates)
```

## Input data

**Included** — `data/natal_origins/2017…2022_Kusko_Natal_Origins_Genetics_CPUE.csv`:
one row per otolith, with the modeled natal isotope ratio (`natal_iso`,
`natal_isose`), day of year (`DOY`), and CPUE weighting (`dailyCPUEprop`, `COratio`).

**You must supply** — the two Kuskokwim spatial layers, placed under
`data/spatial/` using the **same folder structure as the original analysis**
(see `data/spatial/README.txt`):

| Path (under `data/spatial/`) | Contents |
|------|----------|
| `Spatial Data/KuskoUSGS_HUC_joined.shp` | stream edges with `iso_pred`, `isose_pred`, `Str_Order`, `mgmt_river`, `UniPh2oNoE`, `SPAWNING_C`, `Spawner_IP` |
| `Spatial Data/Kusko_basin.shp` | basin outline polygon |

Both shapefiles (with their `.shp`, `.shx`, `.dbf`, `.prj`, … sidecar files) are
in place. They come from the French et al. (2020) Kuskokwim strontium isoscape
and were not redistributed with the original repository. If you move them, point
`KUSKO_EDGES` / `KUSKO_BASIN` in `config.R` at the new location.

## How to run

1. Open `AYK_RunTiming_Publication.Rproj` in RStudio (this sets the working directory
   so all relative paths resolve), or `setwd()` to this folder.
2. Put the two shapefiles in `data/spatial/`.
3. Run the whole pipeline:

   ```r
   source("run_all.R")
   ```

Everything is written under `outputs/`. Stages run in dependency order; later stages
read the CSV intermediates written by earlier stages, so run them through `run_all.R`
(or source the scripts in the numbered order). Individual stages can be re-run after
`source("config.R")`.

## Figure → output map

The pipeline produces the four data figures in the manuscript (Figures 3–6).
**Figures 1, 2A and 2B are GIS-made basemaps (ArcGIS + Esri basemap) and are not
generated here** — they are drawn separately from the exported cluster /
tributary-group assignments.

After a full run, `R/07_export_manuscript_figures.R` copies the final figures into
a clean, paper-ordered **`Figures/`** folder. These copies are byte-identical to
the analysis-script outputs (same data, same rendering):

**`Figures/`** — manuscript figures, named in paper order:

| Manuscript figure | `Figures/` file | Produced by | Source file (in `outputs/`) |
|---|---|---|---|
| Figure 3 — proportional contribution by quartile (boxplots) | `Figure3_proportional_contribution_by_quartile.png` | `R/06` Analysis 1 | `Cluster Analysis Results/Four_Panel_Boxplots_2017-2022.png` |
| Figure 4 — cumulative distribution by regional group | `Figure4_cumulative_distribution_by_group.png` | `R/06` Analysis 2 | `Cluster Analysis Results/cluster_cumulative_distribution_MULTIPANEL.png` |
| Figure 5 — closure protection across timing offsets | `Figure5_closure_protection_by_offset.png` | `R/06` Analysis 5 | `Cluster Analysis Results/cluster_front_end_closure_protection_comparison.png` |
| Figure 6 — protection vs. foregone harvest | `Figure6_protection_vs_foregone_harvest.png` | `R/06` Analysis 6 | `Cluster Analysis Results/cpue_vs_protection_combined.png` |

Figure identities were verified directly against the embedded images in
`CanJFS_AllFigures.docx` (e.g. Figure 4 is the per-year regional-group MULTIPANEL,
**not** the management-unit overview).

**`Figures/Clustering/`** — clustering-results figures (from `R/05`):

| Clustering figure | `Figures/Clustering/` file | Source file (in `outputs/`) |
|---|---|---|
| Cluster-selection diagnostic — 4-panel: WSS elbow, **silhouette score**, variance explained, rate of improvement (k = 2–10) | `cluster_selection_elbow_silhouette.png` | `Management_Distinguishability/elbow_method_analysis_combined_upper_kusko.png` |
| **Boxplots** of isotope values per cluster (k = 6) | `isotope_value_boxplots_by_cluster.png` | `Management_Distinguishability/cluster_isotope_boxplots_k6_combined_upper_kusko.png` |

The numeric backing for the clustering figures (per-k WSS / between-SS /
silhouette in `elbow_method_metrics_combined_upper_kusko.csv`, and the k = 6
`cluster_assignments_…csv`) stays in `outputs/Management_Distinguishability/`.

Everything else the original scripts emitted (per-year/quartile production maps,
run-duration plot, CPUE-curve figure, the k = 3/4/5 maps and histograms, the R
basin map + violin, the management-unit cumulative overview, days-to-50% and
single-window closure plots, and the faceted/annotated/horizontal variants) is
intentionally **not exported** — each call is commented in place with a
`# [trimmed: …]` marker, while the `*.csv` data tables behind the kept figures are
still written (to `outputs/Cluster Analysis Results/CSV_Exports/` and the other
`outputs/` subfolders).

## Key parameters (in `config.R` / `R/00_setup.R`, unchanged from the paper)

- Years: 2017–2022; basin: Kuskokwim.
- Cumulative distribution sampled every **3 days**.
- Front-end closure window: **DOY 153–162 (June 1–11)**; offsets −3/+3/+6 days.
- k-means: `set.seed(42)`, `nstart = 25`, evaluated for k = 2–10, **k = 6 selected**.
- Bayesian assignment error model and stream-order / habitat priors are defined in
  `R/00_setup.R` (`calculate_error`, `setup_priors`, `perform_assignment`).

## Relationship to the original code

The seven scripts in `R/` are derived from the originals with only these changes,
none of which alters any computation, parameter, or random seed:

- **Paths centralized** — every hardcoded absolute path (e.g.
  `/Users/benjaminmakhlouf/...`, and the inconsistent `AYK_RunTiming` vs
  `02_AYK_RunTiming` roots) now resolves through `config.R`.
- **Outputs trimmed** — `ggsave()`/`png()` calls for figures not used in the paper
  are commented out in place (marked `# [trimmed: …]`); the k = 6 map, histogram,
  and violin are gated with `if (k == 6)`. All underlying computations and CSV
  exports are untouched, so the kept figures are bit-for-bit unchanged.
- **Dead code removed** — an RStudio-only `rstudioapi::getActiveDocumentContext()`
  block at the end of `00_setup.R` that sourced a non-existent `00b_cluster_mapping.R`
  was deleted (it errored in non-interactive use and produced no output).
- **Orchestration added** — `run_all.R` sources the stages in dependency order and
  calls the existing `run_*` functions.

The analytical code itself (the Bayesian assignment, k-means clustering, and all six
cluster analyses) is unchanged from the originals.

Scripts not part of the published figures (presentation graphics, the portfolio-effect
figures, the DFA and other archived analyses, and the non-cluster closure/precision
variants) were intentionally left out of this clean repository.

Source scripts these were derived from, in the original `Code/` folder:
`00_setup.R`, `05_visualization_functions.R`, `01_tributary_maps.R`,
`03_cumulative_distribution.R`, `CPUEremainingTiming.R`, `clustering.R`,
`ClusterAnalysisCombining.R`.
