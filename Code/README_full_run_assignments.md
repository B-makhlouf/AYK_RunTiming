# Kuskokwim Full Run Natal Origin Assignments

## Overview

`full_run_kuskokwim_assignments.R` is a single, streamlined script that processes all years of natal origin data (2017-2022) and generates **FULL RUN** assignments for the Kuskokwim watershed. The script exports scaled assignments that sum to 1, representing production proportions across the basin.

## Key Features

- **Single script**: All processing in one well-organized file
- **No functions**: Inline code for maximum readability
- **Full run**: Processes all fish together (not divided by quartiles)
- **Scaled outputs**: Assignments sum to 1 (production proportions)
- **Multi-year**: Loops through all years automatically
- **Two output levels**: Edge-level and management unit aggregations

## What It Does

The script performs the following steps for each year:

1. **Load natal origin data** - Fish isotope measurements with CPUE weights
2. **Bayesian assignment** - Assigns each fish to stream edges using:
   - Isotope likelihood (Gaussian model)
   - Spatial priors (stream order, spawning habitat)
   - Isoscape predictions
3. **Aggregate across fish** - Sums assignments weighted by CPUE
4. **Scale to proportions** - Normalizes so total production = 1
5. **Management unit summary** - Aggregates edges to management units

## Usage

### 1. Update Paths

Before running, open the script and update the `BASE_DIR` variable (line 30):

```r
BASE_DIR <- "/Users/benjaminmakhlouf"  # Change to your path
```

### 2. Run the Script

From RStudio or R console:

```r
source("Code/full_run_kuskokwim_assignments.R")
```

Or from terminal:

```bash
Rscript Code/full_run_kuskokwim_assignments.R
```

### 3. Check Outputs

The script creates two CSV files in `Analysis_Results/Full_Run_Assignments/`:

1. **`kuskokwim_edge_assignments_full_run.csv`**
   - One row per stream edge per year
   - Columns:
     - `year` - Analysis year
     - `edge_id` - Stream edge identifier (1 to N)
     - `Str_Order` - Stream order
     - `mgmt_river` - Management unit name
     - `iso_pred` - Predicted isotope value
     - `assignment_sum` - Raw summed assignments
     - `assignment_proportion` - **Scaled proportion (sums to 1)**
     - `assignment_normalized` - Normalized for visualization (max = 1)

2. **`kuskokwim_management_unit_assignments_full_run.csv`**
   - One row per management unit per year
   - Columns:
     - `year` - Analysis year
     - `mgmt_river` - Management unit name
     - `production_proportion` - **Production proportion (sums to 1)**
     - `total_production` - Raw production value
     - `edge_count` - Number of stream edges in unit

## Key Parameters

These can be modified in Section 2 of the script:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `YEARS` | 2017-2022 | Years to process |
| `SENSITIVITY_THRESHOLD` | 0.7 | Assignment threshold (70%) |
| `MIN_ERROR` | 0.0006 | Minimum isotope error |
| `MIN_STREAM_ORDER` | 3 | Minimum stream order |

## Output Validation

The script prints summary statistics to verify correct scaling:

```
Total production proportion sum: 1.0000
```

This confirms that assignments sum to 1 for each year.

## Management Units

The script processes 19 Kuskokwim management units (upstream to downstream):

1. N. Fork Kusko
2. E. Fork Kuskokwim
3. S. Fork Kusko
4. Upper Kusko Main
5. Big River
6. Takotna and Nixon Fork
7. Tatlawiksuk
8. Swift
9. Stony
10. Holitna
11. Hoholitna
12. Middle Kusko Main
13. George
14. Oskakawlik
15. Holokuk
16. Aniak
17. Tuluksak
18. Kisaralik
19. Kwethluk
20. Lower Kusko (includes consolidated Johnson)

## Notes

- **Johnson consolidation**: The Johnson management unit is consolidated into Lower Kusko (consistent with other analyses)
- **Full run vs. quartiles**: Unlike other scripts that divide by DOY quartiles, this processes all fish together
- **CPUE weighting**: Each fish is weighted by its `COratio` value from the natal data
- **Missing years**: If a year's data file is not found, it is skipped with a warning

## Dependencies

Required R packages (loaded automatically):
- `sf` - Spatial data handling
- `dplyr` - Data manipulation
- `readr` - CSV reading/writing
- `glue` - String formatting

## Script Structure

The script is organized into 8 clear sections:

1. Load libraries
2. Set paths and parameters
3. Load spatial data (once)
4. Setup priors (once)
5. Storage for results
6. Main processing loop (all years)
7. Combine and export results
8. Summary statistics

## Differences from Other Scripts

| Feature | This Script | Tributary Maps (01) |
|---------|-------------|---------------------|
| Temporal division | Full run | DOY quartiles (Q1-Q4) |
| Output | Production proportions | Maps + proportions |
| Functions | Inline code | Function-based |
| Years | All in one run | Looped externally |
| Focus | Data export | Visualization |

## Example Output

**Management unit assignments (2017):**

| year | mgmt_river | production_proportion | total_production | edge_count |
|------|------------|----------------------|------------------|------------|
| 2017 | Stony | 0.0853 | 12.456 | 1450 |
| 2017 | Aniak | 0.0742 | 10.832 | 879 |
| 2017 | Kisaralik | 0.0698 | 10.198 | 1225 |
| ... | ... | ... | ... | ... |

Sum of production_proportion = 1.0000 ✓

## Contact

For questions about this script, refer to the main analysis documentation in `00_setup.R` and related analysis scripts.
