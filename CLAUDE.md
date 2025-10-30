# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R-based data analysis pipeline for analyzing COVID-19 transmission patterns across geographic regions. The pipeline processes sequence metadata and pairwise distance data to calculate relative risk (RR) metrics across geographic divisions, demographics (age, sex), and time periods.

## Core Architecture

### Data Flow
1. **Database Initialization**: Raw metadata and pair data are loaded into DuckDB for efficient querying
2. **Data Cleaning**: Sequences are filtered by coverage, host type, and geographic validity; demographic data is standardized
3. **Pair Exposure Binding**: Pairwise sequences are joined with exposure variables (e.g., state, age_class)
4. **RR Matrix Calculation**: Relative risk matrices are computed across exposure categories
5. **Visualization**: Results are plotted as heatmaps, geographic maps, and line plots

### Key Tables in DuckDB
- `metadata`: Sequence information including strain ID, date, location (region/country/division), demographics (age, sex), viral clade information, and coverage metrics
- `pairs`: Pairwise comparisons of identical or near-identical sequences (strain_1, strain_2)
- `pairs_time`: Time-stratified pair comparisons for temporal analyses

### Critical Functions
The RR calculation pipeline relies on five core functions that work together:

1. **`bind_pairs_exp(db_con, exp_var, ...)`** (scripts/bind_pairs_exp.R:12)
   - Joins the `pairs` table with `metadata` to attach exposure variables (e.g., division, age_class) to both strains in each pair
   - Supports time-bounding and subsampling for bootstrapped confidence intervals
   - Returns a lazy DuckDB table with columns: strain_1, strain_2, x (exposure for strain_1), y (exposure for strain_2)

2. **`calculate_rr_matrix(df_p)`** (scripts/calculate_rr_matrix.R:10)
   - Takes output from `bind_pairs_exp` and computes RR values for all pairs of exposure categories
   - Formula: `RR = ((pair_count * N_total) + 1) / (x_appearances * y_appearances)`
   - Returns a symmetric matrix with columns: x, y, RR, N_pairs, N_x, N_y, N_total

3. **`calculate_rr_ci(db_con, exp_var)`** (scripts/calculate_rr_ci.R)
   - Calculates bootstrap confidence intervals by repeatedly subsampling metadata and recalculating RR matrices
   - Returns quantiles for RR_lower and RR_upper bounds

4. **`normalized_age_rr_diag(df_rr, group_vars = NULL)`** (scripts/age_time_RR_analysis.R:12)
   - Normalizes RR values by the mean of diagonal (within-group) RR values
   - Formula: `nRR = RR(x,y) / mean(RR(x,x), RR(y,y))`
   - `group_vars`: Optional vector of grouping variables (e.g., "date", "state") to normalize within each group separately
   - Use for measuring **preferential mixing** - how much groups mix with themselves vs others
   - Returns original dataframe with added `nRR` column

5. **`normalized_age_rr_fixed(df_rr, baseline_grp, group_vars = NULL)`** (scripts/age_time_RR_analysis.R:47)
   - Normalizes RR values against a fixed baseline group's within-group RR
   - Formula: `nRR_fixed = RR(x,y) / RR(baseline,baseline)`
   - `baseline_grp`: Reference age group (e.g., "Adult", "26-45y")
   - `group_vars`: Optional vector of grouping variables to normalize within each group separately
   - Use for measuring transmission **relative to a stable reference population**
   - Returns original dataframe with added `nRR_fixed` column

### Normalization Strategy Guide

**When to use diagonal normalization (`normalized_age_rr_diag`)**:
- Geographic analyses (state-to-state): No natural baseline state exists
- Research questions about preferential mixing patterns
- Symmetric comparisons where no group is a natural reference

**When to use fixed baseline normalization (`normalized_age_rr_fixed`)**:
- Age analyses: Use "Adult" or "26-45y" as stable working-age baseline
- Time series: Compare changing patterns against a stable reference group
- Policy questions: "How much higher was transmission in schools vs workplaces?"
- Baseline group should be large, stable, and not subject to intervention effects

**Using `group_vars` parameter**:
- `NULL` (default): Single normalization across entire dataset
- `"date"`: Normalize each time point independently (essential for time series)
- `c("date", "state")`: Normalize each date-state combination independently
- Rule: Normalize within the same stratification level to avoid mixing baselines across conditions

### Geographic Hierarchies
The codebase uses multiple geographic classification systems simultaneously:
- **division**: US states/Canadian provinces/Mexico
- **census_div**: US Census divisions (e.g., "Pacific", "Mountain", "New England")
- **census_reg**: US Census regions (West, Midwest, South, Northeast)
- **bea_reg**: Bureau of Economic Analysis regions (8 regions)

All geographic mappings are defined in the `clean_data.R` mutation pipeline (scripts/clean_data.R:68-118).

### Age Binning Strategy
Age aggregation depends on scenario (scripts/clean_data.R:32-49):
- **USA or CAM scenarios**: Individual year bins (00y, 01y, ..., 98y, 99y+)
- **Other scenarios**: 5-year bins by default (00-04y, 05-09y, ..., 90y+)

## Common Development Tasks

### Running the Full Pipeline
```bash
# On SLURM cluster (Fred Hutch)
sbatch batch_analysis.sh

# Locally with Snakemake
snakemake --profile ./profile --cores 32 --group-components duckdb_acc=1
```

The `--group-components duckdb_acc=1` flag ensures DuckDB operations run on a single node to avoid database locking issues.

### Running Individual Analyses

All R scripts use argparse for command-line arguments. Common pattern:
```bash
Rscript ./scripts/<script_name>.R --scenario <scenario> --ci <TRUE|FALSE>
```

Key scripts:
- `clean_data.R`: Filter and standardize metadata (required first step after DB init)
- `state_analysis.R`: Calculate RR matrix for geographic divisions
- `age_analysis.R`: Calculate RR matrix for age classes
- `age_state_analysis.R`: Joint RR matrix for age × state combinations
- `age_time_RR_analysis.R`: Time-stratified age RR analyses with normalization (produces RR, nRR, nRR_fixed)
- `age_time_plot.R`: Visualize time series of age-stratified RR metrics
- `school_share_analysis.R`: Correlate school instruction modality with transmission rates; classifies states into trajectory types (Persistent High, Persistent Low, Rising In-Person)
- `age_school_analysis.R`: Academic year boxplots and ANOVA for school-age transmission

### Modifying Scenarios
Edit `config.yaml` to add/remove analysis scenarios. Scenarios typically correspond to geographic subsets or sampling strategies (e.g., "CAM_100", "CAM_1000", "USA").

### Database Operations
```bash
# Initialize database for a scenario
./scripts/init_db.sh -s <scenario>

# Query database directly
duckdb db_files/db_<scenario>.duckdb
```

### Testing with Smaller Datasets
Use `scripts/init_db_mini.sh` to create a downsampled database for rapid testing without reprocessing full datasets.

## Code Patterns

### Standard R Script Template
```r
library(argparse)
library(dplyr)
library(duckdb)
library(dbplyr)

source("scripts/calculate_rr_matrix.R")
source("scripts/bind_pairs_exp.R")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Scenario name')
  parser$add_argument('--ci', type = 'logical', default = TRUE, help = "Calculate CIs")
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

fn_db <- paste0("db_files/db_", scenario, ".duckdb")
con <- DBI::dbConnect(duckdb(), fn_db)

# Analysis code here...

DBI::dbDisconnect(con, shutdown=TRUE)
```

### Adding a New Stratification Variable

1. Add variable computation in `clean_data.R` (after line 126)
2. Create analysis script following the pattern:
   ```r
   df_rr <- con %>%
     bind_pairs_exp("<new_var>") %>%
     calculate_rr_matrix() %>%
     collect()

   if(ci_flag){
     df_rr_ci <- calculate_rr_ci(con, "<new_var>")
     df_rr <- inner_join(df_rr, df_rr_ci, by=join_by(x,y))
   }
   ```
3. Add corresponding Snakemake rule in `Snakefile`

## File Organization

### Input Data (`data/`)
- `<scenario>/metadata/metadata_<scenario>.tsv.zst`: Sequence metadata (compressed)
- `<scenario>/distance_aggregated/combined_df_identical_pairs_<scenario>.tsv.zst`: Pairwise sequence comparisons
- Geographic reference data (state distances, CBSA distances, travel matrices)

### Outputs
- `results/<scenario>/df_RR_by_<stratification>.tsv`: RR matrices (may include RR, nRR, nRR_fixed columns)
- `results/<scenario>/time_age/`: Time-stratified age RR analyses and school-related outputs
  - `df_RR_by_time_age_series.tsv`: Rolling window age RR time series (2-month windows)
  - `df_RR_by_school_state_time.tsv`: School-age RR by state and time
  - `df_RR_by_school_state_ay.tsv`: School-age RR by state and academic year
  - `state_trajectory_classifications.tsv`: State classification by in-person instruction trajectory
- `results/<scenario>/summary_tables/`: Demographic summary tables (compressed tar.zst)
- `figs/<scenario>/`: Visualizations (JPG/PNG format)
- `db_files/db_<scenario>.duckdb`: DuckDB databases

## Important Notes

- **DuckDB connections**: Always disconnect with `shutdown=TRUE` to prevent file locks
- **Memory management**: Use `dbplyr` lazy evaluation; only `collect()` when necessary
- **Time analyses**: Use `pairs_time` table and pass `time_bounds` parameter to `bind_pairs_exp`
  - **Rolling time windows**: Time-stratified RR analyses use rolling windows with `MONTH_BUFFER = 1` (scripts/age_time_RR_analysis.R:118)
  - This creates 2-month windows (±1 month = 28 days around each mid-point)
  - Mid-points are spaced every 4 weeks from 2020-03-01 to 2024-10-01
  - **Academic year analyses** use fixed periods (Sept 1 to May 31) without buffering
- **NA handling**: Age and sex have explicit cleaning functions; "NA" may appear as character literal in DuckDB
- **Module loading**: Scripts assume `ml fhR` module is loaded (Fred Hutch environment)
- **Color schemes**: ALWAYS use `scripts/color_schemes.R` for visualization colors to maintain consistency across all figures. This file defines standardized color scales for regions (BEA regions), countries, RR gradients, sex, and other categorical variables.

## Dependencies

### R Packages
- tidyverse (dplyr, ggplot2, readr)
- data.table
- dbplyr
- duckdb
- argparse
- RColorBrewer, usmap (visualization)

### System Tools
- Snakemake (7.18.2+)
- DuckDB CLI
- zstd (compression)
