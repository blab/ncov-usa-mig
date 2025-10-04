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
The RR calculation pipeline relies on three core functions that work together:

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
- `results/<scenario>/df_RR_by_<stratification>.tsv`: RR matrices
- `results/<scenario>/summary_tables/`: Demographic summary tables (compressed tar.zst)
- `figs/<scenario>/`: Visualizations (JPG format)
- `db_files/db_<scenario>.duckdb`: DuckDB databases

## Important Notes

- **DuckDB connections**: Always disconnect with `shutdown=TRUE` to prevent file locks
- **Memory management**: Use `dbplyr` lazy evaluation; only `collect()` when necessary
- **Time analyses**: Use `pairs_time` table and pass `time_bounds` parameter to `bind_pairs_exp`
- **NA handling**: Age and sex have explicit cleaning functions; "NA" may appear as character literal in DuckDB
- **Module loading**: Scripts assume `ml fhR` module is loaded (Fred Hutch environment)

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
