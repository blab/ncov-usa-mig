# AI Assistant Instructions for ncov-usa-mig

This document guides AI agents working with the ncov-usa-mig codebase, which analyzes COVID-19 transmission patterns across geographic regions.

## Project Overview

This is an R-based data analysis pipeline that:
1. Processes sequence metadata and pair distances
2. Calculates relative risk (RR) metrics across geographic divisions 
3. Generates visualizations and statistical analyses
4. Manages data through DuckDB for efficient processing

## Key Architecture Components

### Data Pipeline Structure
- `data/` - Raw input data including metadata, distances, and geographic reference data
- `scripts/` - R analysis scripts organized by analysis type
- `results/` - Generated data files and intermediate results 
- `figs/` - Generated plots and visualizations
- `db_files/` - DuckDB databases for efficient data querying

### Core Processing Flow
1. Database initialization (`init_db.sh`) with metadata and pair data
2. Data cleaning and standardization (`clean_data.R`)
3. Analysis across multiple dimensions:
   - Geographic (state, census division)
   - Demographic (age groups)
   - Temporal trends

## Development Workflow

### Environment Setup
```bash
# Required R packages
library(tidyverse)
library(data.table)
library(dbplyr)
library(duckdb)
library(ggplot2)
```

### Running Analyses
The project uses Snakemake for workflow management:
```bash
snakemake --profile ./profile --cores <n> --group-components duckdb_acc=1
```

### Key Configuration
- Analysis scenarios defined in `config.yaml`
- Geographic reference data in `data/regions.csv`
- Database configurations in `db_files/`

## Code Patterns

### Data Loading Pattern
```r
fn_db <- paste0("db_files/db_", scenario, ".duckdb")
con <- DBI::dbConnect(duckdb(), fn_db)
df_meta <- tbl(con, "metadata")
```

### Analysis Functions
Key utility functions in `/scripts/`:
- `calculate_rr_matrix.R` - Core relative risk calculations
- `bind_pairs_exp.R` - Pair data processing
- `calculate_rr_ci.R` - Confidence interval calculations

### File Naming Conventions
- Analysis scripts: `<analysis_type>_<metric>.R`
- Results: `df_RR_by_<stratification>.tsv`
- Figures: `<analysis_type>_<plot_type>.jpg`

## Integration Points

### Database Schema
- metadata: Sequence and demographic information
- pairs: Pairwise sequence comparisons
- Geographical mappings (divisions, regions)

### External Dependencies
- Geographic shapefiles for mapping
- Population data by region
- Census division reference data

## Common Operations

1. Adding a new analysis:
   - Create R script in `scripts/`
   - Add Snakemake rule in `Snakefile`
   - Define outputs in `results/` and/or `figs/`

2. Data processing:
   - Use DuckDB for large dataset operations
   - Follow the bind_pairs_exp → calculate_rr_matrix pattern
   - Save results in standardized TSV format

3. Visualization:
   - Use ggplot2 with consistent themes
   - Save to `figs/<scenario>/` directory
   - Follow existing naming patterns