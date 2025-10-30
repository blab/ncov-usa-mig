# COVID-19 Transmission Pattern Analysis Pipeline

> **⚠️ Status: Under Active Development**
>
> This analysis pipeline is currently being developed and refined. Methods, results, and documentation are subject to change. For detailed technical documentation, see [CLAUDE.md](CLAUDE.md).

## Overview

This repository contains an R-based data analysis pipeline for investigating COVID-19 transmission patterns across geographic regions and demographic groups in North America. The pipeline analyzes SARS-CoV-2 genomic sequences to identify identical or near-identical sequence pairs and calculates relative risk (RR) metrics for transmission between different populations.

## Key Features

### Implemented Analyses

#### Geographic Transmission Patterns
- **State/Division-level RR matrices**: Quantifies transmission risk between US states, Canadian provinces, and Mexican regions
- **Census division and region analyses**: Aggregated transmission patterns across US Census geographic hierarchies
- **Time-stratified geographic RR**: Rolling window analysis (2-month windows, ±1 month) to track changing transmission patterns over time
- **Network visualization**: Identification and visualization of persistent inter-state transmission connections

#### Demographic Transmission Patterns
- **Age-stratified RR matrices**: Fine-grained age bins (individual years or 5-year bins depending on dataset)
- **Age × State joint analyses**: Simultaneous stratification by age and geography
- **Sex-stratified analyses**: Transmission patterns by biological sex
- **Age-time series**: Temporal dynamics of age-specific transmission with normalization

#### School-Related Analyses
- **School instruction modality correlation**: Links school in-person/hybrid/virtual instruction shares with age-specific transmission rates
- **State trajectory classification**: Classifies states into three categories based on in-person instruction patterns during the 2020-2021 academic year:
  - **Persistent High** (6 states): Maintained ≥80% in-person throughout
  - **Persistent Low** (6 states): Maintained ≤20% in-person throughout
  - **Rising In-Person** (3 states): Transitioned from ≤20% to ≥80% during the year
- **Academic year analyses**: School-age transmission RR by state across multiple academic years (2020-2024)

#### Statistical Methods
- **Normalized RR metrics**:
  - `nRR` (diagonal normalization): Measures preferential mixing within groups relative to between-group mixing
  - `nRR_fixed` (baseline normalization): Transmission relative to a stable reference population (e.g., working-age adults)
- **Bootstrap confidence intervals**: Subsampling-based uncertainty quantification
- **Mixed-effects modeling**: Accounts for temporal autocorrelation and state-level random effects in school analyses

## Analysis Approach

### Data Pipeline

1. **Sequence Data Processing**
   - Raw SARS-CoV-2 genomic sequences with metadata (location, date, age, sex, viral clade)
   - Pairwise distance calculations identify identical or near-identical sequences
   - Sequences filtered by coverage (≥90%), host type (human), and geographic validity

2. **Relative Risk Calculation**
   - Pairs of identical sequences represent potential transmission events
   - RR quantifies whether pairs occur more or less frequently than expected between exposure groups
   - Formula: `RR = ((pair_count × N_total) + 1) / (x_appearances × y_appearances)`
   - Pseudocount (+1) prevents division by zero and stabilizes rare category estimates

3. **Temporal Stratification**
   - Rolling time windows: 2-month windows centered on mid-points spaced every 4 weeks
   - Academic year periods: Fixed Sept 1 - May 31 windows for school analyses
   - Time-bounding enables tracking of transmission dynamics across pandemic phases

4. **Normalization**
   - Diagonal normalization controls for overall transmission intensity changes over time
   - Baseline normalization enables comparison against a stable reference group
   - Group-specific normalization (by date, state, etc.) prevents confounding across strata

## Repository Structure

```
ncov-usa-mig/
├── data/                           # Input data (sequence metadata, pairs, reference data)
│   ├── metadata/                   # Compressed sequence metadata
│   ├── distance_aggregated/        # Pairwise sequence comparisons
│   └── state_school_share.csv      # School instruction modality by state/month
├── scripts/                        # Analysis scripts
│   ├── clean_data.R                # Data cleaning and standardization
│   ├── bind_pairs_exp.R            # Pair-exposure variable joining
│   ├── calculate_rr_matrix.R       # Core RR calculation function
│   ├── age_analysis.R              # Age-stratified RR
│   ├── state_analysis.R            # Geographic RR
│   ├── age_time_RR_analysis.R      # Time-stratified age RR with normalization
│   ├── school_share_analysis.R     # School modality correlation & trajectory classification
│   └── *_plot.R / *_heatmap.R      # Visualization scripts
├── results/                        # Analysis outputs
│   ├── df_RR_by_*.tsv              # RR matrices
│   └── time_age/                   # Time-stratified school analyses
├── figs/                           # Visualizations
├── db_files/                       # DuckDB databases for efficient querying
├── Snakefile                       # Workflow automation
├── config.yaml                     # Configuration settings
├── CLAUDE.md                       # Detailed technical documentation
└── README.md                       # This file
```

## Key Outputs

### RR Matrices
- `df_RR_by_age_class.tsv`: Age-stratified relative risk
- `df_RR_by_state.tsv`: State/division-stratified relative risk
- `df_RR_by_age_state.tsv`: Joint age × state relative risk
- `df_RR_by_census_div.tsv`: Census division-level relative risk

### Time-Stratified Analyses
- `time_age/df_RR_by_time_age_series.tsv`: Age RR time series (rolling windows)
- `time_age/df_RR_by_school_state_time.tsv`: School-age RR by state and time
- `time_age/df_RR_by_school_state_ay.tsv`: School-age RR by academic year
- `time_age/state_trajectory_classifications.tsv`: State in-person instruction trajectories

### Visualizations
- Heatmaps: Symmetric RR matrices with dendrograms
- Geographic maps: State-level transmission patterns on US maps
- Time series plots: Temporal dynamics with confidence intervals
- Network visualizations: Persistent inter-state transmission connections
- Trajectory plots: School instruction modality by state classification

## Requirements

### Software Dependencies
- **R (≥4.0)** with packages:
  - tidyverse (dplyr, ggplot2, readr, tidyr, purrr, stringr)
  - DuckDB (efficient database queries)
  - dbplyr (lazy evaluation)
  - data.table (fast data manipulation)
  - lme4/lmerTest (mixed-effects models)
  - argparse (command-line interfaces)
  - RColorBrewer, usmap (visualization)
- **Snakemake (≥7.18)**: Workflow management
- **DuckDB CLI**: Database initialization
- **zstd**: Data compression/decompression

### System Requirements
- SLURM cluster environment (for large-scale analyses)
- Sufficient memory for genomic data processing (recommend ≥32GB RAM)
- Multi-core processor for parallel processing

## Usage

### Quick Start

```bash
# 1. Initialize database
./scripts/init_db.sh

# 2. Clean and standardize metadata
Rscript ./scripts/clean_data.R

# 3. Run a specific analysis
Rscript ./scripts/age_analysis.R --ci TRUE

# 4. Generate visualizations
Rscript ./scripts/age_heatmap.R
```

### Full Pipeline Execution

```bash
# Run complete analysis pipeline via Snakemake
snakemake --profile ./profile --cores 32 --group-components duckdb_acc=1

# Or submit to SLURM cluster
sbatch batch_analysis.sh
```

### School Trajectory Analysis

```bash
# Generate school-age RR data (requires age_time_RR_analysis.R to run first)
Rscript ./scripts/school_share_analysis.R
```

This will:
- Classify states by in-person instruction trajectory
- Generate correlation plots between school modality and transmission
- Create faceted trajectory visualization
- Output state classifications to `time_age/state_trajectory_classifications.tsv`

## Recent Development Highlights

### 2024 - 2025
- ✅ Implemented normalization functions (diagonal and baseline)
- ✅ Added school instruction modality correlation analyses
- ✅ Created state trajectory classification system
- ✅ Refined time windows from 6-month to 2-month rolling windows
- ✅ Added academic year fixed-period analyses
- ✅ Implemented country-stratified age RR time series
- ✅ Enhanced visualization suite with faceted trajectory plots

### Ongoing Work
- 🔄 Validation of normalization approaches for different research questions
- 🔄 Sensitivity analyses for time window parameters
- 🔄 Integration with policy/intervention timelines
- 🔄 Network analysis of persistent transmission routes

## Documentation

- **[CLAUDE.md](CLAUDE.md)**: Comprehensive technical documentation for developers
  - Core architecture and data flow
  - Function reference with line numbers
  - Normalization strategy guide
  - File organization and naming conventions
  - Development patterns and best practices

## Citation

**Note**: Methods and results are preliminary. Formal publication is in preparation.

## Contact

For questions about this analysis pipeline, please contact the repository maintainers.

## License

[License information to be added]

---

**Last Updated**: October 2025
