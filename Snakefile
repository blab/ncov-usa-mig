configfile: "config.yaml"
CI_FLAG = config['ci_flag']
EXCLUDE_DUPLICATES = config['exclude_duplicates']

SCENARIOS = config["scenarios"]

# Parameters for conserved network analysis
SD_THRESHOLD = 2
MIN_OCCURRENCES = 4

# Collect final outputs for all scenarios
all_outputs = expand([
    # Summary tables
    "results/{scenario}/summary_tables/demo_tables.tar.zst",
    # Descriptive figures
    "figs/{scenario}/figure_1.jpg",
    # Age analyses
    "results/{scenario}/df_RR_by_age_class.tsv",
    "figs/{scenario}/age_heatmap.jpg",
    "figs/{scenario}/age_lineplot_all_censor.jpg",
    "figs/{scenario}/age_lineplot_by_state.jpg",
    "figs/{scenario}/age_lineplot_abridged.jpg",
    "figs/{scenario}/age_same_lineplot_by_state.jpg",
    # Age-gender sensitivity analysis
    "figs/{scenario}/age_gender/same_rr_age_gender.png",
    # State analyses
    "results/{scenario}/df_RR_by_state.tsv",
    "figs/{scenario}/state_heatmap.jpg",
    "figs/{scenario}/state_nb_dist_plot.jpg",
    "figs/{scenario}/bea_region_map.png",
    # Age x State analyses
    "results/{scenario}/df_RR_by_age_state.tsv",
    # Census division analyses
    "results/{scenario}/df_RR_by_census_div.tsv",
    "figs/{scenario}/census_divisions_heatmap.jpg",
    # Time series analyses
    "results/{scenario}/time_state/df_state_rr_series.tsv",
    "results/{scenario}/time_age/df_RR_by_time_age_series.tsv",
    # School analyses
    "results/{scenario}/time_age/state_trajectory_classifications.tsv",
    "figs/{scenario}/age_time/school_share_RR_correlation.jpg",
    # Age time visualizations
    "figs/{scenario}/age_time/USA/age_RR_time_faceted.png",
    "figs/{scenario}/age_time/all_countries/vaccine_period_RR_elderly.png",
    # State time visualizations
    "figs/{scenario}/time/state_pair_nRR_heatmap.png",
    "figs/{scenario}/clust/htree.jpg",
    # State regression
    "figs/{scenario}/dist/regression_fit.jpg",
    # Significant connections network analysis
    "results/{scenario}/time_state/df_significant_connections_3sd.tsv",
    "results/{scenario}/time_state/network_plots_complete.txt",
    # Conserved network map
    f"results/{{scenario}}/time_state/df_significant_connections_{SD_THRESHOLD}sd.tsv",
    f"results/{{scenario}}/time_state/df_conserved_connections_{MIN_OCCURRENCES}plus.tsv",
    f"figs/{{scenario}}/time_state_networks/network_conserved_{MIN_OCCURRENCES}plus_map.jpg"
], scenario=SCENARIOS)

shell.prefix("set -euo pipefail; ")

rule all:
    input: all_outputs

rule install_packages:
    output:
        touch(".snakemake/packages_installed.txt")
    shell:
        """
        Rscript scripts/install_packages.R
        """

rule database:
    input:
        ".snakemake/packages_installed.txt"
    output:
        "db_files/db_{scenario}.duckdb"
    group: "duckdb_acc"
    resources:
        duckdb_lock=1
    shell:
        """
        unzstd -f data/{wildcards.scenario}/metadata/metadata_{wildcards.scenario}.tsv.zst
        unzstd -f data/{wildcards.scenario}/distance_aggregated/combined_df_identical_pairs_{wildcards.scenario}.tsv.zst
        ./scripts/init_db.sh -s {wildcards.scenario}
        """

rule clean_data:
    input:
        ancient("db_files/db_{scenario}.duckdb")
    output:
        touch("results/{scenario}/data_cleaned.txt")
    group: "duckdb_acc"
    resources:
        duckdb_lock=1
    shell:
        """
        Rscript ./scripts/clean_data.R --scenario {wildcards.scenario}
        duckdb db_files/db_{wildcards.scenario}.duckdb < scripts/trim_pairs.sql
        Rscript ./scripts/find_duplicates.R --scenario {wildcards.scenario}
        """

rule time_prep:
    input:
        "results/{scenario}/data_cleaned.txt"
    output:
        touch("results/{scenario}/pairs_time_created.txt")
    group: "duckdb_acc"
    resources:
        duckdb_lock=1
    shell:
        """
        Rscript ./scripts/time_prep.R --scenario {wildcards.scenario}
        """

rule demo_tables:
    input:
        "results/{scenario}/data_cleaned.txt"
    output:
        "results/{scenario}/summary_tables/demo_tables.tar.zst"
    group: "duckdb_acc"
    resources:
        duckdb_lock=1
    shell:
        """
        scripts/demographic_tables.sh -s {wildcards.scenario} -c "clade_nextstrain division census_div sex age_class"
        """

rule descriptive_figures:
    input:
        "results/{scenario}/data_cleaned.txt"
    output:
        "figs/{scenario}/bea_region_map.jpg",
        "figs/{scenario}/figure_1.jpg"
    group: "duckdb_acc"
    resources:
        duckdb_lock=1
    shell:
        """
        Rscript ./scripts/descriptive_figures.R --scenario {wildcards.scenario}
        """

# ============================================================================
# Age Analyses
# ============================================================================

rule age_analysis_RR:
    input:
        "results/{scenario}/data_cleaned.txt"
    output:
        "results/{scenario}/df_RR_by_age_class.tsv"
    group: "duckdb_acc"
    resources:
        duckdb_lock=1
    shell:
        """
        Rscript ./scripts/age_analysis.R --scenario {wildcards.scenario} --ci {CI_FLAG} --exclude_duplicates {EXCLUDE_DUPLICATES}
        """

rule age_heatmap:
    input:
        "results/{scenario}/df_RR_by_age_class.tsv"
    output:
        "figs/{scenario}/age_heatmap.jpg"
    shell:
        """
        Rscript ./scripts/age_heatmap.R --scenario {wildcards.scenario}
        """

rule age_gender_sensitivity:
    input:
        "results/{scenario}/data_cleaned.txt"
    output:
        "figs/{scenario}/age_gender/same_rr_age_gender.png",
        "figs/{scenario}/age_gender/same_age_gender_boxplot.png",
        "figs/{scenario}/age_gender/rr_age_gender_corr.png",
        "figs/{scenario}/age_gender/same_gender_age_heatmap.png",
        "figs/{scenario}/age_gender/different_gender_age_heatmap.png",
        "figs/{scenario}/age_gender/gender_rr_over_age.png",
        "figs/{scenario}/age_gender/gender_rr_by_lifestage.png"
    group: "duckdb_acc"
    resources:
        duckdb_lock=1
    shell:
        """
        Rscript ./scripts/age_gender_sensitivity_analysis.R --scenario {wildcards.scenario} --exclude_duplicates {EXCLUDE_DUPLICATES}
        """

# ============================================================================
# State Analyses
# ============================================================================

rule state_analysis_RR:
    input:
        "results/{scenario}/data_cleaned.txt"
    output:
        "results/{scenario}/df_RR_by_state.tsv"
    group: "duckdb_acc"
    resources:
        duckdb_lock=1
    shell:
        """
        Rscript ./scripts/state_analysis.R --scenario {wildcards.scenario} --ci {CI_FLAG} --exclude_duplicates {EXCLUDE_DUPLICATES}
        """

rule state_heatmap:
    input:
        "results/{scenario}/df_RR_by_state.tsv"
    output:
        "figs/{scenario}/state_heatmap.jpg"
    shell:
        """
        Rscript ./scripts/state_heatmap.R --scenario {wildcards.scenario}
        """

rule state_dist:
    input:
        "results/{scenario}/df_RR_by_state.tsv"
    output:
        "figs/{scenario}/state_nb_dist_plot.jpg"
    shell:
        """
        Rscript ./scripts/state_dist_plots.R --scenario {wildcards.scenario}
        """

rule plot_region_map:
    input:
        ".snakemake/packages_installed.txt"
    output:
        "figs/{scenario}/bea_region_map.png"
    shell:
        """
        Rscript ./scripts/plot_region_map.R --scenario {wildcards.scenario}
        """

# ============================================================================
# Age x State Analyses
# ============================================================================

rule age_state_RR:
    input:
        age="results/{scenario}/df_RR_by_age_class.tsv",
        state="results/{scenario}/df_RR_by_state.tsv"
    output:
        "results/{scenario}/df_RR_by_age_state.tsv"
    group: "duckdb_acc"
    resources:
        duckdb_lock=1
    shell:
        """
        Rscript ./scripts/age_state_analysis.R --scenario {wildcards.scenario} --ci {CI_FLAG} --exclude_duplicates {EXCLUDE_DUPLICATES}
        """

rule age_lineplot:
    input:
        "results/{scenario}/df_RR_by_age_state.tsv"
    output:
        "figs/{scenario}/age_lineplot_all_censor.jpg",
        "figs/{scenario}/age_lineplot_by_state.jpg",
        "figs/{scenario}/age_lineplot_abridged.jpg",
        "figs/{scenario}/age_same_lineplot_by_state.jpg"
    shell:
        """
        Rscript scripts/age_lineplot.R --scenario {wildcards.scenario}
        """

# ============================================================================
# Census Division Analyses
# ============================================================================

rule census_div_analysis:
    input:
        "results/{scenario}/data_cleaned.txt"
    output:
        "results/{scenario}/df_RR_by_census_div.tsv"
    group: "duckdb_acc"
    resources:
        duckdb_lock=1
    shell:
        """
        Rscript ./scripts/census_div_analysis.R --scenario {wildcards.scenario} --ci {CI_FLAG} --exclude_duplicates {EXCLUDE_DUPLICATES}
        """

rule census_div_heatmap:
    input:
        "results/{scenario}/df_RR_by_census_div.tsv"
    output:
        "figs/{scenario}/census_divisions_heatmap.jpg"
    shell:
        """
        Rscript ./scripts/census_div_heatmap.R --scenario {wildcards.scenario}
        """

# ============================================================================
# Time Series Analyses
# ============================================================================

rule state_time_rr_analysis:
    input:
        "results/{scenario}/pairs_time_created.txt"
    output:
        "results/{scenario}/time_state/df_state_rr_all.tsv",
        "results/{scenario}/time_state/df_state_rr_snap.tsv",
        "results/{scenario}/time_state/df_state_rr_series.tsv"
    group: "duckdb_acc"
    resources:
        duckdb_lock=1
    shell:
        """
        Rscript ./scripts/state_time_rr_analysis.R --scenario {wildcards.scenario} --exclude_duplicates {EXCLUDE_DUPLICATES}
        """

rule state_time_significant_connections:
    input:
        "results/{scenario}/time_state/df_state_rr_snap.tsv"
    output:
        "results/{scenario}/time_state/df_significant_connections_3sd.tsv",
        "results/{scenario}/time_state/summary_significant_connections_by_time_3sd.tsv"
    shell:
        """
        Rscript ./scripts/state_time_significant_connections.R --scenario {wildcards.scenario} --sd_threshold 3.0
        """

rule state_time_network_viz:
    input:
        "results/{scenario}/time_state/df_significant_connections_3sd.tsv"
    output:
        "results/{scenario}/time_state/network_plots_complete.txt"
    shell:
        """
        Rscript ./scripts/state_time_significant_network_viz.R --scenario {wildcards.scenario} --sd_threshold 3.0
        touch {output}
        """

rule state_time_significant_connections_conserved:
    input:
        "results/{scenario}/time_state/df_state_rr_snap.tsv"
    output:
        f"results/{{scenario}}/time_state/df_significant_connections_{SD_THRESHOLD}sd.tsv",
        f"results/{{scenario}}/time_state/summary_significant_connections_by_time_{SD_THRESHOLD}sd.tsv"
    shell:
        f"""
        Rscript ./scripts/state_time_significant_connections.R --scenario {{wildcards.scenario}} --sd_threshold {SD_THRESHOLD} --min_nb_dist 1
        """

rule state_time_conserved_network_viz:
    input:
        f"results/{{scenario}}/time_state/df_significant_connections_{SD_THRESHOLD}sd.tsv"
    output:
        f"results/{{scenario}}/time_state/df_conserved_connections_{MIN_OCCURRENCES}plus.tsv",
        f"figs/{{scenario}}/time_state_networks/network_conserved_{MIN_OCCURRENCES}plus.jpg",
        f"figs/{{scenario}}/time_state_networks/network_conserved_{MIN_OCCURRENCES}plus_map.jpg"
    shell:
        f"""
        Rscript ./scripts/state_time_conserved_network_viz.R --scenario {{wildcards.scenario}} --sd_threshold {SD_THRESHOLD} --min_occurrences {MIN_OCCURRENCES}
        """

rule age_time_rr_analysis:
    input:
        "results/{scenario}/pairs_time_created.txt"
    output:
        "results/{scenario}/time_age/df_RR_by_time_age_series.tsv",
        "results/{scenario}/time_age/df_RR_by_school_state_time.tsv",
        "results/{scenario}/time_age/df_RR_by_school_state_ay.tsv"
    group: "duckdb_acc"
    resources:
        duckdb_lock=1
    shell:
        """
        Rscript ./scripts/age_time_RR_analysis.R --scenario {wildcards.scenario} --ci FALSE --aggregate TRUE --exclude_duplicates {EXCLUDE_DUPLICATES}
        """

rule school_analyses:
    input:
        "results/{scenario}/time_age/df_RR_by_school_state_time.tsv",
        "results/{scenario}/time_age/df_RR_by_school_state_ay.tsv"
    output:
        # age_school_analysis outputs
        "figs/{scenario}/age_time/age_school_ay_boxplot.jpg",
        "figs/{scenario}/age_time/age_school_ay_nRR_boxplot.jpg",
        "figs/{scenario}/age_time/age_school_ay_nRR_fixed_boxplot.jpg",
        "results/{scenario}/time_age/anova_school_ay_within_group_results.rds",
        "results/{scenario}/time_age/anova_school_ay_nRR_fixed_results.rds",
        # school_share_analysis outputs
        "results/{scenario}/time_age/state_trajectory_classifications.tsv",
        "figs/{scenario}/age_time/school_trajectories_by_type.jpg",
        "figs/{scenario}/age_time/school_share_timeseries.jpg",
        "figs/{scenario}/age_time/school_share_RR_correlation.jpg",
        "figs/{scenario}/age_time/school_share_nRR_fixed_correlation.jpg",
        "results/{scenario}/time_age/df_RR_school_share_combined.tsv"
    shell:
        """
        Rscript ./scripts/age_school_analysis.R --scenario {wildcards.scenario}
        Rscript ./scripts/school_share_analysis.R --scenario {wildcards.scenario}
        """

rule age_time_plots:
    input:
        "results/{scenario}/time_age/df_RR_by_time_age_series.tsv"
    output:
        "figs/{scenario}/age_time/USA/age_RR_time_faceted.png",
        "figs/{scenario}/age_time/USA/age_nRR_time_faceted.png",
        "figs/{scenario}/age_time/USA/age_nRR_fixed_time_faceted_all.png",
        "figs/{scenario}/age_time/all_countries/vaccine_period_RR_elderly.png",
        "figs/{scenario}/age_time/all_countries/vaccine_period_nRR_fixed_elderly.png"
    shell:
        """
        Rscript ./scripts/age_time_plot.R --scenario {wildcards.scenario}
        """

rule state_time_visualizations:
    input:
        "results/{scenario}/time_state/df_state_rr_snap.tsv",
        "results/{scenario}/time_state/df_state_rr_series.tsv"
    output:
        # From state_time_scatter.R
        "figs/{scenario}/time/state_pair_nRR_heatmap.png",
        "figs/{scenario}/time/all_pairs_fold_heatmap.png",
        "figs/{scenario}/time/all_pairs_fold_boxplot_quarters.png",
        "figs/{scenario}/time/rr_series.png",
        # From state_time_dist.R
        "figs/{scenario}/time/rr_boxplot_time.png",
        "figs/{scenario}/time/nb_boxplot_time.png",
        "figs/{scenario}/time/min_cbsa_dist_time.png",
        "figs/{scenario}/time/nhts_time.png",
        "figs/{scenario}/time/safegraph_time.png",
        "figs/{scenario}/time/nhts_time_normalized.png",
        "figs/{scenario}/time/safegraph_time_normalized.png",
        "figs/{scenario}/time/correlations_time_series.png",
        "figs/{scenario}/time/correlations_time_series_split_air.png",
        "figs/{scenario}/time/correlations_time_series_normalized.png"
    shell:
        """
        Rscript ./scripts/state_time_scatter.R --scenario {wildcards.scenario}
        Rscript ./scripts/state_time_dist.R --scenario {wildcards.scenario}
        """

rule state_time_clustering:
    input:
        "results/{scenario}/time_state/df_state_rr_snap.tsv"
    output:
        "figs/{scenario}/clust/htree.jpg",
        "figs/{scenario}/clust/pcoa_variance_explained.jpg",
        "figs/{scenario}/clust/pcoa_V1V2.jpg",
        "figs/{scenario}/clust/pcoa_V1V3.jpg",
        "figs/{scenario}/clust/pcoa_V1V4.jpg",
        "figs/{scenario}/clust/pcoa_V2V3.jpg",
        "figs/{scenario}/clust/pcoa_V2V3_US.jpg"
    shell:
        """
        Rscript ./scripts/state_time_clust.R --scenario {wildcards.scenario}
        """

rule state_regression:
    input:
        "results/{scenario}/df_RR_by_state.tsv"
    output:
        "figs/{scenario}/travel_seq_air_length.jpg",
        "figs/{scenario}/dist/regression_fit.jpg"
    shell:
        """
        Rscript ./scripts/state_regression.R --scenario {wildcards.scenario}
        """
