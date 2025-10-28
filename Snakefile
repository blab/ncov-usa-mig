configfile: "config.yaml"
CI_FLAG = config['ci_flag']

SCENARIOS = config["scenarios"]

# Parameters for conserved network analysis
SD_THRESHOLD = 2
MIN_OCCURRENCES = 4

# Collect final outputs for all scenarios
all_outputs = expand([
    # Summary tables
    "results/{scenario}/summary_tables/tables_{scenario}.tar.zst",
    # Age analyses
    "results/{scenario}/df_RR_by_age_class.tsv",
    "figs/{scenario}/age_heatmap.jpg",
    "figs/{scenario}/age_lineplot_all_censor.jpg",
    "figs/{scenario}/age_lineplot_by_state.jpg",
    "figs/{scenario}/age_lineplot_abridged.jpg",
    "figs/{scenario}/age_same_lineplot_by_state.jpg",
    # State analyses
    "results/{scenario}/df_RR_by_state.tsv",
    "figs/{scenario}/state_heatmap.jpg",
    "figs/{scenario}/state_nb_dist_plot.jpg",
    # Age x State analyses
    "results/{scenario}/df_RR_by_age_state.tsv",
    # Census division analyses
    "results/{scenario}/df_RR_by_census_div.tsv",
    "figs/{scenario}/census_divisions_heatmap.jpg",
    # Time series analyses
    "results/{scenario}/time_state/df_state_rr_series.tsv",
    "results/{scenario}/age_time/df_RR_by_age_time_series.tsv",
    # Significant connections network analysis
    "results/{scenario}/time_state/df_significant_connections_3.0sd.tsv",
    "results/{scenario}/time_state/network_plots_complete.txt",
    # Conserved network map
    f"results/{{scenario}}/time_state/df_significant_connections_{SD_THRESHOLD}sd.tsv",
    f"results/{{scenario}}/time_state/df_conserved_connections_{MIN_OCCURRENCES}plus.tsv",
    f"figs/{{scenario}}/time_state_networks/network_conserved_{MIN_OCCURRENCES}plus_map.jpg"
], scenario=SCENARIOS)

shell.prefix("ml fhR; set -euo pipefail; ")

rule all:
    input: all_outputs

rule database:
    output:
        "db_files/db_{scenario}.duckdb"
    group: "duckdb_acc"
    shell:
        """
        ml zstd
        unzstd -f data/{wildcards.scenario}/metadata/metadata_{wildcards.scenario}.tsv.zst
        unzstd -f data/{wildcards.scenario}/distance_aggregated/combined_df_identical_pairs_{wildcards.scenario}.tsv.zst
        ./scripts/init_db.sh -s {wildcards.scenario}
        """

rule clean_data:
    input:
        "db_files/db_{scenario}.duckdb"
    output:
        touch("results/{scenario}/data_cleaned.txt")
    group: "duckdb_acc"
    shell:
        """
        Rscript ./scripts/clean_data.R --scenario {wildcards.scenario}
        duckdb db_files/db_{wildcards.scenario}.duckdb < scripts/trim_pairs.sql
        """

rule demo_tables:
    input:
        "results/{scenario}/data_cleaned.txt"
    output:
        "results/{scenario}/summary_tables/tables_{scenario}.tar.zst"
    group: "duckdb_acc"
    shell:
        """
        scripts/demographic_tables.sh -s {wildcards.scenario} -c "clade_nextstrain division census_div sex age_class"
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
    shell:
        """
        Rscript ./scripts/age_analysis.R --scenario {wildcards.scenario} --ci {CI_FLAG}
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

# ============================================================================
# State Analyses
# ============================================================================

rule state_analysis_RR:
    input:
        "results/{scenario}/data_cleaned.txt"
    output:
        "results/{scenario}/df_RR_by_state.tsv"
    group: "duckdb_acc"
    shell:
        """
        Rscript ./scripts/state_analysis.R --scenario {wildcards.scenario} --ci {CI_FLAG}
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
    shell:
        """
        Rscript ./scripts/age_state_analysis.R --scenario {wildcards.scenario} --ci {CI_FLAG}
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
    shell:
        """
        Rscript ./scripts/census_div_analysis.R --scenario {wildcards.scenario} --ci {CI_FLAG}
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
        "results/{scenario}/data_cleaned.txt"
    output:
        "results/{scenario}/time_state/df_state_rr_all.tsv",
        "results/{scenario}/time_state/df_state_rr_snap.tsv",
        "results/{scenario}/time_state/df_state_rr_series.tsv"
    group: "duckdb_acc"
    shell:
        """
        Rscript ./scripts/state_time_rr_analysis.R --scenario {wildcards.scenario}
        """

rule state_time_significant_connections:
    input:
        "results/{scenario}/time_state/df_state_rr_series.tsv"
    output:
        "results/{scenario}/time_state/df_significant_connections_3.0sd.tsv",
        "results/{scenario}/time_state/summary_significant_connections_by_time_3.0sd.tsv"
    shell:
        """
        Rscript ./scripts/state_time_significant_connections.R --scenario {wildcards.scenario} --sd_threshold 3.0
        """

rule state_time_network_viz:
    input:
        "results/{scenario}/time_state/df_significant_connections_3.0sd.tsv"
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
        f"results/{{scenario}}/time_state/df_conserved_connections_{MIN_OCCURRENCES}plus.tsv"
    shell:
        f"""
        Rscript ./scripts/state_time_conserved_network_viz.R --scenario {{wildcards.scenario}} --sd_threshold {SD_THRESHOLD} --min_occurrences {MIN_OCCURRENCES}
        """

rule state_time_conserved_network_map:
    input:
        f"results/{{scenario}}/time_state/df_conserved_connections_{MIN_OCCURRENCES}plus.tsv"
    output:
        f"figs/{{scenario}}/time_state_networks/network_conserved_{MIN_OCCURRENCES}plus_map.jpg"
    shell:
        f"""
        Rscript ./scripts/state_time_conserved_network_map.R --scenario {{wildcards.scenario}} --sd_threshold {SD_THRESHOLD} --min_occurrences {MIN_OCCURRENCES}
        """

rule age_time_rr_analysis:
    input:
        "results/{scenario}/data_cleaned.txt"
    output:
        "results/{scenario}/age_time/df_RR_by_age_time_series.tsv",
        "results/{scenario}/age_time/df_RR_by_school_state_time.tsv",
        "results/{scenario}/age_time/df_RR_by_school_state_ay.tsv"
    group: "duckdb_acc"
    shell:
        """
        Rscript ./scripts/age_time_RR_analysis.R --scenario {wildcards.scenario} --ci FALSE --aggregate TRUE
        """
