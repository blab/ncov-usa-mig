configfile: "config.yaml"
CI_FLAG = config['ci_flag']
EXCLUDE_DUPLICATES = config['exclude_duplicates']

SCENARIOS = config["scenarios"]

# Parameters for conserved network analysis
PERCENTILE_THRESHOLD = 0.98
PERCENTILE_LABEL = round(PERCENTILE_THRESHOLD * 100)
MIN_OCCURRENCES = 5
# Extra percentile thresholds used by the conserved-vs-non-conserved boxplots in state_dist_plots.R
AUX_PERCENTILES = ["70", "80", "90", "95"]

# Scenario used for manuscript figures (stitch scripts and scripts that write to manuscript/figures/)
MS = config.get("manuscript_scenario", "CAM_1000")

# Collect final outputs for all scenarios
all_outputs = expand([
    # Summary tables
    "results/{scenario}/summary_tables/demo_tables.tar.zst",
    # Descriptive figures
    "figs/{scenario}/figure_1.jpg",
    # Age analyses
    "results/{scenario}/df_RR_by_age_class.tsv",
    "figs/{scenario}/age_heatmaps/full.jpg",
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
    "figs/{scenario}/census_divisions_heatmap.pdf",
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
    # Significant connections network analysis
    "results/{scenario}/time_state/df_significant_connections_3sd.tsv",
    "results/{scenario}/time_state/network_plots_complete.txt",
    # Conserved network map
    f"results/{{scenario}}/time_state/df_conserved_connections_{PERCENTILE_LABEL}pct_{MIN_OCCURRENCES}plus.tsv"
], scenario=SCENARIOS)

shell.prefix("set -euo pipefail; ")

# Manuscript figures (main + supplementary). Simulation figures (movement_network,
# rt_combined, recovery_by_phase_combined) come from a separate repo and are not built here.
MANUSCRIPT_FIGS = [
    "manuscript/figures/desc.pdf",
    "manuscript/figures/geo_clust.pdf",
    "manuscript/figures/geo_pc.pdf",
    "manuscript/figures/geo_time.pdf",
    "manuscript/figures/geo_dist.pdf",
    "manuscript/figures/age.pdf",
    "manuscript/figures/age_time.pdf",
    "manuscript/figures/supp/ba1_epiweek_seqdist.pdf",
    "manuscript/figures/supp/regions_heatmap.pdf",
    "manuscript/figures/supp/hclust_subset.pdf",
    "manuscript/figures/supp/predictor_cor_heatmap.pdf",
    "manuscript/figures/supp/lm_gam_obs_exp.pdf",
    "manuscript/figures/supp/age_contact_correlation.pdf",
]

rule all:
    input: all_outputs, MANUSCRIPT_FIGS

rule manuscript_figures:
    input: MANUSCRIPT_FIGS

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
        "figs/{scenario}/figure_1.jpg",
        "figs/{scenario}/desc/seq_map.png",
        "figs/{scenario}/desc/natl_effort.svg",
        "figs/{scenario}/desc/age_hist.svg",
        "figs/{scenario}/desc/age_region.svg"
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
        "results/{scenario}/df_RR_by_age_class.tsv",
        "results/{scenario}/df_RR_by_age_state.tsv"
    output:
        "figs/{scenario}/age_heatmaps/full.jpg",
        "figs/{scenario}/age_heatmaps/full.svg",
        "figs/{scenario}/age_heatmaps/young_0_24.svg",
        "figs/{scenario}/age_heatmaps/subsets_row.svg"
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
        "results/{scenario}/df_RR_by_state.tsv",
        "results/{scenario}/time_state/df_significant_connections_98_percentile.tsv",
        expand("results/{{scenario}}/time_state/df_conserved_connections_{pct}pct_5plus.tsv", pct=AUX_PERCENTILES)
    output:
        "figs/{scenario}/state_nb_dist_plot.jpg",
        "figs/{scenario}/dist/state_distance_poster.svg",
        "figs/{scenario}/dist/state_nhts_conserved_boxplot.svg",
        "figs/{scenario}/dist/state_air_db1b_conserved_boxplot.svg",
        "figs/{scenario}/dist/state_safegraph_conserved_boxplot.svg"
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
        "figs/{scenario}/census_divisions_heatmap.pdf"
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

rule state_time_connections_percentile:
    input:
        "results/{scenario}/time_state/df_state_rr_snap.tsv"
    output:
        "results/{scenario}/time_state/df_significant_connections_{pct}_percentile.tsv"
    wildcard_constraints:
        pct="[0-9]+"
    params:
        frac=lambda wc: int(wc.pct) / 100
    shell:
        """
        Rscript ./scripts/state_time_connections_threshold.R --scenario {wildcards.scenario} --percentile_threshold {params.frac}
        """

# Main conserved network (98th percentile) - also generates the network figures used in Fig. geo_time
rule gen_conserved_network:
    input:
        f"results/{{scenario}}/time_state/df_significant_connections_{PERCENTILE_LABEL}_percentile.tsv"
    output:
        f"results/{{scenario}}/time_state/df_conserved_connections_{PERCENTILE_LABEL}pct_{MIN_OCCURRENCES}plus.tsv",
        "figs/{scenario}/time_state_networks/network_conserved_combined.jpg",
        "figs/{scenario}/time_state_networks/network_conserved_combined.svg"
    shell:
        f"""
        Rscript ./scripts/gen_conserved_network.R --scenario {{wildcards.scenario}} --percentile_threshold {PERCENTILE_THRESHOLD} --min_occurrences {MIN_OCCURRENCES} --gen_viz TRUE
        """

# Auxiliary thresholds (data only; viz off so the main network figures are not overwritten)
rule gen_conserved_network_aux:
    input:
        "results/{scenario}/time_state/df_significant_connections_{pct}_percentile.tsv"
    output:
        f"results/{{scenario}}/time_state/df_conserved_connections_{{pct}}pct_{MIN_OCCURRENCES}plus.tsv"
    wildcard_constraints:
        pct="|".join(AUX_PERCENTILES)
    params:
        frac=lambda wc: int(wc.pct) / 100
    shell:
        f"""
        Rscript ./scripts/gen_conserved_network.R --scenario {{wildcards.scenario}} --percentile_threshold {{params.frac}} --min_occurrences {MIN_OCCURRENCES} --gen_viz FALSE
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
        "figs/{scenario}/time/rr_series_poster.svg",
        "figs/{scenario}/time/all_pairs_fold_heatmap.svg",
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
        "figs/{scenario}/clust/pcoa_V1V2.svg",
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
        f"results/{MS}/df_RR_by_state.tsv"
    output:
        f"figs/{MS}/travel_seq_air_length.jpg",
        f"figs/{MS}/dist/regression_fit.jpg",
        "manuscript/figures/supp/predictor_cor_heatmap.pdf"
    shell:
        f"""
        Rscript ./scripts/state_regression.R --scenario {MS}
        """

# ============================================================================
# Manuscript figure panels (scripts previously run by hand)
# ============================================================================

rule state_RR_maps:
    input:
        "results/{scenario}/df_RR_by_state.tsv"
    output:
        "figs/{scenario}/RR_maps/patched_maps.jpg"
    shell:
        """
        Rscript ./scripts/state_RR_maps.R --scenario {wildcards.scenario}
        """

rule state_heatmap_clustered:
    input:
        "results/{scenario}/df_RR_by_state.tsv",
        "results/{scenario}/time_state/df_state_rr_all.tsv"
    output:
        "figs/{scenario}/state_heatmap_clustered.jpg",
        "figs/{scenario}/state_heatmap_clustered.svg",
        "figs/{scenario}/state_dendro_subset.pdf"
    shell:
        """
        Rscript ./scripts/state_heatmap_clustered.R --scenario {wildcards.scenario}
        """

rule age_RR_deviance:
    input:
        "results/{scenario}/df_RR_by_age_class.tsv",
        "results/{scenario}/df_RR_by_age_state.tsv"
    output:
        "figs/{scenario}/age_RR_deviance_geographic_compact.svg"
    shell:
        """
        Rscript ./scripts/age_RR_deviance.R --scenario {wildcards.scenario}
        """

rule age_vaccine_period:
    input:
        "results/{scenario}/pairs_time_created.txt",
        "results/{scenario}/time_age/df_RR_by_time_age_series.tsv"
    output:
        "results/{scenario}/time_age/df_vaccine_period_nRR_ci.tsv",
        "figs/{scenario}/age_time/all_countries/vaccine_period_nRR_fixed_elderly_ci.svg"
    group: "duckdb_acc"
    resources:
        duckdb_lock=1
    shell:
        """
        Rscript ./scripts/age_vaccine_period.R --scenario {wildcards.scenario} --exclude_duplicates {EXCLUDE_DUPLICATES}
        """

rule school_nrr_timeseries:
    input:
        "results/{scenario}/pairs_time_created.txt"
    output:
        "results/{scenario}/time_age/df_nRR_school_season_timeseries.tsv",
        "figs/{scenario}/age_time/school/school_nRR_trace_simplified.svg"
    group: "duckdb_acc"
    resources:
        duckdb_lock=1
    shell:
        """
        Rscript ./scripts/school_nrr_timeseries.R --scenario {wildcards.scenario} --exclude_duplicates {EXCLUDE_DUPLICATES}
        """

# The following scripts write their supplementary PDF straight into manuscript/figures/supp,
# so they run on the manuscript scenario only.

rule state_tensor_plots:
    input:
        f"results/{MS}/df_RR_by_state.tsv"
    output:
        f"figs/{MS}/dist/gam_dev_by_predictor.svg",
        f"figs/{MS}/dist/pred_rr_curves_combined.svg",
        f"figs/{MS}/dist/lm_gam_obs_exp.svg",
        "manuscript/figures/supp/lm_gam_obs_exp.pdf"
    shell:
        f"""
        Rscript ./scripts/state_tensor_plots.R --scenario {MS}
        """

rule ba1_epiweek_supp:
    input:
        f"results/{MS}/data_cleaned.txt"
    output:
        f"figs/{MS}/supp/ba1_epiweek_seqdist.jpg",
        "manuscript/figures/supp/ba1_epiweek_seqdist.pdf"
    group: "duckdb_acc"
    resources:
        duckdb_lock=1
    shell:
        f"""
        Rscript ./scripts/ba1_epiweek_supp.R --scenario {MS}
        """

rule age_contact_correlation:
    input:
        f"results/{MS}/data_cleaned.txt",
        "data/us_contact_matrix.csv",
        "data/canada_contact_matrix.csv"
    output:
        f"results/{MS}/df_RR_by_age_USA.tsv",
        f"results/{MS}/df_RR_by_age_Canada.tsv",
        f"figs/{MS}/age_contact_correlation.svg",
        "manuscript/figures/supp/age_contact_correlation.pdf"
    group: "duckdb_acc"
    resources:
        duckdb_lock=1
    shell:
        f"""
        Rscript ./scripts/age_contact_correlation.R --scenario {MS} --exclude_duplicates {EXCLUDE_DUPLICATES}
        """

rule supp_regions_heatmap:
    input:
        f"figs/{MS}/census_divisions_heatmap.pdf"
    output:
        "manuscript/figures/supp/regions_heatmap.pdf"
    shell:
        "cp {input} {output}"

rule supp_hclust_subset:
    input:
        f"figs/{MS}/state_dendro_subset.pdf"
    output:
        "manuscript/figures/supp/hclust_subset.pdf"
    shell:
        "cp {input} {output}"

# ============================================================================
# Stitched manuscript figures (requires inkscape for PDF export)
# ============================================================================

STITCH_INPUTS = {
    "desc": [f"figs/{MS}/desc/seq_map.png",
             f"figs/{MS}/desc/natl_effort.svg",
             f"figs/{MS}/desc/age_hist.svg",
             f"figs/{MS}/desc/age_region.svg"],
    "geo_clust": [f"figs/{MS}/RR_maps/patched_maps.jpg",
                  f"figs/{MS}/state_heatmap_clustered.svg"],
    "geo_pc": [f"figs/{MS}/bea_region_map.png",
               f"figs/{MS}/clust/pcoa_V1V2.svg"],
    "geo_time": [f"figs/{MS}/time/rr_series_poster.svg",
                 f"figs/{MS}/time/all_pairs_fold_heatmap.svg",
                 f"figs/{MS}/time_state_networks/network_conserved_combined.jpg"],
    "geo_dist": [f"figs/{MS}/dist/state_distance_poster.svg",
                 f"figs/{MS}/dist/gam_dev_by_predictor.svg",
                 f"figs/{MS}/dist/pred_rr_curves_combined.svg"],
    "age": [f"figs/{MS}/age_heatmaps/full.svg",
            f"figs/{MS}/age_heatmaps/young_0_24.svg",
            f"figs/{MS}/age_heatmaps/subsets_row.svg",
            f"figs/{MS}/age_RR_deviance_geographic_compact.svg"],
    "age_time": [f"figs/{MS}/age_time/all_countries/vaccine_period_nRR_fixed_elderly_ci.svg",
                 f"figs/{MS}/age_time/school/school_nRR_trace_simplified.svg"],
}

rule stitch_figure:
    input:
        script="scripts/stitch/stitch_{fig}.py",
        panels=lambda wc: STITCH_INPUTS[wc.fig]
    output:
        svg="manuscript/figures/{fig}.svg",
        pdf="manuscript/figures/{fig}.pdf"
    wildcard_constraints:
        fig="|".join(STITCH_INPUTS.keys())
    shell:
        f"""
        python3 {{input.script}} --scenario {MS} --out {{output.svg}} --pdf --pdf-out {{output.pdf}}
        """
