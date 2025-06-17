configfile: "config.yaml"
CI_FLAG = config['ci_flag']

SCENARIOS = config["scenarios"]  # <-- change this to support more, like ["USA", "CAN", "MEX"]

# Collect final outputs for all scenarios
all_outputs = expand([
#    "figs/{scenario}/age_heatmap.jpg",
#    "figs/{scenario}/age_same_lineplot_by_state.jpg",
#    "results/{scenario}/summary_tables/tables_{scenario}.tar.zst",
    "figs/{scenario}/state_heatmap.jpg"
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

rule age_analysis_RR:
    input:
        "results/{scenario}/data_cleaned.txt"
    output:
        "results/{scenario}/df_RR_by_age_class.tsv"
    group: "duckdb_acc"
    shell:
        """
        Rscript ./scripts/age_analysis.R --scenario {wildcards.scenario} --ci {config[ci_flag]}
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

rule state_analysis_RR:
    input:
        "results/{scenario}/data_cleaned.txt"
    output:
        "results/{scenario}/df_RR_by_state.tsv"
    group: "duckdb_acc"
    shell:
        """
        Rscript ./scripts/state_analysis.R --scenario {wildcards.scenario} --ci {config[ci_flag]}
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

rule age_state_RR:
    input:
        age="results/{scenario}/df_RR_by_age_class.tsv",
        state="results/{scenario}/df_RR_by_state.tsv"
    output:
        "results/{scenario}/df_RR_by_age_state.tsv"
    group: "duckdb_acc"
    shell:
        """
        Rscript ./scripts/age_state_analysis.R --scenario {wildcards.scenario} --ci {config[ci_flag]}
        """

rule age_lineplot:
    input:
        "results/{scenario}/df_RR_by_age_state.tsv"
    output:
        "results/{scenario}/age_same_lineplot_by_state.jpg"
    shell:
        """
        Rscript scripts/age_lineplot.R --scenario {wildcards.scenario}
        """

# rule census_div_analysis:
#     input:
#         "results/{scenario}/age_state_analysis_done.txt"
#     output:
#         touch("results/{scenario}/census_div_analysis_done.txt")
#     group: "duckdb_acc" 
#     shell:
#         """
#         Rscript ./scripts/census_div_analysis.R --scenario {config[scenario]} --ci {config[ci_flag]}
#         Rscript ./scripts/census_div_heatmap.R --scenario {config[scenario]}
#         """