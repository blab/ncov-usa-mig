#! /bin/bash
# Filename: init_db.sh
# Author(s): Amin Bemanian
# Date: 07/08/25
# Description: DuckDB initialization using idseq output

usage()
{
    echo "Usage: [ -s SCENARIO_NAME]"
    exit 2
}

mkdir -p db_files

while getopts 's:?H' flag; do
    case "${flag}" in
        s) SCENARIO=${OPTARG} ;;
        h|?) usage ;; 
    esac
done

DB_FILE="db_files/db_${SCENARIO}.duckdb"
META_FILE="data/${SCENARIO}/metadata/metadata_${SCENARIO}.tsv"
PAIR_FILE="data/${SCENARIO}/distance_aggregated/combined_df_identical_pairs_${SCENARIO}.tsv"
REMOVE_SCRIPT="scripts/remove_duplicates.sql"

echo "Starting metadata table creation..."
START_TIME=$(date +%s)
duckdb ${DB_FILE} -c "CREATE OR REPLACE TABLE metadata AS SELECT strain,date,region,country,division,location,host,age,sex,Nextstrain_clade,pango_lineage,clade_nextstrain,clade_who,coverage,originating_lab,paper_url,sampling_strategy FROM read_csv('${META_FILE}',header = True,sep = '\t',types={'date':'VARCHAR'});"
END_TIME=$(date +%s)
echo "Metadata table creation completed in $(($END_TIME - $START_TIME)) seconds."

echo "Starting pairs table creation..."
START_TIME=$(date +%s)
duckdb ${DB_FILE} -c "CREATE OR REPLACE TABLE pairs AS SELECT * FROM read_csv('${PAIR_FILE}',header = True,sep = '\t');"
END_TIME=$(date +%s)
echo "Pairs table creation completed in $(($END_TIME - $START_TIME)) seconds."

echo "Duplicate pair checking and deletion..."
START_TIME=$(date +%s)
duckdb ${DB_FILE} < scripts/remove_duplicates.sql
END_TIME=$(date +%s)
echo "Pairs table creation completed in $(($END_TIME - $START_TIME)) seconds."

echo "All tasks completed successfully!"
