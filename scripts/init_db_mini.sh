#! /bin/bash

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

OG_DB_FILE="db_files/db_${SCENARIO}.duckdb"
DB_FILE="db_files/db_${SCENARIO}_mini.duckdb"
META_FILE="data/${SCENARIO}/metadata/metadata_${SCENARIO}.tsv"
PAIR_FILE="data/${SCENARIO}/distance_aggregated/combined_df_identical_pairs_${SCENARIO}.tsv"
REMOVE_SCRIPT="scripts/remove_duplicates.sql"

# echo "Starting metadata table creation..."
# START_TIME=$(date +%s)
# duckdb ${DB_FILE} -c "CREATE OR REPLACE TABLE metadata AS SELECT strain,date,region,country,division,location,host,age,sex,Nextstrain_clade,pango_lineage,clade_nextstrain,clade_who,coverage,originating_lab,paper_url,sampling_strategy FROM read_csv('${META_FILE}',header = True,sep = '\t',types={'date':'VARCHAR'});"
# END_TIME=$(date +%s)
# echo "Metadata table creation completed in $(($END_TIME - $START_TIME)) seconds."

# echo "Starting pairs table creation..."
# START_TIME=$(date +%s)
# duckdb ${DB_FILE} -c "CREATE OR REPLACE TABLE pairs AS SELECT * FROM read_csv('${PAIR_FILE}',header = True,sep = '\t');"
# END_TIME=$(date +%s)
# echo "Pairs table creation completed in $(($END_TIME - $START_TIME)) seconds."

# echo "Duplicate pair checking and deletion..."
# START_TIME=$(date +%s)
# duckdb ${DB_FILE} < scripts/remove_duplicates.sql
# END_TIME=$(date +%s)
# echo "Cleaning completed in $(($END_TIME - $START_TIME)) seconds."

echo "Copy db file to make a mini db file"
START_TIME=$(date +%s)
cp ${OG_DB_FILE} ${DB_FILE}
END_TIME=$(date +%s)
echo "Copying completed in $(($END_TIME - $START_TIME)) seconds."

echo "Sampling for the mini set"
START_TIME=$(date +%s)
duckdb ${DB_FILE} < scripts/mini_sample.sql
END_TIME=$(date +%s)
echo "Sampling completed in $(($END_TIME - $START_TIME)) seconds."

echo "All tasks completed successfully!"