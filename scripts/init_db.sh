#! /bin/bash

usage()
{
    echo "Usage: [ -s SCENARIO_NAME]"
    exit 2
}

while getopts 's:?H' flag; do
    case "${flag}" in
        s) SCENARIO=${OPTARG} ;;
        h|?) usage ;; 
    esac
done

DB_FILE="db_files/db_${SCENARIO}.duckdb"
META_FILE="data/${SCENARIO}/metadata/metadata_${SCENARIO}.tsv"
PAIR_FILE="./data/${SCENARIO}/distance_aggregated/combined_df_identical_pairs_${SCENARIO}.tsv"

duckdb ${DB_FILE} -c "CREATE OR REPLACE TABLE metadata AS SELECT strain,date,region,country,division,location,host,age,sex,Nextstrain_clade,pango_lineage,clade_nextstrain,clade_who,coverage FROM read_csv('${META_FILE}',header = True,sep = '\t',types={'date':'VARCHAR'});"
duckdb ${DB_FILE} -c "CREATE OR REPLACE TABLE pairs AS SELECT * FROM read_csv('${PAIR_FILE}',header = True,sep = '\t');"
