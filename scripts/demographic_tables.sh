#!/bin/bash

# Script: demographic_table.sh
# Description: Generates summary tables (count + percentage + total row) for specified columns
#              from a DuckDB database and exports each as a TSV file.
#              Then compresses all TSVs into a single .tar.zst archive.

# Function to show usage
ml zstd

usage() {
  echo "Usage: $0 -s <scenario> -c \"<column1> <column2> ...\""
  echo ""
  echo "Options:"
  echo "  -s    Scenario name (used to locate DB file and name output folder)"
  echo "  -c    Quoted list of column names to summarize (space-separated)"
  echo "  -h    Show this help message"
  echo ""
  echo "Example:"
  echo "  $0 -s USA_mini -c \"clade_nextstrain pango_lineage division\""
  exit 1
}

# Parse arguments
while getopts ":s:c:h" opt; do
  case ${opt} in
    s ) SCENARIO="$OPTARG" ;;
    c ) COLUMN_LIST="$OPTARG" ;;
    h ) usage ;;
    \? )
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
    : )
      echo "Option -$OPTARG requires an argument." >&2
      usage
      ;;
  esac
done

# Validate required args
if [ -z "$SCENARIO" ] || [ -z "$COLUMN_LIST" ]; then
  echo "Error: Both -s (scenario) and -c (columns) are required." >&2
  usage
fi

# Define paths
TABLE_NAME="metadata"
DB_FILE="db_files/db_${SCENARIO}.duckdb"
OUTPUT_DIR="results/${SCENARIO}/summary_tables"
ARCHIVE_FILE="results/${SCENARIO}/summary_tables/demo_tables.tar.zst"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Generate summary tables
for COLUMN_NAME in $COLUMN_LIST; do
  OUTPUT_FILE="${OUTPUT_DIR}/${COLUMN_NAME}_summary.tsv"
  echo "[demographic_table.sh] Generating summary for '${COLUMN_NAME}' → ${OUTPUT_FILE}"

  duckdb "$DB_FILE" <<EOF
COPY (
  SELECT 
    ${COLUMN_NAME} AS category,
    COUNT(*) AS count,
    ROUND(100.0 * COUNT(*) / (SELECT COUNT(*) FROM ${TABLE_NAME}), 1) AS percentage
  FROM ${TABLE_NAME}
  GROUP BY ${COLUMN_NAME}

  UNION ALL

  SELECT 
    'Total' AS category,
    COUNT(*) AS count,
    100.0 AS percentage
  FROM ${TABLE_NAME}
) TO '${OUTPUT_FILE}' (HEADER, DELIMITER '\t');
EOF

done

# Compress all TSVs into a zst archive
echo "[demographic_table.sh] Compressing TSVs into ${ARCHIVE_FILE}"
tar --use-compress-program=zstd -cf "$ARCHIVE_FILE" -C "$OUTPUT_DIR" .


if [ -f "$ARCHIVE_FILE" ]; then
  echo "[demographic_table.sh] Archive verified. Deleting TSVs..."
  find "$OUTPUT_DIR" -type f -name "*.tsv" -delete
  echo "[demographic_table.sh] Done. Archive created at: $ARCHIVE_FILE"
else
  echo "[demographic_table.sh] Warning: Archive not created, skipping cleanup."
fi


