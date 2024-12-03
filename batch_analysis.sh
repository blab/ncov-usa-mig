#! /bin/bash
#SBATCH --job-name=batch_rr
#SBATCH --output=batch_rr.o
#SBATCH --error=batch_rr.e
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=7-0:0:0
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abemania@fredhutch.org

SCENARIO="USA"
THREADS=32
CI_FLAG=TRUE

ml fhR

echo "Starting time: $(date +"%T")"

# ./scripts/init_db.sh -s $SCENARIO
# echo "Database built: $(date +"%T")"

# Rscript ./scripts/clean_data.R --scenario $SCENARIO
# echo "Cleaning done: $(date +"%T")"

# Rscript ./scripts/age_analysis.R --scenario $SCENARIO --ci $CI_FLAG
# Rscript ./scripts/age_heatmap.R --scenario $SCENARIO
# echo "Age analysis done: $(date +"%T")"

# Rscript ./scripts/state_analysis.R --scenario $SCENARIO --ci $CI_FLAG
# Rscript ./scripts/state_heatmap.R --scenario $SCENARIO
# Rscript ./scripts/state_dist_plots.R --scenario $SCENARIO
# echo "State analysis done: $(date +"%T")"

Rscript ./scripts/age_state_analysis.R --scenario $SCENARIO --ci $CI_FLAG
echo "State/Age stratified analysis done: $(date +"%T")"

# Rscript ./scripts/census_div_analysis.R --scenario $SCENARIO --ci $CI_FLAG
# Rscript ./scripts/census_div_heatmap.R --scenario $SCENARIO
# echo "Census Division analysis done: $(date +"%T")"

echo "ANALYSIS COMPLETE!"