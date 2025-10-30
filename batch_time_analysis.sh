#!/bin/bash
#SBATCH --job-name=time_analysis
#SBATCH --output=time_analysis.o
#SBATCH --error=time_analysis.e
#SBATCH --time=14-0:0:0
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abemania@fredhutch.org

ml fhR
ml Pandoc/2.13 #For plotly

SCENARIO="CAM_1000"

# Run time preparation script
echo "Running time_prep.R..."
Rscript scripts/time_prep.R

# Run state-time RR analysis
echo "Running state_time_rr_analysis.R..."
Rscript scripts/state_time_rr_analysis.R --scenario ${SCENARIO}

# Run state-time plotting scripts
echo "Running state_time_clust.R..."
Rscript scripts/state_time_clust.R --scenario ${SCENARIO}

echo "Running state_time_dist.R..."
Rscript scripts/state_time_dist.R --scenario ${SCENARIO}

echo "Running state_time_scatter.R..."
Rscript scripts/state_time_scatter.R --scenario ${SCENARIO}

# Run age-time RR analysis
echo "Running age_time_RR_analysis.R..."
Rscript scripts/age_time_RR_analysis.R --scenario ${SCENARIO} --ci FALSE --aggregate TRUE

# Run age-time plotting script
echo "Running age_time_plot.R..."
Rscript scripts/age_time_plot.R --scenario ${SCENARIO}

# Run age-school analysis
echo "Running age_school_analysis.R..."
Rscript scripts/age_school_analysis.R --scenario ${SCENARIO}

#Run school share analyes
echo "Running school_share_analysis.R..."
Rscript scripts/school_share_analysis.R --scenario ${SCENARIO}


echo "All time-related analyses complete!"
