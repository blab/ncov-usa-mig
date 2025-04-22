#!/bin/bash
#SBATCH --job-name=batch
#SBATCH --output=batch.o
#SBATCH --error=batch.e
#SBATCH --time=14-0:0:0
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abemania@fredhutch.org

ml snakemake/7.18.2-foss-2021b
snakemake --profile ./profile --cores 32 --group-components duckdb_acc=1