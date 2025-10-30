#!/bin/bash
#SBATCH --job-name=compress_cluster
#SBATCH --output=compress_cluster_%j.out
#SBATCH --error=compress_cluster_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=320G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abemania@fredhutch.org

ml zstd

echo "Starting compression of cluster-analysis directory..."
echo "Start time: $(date)"

tar -cf - cluster-analysis/ | zstd -o cluster-analysis.tar.zst

echo "Compression complete!"
echo "End time: $(date)"
echo "Output file: cluster-analysis.tar.zst"
ls -lh cluster-analysis.tar.zst
