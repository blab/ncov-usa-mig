#!/bin/bash
#run from parent repository directory

Rscript idseq-flow/id_seq/scripts/cluster_alloc_from_pairsnp.R --df_id_seq test_data/wa_pairs_abridged.tsv --vec_strain_names test_data/test_pairs.csv --cluster_alloc test_alloc.tsv --pango baby_wa
