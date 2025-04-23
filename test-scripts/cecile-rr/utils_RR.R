## Get distance matrix from a list of sequences
get_dist_mat_from_fasta <- function(sequence_data){
  
  if(length(unique(sapply(sequence_data, length))) > 1){
    stop('Sequences should be aligned and have the same length!')
  }
  
  dist_mat <- dist.dna(sequence_data, model = 'N',
                       as.matrix = T, pairwise.deletion = T)
  
  return(dist_mat)
}

## Get dataframe with genetic distance between sequences in a fasta file
get_df_pairs_from_fasta <- function(sequence_data){
  dist_mat <- get_dist_mat_from_fasta(sequence_data)
  
  df_pairs <- dist_mat %>% 
    as_tibble() %>% 
    mutate(label_tip_1 = rownames(dist_mat)) %>% 
    pivot_longer(cols = -'label_tip_1', names_to = 'label_tip_2', values_to = 'n_mutations') %>% 
    filter(label_tip_1 != label_tip_2)
  
  return(df_pairs)
}

## Get dataframe with pairs of sequences at a given genetic distance from a fasta file
get_df_pairs_given_genetic_distance_from_fasta <- function(sequence_data, n_mut_away){
  
  df_pairs <- get_df_pairs_from_fasta(sequence_data)
  df_pairs_given_mut <- df_pairs %>% 
    filter(n_mutations == n_mut_away)
  
  return(df_pairs_given_mut)
}

## Get relative risk of observing sequences from df_pairs between the groups defined by name_group
get_df_RR <- function(df_pairs, metadata_seq, name_group){
  
  if(! name_group %in% names(metadata_seq)){
    stop(paste0('ERROR: ', name_group, ' is not a column in the metadata file!'))
  }
  
  metadata_of_interest <- metadata_seq %>% 
    rename(group = name_group) %>% 
    select(sequence_name, group)
  
  df_RR <- df_pairs %>% 
    filter(label_tip_1 %in% metadata_of_interest$sequence_name,
           label_tip_2 %in% metadata_of_interest$sequence_name) %>% 
    left_join(metadata_of_interest, by = c('label_tip_1' = 'sequence_name')) %>% 
    rename(group_1 = group) %>% 
    left_join(metadata_of_interest, by = c('label_tip_2' = 'sequence_name')) %>% 
    rename(group_2 = group) %>% 
    group_by(n_mutations, group_1, group_2) %>% 
    summarise(n_pairs = n(), .groups = 'drop') %>% 
    group_by(n_mutations, group_1) %>% 
    mutate(n_pairs_1_x = sum(n_pairs)) %>% 
    group_by(n_mutations, group_2) %>% 
    mutate(n_pairs_x_2 = sum(n_pairs)) %>% 
    ungroup() %>% 
    mutate(n_pairs_x_x = sum(n_pairs),
           RR = n_pairs / n_pairs_1_x / n_pairs_x_2 *n_pairs_x_x)
  print(df_RR)
  return(df_RR)
}

get_df_uncertainty_RR <- function(df_pairs, metadata_seq, name_group, prop_subsample, n_subsamples){
  
  if(! name_group %in% names(metadata_seq)){
    stop(paste0('ERROR: ', name_group, ' is not a column in the metadata file!'))
  }
  
  n_sequences <- nrow(metadata_seq)
  
  df_subsample_RR <- Reduce('bind_rows', lapply(1:n_subsamples, FUN = function(i_subsample){
    id_to_keep <- sample(1:n_sequences, size = round(n_sequences * prop_subsample), replace = F)
    metadata_seq_subsample <- metadata_seq[id_to_keep, ]
    get_df_RR(df_pairs, metadata_seq_subsample, name_group) %>% 
      mutate(i_subsample = i_subsample)
  }))
  
  return(df_subsample_RR)
}
