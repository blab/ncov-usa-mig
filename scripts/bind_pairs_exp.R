#db_con - Connection for dbplyr, typically will be duckdb_connection
#exp_var - Name of exposure variable to bind to the table
#sub_samp - True if you are planning on using this to subsample for confidence interval generation
bind_pairs_exp <- function(db_con,exp_var,sub_samp=FALSE,samp_cov=0.8){
  meta_tbl <- tbl(db_con,'metadata')
  pairs_tbl <- tbl(db_con,'pairs')

  
  X_exp_var <- paste0("X_",exp_var)
  Y_exp_var <- paste0("Y_",exp_var)
  
  if(sub_samp){ #Slice up the metadata for subsampling
    k <- floor(meta_tbl %>% summarize(n()) %>% collect() %>% as.numeric() * samp_cov)  #Yes this is poor grammar
    meta_tbl <- meta_tbl %>% slice_sample(n=k,replace = FALSE)
  }
  exp_dict <- meta_tbl %>% select(c(strain,sym(exp_var))) #Makes a simple dictionary for binding exposures to
  return(pairs_tbl %>%
    inner_join(exp_dict,join_by(strain_1==strain)) %>%
    rename(x = exp_var) %>% 
    inner_join(exp_dict,join_by(strain_2==strain)) %>%
    rename(y = exp_var) %>% 
    filter(!is.na(x) & !is.na(y)) %>%
    filter(x != "NA" & y != "NA") #In cases where DuckDB turned NA into a character literal
  )
}
