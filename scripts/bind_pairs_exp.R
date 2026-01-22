#File: bind_pairs_exp.R
#Author(s): Amin Bemanian
#Date: 07/07/25
#Description: Attach exposure demes to pairs for RR calculation
#Arguments:
#db_con: DuckDB Connection
#exp_var: Name of exposure variable to attach
#sub_samp: Sub-sample (for CI calculation)
#samp_cov: Sample coverage (for CI calculation)
#time_bounds: Date range (as vector of two Date values) for time analyses
#exclude_duplicates: If TRUE, exclude pairs marked as possible_duplicates (default: FALSE)

bind_pairs_exp <- function(db_con,exp_var,sub_samp=FALSE,samp_cov=0.8,time_bounds=NULL,exclude_duplicates=FALSE){
  meta_tbl <- tbl(db_con,'metadata')
  
  #For time analysis use the pairs_time table instead
  if(is.null(time_bounds)){
    pairs_tbl <- tbl(db_con,'pairs')  
  }else{
    pairs_tbl <- tbl(db_con,"pairs_time")
  }

  X_exp_var <- paste0("X_",exp_var)
  Y_exp_var <- paste0("Y_",exp_var)
  
  if(sub_samp){ #Slice up the metadata for subsampling
    k <- floor(meta_tbl %>% summarize(n()) %>% collect() %>% as.numeric() * samp_cov)  #Yes this is poor grammar
    meta_tbl <- meta_tbl %>% slice_sample(n=k,replace = FALSE) %>% compute()  # Materialize to avoid re-sampling on repeated access
  }
  if(!is.null(time_bounds)){
    time_LB <- time_bounds[1]
    time_UB <- time_bounds[2]
    pairs_tbl <- pairs_tbl %>%
      filter(transmit_date > time_LB) %>%
      filter(transmit_date < time_UB)
  }
  if(exclude_duplicates){
    pairs_tbl <- pairs_tbl %>%
      filter(!possible_duplicates | is.na(possible_duplicates))
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
