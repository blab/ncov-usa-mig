#Wrapper function for tsv_join function to wrap whatever exposure of interest to the pairs file
#Saves the data frame for reference and returns the dataframe
bind_pairs_exp <- function(fn_m, fn_p, exp_var,fn_o){
  #Note had to make a temporary file because tsv_join would treat the piped variable as the filter instead of the TSV to join to
  system(
    paste0( #Assumes that the pairs fields are always in the format of strain_1 and strain_2 
      "tsv-join -H --filter-file ",fn_m," --key-fields strain_1 --append-fields ",exp_var," ",fn_p," --prefix X_ > temp"
    )
  )
  system(
    paste0(
      "tsv-join -H --filter-file ",fn_m," --key-fields strain_2 --append-fields ",exp_var," temp --prefix Y_ > ",fn_o  
    )
  )
  system("rm temp") #Clean up the temporary file
  return(fread(fn_o))
}