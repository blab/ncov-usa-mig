T <- 14
exit_prob <- 1/7

transmission_prob <- function(t){
  a = 1
  b = 4
  return(rbeta(1,a,b))
}

num_contacts <- function(){rnbinom(1,size=0.2,mu=2)}

sim_indiv <- function(){
  total_inf <- 0
  
  still_infected <- TRUE
  while(still_infected){
    n <- num_contacts()
    if(n > 0){
      infected <- sum(runif(n) < transmission_prob(t))
      total_inf <- total_inf + infected
    }
    if(runif(1) < exit_prob){
      still_infected <- FALSE
    }
  }
  return(total_inf)
}

sim_pop <- function(N){
  x <- NULL
  for(i in 1:N){
    x <- c(x,sim_indiv())
  }
  return(x)
}

Z <- sim_pop(10000)
hist(Z)
fit <- MASS::fitdistr(Z, "negative binomial")
fit$estimate