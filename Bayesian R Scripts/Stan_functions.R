
# This script contains 2 functions: 
#
#  Beta_Stan - run simulation for posterior 
#
#  Counts_and_ID_Stan - run simulation for percentage ID with diversity taken into account, read counts and counts of species in area

# Author: Ray Stokes

Beta_Stan <- function(prob_array, model, confidence_int = 0.95, show_summary=F, plot_density=F, plot_trace=F, plot_posterior=F, label="", show_time=T)
{
  
  # ARGS:
    # prob_array - an array of probabilities that should be taken into account. Can technically be any size, so can keep feeding more probability values into it as see fit
    # model -  the Stan model to use
    # confidence_int - decimal for desired confidence (credible) interval. Default is 0.95 for a 95% confidence interval
    # show_summary - shows the Bayesian estimates and statistics if TRUE 
    # plot_density - whether you want the density plots for the estimated parameters (alpha, beta) to be shown (TRUE/FALSE)
    # plot_trace - whether you want the trace plots to be shown (TRUE/FALSE)
    # plot_posterior - plot the posterior distribution if TRUE
    # label - label to use for posterior plot heading
    # show_time - If TRUE, print how long it took to run the function. Good for diagnostics and optimization
  
  # RETURNS:
    # A data frame with expected values and confidence intervals for the posterior distribution
  
  start_time <- Sys.time() 
  
  data.in <- list(p=prob_array, n=length(prob_array))
  model.fit <- sampling(model, data=data.in, refresh=0)
  m_alpha = mean(extract(model.fit,pars="alpha")[[1]])
  m_beta = mean(extract(model.fit,pars="beta")[[1]])
  mean_p = m_alpha / (m_alpha+m_beta)
  conf_lower = qbeta((1-confidence_int)/2, m_alpha, m_beta)
  conf_higher = qbeta(0.5+ confidence_int/2, m_alpha, m_beta)
  df <- data.frame(mean_p,conf_lower,conf_higher)
  names(df) <- c("Expected Value","Lower Credible Interval", "Upper Credible Interval")
  
  if (show_summary)
  {
    print(model.fit, pars=c("alpha", "beta"), digits=5)
    check_hmc_diagnostics(model.fit)
  }
  if (plot_density)
  {
    posterior <- as.array(model.fit)
    mcmc_dens(posterior, pars=c("alpha", "beta"))
  }
  if(plot_trace)
  {
    posterior <- as.array(model.fit)
    mcmc_trace(posterior, pars=c("alpha", "beta"))
  }
  if (plot_posterior) # same as BUGS
  {
    x = seq(0,1,0.001)
    b = dbeta(x, m_alpha, m_beta)
    plot(x,b,main=label,xlab="Probability",ylab="Density",type='l',lwd=2,col='coral1')
    lines(c(conf_lower,conf_lower), c(0, dbeta(conf_lower,m_alpha, m_beta)),lty = 2,col='aquamarine4')
    lines(c(conf_higher,conf_higher), c(0,dbeta(conf_higher,m_alpha, m_beta)),lty = 2,col='aquamarine4')
    lines(c(mean_p,mean_p), c(0, dbeta(mean_p,m_alpha, m_beta)),lty = 2,col='dodgerblue')
  }
  if (show_time)
  {
    print(paste("Elapsed time for Beta_Stan function:",  round(Sys.time() - start_time,2),"s"))
  }
  
  return(df)
  
}






Counts_and_ID_Stan <- function(percentage_id, diversity_rate, read_counts, num_obis, num_gbif, model, 
                               show_summary=F, plot_density=F, plot_trace=F, show_time=T)
{
  
  # ARGS:
    # percentage_id - percentage match from BLAST
    # diversity_rate - normalized diversification score between 0 and 1. 0 means species is not diverse, 1 means highly diverse, DNA constantly changing
    # read_counts - the number of counts for the DNA in the sample
    # num_in_area - number of species in the general area
    # num_chains - number of simulations to run at a time, 
    # num_iterations - number of iterations to use for OpenBUGS simulation. The default is 11000
    # num_burnin - number of iterations to discard at the start, to let simulation 'warm up' 
    # show_summary - shows the Bayesian estimates and statistics if TRUE 
    # plot_density - whether you want the density plots for the estimated parameters (alpha, beta) to be shown (TRUE/FALSE)
    # plot_trace - whether you want the trace plots to be shown (TRUE/FALSE)
    # show_time - If TRUE, print how long it took to run the function. Good for diagnostics and optimization
  
  # RETURNS:
    # A data frame with expected values for percentage_id (incorporating diversity) with label "ID"
    #  and probabilities associated with read_counts (RD) and num_in_area (NT)
  
  
  
  

  
  start_time <- Sys.time() 
  
  
  data.in <- list(read_counts = read_counts, num_obis = num_obis, num_gbif=num_gbif, 
                  p_id = percentage_id,diversity_rate=diversity_rate )
  model.fit <- sampling(model, data=data.in, refresh=0)
  
  ID = mean(extract(model.fit,pars="ID")[[1]])
  RD = mean(extract(model.fit,pars="RD")[[1]])
  NTg = mean(extract(model.fit,pars="NTg")[[1]])
  NTo = mean(extract(model.fit,pars="NTo")[[1]])
  
  df <- data.frame(ID,RD,NTo, NTg)
  names(df) <- c("ID","RD", "NTo","NTg")
  
  if (show_summary)
  {
    print(model.fit, pars=c("ID", "RD","NTo","NTg"), digits=5)
    check_hmc_diagnostics(model.fit)
  }
  if (plot_density)
  {
    posterior <- as.array(model.fit)
    mcmc_dens(posterior, pars=c("ID", "RD","NTo","NTg"))
  }
  if(plot_trace)
  {
    posterior <- as.array(model.fit)
    mcmc_trace(posterior, pars=c("ID", "RD","NTo","NTg"))
  }
  if (show_time)
  {
    print(paste("Elapsed time for Counts_and_ID_Stan function:",  round(Sys.time() - start_time,2),"s"))
  }
  
  return(df)
  
}



