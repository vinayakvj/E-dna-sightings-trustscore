
# This script contains 2 functions: 
#
#  Beta_BUGS - run simulation for posterior 
#
#  Counts_and_ID_BUGS - run simulation for percentage ID with diversity taken into account, read counts and counts of species in area

# Author: Ray Stokes

Beta_BUGS <- function(prob_array, confidence_int = 0.95, num_chains = 4, num_iterations=11000, num_burnin=1000, 
                      show_summary=F, plot_density=F, plot_trace=F, plot_posterior=F, label="", show_time=T)
{
  
 # ARGS:
    # prob_array - an array of probabilities that should be taken into account. Can technically be any size, so can keep feeding more probability values into it as see fit
    # confidence_int - decimal for desired confidence (credible) interval. Default is 0.95 for a 95% confidence interval
    # num_chains - number of simulations to run at a time, 
    # num_iterations - number of iterations to use for OpenBUGS simulation. The default is 11000
    # num_burnin - number of iterations to discard at the start, to let simulation 'warm up' 
    # show_summary - shows the Bayesian estimates and statistics if TRUE 
    # plot_density - whether you want the density plots for the estimated parameters (alpha, beta) to be shown (TRUE/FALSE)
    # plot_trace - whether you want the trace plots to be shown (TRUE/FALSE)
    # plot_posterior - plot the posterior distribution if TRUE
    # label - label to use for posterior plot heading
    # show_time - If TRUE, print how long it took to run the function. Good for diagnostics and optimization
  
  # RETURNS:
    # A data frame with expected values and confidence intervals for the posterior distribution
  
  start_time <- Sys.time() 
  res <- bugs(data=list(p = prob_array, n=length(prob_array)),
              inits = NULL, 
              n.chains = num_chains, n.iter = num_iterations, n.burnin = num_burnin,
              parameters.to.save = c("alpha","beta"),
              model.file = "BUGS_beta.txt",
              DIC = FALSE, codaPkg = TRUE)#,debug=T)
  
  codaobj <- read.bugs(res, quiet=TRUE)
  dist <- as.data.frame(as.matrix(codaobj))
  m_alpha = mean(dist$alpha)
  m_beta = mean(dist$beta)
  mean_p = m_alpha / (m_alpha+m_beta)
  conf_lower = qbeta((1-confidence_int)/2, m_alpha, m_beta)
  conf_higher = qbeta(0.5+ confidence_int/2, m_alpha, m_beta)
  df <- data.frame(mean_p,conf_lower,conf_higher)
  names(df) <- c("Expected Value","Lower Credible Interval", "Upper Credible Interval")
  
  if (show_summary)
  {
    print(summary(codaobj)) 
  }
  if (plot_density)
  {
    plot(codaobj, trace=FALSE)
  }
  if(plot_trace)
  {
    plot(codaobj, density=FALSE)
  }
  if (plot_posterior)
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
    print(paste("Elapsed time for Beta_BUGS function:",  round(Sys.time() - start_time,2),"s"))
  }
  
  return(df)
  
}






Counts_and_ID_BUGS <- function(percentage_id, diversity_rate, read_counts, num_in_area, num_chains = 4, num_iterations=11000, num_burnin=1000, 
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
  res <- bugs(data=list(read_counts = read_counts,num_in_area = num_in_area, 
                        p_id = percentage_id,diversity_rate=diversity_rate ),
              inits = NULL, 
              n.chains = num_chains, n.iter = num_iterations, n.burnin = num_burnin,
              parameters.to.save = c("RD","NT","ID"),
              model.file = "Bugs_ID_counts.txt",
              DIC = FALSE, codaPkg = TRUE)#,debug=T)
  
  codaobj <- read.bugs(res, quiet=TRUE)
  dist <- as.data.frame(as.matrix(codaobj))
  df <- data.frame(mean(dist$ID),mean(dist$RD),mean(dist$NT))
  names(df) <- c("ID","RD", "NT")
  
  if (show_summary)
  {
    print(summary(codaobj)) 
  }
  if (plot_density)
  {
    plot(codaobj, trace=FALSE)
  }
  if(plot_trace)
  {
    plot(codaobj, density=FALSE)
  }
  if (show_time)
  {
    print(paste("Elapsed time for Counts_and_ID_BUGS function:",  round(Sys.time() - start_time,2),"s"))
  }
  
  return(df)
  
}