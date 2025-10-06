

# Main code block for Bayesian_score function to calculate the 
# Easier to keep this function in a separate file and call it

# Author: Ray Stokes

Bayesian_score <- function(percentage_id, diversity_rate, read_counts, num_in_area, 
                           nearest_distance, distance_for_50_percent=50,
                           location_prob_array, match_probability_array, confidence_int = 0.95, 
                           calculate_location_too = F, use_stan = F, 
                           num_chains = 4, num_iterations=11000, num_burnin=1000, 
                           stan_beta_model, stan_counts_model,
                           show_summary=F, plot_density=F, plot_trace=F, plot_posterior=F,show_time=F )
{
  # ARGS:
    # percentage_id - percentage match from BLAST
    # diversity_rate - normalized diversification score between 0 and 1. 0 means species is not diverse, 1 means highly diverse, DNA constantly changing
    # read_counts - the number of counts for the DNA in the sample
    # num_in_area - number of species in the general area
    # nearest_distance - closest distance between observation and another observation of the species in the database
    # distance_for_50_percent - distance where we assume there is a 0.5 probability of the species being in location
    # location_prob_array - an array of location based probabilities that should be taken into account. For aquamaps probability, but can be any size so can add more values
    # match_probability_array -  an array of match based probabilities that should be taken into account, for example, probability associated with reliability of assay used, genus coverage in database, etc. The array can be any size
    # confidence_int - decimal for desired confidence (credible) interval. Default is 0.95 for a 95% confidence interval
    # calculate_location_too - Set to TRUE when want to compare posterior against match and location posteriors. Good for determining if species has moved location, but may take more computational time. 
    # use_stan - Set TRUE if want to use Stan to do simulations. Else BUGS will be used with the parameters below:
    # num_chains - number of simulations to run at a time in OpenBUGS (default is 4)
    # num_iterations - number of iterations to use for OpenBUGS simulation. The default is 11000
    # num_burnin - number of iterations to discard at the start, to let simulation 'warm up' (also for OpenBUGS)
    # show_summary - shows the Bayesian estimates and statistics if TRUE 
    # plot_density - whether you want the density plots for the estimated parameters (alpha, beta) to be shown (TRUE/FALSE)
    # plot_trace - whether you want the trace plots to be shown (TRUE/FALSE)
    # show_time - If TRUE, print how long it took to run the function. Good for diagnostics and optimization
  
  # RETURNS:
    # A data frame with expected values and confidence intervals for posterior probability
    # if calculate_location_too is set to TRUE, it will calculate the posteriors for match and location and add them to this dataframe 
    
  
      start_time_overall <- Sys.time() 
    
      # auto correct range for diversity rate
      if (diversity_rate < 0.001)
      {
        diversity_rate=0.001
      }
      
      if (use_stan)
      {
         counts_df = Counts_and_ID_Stan(percentage_id, diversity_rate, read_counts, num_in_area, stan_counts_model, 
                                        show_summary=show_summary, plot_density=plot_density, plot_trace=plot_trace, show_time=show_time)
      }else {
      
        counts_df = Counts_and_ID_BUGS(percentage_id, diversity_rate, read_counts, num_in_area, 
                         num_chains = num_chains, num_iterations=num_iterations, num_burnin=num_burnin, 
                         show_summary=show_summary, plot_density=plot_density, plot_trace=plot_trace, show_time=show_time)
      } 
      
      ID = counts_df$ID
      RD = counts_df$RD
      NT = counts_df$NT
  
      if (nearest_distance==0)
      {
        nearest_distance=0.01
      }
      prob_dist = exp(-nearest_distance*log(2) / distance_for_50_percent)
      
      # combine match-based probabilities, and combine location probabilities
      match_probability_array = c(match_probability_array, ID, RD)
      location_prob_array = c(location_prob_array, NT,prob_dist)
      
      # Simulations don't go so well for values too close to 0 or 1, so restrict the range
      Constrain_Range <- function(array)
      {
        for (i in 1:length(array))
        {
          if (array[i] > 0.9999)
          {
            array[i] = 0.9999
          }
          else if (array[i] < 0.0001)
          {
            array[i] = 0.0001
            
          }
          
        }
        return(array)
      }
      match_probability_array = Constrain_Range(match_probability_array)
      location_prob_array = Constrain_Range(location_prob_array)
      
      # Combine both for the posterior
      total_prob_array = c(match_probability_array,location_prob_array)
      
      if (use_stan) {
        df = Beta_Stan(total_prob_array, stan_beta_model, confidence_int = confidence_int, 
                       show_summary=show_summary, plot_density=plot_density, plot_trace=plot_trace, plot_posterior=plot_posterior, label="Posterior Probability", show_time=show_time)
      }else {  
        df = Beta_BUGS(total_prob_array, confidence_int = confidence_int, num_chains = num_chains, num_iterations=num_iterations, num_burnin=num_burnin, 
                  show_summary=show_summary, plot_density=plot_density, plot_trace=plot_trace, plot_posterior=plot_posterior, label="Posterior Probability", show_time=show_time)  
      }
        
      if (calculate_location_too)
      {
        if (use_stan)
        {
          df_match = Beta_Stan(match_probability_array, stan_beta_model, confidence_int = confidence_int, 
                               show_summary=show_summary, plot_density=plot_density, plot_trace=plot_trace, plot_posterior=plot_posterior, label="Match Probability", show_time=show_time)
          df_location = Beta_Stan(location_prob_array, stan_beta_model, confidence_int = confidence_int, 
                                  show_summary=show_summary, plot_density=plot_density, plot_trace=plot_trace, plot_posterior=plot_posterior, label="Location Probability", show_time=show_time)
        }else  {
        
          df_match = Beta_BUGS(match_probability_array, confidence_int = confidence_int, num_chains = num_chains, num_iterations=num_iterations, num_burnin=num_burnin, 
                                show_summary=show_summary, plot_density=plot_density, plot_trace=plot_trace, plot_posterior=plot_posterior, label="Match Probability", show_time=show_time)
          df_location = Beta_BUGS(location_prob_array, confidence_int = confidence_int, num_chains = num_chains, num_iterations=num_iterations, num_burnin=num_burnin, 
                                show_summary=show_summary, plot_density=plot_density, plot_trace=plot_trace, plot_posterior=plot_posterior, label="Location Probability", show_time=show_time)
        }
        
        df = rbind(df,df_match,df_location)
        df = data.frame(Probability = c("Combined","Match","Location"), df)
        df$Info = paste("Using ", c(confidence_int*100,confidence_int*100,confidence_int*100), "% credible interval")
      }
      
      if (show_time)
      {
        print(paste("Elapsed time for Bayesian_score function:",  round(Sys.time() - start_time_overall,2),"seconds"))
      }
      
      
        
  return(df)
  
}