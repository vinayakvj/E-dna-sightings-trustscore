

# This is the main script to run all the Bayesian code. The code from here can be copied into the main .R file we submit to the client. 
# It calls the other R scripts needed to run the Bayesian Code


# Whether or not to use Stan for simulations. Else use OpenBUGS
use_stan = T



if (use_stan) {
  library(rstan)
  library(bayesplot)
  source("Load_Stan_Code_Script.R")
  source("Stan_functions.R")
  stan_beta_model <- Stan_beta_load()
  stan_counts_model <- Stan_counts_load()
  
}else {
  # Load libraries and scripts so can use their functions
  library(R2OpenBUGS)
  library(coda)
  source("Load_Bugs_Code_Script.R")
  source("Bugs_functions.R")
  
  # Load the Bugs code
  Load_Bugs()
  stan_beta_model <- NULL
  stan_counts_model <- NULL
}


# Load source code for bayesian calculations
source("Bayesian_score.R")



#### Run the above code once, but the code below should be run for each row:



# Say we have this data from databases and all_voyages:

  # compulsory:
  percentage_id = 0.9
  diversity_rate = 0.3
  read_counts = 100
  num_in_area = 4
  nearest_distance = 30

  # extra data
  p_assay = 0.8
  percent_genus_db = 0.6
  p_aquamaps = 0.8
  # can add more here ...

  # put the extra variables into an array (can add even more if needed)
  location_prob_array = c(p_aquamaps) # location-based probabilities 
  match_probability_array = c(p_assay, percent_genus_db) # other probabilities (match probabilities)

  
# Calculate the Bayesian probability
score = Bayesian_score(percentage_id, diversity_rate, read_counts, num_in_area, 
                           nearest_distance, distance_for_50_percent=50,
                           location_prob_array, match_probability_array, confidence_int = 0.95, 
                           calculate_location_too = F, use_stan = use_stan, 
                           num_chains = 4, num_iterations=11000, num_burnin=1000, 
                           stan_beta_model, stan_counts_model,
                           show_summary=F, plot_density=F, plot_trace=F, plot_posterior=F,show_time=F )

score
