# This is the main script to run all the Bayesian code. The code from here can be copied into the main .R file we submit to the client. 
# It calls the other R scripts needed to run the Bayesian Code, then exports a csv file



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


data <- read.csv("data_enriched.csv")
max_index = nrow(data)

# choose indices to start from 
start_index = 1
end_index = 50000
num_rows = end_index-start_index+1
use_random_indices = TRUE # take a random sample instead of sequential rows

if(use_random_indices)
{
  set.seed(300)
  indices = sample(1:max_index, num_rows, replace = FALSE)
  
}else 
{
  indices = 1:max_index
}


# initialize score array (faster than updating each time)
bayes_score = data.frame(matrix(nrow=num_rows, ncol=11)) 
names(bayes_score) <- c("Posterior Mean", "Lower Confidence Interval", "Upper Confidence Interval", 
                        "Match Posterior Mean","Match Lower Confidence Interval","Match Upper Confidence Interval",
                        "Location Posterior Mean","Location Lower Confidence Interval","Location Upper Confidence Interval",
                        "Info","Row Index")

start_time_outside <- Sys.time() 
print(paste("Start time:", start_time_outside))

for (j in start_index:end_index)
{
  
  i = indices[j] # allows for use of random indices 
  
  percentage_id = data$X.ID[i]/100
  diversity_rate = data$dr_scaled[i]
  read_counts = data$count[i]
  # Assuming GBIF_count and OBIS_count are for species in area, not worldwide
  gbif_count = data$GBIF_count[i]
  obis_count = data$OBIS_count[i]
  p_aquamaps = data$am_prob[i]
  nearest_distance_gbif = data$nearest_GBIF_m[i]
  nearest_distance_obis = data$nearest_obis_m[i]
  
  # put probabilities in these arrays, not raw counts or distances
  location_prob_array = c(p_aquamaps) # location-based probabilities 
  match_prob_array = c() # other probabilities (match probabilities)
 
  
  # Calculate the Bayesian probability
  bayes_score[j,1:10] = Bayesian_score(percentage_id, diversity_rate, read_counts, gbif_count, obis_count,
                                   nearest_distance_gbif, nearest_distance_obis, distance_for_50_percent=50000,
                                   location_prob_array, match_prob_array, 
                                   confidence_int = 0.95, calculate_location_too = T, use_stan = use_stan, 
                                   num_chains = 1, num_iterations=3000, num_burnin=1000, 
                                   stan_beta_model, stan_counts_model,
                                   show_summary=F, plot_density=F, plot_trace=F, plot_posterior=F,show_time=F, flat_df=T )
  
  
  # save index for future reference
  bayes_score[j,11] = i
}

# write to uniquely named file based on current time

time = gsub("-","_",Sys.time()) 
time = gsub(":","_",time)
time = gsub(" ","@",time)
time = gsub("[.][0-9]*","",time)
file_name = paste("Bayes_data_",time,".csv",sep="")
write.csv(data.frame(data[indices[start_index:end_index],], bayes_score),file_name , row.names = FALSE)
print(paste("Wrote to file:",file_name))


