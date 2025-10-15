# This is the main script to run all the Bayesian code. 
# It calls the other R scripts needed to run the Bayesian Code, then exports a csv file



# Whether or not to use Stan for simulations. Else use OpenBUGS. Stan seems to be much faster in the long term
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


# The data set is quite large. If you want to run all rows, use start_index = 1 and end_index = max_index
####### choose index range ###################################################################
start_index = 50003
end_index = 50012
use_random_indices = FALSE     # if TRUE, take a random sample instead of sequential rows
##############################################################################################


num_rows = end_index-start_index+1

if(use_random_indices)
{
  set.seed(300)
  indices = sample(1:50000, num_rows, replace = FALSE) #probably should be num_rows
  
}else 
{
  indices = 1:max_index
}


# initialize score array (faster than updating each time)
bayes_score = data.frame(matrix(nrow=num_rows, ncol=11)) 
names(bayes_score) <- c("Combined Probability", "Lower Confidence Interval", "Upper Confidence Interval", 
                        "DNA Match Probability","DNA Match Lower Confidence Interval","DNA Match Upper Confidence Interval",
                        "Location Probability","Location Lower Confidence Interval","Location Upper Confidence Interval",
                        "Info","Row Index")

start_time_outside <- Sys.time() 
print(paste("Start time:", start_time_outside))

for (j in start_index:end_index)
{
  
  i = indices[j] # allows for use of random indices 
  jj = j -start_index+1 # index starting from 1
  
  percentage_id = data$X.ID[i]/100
  diversity_rate = data$dr_scaled[i]
  read_counts = data$count[i]
  # Assuming GBIF_count and OBIS_count are for species in area, not worldwide
  gbif_count = data$gbif_count_AOI[i]
  obis_count = data$OBIS_count[i]
  p_aquamaps = data$am_prob[i]
  nearest_distance_gbif = data$nearest_GBIF_m[i]
  nearest_distance_obis = data$nearest_obis_m[i]
  percent_genus_in_db = data$Pct_GenusDNA_inDB[i]
  
  
  # put probabilities in these arrays, not raw counts or distances
  location_prob_array = c(p_aquamaps) # location-based probabilities 
  match_prob_array = c(percent_genus_in_db) # other probabilities (match probabilities)

  # Calculate the Bayesian probability
  bayes_score[jj,1:10] = Bayesian_score(percentage_id, diversity_rate, read_counts, gbif_count, obis_count,
                                   nearest_distance_gbif, nearest_distance_obis, distance_for_50_percent=50000,
                                   location_prob_array, match_prob_array, 
                                   confidence_int = 0.95, calculate_location_too = T, use_stan = use_stan, 
                                   num_chains = 1, num_iterations=3000, num_burnin=1000, 
                                   stan_beta_model, stan_counts_model,
                                   show_summary=F, plot_density=F, plot_trace=F, plot_posterior=F,show_time=F, flat_df=T )
  
  
  # save index for future reference
  bayes_score[jj,11] = i
  
  #Indicate progress:
  cat("Completed:", jj,"/",num_rows,"\r")
}

# write to uniquely named file based on current time
time = Sys.time()
print(paste("Completed:", time))

time = gsub("-","_",time) 
time = gsub(":","_",time)
time = gsub(" ","@",time)
time = gsub("[.][0-9]*","",time)
file_name = paste("Bayes_data_row_",start_index,"_to_",end_index,"_Date_time_",time,".csv",sep="")
write.csv(data.frame(data[indices[start_index:end_index],], bayes_score),file_name , row.names = FALSE)
print(paste("Wrote to file:",file_name))


