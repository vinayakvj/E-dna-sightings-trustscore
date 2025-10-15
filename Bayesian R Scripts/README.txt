This directory contains files to run Bayesian Code

Running the script for Bayesian_main.R will run the following scripts:

	- Bayesian_score.R
	- Bugs_functions.R
	- Load_Bugs_Code_Script.R
	- Stan_functions
	- Load_Stan_Code_Script.R
	
You will also require the "data_enriched.csv" file

In the section "choose index range", you can choose starting and ending index
For example, if you want to calculate Bayesian Scores for rows 1001 to 2000, set:
	start_index = 1001
	end_index = 2000

If you want to take a random sample, set use_random_indices to TRUE


The function that does most of the calculations is the Bayesian_score function. Some of the parameters are:
	- confidence_int: Specify the confidence (credible) interval, for example 0.95 for a 95% confidence interval
	- calculate_location_too: If TRUE, calculate the posterior probabilities as well. It would be best to leave this as is, 
					that is how the bayes_score data frame is structured
	- plot_posterior: If TRUE, plot the posterior distributions. Best to only use for 1 or 2 rows of data
