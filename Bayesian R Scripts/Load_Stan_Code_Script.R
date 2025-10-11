# The purpose of this script is to write the Stan code to file so it can be used
# Author: Ray Stokes

Stan_beta_load <- function()
{
  
  code <- 
  "data{
    int<lower=1> n;
    vector[n] p;
  }

  parameters{
    real<lower=1,upper=100> alpha;
    real<lower=1,upper=100> beta;
  }


  model{
    // likelihood
    p ~ beta(alpha, beta);
  }"
  
  model <- stan_model(model_code = code)
  return(model)
}


Stan_counts_load <- function()
{
  
  code <- 
    "data{
      real<lower=0.1> read_counts;
      real<lower=0.1> num_obis;
      real<lower=0.1> num_gbif;
      real<lower=0,upper=1> p_id;
      real<lower=0,upper=1> diversity_rate;
    }

    parameters{
      real<lower=0,upper=1> p_rd;
      real<lower=0,upper=1> a_rd;
      real<lower=0,upper=1> p_obis;
      real<lower=0,upper=1> p_gbif;
      real<lower=0,upper=1> k;
    }


    transformed parameters{
      real<lower=0,upper=1> RD;
      real<lower=0,upper=1> NTo;
      real<lower=0,upper=1> NTg;
      real<lower=0,upper=1> ID;
      
      RD = 1 - pow( (1 - p_rd), (a_rd * read_counts) );
      NTo = 1 - pow( (1 - p_obis), num_obis);
      NTg = 1 - pow( (1 - p_gbif), num_gbif);
      ID = pow(p_id, 1-k*diversity_rate);
    }"
  
  model <- stan_model(model_code = code)
  return(model)
}







