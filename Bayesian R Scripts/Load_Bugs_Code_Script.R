
# The purpose of this script is to write the OpenBUGS code to file so it can be used
# Author: Ray Stokes


Load_Bugs <- function()
{
  
  writeLines("
    model{
    
      for (i in 1:n) {
        p[i] ~ dbeta(alpha, beta)
      }
      
      alpha ~ dunif(1,100)
      beta ~ dunif(1,100)
    
    }", con="BUGS_beta.txt")
  
  
  writeLines("
    model{
      # First find probabilities from read counts 
      RD <- 1 - pow( (1 - p_rd), (a_rd * read_counts) )
      p_rd ~ dunif(0, 1)
      a_rd ~ dunif(0, 1)

      # Probability for number of species found in target area
      NT <- 1 - pow( (1 - p_ta), num_in_area)
      p_ta ~ dunif(0, 1)
      
      # Probability of species ID given diversity rate
      ID <- pow(p_id, 1-k*diversity_rate)
      k ~ dunif(0,1)
      
      
    }", con="Bugs_ID_counts.txt")
  

}









