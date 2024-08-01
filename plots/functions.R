
# Create an R script with the following functions to source to plot_results.Rmd:
#  --> summarySE: 
#       - Description: Summarizes and groups variable in w.r.t other chosen variables
#       - Input: measurevar, groupvars
#       - Output: DF w/ Mean, 95% CI, SE
# 
#  --> FST_from_B:
#       - Description: Calculates WF equilibrium FST using inputted B from subpopulation
#       - Input: B, Ne, m, b
#       - Output: FST
# 
#  --> calculate_B_zc:
#       - Description: Calculates B using the method of Hudson & Kaplan (1995) and Nordborg et al. (1996)
#       - Input: mu, s, r_bp, L
#       - Output: B
#
#  --> calculate_predicted_FST:
#       - Description: Calculates FST using B computed from the method of Hudson & Kaplan (1995) & Nordborg et al. (1996)
#       - Input: N, mu, m, s, r_bp, d, L, b
#       - Output: FST



summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N_replicates    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N_replicates)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N_replicates-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


FST_from_B <- function(B, N, m, d, b = TRUE) {
  
  if (b == TRUE) {
    FST = 1 / (1 + ((d-1)*d*(2-m)*B*N*m / (d*(1-m)-1)^2)) 
  }
  else { 
    FST = 1 / (1 + (d^2*(2-m)*B*N*m / (d*(1-m)-1)^2) ) 
  }
  
  return(FST)
}

calculate_B_hk <- function(mu, s, r_bp, L){
  
  U = L*mu
  R = r_bp*L
  
  return(exp(-2*U/(2*s+R)))
}


calculate_predicted_FST <- function(N, mu, m, s, r_bp, d, L, b = TRUE) {
  
  if (b == TRUE) {
    if (s != 0) {
      predicted_FST = 1 / (1 + ((d-1)*d*(2-m)*calculate_B_hk(mu, s, r_bp, L)*N*m / (d*(1-m)-1)^2)) 
    }
    else {
      predicted_FST = 1 / (1 + ((d-1)*d*(2-m)*N*m / (d*(1-m)-1)^2)) 
    }
  }
  else {
    if (s != 0) {
      #predicted_FST =  (1-m)^2 / ((2-m)*(d^2/(d-1)^2)*N*m*calculate_B_hk(mu, s, r_bp, L) + (1-m)^2)
      predicted_FST = 1 / (1 + (d^2*(2-m)*calculate_B_hk(mu, s, r_bp, L)*N*m / (d*(1-m)-1)^2) ) 
    }
    else {
      #predicted_FST = (1-m)^2 / ((2-ml)*(d^2/(d-1)^2)*N*m + (1-m)^2)
      predicted_FST = 1 / (1 + (d^2*(2-m)*N*m / (d*(1-m)-1)^2) ) 
    }
  }

  return(predicted_FST)
}



