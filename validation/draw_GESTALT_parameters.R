## ---------------------------
##
## Script name: draw GESTALT parameters
##
## Purpose of script: draw parameters for the GESTALT mutation process for the validation
##
## Author: Antoine Zwaans
##
## Date Created: 2024-03-28
##
## Copyright (c) Antoine Zwaans, 2024
## Email: antoine.zwaans@bsse.ethz.ch


options(scipen=5)

## specify output location
outputDir = paste0("~/GABIanalysis/validation/simulation_parameters/")

if (! dir.exists(outputDir)){
  dir.create(outputDir)
}

## draw clock rates
clock_rates <- c()
for (seed in 1:650){

  # draw clock on seed
  set.seed(seed)
  clock_rates = c(clock_rates,round(rlnorm(n = 1, meanlog = -4, sdlog = 0.2), digits = 5))

}
clock_rates <- cbind(1:650,clock_rates)
colnames(clock_rates) = c("seed","clock_rate")

write.csv(x = clock_rates, file = paste0(outputDir, "simClockRates.csv"),quote = F, row.names = F)

## draw double cut weight
double_cut_weight <- c()
for (seed in 1:650){
  
  # draw double_cut on seed
  set.seed(seed)
  double_cut_weight = c(double_cut_weight,round(rlnorm(n = 1, meanlog = -3.1, sdlog = 0.2), digits = 5))
}
double_cut_weight <- cbind(1:650,double_cut_weight)

colnames(double_cut_weight) = c("seed","double_cut_rate")

write.csv(x = double_cut_weight, file = paste0(outputDir, "simDoubleWeights.csv"),quote = F, row.names = F)


## draw target cut rates
cut_rates <- c()
for (seed in 1:650){
  
  # draw 4x cut rates on seed
  set.seed(seed)
  cut_rates = rbind(cut_rates,c(c(round(rlnorm(n = 1, meanlog = 0, sdlog = 0.2), digits = 5),round(rlnorm(n = 1, meanlog = 0, sdlog = 0.2), digits = 5),round(rlnorm(n = 1, meanlog = 0, sdlog = 0.2), digits = 5),round(rlnorm(n = 1, meanlog = 0, sdlog = 0.2), digits = 5))))

  
}
cut_rates <- cbind(1:650,cut_rates)
colnames(cut_rates) = c("seed",paste0("cut_rate",0:9))

write.csv(x = cut_rates, file = paste0(outputDir, "simcutRates.csv"),quote = F, row.names = F)


## draw scaling factors on long deletions
longTrimScaling_factor <- c()
for (seed in 1:650){
  
  # draw scaling factors on seed
  set.seed(seed)
  longTrimScaling_factor = rbind(longTrimScaling_factor,round( rbeta(n=2,shape1=5,shape2=80), digits = 5))
  
  
}
longTrimScaling_factor <- cbind(1:650,longTrimScaling_factor)
colnames(longTrimScaling_factor) = c("seed",paste0("trimScalingFactor",1:2))

write.csv(x = longTrimScaling_factor, file = paste0(outputDir, "simlongTrimScaling.csv"),quote = F, row.names = F)

