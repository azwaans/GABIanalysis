## ---------------------------
##
## Script name: get summary data from all log files
##
## Purpose of script: Analyse well-calibrated simulations
##
## Author: Antoine Zwaans (based on S.Seidel's original script)
##
## Date Created: 2024-02-09
##
## Copyright (c) Antoine Zwaans
## Email: antoine.zwaans@bsse.ethz.ch / sophie.seidel@posteo.de
## 
##
## ---------------------------


library(tracerer)
library(HDInterval)
library(ggplot2)
library(tidyr)

## ---------------------------
log_dir = "~/GABIanalysis/validation/validation/log/"

simulation_dir = "~/GABIanalysis/validation/validation/simulation_parameters/"


inc = function(x){
  eval.parent(substitute(x <- x+1))
}

is_parameter_in_hpd = function(hpd_lower, hpd_upper, true_parameter){
  if(true_parameter >= hpd_lower && true_parameter <= hpd_upper)
    return(TRUE)
  else
    return(FALSE)
}

get_bias_rel = function(posterior_median, truth){
  
  if (truth == 0){
    bias = posterior_median - truth
  }else{
    bias = (posterior_median - truth) / truth
  }
  
  return(bias)
}

# Gets the root mean square error relative to the true parameter
get_RMSE_rel = function(posterior_median, truth){
  
  if(truth == 0){
    RMSE = sqrt(mean((posterior_median - truth)^2))
    
  }else{
    RMSE = sqrt(mean((posterior_median - truth)^2))/truth
    
  }
  
  return(RMSE)
}

# Gets the HPD width relative to the true parameter
get_hpd_width_rel = function(hpd_lower, hpd_upper, truth){
  
  if (truth == 0){
    
    width = hpd_upper - hpd_lower
    
  }else{
    width = (hpd_upper - hpd_lower)/truth
  }
  
  return(width)
}

##############################
#  plot the double cut weight coverage results
##############################

nr_converged_chains = 0

cut_rate_inference = data.frame(seed=1:100, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F,cutRate= rep(1,100))

recovered_per_seed = c()

for (seed in 1:100){
  print("Seed: ")
  print(seed)
  # get inference log
  log_file = paste0("validate_GESTALT.", seed, ".log")
  log_data = parse_beast_tracelog_file(paste0(log_dir, log_file))
  log_data_wo_burnin = remove_burn_ins(log_data, burn_in_fraction = 0.3)
  
  # get true parameter from simulation
  simulation_parameter_file = paste0(simulation_dir, "simDoubleWeights.csv")
  simulation_parameters = read.csv(simulation_parameter_file)
  colnames(simulation_parameters) = c("parameter", "value")
  
  #check that chains has converged
  esses = calc_esses(log_data_wo_burnin, sample_interval = 1)
  esses = esses[!colnames(esses) %in% c("birthRate", "deathRate", "prior")]
  esses = esses[! names(esses) %in% c("samplingProportion", "treeHeight.t.alignment", "treeLength.t.alignment")]
  esses
  
  if (any(esses < 200)){
    problematic_ess = esses[which(esses < 200)]
    print(paste("Ess for log file ", log_file, " has ESS < 200: for parameters", paste(names(problematic_ess), collapse = ",")))
    #next
    
  }else{
    inc(nr_converged_chains)
  }
  
  recovered_per_alpha = c()
  
  for (cut_rate_nr in 1:1){
    
    #check HPD recovery rate per credibility level
    for(alpha in seq(5,95,by=5)) {
      
      cut_rate = paste0("doubleCutWeight", "")
      
      cut_rate_hpd = hdi(object = log_data_wo_burnin[, cut_rate], credMass = alpha/100)
      
      median_cut_rate = median(log_data_wo_burnin[, cut_rate])
      
      true_cut_rate = simulation_parameters[seed,cut_rate_nr + 1]
      recovered = is_parameter_in_hpd(hpd_lower = cut_rate_hpd["lower"], hpd_upper = cut_rate_hpd["upper"],true_parameter = true_cut_rate)
      print(recovered)
      recovered_per_alpha <- c(recovered_per_alpha,recovered)
      
      row_index = which(cut_rate_inference$seed == seed & cut_rate_inference$cutRate == cut_rate_nr)
      
      if(alpha/100 == 0.95) {
        cut_rate_inference[row_index, ] = c(seed, as.numeric(cut_rate_hpd), median_cut_rate, true_cut_rate,recovered, cut_rate_nr)
      }
    }
    
  }
  
  recovered_per_seed <- rbind(recovered_per_seed,recovered_per_alpha)
}

#coverage per insert rate
sum(cut_rate_inference$recovered)

# get the coverage per credibility level
calibration <- data.frame(credibility=rep(seq(5,95,by=5),1),recovered=colSums(recovered_per_seed))

# generate CI for the binomial distribution corresponding
binomial_intervals_per_credibility <- c() 
for(i in seq(0.05,0.95,by=0.05)) {
  binomial_intervals_per_credibility <- rbind(binomial_intervals_per_credibility,qbinom(c(0.025,0.975), 100, i, lower.tail = TRUE, log.p = FALSE))
}
binomial_intervals_per_credibility <- cbind(seq(5,95,by=5),binomial_intervals_per_credibility)
colnames(binomial_intervals_per_credibility) <- c("credibility_level","min","max")
binomial_intervals_per_credibility <- data.frame(binomial_intervals_per_credibility)
calibration$title <- as.character(expression(paste(omega, " - Double cut weight")))

g1 = ggplot(data=calibration, aes(x=credibility, y=recovered, fill="black")) + geom_point(size=2) + geom_abline(slope = 1.0, col="darkgreen") + xlab("Credibility level (%)") + ylab("Recovered proportion (%)") +
  geom_line(data= binomial_intervals_per_credibility, aes(x=credibility_level,y=min,linetype="dashed")) + geom_line(data= binomial_intervals_per_credibility, aes(x=credibility_level,y=max),linetype="dashed") + theme_bw()

g1 <- g1 + facet_grid(. ~ title,labeller = label_parsed) + theme(strip.text = element_text(size = 14), legend.position = 'none') 
glegend <- g1 + scale_linetype_manual(name = '', 
                          values =c('dashed'='dashed'), labels = c('95% CI of the Binomial distribution')) + scale_fill_manual(name = '', 
                           values =c('black'='black'), labels = c('Proportion of inferred HPDs containing true value'))

legend <- get_legend(
  glegend) 


ggsave(plot = g1, "~/GABIanalysis/validation/plots/calibration_double_cut_weight.pdf",
       width = 12, height = 12, units = "cm")


#getting the bias and RMSE
biases <- c()
RMSEs <- c()
widths <- c()
for (seed in 1:100){
  biases <- c(biases,get_bias_rel(cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "median"],cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "true_value"]))
  RMSEs <- c(RMSEs,get_RMSE_rel(cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "median"],cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "true_value"]))
  widths <- c(widths,get_hpd_width_rel(cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "hpd_lower"],cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "hpd_upper"],cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "true_value"]))
  
  
}

##############################
#  plot the clock rate coverage results
##############################

nr_converged_chains = 0

cut_rate_inference = data.frame(seed=1:100, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F,cutRate= rep(1,100))

recovered_per_seed = c()

for (seed in  1:100){
  print("Seed: ")
  print(seed)
  # get inference log
  log_file = paste0("validate_GESTALT.", seed, ".log")
  log_data = parse_beast_tracelog_file(paste0(log_dir, log_file))
  log_data_wo_burnin = remove_burn_ins(log_data, burn_in_fraction = 0.3)
  
  # get true parameter from simulation
  simulation_parameter_file = paste0(simulation_dir, "simClockRates.csv")
  simulation_parameters = read.csv(simulation_parameter_file)
  colnames(simulation_parameters) = c("parameter", "value")
  
  #check that chains has converged
  esses = calc_esses(log_data_wo_burnin, sample_interval = 1)
  esses = esses[!colnames(esses) %in% c("birthRate", "deathRate", "prior")]
  esses = esses[! names(esses) %in% c("samplingProportion", "treeHeight.t.alignment", "treeLength.t.alignment")]
  esses
  
  if (any(esses < 200)){
    problematic_ess = esses[which(esses < 200)]
    print(paste("Ess for log file ", log_file, " has ESS < 200: for parameters", paste(names(problematic_ess), collapse = ",")))
    #next
    
  }else{
    inc(nr_converged_chains)
  }
  
  
  recovered_per_alpha = c()
  
  for (cut_rate_nr in 1:1){
    
    #check HPD recovery rate per credibility level
    for(alpha in seq(5,95,by=5)) {
      
      cut_rate = paste0("clockRate", "")
      
      cut_rate_hpd = hdi(object = log_data_wo_burnin[, cut_rate], credMass = alpha/100)
      
      median_cut_rate = median(log_data_wo_burnin[, cut_rate])
      
      true_cut_rate = simulation_parameters[seed,cut_rate_nr + 1]
      
      recovered = is_parameter_in_hpd(hpd_lower = cut_rate_hpd["lower"], hpd_upper = cut_rate_hpd["upper"],true_parameter = true_cut_rate)
      
      recovered_per_alpha <- c(recovered_per_alpha,recovered)
      
      row_index = which(cut_rate_inference$seed == seed & cut_rate_inference$cutRate == cut_rate_nr)
      
      if(alpha/100 == 0.95) {
        cut_rate_inference[row_index, ] = c(seed, as.numeric(cut_rate_hpd), median_cut_rate, true_cut_rate,recovered, cut_rate_nr)
      }
    }
    
  }
  
  recovered_per_seed <- rbind(recovered_per_seed,recovered_per_alpha)
}

#coverage at the 95% level
sum(cut_rate_inference$recovered)

# get coverage per credibility level
calibration <- data.frame(credibility=rep(seq(5,95,by=5),1),recovered=colSums(recovered_per_seed))

binomial_intervals_per_credibility <- c() 
for(i in seq(0.05,0.95,by=0.05)) {
  binomial_intervals_per_credibility <- rbind(binomial_intervals_per_credibility,qbinom(c(0.025,0.975), 100, i, lower.tail = TRUE, log.p = FALSE))
}
binomial_intervals_per_credibility <- cbind(seq(5,95,by=5),binomial_intervals_per_credibility)
colnames(binomial_intervals_per_credibility) <- c("credibility_level","min","max")
binomial_intervals_per_credibility <- data.frame(binomial_intervals_per_credibility)
calibration$title <- as.character(expression(paste(r, " - Molecular clock rate")))

g2 = ggplot(data=calibration, aes(x=credibility, y=recovered)) + geom_point(size=2) + geom_abline(slope = 1.0, col="darkgreen") + xlab("Credibility level (%)") + ylab("Recovered proportion (%)") +
  geom_line(data= binomial_intervals_per_credibility, aes(x=credibility_level,y=min),linetype="dashed") + geom_line(data= binomial_intervals_per_credibility, aes(x=credibility_level,y=max),linetype="dashed") + theme_bw()

g2 <- g2 + facet_grid(. ~ title,labeller = label_parsed) + theme(strip.text = element_text(size = 14))
g2

ggsave(plot = g2, "~/GABIanalysis/validation/plots/calibration_clock.pdf",
       width = 12, height = 12, units = "cm")


##############################
#  plot the clock rate coverage results
##############################

nr_converged_chains = 0

cut_rate_inference = data.frame(seed=rep(1:100, each = 4), hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F,cutRate= rep(1:4,100))

recovered_per_seed = c()

for (seed in 1:100){
  print("Seed: ")
  print(seed)
  # get inference log
  log_file = paste0("validate_GESTALT.", seed, ".log")
  log_data = parse_beast_tracelog_file(paste0(log_dir, log_file))
  log_data_wo_burnin = remove_burn_ins(log_data, burn_in_fraction = 0.3)
  
  # get true parameter from simulation
  simulation_parameter_file = paste0(simulation_dir, "simcutRates.csv")
  simulation_parameters = read.csv(simulation_parameter_file)
  colnames(simulation_parameters) = c("parameter", "value")
  
  #check that chains has converged
  esses = calc_esses(log_data_wo_burnin, sample_interval = 500)
  esses = esses[!colnames(esses) %in% c("birthRate", "deathRate", "prior")]
  esses = esses[! names(esses) %in% c("samplingProportion", "treeHeight.t.alignment", "treeLength.t.alignment")]
  esses
  
  if (any(esses < 200)){
    problematic_ess = esses[which(esses < 200)]
    print(paste("Ess for log file ", log_file, " has ESS < 200: for parameters", paste(names(problematic_ess), collapse = ",")))
    #next
    
  }else{
    inc(nr_converged_chains)
  }
  
  recovered_per_alpha = c()
  
  for (cut_rate_nr in 1:4){
    
    #check HPD recovery rate per credibility level
    for(alpha in seq(5,95,by=5)) {
      
      cut_rate = paste0("cutRates.",cut_rate_nr )
      
      cut_rate_hpd = hdi(object = log_data_wo_burnin[, cut_rate], credMass = alpha/100)
      
      median_cut_rate = median(log_data_wo_burnin[, cut_rate])
      
      true_cut_rate = simulation_parameters[seed,cut_rate_nr + 1]
      
      recovered = is_parameter_in_hpd(hpd_lower = cut_rate_hpd["lower"], hpd_upper = cut_rate_hpd["upper"],true_parameter = true_cut_rate)
      
      recovered_per_alpha <- c(recovered_per_alpha,recovered)
      
      row_index = which(cut_rate_inference$seed == seed & cut_rate_inference$cutRate == cut_rate_nr)
      
      if(alpha/100 == 0.95) {
        cut_rate_inference[row_index, ] = c(seed, as.numeric(cut_rate_hpd), median_cut_rate, true_cut_rate,recovered, cut_rate_nr)
      }
    }
    
  }
  
  recovered_per_seed <- rbind(recovered_per_seed,recovered_per_alpha)
}
#coverages per insert rate at the 95% level
coverages_per_cut <- c()
for(i in 1:4) {
  coverage <- sum(cut_rate_inference[which(cut_rate_inference$cutRate == i),"recovered"])
  coverages_per_cut <- c(coverages_per_cut,coverage)
}
coverages_per_cut

# get coverages per level 
calibration <- data.frame(credibility=rep(seq(5,95,by=5),1),recovered=colSums(recovered_per_seed))

# construct CI on the binomial distribution
binomial_intervals_per_credibility <- c() 
for(i in seq(0.05,0.95,by=0.05)) {
  binomial_intervals_per_credibility <- rbind(binomial_intervals_per_credibility,qbinom(c(0.025,0.975), 100, i, lower.tail = TRUE, log.p = FALSE))
}
binomial_intervals_per_credibility <- cbind(seq(5,95,by=5),binomial_intervals_per_credibility)
colnames(binomial_intervals_per_credibility) <- c("credibility_level","min","max")
binomial_intervals_per_credibility <- data.frame(binomial_intervals_per_credibility)
calibration$title <- as.character(expression(paste(lambda, " - Target cut rates")))

g3 = ggplot(data=calibration, aes(x=credibility, y=recovered)) + geom_point(size=2) + geom_abline(slope = 1.0, col="darkgreen") + xlab("Credibility level (%)") + ylab("Recovered proportion (%)") +
  geom_line(data= binomial_intervals_per_credibility, aes(x=credibility_level,y=min),linetype="dashed") + geom_line(data= binomial_intervals_per_credibility, aes(x=credibility_level,y=max),linetype="dashed") + theme_bw()

g3 <- g3 + facet_grid(. ~ title,labeller = label_parsed) + theme(strip.text = element_text(size = 14))
g3

ggsave(plot = g3, "~/GABIanalysis/validation/plots/calibration_cuts.pdf",
       width = 12, height = 12, units = "cm")


#getting the bias and RMSE
biases <- c()
RMSEs <- c()
widths <- c()
for (seed in 1:100){
  biases <- c(biases,get_bias_rel(cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "median"],cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "true_value"]))
  RMSEs <- c(RMSEs,get_RMSE_rel(cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "median"],cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "true_value"]))
  widths <- c(widths,get_hpd_width_rel(cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "hpd_lower"],cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "hpd_upper"],cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "true_value"]))
  
}

##############################
#  plot coverage results on the long deletion factors
##############################



nr_converged_chains = 0

cut_rate_inference = data.frame(seed=rep(1:100, each = 2), hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F,cutRate= rep(1:2,100))

recovered_per_seed = c()

for (seed in 1:100){
  print("Seed: ")
  print(seed)
  # get inference log
  log_file = paste0("validate_GESTALT.", seed, ".log")
  log_data = parse_beast_tracelog_file(paste0(log_dir, log_file))
  log_data_wo_burnin = remove_burn_ins(log_data )
  
  # get true parameter from simulation
  simulation_parameter_file = paste0(simulation_dir, "simlongTrimScaling.csv")
  simulation_parameters = read.csv(simulation_parameter_file)
  colnames(simulation_parameters) = c("parameter", "value")
  
  #check that chains has converged
  esses = calc_esses(log_data_wo_burnin, sample_interval = 500)
  esses = esses[!colnames(esses) %in% c("birthRate", "deathRate", "prior")]
  esses = esses[! names(esses) %in% c("samplingProportion", "treeHeight.t.alignment", "treeLength.t.alignment")]
  esses
  
  if (any(esses < 200)){
    problematic_ess = esses[which(esses < 200)]
    print(paste("Ess for log file ", log_file, " has ESS < 200: for parameters", paste(names(problematic_ess), collapse = ",")))
    #next
    
  }else{
    inc(nr_converged_chains)
  }
  
  recovered_per_alpha = c()
  
  for (cut_rate_nr in 1:2){
    
    #check HPD recovery rate per credibility level
    for(alpha in seq(5,95,by=5)) {
      
      cut_rate = paste0("longTrimScaling.",cut_rate_nr )
      
      cut_rate_hpd = hdi(object = log_data_wo_burnin[, cut_rate], credMass = alpha/100)
      
      median_cut_rate = median(log_data_wo_burnin[, cut_rate])
      
      true_cut_rate = simulation_parameters[seed,cut_rate_nr + 1]
      
      recovered = is_parameter_in_hpd(hpd_lower = cut_rate_hpd["lower"], hpd_upper = cut_rate_hpd["upper"],true_parameter = true_cut_rate)
      
      recovered_per_alpha <- c(recovered_per_alpha,recovered)
      
      row_index = which(cut_rate_inference$seed == seed & cut_rate_inference$cutRate == cut_rate_nr)
      
      if(alpha/100 == 0.95) {
        cut_rate_inference[row_index, ] = c(seed, as.numeric(cut_rate_hpd), median_cut_rate, true_cut_rate,recovered, cut_rate_nr)
      }
    }
    
  }
  
  recovered_per_seed <- rbind(recovered_per_seed,recovered_per_alpha)
}
#95 % coverages per factor
coverages_per_cut <- c()

for(i in 1:2) {
  coverage <- sum(cut_rate_inference[which(cut_rate_inference$cutRate == i),"recovered"])
  coverages_per_cut <- c(coverages_per_cut,coverage)
}

coverages_per_cut



#get coverage per credibility level
calibration <- data.frame(credibility=rep(seq(5,95,by=5),1),recovered=colSums(recovered_per_seed))

#construct CI for the binonial distribution at all credibility levels
binomial_intervals_per_credibility <- c() 
for(i in seq(0.05,0.95,by=0.05)) {
  binomial_intervals_per_credibility <- rbind(binomial_intervals_per_credibility,qbinom(c(0.025,0.975), 100, i, lower.tail = TRUE, log.p = FALSE))
}
binomial_intervals_per_credibility <- cbind(seq(5,95,by=5),binomial_intervals_per_credibility)
colnames(binomial_intervals_per_credibility) <- c("credibility_level","min","max")
binomial_intervals_per_credibility <- data.frame(binomial_intervals_per_credibility)
calibration$title <- as.character(expression(paste(gamma, " - Long trim scaling factors")))

g4 = ggplot(data=calibration, aes(x=credibility, y=recovered)) + geom_point(size=2) + geom_abline(slope = 1.0, col="darkgreen") + xlab("Credibility level (%)") + ylab("Recovered proportion (%)") +
  geom_line(data= binomial_intervals_per_credibility, aes(x=credibility_level,y=min),linetype="dashed") + geom_line(data= binomial_intervals_per_credibility, aes(x=credibility_level,y=max),linetype="dashed") + theme_bw()

g4 <- g4 + facet_grid(. ~ title,labeller = label_parsed) + theme(strip.text = element_text(size = 14))
g4

ggsave(plot = g4, "~/GABIanalysis/validation/plots/calibration_trims.pdf",
       width = 12, height = 12, units = "cm")

#getting the bias and RMSE
biases <- c()
RMSEs <- c()
widths <- c()

for (seed in 1:100){
  biases <- c(biases,get_bias_rel(cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "median"],cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "true_value"]))
  RMSEs <- c(RMSEs,get_RMSE_rel(cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "median"],cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "true_value"]))
  widths <- c(widths,get_hpd_width_rel(cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "hpd_lower"],cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "hpd_upper"],cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 1, "true_value"]))
  
}


biases <- c()
RMSEs <- c()
widths <- c()

for (seed in 1:100){
  biases <- c(biases,get_bias_rel(cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 2, "median"],cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 2, "true_value"]))
  RMSEs <- c(RMSEs,get_RMSE_rel(cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 2, "median"],cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 2, "true_value"]))
  widths <- c(widths,get_hpd_width_rel(cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 2, "hpd_lower"],cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 2, "hpd_upper"],cut_rate_inference[cut_rate_inference$seed == seed & cut_rate_inference$cutRate == 2, "true_value"]))
  
}


############################
# plot all 4 parameters together
############################

# Generate dummy dataset
foo <- data.frame(x1 = 1, x2 = 2, y1 = 1, y2 = 2,
                  text = paste(letters[1:3], letters[1:3], collapse = "\n"))

# Plot rectangle with text
title <- ggplot(foo) +
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
            color = "black", size = 0.5, fill = "grey") +
  geom_text(aes(x = x1 + (x2 - x1) / 2, y = y1 + (y2 - y1) / 2,
                label = "Inference - 95 percent credibility level"),
            size = 6) +
  theme_classic() +
  theme(axis.line  = element_blank(),
        axis.ticks = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank())

all_plots <- cowplot::plot_grid(g1,g2,g3,g4,legend,nrow=1)
#all_plots_with_title <- cowplot::plot_grid(title,all_plots,nrow=2,rel_heights = c(0.2,1),align = "h")
ggsave(plot = all_plots, "~/GABIanalysis/validation/plots/all_calibrations_aligned.pdf",
       width = 52, height = 12, units = "cm")

