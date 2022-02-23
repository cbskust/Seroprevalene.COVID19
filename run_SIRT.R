#setwd("/Users/cbskust/Dropbox/PAPER/SeroSIR/GitHub")
setwd("~/Documents/SeroSIR/GitHub")


################################################################################
# This script is to accompany the manuscript 
# Analysis of Individual-level Epidemic Data: Study of Louisville, Kentucky.
# 
# Specify the appropriate data paths and source to run.
#
# Written By: Chance Alvarado
# Edited By: GR & Boseung Choi
# Last Edited: Feb 17, 2022
#
################################################################################

### Needed libraries ###

library(rstan)                            # Stan MCMC packaged
require("bayesplot")
require("ggplot2")

set.seed(1098)                           # Fixed seed for repoducibility  
################################################################################

### Specify data paths ###

# Specify path to data: msd_agg, msd_01, msd_02, or msd_0345
zone <- 'msd_agg'
#zone <- 'msd_01'
#zone <- 'msd_02'
#zone <- 'msd_0345'


# Start and end of relevant data
date.start <- 0
date.end <- 300
# Plot beyond last data pt
h.T <- 50

################################################################################
# input data path 
data.path <- paste0("data/",zone,'.csv')

# Save paths
path.I <- paste0('results/prevalence_estimate_',zone,'.csv')
path.parameters <- paste0('results/MCMC_parameters_',zone,'.csv')


# Stan model path
stan.path <- 'stan_SIRT_gammaSSa.stan' 

# Custom ranking function & helper functions
source('helper_functions.R')    


################################################################################

### Read in and format data ###

# Read in current data
data.current_zone.raw <- read.csv(data.path)

# Group by date, not observation
data.current_zone.formatted <- aggregate(antibody_result ~ time_point,
                                         data.current_zone.raw, sum)

# Rename columns
colnames(data.current_zone.formatted) <- c('time_point',
                                           'num_antibody_positive')

# Add column of test date
data.current_zone.formatted$date <- sort(unique(
  as.Date(data.current_zone.raw$test_date)))

origin.date <- data.current_zone.formatted$date[1]

# Add column of total tests by day
tests.by.day <- aggregate(antibody_result ~ time_point,
                          data.current_zone.raw, length)$antibody_result
data.current_zone.formatted$num_antibody_total <- tests.by.day

################################################################################

### Format data and helper variables for Stan ###

# Subset data for specified wave
data.current_zone.subset <- data.current_zone.formatted[
  (data.current_zone.formatted$time_point >= date.start) &
    (data.current_zone.formatted$time_point <= date.end),]

# Number of data points
k <- length(data.current_zone.subset$time_point)                

# Initial date
t.origin <- data.current_zone.subset$time_point[1]             

# Start from day one
data.current_zone.subset$time_point <- data.current_zone.subset$time_point - 
  t.origin + 1

# Data as list for Stan
data.SIRT <- list(k=k,
                  t0=0.5,
                  ti=data.current_zone.subset$time_point,
                  xi=data.current_zone.subset$num_antibody_positive,
                  ni=data.current_zone.subset$num_antibody_total
)

################################################################################

### Run Stan code on data ###
tempdir()
fit <- stan(
  file=stan.path,           # Stan program
  data=data.SIRT,           # named list of data
  chains=2,                 # number of Markov chains
  warmup=1000,              # number of warmup iterations per chain
  iter=3000,                # total number of iterations per chain
  cores=3,                  # number of cores
  refresh=1000,            # show progress every 'refresh' iterations
  control=list(adapt_delta=.9)
)

################################################################################

### output to device 
#print(summary(fit))                # output fit summary status

                                   # output  posterior densities 
#print(plot(fit, show_density=T,pars=c('beta', 'gamma', 'delta'), include=TRUE,fill_color="green")) 

                                   # output  posterior traceplot
#print(traceplot(fit))

### save final samples to file

out.fit = as.data.frame(cbind(extract(fit)$beta, extract(fit)$gamma, extract(fit)$delta,
                              extract(fit)$rho, extract(fit)$eps, extract(fit)$psi,
                              extract(fit)$lp__ ))
colnames(out.fit) = c('beta', 'gamma', 'delta', 'rho', 'eps', 'psi',"lp")


################################################################################
### Solve trajectories ###

# Percent positive column
data.current_zone.subset$percent_positive <- data.current_zone.subset$num_antibody_positive /
  data.current_zone.subset$num_antibody_total

# Rank parameters
parameters <- c('beta', 'gamma', 'delta', 'rho', 'eps', 'psi')
parameters.frame <- as.data.frame(extract(fit)[parameters])
ranked <- rank_by_distance2(parameters.frame, metric='median')
parameters.ranked <- ranked[[1]]
row <- ranked[[2]]

# Get 95% credible parameter draws
ci = 0.90
parameters.ranked.ci <- parameters.ranked[seq(1, as.integer(ci * nrow(parameters.ranked), 1)),]

# Sample trajectories
sample.size <- 500
parameter.draws <- parameters.ranked.ci[sample(nrow(parameters.ranked.ci), sample.size), ]

# Iteration to solve to 
time.point <- data.current_zone.subset$time_point[nrow(data.current_zone.subset)]+h.T
time.step <- 0.5

# Create matrices of trajectories
trajectory.matrix.I <- matrix(nrow=sample.size,
                              ncol=time.point / time.step)

# Index through all parameter draws
for (i in 1:sample.size){
  # Get current row
  row <- parameter.draws[i, ]
  
  # get parameters
  beta <- row$beta
  gamma <- row$gamma
  delta <- row$delta
  rho <- row$rho
  eps <- row$eps
  psi <- row$psi
  
  # Solve DE
  sol <- euler(t=time.point,
               dt=time.step,
               fun=ode.sirt,
               ic=c(1 - rho - eps - psi, rho, eps, psi))
  
  
  # Save I and T only
  trajectory.matrix.I[i, ] <- sol[,2]  # Prevalace 
}

bep=0.0
# Get confidence bounds
top.I <- apply(trajectory.matrix.I, 2, max)*(1+bep)
bottom.I <- apply(trajectory.matrix.I, 2, min)*(1-bep)
median.I <-apply(trajectory.matrix.I, 2, median)


# Format trajectory data
data.I <- data.frame(time=seq(from=0, to=time.point-time.step, by=time.step),
                     top=top.I,
                     med=median.I,
                     bottom=bottom.I)

# Add psuedo median line
data.I$med <- ((data.I$top - data.I$bottom) / 2) + data.I$bottom

# Add date column
data.I$date <- data.I$time + as.Date(data.current_zone.subset$date[1])


############################# Save csvs

write.csv(data.I, path.I,row.names = F)  #MCMC posteior samples of parameter 
write.csv(out.fit, path.parameters, row.names=F,col.names=NULL) # estimated prevalence percent


########################## end of script code ##################################### 
