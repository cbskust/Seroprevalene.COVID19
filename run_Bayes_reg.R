setwd("~/Dropbox/PAPER/SeroSIR/GitHub")



################################################################################
# This script is to accompany the manuscript 
# Analysis of Bayesian regression weekly Epidemic Data: Study of Louisville, Kentucky.
# 
# Specify the appropriate data paths and source to run.
#
# Written By: Boseung Choi
# Edited By: GR & Boseung Choi
# Last Edited: Feb 23, 2022
#
################################################################################


require(dplyr)
require(tibble)
require(MASS)
require(rstan)
require(rstanarm)
require("bayesplot")
require("ggplot2")


### Select one of area 
  site_names=c('MSD01', 'MSD02','MSD03','MSD04','MSD05'); zone = 'msd_agg' #'aggregate'
  #site_names=c('MSD01'); zone = 'msd_01' #zone_1' 
  #site_names=c('MSD02'); zone = 'msd_02' #zone_2'
  #site_names=c('MSD03','MSD04','MSD05'); zone = 'msd_0345' #zone_345'


    
### Save paths
  path.nb.reg <- paste0('results/MCMC_nb.reg_',zone,'.csv')
  path.lm.reg <- paste0('results/MCMC_lm.reg_',zone,'.csv')
  

    
### Setting the first and last date of analysis; depending on Wastewater data
min_date = as.Date("2020-08-17")
max_date = as.Date("2021-02-22")

### Reading wastewater concentration data 
dat_water <-read.csv(file =  paste0("data/WW_",zone,".csv"))

##### Reading population data; calculating the total population   
dat_pop <-read.csv(file =  paste0("data/pop_",zone,".csv")) %>%
  as_tibble() %>%
  transmute(date = as.Date(date),
            site  = as.character(site), 
            pop = as.numeric(pop) ) %>%
  filter(date<=max_date & date >=min_date)
dat_pop<-dat_pop[order(as.numeric(dat_pop$date)),]
site_pop = sum(unique(dat_pop$pop))       

##### prevalence estimation result 
dat_pred.prev <- read.csv(paste0('results/prevalence_estimate_',zone,'.csv')) %>%
  as_tibble() %>%
  transmute(date = as.Date(date),
            prev = as.numeric(med),
            pop = site_pop,
            prev.cnt = round(as.numeric(med) * site_pop)) %>%
  filter(date<=max_date & date >=min_date) %>%
  group_by(date) %>%
  summarize(prev = mean(prev), pop =pop, prev.cnt =round(mean(prev.cnt))) %>%
  distinct() 

  ## converting weekly data 
  time_point=as.numeric(dat_pred.prev$date-min(dat_pred.prev$date))
  dat_pred.prev=add_column(dat_pred.prev,time_point)
  week = round(dat_pred.prev$time_point/7)
  dat_pred.prev = add_column(dat_pred.prev, week)
  
  dat_pred.prev <- dat_pred.prev %>% 
  as_tibble() %>%
  group_by(week) %>%
  summarize(prev = mean(prev), pop =pop, prev.cnt =round(mean(prev.cnt))) %>% 
  distinct() 

###### Merge data for regression analysis 

dat_reg <- merge( x= dat_pred.prev, y = dat_water, by = "week", all = T )
dat_reg <-dat_reg %>% 
  filter(!is.na(avg_n1))


################################################################################
### Run Stan for NB regression ###
nb.fit.Bayes <-stan_glm.nb(prev.cnt ~ avg_n1, # symbolic description of the model  
                  link = "log",           # link function    
                  data= dat_reg,          # name of data 
                  warmup =2000,           # number of warmup iterations per chain
                  iter = 4000,            # total number of iterations per chain
                  chains = 2,             # number of Markov chains
                  refresh = 1000,         # show progress every 'refresh' iterations
                  prior_aux = exponential(.01) # assigning non-inforamtive prior to dispersion         
                  seed=12345
                  )

# output fit summary status
print(summary(nb.fit.Bayes, digits =6, probs=c(0.025, 0.5, 0.975) ))
         

# output  posterior traceplot & histogram 
# mcmc_trace(nb.fit.Bayes, pars = c("avg_n1", "reciprocal_dispersion"), 
#           facet_args = list(ncol = 1, strip.position = "left"))

# mcmc_hist(nb.fit.Bayes, pars = c("avg_n1", "reciprocal_dispersion"))


####### Fitting NB regession 
nb.fit <- glm.nb(prev.cnt ~ avg_n1 , link = log, data= dat_reg, )
print(summary(nb.fit))
#print(prior_summary(nb.fit.Bayes))

################################################################################



####### Fitting linear regession 
lm.fit.Bayes <-stan_glm(prev ~ avg_n1,         # symbolic description of the model
                        family = gaussian(),   # liner regession model 
                        data= dat_reg,        # name of data 
                        warmup =2000,          # number of warmup iterations per chain
                        iter = 4000,           # total number of iterations per chain
                        chains = 2,            # number of Markov chains
                        refresh = 1000,        # show progress every 'refresh' iterations
                        seed=12345
)

# output fit summary status
print(summary(lm.fit.Bayes, digits = 9, probs=c(0.025, 0.5, 0.975)))

# output  posterior traceplot & histogram 
# mcmc_trace(lm.fit.Bayes, pars = c("avg_n1"), 
#           facet_args = list(ncol = 1, strip.position = "left"))

# mcmc_hist(lm.fit.Bayes, pars = c("avg_n1"))


####### Fitting linear regession 

lm.fit <- lm(prev ~ avg_n1 , data= dat_reg)
print(summary(lm.fit))

#print(prior_summary(lm.fit.Bayes))


############################# Save csvs
write.csv(as.data.frame(nb.fit.Bayes), path.nb.reg , row.names = F)  #MCMC posteior samples of parameter 
write.csv(as.data.frame(lm.fit.Bayes), path.lm.reg , row.names = F)  #MCMC posteior samples of parameter 

