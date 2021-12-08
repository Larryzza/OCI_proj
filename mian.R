rm(list = ls())

library(tidyverse)
library(lazyeval)
library(rstan) 
library(shinystan) 
library(purrr)
library(data.table)
options(mc.cores=parallel::detectCores())

OCI_data <- read_csv("OCI_data.csv")
#View(OCI_data)

source("functions.R")
source("set_pars.R")
source("plot.R")

# construct input for stan

indiv_data <- OCI_data %>% 
  mutate(id_clean = rep(1:length(unique(Person.ID)), 
                        table(Person.ID)),
         id = Person.ID,
         t = Date.Index,
         y = ifelse(CT.Mean < global_pars['lod'], CT.Mean, global_pars['lod']),
         special = 1) %>%
  mutate(ref = ifelse(y < global_pars['lod'],1,0),
         ref1 = lag(ref,2),
         ref2 = lead(ref,2),
         ref3 = ref + ref1 + ref2) %>% 
  filter(ref3 > 0) %>% 
  select(id, id_clean, t, y, special)

n_indiv <- length(unique(indiv_data$id))

special <- indiv_data %>%
  group_by(id) %>%
  slice(1) %>%
  select(id, special) %>%
  arrange(id) %>%
  pull(special)

# Useful dataframe for mapping official ID to Stan ID:
id_map <- indiv_data %>% 
  group_by(id) %>%
  summarise(id_clean=first(id_clean)) %>% 
  select(id, id_clean) %>%
  mutate(id_clean=as.character(id_clean))

# Useful dataframe for mapping PersonID to the 'special' entries
special_map <- indiv_data %>% 
  group_by(id) %>%
  summarise(special=first(special))

# set priors
for(i in 1:length(run_pars_list)){
  run_pars <- run_pars_list[[i]]
  
  prior_pars <- list(
    special=special,
    tpsd=run_pars$tpsd,
    dpmin=run_pars$dpmin,
    dpmean_prior=run_pars$dpmean_prior,
    dpsd_prior=run_pars$dpsd_prior,
    wpmin=run_pars$wpmin,
    wpmax=run_pars$wpmax,
    wpmean_prior=run_pars$wpmean_prior,
    wpsd_prior=run_pars$wpsd_prior,
    wrmin=run_pars$wrmin,
    wrmax=run_pars$wrmax,
    wrmean_prior=run_pars$wrmean_prior,
    wrsd_prior=run_pars$wrsd_prior,
    sigma_max=run_pars$sigma_max,
    sigma_prior_scale=run_pars$sigma_prior_scale,
    lambda=run_pars$lambda,
    fpmean=run_pars$fpmean		# so that 90% of mass is <1 and 99% is <2
  )	
  
  ct_model <- stan_model("2fit_posteriors_special.stan") 
  #ct_model <- stan_model("fit_posteriors_special.stan") 
  
  fit_startq <- Sys.time()
  ct_fit <- sampling(ct_model, 
                     data=list(
                       N=nrow(indiv_data), 
                       n_id=length(unique(indiv_data$id_clean)),
                       lod=global_pars[["lod"]], 
                       id=indiv_data$id_clean,
                       t=indiv_data$t, 
                       y=indiv_data$y, 
                       tpsd=as.list(prior_pars)$tpsd,
                       dpmin=as.list(prior_pars)$dpmin,
                       dpmean_prior=as.list(prior_pars)$dpmean_prior,
                       dpsd_prior=as.list(prior_pars)$dpsd_prior,
                       wpmin=as.list(prior_pars)$wpmin,
                       wpmax=as.list(prior_pars)$wpmax,
                       wpmean_prior=as.list(prior_pars)$wpmean_prior,
                       wpsd_prior=as.list(prior_pars)$wpsd_prior,
                       wrmin=as.list(prior_pars)$wrmin,
                       wrmax=as.list(prior_pars)$wrmax,
                       wrmean_prior=as.list(prior_pars)$wrmean_prior,
                       wrsd_prior=as.list(prior_pars)$wrsd_prior,
                       sigma_max=as.list(prior_pars)$sigma_max,
                       sigma_prior_scale=as.list(prior_pars)$sigma_prior_scale,
                       lambda=as.list(prior_pars)$lambda,
                       fpmean=as.list(prior_pars)$fpmean), 
                     iter=10000, chains=4,
                     control = list(adapt_delta=0.90))
  # , control = list(adapt_delta=0.85)
  # , control = list(adapt_delta=0.99)
  # control = list(adapt_delta=0.95, max_treedepth=15)
  fit_endq <- Sys.time()
  print(paste0("Fit time: ",difftime(fit_endq, fit_startq, units="min")," mins"))
  
  # launch_shinystan_nonblocking(ct_fit)
  
  params <- rstan::extract(ct_fit)
  saveRDS(params, paste0(i,"_params.rds"))
}

# prepare results and output
for(i in 1:length(run_pars_list)){
  params <- readRDS(paste0(i, "_params.rds"))
  indiv_params_df <- make_indiv_params_df(params, c("tp","dp","wp","wr"), n_indiv) %>% 
    rename(id_clean=id) %>% 
    left_join(id_map, by="id_clean") %>%
    left_join(special_map, by="id")
  
  shared_params_df <- make_shared_params_df(params, c("dpmeanB","wpmeanB","wrmeanB","dpsd","wpsd","wrsd")) %>% 
    mutate(dpmeanB_trans=truncnormmean(dpmeanB,dpsd,run_pars$dpmin,global_pars[["lod"]])) %>%
    mutate(wpmeanB_trans=truncnormmean(wpmeanB,wpsd,run_pars$wpmin,run_pars$wpmax)) %>%
    mutate(wrmeanB_trans=truncnormmean(wrmeanB,wrsd,run_pars$wrmin,run_pars$wrmax))
  
  params_df <- indiv_params_df %>% 
    left_join(shared_params_df, by="iteration") %>% 
    select(-iteration) 
  
  fig_ct_fit <- plot_ct_fit(params_df, global_pars, indiv_data, 
                            ctalpha=0.01, xlim=c(-20,60), ntraces=100)
  
  meanvalsindiv <- params_df %>% 
    mutate(infdur=wp+wr) %>%
    group_by(id) %>% 
    summarise(tp=mean(tp), dp=mean(dp), wp=mean(wp), wr=mean(wr), 
              infdur=mean(infdur), special=first(special))
  
  # Summarize the posterior distributions for easy conversion into a table: 
  dist_summary <- shared_params_df %>% 
    select(-dpmeanB, -wpmeanB, -wrmeanB) %>%
    rename(dpmeanB=dpmeanB_trans,
           wpmeanB=wpmeanB_trans,
           wrmeanB=wrmeanB_trans) %>% 
    summarise(
      peak.ct.Variant_mean=mean(global_pars[["lod"]]-dpmeanB),
      peak.ct.Variant_lwr=quantile(global_pars[["lod"]]-dpmeanB,0.025),
      peak.ct.Variant_upr=quantile(global_pars[["lod"]]-dpmeanB,0.975),
      peak.geml.Variant_mean=mean(convert_Ct_logGEML(global_pars[["lod"]]-dpmeanB)),
      peak.geml.Variant_lwr=quantile(convert_Ct_logGEML(global_pars[["lod"]]-dpmeanB),0.025),
      peak.geml.Variant_upr=quantile(convert_Ct_logGEML(global_pars[["lod"]]-dpmeanB),0.975),
      proliferation.time.Variant_mean=mean(wpmeanB),
      proliferation.time.Variant_lwr=quantile(wpmeanB,0.025),
      proliferation.time.Variant_upr=quantile(wpmeanB,0.975),
      clearance.time.Variant_mean=mean(wrmeanB),
      clearance.time.Variant_lwr=quantile(wrmeanB,0.025),
      clearance.time.Variant_upr=quantile(wrmeanB,0.975),
      total.duration.Variant_mean=mean(wpmeanB+wrmeanB),
      total.duration.Variant_lwr=quantile(wpmeanB+wrmeanB,0.025),
      total.duration.Variant_upr=quantile(wpmeanB+wrmeanB,0.975)
    ) %>%
    pivot_longer(everything()) %>%
    separate(name, c("parameter", "statistic"), sep="_") %>%
    pivot_wider(names_from=statistic, values_from=value) %>%
    arrange(parameter)
  
  dist_summary$mean <- round(dist_summary$mean, 2)
  dist_summary$lwr <- round(dist_summary$lwr, 2)
  dist_summary$upr <- round(dist_summary$upr, 2)
  dist_summary <- dist_summary[-c(3,5),]
  
  output <- list(plot=fig_ct_fit,
                 summary=dist_summary)
  saveRDS(output, paste0(i, "_output.rds"))
}

# show results
# Uninformative priors
df1 <- readRDS("1_output.rds")
df1[["plot"]]
df1[["summary"]]

# Informative priors 
df2 <- readRDS("2_output.rds")
df2[["plot"]]
df2[["summary"]]

# Biased set of priors
df3 <- readRDS("3_output.rds")
df3[["plot"]]
df3[["summary"]]


