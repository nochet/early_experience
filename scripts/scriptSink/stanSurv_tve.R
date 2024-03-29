# Load packages

library(tidyverse)
library(survival)
library(rstanarm)
library(loo)
library(brms)
options(mc.cores = parallel::detectCores())
library(cowplot)
library(bayesplot)
#theme_set(bayesplot::theme_default())
theme_set(bayesplot::theme_default(base_family = "sans"))
ggplot2::theme_set(bayesplot::theme_default())

eelife <- read.csv("../processed/eeAlldat.csv") %>%
  select(-c(X,larv_adult))


# Bayesian survival models

start_time <- Sys.time()

### M-splines is the default in stan_surv() 


# Investigate association between hazard of death and two longitudinal treatments: larval and adult
mod0 <- stan_surv(formula = Surv(age, status) ~ larvalTreat + 
                  adultTreat + sex, 
                  data    = eelife,
                  chains  = 4, 
                  cores   = 8, 
                  seed    = 29,
                  iter    = 5e4)
print(mod0, digits = 3)
summary(mod0)
mcmc_areas(as.matrix(mod0), prob_outer = .999)
saveRDS(mod0, file="../processed/bhaz_mod0.rds")
mod0 <- readRDS("../processed/bhaz_mod0.rds")

log_lik0 <- log_lik(mod0) 


reff0 <- as.array(mod0)
r_eff <- loo::relative_eff(exp(reff0))
loo_0 <- loo.array(reff0, r_eff = r_eff, 
             cores=getOption("mc.cores",2))
print(loo_0)

saveRDS(loo_0, file = "../processed/loo_0.rds")

# Estimate a time-varying effect (coefficient) to be estimated for each variable
mod1 <- stan_surv(formula = Surv(age, status) ~ tve(larvalTreat) + 
                  tve(adultTreat) + tve(sex), 
                  data    = eelife,
                  chains  = 4, 
                  cores   = 8, 
                  seed    = 4147,
                  iter    = 5e4)
print(mod1, digits = 3)
summary(mod1)
save(mod1, file="../processed/bhaz_mod1.Rda")

reff1 <- as.array(mod1)
r_eff1 <- loo::relative_eff(exp(reff1))
loo_1 <- loo.array(reff1, r_eff = r_eff1, 
                   cores=getOption("mc.cores",4))
print(loo_1)

saveRDS(loo_1, file = "../processed/loo_1.rds")

# Test different parametric hazards
# i.e. fit several models, each with a different baseline hazard

mod1_exp      <- update(mod1, basehaz = "exp") 
mod1_weibull  <- update(mod1, basehaz = "weibull") 
mod1_gompertz <- update(mod1, basehaz = "gompertz") 
mod1_bspline  <- update(mod1, basehaz = "bs") 
mod1_mspline1 <- update(mod1, basehaz = "ms") 
mod1_mspline2 <- update(mod1, basehaz = "ms", 
                        basehaz_ops = list(df = 10))

save(mod1,mod1_exp,mod1_weibull,mod1_gompertz,mod1_bspline,
        mod1_mspline1,mod1_mspline2, file="../processed/bhaz_Compare.Rda")

reff2 <- as.array(mod1_bspline)
r_eff2 <- loo::relative_eff(exp(reff2))
loo_2 <- loo.array(reff2, r_eff = r_eff2, 
                   cores=getOption("mc.cores",4),
                   moment_match=FALSE,
                   k_threshold=0.7)
print(loo_2)

saveRDS(loo_2, file = "../processed/loo_2.rds")


load("../processed/bhaz_Compare.Rda")

loop <- loo(mod1)
loop.ex <- loo(mod1_exp)
loop.wb <- loo(mod1_weibull)
loop.gp <- loo(mod1_gompertz)
loop.bs <- loo(mod1_bspline,k_threshold = 0.7) # advisory from warning msg (see notes)
loop.ms <- loo(mod1_mspline1)
loop.mp <- loo(mod1_mspline2)

saveRDS(loop,loop_ex,loop_wb,loop_gp,loop_bs,
        loop_ms,loop_mp, file="../processed/loop.rds")

# Compare fits using `loo` for `stansurv` objects
loo.comp <- loo_compare(
  loo(mod1),
  loo(mod1_exp),
  loo(mod1_weibull),
  loo(mod1_gompertz),
  loo(mod1_bspline),
  loo(mod1_mspline1),
  loo(mod1_mspline2)) # better of all

saveRDS(loo.comp, file="../processed/loop_compare.rds")



