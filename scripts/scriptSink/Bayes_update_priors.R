
# mod1, update priors

library(tidyverse)
library(rstanarm)
library(loo)
library(brms)
library(cowplot)
library(bayesplot)

load("../processed/bhaz_Compare.Rda")


# declare desired parameters for stan
CHAINS <- 4 
CORES <- 2 
ITER <- 5e3 
SEED <- 39


prior_summary(mod1)

mod1_update_prior <- update(mod1, 
                            prior_intercept = normal(0, 1),
                            prior = normal(0, 0.5),
                            prior_PD = TRUE,
                            chains = CHAINS,
                            cores = CORES,
                            iter = ITER,
                            seed = SEED)

mod1.fit <- update(mod1_update_prior,
                   prior_PD = FALSE,
                   chains = CHAINS,
                   cores = CORES,
                   iter = ITER,
                   seed = SEED)

print(mod1.fit, digits = 3)

saveRDS(mod1_update_prior,mod1.fit, file="../processed/mod1_updated.rds")

# piecewise constant model (with df=5)
mod1.fit.pw5 <- update(mod1.fit,
                   basehaz = "ms",
                   basehaz_ops = list(degree = 0, df = 5),
                   prior_PD = FALSE,
                   chains = CHAINS,
                   cores = CORES,
                   iter = ITER,
                   seed = SEED)

# mod1.fit_ex      <- update(mod1.fit, 
#                            basehaz = "exp",
#                            prior_PD = FALSE,
#                            chains = CHAINS,
#                            cores = CORES,
#                            iter = ITER,
#                            seed = SEED) 

mod1.fit_weib  <- update(mod1.fit, 
                         basehaz = "weibull",
                         prior_PD = FALSE,
                         chains = CHAINS,
                         cores = CORES,
                         iter = ITER,
                         seed = SEED) 
mod1.fit_gomp <- update(mod1.fit, 
                        basehaz = "gompertz",
                        prior_PD = FALSE,
                        chains = CHAINS,
                        cores = CORES,
                        iter = ITER,
                        seed = SEED) 

# mod1.fit_bspl  <- update(mod1.fit, 
#                          basehaz = "bs",
#                          prior_PD = FALSE,
#                          chains = CHAINS,
#                          cores = CORES,
#                          iter = ITER,
#                          seed = SEED) 

mod1.fit_mspl1 <- update(mod1.fit, 
                         basehaz = "ms",
                         prior_PD = FALSE,
                         chains = CHAINS,
                         cores = CORES,
                         iter = ITER,
                         seed = SEED) 

mod1.fit_mspl2 <- update(mod1.fit, 
                         basehaz = "ms", 
                         basehaz_ops = list(df = 10),
                         prior_PD = FALSE,
                         chains = CHAINS,
                         cores = CORES,
                         iter = ITER,
                         seed = SEED)

# piecewise constant model (with df=10)
mod1.fit_pw5 <- update(mod1.fit, 
                       basehaz = "ms",
                       prior_PD = FALSE,
                       chains = CHAINS,
                       cores = CORES,
                       iter = ITER,
                       seed = SEED) 

mod1.fit.pw10 <- update(mod1.fit, 
                        basehaz = "ms",
                        basehaz_ops = list(degree = 0, df = 10),
                        prior_PD = FALSE,
                        chains = CHAINS,
                        cores = CORES,
                        iter = ITER,
                        seed = SEED)

saveRDS(mod1.fit.pw5,mod1.fit_weib,mod1.fit_gomp,mod1.fit_mspl1,
        mod1.fit_mspl2,mod1.fit_pw5,mod1.fit.pw10,
        file="../processed/mod1_updated_dbns.rds")

mod1_fit.pw <- list(
  "basic_updated" = mod1_update_prior,
  "basic_fit" = mod1.fit,
  "Gompertz" = mod1.fit_gomp,
  "Weibull" = mod1.fit_weib,
  "MS1" = mod1.fit_mspl1,
  "MS2" = mod1.fit_mspl2,
  "PW (df = 5)" = mod1.fit.pw5,
  "PW (df = 10)" = mod1.fit.pw10
)

#plots <- map(mod1_fit.pw, plot) 

pw_up <- bayesplot_grid(
  plots = plots,
  ylim = c(0, 0.5),
  titles = names(mod1_fits),
  grid_args = list(ncol = 3))

ggplot2::ggsave(pw_up,filename="../plots/piecewise_mod1_pupdate.pdf",height = 6, width=10)

# LOOs 
loos <- map(mod1_fit.pw, loo, cores=getOption("mc.cores",4)) 
lc_pw <- loo_compare(loos, cores=getOption("mc.cores",4))
saveRDS(loos,lc_pw, file="../processed/loos_pudate.rds")

# Weights
best_pw <- stanreg_list(mod1_fit.pw)
loo_wt_pw <- loo_model_weights(best_pw,cores=getOption("mc.cores",4))
saveRDS(loo_wt_pw, file="../processed/loo_wt_pupdate.rds")

