
# Basehaz = piecewise
###R CMD BATCH --quiet --no-restore --no-save Bayes_Piecewise.R &###

library(tidyverse)
library(rstanarm)
library(loo)
library(brms)
#options(mc.cores = parallel::detectCores())
library(cowplot)
library(bayesplot)

load("../processed/bhaz_Compare.Rda")

eelife <- read.csv("../processed/eeAlldat.csv") %>%
  select(-c(X,larv_adult))


# declare desired parameters for rstan
CHAINS <- 4 
CORES <- 4 
ITER <- 1e4 
SEED <- 42


# piecewise constant model (with df=5)
mod1.pw5 <- update(mod1,
                   basehaz = "ms",
                   basehaz_ops = list(degree = 0, df = 5),
                   chains = CHAINS,
                   cores = CORES,
                   iter = ITER,
                   seed = SEED)

# piecewise constant model (with df=10)
mod1.pw10 <- update(mod1, 
                    basehaz = "ms",
                    basehaz_ops = list(degree = 0, df = 10),
                    chains = CHAINS,
                    cores = CORES,
                    iter = ITER,
                    seed = SEED)

saveRDS(mod1.pw5,mod1.pw10, file="../processed/mod1_pecewise.rds")

mod1_fits <- list(
  "constant" = mod1_exp,
  "Gompertz" = mod1_gompertz,
  "Weibull" = mod1_weibull,
  "MS1" = mod1_mspline1,
  "MS2" = mod1_mspline2,
  "PW (df = 5)" = mod1.pw5,
  "PW (df = 10)" = mod1.pw10
)

# plots <- map(mod1_fits, plot) 

pw <- bayesplot_grid(
                    plots = plots,
                    ylim = c(0, 0.5),
                    titles = names(mod1_fits),
                    grid_args = list(ncol = 3))

ggplot2::ggsave(pw,filename="../plots/piecewise_mod1.pdf",height = 6, width=10)

# LOOs
loos <- map(mod1_fits, loo, cores=getOption("mc.cores",4)) 
lc_pw <- loo_compare(loos, cores=getOption("mc.cores",4))
saveRDS(loos,lc_pw, file="../processed/loos_pw.rds")

# Weights
best_pw <- stanreg_list(mod1_fits)
loo_wt_pw <- loo_model_weights(best_pw,cores=getOption("mc.cores",4))
saveRDS(loo_wt_pw, file="../processed/loo_wt_pw.rds")
