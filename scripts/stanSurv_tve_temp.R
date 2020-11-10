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


load("../processed/bhaz_Compare.Rda")

loop.mp <- loo(mod1_mspline2)

saveRDS(loop.mp, file="../processed/loop_mspline2.rds")






