library(tidyverse)
library(survival)
library(rstanarm)


eelife <- read.csv("../processed/eeAlldat.csv") %>%
  select(-c(X,larv_adult))

load("../processed/BaseHaz_Compare.Rda")


# Compare fits using `loo` for `stansurv` objects

looc <- loo_compare(loo(mod1_exp),
            loo(mod1_weibull),
            loo(mod1_gompertz),
            loo(mod1_bspline),
            loo(mod1_mspline1),
            loo(mod1_mspline2))

saveRDS(looc, file="../processed/loo_compare_bhaz.rds")

