# Load packages

library(tidyverse)
library(survival)
library(rstanarm)
library(cowplot)
library(bayesplot)
theme_set(bayesplot::theme_default())

eelife <- read.csv("../processed/eeAlldat.csv") %>%
  select(-c(X,larv_adult))


# Bayesian survival models

start_time <- Sys.time()

### M-splines is the default in stan_surv() 

# Baseline hazard function - assuming covariate effect of zero
mod1 <- stan_surv(formula = Surv(age, status) ~ tve(larvalTreat) + 
                    tve(adultTreat) + tve(sex), 
                  data    = eelife,
                  chains  = 8, 
                  cores   = 4, 
                  seed    = 4147,
                  iter    = 5e4)
print(mod1, digits = 3)
summary(mod1)

save(mod1, file="../processed/bhaz_mod1.Rda")

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


# Compare fits using `loo` for `stansurv` objects
loo.comp <- loo_compare(loo(mod1_exp),
            loo(mod1_weibull),
            loo(mod1_gompertz),
            loo(mod1_bspline),
            loo(mod1_mspline1),
            loo(mod1_mspline2))

write.table(loo.comp, file="../processed/loo_compare_mod1.rds", 
            sep = "\t",row.names = FALSE)


# Plot the baseline hazards with 95% posterior uncertainty limits
# load("../processed/BaseHaz.Rda")

plotfun <- function(model, title) {
  plot(model, plotfun = "basehaz") +              
    coord_cartesian(ylim = c(0,0.4)) +            
    labs(title = title) +                        
    theme(plot.title = element_text(hjust = 0.5,size = 10),
          axis.title = element_text(size = 9)) + 
    theme_set(theme_cowplot()) +
    theme_half_open() 
}
p_exp      <- plotfun(mod1_exp,      title = "Exponential")
p_weibull  <- plotfun(mod1_weibull,  title = "Weibull")
p_gompertz <- plotfun(mod1_gompertz, title = "Gompertz")
p_bspline  <- plotfun(mod1_bspline,  title = "B-splines with df = 5")
p_mspline1 <- plotfun(mod1_mspline1, title = "M-splines with df = 5")
p_mspline2 <- plotfun(mod1_mspline2, title = "M-splines with df = 10")
cs <-bayesplot::bayesplot_grid(p_exp,
                               p_weibull,
                               p_gompertz,
                               p_bspline,
                               p_mspline1,
                               p_mspline2,
                               grid_args = list(ncol = 3))

ggplot2::ggsave(cs,filename="../plots/baseHaz_tve_fixedEffs_new1.pdf",height = 6, width=10)


## Record execution time and multicore use
end_time <- Sys.time()
diff_time <- end_time - start_time

#readRDS("../processed/loo_mod1.rds")

mm.sum <- summary(mod1_mspline2)
mm.post <- posterior_survfit(mod1_mspline2)
mm.plot <- plot(mmp)

save(mm.sum,mm.post,mm.plot, file="../processed/bhaz_mm.Rda")

