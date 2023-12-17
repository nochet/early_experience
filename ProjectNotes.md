
---
title: "Early Experience Project Record & Notes"
author: "Enoch Ng'oma (EN), Andrew Jones (AJ) and Elizabeth King (EGK)"
date: "7/17/2018"
---


## Project Record

### Summary: This is data from a factorial experiment with Drosophila melanogaster in which larvae were raised on control (C) and dietary restriction (DR) nutrition conditions. Eclosed flies from each treatment were split between C and DR diets followed by measurement of fecundity and lifespan of adult flies. Two replicates were made for each diet treatment resulting into 8 total cages.

`InitialProcess.R`
- `source("PreProcess_lifespan_functions.R")`
- Process raw data for male and female separately.
Input: `/data/processed/ee_lifespan_only.txt`
Output: `/data/processed/Female_events_eelife.txt` and `/data/processed/Male_events_eelife.txt`


`ee_InitialProcess.R`
- attempt to include both sexes. Not usable yet.
Input: `/data/processed/ee_lifespan_only.txt`

STEP 1: Initial processing of raw data

STEP 2: Exploration analyses


STEP 3: Two-way ANOVA

STEP 4: Generalized linear mixed models

STEP 4: A Bayesian approach



## Project Notes 

## 2023-12-16
Made new plots of eggs/female/3hr in `fec_timeSeries`
Changepoint analysis in egg time series
Ref: https://zealous-wiles-e22e83.netlify.app/post/changepoint_detection/
Following parameters used in changeplot():
test.sat="Normal"
method="AMOC" and method="PELT"
pen.value=DEFAULT

## 2023-12-14 Hope this is the final push!!
Remake the design figure

## 2022-08-01
Update time series

## 2022-06-05
`fec_timeSeries.Rmd` - an improved time series analysis of fecundity.

## 2022-05-29
Revise visuals and models for fecundity. Models now integrated in `eefecund_viz.Rmd`.

## 2022-01-19
Completed Leslie matrix analysis to compute a proxy of fitness by combining age-specific lifespan and fecundity into a single metric. The script is `fitnessEstimates.Rmd`.The visualization is `/plots/leslie_fitness.pdf`

# Some time here!
`ee_cph_finalMods.Rmd` Analysis of the Cox model and model selection
`ee_cph_model.Rmd` The selected models above come from a general exploration in this script.
`eefecund_viz.Rmd` Contains most  of the fecundity visualizations for raw data

Reporting Cox model
HR, CI and P-value

http://www.danieldsjoberg.com/gtsummary/
https://rdrr.io/rforge/tab/man/tabcox.html

## 2021-06-8 (EN)
Anova test on coxme model for a factorial design in survival analysis (useful links):
https://stats.stackexchange.com/questions/509905/understanding-coxme-model-results-with-1-or-2-factors

https://towardsdatascience.com/anovas-three-types-of-estimating-sums-of-squares-don-t-make-the-wrong-choice-91107c77a27a

https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/

Use Type I only when there is a serious theoretical reason for it, use Type II when there is no interaction, use Type III when there is interaction.

But, see:
https://www.r-bloggers.com/2014/10/type-iii-tests-and-r/

Our data is balanced, so type of SS may not be relevant.

From:
https://stat.ethz.ch/pipermail/r-help//2013-July/357466.html
GEE = Generalized Estimating Equations
Two choices.
- If this were a linear model, do you like the GEE approach or a mixed effects approach?  Assume that "subject" is a variable containing a per-subject identifier.

1. GEE approach: add "+ cluster(subject)" to the model statement in coxph
2. Mixed models approach: Add " + (1|subject)" to the model statment in coxme.

From:
https://stats.stackexchange.com/questions/178944/precisely-how-does-rs-coxph-handle-repeated-measures
Including cluster(ID) does not change the point estimates of the parameters. It does change the way that the standard errors are computed however.

More details can be found in Therneau & Grambsch's book Extending the Cox Model, chapter 8.2. Note that in their example, they use method = "breslow" as correction for ties, but also with the default (method = "efron") a similar calculation for the se's will be used, and appears in the summary as "robust se".

If cluster(ID) is used, a "robust" estimate of standard errors is imposed and possible dependence between subjects is measured (e.g. by standard errors and variance scores). Not using cluster(ID), on the other hand, imposes independence on each observation and more "information" is assumed in the data. In more technical terms, the score function for the parameters does not change, but the variance of this score does. A more intuitive argument is that 100 observations on 100 individuals provide more information than 100 observations on 10 individuals (or clusters).

Vague indeed. In short, +frailty(ID) in coxph() fits standard frailty models with gamma or log-normal random effects and with non-parametric baseline hazard / intensity. frailtypack uses parametric baseline (also flexible versions with splines or piecewise constant functions) and also fits more complicated models, such as correlated frailty, nested frailty, etc.

Finally, +cluster() is somewhat in the spirit of GEE, in that you take the score equations from a likelihood with independent observations, and use a different "robust" estimator for the standard errors.


## 2021-02-8 (EN)
Consider these two models:
`~ larvalTreat + adultTreat + strata(sex)` : cxph1.2.1
`~ larvalTreat + adultTreat + cageID + strata(sex)` : cxph1.2.2
`~ larvalTreat + adultTreat + strata(sex) + (1 | cageID)` # better, but cox.zph impossible cph3
Last two seem pretty good. See how Bayes analysis compare.
Note: All model summaries indicate larval diet not significant. However, and anova on cph3 has larval diet very significant. Is this the right thing to do?

Next step: Write results in manuscript.

## 2021-01-15 (EN)
`ee_coxme_model.Rmd`
Cox model using coxme (with random effect, cageID - see vignette)
https://cran.r-project.org/web/packages/coxme/vignettes/coxme.pdf


## 2020-12-11 (EN)
Make final plots: 
 - fecundity exploratory
 - lifespan exploratory
Review CoxPH models for the ms


## 2020-11-27 (EN)
Rename `Fec_perFem.Rmd` to `process_fecund.Rmd`. This script is reserved for general processing for the fecundity data. Start a new script for visualization `eefecund_viz.Rmd`.

Visualize.


## 2020-11-25 (EN)
In Andrew initial count (aja1), images labeled 1000, 1001, 1002 don't exist on server. However, set2 on server contains images 0001, 0002, and 0003. The camera numbering had hit the limit of 9999 and restarted from 1 in a new folder called set2. I am therefore renaming 1000, 1001, and 1002 to 0001, 0002, and 0003 in the analysis. After all, tallies seem to agree.

Also, the following images are absent in aj1: IMG_9884.JPG and IMG_9885.JPG. Both have been recounted by both people in recounts. Images 9834 and 9835 are duplicated. The date of flip clearly suggest the duplicates on 4/2/18 really belong to 4/16/18 matching 9834 and 9835. I there resolve the duplicate by transferring duplicates with a later date to 4/16/18:

> aj1[aj1$sampleDate=="4/16/18",]
    collectedEggs startDate sampleDate age        id  deadFemale deadMale censored camera_id numEggs 
185        2/8/18   2/20/18    4/16/18  55 EE_DR_DR1           5        2        0      9888       7       
186        2/8/18   2/20/18    4/16/18  55  EE_DR_C1           2        1        1      9883      19       
187        2/8/18   2/20/18    4/16/18  55   EE_C_C1          11       14        2      9886       7       
188        2/8/18   2/20/18    4/16/18  55  EE_C_DR1           2        3        0      9887      22       
189        2/9/18   2/20/18    4/16/18  55 EE_DR_DR2           4        2        9      9834      33       
190        2/9/18   2/20/18    4/16/18  55  EE_DR_C2           6        2        4      9835      15       
191        2/9/18   2/20/18    4/16/18  55  EE_C_DR2           9        3        3      9889      18       
192        2/9/18   2/20/18    4/16/18  55   EE_C_C2           5        2        1      9882      35       
> aj1[aj1$sampleDate=="4/2/18",]
    collectedEggs startDate sampleDate age        id  deadFemale deadMale censored camera_id numEggs 
137        2/8/18   2/20/18     4/2/18  41 EE_DR_DR1          13        4        0      9834      19      
138        2/8/18   2/20/18     4/2/18  41  EE_DR_C1          10        4        0      9841      14      
139        2/8/18   2/20/18     4/2/18  41   EE_C_C1           4        3        0      9838       5      
140        2/8/18   2/20/18     4/2/18  41  EE_C_DR1           2        1        0      9837       9      
141        2/9/18   2/20/18     4/2/18  41 EE_DR_DR2          23        7        0      9836      17      
142        2/9/18   2/20/18     4/2/18  41  EE_DR_C2          17        4        0      9839      18      
143        2/9/18   2/20/18     4/2/18  41  EE_C_DR2          19        6        1      9835      24      
144        2/9/18   2/20/18     4/2/18  41   EE_C_C2           7        4        0      9840       4      

Flanking 9884 and 9885
> aj1[aj1$camera_id==9883,]
    collectedEggs startDate sampleDate age       id deadFemale deadMale censored camera_id numEggs 
186        2/8/18   2/20/18    4/16/18  55 EE_DR_C1            2        1        1      9883      19       
> aj1[aj1$camera_id==9886,]
    collectedEggs startDate sampleDate age      id  deadFemale deadMale censored camera_id numEggs 
187        2/8/18   2/20/18    4/16/18  55 EE_C_C1            11       14        2      9886       7 


## 2020-11-25 (EN)
Thorough check of the recounts. A few mis-typed image numbers observed and manually rechecked with original files. Excel spreadsheet written `/processed/recounts_checked.xlsx`. The rows with color is where reconciliation was done. The original recounts by Andrew `/original/ee_temp_Values_Sheet_AJ.xlsx` and De'anne `/original/EE_counts_DD.xlsx` may be read to confirm. The images found in the random number file but missing in recounts are non-egg images (markers). The file for further analysis is `/processed/recounts_checked.csv`.

## 2020-11-23 (EN)
Process three independent egg counts for final analysis: `Fec_perFem.Rmd`


## 2020-11-13 (EN)
Piecewise regression: 
 1) 1e4 iter, tve mod1, df=5 - slowest chain took 9.3 hrs
 2) 1e4 iter, tve mod1, df=10 - slowest chain took 9.3 hrs
 
 loos were halted.
 
## 2020-11-03 (EN)
Try Bayesian peicewise regression on mod1 (the favored model), 5e4 iter
Chain 4:  Elapsed Time: 38.51 hrs
Chain 3:  Elapsed Time: 38.95 hrs
Chain 2:  Elapsed Time: 38.98 hrs
Chain 1:  Elapsed Time: 39.16 hrs
This is for one model (tve, df=5)!!!!! Cancel and reduce iter to 1e4

## 2020-11-09 (EN)
Select best model based on leave-one-out (loo)
Plot model results
Script `stanSurvEval.Rmd`


## 2020-11-08 (EN)
stan_surv() for time varying models can accept "delayed entry" data (i.e. left truncated). Our data is 'right-censored', so we do not need the "start-stop" data structure (https://arxiv.org/pdf/2002.09633.pdf, page 11)

The tve() function is used in the model formula to state that we want a time-varying effect (i.e. a time-varying coefficient) to be estimated for a variable. By default, a cubic B-spline basis with 3 degrees of freedom (i.e. two boundary knots placed at the limits of the range of event times, but no internal knots) is used for modeling the time-varying log hazard ratio. To place more knots inside the time interval, use df>4 e.g. Surv(eventtime, status) ~ tve(trt, df = 4, degree = 2).

Hazard ratio = effect of treatment on the hazard of the event

Interpretation of hazard rate vs time plots:
From Figure 4 we can see how the hazard ratio (i.e. the effect of treatment on the hazard of the event) changes as a function of time. The treatment appears to be protective during the first few years following baseline (i.e. HR < 1), and then the treatment appears to become harmful after about 2.5 years post-baseline. (https://arxiv.org/pdf/2002.09633.pdf, page 33)

## 2020-10-28 (EN)
Reporting example:
["Comparing the models on PSIS-LOO reveals an estimated difference in elpd of 10.2 (with a standard error of 5.1) in favor of Model A."]

Do loo (leave one out) cross-validation for each model
loop.bs <- loo(mod1_bspline)
Warning message:
Found 1 observation(s) with a pareto_k > 0.7. We recommend calling 'loo' again with argument 'k_threshold = 0.7' in order to calculate the ELPD without the assumption that these observations are negligible. This will refit the model 1 times to compute the ELPDs for the problematic observations directly.

loo gives us warnings about the Pareto diagnostics, which indicate that for some observations the leave-one-out posteriors are different enough from the full posterior that importance-sampling is not able to correct the difference.https://mc-stan.org/loo/articles/loo2-example.html

Interpretation of loo:

elpd_loo = expected log predictive density 

When a model is well-specified,  we expect the estimated effective number of parameters (p_loo) to be smaller than or similar to the total number of parameters in the model. If not, it is indication of model misspecification.


## 2020-09-12 (EN)
Start final time-varying model in stan_surv

## 2020-09-07 (AJ)
Got three rechecked data table combined and mean taken of them. IMG_9318 DR_C2 has a large number of eggs: 2867. 


### 2020-08-229 (EN)
Move cox models in new script `ee_coxph_models.Rmd` - `eeSurv_viz.Rmd has gotten too long

Model with all terms and a 3-way interaction
Interaction model with a compound variable (2 diets and sex)
Model where sex is fit additively
Model with larval diet fit additively

# Summary from these simple models
1. LOOKS LIKE SEX HAS A MAJOR IMPACT - CONSIDER SEXES SEPARATELY
2. DEAL WITH NON-PH COVARIATES
3. cageID HAS EFFECT, SO SHOULD BE CONTROLLED FOR AS RANDOM EFFECT

### 2020-08-28 (EN)
Suggestions from DataPhiles:

1. cph: Start with a model that fits all terms additively plus a three-way interaction (larval:adult:sex); progressively drop nsig terms

2. Alternatively, explore simplified models - with just one term each on the right side of "~". See how results agrees with KM plots. 
e.g.:
    - split data by sex
    - combine larval_adult diet ids
    - fit model
Then carefully build models that are more complex from intuition gained

3. Before trying cageID as a random effect, fit it as fixed effect first. If sig, then use as random effect, otherwise drop it.

4. Do pp_checks in Bayesian models to discount inadequate models


### 2020-08-23 (EN)

`/scripts/eeSurv_mixedEffects.Rmd`
Base hazards with stan_surv()
Piece-wise regression with stan_surv()
Multi-level models with stan_surv()

### 2020-08-18 (EN)
stan_surv() again:

To install just the function on top of the CRAN package, do `remotes::install_github("stan-dev/rstanarm", ref = "feature/survival") ` instead of the git clone (thanks to Kevin Middleton for identifying this approach)

The documentation module of rstanarm however fails to compile and seems to deactivate package-wide - at least on nivalis.

Useful online guides:

https://github.com/stan-dev/rstanarm/blob/feature/survival/vignettes/surv.Rmd

https://rstudio-pubs-static.s3.amazonaws.com/438966_3b8a25efb9b84454b8d69b7a15e3ebc5.html

This is for a piece-wise exponential model:
https://rpubs.com/kaz_yos/surv_stan_piecewise1

stan_jm:
https://rstudio-pubs-static.s3.amazonaws.com/430561_45236316f46f4edc99d00ca045a17a4a.html

Copies saved to the Documents section of this repo.


### 2020-08-17 (EN)
Do Bayesian models in brms - spooky???

Do Bayesian multilevel models in rstanarm

Clone the development version of rstanarm:
git clone --single-branch --branch feature/frailty-models https://github.com/stan-dev/rstanarm.git

https://discourse.mc-stan.org/t/huge-r-hat-and-low-effective-size-in-piecewise-exponential-model-when-using-more-than-one-predictor-and-interval/10968/8

### 2020-08-14 (EN)
Trying out timeSplit() in Greg package to deal with unpromportional hazards.

### 2020-08-13 (EN)
library(coxphw): Weighted Cox Regression Using the R Package coxphw. Dunkler et al, 2018.

library(Greg). TimeSplitter(), Gordon (2020). https://cran.r-project.org/web/packages/Greg/vignettes/timeSplitter.html

and,

coxph diagnostics. http://www.sthda.com/english/wiki/cox-model-assumptions

Also,
vignette("timedep", package = "survival")

### 2020-08-12 (EN)
Cox models `ee_coxph_models.Rmd`
Models fail the assumption of proportional hazards: that residuals should be time-independent.

Explore approaches to remedy the problem:
John Fox & Sanford Weisberg last revision: 2018-09-28: Cox Proportional-Hazards Regression for Survival Data in R An Appendix to An R Companion to Applied Regression, third edition

 1) Include interaction of non-PH with age. Fails to converge
 2) Stratify the time variable - effects disappear

### 2020-08-11 (EN)
Survival models revisited: subset data by larval diet, then model survival by adult diet for males and females separately. Script: `eeSurv_viz.Rmd`

### 2020-08-10 (EN)
Reaction norms for lifespan and fecundity (baseR matplot()). Script: `eeSurv_viz.Rmd` and `Fec_perFem.Rmd`

### 2020-07-22 (AJ)
Made new plot to show per female fecundity vs age. Also cleaned up some of the older scripts that had gotten clogged up.

### 2020-07-13 (EN)
Fit Cox Proportional Hazard models to lifespan data
Visualize hazard ratios with Forest plots
Next, do reaction norms for Fig. 2 in Monaghan et al, (2008)

### 2020-07-04 (EN)
Identified duplicated image: different flip dates and egg counts
> kk[,c(2:5,11:12)]
    startDate sampleDate age        id camera_id numEggs
137 2/20/2018   4/2/2018  41 EE_DR_DR1      9834      19
189 2/20/2018  4/16/2018  55 EE_DR_DR2      9834      33
143 2/20/2018   4/2/2018  41  EE_C_DR2      9835      24
190 2/20/2018  4/16/2018  55  EE_DR_C2      9835      15

These may be wrongly named - @Andrew, please check and advise!

### 2020-07-02 (EN)
  - bash script `rename_files.sh` adds a random number to the image name
  - A copy of these names kept
  - bash script `simplify.sh` deletes original name leaving just a random number
Made a new data sheet with columns image_id and numEggs.

### 2020-07-01 (EN)
Andrew and De'anne each to recount all images - average out the counts
I replaced image names with random numbers: 

### 2020-06-28 (EN)
Integrate recounted image data
Investigate and compare old and new data - decide approach forward
Differences are large; need for recounting

### 2020-06-27 (EN)
Upgrade figures for descriptive statistics to publication quality

### 2020-06-26 (EN)
Add plotting features to script `Fec_perFem.Rmd`

### 2020-06-04 (EN)
Exploratory analysis of per female fecundity time series.
30 images flagged as outliers. Andrew, pliz blindly recount these images:
Run `Fec_perFem.Rmd`. The last line: spot.check <- rbind(dat[which(dat$eggpFemDay>10),],dat[which(dat$eggpFemDay<1 & dat$age<36),])

### 2020-06-01 (EN)
A bit of cleaning up in the scripts directory

### 2020-05-29
Did more work with the pairwise testing. Used holm and bonf testing to show values < E-16 so they are definitely significantly different. 

### 2020-05-24
Created some more figures to look at our data and made some simple pairwise tests. Signifcant differences between mean eggs per female but unclear about total number of eggs per female. Will investigate.  

### 2020-04-29
Got Per-Female Fecundity script so that it was standalone. 

### 2020-04-18
Using script created by Enoch yesterday, PerFemaleFec was modified and a new data table was created including the ratio of number of eggs per female at each egg collection date. 

### 2020-04-17
Started script `/scripts/Fec_perFem.Rmd` to calculate number of females at each egg quantification time.

Corrections made to the data file read in the script:
cam_id <- c("9318","9389","9397")
new_count <- c(1043,645,709)

### 2020-04-13
Andrew recounted eggs on 3 unreasonably outlier images. On April 7th he reported these counts (below). The original data was copied to a new file called `Early_Experience_Data_rechecked.txt` and saved to `/early_experience/processed`. This data should be used for further analysis. (Note, the original raw data remains unaltered in `/early_experience/original`.) Slack message from Andrew:

Ok Enoch I recounted those 3 images. I got some more reasonable numbers for all three. Originally they were: 9318-2491, 9389-1136, and 9397-1461. The new numbers are 9318-1043, 9389-645, and 9397-709.


### 2019-03-03
status key: 1 = alive, 2 = dead, 3 = censored

### 2018-07- 17(EN)
- Performed a work-around by combining female and male events. This is done at the top of `eeSurv_model_visual.Rmd`. The merged data frame is in the object `eelife`
- Calculated Kaplan-Meier plots various treatment combinations.

### 2018-07-12 (EN)
- Performed initial processing of lifespan data with `InitialProcess.R` modified to calculate `NstartF` and `NstartM` from observed `nDead` separately.
- Started a parallel script `ee_InitialProcess.R` to include both males and females (fixme!)

















