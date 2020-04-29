
---
title: "Early Experience Project Record & Notes"
author: "Enoch Ng'oma (EN), and Elizabeth King (EK)"
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
# 2020-04-29
Got Per-Female Fecundity script so that it was standalone. 


# 2020-04-18
Using script created by Enoch yesterday, PerFemaleFec was modified and a new datatable was created including the ratio of number of eggs per female at each egg collection date. 

# 2020-04-17
Started script `/scripts/Fec_perFem.Rmd` to calculate number of females at each egg quantification time.

## 2020-04-13
Adrew recounted eggs on 3 unreasonably outlier images. On April 7th he reported these counts (below). The original data was copied to a new file called `Early_Experience_Data_rechecked.txt` and saved to `/early_experience/processed`. This data should be used for further analysis. (Note, the original raw data remains unaltered in `/early_experience/original`.) Slack message from Andrew:

Ok Enoch I recounted those 3 images. I got some more reasonable numbers for all three. Originally they were: 9318-2491, 9389-1136, and 9397-1461. The new numbers are 9318-1043, 9389-645, and 9397-709.


### 2019-03-03
status key: 1 = alive, 2 = dead, 3 = censored

### 2018-07- 17(EN)
- Performed a work-around by combining female and male events. This is done at the top of `eeSurv_model_visual.Rmd`. The merged data frame is in the object `eelife`
- Calculated Kaplan-Meier plots various treatment combinations.

### 2018-07-12 (EN)
- Performed initial processing of lifespan data with `InitialProcess.R` modified to calculate `NstartF` and `NstartM` from observed `nDead` separately.
- Started a parallel script `ee_InitialProcess.R` to include both males and females (fixme!)

















