
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

# 2020-07-13 (EN)
Fit Cox Proportional Hazard models to lifespan data
Visualize HRs with Forest plots
Next, do reaction norms for Fig. 2 in Monaghan et al, (2008)

# 2020-07-04 (EN)
Identified duplicated image: different flip dates and egg counts
> kk[,c(2:5,11:12)]
    startDate sampleDate age        id camera_id numEggs
137 2/20/2018   4/2/2018  41 EE_DR_DR1      9834      19
189 2/20/2018  4/16/2018  55 EE_DR_DR2      9834      33
143 2/20/2018   4/2/2018  41  EE_C_DR2      9835      24
190 2/20/2018  4/16/2018  55  EE_DR_C2      9835      15

These may be wrongly named - @Andrew, please check and advise!

# 2020-07-02 (EN)
  - bash script `rename_files.sh` adds a random number to the image name
  - A copy of these names kept
  - bash script `simplify.sh` deletes original name leaving just a random number
Made a new data sheet with columns image_id and numEggs.

# 2020-07-01 (EN)
Andrew and De'anne each to recount all images - average out the counts
I replaced image names with random numbers: 

# 2020-06-28 (EN)
Integrate recounted image data
Investigate and compare old and new data - decide approach forward
Differences are large; need for recounting

# 2020-06-27 (EN)
Upgrade figures for descriptive statistics to publication quality

# 2020-06-26 (EN)
Add plotting features to script `Fec_perFem.Rmd`

# 2020-06-04 (EN)
Exploratory analysis of per female fecundity time series.
30 images flagged as outliers. Andrew, pliz blindly recount these images:
Run `Fec_perFem.Rmd`. The last line: spot.check <- rbind(dat[which(dat$eggpFemDay>10),],dat[which(dat$eggpFemDay<1 & dat$age<36),])

# 2020-06-01 (EN)
A bit of cleaning up in the scripts directory

# 2020-05-29
Did more work with the pairwise testing. Used holm and bonf testing to show values < E-16 so they are defintiely significantly different. 

# 2020-05-24
Created some more figures to look at our data and made some simple pairwise tests. Signifcant differences between mean eggs per female but unclear about total number of eggs per female. Will investigate.  

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

















