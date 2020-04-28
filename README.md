---
title: "Choi_Diane_Capstone"
author: "Diane Choi"
date: "April 24, 2020"
output: html_document
---

## 1. Background and Significance

Frontotemporal dementia (FTD) is a debilitating neurodegenerative disorder that affects, as its name suggests, the frontal and/or temporal lobes of the brain. Although often confused with Alzheimer's Disease, FTD is diagnosed at an earlier age, between 40 to 60 years, and has more characteristic effects on behavior, personality, and difficulties in speech and language comprehension. There has been recent excitement in applying proteomics to the neurodegenerative field in hopes of elucidating up- or down-regulated proteins that may provide some answers to the many questions we have. For FTD, the Hales lab is interested in the protein versican, a chondroitin sulfate proteoglycan, whose V2 isoform is highly expressed in the brain and when dysregulated, may play a role in inhibiting axon regeneration and contribute to neurodegenerative disease progression. Proteomics data indicates versican is upregulated in grey matter and downregulated in white matter. We wish to test the data by performing western blot analyses for both grey and white matter in FTD and control cases.

## 2. Experiment

While proteomics data suggests versican is upregulated in grey matter and downregulated in white matter, it is unknown if it is the V2 isoform that expresses these abnormal protein levels. 

## 3. Prediction

If V2 is found in the grey and white brain matter of FTD cases, then it should express different protein levels compared to control cases. 

## 4. Variables 

The predictor variables will be the brain matter (grey or white) and disease cases (FTD or control) (discrete, factoral variables). The dependent variable will be the V2 isoform protein levels (continuous, measured data). These will be determined via western blot quantifications.

## 5. Statistical Hypothesis 

Null hypothesis: FTD V2 protein levels for both grey and white matter should be the same as control V2 protein levels. **Null Hypothesis:**\[\mu_{V2 protein levels} \ FTD = \mu_{V2 protein levels} \ control\]

Alternate hypothesis: FTD V2 proteins levels for both grey and white matter should be different from control V2 protein levels. **Alternate Hypothesis:** \[\mu_{V2 protein levels} \ present =/= \mu_{V2 protein levels} \ control\]

## 6. Statistical Test

Two-way mixed ANOVA will be used to compare protein levels. There are two predictor variables, sample cases and patient brain matter, with only one outcome variable (protein quantification). Sample case factor is completely randomized, while brain matter factor is related to each sample case as we will be taking both grey and white matter (related within each replicate). This is the best statistical test for my experiment, as it allows for simultaneous analysis of the effects of both predictor variables on V2 protein levels. 

## 7. Procedure

An equal number of FTD and control (non-dementia) samples will be randomly selected. Each sample case will be an independent replicate, as each person is randomly chosen and unique from one another (i.e. a biologically independent and unique). Grey and white brain matter will be taken from each sample case, run on a western blot, and immunoblotted with anti-versican. V2 protein levels will be normalized to GAPDH; experimental units will be arbitrary quantification levels. These normalized values will be compared to average grey matter control value to determine fold change. A fold change greater than or equal to 1.5 will be considered statistically significant. Type 1 error threshold will be set at 5% to a 95% confidence interval (the null hypothesis will be rejected with a p-value less than 0.05). Type 2 error tolerance will be set at 20%, so the sample size will be based on a power of 80%. 

## 8. Plot

```{r}
library(tidyr)
library(tidyverse)
library(datapasta)
library(ez)
library(knitr)

ID <- c(1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3)
Disease <- c(rep("control",6), rep("FTD", 6))
matter <- rep(c("grey", "white"), each=1)
protein_norm <- c(0.7143, 1.1475, 1.1888, 1.8534, 1.0968, 2.1278, 1.5689, 0.4612, 1.3156, 0.6420, 2.2383, 0.6647) #V2 protein levels normalized to average control grey matter

brain <- data.frame(ID, Disease, matter, protein_norm)

protein_plot <- ggplot(brain, aes(x=matter, y=protein_norm, color=Disease))+
 geom_boxplot()+
  geom_jitter(size=3)+
  theme_classic()+
  xlab("Brain Matter")+
  ylab("Average Protein Quantification [arbitrary]")+
  ggtitle("V2 Protein Levels in Grey and White Matter of FTD vs Control Cases")
protein_plot

```

**Pilot data from my Hales lab rotation is plotted above. We see that V2 protein levels are upregulated in grey brain matter and downregulated in white matter of FTD cases, compared to controls.**

## 9. Monte Carlo

```{r warning=FALSE}
#Determining mean and standard deviation for disease and brain matter states
brain %>% group_by(Disease, matter) %>% 
  summarise(m=mean(protein_norm),
            sd=sd(protein_norm))

#initializer data
m1 = 0.9999667 #mean control grey
m2 = 1.7095667 #mean control white	
m3 = 1.7076000 #mean FTD grey
m4 = 0.5893000 #mean FTD white
#f = 1.5 #minimal scientifically relevant fold-to-basal effect of our treatment
sd1 = 0.2516348 #sd control grey
sd2 = 0.5057302 #sd control white
sd3 = 0.4767306 #sd FTD grey
sd4 = 0.1115169 #sd FTD white
n = 3 # number of independent replicates per group
sims = 100 #number of Monte Carlo simulations to run. 

CRdataMaker <- function(m1, m2, m3, m4, sd1, sd2, sd3, sd4, n) { 
  
  
  a1 <- rnorm(n, m1, sd1) #control grey
  a2 <- rnorm(n, m2, sd2)#control white
  a3 <- rnorm(n, m3, sd3) #FTD grey
  a4 <- rnorm(n, m4, sd4) #FTD white
    
    Disease <- c(rep(c("control", "control","FTD","FTD"), each=n))
    matter <- c(rep(c("grey","white","grey","white"), each=n))
    Outcome <- c(a1, a2, a3, a4)
    id <- as.factor(rep(c(1:n)))
    df <-data.frame(Disease, matter, Outcome, id)
    }

dat <- CRdataMaker(m1, m2, m3, m4, sd1, sd2, sd3, sd4, n)

pval <- replicate(
  sims, {
 
    sample.df <- CRdataMaker(m1, m2, m3, m4, sd1, sd2, sd3, sd4, n)
    
    sim.ezaov <- ezANOVA(
            data = sample.df, 
            within = matter,
            wid = id,
            dv = Outcome,
            between = Disease,
            type = 2
            )
  
  pval <- sim.ezaov$ANOVA[3,5] #p-value of interaction of both Disease:matter on V2 protein levels
    
    }
  )

pwr.pct <- sum(pval<0.05)/sims*100
paste(pwr.pct, sep="", "% power. Change 'n' in your initializer for higher or lower power.")
```

**Monte Carlo power analysis indicates a sample size of n = 3 (where n is the number of independent replicates per group) is necessary to achieve a power of 86%, and thereby meets the type2 error tolerance threshold.**

