# Basic population structure and population genetic statistics for microsatellite data
*Your genotypes should be in structure format. See [discreteRFgenotypes.stru](https://github.com/squisquater/Tutorials/blob/main/01.PopStats-PopStr-msats/discreteRFgenotypes.stru) for an example. This tutorial will use this inupt file and associated parameters. You will therefore need to modify all the scripts to suite your study system/needs.*

## Population Structure
**STEP 0:** \
**STEP 1:** \
**STEP 2:**

## Population Summary Statistics and Ne
These statistics will be calculated for discrete groups. It is therefore suggested that you use the results from the population structure analyses to determine how to discretize your data (e.g., geographic clusters, genetic clusters, etc.)

**STEP 0:** Install/load the required packages in R and set your working directory

```
install.packages("adegenet")
install.packages("hierfstat")

library(adegenet)
library(hierfstat)

setwd("~/path/to/working/directory)")
```

**STEP 1:** Read in your structure file using the 'read.structure' command {adegenet} to convert it to a GENIND object. \

```
data <- read.structure("discreteRFgenotypes.stru", n.ind = 301, n.loc = 31, onerowperind = T, col.lab = 1, col.pop = 2, col.others = 0, row.marknames = 1, NA.char = "-9", ask = TRUE, quiet = T)
```
*You can type 'read.structure' into the 'R' help menu to find out more about each of these options/parameters.* 

**STEP 2:** Calculate basic summmary statistics (i.e., Heterozygosity, Allelic Richness, Inbreeding, Pairwise Fst, etc.) \
\
2a: Write a function to calculate the standard error
```
se <- function(x){
  sd(x)/(length(x)^0.5)
} 
```

2b: Calculate allelic richness in 'hierfstat'
```
ar <- allelic.richness(d_discrete) # per locus

ar$min.all # the number of alleles used for rarefaction (default is min number of individuals in a group * 2)

ARmean <- apply(ar$Ar, 2, mean) # mean per population
#    WAC      ORC      LAS       SN       SV     CANN 
# 4.126517 4.382390 3.483252 4.814377 5.077759 5.990235 

ARse <- apply(ar$Ar, 2, se) #standard error of the mean AR per population
#    WAC       ORC       LAS        SN        SV      CANN 
# 0.2633509 0.2881938 0.2016592 0.3015382 0.2971134 0.4114317 
```

2c: Calculate Observed and Expected Heterozygosity, Fst, and Fis.
```
basicstats <- basic.stats(d_discrete)

Nind <- apply(basicstats$n.ind.samp, 2, max, na.rm=TRUE) # Number of individuals per group
# WAC  ORC  LAS   SN   SV CANN 
#  28   35   27   44   58  109 

Homean <- apply(basicstats$Ho,2, mean) #mean observed heterozygosity per group
#     WAC       ORC       LAS        SN        SV      CANN 
# 0.5838710 0.5315194 0.4832387 0.6701258 0.5942613 0.6302258 

Hose <- apply(basicstats$Ho,2, se) #Std.error of observed heterozygosity per group
#     WAC        ORC        LAS         SN         SV       CANN 
# 0.03725410 0.02489624 0.02897126 0.02607292 0.03542426 0.02609991 

Hsmean <- apply(basicstats$Hs, 2, mean) #mean expected heterozygosity per group
#     WAC       ORC       LAS        SN        SV      CANN 
# 0.5840710 0.5473484 0.5063484 0.6615613 0.6344742 0.6952161 

Hsse <- apply(basicstats$Hs, 2, se) #Std.error of expected heterozygosity per group
#     WAC        ORC        LAS         SN         SV       CANN 
# 0.03379070 0.02555326 0.02582648 0.02238832 0.03355569 0.02709591 
```
2d: Calculate the Weir and Cockerham pairwise-Fst between populations
```
pairwise.WCfst(d_discrete, diploid = TRUE)

#           WAC       ORC       LAS        SN         SV       CANN
# WAC         NA 0.2140956 0.2369856 0.1856994 0.15380814 0.13856409
# ORC  0.2140956        NA 0.2860115 0.1739258 0.19034948 0.13596701
# LAS  0.2369856 0.2860115        NA 0.2071802 0.20901366 0.16926347
# SN   0.1856994 0.1739258 0.2071802        NA 0.13027340 0.10470462
# SV   0.1538081 0.1903495 0.2090137 0.1302734         NA 0.08195674
# CANN 0.1385641 0.1359670 0.1692635 0.1047046 0.08195674         NA
```

**STEP 3:** Calculate effective population size \
\
3a: Load the [helper script](https://github.com/squisquater/Tutorials/blob/main/01.PopStats-PopStr-msats/write_genepop_function.R) that writes genepop file from adegenet object. This is used to create input for NeEstimator2. \
*You will need to download this file to create a copy of it in your working directory.* 

```
source("write_genepop_function.R")
```

convert your GENIND object to a GENEPOP file for use in the NeEstimator GUI.

```
write.genepop(d_discrete, "discretepops", "populations with discrete population structure in the western US")
```
3b: Run NeEstimator2 through GUI, using LD method with random mating and 0.05 pcrit. You can install the program [here](http://www.molecularfisherieslaboratory.com.au/neestimator-software/) (you will need to fill out a brief form with your info to download). Youmay also need to download [java](https://www.java.com/en/) if it is not already installed on your computer.
