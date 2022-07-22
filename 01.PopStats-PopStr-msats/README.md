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
ar$min.all # the number of alleles used for rarefaction
ARmean <- apply(ar$Ar, 2, mean) # mean per population
ARse <- apply(ar$Ar, 2, se) #standard error of the mean AR per population
```

**STEP 3:** Calculate effective population size \
\
3a: Load the [helper script](https://github.com/squisquater/Tutorials/blob/main/01.PopStats-PopStr-msats/write_genepop_function.R) that writes genepop file from adegenet object. This is used to create input for NeEstimator2. *You will need to download this file to create a copy of it in your working directory.* 

```
source("write_genepop_function.R")
```
