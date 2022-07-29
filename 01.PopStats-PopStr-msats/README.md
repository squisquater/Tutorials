# Basic population structure and population genetic statistics for microsatellite data

This tutorial uses data from:

Quinn, C.B., Preckler-Quisquater, S., Akins, J.R. et al. **Contrasting genetic trajectories of endangered and expanding red fox populations in the western U.S.** Heredity (2022). https://doi.org/10.1038/s41437-022-00522-4

Your genotypes should be in the proper format to run a structure analyis. See [discreteRFgenotypes_popmod.txt](https://github.com/squisquater/Tutorials/blob/main/01.PopStats-PopStr-msats/discreteRFgenotypes_popmod.txt) for an example. This tutorial will use this input file and associated parameters. To run these analyses on your own data, you will therefore need to modify all the scripts to align with the files associated with your project and modify parameters accordingly.*

You can clone this repository to your own computer by entering the following line of code in your terminal.
```
git clone https://github.com/squisquater/Tutorials.git
```
Note: if you run into an error regarding the *xcrun: error: invalid active developer path* you may have to install xcode to your computer before you can clone the repository.

## Population Structure
Note that structure will accommodate population data but expects a number as opposed to text. For reference, WAC = 1, ORC = 2, LAS = 3, SN = 4, SV = 5, CANN = 6. \
\
**STEP 0:** Open the Structure program \
\
**STEP 1:** Set up your Structure project
Select **'New Project'** from the **File** dropdown menu which will open the **Project Wizard** 
* Step 1 of 4: Project Information
  *  Name your project (i.e. "WestRF")
  *  Select Directory where your output files will go **'discreteRFgenotypes.stru'** file is located
  *  Choose data file (i.e. **'discreteRFgenotypes.stru'**)
* Step 2 of 4: Information of Input Dataset
  *  Number of individuals: **301**
  *  Ploidy of data (This should be 2 for diploids)
  *  Number of loci: **31**
  *  Missing data value: **-9**
* Step 3 of 4: Format of Input Dataset
  *  Select **'Row of marker names'**
  *  Select **'Data file stores data for individuals in a single line'**
* Step 4 of 4: Format of Input Dataset (continued)
  *  Select **'Individual ID for each individual'**
  *  Select **'Putative population origin for each individual'** ONLY IF YOU HAVE THIS INFORMATION AS A COLUMN
* Click **'Finish'**
* Click **'Proceed'**

This should load your input file and all associated project information

**STEP 2:** Create a new parameter set
* Select **'New...'** from the **'Parameter Set'** dropdown menu
* Under the **'Run Length'** tab choose the:
  * Length of your Burnin Period: 10000 (suggested)
  * Number of MCMC Reps after Burnin: 50000 (suggested)
* Under the **'Ancestry Model'** tab select:
  * **Use Admixture Model** (default)
* Under the **'Allele Frequency Model'** tab select:
  * **Allele Frequencies Correlated** (default)
* Under the **'Advanced'** tab select:
  * **Compute probability of the data (for estimating k)** (default)
  * **Print Q-hat** (this generates a nice tab-delimited text file with ancestry proportions)
* Click **'OK'**
* Name the Parameter Set (i.e. 'WestRF_10k50k')

**STEP 3:** Run Structure
* Select **'Run'** from the **'Parameter Set'** dropdown menu
* Set the number of populations assumed. 
 * You should test out multiple iterations of K.
 * You can also run multiple replicates of the same K-value in order to use structure harvester and determine the K-value with the highest likelihood.

**STEP 4:** View Results
* Select the value of K you want to view from the Results folder.
* Select **Bar Plot** from the Results menu to view the admixture proportions. You can sort by POP if they are not already in order.
* These results are also available in the Results folder and can be plotted in R or excel to customize the figure (i.e. colors, labels, etc)

**STEP 5:** Structure Harvester
* Zip your Results folder
* Navigate to the [Structure Harvester](https://taylor0.biology.ucla.edu/structureHarvester/) website
* Load the .zip file
* Click the **Harvest!** button
* This will give you the likelihood scores associated with each K value.
* This will also show the delta-K values (Evanno method), assuming you ran >1 iteration of each K value. 

## Population Summary Statistics and Ne
These statistics will be calculated for discrete groups. It is therefore suggested that you use the results from the population structure analyses to determine how to discretize your data (e.g., geographic clusters, genetic clusters, etc.)

You will also need to output file from structure (i.e. 'project_data.stru') as an input for these analyses. This should be located in the project directory you specified above. This file is also included [here](https://github.com/squisquater/Tutorials/blob/main/01.PopStats-PopStr-msats/project_data.stru) for reference.

**STEP 0:** Install/load the required packages in R and set your working directory

```
install.packages("adegenet")
install.packages("hierfstat")
install.packages("reshape2")
install.packages("ggplot2")

library(adegenet)
library(hierfstat)
library(reshape2)
library(ggplot2

setwd("~/path/to/working/directory")
```

**STEP 1:** Read in your structure file using the 'read.structure' command {adegenet} to convert it to a GENIND object.

```
d_discrete <- read.structure("project_data.stru", n.ind = 301, n.loc = 31, onerowperind = T, col.lab = 1, col.pop = 2, col.others = 0, row.marknames = 1, NA.char = "-9", ask = TRUE, quiet = T)
```
*You can type 'read.structure' into the 'R' help menu to find out more about each of these options/parameters.* 

**STEP 2:** Calculate basic summmary statistics (i.e., Heterozygosity, Allelic Richness, Pairwise Fst, etc.) \
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
ARmean
#     1        2        3        4        5        6 
# 4.126517 4.382390 3.483252 4.814377 5.077759 5.990235 

ARse <- apply(ar$Ar, 2, se) #standard error of the mean AR per population
ARse
#      1         2         3         4         5         6 
# 0.2633509 0.2881938 0.2016592 0.3015382 0.2971134 0.4114317 
```

2c: Calculate Observed and Expected Heterozygosity, Fst, and Fis.
```
basicstats <- basic.stats(d_discrete)

Nind <- apply(basicstats$n.ind.samp, 2, max, na.rm=TRUE) # Number of individuals per group
Nind
#  1   2   3   4   5   6 
# 28  35  27  44  58 109 

Homean <- apply(basicstats$Ho,2, mean) #mean observed heterozygosity per group
Homean
#      1         2         3         4         5         6 
# 0.5838710 0.5315194 0.4832387 0.6701258 0.5942613 0.6302258 

Hose <- apply(basicstats$Ho,2, se) #Std.error of observed heterozygosity per group
Hose
#      1          2          3          4          5          6 
# 0.03725410 0.02489624 0.02897126 0.02607292 0.03542426 0.02609991 

Hsmean <- apply(basicstats$Hs, 2, mean) #mean expected heterozygosity per group
Hsmean
#     1         2         3         4         5         6 
# 0.5840710 0.5473484 0.5063484 0.6615613 0.6344742 0.6952161 

Hsse <- apply(basicstats$Hs, 2, se) #Std.error of expected heterozygosity per group
Hsse
#      1          2          3          4          5          6 
# 0.03379070 0.02555326 0.02582648 0.02238832 0.03355569 0.02709591 
```
2d: Calculate the Weir and Cockerham pairwise-Fst between populations
```
Fst.mat <- pairwise.WCfst(d_discrete, diploid = TRUE)

#         1         2         3         4          5          6
# 1        NA 0.2140956 0.2369856 0.1856994 0.15380814 0.13856409
# 2 0.2140956        NA 0.2860115 0.1739258 0.19034948 0.13596701
# 3 0.2369856 0.2860115        NA 0.2071802 0.20901366 0.16926347
# 4 0.1856994 0.1739258 0.2071802        NA 0.13027340 0.10470462
# 5 0.1538081 0.1903495 0.2090137 0.1302734         NA 0.08195674
# 6 0.1385641 0.1359670 0.1692635 0.1047046 0.08195674         NA

```
Rename populations so they have the correct pop names

```
colnames(Fst.mat) <- c("WAC", "ORC", "LAS", "SN","SV", "CANN")
rownames(Fst.mat) <- colnames(Fst.mat)

#Set up the matrix to be plotted by ggplot
Fst.mat.tri <- Fst.mat
Fst.mat.tri[lower.tri(Fst.mat, diag=TRUE)] <- NA
#fst.mat = data.matrix(fst.df.tri)
melted <- melt(Fst.mat.tri, na.rm =TRUE)
```
Plot it!
```
Fst.plot <- ggplot(data = melted, aes(Var2, Var1, fill = value))+ 
  geom_tile(color = "white") + 
  scale_fill_gradient(low = "white", high = "red", name="FST")  + 
  ggtitle("Pairwise FST, WC (1984)") +
                       labs( x = "Population", y = "Population") + 
                       theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1), axis.text.y = element_text(size = 12)) + 
            coord_fixed()
            

Fst.plot
```
<img align="center" src="/01.PopStats-PopStr-msats/Pairwise.Fst.png" width="500">  \

Save the plot to your working directory
```
ggsave(plot = Fst.plot, width = 5, height = 4, dpi = 300, filename = "Pairwise.Fst.png")
```

**STEP 3:** Calculate effective population size \
\
3a: Load the [helper script](https://github.com/squisquater/Tutorials/blob/main/01.PopStats-PopStr-msats/write_genepop_function.R) that writes genepop file from adegenet object. This is used to create input for NeEstimator2. \

```
source("write_genepop_function.R")
```

convert your GENIND object to a GENEPOP file for use in the NeEstimator GUI.

```
write.genepop(d_discrete, "discretepops", "populations with discrete population structure in the western US")
```
3b: Run NeEstimator2 through GUI, using LD method with random mating and 0.05 pcrit. You can download the program [here](http://www.molecularfisherieslaboratory.com.au/neestimator-software/) (you will need to fill out a brief form with your info to download). You may also need to download [java](https://www.java.com/en/) if it is not already installed on your computer. Unzip the directory. \
\
Make sure the 'discretepops_genepop.gen' file is located within the directory where your executable is found (NeEstimator2x1.jar). \
\
Click on the NeEstimator2x1.jar to open NeEstimator. \
\
Select the 'discretepops_genepop.gen' file from the **Choose File:** dropdown menu \
\
Select the **Linkage Disequilibrium** Model \
\
Select the **Random Mating** Model
\
De-select the **Heterozygote Excess**, **Molecular Coancestry**, and **Temporal, Plan II** Models \
\
Select a Critical Value of **0.05** - Remove/Delete any other critical values \
\
**>>> Run Ne >>>**

3a: You should end up with an output file that has the Ne and associated 95% confidence intervals for each of your input populations. 
```
...
...
Population     1 [S13-0964]  (Number of Individuals = 28)
****************
--------------------------------------------------
Lowest Allele Frequency Used     0.050         0+
--------------------------------------------------


LINKAGE DISEQUILIBRIUM METHOD

Harmonic Mean Sample Size =       27.3        27.2
Independent Comparisons =       3968        4614
OverAll r^2 =                 0.077276    0.071208
Expected r^2 Sample =         0.040965    0.041128
Estimated Ne^ =                    6.2         8.1

95% CIs for Ne^
* Parametric                       5.2         7.1
                                   7.1         9.2

* JackKnife on Samples             2.9         4.2
                                  10.8        14.0
                                  
...
...
```
*Use the Ne estimate that is associated with the 0.05 critical value (i.e. 6.2 in the above example) and the 95% confidence intervals associated with the jacknifing approach (i.e., 2.9 & 10.8)* \
\
Refer to Table S7 from [Quinn et al. (2022)](https://www.nature.com/articles/s41437-022-00522-4#Sec21) as a reference for reporting results from these analyses.

| Population                                                                                       | n  | H<sub>E</sub> (SE)       | H<sub>O</sub> (SE)       | AR (SE)     | N<sub>e</sub> (95% CI) |
| ------------------------------------------------------------------------------------------------------------ | -- | ------------- | ------------- | ----------- | ----------------------- |
| Washington Cascade                                                                                           | 28 | 0.584 (0.034) | 0.584 (0.037) | 4.13 (0.26) | 6.2 (2.9-10.8)          |
| Oregon Cascade                                                                                               | 35 | 0.547 (0.026) | 0.532 (0.025) | 4.38 (0.29) | 9.3 (3.9-17.3)          |
| Lassen Cascade                                                                                               | 27 | 0.506 (0.026) | 0.483 (0.029) | 3.48 (0.20) | 2.6 (1.7-5.4)           |
| Sierra Nevada                                                                                                | 44 | 0.662 (0.022) | 0.67 (0.026)  | 4.81 (0.30) | 3.5 (2.9-5.4)           |



