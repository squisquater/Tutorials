# Population Structure Tutorial

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

## Take a Break and Learn how to plot these results in ArcGIS

## Another way to look at population structure is using a Principle Components Analysis (PCA)

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


**STEP 2:** Generate a PCA \
\
####PCA
pca_result <- dudi.pca(genind_obj, scannf = FALSE, nf = 3)  # nf is the number of axes
s.class(pca_result$li, pop(genind_obj))

# Save original plot settings
old_par <- par(no.readonly = TRUE)

# Create space for the inset plot
par(fig = c(0.7, 0.9, 0.7, 0.9), new = TRUE, mar = c(1, 1, 1, 1))

# Eigenvalues plot (scree plot)
barplot(pca_result$eig, main = "Eigenvalues", xlab = "Axis", ylab = "Eigenvalue")

# Restore original plot settings
par(old_par)


