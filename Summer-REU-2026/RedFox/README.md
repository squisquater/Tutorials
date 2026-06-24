# Population Structure Tutorial

This tutorial uses data from:

Quinn, C.B., Preckler-Quisquater, S., Akins, J.R. et al. **Contrasting genetic trajectories of endangered and expanding red fox populations in the western U.S.** Heredity (2022). https://doi.org/10.1038/s41437-022-00522-4

## Introduction

Population structure analyses are commonly used to determine whether individuals belong to genetically distinct groups and to identify patterns of admixture among populations. Understanding population structure can provide insight into gene flow, dispersal, historical isolation, and population connectivity.

In this tutorial, we will use two complementary approaches:

* **STRUCTURE** – a model-based clustering method that estimates the proportion of an individual's ancestry derived from different genetic clusters.
* **Principal Components Analysis (PCA)** – an ordination method that summarizes genetic variation into a small number of axes and allows patterns of genetic similarity to be visualized.

These methods answer related but slightly different questions and are often used together to evaluate population structure.

### General Workflow

```text
Genotype Data
      ↓
STRUCTURE Analysis
      ↓
Evaluate Multiple K Values
      ↓
Structure Harvester
      ↓
Visualize Admixture Patterns
      ↓
Principal Components Analysis (PCA)
      ↓
Interpret Population Structure
```

Your genotypes should already be formatted appropriately for STRUCTURE. See [discreteRFgenotypes_popmod.txt](https://github.com/squisquater/Tutorials/blob/main/01.PopStats-PopStr-msats/discreteRFgenotypes_popmod.txt) for an example input file.

To run these analyses on your own dataset, you will need to modify file names and parameter settings to match your project.

You can clone this repository by running:

```bash
git clone https://github.com/squisquater/Tutorials.git
```

**Note:** If you receive the error:

```bash
xcrun: error: invalid active developer path
```

you may need to install Xcode command line tools before cloning the repository.

---

# STRUCTURE

STRUCTURE is a model-based clustering method that attempts to assign individuals to genetic clusters based on allele frequencies.

The user specifies the number of clusters (**K**) to evaluate, and the program estimates ancestry proportions for each individual.

### Population Codes

STRUCTURE accepts numerical population identifiers rather than text labels.

| Population | Code |
| ---------- | ---- |
| WAC        | 1    |
| ORC        | 2    |
| LAS        | 3    |
| SN         | 4    |
| SV         | 5    |
| CANN       | 6    |

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
  * The Q-hat output file contains ancestry proportions for each individual and is useful for creating customized plots in R, ArcGIS, or Excel.
* Click **'OK'**
* Name the Parameter Set (i.e. 'WestRF_10k50k')

**STEP 3:** Run Structure
* Select **'Run'** from the **'Parameter Set'** dropdown menu
* Choose a value of **K**, the number of genetic clusters assumed by the model.
```
### What is K?

K represents the number of genetic clusters STRUCTURE attempts to identify.

Because the true number of clusters is typically unknown, analyses are usually run across a range of K values (e.g., K = 1–10).

Different values of K may reveal biologically meaningful structure at different hierarchical levels.

Recommended approach:
* Run multiple K values (e.g., 1–10)
* Run multiple replicates of each K value (e.g., 3-5)
* Compare likelihood scores across K values
* Evaluate biological plausibility alongside statistical results

```

**STEP 4:** View Results
* Select the value of K you want to view from the Results folder.
* Select **Bar Plot** from the Results menu to view the admixture proportions.
* This plot displays the ancestry proportions estimated for each individual.
* Individuals are represented by vertical bars, and colors represent ancestry assigned to different genetic clusters.
* You may sort individuals by population to improve visualization.You can sort by POP if they are not already in order.
* These results are also available in the Results folder and can be plotted in R or excel to customize the figure (i.e. colors, labels, etc)

### Interpreting STRUCTURE Plots

* Individuals represented by a single color are assigned primarily to one genetic cluster.
* Individuals containing multiple colors show evidence of mixed ancestry (admixture).
* Distinct color patterns among populations suggest genetic differentiation.
* Similar color patterns across populations suggest ongoing gene flow or limited differentiation.

The raw output files can also be exported and visualized using R, ArcGIS, or Excel.

**STEP 5:** Structure Harvester
* Zip your Results folder
* Navigate to the [Structure Harvester](https://www.researchgate.net/deref/https%3A%2F%2Flmme.ac.cn%2FStructureSelector%2F?_tp=eyJjb250ZXh0Ijp7ImZpcnN0UGFnZSI6InF1ZXN0aW9uIiwicGFnZSI6InF1ZXN0aW9uIn19) website
* Load the .zip file
* Click the **Harvest!** button
* Structure Harvester provides:

* Mean likelihood values for each K
* ΔK values using the Evanno method (assuming you ran >1 iteration of each K value)
* Summary plots for model comparison

**Important Note**

The K value with the highest ΔK is not necessarily the biologically correct answer.

Population structure should be interpreted using:

* STRUCTURE results
* PCA results
* Sampling design
* Knowledge of the species' biology and geography
---

# Take a Break and Learn How to Plot These Results in ArcGIS

---

# Principal Components Analysis (PCA)

PCA is a dimension-reduction technique that summarizes genetic variation into a small number of axes.

Unlike STRUCTURE, PCA does not require specifying K and makes relatively few assumptions about population structure.

Individuals that cluster together on a PCA plot are genetically similar.

---

## STEP 0: Install Required Packages

```r
install.packages("adegenet")
install.packages("hierfstat")
install.packages("reshape2")
install.packages("ggplot2")

library(adegenet)
library(hierfstat)
library(reshape2)
library(ggplot2)

setwd("/path/to/Tutorials/Summer-REU-2026/RedFox")
```

Replace the path above with the location where you cloned the repository.

---

## STEP 1: Read the STRUCTURE File

The file `project_data.stru` is generated by STRUCTURE and can be imported directly into R.

```r
genind_obj <- read.structure(
  "project_data.stru",
  n.ind = 301,
  n.loc = 31,
  onerowperind = TRUE,
  col.lab = 1,
  col.pop = 2,
  col.others = 0,
  row.marknames = 1,
  NA.char = "-9",
  ask = TRUE,
  quiet = TRUE
)
```

For additional details about these parameters:

```r
?read.structure
```

---

## STEP 2: Generate a PCA

```r
X <- scaleGen(genind_obj, NA.method = "mean")

pca_result <- dudi.pca(
  X,
  scannf = FALSE,
  nf = 3
)

s.class(
  pca_result$li,
  pop(genind_obj)
)
```

This will generate a PCA plot showing the genetic relationships among individuals.

**Interpreting PCA Results**

* Individuals located close together are genetically similar.
* Individuals located far apart are genetically differentiated.
* Distinct clusters suggest population structure.
* Overlapping clusters suggest ongoing gene flow or weak differentiation.

---

**STEP 3: Examine Variance Explained*

```r
eig_percent <- round(
  100 * pca_result$eig / sum(pca_result$eig),
  2
)

eig_percent
```

The first few principal components typically explain the largest proportion of genetic variation.

Reporting the percentage of variation explained by PC1 and PC2 is recommended when presenting PCA figures.

---

**STEP 4: Add Descriptive Population Labels**

For reference:

| Code | Population |
| ---- | ---------- |
| 1    | WAC        |
| 2    | ORC        |
| 3    | LAS        |
| 4    | SN         |
| 5    | SV         |
| 6    | CANN       |

```r
pop_labels <- c(
  "WAC",
  "ORC",
  "LAS",
  "SN",
  "SV",
  "CANN"
)

pop(genind_obj) <- factor(
  pop_labels[as.numeric(pop(genind_obj))],
  levels = pop_labels
)

table(pop(genind_obj))
```

Re-run the PCA:

```r
X <- scaleGen(genind_obj, NA.method = "mean")

pca_result <- dudi.pca(
  X,
  scannf = FALSE,
  nf = 3
)

s.class(
  pca_result$li,
  pop(genind_obj)
)
```

---

**Optional: Create a Publication-Quality PCA Figure**

```r
eig_percent <- round(
  100 * pca_result$eig / sum(pca_result$eig),
  2
)

pca_df <- data.frame(
  PC1 = pca_result$li[,1],
  PC2 = pca_result$li[,2],
  Population = pop(genind_obj)
)

ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = Population
  )
) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    x = paste0("PC1 (", eig_percent[1], "%)"),
    y = paste0("PC2 (", eig_percent[2], "%)"),
    color = "Population"
  )
```

This approach provides greater flexibility for customizing colors, labels, themes, and figure formatting.

---
