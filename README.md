# PADEconnectivity
Files for investigating larval summer flounder connectivity

## Data
* **masterPADElarvae.txt**: Data of 523 larval summer flounder, including date, location, 10 outlier loci, SL and developmental stage
* **masterPADElarvae_metadata.txt**: Metadata pertaining to masterPADElarve.txt

## Scripts
* **assignment_confusion_matrix.R**: Script to sample n individuals from each of 10 spawning regions, calculate genotype likelihoods of multilocus genotypes, determine most likely population and calculate % assignment and plot confusion matrices
* **conceptual_dispersal_diagram.R**: Script to plot **conceptual_dispersal_diagram.png**
* **GAMS_to_generate_likelihoods.R**: fits GAMS between allele frequencies and distance along coast, calculates observed individual genotype likelihoods across 10 populations, plots **indiv_assignments_10pops.png** and calculates observed cluster genotype likelihoods from individual likelihoods
* **group assignment.R**: Given larval clustering data and adult regional allele frequencies, simulates larval genotypes for clusters of size N, then calculates the distribution of genotype likelihoods for each of five adult populations. Plots observed cluster genotype likelihood on top of emperical distribution and calculates a p-value.
* **group assignment GAMS.R**: same as **group assignment.R** but using 10 populations determined from GAMs in **GAMS_to_generate_likelihoods.R**
* **population power analysis.R**: Exploring how to best do simulations for obtaining genotype likelihood distribution for each potential population of origin. This results in the **group assignment.R** script.
* **otoliths.R**: Using larval otolith microchemistry data for MDS, ANOVAs, LDAs, clustering, clustering pie charts, etc.
* **recentIBD.R**: Fst calculations for larvae captured 2008-2012

## Results
### Otoliths
* **cluster_piecharts_bytimeperiod.png**: pie charts of clusters in each time period, colored by ingress site
* **mds_byperiod_10elements.png**: nmds plots for each time period using all 10 elements
* **lda_barplots.png**: LDA barplots for training (data within 68% confidence ellipses), test & combined datasets
* **lda_barplots_3timeperiods.png**: LDA barplots for each time period separately
* **lda_barplots_late.png**: LDA barplots for training (data within 68% confidence ellipses), test & combined datasets for 2008-2012 larvae only
### Genetics
* **indiv_assignments_10pops.png**: barplots of individual assignment by ingress site
* **outlierloci_lineplot_GAM_10pops.png**: Allele frequencies at 10 loci across space
### Otoliths & genetics
* **conceptual_dispersal_diagram.png**: conceptual diagram depicting where larval clusters originated and where they ingressed to
### Assignment validation
* **assignment_confusion_matrix.png**: confusion matrix of individuals
* **assignment_confusion_matrix_5individuals.png**: confusion matrix of groups of 5 individuals
* **assignment_confusion_matrix_10individuals.png**: confusion matrix of groups of 10 individuals
* **assignment_confusion_matrix_20individuals.png**: confusion matrix of groups of 20 individuals
* **assignment_confusion_matrix_30individuals.png**: confusion matrix of groups of 30 individuals
