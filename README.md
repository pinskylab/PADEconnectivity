# PADEconnectivity
Files for investigating larval summer flounder connectivity

## Data
* **masterPADElarvae.txt**: Data of 523 larval summer flounder, including date, location, 10 outlier loci, SL and developmental stage
* **masterPADElarvae_metadata.txt**: Metadata pertaining to masterPADElarve.txt

## Scripts
* **otoliths.R**: Using larval otolith microchemistry data for MDS, ANOVAs, LDAs, clustering, clustering pie charts, etc.
* **group assignment.R**: Given larval clustering data and adult regional allele frequencies, simulates larval genotypes for clusters of size N, then calculates the distribution of genotype likelihoods for each of five adult populations. Plots observed cluster genotype likelihood on top of emperical distribution and calculates a p-value.
* **population power analysis.R**: Exploring how to best do simulations for obtaining genotype likelihood distribution for each potential population of origin. This results in the **group assignment.R** script.
* **recentIBD.R**: Fst calculations for larvae captured 2008-2012

## Results
* **lda_barplots.R**: LDA barplots for training (data within 68% confidence ellipses), test & combined datasets
* **lda_barplots_3timeperiods.R**: LDA barplots for each time period separately
