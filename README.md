# PADEconnectivity
Files for investigating larval summer flounder connectivity throughout the species range from 1989-2012

_These data, methods and scripts are provided in the interests of open science. If you have questions or find errors, please let me know._

Contact:
Jennifer Hoey
jahoey13@gmail.com

## Data
* **Larvae Sampling Database.csv**: csv of metadata for all sampled larvae, whether sequenced or not
* **masterPADElarvae.txt**: Data of 523 larval summer flounder, including date, location, 10 outlier loci, SL and developmental stage
* **masterPADElarvae_metadata.txt**: Metadata pertaining to masterPADElarve.txt
* **otolith_data.txt**: raw larval otolith microchemistry data
* **otolith_data_nozeros.txt**: larval otolith microchemistry data where zeros have been randomly replaced by a number between zero and the smallest non-zero number for that element
* **pop.allele.counts10.txt**: Adult summer flounder allele counts at outlier loci from Hoey & Pinsky (2018)
* **structure_input_Mar14_2017_528_1904.str**: structure-formatted file of 1904 loci across 293 larvae and 235 adult summer flounder used for many of the genetic analyses

## Methods
* **Bioinformatics.txt**: Bioinformatic steps taken to produce filtered SNP genotypes

## Scripts
* **assignment_confusion_matrix.R**: Script to sample n individuals from each of 10 spawning regions, calculate genotype likelihoods of multilocus genotypes, determine most likely population and calculate % assignment and plot confusion matrices
* **conceptual_dispersal_diagram.R**: Script to plot **conceptual_dispersal_diagram.png**
* **GAMS_to_generate_likelihoods.R**: fits GAMS between allele frequencies and distance along coast, calculates observed individual genotype likelihoods across 10 populations, plots **indiv_assignments_10pops.png** and calculates observed cluster genotype likelihoods from individual likelihoods
* **group assignment GAMS.R**: Given larval clustering data and GAM-determined adult allele frequencies calculated from **GAMS_to_generate_likelihoods.R**, simulates larval genotypes for clusters of size N, then calculates the distribution of genotype likelihoods for each of ten adult populations. Plots observed cluster genotype likelihood on top of emperical distribution and calculates a p-value
* **otoliths.R**: Using larval otolith microchemistry data for NMDS, ANOVAs, LDAs, clustering, clustering pie charts, etc.
* **plot allele freqs.R**: Plots allele frequencies of 10 GAM-determined locations to produce **outlierloci_lineplot_GAM_10pops.png**
* **recentIBD.R**: Fst calculations for larvae captured 2008-2012
* **structurethrutime.R**: bulk of larval genomic analyses, including Fst, AMOVA
* **larvalstructure.R**: spatial and temporal PCAs, generating **larval_pca_regionalpops_1904loci.png** and **larval_pca_temppops_1904loci.png**

## Results
### Otoliths
* **cluster197_piecharts_bytimeperiod.png**: pie charts of larval clusters with otolith-only data for each time period (_Figure S8_)
* **nmds_byperiod_10elements.png**: nmds plots for each time period using all 10 elements (_Figure 3_)
* **lda_barplots.png**: LDA barplots for training (data within 68% confidence ellipses), test & combined datasets (_Figure 4_)
* **lda_barplots_3timeperiods.png**: LDA barplots for each time period separately (_Figure S9_)
* **lda_barplots_late.png**: LDA barplots for training (data within 68% confidence ellipses), test & combined datasets for 2008-2012 larvae only
### Genetics
* **larval_pca_regionalpops_1904loci.png**: regional larval PCA (_Figure S2A_)
* **larval_pca_temppops_1904loci.png**: temporal larval PCA (_Figure S2B_)
* **indiv_assignments_10pops.png**: barplots of individual assignment by ingress site (_Figure S5_)
* **outlierloci_lineplot_GAM_10pops.png**: Adult allele frequencies at 10 outlier loci across 10 GAM-determined locations (_Figure S1_)
### Otoliths & genetics
* **cluster_piecharts_bytimeperiod.png**: pie charts of clusters in each time period, colored by ingress site (_Figure S10_)
* **obs_likelihoods_early_6clusters_10pops.txt**: Likelihoods of larval clusters from the early time period (1989-1993) across all 10 GAM-determined adult spawning locations, with Pop 10 being northern-most and Pop 1 being southern-most
* **obs_likelihoods_middle_2clusters_10pops.txt**: Likelihoods of larval clusters from the middle time period (1998-2002) across all 10 GAM-determined adult spawning locations, with Pop 10 being northern-most and Pop 1 being southern-most
* **obs_likelihoods_late_3clusters_10pops.txt**: Likelihoods of larval clusters from the late time period (2008-2012) across all 10 GAM-determined adult spawning locations, with Pop 10 being northern-most and Pop 1 being southern-most
* **conceptual_dispersal_diagram.png**: conceptual diagram depicting where larval clusters originated and where they ingressed to (_Figure 5_)
* **indiv_and_cluster_assignments_hists_10GAMpops.png**: likelihood differences at the individual versus cluster levels (_Figure S4_)
* **cluster_likelihoods_10GAMSgroups.png**: likelihoods of each cluster plotted across 10 GAM-determined locations for each time period
### Assignment validation
_These make up Figure 2A-E_
* **assignment_confusion_matrix.png**: confusion matrix of individuals
* **assignment_confusion_matrix_5individuals.png**: confusion matrix of groups of 5 individuals
* **assignment_confusion_matrix_10individuals.png**: confusion matrix of groups of 10 individuals
* **assignment_confusion_matrix_20individuals.png**: confusion matrix of groups of 20 individuals
* **assignment_confusion_matrix_30individuals.png**: confusion matrix of groups of 30 individuals
