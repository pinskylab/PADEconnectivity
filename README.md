# PADEconnectivity
Using genetics and otolith microchemistry to investigate larval summer flounder connectivity throughout the species range from 1989-2012: *Using multiple natural tags provides evidence for extensive larval dispersal across space and through time in summer flounder* by [Hoey et al. (2020)](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.15414).

_These data, methods and scripts are provided in the interests of open science. If you have questions or find errors, please let me know._

Citable as [![DOI](https://zenodo.org/badge/129143601.svg)](https://zenodo.org/badge/latestdoi/129143601)

Contact:
Jennifer Hoey
(jahoey13@gmail.com)

## Data
* [Larvae Sampling Database.csv](https://github.com/pinskylab/PADEconnectivity/blob/master/Larvae%20Sampling%20Database.csv): csv of metadata for all sampled larvae, whether sequenced or not
* [masterPADElarvae.txt](https://github.com/pinskylab/PADEconnectivity/blob/master/masterPADElarvae.txt): Data of 523 larval summer flounder, including date, location, 10 outlier loci, SL and developmental stage
* [masterPADElarvae_metadata.txt](https://github.com/pinskylab/PADEconnectivity/blob/master/masterPADElarvae_metadata.txt): Metadata pertaining to masterPADElarve.txt
* [otolith_data.txt](https://github.com/pinskylab/PADEconnectivity/blob/master/otolith_data.txt): raw larval otolith microchemistry data
* [otolith_data_nozeros.txt](https://github.com/pinskylab/PADEconnectivity/blob/master/otolith_data_nozeros.txt): larval otolith microchemistry data where zeros have been randomly replaced by a number between zero and the smallest non-zero number for that element
* [pop.allele.counts10.txt](https://github.com/pinskylab/PADEconnectivity/blob/master/pop.allele.counts10.txt): Adult summer flounder allele counts at outlier loci from Hoey & Pinsky (2018)
* [SNP.DP3g95maf05lm75.FIL.recode.vcf](https://github.com/pinskylab/PADEconnectivity/blob/master/SNP.DP3g95maf05lm75.FIL.recode.vcf): all SNPs output following filtering described in [Bioinformatics.txt](https://github.com/pinskylab/PADEconnectivity/blob/master/Bioinformatics.txt)
* [SNP.DP3g95maf05lm75.FIL.recode.firstsnp.vcf](https://github.com/pinskylab/PADEconnectivity/blob/master/SNP.DP3g95maf05lm75.FIL.recode.firstsnp.vcf): same as [SNP.DP3g95maf05lm75.FIL.recode.vcf](https://github.com/pinskylab/PADEconnectivity/blob/master/SNP.DP3g95maf05lm75.FIL.recode.vcf) but only first SNP per contig containing 1904 loci across 293 larvae and 235 adult summer flounder
* [structure_input_Mar14_2017_528_1904.str](https://github.com/pinskylab/PADEconnectivity/blob/master/structure_input_Mar14_2017_528_1904.str): structure-formatted file of 1904 loci across 293 larvae and 235 adult summer flounder used for many of the genetic analyses

## Methods
* [Bioinformatics.txt](https://github.com/pinskylab/PADEconnectivity/blob/master/Bioinformatics.txt): Bioinformatic steps taken to go from sequencing reads to SNP genotypes, as well as steps taken to produce filtered SNP genotypes [SNP.DP3g95maf05lm75.FIL.recode.vcf](https://github.com/pinskylab/PADEconnectivity/blob/master/SNP.DP3g95maf05lm75.FIL.recode.vcf)

## Scripts
* [assignment_confusion_matrix.R](https://github.com/pinskylab/PADEconnectivity/blob/master/assignment_confusion_matrix.R): Script to sample n individuals from each of 10 spawning regions, calculate genotype likelihoods of multilocus genotypes, determine most likely population and calculate % assignment and plot confusion matrices
* [conceptual_dispersal_diagram.R](https://github.com/pinskylab/PADEconnectivity/blob/master/conceptual_dispersal_diagram.R): Script to plot [conceptual_dispersal_diagram.png](https://github.com/pinskylab/PADEconnectivity/blob/master/conceptual_dispersal_diagram.png)
* [GAMS_to_generate_likelihoods.R](https://github.com/pinskylab/PADEconnectivity/blob/master/GAMS_to_generate_likelihoods.R): fits GAMS between allele frequencies and distance along coast, calculates observed individual genotype likelihoods across 10 populations, plots [indiv_assignments_10pops.png](https://github.com/pinskylab/PADEconnectivity/blob/master/indiv_assignments_10pops.png) and calculates observed cluster genotype likelihoods from individual likelihoods
* [group assignment GAMS.R](https://github.com/pinskylab/PADEconnectivity/blob/master/group%20assignment%20GAMS.R): Given larval clustering data and GAM-determined adult allele frequencies calculated from [GAMS_to_generate_likelihoods.R](https://github.com/pinskylab/PADEconnectivity/blob/master/GAMS_to_generate_likelihoods.R), simulates larval genotypes for clusters of size N, then calculates the distribution of genotype likelihoods for each of ten adult populations. Plots observed cluster genotype likelihood on top of emperical distribution and calculates a p-value
* [otoliths.R](https://github.com/pinskylab/PADEconnectivity/blob/master/otoliths.R): Using larval otolith microchemistry data for NMDS, ANOVAs, LDAs, clustering, clustering pie charts, etc.
* [plot allele freqs.R](https://github.com/pinskylab/PADEconnectivity/blob/master/plot%20allele%20freqs.R): Plots allele frequencies of 10 GAM-determined locations to produce [outlierloci_lineplot_GAM_10pops.png](https://github.com/pinskylab/PADEconnectivity/blob/master/outlierloci_lineplot_GAM_10pops.png)
* [recentIBD.R](https://github.com/pinskylab/PADEconnectivity/blob/master/recentIBD.R): Fst calculations for larvae captured 2008-2012
* [structurethrutime.R](https://github.com/pinskylab/PADEconnectivity/blob/master/structurethrutime.R): bulk of larval genomic analyses, including Fst, AMOVA
* [larvalstructure.R](https://github.com/pinskylab/PADEconnectivity/blob/master/larvalstructure.R): spatial and temporal PCAs, generating [larval_pca_regionalpops_1904loci.png](https://github.com/pinskylab/PADEconnectivity/blob/master/larval_pca_regionalpops_1904loci.png) and [larval_pca_temppops_1904loci.png](https://github.com/pinskylab/PADEconnectivity/blob/master/larval_pca_temppops_1904loci.png)

## Results
### Otoliths
* [cluster197_piecharts_bytimeperiod.png](https://github.com/pinskylab/PADEconnectivity/blob/master/cluster197_piecharts_bytimeperiod.png): pie charts of larval clusters with otolith-only data for each time period (_Figure S8_)
* [nmds_byperiod_10elements.png](https://github.com/pinskylab/PADEconnectivity/blob/master/nmds_byperiod_10elements.png): nmds plots for each time period using all 10 elements (_Figure 3_)
* [lda_barplots.png](https://github.com/pinskylab/PADEconnectivity/blob/master/lda_barplots.png): LDA barplots for training (data within 68% confidence ellipses), test & combined datasets (_Figure S9_)
* [lda_barplots_3timeperiods.png](https://github.com/pinskylab/PADEconnectivity/blob/master/lda_barplots_3timeperiods.png): LDA barplots for each time period separately (_Figure S11_)
* [lda_barplots_late.png](https://github.com/pinskylab/PADEconnectivity/blob/master/lda_barplots_late.png): LDA barplots for training (data within 68% confidence ellipses), test & combined datasets for 2008-2012 larvae only
* [otolith_overtimelogged_boxplots.png](https://github.com/pinskylab/PADEconnectivity/blob/master/otolith_overtimelogged_boxplots.png): box plots of otolith microchemistry by time period (_Figure S7_)
* [otolith_sitelogged_boxplots.png](https://github.com/pinskylab/PADEconnectivity/blob/master/otolith_sitelogged_boxplots.png): box plots of otolith microchemistry by ingress site (_Figure S6_)
* [post_group_membership_prob.png](https://github.com/pinskylab/PADEconnectivity/blob/master/post_group_membership_prob.png): frequency distributions of predicted group membership probabilities from LDA classification results (_Figure S10_)
### Genetics
* [larvae293_1904loci.pdf](https://github.com/pinskylab/PADEconnectivity/blob/master/larvae293_1904loci.pdf): STRUCTURE plot of 293 larvae using 1904 loci (_Figure S3_)
* [larval_pca_regionalpops_1904loci.png](https://github.com/pinskylab/PADEconnectivity/blob/master/larval_pca_regionalpops_1904loci.png): regional larval PCA (_Figure S2A_)
* [larval_pca_temppops_1904loci.png](https://github.com/pinskylab/PADEconnectivity/blob/master/larval_pca_temppops_1904loci.png): temporal larval PCA (_Figure S2B_)
* [indiv_assignments_10pops.png](https://github.com/pinskylab/PADEconnectivity/blob/master/indiv_assignments_10pops.png): barplots of individual assignment by ingress site (_Figure S5_)
* [outlierloci_lineplot_GAM_10pops.png](https://github.com/pinskylab/PADEconnectivity/blob/master/outlierloci_lineplot_GAM_10pops.png): Adult allele frequencies at 10 outlier loci across 10 GAM-determined locations (_Figure S1_)
### Otoliths & genetics
* [cluster_piecharts_bytimeperiod.png](https://github.com/pinskylab/PADEconnectivity/blob/master/cluster_piecharts_bytimeperiod.png): pie charts of clusters in each time period, colored by ingress site (_Figure S12_)
* [obs_likelihoods_early_6clusters_10pops.txt](https://github.com/pinskylab/PADEconnectivity/blob/master/obs_likelihoods_early_6clusters_10pops.txt): Likelihoods of larval clusters from the early time period (1989-1993) across all 10 GAM-determined adult spawning locations, with Pop 10 being northern-most and Pop 1 being southern-most
* [obs_likelihoods_middle_2clusters_10pops.txt](https://github.com/pinskylab/PADEconnectivity/blob/master/obs_likelihoods_middle_2clusters_10pops.txt): Likelihoods of larval clusters from the middle time period (1998-2002) across all 10 GAM-determined adult spawning locations, with Pop 10 being northern-most and Pop 1 being southern-most
* [obs_likelihoods_late_3clusters_10pops.txt](https://github.com/pinskylab/PADEconnectivity/blob/master/obs_likelihoods_late_3clusters_10pops.txt): Likelihoods of larval clusters from the late time period (2008-2012) across all 10 GAM-determined adult spawning locations, with Pop 10 being northern-most and Pop 1 being southern-most
* [conceptual_dispersal_diagram.png](https://github.com/pinskylab/PADEconnectivity/blob/master/conceptual_dispersal_diagram.png): conceptual diagram depicting where larval clusters originated and where they ingressed to (_Figure 4_)
* [indiv_and_cluster_assignments_hists_10GAMpops.png](https://github.com/pinskylab/PADEconnectivity/blob/master/indiv_and_cluster_assignments_hists_10GAMpops.png): likelihood differences at the individual versus cluster levels (_Figure S4_)
* [cluster_likelihoods_10GAMSgroups.png](https://github.com/pinskylab/PADEconnectivity/blob/master/cluster_likelihoods_10GAMSgroups.png): likelihoods of each cluster plotted across 10 GAM-determined locations for each time period
### Assignment validation
_These make up Figure 2A-E_
* [assignment_confusion_matrix.png](https://github.com/pinskylab/PADEconnectivity/blob/master/assignment_confusion_matrix.png): confusion matrix of individuals
* [assignment_confusion_matrix_5individuals.png](https://github.com/pinskylab/PADEconnectivity/blob/master/assignment_confusion_matrix_5indiv.png): confusion matrix of groups of 5 individuals
* [assignment_confusion_matrix_10individuals.png](https://github.com/pinskylab/PADEconnectivity/blob/master/assignment_confusion_matrix_10indiv.png): confusion matrix of groups of 10 individuals
* [assignment_confusion_matrix_20individuals.png](https://github.com/pinskylab/PADEconnectivity/blob/master/assignment_confusion_matrix_20indiv.png): confusion matrix of groups of 20 individuals
* [assignment_confusion_matrix_30individuals.png](https://github.com/pinskylab/PADEconnectivity/blob/master/assignment_confusion_matrix_30indiv.png): confusion matrix of groups of 30 individuals
