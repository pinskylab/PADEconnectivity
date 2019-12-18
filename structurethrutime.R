setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files")

library(ade4)
library(adegenet)
library(pegas)
library(hierfstat)
library(mmod)
library(poppr)
library(diveRsity)

# Reading in SNP data file containing only the first SNP at each locus
# full <- read.structure("structure_input_Sept7_2016.str",
#                          n.ind = 568, n.loc = 1651, col.lab = 1, col.pop = 2, row.marknames = 1,
#                          onerowperind = FALSE)

# full <- read.structure("structure_input_Feb10_2017_temppops.str",
#                        n.ind = 568, n.loc = 1651, col.lab = 1, col.pop = 2, row.marknames = 1,
#                        onerowperind = FALSE)

# full <- read.structure("structure_input_March12_2017_554_1109.str",
#                         n.ind = 554, n.loc = 1109, col.lab = 1, col.pop = 2, row.marknames = 1,
#                         onerowperind = FALSE) # too much missing data

# full <- read.structure("structure_input_Mar13_2017_528_3827.str",
#                        n.ind = 528, n.loc = 3827, col.lab = 1, col.pop = 2, row.marknames = 1,
#                        onerowperind = FALSE)

full <- read.structure("structure_input_Mar14_2017_528_1904.str",
                       n.ind = 528, n.loc = 1904, col.lab = 1, col.pop = 2, row.marknames = 1,
                       onerowperind = FALSE)

full.dataframe <- as.data.frame(full@tab)

# adults <- rownames(full.dataframe)[205:439] #528 fish total
# larvs <- full.dataframe[!(rownames(full.dataframe) %in% adults),]
# dim(larvs) #293 x 3831, 300 x 2229
# larvs <- as.genind(larvs)
# larvs@pop <- full@pop[which(full@pop != 4)]

###############################################################
#### How much missing data per locus in larvae and adults? ####
na.count <- sapply(full.dataframe, function(y) sum(length(which(is.na(y)))))
na.count <- data.frame(na.count)
head(na.count)

sum(is.na(full.dataframe[,2]))

# Histogram of SNPs with missing data
counts <- data.matrix(na.count)
hist(counts)

# Max number of missing individuals for a locus is 26. 26/528 = 5%

#### How much missing data per individual in larvae and adults? ####
loci.na.count <- rowSums(is.na(full.dataframe))
loci.na.count <- data.frame(loci.na.count)
head(loci.na.count)

sum(is.na(full.dataframe[6,]))

# Histogram of how much missing data within an individual
counts.byindiv <- data.matrix(loci.na.count)
hist(counts.byindiv)

# Max number of alleles missing within an individual is 1289. 1289/3831 = 34%
##############################
#### How much missing data per locus in larvae only? ####
na.count.larvs <- sapply(larvs, function(y) sum(length(which(is.na(y)))))
na.count.larvs <- data.frame(na.count.larvs)
head(na.count.larvs)

sum(is.na(larvs[,3]))

# Histogram of SNPs with missing data
counts.larvs <- data.matrix(na.count.larvs)
hist(counts.larvs)

# Max number of missing individuals for a locus is 22. 22/293 = 7.5%
# Max number of missing individuals for a locus is 26 using structure_input_March12_2017_554_1109.str. 26/300 = 8.6%

#### How much missing data per individual in larvae? ####
larvs.loci.na.count <- rowSums(is.na(larvs))
larvs.loci.na.count <- data.frame(larvs.loci.na.count)
head(larvs.loci.na.count)

sum(is.na(larvs[6,]))

# Histogram of how much missing data within an individual
counts.byindiv.larvs <- data.matrix(larvs.loci.na.count)
hist(counts.byindiv.larvs)
# Max number of alleles missing within an individual is 1289. 1289/3831 = 34%
# Max number of alleles missing within an individual is 1672 using structure_input_March12_2017_554_1109.str. 1672/2229 = 75%, this is too much

#######################################
#### Working with the full dataset ####
is.genind(full)
head(indNames(full),10)
locNames(full)
sum <- summary(full)
plot(sum$Hobs ~ sum$Hexp)

# # I want to plot a histogram of the number of fish in each year/How even are the years?
# ind_528 <- data.frame(matrix(nrow=528,ncol=2))
# names(ind_528) <- c("ID", "Year")
# ind_528$ID <- indNames(full)
# write.table(ind_528, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/ind_528.txt", sep="\t")
# ind_528 <- read.table("ind_528.txt")
# ind_528
# 
# # Plotting histogram of number of samples/year
# hist(ind_528$Year, breaks=337)
# table(ind_528$Year)

#### Dealing with loci with multiple alleles ####
multiallelic <- names(which(full@loc.n.all > 2))
multiallelic <- paste(multiallelic, ".", sep = "")

multiallelic_data <- full@tab[, grep(paste(multiallelic, "{3,4}", sep = "", collapse = "|"), colnames(full@tab))] # mostly works, but still includes SNP_145*
multiallelic_data <- multiallelic_data[,-c(59:78)] # gets rid of additional SNP_145* that are actually biallelic but got picked up in the grep above
multiallelic_counts <- colSums(multiallelic_data, na.rm = TRUE)

multiallelic_names <- names(multiallelic_counts) #67
biallelic_only <- full.dataframe[,!colnames(full.dataframe) %in% multiallelic_names] # 528x3764 (before: 3831-67 = 3764)

# Separate genind for adults and larvae
adults <- rownames(biallelic_only)[205:439] #528 fish total
larvs <- biallelic_only[!(rownames(biallelic_only) %in% adults),]
# write.table(larvs, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/larvs_allelecounts.txt", sep="\t")
adult_data <- biallelic_only[(rownames(biallelic_only) %in% adults),]
dim(larvs) #293 x 3764
larvs <- as.genind(larvs)
larvs@pop <- full@pop[which(full@pop != 4)]
larvs.sum <- summary(larvs)
plot(larvs.sum$Hobs ~ larvs.sum$Hexp)
adult_data <- as.genind(adult_data)
adult_data@pop <- full@pop[which(full@pop == 4)]

# # For adults & larvae. Turn it back into a genind
# full_bi <- as.genind(biallelic_only)
# full_bi@pop <- full@pop


########################################################################################################
# The rest of the script explores PCA (line 149), summary statistics (line 244), DAPC (line 400),      #
# combining larval allele frequencies with database data (line 416) for use in AMOVA (line 500) and    #
# larval allele frequency changes over time (line 590)                                                 #
########################################################################################################

#### PCA ####
sum(is.na(larvs$tab)) #9424 
X <- scaleGen(larvs, NA.method = "mean")
dim(X)
class (X)

# make PCA
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

pca1

# Plotting PC1 and PC2
s.label(pca1$li)
title("PCA of summer flounder dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

s.class(pca1$li, pop(larvs))
title("PCA of summer flounder dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

### To make a nice plot of the PCA broken down by time period ###
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/results/larval_pca_temppops_b&w.png", width=8, height=7, res=300, units="in")

par(
  mar=c(9, 8, 3, 8), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA,
  bty = 'n'
)

# col <- azur(3)
col <- greypal(4)
col <- c("#C6C6C6", "#686868", "#000000")
s.class(pca1$li, pop(larvs), xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, xlim = c(-50,30), ylim = c(-30,28), clabel = 0)
axis(1, at=seq(-40,20, by=10), labels=seq(-40,20, by= 10), line = 1.5)
axis(2, at=seq(-20,20, by = 10), labels=seq(-20,20, by= 10), line = 3, las = 2)
mtext("PC1 (0.67%)", side = 1, line = 4)
mtext("PC2 (0.65%)", side = 2, line = 5.5)
# text(-34, 17, "B")

legend(-36, -7,
       legend=c("1989-1993 (n = 60)", "1998-2002 (n = 87)", "2008-2012 (n = 146)"),
       pch=c(19, 19, 19),
       # col=col,
       col = c("#C6C6C6", "#686868", "#000000"),
       bty = "n",
       y.intersp = 1)

legend(-40.7,21,
       legend = expression("B"),
       bty = "n",
       cex = 1.3)
dev.off()

eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent [1:3]

#### Closer look at the monomorphic loci in larvae ####
biallelic_only[,2185:2186] #SNP_1111
biallelic_only[,3405:3406] #SNP_1723

biallelic_only[,c(2185,2186,3405,3406)]


#### Break PCA down by year/period, plotting PC1 vs PC2 528 fish and 1882 alleles####
plot(pca1$li[205:439,1], pca1$li[205:439,2], col = "blue", xlab = "PC1 (1.06%)", ylab = "PC2 (0.58%)", xlim = c(-50,10), ylim = c(-30,28)) # plots 2013-2014 adults
points(pca1$li[440:499,1], pca1$li[440:499,2], col = "tomato") # plot 1989-1993
points(pca1$li[500:528,1], pca1$li[500:528,2], col = "gold") # plot 1998-1999
points(pca1$li[1:58,1], pca1$li[1:58,2], col = "gold") # plot 2000-2002
points(pca1$li[59:204,1], pca1$li[59:204,2], col = "green") # plot 2008-2012

legend("topleft",
       legend=c("Larvae (1989-1993)", "Larvae (1998-2002)", "Larvae (2008-2012)", "Adults (2013-2014)"),
       pch=c(1, 1, 1, 1, 1),
       col=c("tomato", "gold", "green", "blue"))

# Plotting PC1 vs. PC3
plot(pca1$li[205:439,1], pca1$li[205:439,3], col = "blue", xlab = "PC1 (1.06%)", ylab = "PC3 (0.51%)", xlim = c(-50,10), ylim = c(-50,20)) # plots 2013-2014 adults
points(pca1$li[440:499,1], pca1$li[440:499,3], col = "tomato") # plot 1989-1993
points(pca1$li[500:528,1], pca1$li[500:528,3], col = "gold") # plot 1998-1999
points(pca1$li[1:58,1], pca1$li[1:58,3], col = "gold") # plot 2000-2002
points(pca1$li[59:204,1], pca1$li[59:204,3], col = "green") # plot 2008-2012

legend("bottomleft",
       legend=c("Larvae (1989-1993)", "Larvae (1998-2002)", "Larvae (2008-2012)", "Adults (2013-2014)"),
       pch=c(1, 1, 1, 1, 1),
       col=c("tomato", "gold", "green", "blue"))

#### PCA of larvae broken down by region ####
# See larvalstructure.R

#### Calculate FST statistics ####
# Read in a text file containing only larvae but multiallelic loci
larvs293_multiallelic <- read.table("structure_input_Mar14_2017_293_1904_multiallelic.txt", sep="\t", header = TRUE)

# Get rid of multiallelic loci
# Multiallelic loci list is in multiallelic
multiallelic <- names(which(full@loc.n.all > 2))
larvs293_biallelic <- larvs293_multiallelic[,!(colnames(larvs293_multiallelic) %in% multiallelic)]

# Split data into odd and even row dataframes
# Multiallelic calculation followed by biallelic calculation
even_indexes<-seq(2,586,2)
odd_indexes<-seq(1,585,2)

# Multiallelic
odds <- data.frame(larvs293_multiallelic[odd_indexes,]) # 293 x 1906
odds2 <- odds[,-c(1:2)] # 293 x 1904
evens <- data.frame(larvs293_multiallelic[even_indexes,]) # 293 x 1906
evens2 <- evens[,-c(1:2)] # 293 x 1904

# Biallelic
odds <- data.frame(larvs293_biallelic[odd_indexes,]) # 293 x 1884
odds2 <- odds[,-c(1:2)] # 293 x 1882
evens <- data.frame(larvs293_biallelic[even_indexes,]) # 293 x 1884
evens2 <- evens[,-c(1:2)] # 293 x 1882

# Now paste each value in one dataframe with its corresponding in the other
# Multiallelic
s <- 1:length(colnames(evens2))
combo <- data.frame(matrix(nrow = 293, ncol = 1904))
for (i in s){
  combo[,i] <-paste(odds2[,i], evens2[,i], sep = '')
}

dim(combo) # 293 x 1904

combo[] <- lapply(combo, function(x) as.numeric(as.character(x)))# Convert to numeric, gives warning because replaces character 'NANA' with NA

pop.names <- as.numeric(as.character(evens$Pop))

# Combine the population numbers with the allele data
combo2 <- cbind(pop.names, combo)
dim(combo2) # 293 x 1905

pairwise.WCfst(combo2,diploid=TRUE) 
genet.dist(combo2, method = 'WC84')

# Biallelic
s <- 1:length(colnames(evens2))
combo <- data.frame(matrix(nrow = 293, ncol = 1882))
for (i in s){
  combo[,i] <-paste(odds2[,i], evens2[,i], sep = '')
}

dim(combo) # 293 x 1882

combo[] <- lapply(combo, function(x) as.numeric(as.character(x)))# Convert to numeric, gives warning because replaces character 'NANA' with NA

pop.names <- as.numeric(as.character(evens$Pop))

# Combine the population numbers with the allele data
combo2 <- cbind(pop.names, combo)
dim(combo2) # 293 x 1883

pairwise.WCfst(combo2,diploid=TRUE) 
genet.dist(combo2, method = 'WC84')

#### FST calculation with each population representing a site/time period combination ####
larvs <- read.table("structure_input_Jan9_2019_293_1904_9pops.txt", sep="\t", header = TRUE)

even_indexes<-seq(2,586,2)
odd_indexes<-seq(1,585,2)

odds <- data.frame(larvs[odd_indexes,]) # 293 x 1906
odds2 <- odds[,-c(1:2)] # 293 x 1904
evens <- data.frame(larvs[even_indexes,]) # 293 x 1906
evens2 <- evens[,-c(1:2)] # 293 x 1904

s <- 1:length(colnames(evens2))
combo <- data.frame(matrix(nrow = 293, ncol = 1904))
for (i in s){
  combo[,i] <-paste(odds2[,i], evens2[,i], sep = '')
}

dim(combo) # 293 x 1904

combo[] <- lapply(combo, function(x) as.numeric(as.character(x)))# Convert to numeric, gives warning because replaces character 'NANA' with NA

pop.names <- as.numeric(as.character(evens$Pop))

# Combine the population numbers with the allele data
combo2 <- cbind(pop.names, combo)
dim(combo2) # 293 x 1905

pairwise.WCfst(combo2,diploid=TRUE) 
genet.dist(combo2, method = 'WC84')

options("scipen"=100, "digits"=4) # forces output to not be in scientific notation


# #####################################################
# 
# # Pairwise population Fst statistics
# # Calculates Nei's pairwise Fst, apparently this overestimates compared to Bhatia et al (2013)
# pairwise.fst(larvs, pop = larvs@pop) # Nei's pairwise Fst (is this the same thing as Nei's Gst?) - this works! Do I want this or WC Fst?
# 
# # Calculates Weir & Cockerham's Fst
# data("gtrunchier")
# pairwise.WCfst(gtrunchier[,-2],diploid=TRUE) # example using a test dataset
# genet.dist(gtrunchier[,-2], method = 'WC84')
# 
# # Need to get the genind object into the correct format first
# loci <- names(larvs@all.names)
# 
# odds <- data.frame() # odd indices
# for(i in 1:3764){
#   if (i %% 2 ==1) 
#    odds <- rbind(odds, i)
# } 
# 
# evens <- data.frame() # even indices
# for(i in 1:3764){
#   if (i %% 2 ==0) 
#     evens <- rbind(evens, i)
# }
# 
# larvs2 <- as.numeric(as.character(larvs@pop))
# 
# j <- 1:1882
# 
# larvs.loci <- data.frame(matrix(nrow = 293, ncol = 1882)) # Create empty dataframe so that alleles for an individual are in a single column
#   
# for (j in 1:length(loci)){
#       larvs.loci[,j] <- paste(pb[,odds[j,]], pb[,evens[j,]], sep="")
# } # This for loop fills in the dataframe by joining adjacent columns. This only works because all SNPs are biallelic.
# colnames(larvs.loci) <- loci # Fill in column names
# larvs2 <- cbind(larvs2,larvs.loci) # combine locus data with individual names and populations
# rownames(larvs2) <- rownames(larvs@tab) # Name individuals
# # These are allele counts. Need a way to indicate different alleles in an individual: 20  --> 11, 11 --> 12, 02 --> 22
# larvs2[larvs2 == "11"] <- 12
# larvs2[larvs2 == "20"] <- 11
# larvs2[larvs2 == "02"] <- 22
# 
# larvs2[] <- lapply(larvs2, function(x) as.numeric(as.character(x)))# Convert to numeric, gives warning because replaces character 'NANA' with NA
# 
# pairwise.WCfst(larvs2,diploid=TRUE)
# genet.dist(larvs2, method = 'WC84') # gives same answer as above, smaller than Nei's Fst estimates
# 
# # Per locus estimates of Fst (Nei's Fst)
# bs <- basic.stats(larvs2)
# boxplot(bs$perloc$Fst)
# which(bs$perloc$Fst > 0.017)
# 
# # W&C Fst is an estimator of Wright's Fst (takes sample size into account)
# fastDivPart(larvs, fst = TRUE) #input files needs to be genepop and this will take some work

#### DAPC ####
# Don't think DAPC is actually that helpful because the PCs explain so little of the variation
# K-means
grp <- find.clusters(larvs, max.n.clust = 10) # no clear asymptote, so keep all PCs, 2
table.value(table(pop(larvs), grp$grp), col.lab=paste("grp", 1:2))
dapc1 <- dapc(larvs, grp$grp)
dapc1

scatter(dapc1)
test <- optim.a.score(dapc1) #80,1 --> 14

# Choosing the number of PCs to retain
dapc2 <- dapc(larvs, n.da = 1, n.pca = 14)
scatter(dapc2)


#### This bit builds a dataframe with hierarchical population strata for use in AMOVA and looking at changes in allele frequencies over time ####
# Differences between larvae?
dim(larvs@tab)

# Need to be able to categorize by region, so reading in larval database so that I can match larvae names to get location
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Larvae")
larvae_database <- read.csv("Larvae Sampling Database.csv", header=TRUE)

allelefreqs_larvs <- as.data.frame(larvs@tab) # Convert larval genind object to a dataframe

# The names are all messed up because of the bioinformatics process...
names1 <- do.call(rbind, strsplit(as.character(larvae_database$ID..), '_'))
names2 <- do.call(rbind, strsplit(as.character(names1[,1]), 'E'))
correctnames <- paste(names2[,1], names2[,2], sep = 'E_')
correctnames2 <- paste(correctnames, names1[,2], sep = '')

# Now put back the fixed names in the larval database and the split names into the allele frequency database
larvae_database$ID.. <- correctnames2

# Need to modify names in allele count data (genind object) so that I can match them up with names in the larval database
freq_names <- as.vector(rownames(larvs@tab))
freq_names_split <- do.call(rbind, strsplit(as.character(freq_names), 'L'))
rownames(allelefreqs_larvs) <- freq_names_split[,1] # replaces rownames in dataframe allele frequency counts with modified format
rownames(larvs@tab) <- freq_names_split[,1] # replaces rownames in genind object with different formatting

# Subset the larval database to only larvae that are staying in the analysis, and order the names so that they match with the allele frequency matrix
larvae_database_sub <- larvae_database[larvae_database$ID.. %in% freq_names_split[,1],]
ordered_larvae_database_sub <- larvae_database_sub[order(larvae_database_sub$ID..),]

# Make sure name order is good
ordered_larvae_database_sub$ID.. == rownames(allelefreqs_larvs)

time_period <- as.vector(larvs@pop) # Converts population numbers to time periods
time_period <- gsub('2', 'half', time_period)
time_period <- gsub('1', 'early', time_period)
time_period <- gsub('3', 'late', time_period)

geo <- ordered_larvae_database_sub[, "Place"] # Subset the ordered larval database and assign region based on place
geo <- gsub('Beaufort, NC', 'south', geo)
geo <- gsub('North Inlet, SC', 'south', geo)
geo <- gsub('Little Egg Inlet, NJ', 'north', geo)
geo <- gsub('York River, VA', 'north', geo)
geo <- gsub('Roosevelt Inlet, DE', 'north', geo)

pop_strata <- data.frame(cbind(rownames(allelefreqs_larvs), time_period, geo))
strata(larvs) <- pop_strata[,-1]

# Hierarchical AMOVA with location and time period. Time period needs to be coded differently for 'north' & 'early' vs. 'south' & 'early', for example
# This just needs to be done once to get the codes for ingress site/time period combos. Can then proceed with reading in pop.id2.
pop.id <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/293larvs_popIDs_forstructure.txt")

# Read in file with population assignment that denotes ingress site/time period combination
pops <- read.structure("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/structure_input_528_10outliers_10pops2.str", n.ind = 528, n.loc = 10, col.lab = 1, col.pop = 2, row.marknames = 1,
                       onerowperind = FALSE)

# Isolate unique larval IDs
larvs.uniq <- as.character(unique(pop.id$V1))

# Isolate pop identifier & cbind with ID
pop <- cbind(rownames(pops@tab), as.vector(pops@pop)) # this also has adults

# Get rid of adults
pop.larvs <- pop[pop[,1] %in% larvs.uniq,] #293 x2
# write.table(pop.larvs, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/293larvs_popIDs.txt", sep = "\t", row.names = FALSE, col.names = FALSE) # Get rid of "'s and quickly fix location of _ in names

# Read pop IDs/hierarchical levels back in and merge with location and time period data, then run AMOVA
pop.id2 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/293larvs_popIDs.txt", header = TRUE)

pop.id3 <- merge(ordered_larvae_database_sub, pop.id2, by.x = 'ID..', by.y = 'Fish_ID')
pop.id3 <- droplevels(pop.id3) # gets rid of unused factor levels
pop.id3$Locations_all <- as.factor(pop.id3$Locations_all)
pop.id3$Time_within_NSLocations <- as.factor(pop.id3$Time_within_NSLocations)
pop.id3$Time <- as.factor(pop.id3$Time)
pop.id3$NSLocations_within_Time <- as.factor(pop.id3$NSLocations_within_Time)

pop.id3$Place <- gsub('Beaufort, NC', 'south', pop.id3$Place) # instead of each ingress site separately, just do north and south
pop.id3$Place <- gsub('North Inlet, SC', 'south', pop.id3$Place)
pop.id3$Place <- gsub('Little Egg Inlet, NJ', 'north', pop.id3$Place)
pop.id3$Place <- gsub('York River, VA', 'north', pop.id3$Place)
pop.id3$Place <- gsub('Roosevelt Inlet, DE', 'north', pop.id3$Place)
pop.id3$Place <- as.factor(pop.id3$Place)

rownames(larvs@tab) == pop.id3$ID..

larvs_amova <- pegas::amova(larvs_dist ~ Place/Time_within_NSLocations, data = pop.id3, nperm = 1000) # no difference between time periods, but no p-value for place? Lowering nperm seems to help.
# larvs_amova <- pegas::amova(larvs_dist ~ Time/NSLocations_within_Time, data = pop.id3, nperm = 1000)
larvs_amova

# Another way to do hierarchical AMOVA using poppr package
library("poppr")
pop_strata <- data.frame(cbind(pop.id3$Place, pop.id3$Time_within_NSLocations, pop.id3$Time, pop.id3$NSLocations_within_Time))
strata(larvs) <- pop_strata

# AMOVA between locations and time periods within locations
amova.poppr <- poppr.amova(larvs, ~ X1/X2, within = FALSE) # Default is ade4
amova.poppr.test <- randtest(amova.poppr, nrepet = 1000) #tests for significance
plot(amova.poppr.test) # variation between populations is significant, but this variation is small compared to variation within pops

# AMOVA between time periods and locations within locations
amova.poppr2 <- poppr.amova(larvs, ~ X3/X4, within = FALSE)
amova.poppr.test2 <- randtest(amova.poppr2, nrepet = 1000) #tests for significance when ade4 amova is used
plot(amova.poppr.test2)

#### AMOVA ####
larvs_dist <- dist(larvs)
larvs_strata <- strata(larvs)
larvs_amova <- pegas::amova(larvs_dist ~ geo, data = larvs_strata, nperm = 1000) # no difference between regions
larvs_amova

# AMOVA another way using the poppr package
table(strata(larvs, ~time_period))
table(strata(larvs, ~time_period/geo, combine = FALSE))
larvs_amova2 <- poppr.amova(larvs, ~time_period/geo, within = FALSE)
larvs_amova2
larvs_amova2_signif <- randtest(larvs_amova2, nrepet = 1000)
plot(larvs_amova2_signif)

fstat(larvs)

# AMOVA with just northern (or southern) fish broken out by time period
# First need to create different datasets for northern (and southern) fish
larvs_north <- pop_strata[which(pop_strata$geo == 'north'),]
larvs_south <- pop_strata[which(pop_strata$geo == 'south'),]

allelefreq.larvs.north <- allelefreqs_larvs[rownames(allelefreqs_larvs) %in% larvs_north$V1, ]
allelefreq.larvs.south <- allelefreqs_larvs[rownames(allelefreqs_larvs) %in% larvs_south$V1, ]

larvs_northonly <- as.genind(allelefreq.larvs.north)
strata(larvs_northonly) <- larvs_north[-1]
larvs_north_dist <- dist(larvs_northonly)
larvs_north_strata <- strata(larvs_northonly)
larvs_north_amova <- pegas::amova(larvs_north_dist ~ time_period, data = larvs_north_strata, nperm = 1000) # time period is not a significant contributor of explaining the variance in the north
larvs_north_amova

larvs_southonly <- as.genind(allelefreq.larvs.south)
strata(larvs_southonly) <- larvs_south[-1]
larvs_south_dist <- dist(larvs_southonly)
larvs_south_strata <- strata(larvs_southonly)
larvs_south_amova <- pegas::amova(larvs_south_dist ~ time_period, data = larvs_south_strata, nperm = 1000) # time period is not a significant contributor of explaining the variance in the south
larvs_south_amova



# This bit works but I'm not sure it is the AMOVA I want to be doing
full.freqs <- as.data.frame(full@tab)
larvs.dist <- dist(larvs.freqs)
strata(larvs.freqs) <- data.frame(pop(full)) #Setting the strata to populations
nameStrata(full) <- ~Population

larvs_dist <- dist(larvs) #Calculates a matrix of Euclidean distances between individuals using the genetic data, missing data is ignored
larvs_strata <- strata(larvs)
str(larvs_dist)

larvs_amova_geo <- amova(larvs_dist ~ geo, data = larvs_strata, nperm = 100) #Just testing for population differentiation
larvs_amova_geo #Populations can only explain a small amount of variance, most of the variance is within a population --> stability of population structure

# Permutation test of AMOVA
full_amova <- amova(full_dist ~ Population, data = full_stra) # 1000 permutations by default
full_amova # results are robust, populations are significantly different

# Trying another implementation of AMOVA
library("vegan")
amova.vegan <- adonis(larvs_dist ~ geo, data = larvs_strata, permutations = 1000)
amova.vegan #Populations are different

library("poppr")
amova.poppr <- poppr.amova(larvs, ~geo, within = FALSE)
amova.poppr.test <- randtest(amova.poppr, nrepet = 1000) #tests for significance
plot(amova.poppr.test) # variation between populations is significant, but this variation is small compared to variation within pops

# Might make sense to create another hierarchical level. Right now it's just "population", which refers to time period. Maybe time period and geographic location/region.
# Or grouping fish by larvae vs adults too


#### Allele frequecies over time ####
allelefreqs <- as.data.frame(scaleGen(larvs, center = FALSE, scale = FALSE, NA.method = "mean")) # this goes from counts to frequencies
rownames(allelefreqs)
dim(allelefreqs)
rownames(allelefreqs) == ordered_larvae_database_sub$ID..

allelefreq.larvs.mean2 <- aggregate(allelefreqs, list(pop_strata$time_period), mean, na.rm = TRUE)

# This way is not as good as above
# allelefreq.larvs.mean <- aggregate(allelefreqs.larvs.place[,3:3319], list(allelefreqs.larvs.place$period), mean) #way 1, this is bad if the column indexes change
# dim(allelefreq.larvs.mean) # 3 x 3318

# Check that means make sense
mean(allelefreqs[59:204,3], na.rm = T) # mean of snp3 during late period
allelefreq.larvs.mean2[,4] # match!
mean(allelefreqs[205:264,1300], na.rm = T) # early, snpn
allelefreq.larvs.mean2[,1301]

# Calculate allele frequency differences between time periods
allelefreq.larvs.mean2.t <- allelefreq.larvs.mean2[,-1] # first need to get rid of the categories so there are fewer problems down the line (num turning to chr)
allelefreq.larvs.mean2.t <- t(allelefreq.larvs.mean2.t) # columns are "early", "half" and "late"

half2early <- allelefreq.larvs.mean2.t[,2] - allelefreq.larvs.mean2.t[,1] # first half of the time period
late2half <- allelefreq.larvs.mean2.t[,3] - allelefreq.larvs.mean2.t[,2] # second half of the time period
late2early <- allelefreq.larvs.mean2.t[,3] - allelefreq.larvs.mean2.t[,1] # allele frequency difference over full time period

hist(half2early)
hist(late2half)
hist(late2early)

boxplot(half2early, xlab = "half2early")
boxplot(late2half, xlab = "late2half")
boxplot(late2early, xlab = "late2early")

summary(half2early)
summary(late2early) # when NA's are ignored, max difference between early and late time periods is 0.18
summary(late2half)

which(late2early > 0.10)
which(half2early > 0.13)
which(late2half > 0.12)

### To make a nice histogram of allele frequency change between earliest and latest time period ###
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/results/temporalchange_allelefreq.png", width=6, height=5, res=300, units="in")

par(
  mar=c(5, 4, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=15, # point size, which is the font size
  bg=NA
)

hist(late2early, main = NULL, xlab = "Temporal change in allele frequency")

dev.off()
#### What about mean allele frequency changes between north & south and by time period?
# Merge the larval samping database and allele frequency dataframe
allelefreqs$ID <- rownames(allelefreqs) # adding an ID column so the IDs are not just rownames. The names are already sorted and match the place data, so another way to do this would be to cbind.
allelefreqs_database <- merge(pop_strata, allelefreqs, by.x = "V1", by.y = "ID")
dim(allelefreqs_database) # 293 x 3767

# Dividing allele frequency dataset into north and south based on sampling location
allelefreqs.larvs.north <- allelefreqs_database[which(allelefreqs_database$geo == 'north'), ]
dim(allelefreqs.larvs.north) # 146 x 3767

allelefreqs.larvs.south <- allelefreqs_database[which(allelefreqs_database$geo == 'south'), ]
dim(allelefreqs.larvs.south) # 147 x 3767

se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x))) # Function to calculate SE

drops <- c("V1", "time_period", "geo")
allelefreq.larvs.mean.north <- allelefreqs.larvs.north[, !(colnames(allelefreqs.larvs.north) %in% drops)] # removing columns that aren't related to allele frequency
allelefreq.larvs.se.north <- aggregate(allelefreq.larvs.mean.north, list(allelefreqs.larvs.north$time_period), se) #take se by period
allelefreq.larvs.mean.north <- aggregate(allelefreq.larvs.mean.north, list(allelefreqs.larvs.north$time_period), mean, na.rm = TRUE) #take mean by period
allelefreq.larvs.mean.south <- allelefreqs.larvs.south[, !(colnames(allelefreqs.larvs.south) %in% drops)] # removing columns that aren't related to allele frequency
allelefreq.larvs.se.south <- aggregate(allelefreq.larvs.mean.south, list(allelefreqs.larvs.south$time_period), se) #take se by period
allelefreq.larvs.mean.south <- aggregate(allelefreq.larvs.mean.south, list(allelefreqs.larvs.south$time_period), mean, na.rm = TRUE) #take mean by period

# Calculate allele frequency differences between time periods & region
allelefreq.larvs.mean.north.t <- allelefreq.larvs.mean.north[,-1] # first need to get rid of the period categories so there are fewer problems down the line (num turning to chr)
allelefreq.larvs.mean.south.t <- allelefreq.larvs.mean.south[,-1]
allelefreq.larvs.se.north.t <- allelefreq.larvs.se.north[,-1]
allelefreq.larvs.se.south.t <- allelefreq.larvs.se.south[,-1]

allelefreq.larvs.mean.north.t <- t(allelefreq.larvs.mean.north.t) # Transpose: columns are "early", "half" and "late"
allelefreq.larvs.mean.south.t <- t(allelefreq.larvs.mean.south.t)
allelefreq.larvs.se.north.t <- t(allelefreq.larvs.se.north.t)
allelefreq.larvs.se.south.t <- t(allelefreq.larvs.se.south.t)

# Difference between time periods within a region
north.half2north.early <- allelefreq.larvs.mean.north.t[,2] - allelefreq.larvs.mean.north.t[,1]
north.late2north.early <- allelefreq.larvs.mean.north.t[,3] - allelefreq.larvs.mean.north.t[,1]
north.late2north.half <- allelefreq.larvs.mean.north.t[,3] - allelefreq.larvs.mean.north.t[,2]

south.half2south.early <- allelefreq.larvs.mean.south.t[,2] - allelefreq.larvs.mean.south.t[,1]
south.late2south.early <- allelefreq.larvs.mean.south.t[,3] - allelefreq.larvs.mean.south.t[,1]
south.late2south.half <- allelefreq.larvs.mean.south.t[,3] - allelefreq.larvs.mean.south.t[,2]

hist(north.half2north.early)
hist(north.late2north.early, breaks = 20)
hist(north.late2north.half) 

boxplot(north.half2north.early)
boxplot(north.late2north.early)
boxplot(north.late2north.half) 

hist(south.half2south.early)
hist(south.late2south.early)
hist(south.late2south.half)

boxplot(south.half2south.early)
boxplot(south.late2south.early)
boxplot(south.late2south.half)

summary(north.late2north.early) # when NA's ignored, max freq difference is 0.365, slightly higher when NA's not replaced with mean
summary(south.late2south.early) # when NA's ignored, max freq difference is 0.294

par(mfrow = c(1,2))
lmts <- range(north.late2north.early, north.late2north.early)
boxplot((north.late2north.early), ylim = lmts)
boxplot(south.late2south.early, ylim = lmts)

# Grouping data so that I can plot multiple boxplots on a single plot
dat <- data.frame(cbind(north.late2north.early, south.late2south.early))

## Reshape data wide -> long
library(reshape2)
long <- melt(dat)

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/results/temporalchange_allelefreq_byregions.png", width=4, height=5, res=300, units="in")

par(
  mar=c(5, 5, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=15, # point size, which is the font size
  bg=NA
)
plot(value ~ variable, data=long, xaxt = 'n', yaxt = 'n', xlab = NULL, ylab = NULL, col = c("#0072B2", "#D55E00"))
axis(1, at=seq(1,2, by=1), labels=c("North", "South"))
axis(2, at= seq(-0.4, 0.4, by=0.2), labels = seq(-0.4, 0.4, by=0.2), las = 2)
mtext ('Temporal change in allele frequency', 2, 3.1)

dev.off()

# Plot loci whose allele frequencies have changed the most over time
# Allele frequencies that have changed the most over time in the south
# south_temporal_outliers <- which(abs(as.vector(south.late2south.early)) > (sd(south.late2south.early))*3) # Gives both alleles
# south_temporal_outliers_single <- south_temporal_outliers[seq(1, length(south_temporal_outliers), 2)]
south_temporal_outliers <- which(as.vector(south.late2south.early) > (sd(south.late2south.early, na.rm = TRUE))*3)

par(mfrow = c(3,4))
for (i in south_temporal_outliers){
  plot(allelefreq.larvs.mean.north.t[10,], ylim=c(0,1), col = NULL, xlab = 'Time Period', ylab = 'Allele Frequency', xaxt = 'n', main = paste(rownames(allelefreq.larvs.mean.south.t)[i]))
  axis(1, at=1:3, labels = c('1989-1993', '1998-2002', '2008-2012'))
  lines(allelefreq.larvs.mean.north.t[i,])
  lines(allelefreq.larvs.mean.south.t[i,], col = 'tomato')
  arrows(c(1,2,3), allelefreq.larvs.mean.north.t[i,] - allelefreq.larvs.se.north.t[i,], c(1,2,3), allelefreq.larvs.mean.north.t[i,] + allelefreq.larvs.se.north.t[i,], length = 0.05, angle = 90, code = 3)
  arrows(c(1,2,3), allelefreq.larvs.mean.south.t[i,] - allelefreq.larvs.se.south.t[i,], c(1,2,3), allelefreq.larvs.mean.south.t[i,] + allelefreq.larvs.se.south.t[i,], length = 0.05, angle = 90, code = 3, col = 'tomato')
}

# Allele frequencies that have changed the most over time in the north
north_temporal_outliers <- which(as.vector(north.late2north.early) > (sd(north.late2north.early, na.rm = TRUE))*3)

par(mfrow = c(4,4))
for (i in north_temporal_outliers){
  plot(allelefreq.larvs.mean.north.t[10,], ylim=c(0,1), col = NULL, xlab = 'Time Period', ylab = 'Allele Frequency', xaxt = 'n', main = paste(rownames(allelefreq.larvs.mean.north.t)[i]))
  axis(1, at=1:3, labels = c('1989-1993', '1998-2002', '2008-2012'))
  lines(allelefreq.larvs.mean.north.t[i,])
  lines(allelefreq.larvs.mean.south.t[i,], col = 'tomato')
  arrows(c(1,2,3), allelefreq.larvs.mean.north.t[i,] - allelefreq.larvs.se.north.t[i,], c(1,2,3), allelefreq.larvs.mean.north.t[i,] + allelefreq.larvs.se.north.t[i,], length = 0.05, angle = 90, code = 3)
  arrows(c(1,2,3), allelefreq.larvs.mean.south.t[i,] - allelefreq.larvs.se.south.t[i,], c(1,2,3), allelefreq.larvs.mean.south.t[i,] + allelefreq.larvs.se.south.t[i,], length = 0.05, angle = 90, code = 3, col = 'tomato')
}

#### How can I summarize allele frequency plots over time for multiple loci into a single plot? ####
# Difference between northern and southern allele frequency in latest time period/difference between northern and southern allele frequency in earliest time period on one axis
# Northern allele frequency slope over time/southern allele frequency slope over time on other axis
north.late2south.late # these are below
north.early2south.early

# I need a for loop to fit a line to all northern AND all southern allele frequencies
n <- 3764
north.lms <- lapply(1:n, function(x) lm(allelefreq.larvs.mean.north.t[x,] ~ c(1,2,3)))
north.coeff <- sapply(north.lms, coef) # extract coefficients
north.slopes <- north.coeff[2,]

south.lms <- lapply(1:n, function(x) lm(allelefreq.larvs.mean.south.t[x,] ~ c(1,2,3)))
south.coeff <- sapply(south.lms, coef) # extract coefficients
south.slopes <- south.coeff[2,]

axis1 <- abs(north.slopes)-abs(south.slopes)
# axis1 <- axis1[seq(1,length(axis1),2)] # plots just one allele per locus
odd_indexes<-seq(1,3763,2)
axis1.odds <- axis1[odd_indexes] # 1882
save(axis1.odds, file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/axis1.odds.RData")

axis2 <- abs(north.early2south.early)-abs(north.late2south.late)
# axis2 <- axis2[seq(1,length(axis2),2)] # plots just one allele per locus
odd_indexes<-seq(1,3763,2)
axis2.odds <- axis2[odd_indexes] # 1882
save(axis2.odds, file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/axis2.odds.RData")

plot(axis1 ~ axis2, ylab = '|Slope in north| - |Slope in south|', xlab = '|Difference between early NvS| - |Difference between late NvS|')
abline(h = 0, v = 0, col = 'blue')
points(axis1[north_temporal_outliers] ~ axis2[north_temporal_outliers], col = 'red')
points(axis1[south_temporal_outliers] ~ axis2[south_temporal_outliers], col = 'gold')
points(axis1[early_spatial_outliers] ~ axis2[early_spatial_outliers], col = 'green')

# For loop to plot all mean allele frequencies over time & linear regressions for each region
n <- 3764
pdf(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/results/allregionalallelefreqs.pdf")
par(mfrow = c(4,5))
for (i in 1:n){
plot(allelefreq.larvs.mean.north.t[10,], ylim=c(0,1), col = NULL, xlab = 'Time Period', ylab = 'Allele Frequency', xaxt = 'n', main = paste(rownames(allelefreq.larvs.mean.north.t)[i]))
axis(1, at=1:3, labels = c('1989-1993', '1998-2002', '2008-2012'))
points(allelefreq.larvs.mean.north.t[i,] ~ c(1,2,3), ylim = c(0,1))
points(allelefreq.larvs.mean.south.t[i,] ~ c(1,2,3), col = 'red')
abline(north.lms[[i]])
abline(south.lms[[i]], col = 'red')
}
dev.off()

# which(abs(north.late2north.early) > 0.15)
# png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/Larvae/regionalallelefreqthrutime.png", width=6, height=5, res=300, units="in")
# par(mar=c(4.1, 4.1, 2, 9), xpd=TRUE,
#     mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
#     tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
#     cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
#     ps=14, # point size, which is the font size
#     bg=NA)
# plot(allelefreq.larvs.mean.north.t[190,], ylim=c(0,1), col = NULL, xlab = 'Time Period', ylab = 'Allele Frequency', xaxt = 'n')
# axis(1, at=1:3, labels = c('1989-1993', '1998-2002', '2008-2012'))
# legend("topright",
#        c("North", "South"),
#        lty=c(1,1), 
#        lwd=c(2,2),col=c("black", "tomato"),
#        inset = c(-0.45,0))
# lines(allelefreq.larvs.mean.north.t[190,])
# lines(allelefreq.larvs.mean.north.t[289,])
# lines(allelefreq.larvs.mean.north.t[653,])
# lines(allelefreq.larvs.mean.north.t[840,])
# lines(allelefreq.larvs.mean.north.t[893,])
# lines(allelefreq.larvs.mean.north.t[929,])
# lines(allelefreq.larvs.mean.north.t[973,])
# lines(allelefreq.larvs.mean.north.t[1035,])
# lines(allelefreq.larvs.mean.north.t[1179,])
# lines(allelefreq.larvs.mean.north.t[1225,])
# lines(allelefreq.larvs.mean.north.t[1243,])
# lines(allelefreq.larvs.mean.north.t[2217,])
# lines(allelefreq.larvs.mean.north.t[2265,])
# lines(allelefreq.larvs.mean.north.t[2569,])
# lines(allelefreq.larvs.mean.north.t[2832,])
# lines(allelefreq.larvs.mean.north.t[2863,])
# lines(allelefreq.larvs.mean.north.t[3054,])
# lines(allelefreq.larvs.mean.north.t[3249,])
# 
# # Plotting the allele frequency of the same alleles in the south
# lines(allelefreq.larvs.mean.south.t[190,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[289,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[653,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[840,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[893,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[929,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[973,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[1035,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[1179,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[1225,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[1243,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[2217,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[2265,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[2569,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[2832,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[2863,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[3054,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[3249,], col='tomato')
# 
# dev.off()

# # Hard to see patterns when all plotted together, so here's just a few alleles over time so that you can see the convergence
# png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/Larvae/regionalallelefreqthrutime_sample.png", width=6, height=5, res=300, units="in")
# par(mar=c(4.1, 4.1, 2, 9), xpd=TRUE,
#     mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
#     tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
#     cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
#     ps=14, # point size, which is the font size
#     bg=NA)
# plot(allelefreq.larvs.mean.north.t[190,], ylim=c(0,1), col = NULL, xlab = 'Time Period', ylab = 'Allele Frequency', xaxt = 'n')
# axis(1, at=1:3, labels = c('1989-1993', '1998-2002', '2008-2012'))
# legend("topright",
#        c("North", "South"),
#        lty=c(1,1), 
#        lwd=c(2,2),col=c("black", "tomato"),
#        inset = c(-0.45,0))
# lines(allelefreq.larvs.mean.north.t[190,])
# lines(allelefreq.larvs.mean.north.t[289,])
# lines(allelefreq.larvs.mean.north.t[893,])
# lines(allelefreq.larvs.mean.north.t[929,])
# lines(allelefreq.larvs.mean.north.t[2863,])
# 
# lines(allelefreq.larvs.mean.south.t[190,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[289,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[893,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[929,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[2863,], col='tomato')
# 
# dev.off()

# Same time period, different regions
plot(allelefreq.larvs.mean.north.t [,1]~ allelefreq.larvs.mean.south.t[,1])
plot(allelefreq.larvs.mean.north.t [,2]~ allelefreq.larvs.mean.south.t[,2])
plot(allelefreq.larvs.mean.north.t [,3]~ allelefreq.larvs.mean.south.t[,3])

# Same time period, different regions, single allele
plot(allelefreq.larvs.mean.north.t [,1][seq(1, length(allelefreq.larvs.mean.north.t [,1]), 2)]~ allelefreq.larvs.mean.south.t[,1][seq(1, length(allelefreq.larvs.mean.south.t [,1]), 2)])
abline(0,1)
plot(allelefreq.larvs.mean.north.t [,2][seq(1, length(allelefreq.larvs.mean.north.t [,1]), 2)]~ allelefreq.larvs.mean.south.t[,2][seq(1, length(allelefreq.larvs.mean.south.t [,1]), 2)])
abline(0,1)
plot(allelefreq.larvs.mean.north.t [,3][seq(1, length(allelefreq.larvs.mean.north.t [,1]), 2)]~ allelefreq.larvs.mean.south.t[,2][seq(1, length(allelefreq.larvs.mean.south.t [,1]), 2)])
abline(0,1)

# Spatial differences within a time period
north.late2south.late <- allelefreq.larvs.mean.north.t [,3] - allelefreq.larvs.mean.south.t[,3]
north.half2south.half <- allelefreq.larvs.mean.north.t [,2] - allelefreq.larvs.mean.south.t[,2]
north.early2south.early <- allelefreq.larvs.mean.north.t [,1] - allelefreq.larvs.mean.south.t[,1]

hist(north.early2south.early, breaks = 20)
hist(north.half2south.half, breaks = 20)
hist(north.late2south.late, breaks = 20) # the tails of the north/south differences are clearly smaller at the end of the time period

# Which ones are the spatial outliers in each of the time periods?
# Early spatial outliers
early_spatial_outliers <- which(as.vector(north.early2south.early) > (sd(north.early2south.early, na.rm = TRUE))*3)

par(mfrow = c(5,4))
for (i in early_spatial_outliers){
  plot(allelefreq.larvs.mean.north.t[10,], ylim=c(0,1), col = NULL, xlab = 'Time Period', ylab = 'Allele Frequency', xaxt = 'n', main = paste(rownames(allelefreq.larvs.mean.north.t)[i]))
  axis(1, at=1:3, labels = c('1989-1993', '1998-2002', '2008-2012'))
  lines(allelefreq.larvs.mean.north.t[i,])
  lines(allelefreq.larvs.mean.south.t[i,], col = 'tomato')
  arrows(c(1,2,3), allelefreq.larvs.mean.north.t[i,] - allelefreq.larvs.se.north.t[i,], c(1,2,3), allelefreq.larvs.mean.north.t[i,] + allelefreq.larvs.se.north.t[i,], length = 0.05, angle = 90, code = 3)
  arrows(c(1,2,3), allelefreq.larvs.mean.south.t[i,] - allelefreq.larvs.se.south.t[i,], c(1,2,3), allelefreq.larvs.mean.south.t[i,] + allelefreq.larvs.se.south.t[i,], length = 0.05, angle = 90, code = 3, col = 'tomato')
}

# Half spatial outliers
half_spatial_outliers <- which(as.vector(north.half2south.half) > (sd(north.half2south.half, na.rm = TRUE))*3)

par(mfrow = c(3,5))
for (i in half_spatial_outliers){
  plot(allelefreq.larvs.mean.north.t[10,], ylim=c(0,1), col = NULL, xlab = 'Time Period', ylab = 'Allele Frequency', xaxt = 'n', main = paste(rownames(allelefreq.larvs.mean.north.t)[i]))
  axis(1, at=1:3, labels = c('1989-1993', '1998-2002', '2008-2012'))
  lines(allelefreq.larvs.mean.north.t[i,])
  lines(allelefreq.larvs.mean.south.t[i,], col = 'tomato')
  arrows(c(1,2,3), allelefreq.larvs.mean.north.t[i,] - allelefreq.larvs.se.north.t[i,], c(1,2,3), allelefreq.larvs.mean.north.t[i,] + allelefreq.larvs.se.north.t[i,], length = 0.05, angle = 90, code = 3)
  arrows(c(1,2,3), allelefreq.larvs.mean.south.t[i,] - allelefreq.larvs.se.south.t[i,], c(1,2,3), allelefreq.larvs.mean.south.t[i,] + allelefreq.larvs.se.south.t[i,], length = 0.05, angle = 90, code = 3, col = 'tomato')
}

# Late spatial outliers
late_spatial_outliers <- which(as.vector(north.late2south.late) > (sd(north.late2south.late, na.rm = TRUE))*3)

par(mfrow = c(3,5))
for (i in late_spatial_outliers){
  plot(allelefreq.larvs.mean.north.t[10,], ylim=c(0,1), col = NULL, xlab = 'Time Period', ylab = 'Allele Frequency', xaxt = 'n', main = paste(rownames(allelefreq.larvs.mean.north.t)[i]))
  axis(1, at=1:3, labels = c('1989-1993', '1998-2002', '2008-2012'))
  lines(allelefreq.larvs.mean.north.t[i,])
  lines(allelefreq.larvs.mean.south.t[i,], col = 'tomato')
  arrows(c(1,2,3), allelefreq.larvs.mean.north.t[i,] - allelefreq.larvs.se.north.t[i,], c(1,2,3), allelefreq.larvs.mean.north.t[i,] + allelefreq.larvs.se.north.t[i,], length = 0.05, angle = 90, code = 3)
  arrows(c(1,2,3), allelefreq.larvs.mean.south.t[i,] - allelefreq.larvs.se.south.t[i,], c(1,2,3), allelefreq.larvs.mean.south.t[i,] + allelefreq.larvs.se.south.t[i,], length = 0.05, angle = 90, code = 3, col = 'tomato')
}

# Grouping data so that I can plot multiple boxplots on a single plot
regionaldiffs <- data.frame(cbind(north.early2south.early, north.half2south.half, north.late2south.late))

## Reshape data wide -> long
library(reshape2)
library(Hmisc)
regionaldiffs.long <- melt(regionaldiffs)

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/results/regional diffs_bytime.png", width=6, height=5, res=300, units="in")

par(
  mar=c(5, 5, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA
)

plot(value ~ variable, data=regionaldiffs.long, xaxt = 'n', yaxt = 'n', xlab = NULL, ylab = NULL, col = azur(3))
mgp.axis(1, at=seq(1,3, by=1), labels=c("Early\n(1989-1993)", "Middle\n(1998-2002)", "Late\n(2008-2012)"), mgp=c(3, 2, 0))
axis(2, at=seq(-0.4,0.4, by = 0.2), labels=seq(-0.4,0.4, by= 0.2), las = 2)
mtext('Regional difference in allele frequency', side = 2, line = 3.1)
     
dev.off()

summary(north.early2south.early)
summary(north.half2south.half)
summary(north.late2south.late) 

which(abs(north.late2south.late) > 0.14)
which(abs(north.half2south.half) > 0.16)
which(abs(north.early2south.early) > 0.3) # spatial outliers between north vs south

# Visualizing allele frequencies over space/time
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/Larvae/regional_allelefreq_dif_early.png", width=5, height=4.5, res=300, units="in")
hist(north.early2south.early, xlab = 'Regional allele frequency difference (1989-1993)', main = NULL)
dev.off()

# # Across time periods and regions
# plot(allelefreq.larvs.mean.north.t [,3]~ allelefreq.larvs.mean.south.t[,1])
# 
# # Difference in allele frequencies between north and south within a time period
# ns <- allelefreq.larvs.mean.north.t - allelefreq.larvs.mean.south.t
# which((ns[,3]- ns[,1]) > 0.20)
# 
# sn <- allelefreq.larvs.mean.south.t - allelefreq.larvs.mean.north.t
# which((sn[,3]- sn[,1]) > 0.20)

##### Plotting greatest allele frequency changes in the north for the 293 larval dataset #####
late_spatial_outliers <- which(abs(as.vector(north.late2south.late)) > (sd(north.late2south.late, na.rm = TRUE)*3))
late_spatial_outliers_single <- late_spatial_outliers[seq(1, length(late_spatial_outliers), 2)]

half_spatial_outliers <- which(abs(as.vector(north.half2south.half)) > (sd(north.half2south.half, na.rm = TRUE)*3))
half_spatial_outliers_single <- half_spatial_outliers[seq(1, length(half_spatial_outliers), 2)]

early_spatial_outliers <- which(abs(as.vector(north.early2south.early)) > (sd(north.early2south.early, na.rm = TRUE)*3))
early_spatial_outliers_single <- early_spatial_outliers[seq(1, length(early_spatial_outliers), 2)]

par(mfrow = c(5,4))
for (i in late_spatial_outliers_single){
  plot(allelefreq.larvs.mean.north.t[10,], ylim=c(0,1), col = NULL, xlab = 'Time Period', ylab = 'Allele Frequency', xaxt = 'n', main = paste(rownames(allelefreq.larvs.mean.south.t)[i]))
  axis(1, at=1:3, labels = c('1989-1993', '1998-2002', '2008-2012'))
  lines(allelefreq.larvs.mean.north.t[i,], lty = 2, col = 'gray50') # Blue = '#0072b2'
  lines(allelefreq.larvs.mean.south.t[i,]) # Red = '#D55E00'
  arrows(c(1,2,3), allelefreq.larvs.mean.north.t[i,] - allelefreq.larvs.se.north.t[i,], c(1,2,3), allelefreq.larvs.mean.north.t[i,] + allelefreq.larvs.se.north.t[i,], length = 0.05, angle = 90, code = 3, col = 'gray50')
  arrows(c(1,2,3), allelefreq.larvs.mean.south.t[i,] - allelefreq.larvs.se.south.t[i,], c(1,2,3), allelefreq.larvs.mean.south.t[i,] + allelefreq.larvs.se.south.t[i,], length = 0.05, angle = 90, code = 3)
  # legend("bottomright",
         # c("North", "South"),
         # lty=c(1,1),
         # lwd=c(2,2),col=c("#0072b2", "#D55E00"))
}

lines(allelefreq.larvs.mean.north.t[127,])
lines(allelefreq.larvs.mean.north.t[245,])
lines(allelefreq.larvs.mean.north.t[319,])
lines(allelefreq.larvs.mean.north.t[347,])
lines(allelefreq.larvs.mean.north.t[563,])
lines(allelefreq.larvs.mean.north.t[761,])
lines(allelefreq.larvs.mean.north.t[867,])
lines(allelefreq.larvs.mean.north.t[1979,])
lines(allelefreq.larvs.mean.north.t[2045,])
lines(allelefreq.larvs.mean.north.t[2111,])
lines(allelefreq.larvs.mean.north.t[2359,])

lines(allelefreq.larvs.mean.south.t[127,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[245,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[319,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[347,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[563,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[761,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[867,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[1979,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[2045,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[2111,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[2359,], col = 'tomato')

which(abs(south.late2south.early) > 0.15)

#### Nice plot for a single locus ####
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/results/regionalallelefreqthrutime_sample.png", width=6, height=5, res=300, units="in")
par(mar=c(5, 5, 4, 2),
    mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
    tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
    cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
    ps=15, # point size, which is the font size
    bg=NA)
plot(allelefreq.larvs.mean.north.t[10,], ylim=c(0,1), col = NULL, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', main = paste(rownames(allelefreq.larvs.mean.south.t)[79]))
axis(1, at=1:3, labels = c('1989-1993', '1998-2002', '2008-2012'))
axis(2, at = seq(0,1, by = 0.2), las = 2)
lines(allelefreq.larvs.mean.north.t[79,], col = '#0072b2')
lines(allelefreq.larvs.mean.south.t[79,], col = '#D55E00')
mtext('Time Period', 1, 2.7)
mtext('Allele Frequency', 2, 3.1)
legend("bottomright",
       c("North", "South"),
       lty=c(1,1),
       lwd=c(2,2),col=c("#0072b2", "#D55E00"))
       
dev.off()       

# Taking a look at spatial outliers from the latest time period
par(mfrow = c(3,3))
late_spatial_outliers <- which(abs(as.vector(north.late2south.late)) > 0.14)
late_spatial_outliers_single <- late_spatial_outliers[seq(1, length(late_spatial_outliers), 2)]
for (i in late_spatial_outliers_single){
  plot(allelefreq.larvs.mean.north.t[10,], ylim=c(0,1), col = NULL, xlab = 'Time Period', ylab = 'Allele Frequency', xaxt = 'n', main = paste('Locus', i))
  axis(1, at=1:3, labels = c('1989-1993', '1998-2002', '2008-2012'))
  lines(allelefreq.larvs.mean.north.t[i,])
  lines(allelefreq.larvs.mean.south.t[i,], col = 'tomato')
}

####################################################################
#### Plotting allele frequency by latitude for each time period ####
allelefreqs <- as.data.frame(scaleGen(larvs, center = FALSE, scale = FALSE, NA.method = "mean")) # this goes from counts to frequencies
rownames(allelefreqs)

early <- c('1989', '1990', '1991', '1992', '1993')
middle <- c('1998', '1999', '2000', '2001', '2002')
late <- c('2008', '2009', '2010', '2011', '2012')

early_allelefreqs <- allelefreqs[which(ordered_larvae_database_sub$Year %in% early),]
middle_allelefreqs <- allelefreqs[which(ordered_larvae_database_sub$Year %in% middle),]
late_allelefreqs <- allelefreqs[which(ordered_larvae_database_sub$Year %in% late),]

early_database <- ordered_larvae_database_sub[which(ordered_larvae_database_sub$Year %in% early),]
middle_database <- ordered_larvae_database_sub[which(ordered_larvae_database_sub$Year %in% middle),]
late_database <- ordered_larvae_database_sub[which(ordered_larvae_database_sub$Year %in% late),]

early_allelefreqs_mean <- aggregate(early_allelefreqs, list(early_database$Lat), mean, na.rm = TRUE)
dim(early_allelefreqs_mean)
middle_allelefreqs_mean <- aggregate(middle_allelefreqs, list(middle_database$Lat), mean, na.rm = TRUE)
dim(middle_allelefreqs_mean)
late_allelefreqs_mean <- aggregate(late_allelefreqs, list(late_database$Lat), mean, na.rm = TRUE)
dim(late_allelefreqs_mean)

plot(early_allelefreqs_mean$SNP_64.03 ~ early_allelefreqs_mean$Group.1, ylim = c(0, 1), xlab = "Latitude", col = NULL, xlim = c(33,40))
lines(early_allelefreqs_mean$SNP_64.03 ~ early_allelefreqs_mean$Group.1)
lines(middle_allelefreqs_mean$SNP_64.03 ~ middle_allelefreqs_mean$Group.1, col = "blue")
lines(late_allelefreqs_mean$SNP_64.03 ~ late_allelefreqs_mean$Group.1, col = "tomato")
legend("bottomleft",
       c("Early", "Middle", "Late"),
       lty=c(1,1,1),
       lwd=c(2,2, 2),col=c("black", "blue","tomato"))

plot(early_allelefreqs_mean$SNP_69.02 ~ early_allelefreqs_mean$Group.1, ylim = c(0, 1), xlab = "Latitude", col = NULL, xlim = c(33,40))
lines(early_allelefreqs_mean$SNP_69.02 ~ early_allelefreqs_mean$Group.1)
lines(middle_allelefreqs_mean$SNP_69.02 ~ middle_allelefreqs_mean$Group.1, col = "blue")
lines(late_allelefreqs_mean$SNP_69.02 ~ late_allelefreqs_mean$Group.1, col = "tomato")

#### ANCOVA or two way ANOVA on larval spatial and temporal outliers???? ####
pop_strata$V1 == rownames(allelefreqs) # are these datasets in the same order?
# Creating a dataframe of the data to be analyzed using ANCOVA
# Subset allele frequency table to only outliers
south_temporal_outliers <- which(abs(as.vector(south.late2south.early)) > (sd(south.late2south.early, na.rm = TRUE))*3)
south_temporal_outliers # 79   80  147  148  245  246  307  308  727  728 1133 1134 1235 1236 1775 1776 2045 2046 3379 3380 3417 3418
north_temporal_outliers <- which(abs(as.vector(north.late2north.early)) > (sd(north.late2north.early, na.rm = TRUE))*3)
north_temporal_outliers # 137  138  319  320  547  548  563  564  911  912 1001 1002 1381 1382 1619 1620 1725 1726 1851 1852 1929 1930 2003 2004 2187 2188 2269 2270 2719 2720 2903 2904
early_spatial_outliers #  127  128  137  138  245  246  319  320  347  348  547  548  563  564  761  762  867  868 1381 1382 1539 1540 1773 1774 1851 1852 1979 1980 2045 2046 2111 2112 2359 2360 2719 2720 2891 2892 3441 3442

outlier.indices <- unique(sort(c(south_temporal_outliers, north_temporal_outliers, early_spatial_outliers)))
outlier.indices.single <- outlier.indices[seq(1, length(outlier.indices), 2)]# take every other one

outlier.df <- allelefreqs[,outlier.indices.single]
dim(outlier.df) # 293 x 38
outlier.df <- cbind(pop_strata[,-1], outlier.df)
dim(outlier.df) #293 x 40

snp40.01 <- lm(SNP_40.01 ~ geo + time_period + geo:time_period, data = outlier.df)
library(car)
options(contrasts = c("contr.sum", "contr.poly"))
Anova(snp40.01, type = 'III')
summary(snp40.01)
snp64.03 <- lm(SNP_64.03 ~ time_period+geo+time_period:geo, data = outlier.df)
Anova(snp64.03, type = 'III')
summary(snp64.03)
snp69.02 <- aov(SNP_69.02 ~ time_period+geo, data = outlier.df)

sample1 <- outlier.df[which(outlier.df$geo == 'south' & outlier.df$time_period == 'early'),]

permute <- vector()
for (i in 1:100){
mysample <- sample1[sample(1:nrow(sample1), 8, replace = FALSE),]
mean <- mean(mysample$SNP_40.01)
permute <- rbind(permute, mean)
}

hist(permute[,1])
mean(permute[,1])

mean(outlier.df[which(outlier.df$time_period =='early' & outlier.df$geo == 'south'), 'SNP_40.01'])

############################################################
#### Create a str file for 293 larvae using 1904 loci, and where the pop column tells you what ingress site/time period the individual was captured in ####
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files")

library(ade4)
library(adegenet)
library("hierfstat")
library(pegas)

# Read in a text file containing 293 larvae and 1904 loci (some of which are multiallelic)
larvs <- read.table("structure_input_Mar14_2017_293_1904_9pops.txt", sep="\t", header = TRUE)

# Read in file with population assignment that denotes ingress site/time period combination
pops <- read.structure("structure_input_528_10outliers_10pops2.str", n.ind = 528, n.loc = 10, col.lab = 1, col.pop = 2, row.marknames = 1,
  onerowperind = FALSE)

# Isolate unique larval IDs
larvs.uniq <- as.character(unique(larvs$ID))

# Isolate pop identifier & cbind with ID
pop <- cbind(rownames(pops@tab), as.vector(pops@pop)) # this also has adults

# Get rid of adults
pop.larvs <- pop[pop[,1] %in% larvs.uniq,] #293 x2

# Order names in the order they are in the file that will be used for STRUCTURE
data <- pop.larvs[match(larvs.uniq, pop.larvs[,1]),]
data[,1] == larvs.uniq # names and pop identifier are in order

# Repeat ID/pop identifier & write file so that it can be pasted into structure_input_Mar14_2017_293_1904_9pops.txt
data2 <- data[rep(seq_len(nrow(data)), each=2),]
write.table(data2, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/293larvs_popIDs_forstructure.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

####################################################################################################################
#### Removing loci not in HWE for another STRUCTURE run ####
# Read in STRUCTURE formatted file of 293 larvae and 1904 loci
larvs <- read.structure("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/structure_input_Jan9_2019_293_1904_9pops.str",
                       n.ind = 293, n.loc = 1904, col.lab = 1, col.pop = 2, row.marknames = 1,
                       onerowperind = FALSE)

larvs.hwt <- hw.test(larvs, B=10000)
larvs.hwt
hist(larvs.hwt[,4], nclass =30)
pval <- larvs.hwt[larvs.hwt[,"Pr.exact"] < 0.0100,] # p<0.01, exact test; like in Wolf population structure & canidate genes under selection paper
length(pval[,"Pr.exact"]) 

# Now read in txt file version of STRUCTURE file and remove SNPs not in HWE (pval)
larvs2 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/structure_input_Jan9_2019_293_1904_9pops.txt', header = TRUE)

larvs2.hwe <- larvs2[,!colnames(larvs2) %in% rownames(pval)]

write.table(larvs2.hwe, '~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/structure_input_Jan9_2019_293_1646_9pops_hwe.txt', row.names = FALSE, col.names = TRUE, sep = '\t')
# Then change extension name from .txt to .str

#### Plot adult outlier loci in larvae ####
# Prepare data. Read in STRUCTURE formatted file of 293 larvae and 1904 loci
larvs <- read.structure("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/structure_input_Jan9_2019_293_1904_9pops.str",
                        n.ind = 293, n.loc = 1904, col.lab = 1, col.pop = 2, row.marknames = 1,
                        onerowperind = FALSE)
pop <- as.character(larvs@pop)

# Make vector of adult canidate SNPs to subset larval dataset. Remember that SNP numbers are different in adults and larvae
# cans <- c('SNP_35\\.', 'SNP_125\\.', 'SNP_214\\.', 'SNP_291\\.', 'SNP_396\\.', 'SNP_442\\.', 'SNP_499\\.', 'SNP_524\\.', 'SNP_542\\.', 'SNP_609\\.', 'SNP_615\\.', 'SNP_743\\.', 'SNP_825\\.', 'SNP_919\\.', 'SNP_990\\.') # these are the SNP #'s in the adults
cans <- c('SNP_170\\.', 'SNP_458\\.', 'SNP_685\\.', 'SNP_833\\.', 'SNP_1062\\.', 'SNP_1129\\.', 'SNP_1171\\.', 'SNP_1194\\.', 'SNP_1282\\.', 'SNP_1443\\.') # these are the corresponding SNP #'s in the larvae

# Search for SNP numbers
larvs.sub <- larvs@tab[,unique (grep(paste(cans,collapse="|"), colnames(larvs@tab), value=TRUE))]

# cbind the population identifiers to the genotype data
larvs.sub2 <- cbind.data.frame(pop, larvs.sub)
table(larvs.sub2[, 'pop'])

# Divide out fish into desired groupings
# Yearly
one <- larvs.sub2[which(larvs.sub2[, 'pop'] == '1'),]
two <- larvs.sub2[which(larvs.sub2[, 'pop'] == '2'),]
three <- larvs.sub2[which(larvs.sub2[, 'pop'] == '3'),]
four <- larvs.sub2[which(larvs.sub2[, 'pop'] == '4'),]
five <- larvs.sub2[which(larvs.sub2[, 'pop'] == '5'),]
six <- larvs.sub2[which(larvs.sub2[, 'pop'] == '6'),]
seven <- larvs.sub2[which(larvs.sub2[, 'pop'] == '7'),]
eight <- larvs.sub2[which(larvs.sub2[, 'pop'] == '8'),]
nine <- larvs.sub2[which(larvs.sub2[, 'pop'] == '9'),]

# period
early <- rbind.data.frame(one, four)
mid <- rbind.data.frame(two, five)
late <- rbind.data.frame(three, six, seven, eight, nine)

north.late <- rbind(three, seven, eight)
south.late <- rbind(six, nine)

# Calculate allele frequencies
# yearly
one2 <- colSums(one[,-1],na.rm=TRUE)/(2*colSums(!is.na(one[,-1]))) # only calculates frequencies of columns containing alleles
two2 <- colSums(two[,-1],na.rm=TRUE)/(2*colSums(!is.na(two[,-1])))
three2 <- colSums(three[,-1],na.rm=TRUE)/(2*colSums(!is.na(three[,-1])))
four2 <- colSums(four[,-1],na.rm=TRUE)/(2*colSums(!is.na(four[,-1])))
five2 <- colSums(five[,-1],na.rm=TRUE)/(2*colSums(!is.na(five[,-1])))
six2 <- colSums(six[,-1],na.rm=TRUE)/(2*colSums(!is.na(six[,-1])))
seven2 <- colSums(seven[,-1],na.rm=TRUE)/(2*colSums(!is.na(seven[,-1])))
eight2 <- colSums(eight[,-1],na.rm=TRUE)/(2*colSums(!is.na(eight[,-1])))
nine2 <- colSums(nine[,-1],na.rm=TRUE)/(2*colSums(!is.na(nine[,-1])))

# period
early2 <- colSums(early[,-1],na.rm=TRUE)/(2*colSums(!is.na(early[,-1])))
mid2 <- colSums(mid[,-1],na.rm=TRUE)/(2*colSums(!is.na(mid[,-1])))
late2 <- colSums(late[,-1],na.rm=TRUE)/(2*colSums(!is.na(late[,-1])))

north.early <- colSums(one[,-1],na.rm=TRUE)/(2*colSums(!is.na(one[,-1])))
south.early <- colSums(four[,-1],na.rm=TRUE)/(2*colSums(!is.na(four[,-1])))
north.mid <- colSums(two[,-1],na.rm=TRUE)/(2*colSums(!is.na(two[,-1])))
south.mid <- colSums(five[,-1],na.rm=TRUE)/(2*colSums(!is.na(five[,-1])))
north.late2 <- colSums(north.late[,-1],na.rm=TRUE)/(2*colSums(!is.na(north.late[,-1])))
south.late2 <- colSums(south.late[,-1],na.rm=TRUE)/(2*colSums(!is.na(south.late[,-1])))

# Plot
year <- rbind(one2, two2, three2, four2, five2, six2, seven2, eight2, nine2)

for (i in 1:ncol(year)) {
  plot(year[,i])
  lines(year[,i], col = 'tomato')
}

period <- rbind(early2, mid2, late2)

for (i in 1:ncol(period)) {
  plot(period[,i])
  lines(period[,i], col = 'tomato')
}

north <- rbind(north.early, north.mid, north.late2)
south <- rbind(south.early, south.mid, south.late2)

for (i in 1:ncol(north)) {
  plot(north[,i], ylim = c(0,1))
  lines(north[,i], col = 'tomato')
  lines(south[,i], col = 'black')
}




r <- vector()
for (i in 1:ncol(year)) {
  r[i] <- max(year[,i] - min(year[,i]))
}

