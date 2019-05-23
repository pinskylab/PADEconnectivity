#### This script fits GAMS between candidate loci and distance along the coast, which can then be used to figure out allele frequencies for an arbitrary number of 'populations' along the coast ####
#### Reads in a file of 232 fish and their allele counts at the 10 candidates present in larvae. Allele frequencies need to be calculated. ####
adults232 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.counts10.txt", header = TRUE)
rownames(adults232) <- adults232$PinskyID

# Replace NAs with column mean?
for(i in 1:ncol(adults232[,-c(1:2)])){
  adults232[,-c(1:2)][is.na(adults232[,-c(1:2)][,i]), i] <- mean(adults232[,-c(1:2)][,i], na.rm = TRUE)
}

# Or just remove the few individuals with some missing data? Same results as letting GAM run with some missing data
# adults232 <- na.omit(adults232)

colSums(is.na(adults232)) # Check there are no zeros in any of the columns

# Convert from count to frequency for allele columns
# adults232.freqs <- adults232[,-c(1:2)]/2
 
envi <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/232envirowithdist.txt", header = TRUE)
envi.ordered <- envi[with(envi, order(PinskyID)),] 
# envi.ordered <- envi.ordered[envi.ordered$PinskyID %in% rownames(adults232),] # If 4 individuals are removed b/c of missing data; 228 x 11
# as.character(envi.ordered[,2]) == rownames(adults.nonas.232.1137)
# as.character(envi.ordered[,2]) == rownames(adults232.freqs) # environmental data is in same order as genetic data?
# rownames(envi.ordered) <- envi.ordered[,"PinskyID"] # so that environmental variable units in GAM plots won't be standardized

# Set up for GAM with binomial errors
adults232 <- adults232[,-c(1:2)] # get rid of names and regions columns
adults232.split <- tapply(as.list(adults232), gl(ncol(adults232)/2, 2), as.data.frame) # split dataframe so that each locus is 232 x 2; array is 232 rows x 2 columns x 10 deep
as.character(envi.ordered[,2]) == rownames(adults232) # environmental data is in same order as genetic data?
rownames(envi.ordered) <- envi.ordered[,"PinskyID"] # so that environmental variable units in GAM plots won't be standardized

# Look up SNP number
# candidates <- c(35, 125 , 214, 291, 396, 442, 499, 524, 542, 609, 615, 743, 825, 919, 990)

# Run a GAM against bottom temperature? Or the associated environmental variable? Or all environmental variables?
# For loop to fit GAM between all 10 // 15 candidates and each 4 enviromental variables
library(mgcv)

# 10 loci
# for(i in 1:ncol(adults232.freqs)){
  # assign(paste0(colnames(adults232.freqs[i]), '.dist.gam'), mgcv::gam(adults232.freqs[,i] ~ s(dist), data = envi.ordered, method = 'REML', select = TRUE))
  # assign(paste0(colnames(adults232.freqs[i]), '.depth.gam'), mgcv::gam(adults232.freqs[,i] ~ s(depth), data = envi.ordered, method = 'REML', select = TRUE))
  # assign(paste0(colnames(adults232.freqs[i]), '.btemp.gam'), mgcv::gam(adults232.freqs[,i] ~ s(b_temp), data = envi.ordered, method = 'REML', select = TRUE))
  # assign(paste0(colnames(adults232.freqs[i]), '.bsalin.gam'), mgcv::gam(adults232.freqs[,i] ~ s(b_salin), data = envi.ordered, method = 'REML', select = TRUE))
# }

# 10 loci with binomial errors
for(i in 1:length(adults232.split)){
  assign(paste0(colnames(adults232.split[[i]][1]), '.dist.gam'), mgcv::gam(as.matrix(adults232.split[[i]]) ~ s(dist), data = envi.ordered, method = 'REML', select = TRUE, family = binomial))
}

# mgcv::gam(as.matrix(adults232[,5:6]) ~ s(dist), data = envi.ordered, method = 'REML', select = TRUE, family = binomial)
# predict.gam(test, data.frame(dist = dist.new), type = 'response')


# plot(GAM.object) as an object and this will save coordinates
par(mfrow = c(2,2))

gam.names <- mget(ls(pattern = ".dist.gam")) # searches through environment and pulls all objects with '.dist.gam'

# Plots all gams quickly!
for (i in 1:length(gam.names)){
  plot(gam.names[[i]], scheme = 1, seWithMean = TRUE, xlab = 'Distance along coast', main = paste0(names(gam.names)[[i]]))
}

# Define new data for each environmental variable to predict
dist.new <- seq(min(envi.ordered$dist), max(envi.ordered$dist), length = 500)
# depth.new <- seq(min(envi.ordered$depth), max(envi.ordered$depth), length = 500)
# btemp.new <- seq(min(envi.ordered$b_temp), max(envi.ordered$b_temp), length = 500)
# bsalin.new <- seq(min(envi.ordered$b_salin), max(envi.ordered$b_salin), length = 500)

for (i in 1:length(names(gam.names))){
  assign(paste0(names(gam.names)[[i]], '.p'), predict.gam(gam.names[[i]], data.frame(dist = dist.new), type = 'response'))
}

gam.names.p <- mget(ls(pattern = ".gam.p")) # searches through environment and pulls all objects with '.gam.p'

# Round any predicted allele frequency that is above 1 down to 0.999, or below 0 to 0.001. Allele frequencies of exactly 1 or zero will lead to null alleles and values of infinity. Cornuet (1999) says null allele frequencies can be replaced by a small value.
# Alleles with predicted allele frequencies > 1: dDocent_Contig_15075_20.02.dist.gam.p
# Alleles with predicted allele frequencies < 0: dDocent_Contig_15075_20.04.dist.gam.p
# gam.names.p$dDocent_Contig_15075_20.02.dist.gam.p[gam.names.p$dDocent_Contig_15075_20.02.dist.gam.p > 1] <- 0.999
# gam.names.p$dDocent_Contig_15075_20.04.dist.gam.p[gam.names.p$dDocent_Contig_15075_20.04.dist.gam.p < 0] <- 0.001

#### Plot newly predicted allele frequency ####
# Plots all predicted gams quickly!
for (i in 1:length(gam.names.p)){
  plot(gam.names.p[[i]] ~ dist.new, main = paste0(names(gam.names.p)[[i]]), ylab = 'Allele frequency')
}

# Break 500 bottom temps into 10 groups
break.points <- dist.new[c(25, 75, 125, 175, 225, 275, 325, 375, 425, 475)]

# Figure out corresponding allele frequency
# cbind temp and predicted allele frequency
for (i in 1:length(names(gam.names))){
  assign(paste0(names(gam.names)[[i]], '.joined'), cbind(dist.new, gam.names.p[[i]]))
}

gam.joined <- mget(ls(pattern = ".joined")) # searches through environment and pulls all objects with '.joined'

allele.freqs10 <- list()
for (i in 1:length(gam.joined)){
  allele.freqs10[[i]] <- gam.joined[[i]][gam.joined[[i]][,1] %in% break.points,][,2]
}

allele.freqs10.df <- data.frame(t(matrix(unlist(allele.freqs10), nrow=length(allele.freqs10), byrow=T)))

# Fix names for columns
allelenames.split <- as.data.frame(do.call(rbind, strsplit(as.character(names(gam.names)), '.', fixed = TRUE)))
allelenames.split2 <- as.factor(paste(allelenames.split$V1, '.', allelenames.split$V2, sep = ''))

colnames(allele.freqs10.df) <- allelenames.split2
# write.table(allele.freqs10.df, '~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.freqs.10pops.txt', row.names = FALSE, col.names = TRUE)
write.table(allele.freqs10.df, '~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.freqs.10pops.gambinomial.txt', row.names = FALSE, col.names = TRUE) # one allele for each locus

# # Keep only one allele per locus in adults. This is only necessary when the datafile contains both alleles.
# odds <- seq(1,20,2) # odd indicies to keep
# allele.freqs10.df.odds <- allele.freqs10.df[,odds]

#### Use these allele frequencies to calculate larval likelihoods for all 293 fish ####
# Read in dataset containing outlier loci
gen.larvae.outs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/masterPADElarvae.txt', header = TRUE, sep = "\t")

# Keep only larvae with genetic assignments
gen.larvae.outs2 <- gen.larvae.outs[which(is.na(gen.larvae.outs$assignment) == FALSE), ]

indiv.allele.counts <- gen.larvae.outs2[,c(16:35)] # Genetic data for all 293 larvae

# indiv.allele.counts.odds <- indiv.allele.counts[, colnames(indiv.allele.counts) %in% colnames(allele.freqs10.df.odds)] # subset alleles in larvae to those in adults when both alleles are in datafile
indiv.allele.counts.odds <- indiv.allele.counts[, colnames(indiv.allele.counts) %in% colnames(allele.freqs10.df)] # 293 x 10
indiv.allele.counts.odds[is.na(indiv.allele.counts.odds)] <- 9 # replace NA's with 9's to make the ifelse statements easier

# Reorder 10 adult allele frequencies to be in same order as larvae
# allele.freqs10.df.odds.ordered <- allele.freqs10.df.odds[names(indiv.allele.counts.odds)] # subset alleles in larvae to those in adults when both alleles are in datafile
allele.freqs10.df.odds.ordered <- allele.freqs10.df[names(indiv.allele.counts.odds)]
# write.table(allele.freqs10.df.odds.ordered, '~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.freqs.10pops.gambinomial.ordered.txt', row.names = FALSE, col.names = TRUE)

colnames(indiv.allele.counts.odds) == colnames(allele.freqs10.df.odds.ordered) # column names match? yay!

#### Now that I have allele frequencies for 10 adult groups from the GAMS and the allele counts for the larvae, I'm ready to calculate the genotype likelihood of each larval indivdiual coming from the 10 potential popuations of origin ####
# Population likelihoods for 10 groups
pop1.likelihoods <- data.frame()
pop2.likelihoods <- data.frame()
pop3.likelihoods <- data.frame()
pop4.likelihoods <- data.frame()
pop5.likelihoods <- data.frame()
pop6.likelihoods <- data.frame()
pop7.likelihoods <- data.frame()
pop8.likelihoods <- data.frame()
pop9.likelihoods <- data.frame()
pop10.likelihoods <- data.frame()

for (j in 1:length(rownames(indiv.allele.counts.odds))){
  
  for (i in 1:length(colnames(indiv.allele.counts.odds))){
    if(indiv.allele.counts.odds[j,i] == 2) {
      pop1.likelihoods[j,i] <- allele.freqs10.df.odds.ordered[1,i]^2
    } else if (indiv.allele.counts.odds[j,i] == 1) {
      pop1.likelihoods[j,i] <- 2*(allele.freqs10.df.odds.ordered[1,i] * (1-allele.freqs10.df.odds.ordered[1,i]))
    } else if (indiv.allele.counts.odds[j,i] == 0) {
      pop1.likelihoods[j,i] <- (1-allele.freqs10.df.odds.ordered[1,i])^2 
    } else {
      pop1.likelihoods[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(indiv.allele.counts.odds))){
    if(indiv.allele.counts.odds[j,i] == 2){
      pop2.likelihoods[j,i] <- allele.freqs10.df.odds.ordered[2,i]^2
    } else if (indiv.allele.counts.odds[j,i] == 1) {
      pop2.likelihoods[j,i] <- 2*(allele.freqs10.df.odds.ordered[2,i] * (1-allele.freqs10.df.odds.ordered[2,i]))
    } else  if (indiv.allele.counts.odds[j,i] == 0) {
      pop2.likelihoods[j,i] <- (1-allele.freqs10.df.odds.ordered[2,i])^2
    } else {
      pop2.likelihoods[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(indiv.allele.counts.odds))){
    if(indiv.allele.counts.odds[j,i] == 2) {
      pop3.likelihoods[j,i] <- allele.freqs10.df.odds.ordered[3,i]^2
    } else if (indiv.allele.counts.odds[j,i] == 1) {
      pop3.likelihoods[j,i] <- 2*(allele.freqs10.df.odds.ordered[3,i] * (1-allele.freqs10.df.odds.ordered[3,i]))
    } else if (indiv.allele.counts.odds[j,i] == 0) {
      pop3.likelihoods[j,i] <- ( 1-allele.freqs10.df.odds.ordered[3,i])^2 
    } else {
      pop3.likelihoods[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(indiv.allele.counts.odds))){
    if(indiv.allele.counts.odds[j,i] == 2) {
      pop4.likelihoods[j,i] <- allele.freqs10.df.odds.ordered[4,i]^2
    } else if (indiv.allele.counts.odds[j,i] == 1) {
      pop4.likelihoods[j,i] <- 2*(allele.freqs10.df.odds.ordered[4,i] * (1-allele.freqs10.df.odds.ordered[4,i]))
    } else if (indiv.allele.counts.odds[j,i] == 0) {
      pop4.likelihoods[j,i] <- ( 1-allele.freqs10.df.odds.ordered[4,i])^2 
    } else {
      pop4.likelihoods[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(indiv.allele.counts.odds))){
    if(indiv.allele.counts.odds[j,i] == 2) {
      pop5.likelihoods[j,i] <- allele.freqs10.df.odds.ordered[5,i]^2
    } else if (indiv.allele.counts.odds[j,i] == 1) {
      pop5.likelihoods[j,i] <- 2*(allele.freqs10.df.odds.ordered[5,i] * (1-allele.freqs10.df.odds.ordered[5,i]))
    } else if (indiv.allele.counts.odds[j,i] == 0) {
      pop5.likelihoods[j,i] <- ( 1-allele.freqs10.df.odds.ordered[5,i])^2 
    } else {
      pop5.likelihoods[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(indiv.allele.counts.odds))){
    if(indiv.allele.counts.odds[j,i] == 2) {
      pop6.likelihoods[j,i] <- allele.freqs10.df.odds.ordered[6,i]^2
    } else if (indiv.allele.counts.odds[j,i] == 1) {
      pop6.likelihoods[j,i] <- 2*(allele.freqs10.df.odds.ordered[6,i] * (1-allele.freqs10.df.odds.ordered[6,i]))
    } else if (indiv.allele.counts.odds[j,i] == 0) {
      pop6.likelihoods[j,i] <- ( 1-allele.freqs10.df.odds.ordered[6,i])^2 
    } else {
      pop6.likelihoods[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(indiv.allele.counts.odds))){
    if(indiv.allele.counts.odds[j,i] == 2) {
      pop7.likelihoods[j,i] <- allele.freqs10.df.odds.ordered[7,i]^2
    } else if (indiv.allele.counts.odds[j,i] == 1) {
      pop7.likelihoods[j,i] <- 2*(allele.freqs10.df.odds.ordered[7,i] * (1-allele.freqs10.df.odds.ordered[7,i]))
    } else if (indiv.allele.counts.odds[j,i] == 0) {
      pop7.likelihoods[j,i] <- ( 1-allele.freqs10.df.odds.ordered[7,i])^2 
    } else {
      pop7.likelihoods[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(indiv.allele.counts.odds))){
    if(indiv.allele.counts.odds[j,i] == 2) {
      pop8.likelihoods[j,i] <- allele.freqs10.df.odds.ordered[8,i]^2
    } else if (indiv.allele.counts.odds[j,i] == 1) {
      pop8.likelihoods[j,i] <- 2*(allele.freqs10.df.odds.ordered[8,i] * (1-allele.freqs10.df.odds.ordered[8,i]))
    } else if (indiv.allele.counts.odds[j,i] == 0) {
      pop8.likelihoods[j,i] <- ( 1-allele.freqs10.df.odds.ordered[8,i])^2 
    } else {
      pop8.likelihoods[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(indiv.allele.counts.odds))){
    if(indiv.allele.counts.odds[j,i] == 2) {
      pop9.likelihoods[j,i] <- allele.freqs10.df.odds.ordered[9,i]^2
    } else if (indiv.allele.counts.odds[j,i] == 1) {
      pop9.likelihoods[j,i] <- 2*(allele.freqs10.df.odds.ordered[9,i] * (1-allele.freqs10.df.odds.ordered[9,i]))
    } else if (indiv.allele.counts.odds[j,i] == 0) {
      pop9.likelihoods[j,i] <- ( 1-allele.freqs10.df.odds.ordered[9,i])^2 
    } else {
      pop9.likelihoods[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(indiv.allele.counts.odds))){
    if(indiv.allele.counts.odds[j,i] == 2) {
      pop10.likelihoods[j,i] <- allele.freqs10.df.odds.ordered[10,i]^2
    } else if (indiv.allele.counts.odds[j,i] == 1) {
      pop10.likelihoods[j,i] <- 2*(allele.freqs10.df.odds.ordered[10,i] * (1-allele.freqs10.df.odds.ordered[10,i]))
    } else if (indiv.allele.counts.odds[j,i] == 0) {
      pop10.likelihoods[j,i] <- ( 1-allele.freqs10.df.odds.ordered[10,i])^2 
    } else {
      pop10.likelihoods[j,i] <- 1
    }
  }
}

pop1.vector <- apply(pop1.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
pop2.vector <- apply(pop2.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
pop3.vector <- apply(pop3.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
pop4.vector <- apply(pop4.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
pop5.vector <- apply(pop5.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
pop6.vector <- apply(pop6.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
pop7.vector <- apply(pop7.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
pop8.vector <- apply(pop8.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
pop9.vector <- apply(pop9.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
pop10.vector <- apply(pop10.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)

bayenv.likelihoods.293indivs <- data.frame(gen.larvae.outs2$PinskyID, gen.larvae.outs2$Place, pop1.vector, pop2.vector, pop3.vector, pop4.vector, pop5.vector, pop6.vector, pop7.vector, pop8.vector, pop9.vector, pop10.vector) # this is for all 293 fish with genetic data
colnames(bayenv.likelihoods.293indivs) <- c("ID", "Place", "Pop1", "Pop2", "Pop3", "Pop4", "Pop5", "Pop6", "Pop7", "Pop8", "Pop9", "Pop10")
# write.table(bayenv.likelihoods.293indivs, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/bayenv_likelihoods_293indivs_10pops.txt", row.names = FALSE, col.names = TRUE)

# Which BayEnv population is each individual most likely from?
bayenv.likelihoods.293indivs <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/bayenv_likelihoods_293indivs_10pops.txt", header = TRUE)
most.like <- colnames(bayenv.likelihoods.293indivs[,-c(1:2)])[apply(bayenv.likelihoods.293indivs[,-c(1:2)],1, which.max)]
t(table(most.like,bayenv.likelihoods.293indivs$Place))

# Plot to summarize individual assignments
locations <- factor(bayenv.likelihoods.293indivs$Place, c("Little Egg Inlet, NJ", "Roosevelt Inlet, DE", "York River, VA", "Beaufort, NC", "North Inlet, SC"))
table <- t(table(most.like,locations))
table2 <- prop.table(table, 1)
table3 <- cbind(table2[,-2], table2[,2])

# Get colors
library(RColorBrewer)
library(plotly)
cols <- brewer.pal(10, "Paired")

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/indiv_assignments_10pops.png", width=8, height=4.5, res=300, units="in")

par(
  mar=c(5, 4, 2, 8), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14,
  xpd=NA
)
barplot(t(table3), ylab = 'Proportion assignment', col = cols, axisnames = FALSE, xaxt = 'n', las = 1)
axis(1, at=c(0.7, 1.9, 3.1, 4.3, 5.5), labels = FALSE)
mtext("Ingress site", side = 1, line = 3)
text(seq(0.7, 5.5, by=1.2), -0.11, srt = 45, labels = c("NJ", "DE", "VA", "NC", "SC"), xpd = TRUE)
text(seq(0.7, 5.5, by=1.2), 1.03, labels = c("n = 52", "n = 50", "n = 44", "n = 146", "n = 1"), cex = 0.8)

legend(6.1,1,
       #legend = levels(rev(locations)),
       legend = c("Pop A", "Pop B", "Pop C", "Pop D", "Pop E", "Pop F", "Pop G", "Pop H", "Pop I", "Pop J"),
       pch=22,
       col = 'black',
       pt.bg = rev(cols),
       title = 'North to south', text.font = 1)

dev.off()

# Which GAM-determined population is each individual most likely from?
bayenv.likelihoods.indivs <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/bayenv_likelihoods_293indivs_10pops.txt", header = TRUE) # 10 GAM pops

# Create logged dataset
bayenv.likelihoods.indivs2 <- as.data.frame(cbind(bayenv.likelihoods.indivs[,c(1:2)], log10(bayenv.likelihoods.indivs[,-c(1:2)])))

# Use either unlogged or logged version
locations <- factor(bayenv.likelihoods.indivs2$Place, c("Little Egg Inlet, NJ", "Roosevelt Inlet, DE", "York River, VA", "Beaufort, NC", "North Inlet, SC"))
most.like <- colnames(bayenv.likelihoods.indivs2[,-c(1:2)])[apply(bayenv.likelihoods.indivs2[,-c(1:2)],1, which.max)]
table <- t(table(most.like,locations))
most.like.geno <- apply(bayenv.likelihoods.indivs2[,-c(1:2)],1,max)
least.like.geno <- apply(bayenv.likelihoods.indivs2[,-c(1:2)],1,min)

# Find second most likely geno & pop
n <- 10
second.most.like <- vector()
second.most.like.geno <- vector()
for (i in 1:nrow(bayenv.likelihoods.indivs2[,-c(1:2)])){
  second.most.like.geno[i] <- sort(bayenv.likelihoods.indivs2[i,-c(1:2)],partial=n-1)[n-1]
  second.most.like[i] <- colnames(sort(bayenv.likelihoods.indivs2[i,-c(1:2)],partial=n-1)[n-1])
}

second.most.like.geno <- unlist(second.most.like.geno)
dif <- most.like.geno - second.most.like.geno

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/indiv_assignments_hists_10GAMpops.png", width=8, height=4.5, res=300, units="in")

par(
  mar=c(5, 5, 3, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14,
  xpd=NA,
  mfrow = c(1,2)
)

hist(dif, xlab = '', main = '')
mtext("Difference between two most \nlikely genotypes", 1, 3.7)
hist(most.like.geno - least.like.geno, xlab = '', main = "")
mtext("Difference between highest and \nlowest genotype likelihoods", 1, 3.7)

dev.off()


#### Use individual genotype likelihoods to calculate group likelihoods ####
# Read in fish data with clustering assignments
clusters <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/oto.gen.merged151.clusterdbytime.txt", header = TRUE) # clusters for each time period separately using all elements except Sn
clusters.bytime <- split(clusters, clusters$Period)
# clusters.bytime.early <- split(clusters.bytime$Early, clusters.bytime$Early$cluster)
# clusters.bytime.mid <- split(clusters.bytime$Mid, clusters.bytime$Mid$cluster)
# clusters.bytime.late <- split(clusters.bytime$Late, clusters.bytime$Late$cluster)

# Read in individual likelihoods and then subset to fish in early, middle and late time periods
bayenv.likelihoods.indivs <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/bayenv_likelihoods_293indivs_10pops.txt", header = TRUE) # these are individual likelihoods for 293 fish
bayenv.likelihoods.indivs.sub <- bayenv.likelihoods.indivs[bayenv.likelihoods.indivs$ID %in% clusters$PinskyID,] # subset to 151 fish in the clustering data
bayenv.likelihoods.early <- bayenv.likelihoods.indivs.sub[bayenv.likelihoods.indivs.sub$ID %in% clusters.bytime$Early$PinskyID,] # important that the fish are in same order
bayenv.likelihoods.middle <- bayenv.likelihoods.indivs.sub[bayenv.likelihoods.indivs.sub$ID %in% clusters.bytime$Mid$PinskyID,] # important that the fish are in same order
bayenv.likelihoods.late <- bayenv.likelihoods.indivs.sub[bayenv.likelihoods.indivs.sub$ID %in% clusters.bytime$Late$PinskyID,] # important that the fish are in same order

# Multiply across all fish within a cluster within a time period
# Early
bayenv.likelihoods.early$ID == clusters.bytime$Early$PinskyID # super important to check fish are in the same order
bayenv.likelihoods.early <- droplevels(bayenv.likelihoods.early) # drop unused factor levels; unmatched factor levels throw error
bayenv.likelihoods.early <- bayenv.likelihoods.early[order(match(bayenv.likelihoods.early$ID,clusters.bytime$Early$PinskyID)),] # order early fish by the clustered data & check order is the same again
as.character(bayenv.likelihoods.early$ID) == as.character(clusters.bytime$Early$PinskyID)

bayenv1.likes.early <- aggregate(bayenv.likelihoods.early$Pop1, by = list(clusters.bytime$Early$cluster), FUN = prod) # likelihoods for groupings from clustering analysis in early time period
bayenv2.likes.early <- aggregate(bayenv.likelihoods.early$Pop2, by = list(clusters.bytime$Early$cluster), FUN = prod)
bayenv3.likes.early <- aggregate(bayenv.likelihoods.early$Pop3, by = list(clusters.bytime$Early$cluster), FUN = prod)
bayenv4.likes.early <- aggregate(bayenv.likelihoods.early$Pop4, by = list(clusters.bytime$Early$cluster), FUN = prod)
bayenv5.likes.early <- aggregate(bayenv.likelihoods.early$Pop5, by = list(clusters.bytime$Early$cluster), FUN = prod) 
bayenv6.likes.early <- aggregate(bayenv.likelihoods.early$Pop6, by = list(clusters.bytime$Early$cluster), FUN = prod) 
bayenv7.likes.early <- aggregate(bayenv.likelihoods.early$Pop7, by = list(clusters.bytime$Early$cluster), FUN = prod) 
bayenv8.likes.early <- aggregate(bayenv.likelihoods.early$Pop8, by = list(clusters.bytime$Early$cluster), FUN = prod) 
bayenv9.likes.early <- aggregate(bayenv.likelihoods.early$Pop9, by = list(clusters.bytime$Early$cluster), FUN = prod) 
bayenv10.likes.early <- aggregate(bayenv.likelihoods.early$Pop10, by = list(clusters.bytime$Early$cluster), FUN = prod) 

bayenv.likelihoods.early.clustered <- as.data.frame(cbind(log10(as.numeric(bayenv1.likes.early[,2])), log10(as.numeric(bayenv2.likes.early[,2])), log10(as.numeric(bayenv3.likes.early[,2])), log10(as.numeric(bayenv4.likes.early[,2])), log10(as.numeric(bayenv5.likes.early[,2])), log10(as.numeric(bayenv6.likes.early[,2])), log10(as.numeric(bayenv7.likes.early[,2])), log10(as.numeric(bayenv8.likes.early[,2])), log10(as.numeric(bayenv9.likes.early[,2])), log10(as.numeric(bayenv10.likes.early[,2]))))
colnames(bayenv.likelihoods.early.clustered) <- c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5", "Pop6", "Pop7", "Pop8", "Pop9", "Pop10")
rownames(bayenv.likelihoods.early.clustered) <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6")
# write.table(bayenv.likelihoods.early.clustered, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_early_6clusters_10pops.txt", row.names = TRUE, col.names = TRUE)

most.like <- colnames(bayenv.likelihoods.early.clustered[apply(bayenv.likelihoods.early.clustered,1, which.max)]) # pick most likely origin region for each cluster

# Middle
bayenv.likelihoods.middle$ID == clusters.bytime$Mid$PinskyID # super important to check fish are in the same order
bayenv.likelihoods.middle <- droplevels(bayenv.likelihoods.middle) # drop unused factor levels; unmatched factor levels throw error
bayenv.likelihoods.middle <- bayenv.likelihoods.middle[order(match(bayenv.likelihoods.middle$ID,clusters.bytime$Mid$PinskyID)),] # order middle fish by the clustered data & check order is the same again
as.character(bayenv.likelihoods.middle$ID) == as.character(clusters.bytime$Mid$PinskyID)

bayenv1.likes.middle <- aggregate(bayenv.likelihoods.middle$Pop1, by = list(clusters.bytime$Mid$cluster), FUN = prod) # likelihoods for BayEnv groupings - same as above but assumes ~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/bayenv_likelihoods_indivs.txt has been read in
bayenv2.likes.middle <- aggregate(bayenv.likelihoods.middle$Pop2, by = list(clusters.bytime$Mid$cluster), FUN = prod)
bayenv3.likes.middle <- aggregate(bayenv.likelihoods.middle$Pop3, by = list(clusters.bytime$Mid$cluster), FUN = prod)
bayenv4.likes.middle <- aggregate(bayenv.likelihoods.middle$Pop4, by = list(clusters.bytime$Mid$cluster), FUN = prod)
bayenv5.likes.middle <- aggregate(bayenv.likelihoods.middle$Pop5, by = list(clusters.bytime$Mid$cluster), FUN = prod)
bayenv6.likes.middle <- aggregate(bayenv.likelihoods.middle$Pop6, by = list(clusters.bytime$Mid$cluster), FUN = prod) 
bayenv7.likes.middle <- aggregate(bayenv.likelihoods.middle$Pop7, by = list(clusters.bytime$Mid$cluster), FUN = prod)
bayenv8.likes.middle <- aggregate(bayenv.likelihoods.middle$Pop8, by = list(clusters.bytime$Mid$cluster), FUN = prod)
bayenv9.likes.middle <- aggregate(bayenv.likelihoods.middle$Pop9, by = list(clusters.bytime$Mid$cluster), FUN = prod)
bayenv10.likes.middle <- aggregate(bayenv.likelihoods.middle$Pop10, by = list(clusters.bytime$Mid$cluster), FUN = prod)

bayenv.likelihoods.middle.clustered <- as.data.frame(cbind(log10(as.numeric(bayenv1.likes.middle[,2])), log10(as.numeric(bayenv2.likes.middle[,2])), log10(as.numeric(bayenv3.likes.middle[,2])), log10(as.numeric(bayenv4.likes.middle[,2])), log10(as.numeric(bayenv5.likes.middle[,2])), log10(as.numeric(bayenv6.likes.middle[,2])), log10(as.numeric(bayenv7.likes.middle[,2])), log10(as.numeric(bayenv8.likes.middle[,2])), log10(as.numeric(bayenv9.likes.middle[,2])), log10(as.numeric(bayenv10.likes.middle[,2]))))
colnames(bayenv.likelihoods.middle.clustered) <- c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5", "Pop6", "Pop7", "Pop8", "Pop9", "Pop10")
rownames(bayenv.likelihoods.middle.clustered) <- c("cluster1", "cluster2")
# write.table(bayenv.likelihoods.middle.clustered, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_middle_2clusters_10pops.txt", row.names = TRUE, col.names = TRUE)

most.like <- colnames(bayenv.likelihoods.middle.clustered[apply(bayenv.likelihoods.middle.clustered,1, which.max)]) # pick most likely origin region for each cluster

# Late
bayenv.likelihoods.late$ID == clusters.bytime$Late$PinskyID # super important to check fish are in the same order
bayenv.likelihoods.late <- droplevels(bayenv.likelihoods.late) # drop unused factor levels; unmatched factor levels throw error
bayenv.likelihoods.late <- bayenv.likelihoods.late[order(match(bayenv.likelihoods.late$ID,clusters.bytime$Late$PinskyID)),] # order late fish by the clustered data & check order is the same again
as.character(bayenv.likelihoods.late$ID) == as.character(clusters.bytime$Late$PinskyID)

bayenv1.likes.late <- aggregate(bayenv.likelihoods.late$Pop1, by = list(clusters.bytime$Late$cluster), FUN = prod) # likelihoods for BayEnv groupings - same as above but assumes ~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/bayenv_likelihoods_indivs.txt has been read in
bayenv2.likes.late <- aggregate(bayenv.likelihoods.late$Pop2, by = list(clusters.bytime$Late$cluster), FUN = prod)
bayenv3.likes.late <- aggregate(bayenv.likelihoods.late$Pop3, by = list(clusters.bytime$Late$cluster), FUN = prod)
bayenv4.likes.late <- aggregate(bayenv.likelihoods.late$Pop4, by = list(clusters.bytime$Late$cluster), FUN = prod)
bayenv5.likes.late <- aggregate(bayenv.likelihoods.late$Pop5, by = list(clusters.bytime$Late$cluster), FUN = prod)
bayenv6.likes.late <- aggregate(bayenv.likelihoods.late$Pop6, by = list(clusters.bytime$Late$cluster), FUN = prod) 
bayenv7.likes.late <- aggregate(bayenv.likelihoods.late$Pop7, by = list(clusters.bytime$Late$cluster), FUN = prod)
bayenv8.likes.late <- aggregate(bayenv.likelihoods.late$Pop8, by = list(clusters.bytime$Late$cluster), FUN = prod)
bayenv9.likes.late <- aggregate(bayenv.likelihoods.late$Pop9, by = list(clusters.bytime$Late$cluster), FUN = prod)
bayenv10.likes.late <- aggregate(bayenv.likelihoods.late$Pop10, by = list(clusters.bytime$Late$cluster), FUN = prod)

bayenv.likelihoods.late.clustered <- as.data.frame(cbind(log10(as.numeric(bayenv1.likes.late[,2])), log10(as.numeric(bayenv2.likes.late[,2])), log10(as.numeric(bayenv3.likes.late[,2])), log10(as.numeric(bayenv4.likes.late[,2])), log10(as.numeric(bayenv5.likes.late[,2])), log10(as.numeric(bayenv6.likes.late[,2])), log10(as.numeric(bayenv7.likes.late[,2])), log10(as.numeric(bayenv8.likes.late[,2])), log10(as.numeric(bayenv9.likes.late[,2])), log10(as.numeric(bayenv10.likes.late[,2]))))
colnames(bayenv.likelihoods.late.clustered) <- c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5", "Pop6", "Pop7", "Pop8", "Pop9", "Pop10")
rownames(bayenv.likelihoods.late.clustered) <- c("cluster1", "cluster2", "cluster3")
# write.table(bayenv.likelihoods.late.clustered, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_late_3clusters_10pops.txt", row.names = TRUE, col.names = TRUE)

most.like <- colnames(bayenv.likelihoods.late.clustered[apply(bayenv.likelihoods.late.clustered,1, which.max)]) # pick most likely origin region for each cluster

#### Histograms of cluster likelihood differences, as for individuals ####
bayenv.likelihoods.early.clustered <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_early_6clusters_10pops.txt", header = TRUE)
bayenv.likelihoods.middle.clustered <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_middle_2clusters_10pops.txt", header = TRUE)
bayenv.likelihoods.late.clustered <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_late_3clusters_10pops.txt", header = TRUE)

cluster.likes <- rbind(bayenv.likelihoods.early.clustered, bayenv.likelihoods.middle.clustered, bayenv.likelihoods.late.clustered) # make sure colnames are the same and in same order

most.like.geno.clusters <- apply(cluster.likes,1,max)
least.like.geno.clusters <- apply(cluster.likes,1,min)

# Find second most likely geno & pop
n <- 10
second.most.like.clusters <- vector()
second.most.like.geno.clusters <- vector()
for (i in 1:nrow(cluster.likes)){
  second.most.like.geno.clusters[i] <- sort(cluster.likes[i,],partial=n-1)[n-1]
  second.most.like.clusters[i] <- colnames(sort(cluster.likes[i,],partial=n-1)[n-1])
}

second.most.like.geno.clusters <- unlist(second.most.like.geno.clusters)
dif.clusters <- most.like.geno.clusters - second.most.like.geno.clusters

hist(dif.clusters, xlab = '', main = '')
mtext("Difference between two most \nlikely cluster genotypes", 1, 3.7)
hist(most.like.geno.clusters - least.like.geno.clusters, xlab = '', main = "")
mtext("Difference between highest and \nlowest cluster genotype likelihoods", 1, 3.7)

###############################################################################
#### What is the meaning of distance? Can I convert this back to lat/long and plot this on a map? ####
library(geosphere)

envi <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/232envirowithdist.txt", header = TRUE)

# Distances from southern most point: 28.384993, -80.546297
# Can these be r? Then (x-h)^2 + (y-k)^2 = r^2?
# x = r*cosine(angle) --> x = Cx + (r*cosine((degrees * pi/180)))
# y = r*sin(angle) --> y = Cy + (r*sine((degrees * pi/180)))
new_latitude  <- latitude  + (dy / 6378137) * (180 / pi)
new_longitude <- longitude + (dx / 6378137) * (180 / pi) / cos(latitude * pi/180)
break.points <- c(158.5453,  343.7802,  529.0152,  714.2502,  899.4852, 1084.7202, 1269.9552, 1455.1902, 1640.4251, 1825.6601) * 1000 # These are km, multiply by 1000 to get meters

x.add <- vector()
y.add <- vector()

x.degs <- vector()
y.degs <- vector()

for (i in 1:length(break.points)){
  x.add[i] <- (r*cos((degrees * pi/180))) #r is distance
  y.add[i] <- (r*sin((degrees * pi/180)))
  
  x.degs[i] <- 28.384993  + ((y.add[i]) / 6378137) * (180 / pi)
  y.degs[i] <- -80.546297 + ((x.add[i]) / 6378137) * (180 / pi) / cos(28.384993 * pi/180)
}


# Plot circles with radius as distance?
library(geosphere)

distantCircle <- function(x, radius) {
  # Creation de 360 points distincts sur le cercle de centre
  # x et de rayon radius
  resul <- do.call("rbind", lapply(0:360, function(bearing) {
    res <- destPoint(p = x, b = bearing, d = radius)
    rownames(res) <- NULL
    return(data.frame(res))
  }))
  resul$dist <- radius / 1000
  return(resul)
}

distantCircle(x = c(0.0000001,89.9999999), radius = 1500*1000)


data <- t(data.frame(c(1,-80.546297,28.384993)))
circle <- make_circles(data, break.points[1]*1000)

library(sp)
longlat <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/adults_locations2.txt")
coordinates(longlat) <- c(1,2)
proj4string(longlat) <- CRS("+proj=longlat +datum=WGS84")
utm <- spTransform(longlat,CRS("+proj=utm +zone=18 +datum=WGS84"))
utm

#############################################################################
library(sp)
library(rgdal)
library(rgeos)
library(maps)
library(mapdata)
library(mapproj)
library(marmap)
library(geosphere)
library(ade4)

# Get bathymetry data
useast <- getNOAA.bathy(lon1 = -85, lon2 = -60, lat1 = 23, lat2 = 48, resolution = 1)

# Calculate buffer zones
the.points.latlong <- data.frame(long = -80.546297, lat = 28.384993) # southern point from which distance is calculated

the.points.sp <- SpatialPointsDataFrame(the.points.latlong, data.frame(ID=seq(1:nrow(the.points.latlong))), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

the.points.projected <- spTransform(the.points.sp, CRS( "+init=epsg:32618" ))  # projected southern point
# bring in the 10 distances
break.points <- c(158.5453,  343.7802,  529.0152,  714.2502,  899.4852, 1084.7202, 1269.9552, 1455.1902, 1640.4251, 1825.6601) * 1000 # These are km, multiply by 1000 to get meters

# For loop to generate buffer and then project those points
for (i in 1:length(break.points)){
  assign(paste0('the.circles.projected',i), gBuffer(the.points.projected, width=break.points[i], byid=TRUE))
}

projected.circles.names <- mget(ls(pattern = "the.circles.projected")) # searches through environment and pulls all objects with 'the.circles.projected'

for (i in 1:length(break.points)){
  assign(paste0('the.circles.sp',i), spTransform(projected.circles.names[[i]], CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")))
}

# And now plot!
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Maps/connectivity maps/adult_connectivity_regions.png", width=7, height=7.5, res=300, units="in")

par(
  mar=c(5, 7, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=18, # point size, which is the font size
  bg=NA
)

map("worldHires", c("us"), xlim=c(-85,-60), ylim=c(23,47.5), col="gray90", fill=TRUE) #plots the region of the USA that I want
# map("worldHires", c("us", "canada", "bah", "cub"), xlim=c(-85,-60), ylim=c(23,47.5), col="gray90", fill=TRUE, projection = 'albers', parameters = c(30,40), lforce = "s") #plots the region of the USA & canada that I want in albers projection
map("worldHires", c("us", "canada", "bah", "cub"), xlim=c(-85,-60), ylim=c(23,47.5), col="gray90", fill=TRUE) #plots the region of the USA & canada that I want
# map("state", xlim=c(-85,-60), ylim=c(23,47.5), add = TRUE, boundary=FALSE, col = 'black', projection = '') # plots US state boundaries in albers projection
map("state", xlim=c(-85,-60), ylim=c(23,47.5), add = TRUE, boundary=FALSE, col = 'black') # plots US state boundaries
# map.grid(cex=0.1 , col="grey30", labels = T) # for albers projection
title(xlab = "Longitude (°)", ylab = "Latitude (°)")
points(-80.546297, 28.384993, pch = "*", cex = 1.8, col = 'deepskyblue3')
# points(mapproject(-80.546297, 28.384993), col = 'red', pch = 15) # for albers projection
text(-76, 46.5,"CANADA")
text(-67.7, 43, "Gulf of Maine", cex=0.45)
text(-72, 35.236, "Cape Hatteras", cex=0.45)
lines(c(-68,-65.5), c(41, 40.5))
text(-63.2,40.2, "Georges Bank", cex=0.45)
lines(c(-75.54, -74.3), c(35.236, 35.236))
text(-65,33, expression(italic("Atlantic Ocean")))
axis(1, at=seq(-85,-60, by=5), labels=seq(-85,-60, by= 5)) # Plot axes
axis(2, at=seq(20,50, by = 5), labels=seq(20,50, by= 5), las = TRUE)
box()

# Figure out where 200 meter isobath and buffer zones intersect. Not necessary for final map
# plot(useast, deep = -200, shallow = -200, step = 0, lwd = 0.5, add = TRUE) # plot 200 meter isobath

# lines(mapproject(the.circles.sp1@polygons[[1]]@Polygons[[1]]@coords[,1], the.circles.sp1@polygons[[1]]@Polygons[[1]]@coords[,2]), col = 'red') # for albers projection

# lines(the.circles.sp1@polygons[[1]]@Polygons[[1]]@coords[,1], the.circles.sp1@polygons[[1]]@Polygons[[1]]@coords[,2], col = 'tomato') # J
# lines(the.circles.sp2@polygons[[1]]@Polygons[[1]]@coords[,1], the.circles.sp2@polygons[[1]]@Polygons[[1]]@coords[,2], col = 'tomato') # A
# lines(the.circles.sp3@polygons[[1]]@Polygons[[1]]@coords[,1], the.circles.sp3@polygons[[1]]@Polygons[[1]]@coords[,2], col = 'tomato') # I
# lines(the.circles.sp4@polygons[[1]]@Polygons[[1]]@coords[,1], the.circles.sp4@polygons[[1]]@Polygons[[1]]@coords[,2], col = 'tomato') # H
# lines(the.circles.sp5@polygons[[1]]@Polygons[[1]]@coords[,1], the.circles.sp5@polygons[[1]]@Polygons[[1]]@coords[,2], col = 'tomato') # G
# lines(the.circles.sp6@polygons[[1]]@Polygons[[1]]@coords[,1], the.circles.sp6@polygons[[1]]@Polygons[[1]]@coords[,2], col = 'tomato') #F
# lines(the.circles.sp7@polygons[[1]]@Polygons[[1]]@coords[,1], the.circles.sp7@polygons[[1]]@Polygons[[1]]@coords[,2], col = 'tomato') # E
# lines(the.circles.sp8@polygons[[1]]@Polygons[[1]]@coords[,1], the.circles.sp8@polygons[[1]]@Polygons[[1]]@coords[,2], col = 'tomato') # D
# lines(the.circles.sp9@polygons[[1]]@Polygons[[1]]@coords[,1], the.circles.sp9@polygons[[1]]@Polygons[[1]]@coords[,2], col = 'tomato') # C
# lines(the.circles.sp10@polygons[[1]]@Polygons[[1]]@coords[,1], the.circles.sp10@polygons[[1]]@Polygons[[1]]@coords[,2], col = 'tomato') # B

# Points where buffer zones intersect with 200 meter isobath
points(-80.2, 29.76001, pch = 18, col = 'deepskyblue3')
text(-79.6, 29.70, "J", cex=0.8)
points(-79.65, 31.36622, pch = 18, col = 'deepskyblue3')
text(-79.1, 31.2, "I", cex=0.8)
points(-78.25, 32.65, pch = 18, col = 'deepskyblue3')
text(-77.5, 32.2, "H", cex=0.8)
points(-76.5, 33.8, pch = 18, col = 'deepskyblue3')
text(-76.0, 33.4, "G", cex=0.8)
points(-75.16604, 34.95, pch = 18, col = 'deepskyblue3')
text(-74.5, 34.6, "F", cex=0.8)
points(-74.7, 36.7, pch = 18, col = 'deepskyblue3')
text(-74.1, 36.3, "E", cex=0.8)
points(-73.7, 38.2, pch = 18, col = 'deepskyblue3')
text(-72.95, 37.9, "D", cex=0.8)
points(-72.3, 39.5, pch = 18, col = 'deepskyblue3')
text(-71.6, 39.2, "C", cex=0.8)
points(-69.7, 40.05, pch = 18, col = 'deepskyblue3')
text(-69.1, 39.6, "B", cex=0.8)
points(-67.45, 40.5, pch = 18, col = 'deepskyblue3')
text(-66.9, 40.1, "A", cex=0.8)

# Label ingress sites
lines(c(-79.16799, -78), c(33.33229, 31.4)) # North Inlet, SC
text(-75.5, 31, "North Inlet, SC", cex=0.55)

lines(c(-76.67293, -75), c(34.71724, 33.5)) # Beaufort, NC
text(-72.7, 33.3, "Beaufort, NC", cex=0.55)

lines(c(-76.49952, -74.5), c(37.24831, 37.24831)) # York River, VA
lines(c(-74.5, -72.5), c(37.24831, 36.5))
text(-69.8, 36.2, "York River, VA", cex=0.55)

#lines(c(-75.42182, -72), c(37.93211, 37)) # Chincoteague, VA
#text(-68.8, 36.9, "Chincoteague, VA", cex=0.55)

lines(c(-75.15865, -73.8), c(38.79292, 38.79292)) # Roosevelt Inlet, DE
lines(c(-73.8, -72.3), c(38.79292, 38.2))
text(-69, 37.9, "Roosevelt Inlet, DE", cex=0.55)

lines(c(-74.32394, -72), c(39.50938, 38.65)) # Little Egg Inlet, NJ
lines(c(-72, -70.8), c(38.65, 38.65))
text(-67.6, 38.65, "Little Egg Inlet, NJ", cex=0.55)


# And plot larval ingress sites
larvalsamps <- read.csv("~/Documents/Graduate School/Rutgers/Summer Flounder/Maps/LarvalSamplingSites.csv") #my data for sampling sites, contains a column of "lat" and a column of "lon" with GPS points in decimal degrees
points(larvalsamps$lon[-4], larvalsamps$lat[-4], pch=21, col="black", cex=1.2, bg = 'tomato') #without Chincoteague

dev.off()
