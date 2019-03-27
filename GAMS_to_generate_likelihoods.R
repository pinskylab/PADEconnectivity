# Bringing in adult genotype data and environmental data
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis")

library(ade4)
library(adegenet)
library(devtools)
library("hierfstat")
library(pegas)
library(fields)

# Reading in SNP data file containing only the first SNP at each locus
adults <- read.structure("structure_input_Nov_11_2015.str",
                         n.ind = 241, n.loc = 1137, col.lab = 1, col.pop = 2, row.marknames = 1, 
                         onerowperind = FALSE)

which(adults@loc.n.all > 2) # which snps have more than 2 alelles?

# adults.nonas <- scaleGen(adults, center = TRUE, scale = FALSE, NA.method = "mean") # filling in NAs with allele frequencies and centering
adults.nonas <- scaleGen(adults, center = FALSE, scale = FALSE, NA.method = "mean") # based on Meirimans 2015 and others, it's okay to perform rda on allele frequencies
hist(adults.nonas)
sum(is.na(adults.nonas)) #0 All NAs have been replaced with means

exclu_names <- c("PADE_14230L1439",
                 "PADE_14231L1440",
                 "PADE_14232L1529",
                 "PADE_14233L1588",
                 "PADE_14234L1441",
                 "PADE_14235L1442",
                 "PADE_14236L1530",
                 "PADE_14237L1531",
                 "PADE_14238L1532") # Creating a list of IDs to exclude from the genetic matrix

adults.nonas.232 <- adults.nonas[ ! rownames(adults.nonas) %in% exclu_names, ] # This is genotype data for 232 fish where NAs have been replaced by mean allele frequencies
rownames(adults.nonas.232)
adults.nonas.232.2274 <- adults.nonas.232[, -c(56, 96, 734, 2115)] # SNP 28, 47, 366 and 1056 have 3 alleles. All the rest have 2. Need to remove 3rd allele for SNP 28, 47, 366 and 1056. Removing one with lowest count
dim(adults.nonas.232.2274)
adults.nonas.232.1137 <- adults.nonas.232.2274[, seq(1, ncol(adults.nonas.232.2274),by = 2)] # including every other column
# adults.nonas.232.1137 <- adults.nonas.232.2274[, -seq(1, ncol(adults.nonas.232.2274),by = 2)] # excluding every other column
dim(adults.nonas.232.1137)


#### Reads in a file of 232 fish and their allele counts at the 10 candidates present in larvae. Allele frequencies need to be calculated. ####
adults232 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.counts10.txt", header = TRUE)
rownames(adults232) <- adults232$PinskyID

# Replace NAs with column mean
for(i in 1:ncol(adults232[,-c(1:2)])){
  adults232[,-c(1:2)][is.na(adults232[,-c(1:2)][,i]), i] <- mean(adults232[,-c(1:2)][,i], na.rm = TRUE)
}

colSums(is.na(adults232)) # Check there are no zeros in any of the columns

# Convert from count to frequency for allele columns
adults232.freqs <- adults232[,-c(1:2)]/2

envi <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/232envirowithdist.txt", header = TRUE)
envi.ordered <- envi[with(envi, order(PinskyID)),] 
# as.character(envi.ordered[,2]) == rownames(adults.nonas.232.1137)
as.character(envi.ordered[,2]) == rownames(adults232.freqs) # environmental data is in same order as genetic data?
rownames(envi.ordered) <- envi.ordered[,"PinskyID"] # so that environmental variable units in GAM plots won't be standardized

# Look up SNP number
# candidates <- c(35, 125 , 214, 291, 396, 442, 499, 524, 542, 609, 615, 743, 825, 919, 990)

# Run a GAM against bottom temperature? Or the associated environmental variable? Or all environmental variables?
# For loop to fit GAM between all 10 // 15 candidates and each 4 enviromental variables
library(mgcv)

# for(i in candidates){
#   assign(paste0('dist',i,'gam'), mgcv::gam(adults.nonas.232.1137[,i] ~ s(dist), data = envi.ordered, method = 'REML', select = TRUE))
#   assign(paste0('depth',i,'gam'), mgcv::gam(adults.nonas.232.1137[,i] ~ s(depth), data = envi.ordered, method = 'REML', select = TRUE))
#   assign(paste0('btemp',i,'gam'), mgcv::gam(adults.nonas.232.1137[,i] ~ s(b_temp), data = envi.ordered, method = 'REML', select = TRUE))
#   assign(paste0('bsalin',i,'gam'), mgcv::gam(adults.nonas.232.1137[,i] ~ s(b_salin), data = envi.ordered, method = 'REML', select = TRUE))
# }

# 10 loci
for(i in 1:ncol(adults232.freqs)){
  assign(paste0(colnames(adults232.freqs[i]), '.dist.gam'), mgcv::gam(adults232.freqs[,i] ~ s(dist), data = envi.ordered, method = 'REML', select = TRUE))
  # assign(paste0(colnames(adults232.freqs[i]), '.depth.gam'), mgcv::gam(adults232.freqs[,i] ~ s(depth), data = envi.ordered, method = 'REML', select = TRUE))
  # assign(paste0(colnames(adults232.freqs[i]), '.btemp.gam'), mgcv::gam(adults232.freqs[,i] ~ s(b_temp), data = envi.ordered, method = 'REML', select = TRUE))
  # assign(paste0(colnames(adults232.freqs[i]), '.bsalin.gam'), mgcv::gam(adults232.freqs[,i] ~ s(b_salin), data = envi.ordered, method = 'REML', select = TRUE))
}

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
  assign(paste0(names(gam.names)[[i]], '.p'), predict.gam(gam.names[[i]], data.frame(dist = dist.new)))
}

gam.names.p <- mget(ls(pattern = ".gam.p")) # searches through environment and pulls all objects with '.gam.p'

# Round any predicted allele frequency that is above 1 down to 0.999, or below 0 to 0.001. Allele frequencies of exactly 1 or zero will lead to null alleles and values of infinity. Cornuet (1999) says null allele frequencies can be replaced by a small value.
# Alleles with predicted allele frequencies > 1: dDocent_Contig_15075_20.02.dist.gam.p
# Alleles with predicted allele frequencies < 0: dDocent_Contig_15075_20.04.dist.gam.p
gam.names.p$dDocent_Contig_15075_20.02.dist.gam.p[gam.names.p$dDocent_Contig_15075_20.02.dist.gam.p > 1] <- 0.999
gam.names.p$dDocent_Contig_15075_20.04.dist.gam.p[gam.names.p$dDocent_Contig_15075_20.04.dist.gam.p < 0] <- 0.001


# Plot newly predicted allele frequency
# plot(dist35gam.p ~ dist.new)
# plot(dist35gam.p ~ envi.ordered$dist)

# Plots all predicted gams quickly!
for (i in 1:length(gam.names.p)){
  plot(gam.names.p[[i]] ~ dist.new, main = paste0(names(gam.names.p)[[i]]))
}

# Break 500 bottom temps into 10 groups
# break_points <- cut(btemp.new, breaks =10)
# break.points <- btemp.new[c(25, 75, 125, 175, 225, 275, 325, 375, 425, 475)]

# break_points <- cut(dist.new, breaks =10)
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
write.table(allele.freqs10.df, '~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.freqs.10pops.txt', row.names = FALSE, col.names = TRUE)

# Fix names for columns
allelenames.split <- as.data.frame(do.call(rbind, strsplit(as.character(names(gam.names)), '.', fixed = TRUE)))
allelenames.split2 <- as.factor(paste(allelenames.split$V1, '.', allelenames.split$V2, sep = ''))

colnames(allele.freqs10.df) <- allelenames.split2

# # Keep only one allele per locus in adults
odds <- seq(1,20,2) # odd indicies to keep
allele.freqs10.df.odds <- allele.freqs10.df[,odds]

#### Use these allele frequencies to calculate larval likelihoods for all 293 fish ####
# Read in dataset containing outlier loci
gen.larvae.outs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/masterPADElarvae.txt', header = TRUE, sep = "\t")

# Keep only larvae with genetic assignments
gen.larvae.outs2 <- gen.larvae.outs[which(is.na(gen.larvae.outs$assignment) == FALSE), ]

indiv.allele.counts <- gen.larvae.outs2[,c(16:35)] # Genetic data for all 293 larvae

indiv.allele.counts.odds <- indiv.allele.counts[, colnames(indiv.allele.counts) %in% colnames(allele.freqs10.df.odds)] # subset alleles in larvae to those in adults
indiv.allele.counts.odds[is.na(indiv.allele.counts.odds)] <- 9 # replace NA's with 9's to make the ifelse statements easier

# Reorder 10 adult allele frequencies to be in same order as larvae
allele.freqs10.df.odds.ordered <- allele.freqs10.df.odds[names(indiv.allele.counts.odds)]
write.table(allele.freqs10.df.odds.ordered, '~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.freqs.10pops.odd.ordered.txt', row.names = FALSE, col.names = TRUE)

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
barplot(t(table3), ylab = 'Proportion assignment', col = cols, axisnames = FALSE, xaxt = 'n')
axis(1, at=c(0.7, 1.9, 3.1, 4.3, 5.5), labels = FALSE)
mtext("Ingress site", side = 1, line = 3.7)
text(seq(0.7, 5.5, by=1.2), -0.135, srt = 45, labels = c("RUMFS", "Roosevelt", "York", "NC", "SC"), xpd = TRUE)
text(seq(0.7, 5.5, by=1.2), 1.03, labels = c("n = 52", "n = 50", "n = 44", "n = 146", "n = 1"), cex = 0.8)

legend(6.1,1,
       #legend = levels(rev(locations)),
       legend = c("Pop A", "Pop B", "Pop C", "Pop D", "Pop E", "Pop F", "Pop G", "Pop H", "Pop I", "Pop J"),
       pch=22,
       col = 'black',
       pt.bg = rev(cols),
       title = 'North to south', text.font = 1)

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


