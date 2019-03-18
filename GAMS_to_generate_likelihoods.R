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
  assign(paste0(colnames(adults232.freqs[i]), '.depth.gam'), mgcv::gam(adults232.freqs[,i] ~ s(depth), data = envi.ordered, method = 'REML', select = TRUE))
  assign(paste0(colnames(adults232.freqs[i]), '.btemp.gam'), mgcv::gam(adults232.freqs[,i] ~ s(b_temp), data = envi.ordered, method = 'REML', select = TRUE))
  assign(paste0(colnames(adults232.freqs[i]), '.bsalin.gam'), mgcv::gam(adults232.freqs[,i] ~ s(b_salin), data = envi.ordered, method = 'REML', select = TRUE))
}

# plot(GAM.object) as an object and this will save coordinates
par(mfrow = c(2,2))

gam.names <- mget(ls(pattern = ".btemp.gam")) # searches through environment and pulls all objects with '.btemp.gam'

# Plots all gams quickly!
for (i in 1:length(gam.names)){
  plot(gam.names[[i]], scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature', main = paste0(names(gam.names)[[i]]))
}

# plot(dDocent_Contig_2558_83.01.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_2558_83.03.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_8420_14.02.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_8420_14.01.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_15075_20.02.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_15075_20.04.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_19728_116.04.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_19728_116.02.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_31661_32.01.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_31661_32.02.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_35399_131.04.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_35399_131.02.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_37577_42.04.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_37577_42.02.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_38595_30.04.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_38595_30.02.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_42614_134.04.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_42614_134.01.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_54409_32.03.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
# plot(dDocent_Contig_54409_32.01.btemp.gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')

# Define new data for each environmental variable to predict
dist.new <- seq(min(envi.ordered$dist), max(envi.ordered$dist), length = 500)
depth.new <- seq(min(envi.ordered$depth), max(envi.ordered$depth), length = 500)
btemp.new <- seq(min(envi.ordered$b_temp), max(envi.ordered$b_temp), length = 500)
bsalin.new <- seq(min(envi.ordered$b_salin), max(envi.ordered$b_salin), length = 500)

# dist35gam.p <- predict(dist35gam, data.frame(dist = dist.new))
# depth35gam.p <- predict(depth35gam, data.frame(depth = depth.new))
# btemp35gam.p <- predict(btemp35gam, data.frame(b_temp = btemp.new))
# bsalin35gam.p <- predict(bsalin35gam, data.frame(b_salin = bsalin.new))

for (i in 1:length(names(gam.names))){
  assign(paste0(names(gam.names)[[i]], '.p'), predict.gam(gam.names[[i]], data.frame(b_temp = btemp.new)))
}

gam.names.p <- mget(ls(pattern = ".gam.p")) # searches through environment and pulls all objects with '.gam.p'

# dDocent_Contig_2558_83.01.btemp.gam.p <- predict(dDocent_Contig_2558_83.01.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_2558_83.03.btemp.gam.p <- predict(dDocent_Contig_2558_83.03.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_8420_14.02.btemp.gam.p <- predict(dDocent_Contig_8420_14.02.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_8420_14.01.btemp.gam.p <- predict(dDocent_Contig_8420_14.01.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_15075_20.02.btemp.gam.p <- predict(dDocent_Contig_15075_20.02.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_15075_20.04.btemp.gam.p <- predict(dDocent_Contig_15075_20.04.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_19728_116.04.btemp.gam.p <- predict(dDocent_Contig_19728_116.04.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_19728_116.02.btemp.gam.p <- predict(dDocent_Contig_19728_116.02.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_31661_32.01.btemp.gam.p <- predict(dDocent_Contig_31661_32.01.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_31661_32.02.btemp.gam.p <- predict(dDocent_Contig_31661_32.02.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_35399_131.04.btemp.gam.p <- predict(dDocent_Contig_35399_131.04.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_35399_131.02.btemp.gam.p <- predict(dDocent_Contig_35399_131.02.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_37577_42.04.btemp.gam.p <- predict(dDocent_Contig_37577_42.04.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_37577_42.02.btemp.gam.p <- predict(dDocent_Contig_37577_42.02.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_38595_30.04.btemp.gam.p <- predict(dDocent_Contig_38595_30.04.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_38595_30.02.btemp.gam.p <- predict(dDocent_Contig_38595_30.02.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_42614_134.04.btemp.gam.p <- predict(dDocent_Contig_42614_134.04.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_42614_134.01.btemp.gam.p <- predict(dDocent_Contig_42614_134.01.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_54409_32.03.btemp.gam.p <- predict(dDocent_Contig_54409_32.03.btemp.gam, data.frame(b_temp = btemp.new))
# dDocent_Contig_54409_32.01.btemp.gam.p <- predict(dDocent_Contig_54409_32.01.btemp.gam, data.frame(b_temp = btemp.new))

# Plot newly predicted allele frequency
# plot(dist35gam.p ~ dist.new)
# plot(dist35gam.p ~ envi.ordered$dist)

# Plots all predicted gams quickly!
for (i in 1:length(gam.names.p)){
  plot(gam.names.p[[i]] ~ btemp.new, main = paste0(names(gam.names.p)[[i]]))
}

# plot(dDocent_Contig_2558_83.01.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_2558_83.03.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_8420_14.02.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_8420_14.01.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_15075_20.02.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_15075_20.04.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_19728_116.04.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_19728_116.02.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_31661_32.01.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_31661_32.02.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_35399_131.04.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_35399_131.02.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_37577_42.04.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_37577_42.02.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_38595_30.04.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_38595_30.02.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_42614_134.04.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_42614_134.01.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_54409_32.03.btemp.gam.p ~ btemp.new)
# plot(dDocent_Contig_54409_32.01.btemp.gam.p ~ btemp.new)

# Break 500 bottom temps into 10 groups
# break_points <- cut(btemp.new, breaks =10)

break.points <- btemp.new[c(25, 75, 125, 175, 225, 275, 325, 375, 425, 475)]

# Figure out corresponding allele frequency
# cbind temp and predicted allele frequency
for (i in 1:length(names(gam.names))){
  assign(paste0(names(gam.names)[[i]], '.joined'), cbind(btemp.new, gam.names.p[[i]]))
}

gam.joined <- mget(ls(pattern = ".joined")) # searches through environment and pulls all objects with '.joined'


# dDocent_Contig_2558_83.01.btemp.gam10 <- as.data.frame(cbind(btemp.new, dDocent_Contig_2558_83.01.btemp.gam.p))
# dDocent_Contig_2558_83.03.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_2558_83.03.btemp.gam.p)
# dDocent_Contig_8420_14.02.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_8420_14.02.btemp.gam.p)
# dDocent_Contig_8420_14.01.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_8420_14.01.btemp.gam.p)
# dDocent_Contig_15075_20.02.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_15075_20.02.btemp.gam.p)
# dDocent_Contig_15075_20.04.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_15075_20.04.btemp.gam.p)
# dDocent_Contig_19728_116.04.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_19728_116.04.btemp.gam.p)
# dDocent_Contig_19728_116.02.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_19728_116.02.btemp.gam.p)
# dDocent_Contig_31661_32.01.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_31661_32.01.btemp.gam.p)
# dDocent_Contig_31661_32.02.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_31661_32.02.btemp.gam.p)
# dDocent_Contig_35399_131.04.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_35399_131.04.btemp.gam.p)
# dDocent_Contig_35399_131.02.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_35399_131.02.btemp.gam.p)
# dDocent_Contig_37577_42.04.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_37577_42.04.btemp.gam.p)
# dDocent_Contig_37577_42.02.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_37577_42.02.btemp.gam.p)
# dDocent_Contig_38595_30.04.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_38595_30.04.btemp.gam.p)
# dDocent_Contig_38595_30.02.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_38595_30.02.btemp.gam.p)
# dDocent_Contig_42614_134.04.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_42614_134.04.btemp.gam.p)
# dDocent_Contig_42614_134.01.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_42614_134.01.btemp.gam.p)
# dDocent_Contig_54409_32.03.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_54409_32.03.btemp.gam.p)
# dDocent_Contig_54409_32.01.btemp.gam10 <- cbind(btemp.new, dDocent_Contig_54409_32.01.btemp.gam.p)

allele.freqs10 <- list()
for (i in 1:length(gam.joined)){
  allele.freqs10[[i]] <- gam.joined[[i]][gam.joined[[i]][,1] %in% break.points,][,2]
}

allele.freqs10.df <- data.frame(t(matrix(unlist(allele.freqs10), nrow=length(allele.freqs10), byrow=T)))
colnames(allele.freqs10.df) <- names(gam.names.p)

# Reorder 10 allele frequencies to be in same order as larvae

# Use these allele frequencies to calculate larval likelihoods
