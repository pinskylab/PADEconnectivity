#### This script aims to calculate and plot IBD for larval summer flounder captured between 2008-2012 ####

library(tidyr)
library(plyr)
library("hierfstat")
library(geosphere)

#### Read in the data ####
# Read in larval genetic data of 293 larvae/1882 biallelic loci spanning 1989-2012
geno <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/larvs_allelecounts.txt", sep="\t")

# Fix genetic rownames so that they match w/ metadata
fishnames.split <- as.data.frame(do.call(rbind, strsplit(as.character(rownames(geno)), 'L', fixed = TRUE)))
fishnames.split2 <- separate(fishnames.split, V1, c("name1", "name2", "name3", "name4"),sep = c(4,5,7))
newnames <- paste(fishnames.split2$name1, fishnames.split2$name3, fishnames.split2$name2, fishnames.split2$name4, sep = '')

rownames(geno) <- newnames

# Read in larval metadata
meta <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/masterPADElarvae.txt", sep="\t", header = TRUE)

#### Subset genetic data to those from 2008-2012 ####
late <- c('2008', '2009', '2010', '2011', '2012')
meta.late <- meta[meta$Year %in% late,]

geno.late <- geno[rownames(geno) %in% meta.late$PinskyID,] # Subset genetic data
meta.late.sub <- meta.late[meta.late$PinskyID %in% rownames(geno.late),] # metadata for 2008-2012 fish having genetic data

meta.late.sub$PinskyID == row.names(geno.late) # they are in the same order

#### Need to convert from allele counts to alleles on my own ####
alleles <- colnames(geno.late)

# Isolate alleles
alleles.split <- as.data.frame(do.call(rbind, strsplit(as.character(alleles), '.', fixed = TRUE)))
alleles.split2 <- separate(alleles.split, V2, c("name1", "name2"),sep = 1)

# Split data into odd and even dataframes
even_indexes<-seq(2,3764,2)
odd_indexes<-seq(1,3763,2)

odds <- data.frame(geno.late[,odd_indexes]) # 146 x 1882
evens <- data.frame(geno.late[,even_indexes]) # 146 x 1882

# Split the alleles up too
odd.alleles <- alleles.split2$name2[odd_indexes]
even.alleles <- alleles.split2$name2[even_indexes]
alleles1882 <- t(cbind(odd.alleles, even.alleles))
colnames(alleles1882) <- as.character(alleles.split2$V1[odd_indexes])

# Build matrix of alleles
odd.allele.table <- matrix(nrow = 146, ncol = 1882, dimnames = list(rownames(odds), colnames(alleles1882)))
even.allele.table <- matrix(nrow = 146, ncol = 1882, dimnames = list(rownames(odds), colnames(alleles1882)))

for(j in 1:length(rownames(odds))){

for (i in 1:length(colnames(odds))){
  if (is.na(odds[j,i])) {
    odd.allele.table[j,i] <- NA
    even.allele.table[j,i] <- NA
  } else if (odds[j,i] == 1) {
    odd.allele.table[j,i] <- paste(alleles1882['odd.alleles', i])
    even.allele.table[j,i] <- paste(alleles1882['even.alleles', i])
  } else if (odds[j,i] == 2) {
    odd.allele.table[j,i] <- paste(alleles1882['odd.alleles',i])
    even.allele.table[j,i] <- paste(alleles1882['odd.alleles',i])
  } else {
    odd.allele.table[j,i] <- paste(alleles1882['even.alleles', i])
    even.allele.table[j,i] <- paste(alleles1882['even.alleles', i])
  }
}

}

# Test to make sure this code converted between allele counts and alleles correctly
allele.counts.after <- matrix(nrow = 146, ncol = 2)
allele.counts.before <- matrix(nrow = 146, ncol = 2)

for (j in 1:length(rownames(odds))){
allele.counts.after[j,] <- sort(as.numeric(t(count(c(odd.allele.table[,j], even.allele.table[,j])))['freq',1:2]), decreasing = TRUE)

allele.counts.before[j,] <- sort(cbind(sum(evens[,j], na.rm = TRUE),sum(odds[,j], na.rm = TRUE)), decreasing = TRUE)

}

allele.counts.before == allele.counts.after # yay!

# Now paste each value in one dataframe with its corresponding in the other
s <- 1:length(colnames(odd.allele.table))
combo <- data.frame(matrix(nrow = 146, ncol = 1882))
for (i in s){
  combo[,i] <-paste(odd.allele.table[,i], even.allele.table[,i], sep = '')
}

dim(combo) #146 x 1882

combo[] <- lapply(combo, function(x) as.numeric(as.character(x)))# Convert to numeric, gives warning because replaces character 'NANA' with NA
colnames(combo) <- colnames(odd.allele.table)

#### Create a vector with integers representing capture location ####
table(meta.late.sub$Place)
population <- gsub('Little Egg Inlet, NJ', 1, meta.late.sub$Place)
population <- gsub('Roosevelt Inlet, DE', 2, population)
population <- gsub('York River, VA', 3, population)
population <- gsub('Beaufort, NC', 4, population)
population <- gsub('North Inlet, SC', 5, population)

#### Combine the population numbers with the allele data ####
combo2 <- cbind(population, combo)
dim(combo2) #146 x 1883

#### Calculate geographic distance ####
geo <- unique(meta.late.sub[c('Lon', 'Lat', 'Place')])
geo.sorted <- geo[order(geo$Lat, decreasing = TRUE),]

geodist <- distm(geo.sorted[,1:2], fun = distCosine)
geodistkm <- geodist * 0.001

#### FST calculation ####
gen <- pairwise.WCfst(combo2,diploid=TRUE)
genet.dist(combo2, method = 'WC84')
fst <- gen/(1-gen)

# Plot
plot(fst[lower.tri(fst)] ~ geodistkm[lower.tri(geodistkm)], xlab = 'Geographic distance (km)', ylab = 'Fst/(1-Fst)')
linreg <- lm(fst[lower.tri(fst)] ~ geodistkm[lower.tri(geodistkm)])
summary(linreg)
abline(linreg, col = "gray65")





#### Build the file in STRUCTURE format so that this can be converted to Geneland format (which is closer to what I want for FST calculations) ####
geno.late.pop <- cbind(rownames(geno.late), population, geno.late)
write.table(geno.late.pop, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/gen_input_recentIBD.txt", sep="\t", col.names = FALSE, row.names = FALSE)
# Get rid of ""'s and change NA to -9

