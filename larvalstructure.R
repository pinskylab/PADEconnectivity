setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files")

library(ade4)
library(adegenet)
library("hierfstat")
library(pegas)


# Reading in the allele frequency data for adults and larvae
full <- read.structure("structure_input_Mar14_2017_528_1904.str",
                       n.ind = 528, n.loc = 1904, col.lab = 1, col.pop = 2, row.marknames = 1,
                       onerowperind = FALSE)

# How much missing data for each allele using full dataset with larvae and adults?
result_full <- vector("integer", length(full$tab[1,]))
for (i in 1:length(full$tab[1,])){
  result_full[i] <- sum(is.na(full$tab[,i]))
}
hist(result_full)
boxplot(result_full)
summary(result_full) # worst allele has 26 individuals missing data

# How much is missing for each fish?
result_fish <- vector("integer", length(full$tab[,1]))
for (i in 1:length(full$tab[,1])){
  result_fish[i] <- sum(is.na(full$tab[i,]))
}

hist(result_fish)
summary(result_fish) # worst fish is missing 33.6% of the data
boxplot(result_fish)

which(result_fish > 500) #row index of fish (most are adults): 29  94 232 234 244 245 259 269 273 275 276 277 278 279 280 281 282 283 284 285 286 290 291 479

#### Multiple alleles ####
multiallelic <- names(which(full@loc.n.all > 2))
multiallelic <- paste(multiallelic, ".", sep = "")
# multiallelic_indices <- c(289:291,384:386,393:395)

multiallelic_data <- full@tab[, grep(paste(multiallelic, "{3,4}", sep = "", collapse = "|"), colnames(full@tab))] # mostly works, but still includes SNP_145*
multiallelic_data <- multiallelic_data[,-c(59:78)]
multiallelic_counts <- colSums(multiallelic_data, na.rm = TRUE)

multiallelic_names <- names(multiallelic_counts) #67
biallelic_only <- full.dataframe[,!colnames(full.dataframe) %in% multiallelic_names] # 528x3764 (before: 3831-67 = 3764)

# Turn it back into a genind
full_bi <- as.genind(biallelic_only) # both adults and larvae
full_bi@pop <- full@pop

# Need to be able to categorize by region, so reading in larval database so that I can match larvae names to get location
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Larvae")
larvs <- read.csv("Larvae Sampling Database.csv", header=TRUE)

#### Getting rid of the adult allele frequencies and fixing the larval names####
# Need to get rid of the adults so I can only look at the larvae
dim(full_bi@tab)
full.freqs <- as.data.frame(full_bi@tab) #biallelic
full.freqs <- as.data.frame(full@tab) #multiallelic
dim(full.freqs) #528 x 3764; 528 x 3831

adults <- rownames(full.freqs)[205:439] #528 fish total

# Biallelic
larvs.freqs <- full.freqs[!(rownames(full.freqs) %in% adults),]
dim(larvs.freqs) #293 x 3764
larvs.freqs2 <- as.genind(larvs.freqs)

# Multiallelic
larvs.freqs <- full.freqs[!(rownames(full.freqs) %in% adults),]
dim(larvs.freqs) #293 x 3831
larvs.freqs2 <- as.genind(larvs.freqs)

# Removing the ligation ID in the names and sticking the names back
larv.names.split <- do.call(rbind, strsplit(rownames(larvs.freqs2@tab), 'L')) # Gets rid of the ligation ID so that I can match fish ID
rownames(larvs.freqs2@tab) <- larv.names.split[,1]

# Removing the adult populations and putting pop back in
nope <- pop(full_bi)[pop(full)==4] #biallelic
pop <- pop(full_bi)[!as.vector(pop(full_bi)) %in% nope]
nope <- pop(full)[pop(full)==4] #multiallelic
pop <- pop(full)[!as.vector(pop(full)) %in% nope]
population <- as.vector(pop)
larvs.freqs2@pop <- as.factor(population)

#### The names for the larval database are all messed up ####
# The names are all messed up because of the bioinformatics process...
names1 <- do.call(rbind, strsplit(as.character(larvs$ID..), '_'))
names2 <- do.call(rbind, strsplit(as.character(names1[,1]), 'E'))
correctnames <- paste(names2[,1], names2[,2], sep = 'E_')
correctnames2 <- paste(correctnames, names1[,2], sep = '')

# Now put back the fixed names in the larval database and the split names into the allele frequency database
larvs$ID.. <- correctnames2

#### Subsetting the larval database to only those who have genetic data ####
larvs.sub <- larvs[larvs$ID.. %in% larv.names.split[,1],]
dim(larvs.sub) #293 x 20

# Extracting the place for all the larvae with genetic data
place <- larvs.sub[,c("ID..", "Place")]
ordered.place <- place[order(place[,1]),]
dim(ordered.place) # 293 x 2

#### Creating a vector of north vs south based on capture location ####
place_only <- ordered.place[,-1]
place_only <- gsub('Beaufort, NC', 'south', place_only)
place_only <- gsub('North Inlet, SC', 'south', place_only)
place_only <- gsub('Little Egg Inlet, NJ', 'north', place_only)
place_only <- gsub('York River, VA', 'north', place_only)
place_only <- gsub('Roosevelt Inlet, DE', 'north', place_only)

# Put into genind object for a PCA
pop(larvs.freqs2) <- place_only

#### Calculate missing data on individual and allele level ####
# How much missing data for each allele?
result <- vector("integer", length(larvs.freqs2$tab[1,]))
for (i in 1:length(larvs.freqs2$tab[1,])){
  result[i] <- sum(is.na(larvs.freqs2$tab[,i]))
}
hist(result)
boxplot(result)
summary(result) # max number of individuals missing data is 22
which(result > 10)

# How much is missing for each larvae?
result_indiv <- vector("integer", length(larvs.freqs2$tab[,1]))
for (i in 1:length(larvs.freqs2$tab[,1])){
  result_indiv[i] <- sum(is.na(larvs.freqs2$tab[i,]))
}

hist(result_indiv)
summary(result_indiv) # max % missing data for an individual is 33.2%
boxplot(result_indiv)
which(result_indiv > 500) #3 larvae missing 500 alleles or more
larvs.freqs2@tab[c(9,94,244),1:2] # which fish are missing these alleles?

#### PCA on larvae only####
sum(is.na(larvs.freqs2$tab)) #9424; try 3: 9424/(293*3764); multiallelic = 9554/(293*3831)
X <- scaleGen(larvs.freqs2, NA.method = "mean") # 4 alleles that aren't variable are being cut out; variance = 0 makes it hard to scale
dim(X) #293 x 3760 = biallelic; 293 x 3827 = multiallelic
class (X)

# make PCA
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

pca1

# Plotting PC1 and PC2
s.label(pca1$li)
title("PCA of summer flounder dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

s.class(pca1$li, pop(larvs.freqs2))
title("PCA of summer flounder dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

col <- azur(2)
s.class(pca1$li, pop(larvs.freqs2), xax=1,yax=2, col = transp(col,0.6), axesell=FALSE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE)

eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent [1:3]

### To make a nice plot of the PCA broken down by region ###
pop(larvs.freqs2) <- place_only
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/results/larval_pca_regionalpops_1904loci.png", width=8, height=7, res=300, units="in")

par(
  mar=c(9, 8, 3, 8), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA,
  bty = 'n'
)

col <- c("#0072B2", "#D55E00")
# col <- c("#969696", "#000000")
s.class(pca1$li, pop(larvs.freqs2), xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, xlim = c(-50,30), ylim = c(-30,28), clabel = 0)
axis(1, at=seq(-40,20, by=10), labels=seq(-40,20, by= 10), line = 1.5)
axis(2, at=seq(-20,20, by = 10), labels=seq(-20,20, by= 10), line = 3, las = 2)
mtext("PC1 (0.66%)", side = 1, line = 4)
mtext("PC2 (0.65%)", side = 2, line = 5.5)

legend(-36, -8,
       legend=c("North (n = 146)", "South (n = 147)"),
       pch=c(19, 19),
       col=col,
       bty = "n",
       y.intersp = 1)

legend(-40.7,21,
       legend = expression("A"),
       bty = "n",
       cex = 1.3)

dev.off()

### To make a nice plot of the PCA broken down by time period ###
larvs.freqs2@pop <- as.factor(population) # Make sure the pop content refers to time periods and not regions
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/results/larval_pca_temppops_1904loci.png", width=8, height=7, res=300, units="in")

par(
  mar=c(9, 8, 3, 8), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA,
  bty = 'n'
)

col <- azur(3)
# col <- greypal(4)
# col <- c("#C6C6C6", "#686868", "#000000")
s.class(pca1$li, pop(larvs.freqs2), xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, xlim = c(-50,30), ylim = c(-30,28), clabel = 0)
axis(1, at=seq(-40,20, by=10), labels=seq(-40,20, by= 10), line = 1.5)
axis(2, at=seq(-20,20, by = 10), labels=seq(-20,20, by= 10), line = 3, las = 2)
mtext("PC1 (0.66%)", side = 1, line = 4)
mtext("PC2 (0.65%)", side = 2, line = 5.5)

legend(-36, -7,
       legend=c("1989-1993 (n = 60)", "1998-2002 (n = 87)", "2008-2012 (n = 146)"),
       pch=c(19, 19, 19),
       col=col,
       # col = c("#C6C6C6", "#686868", "#000000"),
       bty = "n",
       y.intersp = 1)

legend(-40.7,21,
       legend = expression("B"),
       bty = "n",
       cex = 1.3)
dev.off()

# Need to match PC coordinates with the location and period so that I can plot in different colors/draw ellipses
rownames(pca1$li) == ordered.place[,1] #they match
larvalpca.place <- cbind(pca1$li, ordered.place)
write.table(larvalpca.place, file = '~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/larvalpca.place3.txt', col.names = FALSE)

### To make a nice plot of the PCA using larvae###
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/Larvae/larval_pca.png", width=8, height=7, res=300, units="in")

par(
  mar=c(5, 4, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA
)

# Break PCA down by year/period, plotting PC1 vs PC2 (293 larvae, 3760 alleles)
par(mar=c(4.1, 4.1, 10, 2), xpd=TRUE,
    mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
    tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
    cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
    ps=14, # point size, which is the font size
    bg=NA)

plot(pca1$li[c(1:6, 15:28, 40:48, 265:271, 290:291),1], pca1$li[c(1:6, 15:28, 40:48, 265:271, 290:291),2], col = "steelblue3", pch = 6, xlab = "PC1 (0.67%)", ylab = "PC2 (0.65%)", xlim = c(-40,20), ylim = c(-20,15)) # 38 fish (mid period/north)
points(pca1$li[c(7:14, 29:39, 49:58, 272:289, 292:293),1], pca1$li[c(7:14, 29:39, 49:58, 272:289, 292:293),2], col = "lightcoral", pch = 6) # 49 fish (mid period/south)
points(pca1$li[c(59:83, 95:130, 143:154, 178:204),1], pca1$li[c(59:83, 95:130, 143:154, 178:204),2], col = "dodgerblue4", pch = 0) # 100 fish (late period/north)
points(pca1$li[c(84:94, 131:142, 155:177),1], pca1$li[c(84:94, 131:142, 155:177),2], col = "red", pch = 0) # 46 fish (late period/south)
points(pca1$li[c(224:225, 240:244, 264),1], pca1$li[c(224:225, 240:244, 264),2], col = "blue") # 8 fish (early period/north)
points(pca1$li[c(205:223, 226:239, 245:263),1], pca1$li[c(205:223, 226:239, 245:263),2], col = "hotpink") # 52 (early period/south)

legend("bottomleft",
       legend=c("North (1989-1993) n=8", "South (1989-1993) n=52", "North (1998-2002) n=38", "South (1998-2002) n=49", "North (2008-2012) n=100", "South (2008-2010) n=46" ),
       pch=c(1, 1, 6, 6, 0, 0),
       col=c("blue", "hotpink","steelblue3", "lightcoral", "dodgerblue4", "red"))
       # inset = c(0,-0.5),
       # border = NULL)

dev.off()

#### What about a PCA using only spatially divergent alleles?
# Read in dataset containing outlier loci
gen.larvae.outs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/masterPADElarvae.txt', header = TRUE, sep = "\t")

# Keep only larvae with genetic assignments
gen.larvae.outs2 <- gen.larvae.outs[which(is.na(gen.larvae.outs$assignment) == FALSE), ]

library(adegenet)
library(wesanderson)

gen.larvae.outs3 <- as.genind(gen.larvae.outs2[,16:35])

# PCA
sum(is.na(gen.larvae.outs3$tab)) #46
X <- scaleGen(gen.larvae.outs3, NA.method = "mean")
dim(X)
class (X)

col.palette <- wes_palette("Darjeeling1", 5, type = "discrete")
palette(col.palette)
loc.cols <- col.palette[as.numeric(factor(gen.larvae.outs2$Place))]

# make PCA
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig,main="PCA eigenvalues", col=heat.colors(50))

pca1

# Plotting PC1 and PC2
eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent [1:3]

col <- col.palette
s.class(pca1$li, factor(gen.larvae.outs2$Place), xax=1,yax=2, col = transp(col,0.9), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, clabel = 0, xlim = c(-11,5), ylim = c(-9,6))
axis(1, at=seq(-8,4, by=1), labels=seq(-8,4, by= 1), line = -2)
axis(2, at=seq(-6,6, by = 1), labels=seq(-6,6, by= 1), line = 1, las = 2)
mtext("PC1 (12.7%)", side = 1, line = 0.5)
mtext("PC2 (12.3%)", side = 2, line = 3)
title("PCA of larval summer flounder \nusing adult outliers")

legend(-9, 6,
       legend=levels(factor(gen.larvae.outs2$Place)),
       pch=19,
       col=col,
       bty = "n")
