#### Read in data ####
# Read in data for all available larvae
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/")
data <- read.table('masterPADElarvae.txt', header = TRUE)

# Read in otolith microchemistry data
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Otolith microchemistry")
otoliths <- read.table('otolith_data.txt', header = TRUE, sep = "\t")

# Read in available genetic data
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files")

library(ade4)
library(adegenet)
library(pegas)
library(hierfstat)
library(mmod)
library(poppr)
library(diveRsity)
library(tidyr)

full <- read.structure("structure_input_Mar14_2017_528_1904.str",
                       n.ind = 528, n.loc = 1904, col.lab = 1, col.pop = 2, row.marknames = 1,
                       onerowperind = FALSE)

# Which fish were successfully sequenced and have otolith data?
# Merge datasets, but first remove adults from genetic dataset and fix names
genetic <- as.data.frame(full@tab)[-c(205:439),] # remove adults
gen.fishnames.split <- as.data.frame(do.call(rbind, strsplit(as.character(rownames(genetic)), 'L', fixed = TRUE))) # fix names
gen.fishnames.split2 <- separate(gen.fishnames.split, V1, c("name1", "name2", "name3", "name4"),sep = c(4,5,7))
gen.newnames <- paste(gen.fishnames.split2$name1, gen.fishnames.split2$name3, gen.fishnames.split2$name2, gen.fishnames.split2$name4, sep = '')

gen.data <- data[data$PinskyID %in% gen.newnames,] # subset full dataset to only those fish in genetic dataset

oto.gen.data <- gen.data[gen.data$PicID %in% otoliths$Fish.ID,]# subset genetic data set to those in otolith data set
oto.gen.merge <- merge(gen.data, otoliths, by.x = 'PicID', by.y = 'Fish.ID', all = FALSE) # merged otolith and genetic data set

# data.sub <- data[data$PinskyID %in% names$V1,] # genetic data
# data.sub.dna.missing <- data[!data$PinskyID %in% names$V1,] # missing genetic data 230 fish
# data.sub.oto.missing <- data.sub[!data.sub$PicID %in% otoliths$Fish.ID,] # genetic but no otolith 142 fish
# oto.dna.missing <- otoliths[!otoliths$Fish.ID %in% data.sub$PicID,] # otolith but no genetics 42 fish
# oto.dna <- otoliths[otoliths$Fish.ID %in% data.sub$PicID,] # otolith & genetic 151 fish

table(oto.dna$Location)
table(oto.dna$Period)

# Load libraries
library(vegan)
library(factoextra)
library(MASS)
library(car)
library(plyr)
library(dendextend)
library(cluster)
library(purrr)

# Have a look at the genetic data
gen.counts <- ddply(gen.data, .(gen.data$Place, gen.data$Year), nrow)

# Have a look at the otolith data
table(otoliths$Year)
table(otoliths$Location)
table(otoliths$Period)

oto.counts <- ddply(otoliths, .(otoliths$Location, otoliths$Period), nrow)

# Have a look at the genetic & otolith data
oto.gen.counts <- ddply(oto.gen.merge, .(oto.gen.merge$Period, oto.gen.merge$Location), nrow)

####################################################################################
# Microchemistry by location
otoliths.sub <- otoliths[,c(1,4,5,12:21)]
rownames(otoliths.sub) <- otoliths.sub[,1]
otoliths.sub.log <- cbind(otoliths.sub[,c("Fish.ID", "Location", "Period", "Sr")], log10(otoliths.sub[,c('Mg', 'Mn', 'Fe', 'Cu', 'Cd', 'Ba', 'Sn', 'Pb', 'U')] +1) ) # log transform
# chem1 <- aggregate(otoliths.sub[,-1], by = list(otoliths.sub$Location), FUN = mean)
# chem1$Location <- c(4,2,1,3)
# boxplot(chem1$Mg ~ chem1$Location, ylab = "Mg", xaxt = "n")
# axis(1, at = 1:4, labels = c('RUMFS', 'Roosevelt Inlet', 'York River', 'Beaufort'))

locations <- factor(otoliths.sub$Location, c("RUMFS", "Roosevelt", "York", "NC"))

boxplot(otoliths.sub$Mg ~ locations, ylab = "Mg")
boxplot(otoliths.sub$Mn ~ locations, ylab = "Mn")
boxplot(otoliths.sub$Fe ~ locations, ylab = "Fe")
boxplot(otoliths.sub$Cu ~ locations, ylab = "Cu")
boxplot(otoliths.sub$Sr ~ locations, ylab = "Sr")
boxplot(otoliths.sub$Cd ~ locations, ylab = "Cd")
boxplot(otoliths.sub$Ba ~ locations, ylab = "Ba")
boxplot(otoliths.sub$Sn ~ locations, ylab = "Sn")
boxplot(otoliths.sub$Pb ~ locations, ylab = "Pb")
boxplot(otoliths.sub$U ~ locations, ylab = "U")

# Microchemistry by time period
period <- factor(otoliths.sub$Period, c("Early", "Mid", "Late"))

boxplot(otoliths.sub$Mg ~ period, ylab = "Mg")
boxplot(otoliths.sub$Mn ~ period, ylab = "Mn")
boxplot(otoliths.sub$Fe ~ period, ylab = "Fe")
boxplot(otoliths.sub$Cu ~ period, ylab = "Cu")
boxplot(otoliths.sub$Sr ~ period, ylab = "Sr")
boxplot(otoliths.sub$Cd ~ period, ylab = "Cd")
boxplot(otoliths.sub$Ba ~ period, ylab = "Ba")
boxplot(otoliths.sub$Sn ~ period, ylab = "Sn")
boxplot(otoliths.sub$Pb ~ period, ylab = "Pb")
boxplot(otoliths.sub$U ~ period, ylab = "U")

#### Multidimensional scaling ####
oto.chem <- otoliths[-28,c("Fish.ID", "Mg", "Mn", "Fe", "Cu", "Sr", "Cd", "Ba", "Sn", "Pb", "U")] # without PADE12_014 (194)? Scaling might take care of this. Without NCPD 074 (28)
rownames(oto.chem) <- oto.chem[, "Fish.ID"] # Make fish IDs as rownames

# Scatter plot matrix
loc.cols <- col.palette[as.numeric(otoliths$Location)]
pairs(oto.chem[,-1], col = loc.cols, lower.panel = NULL)
par(xpd = TRUE)
legend(x = 0.05, y = 0.35, cex = 2,
       legend = as.character(levels(otoliths$Location)),
       fill = col.palette)
par(xpd = NA)

oto.chem <- scale(oto.chem[, -1]) # scaling probably a good idea since there is a big range of values for each elemental ratio

# Euclidian distance between rows
oto.dist <- dist(oto.chem)

# Fit model
oto.fit <- cmdscale(oto.dist, eig = TRUE, k = 3, add = FALSE) # How many dimensions?

plot(oto.fit$eig[1:10]) # scree plot

# Plot
x <- oto.fit$points[,1]
y <- oto.fit$points[,2]
z <- oto.fit$points[,3]
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Otolith microchemistry")
text(x, y, labels = row.names(oto.chem), pos = 3)
plot(x, z, xlab = "Coordinate 1", ylab = "Coordinate 3", main = "Otolith microchemistry")
text(x, z, labels = row.names(oto.chem), pos = 2)

# Plot locations in color
library(wesanderson)

par(mfrow = c(1,2))
col.palette <- wes_palette("FantasticFox1", 5, type = "discrete")[-1]
palette(col.palette)
plot(x, y, xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry", col = otoliths$Location[-28], pch = 19) # without NCPD 074 (28)
legend("bottomleft",
       legend = levels(otoliths$Location),
       pch=19,
       col = col.palette)

plot(x, z, xlab = "MDS1", ylab = "MDS3", main = "Otolith microchemistry", col = otoliths$Location[-28], pch = 19) # without NCPD 074 (28)
legend("bottomleft",
       legend = levels(otoliths$Location),
       pch=19,
       col = col.palette)

# Plot points by time period
plot(x, y, xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry", col = otoliths$Period[-28], pch = 19) # without NCPD 074 (28)
legend("bottomleft",
       legend = levels(otoliths$Period),
       pch=19,
       col = col.palette)

plot(x, z, xlab = "MDS1", ylab = "MDS3", main = "Otolith microchemistry", col = otoliths$Period[-28], pch = 19) # without NCPD 074 (28)
legend("topleft",
       legend = levels(otoliths$Period),
       pch=19,
       col = col.palette)

# MDMS by ingress site, colored by stage/size
# Read in dataset containing outlier loci
gen.larvae.outs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/masterPADElarvae.txt', header = TRUE, sep = "\t")

# Merge otolith data with predicted sites with genetic data containing outliers
oto.gen.merge2 <- merge(gen.larvae.outs, otoliths, by.x = 'PicID', by.y = 'Fish.ID', all = FALSE) # merged otolith and genetic data set

NC <- oto.gen.merge2[which(oto.gen.merge2$Location == 'NC'), c("PinskyID", "Developmental.Stage", "Average.SL..mm.", "Mg", "Mn", "Fe", "Cu", "Sr", "Cd", "Ba", "Sn", "Pb", "U")]
York <- oto.gen.merge2[which(oto.gen.merge2$Location == 'York'), c("PinskyID", "Developmental.Stage", "Average.SL..mm.", "Mg", "Mn", "Fe", "Cu", "Sr", "Cd", "Ba", "Sn", "Pb", "U")]
Roosevelt <- oto.gen.merge2[which(oto.gen.merge2$Location == 'Roosevelt'), c("PinskyID", "Developmental.Stage", "Average.SL..mm.", "Mg", "Mn", "Fe", "Cu", "Sr", "Cd", "Ba", "Sn", "Pb", "U")]
RUMFS <- oto.gen.merge2[which(oto.gen.merge2$Location == 'RUMFS'), c("PinskyID", "Developmental.Stage", "Average.SL..mm.", "Mg", "Mn", "Fe", "Cu", "Sr", "Cd", "Ba", "Sn", "Pb", "U")]

NC.chem <- scale(NC[, -c(1:3)]) # scaling probably a good idea since there is a big range of values for each elemental ratio
York.chem <- scale(York[, -c(1:3)]) 
Roosevelt.chem <- scale(Roosevelt[, -c(1:3)]) 
RUMFS.chem <- scale(RUMFS[, -c(1:3)]) 

# Euclidian distance between rows
NC.dist <- dist(NC.chem)
York.dist <- dist(York.chem)
Roosevelt.dist <- dist(Roosevelt.chem)
RUMFS.dist <- dist(RUMFS.chem)

# Fit model
NC.fit <- cmdscale(NC.dist, eig = TRUE, k = 3, add = FALSE) 
York.fit <- cmdscale(York.dist, eig = TRUE, k = 3, add = FALSE) 
Roosevelt.fit <- cmdscale(Roosevelt.dist, eig = TRUE, k = 3, add = FALSE) 
RUMFS.fit <- cmdscale(RUMFS.dist, eig = TRUE, k = 3, add = FALSE) 
plot(RUMFS.fit$eig[1:10])

col.palette <- rainbow(6)
palette(col.palette)

# Plot NC
plot(NC.fit$points[,1], NC.fit$points[,2], xlab = "MDS1", ylab = "MDS2", main = "NC otolith microchemistry", col = NC$Developmental.Stage, pch = 19) 
legend("bottomright",
       legend = levels(NC$Developmental.Stage)[-7],
       pch=19,
       col = col.palette)

# Plot York
plot(York.fit$points[,1], York.fit$points[,2], xlab = "MDS1", ylab = "MDS2", main = "York otolith microchemistry", col = York$Developmental.Stage, pch = 19) 
legend("bottomright",
       legend = levels(York$Developmental.Stage)[-7],
       pch=19,
       col = col.palette)

# Plot Roosevelt
plot(Roosevelt.fit$points[,1], Roosevelt.fit$points[,2], xlab = "MDS1", ylab = "MDS2", main = "Roosevelt otolith microchemistry", col = Roosevelt$Developmental.Stage, pch = 19) 
legend("topright",
       legend = levels(Roosevelt$Developmental.Stage)[-7],
       pch=19,
       col = col.palette)

# Plot RUMFS
plot(RUMFS.fit$points[,1], RUMFS.fit$points[,2], xlab = "MDS1", ylab = "MDS2", main = "RUMFS otolith microchemistry", col = RUMFS$Developmental.Stage, pch = 19) 
legend("bottomright",
       legend = levels(RUMFS$Developmental.Stage)[-7],
       pch=19,
       col = col.palette)

# Plot by size
# Color vecotr
zCol <- function(nCols, Z){
  cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(nCols)
  colVec_ind <- cut(Z, breaks=nCols)
  colVec <- cols[colVec_ind]
}

# North Carolina
Z <- NC$Average.SL..mm.
plot(NC.fit$points[,1], NC.fit$points[,2], xlab = "MDS1", ylab = "MDS2", main = "NC otolith microchemistry by size", col = zCol(6,Z), pch = 19) 

my.colors = colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))
z=matrix(1:6,nrow=1)
x=0.25
y=seq(11.402, 15.040,len=7) # range of your data
image(x,y,z,col=my.colors(6),axes=FALSE,xlab='', ylab='',main='')
mtext("Standard length", side=4, line=2.5)
axis(4)

# York
Z <- York$Average.SL..mm.
plot(York.fit$points[,1], York.fit$points[,2], xlab = "MDS1", ylab = "MDS2", main = "York otolith microchemistry by size", col = zCol(6,Z), pch = 19) 

my.colors = colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))
z=matrix(1:6,nrow=1)
x=0.25
y=seq(9.4430, 16.2505,len=7) # range of your data
image(x,y,z,col=my.colors(6),axes=FALSE,xlab='', ylab='',main='')
mtext("Standard length", side=4, line=2.5)
axis(4)

# Roosevelt
Z <- Roosevelt$Average.SL..mm.
plot(Roosevelt.fit$points[,1], Roosevelt.fit$points[,2], xlab = "MDS1", ylab = "MDS2", main = "Roosevelt otolith microchemistry by size", col = zCol(6,Z), pch = 19) 

my.colors = colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))
z=matrix(1:6,nrow=1)
x=0.25
y=seq(10.242, 15.281,len=7) # range of your data
image(x,y,z,col=my.colors(6),axes=FALSE,xlab='', ylab='',main='')
mtext("Standard length", side=4, line=2.5)
axis(4)

# RUMFS
Z <- RUMFS$Average.SL..mm.
plot(RUMFS.fit$points[,1], RUMFS.fit$points[,2], xlab = "MDS1", ylab = "MDS2", main = "RUMFS otolith microchemistry by size", col = zCol(6,Z), pch = 19) 

my.colors = colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))
z=matrix(1:6,nrow=1)
x=0.25
y=seq(10.242, 15.281,len=7) # range of your data
image(x,y,z,col=my.colors(6),axes=FALSE,xlab='', ylab='',main='')
mtext("Standard length", side=4, line=2.5)
axis(4)


#### Clustering ####
# Hierarchical clustering
# fit <- hclust(oto.dist, method = "complete") # default
# fit <- hclust((oto.dist)^2, method = "ward.D")
fit <- hclust(oto.dist, method = "ward.D2") # seems this is the correct method to use
plot(fit, cex = 0.5)
plot(fit, cex = 0.5, labels = otoliths$Location[-28])
rect.hclust(fit, k=4)

# Nicer plot
dend <- as.dendrogram(fit)
dend <- rotate(dend, 1:197)

# Color branches based on clusters
dend <- color_branches(dend, k = 2)
labels_colors(dend) <- col.palette[sort_levels_values(as.numeric(oto.chem)[order.dendrogram(dend)])]


plot(dend)
circlize_dendrogram(dend)


# Assess different methods
m <- c("ward", "single", "complete", "average")
names(m) <- c("ward", "single", "complete", "average")

ac <- function(x) {
  agnes(oto.chem, method = x)$ac
}

map_dbl(m, ac)

hc.agnes <- agnes(oto.chem, method = "ward")
pltree(hc.agnes, cex = 0.6, hang = -1, main = "Dendrogram of agnes")








# Cut dendrogram based on number of sites
nclust <- 1
clusternum <- cutree(fit, k=nclust)
# clusternum <- cutree(as.hclust(hc.agnes), k=nclust)
table(clusternum)
class.table.hier <- table(otoliths$Location[-28], clusternum)
mosaicplot(class.table.hier, color = col.palette)

fviz_cluster(list(data=oto.chem, cluster = clusternum))


# hierarchical with bootstrapped p-values
library(pvclust)
fit2 <- pvclust(t(oto.chem), method.hclust = "ward.D2", method.dist = "euclidean")
plot(fit2)
pvrect(fit2, alpha = 0.95)


# k-means clustering
library(cluster)
kmean.cls <- kmeans(oto.chem, centers = 2)
class.table.km <- table(otoliths$Location, kmean.cls$cluster)
mosaicplot(class.table.km, color = col.palette)

wss <- (nrow(oto.chem)-1)*sum(apply(oto.chem,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(oto.chem, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")


#### Determining optimal number of clusters ####
fviz_nbclust(oto.chem, FUN = hcut, method = "wss")
fviz_nbclust(oto.chem, FUN = hcut, method = "silhouette")
gap_stat <- clusGap(oto.chem, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

################################################################################
#### Discriminant Function Analysis ####
library(MASS)
library(car)

# transform to normal if necessary and standardize
otoliths.sub <- otoliths[,c(1,4,5,12:21)]
rownames(otoliths.sub) <- otoliths.sub[,1]
otoliths.sub.log <- cbind(otoliths.sub[,c("Fish.ID", "Location", "Period", "Sr")], log10(otoliths.sub[,c('Mg', 'Mn', 'Fe', 'Cu', 'Cd', 'Ba', 'Sn', 'Pb', 'U')]+0.00000001) ) # log transform
# otoliths.sub.log <- cbind(otoliths.sub[,c("Fish.ID", "Location", "Period", "Sr")], log10(otoliths.sub[,c('Mg', 'Mn', 'Ba', 'Sn')]), log10(otoliths.sub[,c('Fe', 'Cu', 'Cd', 'Pb', 'U')]+1) ) # log transform, adding 1 to elements with zero values
otoliths.sub.log.trans <- as.data.frame(scale(otoliths.sub.log[4:13]))
pairs(otoliths.sub.log.trans)
scatterplotMatrix(otoliths.sub.log.trans)
Location <- otoliths.sub.log$Location
otoliths.sub.log.trans2 <- cbind(Location, otoliths.sub.log.trans)

dfa1 <- lda(Location ~ Mg + Mn + Fe + Cu + Cd + Sr + Ba + Sn + Pb + U, data = otoliths.sub.log.trans2, na.action = "na.omit", CV = TRUE) #Cd, Pb and U apprently too similar between groups when data not transformed and standardized; but not so when they are

# Assess accuracy of the prediction
# percent correct for each category of Location
ct1 <- table(otoliths.sub.log.trans2$Location, dfa1$class)
diag(prop.table(ct1,1))

# total percent correct
sum(diag(prop.table(ct1)))

#### Add predicted sites to the otolith data ####
otoliths$predicted <- dfa1$class

# Read in dataset containing outlier loci
gen.larvae.outs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/masterPADElarvae.txt', header = TRUE, sep = "\t")

# Keep only larvae with genetic assignments
gen.larvae.outs2 <- gen.larvae.outs[which(is.na(gen.larvae.outs$assignment) == FALSE), ]

# Merge otolith data with predicted sites with genetic data containing outliers
oto.gen.merge2 <- merge(gen.larvae.outs2, otoliths[,-11], by.x = 'PicID', by.y = 'Fish.ID', all = FALSE) # merged otolith and genetic data set; remove column of NAs in otolith data

# Make some plots
plot(oto.gen.merge2$assignment ~ oto.gen.merge2$predicted, xlab = "Predicted population based on otolith microchemistry", ylab = "Population assignment based on genetic likelihood")


#### Validating regional adult allele frequencies ####
# Read in allele counts of 15 spatial outliers in adults
pop.allele.counts <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.counts15.txt')
pop.allele.counts2 <- pop.allele.counts[colnames(pop.allele.counts) %in% colnames(oto.gen.merge2[,16:35])] # subset adult counts to only 10 loci occuring in larvae: 232 x 20

north.counts <- pop.allele.counts2[which(pop.allele.counts$regions == "north"),]
south.counts <- pop.allele.counts2[which(pop.allele.counts$regions == "south"),]

north.allele.freqs <- colSums(north.counts, na.rm = TRUE)/(2*colSums(!is.na(north.counts)))
south.allele.freqs <- colSums(south.counts, na.rm = TRUE)/(2*colSums(!is.na(south.counts)))

regional.allele.freqs <- rbind(north.allele.freqs, south.allele.freqs)

# Read in allele frequencies of 10 adult spatial outliers that exist in the larvae. These will be used to calculate population likelihoods of the otolith populations.
pop.allele.freqs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/regional.adult.outs.freqs10.txt')

# Do these calculated frequencies match?
colnames(regional.allele.freqs) == colnames(pop.allele.freqs) # alleles for one locus are flipped, numbers seem same, but not TRUE for some reason....?
as.data.frame(regional.allele.freqs) == pop.allele.freqs

# Pull out odd indicies
pop.allele.freqs
odds <- seq(1,20,2) # odd indicies to keep
pop.allele.freqs.odds <- pop.allele.freqs[,odds]




# Calculate population allele frequencies based on predicted otolith populations
# pop.allele.counts <- aggregate(oto.gen.merge2[,c(16:35)], by = list(oto.gen.merge2$predicted), FUN = sum, na.rm = TRUE)

# pull out even and odd indicies of allele counts
# odd.indices <- seq(2,20,2)
# even.indices <- seq(3,21,2)
# 
# odds <- pop.allele.counts[,odd.indices]
# evens <- pop.allele.counts[,even.indices]
# 
# total.allele.counts <- odds + evens
# odds.pop.allele.frequencies <- odds/total.allele.counts
# evens.pop.allele.frequencies <- evens/total.allele.counts # sums to 1
# 
# north.odds.pop.allele.frequencies <- colSums(odds[c(2,3),])/colSums(total.allele.counts[c(2,3),])
# south.odds.pop.allele.frequencies <- colSums(odds[c(1,4),])/colSums(total.allele.counts[c(1,4),])
# regional.freqs.odds <- rbind (north.odds.pop.allele.frequencies, south.odds.pop.allele.frequencies)

#### Population assignment using allele frequencies of otolith populations ####
indiv.allele.counts <- oto.gen.merge2[,c(16:35)] # Just genetic data

indiv.allele.counts.odds <- indiv.allele.counts[,odds]
indiv.allele.counts.odds[is.na(indiv.allele.counts.odds)] <- 9 # replace NA's with 9's to make the ifelse statements easier

colnames(indiv.allele.counts.odds) == colnames(pop.allele.freqs.odds) # column names match? yay!

# For loop to loop through each locus & multiply by the otolith allele frequency, and then do this for all 151 larvae. NA's/9's get coded as 1's so they don't make a difference when each row's product is taken
north.likelihoods <- data.frame()
south.likelihoods <- data.frame()

for (j in 1:length(rownames(indiv.allele.counts.odds))){
  
  for (i in 1:length(colnames(indiv.allele.counts.odds))){
    if(indiv.allele.counts.odds[j,i] == 2) {
      north.likelihoods[j,i] <- pop.allele.freqs.odds[1,i]^2
    } else if (indiv.allele.counts.odds[j,i] == 1) {
      north.likelihoods[j,i] <- 2*(pop.allele.freqs.odds[1,i] * (1-pop.allele.freqs.odds[1,i]))
    } else if (indiv.allele.counts.odds[j,i] == 0) {
      north.likelihoods[j,i] <- ( 1-pop.allele.freqs.odds[1,i])^2 
    } else {
      north.likelihoods[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(indiv.allele.counts.odds))){
    if(indiv.allele.counts.odds[j,i] == 2){
      south.likelihoods[j,i] <- pop.allele.freqs.odds[2,i]^2
    } else if (indiv.allele.counts.odds[j,i] == 1) {
      south.likelihoods[j,i] <- 2*(pop.allele.freqs.odds[2,i] * (1-pop.allele.freqs.odds[2,i]))
    } else  if (indiv.allele.counts.odds[j,i] == 0) {
      south.likelihoods[j,i] <- (1-pop.allele.freqs.odds[2,i])^2
    } else {
      south.likelihoods[j,i] <- 1
    }
  }
}

# Multiply everything together
north.vector <- vector()
south.vector <- vector()

for (k in 1:length(north.likelihoods[,1])){
  north.vector[k] <- prod(north.likelihoods[k,])
}

for (l in 1:length(south.likelihoods[,1])){
  south.vector[l] <- prod(south.likelihoods[l,])
}

# Check a few fish by hand
(0.8092105^2)*(2*0.8421053*(1-0.8421053))*(0.9078947^2)*(2*0.4934211*(1-0.4934211))*(2*0.7697368*(1-0.7697368))*(0.9276316^2)*((1-0.6447368)^2)*(0.9539474^2)*(0.9276316^2)*(0.9637681^2) #north fish1
(0.7933333^2)*(2*0.8000000*(1-0.8000000))*(0.9133333^2)*(2*0.5533333*(1-0.5533333))*(2*0.7666667*(1-0.7666667))*(0.9533333^2)*((1-0.5866667)^2)*(0.9533333^2)*(0.9466667^2)*(0.9657534^2) #south

(0.8092105^2)*(0.8421053^2)*(0.9078947^2)*(0.4934211^2)*(0.7697368^2)*(0.9276316^2)*((1-0.6447368)^2)*(0.9539474^2)*(0.9276316^2)*(0.9637681^2) # fish 20, north
(0.7933333^2)*(0.8000000^2)*(0.9133333^2)*(0.5533333^2)*(0.7666667^2)*(0.9533333^2)*((1-0.5866667)^2)*(0.9533333^2)*(0.9466667^2)*(0.9657534^2) #south

(2*0.8092105*(1-0.8092105))*(0.8421053^2)*(0.9078947^2)*(0.4934211^2)*(0.7697368^2)*(0.9276316^2)*(2*0.6447368*(1-0.6447368))*(0.9539474^2)*(2*0.9276316*(1-0.9276316))*(0.9637681^2) # fish 100, north
(2*0.7933333*(1-0.7933333))*(0.8000000^2)*(0.9133333^2)*(0.5533333^2)*(0.7666667^2)*(0.9533333^2)*(2*0.5866667*(1-0.5866667))*(0.9533333^2)*(2*0.9466667*(1-0.9466667))*(0.9657534^2) # south

# Determine if each fish has a higher likelihood of coming from the north or the south
# ratio <- north.vector/south.vector
# 
# assignment <- vector()
# for (m in 1:length(ratio)){
#   if (ratio[m] > 1){
#     assignment[m] <- "north"
#   } else {
#     assignment[m] <- "south"
#   }
# }
# 
# table(assignment) # removed fish with no genetic data

# Plot log likelihood for each population using all fish
# north.log <- log10(north.vector[-c(which(north.vector == 1))]) # Remove individuals who are missing genetic data
# south.log <- log10(south.vector[-c(which(south.vector == 1))])

region1 <- gsub("NC", "South", oto.gen.merge2$predicted) # create otolith population with region (north or south)
region2 <- gsub("York", "South", region1)
region3 <- gsub("RUMFS", "North", region2)
region4 <- gsub("Roosevelt", "North", region3)
oto.gen.merge2$predicted.region <- as.factor(region4)
# geo.color <- as.numeric(oto.gen.merge2$predicted.region)

# plot(south.log ~ north.log, xlab = 'log likelihood (north)', ylab = 'log likelihood (south)', col=ifelse(geo.color == 1, 'blue', 'tomato')) # blue is north, red is south
# abline(a = 0,b=1)
# legend("topleft",
       # legend=c("Northern otolith populations", "Southern otolith populations"),
       # pch=c(1, 1),
       # col=c("blue", "tomato"))

# plot(as.factor(assignment[-c(which(south.vector == 1))]) ~ oto.gen.merge2$predicted[-c(which(south.vector == 1))], xlab = "Predicted population based on otolith microchemistry", ylab = "Population assignment based on genetic likelihood")

#### Group likelihoods by otolith population & take the product ####
north.pop <- aggregate(north.vector, by = list(oto.gen.merge2$predicted.region), FUN = prod) # northern likelihoods for north and south otolith pops
south.pop <- aggregate(south.vector, by = list(oto.gen.merge2$predicted.region), FUN = prod) # southern likelihoods for north and south otolith pops

pops.likelihood <- as.data.frame(cbind(log10(as.numeric(north.pop[,2])), log10(as.numeric(south.pop[,2]))))
rownames(pops.likelihood) <- c("n.oto", "s.oto")
colnames(pops.likelihood) <- c("n.likes", "s.likes")
pops.likelihood$oto.pop <- c('north', 'south')

# Plot population likelihoods of otolith populations
plot(pops.likelihood[,'s.likes'] ~ pops.likelihood[,'n.likes'], ylab = 'Southern likelihood', xlab = 'Northern likelihood', col = "blue", xlim = c(-230,-200), ylim = c(-230,-200))
points(pops.likelihood[2,'n.likes'], pops.likelihood[2,'s.likes'], col = 'tomato')
abline(a=0,b=1)
legend("bottomright",
       legend=c("Northern otolith pop", "Southern otolith pop"),
       pch=c(1, 1),
       col=c("blue", "tomato"))

#########
library(tidyverse)
library(caret)
library(MASS)
library(ggplot2)

# Split data into training (80%) and test sets (20%)
training.samples <- otoliths$Location %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data <- otoliths[training.samples, c(4,12:21)]
test.data <- otoliths[-training.samples, c(4,12:21)]

# Normalize the data
# Estimate preprocessing parameters
preproc.param <- train.data %>% 
  preProcess(method = c("center", "scale"))
# Transform data using estimated parameters
train.transformed <- preproc.param %>% predict(train.data)
test.transformed <- preproc.param %>% predict(test.data)

# Fit the model
model <- lda(Location~., data = train.transformed)
plot(model)
# Make predictions
predictions <- model %>% predict(test.transformed)
# Model accuracy
mean(predictions$class==test.transformed$Location)

lda.data <- cbind(train.transformed, predict(model)$x)
ggplot(lda.data, aes(LD1, LD2)) +
  geom_point(aes(color = Location))


#########################################################
#### MDMS using genetic assignments ####
# Euclidian distance between rows
oto.dist2 <- dist(scale(oto.gen.merge2[,50:59])) # 151 fish with genetic and otolith data

# Fit model
oto.fit2 <- cmdscale(oto.dist2, eig = TRUE, k = 3, add = FALSE) # How many dimensions?

plot(oto.fit2$eig[1:10]) # scree plot

# Plot locations in color
library(wesanderson)

par(mfrow = c(1,2))
col.palette <- wes_palette("FantasticFox1", 5, type = "discrete")[-c(1:2)] # 3 colors
palette(col.palette)
plot(oto.fit2$points[-which(is.na(oto.gen.merge2$assignment)),1], oto.fit2$points[-which(is.na(oto.gen.merge2$assignment)),2], xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry by \ngenetic assignment", col = oto.gen.merge2$assignment[-which(is.na(oto.gen.merge2$assignment))], pch = 19)
legend("topleft",
       legend = levels(oto.gen.merge2$assignment),
       pch=19,
       col = col.palette)

plot(oto.fit2$points[-which(is.na(oto.gen.merge2$assignment)),1], oto.fit2$points[-which(is.na(oto.gen.merge2$assignment)),3], xlab = "MDS1", ylab = "MDS3", main = "Otolith microchemistry by \ngenetic assignment", col = oto.gen.merge2$assignment[-which(is.na(oto.gen.merge2$assignment))], pch = 19)
legend("topleft",
       legend = levels(oto.gen.merge2$assignment),
       pch=19,
       col = col.palette)

