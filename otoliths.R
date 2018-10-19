#### Read in data ####
# Read in data for all available larvae
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/")
data <- read.table('masterPADElarvae.txt', header = TRUE)

# Read in otolith microchemistry data
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/")
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

boxplot(otoliths.sub$Mg ~ locations, ylab = "Mg:Ca")
boxplot(otoliths.sub$Mn ~ locations, ylab = "Mn:Ca")
boxplot(otoliths.sub$Fe ~ locations, ylab = "Fe:Ca")
boxplot(otoliths.sub$Cu ~ locations, ylab = "Cu:Ca")
boxplot(otoliths.sub$Sr ~ locations, ylab = "Sr:Ca")
boxplot(otoliths.sub$Cd ~ locations, ylab = "Cd:Ca")
boxplot(otoliths.sub$Ba ~ locations, ylab = "Ba:Ca")
boxplot(otoliths.sub$Sn ~ locations, ylab = "Sn:Ca")
boxplot(otoliths.sub$Pb ~ locations, ylab = "Pb:Ca")
boxplot(otoliths.sub$U ~ locations, ylab = "U:Ca")

#### Plotting & ANOVAs for each element ####
library(lsmeans)
library(multcompView)
lsmeans = lsmeans::lsmeans

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/otolith_boxplots.png", width=8.5, height=11, res=300, units="in")
par(
  mar=c(4, 4, 1, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  mfrow = c(5,2)
)

boxplot(otoliths.sub.log$Mg ~ locations, ylab = "Mg:Ca", col = c("#B40F20", "#E58601", "#E2D200", "#46ACC8"), outpch = 19)
mg.anova <- aov(otoliths.sub.log$Mg ~ locations)
leastsquare <- lsmeans(mg.anova, pairwise ~ locations, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

boxplot(otoliths.sub.log$Mn ~ locations, ylab = "Mn:Ca", col = c("#B40F20", "#E58601", "#E2D200", "#46ACC8"), outpch = 19)
mn.anova <- aov(otoliths.sub$Mn ~ locations)
leastsquare <- lsmeans(mn.anova, pairwise ~ locations, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

boxplot(otoliths.sub.log$Fe ~ locations, ylab = "Fe:Ca", col = c("#B40F20", "#E58601", "#E2D200", "#46ACC8"), outpch = 19)
fe.anova <- aov(otoliths.sub$Fe ~ locations)
leastsquare <- lsmeans(fe.anova, pairwise ~ locations, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

boxplot(otoliths.sub.log$Cu ~ locations, ylab = "Cu:Ca", col = c("#B40F20", "#E58601", "#E2D200", "#46ACC8"), outpch = 19)
cu.anova <- aov(otoliths.sub$Cu ~ locations)
leastsquare <- lsmeans(cu.anova, pairwise ~ locations, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

boxplot(otoliths.sub.log$Sr ~ locations, ylab = "Sr:Ca", col = c("#B40F20", "#E58601", "#E2D200", "#46ACC8"), outpch = 19)
sr.anova <- aov(otoliths.sub$Sr ~ locations)
leastsquare <- lsmeans(sr.anova, pairwise ~ locations, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

boxplot(otoliths.sub.log$Cd ~ locations, ylab = "Cd:Ca", col = c("#B40F20", "#E58601", "#E2D200", "#46ACC8"), outpch = 19)
cd.anova <- aov(otoliths.sub$Cd ~ locations)
leastsquare <- lsmeans(cd.anova, pairwise ~ locations, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

boxplot(otoliths.sub.log$Ba ~ locations, ylab = "Ba:Ca", col = c("#B40F20", "#E58601", "#E2D200", "#46ACC8"), outpch = 19)
ba.anova <- aov(otoliths.sub$Ba ~ locations)
leastsquare <- lsmeans(ba.anova, pairwise ~ locations, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

boxplot(otoliths.sub.log$Sn ~ locations, ylab = "Sn:Ca", col = c("#B40F20", "#E58601", "#E2D200", "#46ACC8"), outpch = 19)
sn.anova <- aov(otoliths.sub$Sn ~ locations)
leastsquare <- lsmeans(sn.anova, pairwise ~ locations, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

boxplot(otoliths.sub.log$Pb ~ locations, ylab = "Pb:Ca", col = c("#B40F20", "#E58601", "#E2D200", "#46ACC8"), outpch = 19)
mtext("Ingress site",1, line = 2.5)
pb.anova <- aov(otoliths.sub$Pb ~ locations)
leastsquare <- lsmeans(pb.anova, pairwise ~ locations, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

boxplot(otoliths.sub.log$U ~ locations, ylab = "U:Ca", col = c("#B40F20", "#E58601", "#E2D200", "#46ACC8"), outpch = 19)
mtext("Ingress site",1, line = 2.5)
u.anova <- aov(otoliths.sub$U ~ locations)
leastsquare <- lsmeans(u.anova, pairwise ~ locations, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

dev.off()
    
# Microchemistry by time period
period <- factor(otoliths.sub$Period, c("Early", "Mid", "Late"))

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/otolith_overtime_boxplots.png", width=8.5, height=11, res=300, units="in")
par(
  mar=c(4, 4, 1, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  mfrow = c(5,2)
)

boxplot(otoliths.sub.log$Mg ~ period, ylab = "Mg:Ca", col = c("#E58601", "#E2D200", "#46ACC8"), outpch = 19)
mg.anova <- aov(otoliths.sub.log$Mg ~ period)
leastsquare <- lsmeans(mg.anova, pairwise ~ period, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

boxplot(otoliths.sub.log$Mn ~ period, ylab = "Mn:Ca", col = c("#E58601", "#E2D200", "#46ACC8"), outpch = 19)
mg.anova <- aov(otoliths.sub.log$Mn ~ period)
leastsquare <- lsmeans(mg.anova, pairwise ~ period, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

boxplot(otoliths.sub.log$Fe ~ period, ylab = "Fe:Ca", col = c("#E58601", "#E2D200", "#46ACC8"), outpch = 19)
mg.anova <- aov(otoliths.sub.log$Fe ~ period)
leastsquare <- lsmeans(mg.anova, pairwise ~ period, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

boxplot(otoliths.sub.log$Cu ~ period, ylab = "Cu:Ca", col = c("#E58601", "#E2D200", "#46ACC8"), outpch = 19)
mg.anova <- aov(otoliths.sub.log$Cu ~ period)
leastsquare <- lsmeans(mg.anova, pairwise ~ period, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

boxplot(otoliths.sub.log$Sr ~ period, ylab = "Sr:Ca", col = c("#E58601", "#E2D200", "#46ACC8"), outpch = 19)
mg.anova <- aov(otoliths.sub.log$Sr ~ period)
leastsquare <- lsmeans(mg.anova, pairwise ~ period, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

boxplot(otoliths.sub.log$Cd ~ period, ylab = "Cd:Ca", col = c("#E58601", "#E2D200", "#46ACC8"), outpch = 19)
mg.anova <- aov(otoliths.sub.log$Cd ~ period)
leastsquare <- lsmeans(mg.anova, pairwise ~ period, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

boxplot(otoliths.sub.log$Ba ~ period, ylab = "Ba:Ca", col = c("#E58601", "#E2D200", "#46ACC8"), outpch = 19)
mg.anova <- aov(otoliths.sub.log$Ba ~ period)
leastsquare <- lsmeans(mg.anova, pairwise ~ period, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

boxplot(otoliths.sub.log$Sn ~ period, ylab = "Sn:Ca", col = c("#E58601", "#E2D200", "#46ACC8"), outpch = 19)
mg.anova <- aov(otoliths.sub.log$Sn ~ period)
leastsquare <- lsmeans(mg.anova, pairwise ~ period, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

boxplot(otoliths.sub.log$Pb ~ period, ylab = "Pb:Ca", col = c("#E58601", "#E2D200", "#46ACC8"), outpch = 19)
mtext("Time period",1, line = 2.5)
mg.anova <- aov(otoliths.sub.log$Pb ~ period)
leastsquare <- lsmeans(mg.anova, pairwise ~ period, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

boxplot(otoliths.sub.log$U ~ period, ylab = "U:Ca", col = c("#E58601", "#E2D200", "#46ACC8"), outpch = 19)
mtext("Time period",1, line = 2.5)
mg.anova <- aov(otoliths.sub.log$U ~ period)
leastsquare <- lsmeans(mg.anova, pairwise ~ period, adjust = "tukey")
cld(leastsquare, alpha = 0.05, Letters = letters)

dev.off()

#### Multidimensional scaling ####
# oto.chem <- otoliths[-28,c("Fish.ID", "Mg", "Mn", "Fe", "Cu", "Sr", "Cd", "Ba", "Sn", "Pb", "U")] # without PADE12_014 (194)? Scaling might take care of this. Without NCPD 074 (28)
oto.chem <- otoliths[,c("Fish.ID", "Mg", "Mn", "Fe", "Sn")] # why use all the elements when only 4 are different between regions?
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

# Histograms of MDS axis - unimodal or multimodal?
hist(x, xlab = "MDS 1", main = "") # all unimodal
hist(y, xlab = "MDS 2", main = "")
hist(z, xlab = "MDS 3", main = "")

# By ingress site?
hist(x[which(otoliths$Location == "NC")], xlab = "MDS 1", main = "NC")
hist(x[which(otoliths$Location == "York")], xlab = "MDS 1", main = "York")
hist(x[which(otoliths$Location == "Roosevelt")], xlab = "MDS 1", main = "Roosevelt")
hist(x[which(otoliths$Location == "RUMFS")], xlab = "MDS 1", main = "RUMFS")

hist(y[which(otoliths$Location == "NC")], xlab = "MDS 2", main = "NC")
hist(y[which(otoliths$Location == "York")], xlab = "MDS 2", main = "York")
hist(y[which(otoliths$Location == "Roosevelt")], xlab = "MDS 2", main = "Roosevelt")
hist(y[which(otoliths$Location == "RUMFS")], xlab = "MDS 2", main = "RUMFS")

# Plot MDS
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
# plot(x, y, xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry", col = otoliths$Location[-28], pch = 19) # without NCPD 074 (28)
plot(x, y, xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry", col = otoliths$Location, pch = 19) # mds using Mg, Mn, Fe & Sn
legend("bottomleft",
       legend = levels(otoliths$Location),
       pch=19,
       col = col.palette)

# plot(x, z, xlab = "MDS1", ylab = "MDS3", main = "Otolith microchemistry", col = otoliths$Location[-28], pch = 19) # without NCPD 074 (28)
plot(x, z, xlab = "MDS1", ylab = "MDS3", main = "Otolith microchemistry", col = otoliths$Location, pch = 19) # mds using Mg, Mn, Fe & Sn
legend("topright",
       legend = levels(otoliths$Location),
       pch=19,
       col = col.palette)

# Plot points by time period
# plot(x, y, xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry", col = otoliths$Period[-28], pch = 19) # without NCPD 074 (28)
plot(x, y, xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry", col = otoliths$Period, pch = 19) # mds using Mg, Mn, Fe & Sn
legend("bottomleft",
       legend = levels(otoliths$Period),
       pch=19,
       col = col.palette)

# plot(x, z, xlab = "MDS1", ylab = "MDS3", main = "Otolith microchemistry", col = otoliths$Period[-28], pch = 19) # without NCPD 074 (28)
plot(x, z, xlab = "MDS1", ylab = "MDS3", main = "Otolith microchemistry", col = otoliths$Period, pch = 19) # mds using Mg, Mn, Fe & Sn
legend("topright",
       legend = levels(otoliths$Period),
       pch=19,
       col = col.palette)

# Plot by year
col.palette <- rainbow(14)
palette(col.palette)
# plot(x, y, xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry", col = as.factor(otoliths$Year[-28]), pch = 19, xlim = c(-4,6), ylim = c(-5,4)) # without NCPD 074 (28); Year needs to be a factor in order for colors to match properly!
plot(x, y, xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry", col = as.factor(otoliths$Year), pch = 19, xlim = c(-4,6), ylim = c(-5,4)) # # mds using Mg, Mn, Fe & Sn
legend("topright",
        legend = levels(as.factor(otoliths$Year)),
        pch=19,
        col = col.palette,
       cex = 0.8)

# cbind(x,y,otoliths$Year[-28])
cbind(x,y,otoliths$Year)
plot(x[1:6], y[1:6], xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry", col = '#FF0000FF', pch = 19, xlim = c(-4,10), ylim = c(-5,4)) #1989
points(x[124:125], y[124:125], col = "#FF6D00FF", pch = 19) #1990
points(x[c(7:15,126)], y[c(7:15,126)], col = "#FFDB00FF", pch = 19) #1991
points(x[c(16:17,127:130)], y[c(16:17,127:130)], col = "#B6FF00FF", pch = 19) #1992
points(x[c(18:26,131:137)], y[c(18:26,131:137)], col = "#49FF00FF", pch = 19) #1998
points(x[c(27,138:139)], y[c(27,138:139)], col = "#00FF24FF", pch = 19) #1999
points(x[c(28:30,140:145)], y[c(28:30,140:145)], col = "#00FF92FF", pch = 19) #2000
points(x[c(31:41,146:158)], y[c(31:41,146:158)], col = "#00FFFFFF", pch = 19) #2001
points(x[159:162], y[159:162], col = "#0092FFFF", pch = 19) #2002
points(x[c(42:47,58:59, 83:99, 163:167)], y[c(42:47,58:59, 83:99, 163:167)], col = "#0024FFFF", pch = 19) #2008
points(x[c(48:55,60:64,100:118,168:179)], y[c(48:55,60:64,100:118,168:179)], col = "#4900FFFF", pch = 19) #2009
points(x[c(56:57,65:72,119:123,180:185)], y[c(56:57,65:72,119:123,180:185)], col = "#B600FFFF", pch = 19) #2010
points(x[c(73:77,186:188)], y[c(73:77,186:188)], col = "#FF00DBFF", pch = 19) #2011
points(x[c(78:82,189:196)], y[c(78:82,189:196)], col = "#FF006DFF", pch = 19) #2012

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
# Code to calculate and plot elbow method of deciding how many K's
wss <- (nrow(oto.chem)-1)*sum(apply(oto.chem,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(oto.chem, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#### Determining optimal number of clusters ####
library(factoextra)
library(cluster)
fviz_nbclust(oto.chem, FUN = hcut, method = "wss")
fviz_nbclust(oto.chem, FUN = hcut, method = "silhouette")
gap_stat <- clusGap(oto.chem, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

library(NbClust)
nb <- NbClust(oto.chem, distance = 'euclidean', min.nc = 2, max.nc = 10, method = 'kmeans')
fviz_nbclust(nb)
dev.off()

library(mclust) 
mclust.fit <- Mclust(oto.chem)
plot(mclust.fit)
summary(mclust.fit)

# k-means clustering
col.palette <- wes_palette("Darjeeling1", 5, type = "discrete")
palette(col.palette)
kmean.cls <- kmeans(oto.chem, centers = 5, nstart = 50, iter.max = 10)
class.table.km <- table(otoliths$Location, kmean.cls$cluster)
mosaicplot(class.table.km, color = col.palette)

fviz_cluster(kmean.cls, oto.chem, 
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07", "#D8BFD8"), 
             ellipse.type = "euclid", # Concentration ellipse
             star.plot = TRUE, # Add segments from centroids to items
             repel = TRUE, # Avoid label overplotting (slow)
             ggtheme = theme_minimal(),
             geom = "point"
             )

#### Plotting ####
# Hierarchical clustering
# fit <- hclust(oto.dist, method = "complete") # default
# fit <- hclust((oto.dist)^2, method = "ward.D")
fit <- hclust(oto.dist, method = "ward.D2") # seems this is the correct method to use
plot(fit, cex = 0.5)
plot(fit, cex = 0.5, labels = otoliths$Location)
rect.hclust(fit, k=5)

# Nicer plot
library(ape)
library(dendextend)
dend <- as.dendrogram(fit)
dend <- rotate(dend, 1:197)

# Color branches based on clusters
dend <- color_branches(dend, k = 5)
labels_colors(dend) <- col.palette[sort_levels_values(as.numeric(oto.chem)[order.dendrogram(dend)])]
plot(dend)
circlize_dendrogram(dend)

# Cut dendrogram based on number of clusters from heirarchical clustering
nclust <- 5
clusternum <- cutree(fit, k=nclust)
# clusternum <- cutree(as.hclust(hc.agnes), k=nclust)
table(clusternum)
# class.table.hier <- table(otoliths$Location[-28], clusternum)
class.table.hier <- table(otoliths$Location, clusternum)
mosaicplot(class.table.hier, color = col.palette)

fviz_cluster(list(data=oto.chem, cluster = clusternum))

# hierarchical with bootstrapped p-values
# library(pvclust)
# fit2 <- pvclust(t(oto.chem), method.hclust = "ward.D2", method.dist = "euclidean")
# plot(fit2)
# pvrect(fit2, alpha = 0.95) # something not working


################################################################################
#### Discriminant Function Analysis ####
library(MASS)
library(car)

# transform to normal if necessary and standardize
# otoliths.sub <- otoliths[,c(1,4,5,12:21)]
otoliths.sub <- otoliths[,c("Fish.ID","Location","Period","Mg","Mn","Fe","Sn")]
rownames(otoliths.sub) <- otoliths.sub[,1]
# otoliths.sub.log <- cbind(otoliths.sub[,c("Fish.ID", "Location", "Period", "Sr")], log10(otoliths.sub[,c('Mg', 'Mn', 'Fe', 'Cu', 'Cd', 'Ba', 'Sn', 'Pb', 'U')]+0.00000001) ) # log transform
# otoliths.sub.log <- cbind(otoliths.sub[,c("Fish.ID", "Location", "Period", "Sr")], log10(otoliths.sub[,c('Mg', 'Mn', 'Ba', 'Sn')]), log10(otoliths.sub[,c('Fe', 'Cu', 'Cd', 'Pb', 'U')]+1) ) # log transform, adding 1 to elements with zero values
otoliths.sub.log <- cbind(otoliths.sub[,c("Fish.ID", "Location", "Period")], log10(otoliths.sub[,c('Mg', 'Mn', 'Fe', 'Sn')]+0.00000001) ) # log transform
# otoliths.sub.log.trans <- as.data.frame(scale(otoliths.sub.log[4:13]))
otoliths.sub.log.trans <- as.data.frame(scale(otoliths.sub.log[4:7]))
pairs(otoliths.sub.log.trans)
scatterplotMatrix(otoliths.sub.log.trans)
# Location <- otoliths.sub.log$Location
# otoliths.sub.log.trans2 <- cbind(Location, otoliths.sub.log.trans)

# dfa1 <- lda(Location ~ Mg + Mn + Fe + Cu + Cd + Sr + Ba + Sn + Pb + U, data = otoliths.sub.log.trans2, na.action = "na.omit", CV = TRUE) #Cd, Pb and U apprently too similar between groups when data not transformed and standardized; but not so when they are

Cluster <- as.factor(kmean.cls$cluster)
names(kmean.cls$cluster) == rownames(otoliths.sub.log.trans) # make sure individuals are in the right order
otoliths.sub.log.trans2 <- cbind(Cluster, otoliths.sub.log.trans)
dfa1 <- lda(Cluster ~ Mg + Mn + Fe + Sn, data = otoliths.sub.log.trans2, na.action = "na.omit", CV = TRUE) 

# Assess accuracy of the prediction
# percent correct for each category of Location
# ct1 <- table(otoliths.sub.log.trans2$Location, dfa1$class)
ct1 <- table(otoliths.sub.log.trans2$Cluster, dfa1$class)
diag(prop.table(ct1,1))

# total percent correct
sum(diag(prop.table(ct1)))

#### Add predicted sites to the otolith data ####
# otoliths$predicted <- dfa1$class # using DFA classes
otoliths$cluster <- kmean.cls$cluster # using heirarchical clustering

names(kmean.cls$cluster) == otoliths$Fish.ID
# rownames(otoliths)==lda.class.ordered[,1] # should be in same order
# otoliths$predicted <- lda.class.ordered[,2] # from LDA using random 40% from each ingress site as test dataset

# Read in dataset containing outlier loci
gen.larvae.outs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/masterPADElarvae.txt', header = TRUE, sep = "\t")

# Keep only larvae with genetic assignments
gen.larvae.outs2 <- gen.larvae.outs[which(is.na(gen.larvae.outs$assignment) == FALSE), ]

# Merge otolith data with predicted sites with genetic data containing outliers
oto.gen.merge2 <- merge(gen.larvae.outs2, otoliths[,-11], by.x = 'PicID', by.y = 'Fish.ID', all = FALSE) # merged otolith and genetic data set; remove column of NAs in otolith data
write.table(oto.gen.merge2, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/oto.gen.merged151.txt", col.names = TRUE, row.names = FALSE)

# Make some plots
# plot(oto.gen.merge2$assignment ~ oto.gen.merge2$predicted, xlab = "Predicted population based on otolith microchemistry", ylab = "Population assignment based on genetic likelihood")
plot(oto.gen.merge2$assignment ~ oto.gen.merge2$cluster, xlab = "Cluster based on otolith microchemistry", ylab = "Population assignment based on genetic likelihood")

#### Validating regional adult allele frequencies ####
# Read in allele counts of 15 spatial outliers in adults
pop.allele.counts <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.counts15.txt')
pop.allele.counts2 <- pop.allele.counts[colnames(pop.allele.counts) %in% colnames(oto.gen.merge2[,16:35])] # subset adult counts to only 10 loci occuring in larvae: 232 x 20

north.counts <- pop.allele.counts2[which(pop.allele.counts$regions == "north"),]
south.counts <- pop.allele.counts2[which(pop.allele.counts$regions == "south"),]

north.allele.freqs <- colSums(north.counts, na.rm = TRUE)/(2*colSums(!is.na(north.counts)))
south.allele.freqs <- colSums(south.counts, na.rm = TRUE)/(2*colSums(!is.na(south.counts)))

regional.allele.freqs <- rbind(north.allele.freqs, south.allele.freqs)

#### Read in allele frequencies of 10 adult spatial outliers that exist in the larvae. These will be used to calculate population likelihoods of the otolith populations ####
pop.allele.freqs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/regional.adult.outs.freqs10.txt')

# Do these calculated frequencies match?
colnames(regional.allele.freqs) == colnames(pop.allele.freqs) # alleles for one locus are flipped, numbers seem same, but not TRUE for some reason....?
as.data.frame(regional.allele.freqs) == pop.allele.freqs

# Pull out odd indicies
pop.allele.freqs
odds <- seq(1,20,2) # odd indicies to keep
pop.allele.freqs.odds <- pop.allele.freqs[,odds]

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
# (0.8092105^2)*(2*0.8421053*(1-0.8421053))*(0.9078947^2)*(2*0.4934211*(1-0.4934211))*(2*0.7697368*(1-0.7697368))*(0.9276316^2)*((1-0.6447368)^2)*(0.9539474^2)*(0.9276316^2)*(0.9637681^2) #north fish1
# (0.7933333^2)*(2*0.8000000*(1-0.8000000))*(0.9133333^2)*(2*0.5533333*(1-0.5533333))*(2*0.7666667*(1-0.7666667))*(0.9533333^2)*((1-0.5866667)^2)*(0.9533333^2)*(0.9466667^2)*(0.9657534^2) #south
# 
# (0.8092105^2)*(0.8421053^2)*(0.9078947^2)*(0.4934211^2)*(0.7697368^2)*(0.9276316^2)*((1-0.6447368)^2)*(0.9539474^2)*(0.9276316^2)*(0.9637681^2) # fish 20, north
# (0.7933333^2)*(0.8000000^2)*(0.9133333^2)*(0.5533333^2)*(0.7666667^2)*(0.9533333^2)*((1-0.5866667)^2)*(0.9533333^2)*(0.9466667^2)*(0.9657534^2) #south
# 
# (2*0.8092105*(1-0.8092105))*(0.8421053^2)*(0.9078947^2)*(0.4934211^2)*(0.7697368^2)*(0.9276316^2)*(2*0.6447368*(1-0.6447368))*(0.9539474^2)*(2*0.9276316*(1-0.9276316))*(0.9637681^2) # fish 100, north
# (2*0.7933333*(1-0.7933333))*(0.8000000^2)*(0.9133333^2)*(0.5533333^2)*(0.7666667^2)*(0.9533333^2)*(2*0.5866667*(1-0.5866667))*(0.9533333^2)*(2*0.9466667*(1-0.9466667))*(0.9657534^2) # south

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

region1 <- gsub("NC", "South", oto.gen.merge2$predicted) # create otolith population with region (north or south); this is likely wrong approach because based on LDA using ingress sites rather than clustering
region2 <- gsub("York", "South", region1)
region3 <- gsub("RUMFS", "North", region2)
region4 <- gsub("Roosevelt", "North", region3)
oto.gen.merge2$predicted.region <- as.factor(region4)
# geo.color <- as.numeric(oto.gen.merge2$predicted.region)

#### Group likelihoods by otolith population/cluster & take the product ####
# north.pop <- aggregate(north.vector, by = list(oto.gen.merge2$predicted.region), FUN = prod) # northern likelihoods for north and south otolith pops
# south.pop <- aggregate(south.vector, by = list(oto.gen.merge2$predicted.region), FUN = prod) # southern likelihoods for north and south otolith pops

north.pop <- aggregate(north.vector, by = list(oto.gen.merge2$cluster), FUN = prod) # northern likelihoods for 5 otolith clusters
south.pop <- aggregate(south.vector, by = list(oto.gen.merge2$cluster), FUN = prod) # southern likelihoods for 5 otolith clusters


pops.likelihood <- as.data.frame(cbind(log10(as.numeric(north.pop[,2])), log10(as.numeric(south.pop[,2]))))
colnames(pops.likelihood) <- c("n.likes", "s.likes")
rownames(pops.likelihood) <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5")

# pops.likelihood$oto.pop <- c('north', 'south')

# Plot population likelihoods of otolith populations/clusters
# plot(pops.likelihood[,'s.likes'] ~ pops.likelihood[,'n.likes'], ylab = 'Southern likelihood', xlab = 'Northern likelihood', col = "blue", xlim = c(-280,-150), ylim = c(-280,-150))
# points(pops.likelihood[2,'n.likes'], pops.likelihood[2,'s.likes'], col = 'tomato')
# abline(a=0,b=1)
# legend("bottomright",
#        legend=c("Northern otolith pop", "Southern otolith pop"),
#        pch=c(1, 1),
#        col=c("blue", "tomato"))

plot(pops.likelihood[,'s.likes'] ~ pops.likelihood[,'n.likes'], ylab = 'Southern likelihood', xlab = 'Northern likelihood', col = as.factor(rownames(pops.likelihood)), xlim = c(-135,-20), ylim = c(-135,-20), pch = 19)
abline(a=0,b=1)
legend("bottomright",
       legend=c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"),
       pch=c(19, 19, 19, 19, 19),
       col=as.factor(rownames(pops.likelihood)))

# PCA using 10 loci in larvae, colored by otolith cluster
gen.data <- as.genind(indiv.allele.counts)
pop(gen.data) <- oto.gen.merge2$cluster

sum(is.na(gen.data$tab)) #18
X <- scaleGen(gen.data, NA.method = "mean")
dim(X)
class (X)

# make PCA
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

pca1

# Plotting PC1 and PC2
eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent [1:3]

col <- col.palette
s.class(pca1$li, pop(gen.data), xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, xlim = c(-11,6), ylim = c(-7,7), clabel = 0)
axis(1, at=seq(-9,5, by=1), labels=seq(-9,5, by= 1), line = 1)
axis(2, at=seq(-6,6, by = 1), labels=seq(-6,6, by= 1), line = 0.5, las = 2)
mtext("PC1 (14.6%)", side = 1, line = 3.3)
mtext("PC2 (13.3%)", side = 2, line = 3)
title("PCA of summer flounder dataset\naxes 1-2")

legend(-8, -3,
       legend=c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"),
       pch=c(19, 19, 19, 19, 19),
       col=col,
       bty = "n",
       y.intersp = 1)

# Plotting PC1 and PC3
s.class(pca1$li, pop(gen.data), xax=1,yax=3, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, xlim = c(-11,6), ylim = c(-5,4), clabel = 0)
axis(1, at=seq(-9,5, by=1), labels=seq(-9,5, by= 1), line = 0)
axis(2, at=seq(-6,5, by = 1), labels=seq(-6,5, by= 1), line = 0.5, las = 2)
mtext("PC1 (14.6%)", side = 1, line = 3)
mtext("PC2 (12.1%)", side = 2, line = 3)
title("PCA of summer flounder dataset\naxes 1-2")

legend(-8, -3,
       legend=c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"),
       pch=c(19, 19, 19, 19, 19),
       col=col,
       bty = "n",
       y.intersp = 1)

#### Trying to calculate confidence intervals ####
# adults10_counts <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.counts10.txt', header = TRUE)
# 
# adults10_counts_north <- adults10_counts[which(adults10_counts$regions == 'north'), -c(1:2)]
# adults10_counts_south <- adults10_counts[which(adults10_counts$regions == 'south'), -c(1:2)]
# 
# colSums(adults10_counts_north, na.rm = TRUE)/(2*colSums(!is.na(adults10_counts_north)))
# colSums(adults10_counts_south, na.rm = TRUE)/(2*colSums(!is.na(adults10_counts_south)))

# Read in adult outlier allele frequencies
pop.allele.freqs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/regional.adult.outs.freqs10.txt')

# Keep only one allele per locus
odds <- seq(1,20,2) # odd indicies to keep
pop.allele.freqs.odds <- pop.allele.freqs[,odds]
pop.allele.freqs.odds.north <- pop.allele.freqs.odds[1,] # north allele frequency
pop.allele.freqs.odds.south <- pop.allele.freqs.odds[2,] # south allele frequency

# Sample one allele per locus many times
n.alleles <- sapply(pop.allele.freqs.odds.north,function(z){rbinom(2000,1,z)})

# Do this for the southern distribution too
s.alleles <- sapply(pop.allele.freqs.odds.south,function(z){rbinom(2000,1,z)})

# Then create the genotype for each individual by separating even and odd rows, then adding them together to get allele count
odds <- seq(1,2000,2) # odd indicies
evens <- seq(2,2000,2) # even indicies

n.alleles.odds <- n.alleles[odds,]
n.alleles.evens <- n.alleles[evens,]
n.alleles.sum <- n.alleles[odds,] + n.alleles[evens,]

s.alleles.odds <- s.alleles[odds,]
s.alleles.evens <- s.alleles[evens,]
s.alleles.sum <- s.alleles[odds,] + s.alleles[evens,]

# Use adult allele frequencies to calculate likelihoods
colnames(pop.allele.freqs.odds) == colnames(n.alleles.sum) # check to make sure column names are the same
colnames(pop.allele.freqs.odds) == colnames(s.alleles.sum)

# For loop to loop through each of 10 loci & multiply by the adult allele frequency, and then do this for all 1000 offspring. NA's/9's get coded as 1's so they don't make a difference when each row's product is taken
north.likelihoods.ndist <- data.frame()
south.likelihoods.ndist <- data.frame()

for (j in 1:nrow(n.alleles.sum)){
  
  for (i in 1:length(colnames(n.alleles.sum))){
    if(n.alleles.sum[j,i] == 2) {
      north.likelihoods.ndist[j,i] <- pop.allele.freqs.odds[1,i]^2
    } else if (n.alleles.sum[j,i] == 1) {
      north.likelihoods.ndist[j,i] <- 2*(pop.allele.freqs.odds[1,i] * (1-pop.allele.freqs.odds[1,i]))
    } else if (n.alleles.sum[j,i] == 0) {
      north.likelihoods.ndist[j,i] <- ( 1-pop.allele.freqs.odds[1,i])^2 
    } else {
      north.likelihoods.ndist[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(n.alleles.sum))){
    if(n.alleles.sum[j,i] == 2){
      south.likelihoods.ndist[j,i] <- pop.allele.freqs.odds[2,i]^2
    } else if (n.alleles.sum[j,i] == 1) {
      south.likelihoods.ndist[j,i] <- 2*(pop.allele.freqs.odds[2,i] * (1-pop.allele.freqs.odds[2,i]))
    } else  if (n.alleles.sum[j,i] == 0) {
      south.likelihoods.ndist[j,i] <- (1-pop.allele.freqs.odds[2,i])^2
    } else {
      south.likelihoods.ndist[j,i] <- 1
    }
  }
}

north.likelihoods.sdist <- data.frame()
south.likelihoods.sdist <- data.frame()

for (j in 1:nrow(s.alleles.sum)){
  
  for (i in 1:length(colnames(s.alleles.sum))){
    if(s.alleles.sum[j,i] == 2) {
      north.likelihoods.sdist[j,i] <- pop.allele.freqs.odds[1,i]^2
    } else if (s.alleles.sum[j,i] == 1) {
      north.likelihoods.sdist[j,i] <- 2*(pop.allele.freqs.odds[1,i] * (1-pop.allele.freqs.odds[1,i]))
    } else if (s.alleles.sum[j,i] == 0) {
      north.likelihoods.sdist[j,i] <- ( 1-pop.allele.freqs.odds[1,i])^2 
    } else {
      north.likelihoods.sdist[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(s.alleles.sum))){
    if(s.alleles.sum[j,i] == 2){
      south.likelihoods.sdist[j,i] <- pop.allele.freqs.odds[2,i]^2
    } else if (s.alleles.sum[j,i] == 1) {
      south.likelihoods.sdist[j,i] <- 2*(pop.allele.freqs.odds[2,i] * (1-pop.allele.freqs.odds[2,i]))
    } else  if (s.alleles.sum[j,i] == 0) {
      south.likelihoods.sdist[j,i] <- (1-pop.allele.freqs.odds[2,i])^2
    } else {
      south.likelihoods.sdist[j,i] <- 1
    }
  }
}

# Multiply everything together
north.vector.ndist <- vector()
south.vector.ndist <- vector()

for (k in 1:length(north.likelihoods.ndist[,1])){
  north.vector.ndist[k] <- prod(north.likelihoods.ndist[k,])
}

for (l in 1:length(south.likelihoods.ndist[,1])){
  south.vector.ndist[l] <- prod(south.likelihoods.ndist[l,])
}

north.vector.sdist <- vector()
south.vector.sdist <- vector()

for (k in 1:length(north.likelihoods.sdist[,1])){
  north.vector.sdist[k] <- prod(north.likelihoods.sdist[k,])
}

for (l in 1:length(south.likelihoods.sdist[,1])){
  south.vector.sdist[l] <- prod(south.likelihoods.sdist[l,])
}

# Create ratio & plot
hist(log10(north.vector.ndist/south.vector.ndist), xlab = "log10(north likelihood/south likelihood)", main = "", col = rgb(1,0,0,0.5), xlim = c(-4,4), ylim = c(0,400))
hist(log10(north.vector.sdist/south.vector.sdist), col = rgb(0,0,1,0.5), add = TRUE)
legend("topright", c("North", "South"), col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch = 15)

# Read in assignments
larvs.assign <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/results/larval_assignment.txt", header = TRUE)

# Calculate likelihood ratio for each of 293 individuals
ndist.ratio <- log10(north.vector.ndist/south.vector.ndist)
sdist.ratio <- log10(north.vector.sdist/south.vector.sdist)
ass.ratio <- log10(larvs.assign$north.vector/larvs.assign$south.vector)

sapply(ass.ratio,function(z){abline(v=z)})

# P-values given northern distribution
ndist.pvalues <- matrix(nrow = 293, ncol = 1)
for (b in 1:length(ass.ratio)){
  if(ass.ratio[b] > median(ndist.ratio)){
    ndist.pvalues[b] <- round((length(which(ndist.ratio > ass.ratio[b]))+1)/(length(ndist.ratio)+1), 5)
  } else{
    ndist.pvalues[b] <- round((length(which(ndist.ratio < ass.ratio[b]))+1)/(length(ndist.ratio)+1), 5)
  }
}

# How many p-values are < 0.05? Expect 0.05 * 293 = 15
length(which(ndist.pvalues < 0.05)) # 66 for 1000 simulations
larvs.assign[which(ndist.pvalues < 0.05),] # which larvae are these that are unlikely to come from 'northern' parents?

# P-values given southern distribution
sdist.pvalues <- matrix(nrow = 293, ncol = 1)
for (b in 1:length(ass.ratio)){
  if(ass.ratio[b] > median(sdist.ratio)){
    sdist.pvalues[b] <- round((length(which(sdist.ratio > ass.ratio[b]))+1)/(length(sdist.ratio)+1), 5)
  } else {
    sdist.pvalues[b] <- round((length(which(sdist.ratio < ass.ratio[b]))+1)/(length(sdist.ratio)+1), 5)
  }
}

# How many p-values are < 0.05? Expect 0.05 * 293 = 15
length(which(sdist.pvalues < 0.05)) # 36 for 1000 simulations
larvs.assign[which(sdist.pvalues < 0.05),] # which larvae are these that are unlikely to come from 'southern' parents?

# Adjusted p-values: none after adjustment at alpha = 0.05!! Some at alpha = 0.2 though.
ndist.pvalues.adj <- round(p.adjust(ndist.pvalues, method = "BH"), 3)
which(ndist.pvalues.adj < 0.20)

sdist.pvalues.adj <- round(p.adjust(sdist.pvalues, method = "BH"), 3)
which(sdist.pvalues.adj < 0.05)



# Mean absolute difference between 10 outlier loci in the north and south is ~10%
n.simulated <- rbinom(n = 90, size = 135, 0.8)/135 # north adults
s.simulated <- rbinom(n = 90, size = 97, 0.7)/97 # south adults
both.simulated <- rbind(n.simulated, s.simulated)
loci100 <- cbind(pop.allele.freqs.odds, both.simulated) # cbind these simulated allele frequencies to the 10 real adult ones

# Sample from each regional allele frequency distribution to create larval allele counts


###############################################################################################
#### LDA based on ingress site ####
library(tidyverse)
library(caret)
library(MASS)
library(ggplot2)

otoliths$Index <- rownames(otoliths)

# Split data into training (40%) and test sets (60%) 
nc.index <- which(otoliths$Location == "NC") %>% createDataPartition(p=0.4, list = FALSE)
roosevelt.index <- which(otoliths$Location == "Roosevelt") %>% createDataPartition(p=0.4, list = FALSE)
rumfs.index <- which(otoliths$Location == "RUMFS") %>% createDataPartition(p=0.4, list = FALSE)
york.index <- which(otoliths$Location == "York") %>% createDataPartition(p=0.4, list = FALSE)

nc.training <- otoliths[which(otoliths$Location == "NC"), "Index"][nc.index]
roosevelt.training <- otoliths[which(otoliths$Location == "Roosevelt"), "Index"][roosevelt.index]
rumfs.training <- otoliths[which(otoliths$Location == "RUMFS"), "Index"][rumfs.index]
york.training <- otoliths[which(otoliths$Location == "York"), "Index"][roosevelt.index]

training.samples <- sort(as.numeric(c(nc.training, roosevelt.training, rumfs.training, york.training)), decreasing = FALSE)

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
mean(predictions$class==test.transformed$Location) # Assignment accuracy to ingress site: LDA (64%) seems to be slightly more accurate than QDA (57%)
table(predictions$class, test.transformed$Location)

lda.data <- cbind(train.transformed, predict(model)$x)
ggplot(lda.data, aes(LD1, LD2)) +
  geom_point(aes(color = Location))

# Match predicted class with index
cbind(otoliths$Index[-training.samples], as.character(predictions$class))

#### LDA based on ingress region ####
region1 <- gsub("NC", "South", otoliths$Location) # create region (north or south) column based on ingress site
region2 <- gsub("York", "South", region1)
region3 <- gsub("RUMFS", "North", region2)
region4 <- gsub("Roosevelt", "North", region3)
otoliths$region <- as.factor(region4)

train.data2 <- otoliths[training.samples, c(23,12:21)]
test.data2 <- otoliths[-training.samples, c(23,12:21)]

# Normalize the data
# Estimate preprocessing parameters
preproc.param <- train.data2 %>% 
  preProcess(method = c("center", "scale"))
# Transform data using estimated parameters
train.transformed2 <- preproc.param %>% predict(train.data2)
test.transformed2 <- preproc.param %>% predict(test.data2)

# Fit the model
model2 <- lda(region~., data = train.transformed2)
plot(model2)
# Make predictions
predictions2 <- model2 %>% predict(test.transformed2)
# Model accuracy
mean(predictions2$class==test.transformed2$region) # Assignment accuracy to ingress site: LDA (77%)
table(predictions2$class,test.transformed2$region) 

# Match predicted class with index
test.predicted.regions <- cbind(otoliths$Index[-training.samples], as.character(predictions2$class)) # test data, regions
training.regions <- cbind(otoliths$Index[training.samples],as.character(otoliths$region[training.samples])) # training data, regions

lda.class <- as.data.frame(rbind(test.predicted.regions, training.regions))
lda.class[,1] <- as.numeric(lda.class[,1])
lda.class.ordered <- lda.class[order(lda.class[,1]),]


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

#### Back of the envelope calculations regarding power (# of loci) ####
# Read in dataset containing outlier loci and otolith data
oto.gen.merge2 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/oto.gen.merged151.txt", header = TRUE)

oto.gen.merge2

colSums(oto.gen.merge2[,16:35], na.rm = TRUE)/(2*colSums(!is.na(oto.gen.merge2[,16:35]))) 

# Sample 90 more adult allele frequencies based on the regional allele frequencies
# Mean absolute difference between 10 outlier loci in the north and south is ~10%
n.simulated <- rbinom(n = 90, size = 135, 0.8)/135 # north adults
s.simulated <- rbinom(n = 90, size = 97, 0.7)/97 # south adults
both.simulated <- rbind(n.simulated, s.simulated)
loci100 <- cbind(pop.allele.freqs.odds, both.simulated) # cbind these simulated allele frequencies to the 10 real adult ones

# Sample from each regional allele frequency distribution to create larval allele counts



# Then randomly draw 0, 1 or 2 for each of 10 loci for each larvae in each cluster
simulated.counts <- sample(0:2,13590, replace = TRUE) # 90*151
simulated.counts.matrix <- matrix(simulated.counts, nrow = 151, ncol = 90)
counts100 <- cbind(indiv.allele.counts.odds, simulated.counts.matrix) # cbind these simulated counts to the 10 real ones for 151 larval fish

# For loop to loop through each of 100 loci & multiply by the adult allele frequency, and then do this for all 151 larvae. NA's/9's get coded as 1's so they don't make a difference when each row's product is taken
north.likelihoods.100loci <- data.frame()
south.likelihoods.100loci <- data.frame()

for (j in 1:length(rownames(counts100))){
  
  for (i in 1:length(colnames(counts100))){
    if(counts100[j,i] == 2) {
      north.likelihoods.100loci[j,i] <- loci100[1,i]^2
    } else if (counts100[j,i] == 1) {
      north.likelihoods.100loci[j,i] <- 2*(loci100[1,i] * (1-loci100[1,i]))
    } else if (counts100[j,i] == 0) {
      north.likelihoods.100loci[j,i] <- ( 1-loci100[1,i])^2 
    } else {
      north.likelihoods.100loci[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(counts100))){
    if(counts100[j,i] == 2){
      south.likelihoods.100loci[j,i] <- loci100[2,i]^2
    } else if (counts100[j,i] == 1) {
      south.likelihoods.100loci[j,i] <- 2*(loci100[2,i] * (1-loci100[2,i]))
    } else  if (counts100[j,i] == 0) {
      south.likelihoods.100loci[j,i] <- (1-loci100[2,i])^2
    } else {
      south.likelihoods.100loci[j,i] <- 1
    }
  }
}

# Multiply everything together
north.vector.100loci <- vector()
south.vector.100loci <- vector()

for (k in 1:length(north.likelihoods.100loci[,1])){
  north.vector.100loci[k] <- prod(north.likelihoods.100loci[k,])
}

for (l in 1:length(south.likelihoods.100loci[,1])){
  south.vector.100loci[l] <- prod(south.likelihoods.100loci[l,])
}

# Aggregate by k-means clustering groups
# Remember log10(M*N) = log10(M) + log10(N)
north.vector.100loci.log <- log10(north.vector.100loci)
south.vector.100loci.log <- log10(south.vector.100loci)

north.pop.100loci <- aggregate(north.vector.100loci.log, by = list(oto.gen.merge2$cluster), FUN = sum) # northern likelihoods for 5 otolith clusters using 100 loci
south.pop.100loci <- aggregate(south.vector.100loci.log, by = list(oto.gen.merge2$cluster), FUN = sum) # southern likelihoods for 5 otolith clusters using 100 loci

pops.likelihood.100loci <- as.data.frame(cbind(as.numeric(north.pop.100loci[,2]), as.numeric(south.pop.100loci[,2])))
colnames(pops.likelihood.100loci) <- c("n.likes", "s.likes")
rownames(pops.likelihood.100loci) <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5")

# Plot likelihoods for region/season combos
library(wesanderson)
col.palette <- wes_palette("Darjeeling1", 5, type = "discrete")
palette(col.palette)
plot(pops.likelihood.100loci[,'s.likes'] ~ pops.likelihood.100loci[,'n.likes'], ylab = 'Southern log likelihood', xlab = 'Northern log likelihood', col = as.factor(rownames(pops.likelihood.100loci)), pch = 19)
abline(a=0,b=1)
legend("bottomright",
       legend=c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5"),
       pch=c(19),
       col=as.factor(rownames(pops.likelihood.100loci)))

#### Code I previously tried, but now is probably the wrong approach/obscure ####
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

# plot(south.log ~ north.log, xlab = 'log likelihood (north)', ylab = 'log likelihood (south)', col=ifelse(geo.color == 1, 'blue', 'tomato')) # blue is north, red is south
# abline(a = 0,b=1)
# legend("topleft",
# legend=c("Northern otolith populations", "Southern otolith populations"),
# pch=c(1, 1),
# col=c("blue", "tomato"))

# plot(as.factor(assignment[-c(which(south.vector == 1))]) ~ oto.gen.merge2$predicted[-c(which(south.vector == 1))], xlab = "Predicted population based on otolith microchemistry", ylab = "Population assignment based on genetic likelihood")
