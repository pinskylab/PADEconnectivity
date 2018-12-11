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

# MANOVA
otoliths.sub.log
oto10 <- cbind(locations, otoliths.sub.log[,-c(1:3)])
fit <- manova(as.matrix(oto10[,-1]) ~ oto10$locations)
summary(fit)
summary.aov(fit) # Mg, Mn, Fe & Sn are significantly different between ingress sites


#### Multidimensional scaling ####
# oto.chem <- otoliths[-28,c("Fish.ID", "Mg", "Mn", "Fe", "Cu", "Sr", "Cd", "Ba", "Sn", "Pb", "U")] # without PADE12_014 (194)? Scaling might take care of this. Without NCPD 074 (28)
oto.chem <- otoliths[,c("Fish.ID", "Mg", "Mn", "Fe", "Sn")] # why use all the elements when only 4 are different between regions?
rownames(oto.chem) <- oto.chem[, "Fish.ID"] # Make fish IDs as rownames
oto.chem2 <- otoliths[,c("Fish.ID", "Location", "Period", "Mg", "Mn", "Fe", "Sn")] # including time period as a variable
rownames(oto.chem2) <- oto.chem2[, "Fish.ID"] # Make fish IDs as rownames
oto.chem2$Period <- gsub("Early", "1", oto.chem2$Period) # Assign numbers to time periods: 1 = Early, 2 = Mid, 3 = Late
oto.chem2$Period <- gsub("Mid", "2", oto.chem2$Period)
oto.chem2$Period <- gsub("Late", "3", oto.chem2$Period)
oto.chem2$Period <- as.numeric(oto.chem2$Period) # convert time period to numeric for scaling later

# Separate MDS for each time period
oto.chem.early <- oto.chem2[which(oto.chem2$Period == 1),]
oto.chem.middle <- oto.chem2[which(oto.chem2$Period == 2),]
oto.chem.late <- oto.chem2[which(oto.chem2$Period == 3),]

# Scatter plot matrix
loc.cols <- col.palette[as.numeric(otoliths$Location)]
pairs(oto.chem[,-1], col = loc.cols, lower.panel = NULL)
par(xpd = TRUE)
legend(x = 0.05, y = 0.35, cex = 2,
       legend = as.character(levels(otoliths$Location)),
       fill = col.palette)
par(xpd = NA)

oto.chem <- scale(oto.chem[, -1]) # scaling probably a good idea since there is a big range of values for each elemental ratio
oto.chem2 <- scale(oto.chem2[,-1])

oto.chem.early2 <- scale(oto.chem.early[,-c(1:3)])
oto.chem.middle2 <- scale(oto.chem.middle[,-c(1:3)])
oto.chem.late2 <- scale(oto.chem.late[,-c(1:3)])

# Euclidian distance between rows
oto.dist <- dist(oto.chem)
oto.dist2 <- dist(oto.chem2)

oto.chem.early2 <- dist(oto.chem.early2)
oto.chem.middle2 <- dist(oto.chem.middle2)
oto.chem.late2 <- dist(oto.chem.late2)

# Fit model
oto.fit <- cmdscale(oto.dist, eig = TRUE, k = 3, add = FALSE) # How many dimensions?
oto.fit2 <- cmdscale(oto.dist2, eig = TRUE, k = 3, add = FALSE)

oto.fit.early <- cmdscale(oto.chem.early2, eig = TRUE, k = 3, add = FALSE)
oto.fit.middle <- cmdscale(oto.chem.middle2, eig = TRUE, k = 3, add = FALSE)
oto.fit.late <- cmdscale(oto.chem.late2, eig = TRUE, k = 3, add = FALSE)

plot(oto.fit$eig[1:10]) # scree plot
plot(oto.fit2$eig[1:10])

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
x <- oto.fit2$points[,1]
y <- oto.fit2$points[,2]
z <- oto.fit2$points[,3]

plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Otolith microchemistry")
text(x, y, labels = row.names(oto.chem2), pos = 3)
plot(x, z, xlab = "Coordinate 1", ylab = "Coordinate 3", main = "Otolith microchemistry")
text(x, z, labels = row.names(oto.chem2), pos = 2)

# Plot locations in color
library(wesanderson)

par(mfrow = c(1,2))
col.palette <- wes_palette("FantasticFox1", 5, type = "discrete")[-1]
palette(col.palette)
# plot(x, y, xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry", col = otoliths$Location[-28], pch = 19) # without NCPD 074 (28)
plot(x, y, xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry", col = otoliths$Location, pch = 19) # mds using Mg, Mn, Fe & Sn
legend("topleft",
       legend = levels(otoliths$Location),
       pch=19,
       col = col.palette)

# plot(x, z, xlab = "MDS1", ylab = "MDS3", main = "Otolith microchemistry", col = otoliths$Location[-28], pch = 19) # without NCPD 074 (28)
plot(x, z, xlab = "MDS1", ylab = "MDS3", main = "Otolith microchemistry", col = otoliths$Location, pch = 19) # mds using Mg, Mn, Fe & Sn
legend("topleft",
       legend = levels(otoliths$Location),
       pch=19,
       col = col.palette)

# Plot points by time period
# plot(x, y, xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry", col = otoliths$Period[-28], pch = 19) # without NCPD 074 (28)
plot(x, y, xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry", col = otoliths$Period, pch = 19) # mds using Mg, Mn, Fe & Sn
legend("topleft",
       legend = levels(otoliths$Period),
       pch=19,
       col = col.palette)

# plot(x, z, xlab = "MDS1", ylab = "MDS3", main = "Otolith microchemistry", col = otoliths$Period[-28], pch = 19) # without NCPD 074 (28)
plot(x, z, xlab = "MDS1", ylab = "MDS3", main = "Otolith microchemistry", col = otoliths$Period, pch = 19) # mds using Mg, Mn, Fe & Sn
legend("topleft",
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

# Plot separate MDS for each time period
library(wesanderson)

par(mfrow = c(1,2))
col.palette <- wes_palette("FantasticFox1", 5, type = "discrete")[-1]
palette(col.palette)

# Early time period
x <- oto.fit.early$points[,1]
y <- oto.fit.early$points[,2]
z <- oto.fit.early$points[,3]

plot(x, y, xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry: early", col = oto.chem.early$Location, pch = 19) # mds using Mg, Mn, Fe & Sn
legend("topleft",
       legend = levels(oto.chem.early$Location),
       pch=19,
       col = col.palette)
plot(x, z, xlab = "MDS1", ylab = "MDS3", main = "Otolith microchemistry: early", col = oto.chem.early$Location, pch = 19) # mds using Mg, Mn, Fe & Sn

# Middle time period
x <- oto.fit.middle$points[,1]
y <- oto.fit.middle$points[,2]
z <- oto.fit.middle$points[,3]

plot(x, y, xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry: middle", col = oto.chem.middle$Location, pch = 19) # mds using Mg, Mn, Fe & Sn
legend("topleft",
       legend = levels(oto.chem.middle$Location),
       pch=19,
       col = col.palette)
plot(x, z, xlab = "MDS1", ylab = "MDS3", main = "Otolith microchemistry: middle", col = oto.chem.middle$Location, pch = 19) # mds using Mg, Mn, Fe & Sn

# Late time period
x <- oto.fit.late$points[,1]
y <- oto.fit.late$points[,2]
z <- oto.fit.late$points[,3]

plot(x, y, xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry: late", col = oto.chem.late$Location, pch = 19) # mds using Mg, Mn, Fe & Sn
legend("topleft",
       legend = levels(oto.chem.late$Location),
       pch=19,
       col = col.palette)
plot(x, z, xlab = "MDS1", ylab = "MDS3", main = "Otolith microchemistry: late", col = oto.chem.late$Location, pch = 19) # mds using Mg, Mn, Fe & Sn

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

# MDS of 4 elements (Mg, Mn, Sn and Fe) & time period and colored by cluster
oto.gen.merge4 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/oto.gen.merged151.6clusters.txt", header = TRUE)

oto.chem <- oto.gen.merge4[,c("PinskyID", "Period", "Mg", "Mn", "Fe", "Sn")] # including time period as a variable
rownames(oto.chem) <- oto.chem[, "PinskyID"] # Make fish IDs as rownames
oto.chem$Period <- gsub("Early", "1", oto.chem$Period) # Assign numbers to time periods: 1 = Early, 2 = Mid, 3 = Late
oto.chem$Period <- gsub("Mid", "2", oto.chem$Period)
oto.chem$Period <- gsub("Late", "3", oto.chem$Period)
oto.chem$Period <- as.numeric(oto.chem$Period) # convert time period to numeric for scaling later

oto.chem <- scale(oto.chem[,-1]) # Scale
oto.dist <- dist(oto.chem) # Euclidian distance between rows

# Fit model
oto.fit <- cmdscale(oto.dist, eig = TRUE, k = 3, add = FALSE)

# Plot MDS
x <- oto.fit$points[,1]
y <- oto.fit$points[,2]
z <- oto.fit$points[,3]

# Plot clusters in color
par(mfrow = c(1,2))
col.palette <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F")
palette(col.palette)
plot(x, y, xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry: ingress site", col = oto.gen.merge4$Location, pch = 19) # mds using Mg, Mn, Fe & Sn and time period
legend("topleft",
       legend = levels(oto.gen.merge4$Location),
       pch=19,
       col = col.palette)

plot(x, z, xlab = "MDS1", ylab = "MDS3", main = "Otolith microchemistry: ingress site", col = oto.gen.merge4$Location, pch = 19) # mds using Mg, Mn, Fe & Sn and time period
legend("bottomright",
       legend = levels(oto.gen.merge4$Location),
       pch=19,
       col = col.palette)

plot(x, y, xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry: time period", col = oto.gen.merge4$Period, pch = 19) # mds using Mg, Mn, Fe & Sn and time period
legend("topleft",
       legend = levels(oto.gen.merge4$Period),
       pch=19,
       col = col.palette)

plot(x, z, xlab = "MDS1", ylab = "MDS3", main = "Otolith microchemistry: time period", col = oto.gen.merge4$Period, pch = 19) # mds using Mg, Mn, Fe & Sn and time period
legend("bottomright",
       legend = levels(oto.gen.merge4$Period),
       pch=19,
       col = col.palette)

plot(x, y, xlab = "MDS1", ylab = "MDS2", main = "Otolith microchemistry: cluster", col = oto.gen.merge4$cluster6, pch = 19) # mds using Mg, Mn, Fe & Sn and time period
legend("topleft",
       legend = levels(as.factor(oto.gen.merge4$cluster6)),
       pch=19,
       col = col.palette)

plot(x, z, xlab = "MDS1", ylab = "MDS3", main = "Otolith microchemistry: cluster", col = oto.gen.merge4$cluster6, pch = 19) # mds using Mg, Mn, Fe & Sn and time period
legend("bottomright",
       legend = levels(as.factor(oto.gen.merge4$cluster6)),
       pch=19,
       col = col.palette)


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
fviz_nbclust(oto.chem2, FUN = hcut, method = "wss")
fviz_nbclust(oto.chem2, FUN = hcut, method = "silhouette")
gap_stat <- clusGap(oto.chem2, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

library(NbClust)
nb <- NbClust(oto.chem2, distance = 'euclidean', min.nc = 2, max.nc = 10, method = 'kmeans')
fviz_nbclust(nb)
dev.off()

library(mclust) 
mclust.fit <- Mclust(oto.chem2)
plot(mclust.fit)
summary(mclust.fit)

# k-means clustering
col.palette <- wes_palette("Darjeeling1", 5, type = "discrete")
palette(col.palette)
kmean.cls <- kmeans(oto.chem2, centers = 6, nstart = 50, iter.max = 10)
class.table.km <- table(otoliths$Location, kmean.cls$cluster)
class.table.km2 <- table(otoliths$Period, kmean.cls$cluster)
mosaicplot(class.table.km, color = col.palette)

fviz_cluster(kmean.cls, oto.chem2, 
             # palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07", "#D8BFD8"), # for oto.chem
             palette = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F"), # for oto.chem2
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
rect.hclust(fit, k=6)

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
Location <- otoliths.sub.log$Location
Location.ordered <- factor(Location, levels = c("NC", "York", "Roosevelt", "RUMFS"))
otoliths.sub.log.trans2 <- cbind(Location.ordered, otoliths.sub.log.trans)

# dfa1 <- lda(Location ~ Mg + Mn + Fe + Cu + Cd + Sr + Ba + Sn + Pb + U, data = otoliths.sub.log.trans, na.action = "na.omit", CV = TRUE, prior = c(1,1,1,1)/4) #Cd, Pb and U apprently too similar between groups when data not transformed and standardized; but not so when they are

Cluster <- as.factor(kmean.cls$cluster)
names(kmean.cls$cluster) == rownames(otoliths.sub.log.trans) # make sure individuals are in the right order
otoliths.sub.log.trans2 <- cbind(Cluster, otoliths.sub.log.trans)
# dfa1 <- lda(Cluster ~ Mg + Mn + Fe + Sn, data = otoliths.sub.log.trans2, na.action = "na.omit", CV = TRUE, prior = c(1,1,1,1)/4) 
dfa1 <- lda(Location.ordered ~ Mg + Mn + Fe + Sn, data = otoliths.sub.log.trans2, na.action = "na.omit", CV = TRUE, prior = c(1,1,1,1)/4)  # the jack-knifing doesn't result in coordinates for plotting

# Assess accuracy of the prediction
# percent correct for each category of Location
ct1 <- table(otoliths.sub.log.trans2$Location.ordered, dfa1$class)
# ct1 <- table(otoliths.sub.log.trans2$Cluster, dfa1$class)
diag(prop.table(ct1,1))

# total percent correct
sum(diag(prop.table(ct1)))

# barplot of assignments
props <- prop.table(ct1,1) # 'origin' signature/total number of fished that ingressed at a site
library(wesanderson)
col.palette <- wes_palette("FantasticFox1", 5, type = "discrete")[-1]
barplot(props, horiz = TRUE, beside = TRUE, xlim = c(0,1), col = col.palette, xlab = "Assignment proportion", ylab = "Predicted 'origin' signature")
legend("bottomright",
       legend=rev(levels(otoliths.sub.log.trans2$Location)),
       pch=22,
       col = 'black',
       pt.bg= rev(col.palette),
       title = expression(bold('Collection location')), 
       bty = "n")

#### DFA using 67% of core points to 'train' DFA, then rerun with only these individuals, and then assign everybody ####
dfa2 <- lda(Location.ordered ~ Mg + Mn + Fe + Sn, data = otoliths.sub.log.trans2, prior = c(1,1,1,1)/4)
plot(dfa2)
# ld1 <- dfa2$scaling[1,1]*otoliths.sub.log.trans2$Mg + dfa2$scaling[2,1]*otoliths.sub.log.trans2$Mn + dfa2$scaling[3,1]*otoliths.sub.log.trans2$Fe + dfa2$scaling[4,1]*otoliths.sub.log.trans2$Sn
# ld2 <- dfa2$scaling[1,2]*otoliths.sub.log.trans2$Mg + dfa2$scaling[2,2]*otoliths.sub.log.trans2$Mn + dfa2$scaling[3,2]*otoliths.sub.log.trans2$Fe + dfa2$scaling[4,2]*otoliths.sub.log.trans2$Sn

dfa.values <- predict(dfa2) # Calculates linear discriminants, as above
dfa.values2 <- cbind.data.frame(dfa.values$x, Location)

plot(dfa.values2$LD1, dfa.values2$LD2, col = dfa.values2$Location)
legend("topleft",
       legend=levels(Location),
       pch=19,
       col = col.palette)

# Calculate center for each ingress site
center <- aggregate(dfa.values2[,-4], by = list(dfa.values2$Location), FUN = mean) # similar but not exactly the same as center calculation below

# Draw data ellipses for each ingress site
loc.groups <- split(dfa.values2, dfa.values2$Location)
dataEllipse(loc.groups$NC[,"LD1"], loc.groups$NC[,"LD2"], levels = 0.67, xlim = c(-3,3)) # NC
dataEllipse(loc.groups$York[,"LD1"], loc.groups$York[,"LD2"], levels = 0.67) # York
dataEllipse(loc.groups$Roosevelt[,"LD1"], loc.groups$Roosevelt[,"LD2"], levels = 0.67) # Roosevelt
dataEllipse(loc.groups$RUMFS[,"LD1"], loc.groups$RUMFS[,"LD2"], levels = 0.67) # RUMFS

# Fit ellipse for each ingress location
eli.nc <- ellipse(cor(loc.groups$NC[,"LD1"], loc.groups$NC[,"LD2"]), scale=c(sd(loc.groups$NC[,"LD1"]),sd(loc.groups$NC[,"LD2"])), centre=c(mean(loc.groups$NC[,"LD1"]), mean(loc.groups$NC[,"LD2"])), level = 0.67, npoints = 250)
eli.york <- ellipse(cor(loc.groups$York[,"LD1"], loc.groups$York[,"LD2"]), scale=c(sd(loc.groups$York[,"LD1"]),sd(loc.groups$York[,"LD2"])), centre=c(mean(loc.groups$York[,"LD1"]), mean(loc.groups$York[,"LD2"])), level = 0.67, npoints = 250)
eli.roosevelt <- ellipse(cor(loc.groups$Roosevelt[,"LD1"], loc.groups$Roosevelt[,"LD2"]), scale=c(sd(loc.groups$Roosevelt[,"LD1"]),sd(loc.groups$Roosevelt[,"LD2"])), centre=c(mean(loc.groups$Roosevelt[,"LD1"]), mean(loc.groups$Roosevelt[,"LD2"])), level = 0.67, npoints = 250)
eli.rumfs <- ellipse(cor(loc.groups$RUMFS[,"LD1"], loc.groups$RUMFS[,"LD2"]), scale=c(sd(loc.groups$RUMFS[,"LD1"]),sd(loc.groups$RUMFS[,"LD2"])), centre=c(mean(loc.groups$RUMFS[,"LD1"]), mean(loc.groups$RUMFS[,"LD2"])), level = 0.67, npoints = 250)

#Calculate the center of ellipse for each location
eli_center_nc = c(mean(eli.nc[,1]), mean(eli.nc[,2]))
eli_center_york = c(mean(eli.york[,1]), mean(eli.york[,2]))
eli_center_roosevelt = c(mean(eli.roosevelt[,1]), mean(eli.roosevelt[,2]))
eli_center_rumfs = c(mean(eli.rumfs[,1]), mean(eli.rumfs[,2]))

#A function to calculate distance between points 'x1' and 'x2'
dist_2_points <- function(x1, x2) {
  return(sqrt(sum((x1 - x2)^2)))    
}

#Compute distance of each point in ellipse from eli_center
# NC
distance.nc <- numeric(0)
for (i in 1:nrow(eli.nc)){ 
  distance.nc[i] = dist_2_points(eli_center_nc, eli.nc[i,])
}

# York
distance.york <- numeric(0)
for (i in 1:nrow(eli.york)){ 
  distance.york[i] = dist_2_points(eli_center_york, eli.york[i,])
}

# Roosevelt
distance.roosevelt <- numeric(0)
for (i in 1:nrow(eli.roosevelt)){ 
  distance.roosevelt[i] = dist_2_points(eli_center_roosevelt, eli.roosevelt[i,])
}

# RUMFS
distance.rumfs <- numeric(0)
for (i in 1:nrow(eli.rumfs)){ 
  distance.rumfs[i] = dist_2_points(eli_center_rumfs, eli.rumfs[i,])
}

#The maximum distance from eli_center is 'a'
a.nc <- distance.nc[which.max(distance.nc)]
a.york <- distance.york[which.max(distance.york)]
a.roosevelt <- distance.roosevelt[which.max(distance.roosevelt)]
a.rumfs <- distance.rumfs[which.max(distance.rumfs)]

#The minimum distance from eli_center is 'b'
b.nc <- distance.nc[which.min(distance.nc)]
b.york <- distance.york[which.min(distance.york)]
b.roosevelt <- distance.roosevelt[which.min(distance.roosevelt)]
b.rumfs <- distance.rumfs[which.min(distance.rumfs)]

# Figure out if points are inside or outside ellipses
NC.out <- loc.groups$NC[which(round((((loc.groups$NC[,"LD1"] - eli_center_nc[1])^2)/(b.nc)^2) + (((loc.groups$NC[,"LD2"] - eli_center_nc[2])^2)/(a.nc)^2),3) >= 1.000),]
NC.in <- loc.groups$NC[-which(round((((loc.groups$NC[,"LD1"] - eli_center_nc[1])^2)/(b.nc)^2) + (((loc.groups$NC[,"LD2"] - eli_center_nc[2])^2)/(a.nc)^2),3) >= 1.000),]

York.out <- loc.groups$York[which(round((((loc.groups$York[,"LD1"] - eli_center_york[1])^2)/(a.york)^2) + (((loc.groups$York[,"LD2"] - eli_center_york[2])^2)/(b.york)^2),3) >= 1.000),]
York.in <- loc.groups$York[-which(round((((loc.groups$York[,"LD1"] - eli_center_york[1])^2)/(a.york)^2) + (((loc.groups$York[,"LD2"] - eli_center_york[2])^2)/(b.york)^2),3) >= 1.000),]

Roosevelt.out <- loc.groups$Roosevelt[which(round((((loc.groups$Roosevelt[,"LD1"] - eli_center_roosevelt[1])^2)/(b.roosevelt)^2) + (((loc.groups$Roosevelt[,"LD2"] - eli_center_roosevelt[2])^2)/(a.roosevelt)^2),3) >= 1.000),]
Roosevelt.in <- loc.groups$Roosevelt[-which(round((((loc.groups$Roosevelt[,"LD1"] - eli_center_roosevelt[1])^2)/(b.roosevelt)^2) + (((loc.groups$Roosevelt[,"LD2"] - eli_center_roosevelt[2])^2)/(a.roosevelt)^2),3) >= 1.000),]

RUMFS.out <- loc.groups$RUMFS[which(round((((loc.groups$RUMFS[,"LD1"] - eli_center_rumfs[1])^2)/(a.rumfs)^2) + (((loc.groups$RUMFS[,"LD2"] - eli_center_rumfs[2])^2)/(b.rumfs)^2),3) >= 1.000),]
RUMFS.in <- loc.groups$RUMFS[-which(round((((loc.groups$RUMFS[,"LD1"] - eli_center_rumfs[1])^2)/(a.rumfs)^2) + (((loc.groups$RUMFS[,"LD2"] - eli_center_rumfs[2])^2)/(b.rumfs)^2),3) >= 1.000),]

# Combine Fish IDs of data that was within 67% confidence ellipses
in.67 <- c(rownames(NC.in), rownames(York.in), rownames(Roosevelt.in), rownames(RUMFS.in))

# Now subset elemental data to only these fish
sub.67 <- otoliths.sub.log.trans2[rownames(otoliths.sub.log.trans2) %in% in.67,]


#### Add predicted sites to the otolith data ####
# otoliths$predicted <- dfa1$class # using DFA classes
otoliths$cluster <- kmean.cls$cluster # using heirarchical clustering
otoliths$cluster8 <- kmean.cls$cluster
otoliths$cluster6 <- kmean.cls$cluster

names(kmean.cls$cluster) == otoliths$Fish.ID
# rownames(otoliths)==lda.class.ordered[,1] # should be in same order
# otoliths$predicted <- lda.class.ordered[,2] # from LDA using random 40% from each ingress site as test dataset

# Read in dataset containing outlier loci
gen.larvae.outs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/masterPADElarvae.txt', header = TRUE, sep = "\t")

# Keep only larvae with genetic assignments
gen.larvae.outs2 <- gen.larvae.outs[which(is.na(gen.larvae.outs$assignment) == FALSE), ]

# Merge otolith data with predicted sites with genetic data containing outliers
oto.gen.merge2 <- merge(gen.larvae.outs2, otoliths[,-11], by.x = 'PicID', by.y = 'Fish.ID', all = FALSE) # merged otolith and genetic data set; remove column of NAs in otolith data
# write.table(oto.gen.merge2, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/oto.gen.merged151.txt", col.names = TRUE, row.names = FALSE)
oto.gen.merge2 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/oto.gen.merged151.txt", header = TRUE)

oto.gen.merge3 <- merge(gen.larvae.outs2, otoliths[,-11], by.x = 'PicID', by.y = 'Fish.ID', all = FALSE) # merged otolith and genetic data set; remove column of NAs in otolith data
# write.table(oto.gen.merge3, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/oto.gen.merged151.8clusters.txt", col.names = TRUE, row.names = FALSE)
oto.gen.merge3 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/oto.gen.merged151.8clusters.txt", header = TRUE)

oto.gen.merge4 <- merge(gen.larvae.outs2, otoliths[,-11], by.x = 'PicID', by.y = 'Fish.ID', all = FALSE) # merged otolith and genetic data set; remove column of NAs in otolith data
# write.table(oto.gen.merge4, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/oto.gen.merged151.6clusters.txt", col.names = TRUE, row.names = FALSE)
oto.gen.merge4 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/oto.gen.merged151.6clusters.txt", header = TRUE)

# Pie charts of six clusters based on otolith elemental cores
oto.gen.merge4 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/oto.gen.merged151.6clusters.txt", header = TRUE)
library(wesanderson)
col.palette <- wes_palette("FantasticFox1", 5, type = "discrete")[-1]
palette(col.palette)

clusters <- split(oto.gen.merge4, oto.gen.merge4$cluster6)

# Plot together
par(mfrow = c(3,2),
    oma = c(5,4,0,1) +0.1,
    mar = c(0,0,1,1) + 0.1)
pie(table(clusters[[1]]$Location), col = col.palette, labels = '', main = 'cluster 1')
pie(table(clusters[[2]]$Location), col = col.palette, labels = '', main = 'cluster 2')
pie(table(clusters[[3]]$Location), col = col.palette, labels = '', main = 'cluster 3')
pie(table(clusters[[4]]$Location), col = col.palette, labels = '', main = 'cluster 4')
pie(table(clusters[[5]]$Location), col = col.palette, labels = '', main = 'cluster 5')
pie(table(clusters[[6]]$Location), col = col.palette, labels = '', main = 'cluster 6')
legend("bottomright",
       legend=levels(oto.gen.merge4$Location),
       pch=22,
       col = 'black',
       pt.bg= col.palette)

# Plot separately
# Cluster 1
pie(table(clusters[[1]]$Location), col = col.palette, labels = '')
mtext('cluster 1 (n = 23)', side = 3, line = -4.5)
text(0,0.4, "NC")
text(0.3,-0.4, "York")
text(-0.55,-0.75, "Roosevelt")

# Cluster 2
pie(table(clusters[[2]]$Location), col = col.palette, labels = '')
mtext('cluster 2 (n = 32)', side = 3, line = -4.5)
text(0.1,0.4, "NC")
text(-0.45, 0.05, "Roosevelt")
text(-0.79,-0.5, "RUMFS")
text(0.2,-0.35, "York")

# Cluster 3
pie(table(clusters[[3]]$Location), col = col.palette, labels = '')
mtext('cluster 3 (n = 8)', side = 3, line = -4.5)
text(0,0.4, "NC")
text(0.1,-0.35, "York")

# Cluster 4
pie(table(clusters[[4]]$Location), col = col.palette, labels = '')
mtext('cluster 4 (n = 46)', side = 3, line = -4.5)
text(0.1,0.4, "NC")
text(-0.1,-0.35, "RUMFS")

# Cluster 5
pie(table(clusters[[5]]$Location), col = col.palette, labels = '')
mtext('cluster 5 (n = 23)', side = 3, line = -4.5)
text(0.9,0.13, "NC")
text(-0.3, 0.05, "Roosevelt")
text(0.4,-0.3, "RUMFS")
text(0.9,-0.13, "York")

# Cluster 6
pie(table(clusters[[6]]$Location), col = col.palette, labels = '')
mtext('cluster 6 (n = 19)', side = 3, line = -4.5)
text(-0.3, 0.1, "Roosevelt")
text(0.4,-0.2, "RUMFS")




# Make some plots
# plot(oto.gen.merge2$assignment ~ oto.gen.merge2$predicted, xlab = "Predicted population based on otolith microchemistry", ylab = "Population assignment based on genetic likelihood")
plot(oto.gen.merge2$assignment ~ oto.gen.merge2$cluster, xlab = "Cluster based on otolith microchemistry", ylab = "Population assignment based on genetic likelihood")
plot(oto.gen.merge3$assignment ~ oto.gen.merge3$cluster8, xlab = "Cluster based on otolith microchemistry", ylab = "Population assignment based on genetic likelihood")

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
# pop.allele.freqs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/regional.adult.outs.freqs10.txt')

pop.allele.freqs5 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.freqs.5pops.txt')

# Do these calculated frequencies match?
colnames(regional.allele.freqs) == colnames(pop.allele.freqs) # alleles for one locus are flipped, numbers seem same, but not TRUE for some reason....?
as.data.frame(regional.allele.freqs) == pop.allele.freqs

# Pull out odd indicies
pop.allele.freqs
pop.allele.freqs5
odds <- seq(1,20,2) # odd indicies to keep
# pop.allele.freqs.odds <- pop.allele.freqs[,odds] # north & south allele frequencies
pop.allele.freqs5.odds <- pop.allele.freqs5[,odds] # allele frequencies of 5 bayenv populations

#### Population assignment using allele frequencies of otolith populations ####
indiv.allele.counts <- oto.gen.merge2[,c(16:35)] # Just genetic data
indiv.allele.counts <- oto.gen.merge3[,c(16:35)] # Just genetic data
indiv.allele.counts <- oto.gen.merge4[,c(16:35)] # Just genetic data

# Two ways to do this:
indiv.allele.counts.odds <- indiv.allele.counts[,odds] # first way, but this requires the column names in the count data and the allele frequency data to be the same
indiv.allele.counts.odds[is.na(indiv.allele.counts.odds)] <- 9 # replace NA's with 9's to make the ifelse statements easier

indiv.allele.counts.odds <- indiv.allele.counts[, colnames(indiv.allele.counts) %in% colnames(pop.allele.freqs5.odds)] # probably the better way to do this
indiv.allele.counts.odds[is.na(indiv.allele.counts.odds)] <- 9 # replace NA's with 9's to make the ifelse statements easier

colnames(indiv.allele.counts.odds) == colnames(pop.allele.freqs5.odds) # column names match? yay!

#### Calculating likelihoods for each individual. First part assumes two populations (north and south). Second assumes five BayEnv populations.
# For loop to loop through each locus & multiply by the otolith allele frequency, and then do this for all 151 larvae. NA's/9's get coded as 1's so they don't make a difference when each row's product is taken
# Just north and south likelihood calculations
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

# Five BayEnv population likelihoods
pop1.likelihoods <- data.frame()
pop2.likelihoods <- data.frame()
pop3.likelihoods <- data.frame()
pop4.likelihoods <- data.frame()
pop5.likelihoods <- data.frame()

for (j in 1:length(rownames(indiv.allele.counts.odds))){
  
  for (i in 1:length(colnames(indiv.allele.counts.odds))){
    if(indiv.allele.counts.odds[j,i] == 2) {
      pop1.likelihoods[j,i] <- pop.allele.freqs5.odds["Pop1",i]^2
    } else if (indiv.allele.counts.odds[j,i] == 1) {
      pop1.likelihoods[j,i] <- 2*(pop.allele.freqs5.odds["Pop1",i] * (1-pop.allele.freqs5.odds["Pop1",i]))
    } else if (indiv.allele.counts.odds[j,i] == 0) {
      pop1.likelihoods[j,i] <- ( 1-pop.allele.freqs5.odds["Pop1",i])^2 
    } else {
      pop1.likelihoods[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(indiv.allele.counts.odds))){
    if(indiv.allele.counts.odds[j,i] == 2){
      pop2.likelihoods[j,i] <- pop.allele.freqs5.odds["Pop2",i]^2
    } else if (indiv.allele.counts.odds[j,i] == 1) {
      pop2.likelihoods[j,i] <- 2*(pop.allele.freqs5.odds["Pop2",i] * (1-pop.allele.freqs5.odds["Pop2",i]))
    } else  if (indiv.allele.counts.odds[j,i] == 0) {
      pop2.likelihoods[j,i] <- (1-pop.allele.freqs5.odds["Pop2",i])^2
    } else {
      pop2.likelihoods[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(indiv.allele.counts.odds))){
    if(indiv.allele.counts.odds[j,i] == 2) {
      pop3.likelihoods[j,i] <- pop.allele.freqs5.odds["Pop3",i]^2
    } else if (indiv.allele.counts.odds[j,i] == 1) {
      pop3.likelihoods[j,i] <- 2*(pop.allele.freqs5.odds["Pop3",i] * (1-pop.allele.freqs5.odds["Pop3",i]))
    } else if (indiv.allele.counts.odds[j,i] == 0) {
      pop3.likelihoods[j,i] <- ( 1-pop.allele.freqs5.odds["Pop3",i])^2 
    } else {
      pop3.likelihoods[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(indiv.allele.counts.odds))){
    if(indiv.allele.counts.odds[j,i] == 2) {
      pop4.likelihoods[j,i] <- pop.allele.freqs5.odds["Pop4",i]^2
    } else if (indiv.allele.counts.odds[j,i] == 1) {
      pop4.likelihoods[j,i] <- 2*(pop.allele.freqs5.odds["Pop4",i] * (1-pop.allele.freqs5.odds["Pop4",i]))
    } else if (indiv.allele.counts.odds[j,i] == 0) {
      pop4.likelihoods[j,i] <- ( 1-pop.allele.freqs5.odds["Pop4",i])^2 
    } else {
      pop4.likelihoods[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(indiv.allele.counts.odds))){
    if(indiv.allele.counts.odds[j,i] == 2) {
      pop5.likelihoods[j,i] <- pop.allele.freqs5.odds["Pop5",i]^2
    } else if (indiv.allele.counts.odds[j,i] == 1) {
      pop5.likelihoods[j,i] <- 2*(pop.allele.freqs5.odds["Pop5",i] * (1-pop.allele.freqs5.odds["Pop5",i]))
    } else if (indiv.allele.counts.odds[j,i] == 0) {
      pop5.likelihoods[j,i] <- ( 1-pop.allele.freqs5.odds["Pop5",i])^2 
    } else {
      pop5.likelihoods[j,i] <- 1
    }
  }
}

#### Multiply everything together now. First part is assuming 2 regions. Second is assuming five BayEnv regions ###
# Two regions
north.vector <- vector()
south.vector <- vector()

for (k in 1:length(north.likelihoods[,1])){
  north.vector[k] <- prod(north.likelihoods[k,])
}

for (l in 1:length(south.likelihoods[,1])){
  south.vector[l] <- prod(south.likelihoods[l,])
}

# Five BayEnv regions
pop1.vector <- apply(pop1.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
pop2.vector <- apply(pop2.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
pop3.vector <- apply(pop3.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
pop4.vector <- apply(pop4.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
pop5.vector <- apply(pop5.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)

# Individual likelihoods for each BayEnv region
bayenv.likelihoods.indivs <- data.frame(oto.gen.merge4$PinskyID, oto.gen.merge4$Place, pop1.vector, pop2.vector, pop3.vector, pop4.vector, pop5.vector)
colnames(bayenv.likelihoods.indivs) <- c("ID", "Place", "Pop1", "Pop2", "Pop3", "Pop4", "Pop5")
write.table(bayenv.likelihoods.indivs, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/bayenv_likelihoods_indivs.txt", row.names = FALSE, col.names = TRUE)

# Which BayEnv population is each individual most likely from?
bayenv.likelihoods.indivs <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/bayenv_likelihoods_indivs.txt", header = TRUE)

most.like <- colnames(bayenv.likelihoods.indivs[,-c(1:2)])[apply(bayenv.likelihoods.indivs[,-c(1:2)],1, which.max)]
t(table(most.like,bayenv.likelihoods.indivs$Place))

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

north.pop <- aggregate(oto.gen.merge2$north.vector, by = list(oto.gen.merge2$cluster), FUN = prod) # northern likelihoods for 5 otolith clusters
south.pop <- aggregate(oto.gen.merge2$south.vector, by = list(oto.gen.merge2$cluster), FUN = prod) # southern likelihoods for 5 otolith clusters

north.pop <- aggregate(oto.gen.merge2$north.vector, by = list(oto.gen.merge2$cluster), FUN = prod) # same as above, but data may be coming from different places (read in vs. already stored as an object)
south.pop <- aggregate(oto.gen.merge2$south.vector, by = list(oto.gen.merge2$cluster), FUN = prod) 

bayenv1.likes <- aggregate(pop1.vector, by = list(oto.gen.merge4$cluster6), FUN = prod) # likelihoods for BayEnv groupings
bayenv2.likes <- aggregate(pop2.vector, by = list(oto.gen.merge4$cluster6), FUN = prod)
bayenv3.likes <- aggregate(pop3.vector, by = list(oto.gen.merge4$cluster6), FUN = prod)
bayenv4.likes <- aggregate(pop4.vector, by = list(oto.gen.merge4$cluster6), FUN = prod)
bayenv5.likes <- aggregate(pop5.vector, by = list(oto.gen.merge4$cluster6), FUN = prod) 

pops.likelihood <- as.data.frame(cbind(log10(as.numeric(north.pop[,2])), log10(as.numeric(south.pop[,2]))))
colnames(pops.likelihood) <- c("n.likes", "s.likes")
rownames(pops.likelihood) <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5")

bayenv.likelihoods <- as.data.frame(cbind(log10(as.numeric(bayenv1.likes[,2])), log10(as.numeric(bayenv2.likes[,2])), log10(as.numeric(bayenv3.likes[,2])), log10(as.numeric(bayenv4.likes[,2])), log10(as.numeric(bayenv5.likes[,2]))))
colnames(bayenv.likelihoods) <- c("bayenv1.likes", "bayenv2.likes", "bayenv3.likes", "bayenv4.likes", "bayenv5.likes")
rownames(bayenv.likelihoods) <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6")
write.table(bayenv.likelihoods, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/bayenv_likelihoods_6clusters.txt", row.names = TRUE, col.names = TRUE)

most.like <- colnames(bayenv.likelihoods[apply(bayenv.likelihoods,1, which.max)])

#### Plot population likelihoods of otolith populations/clusters ####
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

# Power test
log10(north.vector.ndist/south.vector.ndist)[order(log10(north.vector.ndist/south.vector.ndist))][50] #5% of 1000 is 50
p1 <- length(which(log10(north.vector.sdist/south.vector.sdist) < -0.5864056))/1000

library(pwr)
library(lsr)
pwr.2p.test(abs(mean(log10(north.vector.ndist/south.vector.ndist))-mean(log10(north.vector.sdist/south.vector.sdist))), 1000, 0.05)

cohen <- cohensD(log10(north.vector.ndist/south.vector.ndist), log10(north.vector.sdist/south.vector.sdist))
cohen <- cohensD(north.vector.ndist/south.vector.ndist, north.vector.sdist/south.vector.sdist)
cohen.larvs <- cohensD(ndist.ratio, sdist.ratio)
pwr.t.test(1000, d = cohen, type = "two.sample")
pwr.t2n.test(n1 = 97, n2 = 135, d = cohen.larvs)

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
# Read in adult outlier allele frequencies
pop.allele.freqs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/regional.adult.outs.freqs10.txt')

# Keep only one allele per locus
odds <- seq(1,20,2) # odd indicies to keep
pop.allele.freqs.odds <- pop.allele.freqs[,odds]
pop.allele.freqs.odds.100 <- do.call("cbind", replicate(10, pop.allele.freqs.odds, simplify = FALSE)) # replicate these allele frequencies 10 times
colnames(pop.allele.freqs.odds.100) <- 1:100 # rename the columns so it doesn't get super confusing
pop.allele.freqs.odds.100.north <- pop.allele.freqs.odds.100[1,] # north allele frequency
pop.allele.freqs.odds.100.south <- pop.allele.freqs.odds.100[2,] # south allele frequency

# Sample one allele per locus many times
n.alleles <- sapply(pop.allele.freqs.odds.100.north,function(z){rbinom(2000,1,z)})

# Do this for the southern distribution too
s.alleles <- sapply(pop.allele.freqs.odds.100.south,function(z){rbinom(2000,1,z)})

# Then create the genotype for each individual by separating even and odd rows, then adding them together to get allele count
odds <- seq(1,2000,2) # odd indicies
evens <- seq(2,2000,2) # even indicies

n.alleles.odds <- n.alleles[odds,]
n.alleles.evens <- n.alleles[evens,]
n.alleles.sum.100 <- n.alleles[odds,] + n.alleles[evens,]

s.alleles.odds <- s.alleles[odds,]
s.alleles.evens <- s.alleles[evens,]
s.alleles.sum.100 <- s.alleles[odds,] + s.alleles[evens,]

# Use adult allele frequencies to calculate likelihoods
colnames(pop.allele.freqs.odds.100) == colnames(n.alleles.sum.100) # check to make sure column names are the same
colnames(pop.allele.freqs.odds.100) == colnames(s.alleles.sum.100)

# For loop to loop through each of 10 loci & multiply by the adult allele frequency, and then do this for all 1000 offspring. NA's/9's get coded as 1's so they don't make a difference when each row's product is taken
north.likelihoods.ndist.100 <- data.frame()
south.likelihoods.ndist.100 <- data.frame()

for (j in 1:nrow(n.alleles.sum.100)){
  
  for (i in 1:length(colnames(n.alleles.sum.100))){
    if(n.alleles.sum.100[j,i] == 2) {
      north.likelihoods.ndist.100[j,i] <- pop.allele.freqs.odds.100[1,i]^2
    } else if (n.alleles.sum.100[j,i] == 1) {
      north.likelihoods.ndist.100[j,i] <- 2*(pop.allele.freqs.odds.100[1,i] * (1-pop.allele.freqs.odds.100[1,i]))
    } else if (n.alleles.sum.100[j,i] == 0) {
      north.likelihoods.ndist.100[j,i] <- ( 1-pop.allele.freqs.odds.100[1,i])^2 
    } else {
      north.likelihoods.ndist.100[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(n.alleles.sum.100))){
    if(n.alleles.sum.100[j,i] == 2){
      south.likelihoods.ndist.100[j,i] <- pop.allele.freqs.odds.100[2,i]^2
    } else if (n.alleles.sum.100[j,i] == 1) {
      south.likelihoods.ndist.100[j,i] <- 2*(pop.allele.freqs.odds.100[2,i] * (1-pop.allele.freqs.odds.100[2,i]))
    } else  if (n.alleles.sum.100[j,i] == 0) {
      south.likelihoods.ndist.100[j,i] <- (1-pop.allele.freqs.odds.100[2,i])^2
    } else {
      south.likelihoods.ndist.100[j,i] <- 1
    }
  }
}

north.likelihoods.sdist.100 <- data.frame()
south.likelihoods.sdist.100 <- data.frame()

for (j in 1:nrow(s.alleles.sum.100)){
  
  for (i in 1:length(colnames(s.alleles.sum.100))){
    if(s.alleles.sum.100[j,i] == 2) {
      north.likelihoods.sdist.100[j,i] <- pop.allele.freqs.odds.100[1,i]^2
    } else if (s.alleles.sum.100[j,i] == 1) {
      north.likelihoods.sdist.100[j,i] <- 2*(pop.allele.freqs.odds.100[1,i] * (1-pop.allele.freqs.odds.100[1,i]))
    } else if (s.alleles.sum.100[j,i] == 0) {
      north.likelihoods.sdist.100[j,i] <- ( 1-pop.allele.freqs.odds.100[1,i])^2 
    } else {
      north.likelihoods.sdist.100[j,i] <- 1
    }
  }
  
  for (i in 1:length(colnames(s.alleles.sum.100))){
    if(s.alleles.sum.100[j,i] == 2){
      south.likelihoods.sdist.100[j,i] <- pop.allele.freqs.odds.100[2,i]^2
    } else if (s.alleles.sum.100[j,i] == 1) {
      south.likelihoods.sdist.100[j,i] <- 2*(pop.allele.freqs.odds.100[2,i] * (1-pop.allele.freqs.odds.100[2,i]))
    } else  if (s.alleles.sum.100[j,i] == 0) {
      south.likelihoods.sdist.100[j,i] <- (1-pop.allele.freqs.odds.100[2,i])^2
    } else {
      south.likelihoods.sdist.100[j,i] <- 1
    }
  }
}

# Multiply everything together
north.vector.ndist.100 <- vector()
south.vector.ndist.100 <- vector()

for (k in 1:length(north.likelihoods.ndist.100[,1])){
  north.vector.ndist.100[k] <- prod(north.likelihoods.ndist.100[k,])
}

for (l in 1:length(south.likelihoods.ndist.100[,1])){
  south.vector.ndist.100[l] <- prod(south.likelihoods.ndist.100[l,])
}

north.vector.sdist.100 <- vector()
south.vector.sdist.100 <- vector()

for (k in 1:length(north.likelihoods.sdist.100[,1])){
  north.vector.sdist.100[k] <- prod(north.likelihoods.sdist.100[k,])
}

for (l in 1:length(south.likelihoods.sdist.100[,1])){
  south.vector.sdist.100[l] <- prod(south.likelihoods.sdist.100[l,])
}

# Create ratio & plot
hist(log10(north.vector.ndist.100/south.vector.ndist.100), xlab = "log10(north likelihood/south likelihood)", col = rgb(1,0,0,0.5), main = "100 loci", xlim = c(-12,12), ylim = c(0,250))
hist(log10(north.vector.sdist.100/south.vector.sdist.100), col = rgb(0,0,1,0.5), add = TRUE)
legend("topright", c("North", "South"), col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch = 15)

# Power test
log10(north.vector.ndist.100/south.vector.ndist.100)[order(log10(north.vector.ndist.100/south.vector.ndist.100))][50] # 5% of 1000 is 50
length(which(log10(north.vector.sdist.100/south.vector.sdist.100) < 0.39798))

library(pwr)
library(lsr)
pwr.2p.test(mean(log10(north.vector.ndist.100/south.vector.ndist.100))-mean(log10(north.vector.sdist.100/south.vector.sdist.100)), 1000, 0.05)

cohen <- cohensD(log10(north.vector.ndist.100/south.vector.ndist.100), log10(north.vector.sdist.100/south.vector.sdist.100))
pwr.t.test(1000, d = cohen, type = "two.sample")

#### Likelihood ratios using L(home)/L(max) or L(home)/L(max not home)
bayenv.likelihoods.indivs <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/bayenv_likelihoods_indivs.txt", header = TRUE)

# Create a column for the genotype likelihood of the source population
for (i in 1:length(bayenv.likelihoods.indivs$Place)){
  if(bayenv.likelihoods.indivs$Place[i] == "Little Egg Inlet, NJ"){
    bayenv.likelihoods.indivs$L_home[i] <- bayenv.likelihoods.indivs$Pop3[i]
  } else if (bayenv.likelihoods.indivs$Place[i] == "Roosevelt Inlet, DE") {
    bayenv.likelihoods.indivs$L_home[i] <- bayenv.likelihoods.indivs$Pop3[i]
  } else  if (bayenv.likelihoods.indivs$Place[i] == "York River, VA") {
    bayenv.likelihoods.indivs$L_home[i] <- bayenv.likelihoods.indivs$Pop3[i]
  } else {
    bayenv.likelihoods.indivs$L_home[i] <- bayenv.likelihoods.indivs$Pop4[i]
  }
}

for (i in 1:length(bayenv.likelihoods.indivs$Place)){
  if(bayenv.likelihoods.indivs$Place[i] == "Little Egg Inlet, NJ"){
    bayenv.likelihoods.indivs$L_home_pop[i] <- "3"
  } else if (bayenv.likelihoods.indivs$Place[i] == "Roosevelt Inlet, DE") {
    bayenv.likelihoods.indivs$L_home_pop[i] <- "3"
  } else  if (bayenv.likelihoods.indivs$Place[i] == "York River, VA") {
    bayenv.likelihoods.indivs$L_home_pop[i] <- "3"
  } else {
    bayenv.likelihoods.indivs$L_home_pop[i] <- "4"
  }
}

# Check this is correct
bayenv.likelihoods.indivs$Pop4[which(bayenv.likelihoods.indivs$Place == "Beaufort, NC")]==bayenv.likelihoods.indivs$L_home[which(bayenv.likelihoods.indivs$Place == "Beaufort, NC")]

# Now create a column for the highest genotype likelihood
bayenv.likelihoods.indivs$L_max <- apply(bayenv.likelihoods.indivs[,c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5")],1,max)
bayenv.likelihoods.indivs$L_max_pop <- apply(bayenv.likelihoods.indivs[,c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5")],1,which.max)

# Calculate L_home/L_max
bayenv.likelihoods.indivs$L_home_L_max <- log10(bayenv.likelihoods.indivs$L_home/bayenv.likelihoods.indivs$L_max)

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
