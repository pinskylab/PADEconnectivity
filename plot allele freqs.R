# Read in adult outlier allele frequencies
pop.allele.freqs5 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.freqs.5pops.txt')

# Keep only one allele per locus
odds <- seq(1,20,2) # odd indicies to keep
pop.allele.freqs5.odds <- pop.allele.freqs5[,odds]

# Transpose the data
pop.allele.freqs5.odds2 <- t(pop.allele.freqs5.odds)

# Get colors
library(RColorBrewer)
library(plotly)
cols <- brewer.pal(10, "Paired")

# Plot
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/outlierloci_lineplot.png", width=9, height=6, res=300, units="in")

par(
  mar=c(5, 5, 2, 15), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n',
  xpd = NA
)

plot(pop.allele.freqs5.odds2[1,], type = "l", ylim = c(0.40, 1), xlab = "", ylab = "", xaxt = 'n', yaxt = 'n', col = cols[10])
lines(pop.allele.freqs5.odds2[2,], type = "l", col = cols[9], lty = 2)
lines(pop.allele.freqs5.odds2[3,], type = "l", col = cols[8])
lines(pop.allele.freqs5.odds2[4,], type = "l", col = cols[7], lty = 2)
lines(pop.allele.freqs5.odds2[5,], type = "l", col = cols[6])
lines(pop.allele.freqs5.odds2[6,], type = "l", col = cols[5], lty = 2)
lines(pop.allele.freqs5.odds2[7,], type = "l", col = cols[4])
lines(pop.allele.freqs5.odds2[8,], type = "l", col = cols[3], lty = 2)
lines(pop.allele.freqs5.odds2[9,], type = "l", col = cols[2])
lines(pop.allele.freqs5.odds2[10,], type = "l", col = cols[1], lty = 2)

axis(1, at=seq(1,5, by=1), labels=c('A', 'B', 'C', 'D', 'E'), line = 0)
mtext("Adult region (north to south)", side = 1, line = 2.5)
axis(2, at=seq(0.4,1, by=0.1), labels=seq(0.4,1, by= 0.1), line = 0, las = 1)
mtext("Allele frequency", side = 2, line = 3)

legend(5.2,1, legend = row.names(pop.allele.freqs5.odds2), col = rev(cols), lty = rep(1:2, times = 5), cex = 0.8)

dev.off()

########################################################
# Plot allele frequencies of 10 outlier loci in larvae
oto.gen.merged <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/oto.gen.merged151.txt", header = TRUE)

oto.gen.merged.freqs <- oto.gen.merged[, c(1,2, 42,43,16:35)]# Just allele frequencies
oto.gen.merged.freqs.bytime <- split(oto.gen.merged.freqs, oto.gen.merged.freqs$Period)

# Allele counts for each time period, split into ingress sites
larvs.bytime.early <- split(oto.gen.merged.freqs.bytime$Early, oto.gen.merged.freqs.bytime$Early$Location)
larvs.bytime.mid <- split(oto.gen.merged.freqs.bytime$Mid, oto.gen.merged.freqs.bytime$Mid$Location)
larvs.bytime.late <- split(oto.gen.merged.freqs.bytime$Late, oto.gen.merged.freqs.bytime$Late$Location)

############################################
# Allele frequencies of early ingress sites
early.rumfs <- colSums(larvs.bytime.early$RUMFS[5:24], na.rm = TRUE)/(2*colSums(!is.na(larvs.bytime.early$RUMFS[5:24])))
early.nc <- colSums(larvs.bytime.early$NC[5:24], na.rm = TRUE)/(2*colSums(!is.na(larvs.bytime.early$NC[5:24])))

early <- rbind(early.rumfs, early.nc)

# Keep only one allele per locus
odds <- seq(1,20,2) # odd indicies to keep
early.odds <- early[,odds]

# Transpose the data
early.odds2 <- t(early.odds)

# Plot
par(
  mar=c(5, 5, 2, 15), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n',
  xpd = NA
)
plot(early.odds2[1,], type = "l", ylim = c(0.3, 1), xlab = "", ylab = "", xaxt = 'n', yaxt = 'n', col = cols[10])
lines(early.odds2[2,], type = "l", col = cols[9], lty = 2)
lines(early.odds2[3,], type = "l", col = cols[8])
lines(early.odds2[4,], type = "l", col = cols[7], lty = 2)
lines(early.odds2[5,], type = "l", col = cols[6])
lines(early.odds2[6,], type = "l", col = cols[5], lty = 2)
lines(early.odds2[7,], type = "l", col = cols[4])
lines(early.odds2[8,], type = "l", col = cols[3], lty = 2)
lines(early.odds2[9,], type = "l", col = cols[2])
lines(early.odds2[10,], type = "l", col = cols[1], lty = 2)

axis(1, at=seq(1,2, by=1), labels=c('RUMFS', 'NC'), line = 0)
mtext("Ingress site (north to south)", side = 1, line = 2.5)
axis(2, at=seq(0.3,1, by=0.1), labels=seq(0.3,1, by= 0.1), line = 0, las = 1)
mtext("Allele frequency", side = 2, line = 3)

legend(2.2,1, legend = row.names(early.odds2), col = rev(cols), lty = rep(1:2, times = 5), cex = 0.8)

dev.off()

#### Plotting allele frequencies of clusters ####
clusters <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/oto.gen.merged151.clusterdbytime.txt", header = TRUE) # clusters for each time period separately using all elements except Sn

# Separate data by time period (clustering was done for each time period separately)
clusters.bytime <- split(clusters, clusters$Period)
clusters.bytime.early <- split(clusters.bytime$Early, clusters.bytime$Early$cluster)
clusters.bytime.mid <- split(clusters.bytime$Mid, clusters.bytime$Mid$cluster)
clusters.bytime.late <- split(clusters.bytime$Late, clusters.bytime$Late$cluster)

# Calculate allele frequencies
# Early
early1 <- colSums(clusters.bytime.early[[1]][16:35], na.rm = TRUE)/(2*colSums(!is.na(clusters.bytime.early[[1]][16:35])))
early2 <- colSums(clusters.bytime.early[[2]][16:35], na.rm = TRUE)/(2*colSums(!is.na(clusters.bytime.early[[2]][16:35])))
early3 <- colSums(clusters.bytime.early[[3]][16:35], na.rm = TRUE)/(2*colSums(!is.na(clusters.bytime.early[[3]][16:35])))

# Mid
mid1 <- colSums(clusters.bytime.mid[[1]][16:35], na.rm = TRUE)/(2*colSums(!is.na(clusters.bytime.mid[[1]][16:35])))
mid2 <- colSums(clusters.bytime.mid[[2]][16:35], na.rm = TRUE)/(2*colSums(!is.na(clusters.bytime.mid[[2]][16:35])))

# Late
late1 <- colSums(clusters.bytime.late[[1]][16:35], na.rm = TRUE)/(2*colSums(!is.na(clusters.bytime.late[[1]][16:35])))
late2 <- colSums(clusters.bytime.late[[2]][16:35], na.rm = TRUE)/(2*colSums(!is.na(clusters.bytime.late[[2]][16:35])))
late3 <- colSums(clusters.bytime.late[[3]][16:35], na.rm = TRUE)/(2*colSums(!is.na(clusters.bytime.late[[3]][16:35])))

# Keep only one allele per locus. Must match alleles in candidate adult populations
early1.odds <- early1[names(early1) %in% rownames(pop.allele.freqs5.odds2)]
early2.odds <- early2[names(early2) %in% rownames(pop.allele.freqs5.odds2)]
early3.odds <- early3[names(early3) %in% rownames(pop.allele.freqs5.odds2)]

mid1.odds <- mid1[names(mid1) %in% rownames(pop.allele.freqs5.odds2)]
mid2.odds <- mid2[names(mid2) %in% rownames(pop.allele.freqs5.odds2)]

late1.odds <- late1[names(late1) %in% rownames(pop.allele.freqs5.odds2)]
late2.odds <- late2[names(late2) %in% rownames(pop.allele.freqs5.odds2)]
late3.odds <- late3[names(late3) %in% rownames(pop.allele.freqs5.odds2)]

# Plot points for early1 on top of adult allele frequency graph
points(4, early1.odds[1], col = cols[10])
points(4, early1.odds[2], col = cols[9])
points(4, early1.odds[3], col = cols[8])
points(4, early1.odds[4], col = cols[7])
points(4, early1.odds[5], col = cols[6])
points(jitter(4), early1.odds[6], col = cols[5])
points(4, early1.odds[7], col = cols[4])
points(jitter(4), early1.odds[8], col = cols[3])
points(jitter(4), early1.odds[9], col = cols[2])
points(jitter(4), early1.odds[10], col = cols[1])

# Plot points for early2 on top of adult allele frequency graph
points(2, early2.odds[1], col = cols[10])
points(2, early2.odds[2], col = cols[9])
points(jitter(2), early2.odds[3], col = cols[8])
points(2, early2.odds[4], col = cols[7])
points(2, early2.odds[5], col = cols[6])
points(jitter(2), early2.odds[6], col = cols[5])
points(2, early2.odds[7], col = cols[4])
points(2, early2.odds[8], col = cols[3])
points(jitter(2), early2.odds[9], col = cols[2])
points(2, early2.odds[10], col = cols[1])

# Plot points for early3 on top of adult allele frequency graph
points(4, early3.odds[1], col = cols[10])
points(4, early3.odds[2], col = cols[9])
points(4, early3.odds[3], col = cols[8])
points(4, early3.odds[4], col = cols[7])
points(jitter(4), early3.odds[5], col = cols[6])
points(4, early3.odds[6], col = cols[5])
points(4, early3.odds[7], col = cols[4])
points(4, early3.odds[8], col = cols[3])
points(jitter(4), early3.odds[9], col = cols[2])
points(jitter(4), early3.odds[10], col = cols[1])

# Plot points for mid1 on top of adult allele frequency graph
points(4, mid1.odds[1], col = cols[10])
points(4, mid1.odds[2], col = cols[9])
points(4, mid1.odds[3], col = cols[8])
points(4, mid1.odds[4], col = cols[7])
points(4, mid1.odds[5], col = cols[6])
points(jitter(4), mid1.odds[6], col = cols[5])
points(4, mid1.odds[7], col = cols[4])
points(4, mid1.odds[8], col = cols[3])
points(jitter(4), mid1.odds[9], col = cols[2])
points(4, mid1.odds[10], col = cols[1])

# Plot points for mid2 on top of adult allele frequency graph
points(4, mid2.odds[1], col = cols[10])
points(4, mid2.odds[2], col = cols[9])
points(4, mid2.odds[3], col = cols[8])
points(4, mid2.odds[4], col = cols[7])
points(4, mid2.odds[5], col = cols[6])
points(4, mid2.odds[6], col = cols[5])
points(4, mid2.odds[7], col = cols[4])
points(4, mid2.odds[8], col = cols[3])
points(4, mid2.odds[9], col = cols[2])
points(jitter(4), mid2.odds[10], col = cols[1])

# Plot points for late1 on top of adult allele frequency graph
points(4, late1.odds[1], col = cols[10])
points(4, late1.odds[2], col = cols[9])
points(4, late1.odds[3], col = cols[8])
points(4, late1.odds[4], col = cols[7])
points(4, late1.odds[5], col = cols[6])
points(4, late1.odds[6], col = cols[5])
points(4, late1.odds[7], col = cols[4])
points(4, late1.odds[8], col = cols[3])
points(4, late1.odds[9], col = cols[2])
points(4, late1.odds[10], col = cols[1])

# Plot points for late2 on top of adult allele frequency graph
points(4, late2.odds[1], col = cols[10])
points(4, late2.odds[2], col = cols[9])
points(4, late2.odds[3], col = cols[8])
points(4, late2.odds[4], col = cols[7])
points(4, late2.odds[5], col = cols[6])
points(4, late2.odds[6], col = cols[5])
points(4, late2.odds[7], col = cols[4])
points(4, late2.odds[8], col = cols[3])
points(jitter(4), late2.odds[9], col = cols[2])
points(4, late2.odds[10], col = cols[1])

# Plot points for late3 on top of adult allele frequency graph
points(4, late3.odds[1], col = cols[10])
points(4, late3.odds[2], col = cols[9])
points(4, late3.odds[3], col = cols[8])
points(4, late3.odds[4], col = cols[7])
points(4, late3.odds[5], col = cols[6])
points(4, late3.odds[6], col = cols[5])
points(4, late3.odds[7], col = cols[4])
points(jitter(4), late3.odds[8], col = cols[3])
points(4, late3.odds[9], col = cols[2])
points(jitter(4), late3.odds[10], col = cols[1])

#########################################################################################################################################
# Difference between adult and larval allele frequencies. Smallest difference seems to correspond to ML origin population, except for M2
sum(abs(pop.allele.freqs5.odds2[,'Pop2'] - mid2.odds))
sum(abs(pop.allele.freqs5.odds2[,'Pop4'] - mid2.odds))

diffs <- rbind(
colSums(abs(pop.allele.freqs5.odds2 - early1.odds)),
colSums(abs(pop.allele.freqs5.odds2 - early2.odds)),
colSums(abs(pop.allele.freqs5.odds2 - early3.odds)),
colSums(abs(pop.allele.freqs5.odds2 - mid1.odds)),
colSums(abs(pop.allele.freqs5.odds2 - mid2.odds)),
colSums(abs(pop.allele.freqs5.odds2 - late1.odds)),
colSums(abs(pop.allele.freqs5.odds2 - late2.odds)),
colSums(abs(pop.allele.freqs5.odds2 - late3.odds))
)

rownames(diffs) <- c('E1', 'E2', 'E3', 'M1', 'M2', 'L1', 'L2', 'L3')
