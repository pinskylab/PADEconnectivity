#### Group assignment of clustered individuals ###
#### 
#### This scripts builds upon the 'population power analysis.R' script where I initially explored the best way to do group assignment with 2 and 5 regions
####
#### This script simulates allele frequencies at 10 outlier adult summer flounder loci of 5 BayEnv populations. The likelihood of each genotype 
#### is calculated using the observed allele frequency of the BayEnv population the alleles are drawn from. These genotype likelihoods are then multipled across the number  
#### of individuals in a group, which is the depth of the simulated array

# Read in cluster sizes & set a vector of cluster sizes of which to resample
# oto.gen.merge6 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/oto.gen.merged151.4clusters.txt", header = TRUE)
# table(oto.gen.merge6$cluster4)
# cluster.sizes <- as.numeric(table(oto.gen.merge6$cluster4))

clusters <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/oto.gen.merged151.clusterdbytime.txt", header = TRUE) # clusters for each time period separately using all elements except Sn
clusters.bytime <- split(clusters, clusters$Period)
cluster.sizes <- as.numeric(table(clusters.bytime$Late$cluster)) # this will need to be set for early, middle and late separately

# Read in the likelihoods for each cluster
# obs.likes <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/bayenv_likelihoods_4clusters.txt") # when clustering is done across all time periods

# obs.likes <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_early_6clusters.txt', header = TRUE) # when clustering is done for each time period separately; early clusters
# obs.likes <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_middle_2clusters.txt', header = TRUE) # middle clusters
obs.likes <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_late_3clusters.txt', header = TRUE) # late clusters

# Plot likelihoods across region for each cluster
obs.likes.early <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_early_6clusters.txt', header = TRUE) # when clustering is done for each time period separately; early clusters
obs.likes.mid <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_middle_2clusters.txt', header = TRUE) # middle clusters
obs.likes.late <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_late_3clusters.txt', header = TRUE) # late clusters

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/cluster_likelihoods.png", width=8, height=9.5, res=300, units="in")

par(
  mar=c(2, 2, 1.2, 1), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n',
  par(mfrow = c(4,3)),
  oma = c(0,2,0,0)
)

plot(c(1:5), obs.likes.early[1,], type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'E1')
axis(1, at=seq(1,5, by=1), labels=F)
plot(c(1:5), obs.likes.early[2,], type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'E2')
axis(1, at=seq(1,5, by=1), labels=F)
plot(c(1:5), obs.likes.early[3,], type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'E3')
axis(1, at=seq(1,5, by=1), labels=F)
plot(c(1:5), obs.likes.early[4,], type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'E4')
axis(1, at=seq(1,5, by=1), labels=F)
plot(c(1:5), obs.likes.early[5,], type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'E5')
axis(1, at=seq(1,5, by=1), labels=F)
plot(c(1:5), obs.likes.early[6,], type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'E6')
axis(1, at=seq(1,5, by=1), labels=F)
plot(c(1:5), obs.likes.mid[1,], type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'M1')
axis(1, at=seq(1,5, by=1), labels=F)
plot(c(1:5), obs.likes.mid[2,], type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'M2')
axis(1, at=seq(1,5, by=1), labels=F)
plot(c(1:5), obs.likes.mid[2,], type = 'l', col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = '')
plot(c(1:5), obs.likes.late[1,], type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'L1')
axis(1, at=seq(1,5, by=1), labels=c('A', 'B', 'C', 'D', 'E'), line = 0)
plot(c(1:5), obs.likes.late[2,], type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'L2')
axis(1, at=seq(1,5, by=1), labels=c('A', 'B', 'C', 'D', 'E'), line = 0)
plot(c(1:5), obs.likes.late[3,], type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'L3')
axis(1, at=seq(1,5, by=1), labels=c('A', 'B', 'C', 'D', 'E'), line = 0)
mtext(expression('log'[10]*' (genotype likelihood)'), side = 2, line = 0.7, outer = TRUE)

dev.off()

# Read in adult outlier allele frequencies
pop.allele.freqs5 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.freqs.5pops.txt')

# Keep only one allele per locus
odds <- seq(1,20,2) # odd indicies to keep
pop.allele.freqs5.odds <- pop.allele.freqs5[,odds]

# Sample one allele per locus many times from each BayEnv population, with the z dimension being the number of individuals in a 'cluster'
# Set group size (n) for testing one distribution at a time
cluster.sizes <- 20
obs.likes <- obs.likes[3,]

dist.likes <- data.frame()

for(n in 1:length(cluster.sizes)){

Pop1.alleles <- array(dim = c(10000,10,paste(cluster.sizes[n])))
for (i in 1:paste(cluster.sizes[n])){
  Pop1.alleles[,,i] <- sapply(pop.allele.freqs5.odds[1,],function(z){rbinom(10000,2,z)})
}

# Do this for the rest of the BayEnv populations too
Pop2.alleles <- array(dim = c(10000,10,paste(cluster.sizes[n])))
for (i in 1:paste(cluster.sizes[n])){
  Pop2.alleles[,,i] <- sapply(pop.allele.freqs5.odds[2,],function(z){rbinom(10000,2,z)})
}

Pop3.alleles <- array(dim = c(10000,10,paste(cluster.sizes[n])))
for (i in 1:paste(cluster.sizes[n])){
  Pop3.alleles[,,i] <- sapply(pop.allele.freqs5.odds[3,],function(z){rbinom(10000,2,z)})
}

Pop4.alleles <- array(dim = c(10000,10,paste(cluster.sizes[n])))
for (i in 1:paste(cluster.sizes[n])){
  Pop4.alleles[,,i] <- sapply(pop.allele.freqs5.odds[4,],function(z){rbinom(10000,2,z)})
}

Pop5.alleles <- array(dim = c(10000,10,paste(cluster.sizes[n])))
for (i in 1:paste(cluster.sizes[n])){
  Pop5.alleles[,,i] <- sapply(pop.allele.freqs5.odds[5,],function(z){rbinom(10000,2,z)})
}

# Use adult allele frequencies to calculate likelihoods
# For loop to loop through each of 10 loci & multiply by the corresponding Bayenv pop adult allele frequency, and then do this for all 10000 offspring. NA's/9's get coded as 1's so they don't make a difference when each row's product is taken
Pop1.likelihoods.pop1dist <- array(dim = c(10000,10,paste(cluster.sizes[n])))
Pop2.likelihoods.pop2dist <- array(dim = c(10000,10,paste(cluster.sizes[n])))
Pop3.likelihoods.pop3dist <- array(dim = c(10000,10,paste(cluster.sizes[n])))
Pop4.likelihoods.pop4dist <- array(dim = c(10000,10,paste(cluster.sizes[n])))
Pop5.likelihoods.pop5dist <- array(dim = c(10000,10,paste(cluster.sizes[n])))

for (h in 1:paste(cluster.sizes[n])){
  
  for (j in 1:nrow(Pop1.alleles)){
    
    for (i in 1:ncol(Pop1.alleles)){
      if(Pop1.alleles[j,i,h] == 2) {
        Pop1.likelihoods.pop1dist[j,i,h] <- pop.allele.freqs5.odds['Pop1',i]^2
      } else if (Pop1.alleles[j,i,h] == 1) {
        Pop1.likelihoods.pop1dist[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop1',i] * (1-pop.allele.freqs5.odds['Pop1',i]))
      } else if (Pop1.alleles[j,i,h] == 0) {
        Pop1.likelihoods.pop1dist[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop1',i])^2 
      } else {
        Pop1.likelihoods.pop1dist[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop2.alleles)){
    
    for (i in 1:ncol(Pop2.alleles)){
      if(Pop2.alleles[j,i,h] == 2) {
        Pop2.likelihoods.pop2dist[j,i,h] <- pop.allele.freqs5.odds['Pop2',i]^2
      } else if (Pop2.alleles[j,i,h] == 1) {
        Pop2.likelihoods.pop2dist[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop2',i] * (1-pop.allele.freqs5.odds['Pop2',i]))
      } else if (Pop2.alleles[j,i,h] == 0) {
        Pop2.likelihoods.pop2dist[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop2',i])^2 
      } else {
        Pop2.likelihoods.pop2dist[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop3.alleles)){
  
    for (i in 1:ncol(Pop3.alleles)){
      if(Pop3.alleles[j,i,h] == 2) {
        Pop3.likelihoods.pop3dist[j,i,h] <- pop.allele.freqs5.odds['Pop3',i]^2
      } else if (Pop3.alleles[j,i,h] == 1) {
        Pop3.likelihoods.pop3dist[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop3',i] * (1-pop.allele.freqs5.odds['Pop3',i]))
      } else if (Pop3.alleles[j,i,h] == 0) {
        Pop3.likelihoods.pop3dist[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop3',i])^2 
      } else {
        Pop3.likelihoods.pop3dist[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop4.alleles)){
    
    for (i in 1:ncol(Pop4.alleles)){
      if(Pop4.alleles[j,i,h] == 2) {
        Pop4.likelihoods.pop4dist[j,i,h] <- pop.allele.freqs5.odds['Pop4',i]^2
      } else if (Pop4.alleles[j,i,h] == 1) {
        Pop4.likelihoods.pop4dist[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop4',i] * (1-pop.allele.freqs5.odds['Pop4',i]))
      } else if (Pop4.alleles[j,i,h] == 0) {
        Pop4.likelihoods.pop4dist[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop4',i])^2 
      } else {
        Pop4.likelihoods.pop4dist[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop5.alleles)){
    
    for (i in 1:ncol(Pop5.alleles)){
      if(Pop5.alleles[j,i,h] == 2) {
        Pop5.likelihoods.pop5dist[j,i,h] <- pop.allele.freqs5.odds['Pop5',i]^2
      } else if (Pop5.alleles[j,i,h] == 1) {
        Pop5.likelihoods.pop5dist[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop5',i] * (1-pop.allele.freqs5.odds['Pop5',i]))
      } else if (Pop5.alleles[j,i,h] == 0) {
        Pop5.likelihoods.pop5dist[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop5',i])^2 
      } else {
        Pop5.likelihoods.pop5dist[j,i,h] <- 1
      }
    }
  }
}

# Multiply everything together
Pop1.pop1dist <- apply(Pop1.likelihoods.pop1dist, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop2.pop2dist <- apply(Pop2.likelihoods.pop2dist, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop3.pop3dist <- apply(Pop3.likelihoods.pop3dist, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop4.pop4dist <- apply(Pop4.likelihoods.pop4dist, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop5.pop5dist <- apply(Pop5.likelihoods.pop5dist, FUN = prod, MARGIN = 1, na.rm = TRUE)

distributions <- list(Pop1.pop1dist, Pop2.pop2dist, Pop3.pop3dist, Pop4.pop4dist, Pop5.pop5dist)

# Calculate p-values
for (d in 1:length(distributions)){
  dist.likes[n,d] <- length(which(log10(distributions[[d]]) < obs.likes[n,d]))/length(distributions[[d]])
}
}

rownames(dist.likes) <- rownames(obs.likes)
colnames(dist.likes) <- c('prob_bayenv1', 'prob_bayenv2', 'prob_bayenv3', 'prob_bayenv4', 'prob_bayenv5')


#############


hist(log10(Pop1.pop1dist), col = rgb(1,0,0,0.5), xlim = c(-35,0), ylim = c(0, 2500), xlab = "log10(genotype likelihood)", main = "")
hist(log10(Pop2.pop2dist), add = TRUE, col = rgb(0,0,1,0.5))
hist(log10(Pop3.pop3dist), add = TRUE, col = rgb(1,0,1,0.5))
hist(log10(Pop4.pop4dist), add = TRUE, col = rgb(1,1,0.5,0.5))
hist(log10(Pop5.pop5dist), add = TRUE, col = rgb(0,1,0,0.5))
legend("topright", c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5"), pch = 22, col = "black", pt.bg  = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5), rgb(1,0,1,0.5), rgb(1,1,0.5,0.5), rgb(0,1,0,0.5)))

#### Calculate p-values ####
length(which(log10(Pop1.pop1dist) < obs.likes[2,1]))/length(Pop1.pop1dist) #0.9364, this is the probability that this cluster comes from BayEnv pop1
length(which(log10(Pop2.pop2dist) < obs.likes[2,2]))/length(Pop2.pop2dist) #0.7807, this is the probability that this cluster comes from BayEnv pop2, etc.
length(which(log10(Pop3.pop3dist) < obs.likes[2,3]))/length(Pop3.pop3dist) #0.8173
length(which(log10(Pop4.pop4dist) < obs.likes[2,4]))/length(Pop4.pop4dist) #0.2919
length(which(log10(Pop5.pop5dist) < obs.likes[2,5]))/length(Pop5.pop5dist) #0.0535

# For look to calculate the probability a cluster belogs to a particular population (need to change which cluster I'm referring to in obs.likes)
cluster.sizes
dist.likes <- vector()

for (d in 1:length(distributions)){
  dist.likes[d] <- length(which(log10(distributions[[d]]) < obs.likes[n,d]))/length(distributions[[d]])
}

# Plot histogram with observed value and calculated p-value
# Row index of observed cluster size needs to be correct
par(mfrow = c(3,2))
hist(log10(Pop1.pop1dist), col = rgb(1,0,0,0.5), xlab = "log10(genotype likelihood)", main = "Pop1")
abline(v=obs.likes[n,1])
text(locator(1), paste0("p =", length(which(log10(Pop1.pop1dist) < obs.likes[n,1]))/length(Pop1.pop1dist)))

hist(log10(Pop2.pop2dist), col = rgb(0,0,1,0.5), xlab = "log10(genotype likelihood)", main = "Pop2")
abline(v=obs.likes[n,2])
text(locator(1), paste0("p =", length(which(log10(Pop2.pop2dist) < obs.likes[n,2]))/length(Pop2.pop2dist)))

hist(log10(Pop3.pop3dist), col = rgb(1,0,1,0.5), xlab = "log10(genotype likelihood)", main = "Pop3")
abline(v=obs.likes[n,3])
text(locator(1), paste0("p =", length(which(log10(Pop3.pop3dist) < obs.likes[n,3]))/length(Pop3.pop3dist)))

hist(log10(Pop4.pop4dist), col = rgb(1,1,0.5,0.5), xlab = "log10(genotype likelihood)", main = "Pop4")
abline(v=obs.likes[n,4])
text(locator(1), paste0("p =", length(which(log10(Pop4.pop4dist) < obs.likes[n,4]))/length(Pop4.pop4dist)))

hist(log10(Pop5.pop5dist), col = rgb(0,1,0,0.5), xlab = "log10(genotype likelihood)", main = "Pop5")
abline(v=obs.likes[n,5])
text(locator(1), paste0("p =", length(which(log10(Pop5.pop5dist) < obs.likes[n,5]))/length(Pop5.pop5dist)))

mtext(paste0("Cluster size = ", cluster.sizes[n]), outer = TRUE, line = -2)
