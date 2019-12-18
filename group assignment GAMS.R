#### Group assignment of clustered individuals based on allele frequencies from GAMS ###
#### 
#### This script simulates allele frequencies at 10 outlier adult summer flounder loci of 10 populations based on GAMS with distance. The likelihood of each genotype 
#### is calculated using the GAM allele frequency for each of 10 groups. These genotype likelihoods are then multipled across the number  
#### of individuals in a group, which is the depth of the simulated array

# Read in cluster sizes & set a vector of cluster sizes of which to resample
clusters <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/oto.gen.merged151.clusterdbytime.txt", header = TRUE) # clusters for each time period separately using all elements except Sn
clusters.bytime <- split(clusters, clusters$Period)
cluster.sizes <- as.numeric(table(clusters.bytime$Late$cluster)) # this will need to be set for early, middle and late separately

#### Read in the likelihoods for each cluster ####
# 10 populations from distance GAMS
# obs.likes <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_early_6clusters_10pops.txt', header = TRUE) # when clustering is done for each time period separately; early clusters
# obs.likes <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_middle_2clusters_10pops.txt', header = TRUE)
obs.likes <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_late_3clusters_10pops.txt', header = TRUE)

# Most likely origin population?
most.like <- colnames(obs.likes[apply(obs.likes,1, which.max)]) # pick most likely origin region for each cluster

##########################################################################
# Read in adult outlier allele frequencies based on 10 populations
# Keep only one allele per locus
# odds <- seq(1,20,2) # odd indicies to keep
# pop.allele.freqs5.odds <- pop.allele.freqs5[,odds]

# pop.allele.freqs10.20alleles <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.freqs.10pops.txt', header = TRUE)
# pop.allele.freqs5.odds <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.freqs.10pops.odd.ordered.txt', header = TRUE) # GAM with distance: these are the alleles used to calculate genotype likelihoods of clusters
pop.allele.freqs5.odds <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.freqs.10pops.gambinomial.ordered.txt', header = TRUE) # GAM with distance and binomial errors: these are the alleles used to calculate genotype likelihoods of clusters; 10 x 10, one allele per locus

################################################################################################
# Sample one allele per locus many times from each BayEnv population, with the z dimension being the number of individuals in a 'cluster'
# Set group size (n) for testing one distribution at a time
# cluster.sizes <- 20
# obs.likes <- obs.likes[3,]

dist.likes <- data.frame()

for(n in 1:length(cluster.sizes)){
  
  Pop1.alleles <- array(dim = c(10000,10,paste(cluster.sizes[n])))
  for (i in 1:paste(cluster.sizes[n])){
    Pop1.alleles[,,i] <- sapply(pop.allele.freqs5.odds[1,],function(z){rbinom(10000,2,z)})
  }
  
  # Do this for the rest of the populations too
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
  
  Pop6.alleles <- array(dim = c(10000,10,paste(cluster.sizes[n])))
  for (i in 1:paste(cluster.sizes[n])){
    Pop6.alleles[,,i] <- sapply(pop.allele.freqs5.odds[6,],function(z){rbinom(10000,2,z)})
  }
  
  Pop7.alleles <- array(dim = c(10000,10,paste(cluster.sizes[n])))
  for (i in 1:paste(cluster.sizes[n])){
    Pop7.alleles[,,i] <- sapply(pop.allele.freqs5.odds[7,],function(z){rbinom(10000,2,z)})
  }
  
  Pop8.alleles <- array(dim = c(10000,10,paste(cluster.sizes[n])))
  for (i in 1:paste(cluster.sizes[n])){
    Pop8.alleles[,,i] <- sapply(pop.allele.freqs5.odds[8,],function(z){rbinom(10000,2,z)})
  }
  
  Pop9.alleles <- array(dim = c(10000,10,paste(cluster.sizes[n])))
  for (i in 1:paste(cluster.sizes[n])){
    Pop9.alleles[,,i] <- sapply(pop.allele.freqs5.odds[9,],function(z){rbinom(10000,2,z)})
  }
  
  Pop10.alleles <- array(dim = c(10000,10,paste(cluster.sizes[n])))
  for (i in 1:paste(cluster.sizes[n])){
    Pop10.alleles[,,i] <- sapply(pop.allele.freqs5.odds[10,],function(z){rbinom(10000,2,z)})
  }
  
  # Use adult allele frequencies to calculate likelihoods
  # For loop to loop through each of 10 loci & multiply by the corresponding Bayenv pop adult allele frequency, and then do this for all 10000 offspring. NA's/9's get coded as 1's so they don't make a difference when each row's product is taken
  Pop1.likelihoods.pop1dist <- array(dim = c(10000,10,paste(cluster.sizes[n])))
  Pop2.likelihoods.pop2dist <- array(dim = c(10000,10,paste(cluster.sizes[n])))
  Pop3.likelihoods.pop3dist <- array(dim = c(10000,10,paste(cluster.sizes[n])))
  Pop4.likelihoods.pop4dist <- array(dim = c(10000,10,paste(cluster.sizes[n])))
  Pop5.likelihoods.pop5dist <- array(dim = c(10000,10,paste(cluster.sizes[n])))
  Pop6.likelihoods.pop6dist <- array(dim = c(10000,10,paste(cluster.sizes[n])))
  Pop7.likelihoods.pop7dist <- array(dim = c(10000,10,paste(cluster.sizes[n])))
  Pop8.likelihoods.pop8dist <- array(dim = c(10000,10,paste(cluster.sizes[n])))
  Pop9.likelihoods.pop9dist <- array(dim = c(10000,10,paste(cluster.sizes[n])))
  Pop10.likelihoods.pop10dist <- array(dim = c(10000,10,paste(cluster.sizes[n])))
  
  for (h in 1:paste(cluster.sizes[n])){
    
    for (j in 1:nrow(Pop1.alleles)){
      
      for (i in 1:ncol(Pop1.alleles)){
        if(Pop1.alleles[j,i,h] == 2) {
          Pop1.likelihoods.pop1dist[j,i,h] <- pop.allele.freqs5.odds[1,i]^2
        } else if (Pop1.alleles[j,i,h] == 1) {
          Pop1.likelihoods.pop1dist[j,i,h] <- 2*(pop.allele.freqs5.odds[1,i] * (1-pop.allele.freqs5.odds[1,i]))
        } else if (Pop1.alleles[j,i,h] == 0) {
          Pop1.likelihoods.pop1dist[j,i,h] <- ( 1-pop.allele.freqs5.odds[1,i])^2 
        } else {
          Pop1.likelihoods.pop1dist[j,i,h] <- 1
        }
      }
    }
    
    for (j in 1:nrow(Pop2.alleles)){
      
      for (i in 1:ncol(Pop2.alleles)){
        if(Pop2.alleles[j,i,h] == 2) {
          Pop2.likelihoods.pop2dist[j,i,h] <- pop.allele.freqs5.odds[2,i]^2
        } else if (Pop2.alleles[j,i,h] == 1) {
          Pop2.likelihoods.pop2dist[j,i,h] <- 2*(pop.allele.freqs5.odds[2,i] * (1-pop.allele.freqs5.odds[2,i]))
        } else if (Pop2.alleles[j,i,h] == 0) {
          Pop2.likelihoods.pop2dist[j,i,h] <- ( 1-pop.allele.freqs5.odds[2,i])^2 
        } else {
          Pop2.likelihoods.pop2dist[j,i,h] <- 1
        }
      }
    }
    
    for (j in 1:nrow(Pop3.alleles)){
      
      for (i in 1:ncol(Pop3.alleles)){
        if(Pop3.alleles[j,i,h] == 2) {
          Pop3.likelihoods.pop3dist[j,i,h] <- pop.allele.freqs5.odds[3,i]^2
        } else if (Pop3.alleles[j,i,h] == 1) {
          Pop3.likelihoods.pop3dist[j,i,h] <- 2*(pop.allele.freqs5.odds[3,i] * (1-pop.allele.freqs5.odds[3,i]))
        } else if (Pop3.alleles[j,i,h] == 0) {
          Pop3.likelihoods.pop3dist[j,i,h] <- ( 1-pop.allele.freqs5.odds[3,i])^2 
        } else {
          Pop3.likelihoods.pop3dist[j,i,h] <- 1
        }
      }
    }
    
    for (j in 1:nrow(Pop4.alleles)){
      
      for (i in 1:ncol(Pop4.alleles)){
        if(Pop4.alleles[j,i,h] == 2) {
          Pop4.likelihoods.pop4dist[j,i,h] <- pop.allele.freqs5.odds[4,i]^2
        } else if (Pop4.alleles[j,i,h] == 1) {
          Pop4.likelihoods.pop4dist[j,i,h] <- 2*(pop.allele.freqs5.odds[4,i] * (1-pop.allele.freqs5.odds[4,i]))
        } else if (Pop4.alleles[j,i,h] == 0) {
          Pop4.likelihoods.pop4dist[j,i,h] <- ( 1-pop.allele.freqs5.odds[4,i])^2 
        } else {
          Pop4.likelihoods.pop4dist[j,i,h] <- 1
        }
      }
    }
    
    for (j in 1:nrow(Pop5.alleles)){
      
      for (i in 1:ncol(Pop5.alleles)){
        if(Pop5.alleles[j,i,h] == 2) {
          Pop5.likelihoods.pop5dist[j,i,h] <- pop.allele.freqs5.odds[5,i]^2
        } else if (Pop5.alleles[j,i,h] == 1) {
          Pop5.likelihoods.pop5dist[j,i,h] <- 2*(pop.allele.freqs5.odds[5,i] * (1-pop.allele.freqs5.odds[5,i]))
        } else if (Pop5.alleles[j,i,h] == 0) {
          Pop5.likelihoods.pop5dist[j,i,h] <- ( 1-pop.allele.freqs5.odds[5,i])^2 
        } else {
          Pop5.likelihoods.pop5dist[j,i,h] <- 1
        }
      }
    }
    
    for (j in 1:nrow(Pop6.alleles)){
      
      for (i in 1:ncol(Pop6.alleles)){
        if(Pop6.alleles[j,i,h] == 2) {
          Pop6.likelihoods.pop6dist[j,i,h] <- pop.allele.freqs5.odds[6,i]^2
        } else if (Pop6.alleles[j,i,h] == 1) {
          Pop6.likelihoods.pop6dist[j,i,h] <- 2*(pop.allele.freqs5.odds[6,i] * (1-pop.allele.freqs5.odds[6,i]))
        } else if (Pop6.alleles[j,i,h] == 0) {
          Pop6.likelihoods.pop6dist[j,i,h] <- ( 1-pop.allele.freqs5.odds[6,i])^2 
        } else {
          Pop6.likelihoods.pop6dist[j,i,h] <- 1
        }
      }
    }
    
    for (j in 1:nrow(Pop7.alleles)){
      
      for (i in 1:ncol(Pop7.alleles)){
        if(Pop7.alleles[j,i,h] == 2) {
          Pop7.likelihoods.pop7dist[j,i,h] <- pop.allele.freqs5.odds[7,i]^2
        } else if (Pop7.alleles[j,i,h] == 1) {
          Pop7.likelihoods.pop7dist[j,i,h] <- 2*(pop.allele.freqs5.odds[7,i] * (1-pop.allele.freqs5.odds[7,i]))
        } else if (Pop7.alleles[j,i,h] == 0) {
          Pop7.likelihoods.pop7dist[j,i,h] <- ( 1-pop.allele.freqs5.odds[7,i])^2 
        } else {
          Pop7.likelihoods.pop7dist[j,i,h] <- 1
        }
      }
    }
    
    for (j in 1:nrow(Pop8.alleles)){
      
      for (i in 1:ncol(Pop8.alleles)){
        if(Pop8.alleles[j,i,h] == 2) {
          Pop8.likelihoods.pop8dist[j,i,h] <- pop.allele.freqs5.odds[8,i]^2
        } else if (Pop8.alleles[j,i,h] == 1) {
          Pop8.likelihoods.pop8dist[j,i,h] <- 2*(pop.allele.freqs5.odds[8,i] * (1-pop.allele.freqs5.odds[8,i]))
        } else if (Pop8.alleles[j,i,h] == 0) {
          Pop8.likelihoods.pop8dist[j,i,h] <- ( 1-pop.allele.freqs5.odds[8,i])^2 
        } else {
          Pop8.likelihoods.pop8dist[j,i,h] <- 1
        }
      }
    }
    
    for (j in 1:nrow(Pop9.alleles)){
      
      for (i in 1:ncol(Pop9.alleles)){
        if(Pop9.alleles[j,i,h] == 2) {
          Pop9.likelihoods.pop9dist[j,i,h] <- pop.allele.freqs5.odds[9,i]^2
        } else if (Pop9.alleles[j,i,h] == 1) {
          Pop9.likelihoods.pop9dist[j,i,h] <- 2*(pop.allele.freqs5.odds[9,i] * (1-pop.allele.freqs5.odds[9,i]))
        } else if (Pop9.alleles[j,i,h] == 0) {
          Pop9.likelihoods.pop9dist[j,i,h] <- ( 1-pop.allele.freqs5.odds[9,i])^2 
        } else {
          Pop9.likelihoods.pop9dist[j,i,h] <- 1
        }
      }
    }
    
    for (j in 1:nrow(Pop10.alleles)){
      
      for (i in 1:ncol(Pop10.alleles)){
        if(Pop10.alleles[j,i,h] == 2) {
          Pop10.likelihoods.pop10dist[j,i,h] <- pop.allele.freqs5.odds[10,i]^2
        } else if (Pop10.alleles[j,i,h] == 1) {
          Pop10.likelihoods.pop10dist[j,i,h] <- 2*(pop.allele.freqs5.odds[10,i] * (1-pop.allele.freqs5.odds[10,i]))
        } else if (Pop10.alleles[j,i,h] == 0) {
          Pop10.likelihoods.pop10dist[j,i,h] <- ( 1-pop.allele.freqs5.odds[10,i])^2 
        } else {
          Pop10.likelihoods.pop10dist[j,i,h] <- 1
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
  Pop6.pop6dist <- apply(Pop6.likelihoods.pop6dist, FUN = prod, MARGIN = 1, na.rm = TRUE)
  Pop7.pop7dist <- apply(Pop7.likelihoods.pop7dist, FUN = prod, MARGIN = 1, na.rm = TRUE)
  Pop8.pop8dist <- apply(Pop8.likelihoods.pop8dist, FUN = prod, MARGIN = 1, na.rm = TRUE)
  Pop9.pop9dist <- apply(Pop9.likelihoods.pop9dist, FUN = prod, MARGIN = 1, na.rm = TRUE)
  Pop10.pop10dist <- apply(Pop10.likelihoods.pop10dist, FUN = prod, MARGIN = 1, na.rm = TRUE)
  
  distributions <- list(Pop1.pop1dist, Pop2.pop2dist, Pop3.pop3dist, Pop4.pop4dist, Pop5.pop5dist, Pop6.pop6dist, Pop7.pop7dist, Pop8.pop8dist, Pop9.pop9dist, Pop10.pop10dist)
  
  # Calculate p-values
  for (d in 1:length(distributions)){
    dist.likes[n,d] <- length(which(log10(distributions[[d]]) < obs.likes[n,d]))/length(distributions[[d]])
  }
}

rownames(dist.likes) <- rownames(obs.likes)
colnames(dist.likes) <- c('prob_pop1', 'prob_pop2', 'prob_pop3', 'prob_pop4', 'prob_pop5', 'prob_pop6', 'prob_pop7', 'prob_pop8', 'prob_pop9', 'prob_pop10')


##################################################
# Plot likelihoods across region for each cluster
# 10 GAM-generated groups
obs.likes.early <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_early_6clusters_10pops.txt', header = TRUE) # when clustering is done for each time period separately; early clusters
obs.likes.mid <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_middle_2clusters_10pops.txt', header = TRUE) # middle clusters
obs.likes.late <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_late_3clusters_10pops.txt', header = TRUE) # late clusters


png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/cluster_likelihoods_10GAMSgroups.png", width=12, height=9.5, res=300, units="in")

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

plot(c(1:10), rev(obs.likes.early[1,]), type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'E1')
axis(1, at=seq(1,10, by=1), labels=F)
abline(v = 6)
plot(c(1:10), rev(obs.likes.early[2,]), type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'E2')
axis(1, at=seq(1,10, by=1), labels=F)
abline(v = 5)
plot(c(1:10), rev(obs.likes.early[3,]), type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'E3')
axis(1, at=seq(1,10, by=1), labels=F)
abline(v = 5)
plot(c(1:10), rev(obs.likes.early[4,]), type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'E4')
axis(1, at=seq(1,10, by=1), labels=F)
abline(v = 3)
plot(c(1:10), rev(obs.likes.early[5,]), type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'E5')
axis(1, at=seq(1,10, by=1), labels=F)
abline(v=7)
plot(c(1:10), rev(obs.likes.early[6,]), type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'E6')
axis(1, at=seq(1,10, by=1), labels=F)
abline(v=7)
plot(c(1:10), rev(obs.likes.mid[1,]), type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'M1')
axis(1, at=seq(1,10, by=1), labels=F)
abline(v=5)
plot(c(1:10), rev(obs.likes.mid[2,]), type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'M2')
axis(1, at=seq(1,10, by=1), labels=F)
abline(v=5)
plot(c(1:10), rev(obs.likes.mid[2,]), type = 'l', col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = '')
plot(c(1:10), rev(obs.likes.late[1,]), type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'L1')
axis(1, at=seq(1,10, by=1), labels=c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'), line = 0)
abline(v=6)
plot(c(1:10), rev(obs.likes.late[2,]), type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'L2')
axis(1, at=seq(1,10, by=1), labels=c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'), line = 0)
abline(v=5)
plot(c(1:10), rev(obs.likes.late[3,]), type = 'l', col = 'black', xaxt = 'n', xlab = '', ylab = '', main = 'L3')
axis(1, at=seq(1,10, by=1), labels=c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'), line = 0)
abline(v=6)
mtext(expression('log'[10]*' (genotype likelihood)'), side = 2, line = 0.7, outer = TRUE)

dev.off()

#### Calculate most likely, 2nd most likely and least likely for each cluster. Then plot histograms of these differences ####
# 10 GAM-generated groups. These cluster likelihoods are logged-10 likelihoods. Individual likelihoods are not logged.
obs.likes.early <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_early_6clusters_10pops.txt', header = TRUE) # when clustering is done for each time period separately; early clusters
obs.likes.mid <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_middle_2clusters_10pops.txt', header = TRUE) # middle clusters
obs.likes.late <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_late_3clusters_10pops.txt', header = TRUE) # late clusters

# Early
most.like.early <- colnames(obs.likes.early[apply(obs.likes.early,1, which.max)])
most.like.geno.early <- apply(obs.likes.early,1,max)
least.like.geno.early <- apply(obs.likes.early,1,min)

# Find second most likely geno & pop
n <- 10
second.most.like.early <- vector()
second.most.like.geno.early <- vector()
for (i in 1:nrow(obs.likes.early)){
  second.most.like.geno.early[i] <- sort(obs.likes.early[i,],partial=n-1)[n-1]
  second.most.like.early[i] <- colnames(sort(obs.likes.early[i,],partial=n-1)[n-1])
}

second.most.like.geno.early <- unlist(second.most.like.geno.early)
dif <- 10^(most.like.geno.early) - 10^(second.most.like.geno.early) # these are log rules log(M/N) = log(M) - log(N). Undo the log.

hist(dif, xlab = '', main = '')
mtext("Difference between two most \nlikely genotypes", 1, 3.7)
hist(10^(most.like.geno.early) - 10^(least.like.geno.early), xlab = '', main = "")
mtext("Difference between highest and \nlowest genotype likelihoods", 1, 3.7)

# Middle
most.like.mid <- colnames(obs.likes.mid[apply(obs.likes.mid,1, which.max)])
most.like.geno.mid <- apply(obs.likes.mid,1,max)
least.like.geno.mid <- apply(obs.likes.mid,1,min)

# Find second most likely geno & pop
n <- 10
second.most.like.mid <- vector()
second.most.like.geno.mid <- vector()
for (i in 1:nrow(obs.likes.mid)){
  second.most.like.geno.mid[i] <- sort(obs.likes.mid[i,],partial=n-1)[n-1]
  second.most.like.mid[i] <- colnames(sort(obs.likes.mid[i,],partial=n-1)[n-1])
}

second.most.like.geno.mid <- unlist(second.most.like.geno.mid)
dif <- 10^(most.like.geno.mid) - 10^(second.most.like.geno.mid) # these are log rules log(M/N) = log(M) - log(N). Undo the log.

hist(dif, xlab = '', main = '')
mtext("Difference between two most \nlikely genotypes", 1, 3.7)
hist(10^(most.like.geno.mid) - 10^(least.like.geno.mid), xlab = '', main = "")
mtext("Difference between highest and \nlowest genotype likelihoods", 1, 3.7)

# Late
most.like.late <- colnames(obs.likes.late[apply(obs.likes.late,1, which.max)])
most.like.geno.late <- apply(obs.likes.late,1,max)
least.like.geno.late <- apply(obs.likes.late,1,min)

# Find second most likely geno & pop
n <- 10
second.most.like.late <- vector()
second.most.like.geno.late <- vector()
for (i in 1:nrow(obs.likes.late)){
  second.most.like.geno.late[i] <- sort(obs.likes.late[i,],partial=n-1)[n-1]
  second.most.like.late[i] <- colnames(sort(obs.likes.late[i,],partial=n-1)[n-1])
}

second.most.like.geno.late <- unlist(second.most.like.geno.late)
dif <- 10^(most.like.geno.late) - 10^(second.most.like.geno.late) # these are log rules log(M/N) = log(M) - log(N). Undo the log.

hist(dif, xlab = '', main = '')
mtext("Difference between two most \nlikely genotypes", 1, 3.7)
hist(10^(most.like.geno.late) - 10^(least.like.geno.late), xlab = '', main = "")
mtext("Difference between highest and \nlowest genotype likelihoods", 1, 3.7)

