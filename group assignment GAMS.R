#### Group assignment of clustered individuals based on allele frequencies from GAMS ###
#### 
#### This script simulates allele frequencies at 10 outlier adult summer flounder loci of 10 populations based on GAMS with distance. The likelihood of each genotype 
#### is calculated using the GAM allele frequency for each of 10 groups. These genotype likelihoods are then multipled across the number  
#### of individuals in a group, which is the depth of the simulated array

# Read in cluster sizes & set a vector of cluster sizes of which to resample
clusters <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/oto.gen.merged151.clusterdbytime.txt", header = TRUE) # clusters for each time period separately using all elements except Sn
clusters.bytime <- split(clusters, clusters$Period)
cluster.sizes <- as.numeric(table(clusters.bytime$Early$cluster)) # this will need to be set for early, middle and late separately

# 10 populations from distance GAMS
obs.likes <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/obs_likelihoods_early_6clusters_10pops.txt', header = TRUE) # when clustering is done for each time period separately; early clusters

##########################################################################
# Read in adult outlier allele frequencies based on 10 populations
# Keep only one allele per locus
# odds <- seq(1,20,2) # odd indicies to keep
# pop.allele.freqs5.odds <- pop.allele.freqs5[,odds]

# pop.allele.freqs10.20alleles <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.freqs.10pops.txt', header = TRUE)
pop.allele.freqs5.odds <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.freqs.10pops.odd.ordered.txt', header = TRUE) # these are the alleles used to calculate genotype likelihoods of clusters

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
    
    for (j in 1:nrow(Pop6.alleles)){
      
      for (i in 1:ncol(Pop6.alleles)){
        if(Pop6.alleles[j,i,h] == 2) {
          Pop6.likelihoods.pop6dist[j,i,h] <- pop.allele.freqs5.odds['Pop6',i]^2
        } else if (Pop6.alleles[j,i,h] == 1) {
          Pop6.likelihoods.pop6dist[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop6',i] * (1-pop.allele.freqs5.odds['Pop6',i]))
        } else if (Pop6.alleles[j,i,h] == 0) {
          Pop6.likelihoods.pop6dist[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop6',i])^2 
        } else {
          Pop6.likelihoods.pop6dist[j,i,h] <- 1
        }
      }
    }
    
    for (j in 1:nrow(Pop7.alleles)){
      
      for (i in 1:ncol(Pop7.alleles)){
        if(Pop7.alleles[j,i,h] == 2) {
          Pop7.likelihoods.pop7dist[j,i,h] <- pop.allele.freqs5.odds['Pop7',i]^2
        } else if (Pop7.alleles[j,i,h] == 1) {
          Pop7.likelihoods.pop7dist[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop7',i] * (1-pop.allele.freqs5.odds['Pop7',i]))
        } else if (Pop7.alleles[j,i,h] == 0) {
          Pop7.likelihoods.pop7dist[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop7',i])^2 
        } else {
          Pop7.likelihoods.pop7dist[j,i,h] <- 1
        }
      }
    }
    
    for (j in 1:nrow(Pop8.alleles)){
      
      for (i in 1:ncol(Pop8.alleles)){
        if(Pop8.alleles[j,i,h] == 2) {
          Pop8.likelihoods.pop8dist[j,i,h] <- pop.allele.freqs5.odds['Pop8',i]^2
        } else if (Pop8.alleles[j,i,h] == 1) {
          Pop8.likelihoods.pop8dist[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop8',i] * (1-pop.allele.freqs5.odds['Pop8',i]))
        } else if (Pop8.alleles[j,i,h] == 0) {
          Pop8.likelihoods.pop8dist[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop8',i])^2 
        } else {
          Pop8.likelihoods.pop8dist[j,i,h] <- 1
        }
      }
    }
    
    for (j in 1:nrow(Pop9.alleles)){
      
      for (i in 1:ncol(Pop9.alleles)){
        if(Pop9.alleles[j,i,h] == 2) {
          Pop9.likelihoods.pop9dist[j,i,h] <- pop.allele.freqs5.odds['Pop9',i]^2
        } else if (Pop9.alleles[j,i,h] == 1) {
          Pop9.likelihoods.pop9dist[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop9',i] * (1-pop.allele.freqs5.odds['Pop9',i]))
        } else if (Pop9.alleles[j,i,h] == 0) {
          Pop9.likelihoods.pop9dist[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop9',i])^2 
        } else {
          Pop9.likelihoods.pop9dist[j,i,h] <- 1
        }
      }
    }
    
    for (j in 1:nrow(Pop10.alleles)){
      
      for (i in 1:ncol(Pop10.alleles)){
        if(Pop10.alleles[j,i,h] == 2) {
          Pop10.likelihoods.pop10dist[j,i,h] <- pop.allele.freqs5.odds['Pop10',i]^2
        } else if (Pop10.alleles[j,i,h] == 1) {
          Pop10.likelihoods.pop10dist[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop10',i] * (1-pop.allele.freqs5.odds['Pop10',i]))
        } else if (Pop10.alleles[j,i,h] == 0) {
          Pop10.likelihoods.pop10dist[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop10',i])^2 
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

