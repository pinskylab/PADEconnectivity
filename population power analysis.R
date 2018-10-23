#### Population power analysis using different sized populations ###
#### 
#### This script simulates allele frequencies at 10 outlier adult summer flounder loci north and south of Cape Hatteras. The likelihood of each genotype 
#### is calculated using the observed northern and southern allele frequencies. These genotype likelihoods are then multipled across the number of individuals 
#### in a population, which is the depth of the simulated array

# Read in adult outlier allele frequencies
pop.allele.freqs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/regional.adult.outs.freqs10.txt')

# Keep only one allele per locus
odds <- seq(1,20,2) # odd indicies to keep
pop.allele.freqs.odds <- pop.allele.freqs[,odds]
pop.allele.freqs.odds.north <- pop.allele.freqs.odds[1,] # north allele frequency
pop.allele.freqs.odds.south <- pop.allele.freqs.odds[2,] # south allele frequency

# Sample one allele per locus many times, with the z dimension being the number of individuals in a population
# Population = 5
n.alleles.5 <- array(dim = c(1000,10,5))
for (i in 1:5){
  n.alleles.5[,,i] <- sapply(pop.allele.freqs.odds.north,function(z){rbinom(1000,2,z)})
}

# Do this for the southern distribution too
s.alleles.5 <- array(dim = c(1000,10,5))
for (i in 1:5){
  s.alleles.5[,,i] <- sapply(pop.allele.freqs.odds.south,function(z){rbinom(1000,2,z)})
}

# Use adult allele frequencies to calculate likelihoods
# For loop to loop through each of 10 loci & multiply by the adult allele frequency, and then do this for all 1000 offspring. NA's/9's get coded as 1's so they don't make a difference when each row's product is taken
north.likelihoods.ndist.5 <- array(dim = c(1000,10,5))
south.likelihoods.ndist.5 <- array(dim = c(1000,10,5))

north.likelihoods.sdist.5 <- array(dim = c(1000,10,5))
south.likelihoods.sdist.5 <- array(dim = c(1000,10,5))

for (h in 1:5){
  
  for (j in 1:nrow(n.alleles.5)){
    
    for (i in 1:ncol(n.alleles.5)){
      if(n.alleles.5[j,i,h] == 2) {
        north.likelihoods.ndist.5[j,i,h] <- pop.allele.freqs.odds[1,i]^2
      } else if (n.alleles.5[j,i,h] == 1) {
        north.likelihoods.ndist.5[j,i,h] <- 2*(pop.allele.freqs.odds[1,i] * (1-pop.allele.freqs.odds[1,i]))
      } else if (n.alleles.5[j,i,h] == 0) {
        north.likelihoods.ndist.5[j,i,h] <- ( 1-pop.allele.freqs.odds[1,i])^2 
      } else {
        north.likelihoods.ndist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(n.alleles.5)){
      if(n.alleles.5[j,i,h] == 2){
        south.likelihoods.ndist.5[j,i,h] <- pop.allele.freqs.odds[2,i]^2
      } else if (n.alleles.5[j,i,h] == 1) {
        south.likelihoods.ndist.5[j,i,h] <- 2*(pop.allele.freqs.odds[2,i] * (1-pop.allele.freqs.odds[2,i]))
      } else  if (n.alleles.5[j,i,h] == 0) {
        south.likelihoods.ndist.5[j,i,h] <- (1-pop.allele.freqs.odds[2,i])^2
      } else {
        south.likelihoods.ndist.5[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(s.alleles.5)){
    
    for (i in 1:ncol(s.alleles.5)){
      if(s.alleles.5[j,i,h] == 2) {
        north.likelihoods.sdist.5[j,i,h] <- pop.allele.freqs.odds[1,i]^2
      } else if (s.alleles.5[j,i,h] == 1) {
        north.likelihoods.sdist.5[j,i,h] <- 2*(pop.allele.freqs.odds[1,i] * (1-pop.allele.freqs.odds[1,i]))
      } else if (s.alleles.5[j,i,h] == 0) {
        north.likelihoods.sdist.5[j,i,h] <- (1-pop.allele.freqs.odds[1,i])^2 
      } else {
        north.likelihoods.sdist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(s.alleles.5)){
      if(s.alleles.5[j,i,h] == 2){
        south.likelihoods.sdist.5[j,i,h] <- pop.allele.freqs.odds[2,i]^2
      } else if (s.alleles.5[j,i,h] == 1) {
        south.likelihoods.sdist.5[j,i,h] <- 2*(pop.allele.freqs.odds[2,i] * (1-pop.allele.freqs.odds[2,i]))
      } else  if (s.alleles.5[j,i,h] == 0) {
        south.likelihoods.sdist.5[j,i,h] <- (1-pop.allele.freqs.odds[2,i])^2
      } else {
        south.likelihoods.sdist.5[j,i,h] <- 1
      }
    }
  }
}

# Multiply everything together
north5.ndist <- apply(north.likelihoods.ndist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
south5.ndist <- apply(south.likelihoods.ndist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)

north5.sdist <- apply(north.likelihoods.sdist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
south5.sdist <- apply(south.likelihoods.sdist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)

# Create ratio & plot
hist(log10(north5.ndist/south5.ndist), xlab = "log10(north likelihood/south likelihood)", main = "Population = 5", col = rgb(1,0,0,0.5), xlim = c(-10,10), ylim = c(0,325))
hist(log10(north5.sdist/south5.sdist), col = rgb(0,0,1,0.5), add = TRUE)
legend("topright", c("North", "South"), col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch = 15)

###############################################################################################################
# Sample one allele per locus many times, with the z dimension being the number of individuals in a population
# Population = 10
n.alleles.10 <- array(dim = c(1000,10,10))
for (i in 1:10){
  n.alleles.10[,,i] <- sapply(pop.allele.freqs.odds.north,function(z){rbinom(1000,2,z)})
}

# Do this for the southern distribution too
s.alleles.10 <- array(dim = c(1000,10,10))
for (i in 1:10){
  s.alleles.10[,,i] <- sapply(pop.allele.freqs.odds.south,function(z){rbinom(1000,2,z)})
}

# Use adult allele frequencies to calculate likelihoods
# For loop to loop through each of 10 loci & multiply by the adult allele frequency, and then do this for all 1000 offspring. NA's/9's get coded as 1's so they don't make a difference when each row's product is taken
north.likelihoods.ndist.10 <- array(dim = c(1000,10,10))
south.likelihoods.ndist.10 <- array(dim = c(1000,10,10))

north.likelihoods.sdist.10 <- array(dim = c(1000,10,10))
south.likelihoods.sdist.10 <- array(dim = c(1000,10,10))

for (h in 1:10){

for (j in 1:nrow(n.alleles.10)){
  
  for (i in 1:ncol(n.alleles.10)){
    if(n.alleles.10[j,i,h] == 2) {
      north.likelihoods.ndist.10[j,i,h] <- pop.allele.freqs.odds[1,i]^2
    } else if (n.alleles.10[j,i,h] == 1) {
      north.likelihoods.ndist.10[j,i,h] <- 2*(pop.allele.freqs.odds[1,i] * (1-pop.allele.freqs.odds[1,i]))
    } else if (n.alleles.10[j,i,h] == 0) {
      north.likelihoods.ndist.10[j,i,h] <- ( 1-pop.allele.freqs.odds[1,i])^2 
    } else {
      north.likelihoods.ndist.10[j,i,h] <- 1
    }
  }
  
  for (i in 1:ncol(n.alleles.10)){
    if(n.alleles.10[j,i,h] == 2){
      south.likelihoods.ndist.10[j,i,h] <- pop.allele.freqs.odds[2,i]^2
    } else if (n.alleles.10[j,i,h] == 1) {
      south.likelihoods.ndist.10[j,i,h] <- 2*(pop.allele.freqs.odds[2,i] * (1-pop.allele.freqs.odds[2,i]))
    } else  if (n.alleles.10[j,i,h] == 0) {
      south.likelihoods.ndist.10[j,i,h] <- (1-pop.allele.freqs.odds[2,i])^2
    } else {
      south.likelihoods.ndist.10[j,i,h] <- 1
    }
  }
}
  
for (j in 1:nrow(s.alleles.10)){
  
  for (i in 1:ncol(s.alleles.10)){
    if(s.alleles.10[j,i,h] == 2) {
      north.likelihoods.sdist.10[j,i,h] <- pop.allele.freqs.odds[1,i]^2
    } else if (s.alleles.10[j,i,h] == 1) {
      north.likelihoods.sdist.10[j,i,h] <- 2*(pop.allele.freqs.odds[1,i] * (1-pop.allele.freqs.odds[1,i]))
    } else if (s.alleles.10[j,i,h] == 0) {
      north.likelihoods.sdist.10[j,i,h] <- (1-pop.allele.freqs.odds[1,i])^2 
    } else {
      north.likelihoods.sdist.10[j,i,h] <- 1
    }
  }
  
  for (i in 1:ncol(s.alleles.10)){
    if(s.alleles.10[j,i,h] == 2){
      south.likelihoods.sdist.10[j,i,h] <- pop.allele.freqs.odds[2,i]^2
    } else if (s.alleles.10[j,i,h] == 1) {
      south.likelihoods.sdist.10[j,i,h] <- 2*(pop.allele.freqs.odds[2,i] * (1-pop.allele.freqs.odds[2,i]))
    } else  if (s.alleles.10[j,i,h] == 0) {
      south.likelihoods.sdist.10[j,i,h] <- (1-pop.allele.freqs.odds[2,i])^2
    } else {
      south.likelihoods.sdist.10[j,i,h] <- 1
    }
  }
}
}

# Multiply everything together
north10.ndist <- apply(north.likelihoods.ndist.10, FUN = prod, MARGIN = 1, na.rm = TRUE)
south10.ndist <- apply(south.likelihoods.ndist.10, FUN = prod, MARGIN = 1, na.rm = TRUE)

north10.sdist <- apply(north.likelihoods.sdist.10, FUN = prod, MARGIN = 1, na.rm = TRUE)
south10.sdist <- apply(south.likelihoods.sdist.10, FUN = prod, MARGIN = 1, na.rm = TRUE)

# Create ratio & plot
hist(log10(north10.ndist/south10.ndist), xlab = "log10(north likelihood/south likelihood)", main = "Population = 10", col = rgb(1,0,0,0.5), xlim = c(-12,12), ylim = c(0,250))
hist(log10(north10.sdist/south10.sdist), col = rgb(0,0,1,0.5), add = TRUE)
legend("topright", c("North", "South"), col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch = 15)

# Power test
# log10(north.vector.ndist/south.vector.ndist)[order(log10(north.vector.ndist/south.vector.ndist))][50] #5% of 1000 is 50
# length(which(log10(north.vector.sdist/south.vector.sdist) < -0.5864056))

###############################################################################################################
# Sample one allele per locus many times, with the z dimension being the number of individuals in a population
# Population = 30
n.alleles.30 <- array(dim = c(1000,10,30))
for (i in 1:30){
  n.alleles.30[,,i] <- sapply(pop.allele.freqs.odds.north,function(z){rbinom(1000,2,z)})
}

# Do this for the southern distribution too
s.alleles.30 <- array(dim = c(1000,10,30))
for (i in 1:30){
  s.alleles.30[,,i] <- sapply(pop.allele.freqs.odds.south,function(z){rbinom(1000,2,z)})
}

# Use adult allele frequencies to calculate likelihoods
# For loop to loop through each of 10 loci & multiply by the adult allele frequency, and then do this for all 1000 offspring. NA's/9's get coded as 1's so they don't make a difference when each row's product is taken
north.likelihoods.ndist.30 <- array(dim = c(1000,10,30))
south.likelihoods.ndist.30 <- array(dim = c(1000,10,30))

north.likelihoods.sdist.30 <- array(dim = c(1000,10,30))
south.likelihoods.sdist.30 <- array(dim = c(1000,10,30))

for (h in 1:30){
  
  for (j in 1:nrow(n.alleles.30)){
    
    for (i in 1:ncol(n.alleles.30)){
      if(n.alleles.30[j,i,h] == 2) {
        north.likelihoods.ndist.30[j,i,h] <- pop.allele.freqs.odds[1,i]^2
      } else if (n.alleles.30[j,i,h] == 1) {
        north.likelihoods.ndist.30[j,i,h] <- 2*(pop.allele.freqs.odds[1,i] * (1-pop.allele.freqs.odds[1,i]))
      } else if (n.alleles.30[j,i,h] == 0) {
        north.likelihoods.ndist.30[j,i,h] <- ( 1-pop.allele.freqs.odds[1,i])^2 
      } else {
        north.likelihoods.ndist.30[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(n.alleles.30)){
      if(n.alleles.30[j,i,h] == 2){
        south.likelihoods.ndist.30[j,i,h] <- pop.allele.freqs.odds[2,i]^2
      } else if (n.alleles.30[j,i,h] == 1) {
        south.likelihoods.ndist.30[j,i,h] <- 2*(pop.allele.freqs.odds[2,i] * (1-pop.allele.freqs.odds[2,i]))
      } else  if (n.alleles.30[j,i,h] == 0) {
        south.likelihoods.ndist.30[j,i,h] <- (1-pop.allele.freqs.odds[2,i])^2
      } else {
        south.likelihoods.ndist.30[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(s.alleles.30)){
    
    for (i in 1:ncol(s.alleles.30)){
      if(s.alleles.30[j,i,h] == 2) {
        north.likelihoods.sdist.30[j,i,h] <- pop.allele.freqs.odds[1,i]^2
      } else if (s.alleles.30[j,i,h] == 1) {
        north.likelihoods.sdist.30[j,i,h] <- 2*(pop.allele.freqs.odds[1,i] * (1-pop.allele.freqs.odds[1,i]))
      } else if (s.alleles.30[j,i,h] == 0) {
        north.likelihoods.sdist.30[j,i,h] <- (1-pop.allele.freqs.odds[1,i])^2 
      } else {
        north.likelihoods.sdist.30[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(s.alleles.30)){
      if(s.alleles.30[j,i,h] == 2){
        south.likelihoods.sdist.30[j,i,h] <- pop.allele.freqs.odds[2,i]^2
      } else if (s.alleles.30[j,i,h] == 1) {
        south.likelihoods.sdist.30[j,i,h] <- 2*(pop.allele.freqs.odds[2,i] * (1-pop.allele.freqs.odds[2,i]))
      } else  if (s.alleles.30[j,i,h] == 0) {
        south.likelihoods.sdist.30[j,i,h] <- (1-pop.allele.freqs.odds[2,i])^2
      } else {
        south.likelihoods.sdist.30[j,i,h] <- 1
      }
    }
  }
}

# Multiply everything together
north30.ndist <- apply(north.likelihoods.ndist.30, FUN = prod, MARGIN = 1, na.rm = TRUE)
south30.ndist <- apply(south.likelihoods.ndist.30, FUN = prod, MARGIN = 1, na.rm = TRUE)

north30.sdist <- apply(north.likelihoods.sdist.30, FUN = prod, MARGIN = 1, na.rm = TRUE)
south30.sdist <- apply(south.likelihoods.sdist.30, FUN = prod, MARGIN = 1, na.rm = TRUE)

# Create ratio & plot
hist(log10(north30.ndist/south30.ndist), xlab = "log10(north likelihood/south likelihood)", main = "Population = 30", col = rgb(1,0,0,0.5), xlim = c(-22,22), ylim = c(0,250))
hist(log10(north30.sdist/south30.sdist), col = rgb(0,0,1,0.5), add = TRUE)
legend("topright", c("North", "South"), col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch = 15)

