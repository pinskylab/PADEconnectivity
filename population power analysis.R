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

# Power test
log10(north5.ndist/south5.ndist)[order(log10(north5.ndist/south5.ndist))][50] #5% of 1000 is 50
p5 <- length(which(log10(north5.sdist/south5.sdist) < -0.3087461))/1000

###############################################################################################################
# Sample one allele per locus many times, with the z dimension being the number of individuals in a population
# Population = 8
n.alleles.8 <- array(dim = c(1000,10,8))
for (i in 1:8){
  n.alleles.8[,,i] <- sapply(pop.allele.freqs.odds.north,function(z){rbinom(1000,2,z)})
}

# Do this for the southern distribution too
s.alleles.8 <- array(dim = c(1000,10,8))
for (i in 1:8){
  s.alleles.8[,,i] <- sapply(pop.allele.freqs.odds.south,function(z){rbinom(1000,2,z)})
}

# Use adult allele frequencies to calculate likelihoods
# For loop to loop through each of 10 loci & multiply by the adult allele frequency, and then do this for all 1000 offspring. NA's/9's get coded as 1's so they don't make a difference when each row's product is taken
north.likelihoods.ndist.8 <- array(dim = c(1000,10,8))
south.likelihoods.ndist.8 <- array(dim = c(1000,10,8))

north.likelihoods.sdist.8 <- array(dim = c(1000,10,8))
south.likelihoods.sdist.8 <- array(dim = c(1000,10,8))

for (h in 1:8){
  
  for (j in 1:nrow(n.alleles.8)){
    
    for (i in 1:ncol(n.alleles.8)){
      if(n.alleles.8[j,i,h] == 2) {
        north.likelihoods.ndist.8[j,i,h] <- pop.allele.freqs.odds[1,i]^2
      } else if (n.alleles.8[j,i,h] == 1) {
        north.likelihoods.ndist.8[j,i,h] <- 2*(pop.allele.freqs.odds[1,i] * (1-pop.allele.freqs.odds[1,i]))
      } else if (n.alleles.8[j,i,h] == 0) {
        north.likelihoods.ndist.8[j,i,h] <- ( 1-pop.allele.freqs.odds[1,i])^2 
      } else {
        north.likelihoods.ndist.8[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(n.alleles.8)){
      if(n.alleles.8[j,i,h] == 2){
        south.likelihoods.ndist.8[j,i,h] <- pop.allele.freqs.odds[2,i]^2
      } else if (n.alleles.8[j,i,h] == 1) {
        south.likelihoods.ndist.8[j,i,h] <- 2*(pop.allele.freqs.odds[2,i] * (1-pop.allele.freqs.odds[2,i]))
      } else  if (n.alleles.8[j,i,h] == 0) {
        south.likelihoods.ndist.8[j,i,h] <- (1-pop.allele.freqs.odds[2,i])^2
      } else {
        south.likelihoods.ndist.8[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(s.alleles.8)){
    
    for (i in 1:ncol(s.alleles.8)){
      if(s.alleles.8[j,i,h] == 2) {
        north.likelihoods.sdist.8[j,i,h] <- pop.allele.freqs.odds[1,i]^2
      } else if (s.alleles.8[j,i,h] == 1) {
        north.likelihoods.sdist.8[j,i,h] <- 2*(pop.allele.freqs.odds[1,i] * (1-pop.allele.freqs.odds[1,i]))
      } else if (s.alleles.8[j,i,h] == 0) {
        north.likelihoods.sdist.8[j,i,h] <- (1-pop.allele.freqs.odds[1,i])^2 
      } else {
        north.likelihoods.sdist.8[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(s.alleles.8)){
      if(s.alleles.8[j,i,h] == 2){
        south.likelihoods.sdist.8[j,i,h] <- pop.allele.freqs.odds[2,i]^2
      } else if (s.alleles.8[j,i,h] == 1) {
        south.likelihoods.sdist.8[j,i,h] <- 2*(pop.allele.freqs.odds[2,i] * (1-pop.allele.freqs.odds[2,i]))
      } else  if (s.alleles.8[j,i,h] == 0) {
        south.likelihoods.sdist.8[j,i,h] <- (1-pop.allele.freqs.odds[2,i])^2
      } else {
        south.likelihoods.sdist.8[j,i,h] <- 1
      }
    }
  }
}

# Multiply everything together
north8.ndist <- apply(north.likelihoods.ndist.8, FUN = prod, MARGIN = 1, na.rm = TRUE)
south8.ndist <- apply(south.likelihoods.ndist.8, FUN = prod, MARGIN = 1, na.rm = TRUE)

north8.sdist <- apply(north.likelihoods.sdist.8, FUN = prod, MARGIN = 1, na.rm = TRUE)
south8.sdist <- apply(south.likelihoods.sdist.8, FUN = prod, MARGIN = 1, na.rm = TRUE)

# Create ratio & plot
hist(log10(north8.ndist/south8.ndist), xlab = "log10(north likelihood/south likelihood)", main = "Population = 8", col = rgb(1,0,0,0.5), xlim = c(-10,10), ylim = c(0,270))
hist(log10(north8.sdist/south8.sdist), col = rgb(0,0,1,0.5), add = TRUE)
legend("topright", c("North", "South"), col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch = 15)

# Power test
log10(north8.ndist/south8.ndist)[order(log10(north8.ndist/south8.ndist))][50] #5% of 1000 is 50
p8 <- length(which(log10(north8.sdist/south8.sdist) < 0.2317728))/1000

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
log10(north10.ndist/south10.ndist)[order(log10(north10.ndist/south10.ndist))][50] #5% of 1000 is 50
p10 <- length(which(log10(north10.sdist/south10.sdist) < 0.7742818))/1000

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

# Power test
log10(north30.ndist/south30.ndist)[order(log10(north30.ndist/south30.ndist))][50] #5% of 1000 is 50
p30 <- length(which(log10(north30.sdist/south30.sdist) < 5.351101))/1000

###############################################################################################################
# Sample one allele per locus many times, with the z dimension being the number of individuals in a population
# Population = 45
n.alleles.45 <- array(dim = c(1000,10,45))
for (i in 1:45){
  n.alleles.45[,,i] <- sapply(pop.allele.freqs.odds.north,function(z){rbinom(1000,2,z)})
}

# Do this for the southern distribution too
s.alleles.45 <- array(dim = c(1000,10,45))
for (i in 1:45){
  s.alleles.45[,,i] <- sapply(pop.allele.freqs.odds.south,function(z){rbinom(1000,2,z)})
}

# Use adult allele frequencies to calculate likelihoods
# For loop to loop through each of 10 loci & multiply by the adult allele frequency, and then do this for all 1000 offspring. NA's/9's get coded as 1's so they don't make a difference when each row's product is taken
north.likelihoods.ndist.45 <- array(dim = c(1000,10,45))
south.likelihoods.ndist.45 <- array(dim = c(1000,10,45))

north.likelihoods.sdist.45 <- array(dim = c(1000,10,45))
south.likelihoods.sdist.45 <- array(dim = c(1000,10,45))

for (h in 1:45){
  
  for (j in 1:nrow(n.alleles.45)){
    
    for (i in 1:ncol(n.alleles.45)){
      if(n.alleles.45[j,i,h] == 2) {
        north.likelihoods.ndist.45[j,i,h] <- pop.allele.freqs.odds[1,i]^2
      } else if (n.alleles.45[j,i,h] == 1) {
        north.likelihoods.ndist.45[j,i,h] <- 2*(pop.allele.freqs.odds[1,i] * (1-pop.allele.freqs.odds[1,i]))
      } else if (n.alleles.45[j,i,h] == 0) {
        north.likelihoods.ndist.45[j,i,h] <- ( 1-pop.allele.freqs.odds[1,i])^2 
      } else {
        north.likelihoods.ndist.45[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(n.alleles.45)){
      if(n.alleles.45[j,i,h] == 2){
        south.likelihoods.ndist.45[j,i,h] <- pop.allele.freqs.odds[2,i]^2
      } else if (n.alleles.45[j,i,h] == 1) {
        south.likelihoods.ndist.45[j,i,h] <- 2*(pop.allele.freqs.odds[2,i] * (1-pop.allele.freqs.odds[2,i]))
      } else  if (n.alleles.45[j,i,h] == 0) {
        south.likelihoods.ndist.45[j,i,h] <- (1-pop.allele.freqs.odds[2,i])^2
      } else {
        south.likelihoods.ndist.45[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(s.alleles.45)){
    
    for (i in 1:ncol(s.alleles.45)){
      if(s.alleles.45[j,i,h] == 2) {
        north.likelihoods.sdist.45[j,i,h] <- pop.allele.freqs.odds[1,i]^2
      } else if (s.alleles.45[j,i,h] == 1) {
        north.likelihoods.sdist.45[j,i,h] <- 2*(pop.allele.freqs.odds[1,i] * (1-pop.allele.freqs.odds[1,i]))
      } else if (s.alleles.45[j,i,h] == 0) {
        north.likelihoods.sdist.45[j,i,h] <- (1-pop.allele.freqs.odds[1,i])^2 
      } else {
        north.likelihoods.sdist.45[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(s.alleles.45)){
      if(s.alleles.45[j,i,h] == 2){
        south.likelihoods.sdist.45[j,i,h] <- pop.allele.freqs.odds[2,i]^2
      } else if (s.alleles.45[j,i,h] == 1) {
        south.likelihoods.sdist.45[j,i,h] <- 2*(pop.allele.freqs.odds[2,i] * (1-pop.allele.freqs.odds[2,i]))
      } else  if (s.alleles.45[j,i,h] == 0) {
        south.likelihoods.sdist.45[j,i,h] <- (1-pop.allele.freqs.odds[2,i])^2
      } else {
        south.likelihoods.sdist.45[j,i,h] <- 1
      }
    }
  }
}

# Multiply everything together
north45.ndist <- apply(north.likelihoods.ndist.45, FUN = prod, MARGIN = 1, na.rm = TRUE)
south45.ndist <- apply(south.likelihoods.ndist.45, FUN = prod, MARGIN = 1, na.rm = TRUE)

north45.sdist <- apply(north.likelihoods.sdist.45, FUN = prod, MARGIN = 1, na.rm = TRUE)
south45.sdist <- apply(south.likelihoods.sdist.45, FUN = prod, MARGIN = 1, na.rm = TRUE)

# Create ratio & plot
hist(log10(north45.ndist/south45.ndist), xlab = "log10(north likelihood/south likelihood)", main = "Population = 45", col = rgb(1,0,0,0.5), xlim = c(-30,30), ylim = c(0,220))
hist(log10(north45.sdist/south45.sdist), col = rgb(0,0,1,0.5), add = TRUE)
legend("topright", c("North", "South"), col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch = 15)

# Power test
log10(north45.ndist/south45.ndist)[order(log10(north45.ndist/south45.ndist))][50] #5% of 1000 is 50
p45 <- length(which(log10(north45.sdist/south45.sdist) < 9.584772))/1000

###############################################################################################################
# Sample one allele per locus many times, with the z dimension being the number of individuals in a population
# Population = 80
n.alleles.80 <- array(dim = c(1000,10,80))
for (i in 1:80){
  n.alleles.80[,,i] <- sapply(pop.allele.freqs.odds.north,function(z){rbinom(1000,2,z)})
}

# Do this for the southern distribution too
s.alleles.80 <- array(dim = c(1000,10,80))
for (i in 1:80){
  s.alleles.80[,,i] <- sapply(pop.allele.freqs.odds.south,function(z){rbinom(1000,2,z)})
}

# Use adult allele frequencies to calculate likelihoods
# For loop to loop through each of 10 loci & multiply by the adult allele frequency, and then do this for all 1000 offspring. NA's/9's get coded as 1's so they don't make a difference when each row's product is taken
north.likelihoods.ndist.80 <- array(dim = c(1000,10,80))
south.likelihoods.ndist.80 <- array(dim = c(1000,10,80))

north.likelihoods.sdist.80 <- array(dim = c(1000,10,80))
south.likelihoods.sdist.80 <- array(dim = c(1000,10,80))

for (h in 1:80){
  
  for (j in 1:nrow(n.alleles.80)){
    
    for (i in 1:ncol(n.alleles.80)){
      if(n.alleles.80[j,i,h] == 2) {
        north.likelihoods.ndist.80[j,i,h] <- pop.allele.freqs.odds[1,i]^2
      } else if (n.alleles.80[j,i,h] == 1) {
        north.likelihoods.ndist.80[j,i,h] <- 2*(pop.allele.freqs.odds[1,i] * (1-pop.allele.freqs.odds[1,i]))
      } else if (n.alleles.80[j,i,h] == 0) {
        north.likelihoods.ndist.80[j,i,h] <- ( 1-pop.allele.freqs.odds[1,i])^2 
      } else {
        north.likelihoods.ndist.80[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(n.alleles.80)){
      if(n.alleles.80[j,i,h] == 2){
        south.likelihoods.ndist.80[j,i,h] <- pop.allele.freqs.odds[2,i]^2
      } else if (n.alleles.80[j,i,h] == 1) {
        south.likelihoods.ndist.80[j,i,h] <- 2*(pop.allele.freqs.odds[2,i] * (1-pop.allele.freqs.odds[2,i]))
      } else  if (n.alleles.80[j,i,h] == 0) {
        south.likelihoods.ndist.80[j,i,h] <- (1-pop.allele.freqs.odds[2,i])^2
      } else {
        south.likelihoods.ndist.80[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(s.alleles.80)){
    
    for (i in 1:ncol(s.alleles.80)){
      if(s.alleles.80[j,i,h] == 2) {
        north.likelihoods.sdist.80[j,i,h] <- pop.allele.freqs.odds[1,i]^2
      } else if (s.alleles.80[j,i,h] == 1) {
        north.likelihoods.sdist.80[j,i,h] <- 2*(pop.allele.freqs.odds[1,i] * (1-pop.allele.freqs.odds[1,i]))
      } else if (s.alleles.80[j,i,h] == 0) {
        north.likelihoods.sdist.80[j,i,h] <- (1-pop.allele.freqs.odds[1,i])^2 
      } else {
        north.likelihoods.sdist.80[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(s.alleles.80)){
      if(s.alleles.80[j,i,h] == 2){
        south.likelihoods.sdist.80[j,i,h] <- pop.allele.freqs.odds[2,i]^2
      } else if (s.alleles.80[j,i,h] == 1) {
        south.likelihoods.sdist.80[j,i,h] <- 2*(pop.allele.freqs.odds[2,i] * (1-pop.allele.freqs.odds[2,i]))
      } else  if (s.alleles.80[j,i,h] == 0) {
        south.likelihoods.sdist.80[j,i,h] <- (1-pop.allele.freqs.odds[2,i])^2
      } else {
        south.likelihoods.sdist.80[j,i,h] <- 1
      }
    }
  }
}

# Multiply everything together
north80.ndist <- apply(north.likelihoods.ndist.80, FUN = prod, MARGIN = 1, na.rm = TRUE)
south80.ndist <- apply(south.likelihoods.ndist.80, FUN = prod, MARGIN = 1, na.rm = TRUE)

north80.sdist <- apply(north.likelihoods.sdist.80, FUN = prod, MARGIN = 1, na.rm = TRUE)
south80.sdist <- apply(south.likelihoods.sdist.80, FUN = prod, MARGIN = 1, na.rm = TRUE)

# Create ratio & plot
hist(log10(north80.ndist/south80.ndist), xlab = "log10(north likelihood/south likelihood)", main = "Population = 80", col = rgb(1,0,0,0.5), xlim = c(-50,45), ylim = c(0,400))
hist(log10(north80.sdist/south80.sdist), col = rgb(0,0,1,0.5), add = TRUE)
legend("topright", c("North", "South"), col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch = 15)

# Power test
log10(north80.ndist/south80.ndist)[order(log10(north80.ndist/south80.ndist))][50] #5% of 1000 is 50
p80 <- length(which(log10(north80.sdist/south80.sdist) < 19.88935))/1000

# Plot power curve
# code to generate distributions of population size 1 is in otoliths.R
plot(c(p1, p5, p8, p10, p30, p45, p80) ~ c(1, 5,8,10,30,45,80), xlab = 'Population size', ylab = 'Power', pch = 19, col = 'darkslategray4')

#######################################################################################################
#### Power analysis using the five BayEnv populations ####

# Read in adult outlier allele frequencies
pop.allele.freqs5 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.freqs.5pops.txt')

# Keep only one allele per locus
odds <- seq(1,20,2) # odd indicies to keep
pop.allele.freqs5.odds <- pop.allele.freqs5[,odds]

# Sample one allele per locus many times from each BayEnv population, with the z dimension being the number of individuals in a 'cluster'
# Population = 5
Pop1.alleles.5 <- array(dim = c(10000,10,5))
for (i in 1:5){
  Pop1.alleles.5[,,i] <- sapply(pop.allele.freqs5.odds[1,],function(z){rbinom(10000,2,z)})
}

# Do this for the rest of the BayEnv populations too
Pop2.alleles.5 <- array(dim = c(10000,10,5))
for (i in 1:5){
  Pop2.alleles.5[,,i] <- sapply(pop.allele.freqs5.odds[2,],function(z){rbinom(10000,2,z)})
}

Pop3.alleles.5 <- array(dim = c(10000,10,5))
for (i in 1:5){
  Pop3.alleles.5[,,i] <- sapply(pop.allele.freqs5.odds[3,],function(z){rbinom(10000,2,z)})
}

Pop4.alleles.5 <- array(dim = c(10000,10,5))
for (i in 1:5){
  Pop4.alleles.5[,,i] <- sapply(pop.allele.freqs5.odds[4,],function(z){rbinom(10000,2,z)})
}

Pop5.alleles.5 <- array(dim = c(10000,10,5))
for (i in 1:5){
  Pop5.alleles.5[,,i] <- sapply(pop.allele.freqs5.odds[5,],function(z){rbinom(10000,2,z)})
}

# Use adult allele frequencies to calculate likelihoods
# For loop to loop through each of 10 loci & multiply by the adult allele frequency, and then do this for all 1000 offspring. NA's/9's get coded as 1's so they don't make a difference when each row's product is taken
Pop1.likelihoods.pop1dist.5 <- array(dim = c(10000,10,5))
Pop2.likelihoods.pop1dist.5 <- array(dim = c(10000,10,5))
Pop3.likelihoods.pop1dist.5 <- array(dim = c(10000,10,5))
Pop4.likelihoods.pop1dist.5 <- array(dim = c(10000,10,5))
Pop5.likelihoods.pop1dist.5 <- array(dim = c(10000,10,5))

Pop1.likelihoods.pop2dist.5 <- array(dim = c(10000,10,5))
Pop2.likelihoods.pop2dist.5 <- array(dim = c(10000,10,5))
Pop3.likelihoods.pop2dist.5 <- array(dim = c(10000,10,5))
Pop4.likelihoods.pop2dist.5 <- array(dim = c(10000,10,5))
Pop5.likelihoods.pop2dist.5 <- array(dim = c(10000,10,5))

Pop1.likelihoods.pop3dist.5 <- array(dim = c(10000,10,5))
Pop2.likelihoods.pop3dist.5 <- array(dim = c(10000,10,5))
Pop3.likelihoods.pop3dist.5 <- array(dim = c(10000,10,5))
Pop4.likelihoods.pop3dist.5 <- array(dim = c(10000,10,5))
Pop5.likelihoods.pop3dist.5 <- array(dim = c(10000,10,5))

Pop1.likelihoods.pop4dist.5 <- array(dim = c(10000,10,5))
Pop2.likelihoods.pop4dist.5 <- array(dim = c(10000,10,5))
Pop3.likelihoods.pop4dist.5 <- array(dim = c(10000,10,5))
Pop4.likelihoods.pop4dist.5 <- array(dim = c(10000,10,5))
Pop5.likelihoods.pop4dist.5 <- array(dim = c(10000,10,5))

Pop1.likelihoods.pop5dist.5 <- array(dim = c(10000,10,5))
Pop2.likelihoods.pop5dist.5 <- array(dim = c(10000,10,5))
Pop3.likelihoods.pop5dist.5 <- array(dim = c(10000,10,5))
Pop4.likelihoods.pop5dist.5 <- array(dim = c(10000,10,5))
Pop5.likelihoods.pop5dist.5 <- array(dim = c(10000,10,5))

for (h in 1:5){
  
  for (j in 1:nrow(Pop1.alleles.5)){
    
    for (i in 1:ncol(Pop1.alleles.5)){
      if(Pop1.alleles.5[j,i,h] == 2) {
        Pop1.likelihoods.pop1dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop1',i]^2
      } else if (Pop1.alleles.5[j,i,h] == 1) {
        Pop1.likelihoods.pop1dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop1',i] * (1-pop.allele.freqs5.odds['Pop1',i]))
      } else if (Pop1.alleles.5[j,i,h] == 0) {
        Pop1.likelihoods.pop1dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop1',i])^2 
      } else {
        Pop1.likelihoods.pop1dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop1.alleles.5)){
      if(Pop1.alleles.5[j,i,h] == 2) {
        Pop2.likelihoods.pop1dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop2',i]^2
      } else if (Pop1.alleles.5[j,i,h] == 1) {
        Pop2.likelihoods.pop1dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop2',i] * (1-pop.allele.freqs5.odds['Pop2',i]))
      } else if (Pop1.alleles.5[j,i,h] == 0) {
        Pop2.likelihoods.pop1dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop2',i])^2 
      } else {
        Pop2.likelihoods.pop1dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop1.alleles.5)){
      if(Pop1.alleles.5[j,i,h] == 2) {
        Pop3.likelihoods.pop1dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop3',i]^2
      } else if (Pop1.alleles.5[j,i,h] == 1) {
        Pop3.likelihoods.pop1dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop3',i] * (1-pop.allele.freqs5.odds['Pop3',i]))
      } else if (Pop1.alleles.5[j,i,h] == 0) {
        Pop3.likelihoods.pop1dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop3',i])^2 
      } else {
        Pop3.likelihoods.pop1dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop1.alleles.5)){
      if(Pop1.alleles.5[j,i,h] == 2) {
        Pop4.likelihoods.pop1dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop4',i]^2
      } else if (Pop1.alleles.5[j,i,h] == 1) {
        Pop4.likelihoods.pop1dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop4',i] * (1-pop.allele.freqs5.odds['Pop4',i]))
      } else if (Pop1.alleles.5[j,i,h] == 0) {
        Pop4.likelihoods.pop1dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop4',i])^2 
      } else {
        Pop4.likelihoods.pop1dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop1.alleles.5)){
      if(Pop1.alleles.5[j,i,h] == 2) {
        Pop5.likelihoods.pop1dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop5',i]^2
      } else if (Pop1.alleles.5[j,i,h] == 1) {
        Pop5.likelihoods.pop1dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop5',i] * (1-pop.allele.freqs5.odds['Pop5',i]))
      } else if (Pop1.alleles.5[j,i,h] == 0) {
        Pop5.likelihoods.pop1dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop5',i])^2 
      } else {
        Pop5.likelihoods.pop1dist.5[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop2.alleles.5)){
    
    for (i in 1:ncol(Pop2.alleles.5)){
      if(Pop2.alleles.5[j,i,h] == 2) {
        Pop1.likelihoods.pop2dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop1',i]^2
      } else if (Pop2.alleles.5[j,i,h] == 1) {
        Pop1.likelihoods.pop2dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop1',i] * (1-pop.allele.freqs5.odds['Pop1',i]))
      } else if (Pop2.alleles.5[j,i,h] == 0) {
        Pop1.likelihoods.pop2dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop1',i])^2 
      } else {
        Pop1.likelihoods.pop2dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop2.alleles.5)){
      if(Pop2.alleles.5[j,i,h] == 2) {
        Pop2.likelihoods.pop2dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop2',i]^2
      } else if (Pop2.alleles.5[j,i,h] == 1) {
        Pop2.likelihoods.pop2dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop2',i] * (1-pop.allele.freqs5.odds['Pop2',i]))
      } else if (Pop2.alleles.5[j,i,h] == 0) {
        Pop2.likelihoods.pop2dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop2',i])^2 
      } else {
        Pop2.likelihoods.pop2dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop2.alleles.5)){
      if(Pop2.alleles.5[j,i,h] == 2) {
        Pop3.likelihoods.pop2dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop3',i]^2
      } else if (Pop2.alleles.5[j,i,h] == 1) {
        Pop3.likelihoods.pop2dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop3',i] * (1-pop.allele.freqs5.odds['Pop3',i]))
      } else if (Pop2.alleles.5[j,i,h] == 0) {
        Pop3.likelihoods.pop2dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop3',i])^2 
      } else {
        Pop3.likelihoods.pop2dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop2.alleles.5)){
      if(Pop2.alleles.5[j,i,h] == 2) {
        Pop4.likelihoods.pop2dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop4',i]^2
      } else if (Pop2.alleles.5[j,i,h] == 1) {
        Pop4.likelihoods.pop2dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop4',i] * (1-pop.allele.freqs5.odds['Pop4',i]))
      } else if (Pop2.alleles.5[j,i,h] == 0) {
        Pop4.likelihoods.pop2dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop4',i])^2 
      } else {
        Pop4.likelihoods.pop2dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop2.alleles.5)){
      if(Pop2.alleles.5[j,i,h] == 2) {
        Pop5.likelihoods.pop2dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop5',i]^2
      } else if (Pop2.alleles.5[j,i,h] == 1) {
        Pop5.likelihoods.pop2dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop5',i] * (1-pop.allele.freqs5.odds['Pop5',i]))
      } else if (Pop2.alleles.5[j,i,h] == 0) {
        Pop5.likelihoods.pop2dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop5',i])^2 
      } else {
        Pop5.likelihoods.pop2dist.5[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop3.alleles.5)){
    
    for (i in 1:ncol(Pop3.alleles.5)){
      if(Pop3.alleles.5[j,i,h] == 2) {
        Pop1.likelihoods.pop3dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop1',i]^2
      } else if (Pop3.alleles.5[j,i,h] == 1) {
        Pop1.likelihoods.pop3dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop1',i] * (1-pop.allele.freqs5.odds['Pop1',i]))
      } else if (Pop3.alleles.5[j,i,h] == 0) {
        Pop1.likelihoods.pop3dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop1',i])^2 
      } else {
        Pop1.likelihoods.pop3dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop3.alleles.5)){
      if(Pop3.alleles.5[j,i,h] == 2) {
        Pop2.likelihoods.pop3dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop2',i]^2
      } else if (Pop3.alleles.5[j,i,h] == 1) {
        Pop2.likelihoods.pop3dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop2',i] * (1-pop.allele.freqs5.odds['Pop2',i]))
      } else if (Pop3.alleles.5[j,i,h] == 0) {
        Pop2.likelihoods.pop3dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop2',i])^2 
      } else {
        Pop2.likelihoods.pop3dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop3.alleles.5)){
      if(Pop3.alleles.5[j,i,h] == 2) {
        Pop3.likelihoods.pop3dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop3',i]^2
      } else if (Pop3.alleles.5[j,i,h] == 1) {
        Pop3.likelihoods.pop3dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop3',i] * (1-pop.allele.freqs5.odds['Pop3',i]))
      } else if (Pop3.alleles.5[j,i,h] == 0) {
        Pop3.likelihoods.pop3dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop3',i])^2 
      } else {
        Pop3.likelihoods.pop3dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop3.alleles.5)){
      if(Pop3.alleles.5[j,i,h] == 2) {
        Pop4.likelihoods.pop3dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop4',i]^2
      } else if (Pop3.alleles.5[j,i,h] == 1) {
        Pop4.likelihoods.pop3dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop4',i] * (1-pop.allele.freqs5.odds['Pop4',i]))
      } else if (Pop3.alleles.5[j,i,h] == 0) {
        Pop4.likelihoods.pop3dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop4',i])^2 
      } else {
        Pop4.likelihoods.pop3dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop3.alleles.5)){
      if(Pop3.alleles.5[j,i,h] == 2) {
        Pop5.likelihoods.pop3dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop5',i]^2
      } else if (Pop3.alleles.5[j,i,h] == 1) {
        Pop5.likelihoods.pop3dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop5',i] * (1-pop.allele.freqs5.odds['Pop5',i]))
      } else if (Pop3.alleles.5[j,i,h] == 0) {
        Pop5.likelihoods.pop3dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop5',i])^2 
      } else {
        Pop5.likelihoods.pop3dist.5[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop4.alleles.5)){
    
    for (i in 1:ncol(Pop4.alleles.5)){
      if(Pop4.alleles.5[j,i,h] == 2) {
        Pop1.likelihoods.pop4dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop1',i]^2
      } else if (Pop4.alleles.5[j,i,h] == 1) {
        Pop1.likelihoods.pop4dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop1',i] * (1-pop.allele.freqs5.odds['Pop1',i]))
      } else if (Pop4.alleles.5[j,i,h] == 0) {
        Pop1.likelihoods.pop4dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop1',i])^2 
      } else {
        Pop1.likelihoods.pop4dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop4.alleles.5)){
      if(Pop4.alleles.5[j,i,h] == 2) {
        Pop2.likelihoods.pop4dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop2',i]^2
      } else if (Pop4.alleles.5[j,i,h] == 1) {
        Pop2.likelihoods.pop4dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop2',i] * (1-pop.allele.freqs5.odds['Pop2',i]))
      } else if (Pop4.alleles.5[j,i,h] == 0) {
        Pop2.likelihoods.pop4dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop2',i])^2 
      } else {
        Pop2.likelihoods.pop4dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop4.alleles.5)){
      if(Pop4.alleles.5[j,i,h] == 2) {
        Pop3.likelihoods.pop4dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop3',i]^2
      } else if (Pop4.alleles.5[j,i,h] == 1) {
        Pop3.likelihoods.pop4dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop3',i] * (1-pop.allele.freqs5.odds['Pop3',i]))
      } else if (Pop4.alleles.5[j,i,h] == 0) {
        Pop3.likelihoods.pop4dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop3',i])^2 
      } else {
        Pop3.likelihoods.pop4dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop4.alleles.5)){
      if(Pop4.alleles.5[j,i,h] == 2) {
        Pop4.likelihoods.pop4dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop4',i]^2
      } else if (Pop4.alleles.5[j,i,h] == 1) {
        Pop4.likelihoods.pop4dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop4',i] * (1-pop.allele.freqs5.odds['Pop4',i]))
      } else if (Pop4.alleles.5[j,i,h] == 0) {
        Pop4.likelihoods.pop4dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop4',i])^2 
      } else {
        Pop4.likelihoods.pop4dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop4.alleles.5)){
      if(Pop4.alleles.5[j,i,h] == 2) {
        Pop5.likelihoods.pop4dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop5',i]^2
      } else if (Pop4.alleles.5[j,i,h] == 1) {
        Pop5.likelihoods.pop4dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop5',i] * (1-pop.allele.freqs5.odds['Pop5',i]))
      } else if (Pop4.alleles.5[j,i,h] == 0) {
        Pop5.likelihoods.pop4dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop5',i])^2 
      } else {
        Pop5.likelihoods.pop4dist.5[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop5.alleles.5)){
    
    for (i in 1:ncol(Pop5.alleles.5)){
      if(Pop5.alleles.5[j,i,h] == 2) {
        Pop1.likelihoods.pop5dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop1',i]^2
      } else if (Pop5.alleles.5[j,i,h] == 1) {
        Pop1.likelihoods.pop5dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop1',i] * (1-pop.allele.freqs5.odds['Pop1',i]))
      } else if (Pop5.alleles.5[j,i,h] == 0) {
        Pop1.likelihoods.pop5dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop1',i])^2 
      } else {
        Pop1.likelihoods.pop5dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop5.alleles.5)){
      if(Pop5.alleles.5[j,i,h] == 2) {
        Pop2.likelihoods.pop5dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop2',i]^2
      } else if (Pop5.alleles.5[j,i,h] == 1) {
        Pop2.likelihoods.pop5dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop2',i] * (1-pop.allele.freqs5.odds['Pop2',i]))
      } else if (Pop5.alleles.5[j,i,h] == 0) {
        Pop2.likelihoods.pop5dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop2',i])^2 
      } else {
        Pop2.likelihoods.pop5dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop5.alleles.5)){
      if(Pop5.alleles.5[j,i,h] == 2) {
        Pop3.likelihoods.pop5dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop3',i]^2
      } else if (Pop5.alleles.5[j,i,h] == 1) {
        Pop3.likelihoods.pop5dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop3',i] * (1-pop.allele.freqs5.odds['Pop3',i]))
      } else if (Pop5.alleles.5[j,i,h] == 0) {
        Pop3.likelihoods.pop5dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop3',i])^2 
      } else {
        Pop3.likelihoods.pop5dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop5.alleles.5)){
      if(Pop5.alleles.5[j,i,h] == 2) {
        Pop4.likelihoods.pop5dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop4',i]^2
      } else if (Pop5.alleles.5[j,i,h] == 1) {
        Pop4.likelihoods.pop5dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop4',i] * (1-pop.allele.freqs5.odds['Pop4',i]))
      } else if (Pop5.alleles.5[j,i,h] == 0) {
        Pop4.likelihoods.pop5dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop4',i])^2 
      } else {
        Pop4.likelihoods.pop5dist.5[j,i,h] <- 1
      }
    }
    
    for (i in 1:ncol(Pop5.alleles.5)){
      if(Pop5.alleles.5[j,i,h] == 2) {
        Pop5.likelihoods.pop5dist.5[j,i,h] <- pop.allele.freqs5.odds['Pop5',i]^2
      } else if (Pop5.alleles.5[j,i,h] == 1) {
        Pop5.likelihoods.pop5dist.5[j,i,h] <- 2*(pop.allele.freqs5.odds['Pop5',i] * (1-pop.allele.freqs5.odds['Pop5',i]))
      } else if (Pop5.alleles.5[j,i,h] == 0) {
        Pop5.likelihoods.pop5dist.5[j,i,h] <- ( 1-pop.allele.freqs5.odds['Pop5',i])^2 
      } else {
        Pop5.likelihoods.pop5dist.5[j,i,h] <- 1
      }
    }
  }
}












# Multiply everything together
Pop1.pop1dist5 <- apply(Pop1.likelihoods.pop1dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop2.pop1dist5 <- apply(Pop2.likelihoods.pop1dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop3.pop1dist5 <- apply(Pop3.likelihoods.pop1dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop4.pop1dist5 <- apply(Pop4.likelihoods.pop1dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop5.pop1dist5 <- apply(Pop5.likelihoods.pop1dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)

Pop1.pop2dist5 <- apply(Pop1.likelihoods.pop2dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop2.pop2dist5 <- apply(Pop2.likelihoods.pop2dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop3.pop2dist5 <- apply(Pop3.likelihoods.pop2dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop4.pop2dist5 <- apply(Pop4.likelihoods.pop2dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop5.pop2dist5 <- apply(Pop5.likelihoods.pop2dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)

Pop1.pop3dist5 <- apply(Pop1.likelihoods.pop3dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop2.pop3dist5 <- apply(Pop2.likelihoods.pop3dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop3.pop3dist5 <- apply(Pop3.likelihoods.pop3dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop4.pop3dist5 <- apply(Pop4.likelihoods.pop3dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop5.pop3dist5 <- apply(Pop5.likelihoods.pop3dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)

Pop1.pop4dist5 <- apply(Pop1.likelihoods.pop4dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop2.pop4dist5 <- apply(Pop2.likelihoods.pop4dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop3.pop4dist5 <- apply(Pop3.likelihoods.pop4dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop4.pop4dist5 <- apply(Pop4.likelihoods.pop4dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop5.pop4dist5 <- apply(Pop5.likelihoods.pop4dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)

Pop1.pop5dist5 <- apply(Pop1.likelihoods.pop5dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop2.pop5dist5 <- apply(Pop2.likelihoods.pop5dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop3.pop5dist5 <- apply(Pop3.likelihoods.pop5dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop4.pop5dist5 <- apply(Pop4.likelihoods.pop5dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop5.pop5dist5 <- apply(Pop5.likelihoods.pop5dist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)



north5.sdist <- apply(north.likelihoods.sdist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)
south5.sdist <- apply(south.likelihoods.sdist.5, FUN = prod, MARGIN = 1, na.rm = TRUE)







# Create ratio & plot
hist(log10(Pop1.pop1dist5/Pop1.pop2dist5), xlab = "log10(Pop1 likelihood/Pop2 likelihood)", main = "Cluster size = 5", col = rgb(1,0,0,0.5), xlim = c(-10,10), ylim = c(0,3000))
hist(log10(Pop2.pop1dist5/Pop2.pop2dist5), col = rgb(0,0,1,0.5), add = TRUE)


hist(log10(Pop1.pop1dist5/Pop2.pop1dist5), xlab = "log10(Pop1 likelihood/Pop2 likelihood)", main = "Cluster size = 5", col = rgb(1,0,0,0.5), xlim = c(-5,5), ylim = c(0,3000)) # I think this is the correct way to calculate likelihood ratio
hist(log10(Pop1.pop2dist5/Pop2.pop2dist5), col = rgb(0,0,1,0.5), add = TRUE)
hist(log10(Pop1.pop3dist5/Pop2.pop3dist5), col = rgb(1,0,1,0.5), add = TRUE)
hist(log10(Pop1.pop4dist5/Pop2.pop4dist5), col = rgb(1,1,0.5,0.5), add = TRUE)
hist(log10(Pop1.pop5dist5/Pop2.pop5dist5), col = rgb(0,1,0,0.5), add = TRUE)
legend("topright", c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5"), pch = 22, col = "black", pt.bg = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5), rgb(1,0,1,0.5), rgb(1,1,0.5,0.5), rgb(0,1,0,0.5)))

hist(log10(Pop1.pop1dist5/Pop3.pop1dist5), xlab = "log10(Pop1 likelihood/Pop3 likelihood)", main = "Cluster size = 5", col = rgb(1,0,0,0.5), xlim = c(-5,5), ylim = c(0,3250)) 
hist(log10(Pop1.pop2dist5/Pop3.pop2dist5), col = rgb(0,0,1,0.5), add = TRUE)
hist(log10(Pop1.pop3dist5/Pop3.pop3dist5), col = rgb(1,0,1,0.5), add = TRUE)
hist(log10(Pop1.pop4dist5/Pop3.pop4dist5), col = rgb(1,1,0.5,0.5), add = TRUE)
hist(log10(Pop1.pop5dist5/Pop3.pop5dist5), col = rgb(0,1,0,0.5), add = TRUE)
legend("topright", c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5"), pch = 22, col = "black", pt.bg = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5), rgb(1,0,1,0.5), rgb(1,1,0.5,0.5), rgb(0,1,0,0.5)))

hist(log10(Pop1.pop1dist5/Pop4.pop1dist5), xlab = "log10(Pop1 likelihood/Pop4 likelihood)", main = "Cluster size = 5", col = rgb(1,0,0,0.5), xlim = c(-7,8), ylim = c(0,2250)) 
hist(log10(Pop1.pop2dist5/Pop4.pop2dist5), col = rgb(0,0,1,0.5), add = TRUE)
hist(log10(Pop1.pop3dist5/Pop4.pop3dist5), col = rgb(1,0,1,0.5), add = TRUE)
hist(log10(Pop1.pop4dist5/Pop4.pop4dist5), col = rgb(1,1,0.5,0.5), add = TRUE)
hist(log10(Pop1.pop5dist5/Pop4.pop5dist5), col = rgb(0,1,0,0.5), add = TRUE)
legend("topright", c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5"), pch = 22, col = "black", pt.bg  = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5), rgb(1,0,1,0.5), rgb(1,1,0.5,0.5), rgb(0,1,0,0.5)))

hist(log10(Pop1.pop1dist5/Pop5.pop1dist5), xlab = "log10(Pop1 likelihood/Pop5 likelihood)", main = "Cluster size = 5", col = rgb(1,0,0,0.5), xlim = c(-12,12), ylim = c(0,2500)) 
hist(log10(Pop1.pop2dist5/Pop5.pop2dist5), col = rgb(0,0,1,0.5), add = TRUE)
hist(log10(Pop1.pop3dist5/Pop5.pop3dist5), col = rgb(1,0,1,0.5), add = TRUE)
hist(log10(Pop1.pop4dist5/Pop5.pop4dist5), col = rgb(1,1,0.5,0.5), add = TRUE)
hist(log10(Pop1.pop5dist5/Pop5.pop5dist5), col = rgb(0,1,0,0.5), add = TRUE)
legend("topright", c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5"), pch = 22, col = "black", pt.bg  = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5), rgb(1,0,1,0.5), rgb(1,1,0.5,0.5), rgb(0,1,0,0.5)))

# Test statistic as -log10(genotype likelihood) vs a ratio
hist(-log10(Pop1.pop2dist5), col = rgb(1,0,0,0.5), xlim = c(2,34), ylim = c(0, 3000), xlab = "-log10(genotype likelihood)", main = "")
hist(-log10(Pop2.pop2dist5), add = TRUE, col = rgb(0,0,1,0.5))
hist(-log10(Pop3.pop2dist5), add = TRUE, col = rgb(1,0,1,0.5))
hist(-log10(Pop4.pop2dist5), add = TRUE, col = rgb(1,1,0.5,0.5))
hist(-log10(Pop5.pop2dist5), add = TRUE, col = rgb(0,1,0,0.5))

hist(-log10(Pop1.pop1dist5), col = rgb(1,0,0,0.5), xlim = c(2,34), ylim = c(0, 3000), xlab = "-log10(genotype likelihood)", main = "")
hist(-log10(Pop1.pop2dist5), add = TRUE, col = rgb(0,0,1,0.5))
hist(-log10(Pop1.pop3dist5), add = TRUE, col = rgb(1,0,1,0.5))
hist(-log10(Pop1.pop4dist5), add = TRUE, col = rgb(1,1,0.5,0.5))
hist(-log10(Pop1.pop5dist5), add = TRUE, col = rgb(0,1,0,0.5))





hist(log10(north5.ndist/south5.ndist), xlab = "log10(north likelihood/south likelihood)", main = "Population = 5", col = rgb(1,0,0,0.5), xlim = c(-10,10), ylim = c(0,325))
hist(log10(north5.sdist/south5.sdist), col = rgb(0,0,1,0.5), add = TRUE)
legend("topright", c("North", "South"), col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch = 15)

# Power test
log10(north5.ndist/south5.ndist)[order(log10(north5.ndist/south5.ndist))][50] #5% of 1000 is 50
p5 <- length(which(log10(north5.sdist/south5.sdist) < -0.3087461))/1000

