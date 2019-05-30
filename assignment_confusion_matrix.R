#### Generating a confusion matrix of individual assignment to see how accurate individual assignment is ###
#### 
#### This script simulates allele frequencies at 10 outlier adult summer flounder loci of 10 populations based on GAMS with distance. The likelihood of each genotype 
#### is calculated using the GAM allele frequency for each of 10 groups. The number of correct and incorrect assignments is then calculated.

# Read in adult outlier allele frequencies based on 10 populations
adult10.freqs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/pop.allele.freqs.10pops.gambinomial.ordered.txt', header = TRUE) # GAM with distance and binomial errors: these are the alleles used to calculate genotype likelihoods of clusters; 10 x 10, one allele per locus

################################################################################################
# Sample one allele per locus 1000 times from each of 10 GAM populations
Pop1.alleles <- data.frame() # 1//J is southern-most, 10//A is northern-most
for (i in 1:ncol(adult10.freqs[1,])){
  Pop1.alleles <- sapply(adult10.freqs[1,],function(z){rbinom(1000,2,z)})
}

# Do this for the rest of the populations too
Pop2.alleles <- data.frame() # 1//J is southern-most, 10//A is northern-most
for (i in 1:ncol(adult10.freqs[2,])){
  Pop2.alleles <- sapply(adult10.freqs[2,],function(z){rbinom(1000,2,z)})
}

Pop3.alleles <- data.frame() # 1//J is southern-most, 10//A is northern-most
for (i in 1:ncol(adult10.freqs[3,])){
  Pop3.alleles <- sapply(adult10.freqs[3,],function(z){rbinom(1000,2,z)})
}

Pop4.alleles <- data.frame() # 1//J is southern-most, 10//A is northern-most
for (i in 1:ncol(adult10.freqs[4,])){
  Pop4.alleles <- sapply(adult10.freqs[4,],function(z){rbinom(1000,2,z)})
}

Pop5.alleles <- data.frame() # 1//J is southern-most, 10//A is northern-most
for (i in 1:ncol(adult10.freqs[5,])){
  Pop5.alleles <- sapply(adult10.freqs[5,],function(z){rbinom(1000,2,z)})
}

Pop6.alleles <- data.frame() # 1//J is southern-most, 10//A is northern-most
for (i in 1:ncol(adult10.freqs[6,])){
  Pop6.alleles <- sapply(adult10.freqs[6,],function(z){rbinom(1000,2,z)})
}

Pop7.alleles <- data.frame() # 1//J is southern-most, 10//A is northern-most
for (i in 1:ncol(adult10.freqs[7,])){
  Pop7.alleles <- sapply(adult10.freqs[7,],function(z){rbinom(1000,2,z)})
}

Pop8.alleles <- data.frame() # 1//J is southern-most, 10//A is northern-most
for (i in 1:ncol(adult10.freqs[8,])){
  Pop8.alleles <- sapply(adult10.freqs[8,],function(z){rbinom(1000,2,z)})
}

Pop9.alleles <- data.frame() # 1//J is southern-most, 10//A is northern-most
for (i in 1:ncol(adult10.freqs[9,])){
  Pop9.alleles <- sapply(adult10.freqs[9,],function(z){rbinom(1000,2,z)})
}

Pop10.alleles <- data.frame() # 1//J is southern-most, 10//A is northern-most
for (i in 1:ncol(adult10.freqs[10,])){
  Pop10.alleles <- sapply(adult10.freqs[10,],function(z){rbinom(1000,2,z)})
}

####################################################################################################################################
#### Now that I have allele frequencies for 10 adult groups from the GAMS and the allele counts for the simulated larvae from 10 GAM pops, 
#### I'm ready to calculate the genotype likelihood of each larval indivdiual coming from the 10 potential popuations of origin. This means 
#### that for each set of individuals similated from one GAM pop, I will be calculating geno likelihoods across all 10 pops (1000 x 10 x 10).
#### This will then be done for each of the GAM populations (10)

# Population likelihoods for 10 groups
pop1.likelihoods <- array(dim = c(1000,10,10))
pop2.likelihoods <- array(dim = c(1000,10,10))
# pop3.likelihoods <- array(dim = c(1000,10,10))
# pop4.likelihoods <- array(dim = c(1000,10,10))
# pop5.likelihoods <- array(dim = c(1000,10,10))
# pop6.likelihoods <- array(dim = c(1000,10,10))
# pop7.likelihoods <- array(dim = c(1000,10,10))
# pop8.likelihoods <- array(dim = c(1000,10,10))
# pop9.likelihoods <- array(dim = c(1000,10,10))
# pop10.likelihoods <- array(dim = c(1000,10,10))

# For larvae simulated using GAM pop 1 allele freqs. Then geno likelihoods calculated for all 10 GAM pops
for (h in 1:nrow(adult10.freqs)){
  
  for (j in 1:nrow(Pop1.alleles)){
    
    for (i in 1:ncol(Pop1.alleles)){
      if(Pop1.alleles[j,i] == 2) {
        pop1.likelihoods[j,i,h] <- adult10.freqs[h,i]^2
      } else if (Pop1.alleles[j,i] == 1) {
        pop1.likelihoods[j,i,h] <- 2*(adult10.freqs[h,i] * (1-adult10.freqs[h,i]))
      } else if (Pop1.alleles[j,i] == 0) {
        pop1.likelihoods[j,i,h] <- ( 1-adult10.freqs[h,i])^2 
      } else {
        pop1.likelihoods[j,i,h] <- 1
      }
    }
  }
}

# For larvae simulated using GAM pop 2 allele freqs. 
for (h in 1:nrow(adult10.freqs)){
  
  for (j in 1:nrow(Pop2.alleles)){
    
    for (i in 1:ncol(Pop2.alleles)){
      if(Pop2.alleles[j,i] == 2) {
        pop2.likelihoods[j,i,h] <- adult10.freqs[h,i]^2
      } else if (Pop2.alleles[j,i] == 1) {
        pop2.likelihoods[j,i,h] <- 2*(adult10.freqs[h,i] * (1-adult10.freqs[h,i]))
      } else if (Pop2.alleles[j,i] == 0) {
        pop2.likelihoods[j,i,h] <- ( 1-adult10.freqs[h,i])^2 
      } else {
        pop2.likelihoods[j,i,h] <- 1
      }
    }
  }
}

# GAM pop 3 and so on......
for (h in 1:nrow(adult10.freqs)){
  
  for (j in 1:nrow(Pop3.alleles)){
    
    for (i in 1:ncol(Pop3.alleles)){
      if(Pop3.alleles[j,i] == 2) {
        pop3.likelihoods[j,i,h] <- adult10.freqs[h,i]^2
      } else if (Pop3.alleles[j,i] == 1) {
        pop3.likelihoods[j,i,h] <- 2*(adult10.freqs[h,i] * (1-adult10.freqs[h,i]))
      } else if (Pop3.alleles[j,i] == 0) {
        pop3.likelihoods[j,i,h] <- ( 1-adult10.freqs[h,i])^2 
      } else {
        pop3.likelihoods[j,i,h] <- 1
      }
    }
  }
}





#### Multiply rows together and then assign each individual to a GAM population ####
Pop1.ass <- list()
for (l in 1:10){
Pop1.ass[[l]] <- apply(pop1.likelihoods[,,l], FUN = prod, MARGIN = 1, na.rm = TRUE)
  }

Pop1.ass.df <- data.frame(matrix(unlist(Pop1.ass), nrow=nrow(pop1.likelihoods), byrow=F))
colnames(Pop1.ass.df) <- c('Pop1', 'Pop2', 'Pop3', 'Pop4', 'Pop5', 'Pop6', 'Pop7', 'Pop8', 'Pop9', 'Pop10')
Pop1.most.like <- colnames(Pop1.ass.df)[apply(Pop1.ass.df,1, which.max)]

Pop2.ass <- list()
for (l in 1:10){
  Pop2.ass[[l]] <- apply(pop2.likelihoods[,,l], FUN = prod, MARGIN = 1, na.rm = TRUE)
}
  
Pop2.ass.df <- data.frame(matrix(unlist(Pop2.ass), nrow=nrow(pop2.likelihoods), byrow=F))
colnames(Pop2.ass.df) <- c('Pop1', 'Pop2', 'Pop3', 'Pop4', 'Pop5', 'Pop6', 'Pop7', 'Pop8', 'Pop9', 'Pop10')
Pop2.most.like <- colnames(Pop2.ass.df)[apply(Pop2.ass.df,1, which.max)]





Pop3.pop3ass <- apply(pop3.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop4.pop4ass <- apply(pop4.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop5.pop5ass <- apply(pop5.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop6.pop6ass <- apply(pop6.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop7.pop7ass <- apply(pop7.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop8.pop8ass <- apply(pop8.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop9.pop9ass <- apply(pop9.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop10.pop10ass <- apply(pop10.likelihoods, FUN = prod, MARGIN = 1, na.rm = TRUE)

#### Rbind & organize the assignment results for individuals simulated from each of 10 GAM populations ####
# Row as reference? Column as predicted?
table(Pop1.most.like)

#### Plot confusion matrix ####
