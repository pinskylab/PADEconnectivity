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
pop3.likelihoods <- array(dim = c(1000,10,10))
pop4.likelihoods <- array(dim = c(1000,10,10))
pop5.likelihoods <- array(dim = c(1000,10,10))
pop6.likelihoods <- array(dim = c(1000,10,10))
pop7.likelihoods <- array(dim = c(1000,10,10))
pop8.likelihoods <- array(dim = c(1000,10,10))
pop9.likelihoods <- array(dim = c(1000,10,10))
pop10.likelihoods <- array(dim = c(1000,10,10))

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

# GAM pop 4
for (h in 1:nrow(adult10.freqs)){
  
  for (j in 1:nrow(Pop4.alleles)){
    
    for (i in 1:ncol(Pop4.alleles)){
      if(Pop4.alleles[j,i] == 2) {
        pop4.likelihoods[j,i,h] <- adult10.freqs[h,i]^2
      } else if (Pop4.alleles[j,i] == 1) {
        pop4.likelihoods[j,i,h] <- 2*(adult10.freqs[h,i] * (1-adult10.freqs[h,i]))
      } else if (Pop4.alleles[j,i] == 0) {
        pop4.likelihoods[j,i,h] <- ( 1-adult10.freqs[h,i])^2 
      } else {
        pop4.likelihoods[j,i,h] <- 1
      }
    }
  }
}


# GAM pop 5
for (h in 1:nrow(adult10.freqs)){
  
  for (j in 1:nrow(Pop5.alleles)){
    
    for (i in 1:ncol(Pop5.alleles)){
      if(Pop5.alleles[j,i] == 2) {
        pop5.likelihoods[j,i,h] <- adult10.freqs[h,i]^2
      } else if (Pop5.alleles[j,i] == 1) {
        pop5.likelihoods[j,i,h] <- 2*(adult10.freqs[h,i] * (1-adult10.freqs[h,i]))
      } else if (Pop5.alleles[j,i] == 0) {
        pop5.likelihoods[j,i,h] <- ( 1-adult10.freqs[h,i])^2 
      } else {
        pop5.likelihoods[j,i,h] <- 1
      }
    }
  }
}

# GAM pop 6
for (h in 1:nrow(adult10.freqs)){
  
  for (j in 1:nrow(Pop6.alleles)){
    
    for (i in 1:ncol(Pop6.alleles)){
      if(Pop6.alleles[j,i] == 2) {
        pop6.likelihoods[j,i,h] <- adult10.freqs[h,i]^2
      } else if (Pop6.alleles[j,i] == 1) {
        pop6.likelihoods[j,i,h] <- 2*(adult10.freqs[h,i] * (1-adult10.freqs[h,i]))
      } else if (Pop6.alleles[j,i] == 0) {
        pop6.likelihoods[j,i,h] <- ( 1-adult10.freqs[h,i])^2 
      } else {
        pop6.likelihoods[j,i,h] <- 1
      }
    }
  }
}

# GAM pop 7
for (h in 1:nrow(adult10.freqs)){
  
  for (j in 1:nrow(Pop7.alleles)){
    
    for (i in 1:ncol(Pop7.alleles)){
      if(Pop7.alleles[j,i] == 2) {
        pop7.likelihoods[j,i,h] <- adult10.freqs[h,i]^2
      } else if (Pop7.alleles[j,i] == 1) {
        pop7.likelihoods[j,i,h] <- 2*(adult10.freqs[h,i] * (1-adult10.freqs[h,i]))
      } else if (Pop7.alleles[j,i] == 0) {
        pop7.likelihoods[j,i,h] <- ( 1-adult10.freqs[h,i])^2 
      } else {
        pop7.likelihoods[j,i,h] <- 1
      }
    }
  }
}

# GAM pop 8
for (h in 1:nrow(adult10.freqs)){
  
  for (j in 1:nrow(Pop8.alleles)){
    
    for (i in 1:ncol(Pop8.alleles)){
      if(Pop8.alleles[j,i] == 2) {
        pop8.likelihoods[j,i,h] <- adult10.freqs[h,i]^2
      } else if (Pop8.alleles[j,i] == 1) {
        pop8.likelihoods[j,i,h] <- 2*(adult10.freqs[h,i] * (1-adult10.freqs[h,i]))
      } else if (Pop8.alleles[j,i] == 0) {
        pop8.likelihoods[j,i,h] <- ( 1-adult10.freqs[h,i])^2 
      } else {
        pop8.likelihoods[j,i,h] <- 1
      }
    }
  }
}

# GAM pop 9
for (h in 1:nrow(adult10.freqs)){
  
  for (j in 1:nrow(Pop9.alleles)){
    
    for (i in 1:ncol(Pop9.alleles)){
      if(Pop9.alleles[j,i] == 2) {
        pop9.likelihoods[j,i,h] <- adult10.freqs[h,i]^2
      } else if (Pop9.alleles[j,i] == 1) {
        pop9.likelihoods[j,i,h] <- 2*(adult10.freqs[h,i] * (1-adult10.freqs[h,i]))
      } else if (Pop9.alleles[j,i] == 0) {
        pop9.likelihoods[j,i,h] <- ( 1-adult10.freqs[h,i])^2 
      } else {
        pop9.likelihoods[j,i,h] <- 1
      }
    }
  }
}

# GAM pop 10 is northern-most
for (h in 1:nrow(adult10.freqs)){
  
  for (j in 1:nrow(Pop10.alleles)){
    
    for (i in 1:ncol(Pop10.alleles)){
      if(Pop10.alleles[j,i] == 2) {
        pop10.likelihoods[j,i,h] <- adult10.freqs[h,i]^2
      } else if (Pop10.alleles[j,i] == 1) {
        pop10.likelihoods[j,i,h] <- 2*(adult10.freqs[h,i] * (1-adult10.freqs[h,i]))
      } else if (Pop10.alleles[j,i] == 0) {
        pop10.likelihoods[j,i,h] <- ( 1-adult10.freqs[h,i])^2 
      } else {
        pop10.likelihoods[j,i,h] <- 1
      }
    }
  }
}


#### Multiply rows together and then assign each individual to a GAM population ####
# Pop 1
Pop1.ass <- list()
for (l in 1:10){
Pop1.ass[[l]] <- apply(pop1.likelihoods[,,l], FUN = prod, MARGIN = 1, na.rm = TRUE)
  }

Pop1.ass.df <- data.frame(matrix(unlist(Pop1.ass), nrow=nrow(pop1.likelihoods), byrow=F))
colnames(Pop1.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop1.most.like <- colnames(Pop1.ass.df)[apply(Pop1.ass.df,1, which.max)]

# Pop 2
Pop2.ass <- list()
for (l in 1:10){
  Pop2.ass[[l]] <- apply(pop2.likelihoods[,,l], FUN = prod, MARGIN = 1, na.rm = TRUE)
}
  
Pop2.ass.df <- data.frame(matrix(unlist(Pop2.ass), nrow=nrow(pop2.likelihoods), byrow=F))
colnames(Pop2.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop2.most.like <- colnames(Pop2.ass.df)[apply(Pop2.ass.df,1, which.max)]

# Pop 3
Pop3.ass <- list()
for (l in 1:10){
  Pop3.ass[[l]] <- apply(pop3.likelihoods[,,l], FUN = prod, MARGIN = 1, na.rm = TRUE)
}

Pop3.ass.df <- data.frame(matrix(unlist(Pop3.ass), nrow=nrow(pop3.likelihoods), byrow=F))
colnames(Pop3.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop3.most.like <- colnames(Pop3.ass.df)[apply(Pop3.ass.df,1, which.max)]

# Pop 4
Pop4.ass <- list()
for (l in 1:10){
  Pop4.ass[[l]] <- apply(pop4.likelihoods[,,l], FUN = prod, MARGIN = 1, na.rm = TRUE)
}

Pop4.ass.df <- data.frame(matrix(unlist(Pop4.ass), nrow=nrow(pop4.likelihoods), byrow=F))
colnames(Pop4.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop4.most.like <- colnames(Pop4.ass.df)[apply(Pop4.ass.df,1, which.max)]

# Pop 5
Pop5.ass <- list()
for (l in 1:10){
  Pop5.ass[[l]] <- apply(pop5.likelihoods[,,l], FUN = prod, MARGIN = 1, na.rm = TRUE)
}

Pop5.ass.df <- data.frame(matrix(unlist(Pop5.ass), nrow=nrow(pop5.likelihoods), byrow=F))
colnames(Pop5.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop5.most.like <- colnames(Pop5.ass.df)[apply(Pop5.ass.df,1, which.max)]

# Pop 6
Pop6.ass <- list()
for (l in 1:10){
  Pop6.ass[[l]] <- apply(pop6.likelihoods[,,l], FUN = prod, MARGIN = 1, na.rm = TRUE)
}

Pop6.ass.df <- data.frame(matrix(unlist(Pop6.ass), nrow=nrow(pop6.likelihoods), byrow=F))
colnames(Pop6.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop6.most.like <- colnames(Pop6.ass.df)[apply(Pop6.ass.df,1, which.max)]

# Pop 7
Pop7.ass <- list()
for (l in 1:10){
  Pop7.ass[[l]] <- apply(pop7.likelihoods[,,l], FUN = prod, MARGIN = 1, na.rm = TRUE)
}

Pop7.ass.df <- data.frame(matrix(unlist(Pop7.ass), nrow=nrow(pop7.likelihoods), byrow=F))
colnames(Pop7.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop7.most.like <- colnames(Pop7.ass.df)[apply(Pop7.ass.df,1, which.max)]

# Pop 8
Pop8.ass <- list()
for (l in 1:10){
  Pop8.ass[[l]] <- apply(pop8.likelihoods[,,l], FUN = prod, MARGIN = 1, na.rm = TRUE)
}

Pop8.ass.df <- data.frame(matrix(unlist(Pop8.ass), nrow=nrow(pop8.likelihoods), byrow=F))
colnames(Pop8.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop8.most.like <- colnames(Pop8.ass.df)[apply(Pop8.ass.df,1, which.max)]

# Pop 9
Pop9.ass <- list()
for (l in 1:10){
  Pop9.ass[[l]] <- apply(pop9.likelihoods[,,l], FUN = prod, MARGIN = 1, na.rm = TRUE)
}

Pop9.ass.df <- data.frame(matrix(unlist(Pop9.ass), nrow=nrow(pop9.likelihoods), byrow=F))
colnames(Pop9.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop9.most.like <- colnames(Pop9.ass.df)[apply(Pop9.ass.df,1, which.max)]

# Pop 10
Pop10.ass <- list()
for (l in 1:10){
  Pop10.ass[[l]] <- apply(pop10.likelihoods[,,l], FUN = prod, MARGIN = 1, na.rm = TRUE)
}

Pop10.ass.df <- data.frame(matrix(unlist(Pop10.ass), nrow=nrow(pop10.likelihoods), byrow=F))
colnames(Pop10.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop10.most.like <- colnames(Pop10.ass.df)[apply(Pop10.ass.df,1, which.max)]


#### Rbind & organize the assignment results for individuals simulated from each of 10 GAM populations ####
# Row as reference? Column as predicted?
mat <- rbind(table(Pop10.most.like)/1000, table(Pop9.most.like)/1000, table(Pop8.most.like)/1000, table(Pop7.most.like)/1000, table(Pop6.most.like)/1000, table(Pop5.most.like)/1000, table(Pop4.most.like)/1000, table(Pop3.most.like)/1000, table(Pop2.most.like)/1000, table(Pop1.most.like)/1000) # Each row was drawn from a known population. This should be the reference. Columns are predicted.
rownames(mat) <- c('Pop A', 'Pop B', 'Pop C', 'Pop D', 'Pop E', 'Pop F', 'Pop G', 'Pop H', 'Pop I', 'Pop J')

#### Plot confusion matrix ####
# image(1:ncol(mat), 1:nrow(mat), t(mat), col = rainbow(5), xlab = 'Predicted population', ylab = 'True population', axes = FALSE)
# axis(1, 1:ncol(mat), colnames(mat))
# axis(2, 1:nrow(mat), rownames(mat))

# Nice plot
# png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/assignment_confusion_matrix.png", width=7, height=5, res=300, units="in")

par(
  mar=c(5, 4, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

mat2 <- data.matrix(mat)
# min <- min(mat, na.rm = TRUE)
# max <- max(mat,na.rm = TRUE)

rbPal <- colorRampPalette(c('blue', 'red'))
nocol <- 6
color <- rbPal(nocol)
ColorLevels <- seq(0, 0.5, length=length(color)) # this makes the color bar nicer later on
# ColorLevels <- seq(min, max, length=length(color))

# Set layout.  We are going to include a colorbar next to plot.
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1),
       heights=c(1,1))

#plotting margins.  These seem to work well for me.
par(mar = c(5,5,1.5,1), font = 2)

# Plot it up!
image(1:ncol(mat2), 1:nrow(mat2), t(mat2),
      col=color, xlab="Predicted population", ylab="True population",
      axes=FALSE, zlim = c(0, 0.5),
      main= NA, xlim = c(0.5, 11), ylim = c(0.5, 11))
# image(1:ncol(mat2), 1:nrow(mat2), t(mat2),
#       col=color, xlab="Predicted population", ylab="True population",
#       axes=FALSE, zlim = c(min, max),
#       main= NA, xlim = c(0.5, 11), ylim = c(0.5, 11))
abline(0,1)

# Now annotate the plot
axis(side = 1, at=seq(1,10,1), labels=LETTERS[1:10],
     cex.axis=1.0)
axis(side = 2, at=seq(1,10,1), labels=LETTERS[1:10], las= 1,
     cex.axis=1)

# Add colorbar to second plot region
par(mar = c(3.5,2.5,3.5,2))
image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=color,xlab="",ylab="",xaxt="n", las = 1)
mtext('% Assignment', side = 3, cex = 0.7, line = 0.5, adj = 1)

dev.off()

################################################################################################################################################################

#####################################################################
#### Simulations to generate a confusion matrix of n individuals ####
#####################################################################
# Set n, the number of individuals in the group
n <- 20

# Sample one allele per locus N times (depth of array) from each of 10 GAM populations, and repeat this 1000 times
Pop1.alleles <- array(dim = c(1000,10,n)) # 1//J is southern-most, 10//A is northern-most
for (i in 1:n){
  Pop1.alleles[,,i] <- sapply(adult10.freqs[1,],function(z){rbinom(1000,2,z)})
}

# Do this for the rest of the populations too
Pop2.alleles <- array(dim = c(1000,10,n)) # 1//J is southern-most, 10//A is northern-most
for (i in 1:n){
  Pop2.alleles[,,i] <- sapply(adult10.freqs[2,],function(z){rbinom(1000,2,z)})
}

Pop3.alleles <- array(dim = c(1000,10,n)) # 1//J is southern-most, 10//A is northern-most
for (i in 1:n){
  Pop3.alleles[,,i] <- sapply(adult10.freqs[3,],function(z){rbinom(1000,2,z)})
}

Pop4.alleles <- array(dim = c(1000,10,n)) # 1//J is southern-most, 10//A is northern-most
for (i in 1:n){
  Pop4.alleles[,,i] <- sapply(adult10.freqs[4,],function(z){rbinom(1000,2,z)})
}

Pop5.alleles <- array(dim = c(1000,10,n)) # 1//J is southern-most, 10//A is northern-most
for (i in 1:n){
  Pop5.alleles[,,i] <- sapply(adult10.freqs[5,],function(z){rbinom(1000,2,z)})
}

Pop6.alleles <- array(dim = c(1000,10,n)) # 1//J is southern-most, 10//A is northern-most
for (i in 1:n){
  Pop6.alleles[,,i] <- sapply(adult10.freqs[6,],function(z){rbinom(1000,2,z)})
}

Pop7.alleles <- array(dim = c(1000,10,n)) # 1//J is southern-most, 10//A is northern-most
for (i in 1:n){
  Pop7.alleles[,,i] <- sapply(adult10.freqs[7,],function(z){rbinom(1000,2,z)})
}

Pop8.alleles <- array(dim = c(1000,10,n)) # 1//J is southern-most, 10//A is northern-most
for (i in 1:n){
  Pop8.alleles[,,i] <- sapply(adult10.freqs[8,],function(z){rbinom(1000,2,z)})
}

Pop9.alleles <- array(dim = c(1000,10,n)) # 1//J is southern-most, 10//A is northern-most
for (i in 1:n){
  Pop9.alleles[,,i] <- sapply(adult10.freqs[9,],function(z){rbinom(1000,2,z)})
}

Pop10.alleles <- array(dim = c(1000,10,n)) # 1//J is southern-most, 10//A is northern-most
for (i in 1:n){
  Pop10.alleles[,,i] <- sapply(adult10.freqs[10,],function(z){rbinom(1000,2,z)})
}

#### Population likelihoods for 10 groups ####
#### Alleles drawn from Pop1, genotype likelihoods calculated for all 10 GAM populations ####
Pop1.likelihoods1 <- array(dim = c(1000,10,n))
Pop1.likelihoods2 <- array(dim = c(1000,10,n))
Pop1.likelihoods3 <- array(dim = c(1000,10,n))
Pop1.likelihoods4 <- array(dim = c(1000,10,n))
Pop1.likelihoods5 <- array(dim = c(1000,10,n))
Pop1.likelihoods6 <- array(dim = c(1000,10,n))
Pop1.likelihoods7 <- array(dim = c(1000,10,n))
Pop1.likelihoods8 <- array(dim = c(1000,10,n))
Pop1.likelihoods9 <- array(dim = c(1000,10,n))
Pop1.likelihoods10 <- array(dim = c(1000,10,n))

# For larvae simulated using GAM pop 1 allele freqs. Then geno likelihoods calculated for all 10 GAM pops
# Calculations for genotype likelihoods from 10 GAM populations for groups of n individuals sampled from GAM 1
for (h in 1:n){
  
  for (j in 1:nrow(Pop1.alleles)){
    
    for (i in 1:ncol(Pop1.alleles)){
      if(Pop1.alleles[j,i,h] == 2) {
        Pop1.likelihoods1[j,i,h] <- adult10.freqs[1,i]^2
      } else if (Pop1.alleles[j,i,h] == 1) {
        Pop1.likelihoods1[j,i,h] <- 2*(adult10.freqs[1,i] * (1-adult10.freqs[1,i]))
      } else if (Pop1.alleles[j,i,h] == 0) {
        Pop1.likelihoods1[j,i,h] <- ( 1-adult10.freqs[1,i])^2 
      } else {
        Pop1.likelihoods1[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop1.alleles)){
    
    for (i in 1:ncol(Pop1.alleles)){
      if(Pop1.alleles[j,i,h] == 2) {
        Pop1.likelihoods2[j,i,h] <- adult10.freqs[2,i]^2
      } else if (Pop1.alleles[j,i,h] == 1) {
        Pop1.likelihoods2[j,i,h] <- 2*(adult10.freqs[2,i] * (1-adult10.freqs[2,i]))
      } else if (Pop1.alleles[j,i,h] == 0) {
        Pop1.likelihoods2[j,i,h] <- ( 1-adult10.freqs[2,i])^2 
      } else {
        Pop1.likelihoods2[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop1.alleles)){
    
    for (i in 1:ncol(Pop1.alleles)){
      if(Pop1.alleles[j,i,h] == 2) {
        Pop1.likelihoods3[j,i,h] <- adult10.freqs[3,i]^2
      } else if (Pop1.alleles[j,i,h] == 1) {
        Pop1.likelihoods3[j,i,h] <- 2*(adult10.freqs[3,i] * (1-adult10.freqs[3,i]))
      } else if (Pop1.alleles[j,i,h] == 0) {
        Pop1.likelihoods3[j,i,h] <- ( 1-adult10.freqs[3,i])^2 
      } else {
        Pop1.likelihoods3[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop1.alleles)){
    
    for (i in 1:ncol(Pop1.alleles)){
      if(Pop1.alleles[j,i,h] == 2) {
        Pop1.likelihoods4[j,i,h] <- adult10.freqs[4,i]^2
      } else if (Pop1.alleles[j,i,h] == 1) {
        Pop1.likelihoods4[j,i,h] <- 2*(adult10.freqs[4,i] * (1-adult10.freqs[4,i]))
      } else if (Pop1.alleles[j,i,h] == 0) {
        Pop1.likelihoods4[j,i,h] <- ( 1-adult10.freqs[4,i])^2 
      } else {
        Pop1.likelihoods4[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop1.alleles)){
    
    for (i in 1:ncol(Pop1.alleles)){
      if(Pop1.alleles[j,i,h] == 2) {
        Pop1.likelihoods5[j,i,h] <- adult10.freqs[5,i]^2
      } else if (Pop1.alleles[j,i,h] == 1) {
        Pop1.likelihoods5[j,i,h] <- 2*(adult10.freqs[5,i] * (1-adult10.freqs[5,i]))
      } else if (Pop1.alleles[j,i,h] == 0) {
        Pop1.likelihoods5[j,i,h] <- ( 1-adult10.freqs[5,i])^2 
      } else {
        Pop1.likelihoods5[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop1.alleles)){
    
    for (i in 1:ncol(Pop1.alleles)){
      if(Pop1.alleles[j,i,h] == 2) {
        Pop1.likelihoods6[j,i,h] <- adult10.freqs[6,i]^2
      } else if (Pop1.alleles[j,i,h] == 1) {
        Pop1.likelihoods6[j,i,h] <- 2*(adult10.freqs[6,i] * (1-adult10.freqs[6,i]))
      } else if (Pop1.alleles[j,i,h] == 0) {
        Pop1.likelihoods6[j,i,h] <- ( 1-adult10.freqs[6,i])^2 
      } else {
        Pop1.likelihoods6[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop1.alleles)){
    
    for (i in 1:ncol(Pop1.alleles)){
      if(Pop1.alleles[j,i,h] == 2) {
        Pop1.likelihoods7[j,i,h] <- adult10.freqs[7,i]^2
      } else if (Pop1.alleles[j,i,h] == 1) {
        Pop1.likelihoods7[j,i,h] <- 2*(adult10.freqs[7,i] * (1-adult10.freqs[7,i]))
      } else if (Pop1.alleles[j,i,h] == 0) {
        Pop1.likelihoods7[j,i,h] <- ( 1-adult10.freqs[7,i])^2 
      } else {
        Pop1.likelihoods7[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop1.alleles)){
    
    for (i in 1:ncol(Pop1.alleles)){
      if(Pop1.alleles[j,i,h] == 2) {
        Pop1.likelihoods8[j,i,h] <- adult10.freqs[8,i]^2
      } else if (Pop1.alleles[j,i,h] == 1) {
        Pop1.likelihoods8[j,i,h] <- 2*(adult10.freqs[8,i] * (1-adult10.freqs[8,i]))
      } else if (Pop1.alleles[j,i,h] == 0) {
        Pop1.likelihoods8[j,i,h] <- ( 1-adult10.freqs[8,i])^2 
      } else {
        Pop1.likelihoods8[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop1.alleles)){
    
    for (i in 1:ncol(Pop1.alleles)){
      if(Pop1.alleles[j,i,h] == 2) {
        Pop1.likelihoods9[j,i,h] <- adult10.freqs[9,i]^2
      } else if (Pop1.alleles[j,i,h] == 1) {
        Pop1.likelihoods9[j,i,h] <- 2*(adult10.freqs[9,i] * (1-adult10.freqs[9,i]))
      } else if (Pop1.alleles[j,i,h] == 0) {
        Pop1.likelihoods9[j,i,h] <- ( 1-adult10.freqs[9,i])^2 
      } else {
        Pop1.likelihoods9[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop1.alleles)){
    
    for (i in 1:ncol(Pop1.alleles)){
      if(Pop1.alleles[j,i,h] == 2) {
        Pop1.likelihoods10[j,i,h] <- adult10.freqs[10,i]^2
      } else if (Pop1.alleles[j,i,h] == 1) {
        Pop1.likelihoods10[j,i,h] <- 2*(adult10.freqs[10,i] * (1-adult10.freqs[10,i]))
      } else if (Pop1.alleles[j,i,h] == 0) {
        Pop1.likelihoods10[j,i,h] <- ( 1-adult10.freqs[10,i])^2 
      } else {
        Pop1.likelihoods10[j,i,h] <- 1
      }
    }
  }
}

# Multiply everything together from groups of n individuals originating from GAM 1
Pop1.geno1 <- apply(Pop1.likelihoods1, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop1.geno2 <- apply(Pop1.likelihoods2, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop1.geno3 <- apply(Pop1.likelihoods3, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop1.geno4 <- apply(Pop1.likelihoods4, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop1.geno5 <- apply(Pop1.likelihoods5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop1.geno6 <- apply(Pop1.likelihoods6, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop1.geno7 <- apply(Pop1.likelihoods7, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop1.geno8 <- apply(Pop1.likelihoods8, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop1.geno9 <- apply(Pop1.likelihoods9, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop1.geno10 <- apply(Pop1.likelihoods10, FUN = prod, MARGIN = 1, na.rm = TRUE)

# Assign groups of n individuals originating from GAM 1 to most likely pop
Pop1.ass.df <- data.frame(cbind(Pop1.geno1, Pop1.geno2, Pop1.geno3, Pop1.geno4, Pop1.geno5, Pop1.geno6, Pop1.geno7, Pop1.geno8, Pop1.geno9, Pop1.geno10))
colnames(Pop1.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop1.most.like <- colnames(Pop1.ass.df)[apply(Pop1.ass.df,1, which.max)]

#### Alleles drawn from Pop2, genotype likelihoods calculated for all 10 GAM populations ####
Pop2.likelihoods1 <- array(dim = c(1000,10,n))
Pop2.likelihoods2 <- array(dim = c(1000,10,n))
Pop2.likelihoods3 <- array(dim = c(1000,10,n))
Pop2.likelihoods4 <- array(dim = c(1000,10,n))
Pop2.likelihoods5 <- array(dim = c(1000,10,n))
Pop2.likelihoods6 <- array(dim = c(1000,10,n))
Pop2.likelihoods7 <- array(dim = c(1000,10,n))
Pop2.likelihoods8 <- array(dim = c(1000,10,n))
Pop2.likelihoods9 <- array(dim = c(1000,10,n))
Pop2.likelihoods10 <- array(dim = c(1000,10,n))

# For larvae simulated using GAM pop 2 allele freqs. Then geno likelihoods calculated for all 10 GAM pops
# Calculations for genotype likelihoods from 10 GAM populations for groups of n individuals sampled from GAM 2
for (h in 1:n){
  
  for (j in 1:nrow(Pop2.alleles)){
    
    for (i in 1:ncol(Pop2.alleles)){
      if(Pop2.alleles[j,i,h] == 2) {
        Pop2.likelihoods1[j,i,h] <- adult10.freqs[1,i]^2
      } else if (Pop2.alleles[j,i,h] == 1) {
        Pop2.likelihoods1[j,i,h] <- 2*(adult10.freqs[1,i] * (1-adult10.freqs[1,i]))
      } else if (Pop2.alleles[j,i,h] == 0) {
        Pop2.likelihoods1[j,i,h] <- ( 1-adult10.freqs[1,i])^2 
      } else {
        Pop2.likelihoods1[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop2.alleles)){
    
    for (i in 1:ncol(Pop2.alleles)){
      if(Pop2.alleles[j,i,h] == 2) {
        Pop2.likelihoods2[j,i,h] <- adult10.freqs[2,i]^2
      } else if (Pop2.alleles[j,i,h] == 1) {
        Pop2.likelihoods2[j,i,h] <- 2*(adult10.freqs[2,i] * (1-adult10.freqs[2,i]))
      } else if (Pop2.alleles[j,i,h] == 0) {
        Pop2.likelihoods2[j,i,h] <- ( 1-adult10.freqs[2,i])^2 
      } else {
        Pop2.likelihoods2[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop2.alleles)){
    
    for (i in 1:ncol(Pop2.alleles)){
      if(Pop2.alleles[j,i,h] == 2) {
        Pop2.likelihoods3[j,i,h] <- adult10.freqs[3,i]^2
      } else if (Pop2.alleles[j,i,h] == 1) {
        Pop2.likelihoods3[j,i,h] <- 2*(adult10.freqs[3,i] * (1-adult10.freqs[3,i]))
      } else if (Pop2.alleles[j,i,h] == 0) {
        Pop2.likelihoods3[j,i,h] <- ( 1-adult10.freqs[3,i])^2 
      } else {
        Pop2.likelihoods3[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop2.alleles)){
    
    for (i in 1:ncol(Pop2.alleles)){
      if(Pop2.alleles[j,i,h] == 2) {
        Pop2.likelihoods4[j,i,h] <- adult10.freqs[4,i]^2
      } else if (Pop2.alleles[j,i,h] == 1) {
        Pop2.likelihoods4[j,i,h] <- 2*(adult10.freqs[4,i] * (1-adult10.freqs[4,i]))
      } else if (Pop2.alleles[j,i,h] == 0) {
        Pop2.likelihoods4[j,i,h] <- ( 1-adult10.freqs[4,i])^2 
      } else {
        Pop2.likelihoods4[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop2.alleles)){
    
    for (i in 1:ncol(Pop2.alleles)){
      if(Pop2.alleles[j,i,h] == 2) {
        Pop2.likelihoods5[j,i,h] <- adult10.freqs[5,i]^2
      } else if (Pop2.alleles[j,i,h] == 1) {
        Pop2.likelihoods5[j,i,h] <- 2*(adult10.freqs[5,i] * (1-adult10.freqs[5,i]))
      } else if (Pop2.alleles[j,i,h] == 0) {
        Pop2.likelihoods5[j,i,h] <- ( 1-adult10.freqs[5,i])^2 
      } else {
        Pop2.likelihoods5[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop2.alleles)){
    
    for (i in 1:ncol(Pop2.alleles)){
      if(Pop2.alleles[j,i,h] == 2) {
        Pop2.likelihoods6[j,i,h] <- adult10.freqs[6,i]^2
      } else if (Pop2.alleles[j,i,h] == 1) {
        Pop2.likelihoods6[j,i,h] <- 2*(adult10.freqs[6,i] * (1-adult10.freqs[6,i]))
      } else if (Pop2.alleles[j,i,h] == 0) {
        Pop2.likelihoods6[j,i,h] <- ( 1-adult10.freqs[6,i])^2 
      } else {
        Pop2.likelihoods6[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop2.alleles)){
    
    for (i in 1:ncol(Pop2.alleles)){
      if(Pop2.alleles[j,i,h] == 2) {
        Pop2.likelihoods7[j,i,h] <- adult10.freqs[7,i]^2
      } else if (Pop2.alleles[j,i,h] == 1) {
        Pop2.likelihoods7[j,i,h] <- 2*(adult10.freqs[7,i] * (1-adult10.freqs[7,i]))
      } else if (Pop2.alleles[j,i,h] == 0) {
        Pop2.likelihoods7[j,i,h] <- ( 1-adult10.freqs[7,i])^2 
      } else {
        Pop2.likelihoods7[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop2.alleles)){
    
    for (i in 1:ncol(Pop2.alleles)){
      if(Pop2.alleles[j,i,h] == 2) {
        Pop2.likelihoods8[j,i,h] <- adult10.freqs[8,i]^2
      } else if (Pop2.alleles[j,i,h] == 1) {
        Pop2.likelihoods8[j,i,h] <- 2*(adult10.freqs[8,i] * (1-adult10.freqs[8,i]))
      } else if (Pop2.alleles[j,i,h] == 0) {
        Pop2.likelihoods8[j,i,h] <- ( 1-adult10.freqs[8,i])^2 
      } else {
        Pop2.likelihoods8[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop2.alleles)){
    
    for (i in 1:ncol(Pop2.alleles)){
      if(Pop2.alleles[j,i,h] == 2) {
        Pop2.likelihoods9[j,i,h] <- adult10.freqs[9,i]^2
      } else if (Pop2.alleles[j,i,h] == 1) {
        Pop2.likelihoods9[j,i,h] <- 2*(adult10.freqs[9,i] * (1-adult10.freqs[9,i]))
      } else if (Pop2.alleles[j,i,h] == 0) {
        Pop2.likelihoods9[j,i,h] <- ( 1-adult10.freqs[9,i])^2 
      } else {
        Pop2.likelihoods9[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop2.alleles)){
    
    for (i in 1:ncol(Pop2.alleles)){
      if(Pop2.alleles[j,i,h] == 2) {
        Pop2.likelihoods10[j,i,h] <- adult10.freqs[10,i]^2
      } else if (Pop2.alleles[j,i,h] == 1) {
        Pop2.likelihoods10[j,i,h] <- 2*(adult10.freqs[10,i] * (1-adult10.freqs[10,i]))
      } else if (Pop2.alleles[j,i,h] == 0) {
        Pop2.likelihoods10[j,i,h] <- ( 1-adult10.freqs[10,i])^2 
      } else {
        Pop2.likelihoods10[j,i,h] <- 1
      }
    }
  }
}

# Multiply everything together from groups of n individuals originating from GAM 1
Pop2.geno1 <- apply(Pop2.likelihoods1, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop2.geno2 <- apply(Pop2.likelihoods2, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop2.geno3 <- apply(Pop2.likelihoods3, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop2.geno4 <- apply(Pop2.likelihoods4, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop2.geno5 <- apply(Pop2.likelihoods5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop2.geno6 <- apply(Pop2.likelihoods6, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop2.geno7 <- apply(Pop2.likelihoods7, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop2.geno8 <- apply(Pop2.likelihoods8, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop2.geno9 <- apply(Pop2.likelihoods9, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop2.geno10 <- apply(Pop2.likelihoods10, FUN = prod, MARGIN = 1, na.rm = TRUE)

# Assign groups of n individuals originating from GAM 1 to most likely pop
Pop2.ass.df <- data.frame(cbind(Pop2.geno1, Pop2.geno2, Pop2.geno3, Pop2.geno4, Pop2.geno5, Pop2.geno6, Pop2.geno7, Pop2.geno8, Pop2.geno9, Pop2.geno10))
colnames(Pop2.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop2.most.like <- colnames(Pop2.ass.df)[apply(Pop2.ass.df,1, which.max)]

#### Alleles drawn from Pop3, genotype likelihoods calculated for all 10 GAM populations ####
Pop3.likelihoods1 <- array(dim = c(1000,10,n))
Pop3.likelihoods2 <- array(dim = c(1000,10,n))
Pop3.likelihoods3 <- array(dim = c(1000,10,n))
Pop3.likelihoods4 <- array(dim = c(1000,10,n))
Pop3.likelihoods5 <- array(dim = c(1000,10,n))
Pop3.likelihoods6 <- array(dim = c(1000,10,n))
Pop3.likelihoods7 <- array(dim = c(1000,10,n))
Pop3.likelihoods8 <- array(dim = c(1000,10,n))
Pop3.likelihoods9 <- array(dim = c(1000,10,n))
Pop3.likelihoods10 <- array(dim = c(1000,10,n))

# For larvae simulated using GAM pop 3 allele freqs. Then geno likelihoods calculated for all 10 GAM pops
# Calculations for genotype likelihoods from 10 GAM populations for groups of n individuals sampled from GAM 3
for (h in 1:n){
  
  for (j in 1:nrow(Pop3.alleles)){
    
    for (i in 1:ncol(Pop3.alleles)){
      if(Pop3.alleles[j,i,h] == 2) {
        Pop3.likelihoods1[j,i,h] <- adult10.freqs[1,i]^2
      } else if (Pop3.alleles[j,i,h] == 1) {
        Pop3.likelihoods1[j,i,h] <- 2*(adult10.freqs[1,i] * (1-adult10.freqs[1,i]))
      } else if (Pop3.alleles[j,i,h] == 0) {
        Pop3.likelihoods1[j,i,h] <- ( 1-adult10.freqs[1,i])^2 
      } else {
        Pop3.likelihoods1[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop3.alleles)){
    
    for (i in 1:ncol(Pop3.alleles)){
      if(Pop3.alleles[j,i,h] == 2) {
        Pop3.likelihoods2[j,i,h] <- adult10.freqs[2,i]^2
      } else if (Pop3.alleles[j,i,h] == 1) {
        Pop3.likelihoods2[j,i,h] <- 2*(adult10.freqs[2,i] * (1-adult10.freqs[2,i]))
      } else if (Pop3.alleles[j,i,h] == 0) {
        Pop3.likelihoods2[j,i,h] <- ( 1-adult10.freqs[2,i])^2 
      } else {
        Pop3.likelihoods2[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop3.alleles)){
    
    for (i in 1:ncol(Pop3.alleles)){
      if(Pop3.alleles[j,i,h] == 2) {
        Pop3.likelihoods3[j,i,h] <- adult10.freqs[3,i]^2
      } else if (Pop3.alleles[j,i,h] == 1) {
        Pop3.likelihoods3[j,i,h] <- 2*(adult10.freqs[3,i] * (1-adult10.freqs[3,i]))
      } else if (Pop3.alleles[j,i,h] == 0) {
        Pop3.likelihoods3[j,i,h] <- ( 1-adult10.freqs[3,i])^2 
      } else {
        Pop3.likelihoods3[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop3.alleles)){
    
    for (i in 1:ncol(Pop3.alleles)){
      if(Pop3.alleles[j,i,h] == 2) {
        Pop3.likelihoods4[j,i,h] <- adult10.freqs[4,i]^2
      } else if (Pop3.alleles[j,i,h] == 1) {
        Pop3.likelihoods4[j,i,h] <- 2*(adult10.freqs[4,i] * (1-adult10.freqs[4,i]))
      } else if (Pop3.alleles[j,i,h] == 0) {
        Pop3.likelihoods4[j,i,h] <- ( 1-adult10.freqs[4,i])^2 
      } else {
        Pop3.likelihoods4[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop3.alleles)){
    
    for (i in 1:ncol(Pop3.alleles)){
      if(Pop3.alleles[j,i,h] == 2) {
        Pop3.likelihoods5[j,i,h] <- adult10.freqs[5,i]^2
      } else if (Pop3.alleles[j,i,h] == 1) {
        Pop3.likelihoods5[j,i,h] <- 2*(adult10.freqs[5,i] * (1-adult10.freqs[5,i]))
      } else if (Pop3.alleles[j,i,h] == 0) {
        Pop3.likelihoods5[j,i,h] <- ( 1-adult10.freqs[5,i])^2 
      } else {
        Pop3.likelihoods5[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop3.alleles)){
    
    for (i in 1:ncol(Pop3.alleles)){
      if(Pop3.alleles[j,i,h] == 2) {
        Pop3.likelihoods6[j,i,h] <- adult10.freqs[6,i]^2
      } else if (Pop3.alleles[j,i,h] == 1) {
        Pop3.likelihoods6[j,i,h] <- 2*(adult10.freqs[6,i] * (1-adult10.freqs[6,i]))
      } else if (Pop3.alleles[j,i,h] == 0) {
        Pop3.likelihoods6[j,i,h] <- ( 1-adult10.freqs[6,i])^2 
      } else {
        Pop3.likelihoods6[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop3.alleles)){
    
    for (i in 1:ncol(Pop3.alleles)){
      if(Pop3.alleles[j,i,h] == 2) {
        Pop3.likelihoods7[j,i,h] <- adult10.freqs[7,i]^2
      } else if (Pop3.alleles[j,i,h] == 1) {
        Pop3.likelihoods7[j,i,h] <- 2*(adult10.freqs[7,i] * (1-adult10.freqs[7,i]))
      } else if (Pop3.alleles[j,i,h] == 0) {
        Pop3.likelihoods7[j,i,h] <- ( 1-adult10.freqs[7,i])^2 
      } else {
        Pop3.likelihoods7[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop3.alleles)){
    
    for (i in 1:ncol(Pop3.alleles)){
      if(Pop3.alleles[j,i,h] == 2) {
        Pop3.likelihoods8[j,i,h] <- adult10.freqs[8,i]^2
      } else if (Pop3.alleles[j,i,h] == 1) {
        Pop3.likelihoods8[j,i,h] <- 2*(adult10.freqs[8,i] * (1-adult10.freqs[8,i]))
      } else if (Pop3.alleles[j,i,h] == 0) {
        Pop3.likelihoods8[j,i,h] <- ( 1-adult10.freqs[8,i])^2 
      } else {
        Pop3.likelihoods8[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop3.alleles)){
    
    for (i in 1:ncol(Pop3.alleles)){
      if(Pop3.alleles[j,i,h] == 2) {
        Pop3.likelihoods9[j,i,h] <- adult10.freqs[9,i]^2
      } else if (Pop3.alleles[j,i,h] == 1) {
        Pop3.likelihoods9[j,i,h] <- 2*(adult10.freqs[9,i] * (1-adult10.freqs[9,i]))
      } else if (Pop3.alleles[j,i,h] == 0) {
        Pop3.likelihoods9[j,i,h] <- ( 1-adult10.freqs[9,i])^2 
      } else {
        Pop3.likelihoods9[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop3.alleles)){
    
    for (i in 1:ncol(Pop3.alleles)){
      if(Pop3.alleles[j,i,h] == 2) {
        Pop3.likelihoods10[j,i,h] <- adult10.freqs[10,i]^2
      } else if (Pop3.alleles[j,i,h] == 1) {
        Pop3.likelihoods10[j,i,h] <- 2*(adult10.freqs[10,i] * (1-adult10.freqs[10,i]))
      } else if (Pop3.alleles[j,i,h] == 0) {
        Pop3.likelihoods10[j,i,h] <- ( 1-adult10.freqs[10,i])^2 
      } else {
        Pop3.likelihoods10[j,i,h] <- 1
      }
    }
  }
}

# Multiply everything together from groups of n individuals originating from GAM 1
Pop3.geno1 <- apply(Pop3.likelihoods1, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop3.geno2 <- apply(Pop3.likelihoods2, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop3.geno3 <- apply(Pop3.likelihoods3, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop3.geno4 <- apply(Pop3.likelihoods4, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop3.geno5 <- apply(Pop3.likelihoods5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop3.geno6 <- apply(Pop3.likelihoods6, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop3.geno7 <- apply(Pop3.likelihoods7, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop3.geno8 <- apply(Pop3.likelihoods8, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop3.geno9 <- apply(Pop3.likelihoods9, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop3.geno10 <- apply(Pop3.likelihoods10, FUN = prod, MARGIN = 1, na.rm = TRUE)

# Assign groups of n individuals originating from GAM 1 to most likely pop
Pop3.ass.df <- data.frame(cbind(Pop3.geno1, Pop3.geno2, Pop3.geno3, Pop3.geno4, Pop3.geno5, Pop3.geno6, Pop3.geno7, Pop3.geno8, Pop3.geno9, Pop3.geno10))
colnames(Pop3.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop3.most.like <- colnames(Pop3.ass.df)[apply(Pop3.ass.df,1, which.max)]

#### Alleles drawn from Pop4, genotype likelihoods calculated for all 10 GAM populations ####
Pop4.likelihoods1 <- array(dim = c(1000,10,n))
Pop4.likelihoods2 <- array(dim = c(1000,10,n))
Pop4.likelihoods3 <- array(dim = c(1000,10,n))
Pop4.likelihoods4 <- array(dim = c(1000,10,n))
Pop4.likelihoods5 <- array(dim = c(1000,10,n))
Pop4.likelihoods6 <- array(dim = c(1000,10,n))
Pop4.likelihoods7 <- array(dim = c(1000,10,n))
Pop4.likelihoods8 <- array(dim = c(1000,10,n))
Pop4.likelihoods9 <- array(dim = c(1000,10,n))
Pop4.likelihoods10 <- array(dim = c(1000,10,n))

# For larvae simulated using GAM pop 4 allele freqs. Then geno likelihoods calculated for all 10 GAM pops
# Calculations for genotype likelihoods from 10 GAM populations for groups of n individuals sampled from GAM 4
for (h in 1:n){
  
  for (j in 1:nrow(Pop4.alleles)){
    
    for (i in 1:ncol(Pop4.alleles)){
      if(Pop4.alleles[j,i,h] == 2) {
        Pop4.likelihoods1[j,i,h] <- adult10.freqs[1,i]^2
      } else if (Pop4.alleles[j,i,h] == 1) {
        Pop4.likelihoods1[j,i,h] <- 2*(adult10.freqs[1,i] * (1-adult10.freqs[1,i]))
      } else if (Pop4.alleles[j,i,h] == 0) {
        Pop4.likelihoods1[j,i,h] <- ( 1-adult10.freqs[1,i])^2 
      } else {
        Pop4.likelihoods1[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop4.alleles)){
    
    for (i in 1:ncol(Pop4.alleles)){
      if(Pop4.alleles[j,i,h] == 2) {
        Pop4.likelihoods2[j,i,h] <- adult10.freqs[2,i]^2
      } else if (Pop4.alleles[j,i,h] == 1) {
        Pop4.likelihoods2[j,i,h] <- 2*(adult10.freqs[2,i] * (1-adult10.freqs[2,i]))
      } else if (Pop4.alleles[j,i,h] == 0) {
        Pop4.likelihoods2[j,i,h] <- ( 1-adult10.freqs[2,i])^2 
      } else {
        Pop4.likelihoods2[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop4.alleles)){
    
    for (i in 1:ncol(Pop4.alleles)){
      if(Pop4.alleles[j,i,h] == 2) {
        Pop4.likelihoods3[j,i,h] <- adult10.freqs[3,i]^2
      } else if (Pop4.alleles[j,i,h] == 1) {
        Pop4.likelihoods3[j,i,h] <- 2*(adult10.freqs[3,i] * (1-adult10.freqs[3,i]))
      } else if (Pop4.alleles[j,i,h] == 0) {
        Pop4.likelihoods3[j,i,h] <- ( 1-adult10.freqs[3,i])^2 
      } else {
        Pop4.likelihoods3[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop4.alleles)){
    
    for (i in 1:ncol(Pop4.alleles)){
      if(Pop4.alleles[j,i,h] == 2) {
        Pop4.likelihoods4[j,i,h] <- adult10.freqs[4,i]^2
      } else if (Pop4.alleles[j,i,h] == 1) {
        Pop4.likelihoods4[j,i,h] <- 2*(adult10.freqs[4,i] * (1-adult10.freqs[4,i]))
      } else if (Pop4.alleles[j,i,h] == 0) {
        Pop4.likelihoods4[j,i,h] <- ( 1-adult10.freqs[4,i])^2 
      } else {
        Pop4.likelihoods4[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop4.alleles)){
    
    for (i in 1:ncol(Pop4.alleles)){
      if(Pop4.alleles[j,i,h] == 2) {
        Pop4.likelihoods5[j,i,h] <- adult10.freqs[5,i]^2
      } else if (Pop4.alleles[j,i,h] == 1) {
        Pop4.likelihoods5[j,i,h] <- 2*(adult10.freqs[5,i] * (1-adult10.freqs[5,i]))
      } else if (Pop4.alleles[j,i,h] == 0) {
        Pop4.likelihoods5[j,i,h] <- ( 1-adult10.freqs[5,i])^2 
      } else {
        Pop4.likelihoods5[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop4.alleles)){
    
    for (i in 1:ncol(Pop4.alleles)){
      if(Pop4.alleles[j,i,h] == 2) {
        Pop4.likelihoods6[j,i,h] <- adult10.freqs[6,i]^2
      } else if (Pop4.alleles[j,i,h] == 1) {
        Pop4.likelihoods6[j,i,h] <- 2*(adult10.freqs[6,i] * (1-adult10.freqs[6,i]))
      } else if (Pop4.alleles[j,i,h] == 0) {
        Pop4.likelihoods6[j,i,h] <- ( 1-adult10.freqs[6,i])^2 
      } else {
        Pop4.likelihoods6[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop4.alleles)){
    
    for (i in 1:ncol(Pop4.alleles)){
      if(Pop4.alleles[j,i,h] == 2) {
        Pop4.likelihoods7[j,i,h] <- adult10.freqs[7,i]^2
      } else if (Pop4.alleles[j,i,h] == 1) {
        Pop4.likelihoods7[j,i,h] <- 2*(adult10.freqs[7,i] * (1-adult10.freqs[7,i]))
      } else if (Pop4.alleles[j,i,h] == 0) {
        Pop4.likelihoods7[j,i,h] <- ( 1-adult10.freqs[7,i])^2 
      } else {
        Pop4.likelihoods7[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop4.alleles)){
    
    for (i in 1:ncol(Pop4.alleles)){
      if(Pop4.alleles[j,i,h] == 2) {
        Pop4.likelihoods8[j,i,h] <- adult10.freqs[8,i]^2
      } else if (Pop4.alleles[j,i,h] == 1) {
        Pop4.likelihoods8[j,i,h] <- 2*(adult10.freqs[8,i] * (1-adult10.freqs[8,i]))
      } else if (Pop4.alleles[j,i,h] == 0) {
        Pop4.likelihoods8[j,i,h] <- ( 1-adult10.freqs[8,i])^2 
      } else {
        Pop4.likelihoods8[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop4.alleles)){
    
    for (i in 1:ncol(Pop4.alleles)){
      if(Pop4.alleles[j,i,h] == 2) {
        Pop4.likelihoods9[j,i,h] <- adult10.freqs[9,i]^2
      } else if (Pop4.alleles[j,i,h] == 1) {
        Pop4.likelihoods9[j,i,h] <- 2*(adult10.freqs[9,i] * (1-adult10.freqs[9,i]))
      } else if (Pop4.alleles[j,i,h] == 0) {
        Pop4.likelihoods9[j,i,h] <- ( 1-adult10.freqs[9,i])^2 
      } else {
        Pop4.likelihoods9[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop4.alleles)){
    
    for (i in 1:ncol(Pop4.alleles)){
      if(Pop4.alleles[j,i,h] == 2) {
        Pop4.likelihoods10[j,i,h] <- adult10.freqs[10,i]^2
      } else if (Pop4.alleles[j,i,h] == 1) {
        Pop4.likelihoods10[j,i,h] <- 2*(adult10.freqs[10,i] * (1-adult10.freqs[10,i]))
      } else if (Pop4.alleles[j,i,h] == 0) {
        Pop4.likelihoods10[j,i,h] <- ( 1-adult10.freqs[10,i])^2 
      } else {
        Pop4.likelihoods10[j,i,h] <- 1
      }
    }
  }
}

# Multiply everything together from groups of n individuals originating from GAM 1
Pop4.geno1 <- apply(Pop4.likelihoods1, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop4.geno2 <- apply(Pop4.likelihoods2, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop4.geno3 <- apply(Pop4.likelihoods3, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop4.geno4 <- apply(Pop4.likelihoods4, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop4.geno5 <- apply(Pop4.likelihoods5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop4.geno6 <- apply(Pop4.likelihoods6, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop4.geno7 <- apply(Pop4.likelihoods7, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop4.geno8 <- apply(Pop4.likelihoods8, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop4.geno9 <- apply(Pop4.likelihoods9, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop4.geno10 <- apply(Pop4.likelihoods10, FUN = prod, MARGIN = 1, na.rm = TRUE)

# Assign groups of n individuals originating from GAM 1 to most likely pop
Pop4.ass.df <- data.frame(cbind(Pop4.geno1, Pop4.geno2, Pop4.geno3, Pop4.geno4, Pop4.geno5, Pop4.geno6, Pop4.geno7, Pop4.geno8, Pop4.geno9, Pop4.geno10))
colnames(Pop4.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop4.most.like <- colnames(Pop4.ass.df)[apply(Pop4.ass.df,1, which.max)]

#### Alleles drawn from Pop5, genotype likelihoods calculated for all 10 GAM populations ####
Pop5.likelihoods1 <- array(dim = c(1000,10,n))
Pop5.likelihoods2 <- array(dim = c(1000,10,n))
Pop5.likelihoods3 <- array(dim = c(1000,10,n))
Pop5.likelihoods4 <- array(dim = c(1000,10,n))
Pop5.likelihoods5 <- array(dim = c(1000,10,n))
Pop5.likelihoods6 <- array(dim = c(1000,10,n))
Pop5.likelihoods7 <- array(dim = c(1000,10,n))
Pop5.likelihoods8 <- array(dim = c(1000,10,n))
Pop5.likelihoods9 <- array(dim = c(1000,10,n))
Pop5.likelihoods10 <- array(dim = c(1000,10,n))

# For larvae simulated using GAM pop 5 allele freqs. Then geno likelihoods calculated for all 10 GAM pops
# Calculations for genotype likelihoods from 10 GAM populations for groups of n individuals sampled from GAM 5
for (h in 1:n){
  
  for (j in 1:nrow(Pop5.alleles)){
    
    for (i in 1:ncol(Pop5.alleles)){
      if(Pop5.alleles[j,i,h] == 2) {
        Pop5.likelihoods1[j,i,h] <- adult10.freqs[1,i]^2
      } else if (Pop5.alleles[j,i,h] == 1) {
        Pop5.likelihoods1[j,i,h] <- 2*(adult10.freqs[1,i] * (1-adult10.freqs[1,i]))
      } else if (Pop5.alleles[j,i,h] == 0) {
        Pop5.likelihoods1[j,i,h] <- ( 1-adult10.freqs[1,i])^2 
      } else {
        Pop5.likelihoods1[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop5.alleles)){
    
    for (i in 1:ncol(Pop5.alleles)){
      if(Pop5.alleles[j,i,h] == 2) {
        Pop5.likelihoods2[j,i,h] <- adult10.freqs[2,i]^2
      } else if (Pop5.alleles[j,i,h] == 1) {
        Pop5.likelihoods2[j,i,h] <- 2*(adult10.freqs[2,i] * (1-adult10.freqs[2,i]))
      } else if (Pop5.alleles[j,i,h] == 0) {
        Pop5.likelihoods2[j,i,h] <- ( 1-adult10.freqs[2,i])^2 
      } else {
        Pop5.likelihoods2[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop5.alleles)){
    
    for (i in 1:ncol(Pop5.alleles)){
      if(Pop5.alleles[j,i,h] == 2) {
        Pop5.likelihoods3[j,i,h] <- adult10.freqs[3,i]^2
      } else if (Pop5.alleles[j,i,h] == 1) {
        Pop5.likelihoods3[j,i,h] <- 2*(adult10.freqs[3,i] * (1-adult10.freqs[3,i]))
      } else if (Pop5.alleles[j,i,h] == 0) {
        Pop5.likelihoods3[j,i,h] <- ( 1-adult10.freqs[3,i])^2 
      } else {
        Pop5.likelihoods3[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop5.alleles)){
    
    for (i in 1:ncol(Pop5.alleles)){
      if(Pop5.alleles[j,i,h] == 2) {
        Pop5.likelihoods4[j,i,h] <- adult10.freqs[4,i]^2
      } else if (Pop5.alleles[j,i,h] == 1) {
        Pop5.likelihoods4[j,i,h] <- 2*(adult10.freqs[4,i] * (1-adult10.freqs[4,i]))
      } else if (Pop5.alleles[j,i,h] == 0) {
        Pop5.likelihoods4[j,i,h] <- ( 1-adult10.freqs[4,i])^2 
      } else {
        Pop5.likelihoods4[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop5.alleles)){
    
    for (i in 1:ncol(Pop5.alleles)){
      if(Pop5.alleles[j,i,h] == 2) {
        Pop5.likelihoods5[j,i,h] <- adult10.freqs[5,i]^2
      } else if (Pop5.alleles[j,i,h] == 1) {
        Pop5.likelihoods5[j,i,h] <- 2*(adult10.freqs[5,i] * (1-adult10.freqs[5,i]))
      } else if (Pop5.alleles[j,i,h] == 0) {
        Pop5.likelihoods5[j,i,h] <- ( 1-adult10.freqs[5,i])^2 
      } else {
        Pop5.likelihoods5[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop5.alleles)){
    
    for (i in 1:ncol(Pop5.alleles)){
      if(Pop5.alleles[j,i,h] == 2) {
        Pop5.likelihoods6[j,i,h] <- adult10.freqs[6,i]^2
      } else if (Pop5.alleles[j,i,h] == 1) {
        Pop5.likelihoods6[j,i,h] <- 2*(adult10.freqs[6,i] * (1-adult10.freqs[6,i]))
      } else if (Pop5.alleles[j,i,h] == 0) {
        Pop5.likelihoods6[j,i,h] <- ( 1-adult10.freqs[6,i])^2 
      } else {
        Pop5.likelihoods6[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop5.alleles)){
    
    for (i in 1:ncol(Pop5.alleles)){
      if(Pop5.alleles[j,i,h] == 2) {
        Pop5.likelihoods7[j,i,h] <- adult10.freqs[7,i]^2
      } else if (Pop5.alleles[j,i,h] == 1) {
        Pop5.likelihoods7[j,i,h] <- 2*(adult10.freqs[7,i] * (1-adult10.freqs[7,i]))
      } else if (Pop5.alleles[j,i,h] == 0) {
        Pop5.likelihoods7[j,i,h] <- ( 1-adult10.freqs[7,i])^2 
      } else {
        Pop5.likelihoods7[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop5.alleles)){
    
    for (i in 1:ncol(Pop5.alleles)){
      if(Pop5.alleles[j,i,h] == 2) {
        Pop5.likelihoods8[j,i,h] <- adult10.freqs[8,i]^2
      } else if (Pop5.alleles[j,i,h] == 1) {
        Pop5.likelihoods8[j,i,h] <- 2*(adult10.freqs[8,i] * (1-adult10.freqs[8,i]))
      } else if (Pop5.alleles[j,i,h] == 0) {
        Pop5.likelihoods8[j,i,h] <- ( 1-adult10.freqs[8,i])^2 
      } else {
        Pop5.likelihoods8[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop5.alleles)){
    
    for (i in 1:ncol(Pop5.alleles)){
      if(Pop5.alleles[j,i,h] == 2) {
        Pop5.likelihoods9[j,i,h] <- adult10.freqs[9,i]^2
      } else if (Pop5.alleles[j,i,h] == 1) {
        Pop5.likelihoods9[j,i,h] <- 2*(adult10.freqs[9,i] * (1-adult10.freqs[9,i]))
      } else if (Pop5.alleles[j,i,h] == 0) {
        Pop5.likelihoods9[j,i,h] <- ( 1-adult10.freqs[9,i])^2 
      } else {
        Pop5.likelihoods9[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop5.alleles)){
    
    for (i in 1:ncol(Pop5.alleles)){
      if(Pop5.alleles[j,i,h] == 2) {
        Pop5.likelihoods10[j,i,h] <- adult10.freqs[10,i]^2
      } else if (Pop5.alleles[j,i,h] == 1) {
        Pop5.likelihoods10[j,i,h] <- 2*(adult10.freqs[10,i] * (1-adult10.freqs[10,i]))
      } else if (Pop5.alleles[j,i,h] == 0) {
        Pop5.likelihoods10[j,i,h] <- ( 1-adult10.freqs[10,i])^2 
      } else {
        Pop5.likelihoods10[j,i,h] <- 1
      }
    }
  }
}

# Multiply everything together from groups of n individuals originating from GAM 1
Pop5.geno1 <- apply(Pop5.likelihoods1, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop5.geno2 <- apply(Pop5.likelihoods2, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop5.geno3 <- apply(Pop5.likelihoods3, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop5.geno4 <- apply(Pop5.likelihoods4, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop5.geno5 <- apply(Pop5.likelihoods5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop5.geno6 <- apply(Pop5.likelihoods6, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop5.geno7 <- apply(Pop5.likelihoods7, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop5.geno8 <- apply(Pop5.likelihoods8, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop5.geno9 <- apply(Pop5.likelihoods9, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop5.geno10 <- apply(Pop5.likelihoods10, FUN = prod, MARGIN = 1, na.rm = TRUE)

# Assign groups of n individuals originating from GAM 1 to most likely pop
Pop5.ass.df <- data.frame(cbind(Pop5.geno1, Pop5.geno2, Pop5.geno3, Pop5.geno4, Pop5.geno5, Pop5.geno6, Pop5.geno7, Pop5.geno8, Pop5.geno9, Pop5.geno10))
colnames(Pop5.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop5.most.like <- colnames(Pop5.ass.df)[apply(Pop5.ass.df,1, which.max)]

#### Alleles drawn from Pop6, genotype likelihoods calculated for all 10 GAM populations ####
Pop6.likelihoods1 <- array(dim = c(1000,10,n))
Pop6.likelihoods2 <- array(dim = c(1000,10,n))
Pop6.likelihoods3 <- array(dim = c(1000,10,n))
Pop6.likelihoods4 <- array(dim = c(1000,10,n))
Pop6.likelihoods5 <- array(dim = c(1000,10,n))
Pop6.likelihoods6 <- array(dim = c(1000,10,n))
Pop6.likelihoods7 <- array(dim = c(1000,10,n))
Pop6.likelihoods8 <- array(dim = c(1000,10,n))
Pop6.likelihoods9 <- array(dim = c(1000,10,n))
Pop6.likelihoods10 <- array(dim = c(1000,10,n))

# For larvae simulated using GAM pop 6 allele freqs. Then geno likelihoods calculated for all 10 GAM pops
# Calculations for genotype likelihoods from 10 GAM populations for groups of n individuals sampled from GAM 6
for (h in 1:n){
  
  for (j in 1:nrow(Pop6.alleles)){
    
    for (i in 1:ncol(Pop6.alleles)){
      if(Pop6.alleles[j,i,h] == 2) {
        Pop6.likelihoods1[j,i,h] <- adult10.freqs[1,i]^2
      } else if (Pop6.alleles[j,i,h] == 1) {
        Pop6.likelihoods1[j,i,h] <- 2*(adult10.freqs[1,i] * (1-adult10.freqs[1,i]))
      } else if (Pop6.alleles[j,i,h] == 0) {
        Pop6.likelihoods1[j,i,h] <- ( 1-adult10.freqs[1,i])^2 
      } else {
        Pop6.likelihoods1[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop6.alleles)){
    
    for (i in 1:ncol(Pop6.alleles)){
      if(Pop6.alleles[j,i,h] == 2) {
        Pop6.likelihoods2[j,i,h] <- adult10.freqs[2,i]^2
      } else if (Pop6.alleles[j,i,h] == 1) {
        Pop6.likelihoods2[j,i,h] <- 2*(adult10.freqs[2,i] * (1-adult10.freqs[2,i]))
      } else if (Pop6.alleles[j,i,h] == 0) {
        Pop6.likelihoods2[j,i,h] <- ( 1-adult10.freqs[2,i])^2 
      } else {
        Pop6.likelihoods2[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop6.alleles)){
    
    for (i in 1:ncol(Pop6.alleles)){
      if(Pop6.alleles[j,i,h] == 2) {
        Pop6.likelihoods3[j,i,h] <- adult10.freqs[3,i]^2
      } else if (Pop6.alleles[j,i,h] == 1) {
        Pop6.likelihoods3[j,i,h] <- 2*(adult10.freqs[3,i] * (1-adult10.freqs[3,i]))
      } else if (Pop6.alleles[j,i,h] == 0) {
        Pop6.likelihoods3[j,i,h] <- ( 1-adult10.freqs[3,i])^2 
      } else {
        Pop6.likelihoods3[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop6.alleles)){
    
    for (i in 1:ncol(Pop6.alleles)){
      if(Pop6.alleles[j,i,h] == 2) {
        Pop6.likelihoods4[j,i,h] <- adult10.freqs[4,i]^2
      } else if (Pop6.alleles[j,i,h] == 1) {
        Pop6.likelihoods4[j,i,h] <- 2*(adult10.freqs[4,i] * (1-adult10.freqs[4,i]))
      } else if (Pop6.alleles[j,i,h] == 0) {
        Pop6.likelihoods4[j,i,h] <- ( 1-adult10.freqs[4,i])^2 
      } else {
        Pop6.likelihoods4[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop6.alleles)){
    
    for (i in 1:ncol(Pop6.alleles)){
      if(Pop6.alleles[j,i,h] == 2) {
        Pop6.likelihoods5[j,i,h] <- adult10.freqs[5,i]^2
      } else if (Pop6.alleles[j,i,h] == 1) {
        Pop6.likelihoods5[j,i,h] <- 2*(adult10.freqs[5,i] * (1-adult10.freqs[5,i]))
      } else if (Pop6.alleles[j,i,h] == 0) {
        Pop6.likelihoods5[j,i,h] <- ( 1-adult10.freqs[5,i])^2 
      } else {
        Pop6.likelihoods5[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop6.alleles)){
    
    for (i in 1:ncol(Pop6.alleles)){
      if(Pop6.alleles[j,i,h] == 2) {
        Pop6.likelihoods6[j,i,h] <- adult10.freqs[6,i]^2
      } else if (Pop6.alleles[j,i,h] == 1) {
        Pop6.likelihoods6[j,i,h] <- 2*(adult10.freqs[6,i] * (1-adult10.freqs[6,i]))
      } else if (Pop6.alleles[j,i,h] == 0) {
        Pop6.likelihoods6[j,i,h] <- ( 1-adult10.freqs[6,i])^2 
      } else {
        Pop6.likelihoods6[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop6.alleles)){
    
    for (i in 1:ncol(Pop6.alleles)){
      if(Pop6.alleles[j,i,h] == 2) {
        Pop6.likelihoods7[j,i,h] <- adult10.freqs[7,i]^2
      } else if (Pop6.alleles[j,i,h] == 1) {
        Pop6.likelihoods7[j,i,h] <- 2*(adult10.freqs[7,i] * (1-adult10.freqs[7,i]))
      } else if (Pop6.alleles[j,i,h] == 0) {
        Pop6.likelihoods7[j,i,h] <- ( 1-adult10.freqs[7,i])^2 
      } else {
        Pop6.likelihoods7[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop6.alleles)){
    
    for (i in 1:ncol(Pop6.alleles)){
      if(Pop6.alleles[j,i,h] == 2) {
        Pop6.likelihoods8[j,i,h] <- adult10.freqs[8,i]^2
      } else if (Pop6.alleles[j,i,h] == 1) {
        Pop6.likelihoods8[j,i,h] <- 2*(adult10.freqs[8,i] * (1-adult10.freqs[8,i]))
      } else if (Pop6.alleles[j,i,h] == 0) {
        Pop6.likelihoods8[j,i,h] <- ( 1-adult10.freqs[8,i])^2 
      } else {
        Pop6.likelihoods8[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop6.alleles)){
    
    for (i in 1:ncol(Pop6.alleles)){
      if(Pop6.alleles[j,i,h] == 2) {
        Pop6.likelihoods9[j,i,h] <- adult10.freqs[9,i]^2
      } else if (Pop6.alleles[j,i,h] == 1) {
        Pop6.likelihoods9[j,i,h] <- 2*(adult10.freqs[9,i] * (1-adult10.freqs[9,i]))
      } else if (Pop6.alleles[j,i,h] == 0) {
        Pop6.likelihoods9[j,i,h] <- ( 1-adult10.freqs[9,i])^2 
      } else {
        Pop6.likelihoods9[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop6.alleles)){
    
    for (i in 1:ncol(Pop6.alleles)){
      if(Pop6.alleles[j,i,h] == 2) {
        Pop6.likelihoods10[j,i,h] <- adult10.freqs[10,i]^2
      } else if (Pop6.alleles[j,i,h] == 1) {
        Pop6.likelihoods10[j,i,h] <- 2*(adult10.freqs[10,i] * (1-adult10.freqs[10,i]))
      } else if (Pop6.alleles[j,i,h] == 0) {
        Pop6.likelihoods10[j,i,h] <- ( 1-adult10.freqs[10,i])^2 
      } else {
        Pop6.likelihoods10[j,i,h] <- 1
      }
    }
  }
}

# Multiply everything together from groups of n individuals originating from GAM 1
Pop6.geno1 <- apply(Pop6.likelihoods1, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop6.geno2 <- apply(Pop6.likelihoods2, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop6.geno3 <- apply(Pop6.likelihoods3, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop6.geno4 <- apply(Pop6.likelihoods4, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop6.geno5 <- apply(Pop6.likelihoods5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop6.geno6 <- apply(Pop6.likelihoods6, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop6.geno7 <- apply(Pop6.likelihoods7, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop6.geno8 <- apply(Pop6.likelihoods8, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop6.geno9 <- apply(Pop6.likelihoods9, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop6.geno10 <- apply(Pop6.likelihoods10, FUN = prod, MARGIN = 1, na.rm = TRUE)

# Assign groups of n individuals originating from GAM 1 to most likely pop
Pop6.ass.df <- data.frame(cbind(Pop6.geno1, Pop6.geno2, Pop6.geno3, Pop6.geno4, Pop6.geno5, Pop6.geno6, Pop6.geno7, Pop6.geno8, Pop6.geno9, Pop6.geno10))
colnames(Pop6.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop6.most.like <- colnames(Pop6.ass.df)[apply(Pop6.ass.df,1, which.max)]

#### Alleles drawn from Pop7, genotype likelihoods calculated for all 10 GAM populations ####
Pop7.likelihoods1 <- array(dim = c(1000,10,n))
Pop7.likelihoods2 <- array(dim = c(1000,10,n))
Pop7.likelihoods3 <- array(dim = c(1000,10,n))
Pop7.likelihoods4 <- array(dim = c(1000,10,n))
Pop7.likelihoods5 <- array(dim = c(1000,10,n))
Pop7.likelihoods6 <- array(dim = c(1000,10,n))
Pop7.likelihoods7 <- array(dim = c(1000,10,n))
Pop7.likelihoods8 <- array(dim = c(1000,10,n))
Pop7.likelihoods9 <- array(dim = c(1000,10,n))
Pop7.likelihoods10 <- array(dim = c(1000,10,n))

# For larvae simulated using GAM pop 7 allele freqs. Then geno likelihoods calculated for all 10 GAM pops
# Calculations for genotype likelihoods from 10 GAM populations for groups of n individuals sampled from GAM 7
for (h in 1:n){
  
  for (j in 1:nrow(Pop7.alleles)){
    
    for (i in 1:ncol(Pop7.alleles)){
      if(Pop7.alleles[j,i,h] == 2) {
        Pop7.likelihoods1[j,i,h] <- adult10.freqs[1,i]^2
      } else if (Pop7.alleles[j,i,h] == 1) {
        Pop7.likelihoods1[j,i,h] <- 2*(adult10.freqs[1,i] * (1-adult10.freqs[1,i]))
      } else if (Pop7.alleles[j,i,h] == 0) {
        Pop7.likelihoods1[j,i,h] <- ( 1-adult10.freqs[1,i])^2 
      } else {
        Pop7.likelihoods1[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop7.alleles)){
    
    for (i in 1:ncol(Pop7.alleles)){
      if(Pop7.alleles[j,i,h] == 2) {
        Pop7.likelihoods2[j,i,h] <- adult10.freqs[2,i]^2
      } else if (Pop7.alleles[j,i,h] == 1) {
        Pop7.likelihoods2[j,i,h] <- 2*(adult10.freqs[2,i] * (1-adult10.freqs[2,i]))
      } else if (Pop7.alleles[j,i,h] == 0) {
        Pop7.likelihoods2[j,i,h] <- ( 1-adult10.freqs[2,i])^2 
      } else {
        Pop7.likelihoods2[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop7.alleles)){
    
    for (i in 1:ncol(Pop7.alleles)){
      if(Pop7.alleles[j,i,h] == 2) {
        Pop7.likelihoods3[j,i,h] <- adult10.freqs[3,i]^2
      } else if (Pop7.alleles[j,i,h] == 1) {
        Pop7.likelihoods3[j,i,h] <- 2*(adult10.freqs[3,i] * (1-adult10.freqs[3,i]))
      } else if (Pop7.alleles[j,i,h] == 0) {
        Pop7.likelihoods3[j,i,h] <- ( 1-adult10.freqs[3,i])^2 
      } else {
        Pop7.likelihoods3[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop7.alleles)){
    
    for (i in 1:ncol(Pop7.alleles)){
      if(Pop7.alleles[j,i,h] == 2) {
        Pop7.likelihoods4[j,i,h] <- adult10.freqs[4,i]^2
      } else if (Pop7.alleles[j,i,h] == 1) {
        Pop7.likelihoods4[j,i,h] <- 2*(adult10.freqs[4,i] * (1-adult10.freqs[4,i]))
      } else if (Pop7.alleles[j,i,h] == 0) {
        Pop7.likelihoods4[j,i,h] <- ( 1-adult10.freqs[4,i])^2 
      } else {
        Pop7.likelihoods4[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop7.alleles)){
    
    for (i in 1:ncol(Pop7.alleles)){
      if(Pop7.alleles[j,i,h] == 2) {
        Pop7.likelihoods5[j,i,h] <- adult10.freqs[5,i]^2
      } else if (Pop7.alleles[j,i,h] == 1) {
        Pop7.likelihoods5[j,i,h] <- 2*(adult10.freqs[5,i] * (1-adult10.freqs[5,i]))
      } else if (Pop7.alleles[j,i,h] == 0) {
        Pop7.likelihoods5[j,i,h] <- ( 1-adult10.freqs[5,i])^2 
      } else {
        Pop7.likelihoods5[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop7.alleles)){
    
    for (i in 1:ncol(Pop7.alleles)){
      if(Pop7.alleles[j,i,h] == 2) {
        Pop7.likelihoods6[j,i,h] <- adult10.freqs[6,i]^2
      } else if (Pop7.alleles[j,i,h] == 1) {
        Pop7.likelihoods6[j,i,h] <- 2*(adult10.freqs[6,i] * (1-adult10.freqs[6,i]))
      } else if (Pop7.alleles[j,i,h] == 0) {
        Pop7.likelihoods6[j,i,h] <- ( 1-adult10.freqs[6,i])^2 
      } else {
        Pop7.likelihoods6[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop7.alleles)){
    
    for (i in 1:ncol(Pop7.alleles)){
      if(Pop7.alleles[j,i,h] == 2) {
        Pop7.likelihoods7[j,i,h] <- adult10.freqs[7,i]^2
      } else if (Pop7.alleles[j,i,h] == 1) {
        Pop7.likelihoods7[j,i,h] <- 2*(adult10.freqs[7,i] * (1-adult10.freqs[7,i]))
      } else if (Pop7.alleles[j,i,h] == 0) {
        Pop7.likelihoods7[j,i,h] <- ( 1-adult10.freqs[7,i])^2 
      } else {
        Pop7.likelihoods7[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop7.alleles)){
    
    for (i in 1:ncol(Pop7.alleles)){
      if(Pop7.alleles[j,i,h] == 2) {
        Pop7.likelihoods8[j,i,h] <- adult10.freqs[8,i]^2
      } else if (Pop7.alleles[j,i,h] == 1) {
        Pop7.likelihoods8[j,i,h] <- 2*(adult10.freqs[8,i] * (1-adult10.freqs[8,i]))
      } else if (Pop7.alleles[j,i,h] == 0) {
        Pop7.likelihoods8[j,i,h] <- ( 1-adult10.freqs[8,i])^2 
      } else {
        Pop7.likelihoods8[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop7.alleles)){
    
    for (i in 1:ncol(Pop7.alleles)){
      if(Pop7.alleles[j,i,h] == 2) {
        Pop7.likelihoods9[j,i,h] <- adult10.freqs[9,i]^2
      } else if (Pop7.alleles[j,i,h] == 1) {
        Pop7.likelihoods9[j,i,h] <- 2*(adult10.freqs[9,i] * (1-adult10.freqs[9,i]))
      } else if (Pop7.alleles[j,i,h] == 0) {
        Pop7.likelihoods9[j,i,h] <- ( 1-adult10.freqs[9,i])^2 
      } else {
        Pop7.likelihoods9[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop7.alleles)){
    
    for (i in 1:ncol(Pop7.alleles)){
      if(Pop7.alleles[j,i,h] == 2) {
        Pop7.likelihoods10[j,i,h] <- adult10.freqs[10,i]^2
      } else if (Pop7.alleles[j,i,h] == 1) {
        Pop7.likelihoods10[j,i,h] <- 2*(adult10.freqs[10,i] * (1-adult10.freqs[10,i]))
      } else if (Pop7.alleles[j,i,h] == 0) {
        Pop7.likelihoods10[j,i,h] <- ( 1-adult10.freqs[10,i])^2 
      } else {
        Pop7.likelihoods10[j,i,h] <- 1
      }
    }
  }
}

# Multiply everything together from groups of n individuals originating from GAM 1
Pop7.geno1 <- apply(Pop7.likelihoods1, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop7.geno2 <- apply(Pop7.likelihoods2, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop7.geno3 <- apply(Pop7.likelihoods3, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop7.geno4 <- apply(Pop7.likelihoods4, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop7.geno5 <- apply(Pop7.likelihoods5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop7.geno6 <- apply(Pop7.likelihoods6, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop7.geno7 <- apply(Pop7.likelihoods7, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop7.geno8 <- apply(Pop7.likelihoods8, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop7.geno9 <- apply(Pop7.likelihoods9, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop7.geno10 <- apply(Pop7.likelihoods10, FUN = prod, MARGIN = 1, na.rm = TRUE)

# Assign groups of n individuals originating from GAM 1 to most likely pop
Pop7.ass.df <- data.frame(cbind(Pop7.geno1, Pop7.geno2, Pop7.geno3, Pop7.geno4, Pop7.geno5, Pop7.geno6, Pop7.geno7, Pop7.geno8, Pop7.geno9, Pop7.geno10))
colnames(Pop7.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop7.most.like <- colnames(Pop7.ass.df)[apply(Pop7.ass.df,1, which.max)]

#### Alleles drawn from Pop8, genotype likelihoods calculated for all 10 GAM populations ####
Pop8.likelihoods1 <- array(dim = c(1000,10,n))
Pop8.likelihoods2 <- array(dim = c(1000,10,n))
Pop8.likelihoods3 <- array(dim = c(1000,10,n))
Pop8.likelihoods4 <- array(dim = c(1000,10,n))
Pop8.likelihoods5 <- array(dim = c(1000,10,n))
Pop8.likelihoods6 <- array(dim = c(1000,10,n))
Pop8.likelihoods7 <- array(dim = c(1000,10,n))
Pop8.likelihoods8 <- array(dim = c(1000,10,n))
Pop8.likelihoods9 <- array(dim = c(1000,10,n))
Pop8.likelihoods10 <- array(dim = c(1000,10,n))

# For larvae simulated using GAM pop 8 allele freqs. Then geno likelihoods calculated for all 10 GAM pops
# Calculations for genotype likelihoods from 10 GAM populations for groups of n individuals sampled from GAM 8
for (h in 1:n){
  
  for (j in 1:nrow(Pop8.alleles)){
    
    for (i in 1:ncol(Pop8.alleles)){
      if(Pop8.alleles[j,i,h] == 2) {
        Pop8.likelihoods1[j,i,h] <- adult10.freqs[1,i]^2
      } else if (Pop8.alleles[j,i,h] == 1) {
        Pop8.likelihoods1[j,i,h] <- 2*(adult10.freqs[1,i] * (1-adult10.freqs[1,i]))
      } else if (Pop8.alleles[j,i,h] == 0) {
        Pop8.likelihoods1[j,i,h] <- ( 1-adult10.freqs[1,i])^2 
      } else {
        Pop8.likelihoods1[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop8.alleles)){
    
    for (i in 1:ncol(Pop8.alleles)){
      if(Pop8.alleles[j,i,h] == 2) {
        Pop8.likelihoods2[j,i,h] <- adult10.freqs[2,i]^2
      } else if (Pop8.alleles[j,i,h] == 1) {
        Pop8.likelihoods2[j,i,h] <- 2*(adult10.freqs[2,i] * (1-adult10.freqs[2,i]))
      } else if (Pop8.alleles[j,i,h] == 0) {
        Pop8.likelihoods2[j,i,h] <- ( 1-adult10.freqs[2,i])^2 
      } else {
        Pop8.likelihoods2[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop8.alleles)){
    
    for (i in 1:ncol(Pop8.alleles)){
      if(Pop8.alleles[j,i,h] == 2) {
        Pop8.likelihoods3[j,i,h] <- adult10.freqs[3,i]^2
      } else if (Pop8.alleles[j,i,h] == 1) {
        Pop8.likelihoods3[j,i,h] <- 2*(adult10.freqs[3,i] * (1-adult10.freqs[3,i]))
      } else if (Pop8.alleles[j,i,h] == 0) {
        Pop8.likelihoods3[j,i,h] <- ( 1-adult10.freqs[3,i])^2 
      } else {
        Pop8.likelihoods3[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop8.alleles)){
    
    for (i in 1:ncol(Pop8.alleles)){
      if(Pop8.alleles[j,i,h] == 2) {
        Pop8.likelihoods4[j,i,h] <- adult10.freqs[4,i]^2
      } else if (Pop8.alleles[j,i,h] == 1) {
        Pop8.likelihoods4[j,i,h] <- 2*(adult10.freqs[4,i] * (1-adult10.freqs[4,i]))
      } else if (Pop8.alleles[j,i,h] == 0) {
        Pop8.likelihoods4[j,i,h] <- ( 1-adult10.freqs[4,i])^2 
      } else {
        Pop8.likelihoods4[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop8.alleles)){
    
    for (i in 1:ncol(Pop8.alleles)){
      if(Pop8.alleles[j,i,h] == 2) {
        Pop8.likelihoods5[j,i,h] <- adult10.freqs[5,i]^2
      } else if (Pop8.alleles[j,i,h] == 1) {
        Pop8.likelihoods5[j,i,h] <- 2*(adult10.freqs[5,i] * (1-adult10.freqs[5,i]))
      } else if (Pop8.alleles[j,i,h] == 0) {
        Pop8.likelihoods5[j,i,h] <- ( 1-adult10.freqs[5,i])^2 
      } else {
        Pop8.likelihoods5[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop8.alleles)){
    
    for (i in 1:ncol(Pop8.alleles)){
      if(Pop8.alleles[j,i,h] == 2) {
        Pop8.likelihoods6[j,i,h] <- adult10.freqs[6,i]^2
      } else if (Pop8.alleles[j,i,h] == 1) {
        Pop8.likelihoods6[j,i,h] <- 2*(adult10.freqs[6,i] * (1-adult10.freqs[6,i]))
      } else if (Pop8.alleles[j,i,h] == 0) {
        Pop8.likelihoods6[j,i,h] <- ( 1-adult10.freqs[6,i])^2 
      } else {
        Pop8.likelihoods6[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop8.alleles)){
    
    for (i in 1:ncol(Pop8.alleles)){
      if(Pop8.alleles[j,i,h] == 2) {
        Pop8.likelihoods7[j,i,h] <- adult10.freqs[7,i]^2
      } else if (Pop8.alleles[j,i,h] == 1) {
        Pop8.likelihoods7[j,i,h] <- 2*(adult10.freqs[7,i] * (1-adult10.freqs[7,i]))
      } else if (Pop8.alleles[j,i,h] == 0) {
        Pop8.likelihoods7[j,i,h] <- ( 1-adult10.freqs[7,i])^2 
      } else {
        Pop8.likelihoods7[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop8.alleles)){
    
    for (i in 1:ncol(Pop8.alleles)){
      if(Pop8.alleles[j,i,h] == 2) {
        Pop8.likelihoods8[j,i,h] <- adult10.freqs[8,i]^2
      } else if (Pop8.alleles[j,i,h] == 1) {
        Pop8.likelihoods8[j,i,h] <- 2*(adult10.freqs[8,i] * (1-adult10.freqs[8,i]))
      } else if (Pop8.alleles[j,i,h] == 0) {
        Pop8.likelihoods8[j,i,h] <- ( 1-adult10.freqs[8,i])^2 
      } else {
        Pop8.likelihoods8[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop8.alleles)){
    
    for (i in 1:ncol(Pop8.alleles)){
      if(Pop8.alleles[j,i,h] == 2) {
        Pop8.likelihoods9[j,i,h] <- adult10.freqs[9,i]^2
      } else if (Pop8.alleles[j,i,h] == 1) {
        Pop8.likelihoods9[j,i,h] <- 2*(adult10.freqs[9,i] * (1-adult10.freqs[9,i]))
      } else if (Pop8.alleles[j,i,h] == 0) {
        Pop8.likelihoods9[j,i,h] <- ( 1-adult10.freqs[9,i])^2 
      } else {
        Pop8.likelihoods9[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop8.alleles)){
    
    for (i in 1:ncol(Pop8.alleles)){
      if(Pop8.alleles[j,i,h] == 2) {
        Pop8.likelihoods10[j,i,h] <- adult10.freqs[10,i]^2
      } else if (Pop8.alleles[j,i,h] == 1) {
        Pop8.likelihoods10[j,i,h] <- 2*(adult10.freqs[10,i] * (1-adult10.freqs[10,i]))
      } else if (Pop8.alleles[j,i,h] == 0) {
        Pop8.likelihoods10[j,i,h] <- ( 1-adult10.freqs[10,i])^2 
      } else {
        Pop8.likelihoods10[j,i,h] <- 1
      }
    }
  }
}

# Multiply everything together from groups of n individuals originating from GAM 1
Pop8.geno1 <- apply(Pop8.likelihoods1, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop8.geno2 <- apply(Pop8.likelihoods2, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop8.geno3 <- apply(Pop8.likelihoods3, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop8.geno4 <- apply(Pop8.likelihoods4, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop8.geno5 <- apply(Pop8.likelihoods5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop8.geno6 <- apply(Pop8.likelihoods6, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop8.geno7 <- apply(Pop8.likelihoods7, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop8.geno8 <- apply(Pop8.likelihoods8, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop8.geno9 <- apply(Pop8.likelihoods9, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop8.geno10 <- apply(Pop8.likelihoods10, FUN = prod, MARGIN = 1, na.rm = TRUE)

# Assign groups of n individuals originating from GAM 1 to most likely pop
Pop8.ass.df <- data.frame(cbind(Pop8.geno1, Pop8.geno2, Pop8.geno3, Pop8.geno4, Pop8.geno5, Pop8.geno6, Pop8.geno7, Pop8.geno8, Pop8.geno9, Pop8.geno10))
colnames(Pop8.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop8.most.like <- colnames(Pop8.ass.df)[apply(Pop8.ass.df,1, which.max)]

#### Alleles drawn from Pop9, genotype likelihoods calculated for all 10 GAM populations ####
Pop9.likelihoods1 <- array(dim = c(1000,10,n))
Pop9.likelihoods2 <- array(dim = c(1000,10,n))
Pop9.likelihoods3 <- array(dim = c(1000,10,n))
Pop9.likelihoods4 <- array(dim = c(1000,10,n))
Pop9.likelihoods5 <- array(dim = c(1000,10,n))
Pop9.likelihoods6 <- array(dim = c(1000,10,n))
Pop9.likelihoods7 <- array(dim = c(1000,10,n))
Pop9.likelihoods8 <- array(dim = c(1000,10,n))
Pop9.likelihoods9 <- array(dim = c(1000,10,n))
Pop9.likelihoods10 <- array(dim = c(1000,10,n))

# For larvae simulated using GAM pop 9 allele freqs. Then geno likelihoods calculated for all 10 GAM pops
# Calculations for genotype likelihoods from 10 GAM populations for groups of n individuals sampled from GAM 9
for (h in 1:n){
  
  for (j in 1:nrow(Pop9.alleles)){
    
    for (i in 1:ncol(Pop9.alleles)){
      if(Pop9.alleles[j,i,h] == 2) {
        Pop9.likelihoods1[j,i,h] <- adult10.freqs[1,i]^2
      } else if (Pop9.alleles[j,i,h] == 1) {
        Pop9.likelihoods1[j,i,h] <- 2*(adult10.freqs[1,i] * (1-adult10.freqs[1,i]))
      } else if (Pop9.alleles[j,i,h] == 0) {
        Pop9.likelihoods1[j,i,h] <- ( 1-adult10.freqs[1,i])^2 
      } else {
        Pop9.likelihoods1[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop9.alleles)){
    
    for (i in 1:ncol(Pop9.alleles)){
      if(Pop9.alleles[j,i,h] == 2) {
        Pop9.likelihoods2[j,i,h] <- adult10.freqs[2,i]^2
      } else if (Pop9.alleles[j,i,h] == 1) {
        Pop9.likelihoods2[j,i,h] <- 2*(adult10.freqs[2,i] * (1-adult10.freqs[2,i]))
      } else if (Pop9.alleles[j,i,h] == 0) {
        Pop9.likelihoods2[j,i,h] <- ( 1-adult10.freqs[2,i])^2 
      } else {
        Pop9.likelihoods2[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop9.alleles)){
    
    for (i in 1:ncol(Pop9.alleles)){
      if(Pop9.alleles[j,i,h] == 2) {
        Pop9.likelihoods3[j,i,h] <- adult10.freqs[3,i]^2
      } else if (Pop9.alleles[j,i,h] == 1) {
        Pop9.likelihoods3[j,i,h] <- 2*(adult10.freqs[3,i] * (1-adult10.freqs[3,i]))
      } else if (Pop9.alleles[j,i,h] == 0) {
        Pop9.likelihoods3[j,i,h] <- ( 1-adult10.freqs[3,i])^2 
      } else {
        Pop9.likelihoods3[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop9.alleles)){
    
    for (i in 1:ncol(Pop9.alleles)){
      if(Pop9.alleles[j,i,h] == 2) {
        Pop9.likelihoods4[j,i,h] <- adult10.freqs[4,i]^2
      } else if (Pop9.alleles[j,i,h] == 1) {
        Pop9.likelihoods4[j,i,h] <- 2*(adult10.freqs[4,i] * (1-adult10.freqs[4,i]))
      } else if (Pop9.alleles[j,i,h] == 0) {
        Pop9.likelihoods4[j,i,h] <- ( 1-adult10.freqs[4,i])^2 
      } else {
        Pop9.likelihoods4[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop9.alleles)){
    
    for (i in 1:ncol(Pop9.alleles)){
      if(Pop9.alleles[j,i,h] == 2) {
        Pop9.likelihoods5[j,i,h] <- adult10.freqs[5,i]^2
      } else if (Pop9.alleles[j,i,h] == 1) {
        Pop9.likelihoods5[j,i,h] <- 2*(adult10.freqs[5,i] * (1-adult10.freqs[5,i]))
      } else if (Pop9.alleles[j,i,h] == 0) {
        Pop9.likelihoods5[j,i,h] <- ( 1-adult10.freqs[5,i])^2 
      } else {
        Pop9.likelihoods5[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop9.alleles)){
    
    for (i in 1:ncol(Pop9.alleles)){
      if(Pop9.alleles[j,i,h] == 2) {
        Pop9.likelihoods6[j,i,h] <- adult10.freqs[6,i]^2
      } else if (Pop9.alleles[j,i,h] == 1) {
        Pop9.likelihoods6[j,i,h] <- 2*(adult10.freqs[6,i] * (1-adult10.freqs[6,i]))
      } else if (Pop9.alleles[j,i,h] == 0) {
        Pop9.likelihoods6[j,i,h] <- ( 1-adult10.freqs[6,i])^2 
      } else {
        Pop9.likelihoods6[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop9.alleles)){
    
    for (i in 1:ncol(Pop9.alleles)){
      if(Pop9.alleles[j,i,h] == 2) {
        Pop9.likelihoods7[j,i,h] <- adult10.freqs[7,i]^2
      } else if (Pop9.alleles[j,i,h] == 1) {
        Pop9.likelihoods7[j,i,h] <- 2*(adult10.freqs[7,i] * (1-adult10.freqs[7,i]))
      } else if (Pop9.alleles[j,i,h] == 0) {
        Pop9.likelihoods7[j,i,h] <- ( 1-adult10.freqs[7,i])^2 
      } else {
        Pop9.likelihoods7[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop9.alleles)){
    
    for (i in 1:ncol(Pop9.alleles)){
      if(Pop9.alleles[j,i,h] == 2) {
        Pop9.likelihoods8[j,i,h] <- adult10.freqs[8,i]^2
      } else if (Pop9.alleles[j,i,h] == 1) {
        Pop9.likelihoods8[j,i,h] <- 2*(adult10.freqs[8,i] * (1-adult10.freqs[8,i]))
      } else if (Pop9.alleles[j,i,h] == 0) {
        Pop9.likelihoods8[j,i,h] <- ( 1-adult10.freqs[8,i])^2 
      } else {
        Pop9.likelihoods8[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop9.alleles)){
    
    for (i in 1:ncol(Pop9.alleles)){
      if(Pop9.alleles[j,i,h] == 2) {
        Pop9.likelihoods9[j,i,h] <- adult10.freqs[9,i]^2
      } else if (Pop9.alleles[j,i,h] == 1) {
        Pop9.likelihoods9[j,i,h] <- 2*(adult10.freqs[9,i] * (1-adult10.freqs[9,i]))
      } else if (Pop9.alleles[j,i,h] == 0) {
        Pop9.likelihoods9[j,i,h] <- ( 1-adult10.freqs[9,i])^2 
      } else {
        Pop9.likelihoods9[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop9.alleles)){
    
    for (i in 1:ncol(Pop9.alleles)){
      if(Pop9.alleles[j,i,h] == 2) {
        Pop9.likelihoods10[j,i,h] <- adult10.freqs[10,i]^2
      } else if (Pop9.alleles[j,i,h] == 1) {
        Pop9.likelihoods10[j,i,h] <- 2*(adult10.freqs[10,i] * (1-adult10.freqs[10,i]))
      } else if (Pop9.alleles[j,i,h] == 0) {
        Pop9.likelihoods10[j,i,h] <- ( 1-adult10.freqs[10,i])^2 
      } else {
        Pop9.likelihoods10[j,i,h] <- 1
      }
    }
  }
}

# Multiply everything together from groups of n individuals originating from GAM 1
Pop9.geno1 <- apply(Pop9.likelihoods1, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop9.geno2 <- apply(Pop9.likelihoods2, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop9.geno3 <- apply(Pop9.likelihoods3, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop9.geno4 <- apply(Pop9.likelihoods4, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop9.geno5 <- apply(Pop9.likelihoods5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop9.geno6 <- apply(Pop9.likelihoods6, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop9.geno7 <- apply(Pop9.likelihoods7, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop9.geno8 <- apply(Pop9.likelihoods8, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop9.geno9 <- apply(Pop9.likelihoods9, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop9.geno10 <- apply(Pop9.likelihoods10, FUN = prod, MARGIN = 1, na.rm = TRUE)

# Assign groups of n individuals originating from GAM 1 to most likely pop
Pop9.ass.df <- data.frame(cbind(Pop9.geno1, Pop9.geno2, Pop9.geno3, Pop9.geno4, Pop9.geno5, Pop9.geno6, Pop9.geno7, Pop9.geno8, Pop9.geno9, Pop9.geno10))
colnames(Pop9.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop9.most.like <- colnames(Pop9.ass.df)[apply(Pop9.ass.df,1, which.max)]

#### Alleles drawn from Pop10, genotype likelihoods calculated for all 10 GAM populations ####
Pop10.likelihoods1 <- array(dim = c(1000,10,n))
Pop10.likelihoods2 <- array(dim = c(1000,10,n))
Pop10.likelihoods3 <- array(dim = c(1000,10,n))
Pop10.likelihoods4 <- array(dim = c(1000,10,n))
Pop10.likelihoods5 <- array(dim = c(1000,10,n))
Pop10.likelihoods6 <- array(dim = c(1000,10,n))
Pop10.likelihoods7 <- array(dim = c(1000,10,n))
Pop10.likelihoods8 <- array(dim = c(1000,10,n))
Pop10.likelihoods9 <- array(dim = c(1000,10,n))
Pop10.likelihoods10 <- array(dim = c(1000,10,n))

# For larvae simulated using GAM pop 10 allele freqs. Then geno likelihoods calculated for all 10 GAM pops
# Calculations for genotype likelihoods from 10 GAM populations for groups of n individuals sampled from GAM 10
for (h in 1:n){
  
  for (j in 1:nrow(Pop10.alleles)){
    
    for (i in 1:ncol(Pop10.alleles)){
      if(Pop10.alleles[j,i,h] == 2) {
        Pop10.likelihoods1[j,i,h] <- adult10.freqs[1,i]^2
      } else if (Pop10.alleles[j,i,h] == 1) {
        Pop10.likelihoods1[j,i,h] <- 2*(adult10.freqs[1,i] * (1-adult10.freqs[1,i]))
      } else if (Pop10.alleles[j,i,h] == 0) {
        Pop10.likelihoods1[j,i,h] <- ( 1-adult10.freqs[1,i])^2 
      } else {
        Pop10.likelihoods1[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop10.alleles)){
    
    for (i in 1:ncol(Pop10.alleles)){
      if(Pop10.alleles[j,i,h] == 2) {
        Pop10.likelihoods2[j,i,h] <- adult10.freqs[2,i]^2
      } else if (Pop10.alleles[j,i,h] == 1) {
        Pop10.likelihoods2[j,i,h] <- 2*(adult10.freqs[2,i] * (1-adult10.freqs[2,i]))
      } else if (Pop10.alleles[j,i,h] == 0) {
        Pop10.likelihoods2[j,i,h] <- ( 1-adult10.freqs[2,i])^2 
      } else {
        Pop10.likelihoods2[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop10.alleles)){
    
    for (i in 1:ncol(Pop10.alleles)){
      if(Pop10.alleles[j,i,h] == 2) {
        Pop10.likelihoods3[j,i,h] <- adult10.freqs[3,i]^2
      } else if (Pop10.alleles[j,i,h] == 1) {
        Pop10.likelihoods3[j,i,h] <- 2*(adult10.freqs[3,i] * (1-adult10.freqs[3,i]))
      } else if (Pop10.alleles[j,i,h] == 0) {
        Pop10.likelihoods3[j,i,h] <- ( 1-adult10.freqs[3,i])^2 
      } else {
        Pop10.likelihoods3[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop10.alleles)){
    
    for (i in 1:ncol(Pop10.alleles)){
      if(Pop10.alleles[j,i,h] == 2) {
        Pop10.likelihoods4[j,i,h] <- adult10.freqs[4,i]^2
      } else if (Pop10.alleles[j,i,h] == 1) {
        Pop10.likelihoods4[j,i,h] <- 2*(adult10.freqs[4,i] * (1-adult10.freqs[4,i]))
      } else if (Pop10.alleles[j,i,h] == 0) {
        Pop10.likelihoods4[j,i,h] <- ( 1-adult10.freqs[4,i])^2 
      } else {
        Pop10.likelihoods4[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop10.alleles)){
    
    for (i in 1:ncol(Pop10.alleles)){
      if(Pop10.alleles[j,i,h] == 2) {
        Pop10.likelihoods5[j,i,h] <- adult10.freqs[5,i]^2
      } else if (Pop10.alleles[j,i,h] == 1) {
        Pop10.likelihoods5[j,i,h] <- 2*(adult10.freqs[5,i] * (1-adult10.freqs[5,i]))
      } else if (Pop10.alleles[j,i,h] == 0) {
        Pop10.likelihoods5[j,i,h] <- ( 1-adult10.freqs[5,i])^2 
      } else {
        Pop10.likelihoods5[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop10.alleles)){
    
    for (i in 1:ncol(Pop10.alleles)){
      if(Pop10.alleles[j,i,h] == 2) {
        Pop10.likelihoods6[j,i,h] <- adult10.freqs[6,i]^2
      } else if (Pop10.alleles[j,i,h] == 1) {
        Pop10.likelihoods6[j,i,h] <- 2*(adult10.freqs[6,i] * (1-adult10.freqs[6,i]))
      } else if (Pop10.alleles[j,i,h] == 0) {
        Pop10.likelihoods6[j,i,h] <- ( 1-adult10.freqs[6,i])^2 
      } else {
        Pop10.likelihoods6[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop10.alleles)){
    
    for (i in 1:ncol(Pop10.alleles)){
      if(Pop10.alleles[j,i,h] == 2) {
        Pop10.likelihoods7[j,i,h] <- adult10.freqs[7,i]^2
      } else if (Pop10.alleles[j,i,h] == 1) {
        Pop10.likelihoods7[j,i,h] <- 2*(adult10.freqs[7,i] * (1-adult10.freqs[7,i]))
      } else if (Pop10.alleles[j,i,h] == 0) {
        Pop10.likelihoods7[j,i,h] <- ( 1-adult10.freqs[7,i])^2 
      } else {
        Pop10.likelihoods7[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop10.alleles)){
    
    for (i in 1:ncol(Pop10.alleles)){
      if(Pop10.alleles[j,i,h] == 2) {
        Pop10.likelihoods8[j,i,h] <- adult10.freqs[8,i]^2
      } else if (Pop10.alleles[j,i,h] == 1) {
        Pop10.likelihoods8[j,i,h] <- 2*(adult10.freqs[8,i] * (1-adult10.freqs[8,i]))
      } else if (Pop10.alleles[j,i,h] == 0) {
        Pop10.likelihoods8[j,i,h] <- ( 1-adult10.freqs[8,i])^2 
      } else {
        Pop10.likelihoods8[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop10.alleles)){
    
    for (i in 1:ncol(Pop10.alleles)){
      if(Pop10.alleles[j,i,h] == 2) {
        Pop10.likelihoods9[j,i,h] <- adult10.freqs[9,i]^2
      } else if (Pop10.alleles[j,i,h] == 1) {
        Pop10.likelihoods9[j,i,h] <- 2*(adult10.freqs[9,i] * (1-adult10.freqs[9,i]))
      } else if (Pop10.alleles[j,i,h] == 0) {
        Pop10.likelihoods9[j,i,h] <- ( 1-adult10.freqs[9,i])^2 
      } else {
        Pop10.likelihoods9[j,i,h] <- 1
      }
    }
  }
  
  for (j in 1:nrow(Pop10.alleles)){
    
    for (i in 1:ncol(Pop10.alleles)){
      if(Pop10.alleles[j,i,h] == 2) {
        Pop10.likelihoods10[j,i,h] <- adult10.freqs[10,i]^2
      } else if (Pop10.alleles[j,i,h] == 1) {
        Pop10.likelihoods10[j,i,h] <- 2*(adult10.freqs[10,i] * (1-adult10.freqs[10,i]))
      } else if (Pop10.alleles[j,i,h] == 0) {
        Pop10.likelihoods10[j,i,h] <- ( 1-adult10.freqs[10,i])^2 
      } else {
        Pop10.likelihoods10[j,i,h] <- 1
      }
    }
  }
}

# Multiply everything together from groups of n individuals originating from GAM 1
Pop10.geno1 <- apply(Pop10.likelihoods1, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop10.geno2 <- apply(Pop10.likelihoods2, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop10.geno3 <- apply(Pop10.likelihoods3, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop10.geno4 <- apply(Pop10.likelihoods4, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop10.geno5 <- apply(Pop10.likelihoods5, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop10.geno6 <- apply(Pop10.likelihoods6, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop10.geno7 <- apply(Pop10.likelihoods7, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop10.geno8 <- apply(Pop10.likelihoods8, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop10.geno9 <- apply(Pop10.likelihoods9, FUN = prod, MARGIN = 1, na.rm = TRUE)
Pop10.geno10 <- apply(Pop10.likelihoods10, FUN = prod, MARGIN = 1, na.rm = TRUE)

# Assign groups of n individuals originating from GAM 1 to most likely pop
Pop10.ass.df <- data.frame(cbind(Pop10.geno1, Pop10.geno2, Pop10.geno3, Pop10.geno4, Pop10.geno5, Pop10.geno6, Pop10.geno7, Pop10.geno8, Pop10.geno9, Pop10.geno10))
colnames(Pop10.ass.df) <- c('Pop J', 'Pop I', 'Pop H', 'Pop G', 'Pop F', 'Pop E', 'Pop D', 'Pop C', 'Pop B', 'Pop A')
Pop10.most.like <- colnames(Pop10.ass.df)[apply(Pop10.ass.df,1, which.max)]

#### Combine all the % assignments into a matrix ####
# Rbind & organize the assignment results for groups of individuals simulated from each of 10 GAM populations
# Each row was drawn from a known population. This should be the reference. Columns are predicted.
mat <- rbind(table(factor(Pop10.most.like, levels = c('Pop A', 'Pop B', 'Pop C', 'Pop D', 'Pop E', 'Pop F', 'Pop G', 'Pop H', 'Pop I', 'Pop J')))/1000,
             table(factor(Pop9.most.like, levels = c('Pop A', 'Pop B', 'Pop C', 'Pop D', 'Pop E', 'Pop F', 'Pop G', 'Pop H', 'Pop I', 'Pop J')))/1000,
             table(factor(Pop8.most.like, levels = c('Pop A', 'Pop B', 'Pop C', 'Pop D', 'Pop E', 'Pop F', 'Pop G', 'Pop H', 'Pop I', 'Pop J')))/1000,
             table(factor(Pop7.most.like, levels = c('Pop A', 'Pop B', 'Pop C', 'Pop D', 'Pop E', 'Pop F', 'Pop G', 'Pop H', 'Pop I', 'Pop J')))/1000,
             table(factor(Pop6.most.like, levels = c('Pop A', 'Pop B', 'Pop C', 'Pop D', 'Pop E', 'Pop F', 'Pop G', 'Pop H', 'Pop I', 'Pop J')))/1000,
             table(factor(Pop5.most.like, levels = c('Pop A', 'Pop B', 'Pop C', 'Pop D', 'Pop E', 'Pop F', 'Pop G', 'Pop H', 'Pop I', 'Pop J')))/1000,
             table(factor(Pop4.most.like, levels = c('Pop A', 'Pop B', 'Pop C', 'Pop D', 'Pop E', 'Pop F', 'Pop G', 'Pop H', 'Pop I', 'Pop J')))/1000,
             table(factor(Pop3.most.like, levels = c('Pop A', 'Pop B', 'Pop C', 'Pop D', 'Pop E', 'Pop F', 'Pop G', 'Pop H', 'Pop I', 'Pop J')))/1000,
             table(factor(Pop2.most.like, levels = c('Pop A', 'Pop B', 'Pop C', 'Pop D', 'Pop E', 'Pop F', 'Pop G', 'Pop H', 'Pop I', 'Pop J')))/1000,
             table(factor(Pop1.most.like, levels = c('Pop A', 'Pop B', 'Pop C', 'Pop D', 'Pop E', 'Pop F', 'Pop G', 'Pop H', 'Pop I', 'Pop J')))/1000)

rownames(mat) <- c('Pop A', 'Pop B', 'Pop C', 'Pop D', 'Pop E', 'Pop F', 'Pop G', 'Pop H', 'Pop I', 'Pop J')

#### Plot confusion matrix of groups of n individuals ####
# Nice plot
# png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/assignment_confusion_matrix_groups.png", width=7, height=5, res=300, units="in")

par(
  mar=c(5, 4, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

mat2 <- data.matrix(mat)
# min <- min(mat2, na.rm = TRUE)
# max <- max(mat2,na.rm = TRUE)

rbPal <- colorRampPalette(c('blue', 'red'))
nocol <- 6
color <- rbPal(nocol)
ColorLevels <- seq(0, 0.65, length=length(color)) # this makes the color bar nicer later on
# ColorLevels <- seq(min, max, length=length(color))

# Set layout.  We are going to include a colorbar next to plot.
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1),
       heights=c(1,1))

#plotting margins.  These seem to work well for me.
par(mar = c(5,5,1.5,1), font = 2)

# Plot it up!
image(1:ncol(mat2), 1:nrow(mat2), t(mat2),
      col=color, xlab="Predicted population", ylab="True population",
      axes=FALSE, zlim = c(0, 0.65),
      main= NA, xlim = c(0.5, 11), ylim = c(0.5, 11))
# image(1:ncol(mat2), 1:nrow(mat2), t(mat2),
#       col=color, xlab="Predicted population", ylab="True population",
#       axes=FALSE, zlim = c(min, max),
#       main= NA, xlim = c(0.5, 11), ylim = c(0.5, 11))
abline(0,1)

# Now annotate the plot
axis(side = 1, at=seq(1,10,1), labels=LETTERS[1:10],
     cex.axis=1.0)
axis(side = 2, at=seq(1,10,1), labels=LETTERS[1:10], las= 1,
     cex.axis=1)

# Add colorbar to second plot region
par(mar = c(3.5,2.5,3.5,2))
image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=color,xlab="",ylab="",xaxt="n", las = 1)
mtext('% Assignment', side = 3, cex = 0.7, line = 0.5, adj = 1)

dev.off()
