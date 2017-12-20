library(seqinr)
library(phangorn)
#set working directory
setwd("~/Documents/GitHub/BioinformaticsFinalProject")
#reads in the trees packages from phangorn
fdir <- system.file("extdata/trees", package = "phangorn")
#creates variable dpp that contains my aligned fasta files for the dpp gene
dpp <- read.phyDat(file = "clustalw.fasta", format = "fasta")
dm <- dist.ml(dpp)

#creats a UPGMA tree of this file
treeUPGMA <- upgma(dm)
#creates a neighbor jointing tree for this file 
treeNJ <- NJ(dm)
#the parameters for the layout and appearence of the trees
layout(matrix(c(1,2), 2, 1), height=c(1,1))
par(mar = c(0,0,2,0)+ 0.1)
#plots the trees using these parameters
plot(treeUPGMA, main="Dpp homologs")
plot(treeNJ, "unrooted", main="NJ")
fit = pml(treeNJ, data=dpp)
fitJC <- optim.pml(fit, TRUE)
#Calculates the log likelihood using the best model 
logLik(fitJC)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, 
                    rearrangement = "NNI", control = pml.control(trace = 0))
fitGTR
AIC(fitGTR)
AIC(fitJC)
#calculate bootstrap values and creates a tree based off these values
bs = bootstrap.pml(fitJC, bs=100, optNni=TRUE, 
                   control = pml.control(trace = 0))
par(mfrow=c(1.5,1))
par(mar=c(1,1,3,1))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")
title("Dpp homologs")

#creates variable dpp that contains my aligned fasta files for the hh gene
hh <- read.phyDat(file = "hedgehogAlign.fasta", format = "fasta")
dm <- dist.ml(hh)
#creats a UPGMA tree of this file
treeUPGMA <- upgma(dm)
#creates a neighbor jointing tree for this file 
treeNJ <- NJ(dm)
#the parameters for the layout and appearence of the trees
layout(matrix(c(1,2), 2, 1), height=c(1,1))
par(mar = c(0,0,2,0)+ 0.1)
#plots the trees using these parameters
plot(treeUPGMA, main="Hedgehog homologs")
plot(treeNJ, "unrooted", main="NJ")
fit = pml(treeNJ, data=hh)
fitJC <- optim.pml(fit, TRUE)
#Calculates the log likelihood using the best model 
logLik(fitJC)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, 
                    rearrangement = "Stochastic", control = pml.control(trace = 0))
fitGTR
AIC(fitGTR)
AIC(fitJC)
#calculate bootstrap values and creates a tree based off these values
bs = bootstrap.pml(fitJC, bs=100, optNni=TRUE, 
                   control = pml.control(trace = 0))
par(mfrow=c(1.5,1))
par(mar=c(1,1,3,1))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")
title("Hedgehog Homologs")

#creates variable dpp that contains my aligned fasta files for the punt gene
punt <- read.phyDat(file = "punt.fasta", format = "fasta" )
dm <- dist.ml(punt)
#creats a UPGMA tree of this file
treeUPGMA <- upgma(dm)
#creates a neighbor jointing tree for this file 
treeNJ <- NJ(dm)
#the parameters for the layout and appearence of the trees
layout(matrix(c(1,2), 2, 1), height=c(1,1))
par(mar = c(0,0,2,0)+ 0.1)
#plots the trees using these parameters
plot(treeUPGMA, main="Punt homologs")
plot(treeNJ, "unrooted", main="NJ")
fit = pml(treeNJ, data=punt)
fitJC <- optim.pml(fit, TRUE)
#Calculates the log likelihood using the best model 
logLik(fitJC)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, 
                    rearrangement = "Ratchet", control = pml.control(trace = 0))
fitGTR
AIC(fitGTR)
AIC(fitJC)
#calculate bootstrap values and creates a tree based off these values
bs = bootstrap.pml(fitJC, bs=100, optNni=TRUE, 
                   control = pml.control(trace = 0))
par(mfrow=c(1.5,1))
par(mar=c(1,1,3,1))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")
title("Punt Homologs")




