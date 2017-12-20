library(seqinr)
library(phangorn)
setwd("~/Desktop/Bioinformatics Final Project")
fdir <- system.file("extdata/trees", package = "phangorn")

dpp <- read.phyDat(file = "clustalw.fasta", format = "fasta")
dm <- dist.ml(dpp)
treeUPGMA <- upgma(dm)
treeNJ <- NJ(dm)
layout(matrix(c(1,2), 2, 1), height=c(1,1))
par(mar = c(0,0,2,0)+ 0.1)
plot(treeUPGMA, main="Dpp homologs")
plot(treeNJ, "unrooted", main="NJ")
fit = pml(treeNJ, data=dpp)
fitJC <- optim.pml(fit, TRUE)
logLik(fitJC)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, 
                    rearrangement = "NNI", control = pml.control(trace = 0))
fitGTR
AIC(fitGTR)
AIC(fitJC)
bs = bootstrap.pml(fitJC, bs=100, optNni=TRUE, 
                   control = pml.control(trace = 0))
par(mfrow=c(1.5,1))
par(mar=c(1,1,3,1))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")
title("Dpp homologs")

hh <- read.phyDat(file = "hedgehogAlign.fasta", format = "fasta")
dm <- dist.ml(hh)
treeUPGMA <- upgma(dm)
treeNJ <- NJ(dm)
layout(matrix(c(1,2), 2, 1), height=c(1,1))
par(mar = c(0,0,2,0)+ 0.1)
plot(treeUPGMA, main="Hedgehog homologs")
plot(treeNJ, "unrooted", main="NJ")
fit = pml(treeNJ, data=hh)
fitJC <- optim.pml(fit, TRUE)
logLik(fitJC)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, 
                    rearrangement = "Stochastic", control = pml.control(trace = 0))
fitGTR
AIC(fitGTR)
AIC(fitJC)
bs = bootstrap.pml(fitJC, bs=100, optNni=TRUE, 
                   control = pml.control(trace = 0))
par(mfrow=c(1.5,1))
par(mar=c(1,1,3,1))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")
title("Hedgehog Homologs")


punt <- read.phyDat(file = "punt.fasta", format = "fasta" )
dm <- dist.ml(punt)
treeUPGMA <- upgma(dm)
treeNJ <- NJ(dm)
layout(matrix(c(1,2), 2, 1), height=c(1,1))
par(mar = c(0,0,2,0)+ 0.1)
plot(treeUPGMA, main="Punt homologs")
plot(treeNJ, "unrooted", main="NJ")
fit = pml(treeNJ, data=punt)
fitJC <- optim.pml(fit, TRUE)
logLik(fitJC)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, 
                    rearrangement = "Ratchet", control = pml.control(trace = 0))
fitGTR
AIC(fitGTR)
AIC(fitJC)
bs = bootstrap.pml(fitJC, bs=100, optNni=TRUE, 
                   control = pml.control(trace = 0))
par(mfrow=c(1.5,1))
par(mar=c(1,1,3,1))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")
title("Punt Homologs")


