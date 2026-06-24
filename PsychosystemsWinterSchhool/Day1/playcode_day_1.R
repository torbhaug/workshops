#Installign packages
library(gtools)
library(qgraph)
library(mgm)
library(igraph)
library(pcalg)

# Data is already loaded from file so that is all gravy

# Calculating correlations of graph
graph.cor <- cor(graph)

# Calculatign covariance of graph
graph.pcor <- cor2pcor(graph.cor)

# Calculating som significant like stuff
sin.ag(graph.pcor, n = 100)

# Use the mgm package to extract network graph
dat <- as.matrix(graph)
fit.all <- mgm(dat, type=rep("g", 5), lambdaSeq = 0, lambdaSel = "EBIC")
qgraph(fit.all$pairwise$wadj, threshold=0.15)

# To investigate further we try to make a smaller model
fit.345 <- mgm(dat[, c(3,4,5)], type=rep("g", 3), lambdaSeq = 0, lambdaSel = "EBIC")
qgraph(fit.345$pairwise$wadj, threshold=0.15)

# To investigate further we try to make a smaller model
fit.123 <- mgm(dat[, c(1,2,3)], type=rep("g", 3), lambdaSeq = 0, lambdaSel = "EBIC")
qgraph(fit.345$pairwise$wadj, threshold=0.15)




# Visualizing
qgraph(graph.pcor)
