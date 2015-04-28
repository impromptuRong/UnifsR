library(MCMCpack)
library(dirmult)
library(vegan)
source("C:/Users/impromptu/Desktop/Unifs/R_code/z.source.r")

#######   Simulation based on phylogenetic tree   #######
# setwd("C:/Users/impromptu/Desktop/Unifs")
# sim.tree <- read.tree("./0.raw/sim.tre")
# sim.pack <- sim_group(tree=sim.tree, Nfeature=5, Nnoise=2)
# par(mfrow=c(1,2))
# plot(sim.tree, direction="downwards")
# plot(sim.pack$tree, direction="downwards")
# sim.tree <- sim.pack$tree
# nAbu <- sim.pack$ngroup
# nOTU <- sim.tree$Nnode + 1

#######   Simulate Species Abundance   #######
setwd("C:/Users/impromptu/Desktop/Unifs")
info <- read.csv("./0.raw/meta.info.csv", header=TRUE, row.names=1)
raw.data <- read.csv("./0.raw/oral.summary.genus.unrounded.csv", header=TRUE, row.names=1)
colnames(raw.data) <- sub("(^X)", "", colnames(raw.data), perl=TRUE)
abu_mod <- sim_abufit(raw.data, info$Periodontitis, plot=TRUE)
names(abu_mod) <- c("A", "B")
dev.off()
sim.pack <- sim_matrix(sim_group(tree=100, Nfeature=5, Nnoise=3, srange=c(4,10)), abu_mod, N=c(20,20))
# load(file="./0.raw/sim.100.RData")
# sim.pack <- sim_matrix(sim.pack, abu_mod, N=c(20,20))
save(sim.pack, file="./0.raw/sim.100.RData")
write.csv(sim.pack$prob.mat, "./0.raw/sim.prob.100.csv")
write.tree(sim.pack$tree, "./0.raw/sim.100.tre")

###################################################
library(MCMCpack)
library(dirmult)
library(phyloseq)
source("C:/Users/impromptu/Desktop/Unifs/R_code/z.source.r")

#######   Simulate Species Abundance   #######
setwd("C:/Users/impromptu/Desktop/Unifs")
info <- read.csv("./0.raw/meta.info.csv", header=TRUE, row.names=1)
raw.data <- read.csv("./0.raw/oral.summary.otu.unrounded.csv", header=TRUE, row.names=1)
colnames(raw.data) <- sub("(^X)", "", colnames(raw.data), perl=TRUE)
abu_mod <- sim_abufit(raw.data, info$Periodontitis, plot=TRUE)
names(abu_mod) <- c("A", "B")
abu_mod[[1]]["k"] <- 2
dev.off()
sim.pack <- sim_matrix(sim_group(tree=1000, Nfeature=10, Nnoise=5, srange=c(15,25)), abu_mod, N=c(50,50), cutoff=c(0.025,0.025,0.04,0.045))

save(sim.pack, file="./0.raw/sim.1000.RData")
write.csv(sim.pack$prob.mat, "./0.raw/sim.prob.1000.csv")
write.tree(sim.pack$tree, "./0.raw/sim.1000.tre")


a1, a2, a3, a4, a5, a6
b1, b2, b3, b4, b5, b6, b7, b8




