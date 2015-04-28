setwd("C:/Users/impromptu/Desktop/Unifs")
library(phyloseq)
source("./R_code/z.source.r")

##############     Simulation Small     ###############
setwd("C:/Users/impromptu/Desktop/Unifs")
## Meta Data ##
metainfo <- read.csv("./0.raw/sim.info.1000.csv", header=TRUE, row.names=1)
## Phylogenetic Tree ##
raw.tree <- read.tree("./0.matlab.out/sim_100.ref.tre")
raw.tree$tip.label <- gsub("[^a-zA-Z0-9]", "_", raw.tree$tip.label)
is.rooted(raw.tree)
# raw.tree$Nnode
# plot(raw.tree, direction="downwards")
## Abundance Table ##
raw.data <- read.csv("./0.matlab.out/sim_100_sp.summary.csv", header=TRUE, row.names=1)
colnames(raw.data) <- gsub("[^a-zA-Z0-9]", "_", colnames(raw.data))
## Phyloseq Class and Unifrac ##
raw.physeq <- phyloseq(otu_table(raw.data, taxa_are_rows=FALSE), phy_tree(raw.tree), sample_data(metainfo))
unilist <- fastUnifrac(raw.physeq, seq.depth=TRUE)
write.csv(unilist$edge_bool, "./0.r.out/plot.sim_100_sp.edge_bool.csv")
write.csv(unilist$edge_matrix, "./0.r.out/plot.sim_100_sp.edge_matrix.csv")
write.csv(as.matrix(unilist$dist$dPCoA), "./0.r.out/plot.sim_100_sp.dist.dPCoA.csv")
write.csv(as.matrix(unilist$dist$non_w), "./0.r.out/plot.sim_100_sp.dist.non_w.csv")
write.csv(as.matrix(unilist$dist$w.non_nor), "./0.r.out/plot.sim_100_sp.dist.w.non_nor.csv")
write.csv(as.matrix(unilist$dist$w.nor), "./0.r.out/plot.sim_100_sp.dist.w.nor.csv")
#######################################################

##############     Simulation Large     ###############
setwd("C:/Users/impromptu/Desktop/Unifs")
## Meta Data ##
metainfo <- read.csv("./0.raw/sim.info.1000.csv", header=TRUE, row.names=1)
## Phylogenetic Tree ##
raw.tree <- read.tree("./0.matlab.out/sim_1000.ref.tre")
raw.tree$tip.label <- gsub("[^a-zA-Z0-9]", "_", raw.tree$tip.label)
is.rooted(raw.tree)
# raw.tree$Nnode
# plot(raw.tree, direction="downwards")
## Abundance Table ##
raw.data <- read.csv("./0.matlab.out/sim_1000_sp.summary.csv", header=TRUE, row.names=1)
colnames(raw.data) <- gsub("[^a-zA-Z0-9]", "_", colnames(raw.data))
## Phyloseq Class and Unifrac ##
raw.physeq <- phyloseq(otu_table(raw.data, taxa_are_rows=FALSE), phy_tree(raw.tree), sample_data(metainfo))
unilist <- fastUnifrac(raw.physeq, seq.depth=TRUE)
write.csv(unilist$edge_bool, "./0.r.out/plot.sim_1000_sp.edge_bool.csv")
write.csv(unilist$edge_matrix, "./0.r.out/plot.sim_1000_sp.edge_matrix.csv")
write.csv(as.matrix(unilist$dist$dPCoA), "./0.r.out/plot.sim_1000_sp.dist.dPCoA.csv")
write.csv(as.matrix(unilist$dist$non_w), "./0.r.out/plot.sim_1000_sp.dist.non_w.csv")
write.csv(as.matrix(unilist$dist$w.non_nor), "./0.r.out/plot.sim_1000_sp.dist.w.non_nor.csv")
write.csv(as.matrix(unilist$dist$w.nor), "./0.r.out/plot.sim_1000_sp.dist.w.nor.csv")
#######################################################

##############      HMP Plot Output     ###############
setwd("C:/Users/impromptu/Desktop/Unifs")
## Meta Data ##
metainfo <- read.table('./0.raw/v13_map_uniquebyPSN.txt', header=TRUE, row.names=1)
## Phylogenetic Tree ##
raw.tree <- read.tree("./0.matlab.out/HMP_100_Att_Sup.ref.tre")
raw.tree$tip.label <- gsub("[^a-zA-Z0-9]", "_", raw.tree$tip.label)
is.rooted(raw.tree)
# raw.tree$Nnode
# plot(raw.tree, direction="downwards")
## Taxonomy File ##
taxonomy <- read.csv("./0.raw/otu_taxa_psn_v13.csv", header=TRUE, row.names=1)
rownames(taxonomy) <- gsub("[^a-zA-Z0-9]", "_", rownames(taxonomy))
## Abundance Table ##
raw.data <- read.csv("./0.matlab.out/HMP_100_Att_Sup.summary.csv", header=TRUE, row.names=1)
colnames(raw.data) <- gsub("[^a-zA-Z0-9]", "_", colnames(raw.data))
## Phyloseq Class and Unifrac ##
raw.physeq <- phyloseq(otu_table(raw.data, taxa_are_rows=FALSE), phy_tree(raw.tree), 
                       tax_table(as.matrix(taxonomy)), sample_data(metainfo))
unilist <- fastUnifrac(raw.physeq, seq.depth=TRUE)
write.csv(unilist$edge_bool, "./0.r.out/plot.HMP.edge_bool.csv")
write.csv(unilist$edge_matrix, "./0.r.out/plot.HMP.edge_matrix.csv")
write.csv(as.matrix(unilist$dist$dPCoA), "./0.r.out/plot.HMP.dist.dPCoA.csv")
write.csv(as.matrix(unilist$dist$non_w), "./0.r.out/plot.HMP.dist.non_w.csv")
write.csv(as.matrix(unilist$dist$w.non_nor), "./0.r.out/plot.HMP.dist.w.non_nor.csv")
write.csv(as.matrix(unilist$dist$w.nor), "./0.r.out/plot.HMP.dist.w.nor.csv")
#######################################################

##############   Oral OTU Plot Output   ###############
setwd("C:/Users/impromptu/Desktop/Unifs")
## Meta Data ##
metainfo <- read.csv("./0.raw/meta.info.csv", header=TRUE, row.names=1)
## Phylogenetic Tree ##
raw.tree <- read.tree("./0.matlab.out/oral_otu.ref.tre")
raw.tree$tip.label <- gsub("[^a-zA-Z0-9]", "_", raw.tree$tip.label)
is.rooted(raw.tree)
# raw.tree$Nnode
# plot(raw.tree, direction="downwards")
## Taxonomy File ##
taxonomy <- read.csv("./0.raw/dic.otu.csv", header=TRUE, row.names=1)
rownames(taxonomy) <- gsub("[^a-zA-Z0-9]", "_", rownames(taxonomy))
## Abundance Table ##
raw.data <- read.csv("./0.raw/oral.summary.otu.unrounded.csv", header=TRUE, row.names=1)
colnames(raw.data) <- gsub("[^a-zA-Z0-9]", "_", colnames(raw.data))
## Phyloseq Class and Unifrac ##
raw.physeq <- phyloseq(otu_table(raw.data, taxa_are_rows=FALSE), phy_tree(raw.tree), 
                       tax_table(as.matrix(taxonomy)), sample_data(metainfo))
unilist <- fastUnifrac(raw.physeq, seq.depth=1000)
write.csv(unilist$edge_bool, "./0.r.out/plot.oral_otu.edge_bool.csv")
write.csv(unilist$edge_matrix, "./0.r.out/plot.oral_otu.edge_matrix.csv")
write.csv(as.matrix(unilist$dist$dPCoA), "./0.r.out/plot.oral_otu.dist.dPCoA.csv")
write.csv(as.matrix(unilist$dist$non_w), "./0.r.out/plot.oral_otu.dist.non_w.csv")
write.csv(as.matrix(unilist$dist$w.non_nor), "./0.r.out/plot.oral_otu.dist.w.non_nor.csv")
write.csv(as.matrix(unilist$dist$w.nor), "./0.r.out/plot.oral_otu.dist.w.nor.csv")
#######################################################

##############  Oral Genus Plot Output  ###############
setwd("C:/Users/impromptu/Desktop/Unifs")
## Meta Data ##
metainfo <- read.csv("./0.raw/meta.info.csv", header=TRUE, row.names=1)
## Phylogenetic Tree ##
raw.tree <- read.tree("./0.matlab.out/oral_genus.ref.tre")
raw.tree$tip.label <- gsub("[^a-zA-Z0-9]", "_", raw.tree$tip.label)
is.rooted(raw.tree)
# raw.tree$Nnode
# plot(raw.tree, direction="downwards")
## Taxonomy File ##
taxonomy <- read.csv("./0.raw/dic.genus.csv", header=TRUE, row.names=1)
rownames(taxonomy) <- gsub("[^a-zA-Z0-9]", "_", rownames(taxonomy))
## Abundance Table ##
raw.data <- read.csv("./0.raw/oral.summary.genus.unrounded.csv", header=TRUE, row.names=1)
colnames(raw.data) <- gsub("[^a-zA-Z0-9]", "_", colnames(raw.data))
## Phyloseq Class and Unifrac ##
raw.physeq <- phyloseq(otu_table(raw.data, taxa_are_rows=FALSE), phy_tree(raw.tree), 
                       tax_table(as.matrix(taxonomy)), sample_data(metainfo))
unilist <- fastUnifrac(raw.physeq, seq.depth=1000)
write.csv(unilist$edge_bool, "./0.r.out/plot.oral_genus.edge_bool.csv")
write.csv(unilist$edge_matrix, "./0.r.out/plot.oral_genus.edge_matrix.csv")
write.csv(as.matrix(unilist$dist$dPCoA), "./0.r.out/plot.oral_genus.dist.dPCoA.csv")
write.csv(as.matrix(unilist$dist$non_w), "./0.r.out/plot.oral_genus.dist.non_w.csv")
write.csv(as.matrix(unilist$dist$w.non_nor), "./0.r.out/plot.oral_genus.dist.w.non_nor.csv")
write.csv(as.matrix(unilist$dist$w.nor), "./0.r.out/plot.oral_genus.dist.w.nor.csv")
#######################################################

##############  Oral Phylum Plot Output  ##############
setwd("C:/Users/impromptu/Desktop/Unifs")
## Meta Data ##
metainfo <- read.csv("./0.raw/meta.info.csv", header=TRUE, row.names=1)
## Phylogenetic Tree ##
raw.tree <- read.tree("./0.matlab.out/oral_phylum.ref.tre")
raw.tree$tip.label <- gsub("[^a-zA-Z0-9]", "_", raw.tree$tip.label)
is.rooted(raw.tree)
# raw.tree$Nnode
# plot(raw.tree, direction="downwards")
## Taxonomy File ##
taxonomy <- read.csv("./0.raw/dic.phylum.csv", header=TRUE, row.names=1)
rownames(taxonomy) <- gsub("[^a-zA-Z0-9]", "_", rownames(taxonomy))
## Abundance Table ##
raw.data <- read.csv("./0.raw/oral.summary.phylum.unrounded.csv", header=TRUE, row.names=1)
colnames(raw.data) <- gsub("[^a-zA-Z0-9]", "_", colnames(raw.data))
## Phyloseq Class and Unifrac ##
raw.physeq <- phyloseq(otu_table(raw.data, taxa_are_rows=FALSE), phy_tree(raw.tree), 
                       tax_table(as.matrix(taxonomy)), sample_data(metainfo))
unilist <- fastUnifrac(raw.physeq, seq.depth=1000)
write.csv(unilist$edge_bool, "./0.r.out/plot.oral_phylum.edge_bool.csv")
write.csv(unilist$edge_matrix, "./0.r.out/plot.oral_phylum.edge_matrix.csv")
write.csv(as.matrix(unilist$dist$dPCoA), "./0.r.out/plot.oral_phylum.dist.dPCoA.csv")
write.csv(as.matrix(unilist$dist$non_w), "./0.r.out/plot.oral_phylum.dist.non_w.csv")
write.csv(as.matrix(unilist$dist$w.non_nor), "./0.r.out/plot.oral_phylum.dist.w.non_nor.csv")
write.csv(as.matrix(unilist$dist$w.nor), "./0.r.out/plot.oral_phylum.dist.w.nor.csv")
#######################################################

