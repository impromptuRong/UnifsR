setwd("C:/Users/impromptu/Desktop/Unifsnew")
library(phyloseq)
library(ggplot2)
library(vrmlgen)
library(glmnet)
source("./R_code/z.source.r")

##############    OTU Input   ###############
setwd("C:/Users/impromptu/Desktop/Unifsnew")
phy_tree <- read.tree("./0.raw/HMP/HMPv13.ref.tre")
Leafs <- phy_tree$tip.label
Edges <- c(phy_tree$tip.label, phy_tree$node.label)[phy_tree$edge[,2]]
## OTU Table ##
otu_table <- readMM("./0.raw/HMP/HMPv13.otu.mtx")
rownames <- as.character(as.vector(t(read.csv("./0.raw/HMP/HMPv13.otu.rn", header=FALSE))))
colnames <- as.character(as.vector(t(read.csv("./0.raw/HMP/HMPv13.otu.cn", header=FALSE))))
otu_table@Dimnames <- list(rownames, colnames)
otu_table <- otu_table[,Leafs]
Samps <- rownames(otu_table)
## Edge Table ##
edgematrix <- readMM("./0.raw/HMP/HMPv13.edge.mtx")
rownames <- as.character(as.vector(t(read.csv("./0.raw/HMP/HMPv13.edge.rn", header=FALSE))))
colnames <- as.character(as.vector(t(read.csv("./0.raw/HMP/HMPv13.edge.cn", header=FALSE))))
edgematrix@Dimnames <- list(rownames, colnames)
edgematrix <- edgematrix[Samps,c(Edges,"Root")]
## Meta Info ##
sam_table <- read.table("./0.raw/HMP/HMPv13.meta.txt", header=TRUE, row.names=1)[Samps,]
tax_table <- read.csv("./0.raw/HMP/HMPv13.taxa.csv", header=TRUE, row.names=1)[Leafs,]
## Leaf color ##
species.type <- tax_table$phylum
n <- nlevels(species.type)
hues <- rainbow(n+1)
hues <- seq(15, 375, length=n+1)
col.spe <- c(hcl(h=hues, l=65, c=100)[1:n-1], "#BEBEBE")
names(col.spe) <- c(levels(factor(species.type[species.type!="unclassified"])), "unclassified")
barplot(rep(1,n), yaxt="n", col=col.spe)
tax_table$color <- as.character(col.spe[species.type])
## Distance Matrix ##
phydist <- list()
dist <- read.csv("./0.raw/HMP/HMPv13.dist.dPCoA.csv", header=TRUE, row.names=1)
colnames(dist) <- gsub("^X", "", colnames(dist))
phydist[["dPCoA"]] <- as.dist(dist[Samps,Samps])
dist <- read.csv("./0.raw/HMP/HMPv13.dist.non_w.csv", header=TRUE, row.names=1)
colnames(dist) <- gsub("^X", "", colnames(dist))
phydist[["non_w"]] <- as.dist(dist[Samps,Samps])
dist <- read.csv("./0.raw/HMP/HMPv13.dist.w.non_nor.csv", header=TRUE, row.names=1)
colnames(dist) <- gsub("^X", "", colnames(dist))
phydist[["w.non_nor"]] <- as.dist(dist[Samps,Samps])
dist <- read.csv("./0.raw/HMP/HMPv13.dist.w.nor.csv", header=TRUE, row.names=1)
colnames(dist) <- gsub("^X", "", colnames(dist))
phydist[["w.nor"]] <- as.dist(dist[Samps,Samps])
## Save to File ##
save(otu_table, sam_table, tax_table, col.spe, edgematrix, phydist, file="./0.raw/HMPv13.RData", compress="gzip")
#######################################################

#######    Edge Color Style    #######
# jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
# topo.colors
# heat.colors
# terrain.colors/
# cm.colors
# rainbow
plot.colors <- topo.colors
barplot(rep(1,30), yaxt="n", col=plot.colors(30))
#######################################################


#########   Selection Tree, PCoA/NMDS Figure  #########
#########   Visualization Feature Selection   #########
setwd("C:/Users/impromptu/Desktop/Unifsnew")
######    Read in Coefficient Matrix    ######
name <- "HMPv13.c19"
filename <- list.files(path="./0.hmp.table/", pattern=paste(name, ".fs.table.*", sep=""))[1]
fstable <- data.matrix(read.csv(paste("./0.hmp.table/",filename,sep=""), header=TRUE, row.names=1))
subtaxa <- colnames(fstable)[-1]

weight <- rep(0,length(Edges))
names(weight) <- Edges
weight[subtaxa] <- abs(fstable[1,-1])
weight <- log10(1+weight/max(weight)*99)/2
col.edge <- rev(topo.colors(30))[as.numeric(cut(weight, breaks=30))]
png(paste("./0.hmp.table/",gsub("csv$", "tree.png", filename),sep=""), width=8500, height=3000, res=300)
    plot(phy_tree, edge.color=col.edge, edge.width=as.numeric(cut(weight,4))/2, tip.color=tax_table[Leafs,]$color, font=1, cex=0.2, direction="downwards")
dev.off()

########    Try Using subtree for ploting   ########
subtree <- read.tree("./0.hmp.table/HMPv13.c19.itol.tre")
fultaxa <- read.csv("./0.raw/HMP/HMPv13.taxa.full.csv", header=TRUE, row.names=1)
fultaxa$color <- as.character(col.spe[fultaxa$phylum])
subLeafs <- subtree$tip.label
subEdges <- c(subtree$tip.label, subtree$node.label)[subtree$edge[,2]]

weight <- rep(0,length(subEdges))
names(weight) <- subEdges
weight[subtaxa] <- abs(fstable[1,-1])
weight <- log10(1+weight/max(weight)*99)/2
col.edge <- rev(topo.colors(30))[as.numeric(cut(weight, breaks=30))]

png(paste("./0.hmp.table/",gsub("csv$", "tree.itol.png", filename),sep=""), width=6500, height=6000, res=300)
    plot(subtree, edge.color=col.edge, edge.width=as.numeric(cut(weight,4))*2, tip.color=fultaxa[subLeafs,]$color, font=1, cex=0.2, direction="downwards")
dev.off()

itol <- data.frame(ID=c(subEdges,subtree$tip.label), type=c(rep("clade", length(edgename)),rep("range", length(subtree$tip.label))), 
                   color=c(sub("FF$", "", col.edge),fultaxa[subLeafs,]$color))
write.table(itol[-1,], paste("./0.hmp.table/",gsub("csv$", "itol.label", filename),sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(itol[itol$type=="clade",], "./0.hmp.table/try.label", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

#######   Distance Analysis   #######
D <- read.csv(paste("./0.hmp.table/", name,".fs.dist.dPCoA.csv",sep=""), header=TRUE, row.names=1)
subSamples <- colnames(D) <- rownames(D)
PCoA <- pcoa(as.dist(D))
plotdata <- data.frame(PCoA$vectors[,c(1,2)], sam_table[subSamples,])
ppcoa <- ggplot(plotdata, aes_string(x="Axis.1",y="Axis.2",color="HMPbodysubsite")) + scale_color_manual(values=c("blue","red")) + 
    geom_point(size=5) + theme_bw() + scale_x_continuous(name="PCoA1") + scale_y_continuous(name="PCoA2")
png(paste("./0.hmp.table/",gsub("csv$", "pcoa2d.png", filename),sep=""), width=3200, height=3000, res=300)
    print(ppcoa)
dev.off()

NMDS <- metaMDS(as.dist(D), k=2, trymax=100)
plotdata <- data.frame(NMDS$points, sam_table[subSamples,])
pnmds <- ggplot(plotdata, aes_string(x="MDS1",y="MDS2",color="HMPbodysubsite")) + scale_color_manual(values=c("blue","red")) + 
    geom_point(size=5) + theme_bw() + scale_x_continuous(name="NMDS1") + scale_y_continuous(name="NMDS2")
png(paste("./0.hmp.table/",gsub("csv$", "nmds2d.png", filename),sep=""), width=3200, height=3000, res=300)
    print(pnmds)
dev.off()


alpha = 0.9604
group <- sam_table[subSamples,]$HMPbodysubsite
td.axes <- c(1,2,3)
plotdata <- PCoA$vectors[, td.axes]
cloud3d(plotdata, labels=group, filename=paste("./0.hmp.table/",gsub("csv$", "pcoa3d.wrl", filename),sep=""), type="vrml", pointstyle="s", 
        metalabels=rownames(plotdata), cols=c("blue","red"), scalefac=4, autoscale="independent", lab.axis=paste("PC",td.axes,sep=""), 
        col.axis="darkblue", showaxis=TRUE, col.lab="darkgray", col.bg="white", cex.lab=1, htmlout=NULL, hwidth=1200, hheight=800, showlegend=TRUE)

NMDS.3d <- metaMDS(as.dist(D), k=3, trymax=100)
plotdata <- NMDS.3d$points[, td.axes]
cloud3d(plotdata, labels=group, filename=paste("./0.hmp.table/",gsub("csv$", "nmds3d.wrl", filename),sep=""), type="vrml", pointstyle="s", 
        metalabels=rownames(plotdata), cols=c("blue","red"), scalefac=4, autoscale="independent", lab.axis=paste("MDS",td.axes,sep=""), 
        col.axis="darkblue", showaxis=TRUE, col.lab="darkgray", col.bg="white", cex.lab=1, htmlout=NULL, hwidth=1200, hheight=800, showlegend=TRUE)

#######################################################

#########     NMDS/PCoA Without Selection     #########
for(method in names(unilist$dist)){
    D <- as.matrix(unilist$dist[[method]])
    write.csv(D, paste(name, "ori", method, "dist.csv", sep="."))
    PCoA <- ordinate(raw.physeq, method="PCoA", distance=as.dist(D))
    NMDS <- ordinate(raw.physeq, method="NMDS", distance=as.dist(D))
    NMDS.3d <- metaMDS(as.dist(D), k=3)
    ppcoa <- plot_ordination(raw.physeq, PCoA, color="HMPbodysubsite") + scale_color_manual(values=c("blue","red")) + 
        geom_point(size=5) + theme_bw() + scale_x_continuous(name="PCoA1") + scale_y_continuous(name="PCoA2")
    pnmds <- plot_ordination(raw.physeq, NMDS, color="HMPbodysubsite") + scale_color_manual(values=c("blue","red")) + 
        geom_point(size=5) + theme_bw() + scale_x_continuous(name="NMDS1") + scale_y_continuous(name="NMDS2")
    
    png(filename=paste(name, "ori", method, "pcoa2d.png", sep="."), width=3200, height=3000, res=300)
    print(ppcoa)
    dev.off()
    png(filename=paste(name, "ori", method, "nmds2d.png", sep="."), width=3200, height=3000, res=300)
    print(pnmds)
    dev.off()
    
    group <- sample_data(raw.physeq)$HMPbodysubsite
    td.axes <- c(1,2,3)
    plotdata <- PCoA$vectors[, td.axes]
    cloud3d(plotdata, labels=group, filename=paste(name, "ori", method, "pcoa3d.wrl", sep="."), type="vrml", pointstyle="s", 
            metalabels=rownames(plotdata), cols=c("blue","red"), scalefac=4, autoscale="independent", lab.axis=paste("PC",td.axes,sep=""), 
            col.axis="darkblue", showaxis=TRUE, col.lab="darkgray", col.bg="white", cex.lab=1, htmlout=NULL, hwidth=1200, hheight=800, showlegend=TRUE)
    plotdata <- NMDS.3d$points[, td.axes]
    cloud3d(plotdata, labels=group, filename=paste(name, "ori", method, "nmds3d.wrl", sep="."), type="vrml", pointstyle="s", 
            metalabels=rownames(plotdata), cols=c("blue","red"), scalefac=4, autoscale="independent", lab.axis=paste("MDS",td.axes,sep=""), 
            col.axis="darkblue", showaxis=TRUE, col.lab="darkgray", col.bg="white", cex.lab=1, htmlout=NULL, hwidth=1200, hheight=800, showlegend=TRUE)
}




