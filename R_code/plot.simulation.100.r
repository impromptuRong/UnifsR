setwd("C:/Users/impromptu/Desktop/Unifsnew")
library(phyloseq)
library(ggplot2)
library(vrmlgen)
library(doBy)
source("./R_code/z.source.r")
source("./R_code/sta_test.yq.20130712.r")
source("./R_code/corrplot.rc.20120531.r")

##############     Simulation Small     ###############
setwd("C:/Users/impromptu/Desktop/Unifs")
## Meta Data ##
metainfo <- read.csv("./0.raw/sim.info.100.csv", header=TRUE, row.names=1)
## Phylogenetic Tree ##
raw.tree <- read.tree("./0.matlab.out/sim_100.ref.tre")
raw.tree$tip.label <- gsub("[^a-zA-Z0-9]", "_", raw.tree$tip.label)
is.rooted(raw.tree)
raw.tree$Nnode
# plot(raw.tree, direction="downwards")
## Abundance Table ##
raw.data <- read.csv("./0.matlab.out/sim_100_sp.summary.csv", header=TRUE, row.names=1)
colnames(raw.data) <- gsub("[^a-zA-Z0-9]", "_", colnames(raw.data))
## Phyloseq Class and Unifrac ##
raw.physeq <- phyloseq(otu_table(raw.data, taxa_are_rows=FALSE), phy_tree(raw.tree), sample_data(metainfo))
unilist <- fastUnifrac(raw.physeq, seq.depth=TRUE)
otumatrix <- data.matrix(raw.data/rowSums(raw.data))
edgematrix <- unilist$edge_matrix
edgelength <- c(phy_tree(raw.physeq)$edge.length, 0)
Nsamp <- nsamples(raw.physeq)
Ntaxa <- ntaxa(raw.physeq)
Nedge <- ncol(unilist$edge_matrix)
#######################################################

#########   Visualization Feature Selection   #########
setwd("C:/Users/impromptu/Desktop/Unifs")
load(file="./0.raw/sim.100.RData")
#######    Tip Label color     #######
digit <- floor(log10(Ntaxa))+1
gfeature <- sprintf(paste("Otu%0",digit,"d",sep=""), unlist(sim.pack$feature))
gnoise <- sprintf(paste("Otu%0",digit,"d",sep=""), unlist(sim.pack$noise))
col.label <- rep("black", Nedge)
names(col.label) <- phy_tree(raw.physeq)$tip.label
col.label[gfeature] <- "green"
col.label[gnoise] <- "red"
#######    Edge Color Style    #######
# jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
# topo.colors
# heat.colors
# terrain.colors
# cm.colors
# rainbow
plot.colors <- topo.colors
barplot(rep(1,30), yaxt="n", col=plot.colors(30))
#######################################################

#########   Selection Tree, PCoA/NMDS Figure  #########
setwd("C:/Users/impromptu/Desktop/Unifs")
name <- "sim_100_sp"
unifs_cv.w <- read.csv(paste("./0.matlab.out/",name,".unifs_cv2.weight.csv",sep=""), header=TRUE, row.names=1)[colnames(unilist$edge_bool),]
# weight <- apply(abs(unifs_cv.w), 2, function(x){x/max(x)})
setwd("C:/Users/impromptu/Desktop/Unifs/0.r.out")
write.csv(unilist$edge_matrix, paste(name, "ori.edge_matrix.csv", sep="."))
write.csv(unilist$edge_bool, paste(name, "ori.edge_bool.csv", sep="."))

for(k in 1:ncol(unifs_cv.w)){
    alpha <- as.numeric(sub("^X","", colnames(unifs_cv.w)[k]))
    weight <- unifs_cv.w[,k]
    weight <- abs(weight/max(weight))
    weight.adj <- log10(1+weight*99)/2
    coref <- weight>10^-4
    
    #######   unifs_cv2 Model   #######
#    col.edge <- factor(weight)
#    levels(col.edge) <- c(0, cumsum(table(col.edge))[1:nlevels(col.edge)-1])+1
#    col.edge <- rev(plot.colors(length(weight)))[as.numeric(as.character(col.edge))]
#    col.edge <- rev(plot.colors(Nedge+round(Nedge/10)))[as.numeric(as.character(col.edge))+round(Nedge/10)]
#    col.edge <- rgb(t(col2rgb(col.edge)), alpha=(1+weight)*128-1, maxColorValue=255)
    
    col.edge <- rev(topo.colors(30))[as.numeric(cut(weight.adj, breaks=30))]
    png(filename=paste(name, sprintf("%.1f",alpha), "fs.tree.png", sep="."), width=3500, height=1500, res=300)
        plot(phy_tree(raw.physeq), edge.color=col.edge, edge.width=as.numeric(cut(weight.adj,4)), tip.color=col.label, font=1, cex=0.5, direction="downwards")
    dev.off()
    
    #######   Distance Matrix   #######
    D <- matrix(0, Nsamp, Nsamp, dimnames=list(rownames(edgematrix), rownames(edgematrix)))
    for(i in 1:nrow(edgematrix)){
        for(j in 1:nrow(edgematrix)){
            a1 <- abs(edgematrix[i,coref]-edgematrix[j,coref]) %*% edgelength[coref]
            a2 <- (edgematrix[i,coref]-edgematrix[j,coref])^2 %*% edgelength[coref]
            D[i,j] <- a1*alpha/2 + sqrt(a1^2*alpha^2/4+(1-alpha)*a2)
        }
    }
    write.csv(D, paste(name, sprintf("%.1f",alpha), "fs.dist.csv", sep="."))
    PCoA <- ordinate(raw.physeq, method="PCoA", distance=as.dist(D))
    NMDS <- ordinate(raw.physeq, method="NMDS", distance=as.dist(D))
    NMDS.3d <- metaMDS(as.dist(D), k=3)
    ppcoa <- plot_ordination(raw.physeq, PCoA, color="Group") + scale_color_manual(values=c("blue","red")) + 
        geom_point(size=5) + theme_bw() + scale_x_continuous(name="PCoA1") + scale_y_continuous(name="PCoA2")
    pnmds <- plot_ordination(raw.physeq, NMDS, color="Group") + scale_color_manual(values=c("blue","red")) + 
        geom_point(size=5) + theme_bw() + scale_x_continuous(name="NMDS1") + scale_y_continuous(name="NMDS2")
    png(filename=paste(name, sprintf("%.1f",alpha), "fs.pcoa2d.png", sep="."), width=3200, height=3000, res=300)
        print(ppcoa)
    dev.off()
    png(filename=paste(name, sprintf("%.1f",alpha), "fs.nmds2d.png", sep="."), width=3200, height=3000, res=300)
        print(pnmds)
    dev.off()
    
    group <- sample_data(raw.physeq)$Group
    td.axes <- c(1,2,3)
    plotdata <- PCoA$vectors[, td.axes]
    cloud3d(plotdata, labels=group, filename=paste(name, sprintf("%.1f",alpha), "fs.pcoa3d.wrl", sep="."), type="vrml", pointstyle="s", 
            metalabels=rownames(plotdata), cols=c("blue","red"), scalefac=4, autoscale="independent", lab.axis=paste("PC",td.axes,sep=""), 
            col.axis="darkblue", showaxis=TRUE, col.lab="darkgray", col.bg="white", cex.lab=1, htmlout=NULL, hwidth=1200, hheight=800, showlegend=TRUE)
    plotdata <- NMDS.3d$points[, td.axes]
    cloud3d(plotdata, labels=group, filename=paste(name, sprintf("%.1f",alpha), "fs.nmds3d.wrl", sep="."), type="vrml", pointstyle="s", 
            metalabels=rownames(plotdata), cols=c("blue","red"), scalefac=4, autoscale="independent", lab.axis=paste("MDS",td.axes,sep=""), 
            col.axis="darkblue", showaxis=TRUE, col.lab="darkgray", col.bg="white", cex.lab=1, htmlout=NULL, hwidth=1200, hheight=800, showlegend=TRUE)

}
#######################################################

#########     NMDS/PCoA Without Selection     #########
for(method in names(unilist$dist)){
    D <- as.matrix(unilist$dist[[method]])
    write.csv(D, paste(name, "ori", method, "dist.csv", sep="."))
    PCoA <- ordinate(raw.physeq, method="PCoA", distance=as.dist(D))
    NMDS <- ordinate(raw.physeq, method="NMDS", distance=as.dist(D))
    NMDS.3d <- metaMDS(as.dist(D), k=3)
    ppcoa <- plot_ordination(raw.physeq, PCoA, color="Group") + scale_color_manual(values=c("blue","red")) + 
        geom_point(size=5) + theme_bw() + scale_x_continuous(name="PCoA1") + scale_y_continuous(name="PCoA2")
    pnmds <- plot_ordination(raw.physeq, NMDS, color="Group") + scale_color_manual(values=c("blue","red")) + 
        geom_point(size=5) + theme_bw() + scale_x_continuous(name="NMDS1") + scale_y_continuous(name="NMDS2")
    
    png(filename=paste(name, "ori", method, "pcoa2d.png", sep="."), width=3200, height=3000, res=300)
        print(ppcoa)
    dev.off()
    png(filename=paste(name, "ori", method, "nmds2d.png", sep="."), width=3200, height=3000, res=300)
        print(pnmds)
    dev.off()
    
    group <- sample_data(raw.physeq)$Group
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
#######################################################

#########     Feature Selection Comparison    #########
setwd("C:/Users/impromptu/Desktop/Unifs")
name <- "sim_100_sp"
unifs_cv.w <- read.csv(paste("./0.matlab.out/",name,".unifs_cv2.weight.csv",sep=""), header=TRUE, row.names=1)[colnames(unilist$edge_bool),]
# weight <- apply(abs(unifs_cv.w), 2, function(x){x/max(x)})
setwd("C:/Users/impromptu/Desktop/Unifs/0.r.out")
otu.fs.result <- FS_basic(otumatrix,  metainfo$Group, ID=rownames(metainfo), name=paste(name,"otu",sep="."), method=c("StaTest","RF","ENet","ISA","SVMRFE"))
ext.fs.result <- FS_basic(edgematrix, metainfo$Group, ID=rownames(metainfo), name=paste(name,"ext",sep="."), method=c("StaTest","RF","ENet","ISA","SVMRFE"))










