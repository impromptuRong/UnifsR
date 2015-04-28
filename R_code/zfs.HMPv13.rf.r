# setwd("C:/Users/impromptu/Desktop/Unifsnew")
library(caret)
library(phyloseq)
library(randomForest)
library(Matrix)
library(glmnet)
# source("./R_code/z.source.r")
source("./rfRFE_fold.r")

# load("./0.raw/HMPv13.RData")
##############    OTU Input   ###############
# setwd("C:/Users/impromptu/Desktop/Unifsnew")
phy_tree <- read.tree("./0.raw/HMP/HMPv13.ref.tre")
Leafs <- phy_tree$tip.label
Edges <- c(phy_tree$tip.label, phy_tree$node.label)[phy_tree$edge[,2]]
## OTU Table ##
otu_table <- readMM("./0.raw/HMP/HMPv13.otu.mtx")
rownames <- as.character(as.vector(t(read.csv("./0.raw/HMP/HMPv13.otu.rn", header=FALSE))))
colnames <- as.character(as.vector(t(read.csv("./0.raw/HMP/HMPv13.otu.cn", header=FALSE))))
otu_table@Dimnames <- list(rownames, colnames)
# otu_table <- otu_table[,Leafs]
## Edge Table ##
edgematrix <- readMM("./0.raw/HMP/HMPv13.edge.mtx")
rownames <- as.character(as.vector(t(read.csv("./0.raw/HMP/HMPv13.edge.rn", header=FALSE))))
colnames <- as.character(as.vector(t(read.csv("./0.raw/HMP/HMPv13.edge.cn", header=FALSE))))
edgematrix@Dimnames <- list(rownames, colnames)
# edgematrix <- edgematrix[Samps,c(Edges,"Root")]
## Meta Info ##
Samps <- rownames(otu_table)
sam_table <- read.table("./0.raw/HMP/HMPv13.meta.txt", header=TRUE, row.names=1)[Samps,]
save(edgematrix, otu_table, sam_table, phy_tree, file="./0.hmp.c19.fs.R/HMPv13.fs.RData", compress="gzip")

load("./0.hmp.c19.fs.R/HMPv13.fs.RData")
for (k in 1:50){
##############    EdgeTable    ##############
sprintf("*************************************************")
sprintf("HMPv13.c19.list.R0.%02d.txt", k)
fc <- file(sprintf("./0.hmp.c19.list/HMPv13.c19.list.R0.%02d.txt", k), "r")
mylist <- strsplit(readLines(fc), "\t")
close(fc)
names(mylist) <- sapply(mylist, function(x) x[1])
mylist <- lapply(mylist, function(x) as.numeric(x[-1]))

##############    Generate CV_fold    ##############
trainX <- edgematrix[mylist$train, mylist$edgetaxa]
trainy <- factor(sam_table$HMPbodysubsite[mylist$train])
testX <- edgematrix[mylist$test, mylist$edgetaxa]
testy <- factor(sam_table$HMPbodysubsite[mylist$test])

cv_fold <- list()
for (i in 1:10){
	fold <- list()
	fold[[1]] <- as.matrix(trainX[mylist$foldid!=i,])
	fold[[2]] <- factor(trainy[mylist$foldid!=i])
	fold[[3]] <- as.matrix(trainX[mylist$foldid==i,])
	fold[[4]] <- factor(trainy[mylist$foldid==i])
	cv_fold[[i]] <- fold
}
cv_fold[[11]] <- list()
cv_fold[[11]][[1]] <- cv_fold[[11]][[3]] <- as.matrix(trainX)
cv_fold[[11]][[2]] <- cv_fold[[11]][[4]] <- factor(trainy)

result1 <- rfRFE_fold(cv_fold, step=0.05, recursive=TRUE)
index <- which.max(rowSums(result1$acc_tab[,2:4]))
best1 <- list()
best1$acc <- result1$acc_tab[index,]
best1$mod <- result1$model[[index]]
best1$fea <- result1$feature[[index]]
best1$confmatC <- confusionMatrix(predict(best1$mod, trainX[,best1$fea]), trainy)
best1$confmatT <- confusionMatrix(predict(best1$mod, testX[,best1$fea]), testy)
best1$confmatA <- confusionMatrix(predict(best1$mod, rBind(trainX[,best1$fea], testX[,best1$fea])), c(as.vector(trainy),as.vector(testy)))

##############    Generate CV_fold    ##############
trainX <- otu_table[mylist$train, mylist$otutaxa]
trainy <- factor(sam_table$HMPbodysubsite[mylist$train])
testX <- otu_table[mylist$test, mylist$otutaxa]
testy <- factor(sam_table$HMPbodysubsite[mylist$test])

cv_fold <- list()
for (i in 1:10){
        fold <- list()
        fold[[1]] <- as.matrix(trainX[mylist$foldid!=i,])
        fold[[2]] <- factor(trainy[mylist$foldid!=i])
        fold[[3]] <- as.matrix(trainX[mylist$foldid==i,])
        fold[[4]] <- factor(trainy[mylist$foldid==i])
        cv_fold[[i]] <- fold
}
cv_fold[[11]] <- list()
cv_fold[[11]][[1]] <- cv_fold[[11]][[3]] <- as.matrix(trainX)
cv_fold[[11]][[2]] <- cv_fold[[11]][[4]] <- factor(trainy)

result2 <- rfRFE_fold(cv_fold, step=0.05, recursive=TRUE)
index <- which.max(rowSums(result2$acc_tab[,2:4]))
best2 <- list()
best2$acc <- result2$acc_tab[index,]
best2$mod <- result2$model[[index]]
best2$fea <- result2$feature[[index]]
best2$confmatC <- confusionMatrix(predict(best2$mod, trainX[,best2$fea]), trainy)
best2$confmatT <- confusionMatrix(predict(best2$mod, testX[,best2$fea]), testy)
best2$confmatA <- confusionMatrix(predict(best2$mod, rBind(trainX[,best2$fea], testX[,best2$fea])), c(as.vector(trainy),as.vector(testy)))

save(mylist, result1, best1, result2, best2, file=sprintf("./0.hmp.c19.fs.R/HMPv13.c19.fs.rf.R0.%02d.RData", k), compress="gzip")
sprintf("HMPv13.c19.fs.rf.R0.%02d ... Done!!", k)
sprintf("*************************************************")
}


otu_acc <- edgeacc <- matrix(0, 50, 11)
for (k in 1:50){
load(sprintf("./0.hmp.c19.fs.R/HMPv13.c19.fs.rf.R0.%02d.RData", k))
edgeacc[k, ] <- c(best1$acc[1], sprintf("%.4f",c(best1$acc[2]/100,best1$confmatC$overall[1],best1$confmatC$byClass[1:2], best1$confmatT$overall[1],best1$confmatT$byClass[1:2], best1$confmatA$overall[1],best1$confmatA$byClass[1:2])*100));
otu_acc[k, ] <- c(best2$acc[1], sprintf("%.4f",c(best2$acc[2]/100,best2$confmatC$overall[1],best2$confmatC$byClass[1:2], best2$confmatT$overall[1],best2$confmatT$byClass[1:2], best2$confmatA$overall[1],best2$confmatA$byClass[1:2])*100));
write(c("RandomForest",edgeacc[k,],"RandomForest",otu_acc[k,]), file=sprintf("./0.hmp.c19.acc/HMPv13.c19.%02d.acc.xls", k), ncolumns=24, append=TRUE, sep="\t")
}
colname <- c("No_fea","CVacc","Train_acc","Train_sen","Train_spe","Test_acc","Test_sen","Test_spe","All_acc","All_sen","All_spe")
write.table(edgeacc, file="HMPv13.c19.fs.rf.R0.edgeacc.csv", sep=",", row.names=FALSE, col.names=colname)
write.table(otu_acc, file="HMPv13.c19.fs.rf.R0.otu_acc.csv", sep=",", row.names=FALSE, col.names=colname)

