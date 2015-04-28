setwd("C:/Users/impromptu/Desktop/Unifsnew")
library(phyloseq)
library(ggplot2)
library(vrmlgen)
library(glmnet)
source("./R_code/z.source.r")
load("./0.raw/HMPv13.RData")

##############    EdgeTable    ##############
fc <- file("./HMPv13.c19.list.txt", "r")
mylist <- strsplit(readLines(fc), "\t")
close(fc)
names(mylist) <- sapply(mylist, function(x) x[1])
mylist <- lapply(mylist, function(x) as.numeric(x[-1]))

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
# otu_table <- otu_table[,Leafs]
## Edge Table ##
edgematrix <- readMM("./0.raw/HMP/HMPv13.edge.mtx")
rownames <- as.character(as.vector(t(read.csv("./0.raw/HMP/HMPv13.edge.rn", header=FALSE))))
colnames <- as.character(as.vector(t(read.csv("./0.raw/HMP/HMPv13.edge.cn", header=FALSE))))
edgematrix@Dimnames <- list(rownames, colnames)
# edgematrix <- edgematrix[Samps,c(Edges,"Root")]

##############    Generate FoldId    ##############
trainX <- edgematrix[mylist$train, mylist$edgetaxa]
trainy <- factor(sam_table$HMPbodysubsite[mylist$train])
testX <- edgematrix[mylist$test, mylist$edgetaxa]
testy <- factor(sam_table$HMPbodysubsite[mylist$test])






setwd("C:/Users/impromptu/Desktop/Unifsnew")
trainX <- data.matrix(read.csv('./train.csv', header=TRUE, row.names=1))
testX <- data.matrix(read.csv('./test.csv', header=TRUE, row.names=1))
group <- read.table('./group.txt', header=FALSE)
trainy = factor(group[1:300,])
testy = factor(group[301:371,])

name <- "c19.edge"
setwd("C:/Users/impromptu/Desktop/Unifsnew/classic")


library(randomForest)
oob <- tuneRF(trainX, trainy, cv.fold=10, ntreeTry=1000, stepFactor=1)
mod0 <- randomForest(trainX, trainy, ntree=1000, mtry=oob[which.min(oob[,2]),1], proximity=TRUE)
acc <- 100 - trf$err.rate[nrow(mod0$err.rate),1]*100
confusionMatrix(trainy, predict(mod0, trainX))
confusionMatrix(traty, predict(mod0, testX))

fsrf <- rfcv(trainX, trainy, cv.fold=10, ntree=1000, step=0.05)
acc_tab <- cbind(100-fsrf$error.cv*100, fsrf$n.var)
mod1 <- fsrf


library(party)
cforest <- cforest(group ~ ., data=data.frame(group=factor(group[1:300,]), train), controls=cforest_unbiased(ntree=100, mtry=3))
data.cforest.varimp <- varimpAUC(cforest, conditional=TRUE)








svm <- svm(train, factor(group[1:300,]), cross=10, kernel="radial")
svm.hat <- predict(svm, decision.values=TRUE)
tabTrain <- confusionMatrix(factor(group[1:300,]), svm.hat)
cof3 <- summary(svm)
print(cof3)
result[3] <- cof3$tot.accuracy


cforest <- cforest(group ~ ., data=data.frame(group=factor(group[1:300,]), train), controls=cforest_unbiased(ntree=100, mtry=3))
data.cforest.varimp <- varimpAUC(cforest, conditional=TRUE)


rf <- randomForest(group~., data=data.frame(group=factor(group[1:300,]), train), proximity=TRUE)
acc <- 100-trf$err.rate[nrow(trf$err.rate),1]*100

require(glmnet)
alpha <- seq(0,1,0.1)
cv.enet <- lapply(alpha, function(x){cv.glmnet(train, factor(group[1:300,]), alpha=x, standardize=FALSE, family="binomial", type.measure="class")})
rbind(alpha, sapply(cv.enet, function(x){c(min(x$cvm),x$lambda.min,x$lambda.1se)}))
enet <- cv.enet[[1]]
confusionMatrix(factor(group[1:300,]), predict(enet, train, type="class"))
confusionMatrix(factor(group[301:371,]), predict(enet, test, type="class"))

rownames(result[["cv.enet.par"]]) <- c("alpha", cv.enet[[1]]$name, "lambda.1se", "lambda.min")



load("./0.raw/HMPv13.RData")
setwd("C:/Users/impromptu/Desktop/Unifsnew/classic")
train <- data.matrix(read.csv('./train.csv', header=TRUE, row.names=1))
test <- data.matrix(read.csv('./test.csv', header=TRUE, row.names=1))
group <- read.table('./group.txt', header=FALSE)

library(caret)
train <- as.matrix(edgematrix[rownames(train),])
test <- as.matrix(edgematrix[rownames(test),])

library(e1071)
svm <- svm(train, factor(group[1:300,]), cross=10, kernel="linear")
confusionMatrix(factor(group[1:300,]), predict(svm, decision.values=TRUE))
summary(svm)



require(glmnet)
alpha <- seq(0,1,0.1)
cv.enet <- lapply(alpha, function(x){cv.glmnet(train, factor(group[1:300,]), alpha=x, standardize=FALSE, family="binomial", type.measure="class")})
rbind(alpha, sapply(cv.enet, function(x){c(min(x$cvm),x$lambda.min,x$lambda.1se)}))
enet <- cv.enet[[1]]
confusionMatrix(factor(group[1:300,]), predict(enet, train, type="class"))
confusionMatrix(factor(group[301:371,]), predict(enet, test, type="class"))





