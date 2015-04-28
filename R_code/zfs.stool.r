setwd("C:/Users/impromptu/Desktop/Unifsnew")
library(phyloseq)
library(ggplot2)
library(vrmlgen)
library(glmnet)
library(caret)
library(randomForest)
library(Matrix)
library(e1071)
source("./R_code/sta_test.yq.20130712.r")
source("./R_code/corrplot.rc.20120531.r")

##############    OTU Input   ###############
setwd("C:/Users/impromptu/Desktop/Unifsnew")
foldid <- read.csv("./HMP.stool/KEGG.stool.foldid.csv", header=FALSE)
## OTU Table ##
otu_table <- read.csv("./HMP.stool/KEGG_mac.stool.table.csv", header=TRUE, row.names=1)
sam_table <- read.csv("./HMP.stool/KEGG_mac.stool.meta.csv", header=TRUE, row.names=1)
##############    Generate FoldId    ##############
# trainX <- edgematrix[mylist$train, mylist$edgetaxa]
# trainy <- factor(sam_table$HMPbodysubsite[mylist$train])
# testX <- edgematrix[mylist$test, mylist$edgetaxa]
# testy <- factor(sam_table$HMPbodysubsite[mylist$test])
trainX <- otu_table
trainy <- sam_table$group
testX <- otu_table
testy <- sam_table$group

library(randomForest)
oob <- tuneRF(trainX, trainy, cv.fold=10, ntreeTry=1000, stepFactor=1)
mod0 <- randomForest(trainX, trainy, ntree=1000, mtry=oob[which.min(oob[,2]),1], proximity=TRUE)
acc <- 100 - trf$err.rate[nrow(mod0$err.rate),1]*100
confusionMatrix(trainy, predict(mod0, trainX))
confusionMatrix(testy, predict(mod0, testX))

fsrf <- rfcv(trainX, trainy, cv.fold=10, ntree=1000, step=0.05)
acc_tab <- cbind(100-fsrf$error.cv*100, fsrf$n.var)
mod1 <- fsrf

library(party)
cforest <- cforest(group ~ ., data=data.frame(group=trainy, trainX), controls=cforest_unbiased(ntree=1000, mtry=3), OOB=true)
confusionMatrix(testy, predict(cforest))
data.cforest.varimp <- varimpAUC(cforest, conditional=TRUE)

require(e1071)
c <- 2^(-5:15)
g <- 2^(-15:3)

svm_lin_list <- lapply(c, function(x) svm(trainX, trainy, cross=10, kernel="linear", cost=x))
svm_lin_CVacc <- sapply(svm_lin_list, function(x) x$tot.accuracy)
index <- which.max(svm_lin_CVacc)
svm_lin <- svm_lin_list[[index]]
svm_lin_CVacc[index]
confusionMatrix(testy, predict(svm_lin, decision.values=TRUE))
summary(svm_lin)

c <- 2^(-5:15)
g <- 2^(-15:3)
parameter <- expand.grid(c,g)

svm_rbf_list <- apply(parameter, 1, function(x) svm(trainX, trainy, cross=10, kernel="radial", cost=x[1], gamma=x[2])
svm_rbf_CVacc <- sapply(svm_rbf_list, function(x) x$tot.accuracy)
index <- which.max(svm_rbf_CVacc)
svm_rbf <- svm_rbf_list[[index]]
svm_rbf_CVacc[index]
confusionMatrix(testy, predict(svm_rbf, decision.values=TRUE))
summary(svm_rbf)

require(glmnet)
alpha <- seq(0,1,0.01)
cv.enet <- lapply(alpha, function(x){cv.glmnet(data.matrix(trainX), trainy, alpha=x, standardize=FALSE, family="multinomial", type.measure="class")})
result <- rbind(alpha, sapply(cv.enet, function(x){c(min(x$cvm),min(x$cvup),x$lambda.min,x$lambda.1se)}))
rownames(result) <- c("alpha", "cvErrm", "cvErrup", "lambda_min", "lambda_1se")
index <- rev(which(result[2,]==min(result[2,])))[1]
enet <- cv.enet[[index]]
c(alpha[index], enet$lambda.min)
hat <- as.vector(predict(enet, data.matrix(trainX), s="lambda.min", type="class"))
hat <- factor(hat, level=levels(trainy))
confusionMatrix(as.vector(trainy), hat)

rownames(result[["cv.enet.par"]]) <- c("alpha", cv.enet[[1]]$name, "lambda.1se", "lambda.min")




