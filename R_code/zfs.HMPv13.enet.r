# setwd("C:/Users/impromptu/Desktop/Unifsnew")
library(glmnet)
library(ape)
library(caret)

# load("./0.raw/HMPv13.RData")
##############    OTU Input   ###############
# setwd("C:/Users/impromptu/Desktop/Unifsnew")
# phy_tree <- read.tree("./0.raw/HMP/HMPv13.ref.tre")
# Leafs <- phy_tree$tip.label
# Edges <- c(phy_tree$tip.label, phy_tree$node.label)[phy_tree$edge[,2]]
## OTU Table ##
# otu_table <- readMM("./0.raw/HMP/HMPv13.otu.mtx")
# rownames <- as.character(as.vector(t(read.csv("./0.raw/HMP/HMPv13.otu.rn", header=FALSE))))
# colnames <- as.character(as.vector(t(read.csv("./0.raw/HMP/HMPv13.otu.cn", header=FALSE))))
# otu_table@Dimnames <- list(rownames, colnames)
# otu_table <- otu_table[,Leafs]
## Edge Table ##
# edgematrix <- readMM("./0.raw/HMP/HMPv13.edge.mtx")
# rownames <- as.character(as.vector(t(read.csv("./0.raw/HMP/HMPv13.edge.rn", header=FALSE))))
# colnames <- as.character(as.vector(t(read.csv("./0.raw/HMP/HMPv13.edge.cn", header=FALSE))))
# edgematrix@Dimnames <- list(rownames, colnames)
# edgematrix <- edgematrix[Samps,c(Edges,"Root")]
## Meta Info ##
# Samps <- rownames(otu_table)
# sam_table <- read.table("./0.raw/HMP/HMPv13.meta.txt", header=TRUE, row.names=1)[Samps,]
# save(edgematrix, otu_table, sam_table, phy_tree, file="./HMPv13.fs.rf.meta.RData", compress="gzip")

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
foldid <- mylist$foldid

##############    Remove 0 branch length   #############
Edges <- c(c(phy_tree$tip.label, phy_tree$node.label)[phy_tree$edge[,2]],"Root")
branch <- c(phy_tree$edge.length,0)
branch[branch<1e-3] <- 0;
names(branch) <- Edges
penalty <- 1/branch[colnames(trainX)]
exclude <- which(is.infinite(penalty))
mylist$exclude <- exclude

##############    Cross-Validation Test    ##############
# require(doMC)
# registerDoMC(cores=16)

result <- list()
alpha <- seq(0, 1, 0.01)
family <- ifelse(nlevels(trainy)>2, "multinomial", "binomial")
measure <- ifelse(nlevels(trainy)==2&&length(trainy)>=100, "auc", "class")
measure <- "class"
# cv_mod <- lapply(alpha, function(x){cv.glmnet(trainX, trainy, type.measure=measure, nfolds=max(foldid), foldid=foldid, family=family, alpha=x, standardize=FALSE, exclude=exclude, penalty.factor=penalty)})
cv.enet <- lapply(alpha, function(x){cv.glmnet(trainX[,-exclude], trainy, type.measure=measure, nfolds=max(foldid), foldid=foldid, family=family, alpha=x, standardize=FALSE, penalty.factor=penalty[-exclude])})
result[["cv_mod"]] <- cv.enet

# Get min cv_1se and cv_min 
cv_1se <- sapply(cv.enet, function(x){x$cvm[x$lambda==x$lambda.1se]})
cv_min <- sapply(cv.enet, function(x){x$cvm[x$lambda==x$lambda.min]})
idx_1se <- rev(which(cv_1se==min(cv_1se)))[1]
idx_min <- rev(which(cv_min==min(cv_min)))[1]
coef_1se <- rep(0, length(mylist$edgetaxa)+1)
coef_1se[-mylist$exclude-1] <- as.vector(coef(cv.enet[[idx_1se]], s="lambda.1se"))
coef_min <- rep(0, length(mylist$edgetaxa)+1)
coef_min[-mylist$exclude-1] <- as.vector(coef(cv.enet[[idx_min]], s="lambda.min"))

parlist <- matrix(c(alpha[idx_1se],alpha[idx_min],cv_1se[idx_1se],cv_min[idx_min],cv.enet[[idx_1se]]$lambda.1se,cv.enet[[idx_min]]$lambda.min,sum(abs(coef_1se)>0)-1,sum(abs(coef_min)>0)-1,idx_1se,idx_min), nrow=2, dimnames=list(c("cv_1se","cv_min"), c("alpha",cv.enet[[idx_min]]$name,"lambda","No. feature","index")))
result[["parlist"]] <- parlist
result[["coef.1se"]] <- coef_1se
result[["coef.min"]] <- coef_min

# Calculate CV accuracy, Classification accuracy and Test accuracy
lb <- levels(trainy)
pred <- factor(sign(predict(cv.enet[[idx_1se]], trainX[,-exclude], s="lambda.1se")), label=lb)
confmatC_1se <- confusionMatrix(pred, trainy)
pred <- factor(sign(predict(cv.enet[[idx_1se]], testX[,-exclude], s="lambda.1se")), label=lb)
confmatT_1se <- confusionMatrix(pred, testy)
pred <- factor(sign(predict(cv.enet[[idx_1se]], rBind(trainX,testX)[,-exclude], s="lambda.1se")), label=lb)
confmatA_1se <- confusionMatrix(pred, factor(c(as.vector(trainy),as.vector(testy))))
pred <- factor(sign(predict(cv.enet[[idx_min]], trainX[,-exclude], s="lambda.min")), label=lb)
confmatC_min <- confusionMatrix(pred, trainy)
pred <- factor(sign(predict(cv.enet[[idx_min]], testX[,-exclude], s="lambda.min")), label=lb)
confmatT_min <- confusionMatrix(pred, testy)
pred <- factor(sign(predict(cv.enet[[idx_min]], rBind(trainX,testX)[,-exclude], s="lambda.min")), label=lb)
confmatA_min <- confusionMatrix(pred, factor(c(as.vector(trainy),as.vector(testy))))
result[["confmat.1se"]] <- list(confmatC=confmatC_1se, contmatT=confmatT_1se, confmatA=confmatA_1se)
result[["confmat.min"]] <- list(confmatC=confmatC_min, contmatT=confmatT_min, confmatA=confmatA_min)

# Generate Output Table #
output <- matrix(c(parlist[1,4], 1-parlist[1,2], confmatC_1se$overall[1], confmatC_1se$byClass[1:2], confmatT_1se$overall[1], confmatT_1se$byClass[1:2], confmatA_1se$overall[1], confmatA_1se$byClass[1:2], parlist[2,4], 1-parlist[2,2], confmatC_min$overall[1], confmatC_min$byClass[1:2], confmatT_min$overall[1], confmatT_min$byClass[1:2], confmatA_min$overall[1], confmatA_min$byClass[1:2]), nrow=2, byrow=TRUE, dimnames=list(c("cv_1se","cv_min"), c("No_fea","CVacc","Train_acc","Train_sen","Train_spe","Test_acc","Test_sen","Test_spe","All_acc","All_sen","All_spe")))
output[,2:11] <- output[,2:11]*100
result[["output"]] <- output
output[,2:11] <- sprintf("%.4f",output[,2:11])

save(mylist, result, file=sprintf("./0.hmp.c19.fs.R/HMPv13.c19.fs.enet.w.R0.%02d.RData", k), compress="gzip")

}

load("./0.hmp.c19.fs.R/HMPv13.fs.RData")
for (k in 1:50){
output <- result$output
output[,2:11] <- sprintf("%.4f",output[,2:11])
write(c("Glmnet_1se_w", output[1,], "Glmnet_1se_w", output[1,]), file=sprintf("./0.hmp.c19.acc/HMPv13.c19.%02d.acc.xls", k), ncolumns=24, append=TRUE, sep="\t")
write(c("Glmnet_min_w", output[2,], "Glmnet_min_w", otuput[2,]), file=sprintf("./0.hmp.c19.acc/HMPv13.c19.%02d.acc.xls", k), ncolumns=24, append=TRUE, sep="\t")
}

