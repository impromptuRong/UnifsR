rfRFE_fold <- function(cv_fold, step=0.05, recursive=FALSE, ...){
	require(randomForest)
	require(caret)
	kfold <- length(cv_fold) - 1

	acck <- rep(0, kfold)
	for (m in 1:kfold) {
		foldm <- cv_fold[[m]]
		modk <- randomForest(foldm[[1]], foldm[[2]], foldm[[3]], foldm[[4]], ...)
		acck[m] <- 100 - modk$test$err.rate[nrow(modk$test$err.rate),1]*100
	}
	CVacc <- mean(acck)
	mod <- randomForest(cv_fold[[kfold+1]][[1]], cv_fold[[kfold+1]][[2]], importance=TRUE, ...)
	oob <- 100 - mod$err.rate[nrow(mod$err.rate),1] * 100
	confmat <- confusionMatrix(predict(mod, cv_fold[[kfold+1]][[3]]), cv_fold[[kfold+1]][[4]])
	TCacc <- confmat$overall[1] * 100
	rankingCriteria <- mod$importance[,1]
	nfea <- length(rankingCriteria)
	
	k <- 1
	feature <- list()
	model <- list()
	acc_tab <- matrix(0, nfea, 6, dimnames=list(NULL,c("Nfea","CVacc","TCacc","OOB","Ntree","mtry")))
	acc_tab[k,] <- c(length(rankingCriteria), CVacc, TCacc, oob, modk$ntree, modk$mtry)
	feature[[k]] <- 1:nfea
	model[[k]] <- mod
	print(acc_tab[k,])

	while (acc_tab[k,2] > max(acc_tab[,2])-10 && length(rankingCriteria) > 2){
		k <- k + 1;
		sub_fea <- rankingCriteria > quantile(rankingCriteria, step)
		feature[[k]] <- feature[[k-1]][sub_fea]
		rankingCriteria <- rankingCriteria[sub_fea]

		acck <- rep(0, kfold)
		for (m in 1:kfold) {
			cv_fold[[m]][[1]] <- cv_fold[[m]][[1]][,sub_fea]
			cv_fold[[m]][[3]] <- cv_fold[[m]][[3]][,sub_fea]
			foldm <- cv_fold[[m]]
			modk <- randomForest(foldm[[1]], foldm[[2]], foldm[[3]], foldm[[4]], ...)
			acck[m] <- 100 - modk$test$err.rate[nrow(modk$test$err.rate),1]*100
		}
		CVacc <- mean(acck)

		cv_fold[[kfold+1]][[1]] <- cv_fold[[kfold+1]][[1]][,sub_fea]
		cv_fold[[kfold+1]][[3]] <- cv_fold[[kfold+1]][[3]][,sub_fea]

		mod <- randomForest(cv_fold[[kfold+1]][[1]], cv_fold[[kfold+1]][[2]], importance=recursive, ...)
		oob <- 100 - mod$err.rate[nrow(mod$err.rate),1] * 100
		confmat <- confusionMatrix(predict(mod, cv_fold[[kfold+1]][[3]]), cv_fold[[kfold+1]][[4]])
		TCacc <- confmat$overall[1] * 100
#		TCacc <- (mod$conf[1] + mod$conf[4])/sum(mod$conf[1:4]) * 100
		if (recursive){	rankingCriteria <- mod$importance[,1]	}

		acc_tab[k,] <- c(length(rankingCriteria), CVacc, TCacc, oob, modk$ntree, modk$mtry)
		model[[k]] <- mod
		print(acc_tab[k,])
	}

	acc_tab <- acc_tab[1:k,]
	rownames(acc_tab) <- sprintf(paste("c%0",ceiling(log10(k+1)),"d",sep=""), 1:k)
	feature <- feature[1:k];
	model <- model[1:k];

	return(list(acc_tab=acc_tab, feature=feature, model=model))
}
