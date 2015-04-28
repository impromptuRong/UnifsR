cal_ent <- function(X, y){
	y <- as.factor(y)
	S <- cbind(y=table(y), table(y,X))
	H <- rep(0,ncol(S))
	for(i in 1:ncol(S)){
		Si <- S[,i]
		H[i] <- -sum(sapply(Si, function(x){ifelse(x>0,x/sum(Si)*log2(x/sum(Si)),0)}))
	}
	result <- list()
	result[[1]] <- rbind(S,H)
	H[1] <- -H[1]
	result[[2]] <- -H %*% colSums(S)/length(y)
	return(result)
}

-(5/9)*log2(5/9)-(4/9)*log2(4/9) - 1/9 * (-(2/5)*log2(2/5)-(3/5)*log2(3/5)) - 8/9*(-(2/4)*log2(2/4)-(2/4)*log2(2/4))

X <- c(1,6,5,4,7,3,8,7,5)
border <- unique(sort(X))
border <- sapply(1:(length(border)-1), function(x) {(border[x]+border[x+1])/2})

for(x in 1:length(border)){
	calt <- ifelse(X>border[x],1,-1)
	print(cal_ent(calt, y))
}

