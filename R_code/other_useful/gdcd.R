library(lassoshooting)
# make data 
rm(list = ls(all = TRUE)) # make sure previous work is clear
ls()
p = 70; n = 300
x <- matrix(rnorm(n*p), nrow=n, ncol=p)
# create the y-matrix of dependent variables
y <- matrix(rnorm(n), nrow=n)
m <- nrow(y)
# analytical results with matrix algebra
analytical <- solve(t(x)%*%x)%*%t(x)%*%y 
# w/o feature scaling
grad <- function(x, y, theta) {
    gradient <- (1/m)* (t(x) %*% ((x %*% t(theta)) - y))
    return(t(gradient))
}

# define gradient descent update algorithm
grad.descent <- function(x, maxit){
    thetas = matrix(0,ncol=dim(x)[2], nrow=maxit)
    theta <- matrix(0, nrow=1, ncol=dim(x)[2])
    numcoord = dim(x)[2]
    
    alpha = .1 # set learning rate
    
    for (i in 1:maxit) {
        theta <- theta - alpha*grad(x, y, theta) 
        thetas[i, ] = theta 
    }
    return(thetas)
}
coord.descent <- function(x, maxit, lambda=0.1){
    p = dim(x)[2]
    theta <- matrix(0, ncol=p, nrow=maxit) # Init parameters
    
    for(k in 1:maxit){
        for(i in 1:p) {
            iter = ifelse(k>1, yes=k-1, 1)
            theta[k,i] = t(x[,i])%*% (y - x[,-i]%*%theta[iter,-i])/ (t(x[,i])%*%x[,i])
            theta[k,i] = softthresh(theta[k,i], lambda)
        }
    }
    return(theta)
}

#compute 
maxiter = 20
out_gd = grad.descent(x, maxiter)
out_cd = coord.descent(x, maxiter, lambda=0)

#prepare results and plot
library(ggplot2)
out1 = data.frame(iter=1:maxiter ,p2 = out_cd[,2])
out2 = data.frame(iter=1:maxiter ,p2 = out_gd[,2])
anal = data.frame(iter=1:maxiter, p2 = analytical[2,])
ggplot(out1,aes(iter,p2)) + geom_line(aes(color="Coordinate descent"), size=1.5) +
    geom_line(data=out2,aes(color="Gradient descent", size=1)) + labs(color="") +
    geom_line(data=anal, aes(color="Analytical", size=1))
