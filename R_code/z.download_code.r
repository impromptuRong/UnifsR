###########################
#lecture 6 code
#ozone data
#comparison: sparse regression methods

#libraries required
library(glmnet)
library(ncvreg)

#load in ozone data 
ozone = read.csv("ozone.csv")
Y = as.numeric(ozone[,1]); Y = Y - mean(Y)
X = as.matrix(ozone[,-1]); X = scale(X,center=T,scale=F)

#Lasso
lam = 1
fitl = glmnet(x=X,y=Y,family="gaussian",lambda=lam,alpha=1)
cbind(fit0$coef,as.matrix(fitl$beta))

#Lasso paths
fitl = glmnet(x=X,y=Y,family="gaussian",alpha=1)
plot(fitl,col=1:8)
legend(0,19,legend=names(ozone)[2:9],col=1:8,lty=rep(1,8),cex=.8)



###############################
#least squares, lasso, adaptive lasso, SCAD, ridge, elastic net, MC+
lam = 1

betals = solve(t(X)%*%X)%*%t(X)%*%Y
betar = solve(t(X)%*%X + diag(rep(lam/2*nrow(ozone),8)))%*%t(X)%*%Y
fitl = glmnet(x=X,y=Y,family="gaussian",lambda=lam,alpha=1)
fital = glmnet(x=X,y=Y,family="gaussian",lambda=lam,alpha=1,penalty.factor=1/abs(betals))
fitel = glmnet(x=X,y=Y,family="gaussian",lambda=lam,alpha=.5)
fitscad = ncvreg(X,Y,family="gaussian",penalty="SCAD",lambda=lam)
fitmcp = ncvreg(X,Y,family="gaussian",penalty="MCP",lambda=lam)
mat = cbind(betals,betar,as.matrix(fitl$beta),as.matrix(fital$beta),as.matrix(fitel$beta),fitscad$beta[-1],fitmcp$beta[-1])
colnames(mat) = c("LS","Ridge","Lasso","A-Lasso","EL","SCAD","MC+")
mat


#############################
#compare ridge, lasso, elastic net & SCAD regualrization paths

par(mfrow=c(2,3))
par(mar=c(5,4,3,2))
betals = solve(t(X)%*%X)%*%t(X)%*%Y
lambdas = exp(seq(log(.01),log(100*nrow(ozone)),l=100))
betasr = matrix(0,length(lambdas),8)
for(i in 1:length(lambdas))
{
  betasr[i,] = solve(t(X)%*%X + diag(rep(lambdas[i],8)))%*%t(X)%*%Y
}
plot(c(1,length(lambdas)),range(betals),type="n",ylab="Coefficients",xlab="Lambda Index",main="Ridge")
for(j in 1:8)
{
  lines(betasr[length(lambdas):1,j],col=j)
}
legend(0,20,legend=names(ozone)[2:9],col=1:9,lty=rep(1,9),cex=.75)

fitl = glmnet(x=X,y=Y,family="gaussian",alpha=1)
plot(fitl,col=1:8,main="Lasso")
legend(0,20,legend=names(ozone)[2:9],col=1:8,lty=rep(1,8),cex=.75)

fitel = glmnet(x=X,y=Y,family="gaussian",alpha=.5)
plot(fitel,col=1:8,main="EL alpha=.5")
legend(0,20,legend=names(ozone)[2:9],col=1:8,lty=rep(1,8),cex=.75)

fitel = glmnet(x=X,y=Y,family="gaussian",alpha=.25)
plot(fitel,col=1:8,main="EL alpha=.25")
legend(0,20,legend=names(ozone)[2:9],col=1:8,lty=rep(1,8),cex=.75)

fitscad = ncvreg(X,Y,family="gaussian",penalty="SCAD")
plot(fitscad,col=1:8,main="SCAD",shade=F)
legend(6,30,legend=names(ozone)[2:9],col=1:8,lty=rep(1,8),cex=.75)

fitmcp = ncvreg(X,Y,family="gaussian",penalty="MCP")
plot(fitmcp,col=1:8,main="MC+",shade=F)
legend(6,30,legend=names(ozone)[2:9],col=1:8,lty=rep(1,8),cex=.75)



