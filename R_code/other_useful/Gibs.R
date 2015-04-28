#### Markov Chain ####
N = 10000
signal = vector(length = N)
signal[1] = 0
for(i in 2:N){
# random select one offset (from [-1,1]) to signal[i-1]
    signal[i] = signal[i-1] + sample(c(-1,1),1)
}
plot( signal,type = 'l',col = 'red')

#### MCMC ####
N = 10000
x = vector(length = N)
x[1] = 0

# uniform variable: u
u = runif(N)  
m_sd = 5  
freedom = 5  

for (i in 2:N){
    y = rnorm(1,mean = x[i-1],sd = m_sd)
    print(y)
    rt(1,df = freedom)
    p_accept = dnorm(x[i-1],mean = y,sd = abs(2*y+1)) / dnorm(y, mean = x[i-1],sd = abs(2*x[i-1]+1))
    # print (p_accept)
    
    if ((u[i] <= p_accept)){
        x[i] = y
        print("accept")
    } else {
        x[i] = x[i-1]
        print("reject")
    }
}  

plot(x,type = 'l')
dev.new()
hist(x)

#### Gibs Sampler ####
p_ygivenx <- function(x,m1,m2,s1,s2)
{
    return (rnorm(1,m2+rho*s2/s1*(x-m1),sqrt(1-rho^2)*s2 ))
}

p_xgiveny <- function(y,m1,m2,s1,s2)
{
    return 	(rnorm(1,m1+rho*s1/s2*(y-m2),sqrt(1-rho^2)*s1 ))
}


N = 5000
K = 20 #iteration in each sampling
x_res = vector(length = N)
y_res = vector(length = N)
m1 = 10; m2 = -5; s1 = 5; s2 = 2
rho = 0.5
y = m2

for (i in 1:N)
{
    x = p_xgiveny(y, m1,m2,s1,s2)
    y = p_ygivenx(x, m1,m2,s1,s2)
    # print(x)
    x_res[i] = x;
    y_res[i] = y;
}

hist(x_res,freq = 1)
dev.new()
plot(x_res,y_res)
library(MASS)
valid_range = seq(from = N/2, to = N, by = 1)
MVN.kdensity <- kde2d(x_res[valid_range], y_res[valid_range], h = 10)   #estimate kernel density
plot(x_res[valid_range], y_res[valid_range], col = "blue", xlab = "x", ylab = "y")
contour(MVN.kdensity, add = TRUE)   #2d norm distribution contour plot

#real distribution
# real = mvrnorm(N,c(m1,m2),diag(c(s1,s2)))
# dev.new()
# plot(real[1:N,1],real[1:N,2])


