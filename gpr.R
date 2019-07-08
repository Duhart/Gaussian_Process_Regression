# Libraries
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(nloptr)
library(kernlab)
library(plotly)
source('useful_functions.R')

# Standardising parameters
# Inputs: 
# x a vector of length n or a matrix of dimensions N x m
# Output:
# list with the following elements:
# mean: a vector of length m
# sd: a vector of length m
std_params <- function(x){
  if(is.null(dim(x))){
    m <- 1
    n <- length(x)
    res <- list(mean = mean(x), sd = sd(x))
  }else{
    m <- ncol(x)
    n <- nrow(x)
    res <- list(mean = x %>% apply(2,mean),
                sd = x %>% apply(2,sd))
  }
  return(res)
}

# Covariance kernel
# Inputs: x a vector of length p
# y a vector of length p
# l either a positive number or a vector of length p.
# sf a positive number.
# sn a non-negative number.
# Outputs:
# a number with the prior covariance of vectors x and y.
cov_ker <- function(x,y,l=1,sn=0,sf=1){
  return(sf^2*exp(-sum((x-y)/l)^2) + sn^2 * prod(x == y))
}

# matrix k
# Inputs:
# x: either a vector of length N or a matrix of dimensions N x m.
# y: either a vector of length n or a matrix of dimensions n x m.
# l a positive number.
# sf a positive number.
# sn a non-negative number.
# Outputs:
# a matrix of dimensions N times n
matrix_k <- function(x,y,l = 1, sn = 0, sf = 1){
  #browser()
  if(is.null(dim(x))){
    m <- 1
    N <- length(x)
    n <- length(y)
  }else{
    m <- ncol(x)
    N <- nrow(x)
    n <- nrow(y)
  }
  res <- matrix(rep(0,N*n),nrow = N)
  for(k in 1:N){
    for(j in 1:n){
      #browser()
      if(is.null(dim(x))){
        res[k,j] <- cov_ker(x[k],y[j],l,sn,sf)
      }else{
        res[k,j] <- cov_ker(x[k,],y[j,],l,sn,sf)
      }
    }
  }
  return(res)
}


# Inputs:
# x: a vector of length N or a matrix of dimensions N x m.
# X: a vector of length n or a matrix of dimensions n x m.
# y: a vector of length n.
# sigma: a number.
# l a positive number.
# sf a positive number.
# sn a non-negative number.
# prior_fun: a function from R^m to R that serves as the mean of f.
# ...: all parameters needed for the prior function.
# Outputs:
# A vector of length N 
mu <- function(x,X,y,sigma = 1, l = 1, sn = 0, sf = 1, prior_fun = NULL,...){
  if(is.null(prior_fun)){
    prior_fun <- function(p){return(mean(y))}
  }
  if(is.null(dim(X))){
    N <- length(x)
    n <- length(X)
    m <- 1
    prior_x <- 1:N %>% lapply(function(j){return(prior_fun(x[j],...))}) %>% unlist()
    prior_y <- 1:n %>% lapply(function(j){return(prior_fun(X[j],...))}) %>% unlist()
  }else{
    N <- nrow(x)
    n <- nrow(X)
    m <- ncol(X)
    prior_x <- 1:N %>% lapply(function(j){return(prior_fun(x[j,],...))}) %>% unlist()
    prior_y <- 1:n %>% lapply(function(j){return(prior_fun(X[j,],...))}) %>% unlist()
  }
  res <- rep(0,N)
  id <- diag(n)
  #browser()
  A <- solve(matrix_k(X,X,l,sn,sf) + sigma^2 * id,y - prior_y)
  if(is.null(dim(X))){
    res <- 1:N %>% lapply(function(j){
      aux <- matrix_k(x[j],X,l,sn,sf) %*% A
      return(prior_x[j] + aux[1,1])
    }) %>% unlist()
  }else{
    res <- 1:N %>% lapply(function(j){
      aux <- matrix_k(x[j,],X,l,sn,sf) %*% A
      return(prior_x[j] + aux[1,1])
    }) %>% unlist()
  }
  return(res)
}

mod_f <- mu(seq(0,100,length = 1000),x,y,sigma = sd(y), l = xs$sd, sf = sd(y),sn = 0)


# Inputs:
# x: a vector of length N or matrix of dimensions N x m.
# X: a vector of length n or matrix of dimensions n x m.
# sigma: a number.
# l: a number.
# Outputs:
# A vector of length N 
nu <- function(x,X,sigma = 1, l = 1, sn = 0, sf = 1){
  if(is.null(dim(X))){
    N <- length(x)
    n <- length(X)
    m <- 1
  }else{
    N <- nrow(x)
    n <- nrow(X)
    m <- ncol(X)
  }
  id <- diag(n)
  A <- solve(matrix_k(X,X,l,sn,sf) + sigma^2 * id)
  if(is.null(dim(X))){
    res <- 1:N %>% lapply(function(j){
      k <- cov_ker(x[j],x[j],l,sn,sf) # With current kernel this is always 1
      vec_k <- matrix_k(x[j],X,l,sn,sf)
      y <- k - vec_k %*% A %*% t(vec_k)
      return(y[1,1])
    }) %>% unlist()
  }else{
    res <- 1:N %>% lapply(function(j){
      k <- cov_ker(x[j,],x[j,],l,sn,sf) # With current kernel this is always 1
      vec_k <- matrix_k(x[j,],X,l,sn,sf)
      y <- k - vec_k %*% A %*% t(vec_k)
      return(y[1,1])
    }) %>% unlist()
  }
  return(pmax(0,res))
}


# Derivative of Covariance kernel
# Inputs: 
# x a vector of length p
# y a vector of length p
# l a positive number.
# sf a positive number.
# sn a non-negative number.
# Outputs:
# A vector of length p.
dk <- function(x,y,l=1,sn=0,sf=1){
  gk <- -2*cov_ker(x,y,l,sn,sf) * (x-y) / l^2
  return(gk)
}


lapply(0:100,function(x){return(cov_ker(x,50,l = xs$sd))}) %>% unlist() %>% plot(.,type='l')
lapply(0:100,function(x){return(dk(x,50,l = xs$sd))}) %>% unlist() %>% plot(.,type='l')


10*(cov_ker(x[1]+0.1,x[2],l = sd(x), sn = 0 , sf = 50) - cov_ker(x[1],x[2],l = sd(x), sn = 0 , sf = 50))
dk(x[1],x[2],l = sd(x), sn = 0 , sf = 50)


# The derivative of the matrix K
# Inputs:
# x: a vector of length m.
# X: a vector of length n or a matrix of dimensions n x m.
# l: a number.
# sf: a positive number.
# sn: a non-negative number.
# Outputs:
# A matrix of dimensions p x n.
dmatK <- function(x,X,l = 1, sn = 0, sf = 1){
  m <- length(x)
  if(is.null(dim(X))){
    n <- length(X)
  }else{
    n <- nrow(X)
  }
  res <- diag(0,nrow = m, ncol = n)
  X <- matrix(X,ncol = m)
  for(j in 1:n){
    #browser()
    res[,j] <- dk(x,X[j,],l,sn,sf)
  }
  return(res)
}

plot(seq(20,70,length = 100),as.vector(dmatK(50,seq(20,70,length = 100),l = xs$sd)),type='l')

dmatK(50,seq(20,70,by = 1),l = xs$sd)

dmatK(df$adstocks[1],df$adstocks,l=sd(df$adstocks),sn=0,sf=sd(df$y))
dmatK(df$adstocks[1],df$adstocks,l=sd(df$adstocks),sn=0,sf=sd(df$y))

# The derivative of mu
# Inputs:
# x: a vector of length N or a matrix of dimensions N x m.
# X: a vector of length n or a matrix of dimensions n x m.
# y: a vector of length n.
# sigma: a number.
# l: a number.
# sf: a positive number.
# sn: a non-negative number.
# prior_fun: a function from R^m to R that serves as the mean of f.
# dprior_fun: a function with the gradient of prior_fun. If prior_fun is not NULL
#             one has to send dprior_fun.
# ...: all parameters needed for prior_fun and dprior_fun.
# Output:
# A matrix of dimensions N x m
dmu <- function(x,X,y,sigma = 1, l = 1, sn = 0, sf = 1, prior_fun = NULL,dprior_fun = NULL, ...){
  n <- length(y)
  if(is.null(dim(X))){
    m <- 1
    N <- length(x)
  }else{
    m <- ncol(X)
    N <- nrow(X)
  }
  if(is.null(prior_fun)){
    prior_fun <- function(p){return(mean(y))}
    dprior_fun <- function(p){return(rep(0,m))}
  }
  if(is.null(dim(X))){
    #prior_x <- 1:N %>% lapply(function(j){return(prior_fun(x[j],...))}) %>% unlist()
    prior_y <- 1:n %>% lapply(function(j){return(prior_fun(X[j],...))}) %>% unlist()
    dprior_x <- 1:N %>% lapply(function(j){return(dprior_fun(x[j],...))})
  }else{
    #prior_x <- 1:N %>% lapply(function(j){return(prior_fun(x[j,],...))}) %>% unlist()
    prior_y <- 1:n %>% lapply(function(j){return(prior_fun(X[j,],...))}) %>% unlist()
    dprior_x <- 1:N %>% lapply(function(j){return(dprior_fun(x[j,],...))})
  }
  id <- diag(n)
  res <- matrix(rep(0,N * m),nrow = N)
  A <- solve(matrix_k(X,X,l,sn,sf) + sigma^2 * id, matrix(y,ncol = 1) - matrix(prior_y,ncol = 1))
  x <- matrix(x, ncol = m)
  #browser()
  for(j in 1:N){
    #browser()
    res[j,] <- as.numeric(matrix(dprior_x[[j]],ncol = 1) + 
               dmatK(x[j,],X,l,sn,sf) %*% A)
  }
  return(res)
}

mod_df <- dmu(seq(20,70,length = 1000),x,y,sigma = sd(y), l = 1, sf = sd(y),sn = 0)
plot(seq(20,70,length = 1000),mod_df,type='l')
abline(h=actual)



plot(alpha)


#### Test
# Data generation
set.seed(42)
actual <- 5
x <- rnorm(100,50,10)
y <- actual*x +10 + rnorm(100,0,10)
plot(x,y)

xs <- std_params(x)

cov_ker(30,40,l=xs$sd)


plot(seq(0,100,length.out = 1000),
     as.vector(matrix_k(0,seq(0,100,length.out = 1000),l=xs$sd, sf = 1,sn = 0)),
     type='l')

plot(seq(20,70,length.out = 1000),
     sqrt(nu(seq(20,70,length = 1000),x,sigma = 5, l=xs$sd, sf = sd(y),sn = 0)),
     type='l')

lines(seq(0,100,length.out = 1000),
     sqrt(nu(seq(0,100,length = 1000),x,sigma = 1, l=xs$sd, sf = sd(y),sn = 0)))

lines(seq(0,100,length.out = 1000),
      sqrt(nu(seq(0,100,length = 1000),x,sigma = 5, l=xs$sd, sf = 1,sn = 0)))

xs
mod_f <- mu(seq(0,100,length = 1000),x,y,sigma = 50, l=xs$sd, sf = 50,sn = 0,
            prior_fun = function(x){return(actual*x + 10)})
mod_v <- nu(seq(0,100,length = 1000),x,sigma = 50, l=xs$sd, sf = 50,sn = 0)

plot(x,y)
lines(seq(0,100,length = 1000),mod_f,col='red')
lines(seq(0,100,length = 1000),mod_f + 1.96 * sqrt(mod_v),col='red', lty = 2)
lines(seq(0,100,length = 1000),mod_f - 1.96 * sqrt(mod_v),col='red', lty = 2)


pred <- mu(x,x,y,sigma = 50, l=xs$sd, sf = 50,sn = 0,
           prior_fun = function(x){return(actual*x + 10)})
points(x,pred,col = 'blue')

plot(pred,y)
abline(0,1)

hist(((pred - y)-mean(pred - y))/sd(pred - y),prob = TRUE)

# This is showing a problem. Residuals are smaller as y increases. They should be independent.
plot(y,((pred - y)-mean(pred - y))/sd(pred - y))



mod_df <- dmu(seq(0,100,length = 1000),x,y,sigma = 50, l=xs$sd, sf = 50,sn = 0,
              prior_fun = function(x){return(actual*x+10)},
              dprior_fun = function(x){return(actual)})

plot(seq(0,100,length = 1000),mod_df,type='l')
abline(h=actual)

points(x,dmu(x,x,y,sigma = 50, l=xs$sd, sf = 50,sn = 0,
             prior_fun = function(x){return(actual*x+10)},
             dprior_fun = function(x){return(actual)}),col = 'red')


dmu(x[1],x,y,sigma = 50, l=xs$sd, sf = 50,sn = 0,
    prior_fun = function(x){return(actual*x+10)},
    dprior_fun = function(x){return(actual)})

plot(seq(20,70,by=1),
     -1 * (mu(seq(20,70,by=1),x,y,sigma = 50, l=xs$sd, sf = 50,sn = 0,
   prior_fun = function(x){return(actual*x + 10)}) -
  mu(1 + seq(20,70,by=1),x,y,sigma = 50, l=xs$sd, sf = 50,sn = 0,
     prior_fun = function(x){return(actual*x + 10)})),type='l',ylab = 'slope')
abline(h=actual)


### Test with toy model

df <- read.csv('https://raw.githubusercontent.com/Duhart/Gaussian_Process_Regression/master/data.csv')
names(df) <- c('Date','Awareness','GRPs')

params <- finding_params(logit(df$Awareness),df$GRPs,c(0,0.999),max_theta = 6,grid_size = 1000)
df$adstocks <- adstocks_lag(df$GRPs,params$rate,params$theta)
df$y <- logit(df$Awareness)


head(df)

df$gausspr <-  sigmoid(mu(df$adstocks,df$adstocks, df$y,
                          sigma = 1,
                          l = sd(df$adstocks),
                          sf = 1,
                          sn = 0))

df$LB <-  sigmoid(mu(df$adstocks,df$adstocks, df$y,sigma = sd(df$y),
                     l = sd(df$adstocks),
                     sf = sd(df$y),
                     sn = 0) - 1.96 * 
                  sqrt(nu(df$adstocks,df$adstocks,sigma = sd(df$y),
                     l = sd(df$adstocks),
                     sf = sd(df$y),
                     sn = 0)))

df$UB <- sigmoid(mu(df$adstocks,df$adstocks, df$y,sigma = sd(df$y),
                    l = sd(df$adstocks),
                    sf = sd(df$y),
                    sn = 0) + 1.96 * 
                   sqrt(nu(df$adstocks,df$adstocks,sigma = sd(df$y),
                           l = sd(df$adstocks),
                           sf = sd(df$y),
                           sn = 0)))

df$Base <- sigmoid(mu(0,df$adstocks, df$y,sigma = sd(df$y),
                      l = sd(df$adstocks),
                      sf = sd(df$y),
                      sn = 0))


df %>% ggplot(aes(x = Date, y = Awareness)) + geom_line() + 
  geom_line(aes(x = df$Date, y = df$gausspr),col = 'red') +
  geom_line(aes(x = df$Date, y = df$Base),col = 'red', lty = 2) +
  geom_ribbon(aes(x = df$Date, ymin=df$LB,ymax=df$UB), fill="red", alpha="0.2") 



#### Saturation curve

logit_mean <- mu(seq(0,2500, length = 1000),df$adstocks, df$y,sigma = sd(df$y),
                 l = sd(df$adstocks),
                 sf = sd(df$y),
                 sn = 0)
zero_logit <- mu(0,df$adstocks, df$y,sigma = sd(df$y),
                 l = sd(df$adstocks),
                 sf = sd(df$y),
                 sn = 0)
ads_seq <- seq(0,2500, length = 1000)


plot(ads_seq,logit_mean,type= 'l')
abline(h=zero_logit)
plot(ads_seq,logit_mean - zero_logit,type= 'l')


summary(sigmoid(logit_mean - zero_logit))

summary(sigmoid(logit_mean) - sigmoid(zero_logit))

plot(df$adstocks,df$Awareness - sigmoid(zero_logit),
     xlab = 'adstocks',
     ylab = 'Contribution',
     main = 'Saturation curve')
lines(ads_seq,sigmoid(logit_mean) - sigmoid(zero_logit),
     type='l', col = 'red')

points(df$adstocks,df$gausspr - sigmoid(zero_logit))
points(df$adstocks,df$Awareness - sigmoid(zero_logit))


plot(df$adstocks,df$gausspr - sigmoid(zero_logit))


lines(seq(0,2500, length = 1000),
      sigmoid(mu(seq(0,2500, length = 1000),df$adstocks, df$y,sigma = sd(df$y),
                 l = sd(df$adstocks),
                 sf = sd(df$y),
                 sn = 0) - 1.96 * 
              sqrt(nu(seq(0,2500, length = 1000),df$adstocks,sigma = sd(df$y),
                      l = sd(df$adstocks),
                      sf = sd(df$y),
                      sn = 0))), lty=2, col = 'red')

lines(seq(0,2500, length = 1000),
      sigmoid(mu(seq(0,2500, length = 1000),df$adstocks, df$y,sigma = sd(df$y),
                 l = sd(df$adstocks),
                 sf = sd(df$y),
                 sn = 0) + 1.96 * 
                sqrt(nu(seq(0,2500, length = 1000),df$adstocks,sigma = sd(df$y),
                        l = sd(df$adstocks),
                        sf = sd(df$y),
                        sn = 0))), lty=2, col = 'red')


mod_df <- dmu(df$adstocks,df$adstocks, df$y,sigma = sd(df$y),
              l = sd(df$adstocks),
              sf = sd(df$y),
              sn = 0)
plot(df$Date,mod_df,type='l')


plot(df$Date,mod_df * df$adstocks,type='l')


par(mfrow = c(3,1))

plot(df$Date,df$Awareness,type='l', main = 'Raw, model, and base')
lines(df$Date,df$gausspr,col = 'red')
polygon(c(df$Date,rev(df$Date)),c(df$LB,rev(df$UB)),col = rgb(1,0,0,0.2))
lines(df$Date,df$Base,col = 'red',lty = 2)


plot(df$Date,df$adstocks,type='h', main = 'Media and adstocks')
lines(df$Date,df$GRPs,lwd = 2)


plot(df$Date,mod_df,type='l', col = 'red', main = 'Derivative of logit(tbca) wrt to adstocks')




### Simplest model ever!!!
set.seed(42)
x <- 1:10
y <- 2*x + 4 + rnorm(10)
res <- y-2*x-4
hist(res)
sd(res)

mod <- mu(seq(-3,15,length = 100),x,y,sigma = 1, l = sd(x), sf = sd(res),
          prior_fun = function(x){return(2*x + 4)})
mod_v <- nu(seq(-3,15,length = 100),x,sigma = 1, l = sd(x), sf = sd(res))

par(mfrow = c(2,1))
plot(x,y,xlim = c(-3,15),ylim = c(0,30))
lines(seq(-3,15,length = 100),mod,col='red')
lines(seq(-3,15,length = 100),mod + 1.95*sqrt(mod_v),col='red',lty=2)
lines(seq(-3,15,length = 100),mod - 1.95*sqrt(mod_v),col='red',lty=2)

mod_dy <- dmu(seq(0,11,length = 100),x,y,sigma = 1,l = sd(x), sf = sd(res),
              prior_fun = function(x){return(2*x + 4)},
              dprior_fun = function(x){return(2)})
plot(seq(0,11,length = 100),mod_dy,type = 'l',col = 'red')
#abline(h = 2)

# slope test

x1 <- 0
x2 <- 10
y1 <- mu(0,x,y,sigma = 1, l = sd(x), sf = sd(y-2*x-4),
         prior_fun = function(x){return(2*x + 4)})
y2 <- mu(10,x,y,sigma = 1, l = sd(x), sf = sd(y-2*x-4),
         prior_fun = function(x){return(2*x + 4)})

(y2-y1)/(x2-x1)

s <- dmu(c(0,5,10),x,y,sigma = 1,l = sd(x), sf = sd(res),
          prior_fun = function(x){return(2*x + 4)},
          dprior_fun = function(x){return(2)})

# Multivarate regression test
# Data generation
set.seed(42)
N <- 30
x <- rnorm(N,5,1)
y <- rnorm(N,10,2)
plot(x,y,asp=1)
# Very weird function
# Note there is an interaccion term x * y
#z <- 6 + 2*x  
#z <- 6 + 2*(x-5)^2 + rnorm(N) 
z <- 6 + 2*(x-5)^2 + (y-10)^2 + rnorm(N) 
plot(x,z)
plot(y,z)

dat <- data.frame(x = x,y = y, z = z)
head(dat)

# 3D plot with plotly
plot_ly(dat, x = ~ x,
        y = ~ y,
        z = ~ z,
        marker = list(color = ~ z,
                      showscale = TRUE)) %>% 
  layout(scene = list(xaxis = list(title = 'x'),
                      yaxis = list(title = 'y'),
                      zaxis = list(title = 'z'))) %>% 
  add_markers()

X <- dat %>% select(x,y)


mod <- mu(X,X,dat$z,
          sigma = sd(dat$z), sf = sd(dat$z), l = std_params(cbind(dat$x,dat$y))$sd)
modv <- nu(X,X,
           sigma = sd(dat$z), sf = sd(dat$z), l = std_params(cbind(dat$x,dat$y))$sd)

UB <- mod + 1.96 * sqrt(modv)
LB <- mod - 1.96 * sqrt(modv)

mod_df <- data.frame(z = dat$z, mean = mod,LB = LB, UB = UB)
mod_df <- mod_df[order(z),]

plot(dat$z,mod,ylim = c(0,25))
#lines(mod_df$z,mod_df$mean,col = 'red')
polygon(x = c(mod_df$z,rev(mod_df$z)), y = c(mod_df$LB,rev(mod_df$UB)), col = rgb(1,0,0,0.2))
abline(0,1)

# Contour plots of this model.

test <- expand.grid(list(x = seq(min(x),max(x),length = N),
                      y = seq(min(y),max(y),length = N)))
test$z <- mu(test,X,dat$z,
          sigma = 1, sf = sd(dat$z), l = std_params(cbind(dat$x,dat$y))$sd)
head(test)

# Do the plot
image(seq(min(x),max(x),length = N),
      seq(min(y),max(y),length = N),
      matrix(test$z,ncol = N,byrow = FALSE),
      col = rev(heat.colors(15)),
      xlab = 'x',ylab = 'y', main = 'Contour plot of modelled function',
      xlim = c(min(x),max(x)), ylim = c(min(y),max(y)))

points(dat$x,dat$y,pch=15)

contour(x = seq(min(x),max(x),length = N),
        y = seq(min(y),max(y),length = N),
        z = matrix(test$z,ncol = N,byrow = FALSE),
        add = TRUE)

plot_ly(x = seq(min(x),max(x),length = N),
        y = seq(min(y),max(y),length = N),
        z = matrix(test$z,ncol = N,byrow = FALSE),
        type = 'contour')

## Adding a linear model to help

dat$xx <- 1:N %>% lapply(function(k){
  return(sum((dat[k,1:2] - c(5,9))^2))}) %>% unlist()
  
plot(dat$xx,dat$z)

lm2 <- lm(z~xx,dat)
summary(lm2)

predict(lm2,0)

pf <- function(x){
  r <- sum((x - c(5,9))^2)
  return(as.numeric(lm2$coefficients[1] + lm2$coefficients[2]*r))
}

dpf <- function(x){
  return(as.numeric(2*lm2$coefficients[2]*(x - c(5,9))))
}



mod <- mu(X,X,dat$z,
          sigma = 10, sf = 10, l = std_params(cbind(dat$x,dat$y))$sd,
          prior_fun = pf)
modv <- nu(X,X,
           sigma = 10, sf = 10, l = std_params(cbind(dat$x,dat$y))$sd)

UB <- mod + 1.96 * sqrt(modv)
LB <- mod - 1.96 * sqrt(modv)

mod_df <- data.frame(z = dat$z, mean = mod,LB = LB, UB = UB)
mod_df <- mod_df[order(z),]


par(mfrow = c(1,1))
plot(dat$z,mod, main = 'Actual vs model')
#lines(mod_df$z,mod_df$mean,col = 'red')
polygon(x = c(mod_df$z,rev(mod_df$z)), y = c(mod_df$LB,rev(mod_df$UB)), col = rgb(1,0,0,0.2))
abline(0,1)

# Contour plots of this model.

test <- expand.grid(list(x = seq(min(x),max(x),length = N),
                         y = seq(min(y),max(y),length = N)))
test$z <- mu(test[,1:2],X,dat$z,
             sigma = 10, sf = 10, l = std_params(cbind(dat$x,dat$y))$sd,
             prior_fun = pf)
head(test)

# Do the plot
image(seq(min(x),max(x),length = N),
      seq(min(y),max(y),length = N),
      matrix(test$z,ncol = N,byrow = FALSE),
      col = rev(heat.colors(15)),
      xlab = 'x',ylab = 'y', main = 'Contour plot of modelled function',
      xlim = c(min(x),max(x)), ylim = c(min(y),max(y)))

points(dat$x,dat$y,pch=15)

contour(x = seq(min(x),max(x),length = N),
        y = seq(min(y),max(y),length = N),
        z = matrix(test$z,ncol = N,byrow = FALSE),
        add = TRUE)


## USA data

df <- read.csv("C:/Users/gonzalezh/Kantar/Schubert, Jan (MBWAR) - M2M Simulator/ABI data/data/media_tbca_us.csv")
head(df)
df$Date <- as.Date(df$To.Date,'%Y-%m-%d')
df <- df %>% filter(!is.na(Week))
df <- df %>% filter(!is.na(media))

View(df)
plot(df$Date,df$media,type='l')
plot(df$Date,df$tbca_imputed,type='l')
lines(df$Date,df$tbca_smooth,col = 'red')

# Adstocks


params <- finding_params(logit(df$tbca_smooth/100),df$media,c(0,0.999),max_theta = 6,grid_size = 1000)
df$adstocks <- adstocks_lag(df$media,params$rate,params$theta)
df$y <- logit(df$tbca_smooth/100)

head(df)

df$gausspr <-  sigmoid(mu(df$adstocks,df$adstocks, df$y,
                          sigma = sd(df$y),
                          l = sd(df$adstocks),
                          sf = sd(df$y),
                          sn = 0))

df$LB <-  sigmoid(mu(df$adstocks,df$adstocks, df$y,sigma = sd(df$y),
                     l = sd(df$adstocks),
                     sf = sd(df$y),
                     sn = 0) - 1.96 * 
                    sqrt(nu(df$adstocks,df$adstocks,sigma = sd(df$y),
                            l = sd(df$adstocks),
                            sf = sd(df$y),
                            sn = 0)))

df$UB <- sigmoid(mu(df$adstocks,df$adstocks, df$y,sigma = sd(df$y),
                    l = sd(df$adstocks),
                    sf = sd(df$y),
                    sn = 0) + 1.96 * 
                   sqrt(nu(df$adstocks,df$adstocks,sigma = sd(df$y),
                           l = sd(df$adstocks),
                           sf = sd(df$y),
                           sn = 0)))

df$Base <- sigmoid(mu(0,df$adstocks, df$y,sigma = sd(df$y),
                      l = sd(df$adstocks),
                      sf = sd(df$y),
                      sn = 0))

df$slope <- dmu(df$adstocks,df$adstocks,df$y,
                sigma = sd(df$y),
                l = sd(df$adstocks),
                sf = sd(df$y),
                sn = 0)

zero <- sigmoid(mu(0,df$adstocks,df$y,
                   sigma = sd(df$y),
                   l = sd(df$adstocks),
                   sf = sd(df$y),
                   sn = 0))




df %>% ggplot(aes(x = Date, y = tbca_imputed/100)) + geom_line() + 
  geom_line(aes(x = df$Date, y = df$tbca_smooth/100),col = 'blue') +
  geom_line(aes(x = df$Date, y = df$gausspr),col = 'red') +
  geom_line(aes(x = df$Date, y = df$Base),col = 'red', lty = 2) +
  geom_ribbon(aes(x = df$Date, ymin=df$LB,ymax=df$UB), fill="red", alpha="0.2") 


par(mfrow = c(3,1))

plot(df$Date,df$tbca_imputed/100,type='l', main = 'Raw, smooth, model, and base')
lines(df$Date,df$tbca_smooth/100,type='l',lwd = 2)
lines(df$Date,df$gausspr,col = 'red')
polygon(c(df$Date,rev(df$Date)),c(df$LB,rev(df$UB)),col = rgb(1,0,0,0.2))
abline(h=zero,col = 'red',lty = 2)

plot(df$Date,df$adstocks,type='h', main = 'Media and adstocks')
lines(df$Date,df$media,lwd = 2)

plot(df$Date,df$slope,type='l', main = 'Derivative of logit(tbca) wrt to adstocks')




lines(df$Date,1.5e07 + 3e06*(df$slope-min(df$slope))/(max(df$slope)-min(df$slope)),col = 'red')
lines(df$Date,1.5e07 + 3e06*(df$slope-min(df$slope))/(max(df$slope)-min(df$slope)),col = 'blue')

plot(df$Date,df$slope,type='l',col ='red')
plot(df$Date,df$tbca_smooth/100 * (1 - df$tbca_smooth/100) * df$slope,type='l',col ='red')
plot(df$Date,df$adstocks * (1 - df$tbca_smooth/100) * df$slope,type='l',col ='red')

plot(df$Date,scale(df$tbca_smooth),type='l',col = 'blue',ylim = c(-3,5))
lines(df$Date,scale(df$adstocks))
lines(df$Date,scale(df$media),lty=2)
lines(df$Date,scale(df$slope),col = 'red')


df %>% ggplot(aes(x = Date, y = slope)) + geom_line(col = 'red')
df %>% ggplot(aes(x = Date, y = slope)) + geom_line(col = 'red')


x_ads <- seq(0,max(df$adstocks),length.out = 1000)
pred <- sigmoid(mu(x_ads,df$adstocks,df$y,
                   sigma = sd(df$y),
                   l = sd(df$adstocks),
                   sf = sd(df$y),
                   sn = 0))
dpred <- dmu(x_ads,df$adstocks,df$y,
             sigma = sd(df$y),
             l = sd(df$adstocks),
             sf = sd(df$y),
             sn = 0)


plot(df$adstocks,df$tbca_imputed)
plot(df$tbca_smooth/100,df$gausspr)
abline(0,1)
plot(x_ads,pred,type='l',col = 'red')
plot(x_ads,dpred,type='l',col = 'red')


