# This function generates a random walk of length N,
# starting at level, with drift drift and variance sigma^2.
# It simulates a path at integer times of the SDE
# dX_t = level + drift * t + sigma * dB_t.
# Inputs:
# N: length of the path.
# level: Value for X_0.
# drift: Slope of E[X_t].
# sigma: Variance of X_t/sqrt(t)
# Output:
# B: a vector of length N with the observed path.
bm <- function(N,level,drift,sigma){
  B <- rep(level,N)
  for(k in 2:N){
    B[k] <- B[k-1] + drift + sigma*rnorm(1)
  }
  return(B)
}

# This function generates the usual adstocks
# a[t+1] = x[t+1] + r * a[t].
# Inputs:
# x: a non-negative vector with the values of investment or GRPs
# rate: a value in [0,1). This parameter handles how long the effect of media lasts.
# init: a non-negative number. This parameter handles how the adstocks initiate.
#       It's like asssuming there was some investment before for which we don't
#       have the information.
# Outputs:
# a: a vector the same lenght of x with the transformed adstocks.
adstocks_geom <- function(x,rate,init = NULL){
  n <- length(x)
  if(is.null(init)){
    init <- x[1]
  }
  a <- x
  a[1] <- init
  for(k in 2:n){
    a[k] <- x[k] + rate*a[k-1]
  }
  return(a)
}


# This function calculates the lagged adstocks
# Inputs:
# x: a non-negative vector with the values of investment or GRPs
# rate: a value in [0,1). This parameter handles how long the effect of media lasts.
# theta: a non-negative number, most likely an integer. This parameter handles
#        how long one has to wait to see the maximum effect of media.
# Outputs:
# a: a vector the same lenght of x with the transformed adstocks.
adstocks_lag <- function(x,rate,theta){
  n <- length(x)
  a <- rep(0,n)
  w <- rate^((c(0:(n-1)) - theta)^2)
  for(k in 1:n){
    a <- a + c(rep(0,k-1),x[k]*w)[1:n]
  }
  return(a)
}

# This function calculates the best parameters
# according to the correlation with a variable y.
# Inputs:
# y: a numeric vector with which we want to correlate our adstocks.
# x: a non-negative vector with values of investment or GRPs.
# range_rate: a numeric vector of two components in [0,1).
#             The first component should be smaller than the second.
#             The adstocks will be tested in values within these two.
# max_theta: a non-negative integer. We will test adstocks for all lags
#            from 0 to max_theta.
# grid_size: a natural number larger or equal to two.
#            It is the number of values generated within range_rate.
# Outputs: a list with the following three elements.
# rate: The rate of the adstocks of x that maximises the correlation with y.
# theta: The lag of the adstocks of x that maximises the correlation with y.
# cor: The maximised correlation between the adstocks of x and y.
finding_params <- function(y,x,range_rate,max_theta,grid_size = 100){
  n <- length(x)
  rates <- seq(range_rate[1],range_rate[2],length.out = grid_size)
  best_rate <- 0
  best_theta <- 0
  best_cor <- cor(x,y)
  for(th in 0:max_theta){
    cors <- rates %>% lapply(function(r){
      return(cor(y,adstocks_lag(x,r,th)))
    }) %>% unlist()
    if(max(cors) > best_cor){
      best_cor <- max(cors)
      best_rate <- rates[which(cors == max(cors))]
      best_theta <- th
    }
  }
  return(list(rate = best_rate, theta = best_theta, cor = best_cor))
}

# This function calculates the logit of a numeric vector
# Input: x is a numeric vector
# Output: A vector of the same length of x with its logit values.
logit <- function(x){
  y <- log(x/(1-x))
  y[x<0.0001] <- log(0.0001/0.9999)
  y[x>0.9999] <- log(0.9999/0.0001)
  return(y)
}

# This function gives back the anti-logit of a numeric vector
# Input: x is a numeric vector
# Outpout: A vector of the same length of x with its anti-logit values.
sigmoid <- function(x){
  return(1/(1+exp(-x)))
}

# This function calculates independent linear regressions y~x
# in differnt breaks given by the user.
# Inputs:
# y: the dependent variable
# x: the covariate. Right now, it only works with just one non-constant covariate.
# br: vector with the indices where there's a structural change in y.
# Output: a list with the following elements:
# fitted: a vector with the estimation of y.
# parameters: a data frame with the regression coefficients for each break.
ind_reg <- function(y,x,br){
  N <- length(y)
  n <- length(br)+1
  init <- c(1,br)
  final <- c(br,N)
  y_hat <- y
  params <- c(1:n) %>% lapply(function(k){
    mod <- lm(y[init[k]:final[k]] ~ x[init[k]:final[k]])
    y_hat[init[k]:final[k]] <<- mod$fitted.values
    return(data.frame(intercept = mod$coefficients[1],slope = mod$coefficients[2]))
  }) %>% do.call(rbind,.)
  return(list(fitted = y_hat, parameters = params))
}

# This function repeats the values at vals a number of times given by br to get a
# vector of size n
rep_br <- function(vals,br,n){
  m <- length(vals)
  res <- rep(NA,n)
  for(k in 1:m){
    if(k == 1){
      res[1:br[k]] <- vals[1]
    }else{
      if(k == m){
        res[br[m-1]:n] <- vals[m]
      }else{
        res[br[k-1]:br[k]] <- vals[k]
      }
    }
  }
  return(res)
}

# This function does a DLM with fixed slopes between breakpoints and
# smooth intercept along the whole period.
# Inputs: - y: Dependent variable as numeric vector.
#         - x: Independent variable as numeric vector. (Only one variable!)
#         - br: Numeric vector indicating where are the breaks.
#         - inter: Real number indicating the initial value of the intercept.
#         - sl: Real number indicating the initial value of the slope. 
# Output: A data frame with the following columns:
#         - x: The covariate.
#         - y: The actual data to model.
#         - y_hat: The modelled values of y.
#         - residual: The differnce of y_hat - y.
#         - intercept: The values of the intercept.
#         - slope: The values of the slope.
mixed <- function(y,x,br,inter = NULL, sl = NULL){
  N <- length(y)
  n <- length(br) + 1
  tab <- data.frame(segment = 1:n,
                    init = c(1,br+1),
                    final = c(br,N))
  tab$n <- tab$final - tab$init + 1
  tab$intercept <- lapply(1:n,function(k){
    return(lm(y[tab$init[k]:tab$final[k]] ~ x[tab$init[k]:tab$final[k]])$coefficients[1])
  }) %>% unlist()
  tab$slope <- lapply(1:n,function(k){
    return(lm(y[tab$init[k]:tab$final[k]] ~ x[tab$init[k]:tab$final[k]])$coefficients[2])
  }) %>% unlist()
  
  # Set hyper parameters
  ols <- lm(y ~ x)
  m0 <- c(min(tab$intercept),
          max(0,ols$coefficients[2],tab$slope[1]))
  if(!is.null(inter)) m0[1] <- inter
  if(!is.null(sl)) m0[2] <- sl
  
  C0 <- var(y) * solve(t(cbind(1,x)) %*% cbind(1,x))
  
  # Declare outputs
  m <- vector('list',length=N)
  C <- vector('list',length=N)
  a <- vector('list',length=N)
  R <- vector('list',length=N)
  f <- rep(0,N)
  Q <- rep(0,N)
  s <- vector('list',length=N)
  S <- vector('list',length=N)
  d <- list(c(0.9,0.9))
  W <- vector('list',length=N)
  
  # Estimate filter
  for(k in 1:n){
    for(j in tab$init[k]:tab$final[k]){
      if(j == 1){
        a[[j]] <- m0
        W[[j]] <- 0*C0
        R[[j]] <- C0
      }else{
        a[[j]] <- m[[j-1]]
        W[[j]] <- diag(sqrt(1-d[[j]])/d[[j]]) %*% C[[j-1]] %*% diag(sqrt(1-d[[j]])/d[[j]])
        R[[j]] <- C[[j-1]] + W[[j]]
      }
      f[j] <- sum(a[[j]] * c(1,x[j]))
      Q[j] <- as.numeric(cbind(1,x[j]) %*% R[[j]] %*% rbind(1,x[j])) + var(y)
      m[[j]] <- a[[j]] + (y[j] - f[j]) / Q[j] * R[[j]] %*% rbind(1,x[j]) 
      C[[j]] <- R[[j]] - R[[j]] %*% rbind(1,x[j]) %*% cbind(1,x[j]) %*% R[[j]] / Q[j]
      #d[[j + 1]] <- c(0.98,min(1,d[[j]][2] + 0.02))
      d[[j + 1]] <- c(0.98,1)
    }
    d[[j + 1]] <- c(0.9,0.9)
  }
  
  # Estimate smoothing
  s[[N]] <- m[[N]]
  S[[N]] <- C[[N]]
  for(k in c((N-1):1)){
    s[[k]] <- m[[k]] + C[[k]] %*% solve(R[[k+1]],s[[k+1]] - a[[k+1]])
    S[[k]] <- C[[k]] - C[[k]] %*% solve(R[[k+1]]) %*%
      (R[[k+1]] - S[[k+1]]) %*% solve(R[[k+1]]) %*% C[[k]]
  }
  
  intercept <- s %>% lapply(function(x){ return(x[1])}) %>% unlist()
  slope <- s %>% lapply(function(x){ return(x[2])}) %>% unlist()
  y_hat <- intercept + slope * x 
  residual <- y_hat - y
  
  return(data.frame(x = x,
                    y = y,
                    y_hat = y_hat,
                    residual = residual,
                    intercept = intercept,
                    slope = slope))
}

# This function is to find the breaks. With the idea of finding differnt campaigns after following 
# certain specific rules.
# Input: - spend: a numeric non-negative vector.
# Output: a list with 3 elements. break_index is the one with the indices of the breaks.
find_camp_breaks <- function(spend, h_nospend = 3, h_spikes = 4, spike_sensitivity = .95, correction = FALSE) {
  n <- length(spend)
  no_spend <- rep(FALSE,n)
  
  #find length of sequence of 0s 
  r <- rle(spend)
  r.index <- rep(r$lengths >= h_nospend & r$values == 0,r$lengths)
  
  #create TRUE at the end of and before a sequence of 0s longer than 'h'
  for(i in 1:(n-1)){
    if(r.index[i] == TRUE & r.index[i+1] == FALSE) no_spend[i] <- TRUE
    if(r.index[i] == FALSE & r.index[i+1] == TRUE) no_spend[i] <- TRUE
  }
  
  #create index
  no_spend_index <- c(1:n)[no_spend]
  
  #bonferonni correction
  n_comparisons <- length(h_spikes:(n-1))
  a <- 1-spike_sensitivity
  spike_sensitivity_adj <- 1-(a/(n_comparisons-1))
  
  if(correction) spike_sensitivity <- spike_sensitivity_adj
  
  #check for spikes
  spike <- rep(FALSE,n)
  for(i in h_spikes:(n-1)){
    x <- spend[(h_spikes-1):i]
    m <- mean(x)
    se <- sd(x)/sqrt(4)
    t <- qt(spike_sensitivity,df = h_spikes-1)
    
    up <- m + t*se
    spike[i] <- spend[i+1] > up
    if(any(spike[(i-h_spikes):(i-1)])) spike[i] <- FALSE
  }
  spike_index <- c(1:n)[spike]
  
  #check that spike_break is not 1 to 4 behind/forward no_spend spike
  for(i in 1:4){
    spike_index <- spike_index[!((spike_index-i) %in% no_spend_index)]
    spike_index <- spike_index[!((spike_index+i) %in% no_spend_index)]
  }
  
  #compbine no spend and spike index
  break_index <- sort(unique(c(no_spend_index, spike_index)))
  
  return(list(break_index = break_index,
              spike_index = spike_index,
              no_spend_index = no_spend_index))
}


# This function finds the ADBUDG saturation curve
# sigma(x) = L * x^s / (x^s + K)
# Inputs:
# x: the adstock sequence
# b: the adstock parameter
# L: initial asymptote (param[1])
# K: initial value of K (param[2])
# s: initial value of s (param[3])
# maxeval: maximum number of iterations in the optimisation of the function.
# Output: a list giving:
# func: function giving the saturation function
# sol: the parameters of the saturation function
sat_curve <- function(x,b,L=NULL,K=NULL,s=NULL,maxeval = 10000,norm = 'euclidean'){
  if(is.null(L)){
    L <- 2*max(x*b)
  }
  if(is.null(K)){
    K <- max(x)
  }
  if(is.null(s)){
    s <- 1.2
  }
  
  if(norm == 'euclidean'){
    loss <- function(param){
      return(sum((param[1]*x^param[3]/(x^param[3] + param[2]^param[3]) - x*b)^2))
    }
  }else{
    if(norm == 'infinity'){
      loss <- function(param){
        return(max(abs(param[1]*x^param[3]/(x^param[3] + param[2]^param[3]) - x*b)))
      }
    }else{
      loss <- function(param){
        return(sum(abs(param[1]*x^param[3]/(x^param[3] + param[2]^param[3]) - x*b)))
      }
    }
  }
  
  test <- nloptr::nloptr(c(L,K,s),loss,lb = c(0,0,0), 
                         opts = list(algorithm = 'NLOPT_LN_COBYLA',
                                     maxeval = maxeval))
  
  sat <- function(a){
    return(test$solution[1]*a^test$solution[3]/(a^test$solution[3] + test$solution[2]^test$solution[3]))
  }
  
  return(list(func = sat, sol = test$solution))
}  
