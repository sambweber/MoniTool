# --------------------------------------------------------------------------------------
# acf.negbin
# --------------------------------------------------------------------------------------

# Function for generating autocorrelated count data, where N is the length of the series,
# mu is the mean (or a vector of means of length N) generating counts
# size is the dispersion parameter of negative binomial
# alpha is a numeric vector giving the autocorrelation structure at each lag
# max.iter is the number of iterations attempted to approximate alpha and
# tol is Sum Squared Errors when a solution is accepted to be close enough to the desired alpha. 

acf.negbin <- function(N, mu, size, alpha, max.iter = 100, tol = 1e-5) {
  
  m = length(alpha)
  mu1 = max(mu)
  
  iter <- 0L
  
  # This next line line generates a bunch of random numbers and puts them in 
  # AR1 order, finds their relative size and then puts x in the same order
  generate = function(){
    x <- sort(rnbinom(N,size=size,mu=mu1))
    y <- rnorm(length(x))
    x <- x[rank(stats::filter(y, alpha, circular = TRUE))]
    x * (mu/mu1)
  }
  
  a = generate()
  
  # This extracts the ACF coefficients
  ACF <- function(x) acf(x, lag.max = m - 1, plot = FALSE)$acf[1:m]
  SSE <- function(x,alpha) sum((ACF(log(x) - log(mu)) - alpha)^2)
  
  while ((SSE0 <- SSE(a, alpha)) > tol) {
    if ((iter <- iter + 1L) == max.iter) break
    a1 <- generate()
    if(SSE(a1,alpha) < SSE0) a <- a1
  }
  
  # best to round at end so that no zeros in generated counts for SSE function
  # above. Or could round above only if no zeros and then round again at end.
  return(round(a))
  
}


# --------------------------------------------------------------------------------------
# simulate_trend
# --------------------------------------------------------------------------------------

# Functions to model trends in mean where t is time and init.size is number
# at t = 0. Other parameters control shape

exp.decline = function(init.size,s,t){
  init.size * exp(-((t) / s)^2)
}

lin.decline = function(init.size,beta,t){
  exp(log(init.size) + beta*t) 
}

# Here length is length of the series, model is the trend function above to use, theta is 
# negative binomial dispersion parameter, alpha is an optional autocorrelation structure passed
# to acf.negbin and parms is a list of additional parameters passed to model.

simulate_trend = function(init.size, length, model = c('exp.decline'), theta, alpha, parms,
                          ...){
  
  parms$init.size = init.size
  parms$t = 0:length
  model = match.fun(model)
  mu = do.call(model,args = parms)
  
  if(missing(alpha)){
    y = rnbinom(length,size=theta,mu=mu)
  } else{
    y = acf.negbin(length,mu=mu,size=theta,alpha=alpha,...)
  }
  
  return(data.frame(t=parms$t,mu=mu,y=y))
  
}



# --------------------------------------------------------------------------------------
# simulate_season
# --------------------------------------------------------------------------------------

# Simulate counts across nesting season. Where phenology is either an MTfit model object fit by monitool, or a numeric vector of length days giving the relative amount of activity on each day.
# Days is a numeric vector giving the length of the simulated series, N is the total number of counts that occurred across 'days' and theta is the dispersion parameter of the negative binomial generating counts.

simulate_season = function(phenology,days,N,theta){
  
  if(is(phenology,'MTpred')){
      predict(phenology,samples=1,days = days)$mu
      stop("'phenology' should be a numeric vector of weights for each day in 'days' or an 'MTpred' object")
    }
  } else {p = phenology}
  
  t  = length(days)
  w  = rgamma(t,shape=theta,scale=1/theta)
  mu = N*p
  Y = table(factor(sample(t, N, replace=TRUE,prob=p*w), levels=1:t)) 
  Y = as.numeric(Y)

  if(is.data.frame(proportions)){
    proportions$mu = mu; proportions$Y = Y
    return(proportions)
  } else {return(Y)}
  
  }