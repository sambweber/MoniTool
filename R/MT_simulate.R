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

# Question over what phenology curve to simulate from if both nests and activities in MTfit object and have different curves. Default at the moment is just to use first (or only) which always be activities if both
# modelled together - but could make so that it returns esimtates for any present and they could be subset later. All that is being borrowed at the moment is the shape of the curve and the theta value.

simulate_season = function(phenology,days,N,theta){
  
  t  = length(days)
            
  if(is(phenology,'MTfit')){
    pred = predict(phenology,samples=1,days = days)
    pred = subset(pred,y.var == attr(phenology,'y.names')[1]) #Only select first y.var if several in same model
    p = pred$mu
    theta = MT_sample(phenology,1)$phi[1]     # Sample theta from the posterior (phi is theta in this case) - only for first y.var if several present at moment
  } else if (is(phenology,'numeric') & length(phenology) == t){
    p = phenology
  } else stop("'phenology' should be a numeric vector of weights the same length as 'days' or an 'MTpred' object")

w  = rgamma(t,shape=theta,scale=1/theta)
mu = N*p
Y = table(factor(sample(t, N, replace=TRUE,prob=p*w), levels=1:t)) 
Y = as.numeric(Y)

if(is(phenology,'MTfit')){
  pred$sim = Y
  class(pred) <- c('MTsim',class(pred))
  return(pred)
} else {return(Y)}

}


# --------------------------------------------------------------------------------------
# simulate_phenology
# --------------------------------------------------------------------------------------

# An improvement with methods for different kinds of objects

simulate_phenology = function(x,...) UseMethod('simulate_phenology')

simulate_phenology.numeric = function(phenology,total,theta){
  t  = length(phenology)
  w  = rgamma(t,shape=theta,scale=1/theta)
  Y = table(factor(sample(t, total, replace=TRUE,prob=phenology*w), levels=1:t)) 
  as.numeric(Y)
}

simulate_phenology.MTfit = function(phenology,days,total,n.sims){

    stopifnot("'total' should be of length = 1 or n.sims" = length(total) %in% c(1,n.sims))
    y = attr(phenology,'y.names')[1]   #Only select first y.var if several in same model
  
    pred = predict(phenology,samples=n.sims,days = days) %>% 
           rename(.sim = .draw) %>%
           subset(y.var == y) 
    pred = split(pred,pred$.sim)
    
    # At the moment theta and mean predictions come from different posterior draws - this could be improved.
    theta = MT_sample(phenology,n.sims) %>%
            subset(y.var == y) %$% phi

    pred%<>%
    #map2(theta,~mutate(.x,mu = mu/sum(mu)*total,
    #                   sim = simulate_phenology(mu,total,.y))) %>%
    pmap(list(pred,theta,total),~mutate(..1,mu = mu/sum(mu)*..3,
                       sim = simulate_phenology(mu,..3,..2))) %>%
    bind_rows() 
    
    class(pred) <- c('MTsim',class(pred))
    return(pred)
}

simulate_phenology.MT_df = function(phenology,days,total,n.sims){
  
  if(!has_name(phenology,'fit')) stop ('Error')
  models = phenology$fit[sample(1:nrow(phenology),n.sims,replace=T)]
  lapply(models,simulate_phenology,days = days, total = total, n.sims=1) %>%
  bind_rows(.id = '.sim')
  
}


# -------------------------------------------------------------------------------------------------------------
# Plotting method for MTsim objects
# --------------------------------------------------------------------------------------------------------------

# Optionally highlights points that have been retained by a monitoring protocol applied using MT_protocol

plot.MTsim = function(x){
  
p = ggplot(x,aes(x = day, y=sim)) + geom_point() + geom_line(aes(y=mu), colour = 'blue') + facet_wrap(~.sim)

if(has_name(x,'.include')){ p = p + geom_point(data = subset(x,.include),colour='yellow',size=2) }

return(p)
}


# -------------------------------------------------------------------------------------------------------------
# MT_sim2f: Conversion of an MTsim simulation to an MT_df so it can be used to fit models and evaluate outcome
# --------------------------------------------------------------------------------------------------------------

MT_sim2df = function(MTsim){

if(has_name(MTsim,'.include')) MTsim <- subset(MTsim,.include)
x = dplyr::select(MTsim,.sim,day,y.var,sim) %>% 
spread('y.var','sim') %>%
mutate(window = 1) %>%
nest(data = -.sim)

class(x) <- c('MT_df',class(x)[!class(x) %in% c('MTsim','MTpred')])

return(x)
}

# -------------------------------------------------------------------------------------------------------------
# MAPE: Mean average percentage error
# --------------------------------------------------------------------------------------------------------------

MAPE = function(expected,observed){
  mean(abs((expected - observed)/expected))*100
}


