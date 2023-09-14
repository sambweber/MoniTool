simulate_season = function(N,phenology,theta){
  
  if(is(phenology,'MTpred')){
    if(has_name(phenology,'mu')) p = phenology$mu else {
      stop("'phenology' should be a numeric vector of weights for each day in the time series or an 'MTpred' object")
    }
  } else {p = phenology}
  
  t  = length(p)
  w  = rgamma(t,shape=theta,scale=1/theta)
  mu = N*p
  Y = table(factor(sample(t, N, replace=TRUE,prob=p*w), levels=1:t)) 
  Y = as.numeric(Y)

  if(is.data.frame(proportions)){
    proportions$mu = mu; proportions$Y = Y
    return(proportions)
  } else {return(Y)}
  
  }
