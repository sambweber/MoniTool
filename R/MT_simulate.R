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
