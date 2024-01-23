# ---------------------------------------------------------------------------------------
# phenology_ML
# ---------------------------------------------------------------------------------------

# This function fits the Omeyer et al. MeanFnNim model via maximum likelihood. It needs 
# modifying to take different window lengths for grouped counts. Hessian also needs adding?
# For 95% intervals try this.....
# https://stackoverflow.com/questions/41541528/optimization-of-optim-in-r-l-bfgs-b-needs-finite-values-of-fn
# https://stats.stackexchange.com/questions/27033/in-r-given-an-output-from-optim-with-a-hessian-matrix-how-to-calculate-paramet

# This could be used to quicky find good initial values for the Bayesian approach, or for fast fitting. 

# Arguments: t is a vector of days (must be [0,366]) and N is vector of counts.

phenology_ML = function(t,N){
  
  log.lik = function(par,t,N){
    mu = map_dbl(t,~meanFnNim(.x, s1=par[1],s2=par[2], tp=par[3], tf=par[4], alpha=par[5]))
    -sum(dnbinom(N,size=par[6],mu=mu,log=TRUE))
  }
  
  init.alpha = max(N)
  init.pk = t[which.max(N)]
  s.grid = expand_grid(s1=(1:10)^2,s2=(1:10)^2)
  
  error = TRUE; i = 0
  while(error & i<= nrow(s.grid)){
    i = i+1
    tryCatch({
      a <- optim(par = c(s.grid[i,1],s.grid[i,2],tp=init.pk,tf=0,alpha=init.alpha,theta=50),
        fn=log.lik, N=N,t=t,lower=c(rep(0,6)),
        upper = c(Inf,Inf,365,180,Inf,50),
        method='L-BFGS-B')
  error = FALSE}, error = function(e){}
  )
  }
  
  if(error) {
    print('model failed to initialize'); return(NULL)
    } else {return(a)}
}

# example:
# test=
# tibble(t = 58:207) %>% 
# mutate(yhat = map_dbl(x,~meanFnNim(.x, s1=15,s2=25, tp=160, tf=10, alpha=100))) %>%
# mutate(Y = simulate_phenology(yhat, total=sum(yhat),theta=20))
# ggplot(test,aes(x=t,y=yhat)) + geom_line() + geom_point(aes(y=Y))
# phenology_ML(test$t,test$Y)
# It recovers the original parameters quite well with this saturation dataset...

# 95% CI:
# sigma = sqrt(diag(solve(fit$hessian)))
# upper<-fit$par+1.96*sigma
# lower<-fit$par-1.96*sigma
