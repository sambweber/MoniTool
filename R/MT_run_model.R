require(nimble)
require(tidyverse)

# -----------------------------------------------------------------------------
# meanFnNim: code to create a nimbleFunction object containing the mean function
# -----------------------------------------------------------------------------

#' Returns a numbleFunction implementing the Omeyer et al. 2022 model of mean number of counts
#' as a function of days. Called internally by `MT_make_model`.

meanFnNim <- nimbleFunction(
  run = function(t = double(0), alpha = double(0), s1 = double(0),
                 s2 = double(0), tf = double(0), tp = double(0)) {
    if(t < (tp - tf)) {
      mu <- alpha * exp(-((t - tp + tf) / s1)^2)
    } else {
      if(t < (tp + tf)){
        mu <- alpha
      } else {
        mu <- alpha * exp(-((t - tp - tf) / s2)^2)
      }
    }
    return(mu)
    returnType(double(0))
  }
)

# -----------------------------------------------------------------------
# MT_NimModel: Basic code for the NIMBLE model
# -----------------------------------------------------------------------

# Generates `nimbleCode` for fitting the Omeyer et al. 2022 model in `Nimble`

MT_NimModel <- nimbleCode({
  ## likelihood
  
  for(r in 1:R){
    for(i in 1:N) {
     for(j in 1:nt[i]) {
        mu[i, j, r] <- meanFnNim(t[i, j], alpha[r,bch[i]], s1[r,bch[i]], s2[r,bch[i]], tf[r,bch[i]], tp[r,bch[i]])
      }
      Y[i,r] ~ dSnbinomNim(phi[r,bch[i]], mu[i, ,r], nt[i])
    }
  }
  ## priors

  for(r in 1:R){
  
    for(b in 1:B){
      alpha_mn[r,b] ~ dexp(1)
      alpha_rate[r,b] <- 1 / alpha_mn[r,b]
      alpha[r,b] ~ dexp(alpha_rate[r,b])
      
      tf[r,b] ~ dexp(0.2)
      phi[r,b] ~ T(dinvgamma(shape = 1, rate = 0.1), , 50)
      s_rate[r,b] ~ dunif(0.01, 10)
      s1[r,b] ~ dexp(s_rate[r,b])
      s2[r,b] ~ dexp(s_rate[r,b])
      
      #stp[r,b] ~ dunif(0, 10)
      #tp_ast[r,b] ~ dnorm(0, sd = 1)
      #tp[r,b] <- Pk + tp_ast[r,b] * stp[r,b]
      tp[r,b] ~ dunif(140,220)
    }
    
  }  
  
})


# ----------------------------------------------------------------------------------
# MT_diff_matrix: Helper function for making a matrix containing the days observed in each sampling window
# ----------------------------------------------------------------------------------

#' @param t is a vector of survey days from an object created by `MT_prep` 
#' @param w is a vector of sampling windows giving the number of days represented by each survey ub `t`

# ! This may need modifying in the case where w = 1 as single column matrices
# don't seem to be well handled by the custom distribution in NIMBLE

MT_diff_matrix = function(t,w){
  
  tdiffs =  map2(t, w, function(t,w) (t - w + 1):t)
  
  t <- matrix(NA, length(tdiffs), max(c(2,w)))     #Here it forces ncol to be minimum of 2 until issues can be resolved 
  for(i in 1:length(tdiffs)) {
    t[i, 1:length(tdiffs[[i]])] <- tdiffs[[i]]
  }
  t
}

  
# ----------------------------------------------------------------------------------
# MT_make_model: Function for building the nimble model
# ----------------------------------------------------------------------------------

#' @param data is a dataframe created by MT_prep

#' @param response is a character vector giving the response

#' @param Pk is a single numeric giving an approximate peak day for counts which is used as a starting 
#' point for model fits. Omeyer et al. 2022 found that using an approximate peak in priors resulted
#' in better covnvergence.

MT_make_model <- function(data, response = c('activities','nests'), Pk){

  tmat = MT_diff_matrix(data$day,data$window)
  consts = list(t=tmat,nt = rowSums(!is.na(tmat)), N = nrow(tmat),Pk = Pk)
  consts$bch <- if(has_name(data,'beach')) as.numeric(data$beach) else rep(1,consts$N)
  consts$B = max(consts$bch)
  
  # Remove any columns that are all NAs prior to fitting
  Y = as.matrix(data[,response])
  Y = Y[, colSums(is.na(Y)) < nrow(Y), drop=F]
  
  consts$R = ncol(Y)
  consts$y.names = colnames(Y)
  
  nimbleModel(
  code = MT_NimModel,
  constants = consts,
  data = list(Y = Y),
  dimensions = list(mu = c(dim(consts$t),ncol(Y))),
) %>% compileNimble()
  
}


# ----------------------------------------------------------------------------------
# MT_initialize: Function for initializing the nimble model
# ----------------------------------------------------------------------------------

# Generates a valid set of initial parameters for the model by sampling from the prior 
# distributions.

#' @param model is a model structure created by `MT_make_model`

#' @param smin,smax Single numerics giving the lower and upper bounds of a uniform distribution for 
#' the s_rate parameter which are used for adjusting the default prior if needed. See Details.

#' @attemmpts Number of attempts to find a valid set of starting parameters before giving up.

#' In this function, if the prior on s_rate is kept at the default vague priors `runif(0.01,10)` then it can take 
#' an impossibly long number of iterations to find a valid set of initial parameters, because s1 and s2
#' that are sufficiently high and close to real values are very unlikely to be sampled. You can make it fit by narrowing the
#' prior to something like `runif(0.01,0.2)`. Until some way can be found of quickly estimating a sensible
#' starting point for s_rate, options are provided for manually setting the limits of the uniform distribution. Some
#' trial and error may be needed to find values that initialize quickly.

MT_initialize <- function(model, smin = 0.01, smax = 10, attempts = 100) {
  
  valid <- 0
  consts = model$origData
  B = consts$B
  R = consts$R
  tdiffs = consts$t
  bch = consts$bch
  Pk = consts$Pk
  attempt = 1

  while(valid == 0 & attempt < attempts) {
    ## sample initial values from the priors
    
    alpha_mn <- matrix(rexp(B*R,1),nrow=R)
    alpha_rate <- 1 / alpha_mn
    alpha <- matrix(rexp(B*R,alpha_rate),nrow=R)
  
    s_rate <- matrix(runif(R*B, smin, smax),nrow=R)
    s2 <- s1 <- matrix(rexp(R*B,s_rate),nrow=R)
    
    #tp_ast <- matrix(rnorm(B*R, 0, 1),nrow=R)
    #stp <- matrix(runif(B*R, 0, 10),nrow=R)
    #tp <- Pk + tp_ast * stp
    tp <- matrix(runif(B*R, 140, 220),nrow=R)
    
    tf <- matrix(rexp(B*R,0.2),nrow=R)
    
    phi <- matrix(rep(51,R*B),nrow=R)
    while(any(phi > 50)){
    phi <-matrix(rexp(B*R,1),nrow=R)
    }

    mu = array(dim = c(dim(tdiffs),R))
    y.init <- model$Y 
    
    for(r in 1:R){
      for(i in 1:nrow(tdiffs)) {
       for(j in 1:consts$nt[i]) {
         mu[i, j, r] <- meanFnNim(tdiffs[i, j], alpha[r,bch[i]], s1[r,bch[i]], s2[r,bch[i]], tf[r,bch[i]], tp[r,bch[i]])
       } 
        y.init[i,r] <- rnbinom(1,phi,mu=mu[i,1,r])
      }
    }
    
    mu[is.na(mu)] <- 0
    y.init[!is.na(model$Y)] <- NA
    
    inits <- list(
      alpha_mn = alpha_mn,
      alpha_rate = alpha_rate,
      alpha = alpha,
      s_rate = s_rate,
      s1 = s1,
      s2 = s2,
      #stp = stp,
      #tp_ast = tp_ast,
      tp = tp,
      tf = tf,
      phi = phi,
      mu = mu,
      Y = y.init
    )
    ## check validity of initial values
    model$setInits(inits)
    valid <- ifelse(!is.finite(model$calculate()), 0, 1)
    attempt = attempt+1
  }
  if(valid == 0) stop(paste('Model failed to initialize after',attempts,'attempts. Try increasing attempts or adjusting smin and smax'))
  return(inits)
}

# ----------------------------------------------------------------------------------------------------------------
# MT_run_model: Function for compiling and running the nimble model
# ----------------------------------------------------------------------------------------------------------------
                 
init.control = function() list(smin = 0.01,smax=10)

# ----------------------------------------------------------------------------------------------------------------
# MT_run_model: Function for compiling and running the nimble model
# ----------------------------------------------------------------------------------------------------------------

#' @param model A model created by MT_make_model

#' @param nchains,niter,nburnin,thin Single numerics giving the number of MCMC chains to run, the total number of iterations 
#' per chain, the number of burnin samples to discard, and the sample thining rate, respectively

#' @param parameters.to.monitor Names of the model parameters to monitor. Defaults are generally ok.

#' @param init.attempts Number of attempts to make at finding a valid set of starting parameters for the model before 
#' gracefully giving up.

MT_run_model = function(data, Pk, nchains=2, niter=40000, nburnin = 10000, thin = 1,
                     parameters.to.monitor = c("alpha", "s1", "s2", "tf", "tp", "phi"),
                     init.attempts = 100,init.control = list(smin = 0.01,smax=10)) {

assign('dSnbinomNim', dSnbinomNim, envir = .GlobalEnv)
assign('rSnbinomNim', rSnbinomNim, envir = .GlobalEnv)  
  
cat('#Compiling model')
model = MT_make_model(data,Pk=Pk)

cat('#Initialising model')
init.control = modifyList(init.control(),init.control)
inits <- map(1:nchains, ~MT_initialize(model,attempts = init.attempts, smin = init.control$smin,smax=init.control$smax))

cat('#Configuring model')
config <- configureMCMC(model, monitors = parameters.to.monitor, thin = thin)
built  <- buildMCMC(config)
cbuilt <- compileNimble(built)
  
## run the model for multiple chains
fit = runMCMC(cbuilt, nchains = nchains, niter = niter, nburnin = nburnin, inits = inits, thin = thin,
        progressBar = TRUE, samplesAsCodaMCMC = TRUE)

class(fit) <- c('MTfit',class(fit))
attr(fit,'y.names') <- model$origData$y.names

return(fit)

}


# ----------------------------------------------------------------------------------------------------------------
# MT_fit: Function for running a nimble model across multiple sites/seasons in a nested dataste produced by MT_prep
# ----------------------------------------------------------------------------------------------------------------

MT_fit <- function(data, Pk, nchains=2, niter=40000, nburnin = 10000, thin = 1,
                     params = c("alpha", "s1", "s2", "tf", "tp", "phi"),
                     init.attempts = 100,init.control = list(smin = 0.01,smax=10),ncores=1){
 
if(!is(data,'MT_df')) stop('data should be a MT_df object created by MT_prep')
  
if(ncores>1){
  
cat(paste('#Fitting',nrow(data),'models in parallel'))
plan(tweak(multisession,workers = ncores))
data = mutate(data,fit = furrr::future_map(data,MT_run_model,Pk = Pk,nchains=nchains,niter=niter, nburnin = nburnin, thin = thin,
                     parameters.to.monitor = params,init.attempts = init.attempts,init.control = init.control)) 
plan(sequential)
  
} else {
  
data = mutate(data,fit = purrr::map(data,MT_run_model,Pk = Pk,nchains=nchains,niter=niter, nburnin = nburnin, thin = thin,
                     parameters.to.monitor = params,init.attempts = init.attempts,init.control = init.control)) 
  
  }
}  
  
