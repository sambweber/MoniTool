require(tidybayes)

# ----------------------------------------------------------------------------------
# Function for sampling from the final model
# ----------------------------------------------------------------------------------

# Where this fails it has been because the default sep argument in 
# spread draws "[, ]" doesn't work. Dissection of the source code
# has shown this to be the case and that replacing with sep = ', ' fixes it.

# Use pars argument with caution - not yet fully implemented!!!!!

MT_sample = function(model,n = 2000,pars){

data = attr(model,'data')
model %<>% recover_types(data)
samples = spread_draws(model,alpha[Y,beach],s1[Y,beach],s2[Y,beach]
                       ,tf[Y,beach],tp[Y,beach],phi[Y,beach],ndraws=n,sep=', ') %>%
          mutate(.draw = 1:n) %>% # renumber draws sequentially from one to facilitate integration of sites later
          ungroup()
samples$Y <- attr(model,'y.names')[samples$Y]

if(!has_name(data,'beach')) {samples$beach <- NULL}
  
samples = rename(samples,y.var=Y)

# Subsetting of specific parameters could be done much more efficiently in spread draws above 
# but this works as quick fix. It doesn't currently respect different beaches or y.vars
if(!missing(pars)) samples = samples[,pars,drop=TRUE]

return(samples)
  
}

# ----------------------------------------------------------------------------------
# Function for predicting mean and values from fitted model
# ----------------------------------------------------------------------------------

# This currently matches column positions to argument names in meanFnNimm, but it won't work 
# if meanFnNim has a different number of arguments e.g. a different functional form as 
# we expand the code later.... we could this do this by collapsing the args to a list in a do.call maybe
# This is the better way of doing this using column-argument matching:
# MT_sample(model,samples) %>% 
# mutate(.iter = row_number()) %>% 
# expand_grid(t = days) %>%
# mutate(total = pmap_dbl(list(!!!rlang::parse_exprs(formalArgs('meanFnNim'))), meanFnNim))
# Will also enable different mean functions (models) to be passed. 

predict.MTfit = function(model,days,samples = 2000){

  data = attr(model,'data')
  
  samples = MT_sample(model,samples)
  
  args = formalArgs(meanFnNim) %>% subset(.!='t')
  arg.pos = match(args,names(samples))

  if(missing(days)) days = min(data$day):max(data$day)
  
  samples %<>%
  mutate(preds = pmap(.[arg.pos], ~map_dbl(days,function(t) meanFnNim(t,..1,..2,..3,..4,..5)))) %>%
  mutate(preds = map(preds, ~data.frame(day=days,mu = .x))) %>%
  dplyr::select(.draw,y.var,any_of('beach'),phi,preds) %>% unnest(cols=preds) %>%
  mutate(N.pred = rnbinom(day,phi,mu=mu)) %>%
  dplyr::select(-phi)
  
  samples %<>% MT_append_observed(data)
  
  class(samples) <- c('MTpred',class(samples))                                         
  return(samples)
}
                                           
# ----------------------------------------------------------------------------------
# MT_predict: Make predictions from an MT_tbl object
# ----------------------------------------------------------------------------------

MT_predict <- function(object,days,samples = 2000,ncores){
  
  if(!all(c('data','fit') %in% names(object))) stop ("object should be of class MT_tbl with columns 'data' and 'fit'")
  
  if(missing(days)){ 
    days = bind_rows(object$data)$day
    days = min(days):max(days)
   }


  if(ncores>1){
  
    plan(tweak(multisession,workers = ncores))
    object = mutate(object,predict = furrr::future_map(fit,predict.MTfit,days=days,samples=samples))
    plan(sequential)
  
} else {
  
    object = mutate(object,predict = map(fit,predict.MTfit,days=days,samples=samples))
  
  }

    return(object)
  
}
  
                                             
# ----------------------------------------------------------------------------------
# MT_append_observed: Replaces predictions with observed counts, where available
# ----------------------------------------------------------------------------------

MT_append_observed = function(preds,data){

  t = MT_diff_matrix(data$day,data$window)

  obs.days = dplyr::select(data,any_of(c('beach','day',unique(preds$y.var)))) %>%
             gather('y.var','N.obs',-any_of(c('beach','day')))
  
  obs.days = 
    data.frame(t) %>% 
    {if(has_name(data,'beach')) bind_cols(.,beach=data$beach) else .} %>% 
    gather('.col','day',-any_of('beach')) %>%
    dplyr::select(-.col) %>% 
    drop_na %>%
    left_join(obs.days) %>%
    replace_na(list(N.obs = 0)) 
  
  preds %>% left_join(obs.days,by=intersect(c('day','beach','y.var'),names(obs.days))) %>%
  mutate(N = coalesce(N.obs,N.pred))

}
  
# ---------------------------------------------------------------------------------------------------------------------------------------
# plot.MTpred: plot method for MTpred
# --------------------------------------------------------------------------------------------------------------------------------------- 
  
plot.MTpred <- function(obj,by_site=T){
 
orig = na.omit(distinct(dplyr::select(obj,y.var,any_of('beach'),day,N.obs)))
cat('Calculating mean and CI')
mean.line = MT_meanCI(obj,by_site=by_site)

pl = 
  ggplot(mean.line,aes(x = day,y = mu)) + 
  geom_ribbon(aes(ymin = .lower,ymax = .upper,group = y.var),alpha=.3) +
  geom_line(aes(colour=y.var)) + 
  scale_colour_manual(values = c(nests='orange',activities='blue')) + 
  ylab('Count')

if(by_site & length(orig)) pl = pl + geom_point(data=orig,aes(y = N.obs,colour=y.var),shape=21)
  
if(has_name(mean.line,'beach')) pl + facet_wrap(~beach) else pl

}
    
  
# ---------------------------------------------------------------------------------------------------------------------------------------
# summary.MTfit: summary method for MTfit 
# --------------------------------------------------------------------------------------------------------------------------------------- 

#' @param model A model fit using MT_fit  
  
#' @returns A tibble containing the mean, median, sd, credible intervals, potential scale reduction factor and effective sample size for
#' parameters.
  
  summary.MTfit = function(model){
    
  .s = suppressWarnings(spread_draws(model,alpha[Y,beach],
               phi[Y,beach],
               s1[Y,beach],
               s2[Y,beach],
               tf[Y,beach],
               tp[Y,beach])) %>%
  summarise_draws() %>%
  dplyr::select(-mad,-ess_tail) %>%
  rename(n_eff = ess_bulk)
  
  if(all(is.na(.s$beach))) .s$beach <- NULL
  
  .s$Y <- attr(model,'y.names')[.s$Y]
  rename(.s,y.var=Y)
  
}
  
# ----------------------------------------------------------------------------------------------------------
# season_merge: Helper function to combine predictions for all sites within seasons for summarizing
# ----------------------------------------------------------------------------------------------------------
  
merge_seasons = function(obj){
  
  if(!all(has_name(obj,c('season','beach')))) stop("obj should be of class MT_tbl with a 'season' and 'beach' column")
  
  condense = function(.x) list(setNames(.x,obj$beach) %>% bind_rows(.id = 'beach'))
  group_by(obj,season) %>% 
  summarise(across(any_of(c('data','summary','predict')),condense),
            across(any_of('fit'),list)) %>% 
  dplyr::select(any_of(names(data)))
   
}

# ----------------------------------------------------------------------------------------------------------
# MT_meanCI: Returns the mean and credible interval for counts on each day/site for plotting seasonal trends
# ----------------------------------------------------------------------------------------------------------

MT_meanCI = function(preds,by_site=T,interval = 0.95,what = c('counts','proportions')){
  
  what = match.arg(what)
  
  group_by(preds,y.var,day,.draw) %>%
  {if(by_site & has_name(preds,'beach')) group_by(.,beach,.add=T) else .} %>%
  summarise(mu = sum(mu),.groups = 'keep') %>%
  {if(what=='proportions') {ungroup(.,day) %>% mutate(mu = mu/sum(mu)) %>% group_by(day,.add=T)} else .} %>%
  ungroup(.draw) %>% mean_qi(mu,.width = interval) %>% 
  dplyr::select(-(.width:.interval)) %>%
  ungroup()
  
}

# ----------------------------------------------------------------------------------------------------------
# MT_seasonCI: Calculates the start, end and duration of a season based on a given percentage of nesting
# ----------------------------------------------------------------------------------------------------------

MT_season = function(predictions,by_site=T,quantile = c(0.025,0.975),interval=0.95,full.posterior=F){
  
  group_by(predictions,y.var,.draw) %>%
  {if(by_site & has_name(predictions,'beach')) group_by(.,beach,.add=T) else .} %>%
  mutate(s = cumsum(mu)/sum(mu)) %>% 
  summarise(start = first(day[s>=quantile[1]]),end = first(day[s>=quantile[2]]),duration=end-start) %>% 
  {if(!full.posterior){
  ungroup(.,.draw) %>% mean_qi(start,end,duration,.width = interval) %>% 
  dplyr::select(-(.width:.interval))
  } else .} %>%
  ungroup()
  
}

# This version works directly from fitted model and adds the peak too - not currently from the same draws as the
# start and end though which requires improving
  
MT_seasonality = function(model,samples,quantile = c(0.025,0.975),interval=0.95,full.posterior=F){
  
  predict(model,days=1:365,samples=samples) %>% 
  MT_season(quantile=quantile,full.posterior=TRUE) %>%
  left_join(MT_sample(model,n=samples)) %>%
  dplyr::select(y.var:duration,peak=tp) %>%
  {if(!full.posterior){
  group_by(.,across(any_of(c('y.var','beach')))) %>% 
  mean_qi(start,end,duration,peak,.width = interval) %>% 
  dplyr::select(-(.width:.interval))
  } else .} %>%
  ungroup()
  
  }

# ----------------------------------------------------------------------------------------------------------
# MT_totalCI: Calculates the total number of activities per site or overall along with associated credible intervals
# ----------------------------------------------------------------------------------------------------------
 
# We could optionally here allow use of the totals including observations (Y) or only predicted values (y.pred) to
# assess changes in confidence intervals.

MT_totalCI = function(predictions,by_site=T,interval=0.95,full.posterior=F){
  
  predictions = 
  group_by(predictions,y.var,.draw) %>%
  {if(by_site & has_name(predictions,'beach')) group_by(.,beach,.add=T) else .} %>%
  summarise(N = sum(N),.groups = 'keep') 
   
  if(full.posterior) return(ungroup(predictions))
  
   predictions = ungroup(predictions,.draw)
   mean_qi(predictions,N,.width = interval) %>%
   inner_join(summarise(predictions,log.sd = sd(log(N)))) %>%  
   dplyr::select(-(.width:.interval)) %>%
   ungroup()
 
}

  
# ----------------------------------------------------------------------------------------------------------
# MT_propCI: Calculates proportions on beaches within seasons
# ----------------------------------------------------------------------------------------------------------
 
# We could optionally here allow use of the totals including observations (Y) or only predicted values (y.pred) to
# assess changes in confidence intervals.

MT_propCI = function(preds,interval=0.95,full.posterior=F){
  
  if(!has_name(preds,'beach')) stop ("can only calculate proportions where column 'beach' provided")
  
  group_by(preds,y.var,.draw,beach) %>%
  summarise(N = sum(N),.groups = 'keep') %>% 
  ungroup(beach) %>% mutate(proportion = N/sum(N)) %>% group_by(beach,.add=T) %>%  
  {if(!full.posterior){
  ungroup(.,.draw) %>% mean_qi(proportion,.width = interval) %>% 
  dplyr::select(-(.width:.interval))
  } else .} %>%
  ungroup()
  
}
