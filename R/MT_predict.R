require(tidybayes)

# ----------------------------------------------------------------------------------
# Function for sampling from the final model
# ----------------------------------------------------------------------------------

# Where this fails it has been because the default sep argument in 
# spread draws "[, ]" doesn't work. Dissection of the source code
# has shown this to be the case and that replacing with sep = ', ' fixes it.

MT_sample = function(model,data,n = 2000){
model %<>% recover_types(data)
samples = spread_draws(model,alpha[Y,beach],s1[Y,beach],s2[Y,beach]
                       ,tf[Y,beach],tp[Y,beach],phi[Y,beach],ndraws=2000,sep=', ') %>%
          mutate(.draw = 1:n) %>% # renumber draws sequentially from one to facilitate integration of sites later
          ungroup()
samples$Y <- attr(model,'y.names')[samples$Y]

if(!has_name(data,'beach')) {samples$beach <- NULL}

rename(samples,y.var=Y)
}

# ----------------------------------------------------------------------------------
# Function for predicting mean and values from fitted model
# ----------------------------------------------------------------------------------

# This currently matches column positions to argument names in meanFnNimm, but it won't work 
# if meanFnNim has a different number of arguments e.g. a different functional form as 
# we expand the code later

predict.MTfit = function(model,data,days,samples = 2000){
  
  samples = MT_sample(model,data,samples)
  
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
  
  return(samples)
}
                                           
# ----------------------------------------------------------------------------------
# MT_predict: Make predictions from an MT_tbl object
# ----------------------------------------------------------------------------------

MT_predict <- function(object,days,samples = 2000){
  
  if(!all(c('data','fit') %in% names(object))) stop ("object should be of class MT_tbl with columns 'data' and 'fit'")
  
  if(missing(days)){ 
    days = bind_rows(object$data)$day
    days = min(days):max(days)
   }
  
  mutate(object,predict = map2(fit,data,predict.MTfit,days=days,samples=samples))
  
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
  
season_merge = function(obj){
  
  dplyr::select(obj,any_of(c('season','reference_date','beach','predict'))) %>%
  unnest(predict) %>% 
  nest(predict = -(any_of(c('season','reference_date'))))
   
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
# MT_totalCI: Calculates the total number of activities per site or overall along with associated credible intervals
# ----------------------------------------------------------------------------------------------------------
 
# We could optionally here allow use of the totals including observations (Y) or only predicted values (y.pred) to
# assess changes in confidence intervals.

MT_totalCI = function(predictions,by_site=T,interval=0.95,full.posterior=F){
  
  group_by(predictions,y.var,.draw) %>%
  {if(by_site & has_name(predictions,'beach')) group_by(.,beach,.add=T) else .} %>%
  summarise(N = sum(N),.groups = 'keep') %>% 
  {if(!full.posterior){
  ungroup(.,.draw) %>% mean_qi(N,.width = interval) %>% 
  dplyr::select(-(.width:.interval))
  } else .} %>%
  ungroup()
  
}

  
# ----------------------------------------------------------------------------------------------------------
# MT_totalCI: Calculates proportions on beaches within seasons
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
