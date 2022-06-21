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

MT_predict = function(model,data,days,samples = 2000, append.observed = T){
  
  samples = MT_sample(model,data,samples)
  
  args = formalArgs(meanFnNim) %>% subset(.!='t')
  arg.pos = match(args,names(samples))

  if(missing(days)) days = range(data$day)
  
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

# ----------------------------------------------------------------------------------------------------------
# MT_meanCI: Returns the mean and credible interval for counts on each day/site for plotting seasonal trends
# ----------------------------------------------------------------------------------------------------------

MT_meanCI = function(preds,by_site=T,interval = 0.95){
  
  group_by(preds,y.var,day,.draw) %>%
  {if(by_site & has_name(predictions,'beach')) group_by(.,beach,.add=T) else .} %>%
  summarise(Y = sum(Y),.groups = 'keep') %>%
  ungroup(.,.draw) %>% mean_qi(Y,.width = interval) %>% 
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
  summarise(Y = sum(Y),.groups = 'keep') %>% 
  {if(!full.posterior){
  ungroup(.,.draw) %>% mean_qi(Y,.width = interval) %>% 
  dplyr::select(-(.width:.interval))
  } else .} %>%
  ungroup()
  
}
