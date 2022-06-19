
# -------------------------------------------------------------------------------------------------------------------------
# MT_prep 
# -------------------------------------------------------------------------------------------------------------------------

#' Puts data in right format for a MoniTool model run
#
#' @param data is a dataframe which must have a column called 'dateend' which was the date on which 
#' the count was done, along with at least one of columns entitled'nests' or 'activities' (can be both). 
#' May also optionally have columns called 'beach' to differentiate sites, and 'season' to name nesting 
#' seasons. If 'season' is missing, this will be calculated automatically from the dates. 


MT_prep = function(data, season_start, max.window = 1, min.obs = 5, sites_seperate = FALSE){
  
  if(!has_name(data,'dateend')) stop('Data must contain a column called dateend')
  
  data %<>% mutate(reference_date = set_ref_date(dateend,season_start), 
                   day = as.numeric(dateend - reference_date))
  
  if(!has_name(data,'season')) {
    
  data %<>% group_by(reference_date) %>%
  mutate(season = max(year(dateend))) %>% rowwise()  %>%
  mutate(season = factor(paste(unique(c(year(reference_date),season)),collapse='-')))    
      
  }
  
  if(!has_name(data,'datestart')) {
  
  data %<>% group_by(across(any_of(c('season','Beach')))) %>%
  mutate(window = map_dbl(diff(c(-1,day)),~min(.x,max.window)), datestart = dateend-window)
      
  } else data %<>% mutate(window = dateend-datestart)
  
  # Remove beaches that have too few counts to be effectively modelled
  data %<>% group_by(across(any_of(c('season','Beach')))) %>%
  mutate(too_few = n()<=min.obs) %>% ungroup()
  
  if(sum(data$too_few)>1){
    cat('The following seasons/beaches have insufficient data and have been removed:\n\n')
    print(data.frame(distinct(subset(data,too_few),season,Beach)),row.names = F)
  }
  
  # Prep final output
  subset(data,!too_few) %>%
  dplyr::select(any_of(c('season','Beach','day','datestart','dateend','window')),everything(),-too_few) %>%
  nest(data = -c(season,reference_date)) %>%
  mutate(data = map(data, droplevels))

}
