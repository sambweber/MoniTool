# -----------------------------------------------------------------------------------
# MT_aggregate_counts
# -----------------------------------------------------------------------------------

# Function for aggregating counts carried out on the same dates for groups of beaches prior to analysis.
# This can be helpful where counts on individual beaches are very low. 
# This will return NA for any date in which counts are missing for any one of the beaches in each group.

#' @param data A dataframe created by `MT_prep`
#' @param groupings A named list specifying the groups of beaches to be pooled. For example,
#' `groupings = list('A' = c('B','C'))` will make a new group called 'A' containing the summed counts from
#  beaches 'B' and 'C'.

MT_aggregate_counts = function(data,groupings){
  
  data = mutate(data,beach = as.character(beach))
  
  for(g in 1:length(groupings)){
    
    missing = groupings[[g]] %>% subset(.,!. %in% data$beach)
    if(length(missing)){ 
      warning(paste("The following beaches are not present in the data:",paste(missing,collapse=','),
                    "\n Aggregation step skipped for group",names(groupings)[g]),"\n")
      next
    }
    
    data %<>% 
      subset(beach %in% groupings[[g]]) %>%
      complete(beach,nesting(day,datestart,date,window)) %>%
      group_by(day,datestart,date,window) %>% 
      summarise(across(c('activities','nests'),sum)) %>%
      mutate(beach = names(groupings)[g]) %>%
      subset(!(is.na(activities) & is.na(nests))) %>%
      bind_rows(subset(data,!beach %in% groupings[[g]]),.)
  }
  
  mutate(data, beach = as.factor(beach))
  
}
