# ------------------------------------------------------------------------------------------------------
#  Sampling design functions
# ------------------------------------------------------------------------------------------------------

# These functions assume that x is a sequence of days (or 1:x) if x is a single integer
# We could add optionally here the ability to specify a mix of start and end.

bolus = function(x,start,duration){
  if(length(x)==1) x = 1:x
  x %in% unlist(map(start,~x[.x:(.x+duration-1)]))
}

# A special case of bolus to straddle peak date
peak = function(x,peak,duration){
  bolus(x,start=floor(peak[1]-(duration-1)/2),duration) 
}

staccato = function(x,frequency,duration){
  if(length(x)==1) x = 1:x
  x = x+1-min(x)
  x%%frequency %in% 1:duration 
}

# An improved version of bolus that uses the actual numbers in x, not just their position - but this breaks the
# MT_protocol function currently
#bolus = function(x,start,duration){
#x %in% unlist(map2(start,duration,~.x:(.x+.y)))
#}

# ------------------------------------------------------------------------------------------------------
#  MT_protocol: a function for building monitoring protocols that can be applied to data
# ------------------------------------------------------------------------------------------------------

# A function to combine them elegantly - this can either return a subsetting function 
# or a logical vector indicating whether each point should be retained or not. 
# There are two ways of incorporating staccato sampling. If staccato.in.bolus = TRUE and
# there are blocks included using either bolus or peak then the staccato design will be
# applied within blocks. Overlapping blocks are treated as single contiguous blocks for staccato sampling
# to preserve the periodicity.
# If staccato.in.bolus = FALSE, staccato will be applied to the full set and any blocks included will be sampled daily.

# We could add an option here when staccato.in.bolus = TRUE and peak is specified to centre the staccato on the peak.

MT_protocol = function(x,...) UseMethod('MT_protocol')

MT_protocol.numeric = function(x, bolus, staccato, peak, 
                         staccato.in.block = TRUE, .fun = TRUE){
  
  y = rep(TRUE,length(x))
    
  if(!missing(bolus)){
    bolus$x = x
    y = y & do.call('bolus',args = bolus)
  }
  
  if(!missing(peak)){
    peak$x = x
    y1 = do.call('peak',args = peak)
      if(missing(bolus)) y = y&y1 else y = y|y1
  }
  
  if(!missing(staccato)){
    
    if(staccato.in.block & !(missing(bolus) & missing(peak))){
      
      start = c(which(y[1]),which(diff(y)==1)+1)
      end   = c(which(diff(y)==-1),which(y[length(y)])) 
      
      # The indexing in this section requires fixing for the new bolus function above
      for(i in 1:length(start)){
        z = c(list(x=end[i]-start[i]+1),staccato)
        idx = start[i]:end[i]
        y[idx] <- y[idx] & do.call('staccato',z)
      }
    } else {
  
    staccato$x = x
    y1 = do.call('staccato',args = staccato)
    if(missing(bolus) & missing(peak)) y = y&y1 else y = y|y1
    }
  }
  
  if(.fun){
   f = function(x) x[y,]
   class(f) <- c('protocol.fun',class(f))
   return(f)
  } else {
   return(y)
  }
  
}

# Behaviour for MTsim is if the function doesn't return a subsetting function, add a column called 'include' indicating
# whether each observation is retained.

MT_protocol.MTsim = function(MTsim,...){ 
  x = MT_protocol(MTsim$day, ...) 
  if(!inherits(x,'protocol.fun')) {
    MTsim$.include = x
    return(MTsim) } else {
  return(x)
 }   
} 


