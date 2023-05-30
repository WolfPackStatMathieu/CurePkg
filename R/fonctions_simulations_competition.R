library(survival)
library(cmprsk)


######################################################
# simulations exponential hazard rate exponential
######################################################
# pi, between 0 and 1, =  vector of marginal cumulative incidences at the end of the observation window by dose levels
# obswin, strictly positive= the duration of the observation window 
# type= the distribution for the event process: options are "exponential" or "weibull". If weibull, see below argument typ_wb for options. 
# typ_wb: type of Weibull distribution in terms of shape: "constant", "decreasing" or "increasing" hazard over time. Default shape parameter values hav been specified below and could be modified if desired. 

# output = if "exponential", the vector of instantaneous hazards for each dose level
#          if 'weibull', the table of shape (column 1) and scale (column 2) parameter of the Weibull distribution for each dose level (row)

get_alpha <- function(pi, obswin, typ="exponential", typ_wb=NULL){
  
  if(typ %in% c("exponential", "clayton")){
    alpha1 <- (- log(1-pi)/obswin)
  } 
  
  
  if(typ=="weibull"){
    
    if(typ_wb=="decreasing") {a=0.3}
    if(typ_wb=="increasing") {a=3}
    if(typ_wb=="constant") {a=1}
    b = exp(1/a * log(-obswin^a/log(1-pi)))
    alpha1 <- cbind(rep(a, length(pi)), b)
  }
  
  
  return(alpha1)
}


######################################################
# simulations time exponential
######################################################
# n= the desired number of simulated patients
# K= the number of dose levels
# lambdaT= the vector of instantaneous hazards at each dose level for event T
# lambdaP= the vector of instantaneous hazards at each dose level for event P
# u1= vector of profile values for the time to any event
# u2= vector of profile values for the type of event (T or P)

get_expo <- function(n, K, lambdaT, lambdaP, u1, u2){
  
  tt <- qexp(u1,  lambdaT + lambdaP)
  dd <- qbinom(u2, 1, lambdaP / (lambdaT + lambdaP))+1
  
  return(cbind(u1=u1, u2=u2, tt=tt, dd=dd))
}

######################################################
# simulations time Weibull
######################################################
# n= the desired number of simulated patients
# K= the number of dose levels
# lambdaT= the vector of instantaneous hazards at each dose level for event T
# lambdaP= the vector of instantaneous hazards at each dose level for event P
# u1= vector of profile values for the time to any event
# u2= vector of profile values for the type of event (T or P)


temps <- function (x, a1, b1, a2, b2, u){
  
  udiff <- (x/b1)^a1 + (x/b2)^a2 + log(u)
  
  return(udiff)
}


hx <- function(a,b,x) { # compute the Weibull instantaneous hazard, given shape (a) and scale (b) parameters at a given time (x)
  
  hx1 = a/b * (x/b)^(a-1)  
  return(hx1)
  
}

one_weibull_comp <- function(a1, b1, a2, b2, u1, u2){ # obtain one simulated patient
# a1, the vector (length n) of shape parameters for event T for each patient 
# b1, the vector (length n) of scale parameters for event T for each patient
# a2, the vector (length n) of shape parameters for event P for each patient 
# b2, the vector (length n) of scale parameters for event P for each patient
# u1, the vector of profile for the time to any event (length n) 
# u2, the vector of profile for the type of event (length n) 
  
  u_efs <- u1  
  t_efs <- uniroot(temps, a1=a1, b1=b1, a2=a2, b2=b2, u=u_efs, interval= c(1.e-14, 1e04),
                   extendInt="yes")$root
  
  hx1 <- hx(a1, b1, t_efs)
  hx2 <- hx(a2, b2, t_efs)
  
  evt <- qbinom(u2, 1, hx2/(hx1+hx2))+1
  
  return(c(u1=u1, u2=u2, tt=t_efs, dd=evt))
}



get_weibull <- function(n, a1, b1, a2, b2, u1, u2){ # obtain multiple simulated patients
# n the number of desired patients
# a1, the vector (length n) of shape parameters for event T for each patient 
# b1, the vector (length n) of scale parameters for event T for each patient
# a2, the vector (length n) of shape parameters for event P for each patient 
# b2, the vector (length n) of scale parameters for event P for each patient
# u1, the vector of profile for the time to any event (length n) 
# u2, the vector of profile for the type of event (length n) 
  
  
  
  if(n !=length(a1)){
    print("unequal length")
  }else {
    
    baz <- cbind( a1, b1, a2, b2, u1, u2)
    res <- t(apply(baz, 1, function(x){one_weibull_comp(x[1], x[2], x[3], x[4], x[5], x[6])}))
    colnames(res) <- c("u1", "u2", "tt", "dd")
  } 
  return(res)
}


######################################################
# simulations exponential EFS
######################################################
get_dataset0 <- function(n=60, alpha1, alpha2,  tstar, K=5, graine=1234, type='exponential') {
# obtain a full simulated dataset of multiple patients: one line correspond to the simulated outcome under a given dose level. All potential outcomes at all dose level are generated for all patient (number of lines=K*n)
  
# n= number of patients
# alpha1= vector of marginal instantaneous hazards for event T for each dose level (length=K)
# alpha2= vector of marginal instantaneous hazards for event P for each dose level (length=K)
# tstar= duration of the observation window (time of administrative censoring)
# K= number of dose levels
# graine = seed for computations
# type=type of marginal event distribution (applies to both events)
  
  
  
  set.seed(graine)

  if(type=='exponential'){
    
    alpha1 <- rep(alpha1, each=n)
    alpha2 <- rep(alpha2, each=n)
    
    u1 <- rep(runif(n), K) # individual profile for time to any event for each patient (n distinct values), each repeated for each dose level. total length=n*K
    u2 <- rep(runif(n), K) # individual profile for type of event for each patient (n distinct values), each repeated for each dose level. total length=n*K
    
    res <- as.data.frame(get_expo(n, K, alpha1, alpha2, u1, u2))
  }
  

  if(type=="weibull"){
    u1 <- rep(runif(n), K)
    u2 <- rep(runif(n), K)
    
    res <- as.data.frame(get_weibull(n*K, a1=rep(alpha1[,1], each=n), b1=rep(alpha1[,2], each=n), a2=rep(alpha2[,1], each=n), b2=rep(alpha2[,2],each=n), u1, u2))
    
  }
  
  
  dd <- ifelse(res$tt > tstar, 0, res$dd) #administrative censoring
  tt <- ifelse(res$tt > tstar, tstar+1, res$tt) #administrative censoring
  
  baz <- data.frame(id=1:(n*K), u1=res$u1, u2=res$u2, tt=ceiling(tt*7)/7, dd=dd, status=dd, tstar=tstar)
  baz$id <- rep(1:n, K)
  baz$dose <- rep(1:K, each=n)
  baz <- baz[order(baz$id), c('id','u1' ,'u2', 'tt',"status",'dose')] # ordered by patient
  colnames(baz) <- c('patient','u1' ,'u2', 'time',"status",'dose')
  baz <- as.matrix(baz)
  
  T_entrance_basic = c(0,rep(tstar,n-1))# time of inclusion (on the trial time scale), uniform, 1 patient per window
  T_entrance_basic = cumsum(T_entrance_basic)
  
  T_entrance_fixrapid = c(0, rep(tstar/4 ,n-1))# time of inclusion (on the trial time scale), uniform, 4 patients per window
  T_entrance_fixrapid = cumsum(T_entrance_fixrapid)
  
  T_entrance_fixslow= c(0, rep(tstar/2 ,n-1))# time of inclusion (on the trial time scale), uniform, 2 patients per window
  T_entrance_fixslow = cumsum(T_entrance_fixslow)
  
  
  T_entrance <- data.frame(T_entrance_basic=T_entrance_basic,
                           T_entrance_fixslow=T_entrance_fixslow,
                           T_entrance_fixrapid=T_entrance_fixrapid
                           
  )
  
  return(list(data_complete=baz, T_entrance=T_entrance) )
  
}

