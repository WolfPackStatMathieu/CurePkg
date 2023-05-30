tp.surv <- function(obj, times) {

  x <- summary(obj)
  if(is.null(x$strata))
  {
    y <- as.data.frame(cbind(x$time,x$surv,x$lower,x$upper,x$std.err))
    res <- t(sapply(times,function(z,y){tps.surv(y,z)},y=y))
  } else
  {
    res <- NULL
    ld <- c(0,cumsum(table(x$strata)))
    for(i in 1:(length(ld)-1))
    {
      y <- as.data.frame(matrix(cbind(x$time,x$surv,x$lower,x$upper,x$std.err)[(1+ld[i]):ld[i+1],],ncol=5))
      res[[i]] <- t(sapply(times,function(z,y){tps.surv(y,z)},y=y))
    }
  }
  return(res)
}
clep <- function(x,y)
{
  a <- which(y<=x)
  b <- rep(0,length(y))
  b[a[length(a)]] <- 1
  return(b)
}
extps.surv <- function(obj, times)
{
  y <- obj
  names(y) <- c("time","surv","lower","upper","se")
  y$indic <- apply(sapply(times,clep,y=y[,1]),1,sum)
  y <- y[y$indic==1,]
  y <- cbind(times,y)
  names(y)[1:2] <- c("time","lastev.time")
  return(y[,1:6])
}

tps.surv <- function(obj, time)
{
  y <- rbind(c(0,1,NA,NA,NA),obj)
  names(y) <- c("time","surv","lower","upper","se")
  y$indic <- clep(time,y[,1])
  y <- y[y$indic==1,]
  y <- cbind(time,y)
  names(y)[1:2] <- c("time","lastev.time")
  return(y[,1:6])
}
