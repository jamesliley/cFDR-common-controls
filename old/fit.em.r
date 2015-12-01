##' Fit a specific two Guassian mixture distribution
##'
##' Assumes Z ~ N(0,1) with prob pi0,
##' Z ~ N(0,1+sigma^2) with prob 1-pi0
##'
##' @title fit.em
##' @param Z numeric vector of observed data
##' @param s2 initial value for sigma^2
##' @param pi0 initial proportion of samples from f0
##' @param tol how small a change lhood prompts continued optimization
##' @param maxit maximum number of iterations
##' @return a list of two objects giving fitted values and history
##' @export
##' @author Chris Wallace
##' @examples
##' s2 <- 10
##' pi0 <- 0.8
##' n <- 10000
##' Z <- c(rnorm(round(n*pi0),0,1),rnorm(round(n*(1-pi0)),0,sqrt(1+s2)))
##' fit<-fit.em(Z)
##' fit$pars
##' fit$history
fit.em <- function(Z,pi0=0.9,s2=1,tol=1e-4,verbose=TRUE,maxit=1e4) {
    
  ## probabilities of group0, group1 membership
    p <- c(pi0,1-pi0)
    px <- matrix(p,length(Z),2,byrow=TRUE)

  ## parameter vector
  pars <- c(mu=c(0,0), sigma=sqrt(c(1,1+s2)))
  
  ## define lots of functions within this function to use this environment  
  parvec <- function(pars,index=TRUE) {
    mu <- pars[grep("^mu",names(pars))]
    sigma <- pars[grep("^sigma",names(pars))]
    names(mu) <- sub(".*\\.","",names(mu))
    names(sigma) <- sub(".*\\.","",names(sigma))
    return(list(mu=mu,sigma=sigma))
  }

  pars.fail <- function(pars) {
      if(parvec(pars)$sigma[2] <= 0)
          return(TRUE)
      return(FALSE)
  }
    
  ## likelihood for a single group
  lhood.single <- function(mu, sigma, pi) {
    exp(dnorm(Z,mean=mu,sd=sigma,log=TRUE) + log(pi))
  }
  
  
  ## likelihood function to be maximized
  lhood <- function(pars, px, sumlog=TRUE) {
     parv <- parvec(pars)
     if(pars.fail(pars))
         return(NA)
     ngroup <- ncol(px)
     e <- numeric(length(Z))
     for(i in 1:ngroup) {
      e <- e + lhood.single(mu=parv$mu[i], sigma=parv$sigma[i], pi=px[,i])
    }
     if(!any(is.na(e)) & any(e==0)) {
         wh <- which(e==0)
         e[wh] <- 1e-64
     }
     if(!any(is.na(e)) & any(is.infinite(e))) {
         wh <- which(is.infinite(e))
         e[wh] <- 1e64
     }
     if(sumlog) {
         return(-sum(log(e)))
     } else {
         return(-e)
     }
 }
    
    nit <- 1
    df <- 1
    value <- matrix(NA,maxit,3,dimnames=list(NULL,c("pi0","sigma2","lhood")))
    value[nit,] <- c(p[1], pars[4]^2 - 1, lhood(pars, px))
    while(df>tol & nit<maxit) {
        
        nit <- nit+1    
        parv <- parvec(pars)
    
        ## E step
        px.old <- px
        p <- colMeans(px)
        for(i in 1:2) {
            px[,i] <- lhood.single(parv$mu[i], parv$sigma[i], p[i])
        }
        px <- px/rowSums(px) ## normalise
        if(any(is.nan(px)))
            px[is.nan(px)] <- 0
        
        ## M step
        p <- colMeans(px)
        pars[4] <- sqrt(max(1, sum( px[,2] * Z^2 ) / sum(px[,2])))
        value[nit,] <- c(p[1], pars[4]^2 - 1, lhood(pars, px))
        df <- abs(value[nit,3] - value[nit-1,3])
    
    }
    return(list(pars=c(pi0=p[1], sigma2=pars[4]^2 - 1),
              history=value[1:nit,]))
}
