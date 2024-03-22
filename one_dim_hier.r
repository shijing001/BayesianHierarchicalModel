##use metropolis-hasting algorithm to correct the Gibbs sampler
app.gibbs <-function(size, burnin, indv.sam){
  # fit a -Normal hierarchical model with a Gibbs sampler
  # this Gibbs sampler takes advantage of individual fits of mu
  #  and the individual prior is p(mu_i) \propto 1
  # then draw samples for latent parameters mu_{1},...,mu_{I}
  #
  #    data_{i} ~ Normal(mu_{i}, sigma_{i}^2)
  #    mu_{i} ~ Normal(gamma, tau^2)
  #    p(gamma | tau) \propto 1
  #    p(tau) \propto 1
  #
  #    variance is parametrized in tau^2, not in tau
  # Args:
  #  size: sample size from the posterior distribution
  #  burnin: the number of iterations in burn-in period
  #  indv.sam, a 2-dimensional array, with each column represents the samples of mu
  # 
  # Returns:
  #   the samples for latent parameters and point estimates for hyper-parameters
  #   in the normal-normal hierarchical model
  library(pscl)
  start.time <- Sys.time()
  I <-  dim(indv.sam)[2] ##number of columns, i.e., groups
  TSS <- size + burnin
  indv.sam <- matrix(0, TSS, I)
  app.mu <- matrix(0, TSS, I) ##sample for mu in this method
  app.mu[1,] <- indv.sam[1,]
  app.gamma <- NULL
  app.tausq <- NULL
  app.tausq.scale <-sum((app.mu[1,]-mean(app.mu[1,]))^2)
  app.tausq[1] <- rigamma(1,(I-2)/2,app.tausq.scale/2)
  app.gamma[1] <- rnorm(1,mean(app.mu[1,]),sqrt(app.tausq[1]/I)) #optimal gamma
  app.acp <- matrix(0,TSS,I) ##acceptance matrix
  app.acp[1,] <- rep(1,I)
  for(aind in 2:TSS){
    ### proposals for different units should be randomly selected
    app.prop <- NULL ##proposals for different objects
    prop.ind <- sample.int(TSS, I, replace = T) # sampled indices
    for (i in 1:I){app.prop[i] <- indv.sam[prop.ind[i], i]}
    app.prob <- dnorm(app.prop,app.gamma[aind-1],sqrt(app.tausq[aind-1]))/
      dnorm(app.mu[aind-1,],app.gamma[aind-1],sqrt(app.tausq[aind-1]))
    app.unif <- runif(I)
    app.acp[aind,] <- I(app.unif<=app.prob)
    app.delta <- app.prop - app.mu[aind-1,] ##difference between proposal and previous
    app.mu[aind,] <- app.mu[aind-1,]+app.delta*app.acp[aind,]
    app.tausq.scale <- sum((app.mu[aind,]-mean(app.mu[aind,]))^2)
    app.tausq[aind] <- rigamma(1,(I-2)/2,app.tausq.scale/2)
    app.gamma[aind] <- rnorm(1,mean(app.mu[aind,]),sqrt(app.tausq[aind]/I))
  }
  acprat <- colMeans(app.acp)
  begn <- burnin + 1 ## remove half -- the burn-in period
  app.mu.est <- app.mu[begn:TSS,] ##estimator of mu
  app.gamma.est <- app.gamma[begn:TSS]
  app.tausq.est <- app.tausq[begn:TSS]
  time.len <- Sys.time() - start.time
  time.num = as.numeric(gsub('.*-([0-9]+).*','\\1',time.len))
  lst <- list("app.mu"=app.mu.est,"time"=time.num,
              "acprat"=acprat,"app.gamma"=app.gamma.est,"app.tausq"=app.tausq.est)
  return(lst)
}

### multivariate Gaussian hierarchical 
MHapproxFBifmr <- function(IndvRes, size, burnin){
  # fit a Bayesian hierarchical model given individual fits
  #
  # Args:
  #  size: sample size from the posterior distribution
  #  burnin: the number of iterations in burn-in period
  #  IndvRes, a 3-dimensional array, (num_objects, sample size, dimensions)
  # 
  # Returns:
  #   the samples for latent parameters and point estimates for hyper-parameters
  
  library(MCMCpack)
  library(tmvtnorm) 
  I <-  dim(IndvRes)[1] ##number of groups
  TSS <- dim(IndvRes)[2] ##sample size in the individual fitting
  dims <- dim(IndvRes)[3] ## dimensions
  app.mu <- array(0,dim=c(TSS,I,dims)) ##sample for mu in this method
  indv.mn <- apply(IndvRes,c(1,3),mean) ##means from individual fitting
  app.mu[1,,] <- indv.mn ##initialization
  app.gamma <- matrix(0,TSS,2) ##population meam samples
  app.Sigma <- array(0,dim=c(TSS,2,2)) ##population variance samples
  app.gamma.mn <- matrix(colMeans(app.mu[1,,]),2,1) ##matrix
  app.Sigma.scl <-t(app.mu[1,,])%*%app.mu[1,,]-
    I*app.gamma.mn%*%t(app.gamma.mn)+10^(-9)*diag(2) ##scale matrix
  app.Sigma[1,,] <- riwish(I,app.Sigma.scl)
  app.gamma[1,] <- rmvnorm(1,c(app.gamma.mn),app.Sigma[1,,]/I)
  acp.idx <-matrix(0,TSS,I)  ##acceptance index matrix
  for(aind in 2:TSS){
    app.prop <-IndvRes[,aind,] ##use individual result as proposal
    app.prob <- dmvnorm(app.prop,app.gamma[aind-1,],app.Sigma[aind-1,,])/
      dmvnorm(app.mu[aind-1,,],app.gamma[aind-1,],app.Sigma[aind-1,,])
    app.unif <- runif(I)
    acp.idx[aind,] <-I(app.unif<=app.prob)
    app.delta <-app.prop - app.mu[aind-1,,] ##difference between proposal and previous
    app.mu[aind,,] <- app.mu[aind-1,,]+diag(acp.idx[aind,])%*%app.delta
    app.gamma.mn1 <- matrix(colMeans(app.mu[aind,,]),2,1)
    app.Sigma.scl1 <-t(app.mu[aind,,])%*%app.mu[aind,,]-
      I*app.gamma.mn1%*%t(app.gamma.mn1)+10^(-9)*diag(2) ##scale matrix
    app.Sigma[aind,,] <- riwish(I,app.Sigma.scl1)
    app.gamma[aind,] <- rmvnorm(1,c(app.gamma.mn1),app.Sigma[aind,,]/I)
  }
  acprat <- sum(acp.idx)/(TSS*I)
  lst <- list("app.mu"=app.mu,
              "acprat"=acprat,"acp.ind"=acp.idx,"app.gamma"=app.gamma,"app.Sigma"=app.Sigma)
  return(lst)
}

