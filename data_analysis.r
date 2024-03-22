### empirical Bayes data analysis
## reference: http://varianceexplained.org/r/empirical_bayes_baseball/

data <- read.table('efron_baseball_1975.txt', header = T)
## type of data is data.frame

## column names are: Id FirstName   LastName Hits At.Bats RemainingAt.Bats 
## RemainingHits
at_bats1 = data$At.Bats  ## first 45 at-bats
at_bats2 = data$At.Bats+data$RemainingAt.Bats ## all at-bats
hits1 = data$Hits   ## hits in first 45 at-bats
hits2 = data$Hits + data$RemainingHits ## all hits

true_avg = hits2/at_bats2 ## treat overall batting average as real

## choose a jeffreys' prior as uninformative prior, beta(.5, .5)
indv_avg = (hits1+.5)/(at_bats1+1.) 
mle_avg = hits1/at_bats1

## number of objects
num_obj = length(hits1)
ss = 10000   ## sample size
indv_samp = matrix(0, ss, num_obj)  ## individual samples

### individual fits:
for (i in 1:num_obj){
  indv_samp[,i] <- rbeta(ss, (hits1[i]+.5), (at_bats1[i]-hits1[i]+.5))
}

## Empirical Bayes methods with MCEM algorithm
eb_hier <- function(samp_mat, iters=100, threshold=1.e-3){
  # individual prior: beta(.5, .5)
  # hierarchical prior: normal(gamma, tausq)
  #:params[in]:samp_mat, the sample matrix of theta_1,...,theta_n

  #:params[out]:next_mean,the next update of mean
  #:params[out]:next_tausq,the next update of variance.
  rows = nrow(samp_mat) # sample size 
  cols = ncol(samp_mat) # number of objects
  prior_mn = NULL
  prior_tausq = NULL
  prior_mn[1] = mean(samp_mat)  ## initialize mean
  prior_tausq[1] = mean((samp_mat-prior_mn[1])^2)  ## initialize tausq
  t=2   ## iteration
  wts0 = matrix(0, rows, cols) ## unnormalized weights
  dist = NULL  ## distance between consecutive iterates
  for (t in 2:iters){
    ## iterations 
    wts0 = dnorm(samp_mat, prior_mn[t-1], sqrt(prior_tausq[t-1]))/
               dbeta(samp_mat, .5, .5)
    ## normalize each column
    wts = sweep(wts0, 2, colSums(wts0), FUN="/")
    ## update mean and variance in prior distribution
    prior_mn[t] = sum(samp_mat*wts)/cols  # iteration of mean
    prior_tausq[t] = sum(((samp_mat-prior_mn[t])^2)*wts)/(cols)
    dist = c(dist, sqrt((prior_mn[t]-prior_mn[t-1])^2+
                          (prior_tausq[t]-prior_tausq[t-1])^2))
  }
  ## sample in the end
  fin_mn = tail(prior_mn, 1)
  fin_std = sqrt(tail(prior_tausq, 1))
  wts0 = dnorm(samp_mat, fin_mn, fin_std)/
    dbeta(samp_mat, .5, .5)
  fin_samp <- matrix(0, rows, cols)
  for (j in 1:cols){
    fin_samp[,j] = sample(samp_mat[,j], rows, replace = TRUE, prob=wts0[,j])
  }
  res <- list('prior_mean'=prior_mn, 'prior_var'=prior_tausq,
              'dist'=dist, 'eb_sample'=fin_samp)
  return(res)
  }
  
 


## fully bayesian data analysis
## reference: https://mc-stan.org/users/documentation/case-studies/pool-binary-trials.html
fb_hier <-function(size, burnin, indv.sam){
  # fit a -Normal hierarchical model with a Gibbs sampler
  # this Gibbs sampler takes advantage of individual fits of mu
  #  and the individual prior is p(mu_i) \propto beta(.5,.5)
  # then draw samples for latent parameters mu_{1},...,mu_{I}
  #
  #    data_{i} ~ binomial(mu_{i})
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
  ss <- dim(indv.sam)[2] ## sample size
  app.mu <- matrix(0, TSS, I) ##sample for mu in this method
  app.mu[1,] <- indv.sam[1,]
  app.gamma <- NULL
  app.tausq <- NULL
  app.tausq.scale <-sum((app.mu[1,]-mean(app.mu[1,]))^2)+1.e-8
  app.tausq[1] <- rigamma(1,(I-2)/2,app.tausq.scale/2)
  app.gamma[1] <- rnorm(1,mean(app.mu[1,]),sqrt(app.tausq[1]/I)) #optimal gamma
  app.acp <- matrix(0,TSS,I) ##acceptance matrix
  app.acp[1,] <- rep(1,I)
  for(aind in 2:TSS){
    ### proposals for different units should be randomly selected
    app.prop <- NULL ##proposals for different objects
    prop.ind <- sample.int(ss, I, replace = T) # sampled indices
    for (i in 1:I){app.prop[i] <- indv.sam[prop.ind[i], i]}
    app.prob <- dnorm(app.prop,app.gamma[aind-1],sqrt(app.tausq[aind-1]))/
      dnorm(app.mu[aind-1,],app.gamma[aind-1],sqrt(app.tausq[aind-1]))*
      dbeta(app.mu[aind-1,], .5,.5)/dbeta(app.prop, .5,.5)
    app.unif <- runif(I)
    app.acp[aind,] <- I(app.unif<=app.prob)
    app.delta <- app.prop - app.mu[aind-1,] ##difference between proposal and previous
    app.mu[aind,] <- app.mu[aind-1,]+app.delta*app.acp[aind,]
    app.tausq.scale <- sum((app.mu[aind,]-mean(app.mu[aind,]))^2)+1.e-8
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

### use both empirical bayes and fully bayesian for data analysis
at_bats1 = data$At.Bats  ## first 45 at-bats
at_bats2 = data$At.Bats+data$RemainingAt.Bats ## all at-bats
hits1 = data$Hits   ## hits in first 45 at-bats
hits2 = data$Hits + data$RemainingHits ## all hits

true_avg = hits2/at_bats2 ## treat overall batting average as real

## choose a jeffreys' prior as uninformative prior, beta(.5, .5)
indv_avg = (hits1+.5)/(at_bats1+1.) 
mle_avg = hits1/at_bats1

## number of objects
num_obj = length(hits1)
ss = 10000   ## sample size
indv_samp = matrix(0, ss, num_obj)  ## individual samples

### individual fits:
for (i in 1:num_obj){
  indv_samp[,i] <- rbeta(ss, (hits1[i]+.5), (at_bats1[i]-hits1[i]+.5))
}

### empirical bayes
eb=eb_hier(indv_samp, iters=100, threshold=1.e-3)
eb$prior_mean
eb$prior_var
fb=fb_hier(size=10000, burnin=1000, indv_samp)
fb$app.gamma
fb$app.tausq
dim(fb$app.mu)

fb_ba = colMeans(fb$app.mu)
eb_ba = colMeans(eb$eb_sample)

fb_mse = mean((fb_ba-true_avg)^2)
eb_mse = mean((eb_ba-true_avg)^2)
indv_mse = mean((indv_avg-true_avg)^2)

### load modules
library(devtools)
library(ggplot2)
library(dplyr)
library(LearnBayes)
library(Lahman)

## combine two columns of strings
data <- within(data,  Name <- paste(FirstName, LastName, sep=" "))

df <- rbind(data.frame(Type="Hierarchical FB", 
                       AVG=fb_ba,
                       Player=data$LastName),
            data.frame(Type="Single-case", 
                       AVG=indv_avg,
                       Player=data$LastName),
            data.frame(Type="Hierarchical EB", 
                       AVG=eb_ba,
                       Player=data$LastName))

p <- ggplot(data=df, aes(x=Type, y=AVG, 
                         color=Player)) +
  geom_line(aes(group=Player)) +
  #ggtitle("Shrinkage of Batting Averages")+
  labs(y= "Batting average", x = "Methods")+
  guides(color = guide_legend(override.aes = list(size=.5)))
print(p)
ggsave("batting_avg.png", plot=p, width = 5, height = 6, units = "in")

### make intervals plot
##by comparing the interval length
##we compute the posterior mean and quantiles
star.names <- data$LastName
fb.mat <- matrix(0,18,3) ##fb matrix
eb.mat <- matrix(0,18,3)
indv.mat <- matrix(0,18,3)
for(i in 1:18){
  fb.mat[i,1] <- mean(fb$app.mu[,i]) # posterior mean
  fb.mat[i,2] <- unname(quantile(fb$app.mu[,i],pnorm(-1))) #lower quantile
  fb.mat[i,3] <- unname(quantile(fb$app.mu[,i],pnorm(1))) #upper quantile
  eb.mat[i,1] <- mean(eb$eb_sample[,i]) # posterior mean
  eb.mat[i,2] <- unname(quantile(eb$eb_sample[,i],pnorm(-1))) #lower quantile
  eb.mat[i,3] <- unname(quantile(eb$eb_sample[,i],pnorm(1))) #upper quantile
  indv.mat[i,1] <- mean(indv_samp[,i]) # posterior mean
  indv.mat[i,2] <- unname(quantile(indv_samp[,i],pnorm(-1))) #lower quantile
  indv.mat[i,3] <- unname(quantile(indv_samp[,i],pnorm(1))) #upper quantile
}
xlim1 <- min(indv.mat)*0.9
xlim2 <- max(indv.mat)*1.1
### round 
indv.mat <- round(indv.mat,digits = 3)
fb.mat <- round(fb.mat, digits = 3)
eb.mat <- round(eb.mat,digits =3)
save(star.names,indv.mat, fb.mat, eb.mat, file="estimates.Rdata")

###plot the intervals
pdf(file="ba_intervals.pdf",width = 4,height = 6)
par(mfrow=c(1,1),mar=c(2.3,5.,.2,.1),cex=0.8,
    oma=rep(0.1,4), mgp=c(1.3,.3,.0),tck=0.02)
plot(0,0,xlim=c(xlim1,xlim2),ylim=c(0.8,18.3),yaxt="n",
     xlab="Batting Average",ylab="",main="")
axis(2,at=c(1:18),labels=star.names,las = 1)
for(i in 1:18){
  lines(c(indv.mat[i,2:3]),rep(i,2),lty=1,lwd=2)
  points(indv.mat[i,1],i,pch=16,cex=1.)
  points(true_avg[i],i,pch=15,cex=1.)
  lines(c(fb.mat[i,2:3]),rep(i+0.25,2),lty=2,lwd=2)
  points(fb.mat[i,1],i+.25,pch=2,cex=1.)
  lines(c(eb.mat[i,2:3]),rep(i-.25,2),lty=3,lwd=2)
  points(eb.mat[i,1],i-.25,pch=4,cex=1.)
}
legend(.35,18.5,c("Single-case","FB","EB",'Full season'),
       lty=c(1,2,3,NA),lwd=c(2,2,2,NA),pch = c(16,2,4,15))
dev.off()

eb.mat <- round(eb.mat, 2)
fb.mat <- round(fb.mat, 2)
indv.mat <- round(indv.mat, 2)
