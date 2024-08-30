#################################################################################################
#################################################################################################

###  Supplementary code for 'Uncertainty Quantification for Agent Based Models: A Tutorial'   ###

#################################################################################################
#################################################################################################

# Required R packages:
library(lhs)
library(ggplot2)
library(latex2exp)
library(invgamma)
library(mvtnorm)
library(rdist)
library(DiceKriging)
library(hetGP) 

#################################################################################################
# Section 2.1: Wolf and Sheep Predator Prey ABM
#################################################################################################


# Latin hypercube for initial design
X <- maximinLHS(30, 2)
Xx <- round(20*X, 2) # Scale to [0,20]

# Run initial data through ABM to generate training data

# Load in data
Data <- read.csv('DataSW.csv')

# Inputs
X <- Data[,2:3]/20
# Outputs
YData <- Data[,4:dim(Data)[2]]
YWolf <- YData[,(seq(2,dim(YData)[2],2))]
YSheep <- YData[,(seq(1,dim(YData)[2],2))]

# Training data -> 0/1 wolf extinction
WData <- cbind(X, YWolf)
XWolf1 <- WData[is.na(WData[,3])==FALSE,(1:2)]
YWolf1 <- WData[is.na(WData[,3])==FALSE,(3:dim(WData)[2])]
XWolf0 <- WData[is.na(WData[,3])==TRUE,(1:2)]

# Plot initial design
colours <- c('Yes'='firebrick3', 'No'='dodgerblue3')
ggplot() +
  geom_point(aes(x=XWolf1[,1], y=XWolf1[,2], col='Yes'), size=3) +
  geom_point(aes(x=XWolf0[,1], y=XWolf0[,2], col='No'), size=3) +
  scale_color_manual(values=colours) +
  xlab(TeX('Sheep Reproduction')) + ylab(TeX('Wolf Reproduction')) + labs(colour=TeX('Wolves Extinct'))




#################################################################################################
# Section 2.2: Classification
#################################################################################################


m <- 2 # Dimension

N1 <- dim(XWolf1)[1] # Number of points in each region
N2 <- dim(XWolf0)[1]

Ddata <- as.matrix(rbind(XWolf1,XWolf0)) # Full data

# Run MCMC to calculate latent values

itno <- 50000 # no. of iterations

count <- 0 
MCMCoutput <- metropnd(iterations=itno, eps=0.8, startvalue=c(0,-1,1,1,1,1), D=Ddata, n1=N1, n2=N2) 

# Acceptance rate
MCMCoutput[[2]]/itno

# Trace plots
par(mfrow=c(2,2))
for ( k in 1:(m + 2)) {
  plot(MCMCoutput[[1]][,k],type='l')
}

# Burn-in
burnno <- 5000
out_new <- MCMCoutput[[1]][(burnno:itno),]
par(mfrow=c(2,2))
for ( k in 1:(m + 2)) {
  plot(out_new[,k],type='l')
}


# EITHER sample directly from posterior distribution
# OR sample conditionally

###
# Option 1: Sample directly (expensive)
###

SampOut <- c()

for (ij in 1:dim(out_new)[1]) {
  #print(ij)
  SO <- Sampfun(out_new[ij,], Nn=N1+N2, Dd=Ddata, N1=N1, N2=N2)
  if (is.na(SO[1]) == FALSE) {
    SampOut <- rbind(SampOut, c(out_new[ij,],SO))
  }
  
}

###
# Option 2: Sample conditionally
###

Ddata01 <- cbind(Ddata, c(rep(0,N1),rep(1,N2)))
outGPs <- apply(out_new, 1, condsamp, Dnew=Ddata01)

totalpoints <- outGPs[[1]][[1]]

for (k in 2:dim(out_new)[1]) {
  totalpoints <- cbind(totalpoints,outGPs[[k]][[1]])
}

dim(totalpoints)

newout <- t(totalpoints)  # All latent points
SampOut <- cbind(newout,out_new) # Points and parameters

for (ij in 1:n) {
  newout <- newout[(is.na(newout[,n]) == FALSE),]
  SampOut <- SampOut[(is.na(SampOut[,n]) == FALSE),]
}



##########

SampPar <- SampOut[,(1:6)]
SampVal <- SampOut[,(7:dim(SampOut)[2])]

post2 <- apply(SampVal,2,sort)

lp <- dim(post2)[1]

# Predict classification across input space

medpoints <- post2[round(lp/2),]
minp <- post2[round((lp*0.05)/2),]
maxp <- post2[round(lp-((lp*0.05)/2)),]

# Training data
fD <- medpoints
D <- Ddata

# Set number of data points and parameters
n <- dim(D)[1]
p <- dim(D)[2]
q <- p+1

# Training data to fit the GP
xx <- seq(0,1,0.01) 
D2 <- expand.grid(xx,xx)
Data2 <- data.frame(
  Sheep=D2[,1],
  Wolves=D2[,2]
)

# Fit GP
library(DiceKriging)
outmodel <- km(formula=~., design=as.matrix(Ddata), response=fD, covtype = 'matern3_2')
predmodel <- predict(outmodel, as.matrix(Data2), type='UK')

PM <- predmodel$mean
PM1 <- PM
PM1[PM < 0] <- 'Yes'
PM1[PM > 0] <- 'No'

# Plot
ggplot() +
  geom_tile(aes(x=D2[,1], y=D2[,2], fill=PM1), alpha=0.6) +
  scale_fill_manual(values=c("dodgerblue3","firebrick3")) +
  geom_point(aes(x=Ddata[,1], y=Ddata[,2]), size=2) +
  xlab(TeX('Sheep Reproduction')) + ylab(TeX('Wolf Reproduction')) + labs(fill=TeX('Wolves Extinct'))


# Find GP samples to find classification uncertainty
PredSamp <- t(apply(SampVal, 1, Krig_samples))
PredSamp2 <- PredSamp

PredSamp2[PredSamp < 0] <- 1
PredSamp2[PredSamp > 0] <- 0

SampPer <- apply(PredSamp2, 2, sum)/dim(PredSamp2)[1]

ggplot() +
  geom_raster(aes(x=D2[,1], y=D2[,2], fill=SampPer), alpha=0.8) +
  scale_fill_gradientn(colours=c('dodgerblue3','firebrick3')) +
  geom_point(aes(x=Ddata[,1], y=Ddata[,2]), size=2) +
  xlab(TeX('Sheep Reproduction')) + ylab(TeX('Wolf Reproduction')) + labs(fill=TeX('Prob of Extinct'))


# Store classified inputs for future code
Data_W1 <- Data2[(PM1 == 'Yes'),]
Data_W0 <- Data2[(PM1 == 'No'),]



#################################################################################################
# Section 2.3: Stochastic Gaussian Process
#################################################################################################

# Load in data again
Data <- read.csv('DataSW.csv')

Ddata <- Data[,2:3]/20
YData <- Data[,4:dim(Data)[2]]
YWolf <- YData[,(seq(2,dim(YData)[2],2))]
YSheep <- YData[,(seq(1,dim(YData)[2],2))]

WData <- cbind(Ddata, YWolf)

XWolf1 <- WData[is.na(WData[,3])==FALSE,(1:2)]
YWolf1 <- WData[is.na(WData[,3])==FALSE,(3:dim(WData)[2])]

XWolf0 <- WData[is.na(WData[,3])==TRUE,(1:2)]

colours <- c('Yes'='firebrick3', 'No'='dodgerblue3')
ggplot() +
  geom_point(aes(x=XWolf1[,1], y=XWolf1[,2], col='Yes'), size=3) +
  geom_point(aes(x=XWolf0[,1], y=XWolf0[,2], col='No'), size=3) +
  scale_color_manual(values=colours) +
  xlab(TeX('Sheep Reproduction')) + ylab(TeX('Wolf Reproduction')) + labs(colour=TeX('Wolves Extinct'))


# Put data into format for HetGP
Xhet <- as.matrix(do.call("rbind", rep(list(XWolf1), 10)))
Yhet <- as.vector(as.matrix(YWolf1))

# Fit HetGP
hetgpmod <- mleHetGP(X=Xhet, Z=Yhet)
hetpred <- predict(hetgpmod, as.matrix(Data_W1), noise.var=TRUE)
predmu <- hetpred$mean
predvar <- hetpred$sd2 + hetpred$nugs

# Plot predicted mean response and variance
a1 <- ggplot() +
  geom_raster(aes(x=Data_W1[,1], y=Data_W1[,2], fill=predmu)) +
  scale_fill_gradientn(colours=c('white','firebrick3')) +
  geom_tile(aes(x=Data_W0[,1], y=Data_W0[,2]), fill='dodgerblue3', alpha=0.6) +
  geom_point(aes(x=Ddata[,1], y=Ddata[,2]), size=2) +
  xlab(TeX('Sheep Reproduction')) + ylab(TeX('Wolf Reproduction')) + labs(fill=TeX('Time to Extinct'))

a2 <- ggplot() +
  geom_raster(aes(x=Data_W1[,1], y=Data_W1[,2], fill=predvar)) +
  scale_fill_gradientn(colours=c('white','firebrick3')) +
  geom_tile(aes(x=Data_W0[,1], y=Data_W0[,2]), fill='dodgerblue3', alpha=0.6) +
  geom_point(aes(x=Ddata[,1], y=Ddata[,2]), size=2) +
  xlab(TeX('Sheep Reproduction')) + ylab(TeX('Wolf Reproduction')) + labs(fill=TeX('Variance'))

grid.arrange(a1,a2, ncol=2)



#################################################################################################
# Section 3: Sequential Design
#################################################################################################


# Design using hetGP and IMSPE

# Put data into format for HetGP
Xhet <- as.matrix(do.call("rbind", rep(list(XWolf1), 10)))
Yhet <- as.vector(as.matrix(YWolf1))

# Fit initial HetGP
hetgpmod <- mleHetGP(X=Xhet, Z=(Yhet))

pred <- predict(hetgpmod, as.matrix(Data_W1), noise.var=TRUE)
predmu <- pred$mean
predvar <- pred$sd2 + pred$nugs

# Candidate points for selection (so that no design points are 
# chosen in region where wolves don't go to extinction)
Xcand1 <- as.matrix(rbind(XWolf1, as.matrix(Data_W1)))

# Select new point
SDhetGP <- IMSPE_optim(hetgpmod, Xcand=Xcand1)
SDhetGP$par
SDhetGP$par*20

# Plot with new point
ggplot() +
  geom_raster(aes(x=Data_W1[,1], y=Data_W1[,2], fill=predmu)) +
  scale_fill_gradientn(colours=c('white','firebrick3')) +
  geom_tile(aes(x=Data_W0[,1], y=Data_W0[,2]), fill='dodgerblue3', alpha=0.6) +
  geom_point(aes(x=Ddata[,1], y=Ddata[,2]), size=2) +
  geom_point(aes(x=SDhetGP$par[1], y=SDhetGP$par[2]), col='forestgreen', size=5, shape=4, stroke=4) +
  xlab(TeX('Sheep Reproduction')) + ylab(TeX('Wolf Reproduction')) + labs(fill=TeX('Time to Extinct'))





#################################################################################################
# Section 4: History Matching
#################################################################################################


# Observation - chosen at [2.5,2.5]
Obs <- c(450, 368, 352, 364, 370, 412, 342, 498, 458, 368)


# Load in data again
Data <- read.csv('DataSW.csv')
X <- Data[,2:3]/20
YData <- Data[,4:dim(Data)[2]]
YWolf <- YData[,(seq(2,dim(YData)[2],2))]
YSheep <- YData[,(seq(1,dim(YData)[2],2))]

WData <- cbind(X, YWolf)

XWolf1 <- WData[is.na(WData[,3])==FALSE,(1:2)]
YWolf1 <- WData[is.na(WData[,3])==FALSE,(3:dim(WData)[2])]

XWolf0 <- WData[is.na(WData[,3])==TRUE,(1:2)]

colours <- c('Yes'='firebrick3', 'No'='dodgerblue3')
ggplot() +
  geom_point(aes(x=XWolf1[,1], y=XWolf1[,2], col='Yes'), size=3) +
  geom_point(aes(x=XWolf0[,1], y=XWolf0[,2], col='No'), size=3) +
  scale_color_manual(values=colours) +
  xlab(TeX('Sheep Reproduction')) + ylab(TeX('Wolf Reproduction')) + labs(colour=TeX('Wolves Extinct'))


# Take the sample mean of the data to use as an approximation to the mean response
YW_Mean <- apply(YWolf1, 1, mean)
# Fit a GP to the sample means (log taken for computational stability)
ModelL1 <- km(design = data.frame(XWolf1), response = log(YW_Mean), covtype = "matern3_2",
              control = list(trace = FALSE, BFGSburnin=2, max.generations = 20))
predL1 <-  predict(ModelL1, newdata=(Data_W1), type="UK", checkNames=F, 
                   light.return = F)
# Plot
ggplot() +
  geom_raster(aes(x=Data_W1[,1], y=Data_W1[,2], fill=exp(predL1$mean))) +
  scale_fill_gradientn(colours=c('white','firebrick3')) +
  geom_tile(aes(x=Data_W0[,1], y=Data_W0[,2]), fill='dodgerblue3', alpha=0.6) +
  geom_point(aes(x=XWolf1[,1], y=XWolf1[,2]), size=2) +
  geom_point(aes(x=0.125, y=0.125), size=5, col='forestgreen', shape=4, stroke=4) +
  xlab(TeX('Sheep Reproduction')) + ylab(TeX('Wolf Reproduction')) + labs(fill=TeX('Time to Extinct'))


# Mean of the observation
MObs <- mean(Obs)
MObs

# Implausibility function
ImpMeasure <- function(MPred, VPred, tDisc, tObsErr, tObs) {
  IMPM <- abs(tObs - MPred) / sqrt(tDisc + tObsErr + VPred)
  return(IMPM)
}

# Apply implausibility metric
# True discrepancy is unknown, estimated to be 0.02
# Zero observation error
IMPs <- ImpMeasure(MPred=predL1$mean, VPred=((predL1$sd)^2), tDisc=0.02, tObsErr=0, tObs=log(MObs))
IMPs

# Plot NROY space
ROPsW1 <- Data_W1[(IMPs > 3),]
NROPsW1 <- Data_W1[(IMPs < 3),]

ggplot() +
  geom_tile(aes(x=ROPsW1[,1], y=ROPsW1[,2]), fill='grey50', alpha=0.6) +
  geom_tile(aes(x=NROPsW1[,1], y=NROPsW1[,2]), fill='firebrick3', alpha=0.6) +
  geom_tile(aes(x=Data_W0[,1], y=Data_W0[,2]), fill='dodgerblue3', alpha=0.6) +
  geom_point(aes(x=XWolf1[,1], y=XWolf1[,2]), size=2) +
  geom_point(aes(x=XWolf0[,1], y=XWolf0[,2]), size=2) +
  geom_point(aes(x=0.125, y=0.125), size=5, col='forestgreen', shape=4, stroke=4) +
  xlab(TeX('Sheep Reproduction')) + ylab(TeX('Wolf Reproduction')) + labs(fill=TeX('Time to Extinct'))











#################################################################################################
# FUNCTIONS
#################################################################################################


# metropolis hastings MCMC
metropnd <- function(iterations=100, eps=0.05, startvalue, D, n1, n2) {
  m <- dim(D)[2]
  n <- n1+n2
  chain <- matrix(0,iterations+1,(2*m+2))
  info <- matrix(0,iterations+1,3)
  xS <- startvalue
  newinfo <- c(0,0,0)
  
  beta <- startvalue[1:(m+1)]
  sigma2 <- startvalue[m+2]
  ndelta <- startvalue[(m+3):(2*m+2)]
  
  # covariance matrix
  A <- K12(X1=D, X2=D, rho=ndelta, alpha=sigma2)
  
  H <- cbind(matrix(1,n,1),D)
  mm2 <- as.vector(H %*% beta)
  
  # likelihood, prior, posterior
  oldll <- log(pmvnorm(lower=c(rep(-Inf,n1),rep(0,n2)),upper=c(rep(0,n1),rep(Inf,n2)),mean=mm2,sigma=A))[1]
  prior <-  log(dinvgamma(sigma2,10,10)) + sum(log(dinvgamma((ndelta),10,40)))  + log(dnorm(beta[1],0,5)) + log(dnorm(beta[2],0,5)) + log(dnorm(beta[3],0,5))
  oldd <- oldll+prior
  
  chain[1,] <- startvalue
  info[1,] <- c(oldd,oldll,prior)
  
  prop <- rep(0,(2*m + 2))

  for (i in 1:iterations) {
    
    prop[1] <- xS[1] + rnorm(1,0,eps) # Add on small quantity
    prop[2] <- xS[2] + rnorm(1,0,eps)
    prop[3] <- xS[3] + rnorm(1,0,eps)
    prop[4] <- xS[4] + rnorm(1,0,eps)
    prop[5] <- xS[5] + rnorm(1,0,eps)
    prop[6] <- xS[6] + rnorm(1,0,eps)
    
    # acceptance function
    a <- MCMC_acceptancend(par=prop, olddist=oldd, dD=D, nN1=n1, nN2=n2)
    
    if (a[1]==1) { 
      count = count + 1
      xS <- prop
      newinfo <- a[2:4]
      oldd <- a[2]  
    }               
    chain[(i+1),] <- xS
    info[(i+1),] <- newinfo
  }
  return(list(chain,count,info))
}




# Acceptance function
MCMC_acceptancend <- function(par,olddist=oldd, dD, nN1, nN2){
  
  mM <- dim(dD)[2]
  nN <- nN1+nN2
  if (all(par[(mM+2):(2*mM+2)] > 0) == FALSE) return(c(0,olddist)) # check positivity

  # parameters
  beta <- par[1:(mM+1)]
  sigma2 <- par[mM+2]
  ndelta <- (par[(mM+3):(2*mM+2)])

  hH <- cbind(matrix(1,nN,1),dD)
  mm <- as.vector(hH %*% beta)
 
  # covariance matrix
  A <- K12(X1=dD, X2=dD, rho=ndelta, alpha=sigma2)
  
  # likelihood, prior, posterior
  loglik <- log(pmvnorm(lower=c(rep(-Inf,nN1),rep(0,nN2)),upper=c(rep(0,nN1),rep(Inf,nN2)),mean=mm,sigma=A))[1]
  prior <-  log(dinvgamma(sigma2,10,10)) + sum(log(dinvgamma((ndelta),10,40)))  + log(dnorm(beta[1],0,5)) + log(dnorm(beta[2],0,5)) + log(dnorm(beta[3],0,5))
  distr <- loglik+prior
  
  ratio <- exp(distr-olddist)
  rstart <- min(c(1,ratio))
  
  # comparison with the observed data
  if(runif(1,0,1) < rstart) return(c(1,distr,loglik,prior)) else return(c(0,distr,loglik,prior))
  
}


# covariance function
K12 <-  function(X1, X2, rho, alpha) { 
  #Calculates squared exponential covariance matrix for given inputs
  D = as.matrix(cdist(scale(X1,center=FALSE,scale=rho),scale(X2,center=FALSE,scale=rho)))
  return((alpha^2) * exp(-(D^2)))
}


# Sample GP draws from posterior
Sampfun <- function(Par, Nn, Dd, N1, N2) {
  
  beta <- Par[1]
  sigma2 <- Par[2]
  delta <- Par[3:4]
  
  H <- rep(1, Nn)
  mm <- rep(beta, Nn)
  
  # covariance matrix
  A <- K12(X1=Dd, X2=Dd, rho=delta, alpha=sigma2)
  
  mvtSamples <- rmvnorm(1, mean=mm, sigma=A)
  if (all(mvtSamples[1:N1] < 0) & all(mvtSamples[(N1+1):(N1+N2)] >= 0)) { 
    SampC <- mvtSamples 
  } else {
    SampC <- NA
  }
  return(SampC)
}


# Conditional sampling
condsamp <- function(parameters,Dnew=Dnew_con,no_iter=100) {
  
  beta <- parameters[1:3]
  sigma2 <- parameters[4]
  ndelta <- parameters[5:6]
  
  Data <- Dnew[,(1:m)]
  posneg <- Dnew[,(m+1)]
  n <- dim(Data)[1]
  
  H <- cbind(rep(1,n),Data)
  mx1 <- as.vector(as.matrix(H) %*% beta)
  
  A1 <- K12(X1=Data, X2=Data, rho=ndelta, alpha=sigma2) + diag(1e-5,n)
  
  y_out2 <- matrix(0,no_iter+1,n)
  y_out2[1,] <- posneg
  y_out1 <- posneg
  for ( k in 1:no_iter ) {
    for ( j in 1:n ) {
      if (all(is.na(y_out1) == FALSE)) {
        sEV <- samp_E_var(num=j,A=A1,mx=mx1,y_out=y_out1)
        yo <- samp_accp(num=j,Mean=sEV[[1]],Var=sEV[[2]],pn=posneg)
        y_out1[j] <- yo
      } else {
        y_out1[j] <- NA
      }
    }
    y_out2[(k+1),] <- y_out1
  }
  return(list(y_out1,y_out2))
}

samp_E_var <- function(num,A=A1,mx=mx1,y_out=y_out1) {
  
  covbit2 <- A[num,-num]
  varbit <- chol(A[-num,-num])
  M1 <- backsolve(varbit, covbit2, transpose=TRUE)
  M2 <- backsolve(varbit, y_out[-num] - mx[-num], transpose=TRUE)
  
  Mean <- mx[num] + crossprod(M1,M2)
  Var <- A[num,num] - crossprod(M1)

  return(list(Mean,Var))
}

# Function returning accepted y value
samp_accp <- function(num,Mean,Var,pn) {
  
  yi <- rnorm(1,Mean,sqrt(abs(Var)))
  countit <- 0 
  while( (sign(yi)!=sign(pn[num])) & (is.na(yi) != TRUE) ) {
    if (countit < 5000) {
      yi <- rnorm(1,Mean,sqrt(abs(Var)))
      countit <- countit + 1
    } else {
      yi <- NA
    }
  }
  return(yi)
}


# GP samples
Krig_samples <- function(fD, D1=Ddata, Data=Data2) {
  outmodel1 <- km(formula=~., design=D1, response=fD, covtype='matern3_2')
  predmodel1 <- predict(outmodel1, as.matrix(Data), type='UK')
  
  mstar <- predmodel1$mean
  return(mstar)
}
